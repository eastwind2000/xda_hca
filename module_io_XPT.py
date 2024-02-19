
import warnings
import datetime
from pathlib import Path

from typing import Union

import numpy as np

import xarray as xr

import bz2

from cinrad.constants import deg2rad, con
from cinrad.projection import get_coordinate, height
from cinrad.error import RadarDecodeError

from cinrad.io.base import RadarBase, prepare_file

from cinrad._typing import Number_T

import pandas as pd

from matplotlib import pyplot as plt

from module_io_defineXPT import *

import cmaps

# import polars as pl

# ====================================== AXPT data class ==============================================

"""

revised by chentao@cma.gov.cn 2024-02-03 21:51:49 for Beijing X-radar, VCP21, 11 elevations for dual-PRF

"""

# =======================================================================================================

class XPT_DATA(RadarBase):

    def __init__(self, file):

        f = prepare_file(file)              # Auto-adapted for .bz2

        filename = Path(file).name if isinstance(file, str) else ""
        
        try:
            self.code = self.name = filename.split("_")[3]
        except IndexError:
            self.code = self.name = "None"        

        self._d = []

        # =================== 2024-01-23 15:21:18 by chentao@cma.gov.cn  =====================================
           
        bytebufr = f.read()

        f.close()

        pos00 =  SIZE_AXPT_HEADER00

        data_header    = np.frombuffer(bytebufr[0:pos00], dtype=np.dtype(AXPT_HEADER00))

        scan_starttime = data_header["taskconfig_header"]["scan_starttime"][0]

        taskname       = data_header["taskconfig_header"]["taskname"][0]
        
        cut_number     = data_header["taskconfig_header"]["cut_number"][0]

        self._d.append(data_header)

        self.scan_type = taskname

        # self.available_product(0) = 
        # breakpoint()
        
        # ============================= cut_config_header part ==========================================
        
        total_cutheader = [ ("tch", cut_number*np.dtype(cut_header) ) ]   

        SIZE_TCH        = np.dtype(total_cutheader).itemsize
                        
        data_tch        = np.frombuffer( bytebufr[pos00:pos00+SIZE_TCH], dtype=np.dtype(total_cutheader) ) 

        self._d.append(data_tch)

        self.scan_config = data_tch
        
        # ================================ radia_header part ============================================

        pos_radial00 = pos00 + SIZE_TCH          #! byte position

        num_radial =  0                          #! radail_num counter

        datalist = []                            #! create dump list

        rpos_start = pos_radial00

        while( rpos_start < len(bytebufr) ):      #! to be revised
            
            data_radial_header = np.frombuffer( bytebufr[rpos_start:rpos_start+SIZE_RADIALHEADER], dtype=np.dtype(radial_header) )

            radial_header_info = compose_radial_header( data_radial_header)

            radial_datalength  = data_radial_header['length_data'][0]       #! radial datapart length
                
            rec_pos     = rpos_start + SIZE_RADIALHEADER                      #! radial data position
            data_radial = bytebufr[rec_pos : rec_pos + radial_datalength]

            rabufr = decode_radial(radial_header_info, data_radial)
            datalist.extend(rabufr)                                         #! code trick for list extend method, each moment data piece as one datarecord

            self._d.append(data_radial_header)
            self._d.append(data_radial)

            num_radial = num_radial + 1                                     #! radial counter
 
            rpos_start = rpos_start + SIZE_RADIALHEADER + radial_datalength #! control loop by rpos_start

        # ===================================================================================================================================
        
        pdbufr  =  pd.DataFrame(datalist)   
        # pdbufr = pl.DataFrame(datalist)    
      
        pdbufr.columns = [ "radial_stat",  "radial_time", "elev_number", "elev", "radial_number", "azimuth", 
                           "moment_name", "scale", "offset",  "num_gates", "moment_data"   ]
        
        # =======================================================================================================================

        self.stationlon  =  data_header['siteconfig_header']["longitude"][0]
        self.stationlat  =  data_header['siteconfig_header']["latitude"][0]
        self.radarheight =  data_header['siteconfig_header']["ground_height"][0]
        
        # self.scantime = datetime.datetime.utcfromtimestamp( data_header["taskconfig_header"]["scan_time"][0]    )
        
        log_reso       = data_tch["tch"][0, 0]["log_resolution"]
        dop_reso       = data_tch["tch"][0, 0]["doppler_resolution"]
        start_range    = data_tch["tch"][0, 0]["start_range"]
        nyquist_speed  = data_tch["tch"][0, :]["nyquist_speed"]    #!vital parameter for velocity

        self.log_reso  = log_reso
        self.dop_reso  = dop_reso

        #============================== revised by Chent@NMC 2024-01-25 11:14:34=======================

        el_num = cut_number       
        
        # az_num = ref.shape[1]                                  #!  unused var to be revsised
        
        elevation = data_tch["tch"]["elevation"][0]

        # ============================================================================================= kernel parameters
        
        # self.scantime        = scantime
        self.reso            = log_reso
        self.first_gate_dist = start_range
        self.el              = elevation        
        self.nyquist_speed   = nyquist_speed

        self.data = dict()
        self.aux = dict()

        # ============================================================================================== 2024-01-25 09:19:42
        # az = data["data"]["azimuth"].astype(int).reshape(az_num, el_num) * 360 / 65535
        
        radialtime_list = []

        for el_idx in np.arange(el_num):

            self.data[el_idx] = dict()

            tref = compose_var_sweep(pdbufr, "DBT", el_idx)
            ref  = compose_var_sweep(pdbufr, "DBZ", el_idx)
            vel  = compose_var_sweep(pdbufr, "vel", el_idx)
            sw   = compose_var_sweep(pdbufr, "sw",  el_idx)
            zdr  = compose_var_sweep(pdbufr, "zdr", el_idx)
            phi  = compose_var_sweep(pdbufr, "phi", el_idx)
            kdp  = compose_var_sweep(pdbufr, "kdp", el_idx)
            rho  = compose_var_sweep(pdbufr, "rho", el_idx)

            if ( len(tref) > 1):
                self.data[el_idx]["TREF"] =  tref

            if ( len(ref) > 1):
                self.data[el_idx]["REF"]  =  ref

            if ( len(vel) > 1):
                self.data[el_idx]["VEL"]  =  vel
            
            if ( len(sw) > 1):
                self.data[el_idx]["SW"]   =  sw
            
            if ( len(zdr) > 1):
                self.data[el_idx]["ZDR"]  =   zdr

            if ( len(phi) > 1):
                self.data[el_idx]["PHI"]  =   phi
      
            if ( len(kdp) > 1):
                self.data[el_idx]["KDP"]  =   kdp
            
            if ( len(rho) > 1):
                self.data[el_idx]["RHO"]  =   rho

            self.aux[el_idx] = dict()

            if( el_idx in [1, 3] ):                                #! code trick for VCP21   
                cond_var = pdbufr["moment_name"]=="vel"  
            else:            
                cond_var = pdbufr["moment_name"]=="DBZ" 
            
            cond_elev = pdbufr["elev_number"]== el_idx+1

            pdsec = pdbufr[ cond_var & cond_elev ]
        
            az  = pdsec["azimuth"].values     
            self.aux[el_idx]["azimuth"] = az
            
            radialtime = pdsec["radial_time"].values   
            radialtime_list.extend(radialtime)
            self.aux[el_idx]["radialtime"] = radialtime

            self.aux[el_idx]["elevation"] =  pdsec["elev"].values   

            # del(tref)
            # del(ref)
            # del(vel)
            # del(sw)
            # del(zdr)
            # del(kdp)
            # del(rho)            

        self._radialtime = radialtime_list
        
        self.scantime = radialtime_list[0]
        
        # breakpoint()

        return

    # ===================================================================================


    def get_raw( self, tilt: int, drange: Number_T, dtype: str) -> Union[np.ndarray, tuple]:   #! 2024-02-05 11:45:21 revised by chentao@cma.gov.cn
        r"""
        Get radar raw data

        Args:
            tilt (int): Index of elevation angle starting from zero.

            drange (float): Radius of data.

            dtype (str): Type of product (REF, VEL, etc.)

        Returns:
            numpy.ndarray or tuple of numpy.ndarray: Raw data
        """
        rf_flag = False

        self.tilt = tilt

        try:
            dmax = np.round(self.data[tilt][dtype][0].shape[0] * self.reso)
        except:
            dmax = self.scan_config['tch']["max_range1"][0,0] / 1000

        if dmax < drange:
            warnings.warn("Requested data range exceed max range in this tilt")

        self.drange = drange

        self.elev = self.el[tilt]
    
        ngates = int(drange // self.reso)

        try:
            data = np.ma.array(self.data[tilt][dtype])
            # The structure of `out_data`:
            # The key of `out_data` is the number of scan counting from zero (int).
            # The value of `out_data` is also dictionary, the key of it are the abbreviations of
            # product name, such as `REF`, `VEL`.
            # The value of this sub-dictionary is the data stored in np.ma.MaskedArray.
        except KeyError:
            raise RadarDecodeError("Invalid product name")     #! revised by 2024-02-05 11:44:50 
            # warnings.warn("Empty data", RuntimeWarning)
            # # Calculate size equivalent
            # nrays = len(self.aux[tilt]["azimuth"])
            # out = np.zeros((nrays, ngates)) * np.ma.masked       #! cod trick, created masked array for invalid moment data
            # return out
        
        if data.size == 0:
            warnings.warn("Empty data", RuntimeWarning)
            # Calculate size equivalent
            nrays = len(self.aux[tilt]["azimuth"])
            out = np.zeros((nrays, ngates)) * np.ma.masked
            return out

        cut = data.T[:ngates]
        shape_diff = ngates - cut.shape[0]
        append = np.zeros((int(np.round(shape_diff)), cut.shape[1])) * np.ma.masked
        if dtype in ["VEL", "SW"]:
            try:
                rf = self.data[tilt]["RF"]
            except KeyError:
                pass
            else:
                rf_flag = True
                rf = rf.T[:ngates]
                rf = np.ma.vstack([rf, append])
        r = np.ma.vstack([cut, append])
        if rf_flag:
            r.mask = np.logical_or(r.mask, ~rf.mask)
            ret = (r.T, rf.T)
        else:
            ret = r.T
        
        return ret
    
    

    def get_data(self, tilt: int, drange: Number_T, dtype: str) -> xr.Dataset:
        r"""
        Get radar data with extra information

        Args:
            tilt (int): Index of elevation angle starting from zero.

            drange (float): Radius of data.

            dtype (str): Type of product (REF, VEL, etc.)

        Returns:
            xarray.Dataset: Data.
        """
        # task = getattr(self, "task_name", None)
        task = ""
        ret = self.get_raw(tilt, drange, dtype)
        rf_flag = (dtype in ["VEL", "SW"]) and ("RF" in self.data[tilt])
        x, y, z, d, a = self.projection(self.reso)
        shape = ret[0].shape[1] if rf_flag else ret.shape[1]
        if rf_flag:
            da = xr.DataArray(ret[0], coords=[a, d], dims=["azimuth", "distance"])
        else:
            da = xr.DataArray(ret, coords=[a, d], dims=["azimuth", "distance"])
        ds = xr.Dataset(
            {dtype: da},
            attrs={
                "elevation": self.elev,
                "range": int(np.round(shape * self.reso)),
                "scan_time": self.scantime.strftime("%Y-%m-%d %H:%M:%S"),
                "site_code": self.code,
                "site_name": self.name,
                "site_longitude": self.stationlon,
                "site_latitude": self.stationlat,
                "tangential_reso": self.reso,
                # "nyquist_vel": self.nyquist_v[tilt],
                "task": task,
            },
        )
        ds["longitude"] = (["azimuth", "distance"], x)
        ds["latitude"] = (["azimuth", "distance"], y)
        ds["height"] = (["azimuth", "distance"], z)
        if rf_flag:
            ds["RF"] = (["azimuth", "distance"], ret[1])
        return ds

    def projection(self, reso: float) -> tuple:
        r = self.get_range(self.drange, reso)
        theta = np.array(self.aux[self.tilt]["azimuth"]) * deg2rad
        lonx, latx = get_coordinate(
            r, theta, self.elev, self.stationlon, self.stationlat
        )
        hght = (
            height(r, self.elev, self.radarheight)
            * np.ones(theta.shape[0])[:, np.newaxis]
        )
        return lonx, latx, hght, r, theta

    # ========================================================================


    def available_tilt(self, product: str):
        r"""Get all available tilts for given product"""
        tilt = list()
        for i in list(self.data.keys()):
            if product in self.data[i].keys():
                tilt.append(i)
        return tilt

    def iter_tilt(self, drange: Number_T, dtype: str):
        for i in self.available_tilt(dtype):
            yield self.get_data(i, drange, dtype)
    
# ======================================================================

def compose_radial_header(radial_header):

    radial_stat        = radial_header["radial_stat"][0]
    radial_elev        = radial_header["elevation"][0]
    radial_azimuth     = radial_header["azimuth"][0]
    radial_number      = radial_header["radial_number"][0]
    radial_datalength  = radial_header['length_data'][0]
    seq_number         = radial_header["seq_number"][0]
    elev_number        = radial_header["elev_number"][0]

    # radial_time       = datetime.datetime.utcfromtimestamp( data_radialheader["seconds"][0] )
    radial_time        = radial_header["seconds"][0]    #  int64,  #! key parameters
    moment_number      = radial_header["moment_number"][0]

    radial_header_info = { "radial_stat":radial_stat, 
                           "radial_datalength":radial_datalength, 
                           "radial_time":radial_time, 
                           "elev_number":elev_number, 
                           "radial_elev":radial_elev, 
                           "radial_number":radial_number, 
                           "radial_azimuth":radial_azimuth, 
                           "moment_number":moment_number  }

    return radial_header_info

# =========================================================================================================================

def decode_radial(radial_header_info, radialbufr):

    datalist = []

    radial_stat        = radial_header_info["radial_stat"]
    radial_datalength  = radial_header_info["radial_datalength"]
    radial_time        = radial_header_info["radial_time"]    #  int64, 
    elev_number        = radial_header_info["elev_number"]
    radial_elev        = radial_header_info["radial_elev"]
    radial_number      = radial_header_info["radial_number"]
    radial_azimuth     = radial_header_info["radial_azimuth"]
    moment_number      = radial_header_info["moment_number"]

    mpos_start = 0

    for imoment in np.arange(moment_number): 

        data_momentheader = np.frombuffer(radialbufr[mpos_start:mpos_start+SIZE_MOMENTHEADER], dtype=np.dtype(MOMENT_HEADER) )

        moment_type      = data_momentheader["data_type"][0]  
        moment_name      = moment_dict[ moment_type ]
        
        moment_scale     = data_momentheader["scale"][0]                
        moment_offset    = data_momentheader["offset"][0]
        
        moment_binlength = data_momentheader["binlength"][0]
        moment_length    = data_momentheader["length"][0]

        num_gates        = np.int32(moment_length/moment_binlength)
        
        mpos_end = mpos_start + SIZE_MOMENTHEADER + moment_length    #! byte position for certrain moment

        if(moment_binlength == 1):
            moment_data = np.frombuffer( radialbufr[ mpos_start+SIZE_MOMENTHEADER:mpos_end], 
                                         dtype = num_gates*np.dtype("B") )[0]
        elif(moment_binlength == 2):
            moment_data = np.frombuffer( radialbufr[ mpos_start+SIZE_MOMENTHEADER:mpos_end], 
                                         dtype = num_gates*np.dtype("H") )[0]
        
        data_rec = [ radial_stat, radial_time, elev_number, radial_elev, radial_number, radial_azimuth, 
                     moment_name, moment_scale, moment_offset, num_gates, moment_data ]

        datalist.append(data_rec)     #! put into dataframe

        print(data_rec[0:10])
        
        mpos_start = mpos_end         #! relocate moment byte  position

    # print(datalist[1])
    # breakpoint()
    return datalist

# ============================================================================================

def compose_var_sweep(pdbufr, varname, cut_index):

    cond_var = pdbufr["moment_name"] == varname
                                                         
    cond_elev = pdbufr["elev_number"] == cut_index + 1

    pdsec = pdbufr[ cond_var & cond_elev ]  #! trick，number of azimuth

    num_azimuth = len(pdsec)
    
    if(  num_azimuth > 0 ):
        
        # databufr = np.zeros( ( num_azimuth, max_numgates), np.float64  )    #! var initialization integet datatype
        
        num_gates = pdsec.iloc[0]["num_gates"]

        print([varname, " ilev: ", cut_index, num_azimuth, num_gates ])

        scale  = pdsec.iloc[0]["scale"]
        offset = pdsec.iloc[0]["offset"]

        # databufr[0:num_azimuth, 0:num_gates] = np.array(pdsec["moment_data"].explode()).reshape(num_azimuth, num_gates)   #! code trick, pay attention to pd.explode()         
        
        moment_bufr = np.array(pdsec["moment_data"].explode()).reshape(num_azimuth, num_gates)   #! code trick, pay attention to pd.explode()         

        mfilter = lambda x: FillValue if x < 5 else (x - offset)/float(scale)                      #! decoding function 
   
        var = np.vectorize(mfilter)(moment_bufr)                 #! to be imcproved
    
    else:

        # databufr[:, :] = FillValue

        var = []

        print(["Warning: ", varname, " ilev: ", cut_index, "NULL"])

    
    # var = databufr    
    # breakpoint()

    return var

# -============================================================================================

# rdafname = "//home/lse/haikui_tc2311/data/xrda/ZG003/Z_RADR_I_ZG003_20230907060004_O_DOR_AXPT0364_CRA_FMT.bin"
# rdafname = "/home/lse/haikui_tc2311/data/xrda/ZG000/Z_RADR_I_ZG000_20230907060007_O_DOR_AXPT0364_CRA_FMT.bin"

# fbufr = AXPT_Data(rdafname)

# breakpoint()
