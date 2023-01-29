# coding = utf-8
# Author: Puyuan Du

from typing import Callable
from functools import wraps
from cinrad.io.level2 import PhasedArrayData

import numpy as np
import pdb

try:
    import pyart

    PYART_INSTALLED = True
except ImportError:
    PYART_INSTALLED = False


def check_pyart_installed(func: Callable) -> Callable:
    @wraps(func)
    def deco(*args, **kwargs):
        if not PYART_INSTALLED:
            raise ImportError("pyart is not installed")
        return func(*args, **kwargs)

    return deco


# ==============================================================================
# ================== adding by chent 202206=====================================
# ==============================================================================
# ==============================================================================


mapping_xpar = {    "TREF": "corrected reflectivity",
                    "REF": "reflectivity",
                    "VEL": "velocity",
                    "SW":  "spectrum_width",
                    "PHI": "differential_phase",
                    "ZDR": "differential_reflectivity",
                    "RHO": "cross_correlation_ratio",
                    "KDP": "specific_differential_phase",  # adding by chentao 2022.05.26 
                }


@check_pyart_installed
def xpar_data_to_pyart(f: PhasedArrayData, radius: int = 60) -> pyart.core.Radar:

    """     
     radial_num = 2000 ? may be num_of_gates
     f.data[0-12] , 12 elvevation angels
     f.data[0]["REF"].shapes :-> (400, 2000)
     means eachs scan : 400 azimuths, 2000 gates     
     raidials : 4800, 400 * 12 sweeps
     f._d["data"]["radial_time"]
     f._d["data"]["az_num"] = 400
     el=12
     f._d["data"]["gate_num"] = [1433, 1433], len=4800
     f._d["data"]["gate_length"] = [1996, 1996,.....], len=4800    \
     f._d["data"]["radial_index"] : 0->4799    
     f._d["data"]["azimuth"] : 0->4799
     REF : 400*12*2000 ,real gate 1996
     f._d["data"]["gate_num"] = [1433]*4800
     ! gate_num*[reso=0.03] = 42990 Radar range about 42.9km
     f._d["data"]["gate_length"] = [1966]*4800
     REF,ZDR,KDP matrix: 4800*2000 

     #============ moved to WSL platform 20220617 ========
     
    """

    filemetadata = pyart.config.FileMetadata("XPAR Radar")
    time = filemetadata("time")
    time["calendar"] = "standard"
    time["units"] = "seconds since 1970-01-01 00:00"
    time["standard_name"] = "time"
    time["long_name"] = "time in seconds since volume start"

    # pdb.set_trace()
    # time["data"] = f._time_radial
    time["data"] =  f._d["data"]["radial_time"]

    # ===============================================================

    _range = filemetadata("range")
    
    # reso = f.scan_config[0].dop_reso
    
    reso = 1000.0 * f.reso                             #! unit in m, generally XPAR ~ 30m in horizontal resolution

    _range["data"] = f.get_range(radius * 1000, reso)  #! unit in m

    _range["meters_to_center_of_first_gate"] = float(reso)
    _range["meters_between_gates"] = reso

    metadata = filemetadata("metadata")
    metadata["original_container"] = "CINRAD_XPAR"
    
    # vcp_pattern = f.scan_type  # revised by chentao 2022/06/10
    vcp_pattern = "XPAR-VCP11"

    metadata["instrument_name"] = f.code

    scan_type = "ppi"

    latitude = filemetadata("latitude")
    longitude = filemetadata("longitude")
    altitude = filemetadata("altitude")

    latitude["data"] = np.array([f.stationlat], dtype="float64")
    longitude["data"] = np.array([f.stationlon], dtype="float64")
    altitude["data"] = np.array([f.radarheight], dtype="float64")

    sweep_number = filemetadata("sweep_number")
    sweep_mode = filemetadata("sweep_mode")
    sweep_start_ray_index = filemetadata("sweep_start_ray_index")
    sweep_end_ray_index = filemetadata("sweep_end_ray_index")

    nsweeps = len(f.el)
    sweep_number["data"] = np.arange(nsweeps, dtype="int32")
    sweep_mode["data"] = np.array(nsweeps * ["azimuth_surveillance"], dtype="S")

    #
    # sweep_end_ray_index["data"] = np.array(f._sweep_end_ray_index)
    # sweep_start_ray_index["data"] = np.array(f._sweep_start_ray_index)
    # sweep_end_ray_index["data"] = np.array( []  )
    # sweep_start_ray_index["data"] = np.array()
    #

    # sweep_end_ray_index["data"] = np.arange(0, 4800, 400)
    # sweep_start_ray_index["data"] = np.arange(399, 4800, 400)

    sweep_end_ray_index["data"] = np.zeros(nsweeps, dtype=int)

    num_radials_iswp =  np.array( [ len(f.aux[i]["azimuth"]) for i in np.arange(nsweeps)    ] )

    sweep_start_ray_index["data"] = np.array([ int(np.sum(num_radials_iswp[0:i])) for i in np.arange(12) ])

    sweep_end_ray_index["data"]  = sweep_start_ray_index["data"] + num_radials_iswp - 1
    
    # pdb.set_trace()


    azimuth = filemetadata("azimuth")
    azimuth["data"] = np.hstack([f.aux[i]["azimuth"] for i in f.aux.keys()])
    
    # elevation["data"] = np.hstack([f.aux[i]["elevation"] for i in f.aux.keys()])
    elevation = filemetadata("elevation")
    elevbufr = np.zeros( (nsweeps*num_radials_iswp[0]))
    for iswp in np.arange(nsweeps):
        elevbufr[ iswp*400:(iswp+1)*400 ] = f.el[iswp] 
    elevation["data"] = elevbufr

    fixed_angle = filemetadata("fixed_angle")
    fixed_angle["data"] = np.array(f.el)

    fields = {}
    nscans = f.get_nscans()

    for mom in mapping_xpar.keys():
        name = mapping_xpar[mom]
        dic = filemetadata(name)
        dic["_FillValue"] = pyart.config.get_fillvalue()
        raw_arr = [f.get_raw(nel, radius, mom) for nel in range(nscans)]
        sel_arr = [i if not isinstance(i, tuple) else i[0] for i in raw_arr]
        moment_data = np.ma.vstack(sel_arr)
        dic["data"] = moment_data
        fields[name] = dic

    
    unambiguous_range = filemetadata("unambiguous_range")
    unambiguous_range["data"] = f._d["data"]["unambiguous_dist"]   # 432
    
    nyquist_velocity = filemetadata("nyquist_velocity")
    nyquist_velocity["data"] = f._d["data"]["unambiguous_dist"]
    
    # nyquist_velocity["data"] = np.array(
    #      [i.nyquist_spd for i in f.scan_config], "float32"
    # # )


    # unambiguous_range["data"] = np.array(
    #     [i.max_range1 for i in f.scan_config], "float32"
    # )

    instrument_parameters = {
        "unambiguous_range": unambiguous_range,
        "nyquist_velocity": nyquist_velocity,
    }

    gate_num = f._d["data"]["gate_num"]

    gate_length = f._d["data"]["gate_length"]
    
    # pdb.set_trace()
    

    radar = pyart.core.Radar(
        time,
        _range,
        fields,
        metadata,
        scan_type,
        latitude,
        longitude,
        altitude,
        sweep_number,
        sweep_mode,
        fixed_angle,
        sweep_start_ray_index,
        sweep_end_ray_index,
        azimuth,
        elevation,
        instrument_parameters=instrument_parameters,
    )
    
    # pdb.set_trace()

    return radar

    