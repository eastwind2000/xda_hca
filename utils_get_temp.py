
import numpy as np

import pandas as pd

import xarray as xr

# from metpy.interpolate import interpolate_to_points
# from scipy.spatial import cKDTree

from scipy.interpolate import interpn

from utils_rda_hca import *

import pdb

# ===============================================================================
# ===============================================================================
# ===============================================================================


def setup_g3km_model(rdatime_utc):    #! TODO: waiting for revising 

    rda_time_hour   = int(rdatime_utc[8:10])

    rda_time_min    = int(rdatime_utc[10:12])

    if (rda_time_hour < 12  ):
        model_inittime = rdatime_utc[0:8] + "00"

        if(rda_time_min < 30):
            fcst_chour = str(rda_time_hour).zfill(3)
        else:
            fcst_chour = str(rda_time_hour+1).zfill(3)
    else:
        model_inittime = rdatime_utc[0:8] + "12"

        fcst_chour = str(rda_time_hour-12).zfill(3)

        if(rda_time_min < 30):
            fcst_chour = str(rda_time_hour).zfill(3)
        else:
            fcst_chour = str(rda_time_hour+1).zfill(3)

    txtbufr = [ "./", rdatime_utc[0:8],  "/" ,
                "Z_NAFP_C_BABJ_", 
                model_inittime + "0000", 
                "_P_NWPC-GRAPES-3KM-ORIG-", 
                fcst_chour + "00", 
                ".grb2"       ]

    model_fname = ''.join(txtbufr)

    return model_fname


# =====================================================================================

def setup_era5_model(rdatime_utc):

    model_fname = "era5"

    return model_fname

# =====================================================================================
# =====================================================================================

def interpolate_temp_to_radar(radar, model_fname):
    """ Takes in  temperature data and interpolates it to every radar gate."""

    radar_z = get_z_from_radar(radar)     # 3994*1840,  narys*ngates

    radar_T = None

    shape = np.shape(radar_z)

    rad_z1d = radar_z.ravel()

    radar_site_lat = radar.latitude['data'][0]
    radar_site_lon = radar.longitude['data'][0]

    # model_name = "g3km"  

    databufr = xr.open_dataset(model_fname, engine="pynio")

    print(["NWP model: ", model_fname])

    temp00 = databufr.variables["TMP_P0_L100_GLL0"]
    hgt00 = databufr.variables["HGT_P0_L100_GLL0"]

    gridlat00 = databufr.variables["lat_0"]
    gridlon00 = databufr.variables["lon_0"]
    gridlev00 = databufr.variables["lv_ISBL0"]

    # ========================================================================================================

    cond_lat = (gridlat00 > radar_site_lat-5) & (gridlat00 < radar_site_lat+5)

    cond_lon = (gridlon00 > radar_site_lon-5) & (gridlon00 < radar_site_lon+5)

    tempbufr = temp00.data[:, cond_lat, :][:, :, cond_lon] - 273.15  # ! slice coding trick

    hgtbufr = hgt00.data[:, cond_lat, :][:, :, cond_lon]

    gridlev = gridlev00.data
    gridlat = gridlat00[cond_lat]
    gridlon = gridlon00[cond_lon]

    nz, ny, nx = tempbufr.shape

    # ================================== Put into xarray DataArray  =================================

    temp3d = xr.DataArray(data=tempbufr,
                          dims=["lev", "lat", "lon"],
                          coords=dict(lev=(["lev"], gridlev),
                                      lat=(["lat"], gridlat),
                                      lon=(["lon"], gridlon),
                                      ),
                          attrs=dict(description="Ambient temperature.",
                                     units="degC",
                                     ),
                          )

    hgt3d = xr.DataArray(data=hgtbufr,
                         dims=["lev", "lat", "lon"],
                         coords=dict(lev=(["lev"], gridlev),
                                     lat=(["lat"], gridlat),
                                     lon=(["lon"], gridlon),
                                     ),
                         attrs=dict(description="geopotential height",
                                    units="gpm",
                                    ),
                         )

    # print(temp3d)
    # print(hgt3d)

    # ! Averaged heigh on pressure levels centered at radar site.
    hgt1d_ave = hgt3d.interp(lat=radar_site_lat, lon=radar_site_lon)

    grid_coord = (hgt1d_ave.data[::-1], gridlat.data[::-1], gridlon.data)

    # ============================================================ adapted by chentao@cma.gov.cn 2023-01-25 00:27:57

    nsweeps = radar.nsweeps

    # pdb.set_trace()
    rtlist = []

    for isweep in np.arange(nsweeps):

        rda_lat, rda_lon, rda_lev = radar.get_gate_lat_lon_alt(sweep=isweep)

        num_rays, num_gates = rda_lat.shape

        print([isweep, num_rays, num_gates, rda_lat.shape])

        rtempbufr = np.zeros((num_rays, num_gates), dtype=float)
        px_coord = np.zeros((num_rays, num_gates, 3), dtype=float)

        # =======================  interpolate ray by ray ==============================

        for iray in np.arange(num_rays):
            # print(["iray: ", iray])
            px_coord[iray, :, 0] = rda_lev[iray, :]
            px_coord[iray, :, 1] = rda_lat[iray, :]
            px_coord[iray, :, 2] = rda_lon[iray, :]

            rtempbufr[iray, :] = interpn( grid_coord, temp3d.data[::-1, ::-1, :], px_coord[iray, :, :], bounds_error=False, fill_value=None)

        rtlist.append(rtempbufr.tolist())

    radar_T = np.vstack(rtlist)

    # pdb.set_trace()

    return radar_T, radar_z


# ===============================================================================
# ===============================================================================
# ===============================================================================

def setup_snd( snd_fname ):
    """
    ! waiting for fix , sounding-radar coehence
    """
   
    txtbufr = pd.read_csv(snd_fname,  header=None, skiprows=1, sep="\\s+")
    
    txtbufr.columns=("pressure", "hght", "temp", "dwpt", "wind_direction", "wind_speed")

    txtbufr['hght'] = txtbufr['hght']*10

    # cond_qc

    print("")
    print(["sounding data: ", snd_fname])
    print(txtbufr)
    
    return txtbufr


def interpolate_sounding_to_radar(sounding, radar):
    """Takes sounding data and interpolates it to every radar gate."""
    radar_z = get_z_from_radar(radar)
    radar_T = None
    snd_T, snd_z = check_sounding_for_montonic(sounding)
    shape = np.shape(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z

# ========================================================================================


def check_sounding_for_montonic(sounding):
    """
    So the sounding interpolation doesn't fail, force the sounding to behave
    monotonically so that z always increases. This eliminates data from
    descending balloons. 
    """
    snd_T = sounding['temp']  # In old SkewT, was sounding.data
    snd_z = sounding['hght']  # In old SkewT, was sounding.data
    dummy_z = []
    dummy_T = []
    #print(snd_T.mask)
    # if not snd_T.mask.all(): #May cause issue for specific soundings
    dummy_z.append(snd_z[0])
    dummy_T.append(snd_T[0])
    for i, height in enumerate(snd_z):
        if i > 0:
            if snd_z[i] > snd_z[i-1]:
                if np.isfinite(snd_z[i]) and np.isfinite(snd_T[i]):
                    #print(snd_T[i])
                    dummy_z.append(snd_z[i])
                    dummy_T.append(snd_T[i])
    snd_z = np.array(dummy_z)
    snd_T = np.array(dummy_T)
    # else:
    #     print('uh-oh. sounding problem')
    # pdb.set_trace()
    return snd_T, snd_z


# def bak_interpolate_temp_to_radar(radar, model_name):
#     """Takes temperature data and interpolates it to every radar gate."""

#     radar_z = get_z_from_radar(radar)

#     radar_T = None

#     shape = np.shape(radar_z)

#     rad_z1d = radar_z.ravel()

#     radar_site_lat = radar.latitude['data'][0]
#     radar_site_lon = radar.longitude['data'][0]

#     model_name = "g3km"

#     fname = "/mnt/i/fcst2023/coldwave/g3km/20230114/Z_NAFP_C_BABJ_20230114000000_P_NWPC-GRAPES-3KM-ORIG-00600.grb2"

#     databufr = xr.open_dataset(fname, engine="pynio")

#     print(["NWP model: ", fname])

#     tempbufr = databufr.variables["TMP_P0_L100_GLL0"]
#     hgtbufr = databufr.variables["HGT_P0_L100_GLL0"]

#     gridlat = databufr.variables["lat_0"]
#     gridlon = databufr.variables["lon_0"]
#     gridlev = databufr.variables["lv_ISBL0"]

#     # ================================== Put into xarray DataArray structure =================================

#     temp3d = xr.DataArray(data=tempbufr,

#                           dims=["lev", "lat", "lon"],

#                           coords=dict(lev=(["lev"], gridlev),
#                                       lat=(["lat"], gridlat),
#                                       lon=(["lon"], gridlon),
#                                       ),
#                           attrs=dict(description="Ambient temperature.", units="degC",  ),
#                           )

#     hgt3d = xr.DataArray(data=hgtbufr,

#                          dims=["lev", "lat", "lon"],

#                          coords=dict(lev=(["lev"], gridlev),
#                                      lat=(["lat"], gridlat),
#                                      lon=(["lon"], gridlon),
#                                      ),
#                          attrs=dict(description="geopotential height",  units="gpm",   ),
#                          )

#     print(temp3d)

#     print(hgt3d)

#     # ! Averaged heigh on pressure levels centered at radar site.
#     hgt1d_ave = hgt3d.interp(lat=radar_site_lat, lon=radar_site_lon)

#     # ==========================================================================================================

#     # TODO:  watiting for optimization,  #### 2023-01-19 21:45:59
#     nsweeps = radar.nsweeps
#     rda_lat = np.zeros((nsweeps), dtype=float)

#     rda_lat00, rda_lon00, rda_lev00 = radar.get_gate_lat_lon_alt(sweep=0)
#     rda_lat01, rda_lon01, rda_lev01 = radar.get_gate_lat_lon_alt(sweep=1)
#     rda_lat02, rda_lon02, rda_lev02 = radar.get_gate_lat_lon_alt(sweep=2)
#     rda_lat03, rda_lon03, rda_lev03 = radar.get_gate_lat_lon_alt(sweep=3)
#     rda_lat04, rda_lon04, rda_lev04 = radar.get_gate_lat_lon_alt(sweep=4)
#     rda_lat05, rda_lon05, rda_lev05 = radar.get_gate_lat_lon_alt(sweep=5)
#     rda_lat06, rda_lon06, rda_lev06 = radar.get_gate_lat_lon_alt(sweep=6)
#     rda_lat07, rda_lon07, rda_lev07 = radar.get_gate_lat_lon_alt(sweep=7)
#     rda_lat08, rda_lon08, rda_lev08 = radar.get_gate_lat_lon_alt(sweep=8)
#     rda_lat09, rda_lon09, rda_lev09 = radar.get_gate_lat_lon_alt(sweep=9)
#     rda_lat10, rda_lon10, rda_lev10 = radar.get_gate_lat_lon_alt(sweep=10)

#     rda_lat = np.vstack((rda_lat00, rda_lat01, rda_lat02, rda_lat03, rda_lat04,
#                         rda_lat05, rda_lat06, rda_lat07, rda_lat08, rda_lat09, rda_lat10))
#     rda_lon = np.vstack((rda_lon00, rda_lon01, rda_lon02, rda_lon03, rda_lon04,
#                         rda_lon05, rda_lon06, rda_lon07, rda_lon08, rda_lon09, rda_lon10))
#     rda_lev = np.vstack((rda_lev00, rda_lev01, rda_lev02, rda_lev03, rda_lev04,
#                         rda_lev05, rda_lev06, rda_lev07, rda_lev08, rda_lev09, rda_lev10))

#     print([rda_lat.shape])

#     rlatbufr = rda_lat.flatten()
#     rlonbufr = rda_lon.flatten()
#     rlevbufr = rda_lev.flatten()

#     # rad_T1d = np.interp(rad_z1d, snd_z, snd_T)

#     num_rec = len(rda_lat.flatten())

#     # p1d     = gridlev.data

#     nrays = radar.nrays
#     ngates = radar.ngates

#     rtempbufr = np.zeros((nrays * ngates), dtype=float)

#     #!  grid_coord = ( gridlev.data, gridlat.data[::-1], gridlon.data   )   # TODO: gridlev on pressure levels ~> approximately on height levels
#     grid_coord = (hgt1d_ave.data[::-1], gridlat.data[::-1], gridlon.data)

#     # px_coord  = np.array( [ rlevbufr[0:10], rlatbufr[0:10], rlonbufr[0:10]   ] ).transpose()
#     # px_coord = np.array( [ [1500, 32, 115],    [1500, 32, 115]  ]        )

#     for irec in np.arange(num_rec):

#         px_coord = np.array([rlevbufr[irec], rlatbufr[irec], rlonbufr[irec]])

#         rtempbufr[irec] = interpolate.interpn(
#             grid_coord, temp3d[::-1, ::-1, :], px_coord, bounds_error=False, fill_value=None)

#         # print([irec,  px_coord, " : ",   rtempbufr[irec]])

#     pdb.set_trace()

#     # xx= interpolate.interpn( grid_coord, hgt3d[::-1, ::-1, :], px_coord.T, fill_value="extrapolate" )

#     pdb.set_trace()

#     #! TODO: accurate interpolaton,  however too slowly ================================================
#     # for irec in np.arange(num_rec):
#     #     print(irec)
#     #     grid_h1d = hgt3d.interp( lat=rlatbufr[irec], lon=rlonbufr[irec]  )
#     #     grid_t1d = temp3d.interp( lat=rlatbufr[irec], lon=rlonbufr[irec] )
#     #     inp1d = interpolate.interp1d(grid_h1d, grid_t1d, bounds_error=False, fill_value="extrapolate")
#     #     rtempbufr[irec] = inp1d( rlevbufr[irec] )
#     #! ==================================================================================================

#     pdb.set_trace()

#     radar_T = rtempbufr.reshape(nrays, ngates)

#     return radar_T, radar_z
