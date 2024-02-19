
"""
2024-01-28 15:39:27

Created by chentao@cma.gov.cn

"""

# from cinrad.io.level2 import PhasedArrayData

import numpy as np

import datetime

import pyart

from cinrad.io.level2 import StandardData

# from cinrad.io

from module_io_XPT import XPT_DATA

from module_io_defineXPT import FillValue

# ==============================================================================
# ==============================================================================

vcp_pattern = {  }

mapping = {
    "REF": "reflectivity",
    "VEL": "velocity",
    "SW": "spectrum_width",
    "PHI": "differential_phase",
    "ZDR": "differential_reflectivity",
    "RHO": "cross_correlation_ratio",
    "KDP": "specific_differential_phase",  # adding by chentao 2022.05.26
    "SNRH": "snrh"
}


mapping_xpar = {    "TREF": "corrected reflectivity",
                    "REF":  "reflectivity",
                    "VEL":  "velocity",
                    "SW":   "spectrum_width",
                    "PHI":  "differential_phase",
                    "ZDR":  "differential_reflectivity",
                    "RHO":  "cross_correlation_ratio",
                    "KDP":  "specific_differential_phase",  # adding by chentao 2022.05.26 
                }

# mom_log_vars = ["TREF", "REF", "PHI", "ZDR", "RHO", "KDP" ]
# mom_dop_vars = ["VEL", "SW"  ]

mom_group = {  "LOG":"REF", 
               "DOP":"VEL",
               "log_vars": ["REF", "PHI", "ZDR", "RHO", "KDP"],
               "dop_vars": ["VEL", "SW"  ]
                }

# ==========================================================================================   


def XPT_to_pyart(f: XPT_DATA, radius: int = 60, mom_kind: str = "LOG") -> pyart.core.Radar:

    filemetadata = pyart.config.FileMetadata("XPAR Radar")

    rx = radius * 1000                                           #! key parameters to control number of valid gates
    
    valid_swplist = f.available_tilt(mom_group[mom_kind])

    if(mom_kind=="LOG"):
        mom_vars  = mom_group["log_vars"]
    
    if(mom_kind=="DOP"):
        mom_vars  = mom_group["dop_vars"]
    
    time                    = filemetadata("time")
    time["calendar"]        = "standard"
    time["units"]           = "seconds since 1970-01-01 00:00"
    time["standard_name"]   = "time"
    time["long_name"]       = "time in seconds since volume start"

    # radialtime = np.array(f.aux[valid_swplist ]["radialtime"]).flatten()
    # time["data"] =  f._radialtime   # Should be toal_raidials, [elev, num_radial]
    
    radialtime = []
    for iswp in valid_swplist:
        radialtime.extend( f.aux[iswp]["radialtime"]  )

    time["data"] = np.array(radialtime)

    print( ["Scan Time: ", datetime.datetime.utcfromtimestamp(f.scantime) ] ) 

    # ===============================================================

    _range = filemetadata("range")
    
    # reso = f.scan_config[0].dop_reso
    
    reso = f.reso                             #! unit in m, generally XPAR ~ 30m in horizontal resolution

    _range["data"] = f.get_range(rx, reso)    #! rx parameter, unit in m

    _range["meters_to_center_of_first_gate"] = f.first_gate_dist

    _range["meters_between_gates"] = reso

    metadata = filemetadata("metadata")

    metadata["original_container"] = "CHINA_AXPT"
    
    scan_type = "other"

    vcp_pattern = "VCP23D"

    metadata["instrument_name"] = f.code

    metadata["vcp_pattern"] = vcp_pattern

    latitude = filemetadata("latitude")
    longitude = filemetadata("longitude")
    altitude = filemetadata("altitude")

    latitude["data"]  = np.array([f.stationlat], dtype="float64")
    longitude["data"] = np.array([f.stationlon], dtype="float64")
    altitude["data"]  = np.array([f.radarheight], dtype="float64")

    nsweeps =  len(valid_swplist)

    sweep_number          = filemetadata("sweep_number")
    sweep_mode            = filemetadata("sweep_mode")
    sweep_start_ray_index = filemetadata("sweep_start_ray_index")
    sweep_end_ray_index   = filemetadata("sweep_end_ray_index")

    sweep_number["data"] = np.arange(nsweeps, dtype="int32")
    sweep_mode["data"]   = np.array(nsweeps * ["azimuth_surveillance"], dtype="S")

    num_radials_iswp =  np.array( [ len(f.aux[i]["azimuth"]) for i in valid_swplist    ] )

    sweep_start_ray_index["data"] = np.array([ int(np.sum(num_radials_iswp[0:i])) for i in np.arange(nsweeps) ])

    sweep_end_ray_index["data"]  = sweep_start_ray_index["data"] + num_radials_iswp -1
    
    # pdb.set_trace()

    azimuth         = filemetadata("azimuth")
    # azimuth["data"] = np.hstack([f.aux[i]["azimuth"] for i in f.aux.keys()])
    azimuth["data"] = np.hstack([f.aux[i]["azimuth"] for i in valid_swplist])
    
    
    elevation = filemetadata("elevation")
    # elevation["data"] = np.hstack([f.aux[i]["elevation"] for i in f.aux.keys()])
    elevation["data"] = np.hstack([f.aux[i]["elevation"] for i in valid_swplist])


    fixed_angle = filemetadata("fixed_angle")
    fixed_angle["data"] = f.el[valid_swplist]

    fields = {}
    nscans = f.get_nscans()

    # for mom in mapping_xpar.keys():
    for mom in mom_vars:
        name = mapping_xpar[mom]
        dic = filemetadata(name)
        dic["_FillValue"] = pyart.config.get_fillvalue()
    
        raw_arr = [f.get_raw(nel, rx, mom) for nel in valid_swplist]
        sel_arr = [i if not isinstance(i, tuple) else i[0] for i in raw_arr]
        
        moment_data = np.ma.vstack(sel_arr)
        moment_data.fill_value = FillValue                                     #! code trick for undef values 2024-01-26 23:01:52 by chenetao@cma.gov.cn
        
        dic["data"] = moment_data
        fields[name] = dic

    unambiguous_range = filemetadata("unambiguous_range")
    # unambiguous_range["data"] = f._d["data"]["unambiguous_dist"]   # 432 :?  #! to be rvised
    unambiguous_range["data"] =  0   # 432 :?  #! to be rvised

    nyquist_velocity         = filemetadata("nyquist_velocity")
    nyquist_velocity["data"] = f.nyquist_speed
    
    # nyquist_velocity["data"] = np.array(
    #      [i.nyquist_spd for i in f.scan_config], "float32"
    #       # )
    # unambiguous_range["data"] = np.array(
    #     [i.max_range1 for i in f.scan_config], "float32"
    # )

    instrument_parameters = { "unambiguous_range": unambiguous_range,
                               "nyquist_velocity": nyquist_velocity,
                            }

    # gate_num = f._d["data"]["gate_num"]
    # gate_length = f._d["data"]["gate_length"]
    # pdb.set_trace()
   
    radar = pyart.core.Radar(   time                = time,
                                _range              = _range,
                                fields              = fields,
                                metadata            = metadata,
                                scan_type           = scan_type,
                                latitude            = latitude,
                                longitude           = longitude,
                                altitude            = altitude,
                                sweep_number        = sweep_number,
                                sweep_mode          = sweep_mode,
                                fixed_angle         = fixed_angle,
                                sweep_start_ray_index   = sweep_start_ray_index,
                                sweep_end_ray_index     = sweep_end_ray_index,
                                azimuth                 = azimuth,
                                elevation               = elevation,
                                instrument_parameters   = instrument_parameters,
                            )
    

    # breakpoint()
    
    return radar


# ==========================================================================================



def standard_data_to_pyart(f: StandardData, radius: int = 460) -> pyart.core.Radar:
    filemetadata = pyart.config.FileMetadata("cinrad standard")
    time = filemetadata("time")
    time["calendar"] = "standard"
    time["units"] = "seconds since 1970-01-01 00:00"
    time["standard_name"] = "time"
    time["long_name"] = "time in seconds since volume start"

    time["data"] = f._time_radial

    # pdb.set_trace()

    _range = filemetadata("range")
    reso = f.scan_config[0].dop_reso
    _range["data"] = f.get_range(radius * 1000, reso)
    _range["meters_to_center_of_first_gate"] = float(reso)
    _range["meters_between_gates"] = float(reso)

    metadata = filemetadata("metadata")
    metadata["original_container"] = "CINRAD"
    vcp_pattern = f.scan_type
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

    sweep_end_ray_index["data"] = np.array(f._sweep_end_ray_index)
    sweep_start_ray_index["data"] = np.array(f._sweep_start_ray_index)

    azimuth = filemetadata("azimuth")
    elevation = filemetadata("elevation")
    fixed_angle = filemetadata("fixed_angle")
    azimuth["data"] = np.hstack([f.aux[i]["azimuth"] for i in f.aux.keys()])
    elevation["data"] = np.hstack([f.aux[i]["elevation"] for i in f.aux.keys()])
    fixed_angle["data"] = np.array(f.el)

    fields = {}
    nscans = f.get_nscans()

    for mom in mapping.keys():
        name = mapping[mom]
        dic = filemetadata(name)
        dic["_FillValue"] = pyart.config.get_fillvalue()
        raw_arr = [f.get_raw(nel, radius, mom) for nel in range(nscans)]
        sel_arr = [i if not isinstance(i, tuple) else i[0] for i in raw_arr]
        moment_data = np.ma.vstack(sel_arr)
        dic["data"] = moment_data
        fields[name] = dic

    nyquist_velocity = filemetadata("nyquist_velocity")
    nyquist_velocity["data"] = np.array(
        [i.nyquist_spd for i in f.scan_config], "float32"
    )
    unambiguous_range = filemetadata("unambiguous_range")
    unambiguous_range["data"] = np.array(
        [i.max_range1 for i in f.scan_config], "float32"
    )

    instrument_parameters = {
        "unambiguous_range": unambiguous_range,
        "nyquist_velocity": nyquist_velocity,
    }

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
    
    return radar

