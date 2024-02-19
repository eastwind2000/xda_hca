
import numpy as np

common_header = [   ("magic_number",   "S4"),
                    ("major_version",  "u2"),
                    ("minor_version",  "u2"),
                    ("generic_type",   "i4"),
                    ("file_type",      "i4"),
                    ("res00",          "16B")                              
                ]

siteconfig_header = [   ("sitecode",        "S8"),
                        ("sitename",        "S32"),
                        ("latitude",        "f4"),
                        ("longitude",       "f4"),
                        ("antenna_height",  "i4"),
                        ("ground_height",   "i4"),
                        ("frequency",       "i4"),
                        ("beam_width_hori",       "f4"),
                        ("beam_width_vert",       "f4"),
                        ("rda_version",              "i4"),
                        ("radar_type",               "H"),
                        ("antenna_gain",             "H"),
                        ("transmitting_feeder_loss", "H"),
                        ("receiving_feeder_loss",    "H"),
                        ("other_loss",               "H"),
                        ("res01",                    "46B"),          
                    ]   


taskconfig_header = [   ("taskname",              "S32"),
                        ("taskdes",              "S128"),
                        ("pol_type",               "i4"),
                        ("scan_type",              "i4"),
                        ("pulse_width",              "i4"),
                        ("scan_starttime",          "i4"),
                        ("cut_number",             "i4"),
                        ("horizontal_noise",       "f4"),
                        ("vertical_noise",             "f4"),
                        ("horizontal_cali",       "f4"),
                        ("vertical_cali",             "f4"),
                        ("HN_temp",           "f4"),
                        ("VN_temp",             "f4"),
                        ("zdr_cali",           "f4"),
                        ("phi_cali",             "f4"),
                        ("ldr_cali",             "f4"),                                         
                        ("res",                   "40B"),                                            
                    ]

# ===================================================================

# beam_header = [     ("beam_index",            "i4"),
#                     ("beam_type",             "i4"),
#                     ("subplse_number",        "i4"),
#                     ("txbeam_direction",      "f4"),
#                     ("txbeam_widthH",         "f4"),
#                     ("txbeam_widthV",         "f4"),
#                     ("txbeam_gain",           "f4"),
#                     ("res",                 "100B"),
#                 ]

# subpulse_header = [   ("sp_stategy",            "i4"),
#                       ("sp_modulation",         "i4"),
#                       ("sp_frequency",          "f4"),
#                       ("sp_bandwidth",          "f4"),
#                       ("sp_width",              "i4"),
#                       ("horizontal_noise",      "f4"),
#                       ("vertical_noise",        "f4"),
#                       ("horizotnal_cali",       "f4"),
#                       ("vertical_cali",         "f4"),
#                       ("noise_temperature",    "2f4"),
#                       ("zdr_cali",              "f4"),
#                       ("phidp_cali",            "f4"),
#                       ("ldr_cali",              "f4"),
#                       ("pulse_points",           "H"),
#                       ("res",                  "70B")
#                      ]

# comp_beamheader = [  ("bh",      np.dtype(beam_header) ),
#                      ("subph", 4*np.dtype(subpulse_header) )
#                      ]

# =====================================================================

cut_header = [     ("process_mode",           "i4"),
                   ("wave_form",              "i4"),
                   ("prf1",                   "f4"),
                   ("prf2",                   "f4"),
                   ("deliasing_mode",         "i4"),
                   ("azimuth",                "f4"),
                   ("elevation",              "f4"),
                   ("start_angle",            "f4"),
                   ("end_angle",              "f4"),
                   ("angular_reso",           "f4"),
                   ("scan_speed",             "f4"),
                   ("log_resolution",         "i4"),
                   ("doppler_resolution",     "i4"),
                   ("max_range1",             "i4"),
                   ("max_range2",             "i4"),
                   ("start_range",            "i4"),
                   ("sample1",                "i4"),
                   ("sample2",                "i4"),
                   ("phase_mode",             "i4"),
                   ("atmos_loss",             "f4"),
                   ("nyquist_speed",          "f4"),
                   ("moment_mask",            "i8"),
                   ("moment_size",            "i8"),
                   ("filter_mask",            "i4"),
                   ("var_threshold",         "8f4"),
                   ("dbt_mask",               "i4"),
                   ("dbz_mask",               "i4"),
                   ("vel_mask",               "i4"),
                   ("sw_mask",                "i4"),
                   ("dp_mask",                "i4"),
                   ("mask_res",              "12B"),
                   ("res00",                  "4B"),
                   ("direction",              "i4"),
                   ("gc_type",                 "H"),
                   ("gc_filter",              "3H"),
                   ("res01",                 "72B"),
                ]


# =====================================================================

axpt_comheader  = np.dtype(common_header)

axpt_siteheader = np.dtype(siteconfig_header)

axpt_taskheader = np.dtype(taskconfig_header)

# =====================================================================

radial_header = [ ("radial_stat",        "i4"),
                  ("spot_blank",         "i4"),
                  ("seq_number",         "i4"),
                  ("radial_number",      "i4"),
                  ("elev_number",        "i4"),
                  ("azimuth",            "f4"),
                  ("elevation",          "f4"),
                  ("seconds",            "i4"),
                  ("mircrosecconds",     "i4"),
                  ("length_data",        "i4"),
                  ("moment_number",      "i4"),
                  ("res00",              "2B"),
                  ("noise",              "2H"),
                  ("res01",             "14B"),
                 ] 

SIZE_RADIALHEADER = np.dtype(radial_header).itemsize

MOMENT_HEADER = [     ("data_type",               "i4"),
                      ("scale",                   "i4"),
                      ("offset",                  "i4"),
                      ("binlength",               "h"),
                      ("flags",                   "h"),
                      ("length",                  "i4"),
                      ("res",                    "12B"),
                    ]
SIZE_MOMENTHEADER = np.dtype(MOMENT_HEADER).itemsize

# =====================================================================

AXPT_HEADER00  = [ ("common_header",     axpt_comheader), 
                   ("siteconfig_header", axpt_siteheader ),
                   ("taskconfig_header", axpt_taskheader )                           
                 ]

SIZE_AXPT_HEADER00 = np.dtype(AXPT_HEADER00).itemsize

moment_dict = dict( [ (1,  "DBT"), 
                      (2,  "DBZ"),
                      (3,  "vel"),
                      (4,  "sw"),
                      (5,  "SQI"),
                      (7,  "zdr"),
                      (9,  "rho"),
                      (10, "phi"),
                      (11, "kdp"),
                      (16, "snrh"),
                      (17, "snrv"),
                        ]   
                   )

FillValue = -9999.0

max_numgates = 3000   #! predefined value

max_numazimuth =  500
