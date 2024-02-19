import ssl
import matplotlib.pyplot as plt
import os
import numpy as np
import argparse

import pyart

import warnings 

import pandas as pd

import bz2

# from sklearn.covariance import log_likelihood

from csu_radartools import ( csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, csu_kdp, csu_misc, fundamentals )

from module_snd import *

from module_io_XPT import XPT_DATA

from module_io_XPT2pyart import  standard_data_to_pyart, XPT_to_pyart

from module_plot_rda  import *

from module_xpt_gridding import xpt_gridding

from cinrad.io import StandardData

# ===========================================================================================

warnings.filterwarnings("ignore")

parser = argparse.ArgumentParser(description='manual to this script, input Standard SPOL-RADAR Level2 Dasta')

parser.add_argument('--rdafname', type=str, default = None)

parser.add_argument('--figroot', type=str, default = None)

args = parser.parse_args()

rdafname = args.rdafname

figroot= args.figroot

# pdb.set_trace()
# ===============================================================================
# figroot = "../remote_figs/"

rband        = "X"

xpar_radius  = 75.0   #! unit in km

mom_log_list = ["reflectivity",  "differential_reflectivity", "specific_differential_phase", "qpe"   ]

mom_dop_list = ["velocity"]

rng_list = [20.0, 40.0, 60.0, 80.0, 100.0]  

plot_radius = xpar_radius  #! units in km

#
#
# rdafname = "//home/lse/haikui_tc2311/data/xrda/ZG003/Z_RADR_I_ZG003_20230907060004_O_DOR_AXPT0364_CRA_FMT.bin"

# ===============================================================================
#

absfname = os.popen("basename " + rdafname).readline().strip("\n")

datapath = os.popen("dirname " + rdafname).readline().strip("\n")

rda_site = absfname[9:14]

ftime_utc = absfname[15:15+14]

print([rda_site, ftime_utc])

#! ===================================================== Reading in AXPT data and transform to ART radar object

decode_method = 0

if (decode_method == 0):

    fbufr = XPT_DATA(rdafname)

    # rdabufr = XPT_data_to_pyart(fbufr, radius=xpar_radius)

    rdabufr_log = XPT_to_pyart(fbufr, radius=xpar_radius, mom_kind="LOG")
    
    rdabufr_dop = XPT_to_pyart(fbufr, radius=xpar_radius, mom_kind="DOP")
    
    rdabufr     = rdabufr_log
    
else:

    fbufr = StandardData(rdafname)

    rdabufr_log = standard_data_to_pyart(fbufr, radius=xpar_radius)  #! also works ! allmost same VolumnScan method

    rdabufr_dop = rdabufr_log

    rdabufr     = rdabufr_log

# breakpoint()
# ==========================================================

print(rdabufr.info())

print(rdabufr.elevation)

print(["nsweeps: ", rdabufr.nsweeps])

rda_lat = rdabufr.latitude['data'][0]

rda_lon = rdabufr.longitude['data'][0]

# ========================================================= setup sounding temperature
# TODO : need to revise temerature algorithm

snd_fname =  "../data/snd54511_2023073000UTC.dat"

sndbufr = setup_snd(snd_fname)

# snddata = SkewT.Sounding(snd_fname).soundingdata

wh0 = np.where(np.isclose(np.abs(sndbufr['temp']), 0.0, atol=1.0))
expected_ML = np.array(sndbufr['hght'])[wh0[0]][0]/1000.

####Interpolate sounding to the radar heights

radar_T, radar_z = interpolate_sounding_to_radar(sndbufr, rdabufr)  # 

rheight = radar_z/1000.0                     #!  rheight in km ?

# rdabufr.shape

rdabufr = add_field_to_radar_object( radar_T, rdabufr, field_name='temperature', units="degree", 
                                     dz_field='reflectivity')
# breakpoint()
# ============================================================= 

# pdb.set_trace()
#
# sn = pyart.retrieve.simple_moment_calculations.calculate_snr_from_reflectivity(rdabufr, refl_field='reflectivity',toa=15000.0)
# sndat= sn['data'][:]
#add SNR to the radar structure
# radar = add_field_to_radar_object(sndat, rdabufr, field_name='SN', dz_field='reflectivity')
# sndat =np.ma.masked_array(sndat)
# pdb.set_trace()
###If you use SNR, select an appropriate threshold for your data.
#

snthresh = -30

###Azimuths are important for the PPI data for the bright band detection.
azimuths = rdabufr.azimuth['data']

##The melting layer algorithm for PPIs takes different sectors. The easiest is to do 10º sectors, and pass nsect=36 (also
# the default). But you could change to any multiple of 360 (or the number of azimuths you have).
nsect = 10

#Radar variables. These should be just the data (not an object), and masked arrays

ref = rdabufr.fields['reflectivity']['data']
vel = rdabufr_dop.fields['velocity']['data']
zdr = rdabufr.fields['differential_reflectivity']['data']
kdp = rdabufr.fields['specific_differential_phase']['data']
rho = rdabufr.fields['cross_correlation_ratio']['data']
phi = rdabufr.fields['differential_phase']['data']

# breakpoint()


##Now what if we use trapazoidal functions and the linear method?
# CSU-HID types:           Species #:
# -------------------------------
# Drizzle                  1
# Rain                     2
# Ice Crystals             3
# Aggregates               4
# Wet Snow                 5
# Vertical Ice             6
# Low-Density Graupel      7   LD-Graupel
# High-Density Graupel     8   HD-Graupel
# Hail                     9
# Big Drops                10
#

csu_hidt = csu_fhc.csu_fhc_summer( dz=ref, zdr=zdr, rho=rho, kdp=kdp, use_temp=True, T=radar_T, band=rband, 
                                   verbose=False, use_trap=False,  method='hybrid ')

rdabufr = add_field_to_radar_object(csu_hidt, rdabufr, field_name='csu_hidt', dz_field='reflectivity')

# breakpoint()
#
#
# hcawinter = csu_fhc.run_winter(dz=dz,zdr=dr,kdp=kd,rho=rh,azimuths=azimuths,sn_thresh=snthresh,
#                               expected_ML=expected_ML,sn=sndat,T=radar_T,heights=rheights,nsect=nsect,
#                               scan_type = scan_type,verbose=verbose,use_temp = use_temp,
#                               band=rband, return_scores=return_scores)
# radar = add_field_to_radar_object(hcawinter, radar,field_name='HCA',dz_field='corrected_reflectivity')
#
#
#
#
# ================================================================================
#
#

qpe = pyart.retrieve.est_rain_rate_z(rdabufr, refl_field="reflectivity" )

rdabufr.add_field("qpe", qpe, replace_existing=False)

#
# kdp = pyart.retrieve.kdp_maesaka(rdabufr, method='cg', phidp_field = "differential_phase", verbose=True)
# kdp = pyart.retrieve.kdp_maesaka(rdabufr, method='cg', kdp_field="kdp", verbose=True)
# kdp = pyart.retrieve.kdp_maesaka(rdabufr, method='cg', verbose=True)
# pyart:sshca
# 0  None
# 1  DS, drysnow
# 2  CR, crystals
# 3  LR, light rain
# 4  GR, ld grapel
# 5  RN, rain, 
# 6  VI, verticle ice
# 7  WS, wet snow
# 8  MH, melting hail
# 9  IH/HDG, high-density grapel, 
#10  HeavyRain, HR

ss_hca = pyart.retrieve.hydroclass_semisupervised(rdabufr, temp_field="temperature" )

rdabufr.add_field("art_hca", ss_hca, replace_existing=False)

# pdb.set_trace()

# ==================== pyart output  ================================

# rdafname = "../data/" +  rda_site + "_hca_" + ftime_utc + ".nc"
# pyart.io.write_cfradial(rdafname, rdabufr,  format="NETCDF4", arm_time_variables= True )
# print(["Write rdafile: ", rdafname])


# ================================================================================= Grdding AXPT variants

gridding = True

if(gridding):
    
    var_list   = [ "reflectivity",  "differential_reflectivity", "specific_differential_phase",  
                       "qpe", "art_hca" ]
    xpt_gridding(rdabufr, xpar_radius, var_list)

    var_list   = [ "velocity" ]
    xpt_gridding(rdabufr_dop, xpar_radius, var_list)

# ========================================================================== Plotting RDA parameter figs



# print("")
# print("Plotting RDA-SECTION Figs")

# fig_path =  figroot+ "/" + rda_site + "/"
# os.system("mkdir " + fig_path)
# fig_prefix = rda_site  + "_HCA_SEC_" + ftime_utc
# plot_rda_sect(rdabufr, rdabufr_dop, fig_path+fig_prefix, [0, 30, 60, 90, 120, 150,  180, 210, 240,  270, 300, 330],  plot_radius=60  )

# breakpoint()


# ==========================================================================

print("")
print("Plotting RDA Parameters Figs")


os.system("mkdir " + figroot + "/" + rda_site)  # make fig dir

fig_path = figroot + "/" + rda_site + "/"
os.system("mkdir " + fig_path)

# ============================================================================

fig_prefix = rda_site  + "_RDA_parameter_" + ftime_utc

plot_rda_parameter(rdabufr, mom_log_list,  fig_path+fig_prefix, plot_radius, rng_list)

plot_rda_parameter(rdabufr_dop, mom_dop_list,  fig_path+fig_prefix, plot_radius, rng_list)
# =========================================================================== Plotting HCA figs

print("")
print("Plotting RDA-HCA Figs")

fig_path = figroot + "/" + rda_site + "/"
fig_prefix = rda_site  + "_HCA_" + ftime_utc

plot_hca_map(rdabufr_log, fig_path+fig_prefix, plot_radius, rng_list )


# ================= Radar Azimuth Angle ================= RDA Section figs
#  
#          North=0 
#
# WEST=270            East=90
#
#          South=180
#
# ========================================================



# pdb.set_trace()

