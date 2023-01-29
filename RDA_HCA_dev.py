# import matplotlib.pyplot as plt
# import pandas as pd
# from sklearn.covariance import log_likelihood
# from cinrad.visualize import PPI

import os

import numpy as np

import argparse

import pyart

import cinrad

from cinrad.io.export import standard_data_to_pyart

from csu_radartools import ( csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, csu_kdp, csu_misc, fundamentals )

# ====================================== local utils  =======================================

from define_hca_params import *

from utils_rda_hca import *
from utils_get_temp import *
from utils_plot_rda import *

import pdb
import warnings 

warnings.filterwarnings("ignore")

# ===========================================================================================
# ===========================================================================================
# ===========================================================================================

parser = argparse.ArgumentParser(description='manual to this script, input Standard SPOL-RADAR Level2 Dasta')

parser.add_argument('--rdafname', type=str, default = None)

args = parser.parse_args()

rdafname = args.rdafname

# ==============================================================================

# rdafname = "./data/Z9662/Z_RADR_I_Z9662_20220510203000_O_DOR_SAD_CAP_FMT.bin.bz2"
# fname = "./data/Z9200/Z_RADR_I_Z9200_20220511000000_O_DOR_SAD_CAP_FMT.bin.bz2"

absfname = os.popen("basename " + rdafname).readline().strip("\n")

datapath = os.popen("dirname " + rdafname).readline().strip("\n")

# pdb.set_trace()

rda_site = absfname[9:14]

ftime_utc = absfname[15:15+14]

print([rda_site, ftime_utc])

# datapath = "./data/Z9755/"
# fname = "Z_RADR_I_Z9755_20220511000000_O_DOR_SAD_CAP_FMT.bin.bz2"



# ===========================================================#! reading in RDA-Standard Level2 data by CINRAD
# # TODO : is there anything code to revise?

fbufr = cinrad.io.StandardData(rdafname)  

# pdb.set_trace()
# ===========================================================

scantime = fbufr.scantime.strftime("%Y-%m-%d %H:%M:%S")

rdabufr = standard_data_to_pyart(fbufr)

print(rdabufr.info())

print(rdabufr.elevation)

print(["nsweeps: ", rdabufr.nsweeps])

rda_lat = rdabufr.latitude['data'][0]
rda_lon = rdabufr.longitude['data'][0]

# pdb.set_trace()

# ================================================================================ #! setup sounding temperature
#  adding  model-forecast temperature   by chentao@cma.gov.cn 2023-01-25 12:30:16

match temp_source:
    case "sounding":
        sndbufr = setup_snd(snd_fname)
        # snddata = SkewT.Sounding(snd_fname).soundingdata       

        wh0 = np.where(np.isclose(np.abs(sndbufr['temp']), 0.0, atol=0.5))

        expected_ML = np.array(sndbufr['hght'])[wh0[0]][0]/1000.

        radar_T, radar_z = interpolate_sounding_to_radar(sndbufr, rdabufr)  # radar_T: interpolated sounding temperature, units in degC

    case "g3km":         # =========  using NWP temperature 

        model_fname = setup_g3km_model(ftime_utc)

        radar_T, radar_z = interpolate_temp_to_radar( rdabufr, g3km_dir + model_fname)  # radar_T: interpolated model temperature, units in degC

# pdb.set_trace()

# ================================================================================ 


rheight = radar_z/1000

rdabufr = add_field_to_radar_object(radar_T, rdabufr, field_name='temperature', units="degree", dz_field='reflectivity')

# pdb.set_trace()
# sn = pyart.retrieve.simple_moment_calculations.calculate_snr_from_reflectivity(rdabufr, refl_field='reflectivity',toa=15000.0)
# sndat= sn['data'][:]
#add SNR to the radar structure
# radar = add_field_to_radar_object(sndat, rdabufr, field_name='SN', dz_field='reflectivity')
# sndat =np.ma.masked_array(sndat)
# pdb.set_trace()
###If you use SNR, select an appropriate threshold for your data.

snthresh = -30

###Azimuths are important for the PPI data for the bright band detection.
azimuths = rdabufr.azimuth['data']

##The melting layer algorithm for PPIs takes different sectors. The easiest is to do 10º sectors, and pass nsect=36 (also
# the default). But you could change to any multiple of 360 (or the number of azimuths you have).
nsect = 10

#Radar variables. These should be just the data (not an object), and masked arrays
dz = np.ma.masked_array(rdabufr.fields['reflectivity']['data'])
dr = np.ma.masked_array(rdabufr.fields['differential_reflectivity']['data'])
kd = np.ma.masked_array(rdabufr.fields['specific_differential_phase']['data'])
rh = np.ma.masked_array(rdabufr.fields['cross_correlation_ratio']['data'])

# pdb.set_trace()

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

csu_hidt = csu_fhc.csu_fhc_summer( dz=dz, zdr=dr, rho=rh, kdp=kd, use_temp=True, T=radar_T, band=rband, 
                                   verbose=False, use_trap=False, method='hybrid ')


rdabufr = add_field_to_radar_object(csu_hidt, rdabufr, field_name='csu_hidt', dz_field='reflectivity')

# hcawinter = csu_fhc.run_winter(dz=dz,zdr=dr,kdp=kd,rho=rh,azimuths=azimuths,sn_thresh=snthresh,
#                               expected_ML=expected_ML,sn=sndat,T=radar_T,heights=rheights,nsect=nsect,
#                               scan_type = scan_type,verbose=verbose,use_temp = use_temp,
#                               band=rband, return_scores=return_scores)
# radar = add_field_to_radar_object(hcawinter, radar,field_name='HCA',dz_field='corrected_reflectivity')
#
#
###Let's see what we got.
# plt.show()
#
# ================================================================================

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

art_hca = pyart.retrieve.hydroclass_semisupervised(rdabufr, temp_field="temperature" )

rdabufr.add_field("art_hca", art_hca, replace_existing=False)

# pdb.set_trace()

# ====================#! pyart output in RDA-Standard format  ===========================================

rdafname = rdagrid_dir +  rda_site + "_hca_" + ftime_utc + ".nc"

pyart.io.write_cfradial(rdafname, rdabufr, format="NETCDF4", arm_time_variables= True )

print(["Write rdafile: ", rdafname])

# ====================#! pyart output in gridded RDA format with [HCA, QPE] parameters   ================

rdagrid = pyart.map.grid_from_radars(rdabufr, grid_shape=(21, 1000, 1000),
                                     grid_limits=((0, 18000), (-250000.0, 250000.0), (-250000.0, 250000.0)),
                                     fields=[ "reflectivity",  "velocity",  "differential_reflectivity", 
                                              "specific_differential_phase", "csu_hidt", "art_hca" ])

rdagrid_fname = rdagrid_dir + rda_site + "_rdagrid_"  + ftime_utc +  ".nc"

pyart.io.write_grid(rdagrid_fname, rdagrid, format="NETCDF4", write_point_lon_lat_alt=True, arm_alt_lat_lon_variables=False )

print(["Write grid_rdafile: ", rdafname])

# ============================================================================== debug code
# perform Cartesian mapping, limit to the reflectivity field.
# fig = plt.figure()
# ax = fig.add_subplot(121)
# ptx = ax.imshow(rdagrid.fields['reflectivity']['data'].data[0, :, :], origin='lower',
#                 vmin=-5, vmax=70,  cmap=cmaps.radar)               
# plt.colorbar(ptx, ax=ax)
# ax = fig.add_subplot(122)
# ptx = ax.imshow(rdagrid.fields['velocity']['data'].data[0, :, :], origin='lower',
#                 vmin=-30, vmax=30,  cmap=plt.cm.bwr)               
# plt.colorbar(ptx, ax=ax)
# plt.show()
# pdb.set_trace()

# =========================================================================== HCA figs

print("")
print("Plotting RDA-HCA Figs")

fig_path = figroot + "/" + rda_site + "/" + ftime_utc[0:8] + "/"

os.system("mkdir " + fig_path)

fig_prefix = rda_site  + "_HCA_" + ftime_utc

# plot_hca(rdabufr, fig_path +fig_prefix)

plot_hca_map(rdabufr, fig_path +fig_prefix, plot_radius)

# ========================================================================= RDA parameter figs

print("")
print("Plotting RDA Parameters Figs")

fig_path = figroot + "/" + rda_site + "/"+ ftime_utc[0:8] + "/"
os.system("mkdir " + fig_path)

fig_prefix = rda_site  + "_RDA_parameter_" + ftime_utc

plot_rda_parameter(rdabufr, fig_path+fig_prefix, plot_radius )

# ====================== Radar Azimuth Angle ================= RDA Section figs
#  
#          North=0 
#
# WEST=270            East=90
#
#          South=180
#
#

print("")
print("Plotting RDA-SECTION Figs")

fig_path =  figroot+ "/" + rda_site + "/" + ftime_utc[0:8] + "/"
fig_prefix = rda_site  + "_HCA_SEC_" + ftime_utc
plot_rda_sect(rdabufr, fig_path+fig_prefix, [0, 30, 60,  90, 120, 150,  180, 210, 240,  270, 300, 330], plot_radius )


# ====================== release memory ====================================

del(fbufr)
del(rdabufr)
del(rdagrid)

# pdb.set_trace()


