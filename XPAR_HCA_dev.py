import matplotlib.pyplot as plt
import os
import numpy as np
import argparse
import pyart
import cinrad

# from cinrad.visualize import PPI
# from cinrad.io.export import standard_data_to_pyart

from csu_radartools import ( csu_fhc, csu_liquid_ice_mass, csu_blended_rain, csu_dsd, csu_kdp, csu_misc, fundamentals )

import warnings 

import pandas as pd

import pdb

from sklearn.covariance import log_likelihood

from util_xpar2pyart import xpar_data_to_pyart

from utils_XPAR_HCA  import *


#
# dtype([('header', 
#         [ ('filetype', '<u2'), 
#           ('longitude', '<i2'), 
#           ('latitude', '<i2'), 
#           ('v_beam_width', 'u1'), 
#           ('h_beam_width', 'u1'), 
#           ('height', '<u2'), 
#           ('reserved', 'S22')]), 
#           ('data', 
#            [('radial_time', '<u4'), 
#             ('radial_date', '<u4'), 
#             ('unambiguous_dist', '<u2'), 
#             ('azimuth', '<u2'), 
#             ('az_num', '<u2'), 
#             ('elevation', '<u2'), 
#             ('el_num', '<u2'), 
#             ('radial_index', '<u2'), 
#             ('first_gate_dist', '<u2'), 
#             ('gate_length', '<u2'), 
#             ('gate_num', '<u2'), 
#             ('reserved', 'S26'), 
#             ('tref_pt', '<u2'), 
#             ('ref_pt', '<u2'), 
#             ('vel_pt', '<u2'), 
#             ('sw_pt', '<u2'), 
#             ('zdr_pt', '<u2'), 
#             ('phi_pt', '<u2'), 
#             ('kdp_pt', '<u2'), 
#             ('rho_pt', '<u2'), 
#             ('tref', 'u1', (2000,)), 
#             ('ref', 'u1', (2000,)), 
#             ('vel', 'u1', (2000,)), 
#             ('sw', 'u1', (2000,)), 
#             ('zdr', 'u1', (2000,)), 
#             ('phi', 'u1', (2000,)), 
#             ('kdp', 'u1', (2000,)), 
#             ('rho', 'u1', (2000,)), 
#             ('res2', 'S4')])])
#
# fbufr = cinrad.io.StandardData(rdafname)
# rdabufr = standard_data_to_pyart(fbufr, radius=460)
# pdb.set_trace()
#

# ===========================================================================================

def setup_snd():

    """
    ! waiting for fix , sounding-radar coehence
    """

    snd_fname =  "../data/snd59280_2022051108.dat"

    txtbufr = pd.read_csv(snd_fname,  header=None, skiprows=1, sep="\\s+")
    
    txtbufr.columns=("pressure", "hght", "temp", "dwpt", "wind_direction", "wind_speed")

    txtbufr['hght'] = txtbufr['hght']*10

    # cond_qc

    print("")
    print(["sounding data: ", snd_fname])
    print(txtbufr)
    
    return txtbufr

# ===========================================================================================
# ===========================================================================================
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

rband = "X"

#
# rdafname = "../data/XPAR_L2_scsmex2022/ZZH01/Z_RADR_I_ZZH01_20220511000023_O_DOR_DXK_CAR.bin.bz2"
# rdafname = "../data/T_RADR_I_ZG040_20220607003023_O_DOR_AXPT0364_CRA_FMT.bin.bz2"
# rband = "S"
# rdafname = "../data/Z9200/Z_RADR_I_Z9200_20220513224800_O_DOR_SAD_CAP_FMT.bin.bz2"
#
# ===============================================================================
#
#

absfname = os.popen("basename " + rdafname).readline().strip("\n")

datapath = os.popen("dirname " + rdafname).readline().strip("\n")

rda_site =absfname[9:14]

ftime_utc = absfname[15:15+14]

print([rda_site, ftime_utc])

os.system("mkdir " + figroot + "/" + rda_site)  # make fig dir


#! ===================================================== Reading in XPAR data by CINRAD and transform to ART radar obj


fbufr = cinrad.io.PhasedArrayData(rdafname)

rdabufr = xpar_data_to_pyart(fbufr, radius=50)

# pdb.set_trace()  
# ==========================================================

print(rdabufr.info())

print(rdabufr.elevation)

print(["nsweeps: ", rdabufr.nsweeps])

rda_lat = rdabufr.latitude['data'][0]
rda_lon = rdabufr.longitude['data'][0]

# pdb.set_trace()

# ========================================================= setup sounding temperature
# TODO : need to revise temerature algorithm

sndbufr = setup_snd()

# snddata = SkewT.Sounding(snd_fname).soundingdata

wh0 = np.where(np.isclose(np.abs(sndbufr['temp']), 0.0, atol=1.0))
expected_ML = np.array(sndbufr['hght'])[wh0[0]][0]/1000.

# pdb.set_trace()

####Interpolate sounding to the radar heights
radar_T, radar_z = interpolate_sounding_to_radar(sndbufr, rdabufr)  # 



rheight = radar_z/1000.0                     #!  unit trans

rdabufr = add_field_to_radar_object( radar_T, rdabufr, field_name='temperature', units="degree", 
                                     dz_field='reflectivity')

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

ref = np.ma.masked_array(rdabufr.fields['reflectivity']['data'])
vel = np.ma.masked_array(rdabufr.fields['velocity']['data'])

zdr = np.ma.masked_array(rdabufr.fields['differential_reflectivity']['data'])
kdp = np.ma.masked_array(rdabufr.fields['specific_differential_phase']['data'])
rho = np.ma.masked_array(rdabufr.fields['cross_correlation_ratio']['data'])


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
#


csu_hidt = csu_fhc.csu_fhc_summer( dz=ref, zdr=zdr, rho=rho, kdp=kdp, use_temp=True, T=radar_T, band=rband, 
                                   verbose=False, use_trap=False,  method='hybrid ')

rdabufr = add_field_to_radar_object(csu_hidt, rdabufr, field_name='csu_hidt', dz_field='reflectivity')

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
###Let's see what we got.
# plt.show()
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

art_hca = pyart.retrieve.hydroclass_semisupervised(rdabufr, temp_field="temperature" )

rdabufr.add_field("art_hca", art_hca, replace_existing=False)

# pdb.set_trace()

# ==================== pyart output  ================================

rdafname = "../data/" +  rda_site + "_hca_" + ftime_utc + ".nc"
pyart.io.write_cfradial(rdafname, rdabufr,  format="NETCDF4", arm_time_variables= True )
print(["Write rdafile: ", rdafname])


# ============================================================================
# perform Cartesian mapping, limit to the reflectivity field.

rdagrid = pyart.map.grid_from_radars(rdabufr, grid_shape=(6, 800, 800),
                                     grid_limits=((0, 3000), (-200000.0, 200000.0), (-200000.0, 200000.0)),
                                     fields=[ "reflectivity",  "velocity",  "differential_reflectivity", 
                                              "specific_differential_phase", "csu_hidt", "art_hca" ])

rdagrid_fname = "../data/" + rda_site + "_rdagrid_"  + ftime_utc +  ".nc"

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

# =========================================================================== Plotting HCA figs

plot_radius = 45  #! units in km

print("")
print("Plotting RDA-HCA Figs")

fig_path = figroot + "/" + rda_site + "/"
fig_prefix = rda_site  + "_HCA_" + ftime_utc

# plot_hca(rdabufr, fig_path +fig_prefix)
plot_hca_map(rdabufr, fig_path+fig_prefix, plot_radius )

# ========================================================================== Plotting RDA parameter figs

print("")
print("Plotting RDA Parameters Figs")


fig_path = figroot + "/" + rda_site + "/"
os.system("mkdir " + fig_path)
fig_prefix = rda_site  + "_RDA_parameter_" + ftime_utc

plot_rda_parameter(rdabufr, fig_path+fig_prefix, plot_radius )



# ================= Radar Azimuth Angle ================= RDA Section figs
#  
#          North=0 
#
# WEST=270            East=90
#
#          South=180
#
# ========================================================

print("")
print("Plotting RDA-SECTION Figs")

fig_path =  figroot+ "/" + rda_site + "/"
os.system("mkdir " + fig_path)
fig_prefix = rda_site  + "_HCA_SEC_" + ftime_utc
plot_rda_sect(rdabufr, fig_path+fig_prefix, [0, 30, 60,  90, 120, 150,  180, 210, 240,  270, 300, 330], 40  )


# pdb.set_trace()


