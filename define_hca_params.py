
import os

workdir = "/home/lse/scsmex2022/radarpro/"

# rda_dir = "/home/lse/scsmex2022/data/radar_L2_scsmex2022/"
rda_dir = "/mnt/i/fcst2023/coldwave/radar/"

# 2023014 coldwave rain->sleet->snow Case / Nanjing-SPOL
# setting up  parallel processes using pool
# rdasite_list = ["Z9250"]                                           #! 2023/01/14 Nanjing rain, sleet to snsow, precipitation type transition

# rdasite_list= [ "Z9758", "Z9762", "Z9200", "Z9662",  "Z9755"  ]  #! 2022 SCSMEX GBA
# rdasite_list= [ "Z9518", "Z9527", "Z9516" ]                      #! 20220720 JiangSu Province Tornadoes Case
# rdasite_list= [ "Z9370" ]                                        #! 2022072200-2300UTC Henan-Shangqiu Tornadoes Case
# rdasite_list= [ "Z9710"  ]                                       #! 20210812 SuiZhou, HuBei Extreme Precipitaion Case
# rdasite_list= [ "Z9971", "Z9970"   ]                             #! 20220818 Qinghai DaTong Extreme precipitation Case
# rdasite_list= [ "Z9371", "Z9379", "Z9375"   ]                    #! 20210719-20 Zhengzhou 7.20 rainstorm

rdasite_list = ["Z9010"]                                           #! 2023/01/12 Beijing snsow,sleet, rain, complex  precipitation type transition

rband = "S"                                   #! "S", "X

nproc = 4

# ====================================================================== fig settings

# figroot = workdir + "./remote_figs/"
figroot = "/mnt/i/fcst2023/coldwave/figs/"
# ====================================================================== fig settings

fig_dpi = 400

plot_radius = 230                             #! fig producsts range, unit in km

shp_path = "/home/lse/res/shp/"

shp_list = ["bou1_4l.shp", "continents_lines.shp", "hyd1_4l.shp", "province_l.shp", "City.shp", "river_c1.shp"]

# ====================================================================== temperature settings

temp_source = "g3km"
# model_fname = model_dir + "./20230114/Z_NAFP_C_BABJ_20230114000000_P_NWPC-GRAPES-3KM-ORIG-00600.grb2"
g3km_dir = "/mnt/i/fcst2023/coldwave/g3km/"

# temp_source = "sounding"
# snd_fname =  "../data/snd52866_2022081712UTC.dat"
# snd_fname =  "/mnt/f/fcst2022/chaba2022/snd57687_2022070400UTC.dat"
# snd_fname =  "/mnt/f/fcst2022/Tornado_IOP2022/snd57083_2022072200UTC.dat"
# snd_fname =  "/mnt/f/fcst2022/Tornado_IOP/snd58027_2022072000UTC.dat"   
# snd_fname =  "/mnt/i/henan2021/data/micaps/snd57083_2021071912utc.dat"    
snd_fname =  "/mnt/i/fcst2023/coldwave/snd58238_2023011412UTC.dat"

# ===============================================================================

rdagrid_dir =  workdir + "./data/interm_nc/"

# ===============================================================================

