import matplotlib.pyplot as plt

import os

import numpy as np

import pyart

import cinrad

from cinrad.visualize import PPI

from cinrad.io.export import standard_data_to_pyart

import cartopy.crs as ccrs

import cmaps

# import warnings 

import pdb

# import warnings

# ====================================================================================
# Z_RADR_I_Z9371_20210720000000_O_DOR_SAD_CAP_FMT.bin.bz2

def radar_to_sweep_netcdf(fname, sweep_dir): 
    
    """ Transeform CINRAD LevelII Standard basefile  to single sweep netcdf. by chentao@cma.gov.cn """
    
    base_fname = os.path.basename(fname)
    
    # base_fname[0:]
    
    print(base_fname[9:29])   
    
    # pdb.set_trace()

    f = cinrad.io.StandardData(fname)

    print( [f.name,  f.stationlon] )

    cent_lat = f.stationlat

    cent_lon = f.stationlon

    dd = 2 

    for k in np.arange(11):        
        print([k, f.el[k],  f.available_product(k)])    
       
    print( ["tilt: ", f.el] )

    print( ["REF: ", f.available_tilt("REF")] )
    
    print( ["VEL: ", f.available_tilt("VEL")] )

    print( ["ZDR: ", f.available_tilt("ZDR")] )

    print( ["KDP: ", f.available_tilt("KDP")] )
    
    print( ["RHO: ", f.available_tilt("RHO")] )
    
    print( ["PHI: ", f.available_tilt("PHI")] )  
    
    # ======================================================
    
    for varname in ["REF", "VEL", "ZDR", "KDP"]:
            
        for k in f.available_tilt(varname):      
            
            sweepbufr = f.get_data(k, 460, varname)        
            
            sweepfname  = sweep_dir + varname  + "_SWEEP_" + str(k).zfill(2) + "_" +  base_fname[9:29]  + ".nc"        
            
            print(sweepfname)
            
            sweepbufr.to_netcdf(sweepfname)   

    # testing plot
    
    # vel_fig = PPI(vel, style="white", extent=[cent_lon-dd, cent_lon+dd, cent_lat-dd, cent_lat+dd])  

    # vel_fig.plot_range_rings()

    # plt.show()

    # pdb.set_trace()
    
    return


# ===========================================================================================
# ===========================================================================================
#
# 洛阳	Z9379	洛阳	河南	SB	112	27	7	34	38	34	244
# 南阳	Z9377	南阳	河南	SB	112	29	34	33	1	15	231.4
# 濮阳	Z9393	濮阳	河南	SB	114	57	28	35	47	49	145
# 三门峡	Z9398	三门峡	河南	SB	111	7	44	34	40	44	760.5q
# 商丘	Z9370	商丘	河南	SB	115	37	47	34	24	25	141.2
# 郑州	Z9371	郑州	河南	SA	113	41	49	34	42	14	202
# 驻马店	Z9396	驻马店	河南	SB	114	1	11	33	0	35	189.1
# 信阳	Z9376	信阳	河南	SA	114	4	11	32	7	23	150
# 平顶山	Z9375	平顶山	河南	SA	113	6	57	33	46	14	190.2
# 广州  Z9200
# 深圳  Z9755
# 阳江  Z9662

# warnings.filterwarnings("ignore", category=DeprecationWarning)

datadir = "./data/Z9662/"

sweep_dir = "./sweep_data/"    

flist = os.popen("ls -1 " + datadir + "*202205* | cut -c 14-" ).readlines()

# fname = datadir + "Z_RADR_I_Z9371_20210720000000_O_DOR_SAD_CAP_FMT.bin.bz2"

for txtbufr in flist:
    
    fname = datadir + txtbufr.strip("\n")
    
    radar_to_sweep_netcdf( fname, sweep_dir )



