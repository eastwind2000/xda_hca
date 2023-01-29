import matplotlib.pyplot as plt

import matplotlib.cm as cm

from  matplotlib.colors import ListedColormap

import os

import numpy as np

# import pyart
# import cinrad
# from cinrad.visualize import PPI
# from cinrad.io.export import standard_data_to_pyart

import cmaps
import pandas as pd
import xarray as xr

from scipy.interpolate import griddata

# import warnings 
import re

import pdb

# ====================================================================================
# Z_RADR_I_Z9371_20210720000000_O_DOR_SAD_CAP_FMT.bin.bz2
# 广州  Z9200
# 深圳  Z9755
# 阳江  Z9662

# warnings.filterwarnings("ignore", category=DeprecationWarning)

# datadir = "./data/Z9662/"

def sweep_preproc(rda_site, rdavar, out_prefix):

    sweep_dir = "./sweep_data/"    

    flist = os.popen("ls -1 " + sweep_dir + rdavar + "_SWEEP_00_*" + rda_site + "*2022051[0-1]* | sort "  ).readlines()

    # fname = datadir + "Z_RADR_I_Z9371_20210720000000_O_DOR_SAD_CAP_FMT.bin.bz2"

    # ./sweep_data/REF_SWEEP_00_Z9200_20220510000001.nc

    # focus on ZhongShan
    focus_lat, focus_lon = 22.3, 113.5

    print(flist)

    # pdb.set_trace()

    htbufr = []

    for txtbufr in flist[:]:
        
        fname_sweep00 = txtbufr.strip("\n")
        print("")
        print("")
        print("")
    
        for isweep in ["00", "02", "04", "05", "06", "07", "08", "09", "10" ]:

            fname = sweep_dir + rdavar + "_SWEEP_" + isweep + fname_sweep00[25:]
            
            print(fname)
            # pdb.set_trace()    

            rdabufr = xr.open_dataset(fname, engine="netcdf4")

            site_lat = rdabufr.attrs["site_latitude"]
            site_lon = rdabufr.attrs["site_longitude"]

            # ================= perform gridding by scipy.interploation  ==============================

            rlat    = rdabufr["latitude"].data
            rlon    = rdabufr["longitude"].data
            rheight = rdabufr["height"].data

            points = np.vstack( (rlat.flatten(), rlon.flatten()) ).transpose()
            
            refbufr  = rdabufr[rdavar].data           

            dd = 2.0

            lat2d, lon2d = np.mgrid[ site_lat-dd:site_lat+dd:400j,  site_lon-dd:site_lon+dd:400j  ]
            
            grid_height  = griddata(points, rheight.flatten(), (lat2d, lon2d), method='nearest' )
            
            grid_ref     = griddata(points, refbufr.flatten(), (lat2d, lon2d), method='nearest' )
            
            # =========================== get radar obs in focus-district ==============================

            ddx = 0.02

            cond_lat = (lat2d >= focus_lat -ddx) &   (lat2d <= focus_lat + ddx) 
            cond_lon = (lon2d >= focus_lon -ddx) &   (lon2d <= focus_lon + ddx)  
            cond_area = cond_lat & cond_lon
            
            seclat     = lat2d[cond_area]
            seclon     = lon2d[cond_area]
            secheight  = grid_height[cond_area]

            secvar     = grid_ref[cond_area]

            # print(["seclat: ",  isweep, seclat])
            # print(["seclon: ",  isweep, seclon])
            # print(["sech:   ",  isweep, secheight])
            print([fname, rdavar,  isweep, secvar])
            # 
            #
            # testing code 
            # plt.contourf(grid_var[::-1, :]);plt.colorbar();plt.show()
            #

            # ===================== save inermediate data ==========================
            htrec = [fname, isweep, seclat, seclon, secheight, secvar]        

            # print(htrec)

            htbufr.append(htrec)

        # pdb.set_trace()
    
    
        # ===================== save inermediate data ==========================

    htdata = pd.DataFrame(htbufr)
    htdata.columns=["fname", "isweep", "seclat", "seclon", "secheight", "sec_" + rdavar]

    htdata.to_csv(out_prefix + ".csv", index=False)
    htdata.to_pickle(out_prefix + ".pickle" )
    print(htdata)

    # pdb.set_trace()

    return  

#
# ================================================================================
# ================================================================================
# ================================================================================

sweep_preproc("Z9200", "KDP", "radar_Z9200_ht_KDP_zhongshan")
sweep_preproc("Z9200", "ZDR", "radar_Z9200_ht_ZDR_zhongshan")
sweep_preproc("Z9200", "REF", "radar_Z9200_ht_REF_zhongshan")

# pdb.set_trace()
# ================================================================================
# ================================================================================


radar_sid = "Z9200"

nlev = 9

# ===========================================================  setting ZDR

dfbufr = pd.read_pickle("radar_" + radar_sid + "_ht_ZDR_zhongshan.pickle")
nrec = len(dfbufr)
nt = int(nrec/nlev)
databufr = dfbufr["sec_ZDR"]

nsample = len(dfbufr["sec_ZDR"].iloc[0])

secheight = [ np.mean(x)  for x in  dfbufr["secheight"].iloc[0:9] ]

fbufr  = np.zeros((nrec, nsample), dtype=float)
for irec in np.arange(nrec):
    fbufr[irec, :] =  databufr[irec]

fdata = fbufr.reshape((nt, nlev, nsample))

# fmean = np.transpose( np.nanmean(fdata, axis=2) )

zdr_htime = np.transpose( np.nanmean(fdata, axis=2) )
zdr_std= np.transpose( np.nanstd(fdata, axis=2) ) 


# =============================== for testting ===============================

# xx = dfbufr.iloc[2196:2196+9]["sec_ZDR"]
# pdb.set_trace()

# ============================================================================= setting REF


dfbufr = pd.read_pickle("radar_" + radar_sid + "_ht_REF_zhongshan.pickle")
nrec = len(dfbufr)
nt = int(nrec/nlev)

databufr = dfbufr["sec_REF"]

nsample = len( dfbufr["sec_REF"].iloc[0] )

secheight = [ np.mean(x)  for x in  dfbufr["secheight"].iloc[0:nlev] ]

fbufr  = np.zeros((nrec, nsample), dtype=float)

for irec in np.arange(nrec):

    fbufr[irec, :] =  databufr[irec]

fdata = fbufr.reshape((nt, nlev, nsample))

# fmean = np.transpose( np.nanmean(fdata, axis=2) )

ref_htime = np.transpose( np.nanmax(fdata, axis=2) )
ref_std= np.transpose( np.nanstd(fdata, axis=2) ) 


# ============================================================  setting REF

dfbufr = pd.read_pickle("radar_" + radar_sid + "_ht_KDP_zhongshan.pickle")
nrec = len(dfbufr)
nt = int(nrec/nlev)
databufr = dfbufr["sec_KDP"]

nsample = len(dfbufr["sec_KDP"].iloc[0])

secheight = [ np.mean(x)  for x in  dfbufr["secheight"].iloc[0:9] ]

fbufr  = np.zeros((nrec, nsample), dtype=float)
for irec in np.arange(nrec):
    fbufr[irec, :] =  databufr[irec]

fdata = fbufr.reshape((nt, nlev, nsample))

# fmean = np.transpose( np.nanmean(fdata, axis=2) )

kdp_htime = np.transpose( np.nanmax(fdata, axis=2) )
kdp_std= np.transpose( np.nanstd(fdata, axis=2) ) 


# pdb.set_trace()

# ================================================================================
# nrec = len(dfbufr)
# nt = int(nrec/nlev)
# fbufr  = np.zeros((nrec, 9), dtype=float)
# for irec in np.arange(nrec):
#     sbufr = re.sub('[\[\]]', '', dfbufr["secref"].iloc[irec])
#     print(sbufr)
#     fbufr[irec, :] = [ float(x) for x in sbufr.split()  ]
# # pdb.set_trace()
# fdata = fbufr.reshape((nt, nlev, 9))
# fmean = np.nanmean(fdata, axis=2)
# ========================================================================================

xtime =  np.array(dfbufr["fname"].iloc[0::9])

# ========================================================================================

fig, axlist = plt.subplots(3, 2, figsize=(18, 10))

ax = axlist.flatten()

reflevels = [0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75]
# plot_obj =  ax[0].imshow(fmean[:, :].transpose(),  cmap=cmaps.radar )
plot_obj =  ax[0].contourf(ref_htime[:, :], levels=reflevels,  cmap=cmaps.radar, extend="both" )

xticks =  np.arange(0, 360, step=30)
ax[0].set_xticks( xticks   )
ax[0].set_xticklabels(  [  txtbufr[36:42] for txtbufr in xtime[xticks]     ]  )

ax[0].set_yticks( [0, 1, 2, 3, 4, 5, 6, 7, 8]  )
ax[0].set_yticklabels( [ str(round(x, 1)) for x in secheight    ]       )

ax[0].set_xlabel("Time")
ax[0].set_ylabel("REF \dBZ")

ax[0].set_xlim(120, 300)
plt.colorbar(plot_obj, ax=ax[0], shrink=1.0)

plot_ref_std =  ax[1].contourf(ref_std,  cmap=cmaps.radar, extend="both" )
ax[1].set_xlim(120, 300)
ax[1].set_ylabel("REF STD")
plt.colorbar(plot_ref_std, ax=ax[1], shrink=1.0)


# =======================================================================================

kdplevels = [-0.4,  -0.2,  -0.1,  0.1,  0.15,  0.22,  0.33,  0.5, 0.75,  1.1,  1.7,  2.4,  3.1,  7,  20]

kdpcolors = [ "#00efef", "#00979a", "#b4b4b4", "#b4b4b4", "#00c027", "#00e80a", "#24ff24", "#ffff1e", 
              "#ffff1e", "#ffe600", "#ffbc00", "#ff9800", "#ff5e00", "#f20f00", "#bb003a", "#ff00ff" ]

kdpcmap = ListedColormap( kdpcolors )

plot_obj =  ax[2].contourf(kdp_htime[:, :], levels=kdplevels, colors=kdpcolors,  extend="both")

xticks =  np.arange(0, 360, step=30)
ax[2].set_xticks( xticks   )
ax[2].set_xticklabels(  [  txtbufr[36:42] for txtbufr in xtime[xticks]     ]  )

ax[2].set_yticks( [0, 1, 2, 3, 4, 5, 6, 7, 8]  )
ax[2].set_yticklabels( [ str(round(x, 1)) for x in secheight    ]       )

ax[2].set_xlabel("Time")
ax[2].set_ylabel("KDP ")

ax[2].set_xlim(120, 300)
plt.colorbar(plot_obj, ax=ax[2], shrink=1.0)

plot_std =  ax[3].contourf(kdp_std,  cmap=cmaps.radar, extend="both" )
ax[3].set_xlim(120, 300)
ax[3].set_ylabel("KDP STD")
plt.colorbar(plot_std, ax=ax[3], shrink=1.0)


# ======================================================================================== draw ZDR

zdrlevels = [ -3, -2, -1, 0, 0.2, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5  ]

zdrcolors = [ "#6b6d6a", "#949693", "#cbcbd0", "#dcf3de", "#01c21f", "#00eb0b", "#22fd22", 
              "#fffe19", "#ffe601", "#ffbc02", "#ff9a00", "#fe5d00", "#f50e00", "#bd003a", "#ff00fd" ]

zdrcmap = ListedColormap( zdrcolors )

plot_obj =  ax[4].contourf(zdr_htime[:, :], levels=zdrlevels, colors=zdrcolors,  extend="both")

xticks =  np.arange(0, 360, step=30)
ax[4].set_xticks( xticks   )
ax[4].set_xticklabels(  [  txtbufr[36:42] for txtbufr in xtime[xticks]     ]  )

ax[4].set_yticks( [0, 1, 2, 3, 4, 5, 6, 7, 8]  )
ax[4].set_yticklabels( [ str(round(x, 1)) for x in secheight    ]       )

ax[4].set_xlabel("Time")
ax[4].set_ylabel("ZDR ")

ax[4].set_xlim(120, 300)
plt.colorbar(plot_obj, ax=ax[4], shrink=1.0)


plot_std =  ax[5].contourf(zdr_std,  cmap=cmaps.radar, extend="both" )
ax[5].set_xlim(120, 300)
ax[5].set_ylabel("ZDR STD")
plt.colorbar(plot_std, ax=ax[5], shrink=1.0)




# ========================================================================================
plt.tight_layout()
plt.savefig( "radar_" + radar_sid + "_htime_zhongshan.eps")
plt.savefig( "radar_" + radar_sid + "_htime_zhongshan.png", dpi=600)

plt.show()



pdb.set_trace()


