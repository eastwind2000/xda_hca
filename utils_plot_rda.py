
import  os
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm

import pyart
import cartopy.crs as ccrs

# import xarray as xr
# from scipy import interpolate

import cmaps
import warnings 

from define_hca_params import (fig_dpi, shp_path, shp_list)

import pdb
import pycwr

# ========================================================================================

warnings.filterwarnings("ignore")


"""

cmap="pyart_NWSVel"


https://arm-doe.github.io/pyart/examples/plotting/plot_choose_a_colormap.html

"""

## setup the summre colors.

hid_colors =  [ 'White',      
                'LightBlue',   # light rain
                'MediumBlue',  # rain
                'Darkorange',  # ice crystals
                'LightPink',   # aggregates
                'Cyan',        # wetsnow
                'DarkGray',    # verticle ice
                'Lime',        # LD graupel
                'green',       # HD graupel 
                'Yellow',      # hail  # red
                'Fuchsia'      # BD 
                ]
        
cmaphid = colors.ListedColormap(hid_colors)

sshca_colors = [ 'white', 	     # 0  None
                 'LightPink',    # 1  dyrsnow
                 'Darkorange',   # 2  CR
                 'LightBlue',    # 3  lightrain
                 'Lime',         # 4  LD graupel
                 'MediumBlue',   # 5  Rain
                 'DarkGray',     # 6  Verticle Ice
                 'Cyan',         # 7  wet snow
                 'Yellow',       # 8  Melting Hail  # red
                 'green',        # 9  HD Graupel 
                 'Fuchsia'       # 10 Heavy Rain
                ]
                

cmap_sshca = colors.ListedColormap(sshca_colors)

##Setup the winter colors.
hid_colors_winter = ['White','Orange', 'Purple', 'Fuchsia', 'Pink', 'Cyan', 'LightBlue', 'Blue']#, 'Lime', 'Yellow', 'Red', 'Fuchsia']
cmaphidwinter = colors.ListedColormap(hid_colors_winter)

kdplevels = np.array([-0.4,  -0.2,  -0.1,  0.1,  0.15,  0.22,  0.33,  0.5, 0.75,  1.1,  1.7,  2.4,  3.1,  7,  20])
kdpcolors = [ "#00efef", "#00979a", "#b4b4b4", "#b4b4b4", "#00c027", "#00e80a", "#24ff24", "#ffff1e", 
              "#ffff1e", "#ffe600", "#ffbc00", "#ff9800", "#ff5e00", "#f20f00", "#bb003a", "#ff00ff" ]
cmapkdp =   colors.ListedColormap(kdpcolors)

zdrlevels = np.array([ -3, -2, -1, 0, 0.2, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5  ])
zdrcolors = [ "#6b6d6a", "#949693", "#cbcbd0", "#dcf3de", "#01c21f", "#00eb0b", "#22fd22", 
              "#fffe19", "#ffe601", "#ffbc02", "#ff9a00", "#fe5d00", "#f50e00", "#bd003a", "#ff00fd" ]
cmapzdr =  colors.ListedColormap(zdrcolors)

vellevels = np.array([-32, -27, -20, -15, -10, -5, -1, 1, 5, 10, 15, 20, 27, 32])
velcolors = [ "#7be3ff", "#00e3ff", "#00b2b5", "#00ff00", "#00c700", "#008200", "#ffffff", 
              "#ff0000", "#ff5a5a", "#ffb5b5", "#ff7d00", "#ffd300", "#ffff00", "#7c007c"   ]
cmapvel = colors.ListedColormap(velcolors)

# cmapvel= "pyart_NWSVel"

# =================================================================================================

def plot_rda_parameter(rdabufr, figprefix, plot_radius):   

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n") 

    angle_list = rdabufr.fixed_angle['data']

    rng_list = np.arange(50, plot_radius, step=50 )

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]
    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)       

    # pdb.set_trace()
    dd = (180/np.pi) *  plot_radius*1e3/6371e3

    xlim = (rda_lon-dd, rda_lon+dd)
    ylim = (rda_lat-dd, rda_lat+dd)

    display = pyart.graph.RadarMapDisplay(rdabufr)      

    fig= plt.figure(figsize=(16, 16))
    
    #!   ======================================================= Plot vel on vel levels 
    for isweep in [1, 3, 4, 5, 6, 7, 8, 9, 10]:   
                
        ppi_angle  = str(round(angle_list[isweep], 2)) 

        velnorm = BoundaryNorm(vellevels, ncolors=len(velcolors), clip=True)

        # display.plot_ppi_map('velocity', sweep=isweep, vmin=-30, vmax=30,  cmap="pyart_NWSVel", 
        #                         title= "VEL " + ppi_angle + '\n ' + scantime,
        #                         min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
        #                         resolution='10m',  projection=projection, fig=fig, colorbar_flag=False,                                
        #                         lat_0=rda_lat, lon_0=rda_lon)


        display.plot_ppi_map('velocity', sweep=isweep, norm=velnorm,  cmap=cmapvel, 
                                title= "VEL " + ppi_angle + '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  projection=projection, fig=fig, colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)
    

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )

        display.plot_colorbar(shrink=0.7, ticks=vellevels )
        # display.set_limits(xlim, ylim)

        print( [ "Fig:  ", figprefix + "_vel_swp" + str(isweep).zfill(2) +  ".png"] )
        plt.savefig(figprefix + "_vel_swp" + str(isweep).zfill(2) + ".png", 
                     dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        
        # plt.show()

        plt.clf()   # clear all axis on this fig

    #! ============================================================ plot REF on REF levels
    
    for isweep in [0, 2, 4, 5, 6, 7, 8, 9, 10]:   
     
        ppi_angle  = str(round(angle_list[isweep], 2)) 
        
        display.plot_ppi_map('reflectivity', sweep=isweep, vmin=-5, vmax=65,  cmap=cmaps.radar, 
                                title= "REF " + ppi_angle + '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  projection=projection, fig=fig, colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue"  )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue"  )

        display.plot_colorbar(shrink=0.7, ticks=np.arange(start=-5, stop=65, step=5))


        # display.set_limits(xlim, ylim)

        print( [ "Fig:  ", figprefix + "_ref_swp" + str(isweep).zfill(2) +  ".png"] )
        
        plt.savefig(figprefix + "_ref_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        # plt.show()
        plt.clf()   # clear all axis on this fig       
        
        #!  ======================================================== Plot ZDR 

        zdrnorm = BoundaryNorm(zdrlevels, ncolors=len(zdrcolors), clip=True)

        display.plot_ppi_map('differential_reflectivity', sweep=isweep, norm=zdrnorm,  cmap=cmapzdr, 
                                title= "ZDR " + ppi_angle + '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  projection=projection, fig=fig, colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_colorbar(shrink=0.7, ticks=zdrlevels)

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )

        print( [ "Fig:  ", figprefix + "_zdr_swp" + str(isweep).zfill(2) +  ".png"] )

        plt.savefig(figprefix + "_zdr_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        plt.clf()

        #!  ======================================================== Plot QPE

        display.plot_ppi_map('qpe', sweep=isweep, vmin=0.5, vmax=50,  cmap=cmaps.precip2_17lev, 
                                title= "QPE " + ppi_angle + '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  projection=projection, fig=fig, colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)
        
        display.plot_colorbar(shrink=0.7)
        
        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )
        

        print( [ "Fig:  ", figprefix + "_qpe_swp" + str(isweep).zfill(2) +  ".png"] )
        
        plt.savefig(figprefix + "_qpe_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        plt.clf()


        #!  ======================================================== Plot KDP

        kdpnorm = BoundaryNorm(kdplevels, ncolors=len(kdpcolors), clip=True)

        display.plot_ppi_map('specific_differential_phase', sweep=isweep, norm=kdpnorm,  cmap=cmapkdp, 
                                title= "KDP " + ppi_angle + '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  projection=projection, fig=fig,  colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_colorbar(shrink=0.7, ticks=kdplevels)
        
        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )
        

        print( [ "Fig:  ", figprefix + "_kdp_swp" + str(isweep).zfill(2) +  ".png"] )
        
        plt.savefig(figprefix + "_kdp_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        plt.clf()


        # # ========================================================= Save figs
        # plt.tight_layout()
      

    return


# ===================================================================
# ===================================================================


def adjust_fhc_colorbar_for_pyart(cb, type='summer'):

    if type == 'summer':
        cb.set_ticks(np.arange(1.4, 10, 0.9))
        cb.ax.set_yticklabels(['Drizzle', 'Rain', 'Ice Crystals', 'Aggregates',
                               'Wet Snow', 'Vertical Ice', 'LD Graupel',
                               'HD Graupel', 'Hail', 'Big Drops'])
    if type == 'winter':
        cb.set_ticks(np.arange(1.4, 7, 0.9))
        cb.ax.set_yticklabels(['IC','Plates', 'Dendrites', 'Aggregates', 'Wet Snow',
                           'LightRain', 'Rain'])

    if type == 'art_ss':   # 9 classes
        cb.set_ticks(np.arange(1.4, 10, 0.9)) 
        cb.ax.set_yticklabels([ 'DrySnow','CRystals', 'LightRain', 'LDGraupel', 'Rain',
                                'VerticalIce', 'WetSnow', "MeltingHail", "HDgraupel", "Heavy Rain"])

    cb.ax.set_ylabel('HCA')
    cb.ax.tick_params(length=0)

    return cb

# =================================================================================

# def three_panel_plot(radar, sweep=0, var1='reflectivity', vmin1=0, vmax1=65,
#                    cmap1='RdYlBu_r', units1='dBZ', var2='VEL',
#                    vmin2=-20, vmax2=20, cmap2='RdYlBu_r', units2='m/s',
#                    var3 = 'HCA',vmin3=0,vmax3=8,cmap3='Reds',units3='',return_flag=False,norm3=None,
#                    xlim=[-150,150], ylim=[-150,150]):
#     display = pyart.graph.RadarDisplay(radar)
#     fig = plt.figure(figsize=(13,8))
#     ax1 = fig.add_subplot(311)
#     if radar.scan_type == 'rhi':
#         display.plot_rhi(var1, sweep=sweep, vmin=vmin1, vmax=vmax1, cmap=cmap1,
#                          colorbar_label=units1, mask_outside=True)
#         display.set_limits(xlim=xlim, ylim=ylim)
#         ax2 = fig.add_subplot(312)
#         display.plot_rhi(var2, sweep=sweep, vmin=vmin2, vmax=vmax2, cmap=cmap2,
#                          colorbar_label=units2, mask_outside=True)
#         display.set_limits(xlim=xlim, ylim=ylim)

#         ax3 = fig.add_subplot(313)
#         if norm3 is not None:
#             display.plot_rhi(var3, sweep=sweep, vmin=vmin3, vmax=vmax3, cmap=cmap3,
#                              colorbar_label=units3, mask_outside=True,norm=norm3)
#         else:
#             display.plot_rhi(var3, sweep=sweep, vmin=vmin3, vmax=vmax3, cmap=cmap3,
#                              colorbar_label=units3, mask_outside=True)

            
#         display.set_limits(xlim=xlim, ylim=ylim)       
    
#     else:
#         display.plot_ppi(var1, sweep=sweep, vmin=vmin1, vmax=vmax1, cmap=cmap1,
#                          colorbar_label=units1, mask_outside=True)
#         display.set_limits(xlim=xlim, ylim=ylim)
#         ax2 = fig.add_subplot(212)
#         display.plot_ppi(var2, sweep=sweep, vmin=vmin2, vmax=vmax2, cmap=cmap2,
#                          colorbar_label=units2, mask_outside=True)
#         display.set_limits(xlim=xlim, ylim=ylim)
#     if return_flag:
#         return fig, ax1, ax2, ax3,display

# =========================================================================
# =========================================================================
# =========================================================================


# =================================================================================================================

def plot_hca_map(rdabufr, figprefix, plot_radius):

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    angle_list = rdabufr.fixed_angle['data']

    rng_list = np.arange(50, plot_radius, step=50 )

    display = pyart.graph.RadarMapDisplay(rdabufr)

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]
    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)       

    # pdb.set_trace()
    dd = (180/np.pi) * plot_radius*1e3/6371e3

    xlim = (rda_lon-dd, rda_lon+dd)
    ylim = (rda_lat-dd, rda_lat+dd)

    fig = plt.figure(figsize=(16, 16))

    for isweep in [0, 2, 4, 5, 6, 7, 8, 9, 10]:
        
        ppi_angle  = str(round(angle_list[isweep], 2)) 

        # ============================================================================= plot csu_hidt

        px = display.plot_ppi_map('csu_hidt', sweep=isweep, vmin=0, vmax=10, cmap=cmaphid, 
                                   title='CSU_HIDT '+ppi_angle + '\n ' + scantime,
                                   min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                   resolution='10m',  projection=projection, colorbar_flag=False )

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )

        # pdb.set_trace()

        display.plot_colorbar(shrink=0.7)
        num_cbs = len(display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[num_cbs-1], 'summer')
        
        # display.set_limits((-200, 200), (-200, 200), ax=ax)

        plt.savefig(figprefix + "_swp" + str(isweep).zfill(2) +  "_CSU_HIDT.png", 
                    dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        print( [ "Fig:  ", figprefix + "_swp" + str(isweep).zfill(2) +  "_CSU_HIDT.png"] )
        plt.clf()

        # ============================================================================= plot art_hca
        
        display.plot_ppi_map('art_hca', sweep=isweep, vmin=0, vmax=10, cmap=cmap_sshca, 
                               title='ART_SS '+ppi_angle + '\n ' + scantime,
                               min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                               resolution='10m',  projection=projection, colorbar_flag=False )        

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )

        display.plot_colorbar(shrink=0.7)
        num_cbs = len(display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[num_cbs-1], 'art_ss')
        # display.set_limits((-200, 200), (-200, 200), ax=ax)
        
        plt.savefig(figprefix + "_swp" + str(isweep).zfill(2) +  "_ART_SS.png", 
                     dpi=fig_dpi, pad_inches=0, bbox_inches="tight")

        print( [ "Fig:  ", figprefix + "_swp" + str(isweep).zfill(2) +  "_ART_SS.png"] )
        plt.clf()
        
        # plt.show()

    return

# ============================================================================================== 
# ============================================================================================== 
# ============================================================================================== 


def plot_rda_sect(rdabufr, figprefix, azimuth_list, plot_radius):

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    xsect = pyart.util.cross_section_ppi(rdabufr, azimuth_list, az_tol=2 )
 
    display = pyart.graph.RadarDisplay(xsect)

    fig = plt.figure(figsize=(16, 20))    

    for ith in np.arange(len(azimuth_list)): 

        # ==============================================================================REF
        
        ax = fig.add_subplot(321)
        display.plot('reflectivity', ith, ax=ax, 
                        title='REF ' + str(azimuth_list[ith]) + '\n ' + scantime,  
                        colorbar_label='REF',  
                        axislabels=('', 'Height /km'), vmin=-5, vmax=70, cmap=cmaps.radar)
        display.set_limits((0, plot_radius), (0, 20), ax=ax)

        # ============================================================================== VEL
        ax = fig.add_subplot(322)
        display.plot('velocity', ith, ax=ax, title='VEL ' + str(azimuth_list[ith]) + '\n ' + scantime, 
                        colorbar_label='VEL',  
                        axislabels=('', 'Height /km'), vmin=-30, vmax=30, cmap=cmaps.testcmap)
        display.set_limits((0,plot_radius), (0, 20), ax=ax)

        # ============================================================================== ZDR
        ax = fig.add_subplot(323)
        zdrnorm = BoundaryNorm(zdrlevels, ncolors=len(zdrcolors), clip=True)
        display.plot('differential_reflectivity', ith, ax=ax, title='ZDR '+ str(azimuth_list[ith]) + '\n ' + scantime, 
                      norm=zdrnorm, cmap=cmapzdr,
                      colorbar_label='ZDR',  axislabels=('', 'Height /km'))
        display.set_limits((0, plot_radius), (0, 20), ax=ax)

        # ============================================================================== KDP
        ax = fig.add_subplot(324)
        # display.plot('specific_differential_phase', ith, ax=ax, title='KDP '+ str(azimuth_list[ith]) , vmin=-5, vmax=20, cmap=cmapkdp,
        #             colorbar_label='KDP',  axislabels=('', 'Height \km'))
        # pdb.set_trace()
        kdpnorm = BoundaryNorm(kdplevels, ncolors=len(kdpcolors), clip=True)
        display.plot('specific_differential_phase', ith, ax=ax, 
                        title='KDP '+ str(azimuth_list[ith]) + '\n ' + scantime, 
                        cmap=cmapkdp, norm=kdpnorm,
                        colorbar_label='KDP',  axislabels=('', 'Height /km'))
        display.set_limits((0, plot_radius), (0, 20), ax=ax)
        
        # =============================================================================== CSU_HIDT
        ax = fig.add_subplot(325)
        display.plot('csu_hidt' , ith, vmin=0, vmax=10, cmap=cmaphid)
        numcbs = len( display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[numcbs-1], 'summer')
        display.set_limits((0, plot_radius), (0, 20), ax=ax)

        # ================================================================================ ART_HCA
        ax = fig.add_subplot(326)
        display.plot('art_hca', ith, vmin=0, vmax=10, cmap=cmap_sshca)
        numcbs = len( display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[numcbs-1], 'art_ss')
        display.set_limits((0, plot_radius), (0, 20), ax=ax)

        # ================================================================================ Save figs
        
        plt.tight_layout()

        print( [ "Fig:  ", figprefix + "_sec_" + str(azimuth_list[ith]).zfill(3) +  ".png"] )
        plt.savefig(figprefix + "_sec_" + str(azimuth_list[ith]).zfill(3) +  ".png", 
                    dpi=fig_dpi, pad_inches=0, bbox_inches="tight")
        

        plt.clf()

    return

