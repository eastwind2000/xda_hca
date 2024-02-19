import os
from sys import breakpointhook
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm

import pyart

import cartopy.crs as ccrs

import cmaps

import warnings 

import bz2

from module_io_defineXPT import FillValue

"""

cmap="pyart_NWSVel"

https://arm-doe.github.io/pyart/examples/plotting/plot_choose_a_colormap.html

"""


# ========================================================================================

warnings.filterwarnings("ignore")

## setup the summre colors.
hid_colors =  [ 'White',      
                'LightBlue',   # light rain
                'MediumBlue',  # rain
                'Darkorange',  # ice crystals
                'LightPink',   # aggregates
                'Cyan',        # wetsnow
                'DarkGray',    # verticle ice
                'Lime',        # LD graupel
                'green',      # HD graupel 
                'Yellow',         # hail  # red
                'Fuchsia'      # BD 
                ]
cmaphid = colors.ListedColormap(hid_colors)

sshca_colors = [ 'white', 
                 'LightPink',    # dyrsnow
                 'Darkorange',   # CR
                 'LightBlue',    # lightrain
                 'Lime',         # LD graupel
                 'MediumBlue',   # Rain
                 'DarkGray',     # Verticle Ice
                 'Cyan',         # wet snow
                 'Yellow',        # Melting Hail  # red
                 'green',       # HD Graupel 
                 'Fuchsia'       # Heavy Rain
                ]

cmap_sshca = colors.ListedColormap(sshca_colors)

##Setup the winter colors.
hid_colors_winter = ['White','Orange', 'Purple', 'Fuchsia', 'Pink', 'Cyan', 'LightBlue', 'Blue'] #, 'Lime', 'Yellow', 'Red', 'Fuchsia']
cmaphidwinter = colors.ListedColormap(hid_colors_winter)

# =========================================================================================================

reflevels = np.arange(-5.0, 70.0, step=5.0)
cmapref   = cmaps.radar
refnorm = BoundaryNorm(reflevels, ncolors=cmapref.N, clip=True)

# =========================================================================================================
vellevels = np.array([-24, -20, -16, -12, -8, -4, -1, 1, 4, 8, 12, 16, 20, 24])
velcolors = [ "#7be3ff", "#00e3ff", "#00b2b5", "#00ff00", "#00c700", "#008200", "#ffffff", 
              "#ff0000", "#ff5a5a", "#ffb5b5", "#ff7d00", "#ffd300", "#ffff00", "#7c007c"   ]
cmapvel = colors.ListedColormap(velcolors)
velnorm = BoundaryNorm(vellevels, ncolors=len(velcolors), clip=True)

# =========================================================================================================
kdplevels = np.array([-0.4,  -0.2,  -0.1,  0.1,  0.15,  0.22,  0.33,  0.5, 0.75,  1.1,  1.7,  2.4,  3.1,  7,  20])
kdpcolors = [ "#00efef", "#00979a", "#b4b4b4", "#b4b4b4", "#00c027", "#00e80a", "#24ff24", "#ffff1e", 
              "#ffff1e", "#ffe600", "#ffbc00", "#ff9800", "#ff5e00", "#f20f00", "#bb003a", "#ff00ff" ]
cmapkdp =   colors.ListedColormap(kdpcolors)
kdpnorm = BoundaryNorm(kdplevels, ncolors=len(kdpcolors), clip=True)

# =========================================================================================================

zdrlevels = np.array([ -3, -2, -1, 0, 0.2, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5  ])
zdrcolors = [ "#6b6d6a", "#949693", "#cbcbd0", "#dcf3de", "#01c21f", "#00eb0b", "#22fd22", 
              "#fffe19", "#ffe601", "#ffbc02", "#ff9a00", "#fe5d00", "#f50e00", "#bd003a", "#ff00fd" ]
cmapzdr =  colors.ListedColormap(zdrcolors)
zdrnorm = BoundaryNorm(zdrlevels, ncolors=len(zdrcolors), clip=True)

# =========================================================================================================

qpelevels = np.arange(5, 80, 5)
cmapqpe   = cmaps.precip2_17lev
qpenorm = BoundaryNorm(qpelevels, ncolors=cmapqpe.N, clip=True)

# =========================================================================================================

norm_dic  = { 
              "reflectivity":refnorm, 
              "velocity":velnorm, 
              "differential_reflectivity":zdrnorm, 
              "specific_differential_phase":kdpnorm, 
              "qpe":qpenorm
            }

cmap_dic  = { 
              "reflectivity":cmapref, 
              "velocity":cmapvel, 
              "differential_reflectivity":cmapzdr, 
              "specific_differential_phase":cmapkdp, 
              "qpe":cmapqpe
            }

varlevels_dic = {      "reflectivity":reflevels, 
                        "velocity":vellevels, 
                        "differential_reflectivity":zdrlevels, 
                        "specific_differential_phase":kdplevels, 
                        "qpe":qpelevels
                    }

var_shortname = {  "reflectivity":"ref", 
                   "velocity":"vel", 
                   "differential_reflectivity":"zdr", 
                   "specific_differential_phase":"kdp", 
                   "qpe":"qpe"                    
                }

# =====================================================================

shp_path = "/home/lse/res/ncarg/shp/"

shp_list = [
            "/home/lse/res/ncarg/shp/bou1_4l.shp", 
            "/home/lse/res/ncarg/shp/continents_lines.shp", 
            "/home/lse/res/ncarg/shp/hyd1_4l.shp", 
            # "/home/lse/res/ncarg/shp/province_l.shp", 
            # "/home/lse/res/ncarg/shp/City_l.shp"            
            ]

plot_dpi = 300

# ======================================================================================

# =========================================================================
#  plotting XPAR parameters [REF,VEL,ZDR,KDP,QPE] on PPI
# =========================================================================

def plot_rda_parameter(rdabufr, varlist,  figprefix, plot_radius, rng_list):    

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]

    angle_list = rdabufr.fixed_angle['data']

    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)   #! key function   

    dd = (180/np.pi)*plot_radius*1e3/6371e3  #!unit in degree

    xscale = 1.2

    xlim = (rda_lon-dd*xscale, rda_lon+dd*xscale)
    ylim = (rda_lat-dd, rda_lat+dd)

    nsweeps = rdabufr.nsweeps

    # var_list = ["velocity", "reflectivity",  "differential_reflectivity", "specific_differential_phase", "qpe"   ]

    # ===========================================================
    
    display = pyart.graph.RadarMapDisplay(rdabufr)        #! key object

    fig = plt.figure(figsize=(16, 16))

    ax = plt.subplot( projection=projection )

    for varname in varlist:
    
        swp_list = np.arange(nsweeps)       
        
        for isweep in swp_list:
        # ======================================================= vel
        
            ppi_angle  = str(round(angle_list[isweep], 2)) 

            var_mask = (varname, FillValue+0.01)

            vsn = var_shortname[varname]  #var shortname

            display.plot_ppi_map(   varname, sweep=isweep, mask_tuple=var_mask, 
                                    cmap=cmap_dic[varname], norm=norm_dic[varname],
                                    title= varname + " " + ppi_angle + '\n ' + scantime,
                                    min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                    resolution='10m',  
                                    colorbar_flag=False,
                                    lat_0=rda_lat, lon_0=rda_lon )

            display.plot_range_rings(rng_list, col="darkblue"  )
            display.plot_cross_hair(rng_list[-1])
            display.plot_line_geo( line_lons=[ xlim[0], xlim[1]  ], line_lats=[rda_lat, rda_lat ])
            display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ] )
            display.plot_colorbar(shrink=0.8, ticks=varlevels_dic[varname]) 

            print( [ "Fig:  ", figprefix + "_" + vsn + "_swp" + str(isweep).zfill(2) +  ".png"] )

            plt.savefig(figprefix + "_" + vsn + "_swp" + str(isweep).zfill(2) + ".png", 
                        dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
                       
            plt.clf()   # clear all axis on this fig
        
            # breakpoint()
    
    plt.close()

    return

# =================================================================================================================
# =================================================================================================================
# =================================================================================================================

def plot_hca_map(rdabufr, figprefix, rda_radius, rng_list):

    # pdb.set_trace()

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    angle_list = rdabufr.fixed_angle['data']

    # =================================================================================

    display = pyart.graph.RadarMapDisplay(rdabufr)

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]
    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)       

    # pdb.set_trace()===============================================

    dd = (180/np.pi) * rda_radius*1e3/6371e3  #! important unit trans
    
    xscale = 1.2
    xlim = (rda_lon-dd*xscale, rda_lon+dd*xscale)
    ylim = (rda_lat-dd, rda_lat+dd)

    fig = plt.figure(figsize=(16, 16))

    # if vcp_pattern in ["VCP21", "VCP23D"]:
    #     swp_list = [ 0, 2, 4, 5, 6, 7, 8, 9, 10 ]
    # else:
    #     swp_list = np.arange(rdabufr.nsweeps)
    
    swp_list = np.arange(rdabufr.nsweeps)

    for isweep in swp_list:
        
        ppi_angle  = str(round(angle_list[isweep], 2)) 

        # ============================================================================= plot csu_hidt

        px=display.plot_ppi_map('csu_hidt', sweep=isweep, vmin=0, vmax=10, cmap=cmaphid, 
                               title='CSU_HIDT ' + ppi_angle + '\n ' + scantime,
                               min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                               resolution='10m',  
                               colorbar_flag=False )

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
                    bbox_inches= "tight", dpi=400, pad_inches=0)
        print( [ "Fig:  ", figprefix + "_swp" + str(isweep).zfill(2) +  "_CSU_HIDT.png"] )
        plt.clf()

        # ============================================================================= plot art_hca
        
        display.plot_ppi_map('art_hca', sweep=isweep, vmin=0, vmax=10, cmap=cmap_sshca, 
                               title='ART_SS ' + ppi_angle  + '\n ' + scantime, 
                               min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                               resolution='10m',  
                               colorbar_flag=False )        

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue"  )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue"  )

        display.plot_colorbar(shrink=0.7)
        num_cbs = len(display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[num_cbs-1], 'art_ss')
        # display.set_limits((-200, 200), (-200, 200), ax=ax)
        
        plt.savefig(figprefix + "_swp" + str(isweep).zfill(2) +  "_ART_SS.png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
        print( [ "Fig:  ", figprefix + "_swp" + str(isweep).zfill(2) +  "_ART_SS.png"] )

        plt.clf()
        
        # plt.show()

    plt.close()

    return


# ========================================================================= 


def plot_rda_sect(rdabufr, rdabufr_dop, figprefix, azimuth_list, plot_radius):

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    xsect_log = pyart.util.cross_section_ppi(rdabufr, azimuth_list, az_tol=0.25 )

    xsect_dop = pyart.util.cross_section_ppi(rdabufr_dop, azimuth_list, az_tol=0.25 )
    
    display = pyart.graph.RadarDisplay(xsect_log)

    # display_dop = pyart.graph.RadarDisplay(xsect_dop)

    for ith in np.arange(len(azimuth_list)): 
        
        fig, axs = plt.subplots(3, 2, figsize=(16, 20))

        ax = axs.flatten()

        # ===================================================REF
        
        display.plot('reflectivity', ith,  ax=ax[0], edges=False,
                      title='REF ' + str(azimuth_list[ith])  + '\n ' + scantime, 
                      colorbar_label='REF',  
                      axislabels=('', 'Height /km'), vmin=-5, vmax=70, cmap=cmaps.radar)
        display.set_limits((0, plot_radius), (0, 20), ax=ax[0])

        # ===================================================VEL

        # display_dop.plot('velocity', ith, ax=ax, title='VEL ' + str(azimuth_list[ith]) + '\n ' + scantime, 
        #               colorbar_label='VEL',  
        #               axislabels=('', 'Height /km'), vmin=-30, vmax=30, cmap=cmaps.testcmap)
        # display_dop.set_limits((0,  plot_radius), (0, 20), ax=ax)

        # ===================================================ZDR

        zdrnorm = BoundaryNorm(zdrlevels, ncolors=len(zdrcolors), clip=True)
        display.plot('differential_reflectivity', ith, ax=ax[2], 
                     title='ZDR '+ str(azimuth_list[ith]) + '\n ' + scantime,  
                     norm=zdrnorm, cmap=cmapzdr,
                     colorbar_label='ZDR',  axislabels=('', 'Height /km'))
        display.set_limits((0,  plot_radius), (0, 20), ax=ax[2])

        # ===================================================KDP

        kdpnorm = BoundaryNorm(kdplevels, ncolors=len(kdpcolors), clip=True)
        display.plot('specific_differential_phase', ith, ax=ax[3], 
                    title='KDP '+ str(azimuth_list[ith]) + '\n ' + scantime, 
                    cmap=cmapkdp, norm=kdpnorm,
                    colorbar_label='KDP',  axislabels=('', 'Height /km'))
        display.set_limits((0,  plot_radius), (0, 20), ax=ax[3])
        
        # ===================================================CSU_HIDT

        display.plot('csu_hidt' , ith, vmin=0, vmax=10, cmap=cmaphid, ax=ax[4],
                      title='CSU HDIT '+ str(azimuth_list[ith]) + '\n ' + scantime   )
        numcbs = len( display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[numcbs-1], 'summer')
        display.set_limits((0,  plot_radius), (0, 20), ax=ax[4])

        # ==================================================ART_HCA

        display.plot('art_hca', ith, vmin=0, vmax=10, cmap=cmap_sshca, ax=ax[5],
                      title='ART SS '+ str(azimuth_list[ith]) + '\n ' + scantime   )
        numcbs = len( display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[numcbs-1], 'art_ss')
        display.set_limits((0,  plot_radius), (0, 20), ax=ax[5])

        # ================================================== save figs
        
        plt.tight_layout()

        # plt.show()

        fig_name = figprefix + "_sec_" + str(azimuth_list[ith]).zfill(3) +  ".png"

        print( [ "Fig:  ", fig_name ] )
        plt.savefig(fig_name, dpi=plot_dpi, pad_inches=0)
        
        breakpoint()

        # plt.close(fig)

    return

# ================================================================================================ auxiliary subs 

def add_field_to_radar_object(field, radar, field_name='FH', units='unitless', 
                              long_name='Hydrometeor ID', standard_name='Hydrometeor ID',
                              dz_field='CZ'):
    """
    Adds a newly created field to the Py-ART radar object. If reflectivity is a masked array,
    make the new field masked the same as reflectivity.
    """
    fill_value = -32768
    masked_field = np.ma.asanyarray(field)
    masked_field.mask = masked_field == fill_value
    if hasattr(radar.fields[dz_field]['data'], 'mask'):
        setattr(masked_field, 'mask', 
                np.logical_or(masked_field.mask, radar.fields[dz_field]['data'].mask))
        fill_value = radar.fields[dz_field]['_FillValue']
    field_dict = {'data': masked_field,
                  'units': units,
                  'long_name': long_name,
                  'standard_name': standard_name,
                  '_FillValue': fill_value}
    radar.add_field(field_name, field_dict, replace_existing=True)
    return radar

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

    if type == 'art_ss':                            #! == HCA in 9 classes, revised by chentao@cma.gov.cn 2022.06
        cb.set_ticks(np.arange(1.4, 10, 0.9)) 
        cb.ax.set_yticklabels([ 'DrySnow','CRystals', 'LightRain', 'LDGraupel', 'Rain',
                                'VerticalIce', 'WetSnow', "MeltingHail", "HDgraupel", "Heavy Rain"])

    cb.ax.set_ylabel('HCA')
    cb.ax.tick_params(length=0)

    return cb

# =============================================================================================================



def bak_plot_rda_parameter(rdabufr, figprefix, plot_radius, rng_list):    

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    angle_list = rdabufr.fixed_angle['data']

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]

    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)   #! key function   

    dd = (180/np.pi)*plot_radius*1e3/6371e3  #!unit in degree

    xscale = 1.1

    xlim = (rda_lon-dd*xscale, rda_lon+dd*xscale)
    
    ylim = (rda_lat-dd, rda_lat+dd)

    nsweeps = rdabufr.nsweeps

    # ===========================================================
    
    display = pyart.graph.RadarMapDisplay(rdabufr)        #! key object

    fig = plt.figure(figsize=(16, 16))

    ax = plt.subplot( projection=projection )
    
    # for isweep in np.arange(rdabufr.nsweeps):   
    for isweep in  np.arange(nsweeps):
        # ======================================================= vel
        
        ppi_angle  = str(round(angle_list[isweep], 2)) 

        var_mask = ("velocity", FillValue+0.01)

        display.plot_ppi_map('velocity', sweep=isweep, mask_tuple=var_mask, 
                              cmap=cmapvel, norm=velnorm, 
                              title= "VEL " + ppi_angle + '\n ' + scantime,
                              min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                              resolution='10m',  
                                # shapefile=[],
                              colorbar_flag=False,
                              lat_0=rda_lat, lon_0=rda_lon )
       # cartopy.crs.PlateCarree()

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ xlim[0], xlim[1]  ], line_lats=[rda_lat, rda_lat ] )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ]  )
        display.plot_colorbar(shrink=0.8, ticks=vellevels) 

        # display.plot_colorbar(shrink=0.8, ticks=np.arange(start=-30, stop=30, step=2) )
        # display.set_limits(xlim, ylim)

        print( [ "Fig:  ", figprefix + "_vel_swp" + str(isweep).zfill(2) +  ".png"] )

        plt.savefig(figprefix + "_vel_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
        
        plt.clf()   # clear all axis on this fig
    
    # breakpoint()
    # ============================================================

    # for isweep in np.arange(rdabufr.nsweeps):   

    for isweep in np.arange(nsweeps):
     
        ppi_angle  = str(round(angle_list[isweep], 2))        

        # ========================================================REF
        
        var_mask =  ("reflectivity", FillValue+0.01)
        display.plot_ppi_map('reflectivity', sweep=isweep, mask_tuple=var_mask, 
                                cmap=cmapref, norm=refnorm, 
                                title= "REF " + ppi_angle + '\n ' + scantime, 
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m', 
                                # shapefile=[],                                 
                                colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[  xlim[0], xlim[1] ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )
        display.plot_colorbar(shrink=0.7, ticks=reflevels)

        print( [ "Fig:  ", figprefix + "_ref_swp" + str(isweep).zfill(2) +  ".png"] )
        plt.savefig(figprefix + "_ref_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
    
        plt.clf()   # clear all axis on this fig       
        
        # ========================================================ZDR

        var_mask =  ("differential_reflectivity", FillValue+0.01)
        display.plot_ppi_map('differential_reflectivity', sweep=isweep, mask_tuple=var_mask,
                                norm=zdrnorm,  cmap=cmapzdr, 
                                title= "ZDR " + ppi_angle + '\n ' + scantime, 
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  
                                colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[  xlim[0], xlim[1]  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue")
        display.plot_colorbar(shrink=0.7, ticks=zdrlevels)

        print( [ "Fig:  ", figprefix + "_zdr_swp" + str(isweep).zfill(2) +  ".png"] )
        plt.savefig(figprefix + "_zdr_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")

        plt.clf()

        # # ========================================================QPE

        display.plot_ppi_map('qpe', sweep=isweep, vmin=0.5, vmax=50,  cmap=cmaps.precip2_17lev, 
                                title= "QPE " + ppi_angle+ '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  
                                colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)
        

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[  xlim[0], xlim[1]  ], line_lats=[rda_lat, rda_lat ], color="darkblue")
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )
        display.plot_colorbar(shrink=0.7)

        print( [ "Fig:  ", figprefix + "_qpe_swp" + str(isweep).zfill(2) +  ".png"] )
        plt.savefig(figprefix + "_qpe_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
        plt.clf()


        # # ========================================================KDP

        var_mask =  ("specific_differential_phase", FillValue+0.01)

        display.plot_ppi_map('specific_differential_phase', mask_tuple=var_mask,
                                sweep=isweep, norm=kdpnorm,  cmap=cmapkdp, 
                                title= "KDP " + ppi_angle + '\n ' + scantime, 
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  
                                colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)
        
        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ xlim[0], xlim[1] ], line_lats=[rda_lat, rda_lat ], color="darkblue", transform=projection )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue", transform=projection )
        display.plot_colorbar(shrink=0.7, ticks=kdplevels)

        print( [ "Fig:  ", figprefix + "_kdp_swp" + str(isweep).zfill(2) +  ".png"] )
        plt.savefig(figprefix + "_kdp_swp" + str(isweep).zfill(2) + ".png", 
                     dpi=plot_dpi, 
                     pad_inches=0,  bbox_inches= "tight")
        plt.clf()

        # ======================================================================================================== Save figs
        # plt.tight_layout()
      
    plt.close()

    return
