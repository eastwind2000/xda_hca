import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
from matplotlib.colors import BoundaryNorm

import pyart

import cartopy.crs as ccrs

import cmaps
import warnings 
import pdb
import gc

# from skewt import SkewT 
# from cinrad.io.export import standard_data_to_pyart
# from export_xpar2pyart import  xpar_data_to_pyart
# # from CSU_RadarTools.csu_radartools import csu_fhc
# from csu_radartools import (csu_fhc, csu_liquid_ice_mass, csu_blended_rain,
#                             csu_dsd, csu_kdp, csu_misc, fundamentals)


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

kdplevels = np.array([-0.4,  -0.2,  -0.1,  0.1,  0.15,  0.22,  0.33,  0.5, 0.75,  1.1,  1.7,  2.4,  3.1,  7,  20])
kdpcolors = [ "#00efef", "#00979a", "#b4b4b4", "#b4b4b4", "#00c027", "#00e80a", "#24ff24", "#ffff1e", 
              "#ffff1e", "#ffe600", "#ffbc00", "#ff9800", "#ff5e00", "#f20f00", "#bb003a", "#ff00ff" ]
cmapkdp =   colors.ListedColormap(kdpcolors)

zdrlevels = np.array([ -3, -2, -1, 0, 0.2, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5  ])
zdrcolors = [ "#6b6d6a", "#949693", "#cbcbd0", "#dcf3de", "#01c21f", "#00eb0b", "#22fd22", 
              "#fffe19", "#ffe601", "#ffbc02", "#ff9a00", "#fe5d00", "#f50e00", "#bd003a", "#ff00fd" ]
cmapzdr =  colors.ListedColormap(zdrcolors)

shp_path = "/home/lse/res/ncarg/shp/"

shp_list = [
            "/home/lse/res/ncarg/shp/bou1_4l.shp", 
            "/home/lse/res/ncarg/shp/continents_lines.shp", 
            "/home/lse/res/ncarg/shp/hyd1_4l.shp", 
            # "/home/lse/res/ncarg/shp/province_l.shp", 
            # "/home/lse/res/ncarg/shp/City_l.shp"            
            ]

plot_dpi = 400

# ======================================================================================

def interpolate_sounding_to_radar(sounding, radar):
    """Takes sounding data and interpolates it to every radar gate."""
    radar_z = get_z_from_radar(radar)
    radar_T = None
    snd_T, snd_z = check_sounding_for_montonic(sounding)
    shape = np.shape(radar_z)
    rad_z1d = radar_z.ravel()
    rad_T1d = np.interp(rad_z1d, snd_z, snd_T)
    return np.reshape(rad_T1d, shape), radar_z

def check_sounding_for_montonic(sounding):
    """
    So the sounding interpolation doesn't fail, force the sounding to behave
    monotonically so that z always increases. This eliminates data from
    descending balloons.
    """
    snd_T = sounding['temp']  # In old SkewT, was sounding.data
    snd_z = sounding['hght']  # In old SkewT, was sounding.data
    dummy_z = []
    dummy_T = []
    #print(snd_T.mask)
    # if not snd_T.mask.all(): #May cause issue for specific soundings
    dummy_z.append(snd_z[0])
    dummy_T.append(snd_T[0])
    for i, height in enumerate(snd_z):
        if i > 0:
            if snd_z[i] > snd_z[i-1]:
                if np.isfinite(snd_z[i]) and np.isfinite(snd_T[i]):
                    #print(snd_T[i])
                    dummy_z.append(snd_z[i])
                    dummy_T.append(snd_T[i])
    snd_z = np.array(dummy_z)
    snd_T = np.array(dummy_T)
    # else:
    #     print('uh-oh. sounding problem')
    # pdb.set_trace()
    return snd_T, snd_z

# ====================================================================

def get_z_from_radar(radar):
    """Input radar object, return z from radar (km, 2D)"""
    azimuth_1D = radar.azimuth['data']
    elevation_1D = radar.elevation['data']
    srange_1D = radar.range['data']
    sr_2d, az_2d = np.meshgrid(srange_1D, azimuth_1D)
    el_2d = np.meshgrid(srange_1D, elevation_1D)[1]
    xx, yy, zz = radar_coords_to_cart(sr_2d/1000.0, az_2d, el_2d)
    #print(np.shape(zz), np.shape(radar.altitude['data']))
    try:
        return zz + radar.altitude['data'][0]
    except:
        return zz + radar.altitude['data']


def radar_coords_to_cart(rng, az, ele, debug=False):
    """
    TJL - taken from old Py-ART version
    Calculate Cartesian coordinate from radar coordinates
    Parameters
    ----------
    rng : array
        Distances to the center of the radar gates (bins) in kilometers.
    az : array
        Azimuth angle of the radar in degrees.
    ele : array
        Elevation angle of the radar in degrees.
    Returns
    -------
    x, y, z : array
        Cartesian coordinates in meters from the radar.
    Notes
    -----
    The calculation for Cartesian coordinate is adapted from equations
    2.28(b) and 2.28(c) of Doviak and Zrnic [1]_ assuming a
    standard atmosphere (4/3 Earth's radius model).
    .. math::
        z = \\sqrt{r^2+R^2+r*R*sin(\\theta_e)} - R
        s = R * arcsin(\\frac{r*cos(\\theta_e)}{R+z})
        x = s * sin(\\theta_a)
        y = s * cos(\\theta_a)
    Where r is the distance from the radar to the center of the gate,
    :math:\\theta_a is the azimuth angle, :math:\\theta_e is the
    elevation angle, s is the arc length, and R is the effective radius
    of the earth, taken to be 4/3 the mean radius of earth (6371 km).
    References
    ----------
    .. [1] Doviak and Zrnic, Doppler Radar and Weather Observations, Second
        Edition, 1993, p. 21.
    """
    theta_e = ele * np.pi / 180.0  # elevation angle in radians.
    theta_a = az * np.pi / 180.0  # azimuth angle in radians.
    R = 6371.0 * 1000.0 * 4.0 / 3.0  # effective radius of earth in meters.
    r = rng * 1000.0  # distances to gates in meters.

    z = (r ** 2 + R ** 2 + 2.0 * r * R * np.sin(theta_e)) ** 0.5 - R
    s = R * np.arcsin(r * np.cos(theta_e) / (R + z))  # arc length in m.
    x = s * np.sin(theta_a)
    y = s * np.cos(theta_a)
    return x, y, z


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

# =========================================================================
# ! plotting XPA-RDA parameters [REF,VEL,ZDR,KDP,QPE] on PPI
# =========================================================================

def plot_rda_parameter(rdabufr, figprefix, plot_radius):    

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    angle_list = rdabufr.fixed_angle['data']

    rng_list = [10.0, 20.0, 30.0, 40.0, 50.0]

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]
    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)       

    # rdabufr.nsweeps
    # pdb.set_trace()

    dd = (180/np.pi)*plot_radius*1e3/6371e3  #!unit in degree

    xlim = (rda_lon-dd, rda_lon+dd)
    
    ylim = (rda_lat-dd, rda_lat+dd)
    
    display = pyart.graph.RadarMapDisplay(rdabufr)      

    fig= plt.figure(figsize=(16, 16))
    
    # for isweep in np.arange(rdabufr.nsweeps):   
    for isweep in [0, 1]:
        # ======================================================= vel
        ppi_angle  = str(round(angle_list[isweep], 2)) 

        display.plot_ppi_map_shplist('velocity', sweep=isweep, vmin=-30, vmax=30,  cmap=cmaps.testcmap, 
                                title= "VEL " + ppi_angle + '\n ' + scantime,
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m',  
                                shapefile=[],
                                projection=projection, fig=fig, colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ] )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ] )

        display.plot_colorbar(shrink=0.7, ticks=np.arange(start=-30, stop=30, step=2) )
        # display.set_limits(xlim, ylim)

        print( [ "Fig:  ", figprefix + "_vel_swp" + str(isweep).zfill(2) +  ".png"] )

        plt.savefig(figprefix + "_vel_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
        
        # plt.show()

        plt.clf()   # clear all axis on this fig

    # ============================================================
    
    
    # for isweep in np.arange(rdabufr.nsweeps):   
    for isweep in [0, 1]:
     
        ppi_angle  = str(round(angle_list[isweep], 2))        

        # ========================================================REF
        
        display.plot_ppi_map_shplist('reflectivity', sweep=isweep, vmin=-5, vmax=65,  cmap=cmaps.radar, 
                                title= "REF " + ppi_angle + '\n ' + scantime, 
                                min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                                resolution='10m', shapefile=[],                                 
                                projection=projection, fig=fig, colorbar_flag=False,
                                lat_0=rda_lat, lon_0=rda_lon)

        display.plot_range_rings(rng_list, col="darkblue"  )
        display.plot_cross_hair(rng_list[-1])
        display.plot_line_geo( line_lons=[ rda_lon-dd, rda_lon+dd  ], line_lats=[rda_lat, rda_lat ], color="darkblue" )
        display.plot_line_geo( line_lons=[ rda_lon, rda_lon  ], line_lats=[rda_lat-dd, rda_lat+dd ], color="darkblue" )

        display.plot_colorbar(shrink=0.7, ticks=np.arange(start=-5, stop=65, step=5))

        # display.set_limits(xlim, ylim)

        print( [ "Fig:  ", figprefix + "_ref_swp" + str(isweep).zfill(2) +  ".png"] )
        plt.savefig(figprefix + "_ref_swp" + str(isweep).zfill(2) + ".png", 
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
        
        # plt.show()

        plt.clf()   # clear all axis on this fig       
        
        # ========================================================ZDR

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
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")

        plt.clf()

        # # ========================================================QPE

        display.plot_ppi_map('qpe', sweep=isweep, vmin=0.5, vmax=50,  cmap=cmaps.precip2_17lev, 
                                title= "QPE " + ppi_angle+ '\n ' + scantime,
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
                    dpi=plot_dpi, pad_inches=0,  bbox_inches= "tight")
        plt.clf()


        # # ========================================================KDP

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
                     dpi=plot_dpi, 
                     pad_inches=0,  bbox_inches= "tight")
        plt.clf()


        # # ========================================================= Save figs
        # plt.tight_layout()
      

    return

# =================================================================================================================
# =================================================================================================================
# =================================================================================================================

def plot_hca_map(rdabufr, figprefix, rda_radius):

    # pdb.set_trace()

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    angle_list = rdabufr.fixed_angle['data']

    rng_list = [10.0, 20.0, 30.0, 40.0, 50.0 ]

    display = pyart.graph.RadarMapDisplay(rdabufr)

    rda_lat = rdabufr.latitude['data'][0]
    rda_lon = rdabufr.longitude['data'][0]
    projection = ccrs.LambertConformal(central_longitude=rda_lon, central_latitude=rda_lat)       

    # pdb.set_trace()

    dd = (180/np.pi) * rda_radius*1e3/6371e3  #! important unit trans
    
    ylim = (rda_lat-dd, rda_lat+dd)
    xlim = (rda_lon-dd, rda_lon+dd)

    fig = plt.figure(figsize=(16, 16))

    for isweep in np.arange(rdabufr.nsweeps):
        
        ppi_angle  = str(round(angle_list[isweep], 2)) 

        # ============================================================================= plot csu_hidt

        px=display.plot_ppi_map('csu_hidt', sweep=isweep, vmin=0, vmax=10, cmap=cmaphid, 
                               title='CSU_HIDT ' + ppi_angle + '\n ' + scantime,
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
                    bbox_inches= "tight", dpi=400, pad_inches=0)
        print( [ "Fig:  ", figprefix + "_swp" + str(isweep).zfill(2) +  "_CSU_HIDT.png"] )
        plt.clf()

        # ============================================================================= plot art_hca
        
        display.plot_ppi_map('art_hca', sweep=isweep, vmin=0, vmax=10, cmap=cmap_sshca, 
                               title='ART_SS ' + ppi_angle  + '\n ' + scantime, 
                               min_lon=xlim[0], max_lon=xlim[1], min_lat=ylim[0], max_lat=ylim[1],
                               resolution='10m',  projection=projection, colorbar_flag=False )        

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

    return


# ========================================================================= 


def plot_rda_sect(rdabufr, figprefix, azimuth_list, plot_radius):

    scantime = os.popen("basename " + figprefix).readlines()[0].strip("\n")

    xsect = pyart.util.cross_section_ppi(rdabufr, azimuth_list, az_tol=0.5 )
 
    display = pyart.graph.RadarDisplay(xsect)

    fig = plt.figure(figsize=(16, 20))
    

    for ith in np.arange(len(azimuth_list)): 

        # ===================================================REF
        
        ax = fig.add_subplot(321)
        display.plot('reflectivity', ith, ax=ax, title='REF ' + str(azimuth_list[ith])  + '\n ' + scantime, 
                      colorbar_label='REF',  
                      axislabels=('', 'Height /km'), vmin=-5, vmax=70, cmap=cmaps.radar)
        display.set_limits((0, plot_radius), (0, 20), ax=ax)

        # ===================================================VEL
        ax = fig.add_subplot(322)
        display.plot('velocity', ith, ax=ax, title='VEL ' + str(azimuth_list[ith]) + '\n ' + scantime, 
                      colorbar_label='VEL',  
                      axislabels=('', 'Height /km'), vmin=-30, vmax=30, cmap=cmaps.testcmap)
        display.set_limits((0,  plot_radius), (0, 20), ax=ax)

        # ===================================================ZDR
        ax = fig.add_subplot(323)
        zdrnorm = BoundaryNorm(zdrlevels, ncolors=len(zdrcolors), clip=True)
        display.plot('differential_reflectivity', ith, ax=ax, 
                     title='ZDR '+ str(azimuth_list[ith]) + '\n ' + scantime,  
                     norm=zdrnorm, cmap=cmapzdr,
                     colorbar_label='ZDR',  axislabels=('', 'Height /km'))
        display.set_limits((0,  plot_radius), (0, 20), ax=ax)

        # ===================================================KDP
        ax = fig.add_subplot(324)
        # display.plot('specific_differential_phase', ith, ax=ax, title='KDP '+ str(azimuth_list[ith]) , vmin=-5, vmax=20, cmap=cmapkdp,
        #             colorbar_label='KDP',  axislabels=('', 'Height \km'))
        # pdb.set_trace()
        kdpnorm = BoundaryNorm(kdplevels, ncolors=len(kdpcolors), clip=True)
        display.plot('specific_differential_phase', ith, ax=ax, 
                    title='KDP '+ str(azimuth_list[ith]) + '\n ' + scantime, 
                    cmap=cmapkdp, norm=kdpnorm,
                    colorbar_label='KDP',  axislabels=('', 'Height /km'))
        display.set_limits((0,  plot_radius), (0, 20), ax=ax)
        
        # ===================================================CSU_HIDT
        ax = fig.add_subplot(325)
        display.plot('csu_hidt' , ith, vmin=0, vmax=10, cmap=cmaphid,
                      title='CSU HDIT '+ str(azimuth_list[ith]) + '\n ' + scantime   )
        numcbs = len( display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[numcbs-1], 'summer')
        display.set_limits((0,  plot_radius), (0, 20), ax=ax)

        # ==================================================ART_HCA
        ax = fig.add_subplot(326)
        display.plot('art_hca', ith, vmin=0, vmax=10, cmap=cmap_sshca, 
                      title='ART SS '+ str(azimuth_list[ith]) + '\n ' + scantime   )
        numcbs = len( display.cbs)
        adjust_fhc_colorbar_for_pyart(display.cbs[numcbs-1], 'art_ss')
        display.set_limits((0,  plot_radius), (0, 20), ax=ax)

        # ================================================== save figs
        
        plt.tight_layout()

        print( [ "Fig:  ", figprefix + "_sec_" + str(azimuth_list[ith]).zfill(3) +  ".png"] )
        plt.savefig(figprefix + "_sec_" + str(azimuth_list[ith]).zfill(3) +  ".png", 
                    dpi=plot_dpi, pad_inches=0)
        
        plt.clf()

    return

