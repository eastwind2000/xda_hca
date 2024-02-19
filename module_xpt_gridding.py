
import pyart 

import os

import numpy as np

import matplotlib.pyplot as plt

import datetime

import cmaps

def xpt_gridding(rdabufr, xpar_radius, varname_list):

    rda_site = rdabufr.metadata["instrument_name"]

    ftime_utc = datetime.datetime.utcfromtimestamp( rdabufr.time['data'][0]).strftime("%Y%m%d%H%M")

    xpargrid_dx  = 50.0
    xpargrid_dy  = 50.0
    xpargrid_dz  = 150

    xpargrid_nx = np.int32(xpar_radius*1000/xpargrid_dx) -1
    xpargrid_ny = np.int32(xpar_radius*1000/xpargrid_dy) -1
    xpargrid_nz = np.int32(6000/xpargrid_dz) -1

    grid_range = ((0,  12000), (-xpar_radius*1000, xpar_radius*1000), (-xpar_radius*1000, xpar_radius*1000))

    grid_dim   = (xpargrid_nz, xpargrid_ny, xpargrid_nx)


    rdagrid = pyart.map.grid_from_radars(rdabufr, grid_shape=grid_dim, grid_limits=grid_range,
                                         fields=varname_list)
    if "velocity" not in varname_list:

        rdagrid_fname = "../data/" + rda_site + "_rdagrid_"  + ftime_utc +  ".nc"
    else: 

        rdagrid_fname = "../data/" + rda_site + "_rdagrid_VEL_"  + ftime_utc +  ".nc"

    pyart.io.write_grid(rdagrid_fname, rdagrid, format="NETCDF4", write_point_lon_lat_alt=True, arm_alt_lat_lon_variables=False )

    print(["Write grid_rdafile: ", rdagrid_fname])

    # ============================================================================== debug code

    debug_fig = True

    if debug_fig:

        fig = plt.figure()
        ax = fig.add_subplot(121)
        ptx = ax.imshow(rdagrid.fields['reflectivity']['data'].data[0, :, :], origin='lower',
                        vmin=-5, vmax=70,  cmap=cmaps.radar)               
        plt.colorbar(ptx, ax=ax)

        # ax = fig.add_subplot(122)
        # ptx = ax.imshow(rdagrid.fields['']['data'].data[10, :, :], origin='lower',
        #                 vmin=-30, vmax=30,  cmap=plt.cm.bwr)               
        # plt.colorbar(ptx, ax=ax)
        
        plt.show()

        # breakpoint()
