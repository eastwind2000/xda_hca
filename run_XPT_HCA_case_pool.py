
import os

import multiprocessing 

import pdb

# =================================================================================
# adding commenets and git version control
# =================================================================================

# ============================================== 2024-01-25 16:23:00 by chentao@cma.gov.cn

# rdafname = "//home/lse/haikui_tc2311/data/xrda/ZG003/Z_RADR_I_ZG003_20230907060004_O_DOR_AXPT0364_CRA_FMT.bin"
# rdafname = "/home/lse/haikui_tc2311/data/xrda/ZG000/Z_RADR_I_ZG000_20230907060007_O_DOR_AXPT0364_CRA_FMT.bin"

# ===============================================================================================

# xrda_datadir  =  "/home/lse/haikui_tc2311/data/xrda/  #! TC2311 Haikui, X-band phased-array radar

# /home/lse/237rainstorm/data/xrdas/ZA001/Z_RADR_I_ZA001_20230730000000_O_DOR_YLD2-D_CAP_FMT.bin.bz2  #! Beijing 23.7 rainstorm, VP23D, X-band

XRDA_DATADIR  =  "/home/lse/237rainstorm/data/xrdas/" 

RDASITE_LIST = [ "ZA001"  ]

NPROC = 3

WORK_DIR = "/mnt/f/fcst2022/scsmex2022/radarpro/XPT_HCA/"

FIG_ROOT_DIR = WORK_DIR + "../figs/"


# ===============================================================================================

def met_cmd_exec( cmd  ):
    print(cmd)
    os.system(cmd)
    return

if __name__ == "__main__":

    pool = multiprocessing.Pool(processes = NPROC)    # setting up  parallel processes using pool  

    for rdasite_id in RDASITE_LIST:

        flist = os.popen("ls -1 "+ XRDA_DATADIR  + rdasite_id + "/*" ).readlines()

        for txtbufr in flist[0:1]:

            rdafname = txtbufr.strip("\n")

            cmdbufr = "python " + WORK_DIR + "./XPT_HCA_dev.py --rdafname=" + rdafname + " --figroot=" + FIG_ROOT_DIR   #! pay attention to last splash "/"

            print(cmdbufr)

            os.system(cmdbufr)            
            # breakpoint()

    #         pool.apply_async( met_cmd_exec, (cmdbufr, ) )  # Adding parallel jobs into pool            #! pay attention (cmdbufr, ) "," can not be ignored
    # pool.close()    
    # pool.join()


