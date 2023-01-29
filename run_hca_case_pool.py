
import os

import multiprocessing

import pdb

from define_hca_params import *

# adding commenets and git version control
# sudo mount.cifs -t cifs   //10.10.31.21/data/rda/figs  ./remote_figdir/  -o user=root,pass="nmc601)(*"

# =============================================== testing git hub  change 2023-01-25 22:43:47

def rda_cmd_exec(cmd):
    print(cmd)
    os.system(cmd)
    return

if __name__ == "__main__":

    pool = multiprocessing.Pool(processes=nproc)

    for rdasite_id in rdasite_list:

        flist = os.popen("ls -1 " + rda_dir + rdasite_id + "/*").readlines()      
        
        # pdb.set_trace()

        for txtbufr in flist[:]:

            rdafname = txtbufr.strip("\n")

            cmdbufr = "python " + srcdir + "RDA_HCA_dev.py --rdafname=" + rdafname  
            # print(cmdbufr)
            # pdb.set_trace()
            # Adding parallel jobs into pool       #! pay attention (cmdbufr, ) "," can not be ignored

            pool.apply_async(rda_cmd_exec, (cmdbufr, ))

    pool.close()

    pool.join()
