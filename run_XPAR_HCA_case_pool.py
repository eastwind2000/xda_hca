
import os

import multiprocessing 

import pdb

# ===============================================
# adding commenets and git version control
# sudo mount.cifs -t cifs  //10.10.31.21/data/rda/figs  ./remote_figdir/  -o user=root,pass="nmc601)(*"
# ===============================================


def met_cmd_exec( cmd  ):
    print(cmd)
    os.system(cmd)
    return

if __name__ == "__main__":

    nproc = 1

    workdir = "/mnt/f/fcst2022/scsmex2022/radarpro/" 

    rdasite_list= [  "ZZH03"  ]

    pool = multiprocessing.Pool(processes = nproc)    # setting up  parallel processes using pool  

    for rdasite_id in rdasite_list:

        flist = os.popen("ls -1 /mnt/f/fcst2022/scsmex2022/radarpro/data/XPAR_L2_scsmex2022/" + rdasite_id+"/*" ).readlines()

        # for txtbufr in flist[360::5]:
        #     print(txtbufr)    
        # pdb.set_trace()

        for txtbufr in flist[960:960+240]:

            rdafname = txtbufr.strip("\n")

            # pdb.set_trace()

            cmdbufr = "python " + workdir + "./src/XPAR_HCA_dev.py --rdafname=" + rdafname + " --figroot=" + workdir+ "./remote_figs/"   #! pay attention to last splash "/"

            # print(cmdbufr)

            pool.apply_async( met_cmd_exec, (cmdbufr, ) )  # Adding parallel jobs into pool            #! pay attention (cmdbufr, ) "," can not be ignored

    pool.close()
    
    pool.join()

