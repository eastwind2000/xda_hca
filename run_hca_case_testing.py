
import os

import multiprocessing 

import pdb

# adding commenets and git version control

# sudo mount.cifs -t cifs   //10.10.31.21/data/rda/figs  ./remote_figdir/  -o user=root,pass="nmc601)(*"

# ===============================================

def met_cmd_exec( cmd  ):
    print(cmd)
    os.system(cmd)
    return

if __name__ == "__main__":

    nproc = 4

    workdir = "/mnt/f/fcst2022/scsmex2022/radarpro/" 

    # rdasite_list= [   "Z9758", "Z9762", "Z9200", "Z9662",  "Z9755"  ]  #! 2022 SCSMEX GBA         
    # rdasite_list= [  "Z9518", "Z9527", "Z9516" ] # 20220720 JiangSu Province Tornadoes Case
    # rdasite_list= [ "Z9370" ] # 2022072200-2300UTC shangqiu Tornadoes Case
    # rdasite_list= [ "Z9710"  ] # 20210812 SuiZhou, HuBei Extreme Precipitaion Case

    # rdasite_list= [ "Z9971"   ] # 20220818 Qinghai DaTong Extreme precipitation Case
    
    rdasite_list= [ "Z9971"   ] # 20220818 Qinghai DaTong Extreme precipitation Case

    pool = multiprocessing.Pool(processes = nproc)    # setting up  parallel processes using pool  

    for rdasite_id in rdasite_list:

        flist = os.popen("ls -1 /mnt/f/fcst2022/scsmex2022/data/radar_L2_scsmex2022/" + rdasite_id + "/*" ).readlines()
        
        # flist = os.popen("ls -1 /mnt/f/fcst2022/Tornado_IOP/rda/" + rdasite_id+"/*" ).readlines()        
        # flist = os.popen("ls -1 /mnt/f/fcst2022/Mesovortex_IOP/data/rda/" + rdasite_id + "/*" ).readlines()
        
        # for txtbufr in flist[360::5]:
        #     print(txtbufr)    
        # pdb.set_trace()

        for txtbufr in flist[131:144]:

            rdafname = txtbufr.strip("\n")

            cmdbufr = "python " + workdir + "./src/RDA_HCA_dev.py --rdafname=" + rdafname + " --figroot=" + workdir + "./remote_figs/"   #! pay attention to last splash "/"

            # print(cmdbufr)
            # pdb.set_trace()

            pool.apply_async( met_cmd_exec, (cmdbufr, ) )  # Adding parallel jobs into pool       #! pay attention (cmdbufr, ) "," can not be ignored

    pool.close()
    
    pool.join()

