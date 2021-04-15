# Clearup version of data generation:
# Install open source for mtpy:  https://mtpy2.readthedocs.io/en/develop/index.html
# Need download open-source package Occam1D and set the running path.
# If found any bugs, please contact jon.liu@canada.ca
#
#-----------------------------------------------------------------------------------------------------------------------
import os
import shutil
import subprocess
import matplotlib.pyplot as plt
import mtpy.modeling.occam1d as occam1d

from tqdm import tqdm_notebook
from mtpy.core.mt import MT
from itools import forWeb_cnn_1D as un
occam1d_path = '~/Occam1DCSEM/Source/OCCAM1DCSEM'  # path to Occam1D open source
#-----------------------------------------------------------------------------------------------------------------------
edi_file = '~/edi/plc002.edi'
edi_path = '~/edi'
savepath = os.getcwd()
fold_train = 'train_unet'
mt_obj = MT(edi_file)

n_layers = 49
input_size = 128

#-----------------------------------------------------------------------------------------------------------------------
savepath = os.path.join(savepath, fold_train)
if True:   #create the training dataset for cnn model
    n_sample = 50000
    if not os.path.exists(savepath):
        os.mkdir(savepath)
    savepath_m = os.path.join(savepath,'mask') #save model
    if not os.path.exists(savepath_m):
        os.mkdir(savepath_m)
    savepath_d = os.path.join(savepath,'data') #save app resisitivity
    if not os.path.exists(savepath_d):
        os.mkdir(savepath_d)

    ocd = occam1d.Data()
    ocd.write_data_file(edi_file=edi_file, mode='det',save_path=savepath, res_errorfloor=4,phase_errorfloor=2, remove_outofquadrant=True) 
    ocm = occam1d.Model(n_layers=n_layers, target_depth=50000,bottom_layer=110000,z1_layer=20)
    ocm.write_model_file(save_path=savepath)
    ocs = occam1d.Startup(data_fn=ocd.data_fn,model_fn=ocm.model_fn,max_iter=1,target_rms=0.0)
    #ocs.write_startup_file(save_path=savepath,)

    for i in range(0,n_sample):   #generate model and forward
        #generate resistivity.
        res_i = un.res_generator_0(ocm.num_params)
        save_startup_fn = f'OccamStartup1D_{i}'
        save_appres_fn = f'occam1D_{i}'
        ocs.write_startup_file_Jon(save_path=savepath, save_fn=save_startup_fn, res_i=res_i)
        #modeling
        subprocess.Popen([occam1d_path, '-F', save_startup_fn, save_appres_fn], cwd=savepath).wait()    #
        #save dataset.
        shutil.move(os.path.join(savepath,save_startup_fn), os.path.join(savepath_m,save_startup_fn))   # move model to subfolder
        shutil.move(os.path.join(savepath,save_appres_fn+'_0.resp'), os.path.join(savepath_d,save_appres_fn+'_0.resp'))  # move data to subfolder
    print(f'check: {i} is done. ')
