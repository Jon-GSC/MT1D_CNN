# this part is an example to load model for real 'edi' data predicting, and plot figures of result.
# If found any bugs, please contact jon.liu@canada.ca
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mtpy.modeling.occam1d as occam1d
from tqdm import tqdm_notebook
from mtpy.core.mt import MT
from scipy.interpolate import interp1d
from sklearn.model_selection import train_test_split
from keras.models import load_model
plt.rcParams.update({'font.size': 17})
from itools import forWeb_cnn_1D as un
current_path = os.getcwd()
fold_train = 'train_unet'
folder_data = 'data'
folder_mask = 'mask'
occam1d_path = os.path.join(current_path,'itools/Occam1DCSEM/Source/OCCAM1DCSEM')  # path to Occam1D
#-----------------------------------------------------------------------------------------------------------------------

n_layers = 49  # number of layers
input_size = 128  # size of input

version = 1   # model version number

edi_folder = 'edi_data/'  # real data file, user can add edi file to list.
save_plot_0 = 'save_plot_0'   # output folder

#-----------------------------------------------------------------------------------------------------------------------
t_start = time.time()
basic_name = f'CNN_model_v{version}'
save_model_name = basic_name + '.model'
train_fold = os.path.join(current_path, fold_train)
slst = [edi for edi in sorted(os.listdir(os.path.join(current_path,edi_folder))) if edi.find('.edi')>0]
output_size = n_layers+1
freq0 = un.read_resp_file_0(input_size, train_fold + '/data/occam1D_0_0.resp')[2]

savepath0 = os.path.join(current_path, save_plot_0)
if not os.path.exists(savepath0):
    os.mkdir(savepath0)
if True:
    train_df = pd.DataFrame(columns=[])
    filelist = []
    for file in sorted([f for f in os.listdir(train_fold + '/' + 'mask' + '/')],
                       key=lambda x: int(''.join(filter(str.isdigit, x)))):
        filelist.append(os.path.splitext(file)[0])
    train_df['mask_name'] = filelist
    filelist = []
    for file in sorted([f for f in os.listdir(train_fold + '/' + 'data' + '/')],
                       key=lambda x: int(''.join(filter(str.isdigit, x)))):
        filelist.append(os.path.splitext(file)[0])
    train_df['data_name'] = filelist
    train_df["data0"] = [un.read_resp_file_0(input_size, train_fold + '/' + 'data' + '/{}.resp'.format(idx))[0] for
                         idx in tqdm_notebook(train_df.data_name)]
    train_df["data1"] = [un.read_resp_file_0(input_size, train_fold + '/' + 'data' + '/{}.resp'.format(idx))[1] for
                         idx in tqdm_notebook(train_df.data_name)]
    train_df["mask"] = [un.read_startup_0(train_fold + '/' + 'mask' + '/{}'.format(idx)) for idx in
                        tqdm_notebook(train_df.mask_name)]
    train_data0 = np.array([np.array(train_df['data0'][i][0:]) for i in range(len(train_df))])
    data_mean0 = np.mean(train_data0)
    data_dev0 = np.std(train_data0)
    train_data0_N = (train_data0 - data_mean0) / data_dev0
    train_data1 = np.array([np.array(train_df['data1'][i][0:]) for i in range(len(train_df))])
    data_mean1 = np.mean(train_data1)
    data_dev1 = np.std(train_data1)
    train_data1_N = (train_data1 - data_mean1) / data_dev1

if True:
    model = load_model(save_model_name, custom_objects={'rmse_m': un.rmse_m})
    slst0 = [s0[:-4] for s0 in slst]
    ocm = occam1d.Model()
    ocm.read_model_file(os.path.join(train_fold, 'Model1D'))
    plot_depth = ocm.model_depth[1:].copy();  plot_depth[plot_depth==0]=1
    res_all = pd.DataFrame(columns=[])
    res_all['depth(m)'] = ocm.model_depth[1:]
    for idx,fn in enumerate(slst):
        mt_obj = MT(os.path.join(current_path,edi_folder+fn))
        freq = mt_obj.Z.freq
        res_det0 = mt_obj.Z.res_det
        phase_det0 = un.phase_shift_1(mt_obj.Z.phase_det)
        f1d0 = interp1d(freq, res_det0, kind='linear', fill_value='extrapolate')
        res_det0 = f1d0(freq0)
        f1d1 = interp1d(freq, phase_det0, kind='linear', fill_value='extrapolate')
        phase_det0 = f1d1(freq0)
        res_det = un.savgol_smooth_Jon_3(res_det0).reshape(-1,input_size,1)
        phase_det = un.savgol_smooth_Jon_3(phase_det0).reshape(-1,input_size,1)
        res_det = (res_det - data_mean0) / data_dev0
        phase_det = (phase_det - data_mean1) / data_dev1
        plot_pred = np.mean(un.predict_result(model,[res_det,phase_det],output_size),axis=0)#.reshape(-1)
        res_all[slst0[idx]] = plot_pred
        fig1, axs1 = plt.subplots(1, 1, figsize=(9, 15))
        line1_1 = axs1.plot(plot_pred[::-1], np.log10(plot_depth[::-1]), color='b', lw=1.5)
        line1_1[0].set_drawstyle('steps')    #
        line1_1[0].set_linestyle('-')
        axs1.set_xlabel(r'$log(\rho)$'+ ' $(\Omega m)$')
        axs1.set_ylabel('log(Depth) (m)')
        axs1.set_title(str(slst0[idx]), color='blue', fontsize=15)
        axs1.invert_yaxis()
        axs1.grid(linestyle='--')
        plt.savefig(os.path.join(savepath0, f'cnn_res_{slst0[idx]}.png'))
        plt.clf()
    res_all.to_csv(savepath0 + '/All_resistivity.csv')
