# Clearup concised version of cnn-inversion for MT data, 2021-04-09
# install mtpy, https://mtpy2.readthedocs.io/en/develop/index.html
# Need download open-source package Occam1D and set the running path.
# Keras2.4.3 is used. Please contact jon.liu@canada.ca if there is any bugs.
#
#-----------------------------------------------------------------------------------------------------------------------
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mtpy.modeling.occam1d as occam1d
from tqdm import tqdm_notebook
from mtpy.core.mt import MT
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping, ModelCheckpoint, ReduceLROnPlateau
from keras.models import load_model
from keras.layers.merge import concatenate
from keras.engine.training import Model
from keras.layers import Input,Flatten,Dense,Reshape
from keras import optimizers
from itools import forWeb_cnn_1D as un

occam1d_path = '~/Occam1DCSEM/Source/OCCAM1DCSEM'  # change to your path of Occam1D open source
savepath = os.getcwd()
#-----------------------------------------------------------------------------------------------------------------------
edi_file = '~/edi/plc002.edi'
edi_path = '~/edi'
fold_train = 'train_unet'
mt_obj = MT(edi_file)

n_layers = 49   #input layers
input_size = 128  #number of frequency
epochs = 10
batch_size = 16
version = 0          # test version
iflag = [ 1, 2, 3  ]

#-----------------------------------------------------------------------------------------------------------------------
t_start = time.time()
basic_name = f'Unet_version_{version}'
save_model_name = basic_name + '.model'
savepath = os.path.join(savepath, fold_train)
output_size = n_layers+1
frequency = mt_obj.Z.freq

if 1 in iflag:  #load data to matrix..
    train_fold = savepath
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
    train_mask = np.array([np.array(train_df['mask'][i][0:]) for i in range(len(train_df))])
    ids_train, ids_valid, x_train0, x_valid0, x_train1, x_valid1, y_train, y_valid = train_test_split(
        train_df.index.values, train_data0_N.reshape(-1, input_size, 1), train_data1_N.reshape(-1, input_size, 1),
        train_mask.reshape(-1, output_size, 1), test_size=0.2, random_state=1588)

if 2 in iflag:    #training.
    # training 1.
    input_layer0 = Input(shape=(input_size, 1))
    input_layer1 = Input(shape=(input_size, 1))
    model0 = un.build_model(input_layer0, output_size, 32, 3, 0.5)
    model1 = un.build_model(input_layer1, output_size, 32, 3, 0.5)
    # model1.summary()
    combined = concatenate([model0.output, model1.output])
    uconv0 = Dense(output_size, activation='relu', name='predictions')(combined)
    uconv0 = Reshape((uconv0.shape[1], 1))(uconv0)
    model = Model([model0.input, model1.input], uconv0)
    c = optimizers.Adam(lr=0.01)
    model.compile(loss=un.rmse_m, optimizer=c, metrics=[un.rmse_m])
    early_stopping = EarlyStopping(monitor='val_loss', mode='min', patience=50, verbose=1)
    model_checkpoint = ModelCheckpoint(save_model_name, monitor='val_loss',
                                       mode='min', save_best_only=True, verbose=1)
    reduce_lr = ReduceLROnPlateau(monitor='val_loss', mode='min', factor=0.8, patience=3, min_lr=0.0001, verbose=1)

    history = model.fit([x_train0, x_train1], y_train,
                        validation_data=([x_valid0, x_valid1], y_valid),
                        epochs=epochs,
                        batch_size=batch_size,
                        callbacks=[model_checkpoint, reduce_lr, early_stopping],
                        verbose=2)
    pd.DataFrame(history.history).to_csv(f'TrainHistory1_{version}.csv')  # save history1
    # Training 2

if 3 in iflag:   # prediction
    model = load_model(save_model_name, custom_objects={'rmse_m':un.rmse_m})  #
    preds_valid = un.predict_result(model,[x_valid0,x_valid1],output_size)
    offset = 0
    max_images = 32
    grid_width = 8
    Depth = np.arange(len(y_valid[0]))
    freq0 = np.logspace(np.log10(frequency[0]), np.log10(frequency[-1]), input_size)
    grid_height = int(max_images / grid_width)

    fig1, axs1 = plt.subplots(grid_height, grid_width, figsize=(2 * grid_width, 2 * grid_height))  #plot model
    for i, idx in enumerate(ids_valid[offset:offset+max_images]):
        plot_model = y_valid[i + offset]
        plot_pred = preds_valid[i + offset]
        ax1 = axs1[int(i / grid_width), i % grid_width]
        ax1.plot(plot_model,Depth, 'g', linewidth=0.8)
        ax1.plot(plot_pred,Depth, 'r', linewidth=0.8)
        ax1.set_title(str(ids_valid[i+offset]), color='red')
        ax1.set_yticklabels([])
        ax1.invert_yaxis()
        ax1.grid(linestyle='-')
    plt.savefig(os.path.join(savepath,f'plotM_{ids_valid[i+offset]}.png'))

    fig2, axs2 = plt.subplots(grid_height, grid_width, figsize=(2 * grid_width, 2 * grid_height))
    for i, idx in enumerate(ids_valid[offset:offset+max_images]):
        ax2 = axs2[int(i / grid_width), i % grid_width]
        ax2.loglog(train_data0[idx],freq0, 'g', linewidth=0.8)
        ax2.set_title(str(idx), color='red')  # .z
        ax2.grid(linestyle='-')
    plt.savefig(os.path.join(savepath,f'plotD_{idx}.png'))

    fig3, axs3 = plt.subplots(grid_height, grid_width, figsize=(2 * grid_width, 2 * grid_height))
    for i, idx in enumerate(ids_valid[offset:offset+max_images]):
        ax3 = axs3[int(i / grid_width), i % grid_width]
        ax3.plot(train_data0_N[idx],freq0, 'g', linewidth=0.8)
        ax3.plot(train_data1_N[idx],freq0, 'b', linewidth=0.8)
        ax3.set_yscale('log')
        ax3.set_title(str(idx), color='blue')
        ax3.grid(linestyle='-')
    plt.savefig(os.path.join(savepath,f'plot_nD_{idx}.png'))

    fig4, axs4 = plt.subplots(grid_height, grid_width, figsize=(2 * grid_width, 2 * grid_height))
    ocm = occam1d.Model()
    ocm.read_model_file(os.path.join(savepath,'Model1D'))
    for i, idx in enumerate(ids_valid[offset:offset+max_images]):
        ax4 = axs4[int(i / grid_width), i % grid_width]
        plot_depth = ocm.model_depth[1:] / 1
        plot_model4 = y_valid[i+offset]
        plot_pred4 =  preds_valid[i+offset]
        line4_1 = ax4.plot(plot_model4[::-1], np.log10(plot_depth[::-1]), color='g', lw=0.8)
        line4_2 = ax4.plot(plot_pred4[::-1], np.log10(plot_depth[::-1]), color='r', lw=0.8)
        line4_1[0].set_drawstyle('steps')
        line4_2[0].set_drawstyle('steps')
        line4_2[0].set_linestyle('--')
        ax4.set_title(str(ids_valid[i+offset]), color='red')
        ax4.invert_yaxis()
        ax4.grid(linestyle='-')
    plt.savefig(os.path.join(savepath,f'plotMzig_{ids_valid[i+offset]}.png'))
    plt.show()
