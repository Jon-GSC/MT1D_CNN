# MT1D_CNN
1D inversion study with CNN model.

# 1D MT forward modeling with Occam1D

### Download link-note:

* The aim of this inversion study is to train a model for retriving subsurface resistivity distribution using deep learning framework.    

### Requirements: 

* Install python 3.6/newer and anaconda packages: conda install keras-2.4.3, tensorflow-2.3.1, scipy, scikit-image, tqdm, pandas, numpy, seaborn, shutil, matplotlib, MTpy etc. libraries.

### Dataset: 

* The original dataset of apparent resistivity and phase are simulated and saved as input data. The loaded data are used for demonstration purpose, the default frequency range and size is 128. During the training, the dataset are split into training and validation sets. 


### Instruction:

   1. The three main python codes and /itools need be save in one folder. Before run the code: a). Install open-source mtpy: https://pypi.org/project/mtpy/. b). Download and setup open-source MT-Occam1D: https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/index.html, the enclosed Occam1DCSEM in /itools has been tested in Ubuntu operating system. User might need compile the Occam Fortran source for different operating environments.

   2. _forWeb1_Data_generate.py_ is used to generate the training dataset, and saved automatically in two folders. Some parameters at the beginning of code are default, user can change it as wish.

   3. _forWeb2_Model_training.py_ is used to create the model, and predict the resistivity of evaluation set. All of the data should be saved in one folder, code will run through each of the data and save all of the results in same folder.

   4. _forWeb3_Prediction.py_ is used for predicting of real '.edi' data, and plot figures of result.
   

### Hardware tested: 

* HP-7920 workstation: 56core CPU; 64G memory; GPU Nvidia Quadro P5000.


### Contact: 

* jon.liu@canada.ca


### Acknowledgments

* MT1D_CNN used open source codes and library from github, google, kaggle, and open-sourced geophysical inversion packages mtpy and occam1d.
 
