# MT1D_CNN
1D inversion study with CNN model.  (To be continued..) 

# 1D MT forward modeling with Occam1D

Download link:

* The aim of this inversion study is to training a model for retrive subsurface resistivity distribution using deep learning.    

### Requirements: 

* Install python 3.6/newer and anaconda packages: conda install keras, tensorflow, scikit-image, tqdm, pandas, numpy, seaborn, MTpy libraries.


### Dataset: 

* The original dataset of apparent resistivity and phase are simulated and saved as input data. The loaded data are used for demonstration purpose, the default frequency range and size is 128. During the training, the dataset are split into training and validation groups. The file structure is:


### Instruction:

   1. The two main python codes and /itools need be save in one folder.

   2. _forWeb_DataGen_1.py_ are used to training dataset generation, and are saved automatically after validation. Some parameters at the beginning of code are default, user can change it as wish.

   3. _forWeb_Training_.py_ are used to training model, and used it to predict the features of resistivity. All of the data should be saved in one folder, code will run through each of the data and save all of the results in new folder individually.


### Hardware used: 

* HP-7920 workstation: 56core CPU; 64G memory; GPU Nvidia Quadro P5000.


### Contact: 

* jon.liu@canada.ca


### Acknowledgments

* MT1D_CNN used open source codes and library from github, google, kaggle, and geophysical inversion packages.
* 
