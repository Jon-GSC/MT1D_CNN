# MT1D inversion study with CNN model.

### Download link-note:

* The aim of this inversion study is to train a neural network model for retriving subsurface resistivity distribution using deep learning framework. The code has been tested without any issue in PyCharm IDE.   

### Requirements: 

* Install python 3.6/newer and anaconda packages: conda install keras-2.4.3, tensorflow-2.3.1, scipy, scikit-image, tqdm, pandas, numpy, seaborn, shutil, matplotlib, MTpy etc. libraries.

### Dataset: 

* The original dataset of apparent resistivity and phase are simulated and saved as input data in folder _train_unet_. The loaded data are used for demonstration purpose, the default frequency range and size is 128. During the training, the dataset are split into training and validation sets. 


### Instruction: 
     
   1. The three main python codes and _/itools_, _/edi_data_ need be save in same folder. Before run the code: a). Install open-source mtpy: https://pypi.org/project/mtpy/. b). Download and setup open-source MT-Occam1D: https://marineemlab.ucsd.edu/Projects/Occam/1DCSEM/index.html (optional), the enclosed Occam1DCSEM EXE in _/itools_ has been tested in Ubuntu operating system. User might need compile the Occam Fortran source for different operating environments.

   2. _forWeb1_Data_generate.py_ is used to generate the training dataset, and saved automatically in two folders: _data_ and _mask_. Some parameters at the beginning of code are default, user can change it as wish.

   3. _forWeb2_Model_training.py_ is used to create the model, and predict the resistivity of evaluation set. All of the data should be saved in one folder, code will run through each of the data and save all of the results in same folder.

   4. _forWeb3_Prediction.py_ is used for predicting of real '.edi' data, and plotting figures. The 'plc002.edi' in _edi_data_ is real data for testing purpose. User can add '.edi' files to the folder _edi_data_ without other change. The output will be save in folder _save_plot_0_, the predicted resistivity is saved in '_All_resistivity.csv_'.

   5. If user wants to validate the code, one can create some synthetic models and save the forward results as '.edi' format, then load into model by operating the step 4.
   
   For pilot running, user does not need change any parameters, just run the code files in steps 2,3,4. The code should run without any error if environment setting is correct. Please contact at following email address if any bugs popup.


### Hardware tested: 

* HP-7920 workstation: 56core CPU; 64G memory; one Nvidia Quadro P5000 GPU.


### Contact: 

* jon.liu@nrcan-rncan.gc.ca


### License

Unless otherwise noted, the source code of this project is covered under Crown Copyright, Government of Canada, and is distributed under the MIT Licence

The Canada wordmark and related graphics associated with this distribution are protected under trademark law and copyright law. No permission is granted to use them outside the parameters of the Government of Canada's corporate identity program. For more information, see Federal identity requirements.


### Acknowledgments:

* MT1D_CNN used open source codes and library from github, google, kaggle, and open-sourced geophysical inversion packages mtpy and occam1d. Please cite the related references in your publications.
 
### Citation: 

* Liu, X.; Craven, J.A.; Tschirhart, V. Retrieval of Subsurface Resistivity from Magnetotelluric Data Using a Deep-Learning-Based Inversion Technique. Minerals 2023, 13, 461. https://doi.org/10.3390/min13040461.

  
![plotMzig_94](https://user-images.githubusercontent.com/39324742/134573289-c75ca6e1-f622-4fe1-8d02-2e4cfbf67c3a.png)

