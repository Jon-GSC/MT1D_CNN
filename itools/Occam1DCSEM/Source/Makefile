# 
# Makefile for Occam1DCSEM
# 
# Kerry Key
# Scripps Institution of Oceanography
# 
# 

#----------------------------------------------------
# Step 1: Specify your Fortran compiler (FC) and the 
# compiler optimization settings (FCFLAGS):
#----------------------------------------------------

# 
# Intel Fortran compiler:
# 
 FC = ifort 

# Use this for faster runtime: 
 #FCFLAGS = -O2 -m64 -mmacosx-version-min=10.4 
 
# use   -mmacosx-version-min=10.4   to create a portable executable
 
# Use this for debugging source code:
# FCFLAGS = -stand f03 -g  -warn all -fstack-security-check -check all

# 
# GFortran Compiler (its free!):
# 
 FC = gfortran
# 
# Use this for faster runtime: 
# FCFLAGS =  -O2  -m64

# Use this for debugging source code:
# FCFLAGS =    -Wall  -pedantic -std=f2003
  
#----------------------------------------------------
# OPTIONAL for compiling the mexDipole1D.f90 Matlab interface:
# 
# MATLIBPATH = Path to the MatLab libraries used for 
#              building the MatLab interfaces.
#       You also need to change the path hardcoded
#       in mexDipole1D.f90!
# 
# MEXBIN,MEXINC =  Path to specific folders in the Matlab program. 
# 
# 
# MATEXT = File extension required by MatLab to
#          recognize the MatLab mex function.
#          Platform dependent.
# 
# Note that the MEXFLAGS variable needs either -cpp (for ifort)
# or  -x f95-cpp-input (for gfortran) to preprocess the include 
# directive in mexDipole1D.f90.
# 
#----------------------------------------------------

MATLABPATH = /Applications/MATLAB_R2010b.app
#MEXBIN = $(MATLABPATH)/bin/maci/         # 32 bit
MEXBIN = $(MATLABPATH)/bin/maci64/		 # 64 bit
MEXINC = $(MATLABPATH)/extern/include/

# don't forget to add the MATLABPATH to the include statement in mexDipole1D.f90


MATEXT = mexmaci64   # 64 bit Mac Matlab
# MATEXT = mexmaci   # 32 bit Mac Matlab
# MATEXT = dll       # use this on PC

# MEXFLAGS =                   #  use nothing if you don't want the mex files.
MEXFLAGS = -x f95-cpp-input  # use this with  gfortran 
#MEXFLAGS = -cpp 	      # use this with ifort

#----------------------------------------------------
# Step 2: Specify the location of your LAPACK and BLAS
# libraries.  These contain the linear algebra routines
# used by CSEM1D and Occam.  If you only want to compile
# Dipole1D, you don't need these libraries.
#----------------------------------------------------

# 
# LAPACK and BLAS codes:
# 
#  Download the latest Lapack library from netlib.org.
# Compile the lapack using this:
# "cd" to ./lapack-3.2.1 and make sure the file 
# make.inc lists the correct Fortran compiler for your machine.
# Then type "make blaslib" and then "make lapacklib".  
# This makes the libraries libblas.a and liblapack.a, 
# which are linked in the following command:
#  
#LIBS =  -L./lapack-3.2.1  -lblas -llapack
 LIBS = -L/home/jliu/My_Projects/Geothermal_MT/myCodes/lapack-3.2.1 -lblas -llapack
 
# 
# On Mac OS X platforms, install Apple's developer tools. Then 
# all you need to get lapack and blas is this
# command:
# 
 #LIBS = -framework veclib
# 
# You're done.  You shouldn't need to change anything below here.
# 
#----------------------------------------------------
# Build commands:
#----------------------------------------------------

TARGETS=  clean OCCAM1DCSEM DIPOLE1D   # mexDipole1D

OBJSOC= FilterModules.o Occam.o  \
		Dipole1D.o  CSEM1D.o   RunOccam.o 
		
OBJSDP= FilterModules.o Dipole1D.o CallDipole1D.o
		
OBJSMEX= FilterModules.o Dipole1D.o mexDipole1D.o
MEXlibs = -L$(MEXBIN) -I$(MEXINC)  -lmx -lmex -lmat \

		
all:  $(TARGETS)
		
clean:	clean_msg
		rm -f *.o *~ core *.mod
		rm -f OCCAM1DCSEM
		rm -f DIPOLE1D
		rm -f *.$(MATEXT)


OCCAM1DCSEM: build_msg_occ $(OBJSOC)
		$(FC) $(FCFLAGS) -o $@ $(OBJSOC)  $(LIBS)

DIPOLE1D:  build_msg_dp $(OBJSDP)
		$(FC) $(FCFLAGS) -o $@ $(OBJSDP) 

mexDipole1D: build_msg_dpmx $(OBJSMEX)
		$(FC) $(FCFLAGS)  $(MEXlibs) -DMATLAB_MEX_FILE -undefined dynamic_lookup -bundle -o $@.$(MATEXT) $(OBJSMEX)
 
	
#		
# Compile rules
#		

# General Fortran compile:
%.o: %.f90 
	$(FC) $(FCFLAGS)    -c -o $@ $^
	
# Fortran compile with Matlab Mex libraries:	
$(OBJSMEX): %.o: %.f90 
	$(FC) $(FCFLAGS) $(MEXFLAGS)  -c -o $@ $^
	
#	
# Build Messages:
#	
clean_msg: 
	@printf "#\n# Cleaning files: \n#\n"
	
build_msg_occ: 
	@printf "#\n# Building OCCAM1DCSEM: \n#\n"
	
build_msg_dp: 
	@printf "#\n# Building DIPOLE1D: \n#\n"

build_msg_dpmx: 
	@printf "#\n# Building mexDipole1D: \n#\n"	
