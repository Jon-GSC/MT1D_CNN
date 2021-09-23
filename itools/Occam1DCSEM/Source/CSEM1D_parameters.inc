!
! CSEM1D_parameters.inc
!
! Stores constants and magic numbers for 1D CSEM inversion
!

    real(8), parameter    :: pi = 3.141592653589793d0
    real(8), parameter    :: rad2deg = 57.295779513082323d0
    real(8), parameter    :: deg2rad = 0.017453292519943d0
!
! SVD Method for OPRA:
!
    character(24), parameter :: svdMethod = 'lapack_svd'   ! options are 'lapack_svd' and 'parker_svd'
!
! Data Type Index Dictionary
!
! These are the codes used in the data file parameter table.
! 
! Format:  EMData_1.0
!
! CSEM Data:
! 
! This format is compatible with 2D and 3D CSEM modeling.
! 
!   CSEM Electric Data, Real/Imag:
    integer, parameter :: indRealEx = 1
    integer, parameter :: indImagEx = 2
    integer, parameter :: indRealEy = 3
    integer, parameter :: indImagEy = 4
    integer, parameter :: indRealEz = 5
    integer, parameter :: indImagEz = 6

    
!   CSEM Magnetic Data, Real/Imag:
    integer, parameter :: indRealBx = 11
    integer, parameter :: indImagBx = 12
    integer, parameter :: indRealBy = 13
    integer, parameter :: indImagBy = 14
    integer, parameter :: indRealBz = 15
    integer, parameter :: indImagBz = 16
    
!   CSEM Electric Data, Amp/Phs:
    integer, parameter :: indAmpEx = 21
    integer, parameter :: indPhsEx = 22
    integer, parameter :: indAmpEy = 23
    integer, parameter :: indPhsEy = 24
    integer, parameter :: indAmpEz = 25
    integer, parameter :: indPhsEz = 26

    
!   CSEM Magnetic Data, Amp/Phs:
    integer, parameter :: indAmpBx = 31
    integer, parameter :: indPhsBx = 32
    integer, parameter :: indAmpBy = 33
    integer, parameter :: indPhsBy = 34
    integer, parameter :: indAmpBz = 35
    integer, parameter :: indPhsBz = 36
    



!   CSEM Horizontal Polarization Ellipses:
    integer, parameter :: iPEmax = 41
    integer, parameter :: iPEmin = 42
    
    integer, parameter :: iPBmax = 43
    integer, parameter :: iPBmin = 44
    


! To do: E/E, E/H transfer functions, etc

!   
! Magnetotelluric Data:  
!
! This format is compatible with 2D and 3D MT modeling codes.
!
!   Apparent Resistivity and Phase:
    integer, parameter :: indRhoZXX  = 101   ! not used for 1D
    integer, parameter :: indPhsZXX  = 102   ! not used for 1D
    integer, parameter :: indRhoZXY  = 103   ! 
    integer, parameter :: indPhsZXY  = 104   !  
    integer, parameter :: indRhoZYX  = 105   !
    integer, parameter :: indPhsZYX  = 106   ! 
    integer, parameter :: indRhoZYY  = 107   ! not used for 1D
    integer, parameter :: indPhsZYY  = 108   ! not used for 1D
    
!   Impedances: 
    integer, parameter :: indRealZXX  = 111   ! not used for 1D
    integer, parameter :: indImagZXX  = 112   ! not used for 1D
    integer, parameter :: indRealZXY  = 113   ! 
    integer, parameter :: indImagZXY  = 114   ! 
    integer, parameter :: indRealZYX  = 115   !
    integer, parameter :: indImagZYX  = 116   ! 
    integer, parameter :: indRealZYY  = 117   ! not used for 1D
    integer, parameter :: indImagZYY  = 118   ! not used for 1D 

! To do:  other transfer functions (B tipper, marine Jz tipper, interstation E-E, H-H and hybrid E-H), etc

