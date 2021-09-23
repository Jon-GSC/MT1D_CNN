!==============================================================================!    
!
! 1D CSEM/MT Occam's Inversion Code
!
! This file contains the interface routines that reside between the Dipole1D and
! Occam's inversion Fortran codes.
!
!------------------------------------------------------------------------------!
! License for Dipole1D and CSEM1D codes:
!------------------------------------------------------------------------------!
!
!    Copyright (C) 2007-2010
!    Kerry Key
!    Scripps Institution of Oceanography
!    kkey@ucsd.edu
!
!    This file is part of Dipole1D.
!
!    Dipole1D is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Dipole1D is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Dipole1D.  If not, see <http://www.gnu.org/licenses/>.
!
!
!------------------------------------------------------------------------------!
! Version numbers for the complete Occam 1D CSEM/MT Inversion Package:
!------------------------------------------------------------------------------!
!
! 3.13, March-Dec, 2010     DGM added several roughness types, merged various codes.
!
! 3.12, March 26, 2010      Minor tweaks for roughness Type parameter, added option for
!                           minimum gradient support regularization.
! 
! 3.11, Feb 05, 2010        David Myer added support for finite dipole length.
! 
! 3.10, January 21, 2010.   Fixed a minor bug affecting only joint CSEM/MT modeling.
!                           Thanks to DGM for identifying the problem.
!
! 3.9, October 20, 2009.    Minor changes, added parameter for selecting
!                           which svd routine to use since svd() occasionally gives wrong svd.
!
! 3.8, Aug 2009.  DGM       Lots of memory was NOT being deallocated before exit.
!                           Added provision in Occam for deallocateFwd subroutine
!                           and coded here.  Also rearranged some bits of code to allow
!                           for a MatLab interface to the forward and inversion codes.
!                           Since the MatLab interface is not a separate executable,
!                           memory must be very carefully managed.
!                                   
! 3.7, May 8, 2009.         A few minor changes to accomodate some Occam.f90 updates.
!
! 3.6, January, 2009.       Added Steve's polarization ellipse code for 
!                           horizontal E and B data.
!
! 3.5  December, 2008.      Fixed a few things that were not compliant with the
!                           Fortran2003 standard.  Updated model and data file 
!                           formats.  Add support for RealPrec parameter defined in 
!                           module Occam_RealPrecision.
! 
! 3.4  November, 2008.      Added in support for ascii data types in data 
!                           file parameter table. Model files now can 
!                           use ? for flagging inversion layers (or use -1). 
!                           Added support for phase data and option 
!                           for specifying Phase Lag or Phase Lead conventions.
! 
! 3.3  September 19, 2008.  Switched to R.L. Parker's svd routine for OPRA.
!
! 3.2  September 17, 2008.  Cleaned up a few things, fixed a bug in the 
!                           rotation code.
!
! 3.1   September 4, 2008.  Added Rx-Tx transformations, cleaned up 
!                           various sections of the code in preparation  
!                           for distribution.
!
! 3.0   August 2008.        Added orthogonal Procrustes method to estimate  
!                           3D receiver orientations, if requested.
!
! 2.0                       Added 1D MT for joint 1D CSEM-MT inversion.
!
! 1.0   March, 2008.        1D CSEM inversion using Dipole1D.f90.
!
!
!------------------------------------------------------------------------------!
! Authorship:
!------------------------------------------------------------------------------!
!
! Please reference both papers listed below when citing this code.  
!
! Occam.f90:
! 
!    The main Occam's inversion was originally written by Steven Constable 
!    and is described in the paper:
!  
!    Constable, S. C., R. L. Parker, and C. G. Constable, 1987, 
!    Occam's inversion - A practical algorithm for generating smooth models 
!    from electromagnetic sounding data: Geophysics, 52, 289 300.
!
!    See the notes at the top of Occam.f90 for more information and  
!    documentation of Occam's revision history.
!
!
! CSEM1D.f90 and Dipole1D.f90:
!
!    These routines perform the 1D CSEM forward and sensitivity computations. 
!    They were written in 2007-2009 by Kerry Key and are described in 
!    this manuscript:
!
!    Key, K., 2009, One-dimensional inversion of multi-component, 
!    multi-frequency marine CSEM data: Methodology and synthetic studies for 
!    resolving thin resistive layers: Geophysics, Vol.74,No.2, March-April 2009,
!    doi:  10.1190/1.3058434 
!
!
! Orthogonal Procrustes Rotation Analysis in subroutine computeCSEM_OPRA:
! 
!    The OPRA method for determining the orientation of marine CSEM receivers
!    is described this manuscript:
!
!    Key, K. and A. Lockwood, 2010, Determining the orientation of marine CSEM
!    receivers using orthogonal Procrustes rotation analysis, accepted 
!    to Geophysics.
!
!
!------------------------------------------------------------------------------!
! Acknowledgments:
!------------------------------------------------------------------------------!
!
! This work was supported by:
!
!    The Seafloor Electromagnetic Methods Consortium   
!       at Scripps Institution of Oceanography  
!    
! See this URL for a list of SEMC funding sponsors:
!
!    http://marineemlab.ucsd.edu/semc.html
!    
! 
!------------------------------------------------------------------------------!
! Description of major subroutines contained in this file:
!------------------------------------------------------------------------------!
!
! The first 7 routines listed here are called directly by Occam.f90, but are specific
! to the geophysical method and dimensionality (1D,2D,3D), so it's best to 
! put them here in an interface external to Occam.  See the notes in Occam.f90 for more
! specific details of what these routines need to do for Occam.f90.
!
! readData     - Reads in CSEM data file, allocates some data sized arrays.
!                Data and std errors are stored in vectors d and sd, listed in module Occam 
!                Data parameters are stored in array dp as a lookup table to other arrays (freq, tx location, etc)
!                Rx,Tx parameters are stored in module csem1d_mod since Occam doesn't need these.
!
! readModel    - Reads in the 1D model structure, allocates layers for fwd comps.
!                Occam doesn't use this stuff, so it is stored in module csem1d_mod
!
! computeFwd   - Sets up the calls to Dipole1D and MT1D for forward responses and model derivatives.
!
! constructPenaltyMat  - Constructs the penalty matrix for the 1D problem
!
! countPenaltyTerms - Counts the total number of penalty terms so that arrays can be allocated
!
! writeResponse -  Occam calls this routine at the end of each iteration.  This is where the model
!                  response are output to a file. If you want your own custom output files, 
!                  this is where to place that code.
!
! deallocateFwd - Deallocate all memory associated with the forward problem.  Useful when
!                 calling Occam from an external interface.
!
!
! These routines are called by computeFwd:
!
! computeFwd_CSEM - Computes 1D CSEM forward responses and model derivatives by calling Dipole1D
!
! computeFwd_MT   - Computes 1D MT forward responses and model derivatives by calling MT1D
!
! computeCSEM_OPRA  - Computes the orthogonal Procrustes estimate of the CSEM receiver's 3D rotation, if requested.
!                     Then rotates the data vector d and standard error vector sd accordingly.
!
! svd - An SVD routine from R.L. Parker that is optionally used for the OPRA method.
! 
!==============================================================================!    
!================================================================! CSEM1D_MOD  ! 
!==============================================================================!
    module csem1d_mod

    use Occam, only : realprec
    
    implicit none    
!
! Module to store various parameters for 1D CSEM inversion
! This is stuff that is not used in Occam.f90 and not Dipole1D.f90
! but is required for interfacing between them.

!   Magic codes and constants in this include
    include 'CSEM1D_parameters.inc'

! While the data parameter array dp is passed in/out of Occam, it is only used
! in the forward routines.  So we set cnDataParams in the forward routines
! and then allocate the dp array in readData.
!        
    integer :: cnDataParams
 
    real(8), dimension(:,:), allocatable :: layer_params 
!
! layer_params is [ top_depth, rho, penalty, preference, prej_penalty ]          
!
! top_depth :: Depth to the top of the layer (ignores first value since top of 1D model is infinity)
! rho       :: Either a fixed value or else -1, the flag for an inversion free parameter.
!              If the resistivity is fixed, then none of the remaining 3 layer parameters are used 

! penalty       :: The roughness penalty to apply across the layer top boundary. Use zero to allow for a jump
! preference     :: A linear resistivity preference value that you can use to push the inversion in a certain direction.  
!                  If you don't want a preference, set this and the prej_penalty to 0.
! prej_penalty  :: The roughness penalty to apply to the difference of the current model parameter value and the Preference.
!                  If you really want to force the inversion towards the preference value, make prej_penalty >> penalty for 
!                  the layer
!

!
! Store a copy of the d and sd vectors to use for the OPRA method:
!
    real(RealPrec), dimension(:), allocatable :: d_copy, sd_copy   ! using RealPrec from module Occam_realprecision
       
!
! Data dependent Tx-Rx variables read from data file in subroutine readData:
!
    integer                               :: nTxIn, nRxIn, nFreqIn
    real(8), dimension(:), allocatable    :: azimuthTxIn,dipTxIn
    real(8), dimension(:), allocatable    :: xTxIn,yTxIn,zTxIn,xRxIn,yRxIn,zRxIn
    real(8), dimension(:), allocatable    :: ThetaRxIn,AlphaRxIn, BetaRxIn
    real(8), dimension(:), allocatable    :: fTxIn, sigsite 
 
!
! Rx orientation solver variables:
!
    integer                               :: nRxRotSolve ! Number of Rx to solve for orientations
    logical, dimension(:), allocatable    :: lRxRotSolve ! [nRxIn] logical array of flags for solving for Rx rotations

!
! Transformed Tx-Rx arrays used to speed up computations (reduces nTx when dip, azimuth and depth are constant): 
!
    integer                               :: nTxTfm
    real(8), dimension(:), allocatable    :: xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm
    integer                               :: nRxTfm
    real(8), dimension(:), allocatable    :: xRxTfm, yRxTfm, zRxTfm
    integer, dimension(:), allocatable    :: iRxTfm ! index pointing from Tfm site number to input site number (RxIn)
    integer, dimension(:), allocatable    :: RxInd  ! index pointing from Tfm site number to site sent to Dipole1D
    integer, dimension(:), allocatable    :: rx1DtoRxIn ! index pointing from Dipole1D site # to RxIn number,used for rotations
    integer, dimension(:,:), allocatable  :: dpTfm  ! data parameter array for transformed Tx to Rx's


    end module csem1d_mod
    
!==============================================================================!    
!==================================================================! MT1D_MOD  ! 
!==============================================================================!
    module mt1d_mod
!
! Module for 1D MT response and layer derivatives:
!
    real(8) :: zRx_MT
    real(8) :: freqMT
    
    complex(8) :: CResp
    complex(8), dimension(:), allocatable :: dCRespdRho  
    
    
    end module mt1d_mod

!==============================================================================!  
!========================================================! dataTransformations !  
!==============================================================================!  
!
! Module for functions that perform various data transformations
!

    module dataTransformations
    
    use csem1d_mod
    
    implicit none
        
    real(8), parameter  :: mu = 4d-7*pi 
    
    contains 

!
! MT Apparent resistivity
!
        pure real(8) function getAppRes(z,freq)
        
          implicit none
          complex(8), intent(in)   :: z    ! impedance Z, not C
          real(8), intent(in)      :: freq ! linear frequency
        
          
          getAppRes = abs(z)**2  / (2*pi*freq*mu)
          
        end function

        elemental real(8) function getAppResDeriv(z,dz,freq)
        
          implicit none
          complex(8), intent(in)   :: z    ! impedance Z, not C
          complex(8), intent(in)   :: dz   ! derivative of impedance Z 
          real(8), intent(in)      :: freq ! linear frequency
          
          getAppResDeriv = 2.d0 /  (2*pi*freq*mu) * (real(z)*real(dz) + aimag(z)*aimag(dz))

        end function

!
! Phase with unwrapping:
!  
        pure real(8) function getPhase(z,dataphi)
        
        ! JH hacky method of getting around lack of modulo subtraction in Occam.
        
          implicit none
          complex(8), intent(in)   :: z
          real(8), intent(in)      :: dataphi
          real(8)                  :: tmp 
    
          tmp = atan2( aimag(z), dble(z) )*rad2deg  ! phase in degrees
          
          getPhase = tmp
          if( abs(tmp + 360.0 - dataphi) .le. abs(tmp - dataphi) ) getPhase = (tmp + 360.0)
          if( abs(tmp - 360.0 - dataphi) .le. abs(tmp - dataphi) ) getPhase = (tmp - 360.0)
          
    
        end function
        
        elemental real(8) function phaseDeriv(z,dz) 
        
          implicit none
          complex(8),intent(in) :: z, dz
          phaseDeriv = rad2deg * (real(z)*aimag(dz) - aimag(z)*real(dz) ) / abs(z)**2
          
        end function          
        
!
! Absolute value:
!
        elemental real(8) function absDeriv(phi,dphi) 
        
          implicit none
          complex(8),intent(in) :: phi, dphi
    
          absDeriv = ( dble(phi)*dble(dphi) + aimag(phi)*aimag(dphi))  / abs(phi)
          
        end function
     
!
! Polarization ellipse parameters:
!

        pure real(8) function getPE(e1,e2,comp)
        
      ! Taken from Steve's old Occam routines - does not match  Smith and Ward 
      ! but probably works OK.
        implicit none
        complex(8), intent(in)   :: e1    ! E/B-field channel 1
        complex(8), intent(in)   :: e2    ! E/B-field channel 2
        character(4), intent(in) :: comp  ! flag  'pmax' or 'pmin' for output 
        
      ! Local:
        real(8)                  :: x1, x2, y1, y2
        real(8)                  :: a, b, phi, s, c, p1, p2, pmin, pmax


      ! bust up the complex components of the electric field:
        x1 = real(e1)
        y1 = aimag(e1)
        x2 = real(e2)
        y2 = aimag(e2)
      
      ! find the critical angle
        a   = -(x1*y1 + x2*y2)
        b   = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
        phi = atan2(a,b)/2.
      
      ! find the two critical axes of the pe 
        s  = sin(phi)
        c  = cos(phi)
        p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p1 = sqrt(abs(p1))
        
        s  = sin(phi+pi/2.d0)
        c  = cos(phi+pi/2.d0)
        p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p2 = sqrt(abs(p2))

        pmax = max(p1,p2)
        pmin = min(p1,p2)
        
        select case (comp)
        
          case('pmin')
            getPE = pmin
            
          case('pmax')
            getPE = pmax
            
        end select
    

        end function  
        
        

        
        elemental real(8) function getPEDeriv(e1,e2, de1, de2,comp)
      ! Taken from Steve's old Occam routines - does not match Smith and Ward 
      ! but probably works OK
        implicit none
        complex(8), intent(in)   :: e1    ! E/B-field channel 1
        complex(8), intent(in)   :: e2    ! E/B-field channel 2
        complex(8), intent(in)   :: de1   ! E/B-field deriv channel 1
        complex(8), intent(in)   :: de2   ! E/B-field deriv channel 2
        character(4), intent(in) :: comp  ! flag  'pmax' or 'pmin' for output 
        
        ! Local:         
        real(8)                  :: x1, x2, y1, y2, dx1, dx2, dy1, dy2
        real(8)                  :: a, b, phi, s, c, p1, d1, p2, d2
        
      ! bust up the complex components of the electric field:
        x1  = real(e1)
        y1  = aimag(e1)
        x2  = real(e2)
        y2  = aimag(e2)
        dx1 = real(de1)
        dy1 = aimag(de1)
        dx2 = real(de2)
        dy2 = aimag(de2)
      
      ! find the critical angle
        a = -(x1*y1 + x2*y2)
        b = 0.5*(y1*y1 + y2*y2 -x1*x1 - x2*x2)
        phi = atan2(a,b)/2.
      
      ! find the two critical axes of the pe and derivatives
        s = sin(phi)
        c = cos(phi)
        p1 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p1 = sqrt(abs(p1))
        d1 = c*c*(x1*dx1+x2*dx2) + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
        d1 = d1/p1
        
        s = sin(phi+pi/2.d0)
        c = cos(phi+pi/2.d0)
        p2 = (x1*x1+x2*x2)*c*c + (y1*y1+y2*y2)*s*s + 2.*s*c*(x1*y1+x2*y2)
        p2 = sqrt(abs(p2))
        d2 = c*c*(x1*dx1+x2*dx2) + s*s*(y1*dy1+y2*dy2) + s*c*(x1*dy1+dx1*y1+x2*dy2+dx2*y2)
        d2 = d2/p2


        if (p1 >= p2) then
        
            select case (comp)
            
                case('pmin')
                getPEDeriv = d2
                
                case('pmax')
                getPEDeriv = d1
            
            end select
        
        else ! ( p1 < p2)
        
            select case (comp)
            
                case('pmin')
                getPEDeriv = d1
                
                case('pmax')
                getPEDeriv = d2
            
            end select
        
        end if
     
        
        end function  

        
    end module dataTransformations

!==============================================================================!  
!===================================================================! readData !  
!==============================================================================!    
    subroutine readData()
!
! Subroutine to read in an Occam1DCSEM data file and allocate data dependent
! arrays.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 2.2  November, 2008     Added in support for ascii data types (RealEx, etc)
!                                 in the input data parameter table.
! Version 2.1  September 24, 2008 Fixed a bug where iRxTfm was undefined.
! Version 2.0  September 2008     Added some error checks. Created Tx-Rx transformations 
!                                 to speed up computations when simple transmitters are modeled.
!                                 The speed of Dipole1D depends on the number of transmitters, so 
!                                 transforming many simple transmitters into a single transmitter 
!                                 can greatly speed up the 1D inversion.
!                                 If only the x,y position of the transmitter varies,
!                                 but not the azimuth, dip or z position, then this can be modeled as a single
!                                 transmitter with many pseudo-receivers added to the Rx array.  There are even
!                                 further transformations that can be added for variable azimuth Tx's, but 
!                                 I'll skip that for now.
!                                   
! Version 1.0  Spring 2008     
!
!
    use Occam
    use csem1d_mod
    use dipole1d
    
    implicit none
      
!
! Local variables:
! 
    integer         :: i, err, iTxRead, iRxRead, iFreqRead, iDataRead, nAllocErr, nRxRots, iRx, iLen, itemp
    real(8)         :: rxtheta,rxalpha,rxbeta
    character(180)  :: sLine, sCode, sValue ! These are for reading the lines of the file
    character(180)  :: sFields(6)
    logical         :: bComment, lerror   
    character(180)  :: sDataFormat  ! use this for noting 1.0 or 1.1 version of EMDATA
    
    integer, parameter :: iof = 15   ! I/O File identifier  
    
    
!
! First step is to initialize Dipole1D defaults.  These can then be overridden
! by values in the model file:
! 
    call init_defaults_Dipole1D         

!
! Display a banner for the 1D CSEM/MT Occam implementation:
! 
    write(*,*) ' '      
    write(*,*) '---------------Dipole1D.f90 and CSEM1D.f90-----------------------------'
    write(*,*) ' '
    write(*,*) '         1D Controlled Source EM Kernel Routines                '   
    write(*,*) ' '  
    write(*,*) '                  Copyright 2007-2010                           '
    write(*,*) '                      Kerry Key                                 '
    write(*,*) '           Scripps Institution of Oceanography                  '
    write(*,*) '                     kkey@ucsd.edu                              '
    write(*,*) ''    
    write(*,*) '    CSEM1D.f90:    Version 3.13, Dec, 2010.               '
    write(*,*) '    Dipole1D.f90:  Version 7.3, Feb 10, 2010.               '
    write(*,*) ' '
    write(*,*) '            Developed under sponsorship from the                '
    write(*,*) '         Seafloor Electromagnetic Methods Consortium            '
    write(*,*) '           at Scripps Institution of Oceanography               '
    write(*,*) ' '
    write(*,*) '        See this URL for a list of SEMC funding sponsors:       '
    write(*,*) '            http://marineemlab.ucsd.edu/semc.html               '
    write(*,*) ' '
    write(*,*) ' '     
    write(*,*) '     Dipole1D is free software: you can redistribute it and/or modify'
    write(*,*) '     it under the terms of the GNU General Public License as published by'
    write(*,*) '     the Free Software Foundation, either version 3 of the License, or'
    write(*,*) '     (at your option) any later version.'
    write(*,*) ' '
    write(*,*) '     Dipole1D is distributed in the hope that it will be useful,'
    write(*,*) '     but WITHOUT ANY WARRANTY; without even the implied warranty of'
    write(*,*) '     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
    write(*,*) '     GNU General Public License for more details.'
    write(*,*) ' '
    write(*,*) '     You should have received a copy of the GNU General Public License'
    write(*,*) '     along with Dipole1D.  If not, see <http://www.gnu.org/licenses/>.'
    write(*,*) ' '      
    write(*,*) '-----------------------------------------------------------------------'
!
! Initialize 
!
    nTxIn  = 0
    nRxIn  = 0
    nRxRotSolve = 0

!
! Open the data file:
!
    open (iof, file=trim(dataFileName), status='old', iostat=err)
    if (err /= 0) then
        write(*,*) ' Error opening data file:', trim(dataFileName)
        stop 
    end if  
    write(*,*) ' '
    write(*,*) 'Reading data from file: ', trim(dataFileName)
!
! Loop through the model file looking for the Format and Number of Layers flags
! and skipping over any comments
!
    do while (.true.) 
    
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( iof, '(A180)', iostat = err ) sLine
        
        if (err /= 0) exit  ! end of file read, escape from while(.true.) loop and 
                            ! proceed to checking that required inputs are defined
        call ParseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
!
! What do we have?
!
        select case (trim(sCode))
       
        case ('format')         ! Set cnDataParams for the input data format
        
            call lower(sValue) 
            sDataFormat = trim(sValue)
            if ( (trim(sDataFormat) == 'emdata_1.0') .or. &
               & (trim(sDataFormat) == 'emdata_1.1') .or. & 
               & (trim(sDataFormat) == 'emdata_1.2'))  then
                cnDataParams = 4    
            else
                write(*,*) ' Error: data format unsupported: ', trim(sDataFormat)
                write(*,*) ' Try using EMDATA_1.0 to EMDATA_1.2'
                close(iof)
                stop
            endif
            
        case ('# transmitters')
        
            read(sValue,*) nTxIn
            write(*,fmt='(a24,i6)') '# Transmitters:',nTxIn
            write(*,*) ''    
            write(*,*) '  Transmitter Positions: '
            write(*,54) 'X','Y','Z','Azimuth','Dip'         
            !
            !  Allocate arrays dependent on nTx
            !
            allocate ( azimuthTxIn(nTxIn),dipTxIn(nTxIn), &
                     & xTxIn(nTxIn),yTxIn(nTxIn),zTxIn(nTxIn), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many transmitters (', nTxIn, ')'
                stop 
            endif        
                
            !
            ! Now read in block of transmitters, skipping any comment lines:
            !
            iTxRead = 0
            do while (iTxRead < nTxIn) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iTxRead  = iTxRead + 1                  
                read(sLine,*)   xTxIn(iTxRead),yTxIn(iTxRead),zTxIn(iTxRead),azimuthTxIn(iTxRead),dipTxIn(iTxRead)
                
                ! DGM 10/2010 Don't spew everything to the screen - only need to
                ! check the first so many.
                if( iTxRead < 21 ) then
                    write(*,55)     xTxIn(iTxRead),yTxIn(iTxRead),zTxIn(iTxRead),azimuthTxIn(iTxRead),dipTxIn(iTxRead)
                elseif( iTxRead == 21 ) then
                    write(*,*) 'More than 20 transmitters. Not writing the rest.'
                endif
        
            enddo
            !
            ! Check and make sure all the transmitters were read in:
            !
            if (iTxRead/=nTxIn) then
                write(*,*) ' Error reading transmitters: iTxRead/=nTxIn:',iTxRead, nTxIn
                close(iof)
                stop            
            
            endif       
            
        
                
        case ('# frequencies','nfreq','#freq','#freqs', '# freq')
            
            read(sValue,*) nFreqIn
            write(*,*)  ' '    
            write(*,fmt='(a24,i6)') '# Frequencies:',nFreqIn   
            write(*,fmt='(a24)') 'Frequencies [Hz]:'
            allocate (fTxIn(nFreqIn), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many frequencies (', nFreqIn, ')'
                stop 
            endif   
            
            !
            ! Now read in block of frequencies, skipping any comment lines:
            !
            iFreqRead = 0
            do while (iFreqRead < nFreqIn) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iFreqRead  = iFreqRead + 1                  
                read(sLine,*)  fTxIn(iFreqRead)
                write(*,'(24x,g12.3)')     fTxIn(iFreqRead)
        
            enddo
            !
            ! Check and make sure all the frequencies were read in:
            !
            if (iFreqRead/=nFreqIn) then
                write(*,*) ' Error reading frequencies: iFreqRead/=nFreqIn:',iFreqRead, nFreqIn
                close(iof)
                stop                        
            endif           

        case ('# receivers','nrx','#receivers')
         
            read(sValue,*) nRxIn
            write(*,*) ' ' 
            write(*,fmt = '(a24,i6)') ' Number of Receivers:',nRxIn
            write(*,'(6(a12,2x))') 'x','y', 'z','Theta','Alpha','Beta' 
            !
            ! Allocate the field arrays: 
            !
            allocate(xRxIn(nRxIn), yRxIn(nRxIn), zRxIn(nRxIn),   &
                      ThetaRxIn(nRxIn),AlphaRxIn(nRxIn),BetaRxIn(nRxIn),  &
                      lRxRotSolve(nRxIn), stat=nAllocErr )   
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many receivers (', nRxIn, ')'
                stop 
            endif         
            
            ! initialize rotations to 0 degrees:
            ThetaRxIn = 0
            AlphaRxIn = 0
            BetaRxIn  = 0
            
            ! Intialize rotation solver flags to none:
            nRxRotSolve = 0
            lRxRotSolve = .false.
         
            !
            ! Now read in block of receiver locations, skipping any comment lines:
            !
            iRxRead = 0
 
            do while (iRxRead < nRxIn) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1        
                
                if (trim(sDataFormat) == 'emdata_1.0') then
                   read(sLine,*)  xRxIn(iRxRead),yRxIn(iRxRead),zRxIn(iRxRead)
                   write(*,'(3(f12.2,2x))')  xRxIn(iRxRead),yRxIn(iRxRead),zRxIn(iRxRead)  
                   
                elseif ( (trim(sDataFormat) == 'emdata_1.1').or.(trim(sDataFormat) == 'emdata_1.2')) then
                    call parseFields( len(sLine), sLine, 6, sFields)
                    read(sFields(1),*) xRxIn(iRxRead)
                    read(sFields(2),*) yRxIn(iRxRead)
                    read(sFields(3),*) zRxIn(iRxRead)  
                    
                    if ( (trim(sFields(4)) == '?').or.(trim(sFields(5)) == '?').or.(trim(sFields(6)) == '?') ) then 
                        nRxRotSolve = nRxRotSolve + 1
                        lRxRotSolve(iRxRead) = .true.
                    else 
                      read(sFields(4),*) ThetaRxIn(iRxRead)
                      read(sFields(5),*) AlphaRxIn(iRxRead)
                      read(sFields(6),*) BetaRxIn(iRxRead)
                    endif
                    write(*,'(3(f12.2,2x),3(a12,2x))')  xRxIn(iRxRead),yRxIn(iRxRead),zRxIn(iRxRead), & 
                          & trim(sFields(4)),trim(sFields(5)),trim(sFields(6))  
                endif
                
            enddo
            !
            ! Check and make sure all the layers were read in:
            !
            if (iRxRead/=nRxIn) then
                write(*,*) ' Error reading receivers: iRxRead/=nRxIn:',iRxRead,nRxIn
                close(iof)
                stop                        
            endif           
    
    !
    ! Optional parameter for receiver rotations.  This allows for modeling data
    ! at the orientation of the receiver, which can be different from x,y,z.
    !  So this is actually the angle the model responses are rotated by. Degrees clockwise from x toward y.
    !
        case ('# receiver orientations')     ! EMDATA_1.0 only...
            read(sValue,*) nRxRots
            write(*,*) ' ' 
            write(*,fmt = '(a42,i6)') ' Number of Modified Receiver Orientations:',nRxRots
            write(*,'(a6,2x,3(a7,2x),2x,a9)') ' Rx#', 'Theta', 'Alpha','Beta','[degrees]'
            iRxRead = 0
            do while (iRxRead < nRxRots) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iRxRead  = iRxRead + 1                  
                read(sLine,*)  irx, rxtheta, rxalpha, rxbeta
                ThetaRxIn(irx) = rxtheta
                AlphaRxIn(irx) = rxalpha
                BetaRxIn(irx)  = rxbeta
                write(*,'(i6,2x,3(f7.2,2x))') irx,ThetaRxIn(irx), AlphaRxIn(irx),BetaRxIn(irx)
                
        
            enddo           
    !
    ! Should we solve for Rx orientation?
    !
        case ('# solve orientations','# unknown orientations','# unknown rotations')  
          
            read(sValue,*) nRxRotSolve  ! EMDATA_1.0 only...
            write(*,*) ' ' 
            write(*,fmt = '(a24,i6)') ' # Unknown Rotations:',nRxRotSolve
!            write(*,*) ' The Rx orientation solutions will be output to the  '
!            write(*,*) ' RESP##.assoc files'
!           allocate( IndRxRotSolve(nRxRotSolve), stat=nAllocErr ) 
            write(*,'(a24)') ' Receiver #:'
            iRxRead = 0
            do while (iRxRead < nRxRotSolve) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Rx index
                iRxRead  = iRxRead + 1                  
                read(sLine,*)  itemp
                lRxRotSolve(itemp) = .true.
                write(*,'(24x,i6)') itemp
        
            enddo       
    
    !
    ! DGM Feb 2010 - add support for finite length dipole
    !
        case ('dipole length','dipole','dipolelength')
            read(sValue,*) lenTx1D 
            write(*,fmt = '(a24,f6.1)') 'Dipole Length:', lenTx1D
    
        case ('# integ pts','# integration points','numIntegPts','number of integration points')
            read(sValue,*) numIntegPts
            write(*,fmt = '(a24,I6)') '# integration points:', numIntegPts
            
    !
    ! Phase Convention:
    !
        case ('phase','phase convention')  
            call lower(sValue)
            
            select case (trim(sValue))
            
            case ('lag')
                    phaseConvention = 'lag'
            case ('lead')
                    phaseConvention = 'lead'
            end select
            
            write(*,*) ' '
            write(*,fmt = '(a24,a4)') 'Phase convention: ', phaseConvention
            
    !
    ! Lastly we need to read in the DATA block.
    !
         case ('# data','#data','ndata')
         
            read(sValue,*) nd
            write(*,*) ' ' 
            write(*,fmt = '(a24,i6)') ' Number of Data:',nd
            
            !
            ! Allocate the nd data arrays: 
            !
            allocate(dp(nd,cnDataParams), d(nd), d_copy(nd), sd(nd), sd_copy(nd), stat=nAllocErr ) !kwk d_copy,sd_copy need dealloc
            
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many data points (', nd, ')'
                stop 
            endif        
         
            iDataRead = 0
            write(*,'(4(a8,1x),a12,2x,a12,2x)') 'Type','Freq #','Tx #','Rx #', 'Data','Std Error'    
            
            do while (iDataRead < nd) 
    
                ! Get the next code/value pair
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
               
               ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Else it is not a comment, so read in the Tx parameters from sLine
                iDataRead  = iDataRead + 1     
                
                ! Read in 6 fields as character strings,
                ! this is clunky but works for mixed character/numeric tables:
                call parseFields( len(sLine), sLine, 6, sFields)    

                ! Test first field to see if it is character data type, 
                ! if so convert it to numeric code given in module csem1d_mod
                iLen = len_trim(sFields(1))
                
                if (iLen < 5) then  ! Field is integer code for data type
                    read(sFields(1),*) dp(iDataRead,1)  
                    
                else ! Field is ascii character string describing data:
                    call getDataCode(sFields(1),dp(iDataRead,1))
                
                endif

                ! Parse the remaining fields for the data parameters and data 
                read(sFields(2),*) dp(iDataRead,2)
                read(sFields(3),*) dp(iDataRead,3)
                read(sFields(4),*) dp(iDataRead,4)
                read(sFields(5),*) d(iDataRead)
                read(sFields(6),*) sd(iDataRead)

                
                ! DGM 10/2010 Don't spew everything to the screen - only need to
                ! check the first so many.
                if( iDataRead < 101 ) then
                    write(*,fmt= '(4(i8,1x),ES12.3,2x,ES12.3,2x)') (dp(iDataRead,i),i=1,4),d(iDataRead),sd(iDataRead)
                elseif( iDataRead == 101 ) then
                    write(*,*) 'More than 101 data. Not writing the rest.'
                endif
 
            enddo        
            !
            ! Check and make sure all the data was read in:
            !
            if (iDataRead/=nd) then
                write(*,*) ' Error reading data: iDataRead/=nd:',iDataRead,nd
                close(iof)
                stop                        
            endif      
            !
            ! Copy data array for OPRA:
            !
            d_copy = d
            sd_copy = sd
            
            
            case default
                write(*,*) 'Error, unknown code in Data file: ',trim(sCode)
                write(*,*) 'Stopping!'
                stop
            
        end select ! case (trim(sCode))
        
        
    enddo ! while (.true.) ! The main data reading loop
!
! Close the data file
!
    close(iof)


!
! Check that parameters have been defined:
!
    lerror = .false.
    if (nTxIn < 1) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of transmitters: nTx = ',nTxIn
    endif
    if (nRxIn < 1) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of receivers: nRx = ',nRxIn
    endif
    if (nFreqIn < 1) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of frequencies: nFreq = ',nFreqIn
    endif
    if (nd < 1) then
        lerror = .true.
        write(*,*) ' Error in Data File, insufficient number of data: # Data = ',nd
    endif   
    
    if (lerror) then
        write(*,*)  ' Stopping ! '
        stop
    endif

!
! Allocate arrays for associated parameters if OPRA rotation solver has been requested:
!
    if (nRxRotSolve > 0) then   
        ! These are stored in pm_assoc  (associated parameters)
        npm_assoc = nRxRotSolve*3 ! solving for 3 angles to get full 3D orientation (tilts set to zero if only horiz data).
        allocate ( pm_assoc (npm_assoc), stat=nAllocErr ) 
    endif
    
 !
 ! Finally, apply reciprocity, if possible:
 !
    call apply1Dreciprocity
!
! Print format statements:
!
54  format(a9,1x,a9,1x,a9,1x,a9,1x,a9)   
55  format(f9.1,1x,f9.1,1x,f9.1,1x,f9.1,1x,f9.1)  

      
    end subroutine readData       
    
!==============================================================================!  
!=========================================================! apply1Dreciprocity !  
!==============================================================================!    
    subroutine apply1Dreciprocity()
!
! Subroutine apply 1D reciprocity to input Rx and Tx arrays.  This will speed 
! up the computations by turning multiple transmitters into multiple receivers,
! if permissible (same zTx, AzimuthTz and DipTx)
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 1.0  Oct 2009.  DGM moved this out of readData for compatibility with 
!                         his external interface.
!
!------------------------------------------------------------------------------
! Create Rx and Tx arrays for fast computations.  This means reducing the 
! number of discrete transmitters since Dipole1D's runtime is proportional 
! to the number of transmitters.  So, we will look through the array of 
! transmitters. If all transmitters have the same zTx,AzimuthTx and DipTx 
! then we can simply reduce them to a single transmitter and 
! many sites (in x and y).
!------------------------------------------------------------------------------

    use Occam
    use csem1d_mod
    use dipole1d
    
    implicit none
    
    integer(4)      :: nAllocErr, iRx, icnt, iTx, i
    logical         :: lTxSame
    
    !
    ! Allocate a dpTfm array that contains transformed data parameter indices:
    !
    
    allocate ( dpTfm(nd,2), stat = nAllocErr ) ! dpTfm(nd,[iTxTfm iRxTfm]) 
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many receivers  (', nRxIn, ')'
        stop 
    endif   
        
    !
    ! Loop through nTx  and check for identical transmitters parameters:
    ! 
    
     lTxSame = .true.
     do iTx = 2,nTxIn
        if ( abs( zTxIn(iTx)       - zTxIn(1)       ) > 1d-2 ) lTxSame = .false.  
        if ( abs( azimuthTxIn(iTx) - azimuthTxIn(1) ) > 1d-2 ) lTxSame = .false.
        if ( abs( dipTxIn(iTx)     - dipTxIn(1)     ) > 1d-2 ) lTxSame = .false.
     enddo

    !
    ! If not similar, then just copy default Rx and Tx params:
    !
    
    if (.not.lTxSame) then  
    
        allocate ( xTxTfm(nTxIn), yTxTfm(nTxIn), zTxTfm(nTxIn), azimuthTxTfm(nTxIn), dipTxTfm(nTxIn), &
                 & stat = nAllocErr)
                 
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Too many transmitters  (', nTxIn, ')'
            stop 
        endif   
        
        nTxTfm       = nTxIn
        xTxTfm       = xTxIn
        yTxTfm       = yTxIn
        zTxTfm       = zTxIn
        azimuthTxTfm = azimuthTxIn
        dipTxTfm     = dipTxIn
        
        allocate ( xRxTfm(nRxIn), yRxTfm(nRxIn), zRxTfm(nRxIn),iRxTfm(nRxIn), stat = nAllocErr)
        
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Too many receivers  (', nRxIn, ')'
            stop 
        endif   
        
        nRxTfm = nRxIn
        xRxTfm = xRxIn
        yRxTfm = yRxIn
        zRxTfm = zRxIn
        do i= 1,nRxIn
            iRxTfm(i) = i
        enddo 
        
        dpTfm(1:nd,1) = dp(1:nd,3)
        dpTfm(1:nd,2) = dp(1:nd,4)  
        
    !
    ! If Tx's are similar, then create transformation arrays:
    !   
    
    else ! lTxSame == .true.
    
       ! Define a single transmitter:
        nTxTfm = 1
        allocate ( xTxTfm(1), yTxTfm(1), zTxTfm(1), azimuthTxTfm(1), dipTxTfm(1), &
                 & stat = nAllocErr)
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Too many transmitters  (', 1, ')'
            stop 
        endif           
       
        xTxTfm       = 0d0
        yTxTfm       = 0d0
        zTxTfm       = zTxIn(1)
        azimuthTxTfm = azimuthTxIn(1)
        dipTxTfm     = dipTxIn(1)
        
      ! Define array of sites:
        nRxTfm = nTxIn*nRxIn
        allocate ( xRxTfm(nRxTfm), yRxTfm(nRxTfm), zRxTfm(nRxTfm),iRxTfm(nRxTfm), stat = nAllocErr)
        
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Too many receivers  (', nRxIn, ')'
            stop 
        endif         
    
        icnt = 0    
        do iRx = 1,nRxIn
           do iTx = 1,nTxIn
               icnt = icnt + 1
                xRxTfm(icnt) = xRxIn(iRx) - xTxIn(iTx)
                yRxTfm(icnt) = yRxIn(iRx) - yTxIn(iTx)
                zRxTfm(icnt) = zRxIn(iRx)
                iRxTfm(icnt) = iRx ! store the original site number for the Tfm site
           enddo
        enddo
        
        ! Now update the dpTfmarray to point to transformed transmitter and sites:
        dpTfm(1:nd,1) = 1  ! only 1 transformed transmitter 
        
        do i = 1,nd  ! 
           iTx = dp(i,3)        
           iRx = dp(i,4)
           dpTfm(i,2) = iTx + (iRx-1)*nTxIn  ! in other words transformed receiver index is now in blocks of iTx
        enddo
                
        write(*,*)  ' '
        write(*,*) ' Transmitters have same depth, azimuth and dip, so I transformed '
        write(*,*) ' the transmitters into many receivers for a single transmitter '
        write(*,*) ' in order to speed up the computations. '
        write(*,*) ' The number of effective receivers is: ',icnt 
        
    endif ! (.not.lTxSame) 
    
    
!   
! Allocate the arrays passed in dipole1d module.  Note that derivative arrays are allocated in readModel.
!
    allocate( x1D(nRxTfm),y1D(nRxTfm),z1D(nRxTfm),    &   
            & ex1D(nRxTfm),ey1D(nRxTfm),jz1D(nRxTfm), &
            & bx1D(nRxTfm),by1D(nRxTfm),bz1D(nRxTfm), &
            & sigsite(nRxTfm),RxInd(nRxTfm), rx1DtoRxIn(nRxTfm), &
            & stat=nAllocErr )  
            
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many receivers  (', nRxTfm, ')'
        stop 
    endif   
    
    return 
    
    end subroutine apply1Dreciprocity   
!==============================================================================!  
!===============================================================! getDataCodes !  
!==============================================================================!  
    subroutine getDataCode( sField, iCode)
!
! Management routine to get data format magic number code from ascii text label.
!
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
!
! Revision      Date            Author          Notes
! 1.0           November 2008   Kerry Key       Created
!
    use occam, only : lower 
    use csem1d_mod  ! contains the magic numbers
    
    implicit none
    
    character(180), intent(inout) :: sField
    integer, intent(out) :: iCode
    
    !
    ! Get lower case string:
    !
    call lower(sField) 
    
    !
    ! Decode string:
    !
    select case (trim(sField))
    
    !
    ! CSEM Data:
    !
        ! Electric Real/Imag Data:
        case ('realex','reex')
            iCode = indRealEx
        case ('imagex','imex','imaginaryex')    
            iCode = indImagEx
        case ('realey','reey')
            iCode = indRealEy
        case ('imagey','imey','imaginaryey')    
            iCode = indImagEy   
        case ('realez','reez')
            iCode = indRealEz   
        case ('imagez','imez','imaginaryez')    
            iCode = indImagEz
        
        ! Magnetic Real/Imag Data:
        case ('realbx','rebx')
            iCode = indRealBx
        case ('imagbx','imbx','imaginarybx')    
            iCode = indImagBx
        case ('realby','reby')
            iCode = indRealBy
        case ('imagby','imby','imaginaryby')    
            iCode = indImagBy       
        case ('realbz','rebz')
            iCode = indRealBz
        case ('imagbz','imbz','imaginarybz')    
            iCode = indImagBz   
        
        ! Electric Amplitude and Phase Data:
        case ('ampex','amplitudeex')
            iCode = indAmpEx
        case ('phsex','phaseex')
            iCode = indPhsEx
        case ('ampey','amplitudeey')
            iCode = indAmpEy
        case ('phsey','phaseey')
            iCode = indPhsEy
        case ('ampez','amplitudeez')
            iCode = indAmpEz                
        case ('phsez','phaseez')
            iCode = indPhsEz        

        ! Magnetic Amplitude and Phase Data:
        case ('ampbx','amplitudebx')
            iCode = indAmpBx        
        case ('phsbx','phasebx')
            iCode = indPhsBx
        case ('ampby','amplitudeby')
            iCode = indAmpBy
        case ('phsby','phaseby')
            iCode = indPhsBy
        case ('ampbz','amplitudebz')
            iCode = indAmpBz        
        case ('phsbz','phasebz')
            iCode = indPhsBz    
            
        ! CSEM polarization ellipse parameters:
        case('pemax')
            icode = iPEmax
        case('pemin')
            icode = iPEmin
        case('pbmax')
            icode = iPBmax
        case('pbmin')
            icode = iPBmin
    !   
    ! Magnetotelluric Data:
    !
        ! Real/Imag
        case ('realzxx')        ! ignored in 1D
            iCode = indRealZXX  
        case ('imagzxx')        ! ignored in 1D
            iCode = indImagZXX      
        case ('realzxy')
            iCode = indRealZXY  
        case ('imagzxy')
            iCode = indImagZXY      
        case ('realzyx')
            iCode = indRealZYX  
        case ('imagzyx')
            iCode = indImagZYX          
        case ('realzyy')        ! ignored in 1D
            iCode = indRealZYY  
        case ('imagzyy')        ! ignored in 1D
            iCode = indImagZYY                  
        
        ! Apparent resistivity and phase:
        case ('rhozxx')         ! ignored in 1D
            iCode = indRhoZXX   
        case ('phszxx')         ! ignored in 1D
            iCode = indPhsZXX       
        case ('rhozxy')         
            iCode = indRhoZXY   
        case ('phszxy')         
            iCode = indPhsZXY       
        case ('rhozyx')         
            iCode = indRhoZYX   
        case ('phszyx')         
            iCode = indPhsZYX       
        case ('rhozyy')         ! ignored in 1D
            iCode = indRhoZYY   
        case ('phszyy')         ! ignored in 1D
            iCode = indPhsZYY       
                
        
        ! Error if unknown type:
        case default
            write(*,*) ' Error decoding data type: ', sField
            write(*,*) ' Stopping!'
            stop
        
    
    end select

    
    end subroutine getDataCode       
!==============================================================================!  
!==================================================================! readModel !  
!==============================================================================!  
    subroutine readModel()
!
! Subroutine to read in an Occam1DCSEM model file and allocate model dependent
! arrays.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!   
    use Occam    ! prewts, premod are initialized to 0, then read in from model file (if used)
                     
    use csem1d_mod   ! stores all the variables needed to go between Occam and Dipole1D 
    
    use mt1D_mod     ! used for passing 1D MT parameters and responses
    
    use dipole1d  ! Model nLay1D and zlay1D set, allocates sig1D(nlay1D), zlay1D(nlay1D)
                     ! note that sig1D is not set since it needs to be updated to 
                     ! the current model on each call to Dipole1D (done in subr. Fwd_Calc)
                     
!
! Layer block in file is
! 
! [top_depth    resistivity  penalty    preference   pref_penalty]
! 
! Use 0 for Preference to not have a Preferenced parameter
! Penalty is penalty to apply at the layer top boundary.  Use 0 for no roughness penalty at layer top
!
    
    implicit none
    
!
! Local variables:
! 
    integer         :: i, j, err, iof, ilayread, ichecknParams, nAllocErr
    character(180)  :: sLine, sCode, sValue ! These are for reading the lines of the file
    logical         :: bComment
    character(180)   :: sFields(5)    
    
    parameter (iof = 15)     ! I/O File identifier
!
! Initialize number of layers in model:
!
    nlay1D = 0
!
! Open the Model File:
!     
    open (iof, file=trim(modelFileName), status='old', iostat=err)
    if (err /= 0) then
        write(*,*) ' Error opening model file:', trim(modelFileName)
        stop 
    end if   
    write(*,*) ' '
    write(*,*) 'Reading model from file: ', trim(modelFileName)
!
! Loop through the model file looking for the Format and Number of Layers flags
! and skipping over any comments
!
    do while (.true.) 
    
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( iof, '(A180)', iostat = err ) sLine
        
        if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
        call ParseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        
        ! What do we have?
        select case (trim(sCode))
       
        case ('format')
            call lower(sValue)
            if ( (trim(sValue) /= 'occam1dmod_1.0') .and.   &  
            &    (trim(sValue) /= 'resistivity1dmod_1.0') ) then
                write(*,*) ' Erorr: model format unsupported: ', trim(sValue)
                write(*,*) ' Try using Resistivity1DMod_1.0'
                close(iof)
                stop
            endif
            
        case ('# layers','#layers', 'number of layers','nlayers')   
            
            ! Get the number of layers value:
            read(sValue,*)  nlay1D
                                    
            ! [top_depth,rho,penalty,preference,prej_penalty]    
            allocate ( layer_params(nlay1D,5) , zlay1D(nlay1D), sig1D(nlay1D), stat=nAllocErr )   

            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many layers in readModel (', nlay1D, ')'
                stop 
            endif        
            
        
            
            ! read in the layering blocks (skipping any comment lines)
            ilayread = 0
            do while (ilayread < nlay1D) 
    
                ! read in  a line:
       
                read( iof, '(A180)', iostat = err ) sLine
        
                if (err /= 0) exit  ! end of file read, escape from while loop and 
                          
                ! Check if line is a comment, if not strip off comments at end of line:
                call ParseLine( len(sLine), sLine, bComment )
                if( bComment ) cycle
                ! Not a comment, so read in the values from sLine:  
                ilayread  = ilayread + 1    
                !read(sLine,*) layer_params(ilayread,1:5)
                
                ! Read in 5 fields as character strings,
                ! this is clunky but works for mixed character/numeric tables:
                call parseFields( len(sLine), sLine, 5, sFields)
                
                if (trim(sFields(2)) == '?') then
                     layer_params(ilayread,2) = -1
                else 
                    read(sFields(2),*) layer_params(ilayread,2)
                endif
                
                read(sFields(1),*) layer_params(ilayread,1)
                read(sFields(3),*) layer_params(ilayread,3)
                read(sFields(4),*) layer_params(ilayread,4)
                read(sFields(5),*) layer_params(ilayread,5)
          
          
            enddo
            !
            ! Check and make sure all the layers were read in:
            !
            if (ilayread/=nlay1D) then
                write(*,*) ' Error: model ilayread/=nlay1D:',ilayread,nlay1D
                close(iof)
                stop            
            
            endif
            
            ! Now print them to the screen as a debug check:
            write(*,*) ''
            write(*,'(5(a12,2x))') 'Top Depth','Linear Rho', 'Penalty', 'Preference', 'Pref Penalty'
            do i = 1,nlay1D
                    write(*,'(F12.3,2x,G12.2,2x,3(f12.1,2x))') (layer_params(i,j),j=1,5)
            enddo
    !
    ! The rest of the cases are for overriding the Dipole1D default setup
    !       
        case ('ht filters')
        
            HTmethod1D      = sValue(1:len(HTmethod1D))
            write(*,*) ' '
            write(*,fmt = '(a,a24)') ' Using Hankel TF Filters: ', HTmethod1D

            
        case('usespline1d')
        
            select case (trim(sValue))
           case ('yes')
               lUseSpline1D = .true.
               write(*,*) 'Using Spline Interpolation for speed'
               write(*,*) ' '
           case default
               lUseSpline1D = .false.
           end select
            
        case default
            write(*,*) 'Error reading Model file!'
            write(*,*) 'Unknown or unsupported code:', sCode
            stop            

        end select ! case (trim(sCode))
        
    enddo
!
! Close the model file:
!
    close(iof)      
!
! Check to make sure layers are defined: 
!
    if (nlay1D <= 0) then   
        write(*,*) ' Error reading model file: number of layers undefined! '    
        stop    
    endif
! 
!  Count the number of free parameters to see if it agrees with the startup iteration file
!
    ichecknParams = 0
    do i = 1,nlay1D
        if (layer_params(i,2) < 0) then ! negative resistivity is free parameter
            ichecknParams = ichecknParams + 1
        endif
    enddo
    if (ichecknParams /= nParams) then
        write(*,*) ' Error number of free parameters in model file disagrees with startup/iteration file! ' 
        write(*,*) ' Fix them so they agree. Stopping!'
        write(*,*) ' Number of free paramters in iteration file: ', nParams
        write(*,*) ' Number of free layers in model file: ',        ichecknParams
        write(*,*) '  Stopping!'        
        stop        
    endif

!
! Now set up a few model related variable:
!
    call setup1Dlayers

     
    end subroutine readModel       
    
!==============================================================================!  
!==============================================================! setup1Dlayers !
!==============================================================================!
! Broken out of readModel for support of MatLab interface.  These don't need 
! the input file.
!
    subroutine setup1Dlayers
    
    use occam
    use csem1d_mod
    use dipole1d
    
    implicit none
    
    integer(4)   :: nAllocErr, ichecknParams, i
    logical      :: lhasPref
    
!
! Assign layer top depths:
!
    zlay1D(1:nlay1D) = layer_params(1:nlay1D,1)
    
!
! Notes: upon each call to dipole1d, sig1d will be filled using layer_params(i,2) 
! or pm where layer_params(i,2) < 0 (ie free parameter)
!
   
!
! Set Preference and weights:
!
    ichecknParams = 0
    premod = 0d0 ! initialize these to 0 (unused)
    prewts = 0d0
    lhasPref = .false.

    do i=1,nlay1D
        if (layer_params(i,2) < 0) then ! negative resistivity is free parameter
            ichecknParams = ichecknParams + 1
            if (abs(layer_params(i,4)) /= 0) then
                premod(ichecknParams) = log10(layer_params(i,4)) ! convert to log10 scale resistivity (only iteration file has 
                                                                 ! log10 resistivity)
            else
                premod(ichecknParams) = 0d0
            endif
            
            prewts(ichecknParams) = layer_params(i,5) 
            if ( prewts(ichecknParams)  /= 0 ) then 
                if (.not.lhasPref) then
                    write(*,*) ' '
                    write(*,*) ' Model Preference (linear ohm-m) and Preference Penalty Weights:'
                    write(*,'(A12,2x,A12,2x,A12)') ' Parameter #',' Preference ', 'Weight'
                endif
                lhasPref = .true.
                write(*,'(I12,2x,G12.2,2x,G12.2)') ichecknParams,10**premod(ichecknParams),prewts(ichecknParams)
            endif
        endif       
    enddo
    if (.not.lhasPref) then
        write(*,*) ' '
        write(*,*) '      No model preference values used.'
    endif
    write(*,*) ' '
    
 
!
! Allocate the jacobian arrays output from dipole1d:
!
    allocate( dexdsig(nRxTfm,nlay1D),deydsig(nRxTfm,nlay1D),djzdsig(nRxTfm,nlay1D), &
            & dbxdsig(nRxTfm,nlay1D),dbydsig(nRxTfm,nlay1D),dbzdsig(nRxTfm,nlay1D), &
            &  stat=nAllocErr )             
            
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many receivers x 1D layers  (', nRxIn,' x ',nlay1D, ')'
        stop 
    endif   

     
    end subroutine setup1Dlayers   
    
!==============================================================================!  
!=================================================================! computeFwd !  
!==============================================================================!  
    subroutine computeFwd( bDoPartials, currentMod )
!
! Routine to compute the 1DCSEM forward response and model Jacobian matrix
! This is called from Occam.f90.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 2.0   August, 2008.  Moved main sections to individual subroutines, cleaned up
! Version 1.0   April, 2008.
!
    use Occam        ! Passes out dm and wj here
    use csem1d_mod
    use dipole1d
    
    implicit none

!
! Arguments
!
    logical, intent(in)        :: bDoPartials           ! logical flag for computing Jacobians
    real(RealPrec), intent(in) :: currentMod(nParams)   ! input model parameters

!
! Local variables:
!
    integer :: i, iparam, nAllocErr
    integer, dimension(:), allocatable :: ipm2ilay ! array index mapping parameters to 1D layers
      
    allocate (ipm2ilay(nParams), stat=nAllocErr)
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many parameters (', nParams, ')'
        stop 
    endif           
! 
! Set up a few things:
! 
!
! 1. Copy the model parameters (log10(resistivity)) in currentMod into sig1d.
! 2. Create ipm2ilay index that points from parameters to 1D layers.
! 3. Check to make sure the conductivity is within reasonable bounds.
!
    iparam = 0
    do i=1,nlay1D
        if (layer_params(i,2) < 0) then                 ! negative resistivity is free parameter    
            iparam = iparam + 1
            ipm2ilay(iparam) = i                        ! store mapping of parameter number to 1D model layer 
            sig1D(i) = 1.d0/ (10**currentMod(iparam))   ! convert to conductivity from the log10(resistivity) parameters
        !
        ! Make sure conductivity isn't off scale: 
        !
            if (log10(sig1D(i)) > 12) then
                sig1d(i) = 1d12
              ! write(*,*) ' warning, resistivity to low, moving back up to 1d-12 for layer: ',i
            endif
            if (log10(sig1D(i)) < -12) then
                sig1d(i) = 1d-12
               !write(*,*) ' warning, resistivity to high, moving back down to 1d12 for layer: ',i
            endif           
        else   ! this is a fixed layer, read it in from the layer_params array
            sig1D(i) = 1.d0/layer_params(i,2)
        endif 
    enddo 

! 
! Compute responses and sensitivities for 1D CSEM data:
! 
    call computeFwd_CSEM( bDoPartials, ipm2ilay )
    
! 
! Solve for CSEM Site orientations using orthogonal Procrustes rotation analysis:
! 
    call computeCSEM_OPRA()
    
! 
! Compute responses and sensitivities for 1D MT data:
!  
    call computeFwd_MT( bDoPartials, ipm2ilay  )
 
! 
! We're done, deallocate and go home, hooray! :)
! 
    deallocate ( ipm2ilay )
    
    ! DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
    ! DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
    ! This bit of code is useful for keeping track of where your
    ! sensitivities are going. Just slurp it into MatLab and use
    ! pcolor or surf to plot it.
!    if(bDoPartials) then
!        open( unit=42, file='LatestWJ.dat' )
!        write(42,*) wj
!        close(42)
!    endif
    ! DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
    ! DEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUGDEBUG
    
    return
    
    end subroutine computeFwd
    
!==============================================================================!  
!============================================================! computeFwd_CSEM !  
!==============================================================================!  
    subroutine computeFwd_CSEM( bDoPartials, ipm2ilay )
!
! Routine to compute the 1DCSEM forward response and model Jacobian matrix
!
! This is structured in a similar manner to David Myer's 
! Fortran90 computeFwd for the 2DMT problem. There's some complicated 
! stuff that deals with finding only the minimum number of Rx's and Tx's so
! that the 1D CSEM comps are done somewhat efficiently.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! 
! Version 2.1   November, 2008.   Added in support for phase data.
! Version 2.0   September, 2008.  Added support for fixed 3D rotation of responses
! Version 1.0   April, 2008.
!
!------------------------------------------------------------------------------

    use Occam
    use csem1d_mod
    use dipole1d
    use dataTransformations
    
    implicit none

!
! Arguments
!
    logical, intent(in) :: bDoPartials         ! logical flag for computing Jacobians  
    integer,intent(in)  :: ipm2ilay(nParams) 
!
! Local variables:
!
    integer    :: i,j, iTx, iFreq,  iRx, nAllocErr    

    
    real(8)                     :: theta, alpha, beta           ! site rotation angles
    real(8)                     :: cct, sst, cca, ssa, ccb, ssb ! cos and sin values for rotation angles
    real(8) , dimension(3,3)    :: RotR,Rotx,Roty,Rotz          ! 3x3 rotation matrices
    complex(8), dimension(3)    :: tempvec                      ! temporary vector of ex,ey,ez for rotations
    
    integer, dimension(:), allocatable :: iRxlayer
!
! Initialize:
!  
    if (bDoPartials) then 
        linversion = .true.
    else
        linversion = .false.
    endif
    
            
!
! Loop over transmitters and frequencies and compute the CSEM responses:
!

!
! Loop over each transmitter:
! 
    do iTx = 1,nTxTfm
!
! Assign Tx parameters used by Dipole1D:
!
        xTx1D       = xTxTfm(iTx)
        yTx1D       = yTxTfm(iTx)
        zTx1D       = zTxTfm(iTx)
        azimuthTx1D = azimuthTxTfm(iTx)
        dipTx1D     = dipTxTfm(iTx)

!
! Inner loop over requested frequencies
!
        do iFreq= 1,nFreqIn
!
! Get the requested frequency:
!
            ftx1D = fTxIn(iFreq) 
            
!
! Loop through the data array and get only the sites used for this iTx and iFreq:
!

! 1. Set x1d = -1 as default (site isn't part of current batch)
            x1d = -1

! 2. Loop over data and find sites that belong in the current batch, also see if B's in the data in addition to E's
            lbcomp = .false.
            do i = 1,nd
              ! if iTx and iFreq, x1d = 1
                if (  ( dp(i,2) == iFreq) .and.  ( dpTfm(i,1) == iTx) .and. (dpTfm(i,2) /= 0)) then
                     
                    x1d( dpTfm(i,2) ) = 1
                    if ( ( (dp(i,1) >  10).and.(dp(i,1) <  20) ) .or. &
                         ( (dp(i,1) >  30).and.(dp(i,1) <  40) ) .or. &
                         ( (dp(i,1) == 43) ) )  then
                       lbcomp =.true.  
                    endif
                endif
            enddo
  
! 3. n1d = 0
            n1d = 0
            RxInd = 0
! 4. Loop over i=1,nRxIn
            do i = 1,nRxTfm
!   a. Add site to current batch if x1d > 0
                if (x1d(i) > 0) then
                    n1d = n1d + 1
                    x1d(n1d)  = xRxTfm(i)
                    y1d(n1d)  = yRxTfm(i)                    
                    z1d(n1d)  = zRxTfm(i)    
                    RxInd(i)  = n1d               ! index from Tfm array site number to current batch site number, 0 if not in batch
                    rx1DtoRxIn(n1d) = iRxTfm(i)   ! index of original input Rx number to use for rotations. iRxTfm points to RxIn 
                                                  ! site number
                endif
            enddo

            
               
! If no sites for iTx and iFreq, cycle the loop over iFreq:
!
            if (n1d == 0) cycle
!       
! Loop over the sites and assign the conductivity at the site (for converting Jz to Ez)
!
            allocate(iRxlayer(n1d), stat = nAllocErr)  ! stores layer each site resides in
            if (nAllocErr .ne. 0) then
                write(*,*) 'Out of memory.  Too many layers (', n1d, ')'
                stop 
            endif        
            
            do i=1,n1d
                sigsite(i)  = sig1D(1) 
                iRxlayer = 1
                do j = 2,nlay1D
                    if (zlay1D(j).lt.z1d(i)) then ! zlay1D contains layer top depths
                       sigsite(i) = sig1D(j)      ! Ends on last layer top being shallower than site location, 
                                                  ! so a site on a boundary uses the upper layer,
                                                  ! (i.e., a seafloor site uses the sea conductivity).
                        iRxlayer(i) = j             
                    endif
                enddo       
            enddo   
                

!
! Compute CSEM responses using Dipole1D:
!           
            call comp_dipole1D
 
!
! Convert Jz and Jz derivatives to Ez and Ez derivatives.  Here the assumption is that a receiver located on a boundary
! measures Ez in the top layer, since present day EM receivers on the seafloor measure Ez in the 
! sea, not the seabed.
!
            do i=1,n1d
                jz1d(i) = jz1d(i)/sigsite(i)  ! converts jz to ez
                  
                ! Convert sensitivity as well:  
                if (bDoPartials) then
                    forall( j=1:nlay1D ) 
                         djzdsig(i,j)  = djzdsig(i,j)/sigsite(i)
                    end forall
                    ! Extra term when i=j:
                    djzdsig(i,iRxlayer(i))  = djzdsig(i,iRxlayer(i))  - jz1d(i)/sigsite(i)  ! note that here jz is ez
                endif
            enddo
     


 
!
! Rotate the model response if the data file gives a modified orientation for the receiver(s):
!
            do i=1,n1d
            
                ! Get original site number:
                iRx = rx1DtoRxIn(i)  !rx1dtoRxIn(i) points to RxIn site number
                
                ! Get site rotation parameters:
                theta = ThetaRxIn(iRx)
                alpha = AlphaRxIn(iRx)
                beta  = BetaRxIn(iRx)
                
                if ( (theta /= 0) .or. (alpha /= 0) .or. ( beta /= 0) ) then
                
                    cct = cos(deg2rad* theta)
                    sst = sin(deg2rad* theta)
                    cca = cos(deg2rad* alpha)
                    ssa = sin(deg2rad* alpha)
                    ccb = cos(deg2rad* beta)
                    ssb = sin(deg2rad* beta)                
                    
                    Rotx = 0d0
                    Rotx(1,1) = 1.d0
                    Rotx(2,2) =  ccb
                    Rotx(2,3) = -ssb
                    Rotx(3,2) =  ssb
                    Rotx(3,3) =  ccb
                    
                    Roty = 0d0
                    Roty(1,1) =  cca  
                    Roty(1,3) = -ssa
                    Roty(2,2) =  1.d0 
                    Roty(3,1) =  ssa
                    Roty(3,3) =  cca
                    
                    
                    Rotz = 0d0
                    Rotz(1,1) =  cct
                    Rotz(1,2) = -sst
                    Rotz(2,1) =  sst
                    Rotz(2,2) =  cct
                    Rotz(3,3) =  1.d0    
                    
                    RotR = matmul(Roty,Rotx)
                    RotR = matmul(Rotz,RotR)
                
                    ! Rotate Responses:
                    tempvec(1) = ex1d(i)
                    tempvec(2) = ey1d(i)
                    tempvec(3) = jz1d(i)
                    tempvec = matmul(tempvec,RotR)
                    ex1d(i) = tempvec(1) 
                    ey1d(i) = tempvec(2) 
                    jz1d(i) = tempvec(3) 
                    
                    tempvec(1) = bx1d(i)
                    tempvec(2) = by1d(i)
                    tempvec(3) = bz1d(i)
                    tempvec = matmul(tempvec,RotR)
                    bx1d(i) = tempvec(1) 
                    by1d(i) = tempvec(2) 
                    bz1d(i) = tempvec(3)                
                    
                    ! Rotate sensitivities as well:
                    if (bDoPartials) then
                        do j=1,nlay1D 
                            
                            tempvec(1)     = dexdsig(i,j)
                            tempvec(2)     = deydsig(i,j)
                            tempvec(3)     = djzdsig(i,j)
                            tempvec        = matmul(tempvec,RotR)
                            dexdsig(i,j) = tempvec(1) 
                            deydsig(i,j) = tempvec(2) 
                            djzdsig(i,j) = tempvec(3) 
                            
                            tempvec(1)     = dbxdsig(i,j)
                            tempvec(2)     = dbydsig(i,j)
                            tempvec(3)     = dbzdsig(i,j)
                            tempvec        = matmul(tempvec,RotR)
                            dbxdsig(i,j) = tempvec(1) 
                            dbydsig(i,j) = tempvec(2) 
                            dbzdsig(i,j) = tempvec(3)                       
                        enddo
                    endif
                    
                endif   ! site has rotations
                
            enddo ! i=1,n1d Rotation loop
                    

!
! Lastly, if computing Jacobians, convert df/dsigma to df/dlog10(rho):
!
            if (bDoPartials) then
                forall( j=1:nlay1D ) 
                    
                    ! Convert E's
                    dexdsig(:,j) = -1.d0*dlog(10.d0)*sig1d(j)*dexdsig(:,j) 
                    deydsig(:,j) = -1.d0*dlog(10.d0)*sig1d(j)*deydsig(:,j) 
                    djzdsig(:,j) = -1.d0*dlog(10.d0)*sig1d(j)*djzdsig(:,j) 
                    
                    ! Convert B's
                    dbxdsig(:,j) = -1.d0*dlog(10.d0)*sig1d(j)*dbxdsig(:,j) 
                    dbydsig(:,j) = -1.d0*dlog(10.d0)*sig1d(j)*dbydsig(:,j) 
                    dbzdsig(:,j) = -1.d0*dlog(10.d0)*sig1d(j)*dbzdsig(:,j) 
                    
                end forall
            endif

!
! Now loop through all the data and pull out model responses computed by the current call to Dipole1D
! and insert them into dm and wj
!
! This is a giant block that handles many different data types:
!
! JH Split into two loops and refactored slightly.

        do i = 1,nd
            
            ! Skip rest of loop if dp(i,:) does not use current iFreq and iTx:
            if (  ( dp(i,2) .ne. iFreq) .or.  ( dpTfm(i,1) .ne. iTx) ) cycle
                   
            ! Get the transformed receiver number,
            ! Then get mapping of that parameter into Rx subset passed to Dipole1D: 
            iRx = RxInd( dpTfm(i,2))                  

            select case (  dp(i,1) ) ! what data type is this?
                                              
            ! Electric: Real and Imaginary Data:
                    
                    case (indRealEx)  
                      dm(i) = dble(ex1d(iRx))
                    case (indImagEx)  
                      dm(i) = aimag(ex1d(iRx))
                    case (indRealEy)  
                      dm(i) = dble(ey1d(iRx))
                    case (indImagEy)  
                      dm(i) = aimag(ey1d(iRx))
                    case (indRealEz)  
                      dm(i) = dble(jz1d(iRx)) !/sigsite(iRx)
                    case (indImagEz)  
                      dm(i) = aimag(jz1d(iRx)) 
                        
            ! Electric: Amplitude and Phase Data:
                               
                    case (indAmpEx)   
                      dm(i) = abs(ex1d(iRx))    !     |Ex|
                    case (indAmpEy)   
                      dm(i) = abs(ey1d(iRx))    !     |Ey|
                    case (indAmpEz)   
                      dm(i) = abs(jz1d(iRx))    !     |Ez|
                    case (indPhsEx) 
                      dm(i) = getPhase(ex1d(iRx),dble(d(i)))
                    case (indPhsEy) 
                      dm(i) = getPhase(ey1d(iRx),dble(d(i)))
                    case (indPhsEz) 
                      dm(i) = getPhase(jz1d(iRx),dble(d(i)))

            ! Magnetic: Real and Imaginary Data:

                    case (indRealBx) 
                      dm(i) = dble(bx1d(iRx))   !        real(Bx)
                    case (indImagBx) 
                      dm(i) = aimag(bx1d(iRx))  !        imag(Bx)
                    case (indRealBy) 
                      dm(i) = dble(by1d(iRx))   !        real(By)
                    case (indImagBy) 
                      dm(i) = aimag(by1d(iRx))  !        imag(By)
                    case (indRealBz) 
                      dm(i) = dble(bz1d(iRx))   !        real(Bz)
                    case (indImagBz) 
                      dm(i) = aimag(bz1d(iRx))  !        imag(Bz)
                        
            ! Magnetic: Amplitude and Phase Data:
                        
                    case (indAmpBx) 
                      dm(i) = abs(bx1d(iRx))   !     |Bx|
                    case (indAmpBy) 
                      dm(i) = abs(by1d(iRx))   !     |By|
                    case (indAmpBz) 
                      dm(i) = abs(bz1d(iRx))   !     |Bz|
                    case (indPhsBx) 
                      dm(i) = getPhase(bx1d(iRx),dble(d(i)))
                    case (indPhsBy) 
                      dm(i) = getPhase(by1d(iRx),dble(d(i)))
                    case (indPhsBz) 
                      dm(i) = getPhase(bz1d(iRx),dble(d(i)))
                      
               
             ! Polarization Ellipse Parameters:
             
                     case (iPEmax)
                       dm(i) = getPE(ex1d(iRx),ey1d(iRx),'pmax') ! E-field max
                     case (iPEmin)
                       dm(i) = getPE(ex1d(iRx),ey1d(iRx),'pmin') ! E-field min
                     case (iPBmax)
                       dm(i) = getPE(bx1d(iRx),by1d(iRx),'pmax') ! B-field max
                     case (iPBmin)
                       dm(i) = getPE(bx1d(iRx),by1d(iRx),'pmin') ! B-field min
                    
                        
            end select  ! case dp(i,1)
        
          enddo  ! loop over nd


          if( bDoPartials) then  ! calculate the partial derivatives as well.
 
            do i=1,nd
              if (  ( dp(i,2) .ne. iFreq) .or.  ( dpTfm(i,1) .ne. iTx) ) cycle

              ! Get the transformed receiver number, mapped into Dipole1D Rx subset.
              iRx = RxInd( dpTfm(i,2))                  

               select case (  dp(i,1) )  ! what data type is this?
                  
               ! Electric: Real and Imaginary Data:
                          
                  case (indRealEx)  
                    wj(i, :) = dble(dexdsig(iRx,ipm2ilay))    !        real(Ex)
                  case (indImagEx)  
                    wj(i, :) = aimag(dexdsig(iRx,ipm2ilay))   !        imag(Ex)
                  case (indRealEy)  
                    wj(i, :) = dble(deydsig(iRx,ipm2ilay))    !        real(Ey)
                  case (indImagEy)  
                    wj(i, :) = aimag(deydsig(iRx,ipm2ilay))   !        imag(Ey)
                  case (indRealEz)  
                    wj(i, :) = dble( djzdsig(iRx,ipm2ilay))   !        real(Ez)
                  case (indImagEz)  
                    wj(i, :) = aimag( djzdsig(iRx,ipm2ilay))  !        imag(Ez)
                                  
               ! Electric: Amplitude and Phase Data:
                                     
                  case (indAmpEx)   
                    wj(i, :) = absDeriv( ex1d(iRx), dexdsig(iRx,ipm2ilay) )   !     |Ex|
                  case (indAmpEy)   
                    wj(i, :) = absDeriv( ey1d(iRx), deydsig(iRx,ipm2ilay) )   !     |Ey|
                  case (indAmpEz)   
                    wj(i, :) = absDeriv( jz1d(iRx), djzdsig(iRx,ipm2ilay) )   !     |Ez|
                  case (indPhsEx) 
                    wj(i, :) = phaseDeriv( ex1d(iRx), dexdsig(iRx,ipm2ilay) )   !     Ex phase 
                  case (indPhsEy)
                    wj(i, :) = phaseDeriv( ey1d(iRx), deydsig(iRx,ipm2ilay) )   !     Ey phase 
                  case (indPhsEz) 
                    wj(i, :) = phaseDeriv( jz1d(iRx), djzdsig(iRx,ipm2ilay) )   !     Ez phase 

               ! Magnetic: Real and Imaginary Data:

                  case (indRealBx) 
                    wj(i, :) = dble(dbxdsig(iRx,ipm2ilay))   !        real(Bx)
                  case (indImagBx) 
                    wj(i, :) = aimag(dbxdsig(iRx,ipm2ilay))  !        imag(Bx)
                  case (indRealBy) 
                    wj(i, :) = dble(dbydsig(iRx,ipm2ilay))   !        real(By)
                  case (indImagBy) 
                    wj(i, :) = aimag(dbydsig(iRx,ipm2ilay))  !        imag(By)
                  case (indRealBz) 
                    wj(i, :) = dble(dbzdsig(iRx,ipm2ilay))   !        real(Bz)
                  case (indImagBz) 
                    wj(i, :) = aimag(dbzdsig(iRx,ipm2ilay))  !        imag(Bz)
                                  
               ! Magnetic: Amplitude and Phase Data:
                              
                  case (indAmpBx) 
                    wj(i, :) = absDeriv( bx1d(iRx), dbxdsig(iRx,ipm2ilay) )   !     |Bx|
                  case (indAmpBy) 
                    wj(i, :) = absDeriv( by1d(iRx), dbydsig(iRx,ipm2ilay) )   !     |By|
                  case (indAmpBz) 
                    wj(i, :) = absDeriv( bz1d(iRx), dbzdsig(iRx,ipm2ilay) )   !     |Bz|
                  case (indPhsBx) 
                    wj(i, :) = phaseDeriv( bx1d(iRx), dbxdsig(iRx,ipm2ilay) )   !     Bx phase 
                  case (indPhsBy) 
                    wj(i, :) = phaseDeriv( by1d(iRx), dbydsig(iRx,ipm2ilay) )   !     By phase 
                  case (indPhsBz) 
                    wj(i, :) = phaseDeriv( bz1d(iRx), dbzdsig(iRx,ipm2ilay) )   !     Bz phase 
             
               ! Polarization Ellipse Parameters:   
               
                  case (iPEmax)
                    wj(i, :) = getPEDeriv(ex1d(iRx),ey1d(iRx),dexdsig(iRx,ipm2ilay),deydsig(iRx,ipm2ilay),'pmax') ! E-field max
                  case (iPEmin)
                    wj(i, :) = getPEDeriv(ex1d(iRx),ey1d(iRx),dexdsig(iRx,ipm2ilay),deydsig(iRx,ipm2ilay),'pmin') ! E-field min
                  case (iPBmax)
                    wj(i, :) = getPEDeriv(bx1d(iRx),by1d(iRx),dbxdsig(iRx,ipm2ilay),dbydsig(iRx,ipm2ilay),'pmax') ! B-field max
                  case (iPBmin)
                    wj(i, :) = getPEDeriv(bx1d(iRx),by1d(iRx),dbxdsig(iRx,ipm2ilay),dbydsig(iRx,ipm2ilay),'pmin') ! B-field min
   
   
               end select  ! case dp(i,1)

            end do

          endif !bDoPartials
          
          deallocate(iRxlayer)
          
        enddo ! loop over nFreq
    enddo ! loop over nTx
    

    end subroutine computeFwd_CSEM
    

!==============================================================================!  
!===========================================================! computeCSEM_OPRA !  
!==============================================================================!     
    subroutine computeCSEM_OPRA( )
!
! Subroutine to compute the orthogonal Procrustes rotation between
! the model vector and data vector for each vector datum.  
! 
! Rotates the original data to the OPRA angles. Model response is unperturbed
!
!  Assumes data ordered as
!  real Ex
!  imag Ex
!  real Ey
!  imag Ey
!  real Ez
!  imag Ez
!
! It searches for real(Ex) and then assumes the others immediately follow.
!
! Assumes that the standard errors are isotopic among the x,y,z components (ie, all the same). 
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 2.0    August-September, 2008.    Added full 3D rotation solver.
! Version 1.0.   April, 2008.               Horizontal rotation solver only.
!
!
    use Occam
    use csem1d_mod
    
    implicit none
    
!
! Local variables:
!
    logical                         :: lTilt            ! this is used to check number of # data for rotation is correct
    integer                         :: i,n,iRx, ndrot, irotcounter, info
    real(8)                            :: stderr, rot, detrxy, alpha, beta, work(100)
    real(8), dimension(3,3)            :: zxy,uxy,vxy,sxy,rxy, ss, RotR,Rotx,Roty,Rotz 
    real(8),dimension(:,:),allocatable :: dxyz, mxyz

    real(8) :: ccr, ssr, cca, ssa, ccb, ssb
    real(8) :: tempr(3)
    
! Assumes data ordered as
!  real Ex
!  imag Ex
!  real Ey
!  imag Ey
! So that rotations can be easily compute without complicated sorting
! Error on data should be:  %error * | (Ex,Ey) | ; i.e., some percentage of the vector amplitude
! Same ordering applies to B field data
!

  !  write(*,*) 'computeCSEM_OPRA :',nRxRotSolve
    irotcounter = 0
    
    if (nRxRotSolve > 0) then
    !
    ! 1. Loop over nRx
    !
        do iRx = 1,nRxIn
            
            if (.not.lRxRotSolve(iRx)) cycle ! skip receiver

            irotcounter  = irotcounter + 1
            lTilt = .false.
    !
    !   2. Count # x,y,z data pairs in dp for current Rx:
    !
            n = 0
            do i = 1,nd
                ! KWK debug, need to put a check in here to make sure data is in correct order
                if ( dp(i,4) == iRx) then  ! Current receiver
                
                    if ( ( dp(i,1) == indRealEx ).or.( dp(i,1) == indImagEx ).or. &       ! Horizontal E's
                       & ( dp(i,1) == indRealEy ).or.( dp(i,1) == indImagEy ) ) then
                       n = n + 1   
                    endif
                    
                    if ( ( dp(i,1) == indRealBx ).or.( dp(i,1) == indImagBx ).or. &     ! Horizontal B's
                       & ( dp(i,1) == indRealBy ).or.( dp(i,1) == indImagBy ) ) then
                       n = n + 1   
                    endif                   
                    
                    if ( ( dp(i,1) == indRealEz).or.( dp(i,1) == indImagEz ) ) then    ! Vertical E's
                        n = n + 1   
                        lTilt = .true.
                    endif 
                    
                    if ( ( dp(i,1) == indRealBz ).or.( dp(i,1) == indImagBz ) ) then    ! Vertical B's
                       n = n + 1   
                       lTilt = .true.
                    endif    

                endif
            enddo
            
            ! 
            ! Now check n: it should be divisible by 4 (x,y) or 6 (x,y,z).  If not, the data array isn't complete
            ! 
            if (lTilt) then
                ndrot = n / 6
                if (mod(n,6) /= 0 ) then
                    write(*,*) ' '
                    write(*,*) ' !!! Error in computeCSEM_OPRA !!!'
                    write(*,*) '  '
                    write(*,*) ' Number of data is not divisible by 6, as required for 3D rotation solution! '
                    write(*,*) ' Make sure the data have real(Ex),Imag(Ex),Real(Ey),Imag(Ey),Real(Ez),Imag(Ez) '
                    write(*,*) ' and similarly if including magnetic data '
                    write(*,*) ' '
                    write(*,*) ' Stopping. '
                    stop
                
                endif
            else
                ndrot = n / 4
                if (mod(n,4) /= 0 ) then
                    write(*,*) ' '
                    write(*,*) ' !!! Error in computeCSEM_OPRA !!!'
                    write(*,*) '  '
                    write(*,*) ' Number of data is not divisible by 4, as required for horizontal rotation solution! '
                    write(*,*) ' Make sure the data have real(Ex),Imag(Ex),Real(Ey),Imag(Ey) '
                    write(*,*) ' and similarly if including magnetic data '
                    write(*,*) ' '
                    write(*,*) ' Stopping. '
                    stop
                    
                endif
            endif

    !
    !     3. Allocate arrays for dxy (n x 3) and dmxy (n x 3)
    !
            if (ndrot>0) then
                
                allocate ( dxyz(3,ndrot*2), mxyz(ndrot*2,3) ) ! *2 for real, imag rws, note dxy is actually transpose(dxy) 
                
                ! initialize them to 0:
                dxyz = 0d0
                mxyz = 0d0
         
    !            
    !         4. Loop through dp and form dxy, dmxy arrays
    !
                n = 1
                do i = 1,nd
                
                  if ( dp(i,4) == iRx) then  ! Current receiver
                
                      
                     if ( ( dp(i,1) == indRealEx).or.( dp(i,1) == indRealBx) ) then ! real Ex | Bx  
                         
                       ! Get all the data and model responses, push it into the arrays:
                         dxyz(1,n)    = d_copy(i)   ! real Ex ! kwk debug, using dcopy og data array
                         dxyz(1,n+1)  = d_copy(i+1) ! imag Ex                         
                         dxyz(2,n)    = d_copy(i+2) ! real Ey  
                         dxyz(2,n+1)  = d_copy(i+3) ! imag Ey        
                         
                         stderr = max(sd(i),sd(i+2)) ! max of x,y errors
                         
                         if (lTilt) then
                            dxyz(3,n)    = d_copy(i+4) ! real Ez 
                            dxyz(3,n+1)  = d_copy(i+5) ! imag Ez 
                            stderr = max(stderr,dble(sd(i+4))) ! max of x,y,z errors
                         endif
                         

                         mxyz(n,1)    = dm(i)   ! real Ex
                         mxyz(n+1,1)  = dm(i+1) ! imag Ex                         
                         mxyz(n,2)    = dm(i+2) ! real Ey  
                         mxyz(n+1,2)  = dm(i+3) ! imag Ey                         
                         
                         if (lTilt) then
                            mxyz(n,3)   = dm(i+4) ! real Ez   
                            mxyz(n+1,3) = dm(i+5) ! imag Ez 
                         endif
                                           

                         ! Normalize by standard error (larger of x,y,z components).  
                         ! Here i've assumed that real/imag errors are based on amplitude of complex field and are the same
                       
                         dxyz( 1:3   , n:n+1 ) = dxyz( 1:3   , n:n+1 ) / stderr
                         mxyz( n:n+1 , 1:3   ) = mxyz( n:n+1 , 1:3   ) / stderr
                          
                         n = n + 2 ! increment n, the vector data counter, adding a row of Real and a Row of Imag data to array
                                            
                    endif  
                  endif
               enddo    
    
    !
    !         5. Compute OPRA
    !

              ! Compute  D'*M:
                zxy = matmul(dxyz,mxyz)
                
                deallocate ( dxyz, mxyz) 
                 
                
                select case (svdMethod)
                
                case('lapack_svd')
                    
                    ! Lapack SVD Routine:
                    call dgesvd('A','A',3,3,zxy,3,sxy,uxy,3,vxy,3,work ,100, info)
                    vxy = transpose(vxy)
                
                case('parker_svd')
                
                    ! Parker's SVD routine:
                    ! Note that this routine doesn't always give the correct SVD:
                    ! 
                    call svd(3,3,3,3, zxy, uxy, vxy, sxy, 1)
                    
                case default
                    write(*,*) 'CSEM1D.f90 error, undefined svdMethod variable! Stopping!'
                    stop
                end select
                
                uxy = transpose(uxy)
                rxy = matmul(vxy,uxy)  
                
                detrxy = ( rxy(1,1)*rxy(2,2)*rxy(3,3) + rxy(1,2)*rxy(2,3)*rxy(3,1) + rxy(1,3)*rxy(2,1)*rxy(3,2) ) - &
                       & ( rxy(3,1)*rxy(2,2)*rxy(1,3) + rxy(3,2)*rxy(2,3)*rxy(1,1) + rxy(3,3)*rxy(2,1)*rxy(1,2) )
                        
                ss = 0
                ss(1,1) = 1
                ss(2,2) = 1
                ss(3,3) = detrxy
               ! write(*,*) detrxy
            
                
                rxy = matmul(ss,uxy)
                rxy = matmul(vxy,rxy)

                rot   = rad2deg*atan2(  rxy(2,1) , rxy(1,1))
                alpha = rad2deg*asin(   rxy(3,1) )
                beta  = rad2deg*atan(  rxy(3,2) / rxy(3,3))  ! atan forces |beta| < 90 degrees
                
                !write(*,*) detrxy, rot, alpha, beta
             !
             ! Store rotations
             !
                
                 pm_assoc (3*(irotcounter-1) + 1)  = rot
                 pm_assoc (3*(irotcounter-1) + 2)  = alpha
                 pm_assoc (3*(irotcounter-1) + 3)  = beta
                 
    !           6. Rotate the data vectors
    !  
    !  Note: In the first version of OPRA, I rotated the model response and sensitivities, but then on forward only calls
    !        by Occam the model response and sensitivities can have different rotations. So now i'm rotating the data
    !        and the inversions converge much better.
    
                     ccr = cos(deg2rad* rot)
                     ssr = sin(deg2rad* rot)
                     cca = cos(deg2rad* alpha)
                     ssa = sin(deg2rad* alpha)
                     ccb = cos(deg2rad* beta)
                     ssb = sin(deg2rad* beta)
                     
                     Rotx = 0
                     Rotx(1,1) = 1
                     Rotx(2,2) =  ccb
                     Rotx(2,3) = -ssb
                     Rotx(3,2) =  ssb
                     Rotx(3,3) =  ccb
                     
                     Roty = 0
                     Roty(1,1) =  cca  
                     Roty(1,3) = -ssa
                     Roty(2,2) = 1 
                     Roty(3,1) =  ssa
                     Roty(3,3) =  cca
                     
                     
                     Rotz = 0
                     Rotz(1,1) =  ccr
                     Rotz(1,2) = -ssr
                     Rotz(2,1) =  ssr
                     Rotz(2,2) =  ccr
                     Rotz(3,3) =  1    
                     
                     RotR = matmul(Roty,Rotx)
                     RotR = matmul(Rotz,RotR)
                     
                     RotR = transpose(RotR) ! transpose for multiplying Data array
                     
                 
                !     
                ! Loop through and extract rotated data:
                !
                    
                    do i = 1,nd
                      if ( dp(i,4) == iRx) then  ! Current receiver
                
                        if ( ( dp(i,1) == indRealEx) .or. ( dp(i,1) == indRealBx)) then ! real Ex | Bx
                      ! Rotate original unrotated data in d_copy:
                            tempr = 0
                            tempr(1) = d_copy(i)
                            tempr(2) = d_copy(i+2)
                            if (lTilt) then
                                tempr(3) = d_copy(i+4)
                            endif
                            tempr = matmul(tempr,RotR)
                            
                            d(i)   = tempr(1)
                            d(i+2) = tempr(2)
                            if (lTilt) then
                                 d(i+4) = tempr(3)
                            endif       
                            
                          ! Now rotate imag Exy data:
                            tempr = 0
                            tempr(1) = d_copy(i+1)
                            tempr(2) = d_copy(i+3)
                            if (lTilt) then
                                tempr(3) = d_copy(i+5)
                            endif
                            
                            tempr = matmul(tempr,RotR)
                            
                            d(i+1) = tempr(1)
                            d(i+3) = tempr(2)
                            if (lTilt) then
                                 d(i+5) = tempr(3)
                            endif          
                            
                         ! KWK debug: need to rotate the standard errors?
                         ! For now let's assume the error bars are isotropic
                            
                 
                        endif ! Ex or Bx
                        
                    endif !  if ( dp(i,4) == iRx) then  ! Current receiver
               
                enddo

            else 
                write(*,*) ' No rotation data for current Rx! iRx = ',iRx
            endif   ! n > 0
            
        enddo ! Loop over nRxIn  
    endif

    end subroutine computeCSEM_OPRA
    
!==============================================================================!  
!==============================================================! computeFwd_MT !  
!==============================================================================! 
    subroutine computeFwd_MT( bDoPartials, ipm2ilay )
!
! Compute responses and sensivities for 1DMT data:
!

    use Occam
    use csem1d_mod  ! layer_params
    use dipole1d
    use mt1d_mod 
    use dataTransformations
    
    implicit none

!
! Arguments
!
    logical, intent(in) :: bDoPartials  ! logical flag for computing Jacobians
    integer,intent(in)  :: ipm2ilay(nParams) 
!
! Local Variables
!
    integer :: i, j,nAllocErr

!
! We assume that pm values have already been inserted into sig1d array the main computeFwd_CSEM subroutine
! These are then used by MT1D
!
    if (bDoPartials) then
        allocate (  dCRespdRho(nlay1D),   stat=nAllocErr )
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Too many 1D layers  (',nlay1D, ')'
            stop 
        endif           
    endif
    
!
! Loop over data points and compute MT responses and sensitivities
! 
    do i = 1,nd
    
        !   TYPE   FREQ#   TX#    RX#     DATA  SD_ERROR
        
        
        if ( dp(i,1) > 100) then
            
            zRx_MT = zRxIn(dp(i,4))  ! get the Rx depth
            freqMT = fTxIn(dp(i,2))  ! get the frequency
            
           ! 
           ! Compute the MT response (note MT1D gives the C response, not the impedance Z):
           ! 
            call MT1D(bDoPartials)  
            
           !
           ! If computing Jacobians, convert dC/rho to dC/dlog10(rho):
           !     
            if (bDoPartials) then
                forall( j=1:nlay1D )   
                    dCRespdRho(j) = dCRespdRho(j) * dlog(10.d0) / sig1d(j)  
                end forall
            endif
            
           !
           ! Convert C response to impedance Z:
           !
            CResp = cmplx(0,1)*7.8957E-06*freqMT*CResp
            if (bDoPartials) then
                dCRespdRho= cmplx(0,1)*7.8957E-06*freqMT*dCRespdRho
            endif       
            
           ! 
           ! Decode the MT response data:
           !
            select case ( dp(i,1) )
                
                case (indRhoZXY,indRhoZYX)   !     Apparent resistivity
                    
                    dm(i) = getAppRes(CResp,freqMT) 
                    
                    if (bDoPartials) then
                        wj(i, :) =  getAppResDeriv(CResp,dCRespdRho(ipm2ilay),freqMT)
                    endif
                    
                case (indPhsZXY,indPhsZYX)   !    Phase
                        
                    dm(i) = getPhase(CResp,dble(d(i)) ) 
                    
                    ! KWK debug: should add check for yx phase at -135 degrees.
                    
                    if (bDoPartials) then
                        wj(i, :) = phaseDeriv(CResp, dCRespdRho(ipm2ilay) )
                    endif   
    
                 case (indRealZXY)          
                 
                    dm(i) = real(CResp)
                    
                    if (bDoPartials) then
                        wj(i, :) = real(dCRespdRho(ipm2ilay)) 
                    endif
                    
                 case (indImagZXY)
                 
                    dm(i) = aimag(CResp)
                    
                    if (bDoPartials) then
                        wj(i, :) = aimag(dCRespdRho(ipm2ilay)) 
                    endif               
                    
                 case (indRealZYX)  ! In 1D, Zyx = -Zxy
                 
                    dm(i) = -real(CResp)
                    
                    if (bDoPartials) then
                        wj(i, :) = -real(dCRespdRho(ipm2ilay)) 
                    endif                    
                    
                 case (indImagZYX)
                
                    dm(i) = -aimag(CResp)
                    
                    if (bDoPartials) then
                        wj(i, :) = -aimag(dCRespdRho(ipm2ilay)) 
                    endif   
                    
            end select  
            
        endif

    enddo
    
    if (bDoPartials) then
        deallocate (  dCRespdRho )  
    endif
     
    end subroutine computeFwd_MT
 
!==============================================================================!  
!=======================================================================! MT1D !  
!==============================================================================!  
    subroutine MT1D(bDoPartials)
!
! MT1D computes the 1D MT response as complex C and optionally the partial 
! derivatives with respect to the RESISTIVITY.
!
! KWK: Modified to Fortran90 from MTLAY and MTDIV routines in 
! the original Occam1DMT.  
!
! MTLAY COMPUTES THE COMPLEX C VALUE OVER A LAYERED EARTH
! S.C. CONSTABLE, S.I.O. LA JOLLA CA92093, AUGUST 1986
!
! MTDIV COMPUTES THE COMPLEX C VALUE OVER A LAYERED EARTH AND DC/DRHO FOR
!   EACH LAYER (I.E. A ROW OF THE JACOBIAN MATRIX)
! S.C. CONSTABLE, S.I.O. LA JOLLA CA92093, AUGUST 1986
!
! Modified for usage in CSEM1D by Kerry Key, 2008.
! 
    use dipole1d   ! this module has the 1D model as nlay1D, siglay1D and zlay1D (model top depth)
    use mt1d_mod      ! inputs zRx_MT, freqMT, and solution is passed out as CResp, dCRespdRho(nlay1d)
    
    implicit none
!
! Arguments
!
!
    logical, intent(in) :: bDoPartials  ! logical flag for computing Jacobians  
    
!
! Local variables
!
    real(8)     :: omega, thick
    complex(8)  :: omi, ck, c
    complex(8)  :: CKM1,CKP1,OPEXP,OMEXP,DCIIP1,THEEXP,RI
    integer     :: iRxLayer, i, j


!
! First get the layer with the MT receiver:
!
                        
    do j= 1,nlay1D
        if (zRx_MT.ge.zlay1D(j))  then   ! Sites on boundaries use the layer below, impedance computed at top of that layer.
            iRxlayer = j
        endif
    enddo
    
 ! kwk debug: need to add in support for MT sites not on a layer boundary...is this ever done?
    
!
! Compute the 1DMT response, assumes receiver is on top boundary of that layer
!

        
    OMEGA = 6.283185308*freqMT                  ! w = 2 PI F
    OMI = OMEGA*1.256637062E-06*(0.,1.)     ! OMI = i w mu
    
    C = 1./SQRT(OMI*sig1D(nlay1d))          ! C = 1 / sqrt( i w mu sig)
    DO  I = nlay1d-1,iRxlayer,-1                    ! run up the stack of layers computing impedance
        CK = SQRT(OMI*sig1D(i))
        thick  =  zlay1d(i+1) - zlay1d(i)
        C = EXP(-2.*CK*THICK)*(CK*C-1.)/(CK*C+1.)
        C = (1.+C)/(1.-C)/CK
    ENDDO
    cResp = c
!
! Compute the partial derivatives with respect to the layer resistivity:
!
    if (bdopartials) then
! Initialize:
    dCRespdRho = 0
    
! FIRST OF ALL DC/DRHO FOR I=1,2,..NLAYER WILL BE STORED IN D()
! THEN PROPAGATED TO THE SURFACE USING DC(I)/DC(I-1) AT EACH LAYER.
!
! ESTABLISH THE VALUES FOR THE TOP OF THE TERMINATING HALF-SPACE
        CK = SQRT(OMI*sig1D(nlay1d) )   
        C = 1./CK
        dCRespdRho(nlay1d) = 1./(CK/sig1d(nlay1d)*2.)
! NOW SWEEP UP THROUGH THE LAYERS
        DO  I = nlay1d-1,iRxlayer,-1
            CK = SQRT(OMI*sig1D(i))
            CKP1 = C*CK + 1.
            CKM1 = C*CK - 1.
            RI = CKM1/CKP1
            thick  =  zlay1d(i+1) - zlay1d(i)
            THEEXP = EXP(-2.*CK*cmplx(thick))
            OPEXP = THEEXP*RI
            OMEXP = 1. - OPEXP
            OPEXP = 1. + OPEXP
! COMPUTE DC/DRHO FOR THE (I)TH LAYER
            dCRespdRho(I) = 2.*THEEXP*(C/CKP1**2 - cmplx(thick)*RI)/OMEXP**2
            dCRespdRho(I) = OPEXP/OMEXP/CK/2. - dCRespdRho(I)
            dCRespdRho(I) = dCRespdRho(I)*sig1D(i)
! COMPUTE DC(I)/DC(I+1) FOR THE (I)TH LAYER
            DCIIP1 = 4.*THEEXP/(CKP1*OMEXP)**2
! NOW COMPUTE THE VALUE OF C AT THE TOP OF THIS LAYER
            C = OPEXP/OMEXP/CK
! MULTIPLY THE PRODUCT SUMS FOR ALL THE LAYERS BELOW THIS ONE, WHICH HAVE
!   ALREADY BEEN PARTIALLY ACCUMULATED.
            DO J = I+1,nlay1d
                dCRespdRho(J) = dCRespdRho(J)*DCIIP1
            enddo
        enddo

    endif
    
    end subroutine MT1D  

!==============================================================================!  
!=============================================================! WRITE_RESPONSE !  
!==============================================================================!  
    subroutine writeResponse( )
!
! Subroutine to output a model response file for the Occam1DCSEM code.
! This routine is placed in the CSEM1D interface to Occam since the 
! response files will have data parameters that Occam knows nothing about. 
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 1.0.   March, 2008.
!
    use Occam 
    use csem1d_mod
    use dipole1d
    
    implicit none
    
    character(50)   :: cNum
    integer :: lerr, i, j, ncolumns, irotcounter
    
!
! Write out Receiver orientation solutions if requested:
!
    if ( nRxRotSolve > 0) then
    
    !
    ! Write them to the screen:
    !
      write(*,*) ' Estimated Receiver Rotations : '
      write(*,'(5(a12,2x))') 'Receiver #', 'theta', 'alpha', 'beta', '[degrees]'
      irotcounter = 0
       do i=1,nRxIn
          if (lRxRotSolve(i)) then
            irotcounter = irotcounter+1
            write(*,'(i12,2x,3(f12.2,2x))') i,pm_assoc ( 3*(irotcounter-1)+1:3*irotcounter) 
         endif
       enddo 
        
    !    
    ! Write them to a file:
    !
        write (cNum,*) nCurrentIter
        !write (cNum,'(I4.3)') nCurrentIter      ! DGM put leading zeros so files sort properly
        open (unit=21, file=trim(cRootName)//'_'//trim(adjustl(cNum))//'.rotations', iostat=lerr)
    
    !
    ! Catch error opening response file:
    !
        if (lerr /= 0) then
            ! we could be in the middle of an iteration, so don't die on error, but
            ! send a message and forge on regardless
            write(*,*) ' Error opening associated parameter response file'
            return
        end if
    !
    ! No error, write out the responses:
    !    
!           do i=1,nRxRotSolve
!              write(21,*) IndRxRotSolve(i),pm_assoc ( 3*(i-1)+1:3*i) 
!           enddo 
        irotcounter = 0
        do i=1,nRxIn
          if (lRxRotSolve(i)) then
            irotcounter = irotcounter+1
            write(21,'(i12,2x,3(f12.2,2x))') i,pm_assoc ( 3*(irotcounter-1)+1:3*irotcounter) 
         endif
        enddo     
        close(21)
    
    
    
    endif
    

!----------------------------
!    
! Open response file
!
    write (cNum,*) nCurrentIter
    !write (cNum,'(I4.3)') nCurrentIter      ! DGM put leading zeros so files sort properly
    open (unit=21, file=trim(cRootName)//'_'//trim(adjustl(cNum))//'.resp', iostat=lerr)

!
! Catch error opening response file:
!
    if (lerr /= 0) then
        ! we could be in the middle of an iteration, so don't die on error, but
        ! send a message and forge on regardless
        write(*,*) ' Error opening response file'
        return
    end if
!
! No error, write out the responses:
!    

! New Feb 2010:  EMResp_1.2, now has dipole length
    write(21,*) 'Format:    EMResp_1.2'
    write(21,'(a18,1x,g13.5)') 'Dipole length: ',  lenTx1D 
    write(21,*) '# Transmitters:', nTxIn
    write(21,'(a1,5(a12,1x))') '!','X','Y','Z','Azimuth','Dip'   
    do i = 1,nTxIn
       write(21,'(1x,5(f12.1,1x))') xTxIn(i),yTxIn(i),zTxIn(i),azimuthTxIn(i),dipTxIn(i)
    enddo
    
    write(21,*) '# Frequencies:',   nFreqIn
    do i = 1,nFreqIn
        write(21,'(g13.5)') fTxIn(i)
    enddo
    
    write(21,*) '# Receivers:',   nRxIn
    write(21,'(a1,6(a12,1x))') '!','X','Y','Z','Theta','Alpha','Beta'
    do i = 1,nRxIn
       write(21,'(1x,6(f12.1,1x))') xRxIn(i),yRxIn(i),zRxIn(i), ThetaRxIn(i),AlphaRxIn(i),BetaRxIn(i)
    enddo
    
    
    write(21,*) '# Data:',   nd
    write(21,'(a1,4(a12,1x),4(a13,1x))') '!','Type','Freq#','Tx#','Rx#','Data','StdError','Response','Residual' 
    ncolumns = size(dp,2)
    do i = 1,nd
       write(21,'(1x,4(i12,1x),3(g13.5,1x),f13.2)') (dp(i,j),j=1,ncolumns), d(i),sd(i), dm(i), (d(i)-dm(i))/sd(i)
    enddo   
    close(21)
    
    end subroutine writeResponse
    
!==============================================================================!  
!========================================================! constructPenaltyMat !  
!==============================================================================!     
    subroutine constructPenaltyMat( del, nPTerms, bHaveLink, linkpr)
!
! Constructs the penalty matrix for the 1D layered model
!
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 2.0    March 2010     Added option for weighted L2 norm (aka minimum gradient support)
! Version 1.0    April, 2008
!
! Input:
!   nPTerms:  (INTEGER) number of penalty terms
!   bHaveLink: Logical flag for forming the linkpr array (carryover from 2D code, not really needed here)
!
! Output:
!   del:      (REAL) penaly matrix (nPTerms x nParams)
!   linkpr:   optional, (2 x nPTerms), brick #s of pairs of bricks linked in the penalty matrix.
!           Used in deltdel to calculate (DEL)Transpose * DEL using
!           essentially sparse matrix methods. Carryover from 2D code.
!
!-----------------------------------------------------------------------
    use Occam 
    use csem1d_mod
    use dipole1d     

    implicit none
    
    integer, intent(in)            :: nPTerms    ! # of penalty terms (for dim of DEL & linkpr)
    real(RealPrec), intent(out)    :: del(nPTerms,nParams)
    logical, intent(in)            :: bHaveLink  ! true if linkpr was passed.  NB: intrinsic fctn present() doesn't appear to work.
    integer, optional, intent(out) :: linkpr(nPTerms,2)        ! DGM Oct 2006 - not always needed.  Quite large!
    
    character(180) :: sValue, sExtra
!
! Local variables:
! 
    integer        :: i, nPen, nPar
    real(RealPrec) :: penalt, delta
    real(8)        :: z0

!------------------------------------------------

!
! Loop over model layers and add penalties where free parameters exist
! 
    del = 0.0     ! Array math
    nPen = 0
    nPar = 0
 
    ! DGM March 2010 - Add roughness type 5: sum of model values (log10 so
    ! this is essentially penalizing distance from 1 Ohm m).  SC's idea.
    if( iRufType == eSumParam ) then
        ! DGM 10/2010 - allow use of the weighting param in the input file.
!        do i=1,nParams
!            del(i,i) = 1
!        enddo
!        ! Don't need the linkpr stuff, but *must* fill it in because
!        ! Occam expects it. So just point to the first pair. It is used
!        ! in forming: (del)T * del
!        if (bHaveLink) then
!            linkpr(:,1) = 1
!            linkpr(:,2) = 2
!        endif
        do i=1,nLay1D
            if( layer_params(i,2) < 0 ) then    ! neg = free param
                nPar = nPar + 1
                del(nPar,nPar) = real(layer_params(i,3))
                if( bHaveLink ) linkpr(nPar,1:2) = nPar
            endif
        enddo
        
        return
    endif
    
    ! DGM Oct 2010 - Roughness type 7: attempt at 1D anisotropy by
    !   having the layer-to-layer penalties alternate layers. This
    !   in effect makes two model spaces interspersed with one
    !   another. For anisotropy, one can be conductive and the other
    !   resistive. For vertical rho, pairs of layers are resistors
    !   in series. For horiz rho, pairs of layers are resistors in
    !   parallel. Hopefully this will develop naturally and not just
    !   be completely unstable....
    if( iRufType == eTwoModel ) then
        do i = 1, nLay1D
            if( layer_params(i,2) < 0 ) then    ! neg = free param
                nPar = nPar + 1
                if( nPar < 3 ) cycle
                nPen = nPen + 1
                del(nPen,nPar)  = -real(layer_params(i,3))
                del(nPen,nPar-2)=  real(layer_params(i,3))
                if( bHaveLink ) then
                    linkpr(nPen,1) = nPar
                    linkpr(nPen,2) = nPar-2
                endif
                
                ! If a linkage has been requested between the two
                ! normally independent models, then setup the penalties.
                ! This is used to place some weight against anistropy.
                if( rRufParam .ne. 0.0 ) then
                    nPen = nPen + 1
                    del(nPen,nPar)  = -rRufParam
                    del(nPen,nPar-1)=  rRufParam
                    if( bHaveLink ) then
                        linkpr(nPen,1) = nPar
                        linkpr(nPen,2) = nPar-1
                    endif
                endif
            endif
        enddo
        
        return
    endif
    
!
! If depth weighted roughness, get top depth of first free layer:
! 
    if( iRufType == eDepthWt ) then
        do i=1,nlay1D
            if (layer_params(i,2) < 0) then
                z0 = layer_params(i,1)
                exit
            endif
        enddo
    endif
    
!
! Loop over 1D model layers and add roughness penalty terms for free layers
!
    do i=1,nlay1D
    
        if (layer_params(i,2) < 0) then ! negative resistivity is free parameter    
            nPar = nPar + 1
            if (nPar == 1 ) cycle 
            
            !
            ! Set the penalty weight:
            !
            select case (iRufType)
            case (e1stDeriv)        ! first difference penalty
                penalt = real(layer_params(i,3))
                
            !case (e2ndDeriv)        ! 2nd difference penalty
            !   penalt = 
            
            case (eDepthWt)         ! log10(depth) weighted penalty
                penalt = real(layer_params(i,3)*log10( layer_params(i,1) - z0 ))
                !   write(*,*) layer_params(i,3),layer_params(i,1),layer_params(1,1)
                
            case (eMinGrad)         ! 1st diff with focusing (minimum gradient support)
                ! from KK email March 2010 - approximate minimum gradient support, using last iteration's
                ! model parameters for weighting:
                penalt =  real(layer_params(i,3))  / sqrt( (pm(npar) - pm(npar-1) )**2 + rRufParam**2 )
                ! KWK: note rRufParam**2 keeps rRufParam in units of resistivity so that its value on input
                ! is easy to understand: i.e. the cutoff resistivity jump for the weights to become smaller than unity.
            case default
                write(*,*) ' Error. Roughness type not supported: ', iRufType
                write(*,*) ' stopping.'
                stop
                
            end select
            
            !
            ! Now insert the penalty across the layer boundary, ie between layers (i) and (i-1) 
            !
            nPen = nPen + 1
            del(nPen, nPar) = -1. * penalt
            del(nPen, nPar-1) =  1. * penalt
            
            if (bHaveLink) then
                linkpr(nPen,1) = nPar
                linkpr(nPen,2) = nPar-1
               ! write(*,*) nPen,linkpr(nPen,1:2),del(nPen, nPar-1:nPar)
            endif
        
        endif ! (layer_params(i,2) < 0)
        
    enddo ! i=1,nlay1D

    return
    
    end subroutine constructPenaltyMat

!==============================================================================!  
!==========================================================! countPenaltyTerms !  
!==============================================================================!     
    integer function countPenaltyTerms(  ) result(nTerms)
!
! Counts the number of penalty terms in the roughness matrix.
! This routine is a carryover from the more complicated 2D Occam.
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!-----------------------------------------------------------------------

    use Occam 
    use csem1d_mod
    use dipole1d 
    
    implicit none
    
! DGM 10/2010 - Is there some reason a loop was used when we already know?
!   For most roughness types, The # of terms is the # of params minus 1.
    select case( iRufType )
    ! DGM March 2010 - Add roughness type 5: sum of model values (log10 so
    ! this is essentially penalizing distance from 1 Ohm m).  SC's idea.
    case( eSumParam )
        nTerms = nParams
    case( eTwoModel )
        nTerms = nParams - 2
        if( rRufParam .ne. 0.0 ) nTerms = nTerms * 2
    case default
        nTerms = nParams - 1
    end select
    return

    end function countPenaltyTerms


!-----------------------------------------------------------------------
    subroutine svd(mdim,ndim,m, n, a, u, v, q, index)
!-----------------------------------------------------------------------
   
    implicit none
!$$$$$ calls no other routines
!  singular value decomposition)  for algol program see wilkinson+reinsch
!  handbook for automatic computation vol 2 - linear algebra, pp140-144
!  translated from algol by r.l.parker
!  the matrix a(m,n) is decomposed.  singular values in q, pre-matrix in u,
!  post-matrix in v.   index may be 1,2,3 or 4.  if 1, find u,v. if 2, find
!  only u. if 3, find only v. if 4, find neither. in all cases, the array  u
!  must be supplied as it is used as working space for the routine.
      integer :: mdim,ndim,m,n,index
      real(8) :: a,u,v,q
      dimension a(mdim,ndim),u(mdim,ndim),v(ndim,ndim),q(ndim)
      
      real(8) :: eps,tol
      
      real(8), dimension(:), allocatable :: e
      
      integer :: i,j,k,l, iback, kback, lback, l1, lplus
      
      real(8) :: c,f,g,h,s,x,y,z
      
      allocate (e(n))
      
!
!
!
      eps=1.0e-12
      tol=1.0e-35
      do 1100 i=1,m
      do 1100 j=1,n
 1100 u(i,j)=a(i,j)
!  householder reduction to bi-diagonal form
      g=0.0
      x=0.0
      do 2900 i=1,n
      e(i)=g
      s=0.0
      l=i+1
      do 2100 j=i,m
 2100 s=u(j,i)**2 + s
      if (s .lt. tol) go to 2500
      f=u(i,i)
!      g=-dsign(dsqrt(s),f)                                              
      g=- sign( sqrt(s),f)                                              
      h=f*g - s
      u(i,i)=f - g
      if (l.gt.n) go to 2501
      do 2400 j=l,n
      s=0.0
      do 2200 k=i,m
 2200 s=u(k,i)*u(k,j) + s
      f=s/h
      do 2300 k=i,m
 2300 u(k,j)=u(k,j) + f*u(k,i)
 2400 continue
      go to 2501
 2500 g=0.0
!
 2501 continue
      q(i)=g
      s=0.0
      if (l.gt.n) go to 2601
      do 2600 j=l,n
 2600 s=u(i,j)**2 + s
 2601 if (s.lt.tol) go to 2800
      f=u(i,i+1)
!      g=-dsign(dsqrt(s),f)                                              
      g=- sign( sqrt(s),f)                                              
      h=f*g - s
      u(i,i+1)=f - g
      if (l.gt.n) go to 2651
      do 2650 j=l,n
 2650 e(j)=u(i,j)/h
 2651 continue
      if (l.gt.m) go to 2850
      do 2700 j=l,m
      s=0.0
      if (l.gt.n) go to 2700
      do 2670 k=l,n
 2670 s=u(j,k)*u(i,k) + s
      do 2690 k=l,n
 2690 u(j,k)=u(j,k) + s*e(k)
 2700 continue
      go to 2850
 2800 g=0.0
! 2850 y=dabs(q(i)) + dabs(e(i))                                         
2850  y= abs(q(i)) +  abs(e(i))                                         
      if (y .gt. x) x=y
 2900 continue
!
!  accumulation of right-hand transforms (v)
!
      go to (3000,3701,3000,3701       ),index
 3000 continue
      do 3700 iback=1,n
      i=n+1-iback
      if (g .eq. 0.0) go to 3500
      h=u(i,i+1)*g
      if (l.gt.n) go to 3500
      do 3100 j=l,n
 3100 v(j,i)=u(i,j)/h
      do 3400 j=l,n
      s=0.0
      do 3200 k=l,n
 3200 s=u(i,k)*v(k,j) + s
      do 3300 k=l,n
 3300 v(k,j)=v(k,j) + s*v(k,i)
 3400 continue
 3500 continue
      if (l.gt.n) go to 3601
      do 3600 j=l,n
      v(j,i)=0.0
 3600 v(i,j)=0.0
 3601 v(i,i)=1.0
      g=e(i)
      l=i
 3700 continue
 3701 continue
!
!  accumulation of left-hand transforms
      go to (4000,4000,4701,4701       ),index
 4000 continue
      do 4700 iback=1,n
      i=n+1-iback
      l=i+1
      g=q(i)
      if (l.gt.n) go to 4101
      do 4100 j=l,n
 4100 u(i,j)=0.0
 4101 if (g.eq. 0.0) go to  4500
      h=u(i,i)*g
      if (l.gt.n) go to 4401
      do 4400 j=l,n
      s=0.0
      do 4200 k=l,m
 4200 s=u(k,i)*u(k,j) + s
      f=s/h
      do 4300 k=i,m
 4300 u(k,j)=u(k,j) + f*u(k,i)
 4400 continue
 4401 continue
      do 4550 j=i,m
 4550 u(j,i)=u(j,i)/g
      go to 4700
 4500 continue
      do 4600 j=i,m
 4600 u(j,i)=0.0
 4700 u(i,i)=u(i,i) + 1.0
!
!  diagonalization of bi-diagonal form
 4701 eps=eps*x
      do 9000 kback=1,n
      k=n+1-kback
!  test f-splitting
 5000 continue
      do 5100 lback=1,k
      l=k+1-lback
!      if (dabs(e(l)).le. eps) go to 6500                                
      if ( abs(e(l)).le. eps) go to 6500                                
!      if (dabs(q(l-1)) .le. eps) go to 6000                             
      if ( abs(q(l-1)) .le. eps) go to 6000                             
 5100 continue
!  cancellation of e(l), if l.gt. 1
 6000 c=0.0
      s=1
      l1=l - 1
      do 6200 i=l,k
      f=s*e(i)
      e(i)=c*e(i)
!      if (dabs(f) .le. eps) go to 6500                                  
       if ( abs(f) .le. eps) go to 6500                                  
      g=q(i)
!      q(i)=dsqrt(f*f + g*g)                                             
      q(i)= sqrt(f*f + g*g)                                             
      h=q(i)
      c=g/h
      s=-f/h
      go to (6050,6050,6200,6200       ),index
 6050 continue
      do 6100 j=1,m
      y=u(j,l1)
      z=u(j,i)
      u(j,l1)=y*c + z*s
      u(j,i)=-y*s + z*c
 6100 continue
 6200 continue
!  test f-convergence
 6500 z=q(k)
      if (l .eq. k) go to  8000
!  shift from bottom 2 x 2 minor
      x=q(l)
      y=q(k-1)
      g=e(k-1)
      h=e(k)
      f=((y-z)*(y+z) + (g-h)*(g+h))/(2.0*h*y)
!      g=dsqrt(f*f + 1.0)                                                
      g= sqrt(f*f + 1.0)                                                
!      f=((x-z)*(x+z) + h*(y/(f + dsign(g,f))-h))/x                      
      f=((x-z)*(x+z) + h*(y/(f +  sign(g,f))-h))/x                      
!  next q-r transformation
      c=1.0
      s=1.0
      lplus=l + 1
      do 7500 i=lplus,k
      g=e(i)
      y=q(i)
      h=s*g
      g=c*g
!      z=dsqrt(f*f + h*h)                                                
      z= sqrt(f*f + h*h)                                                
      e(i-1)=z
      c=f/z
      s=h/z
      f=x*c + g*s
      g=-x*s + g*c
      h=y*s
      y=y*c
      go to (7100,7201,7100,7201       ),index
 7100 do 7200 j=1,n
      x=v(j,i-1)
      z=v(j,i)
      v(j,i-1)=x*c + z*s
      v(j,i)=-x*s + z*c
 7200 continue
! 7201 z=dsqrt(f*f + h*h)                                                
 7201 z= sqrt(f*f + h*h)                                                
      q(i-1)=z
      c=f/z
      s=h/z
      f=c*g + s*y
      x=-s*g + c*y
      go to (7300,7300,7500,7500       ),index
 7300 do 7400 j=1,m
      y=u(j,i-1)
      z=u(j,i)
      u(j,i-1)=y*c + z*s
      u(j,i)=-y*s + z*c
 7400 continue
 7500 continue
      e(l)=0.0
      e(k)=f
      q(k)=x
      go to  5000
!  convergence
 8000 if (z .ge. 0.0) go to 9000
!  q is made non-negative
      q(k)=-z
      go to (8100,9000,8100,9000       ),index
 8100 do 8200 j=1,n
 8200 v(j,k)=-v(j,k)
 9000 continue
 
      deallocate (e) 
       
      return
      end         
      
      
!-----------------------------------------------------------------------
subroutine parseFields( nLen, sLine, nFields, sFields)
!-----------------------------------------------------------------------
!
! Routine to parse out mixed numeric and character fields from a line
! This is useful for mixed format input tables, which are awkward for Fortran
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
! Version 1.0.   November 21, 2008.  kwk debug: Needs some error checking.
!
    implicit none
    
    ! Args
    integer, intent(in)     :: nLen
    character(nLen)         :: sLine
    integer, intent(in)     :: nFields
    character(nLen), intent(out) :: sFields(nFields)
   
    
    ! Local vars
    integer :: iFrom, iTo, i
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    iFrom = 1
    
    ! Loop through the line and get the nFields:
    do i=1,nFields
    
        ! Skip any beginning blanks
        do iFrom = iFrom,nLen
            if (sLine(iFrom:iFrom) .ne. ' ') exit
        enddo
        
        ! Pull out nonblank character string:
        do iTo = iFrom,nLen
             if (sLine(iTo:iTo) .eq. ' ') exit
        enddo
        sFields(i) = trim(sLine(iFrom:iTo-1))
        iFrom = iTo
    enddo  
    
    end subroutine parseFields     
    
!==============================================================================!  
!==============================================================! deallocateFwd !  
!==============================================================================!      

    subroutine deallocateFwd
!
! Routine from DGM for deallocating all arrays in FWD code, useful for when 
! Occam1DCSEM is called many times from an external interface.
! 
    use csem1d_mod
    use dipole1d
    
    implicit none
    
!! readData
    deallocate( fTxIn )                                     ! frequency list
    deallocate( xTxIn, yTxIn, zTxIn, azimuthTxIn, dipTxIn ) ! transmitter lists
    deallocate( xRxIn, yRxIn, zRxIn, ThetaRxIn, AlphaRxIn &
              , BetaRxIn, lRxRotSolve )                     ! receiver lists
    deallocate( d_copy, sd_copy )                           ! data & std copy for OPRA

!! LookForRxTxShortcut
    ! Variables used in Rx/Tx reciprocity transformation
    deallocate( dpTfm, xTxTfm, yTxTfm, zTxTfm, azimuthTxTfm, dipTxTfm &
              , xRxTfm, yRxTfm, zRxTfm, iRxTfm &
              , sigsite, RxInd, rx1DtoRxIn )
    ! Field return variables
    deallocate( x1D, y1D, z1D &             ! from dipole1d
              , ex1D, ey1D, jz1D &
              , bx1D, by1D, bz1D )
    
!! readModel
    deallocate( layer_params )
    deallocate( zlay1D, sig1D )             ! from dipole1d

!! readModel_PartTwo
    deallocate( dexdsig, deydsig, djzdsig & ! from dipole1d
              , dbxdsig, dbydsig, dbzdsig )

end subroutine deallocateFwd
