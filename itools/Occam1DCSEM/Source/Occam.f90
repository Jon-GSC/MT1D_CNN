!=======================================================================
! Occam's Inversion Fortran2003 Code  
!=======================================================================
!
!    Copyright (C) 1986-2011
!    Steven Constable, Kerry Key, David Myer, Catherine deGroot-Hedlin
!    Scripps Institution of Oceanography
!    University of California, San Diego
!
!    This file is part of Occam's Inversion.
!
!    Occam's Inversion is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Occam's Inversion is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with Occam's Inversion.  If not, see <http://www.gnu.org/licenses/>.
!
!-----------------------------------------------------------------------
! Revision History
!-----------------------------------------------------------------------
!
! These revisions refer to the Occam's inversion code, which is distinctly
! separate from the specific problem under consideration (EM, seismic, gravity, etc).
!
! Revision  Date                Authors     Notes 
!
! 3.10      March 22, 2011      KK          Added support for model bounds using nonlinear transformations.  Supports the Habashy 
!                                           and Abubakar (2004) exponential transformation and a new one that I created to have a
!                                           flat sensitivity scaling in the pass-band between the bounds. Both seem to work well but 
!                                           have only been moderately tested here so far.
!
! 3.9       March 2, 2011       DM, KK      Fix to allow Min Norm inversion to have 2nd stage in Occam procedure.
!                                           Fixed Lapack calls to use correct precision routines depending on RealPrec parameter. 
!
! 3.8       March 26, 2010      KK          Changed integer iruf parameter to character(180) roughnessPenalty.
!
! 3.7       Oct 19, 2009        KK          Minor changes, add linearSolver parameter for specifying Lapack or Parker's 
!                                           cholesky codes or Lapack LU.
!
! 3.6       Aug 2009            DGM         Rearranged to allow for MatLab interface.  Also added deallocation 
!                                           routine so that non-occam memory allocated by readData, readModel,
!                                           etc... can be properly deallocated before the process ends.  While
!                                           this is not important for the executable, it is critical for the 
!                                           MatLab interface.  Quite a bit of the Occam1DCSEM memory was being leaked.
! 
! 3.5       May 8, 2009         KK          Fixed a bug that was crashing Occam on all 64-bit systems and some 32 bit compilers. 
!                                           My sloppy programming for version 3.3 left a few functions as single precision, rather 
!                                           than declaring them real(RealPrec).  Also completely rearranged the subroutines
!                                           into a more object oriented fashion:  (1) Created module Occam and made public and
!                                           private variables, (2) moved all subroutines inside the module as contained subroutines
!                                           and set many of them to private.  This helps encapsulate memory and shields 
!                                           from outside prying eyes. Rewrote the instructions below to reflect these new changes.  
!                                           Occam_inteface_module has been eliminated.
!
! 3.4       May 1, 2009         KK, DM      Added some bounds on mu in tofmu.  If mu is outside the bounds,
!                                           the tofmu returns with an artifically large misfit 
!                                           Bounds: ( 1d-6 < mu < 1d8 ).  Upper: tofmu = mu, lower: tofmu = 1/mu
!                                           This helps keep mu away from infinity. 
!
! 3.3       December 2008       KK          Switched to Fortran2003 intrinsic commands for command line arguments
!                                           Added RealPrec variable in module Occam_RealPrecision for so Occam
!                                           can be easily complied in single or double precision.
!
! 3.2       November 2008       KK          Implemented J. Gunning's suggestion to move interface module 
!                                           to a separate file to make it easier for interfacing with 
!                                           other languages.
!
! 3.1       April-Sept 2008     KK          Added pass-through "associated parameters" for the best fitting model. 
!                                           Split occam_mod into internal and interface portions to keep working parameters
!                                           isolated from forward routines. Added module Occam_Timers. Modified datime 
!                                           routine to use Fortran90/95 intrinsic date_and_time command. Renamed some of the 
!                                           subroutines for clarity.
!                                           Added programming notes explaining how to apply the Fortran90/95 Occam for your 
!                                           own inverse problems.
!
! 3.0       Oct/Nov 2006        DM, KK      Big re-write! Moved to F90, allocated memory, cmdline params,
!                                           modules, and general code reorg.  Flexible file input file formats.
!                                           
! 
! 2.1       July 2001           CDH         deltdel for faster roughness matrix multiply
!
! SC: incredibly stable for 9 years!
!
! 2.0       April 1992          SC,CDH      A re-write to make OCCAM independent of the dimensionality
!                                           of the forward problem and make maximum use of dynamic memory allocation.
!                                           The calling program is now an integral part of the package since all the
!                                           model-dependent stuff has been removed.  (MORE THANKS TO CATHERINE AND ALSO GEOTOOLS) 
! 1.5.2     September 1989      SC
! 1.5       January 1989        SC,CDH
! 1.4       August 1986         SC
! 1.3       October 1987        SC  
! 1.2       March 1986          SC
!
! Authors: 
!
! SC    Steve Constable
! CDH   Catherine deGroot-Hedlin
! DM    David Myer
! KK    Kerry Key
!
!-----------------------------------------------------------------------
!  References:
!-----------------------------------------------------------------------
!
!  Constable, S. C., R. L. Parker, and C. G. Constable, 1987, 
!  Occam s inversion - A practical algorithm for generating smooth models 
!  from electromagnetic sounding data: Geophysics, 52, 289 300.
!
!  deGroot-Hedlin, C. and S. Constable, 1990, Occam s inversion to 
!  generate smooth two-dimensional models from magnetotelluric data: Geophysics, 55, 1613 1624.
!
!  Constable, S., 1991, Comment on  Magnetic appraisal using simulated annealing  by S. E. 
!  Dosso and D. W. Oldenburg: Geophys. J. Int., 106, 387 388.
!
!  Constable, S., 1992: Occam Distribution Notes.
!   
!  Myer et al., 2007: OCCAM 3.0 release notes.
!
!  Key, K., 2009, One-dimensional inversion of multi-component, 
!  multi-frequency marine CSEM data: Methodology and synthetic studies for 
!  resolving thin resistive layers: Geophysics, Vol.74,No.2, March-April 2009,
!  doi:  10.1190/1.3058434 
!
!  Key, K., and A. Lockwood, 2010, Determining the orientation of marine 
!  CSEM receivers using orthogonal Procrustes rotation analysis: Geophysics,
!  75(3), F63-F70. doi:10.1190/1.3378765
!
! See http://marineemlab.ucsd.edu/Projects/Occam/
!
!-----------------------------------------------------------------------
!  Usage:
!-----------------------------------------------------------------------
! 
!  Occam2D [-F] [<startupfile>] [<outputprefix>] 
! 
! Optional Parameters:
!   -F              Only compute the forward response of the startupfile.
!   <startupfile>   Name of the startup file. 
!   <outputprefix>  Prefix to apply to all output from Occam2D.  You 
!                   must specify a startupfile parameter in order to 
!                   specify this parameter.  Without this parameter, model 
!                   iteration files will be named iterxx.iter and model 
!                   response files will be named respxx.resp, where xx is the
!                   Occam iteration number. The Occam run log will be called 
!                   logfile.logfile.  If an output prefix is specified, 
!                   then the names will be <outputprefix>xx.iter, 
!                   <outputprefix>xx.resp, and <outputprefix>.logfile.
!
!-----------------------------------------------------------------------
! Programming Notes:
!-----------------------------------------------------------------------
!
! These notes apply to Version 3.7 of the Occam's inversion engine.
!
! Around 1992 Steve Constable re-wrote Occam so that the inversion engine was distinctly separate from the 
! specific geophysical problem being considered. Thus, it became relatively easy to apply Occam to various methods. 
! The Occam code was then migrated to Fortran90/95/2003 in 2006-2009 by David Myer and Kerry Key.  While the inversion 
! algorithm code remains largely unchanged, the implementation of Fortran90/95/2003 constructs allow for much simpler 
! interfacing with Occam and allow for dynamic memory allocation, memory encapsulation, etc.   
!
! Occam.f90 is controlled by RunOccam.f90, which reads in only a single file, the iteration file,
! which by default is called 'startup' on the first iteration. The iteration file lists a data file, a model file, 
! as well as the many other parameters used by Occam,  including the starting values for the inversion model parameters
! of interest.  It is up to the programmer to supply the following external routines listed below, which are specific 
! to the problem at hand.  Data passing from Occam to the forward routines is accomplished entirely through the module 
! 'Occam', which is what this file contains.
! 
!
! Occam Interface:
!
!   module Occam   A module that contains the core Occam routines (private) and public data shared with the forward routine.
!
!   *Note that all real parameters in Occam are determined precision variable RealPrec, defined near the top of module Occam.
!
!
! External subroutines called by Occam.f90:
!
!    The routines should read in problem specific parameters as well, but you must store those in your own external modules.
!    Below I only describe the parameters required by the Occam's inversion engine.  
!
!   readData()          Problem specific data file input routine that reads in the data arrays. Also can read any other 
!                       data dependent parameters but these need to be stored in the Forward code modules.
!
!                       Uses:           module Occam
!
!                       Assigns:        nd          number of data
!                                       d           vector of data
!                                       sd          vector of data errors
!                                       dp          array of data parameters (not used by Occam)
!                                       npm_assoc   number of associated model or data parameters (interesting things associated 
!                                                   with the best fitting model, but not used in the misfit computation). 
!                                                   Set to 0 if your problem doesn't have associated parameters.
!                                       pm_assoc    Vector of associated model parameters. These are optional and just pass through
!                                                   Occam, which keeps track of the pm_assoc vector for the best fitting model.
!                                                   For example, these could be estimated site orientations, static shifts, 
!                                                   temperatures or porosity predicted by conductivity, or anything else 
!                                                   associated with the current model parameters, but not used in the misfit 
!                                                   computation.
!
!                       Allocates:      d(nd), sd(nd), dp(nd, # data parameters), pm_assoc(npm_assoc)
!
!
!   readModel()         Problem specific model file input routine.  This should read in the parameters that describe the 
!                       model layers or blocks, and anything else needed for performing the Forward and Sensitivity 
!                       computations for the specific problem at hand.  The parameters that describe the model geometry
!                       and any fixed structure should be stored in a module outside of Occam.  Occam only needs to
!                       get the preference model and its weights from readModel.  
! 
!
!                       Uses:           module Occam 
!
!                       Assigns:        premod(nParams)     A vector of preference values for each model parameter.
!                                       prewts(nParams)     A vector of weights to apply to the preference model parameters.
!                                                           Set these to 0 if no preferences are desired.  These are allocated
!                                                           within Occam in subroutine readIteration, so you just need to assign 
!                                                           their values externally.
!                                       
!                       Allocates:      Nothing used in Occam.f90.  nParams is read in readIteration so Occam arrays are already 
!                                       allocated within Occam.f90.
!
!
!   computeFwd(bDoPartials, CurrentModParams)     Computes the current model's forward response and optionally the model 
!                                                 sensitivities.      
!
!                       Uses:           model Occam
!
!                       Input:          logical :: bDoPartials      True if requesting Jacobian (sensitivity) matrix computations.
!                                       
!                                       real(RealPrec), dimension(nParams) :: CurrentModParams  The vector of model parameters to 
!                                                                                               compute the forward response for.  
!
!                       Assigns:        dm(nd)              The forward response of the current model
!                                       wj(nd,nParams)      Optional. The Jacobian matrix (if bDoPartials is true).
!                                       pm_assoc(npm_assoc) Optional. The associated model parameters.
!                                                           These are optional and just pass through  Occam, which keeps track of 
!                                                           the pm_assoc array for the best fitting model.
!                                                           For example, these could be estimated site orientations, static shifts,
!                                                           temperatures or porosity predicted by conductivity, or anything else 
!                                                           associated with the current model parameters.   
!
!   constructPenaltyMat()  Constructs the sparse penalty matrix for the 1D problem.  
!
!   countPenaltyTerms()    Function that counts the total number of penalty terms so that arrays can be allocated.
!                         
!
!   writeResponse()     Occam always writes out an iteration file for each iteration, but you need to supply a routine
!                       that outputs the model response.  This writeResponse subroutine is called at the end of each Occam iteration
!                       and can be used to output the current model responses in module Occam,  as well as any other custom output 
!                       you desire.   This can be strongly problem depenendent and so is best left outside of Occam.f90.
!
!   deallocateFwd() subroutine to allow you to clean up your memory allocations before exit.
!
!==============================================================================!    
!===============================================================! OCCAM MODULE ! 
!==============================================================================!
module Occam
    
    implicit none
    
!--------------
! Public data:
!--------------
! These can be accessed outside Occam by placing a "use Occam" statement in your subroutine.
!
    integer, parameter, public :: RealPrec = kind(0d0)   ! Double precision:  RealPrec = kind(0d0) 
                                                         ! Single precision:  RealPrec = kind(0e0)  
    ! This controls the real variable precision used in Occam, which doesn't have to be the same precision used in the forward code
    ! to compute the responses and model derivatives.
    ! Note that for most problems there will probably be little difference in the Occam results for single or double precision.
    ! Single precision will reduce the memory requirements for the arrays below. So if you're close to maxing out your cpu's memory,
    ! and need more unknowns, try using single precision above.
    !

    integer, public       :: nCurrentIter     ! Current Occam iteration number
    character(180), public :: cRootName        ! Iteration fileroot name. Defaults to 'ITER' if not given as a command line argument
    character(180), public :: dataFileName     ! Data file name
    character(180), public :: modelFileName    ! Model file name    
    
    ! DGM Aug 2009, changed to a variable from a parameter so the 
    ! MatLab interface can set to 6 - equivalent to write(*,*) which
    ! MatLab ignores.
    integer, public                                     :: ioUnitOccamLogFile
  
    integer, public                                     :: nd        ! number of data
    real(RealPrec), dimension(:), allocatable, public   :: d         ! data array d(nd) 
    real(RealPrec), dimension(:), allocatable, public   :: sd        ! standard error array sd(nd)      
    integer, dimension(:,:), allocatable, public        :: dp        ! data parameters (nd x nDataParams) lookup table
                                                                     ! of parameter indices (frequncies, positions, etc)
    integer, public                                     :: nParams   ! number of model parameters
    real(RealPrec), dimension(:), allocatable, public   :: pm        ! model parameter array pm(nParams)
    real(RealPrec), dimension(:), allocatable, public   :: dm        ! model response vector dm(nd)
    real(RealPrec), dimension(:,:), allocatable, public :: wj        ! model response Jacobian matrix
    real(RealPrec), dimension(:), allocatable, public   :: premod    ! model preference values premod(nParams)
    real(RealPrec), dimension(:), allocatable, public   :: prewts    ! weights for model preference prewts(nParams). Set to 0 for 
                                                                     ! no preference   
    integer, public                                     :: npm_assoc ! Number of associate (pass through) parameters
    real(RealPrec), dimension(:), allocatable, public   :: pm_assoc  ! Assocated parameter arrays pm_assoc(npm_assoc)
                                                                  
    ! Note: associated data/model parameters for the current iteration  are 
    ! are things that are derived from the current model and/or data
    ! but are not used for the model roughness or misfit computations.
    ! These are only along for the ride in the Occam routines and are entirely
    ! dealt with in the Forward code.  They need to be included here
    ! so Occam can keep track of pm_assoc for the best fitting model 
    ! of each iteration.  Then the associated parameters can be written to file 
    ! at the end of each Occam iteraion with the call to subroutine writeResponse.
    
    real(RealPrec), public                              :: rDeltaMisfit   
    ! DGM Sept 2009: if misfit changes by < rDeltaMisfit, then end so user can 
    ! examine residual & edit data as necessary.
    
    ! DGM 12/1/2010: roughness type has had a make-over by 2 people simultaneously.
    ! I added some new types and the ability to have 2 numbers, whereas KKey changed
    ! it to strings for roughness type. So I've combined the two. Strings are allowed
    ! in the iter file, but are translated into type numbers which are very fast
    ! to determine at runtime and are now easily readable through the enumeration below.
    enum, bind(c)
        enumerator :: e1stDeriv = 1     ! 1st derivative penalty
        enumerator :: e2ndDeriv = 2     ! 2nd deriv -- Not all Occam uses support this
        enumerator :: eDepthWt  = 4     ! 1st deriv * log10(depth)
        enumerator :: eSumParam = 5     ! DGM 3/2010 sum of model params (e.g. log10(rho))
                                        ! In 1D, tends to make thin sheet models (!!)
        enumerator :: eMinGrad  = 6     ! DGM 3/2010 Minimum gradient support (aka weighted L2)
                                        ! In 1D, tends to make really CRAPPY models.
            ! Needs a 2nd number in the iter file - the stablizing value used in the
            ! denominator of: real(layer_params(i,3)) / sqrt ( (pm(npar) - pm(npar-1) )**2 + delta**2)
            ! Use a space or comma like this: 
            !               Roughness Type:   MinimumGradientSupport 0.01
        enumerator :: eTwoModel = 7     ! DGM 10/2010: 1D anistropy attempt 
                            ! Applies 1st deriv on alternating layers.
                            ! Should allow two different models to develop simultaneously.
                            ! For vertical rho, take mean of every pair of layers (resistors in series)
                            ! For horiz rho, treat each layer pair as resistors in parallel.
            ! Optional 2nd number: the penalty between the layers of the two models.
            ! If this is zero or missing, then the odd #'d layers will be penalized 1st deriv
            ! and the even #'d layers will be penalized 1st deriv and "two models" will have
            ! no penalty between them. This optional param is for penalty between the pairs
            ! of even & odd #'d layers (e.g. 1 & 2, 3 & 4, 5 & 6, etc...)
            !   While this alternate layer approach is stable, without the optional penalty
            ! it tends to introduce anisotropy everywhere - even for synthetic non-aniso data.
            ! I've found that setting this value to between 0.1 and 0.3 reduces and often 
            ! eliminates that.
    end enum
    ! Text values allowed (case ignored, but spaces are important):
    ! 1     FirstDiff
    ! 4     DepthWeighted
    ! 6     mgs, weightedL2, MinimumGradientSupport
    ! 7     TwoModel
    
    integer(kind(e1stDeriv)), public                    :: iRufType
    real(RealPrec), public                              :: rRufParam
    ! DGM Oct 2010: some new roughness types want an additional scalar from
    ! the startup/iter file to be used in the forward modules countPenaltyTerms
    ! and constructPenaltyMat. Rather than carrying yet another local variable
    ! all over creation, use this nifty module.

!
!---------------
! Private data:
!---------------
! The private attribute means these can only be used within the contained subroutines below.
!
    integer, parameter, private                          :: ioUnitOccamIterFile = 21    
    real(RealPrec), dimension(:,:), allocatable, private :: ptp   ! product of penalty matrix and its transpose    
    real(RealPrec), dimension(:,:), allocatable, private :: wjtwj ! product of weighted jacobian matrix and its transpose
    real(RealPrec), dimension(:), allocatable, private   :: wjtwd ! product of weighted jacobian matrix and weighted translated data
    real(RealPrec), dimension(:), allocatable, private   :: pwk1, pwk2, pwk3 ! vectors to carry around working copies of models
    real(RealPrec), dimension(:), allocatable, private   :: dwk1, dwk3  ! working copies of the model responses,'dm' is used as dwk2
    real(RealPrec), dimension(:), allocatable, private   :: dwk1_assoc, dwk2_assoc, dwk3_assoc ! working associated parameters arrys
    real(RealPrec), private                              :: frac        ! current reduction in step size
    integer, private                                     :: nFracTimes  ! # of times to cut stepsize before giving up
    integer, private                                     :: idebug      ! current debug level
    integer, private                                     :: nfor        ! tally for total number of forward calculations
    real(RealPrec), dimension(:,:), allocatable, private :: aMat        ! working arrays
    real(RealPrec), dimension(:), allocatable, private   :: pwk4, dwk4  ! working arrays
    
    character(80), parameter, private :: linearSolver = 'lapack_cholesky'  ! Specify the linear solver for Ax=b in tofmu.  
                                                                           ! Options are:
                                                                           ! 'lapack_cholesky', 'lapack_lu', 'parker_cholesky'  
        
    ! DGM Nov 2006 - Force model to have finite discrete steps:    
    logical, private        ::  gbModelStepped    
    real(RealPrec), private ::  gnModelSteps

    ! KWK March 2011
    ! Non-linear transformation for constraining parameters with upper and lower bounds:
    logical, private        :: lConstrainParams
    real(RealPrec), parameter, private ::bandPassfactor = 15.  ! factor in exp(bandPassFactor/(a-b) * x)
    character(80), private  :: cBoundsTransform    ! 'exponential' or 'bandpass'. The default is set in intializeOccamParamLimits,
                                                  ! but can be modified by an optional parameter in the iteration file. 
                                                  ! 'exponential' - from Habashy and Abubakar (2004), Commer and Newman (2008)
                                                  ! 'bandpass' - designed by K. Key in March 2011.
                                                  !
    real(RealPrec), private :: lowerBound,upperBound

!---------------------------------------------------------------------
! Public subroutines called from program runOccam, outside the module:
!---------------------------------------------------------------------

    public :: intializeOccam, intializeOccamParamLimits
    public :: readIteration, writeIteration, compOccamIteration
    public :: allocateNdArrays, allocateAssocArrays, allocatenParamArrays
    public :: getOccamCommandLineArguments, datime, filerr, Lower, parseLine, parseCode, get_time_offset,openOccamLogFile
    public :: deallocateOccam
    
!--------------------------------------------------------
! Private subroutines only called from within this module:
!--------------------------------------------------------

    private :: scanmu, makjtj, makptp, deltdel, fndruf, minbrk, fminocc, froot, cholsl, cholin, trmult, atamul, anorm, tofmu

!
!------------------------
! Contained subroutines:
!------------------------
! Thes variables above have global scope within the contained subroutines below.
!
    contains
    
!-----------------------------------------------------------------------
! Public subroutines
!-----------------------------------------------------------------------
    subroutine openOccamLogFile()
    
    integer :: ierr
    
    ! DGM Aug 2009 - changed from param to variable for MatLab interface
    ioUnitOccamLogFile  = 20
    open (ioUnitOccamLogFile, file=trim(cRootName)//'.logfile', iostat=ierr)
    if (ierr .ne. 0) then
        write(*,*) ' Error opening log file'
        stop 
    end if
    
    end subroutine openOccamLogFile

    !-----------------------------------------------------------------------
    subroutine intializeOccam
    !-----------------------------------------------------------------------
        nCurrentIter    = 0
        nParams         = 0
        nFracTimes      = 8
        rDeltaMisfit    = 0.0
        idebug          = 0     
            ! write() status msgs do NOT work through MatLab! So default idebug
            !   off in this fctn. readIteration will set it to 1 for the 
            !   executable version of this code, so no worries.  The MatLab
            !   interface code can't touch this variable because it is private.
        
        iRufType        = e1stDeriv
        rRufParam       = 0.0
        
    end subroutine intializeOccam

    !-----------------------------------------------------------------------
    subroutine  get_time_offset(timein,timeout)
    !-----------------------------------------------------------------------    
    !    
    ! timein is the clock start time or 0.
    ! timeout is the time in seconds since the input time
    !
    ! Version 2.0  February 25, 2008  Now uses date_and_time Fortran intrinsic
    !
    ! Kerry Key
    ! Scripps Institution of Oceanography
    ! kkey@ucsd.edu
    
    
    integer, dimension(8) :: values
    integer               :: i,j,k,mjd
    
    real(8) :: timein, timeout, fracday
    
    !
    ! Fortran90 Time function:
    !
    call date_and_time(values=values) !this ouputs only values, ignoring other optional arguments
    
    ! Convert year, month day to modified julian day:
      
    i = values(1)
    j = values(2)
    k = values(3)
    mjd = -678927 + k + 1461*(i+(j-14)/12)/4 + 367*(j-2-12 * &
          & ((j-14)/12))/12 + (24002-12*i-j)/1200
    
               
    ! Add on fractional day:
                ! hour            ! minute          ! sec       ! millisec
    fracday = ((values(5)*60.d0 + values(6))*60.d0 + values(7) + values(8)/1000.d0 )/86400.d0
    
    timeout = mjd + fracday
    
    ! Finally, convert timeout to time difference between it and timein:  
    timeout =  timeout*86400.d0  - timein               
           
    end subroutine  get_time_offset
    
    
    !-----------------------------------------------------------------------
    subroutine readIteration (cStartup,descr,maxitr,tolreq,pmu,rlast,tobt,ifftol)      
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine Revision 3.0, November 2006
    ! Subroutine Revision 2.01, 13 Jan 1993
    ! DGM Nov 2006 - new startup file type (OCCAMITER_FLEX) allows
    !   flexible placement & specification of header items.
    !
    ! Reads a STARTUP or ITER file for OCCAM.
    !
    !

    
    ! arguments:
    integer maxitr, ifftol, lerr
    real(RealPrec) tolreq, pmu, rlast, tobt, pmtol
    character(50)  itform, descr, cTemp, cStartup
    ! local variables
    integer i 
  ! 
    
    character(180)  sLine, sCode, sValue, sValueToComma
    logical         bComment
    
    
    ! 
    ! Set default values to the file params so that any not in the file
    !   which are optional don't cause a crash.
    !
    call intializeOccam
    call intializeOccamParamLimits
    
    idebug      = 1     ! print status msgs by default.  Everyone likes them.
    !! return variables
    maxitr           = 100
    tolreq           = 1.0
    pmu              = 5.0
    rlast            = 1.0e+10
    tobt             = 1000
    ifftol           = 0
    
    write(*,*) 'Reading startup parameters from file: ', trim(cStartup)
    open (ioUnitOccamIterFile, file=trim(cStartup), status='old', iostat=lerr)
    if (lerr .ne. 0) then
        write(*,*) ' Error opening startup file'
        stop 
    end if    
    
    ! Read the file header.
    do while (.true.) !KWK nov 16,2006 eof not supported generally (.not. EOF(ioUnitOccamIterFile))
    ! Get the next code/value pair
    ! parseCode forces the code portion to be all lowercase with no
    !   padding and ending colon stripped off.  user comments are 
    !   stripped from the value portion.
    read( ioUnitOccamIterFile, '(A180)', end=198, err=199 ) sLine
    call parseCode( len(sLine), sLine, sCode, sValue, bComment )
    if( bComment ) cycle
    
    ! What do we have?
    select case (trim(sCode))
    
    case ('format')
        itform = trim(sValue)
        
    case ('description')
        descr = trim(sValue) 
        
    case ('model file')
        modelFileName = trim(sValue) 
        
    case ('data file')
        dataFileName = trim(sValue)
        
    case ('date/time')
        ! ignored!
        
    case ('max iter', 'iterations to run')
       read(sValue,*) maxitr 
        
    case ('req tol', 'target misfit')
        read(sValue,*) tolreq
        
    case ('iruf', 'roughness type')
        ! DGM 12/2010 If this is one of the new text entries KKey likes, decode it now.
        ! This is the ONLY PLACE these text values appear, so as new ones are added
        ! just put them here.
 
        ! KWK: fix for text input with comma:
        ! DGM allow space as well, this way if the user uses a numeric roughness
        ! type, Fortran doesn't mistake a comma for a numeric 1000's separator.
        i = index( sValue, ',' )
        if( i == 0 ) i = index( sValue, ' ' )
        if (i /= 0) then
            sValueToComma = sValue(1:i-1)
        else
            sValueToComma  = sValue
        endif
    
    
        call Lower(sValueToComma)
        select case (trim(adjustl(sValueToComma)))
        case( 'firstdiff' )
            iRufType = e1stDeriv
        case( 'depthweighted' )
            iRufType = eDepthWt
        case( 'mgs', 'minimumgradientsupport', 'weightedl2' )
            iRufType = eMinGrad
        case( 'twomodel' )
            iRufType = eTwoModel
        case default
            read(sValue,*) iRufType
        end select
         
        ! DGM 10/2010 Some new roughness types want an additional
        ! value to go along. It should be on the same line as 
        ! the roughness type as an additional number. Read it.
        if( iRufType == eMinGrad .or. iRufType == eTwoModel ) then
!            i = index( sValue, ',' )
!            if( i == 0 ) then
!                i = index( sValue, ' ' )
!            endif
            if( i .ne. 0 ) then
                read(sValue(i+1:),*,err=123,end=123) rRufParam
            endif
        endif
123     continue
        
    case ('diagonal penalties')
      ! removing this from Occam.f90, but leaving here so case
      ! so old files don't crash on read in
      write(*,*) ' WARNING: diagonal penalties not supported here...'
      write(*,*) ''
      !  if( idcode( sValue, ioUnitOccamIterFile ) == 0 ) then
      !      gbAddDiags  = .false.
      !  else
      !      gbAddDiags  = .true.
      !  endif
        
    case ('debug level')
        read(sValue,*) idebug
        
    case ('iteration')
        ! Actually the # of the previous iteration - so that at startup an iter00.iter
        !   file will be written reflecting the startup values.
        read(sValue,*)  nCurrentIter 
        
    case ('pmu', 'lagrange value')
        read(sValue,*) pmu 
        
    case ('rlast', 'roughness value')
        read(sValue,*) rlast
        
    case ('tlast', 'misfit value')
        read(sValue,*) tobt 
        
    case ('ifftol', 'misfit reached')
        read(sValue,*) ifftol
        
    case ('stepsize cut count')
        read(sValue,*) nFracTimes 
        
    case ('delta misfit limit')
    ! DGM Sept 2009 - if misfit changes by < this amount, then "auto converge"
    !   and STOP the inversion without any more iterations.  This allows the 
    !   user to not waste iterations dwelling in a high spot and make a decision
    !   about bad data points (based on residuals) and/or minimum error size.
        read(sValue,*) rDeltaMisfit
    case ('bounds transform')
         cBoundsTransform = trim(sValue)  
         !
         ! Check to make sure input is valid:
         !
         select case (cBoundsTransform)
         case ('exponential','bandpass')
             write(*,*) '  Using model bounds transform:  ', trim(cBoundsTransform) 
         case default 
             write(*,*) '  Error, unrecognized model bounds transform:  ', trim(cBoundsTransform)
             stop
         end select
       
         
         
    case ('model limits','model bounds')   ! KWK March 2011: Now using true bounds transformations...
        ! Limit value string must have two values: min, max
        ! Comma is required.
        i = index( sValue, ',' )
        if (i == 0) goto 196
        read(sValue(:i-1),*) lowerBound
        
        ! Skip space padding
        do i=i+1,len_trim(sValue)
            if (sValue(i:i) .ne. ' ') exit
        enddo
         read(sValue(i:),*) upperBound
        
        ! Min & Max must be reasonable
        if (upperBound < lowerBound .or. upperBound == lowerBound) then
            goto 196
        endif
        
        lConstrainParams = .true.
        write(*,*) ' Using nonlinear transformation to bound parameters:'
        write(*,'(5x,g9.2,a7,g9.2)')  lowerBound, '< m <', upperBound
        
    case ('model value steps')   ! DGM Nov 2006 see module  Occam_internal_mod
        read(sValue,*) gnModelSteps 
        if (gnModelSteps > 0) then
            gbModelStepped  = .true.
            write(*,*) '  Limiting model values to steps of', gnModelSteps
        endif
        
    case ('modify roughness')
        ! Left-over from previous investigation of ways to make 
        ! inversion sharper.  Any startup file with this tag in
        ! it will no longer function properly.  DGM Dec 2006.
        write(*,*) '  Statement: Modify Roughness :no longer supported'
        
        
    case ('no. parms', 'param count')
        read(sValue,*) nParams 
        exit    ! The next thing in the file is the data param list!
        
    case default
        write(*,*) 'Error reading startup/iteration file!'
        write(*,*) 'Unknown or unsupported code:', sCode
        stop
    end select
    enddo
    
    ! Did we end early without a params statement??
    if (nParams <= 0) then
        write(*,*) 'Cannot have zero parameters.'
        write(*,*) 'Is header of startup file formatted properly?'
        stop
    endif
    
    ! write out some values
    write(*,*) ' '
    write(cTemp,*) nParams
    write(*,'(a24,a)') ' Number of Parameters:  ', trim(adjustl(cTemp))
    write(cTemp,*) tolreq
    write(*,'(a24,a)') ' Target Misfit:  ', trim(adjustl(cTemp))
    write(cTemp,*) ifftol
    write(*,'(a24,a)') ' Converged Flag:  ', trim(adjustl(cTemp))
    write(*,'(a24,I1,", ",F8.4)') ' Roughness Type: ', iRufType, rRufParam
    
    ! DGM Oct 2006 - allocate what we need for the rest of this routine.
    ! Most vars allocated in the main program.
    ! If we don't survive allocation, explain nicely & stop the program.
    call allocatenParamArrays
    
    ! Read the model parameters
    read (unit=ioUnitOccamIterFile,fmt= *, end=198, err=199) (pm(i), i=1,nParams)

    close(ioUnitOccamIterFile)
    
    ! Apply non-linear transformation for model bounds if requested:
    if (lConstrainParams) then
               
       ! Make sure pm is inside the bounds (nudged just inside bounds if not):
       pmtol = (upperBound - lowerBound) / 232
        do i = 1,nParams
         if (pm(i) <= lowerBound + pmtol) pm(i) = lowerBound + pmtol
         if (pm(i) >= upperBound - pmtol) pm(i) = upperBound - pmtol
        enddo
        ! now transform it:
        pm = transformToUnbound(pm) 
    endif
    return
    ! Various error conditions
    196 continue
    write(*,*) 'Invalid "Model Bounds" entry in startup file.'
    write(*,*) 'The format for this line is'
    write(*,*) '   Model Bounds: <min>,<max>'
    write(*,*) 'Do not include angle brackets.  Comma required.'
    write(*,*) 'Min must be less than max.'
    close(ioUnitOccamIterFile)
    stop
    
    198 call filerr (' Startup file ended prematurely', ioUnitOccamIterFile)
    199 call filerr (' Error reading startup file', ioUnitOccamIterFile)
    
    end subroutine readIteration


    !-----------------------------------------------------------------------
    subroutine writeIteration (descr,maxitr,tolreq,pmu,rlast,tobt,ifftol)      
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    !
    ! writes an iteration output file for OCCAM.
    ! calls RespOut
    !
    ! Revisions:
    !  3.1, Sept 2008       KWK: applied the clunkly F95 character conversion + adjustl 
    !                       commands to left justify the output numbers for neatness.
    !
    !  3.0, November 2006   DGM updated to new FLEX format & to add roughness setup.
    !  2.0, 13 May 1992
    !
    
    ! arguments:
    integer maxitr, ifftol
    real(RealPrec) tolreq, pmu, rlast, tobt
    character(50)  descr, cTemp
    character(80) datetm
    ! local variables
    integer lerr, i
    !   lerr = error flag on file operation
    !   i = loop index
    character(18) cNum
    
    !-----------------------------------------------------------------------
    ! create file name suffix
    write (cNum,*) nCurrentIter
    !write (cNum,'(I4.3)') nCurrentIter      ! DGM put leading zeros so files sort properly
    
    !if (nCurrentIter <= 9) cNum(1:1) = '0'
    
    ! open iteration file
    open (unit=ioUnitOccamIterFile, file=trim(cRootName)//'_'//trim(adjustl(cNum))//'.iter', iostat=lerr)
    
    if (lerr .ne. 0) then
        ! we could be in the middle of an iteration, so don't die on error, but
        ! send a message and forge on regardless
        write(*,*) ' Error opening iteration file'
        return
    end if
    
    ! Get the current date & time
    call datime(datetm)
    
    ! write the header
    write(ioUnitOccamIterFile,'(a)')       'Format:             OCCAMITER_FLEX'
    write(ioUnitOccamIterFile,'(a,a)')     'Description:        ', trim(descr)
    write(ioUnitOccamIterFile,'(a,a)')     'Model File:         ', trim(modelFileName)
    write(ioUnitOccamIterFile,'(a,a)')     'Data File:          ', trim(dataFileName)
    write(ioUnitOccamIterFile,'(a,a)')     'Date/Time:          ', trim(datetm)
    write(cTemp,*) maxitr
    write(ioUnitOccamIterFile,'(a,a)')     'Iterations to run:  ', trim(adjustl(cTemp))
    write(cTemp,*) tolreq
    write(ioUnitOccamIterFile,'(a,a)')     'Target Misfit:      ', trim(adjustl(cTemp))
    if( rRufParam .ne. 0.0 ) then
        write(cTemp,*) iRufType, ' ', rRufParam
    else
        write(cTemp,*) iRufType
    endif
    write(ioUnitOccamIterFile,'(a,a)')     'Roughness Type:     ', trim(adjustl(cTemp))
    write(cTemp,*) nFracTimes
    write(ioUnitOccamIterFile,'(a,a)')     'Stepsize Cut Count: ', trim(adjustl(cTemp))
    if (lConstrainParams) then
        write(ioUnitOccamIterFile,'(a,a)') &
                            'Bounds Transform:       ', trim(cBoundsTransform) 
        write(ioUnitOccamIterFile,'(a,f10.6,",",f10.6)') &
                            'Model Bounds:       ', lowerBound, upperBound
    else
        write(ioUnitOccamIterFile,'(a)')   '!Bounds Transform:       bandpass'
        write(ioUnitOccamIterFile,'(a)')   '!Model bounds:       min,max'
    endif
    if (gbModelStepped) then
        write(ioUnitOccamIterFile,'(a,f10.6)') &
                            'Model Value Steps:  ', gnModelSteps
    else
        write(ioUnitOccamIterFile,'(a)')   '!Model Value Steps:  stepsize (e.g. 0.2 or 0.1)'
    endif
    write(cTemp,*) rDeltaMisfit
    write(ioUnitOccamIterFile,'(a,a)')    'Delta Misfit Limit: ', trim(adjustl(cTemp))
    
    write(cTemp,*) idebug
    write(ioUnitOccamIterFile,'(a,a)')    'Debug Level:        ', trim(adjustl(cTemp))
    write(cTemp,*) nCurrentIter
    write(ioUnitOccamIterFile,'(a,a)')    'Iteration:          ', trim(adjustl(cTemp))
    write(cTemp,*) pmu
    write(ioUnitOccamIterFile,'(a,a)')    'Lagrange Value:     ', trim(adjustl(cTemp))
    write(cTemp,*) rlast
    write(ioUnitOccamIterFile,'(a,a)')    'Roughness Value:    ', trim(adjustl(cTemp))
    write(cTemp,*) tobt
    write(ioUnitOccamIterFile,'(a,a)')    'Misfit Value:       ', trim(adjustl(cTemp))
    write(cTemp,*) ifftol
    write(ioUnitOccamIterFile,'(a,a)')    'Misfit Reached:     ', trim(adjustl(cTemp))
    
    ! write the end of the header and all of the param data.
    write(cTemp,*) nParams
    write(ioUnitOccamIterFile,'(a,a)')    'Param Count:        ', trim(adjustl(cTemp)) 
    
    if (lConstrainParams) then
       write(ioUnitOccamIterFile,'(4(2x,g15.7))') (transformToBound(pm(i)), i=1,nParams)
    else
       write(ioUnitOccamIterFile,'(4(2x,g15.7))') (pm(i), i=1,nParams)
    endif

    close(ioUnitOccamIterFile)
       
    end subroutine writeIteration

    !-----------------------------------------------------------------------
    subroutine parseCode( nLen, sLine, sCode, sValue, bComment )
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine Revision 3.0, November 2006
    ! DGM Nov 2006 - parse a line read from a file into a code & value.
    ! Force the code to be all lowercase with no ending colon.  Terminate
    ! the line at a '%' or '!' sign (these allow for user comments!)

    ! Args
    integer, intent(in)   :: nLen
    character(nLen)       :: sLine
    character(nLen), intent(out) :: sCode, sValue
    logical, intent(out)    :: bComment
    
    ! Local vars
    integer :: iFrom, iTo
    
    ! Init returns
    bComment = .false.
    sCode = ' '
    sValue = ' '
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    ! Skip any beginning blanks
    do iFrom = 1,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    ! If the first char is a comment char, then the whole line is a comment.
    ! DGM April 2008 Also, if the line is blank, consider it a comment.
    if (iFrom >= nLen) then !KWK may 2009 pulled this out in from since sometimes iFrom > nLen and this kills (iFrom:iFrom) below
        bComment = .true.
        return
    endif
    
    if(  sLine(iFrom:iFrom) == '%' &
        .or. sLine(iFrom:iFrom) == '!' ) then
        bComment = .true.
        return
    endif
    
    ! Pull off the code value. Cvt to lowercase as we go.
    iTo = index(sLine,':') - 1
    if (iTo < iFrom) then
        write(*,*) 'Parsing Error: missing colon in line below:'
        write(*,*) sLine
        return
    endif
    sCode = sLine(iFrom:iTo)
    call Lower(sCode)
    
    ! Skip spaces after the colon
    do iFrom = iTo+2,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! Get the rest, up to any comment
    sValue = sLine(iFrom:)
    iTo = len_trim(sValue)
    
    iFrom = index(sValue,'%')
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    iFrom = index(sValue,'!')
    if (iFrom > 0 .and. iFrom < iTo) then
        sValue(iFrom:iTo) = ' '
    endif
    !call Lower(sValue)   ! No: Some values are filenames which are case-sensitive on UNIX!
    
    end subroutine parseCode
    
    !==============================================================================!  
    !==================================================================! parseLine !  
    !==============================================================================! 
    subroutine parseLine( nLen, sLine, bComment )
    !   
    ! Subroutine to check if the sLine is blank or a comment line, and if it isn't
    ! then any comment at the end of the line is blanked. 
    !
    ! This is styled after D. Myer's parseCode function.
    !
    ! This is useful for reading in data tables that have user comment lines 
    ! or comments at the end of the line, denoted by the ! and % symbols.
    ! 
    ! If the entire line is a comment, bComment = .true.
    !
    ! Kerry Key
    ! Scripps Institution of Oceanography
    ! kkey@ucsd.edu
    !
    ! Version 1.0.   April, 2008.
    !
 
    
    ! Args
    integer, intent(in)     :: nLen
    character(nLen)         :: sLine
    logical, intent(out)    :: bComment
    
    ! Local vars
    integer :: iFrom, iTo
    
    ! Init returns
    bComment = .false.
    
    
    ! Convert all tab characters to spaces
    forall( iTo = 1:nLen, ichar(sLine(iTo:iTo)) == 9 ) sLine(iTo:iTo) = ' '
    
    ! Skip any beginning blanks
    do iFrom = 1,nLen
        if (sLine(iFrom:iFrom) .ne. ' ') exit
    enddo
    
    ! If the first char is a comment char, then the whole line is a comment.
    ! DGM April 2008 Also, if the line is blank, consider it a comment.
    if( iFrom >= nLen .or. sLine(iFrom:iFrom) == '%' &
        .or. sLine(iFrom:iFrom) == '!' ) then
        bComment = .true.
        return
    endif
    
    ! Now trim off any comments at the end of the line  
    iTo = len_trim(sLine)
    iFrom = index(sLine,'%')
    if (iFrom > 0 .and. iFrom < iTo) then
        sLine(iFrom:iTo) = ' '
    endif
    iFrom = index(sLine,'!')
    if (iFrom > 0 .and. iFrom < iTo) then
        sLine(iFrom:iTo) = ' '
    endif
    
    end subroutine parseLine

    !-----------------------------------------------------------------------
    subroutine Lower( s )
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine Revision 3.0, November 2006
    ! DGM Nov 2006 - convert string to lower case
     
    character(*), intent(out)  :: s
    integer i
    do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
        s(i:i) = char(ichar(s(i:i)) + 32)
      endif
    enddo
    
    end subroutine Lower


    !-----------------------------------------------------------------------
    subroutine filerr(mssg, io1)
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine Revision 3.0, November 2006
    ! Subroutine Revision 2.00, 13 May 1992
    !
    ! filerr prints an error message, closes file io1, and stops
    !
    ! on input:
    !    mssg = character string containing error message
    !    i01 = unit number of file to close (0 = no file open)
    ! on output:
    !    outputs nothing
    ! calls:
    !    no other routines
    !
    !
    
    ! input arguments:
      character(*) mssg
      integer io1
    !
    !-----------------------------------------------------------------------
      write(*,*) mssg
      if (io1 > 0) close (io1)
      stop
    end subroutine filerr


    !-----------------------------------------------------------------------
    subroutine datime(datetm)
    !-----------------------------------------------------------------------
    ! OCCAM 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    !
    ! Revision  Date            Comments:
    ! 4.0       April 2008      KWK: entirely rewritten using intrinsic Fortran90 time call
    !                           eliminating the heinous compiler dependent timing functions.
    ! 3.0       November 2006
    ! 1.0       18 May 1992
    !
    ! Date and time utility for OCCAM 2.0.  
    !
    ! on input:
    !   nothing
    ! on output:
    !   datetm = 80 character string containing data and time
    !
    character(80) datetm
    character(10) cdate,ctime
    
    !
    ! Standard intrinsic Fortran90 data_and_time call:
    !
    call date_and_time(cdate,ctime) !this ouputs only values, ignoring other optional arguments
    
    datetm = cdate(5:6)//'/'//cdate(7:8)//'/'//cdate(1:4)//' '// ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:)               
    
    end subroutine datime
    
    !==============================================================================!  
    !================================================! getOccamCommandLineArguments!  
    !==============================================================================!    
    subroutine getOccamCommandLineArguments(cStartup,cRootName,bFwdOnly)
    !
    ! Subroutine to get command line arguments for Occam.
    ! 
    ! Uses the Fortran2003 intrinsics for command line stuff:
    !
    ! command_argument_count()
    ! get_command_argument()
    !
    ! Kerry Key
    ! Scripps Institution of Oceanography
    ! kkey@ucsd.edu
    !
    !
    ! December 13, 2008     Created.
    !
    
    character(50), intent(inout) :: cStartup,cRootName
    logical, intent(inout)       :: bFwdOnly    
    
    character(50) :: arg
    integer :: n
    
    n = command_argument_count()
    
    if (n > 0 ) then  
    
       call get_command_argument(1, arg)
       
       select case (arg)
       
            case ('?')
           !
           ! Help info requested:
           !
            write(*,*) ' '
            write(*,*) 'Occam accepts one flag:'
            write(*,*) '      -F    When given, this flag means produce a forward model only.'
            write(*,*) '            The forward model will have the name <output root name>.fwd.'
            write(*,*) '            If no output root is given, the name is Forward.fwd'
            write(*,*) ' '
            write(*,*) 'Occam has 2 optional parameters:'
            write(*,*) '      <startup file> -- name of the file with the starting model'
            write(*,*) '                    and runtime parameters.  Traditionally this'
            write(*,*) '                    file has been called ''startup''.  To continue'
            write(*,*) '                    a previous inversion, specify an iteration file.'
            write(*,*) '                    If not given, ''startup'' is assumed.'
            write(*,*) '      <output root name> -- This parameter can only be given if '
            write(*,*) '                    a startup filename is given.  It is the root'
            write(*,*) '                    name to use for the output iteration, response,'
            write(*,*) '                    and log files.  If not given, these files will'
            write(*,*) '                    be named ITERxx.iter, RESPxx.resp, LogFile.logfile.'
            write(*,*) '                    If given, the root name will be used in place of'
            write(*,*) '                    ''ITER'', ''RESP'', and ''LogFile'' in the above names.'
            write(*,*) '                    Do NOT include a path or a ''.''.'
            stop
            
            case ('-f','-F','/F','/f',char(92)//'F',char(92)//'f') ! Forward run only
                bFwdOnly = .true.
                
                if ( n  > 1 )  call get_command_argument(2, cStartup)
                if ( n  == 3 ) call get_command_argument(3, cRootName)
                
            case default
                cStartup = arg
                if ( n == 2 ) call get_command_argument(2, cRootName)
                
       end select
       
    endif 
    end subroutine  getOccamCommandLineArguments

!-----------------------------------------------------------------------
    subroutine intializeOccamParamLimits
        lConstrainParams = .false.
        cBoundsTransform = 'bandpass'
        gbModelStepped  = .false.
        gnModelSteps    = 0
    end subroutine intializeOccamParamLimits
    
    !-----------------------------------------------------------------------
    subroutine allocatenParamArrays
    !-----------------------------------------------------------------------
    implicit none
    integer(4)   :: nAllocErr
    
    allocate( pm(nParams)         &
        , ptp(nParams,nParams)    &
        , wjtwj(nParams,nParams)  &
        , wjtwd(nParams)          &
        , pwk1(nParams)           &
        , pwk2(nParams)           &
        , pwk3(nParams)           &
        , prewts(nParams)         &
        , premod(nParams)         &
        , pwk4(nParams)           &   ! work var or TOFMU
        , stat=nAllocErr )   
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif
    
    
    select case (linearSolver)
    
    case ('parker_cholesky')
           allocate( aMat(nParams,nParams+1),stat=nAllocErr  ) ! extra dim for use in cholin() 
    case default 
       allocate( aMat(nParams,nParams), stat=nAllocErr  )
    end select
    
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif

   
    end subroutine allocatenParamArrays
    
    !-----------------------------------------------------------------------
    subroutine allocateNdArrays
    !-----------------------------------------------------------------------
   
    integer  :: nAllocErr
    
    allocate( dm(nd)               &
            , dwk1(nd)             &
            , dwk3(nd)             &
            , dwk4(nd)             &
            , stat=nAllocErr )
    if (nAllocErr .ne. 0) then
        write(*,*) 'Out of memory.  Try reducing the model size'
        stop
    endif
    
    end subroutine allocateNdArrays

    !-----------------------------------------------------------------------
    subroutine allocateAssocArrays
    !-----------------------------------------------------------------------    
    
    integer  :: nAllocErr
    
    if (npm_assoc > 0 ) then
        allocate ( dwk1_assoc(npm_assoc) &
                 , dwk2_assoc(npm_assoc) &
                 , dwk3_assoc(npm_assoc) & 
                 , stat=nAllocErr )  ! KWK May 2008
        if (nAllocErr .ne. 0) then
            write(*,*) 'Out of memory.  Try reducing the model size and number of associated parameters'
            stop
        endif    
    endif
    
    end subroutine allocateAssocArrays
    
    !-----------------------------------------------------------------------
    subroutine deallocateOccam
    !-----------------------------------------------------------------------
    
    integer :: n
    
    !! NB: the variable wj is deallocated after each use to conserve memory.
    
    deallocate( pm, ptp, wjtwj, wjtwd, pwk1, pwk2, pwk3, prewts, premod &
              , aMat, pwk4, d, sd, dp, dm, dwk1, dwk3, dwk4 &
              , stat=n )        ! stat=  doesn't allow process to blow chunks on dealloc error
   
    if (npm_assoc > 0 ) then
        deallocate ( pm_assoc ,dwk1_assoc, dwk2_assoc, dwk3_assoc , stat=n ) 
    endif
    
    ! Tell the forward module (whomever it is) to deallocate its memory
    call deallocateFwd
    
    
    end subroutine deallocateOccam
    
    !-----------------------------------------------------------------------
    subroutine compOccamIteration(tolreq,tobt,ifftol,pmu,rlast,stepsz,konv)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla, CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.01, 20 Jan 1993
    
    ! Occam executes one iteration of a smooth model finder
    
    !  Arguments:
    !   
    !     tolreq = required rms misfit
    !     Tobt = rms misfit from last iteration (large initially)
    !     ifftol = record of whether feasible model has been attained (0=no,
    !         1=Yes).  Initially set 0 and then left alone.
    ! KWK now in module:    nit = the number of previous iterations in this inversion. Set to 0 on first 
    !         call, it will be updated by occam. 
    !     pmu = log10(lagrange) multiplier to start searches with.  Initially 
    !         Set large and then left from last call for efficiency
    !     Rlast = roughness of last model.  Initially set large and then left.        
    
    ! On output:
    
    !  Occam_interface_mod  contains
    !     pm(nParams) = vector of updated model parameters
    !     Dm(ndata+) = vector containing *approximate* response of new model
    !  Arguments:
    !     Tobt = rms misfit of new model response.
    !     Nit = number of previous calls to occam1 during this inversion
    !     Stepsz = rms measure of the changes in the model parameters 
    !     Konv = status flag:
    !        0 = Normal exit for a successful iteration,
    !        1 = A perfectly smooth model has been found for the required tolerance
    !        2 = Convergence problems.  Unable to find a model with an R.M.S.. Misfit
    !          Less than or equal to that of the last iteration's
    !        3 = Convergence problems.  Unable to find a model with an R.M.S.. Misfit
    !          Of tolreq and smoother than the last iteration's
    
    ! Subroutines required and supplied:
    
    !   Makptp, a subroutine which creates the roughening matrix times its trans.
    !   Tofmu(amu), a function which returns the rms misfit of the response
    !      Of the model constructed using the lagrange multiplier amu.
    !      (calls cholin,cholsl,computefwd,anorm)
    !   Trmult, mult, atamult and anorm, subroutines for matrix operations
    !   Fminocc, a subroutine which minimises a univariate function
    !   Froot, a subroutine which finds the root of a univariate function
    !   Minbrk, a subroutine which brackets a univariate minimum
    !   Scanmu, a subroutine which generates the misfit table when ibug > 0
    
    ! Subroutines which must be supplied by the user:
    
    !  Computefwd: computes the forward function for model pm() at
    !    The data parameters dp() and returns it in dm().  Optionally
    !    Returns the matrix of partials, wj(nd,nParams).
    !  constructPenaltyMat(), assembles the penalty matrix (this will be dependent on model 
    !    Type and dimension)
    
    
    ! Arguments:
          integer ifftol,konv
          real(RealPrec) tolreq,tobt,pmu,rlast,stepsz
    ! Local variables:
          real(RealPrec) tol0, tmin, tint
          real(RealPrec) amu1, amu2, amu3, tamu1, tamu2, tamu3
          real(RealPrec) ruf, pmu2
 
    
    ! Local parameters
          real(RealPrec) tolm, toli
    ! Tolm and toli are the tolerances for fminocc and froot respectively
    ! Decreasing them will improve accuracy at the cost of more forward 
    ! Calculations
          parameter (tolm= 0.1, toli = 0.001)
    
    !-----------------------------------------------------------------------
          if (abs(idebug) >= 1) write(*,*) ' Entering Occam iteration...'
    ! Set things up for this iteration:
          konv = 0
    ! Frac controls the step size; normally will remain at 1.0 Unless we have
    ! Convergence problems
          frac = 1.0
    ! Construct wjtwj, wjtwd, and find initial misfit
          if (abs(idebug) >= 1) write(*,*) &
         &  ' Constructing derivative dependent matrices...'
          call makjtj(tol0)
    
          write(ioUnitOccamLogFile,*) 'Starting R.M.S. = ',tol0
          if (nCurrentIter == 0) then
              if ((tol0 <= tolreq) .and. (ifftol == 0)) then
                if (abs(idebug) >= 1) &
                    write(*,*) ' Tolerance met.  This iteration begins smoothing.'
                ifftol = 1
              end if
          end if
          if (abs(idebug) >= 1) then
             write(*,fmt='(a,1x,g12.5)') '  Staring misfit:', tol0 !, ' Mu:', 10**pmu
          endif
          write(ioUnitOccamLogFile,*) '** Iteration ',nCurrentIter+1,' **'
    
    ! Construct the penalty matrix multiplied into its transpose
          if (abs(idebug) >= 1) write(*,*)' Constructing penalty matrix...'
          call makptp() 
    
    
    ! The next block of code controls the selection of the lagrange multiplier,
    ! pmu:
          if (abs(idebug) >= 1)  write(*,*) ' Searching Lagrange multiplier...'
         
    120   nfor = 0
    ! Produce the misfit function if required (this will be used to bracket min)
          if (abs(idebug) >= 1) write(*,*) ' Bracketing minimum...'
          if (abs(idebug) >= 2) then
            call scanmu(amu1,amu2,amu3,tamu2,tolreq)
          else
    ! Bracket the minimum using minbrk and two guesses
            amu1 = pmu - 1.0
            amu2 = pmu
            
            ! In: minbrk( mu guess #1, #2, ..., Functional to calc rms)
            ! Out: amu1,3 - brackets around amu2 - the minimum mu
            !      Tamux = tofmu( amux ) -- the value of the rms misfit
            call minbrk(amu1,amu2,amu3,tamu1,tamu2,tamu3,tofmu)
          end if
          
    ! Find the minimum
          if (tamu2 < tolreq) then
    ! We've been lucky and found an acceptable minimum using minbrk
            tmin = tamu2
            pmu = amu2
            write(ioUnitOccamLogFile,*) 'Minimum tol from minbrk is at mu =',pmu
          else
            if (abs(idebug) >= 1) write(*,*) ' Finding minimum...'
            
            tmin = fminocc(amu1,amu2,amu3,tamu2,tofmu,tolm,pmu)
            
            write(ioUnitOccamLogFile,*) 'Minimum tol from fminocc is at mu =',pmu
          end if
          write(ioUnitOccamLogFile,*) 'and is =',tmin
          write(ioUnitOccamLogFile,*) 'using ',nfor,' evaluations of function'
          
    ! If the new minimum tolerance is greater than the tolerance from the
    ! Previous model, we are having convergence problems. Cut the step size
    ! tol0 is tol of the previous model.
    ! tolreq is the target tol.
          if( (  ((tmin >= tol0)  .and. (ifftol == 0))    & ! Misfit greater than starting value
            .Or. ((tmin > tolreq) .and. (ifftol == 1))    & ! Starting model reached target, but current model doesn't
              ) .and. idebug >= 0 ) then
            !   .Or. ( (Abs(tmin - tol0) < .01) .and. (ifftol == 0)  )    &  ! Mistfit nearly the same as starting value
    
            write(ioUnitOccamLogFile,*) 'Divergence problems, cutting step size'
            if (abs(idebug) >= 1)  write(*,*) ' Divergence problems, cutting stepsize..'
            Frac = frac*2.0
            
            ! We have cut step size a lot to no effect: give up
            if (frac >= 2**nfractimes) then   
              konv = 2
              nCurrentIter = nCurrentIter + 1     ! Dgm nov 2006 - must increment or previous iter file will be overwritten
              return
            else
              goto 120
            end if
          end if
    
          if (tmin < tolreq) then
            if (abs(idebug) >= 1) write(*,*) ' Finding intercept...'
    ! Tolerance is below that required; find intercept.
    ! The lower value of mu bracketing the intercept is easily found: it is
    ! Just the minimum
            amu1 = pmu
            tamu1 = tmin
    ! The upper bound is found by testing successively greater mu
            amu2 = pmu
            tamu2 = tmin
            nfor = 0
            do while (tamu2 < tolreq)
    ! Successively double Mu:
              amu2 = amu2 + 0.30103 ! Log10(2x) = log10(x) + (0.30103)
              Tamu2 = tofmu(amu2)
            enddo
            
            pmu2 = froot(tofmu,amu1,amu2,tamu1,tamu2,tolreq,toli)
            
            tint = tamu2 + tolreq
            write(ioUnitOccamLogFile,*) 'Intercept is at mu = ',pmu2
            write(ioUnitOccamLogFile,*) 'and is = ',tint
            write(ioUnitOccamLogFile,*) 'using ',nfor,' function evaluations'
            
            if (abs(idebug) >= 1) write(*,'(a,g16.8)') '  Optimal Mu: ', 10**pmu2
          else
    ! Tolerance is not yet small enough.  We will keep the minimum.
    ! Since fminocc returns with pwk1() instead of pwk2() we need this copy
            pwk2 = pwk1     ! Array math
            dm = dwk1       ! Array math
            if (npm_assoc > 0) then
                pm_assoc  = dwk1_assoc
            endif
            if (abs(idebug) >= 1) write(*,'(a,g12.5)') '  Optimal Mu: ', 10**pmu
          end if
    !***End lagrange multiplier selection
    
    ! Compute roughness.  We do this by a function call to avoid having the penalty
    ! Matrix hang around taking up memory.
          if (abs(idebug) >= 1) write(*,*) ' Computing roughness...'
          Ruf = fndruf()
    
    ! If we attained the intercept last iteration but the model is getting 
    ! Rougher we have problems with convergence.
! DGM 11/2010 Occam has had a problem for some years where it will reach the 
! target RMS but not converge and the roughness will slowly creep upwards.
! Taking out the 1.01 factor to force the roughness to always decline.
!          If ((ruf>1.01*Rlast) .and. (ifftol==1) .and. (idebug>=0)) then
          if ((ruf >= Rlast) .and. (ifftol==1) .and. (idebug>=0)) then
            write(ioUnitOccamLogFile,*) 'Roughness problems, cutting step size'
            frac = frac*2.0
    ! Check to see if all is hopeless
            if (frac >= 2**nfractimes) then   ! > 1.0E+05) then
                if (abs(idebug) >= 1) &
         &      write(*,*) 'Roughness problem not resolved by small steps'
                konv = 3
                nCurrentIter = nCurrentIter + 1     ! Dgm nov 2006 - must increment or previous iter file will be overwritten
                return
            end if
    ! Otherwise plow on
            if (abs(idebug) >= 1) then
                write(*,*) 'Roughness growing! Cut step size & recalc mu'
                write(*,*) '     roughness =', ruf
            endif
            goto 120
          end if
    
    ! Save new model and compute step size
          if (abs(idebug) >= 1) write(*,*) ' Computing model stepsize...'
          pwk3 = pwk2 - pm      ! Array math -- may have memory allocation drawback!
          pm = pwk2             ! Array math
          
    ! The stepsize is the actual change in the model, normalized by nParams:
          stepsz = sqrt(anorm(nParams,pwk3)/nParams)
          write(ioUnitOccamLogFile,*) 'Stepsize is = ',stepsz
          write(ioUnitOccamLogFile,*) 'Roughness is = ',ruf
          write(ioUnitOccamLogFile,*) ' '
    ! Tidy up
          nCurrentIter = nCurrentIter + 1
          rlast = ruf
          if (tmin < tolreq) then
            if (abs(idebug) >= 1 .and. ifftol == 0) &
         &      write(*,*) ' Tolerance met.  Next iteration begins smoothing.'
            ifftol = 1
            tobt = tint
          else
            tobt = tmin
            
            ! DGM Sept 2009 - not converged normally.  Check the auto-converge
            ! setting and set a konv code if necessary.
            ! NB: If no auto-converge specified, then rDeltaMisfit == 0.0
            if (abs(tobt - tol0) < rDeltaMisfit) then
                konv = 4
            endif
          end if
    ! See if we have a perfectly smooth model that fits data:
          if (ruf < 1.0E-5 .and. Tobt <= 1.01*tolreq) konv = 1
          if (abs(idebug) >= 1) write(*,*) ' Leaving Occam iteration'
          write(*,*) ''
          return
    end subroutine compOccamIteration

    
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! private subroutines:
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
    
 
    
    !-----------------------------------------------------------------------
    subroutine scanmu(am1,am2,am3,t2,t0)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    
    ! Scanmu produces a report of misfit versus lagrange multiplier.
    ! It calls the forward routine a lot so is really only used for debugging
    ! And particularly nasty cases
    ! It looks for the global minimum as follows:
    !   Starting at the largest mu used for the debug report, we look for the 
    !   First minimum.  We also look for the glolbal min.  We want the global 
    !   Min unless the local min at larger mu is lower than the required 
    !   Tolerance. 
    
    ! On input:
    !   T0 is the required tolerance
    !   /Result/ idebug is the print level from occam1
    !   /Result/ ioUnitOccamLogFile is a unit number for output
    ! On output:
    !   Am1, am2, am3 are three lagrange multipliers which bracket a minimum
    !   T2 is the tolerance at am2.
    !   /Result/ ..pmw1().. Has the model associated with am2,t2.
    ! Subroutines called: 
    !   Tofmu

    
    ! Arguments:
          real(RealPrec) am1,am2,am3,t2,t0
    ! Local variables
          logical search
          integer j,k
          real(RealPrec) t1,tg2,t3,tl2,aml3,aml2,aml1,amg3,amg2,amg1,amu,t

          write(ioUnitOccamLogFile,*) 'misfit as a function of Mu:'
          t1 = 1.0E+09
          t2 = 1.1E+09
          tg2 = 0.9E+09
          am2 = 16.
          Am1 = 16.
          Search = .True.
    
    ! Start at the top (smoothest and presumably worst fit), sweeping from
    ! Log10(mu) = 16 to log10(mu) = -4 at two values per decade
          do 100 k = 16,-4,-1
    ! Keep the last two values for minimum search
            t3 = t2
            t2 = t1
            am3 = am2
            am2 = am1
    ! Keep the model associated with the middle value if we are searching
    ! For the minimum
            if (search) then
                pwk1 = pwk2     ! Array math
                dwk1 = dm
                if (npm_assoc > 0) then  ! Pass associated parameters
                   dwk1_assoc = pm_assoc 
                endif              
            end if
    ! Keep the model associated with the middle value anyway in case
    ! Global and local minima are found
            pwk3 = pwk2     ! Array math
            dwk3 = dm
            if (npm_assoc > 0) then  ! Pass associated parameters
               dwk3_assoc = pm_assoc 
            endif          
    ! Get the next mu and associated tolerance
            am1 = float(k)/2.
            T1 = tofmu(am1)
    
            if ((t2 <= T1) .and. (T2 <= T3) .and. Search) then
    ! We have found the first minimum. Keep it and turn off search
    ! (Model is already stashed in pm1())
              search = .False.
              Tl2 = t2
              aml3 = am3
              aml2 = am2
              aml1 = am1
            end if
    ! Keep an eye out for the global min in case it is different
            if (t2 <= Tg2) then
              tg2 = t2
              amg3 = am3
              amg2 = am2
              amg1 = am1
    ! Copy across model associated with global min.
              pwk1 = pwk3       ! Array math
              dwk1 = dwk3
              if (npm_assoc > 0) then  ! Pass associated parameters
                 dwk1_assoc = dwk3_assoc
              endif            
            end if
    ! Record in the log file
    !kwk debug:: not needed since this is written in tofmu:
    ! write(ioUnitOccamLogFile,*) am1,t1
    
    ! Output the model as well if required
            if (abs(idebug) >= 3) then
              do j = 1,nParams
                write(ioUnitOccamLogFile,*) '   ',pwk2(j)
              enddo
            end if
    100   continue
    ! write out the gauss step at highest debug level
          if (abs(idebug) >= 3) then
            amu = -99.0
            T = tofmu(amu)
            write(ioUnitOccamLogFile,*) am1,t1
            do j = 1,nParams
              write(ioUnitOccamLogFile,*) '   ',pwk2(j)
            enddo
          end if
    
    ! Now check out the minima
    !  If the first minimum achieves the tolerance required we want it
          if (tl2 <= T0) then
            t2 = tl2
            am1 = aml1
            am2 = aml2
            am3 = aml3
          else 
    ! Otherwise the global min is just fine
            t2 = tg2
            am1 = amg1
            am2 = amg2
            am3 = amg3
          end if
          return
    end subroutine scanmu
    
    !-----------------------------------------------------------------------
    subroutine makjtj(tol0)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    
    ! Constructs matrices that depend on the jacobian (wjtwj and wjtwd).
    ! The jacobian only needs to be around long enough to do
    ! This, so by stuffing this into a subroutine dynamic memory allocation
    ! Will make wj disappear after use.
    
    ! On input:
    !   Current model and data stored in common /model/ and /data/ 
    ! On output:
    !   /Result/ wjtwj = w.J.Trans(w.J)
    !   /Result/ wjtwd = weighted, translated data premultiplied by trans(w.J)
    !   tol0 = rms misfit associated with current model
    ! calls:
    !   Atamul, mult, trmult, computefwd, anorm

    
    ! Arguments:
          real(RealPrec) tol0
          
    ! Local variables
          real(RealPrec), allocatable :: dhat(:)
          integer i !, J

    
          allocate( wj(nd,nParams), dhat(nd), stat=i )
          if (i .ne. 0) then
            write(*,*) 'out of memory.  Too many free parameters (', nParams, ')'
            stop 
          endif

    ! Calculate the matrix of partials and model response
         ! Convert back to Bound parameters:
        ! write(*,*) 'pm(1:10) b4:',pm(1:10)
         if (lConstrainParams) then
             call computefwd( .True., transformToBound(pm)) 
         else               
             call computefwd( .True., pm)   ! Kwk sept 2008.  Output now done through model occam_interface_mod 
         endif
        !write(*,*) 'pm(1:10) aft:',pm(1:10)
    
    
          ! Calc misfit vector and misfit
          dwk1 = (d - dm) / sd      ! Array math - stack alloc of only nd size
    
          tol0 = sqrt(anorm(nd,dwk1)/nd)
    
          ! Weight the jacobian matrix
          do i = 1,nd
            wj(i,:) = wj(i,:) / sd(i)  
          enddo
    
          ! Scale Jacobian if using nonlinear transformation for constraining model bounds:
           if (lConstrainParams) then
              call transformWJ()
           endif
          
          ! Form w.J.Trans(w.J)
          call atamul(nd,nd,nParams,nParams,wj,wjtwj)  
    
          ! Form the weighted, translated data and premultiply by trans(w.J)
         ! call mult(nd,nd,nParams,nParams,1,wj,pm,dhat)     
          dhat = matmul(wj,pm)
          dhat = dhat + dwk1
    
          call trmult(nd,nd,nParams,nParams,1,wj,dhat,wjtwd)
    
          deallocate( wj, dhat )
          return
    end subroutine makjtj


    !-----------------------------------------------------------------------
    subroutine makptp()
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    
    ! Constructs penalty matrix (del) multiplied by its transpose.
    ! The penalty matrix only needs to be around long enough to do
    ! This, so by stuffing this into a subroutine dynamic memory allocation
    ! Can get rid of it afterwards
    
    ! On input:
    !   nothing anymore
    ! On output:
    !   /Result/ ptp() = del(trans)*del


        ! Local variables
        integer                 :: nPTerms, i
        integer, allocatable    :: linkpr(:,:)
        logical                 :: b
        real(RealPrec), allocatable       :: del(:,:)
        
        ! Functions ! External, must be in forward code!
        integer countPenaltyTerms
        
        ! How many penalty terms are there?
        nPTerms = countPenaltyTerms(  )

        
        allocate( del(nPTerms,nParams), linkpr(nPTerms,2), stat=i )
        if (i .ne. 0) then
            write(*,*) 'out of memory.  Too many free parameters (', nParams, ')'
            stop 
        endif
        ! Make partials matrix:
        b = .true.
        call constructPenaltyMat( del, nPTerms, b, linkpr)
    
    
        ! Form trans(del).del
        call deltdel(nPTerms,nParams,del,ptp,linkpr)
        
        ! Add the weights for the model prejudice to the diagonal
        forall (i = 1:nParams, prewts(i) .ne. 0)        ! Prewts is mostly zero!
            ptp(i,i) = ptp(i,i) + prewts(i)**2
        end forall
    
        deallocate( del, linkpr )
        
        return
    end subroutine makptp


    !-----------------------------------------------------------------------
    subroutine deltdel(ma,na,a,c,lpair)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Catherine DeGroot-Hedlin IGPP/SIO La Jolla Ca 92093-0225
    ! Subroutine revision 3.0, November 2006
    
    ! Multiplies the transpose of del by del 
    ! Takes advantage of sparseness of the del matrix 
    ! (Ie multiply out only connected bricks)
    
          integer ma,na,i,j,kk,k,lpair(ma,2)
          real(RealPrec) a(ma,na),c(na,na), cij
    
          c = 0.0
    ! Connected pairs
          do kk=1,ma
             i=lpair(kk,1)
             j=lpair(kk,2)
             cij = 0.0
             do k=1,ma
                cij = a(k,i)*a(k,j) + cij
             enddo
             c(i,j) = cij
             c(j,i) = cij
          enddo
    ! Diagonal component
          do i=1,na
             cij = 0.0
             do k=1,ma
                cij = a(k,i)*a(k,i) + cij
             enddo
             c(i,i) = cij
          enddo
          return
    end subroutine deltdel


    !-----------------------------------------------------------------------
    real(RealPrec)  function fndruf()      
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    
    ! Computes the roughness of the model. For memory allocation considerations,
    ! The penalty matrix has not been kept around, so we need to reform it.
    
    ! On input:
    !   /Result/ pwk2() = location of the model we want to measure
    ! On output:
    !   Fndruf = roughness of model
    ! calls:
    !   constructPenaltyMat, trmult
    !
    ! Functions
    integer countPenaltyTerms
    !  Real(RealPrec) anorm
    
    ! Local variables
    integer :: i, nPTerms, linkpr(1,1)
    logical :: b
    real(RealPrec), allocatable :: del(:,:), vect(:)
    
    !------------------------------------------------
    ! How many penalty terms are there?
    nPTerms = countPenaltyTerms(  )
    allocate( del(nPTerms,nParams), vect(nPTerms), stat=i )
    if (i .ne. 0) then
        write(*,*) 'out of memory.  Too many free parameters (', nParams, ')'
        stop 
    endif
    
    ! Make partials matrix:
    b = .false.
    call constructPenaltyMat( del, nPTerms, b, linkpr )
        ! DGM Aug 2009 - if the "optional" parameter (linkpr) is not passed
        !   then the value in nPTerms CHANGES!  This affects the MatLab
        !   interface by causing a segmentation fault.
    
    ! Calculate the new roughness (del * model2)
    do i = 1,nPTerms
        vect(i) = sum( del(i,1:nParams)*pwk2(1:nParams) ) ! Array math
    enddo
    fndruf = anorm( nPTerms, vect )
    
    ! DGM 2/2011 in a minimum norm inversion where rather than penalizing
    !   the sqrd sum of the model params, the prejudice is used to alone
    !   to penalize the difference from a background model, the roughness
    !   is calculated to be zero because the penalty matrix is all zeros.
    ! So calculate the roughness using the prejudice instead.
    ! Note: because of floating point cannot just compare fndruf==0. Rats.
    if( fndruf < 0.01 .and. iRufType == eSumParam ) then
        fndruf = sum( (prewts*premod - pwk2)**2 )
    endif
    
    deallocate( del, vect )
    
    return 
    
    end function fndruf


    !-----------------------------------------------------------------------
    subroutine minbrk(ax,bx,cx,fa,fb,fc,func)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    !
    ! Minbrk brackets a univariate minimum of a function. To be used prior to
    ! A univariate minimisation routine.
    ! Based on a routine in numerical recipes by press et al
    ! Modified so that the model associated with the misfits is carried around
    ! For use in the minimisation routines and possibly ultimately kept as the
    ! Result of this iteration;
    !
    ! On input:
    !   Ax, bx = two distinct estimates of the minimum's whereabouts
    !   Func = the function in question
    ! On output:
    !   Ax,bx,cx = three new points which bracket the minimum
    !   Fa,fb,fc = func(ax), func(bx) etc.
    !       Apparently the minimum value is in fb, such that fb <= fa and fb < fc.
    !       No guaranteed relationship between fa and fc.
    !   /Result/ pwk1() = the model associated with bx and fb
    !

    !
    ! Arguments:
          real(RealPrec) ax,bx,cx,fa,fb,fc
    ! Local variables:
    !      Integer i
          real(RealPrec) dum,r,q,ulim,u,fu
    ! Local parameters
          real(RealPrec) gold,glimit,tiny
          parameter (gold=1.618034, Glimit=100., Tiny=1.E-32)
          
    ! Functions:
          real(RealPrec),external:: func
    !
    !-----------------------------------------------------------------------
    
    
          fb = func(bx)
    ! Make pwk1 associated with b
          pwk1 = pwk2   ! Array math
          dwk1 = dm     ! Array math
          if (npm_assoc > 0) then  ! Pass associated parameters
             dwk1_assoc = pm_assoc 
          endif  
    ! pm2 will be associated with a
          fa = func(ax)
          if (fb>Fa) then
    ! Switch a and b
            dum = ax
            ax = bx
            bx = dum
            dum = fb
            fb = fa
            fa = dum
    ! Keep pwk1 associated with b (a will be lost so doesn't matter)
            pwk1 = pwk2     ! Array math
            dwk1 = dm       ! Array math
            if (npm_assoc > 0) then  ! Pass associated parameters
               dwk1_assoc = pm_assoc 
            endif          
          endif
          ! At this point, fb <= fa
          
          cx = bx + gold*(bx - ax)
          fc = func(cx)
    ! Main test. If it fails we exit with b
    ! I.E. If fb <= fa and < fc then exit.
    ! Note: nothing inferred about fc wrt fa.
    10    if (fb>=Fc) then
    !      But in here: fc <= fb <= fa
    ! Keep pwk1 associated with c 
            pwk1 = pwk2     ! Array math
            dwk1 = dm       ! Array math
            if (npm_assoc > 0) then  ! Pass associated parameters
               dwk1_assoc = pm_assoc 
            endif          
            ! Guess a new value (u) near bx for mu to try.
            R = (bx - ax)*(fb - fc)
            q = (bx - cx)*(fb - fa)
            u = bx-((bx-cx)*q-(bx-ax)*r)/(2.*Sign(max(abs(q-r),tiny),q-r))
            ulim = bx + glimit*(cx - bx)
            ! If "delta" (stuff subtracted from bx above) is positive
            ! then this is if u > cx.
            ! If delta is negative then this is if cx > u
            ! I "think" this is trickery to account for ax,bx,cx being either directly
            ! Or inversely related to the output of its functional fa,fb,fc
            if ((bx - u)*(u - cx)>0.) then
              fu = func(u)
              if (fu<Fc) then
              ! Have: fu < fc <= fb <= fa
              ! Make: a=b & b=u (want fb < fc & fb <= fa)
                ax = bx
                fa = fb
                bx = u
                fb = fu
    ! Keep pwk1 associated with b 
                pwk1 = pwk2     ! Array math
                dwk1 = dm       ! Array math
                if (npm_assoc > 0) then  ! Pass associated parameters
                   dwk1_assoc = pm_assoc 
                endif              
                go to 10
              else if (fu>Fb) then
              ! Have: fc <= fb < fu and <= fa (fu may be >,=,< fa)
              ! Make: fc=u so fb < fc and fb <= fa
                cx = u
                fc = fu
                go to 10
              endif
              ! Have: fc <= fu, fb <= fu
              ! So, take the low value (fc) and try again
              u = cx + gold*(cx - bx)
              fu = func(u)
              
            else if ((cx - u)*(u - ulim)>0.) then
              fu = func(u)
              if (fu<Fc) then
                bx = cx
                cx = u
                u = cx + gold*(cx - bx)
                fb = fc
                fc = fu
    ! Keep pwk1 associated with b 
                pwk1 = pwk2     ! Array math
                dwk1 = dm       ! Array math
                if (npm_assoc > 0) then  ! Pass associated parameters
                   dwk1_assoc = pm_assoc 
                endif               
                fu = func(u)
              endif
            else if ((u - ulim)*(ulim - cx)>=0.) then
              u = ulim
              fu = func(u)
            else
              u = cx + gold*(cx - bx)
              fu = func(u)
            endif
            
            ! It appears that a somewhat linear relationship is being assumed
            ! Between successive values of out = func(in).  Above bits that 
            ! Have fallen through have used fc to generate u and then called
            ! Fu = func(u).  And here we presume this generates a lower value
            ! Than fc and "pop" a off the top.
            ax = bx
            bx = cx
            cx = u
            fa = fb
            fb = fc
            fc = fu
            go to 10
          endif
          return
    end subroutine minbrk


    !-----------------------------------------------------------------------
    real(RealPrec) function fminocc(ax,bx,cx,fbx,f,tol,xmin)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    !
    ! Fminocc returns the minimum value of a function
    ! Based on a routine in numerical recipes by press et al 
    ! Modified so that the model associated with the misfits is carried around
    ! For use possible ultimate use as the
    ! Result of this iteration;
    !
    ! On input:
    !  Ax,bx,cx = independent variable which bracket the minimum
    !  Fbx = f(bx) (usually available from the bracketing procedure)
    !  F = function in question
    !  Tol = fractional tolerance required in the independent variable
    !  /Result/ pwk1 = model associated with bx and fbx
    ! On output:
    !  Xmin = abscissa of minimum
    !  Fminocc = f(xmin)
    !  /Result/ pwk1 = model associated with xmin and fminocc
    !
    ! Includes:

          external f
    !
    ! Arguments:
          real(RealPrec) ax,bx,cx,fbx,tol,xmin
    ! Local variables
          real(RealPrec) aa,b,v,w,x,e,fx,fv,fw,tol1,q,r,p,etemp,dd,fu,u,xm,tol2
          integer iter !, I
    ! Local parameters
          integer itmax
          real(RealPrec) cgold,zeps
          parameter (itmax=100,cgold=.3819660,Zeps=1.0E-10)
    ! Function:
          real(RealPrec) f
    !
    !-----------------------------------------------------------------------
          aa = min(ax,cx)
          b = max(ax,cx)
          v = bx
          w = v
          x = v
          e = 0.
          fx = fbx
          fv = fx
          fw = fx
          do 11 iter = 1,itmax
            ! Calc a new value (u) to run the function on.
            xm = 0.5*(Aa + b)
            tol1 = tol*abs(x) + zeps
            tol2 = 2.*Tol1
            if (abs(x - xm)<=(Tol2 - .5*(B - aa))) then
                goto 30     ! Exit!
            End if
            if (abs(e)>Tol1) then
              r = (x - w)*(fx - fv)
              q = (x - v)*(fx - fw)
              p = (x - v)*q - (x - w)*r
              q = 2.*(Q - r)
              if (q>0.) P =  - p
              q = abs(q)
              etemp = e
              e = dd
              if (abs(p)>=Abs(.5*Q*etemp).Or.P<=Q*(aa - x) .Or. &
         &        P>=Q*(b - x)) goto 1
              dd = p/q
              u = x + dd
              if (u - aa<Tol2 .Or. B - u<Tol2) dd = sign(tol1,xm - x)
              goto 2
            endif
    1       if (x>=Xm) then
              e = aa - x
            else
              e = b - x
            endif
            dd = cgold*e
    2       if (abs(dd)>=Tol1) then
              u = x + dd
            else
              u = x + sign(tol1,dd)
            endif
            
            ! Run the function
            fu = f(u)
            
            ! Evaluate the result - looking for a minimum
            if (fu<=Fx) then
              if (u>=X) then
                aa = x
              else
                b = x
              endif
              v = w
              fv = fw
              w = x
              fw = fx
              x = u
              fx = fu
              pwk1 = pwk2   ! Array math
              dwk1 = dm     ! Array math
              if (npm_assoc > 0) then  ! Pass associated parameters
                 dwk1_assoc = pm_assoc 
              endif             
            else
              if(u<X) then
                aa = u
              else
                b = u
              endif
              if (fu<=Fw .Or. W==X) then
                v = w
                fv = fw
                w = u
                fw = fu
              else if (fu<=Fv .Or. V==X .Or. V==W) then
                v = u
                fv = fu
              endif
            endif
    11    continue
          write(*,*) 'maximum iterations exceeded in fminocc'
          
    30    xmin = x
          fminocc = fx
          return
          
    end function fminocc


    !-----------------------------------------------------------------------
    real(RealPrec) function froot(func,x1,x2,fa,fb,off,tol)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.00, 13 May 1992
    !
    ! Froot finds the point at which a univariate function attains a given
    !  Value
    ! Based on a routine in numerical recipes by press et al 
    ! Modified so that the model associated with the misfits is carried around
    ! For use possible ultimate use as the
    ! Result of this iteration;
    !
    ! On input:
    !  func = the function in question
    !  x1,x2 = independent variables bracketing the root
    !  fa,fb = func(x1),func(x2) (usually available from the bracketing)
    !  off = value of the functional at the required abscissa
    !  tol = fractional tolerance required for abscissa
    !  /Result/ pwk2 = model associated with x2 and fb
    !  /Result/ pwk1 = model associated with x1 and fa
    !
    ! On output:
    !  Fb = func(froot) - off
    !  Froot = abscissa required
    !  /Result/ pwk2 = model associated with fb, froot
    !


    !
    ! Arguments:
          real(RealPrec) :: x1,x2,fa,fb,off,tol
    ! Local variables
          real(RealPrec) :: aa,b,fc,c,e,dd,tol1,p,q,r,s,xm
          integer        :: iter 
          integer        :: itmax
          real(RealPrec) :: eps
          parameter (itmax = 100,eps = 3.E-8)
          
    ! Functions:
          real(RealPrec), external  :: func
    !
    !-----------------------------------------------------------------------
          Aa = x1
          b = x2
          fa = fa - off
          fb = fb - off
          if (fb*fa>0.) write(*,*) 'root not bracketed in froot'
          fc = fb
    ! Models will be carried around as:
    !   Aa: pwk1()
    !   B: pwk2()
    !   C: pwk3()
    ! So copy pm2 (with b) into wk (with c)
          pwk3 = pwk2   ! Array math
          dwk3 = dm     ! Array math
          if (npm_assoc > 0) then  ! Pass associated parameters
              dwk3_assoc = pm_assoc 
          endif      
          do 11 iter = 1,itmax
    ! If b and c are not on different sides of root then make them  so
            if (fb*fc>0.) then
              c = aa
              fc = fa
              dd = b - aa
              e = dd
    ! Copy pmf (with aa) into wk (with c)
              pwk3 = pwk1   ! Array math
              dwk3 = dwk1   ! Array math
              if (npm_assoc > 0) then  ! Pass associated parameters
                 dwk3_assoc = dwk1_assoc
              endif               
            endif
            if (abs(fc)<Abs(fb)) then
    ! Dgm oct 2006 - found this code was rotating values
    ! Using the old a=b, b=c, c=a bug.  However, sc confirmed that
    ! This is the way it is in numerical recipes.  Inversion test
    ! Showed that it made no difference.  So leaving it the confusing
    ! Way just in case it does make a difference in some odd case.
    ! If b is not closer than c make it so
              aa = b
              b = c
              c = aa
              fa = fb
              fb = fc
              fc = fa
    ! Switch pm2 (with b) and wk (with c) (aa will be lost if we don't exit)
              pwk1 = pwk2   ! Array math
              pwk2 = pwk3
              pwk3 = pwk1
              dwk1 = dm
              dm   = dwk3
              dwk3 = dwk1
              if (npm_assoc > 0) then  ! Pass associated parameters
                 dwk1_assoc = pm_assoc 
                 pm_assoc    = dwk3_assoc
                 dwk3_assoc = dwk1_assoc
              endif   
              
            endif
            tol1 = 2.*Eps*abs(b) + 0.5*Tol
            xm = .5*(C - b)
            if (abs(xm)<=Tol1 .Or. Fb==0.)  then
              froot = b
              return
            endif
            if (abs(e)>=Tol1 .and. Abs(fa)>Abs(fb)) then
              s = fb/fa
              if (aa==C) then
                p = 2.*Xm*s
                q = 1. - S
              else
                q = fa/fc
                r = fb/fc
                p = s*(2.*Xm*q*(q - r) - (b - aa)*(r - 1.))
                Q = (q - 1.)*(R - 1.)*(S - 1.)
              Endif
              if (p>0.) Q =  - q
              p = abs(p)
              if (2.*P < Min(3.*Xm*q - abs(tol1*q),abs(e*q))) then
                e = dd
                dd = p/q
              else
                dd = xm
                e = dd
              endif
            else
              dd = xm
              e = dd
            endif
            aa = b
            fa = fb
    ! Copy pm2 (with b) into pmf (with aa)
            pwk1 = pwk2     ! Array math
            dwk1 = dm
            if (npm_assoc > 0) then  ! Pass associated parameters
               dwk1_assoc = pm_assoc 
            endif           
            if (abs(dd) > Tol1) then
              b = b + dd
            else
              b = b + sign(tol1,xm)
            endif
    ! The new funcion evaluation automatically uses pm2 for b
            fb = func(b) - off
    11    continue
          write(*,*) 'maximum iterations exceeded in froot'
          froot = b
          return
    end function froot


    !-----------------------------------------------------------------------
     subroutine cholsl(la, n, a, x, b, istage)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! R.L. Parker igpp/sio la jolla ca 92093-0225
    ! Subroutine revision 3.0, November 2006
    !
    !$$$$$ calls no other routines
    !  Routine for back-substitution solution of a linear system that has
    !  Already been factored into cholesky factors by 'cholin' (q.V.)
    !  La      leading dimension of a
    !  N       the order of matrix
    !  A       array containing on and above its diagonal the transpose of
    !          The cholesky factor.  The normal way of finding the factors
    !          Is to call  cholin  as follows
    !     call cholin(n, a, 1, ifbad)
    !     call cholsl(n, a, x, b, 0)
    !  X       the unknown vector.  Solution to equations is returned here
    !  B       the right side of system
    !  Istage  the mode of the solution.  If=1 do not complete back-
    !     Substitution.  Instead return inverse(l-trans)*b .  Other values
    !     Form the complete solution.
    ! R.L. Parker
    !

    ! Arguments
          integer la,n,istage
          real(RealPrec) a(la,n+1),x(n),b(n)
    ! Local variables
          integer i,ii !,K,kk
          real(RealPrec) z
    !
    !-----------------------------------------------------------------------
        do i=1,n
            z=b(i)
            if (i .ne. 1) then
                z = z - sum( a(1:i-1,i+1) * x(1:i-1) )
    !            Do kk=2,i
    !                K=i-kk+1
    !                Z=z - a(k,i+1)*x(k)
    !            Enddo
            endif
            x(i)=z/a(i,i+1)
        enddo
        if (istage==1) Return
        
    !  Complete back-substitution
        do ii=1,n
            i=n-ii+1
            z=x(i)
            if (i .ne. N) then
                z = z - sum( a(i,i+2:n+1) * x(i+1:n) )
    !            Do k=i+1,n
    !                Z= z - a(i,k+1)*x(k)
    !            Enddo
            end if
            x(i)=z/a(i,i+1)
        enddo
        return
    end subroutine cholsl


    !-----------------------------------------------------------------------
    subroutine cholin(la,n, a, jstage, ifbad)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! R.L. Parker igpp/sio la jolla ca 92093-0225
    ! Subroutine revision 3.0, November 2006
    !
    !$$$$$ calls no other routines
    !  Inversion and factorization of positive definite symmetric matrices
    !  From algol routine of wilkinson & reinsch 'linear algebra' pp16-17.
    !  Explanation
    !  Cholin  calculates the factor  l  in the cholesky factorization
    !  A = l l-trans, or the inverse of l , or the inverse of  a.
    !  Transposed storage is for the convenience of fortran users.  However,
    !  It means that l-trans and inverse of l-trans are computed.
    !  Arguments
    !    A     an array dimension n,n+1.  On entry the input matrix should
    !          Should occupy the 1st  n cols, or if desired, just the lower
    !          Lower left triangle.  On exit the upper triangle of the array
    !          Defined by cols 2,3...N+1 contains the desired result.
    !          This can be conveniently accessed as a(1,2), or a(n+1) if
    !          Singly dimensioned in calling program.
    !    N     the order of the matrix
    !    La    leading dimension of a
    !  Jstage  indicates desired calculation - if=1 upper triangle contains
    !          L-trans.  If=2 upper triangle contains inverse of l-trans.
    !          If=0 or 3 contains the upper half of the inverse of a.  Note
    !          Lower triangle is not destroyed. If=-2,-3, the cholesky
    !          Factors are assumed to be in place on entry.
    !  Ifbad   error flag, set to 1 if matrix not positive definite. Set to
    !          0 If ok.
    ! R.L. Parker
    !

        
        ! Pass parameters
        integer la,n,jstage,ifbad
        real(RealPrec) a(la,n+1)
        
        ! Local variables
        integer i,j,k,kback,istage,i1,j1,j2,n1
        real(RealPrec) x,y,z
    !
    !-----------------------------------------------------------------------
    !  Stage 1 - form factors (here l-trans is stored above diagonal)
        ifbad  = 0
        istage = iabs(jstage)
        if (jstage>=0) then
            do i=1,n
                i1=i + 1
                do j=i,n
                    j1=j + 1
                    x=a(j,i)
                    if (i.ne.1) then
                        ! Dgm: array math here is found to be 25% faster in tests!
                        ! Drawback: stack space allocated equal to 3 times a(1:i-1,1)
                        !   - Once for each sub-array, then for the sum.
                        X = x - sum( a(1:i-1,j1)*a(1:i-1,i1) )
    !                    Do kback=2,i
    !                        K=i-kback+1
    !                        X=x - a(k,j1)*a(k,i1)
    !                    End do
                    endif
                    if (i.ne.j) then
                        a(i,j1)=x*y
                    else
                        if (x<=0.0) Go to 4000
                        y=1.0/Sqrt(x)
                        a(i,i1)=y
                    endif
                enddo
            enddo
            if (istage.ne.1) Go to 2000
        endif
        
        do i=1,n
            a(i,i+1)=1.0/A(i,i+1)
        enddo
        if (istage==1) Return
        
        ! Nb: as of oct 2006, rest of routine never used by occam2d.
        ! Always called with jstage==1
        
    !  Stage 2 - inversion of l-trans
     2000 do 2600 i=1,n
          i1=i+1
          if (i1>N) go to 2600
          do 2500 j=i1,n
          z=0.0
          J1=j + 1
          j2=j - 1
          do 2200 kback=i,j2
          k=j2-kback+i
     2200 z=z - a(k,j1)*a(i,k+1)
          a(i,j1)=z*a(j,j1)
     2500 continue
     2600 continue
          if (istage==2) Return
    !  Stage 3  - construction of inverse of  a  above diagonal.
          do 3500 i=1,n
          do 3500 j=i,n
          z=0.0
          N1=n + 1
          j1=j + 1
          do 3200 k=j1,n1
     3200 z=z + a(j,k)*a(i,k)
          a(i,j+1)=z
     3500 continue
          return
    !  Error exit
     4000 ifbad=1
    ! Error message supressed - s.Constable
    !      write(*,400)
    ! 400  Format(44h0cholin fails - matrix not positive definite           )
          return
    end   subroutine cholin

    !-----------------------------------------------------------------------
    subroutine trmult(mad,ma,nad,na,nb,a,b,c)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Multiplies the transpose of a 2d matrix by another matrix


        integer, intent(in)         :: mad, ma, nad, na, nb
        real(RealPrec), intent(in)  :: a(mad,na),b(mad,nb) 
        real(RealPrec), intent(out) :: c(nad,nb)  
        integer :: i, j 
        
          do i = 1,na
            do j = 1,nb
               c(i,j) = sum( a(1:ma,i) * b(1:ma,j) )
            enddo
          enddo
          return
          
    end subroutine trmult

    !-----------------------------------------------------------------------
    real(RealPrec) function anorm(n,d)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! Subroutine revision 3.0, November 2006
    ! Returns the square of the euclidean norm of a vector

          integer        :: n 
          real(RealPrec) :: d(n)

          anorm = sum( d * d )
          
          return
          
    end function anorm

    !-----------------------------------------------------------------------
    subroutine atamul(mad,ma,nad,na,a,c)
    !-----------------------------------------------------------------------    
    ! Occam 3.0 Package
    ! C. Degroot-hedlin igpp/sio la jolla ca 92093-0225
    ! Subroutine revision 3.0, November 2006
    !
    ! Description: 
    !
    ! Multiplies the transpose of a 2d matrix by itself
    ! Takes advantage of symmetry about diagonal
    !
    
    ! Arguments
    integer, intent(in)         :: mad,ma,nad,na
    real(RealPrec), intent(in)  :: a(mad,na)
    real(RealPrec), intent(out) :: c(nad,na)
    
    ! Local:
    integer :: i,j 
 
    !-----------------------------------------------------------------------
    ! Note: cannot simply do the below.  It allocs a-size on the stack 
    !       Twice: once as transpose, then as output of matmul.
    !    C = matmul( transpose(a), a )
    !    Return
    ! This only allocs one a-size, but this is still too much when 
    !   A is very large - which it usually is.
    !    C = transpose(a)
    !    C = matmul( c, a )
    !    Return
        
          do i = 1,na
            do j = i,na
              ! Note assumption of symmetry
              c(i,j) = sum( a(1:ma,i)*a(1:ma,j) )
            enddo
          enddo
          ! Copy into the symmetric part of the matrix
          do i=2,na
            c(i,1:i-1) = c(1:i-1,i)
          enddo
          return
          
    end subroutine atamul
            
    !-----------------------------------------------------------------------
    real(RealPrec) function tofmu(alogmu)
    !-----------------------------------------------------------------------
    ! Occam 3.0 Package
    ! Steven Constable, Kerry Key, David Myer IGPP/SIO La Jolla CA 92093-0225
    ! 
    ! Subroutine revision 3.0, November 2006
    ! Subroutine revision 2.01, 13 Jan 1993
    
    ! Function tofmu returns the rms misfit of the response of a model 
    ! Constructed using the given value of lagrange multiplier.
    
    ! On input:
    !   Alogmu = log10 of the lagrange multiplier 
    !            (-99. If a zero lagrange multiplier needed)
    !   Stuff in occam_interface_mod 
    
    ! On output
    !   Tofmu = R.M.S.. Misfit or 1.0E+10 if cholesky decomposition failed
    !   /Result/ pwk2() contains the model which produces the misfit tofmu
    !   /Model/ dm() contians the response of the model 
    
    ! Subroutines called:
    !   Anorm, computefwd are explained in occam
    
    !  At version 1.4 Log10(mu) used to improve performance of 1d optimisation
    

    
    ! Arguments:
    real(RealPrec) alogmu
    
    ! Local variables
    real(RealPrec) amu, part
    integer istat
    
    !Lapack        
    integer, dimension(:), allocatable :: idummy
    
    ! Functions:
    !Real(RealPrec) anorm
    
    
    !-----------------------------------------------------------------------
    
    
    if (alogmu == -99.) then
        ! The gauss step is required
        amu = 0.0
    else
        ! Dgm 5/2009 - found this had just > and < which is useless since
        !   A floating point cannot be > huge or < tiny by defn. Also found
        !   That when this occurs, the inversion is screwed because the
        !   Misfit surface has flattened out. So end.
        if ( alogmu >= log10(huge(alogmu)) )  then
            write(*,*) 'amu = huge(). Has gone to infinity.'
            Stop 'lagrange multiplier went to infinity. Misfit surface has flattened. Stopping.'
            Alogmu = log10(huge(alogmu)) ! Kwk debug: keep mu from overflowing
        elseif ( alogmu <= log10(tiny(alogmu)) )  then
            write(*,*) 'amu = tiny(). Has gone to 1/infinity.'
            Stop 'lagrange multiplier went to 1/infinity. Misfit surface has flattened. Stopping.'
            Alogmu = log10(tiny(alogmu))
        endif
        amu = 10.**(Alogmu)
        !
        ! KWK: may 1, 2009:
        ! If  amu< min_tol or amu > max_tol, just output an artificially large misfit and return
        ! This should keep the search algorithm away form mu = 0 and mu = infinity
        !
        if (amu > 1d8 ) then   ! Mu tolerance = 1d20?
            tofmu = amu
            return  
            
        elseif (amu < 1d-6) then
            tofmu = 1./amu
            return  
            
        endif
        
    end if
    nfor = nfor + 1
      
    ! Add amu.Trans(partial).Partial to trans(w.J).W.J
    ! Also add model prejudice to wjtwd()
    ! Nb: memory allocaton drawback to the below array math.  Temp
    !     Stack space is allocated for the rhs of real(nParams) size,
    !     Possibly two or three times.
    pwk4  = wjtwd + amu*prewts*premod  ! Array math
    amat  = 0.0
    Amat(:,1:nParams) = amu * ptp + wjtwj
    
    ! Solve the linear system
    ! Eqn: [mu.Dtd + (wj)t.Wj] m = (wj)t.Wd + mu*prejudice-wts
    !    Aka:     amat   *     m = pwk4
    !    Note: prejudice also included in diagonal of dtd (follow variable ptp)
    !         Prejudice almost always zero!
    ! Soln for m returned in pwk4 -- original value is overwritten.
    
    select case (linearSolver)
    
    case('lapack_cholesky')
    !
    ! Try lapack cholesky routines:
    !
       if (RealPrec == kind(0d0) ) then  ! Double
           call dpotrf( 'u', nparams, amat, nparams, istat )
           call dpotrs( 'u', nparams, 1, amat, nparams, pwk4, nparams, istat )
       else ! (RealPrec == kind(0e0) ) Single
           call spotrf( 'u', nparams, amat, nparams, istat )
           call spotrs( 'u', nparams, 1, amat, nparams, pwk4, nparams, istat )       
       endif 
       pwk2 = pwk4
       
    case('lapack_lu')   
    !
    ! Or try lapack lu decomposition:
    !
        Allocate(idummy(nparams))
        Idummy = 0
        if (RealPrec == kind(0d0) ) then  ! Double
            call dgesv( size(amat,1), 1, amat, size(amat,1), idummy &
                      , pwk4, size(amat,1), istat )
        else ! (RealPrec == kind(0e0) ) Single
            call sgesv( size(amat,1), 1, amat, size(amat,1), idummy &
                      , pwk4, size(amat,1), istat )        
        endif
        Deallocate(idummy)
        pwk2 = pwk4    ! Rest of code expects soln to be in pwk2.  So accommodate it.
        
    case('parker_cholesky')   
    !
    ! Bob Parker's chol routines:       
        call cholin(nparams,nparams,amat,1,istat)       ! Comment out if using lapack
        
    end select
    
    if (istat .ne. 0) then
        ! If decomp failed we give warning but might still get by if tofmu
        ! Is set to a large value
        if (amu >= huge(amu)  ) then !Within 1 order magnitude of maximum
            write(*,*) 'Error: lagrange param went to infinity.   Stopping.'
             Stop 'Error: lagrange param went to infinity.  Stopping.'
        Endif
        write(*,*) 'warning:- cholesky decomp. Failed at mu = ',amu
        tofmu = 1.0E+10
    
    else
        select case (linearSolver)
        
        case('parker_cholesky')   
    
            call cholsl(nparams,nparams,amat,pwk2,pwk4,0)   ! Comment out if using lapack
            ! Amat is mtx of cholesky factors - upper right triangle + diagonal
            ! pwk4 is rhs of system of equations
            ! pwk2 is output: model
            ! 0=Do full solution
            
        
        end select

        
        ! Cut step size if necessary
        if (frac > 1.001) then
            part = 1./frac
            pwk2 = (1.-part)*pm + part*pwk2     ! Array math - possible memory allocation problem
        end if
        
        ! Dgm nov 2005
        ! If the system is configured to only allow discrete values
        ! in the model space, limit them now.   
        if (gbModelStepped) then
            pwk2 = anint( pwk2 / gnModelSteps ) * gnModelSteps
        endif
        
        ! Calculate the model response and misfit
        if (lConstrainParams) then
          call computefwd( .False., transformToBound(pwk2)) ! dont forget to convert back to bound parameters...
        else
           call computefwd( .False.,pwk2 )      
        endif

        
        dwk4 = (d - dm)/sd        ! Array math - stack alloc of only nd size
        tofmu = sqrt(anorm(nd,dwk4)/nd)
             
        ! If the misfit is nan then exit early.  Inversion has failed.
        !If (isnan(tofmu)) then
        if (tofmu >= huge(tofmu)) then
            write(*,*) 'Error: nan misfit.  Process aborted.'
            write(ioUnitOccamLogFile,*) 'Error: nan misfit.  Process aborted.'
            Stop
        endif
    endif
    
    if (abs(idebug) >= 1) then
            write(*,'(a,g12.5,A,g16.9,A,f12.5,A,f12.5)') '   Misfit: ', tofmu, '  Mu:', amu, &
                 & '   min/max params:   ',minval(pwk2),', ',maxval(pwk2)
            write(ioUnitOccamLogFile ,'(a,g12.5,A,g16.9,A,f12.5,A,f12.5)') 'tofmu: misfit =', tofmu,' amu =', alogmu, &
                 & '   min/max params:   ',minval(pwk2),', ',maxval(pwk2)
    endif
      
    return
    end function tofmu   
    
!---------------------------------------------------------------------
! Nonlinear transforms to bound model parameters:
!
! Computes x(m) and m(x) where b < m < a
!
! Kerry Key
! Scripps Institution of Oceanography
!
! March 2011    Implemented with exponential and bandpass options.
!
! Usage: 
!
! m = transformToBound(x)
! x = transformToUnbound(m)
! call transformWJ() - converts Jacobian matrix WJ using current x
!
! Uses model variables upperBound, lowerBound, cBoundsTransform
!
! cBoundsTransform must be:
!
! 'exponential' is from Habashy and Abubakar(2004); Commer and Newman (2008)
!
! 'bandpass' is one I designed in March 2011, which has a flat pass 
!            band for the sensitivity transform rather than the peak 
!            of the exponential.
!  
!---------------------------------------------------------------------
    elemental real(RealPrec) function transformToBound(x)
!---------------------------------------------------------------------
!
! Converts the unbound parameter x to the bound parameter b < m < a    
!       
      implicit none
      real(RealPrec), intent(in) :: x     
      real(RealPrec) :: a, b, c, p, q
  
      a = upperBound
      b = lowerBound
  
      select case (cBoundsTransform)
      
      case ('exponential')  
         
         if (x.le.0) then
            transformToBound =  ( a*exp( x) + b )  / ( exp( x) + 1.0 )
         else
            transformToBound =  ( b*exp(-x) + a )  / ( exp(-x) + 1.0 )
         endif
         
      case ('bandpass')
         
         c = bandPassFactor / (a-b)
         p = c*(1.0-exp(c*(b-a)))
         
         if (x.le.0) then
            q =  log( (exp(c*x)+exp(c*b)) / (exp(c*x)+exp(c*a)) )
         else
            q =  log( (1.0+exp(c*(b-x))) / (1.0+exp(c*(a-x))) )
         endif       
         transformToBound = (a*c + q) / p
         
      end select
 
    end function
!---------------------------------------------------------------------
    elemental real(RealPrec) function transformToUnbound(m)
!---------------------------------------------------------------------    
! Converts the bound parameter b < m < a to the unbound parameter x    
!
      implicit none
      real(RealPrec), intent(in) :: m     
      real(RealPrec) :: a, b, c, p
  
      a = upperBound
      b = lowerBound
  
      select case (cBoundsTransform)
      
      case ('exponential')  
        
         transformToUnbound = log(m-b) - log(a-m)
        
      case ('bandpass')
         
         c = bandPassFactor / (a-b)
         p = m*c*(1.0-exp(c*(b-a)))
         
         transformToUnbound =  (log(exp(a*c)*(exp(p)-exp(b*c)))-log(exp(a*c)-exp(p))) / c
         
      end select
 
    
    end function
 !---------------------------------------------------------------------   
    subroutine transformWJ()
 !---------------------------------------------------------------------   
 !
 ! Converts WJ from bound parameter m bound parameter x: dF/dx = dm/dx *dF/dm
 !  
 
    implicit none
    
      integer        :: i
      real(RealPrec) :: a, b, c
      real(RealPrec), allocatable ::  p(:), q(:) 
      
      a = upperBound
      b = lowerBound
      
      allocate (p(size(pm)), q(size(pm)) )
      
      select case (cBoundsTransform)
      
      case ('exponential')  
         where (pm.le.0) p = exp( pm)
         where (pm.gt.0) p = exp(-pm)
         do i = 1,nd
             wj(i,:) = wj(i,:) * (a-b)*p / (1.0 + p)**2
         enddo 
      
      case ('bandpass')
         
         c = bandPassFactor / (a-b)
 
         where (pm.le.0) 
            p = exp( c*(pm-b))
            q = exp( c*(pm-a))
         elsewhere (pm.gt.0)
            p = exp(-c*(pm-b)) 
            q = exp(-c*(pm-a)) 
         end where
         
         do i = 1,nd
             where (pm.le.0)
                 wj(i,:) = wj(i,:) * p / ( (1.0+p)*(1.0+q) )
             elsewhere (pm.gt.0) 
                 wj(i,:) = wj(i,:) * q / ( (1.0+p)*(1.0+q) )
            end where   
         enddo 
 
      end select

      deallocate(p,q)
    
    end subroutine transformWJ
    


end module Occam
