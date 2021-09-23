!==============================================================================!
! nFields = mexDipole1D( nXmtr, nFreq, nModel, nRcvr, nTalk, bPhaseConv, nDipLen, nNumIntegPts )
! 
! MatLab interface for Kerry Key's Dipole 1D CSEM forward model code.
!
! David Myer
! Copyright 2009
! Scripps Institution of Oceanography
! dmyer@ucsd.edu
!
! Updates:
! 
! June 23, 2010   KK: Fixed for compatibility with 64 bit Matlab.  Also added
! a matlab gateway routine Dipole1D.m that does error checking on the input arrays
! instead of here, that way Matlab errors nicely instead of completely crashing if it errors
! in mex.  Also commented out all int32 conversions since it seems like these will never be used.
!
! Dec 1, 2010 DGM: Keep the int32 conversions since these ARE used by some of us.
!
! This function ('mexDipole1D' in MatLab, NOT 'MEXFUNCTION') definition:
! Params:
!       nXmtr     - array with cols: X, Y, Z, Azimuth, Dip
!       nFreq     - vector of frequencies
!       nModel    - array with cols: TopDepth, Resistivity
!       nRcvr     - array with cols: X, Y, Z
!       nTalk     - (optional) 0=say nothing, 1=show fwd calls, 
!                    2=spew lots of text (default)
!       bPhaseConv   - (optional) 0=phase lag (dflt), 1=phase lead
!       nDipLen   - (optional; dflt=0) length of dipole in meters centered
!                    at the transmitter position(s).  Zero = point dipole.
!       nNumIntegPts - (optional; dflt=10) number of integration points to
!                       use along nDipLen in gaussian quadrature
!                       integration.
! Returns:
!       nFields     - array with cols: 
!                       Xmtr #, Freq #, Rcvr #, Ex, Ey, Ez, Bx, By, Bz
!                     Sorted by columns 1,2,3
!
!==============================================================================!

#define MX_COMPAT_32

#IF .NOT. DEFINED(_WIN32) .AND. .NOT. DEFINED(_WIN64)
#include "/Applications/MATLAB_R2010b.app/extern/include/fintrf.h"
#ELSE
#include "c:\Program Files\MATLAB\R2007b\extern\include\fintrf.h"
#ENDIF

subroutine MEXFUNCTION(nLHS, pLHS, nRHS, pRHS)
!DEC$ ATTRIBUTES DLLEXPORT :: MEXFUNCTION
    
    use dipole1d ! module shared with Dipole1D code: model params in, fields out
    
    implicit none
    
    !!!!!!!!!!!!!!!!!!!!
    !!! Declarations !!!
    !!!!!!!!!!!!!!!!!!!!
    ! Parameters - # & pointer to left/right-hand-side of MatLab function call
    integer(4)  :: nLHS, nRHS
    mwPointer   :: pLHS(*), pRHS(*)
    
    ! Matlab functions used by this routine
    mwSize      :: mxGetM, mxGetN, mxGetNumberOfElements
    mwPointer   :: mxGetPR, mxGetPI
    integer(4)  :: mxIsDouble, mxIsInt32, mxIsLogical, mxIsEmpty
    mwPointer   :: mxCreateDoubleMatrix ! (mwSize m, mwSize n, integer*4 ComplexFlag)
    
    mwSize      :: m, n, one
    integer*4   :: ComplexFlag
    
    ! Misc working variables
    integer(4)       :: i, j, iTx, iFreq, nTx, nFreq, nTalk
    integer(1)       :: b    ! for MatLab logical datatype
    real(8)          :: r
    character(1000)  :: c    ! for output printing
    real(8), allocatable     :: rFreq(:), nSigSite(:)
    real(8), allocatable     :: rWork(:,:), rTx(:,:)
    integer(4), allocatable  :: iWork(:,:)
    complex(8), allocatable :: cResults(:,:)
    
    one = 1  ! note that l,m,n,one are type mwSize, which ensures 32/64 bit Matlab compatibility for mxFunctions.
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Check params & set up modules !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (nRHS < 4) then
        call mexErrMsgTxt( 'mexDipole1D requires 4 inputs.  See fortran code.' );
    endif
    if (mxGetN(pRHS(1)) .NE. 5) then
        call mexErrMsgTxt( 'mexDipole1D input 1 requires 5 columns: xmtr X,Y,Z,Azimuth,Dip' );
    endif
    if (mxGetN(pRHS(3)) .NE. 2) then
        call mexErrMsgTxt( 'mexDipole1D input 3 requires 2 columns: Layer TopDepth,Resistivity' );
    endif
    if (mxGetN(pRHS(4)) .NE. 3) then
        call mexErrMsgTxt( 'mexDipole1D input 4 requires 3 columns: rcvr X,Y,Z' );
    endif
    if (nRHS >= 5 .and. .not. mxIsEmpty(pRHS(5))) then
        if (mxIsDouble( pRHS(5) ) .ne. 0) then
            call mxCopyPtrToReal8( mxGetPR( pRHS(5) ), r, one )
            nTalk = int( r, 4 )
        elseif (mxIsInt32( pRHS(5) ) .ne. 0) then
            call mxCopyPtrToInteger4( mxGetPR( pRHS(5) ), nTalk, one )
        elseif (mxIsLogical( pRHS(5) ) .ne. 0) then
            call mxCopyPtrToInteger1( mxGetPR( pRHS(5) ), b, one )
            nTalk = int( b, 4 )
        else
            call mexErrMsgTxt( 'nTalk parameter must be double, integer, or logical' )
        endif
    else
        nTalk = 1
    endif
    
    ! Start talking
    if (nTalk > 0) then
        call mexPrintf( '-----------' // achar(13) )
        call mexPrintf( 'mexDipole1D' // achar(13) )
        call mexPrintf( '-----------' // achar(13) )
        call mexPrintf( 'DMyer''s MatLab interface to KKey''s 1D CSEM forward model code. Aug 2009' // achar(13) // achar(13) )
        call mexDrawNow
    endif
    
    ! Fill out the dipole 1d module's vars
    call init_defaults_Dipole1D
    lbcomp          = .true.            ! compute B fields
    lUseSpline1D    = .true.            ! Use spline. KKey swears this gives exact results in almost every case
 
        
    if (nRHS >= 6) then
        if (mxIsEmpty(pRHS(6))) then
            continue    ! do nothing
        elseif (mxIsDouble( pRHS(6) ) .ne. 0) then
            call mxCopyPtrToReal8( mxGetPR( pRHS(6) ), r, one )
            if (r .ne. 0) phaseConvention = 'lead'
        elseif (mxIsInt32( pRHS(6) ) .ne. 0) then
            call mxCopyPtrToInteger4( mxGetPR( pRHS(6) ), i, one )
            if (i .ne. 0) phaseConvention = 'lead'
        elseif (mxIsLogical( pRHS(6) ) .ne. 0) then
            call mxCopyPtrToInteger1( mxGetPR( pRHS(6) ), b, one )
            if (b .ne. 0) phaseConvention = 'lead'
        else
            call mexErrMsgTxt( 'Input bPhaseConv type must be double, integer, or logical' )
        endif
        
        ! DGM Feb 2010 - finite dipole length now supported
        if (nRHS >= 7) then
            if (mxIsEmpty(pRHS(7))) then
                continue    ! do nothing
            elseif (mxIsDouble( pRHS(7) ) .ne. 0) then
                call mxCopyPtrToReal8( mxGetPR( pRHS(7) ), lenTx1D, one )
            elseif (mxIsInt32( pRHS(7) ) .ne. 0) then
                call mxCopyPtrToInteger4( mxGetPR( pRHS(7) ), i, one )
                lenTx1D = real( i, 8 );
            else
                call mexErrMsgTxt( 'Input nDipoleLen type must be double or integer' )
            endif
 
            ! KWK February 9, 2010:
            ! Is number of integration points for finite dipole approximation given?
            if (nRHS >= 8 .and. .not. mxIsEmpty(pRHS(8))) then
                if (mxIsDouble( pRHS(8) ) .ne. 0) then
                    call mxCopyPtrToReal8( mxGetPR( pRHS(8) ), r, one )
                    numIntegPts= int( r, 4 )
                elseif (mxIsInt32( pRHS(8) ) .ne. 0) then
                    call mxCopyPtrToInteger4( mxGetPR( pRHS(8) ), numIntegPts, one )
                else
                    call mexErrMsgTxt( 'Input numIntegPts type must be double or integer' )
                endif
                      
            endif
        endif
    endif

    ! Copy the transmitter & frequency data
    nTx     = mxGetM( pRHS(1) )
    nFreq   = mxGetM( pRHS(2) ) * mxGetN( pRHS(2) )
    allocate( rTx(nTx,5), rFreq(nFreq) )
    
    if (mxIsDouble( pRHS(1) ) .ne. 0) then
    	n =  nTx * 5
        call mxCopyPtrToReal8( mxGetPR( pRHS(1) ), rTx, n )
    elseif (mxIsInt32( pRHS(1) ) .ne. 0) then
    	n =  nTx * 5
        allocate( iWork(nTx,5) )
        call mxCopyPtrToInteger4( mxGetPR( pRHS(1) ), iWork, n )
        rTx = real( iWork, 8 )
        deallocate( iWork )
    else
        call mexErrMsgTxt( 'Transmitter datatype must be integer or double' )
    endif
    
    if (mxIsDouble( pRHS(2) ) .ne. 0) then
        n = nFreq
        call mxCopyPtrToReal8( mxGetPR( pRHS(2) ), rFreq, n )
    elseif (mxIsInt32( pRHS(2) ) .ne. 0) then
        allocate( iWork(nFreq,1) )
        call mxCopyPtrToInteger4( mxGetPR( pRHS(2) ), iWork, nFreq )
        rFreq = real( iWork(1:nFreq,1), 8 )
        deallocate( iWork )
    else
        call mexErrMsgTxt( 'Frequency datatype must be integer or double' )
    endif
    
    if (nTalk == 2) then
        write(c,'(A,I4,A)') 'Transmitters: ', nTx, achar(13)
        call mexPrintf( c )
        call mexPrintf( '      X         Y         Z   Azimuth   Dip' // achar(13) )
        do i = 1,nTx
            write(c,'(3(F9.1,1X),F7.1,1X,F5.1,A)') rTx(i,1), rTx(i,2), rTx(i,3), rTx(i,4), rTx(i,5), achar(13)
            call mexPrintf( c )
        enddo
        
        write(c,'(A,I4,A)') 'Frequencies: ', nFreq, achar(13)
        call mexPrintf( c )
        do i = 1,nFreq
            write(c,'(F8.3,1X)') rFreq(i)
            call mexPrintf( c )
        enddo
        call mexPrintf( achar(13) )
        
        if (nRHS >= 7) then
            write(c,'(A,F12.3,A)') 'Dipole Length: ', lenTx1D, achar(13)
            call mexPrintf( c )     
            
            if (nRHS >= 8) then
                write(c,'(A,I4,A)') '# integration points: ', numIntegPts, achar(13)
                call mexPrintf( c )

            endif
        endif
        
        call mexDrawNow
        
        
        
    endif
    
    ! Copy the model layer configuration
    nlay1D  = mxGetM(pRHS(3))           ! # of layers
    allocate( sig1D(nlay1D), zlay1D(nlay1D) )
    if (nTalk == 2) then
        write(c,'(A,I4,A)') 'Model Layers: ', nlay1D, achar(13)
        call mexPrintf( c )
        call mexPrintf( '   Top Depth      Sigma'//achar(13) )
    endif
    allocate( rWork(nlay1D,2) )
    if (mxIsDouble( pRHS(3) ) .ne. 0) then
   	    n = nlay1D*2
        call mxCopyPtrToReal8( mxGetPR( pRHS(3) ), rWork, n )
    elseif (mxIsInt32( pRHS(3) ) .ne. 0) then
   	    n = nlay1D*2
        allocate( iWork(nlay1D,2) )
        call mxCopyPtrToInteger4( mxGetPR( pRHS(3) ), iWork, n )
        rWork = real( iWork, 8 )
        deallocate( iWork )
    else
        call mexErrMsgTxt( 'Model datatype must be integer or double' )
    endif
    do i = 1,nlay1D
        zlay1D(i)   = rWork(i,1)
        sig1D(i)    = 1.d0 / rWork(i,2)
        if (nTalk == 2) then
            write(c,'(F12.1,1X,F12.3,A)') zlay1D(i), sig1D(i), achar(13)
            call mexPrintf( c )
        endif
    enddo
    deallocate(rWork)
    if (nTalk >= 2) call mexDrawNow
    
    ! Copy the receiver list & find SIGMA at each site's layer. It is used 
    ! for normalizing Jz in the return.
    n1D     = mxGetM(pRHS(4))           ! # of receivers
    allocate( ex1D(n1D), ey1D(n1D), jz1D(n1D) &     ! output fields
            , bx1D(n1D), by1D(n1D), bz1D(n1D) &
            , x1D(n1D),  y1D(n1D),  z1D(n1D)  &     ! receiver locations
            , nSigSite(n1D) )
    allocate( rWork(n1D,3) )
    if (mxIsDouble( pRHS(4) ) .ne. 0) then
   	    n = n1D*3
        call mxCopyPtrToReal8( mxGetPR( pRHS(4) ), rWork, n )
    elseif (mxIsInt32( pRHS(4) ) .ne. 0) then
   	    n = n1D*3
        allocate( iWork(n1D,3) )
        call mxCopyPtrToInteger4( mxGetPR( pRHS(4) ), iWork, n )
        rWork = real( iWork, 8 )
        deallocate( iWork )
    else
        call mexErrMsgTxt( 'Receiver datatype must be integer or double' )
    endif
    if (nTalk == 2) then
        write(c,'(A,I4,A)') 'Receivers: ', n1D, achar(13)
        call mexPrintf( c )
        call mexPrintf( '      X         Y         Z     LayerSigma' // achar(13) )
    endif
    do i = 1,n1D
        ! Copy the receiver location
        x1D(i)  = rWork(i,1)
        y1D(i)  = rWork(i,2)
        z1D(i)  = rWork(i,3)
        
        ! Figure out which layer the receiver is in. If it is on a boundary
        ! (which is true most of the time) use the layer above.
        nSigSite(i) = sig1D(1)
        do j = 2,nlay1D
            if (zlay1D(j) < z1D(i)) then
                nSigSite(i) = sig1D(j)
                exit
            endif
        enddo
        
        if (nTalk == 2) then
            write(c,'(3(F9.1,1X),F12.3,A)') x1D(i), y1D(i), z1D(i), nSigSite(i), achar(13)
            call mexPrintf( c )
        endif
    enddo
    deallocate(rWork)
    if (nTalk >= 2) call mexDrawNow
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!! Run forward models, collecting results !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    allocate( cResults(nTx * nFreq * n1D, 9) )
    
    ! For each transmitter location...
    i = 1       ! current row in cResults
    do iTx = 1,nTx
        ! Assign transmitter params to Dipole1D's variables
        xTx1D       = rTx(iTx,1)
        yTx1D       = rTx(iTx,2)
        zTx1D       = rTx(iTx,3)
        azimuthTx1D = rTx(iTx,4)
        dipTx1D     = rTx(iTx,5)
        
        ! For each frequency...
        do iFreq = 1,nFreq
            ! Assign frequency to Dipole1D's variables
            ftx1D   = rFreq(iFreq)
            
            ! Talk about what's going on
            if (nTalk > 0) then
                write(c,'(A,I4,A,F8.3,A)') 'Transmitter # ', iTx, '  Frequency ', ftx1D, achar(13)
                call mexPrintf( c )
                call mexDrawNow
            endif
            
            ! Run the forward solution
            call comp_Dipole1D
            
            ! Copy return values from Dipole1D into the MatLab rtn var.
            ! Jz needs to be converted to Ez using J = sigma*E
            do j = 1,n1D
                ! Cols: Xmtr #, Freq #, Rcvr #, Ex, Ey, Ez, Bx, By, Bz
                cResults(i,1) = iTx
                cResults(i,2) = iFreq
                cResults(i,3) = j
                cResults(i,4) = ex1D(j)
                cResults(i,5) = ey1D(j)
                cResults(i,6) = jz1D(j) / nSigSite(j)
                cResults(i,7) = bx1D(j)
                cResults(i,8) = by1D(j)
                cResults(i,9) = bz1D(j)
                
                i = i + 1
            enddo
        enddo
    enddo
    
    !!!!!!!!!!!!
    !!! Done !!!
    !!!!!!!!!!!!
    
    ! Allocate MatLab's return variable memory & copy the solutions to it
    m = nTx * nFreq * n1D
    n = 9
    ComplexFlag = 1
    pLHS(1) = mxCreateDoubleMatrix( m, n, ComplexFlag )
    call mxCopyComplex16ToPtr( cResults, mxGetPR(pLHS(1)), mxGetPI(pLHS(1)) &
                             , mxGetNumberOfElements(pLHS(1)) )
!write(c,*) mxGetNumberOfElements(pLHS(1))
!call mexPrintf( c // achar(13) )
!do i = 1,(mxGetNumberOfElements(pLHS(1)) / 9)
!    write(c,'(9(ES12.5,X),A)') real(cResults(i,1)), real(cResults(i,2)), real(cResults(i,3)),   &
!                               real(cResults(i,4)), real(cResults(i,5)), real(cResults(i,6)),   &
!                               real(cResults(i,7)), real(cResults(i,8)), real(cResults(i,9)), achar(13)
!    call mexPrintf( c )
!    write(c,'(9(ES12.5,X),A)') aimag(cResults(i,1)), aimag(cResults(i,2)), aimag(cResults(i,3)), &
!                               aimag(cResults(i,4)), aimag(cResults(i,5)), aimag(cResults(i,6)), &
!                               aimag(cResults(i,7)), aimag(cResults(i,8)), aimag(cResults(i,9)), achar(13)
!    call mexPrintf( c )
!enddo
    deallocate( cResults )
    
    ! Deallocate working memory
    deallocate( rTx, rFreq )
    deallocate( sig1D, zlay1D )
    deallocate( ex1D, ey1D, jz1D, bx1D, by1D, bz1D )
    deallocate( x1D, y1D, z1D, nSigSite )
    
    return
end subroutine MEXFUNCTION

!---------------------------------------------------------------------------
subroutine mexDrawNow()
    integer(4) i, mexEvalString
    i = mexEvalString( 'drawnow();' )
end subroutine mexDrawNow
