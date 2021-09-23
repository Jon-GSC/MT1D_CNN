!==============================================================================!
!============================================================ CallDipole1D.f90 !
!==============================================================================!
!
!    Copyright 2007,2008
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
!==============================================================================!
!===================================================================== runfile !
!==============================================================================!    
    module runfile
!
! Module for some of the parameters read in from the RUNFILE
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
    implicit none

    integer :: nTx, nFreq
    real(8), dimension(:), allocatable :: azimuthIn,dipIn
    real(8), dimension(:), allocatable :: xTxIn,yTxIn,zTxIn
    real(8), dimension(:), allocatable :: fTxIn, sigsite
    
    integer,parameter :: outfileunit = 16
    
    character(120)  :: outputfilename, sRunfile
    
    integer,parameter :: derivoutfileunit = 17
    
    end module runfile


!==============================================================================!
!========================================================! Program CallDipole1D 
!==============================================================================!
!
! Program that reads in a runfile and then calls Dipole1D subroutine to 
! compute CSEM fields from an arbitrarily oriented electric dipole.
!
! Version 1.1   February 27, 2008   Switched get_time_offset to use get_time_offset subroutine
!                                   Added in flexible format for reading RUNFILE
! Version 1.0   Fall, 2007
!
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
!
!==============================================================================!
    Program CallDipole1D
    
    use dipole1d     ! Model parameters are passed in here, fields passed out here
    use runfile      ! This stores the info read in from the RUNFILE
    
    implicit none
    
    integer               :: i,j, iTx, iFreq, iRxlayer
    real(8)               :: t0, t1 ! timing variables
    

!
! Initialize the default values, can be overidden with values in RUNFILE
!
   call init_defaults_Dipole1D 
   call init_defaults 
   
    ! DGM 10/2010 Allow command line arguments for the in & out filenames
    ! Note that the input filename allows an "Output FileName: xxxx" entry
    ! which will override the command line, if present.
    i = command_argument_count()
    if( i > 0 ) then
        call get_command_argument(1,sRunfile)
        
        if( i > 1 ) then
            call get_command_argument(1,outputfilename)
        else
            ! If an input name was specified, but no output name
            ! then use the input name but change the extension
            ! to '.csem'
            i = index( sRunfile, '.' )
            if( i == 0 ) then
                outputfilename = sRunfile // '.csem'
            else
                outputfilename = sRunfile(1:i-1) // '.csem'
            endif
        endif
    endif

!
! Read in the runfile
!
    call readrunfile_dipole1D   
!
! Store the all the model, Tx and site parameters to the file so that each file
! completely describes the problem solved
!
    call initialize_outputfile      
 
    if (linversion) then
        call initialize_deriv_outputfile
    endif
!
! Start the get_time_offset:
!
   write(*,*) 'Starting 1D computations...'
   write(*,*) ' '
   call get_time_offset(0d0,t0)  ! get_time_offset should work on all F90/95 compilers
   
   
!
! Loop over each transmitter:
! 
    do iTx = 1,nTx
        write(*,'(a24,i6)') 'Transmitter #:',iTx    
!
! Assign Tx parameters:
!
        xTx1D       = xTxIn(iTx)
        yTx1D       = yTxIn(iTx)
        zTx1D       = zTxIn(iTx)
        azimuthTx1D = azimuthIn(iTx)
        dipTx1D     = dipIn(iTx)
!
! Inner loop of requested frequencies
!
        do iFreq= 1,nFreq
!
! Get the requested frequency:
!
            ftx1D = fTxIn(iFreq) 
            write(*,'(a24,i6)') 'Frequency #:',iFreq            
!
! Compute CSEM fields:
!           
            call comp_dipole1D
!   
! Output the response for the current transmitter, note that I convert jz to ez here
! The assumption is that receivers on a boundary are actually in the top layer.
! So a seafloor site would use the sea layer conductivity to convert to vertical
! current to electric field.
! sigsite is defined for each site in subroutine readrunfile

        
            do i=1,n1D
                jz1D(i) = jz1D(i)/ sigsite(i)  
                write(outfileunit,100) &
                   &   real(ex1D(i)),aimag(ex1D(i)),real(ey1D(i)),aimag(ey1D(i)),real(jz1D(i)),aimag(jz1D(i)), &
                   &   real(bx1D(i)),aimag(bx1D(i)),real(by1D(i)),aimag(by1D(i)),real(bz1D(i)),aimag(bz1D(i))
                   
            enddo
            
            if (linversion) then
    
                
                do i = 1,n1D
                ! convert djzdsig to dezdsig:
                    djzdsig(i,:) = djzdsig(i,:)/sigsite(i)
                    
                    ! Get receiver layer
                    iRxlayer = 1
                    do j= 2,nlay1D
                      if (z1d(i).gt.zlay1D(j))  then   ! Sites on boundaries use the top layer 
                        iRxlayer = j
                      endif
                    enddo
                    djzdsig(i,iRxlayer) =   djzdsig(i,iRxlayer) - jz1D(i)/sigsite(i)    
                
                   do j = 1,nlay1d
                    write(derivoutfileunit,fmt='(E18.6)') dble(dexdsig(i,j)), aimag(dexdsig(i,j))
                    write(derivoutfileunit,fmt='(E18.6)') dble(deydsig(i,j)), aimag(deydsig(i,j)) 
                    write(derivoutfileunit,fmt='(E18.6)') dble(djzdsig(i,j)), aimag(djzdsig(i,j))
                    write(derivoutfileunit,fmt='(E18.6)') dble(dbxdsig(i,j)), aimag(dbxdsig(i,j)) 
                    write(derivoutfileunit,fmt='(E18.6)') dble(dbydsig(i,j)), aimag(dbydsig(i,j))
                    write(derivoutfileunit,fmt='(E18.6)') dble(dbzdsig(i,j)), aimag(dbzdsig(i,j))                   
                    enddo
                enddo
        
            endif
            
100    format(12(1x,e15.8)) 

!
! Go on to the next frequency
! 
        enddo
!
! Go on to the next Tx
! 
   enddo 
!
! Print out the run time:
!
    call get_time_offset(t0,t1)  ! get_time_offset should work on all F90/95 compilers 
    write(*,*) ' '     
    write(*,'(a32,2x,g12.4)') 'Total time for computations (s): ',t1
    write(*,*) ' ' 
   
!
! Close the output file
!
    close(outfileunit)   

    if (linversion) close(derivoutfileunit)
            
!
! Deallocate arrays
!
    deallocate ( x1D,y1D,z1D,sigsite )
    deallocate ( ex1D,ey1D,jz1D,bx1D,by1D,bz1D )
    deallocate ( sig1D, zlay1D )
    deallocate ( azimuthIn,dipIn,xTxIn,yTxIn,zTxIn,fTxIn )
    
    if (linversion) then
        deallocate (dexdsig,deydsig,djzdsig,dbxdsig,dbydsig,dbzdsig)
    endif

!
! Hasta la vista 
!
    end program CallDipole1D
  
!==============================================================================!    
!===============================================================! init_defaults 
!==============================================================================!
    subroutine init_defaults 
    
    use runfile
    use dipole1d     ! Model parameters are passed in here, fields passed out here
!
! Specify some parameters required by Dipole1D:
!
    sRunfile        = 'RUNFILE'
    outputfilename  = 'dipole1Doutput.csem'
    HTmethod1D      = 'kk_ht_201' 
    outputdomain1D  = 'spatial'
    lbcomp          = .true.
    sdm1D           = 1.0   ! (Am), dipole moment 
    lUseSpline1D    = .false.
    linversion      = .false. ! Compute Derivatives with respect to sigma(layers)
    
    
    end subroutine init_defaults 
  
!==============================================================================!    
!=============================================================! get_time_offset
!==============================================================================!
    subroutine  get_time_offset(timein,timeout)
!    
! timein is the clock start time or 0.
! timeout is the time in seconds since the input time
!
! Version 2.0  February 25, 2008  Now uses date_and_time Fortran intrinsic
!
! Written by
! Kerry Key
! Scripps Institution of Oceanography
! kkey@ucsd.edu
! ---------------------------------------------------------------------
   implicit none
   
   integer, dimension(8) :: values
   integer               :: i,j,k,mjd
   
   real(8) :: timein, timeout, fracday
    

 !
 ! New standard Fortran90 Time function:
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
               
           
    end
        
!==============================================================================! 
!=================================================================! readrunfile
!==============================================================================! 
      subroutine readrunfile_dipole1d       
 
    use runfile
    use dipole1D  
      
    implicit none
    
    integer         :: err,i,j
    character(56)   :: vers,ctemp
    real(8)         :: rhotemp
    character(180)  :: sLine, sCode, sValue
    logical         :: bComment
    
    write(*,*) ' ' 
    write(6,*) '============================= Dipole1D ==================================='
    write(*,*) ' '  
    write(*,*) '        A program for computing the EM fields from an arbitrarily ' 
    write(*,*) '          oriented electric dipole in an N-layered model using ' 
    write(*,*) '              a Lorenz gauged vector potential formulation. ' 
    write(*,*) ' '      
    write(*,*) '                        Copyright 2007-2010 '
    write(*,*) '                             Kerry Key                     ' 
    write(*,*) '                  Scripps Institution of Oceanography      '  
    write(*,*) '                             kkey@ucsd.edu                 '        
    write(*,*) ' '     
    write(*,*) '                      Version 7.2  Feb 5, 2010.             ' 
    write(*,*) ' '
    write(*,*) '                  This work was supported by: '
    write(*,*) ' '
    write(*,*) '             The Seafloor Electromagnetic Methods Consortium   '
    write(*,*) '                at Scripps Institution of Oceanography  '
    write(*,*) ' '
    write(*,*) '          See this URL for a list of SEMC funding sponsors:'
    write(*,*) ' '
    write(*,*) '                http://marineemlab.ucsd.edu/semc.html'
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
    write(6,*) '======== Reading in "RUNFILE" ============================================'
    write(*,*) ' ' 
    
!  
!  Open RUNFILE and read in modeling info:
!
    open (unit=10,file=trim(sRunfile),status='old',iostat=err)
    if (err .ne. 0) then
        write(*,*) ' Error opening RUNFILE'
        stop
    end if
!
! Feb, 2008  New version uses flexible format that is backwards compatible
! so we don't need the version check anymore:
!
!   This code will only support version 2.0 and newer
!
!    read(unit=10,fmt='(18x,a)')    vers
!    if  (vers(1:16) .ne. 'DIPOLE1D_1.0') then
!        write(*,*) 'Error: Runfile format ',vers
!        write(*,*) 'not supported, stopping.'           
!        write(*,*) 'Try using DIPOLE1D_1.0'
!        stop
!    endif

! Read in the file a line at a time and decode the sCode/sValue pairs that are separated by a semicolon
! Certain sCodes are followed by lists such as model parameters, site locations, etc.
! Parsecode is D Myer's useful code tidbit for separating sCode and sValue pairs

    do while (.true.) 
    
        ! Get the next code/value pair
        ! ParseCode forces the code portion to be all lowercase with no
        !   padding and ending colon stripped off.  User comments are 
        !   stripped from the value portion.
        
        read( 10, '(A)', iostat = err ) sLine
        
        if (err /= 0) exit  ! end of file read, escape from while loop and 
                            ! proceed to checking that required inputs defined
        
        call ParseCode( len(sLine), sLine, sCode, sValue, bComment )
        if( bComment ) cycle
        
        ! What do we have?
        select case (trim(sCode))
       
        case ('version')
       
            vers = sValue(1:len(vers))
            
        case ('dipole length','dipole','dipolelength')
            ! DGM Feb 2010 - add support for finite length dipole
            read(sValue,*) lenTx1D
            write(*,fmt = '(a24,f6.1)') 'Dipole Length:', lenTx1D
            
        case ('# integ pts','# integration points','numIntegPts','number of integration points')
            read(sValue,*) numIntegPts
            write(*,fmt = '(a24,I6)') '# integration points:', numIntegPts
            
        ! DGM 10/2010 - someone forgot this
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
            
        case ('# transmitters')
        
            read(sValue,*) nTx
            write(*,fmt='(a24,i6)') '# Transmitters:',nTx
            read(10,fmt='(a56)',end=198, err=199 ) ctemp ! skip header text
            write(6,*) ''            
            !
            !  Allocate arrays dependent on nTx
            !
                allocate ( azimuthIn(nTx),dipIn(nTx))
                allocate ( xTxIn(nTx),yTxIn(nTx),zTxIn(nTx))  
            !
            ! Read in each Tx position info line:
            !
                write(*,*) '  Transmitter Positions: '
                write(*,*) ''
                write(6,54) '#','X','Y','Z','Azimuth','Dip'
            !
            ! Loop over number of transmitters:
            !
                do i = 1,nTx  
            !
            ! Read in a line of transmitter parameters:
            !
                    read(10,*,end=198, err=199 ) xTxIn(i),yTxIn(i),zTxIn(i), azimuthIn(i), dipIn(i)
                    write(6,55) i,xTxIn(i),yTxIn(i),zTxIn(i), azimuthIn(i), dipIn(i)
            !
            ! Go on to the next transmitter
            !
                enddo      
                
        case ('# frequencies')
            
            read(sValue,*) nFreq
            write(*,*)  ' '    
            write(*,fmt='(a24,i6)') '# Frequencies:',nFreq   
            write(*,*)  ' '
            write(6,fmt='(a24,a12)') 'Frequencies:', '[Hz]:'
            allocate (fTxIn(nfreq))
            do i  = 1,nfreq
                read(10,*,end=198, err=199 ) fTxIn(i)
                write(*,fmt='(i24,e12.4)') i, fTxIn(i) 
            enddo          
            
        case ('# layers')
        
            read(sValue,*)  nlay1D
            write(*,*) ' ' 
            write(*,fmt = '(a24)') ' Background 1D Model:'
            write(*,*) ' ' 
            allocate(sig1D(nlay1D), zlay1D(nlay1D))
            write(*,56) 'Layer #','Top Depth (m)', 'Resistivity (ohm-m)'
            do i = 1,nlay1D
                read(10,*,end=198, err=199) zlay1D(i), rhotemp
                if (rhotemp == 0) then
                    write(*,*) 'Error:  zero resistivity in layer ', i
                    write(*,*) ' Please replace with non-zero value! '
                    write(*,*) ' stopping! '
                    stop
                endif
                sig1D(i) = 1.d0/rhotemp   ! convert resistivity in file to conductivity     
                write(6,57) i, zlay1D(i),  1./sig1D(i)
            enddo
            
        case ('# receivers')
         
            read(sValue,*) n1D
            write(*,*) ' ' 
            write(*,fmt = '(a24,i6)') ' Number of Receivers:',n1D
            write(*,*) ' ' 
            !
            ! Allocate the field arrays: 
            !
            allocate(x1D(n1D), y1D(n1D), z1D(n1D),sigsite(n1D) )
            allocate ( ex1D(n1D),ey1D(n1D),jz1D(n1D),bx1D(n1D),by1D(n1D),bz1D(n1D) )


        !   write(*,58) 'Rx #','x', 'y','z','[m]'
            do i = 1,n1D
                read(10,*,end=198, err=199) x1D(i),y1D(i),z1D(i)   
!               write(6,59) i, x1D(i),y1D(i),z1D(i)   
            enddo   

        case ('ht filters')
        
            HTmethod1D      = sValue(1:len(HTmethod1D))
            write(*,fmt = '(a36,a24)') ' Using Hankel TF Filters: ', HTmethod1D
            write(*,*) ' '
        case ('output filename')
        
            outputfilename = sValue(1:len(outputfilename))
            write(*,fmt = '(a36,a56)') ' Saving output to file: ', outputfilename
            write(*,*) ' '
            
        case('usespline1d')
        
            select case (trim(sValue))
            case ('yes')
                lUseSpline1D = .true.
                write(*,*) 'Using Spline Interpolation for speed'
                write(*,*) ' '
            case default
                lUseSpline1D = .false.
            end select
            
            
        case ('compderivatives')
        
            select case (trim(sValue))
            case ('yes')
                linversion = .true.
                write(*,*) 'Computing derivatives of fields wrt sigma(layers)'
                write(*,*) ' '  
            case default
                linversion = .false.
            end select      
        
             
        case default
            write(*,*) 'Error reading RUNFILE file!'
            write(*,*) 'Unknown or unsupported code:', sCode
            stop
            
        end select
        
    enddo   
        
!
! Now perform a few rudimentary checks to make sure freq, Tx, Rx and model defined:
!
    
    if (nTx < 1 ) then
        write(*,*) ' Error, no transmitters defined in RUNFILE'
        stop
    end if
    
    if (nFreq < 1 ) then
        write(*,*) ' Error, no frequencies defined in RUNFILE'
        stop
    end if

    if (nlay1D < 1 ) then
        write(*,*) ' Error, no layers defined in RUNFILE'
        stop
    end if
    
    if (n1D < 1 ) then
        write(*,*) ' Error, no receivers defined in RUNFILE'
        stop
    end if
        
!
! Allocate the inversion field derivatives
!    
    if (linversion) then
        allocate (dexdsig(n1D,nlay1D), deydsig(n1D,nlay1D), djzdsig(n1D,nlay1D)) 
        allocate (dbxdsig(n1D,nlay1D), dbydsig(n1D,nlay1D), dbzdsig(n1D,nlay1D)) 
    endif
            
!
! Now lets store the conductivity at the site location for later use.  Sites
! The assumption is that receivers on a boundary are actually in the top layer.
! So a seafloor site would use the sea layer conductivity to convert to vertical
! current to electric field.
    do i=1,n1D
        sigsite(i)  = sig1D(1) 
        do j = 2,nlay1D
            if (zlay1D(j).lt.z1D(i)) then
               sigsite(i) = sig1D(j) ! ends on last layer top being shallower than site location
            endif
        enddo       
    enddo

!
! Close RUNFILE
!
    close(10)
    write(*,*) ' '
    write(6,*) '====== Starting Computations ============',&
             & '================================='  
    write(*,*) ' '   
    
!
! Everything was read in okay, proceed with computations
!
    return
    
!
! Format and error statements:
!

54  format(a3,1x,a9,1x,a9,1x,a9,1x,a9,1x,a9)   
55  format(i3,1x,f9.1,1x,f9.1,1x,f9.1,1x,f9.1,1x,f9.1)  
56  format(a18,1x,a18,1x,a19) 
57  format(i18,1x,f18.1,1x,e18.3) 
!58  format(a6,1x,a12,1x,a12,1x,a12,1x,a4) 
!59  format(i6,1x,f12.1,1x,f12.1,1x,f12.1) 

198 write(*,*) ' Error: RUNFILE ended prematurely'
    close(10)   
    stop
        
199 write(*,*) '  Error reading RUNFILE '
    close(10)
    stop

    end subroutine readrunfile_dipole1d     
 
 
!==============================================================================!
!====================================================================! DIPOLE1D 
!==============================================================================!
    subroutine initialize_outputfile
!
! Opens and initializes the formatted output file "dipole1Doutput.csem"
! This will we store all parameters for the problem being solved
! along with the EM fields computed by dipole1D, this way the
! file is a self contained description and solution to the CSEM problem.
!
!
!  Kerry Key 
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!==============================================================================! 
    use runfile
    use dipole1D
    
    implicit none
    
    integer :: i
    
!
! Open the output file:
!
    open(unit=outfileunit,file=outputfilename,status='replace') 

!
! Write out the transmitters:
!
    write(outfileunit,*) 'Dipole1D_1.1'
    write(outfileunit,*) lenTx1D 
    write(outfileunit,*) nTx
    do i = 1,nTx  
        write(outfileunit,10) xTxIn(i),yTxIn(i),zTxIn(i), azimuthIn(i), dipIn(i)
    enddo    
10  format(f14.1,1x,f14.1,1x,f14.1,1x,f6.1,1x,f6.1) 
!
! Write out the frequencies:
!
    write(outfileunit,*) nFreq
    do i  = 1,nfreq
         write(outfileunit,*) fTxIn(i)      
    enddo
    
!
! Write out the 1D model:
!
    write(outfileunit,*) nlay1D
    do i = 1,nlay1D     
        write(outfileunit,11)  zlay1D(i),  1./sig1D(i)
    enddo
11  format(f14.1,1x,e14.6)     
!
! Write out the site locations
!
    write(outfileunit,*) n1D
    do i = 1,n1D 
        write(outfileunit,12)  x1D(i),y1D(i),z1D(i)   
    enddo       
12  format(f14.1,1x,f14.1,1x,f14.1) 
 
    end subroutine initialize_outputfile    
    
!==============================================================================!
!====================================================================!  
!==============================================================================!
    subroutine initialize_deriv_outputfile
!
! Opens and initializes the formatted output file "dipole1Doutput.csem"
! This will we store all parameters for the problem being solved
! along with the EM fields computed by dipole1D, this way the
! file is a self contained description and solution to the CSEM problem.
!
!
!  Kerry Key 
!  Scripps Institution of Oceanography
!  kkey@ucsd.edu
!
!==============================================================================! 
    use runfile
    use dipole1D
    
    implicit none
    
    integer :: i
    
!
! Open the output file:
!
    open(unit=derivoutfileunit,file=trim(outputfilename)//'_deriv',status='replace')    

!
! Write out the transmitters:
!
    write(derivoutfileunit,*) 'Dipole1D_1.0'
    write(derivoutfileunit,*) nTx
    do i = 1,nTx  
        write(derivoutfileunit,10) xTxIn(i),yTxIn(i),zTxIn(i), azimuthIn(i), dipIn(i)
    enddo    
10  format(f14.1,1x,f14.1,1x,f14.1,1x,f6.1,1x,f6.1) 
!
! Write out the frequencies:
!
    write(derivoutfileunit,*) nFreq
    do i  = 1,nfreq
         write(derivoutfileunit,*) fTxIn(i)      
    enddo
    
!
! Write out the 1D model:
!
    write(derivoutfileunit,*) nlay1D
    do i = 1,nlay1D     
        write(derivoutfileunit,11)  zlay1D(i),  1./sig1D(i)
    enddo
11  format(f14.1,1x,e14.6)     
!
! Write out the site locations
!
    write(derivoutfileunit,*) n1D
    do i = 1,n1D 
        write(derivoutfileunit,12)  x1D(i),y1D(i),z1D(i)   
    enddo       
12  format(f14.1,1x,f14.1,1x,f14.1) 
 
    end subroutine initialize_deriv_outputfile      
!==============================================================================!
!===================================================================! ParseCode
!==============================================================================!-
    subroutine parseCode( nLen, sLine, sCode, sValue, bComment )
  
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
!=======================================================================! LOWER
!==============================================================================!
    subroutine Lower( s )

! David Myer IGPP/SIO La Jolla CA 92093-0225
! DGM Nov 2006 - convert string to lower case
    implicit none
    character(*), intent(out)  :: s
    integer i

    do  i=1,len_trim(s)
      if  ( s(i:i) >= 'A' .and. s(i:i) <= 'Z' ) then
        s(i:i) = char(ichar(s(i:i)) + 32)
      endif
    enddo
    
    end subroutine Lower    
    
    