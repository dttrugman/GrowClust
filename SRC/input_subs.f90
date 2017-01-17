!-----------------------------------------------------------------------
! Copyright 2016 Daniel Trugman
!
! This file is part of GrowClust.
!
! GrowClust is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! GrowClust is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with GrowClust.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------



!--------------------------SUBROUTINES used in GROWCLUST---------------------------------!


!------------------------------------------------------------------------------------
! ***** GET_INP: Simple subroutine to get next non-commented line in input file *****
!------------------------------------------------------------------------------------
!  Inputs:  fnum  = file ID number (assumed to be open)
!  Returns: line  = next line (char*100)
!------------------------------------------------------------------------------------
    subroutine GET_INP(fnum, line)
    
    implicit none
    integer :: fnum, i
    character (len=1) :: charin
    character (len=100) :: line
    
    charin = '*'

100 continue 
    if (charin .eq. '*') then 
      read(fnum, '(a)') line
      i = 1
      charin = line(i:i)

200   continue
      if(charin.eq.' ')then
         i=i+1
         charin=line(i:i)
         goto 200
      endif
    goto 100
    
    endif

    return
    end
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------     
! **** INPUT_CHECK: Checks set of input parameters for invalid values 
! ****** (note that not all parameters are checked, only those that may prove problematic
! ****** later in the algorithm AND are not checked elsewhere. 
! ****** E.g., it's okay to have rmin < 0 ....)
!------------------------------------------------------------------------------------
!  Inputs:  plongcutP =  minimum ray parameter, P-phase
!           plongcutP =  minimum ray parameter, P-phase
!         vpvs_factor =  Vp/Vs ratio is Vs not specified
!           tt_dep1   =  min. depth in travel-time table          
!           tt_dep2   =  max. depth in travel-time table
!           tt_dep3   =  depth spacing in travel-time table
!           tt_del1   =  min. del in travel-time table          
!           tt_del2   =  max. del in travel-time table
!           tt_del3   =  del spacing in travel-time table
!           rmsmax    =  max. rms residual allowed to join clusters
!           delmax    =  max. station distance for xcor observation to count as "good"
!           iponly    =  option to use P and S (0 or 2), or P-phase only (1)
!           nboot     =  array(len=nq) of event latitude (from input catalog)
!           maxboot   =  array(len=nq) of event longitude (from input catalog)
!-----------------------------------------------------------------------------------------
   subroutine INPUT_CHECK(plongcutP, plongcutS, vpvs_factor, tt_dep1, tt_dep2, tt_dep3, &
      tt_del1, tt_del2, tt_del3, rmsmax, delmax, iponly, nboot, maxboot)
    
   implicit none
   integer :: iponly, nboot, input_ok, maxboot
   real :: plongcutP, plongcutS, vpvs_factor, tt_dep1, tt_dep2, tt_dep3
   real :: tt_del1, tt_del2, tt_del3, delmax, rmsmax
   
   print *, 'Checking input parameters...'
   input_ok = 1 ! input is ok unless problem is found
   
   ! check min ray parameters plongcutP, plongcutS
   if ((plongcutP < 0) .or. (plongcutS < 0)) then
      print *, 'Input error (velocity model param.): plongcutP, plongcut:'
      print *, plongcutP, plongcutS   
      input_ok = 0
   endif
   
   ! check vp/vs ratio
   if (vpvs_factor < 0.0) then
      print *, 'Input error (velocity model param.): vpvs_factor:'
      print *, vpvs_factor 
      input_ok = 0
   endif
   
   ! check travel-time table: dep and del
   if ((tt_dep1 < 0.0) .or. (tt_dep2 < tt_dep1) .or. (tt_dep3 > tt_dep2)) then
     print *, 'Input error (travel-time table param.) : min_dep, max_dep, d_dep'
     print *, tt_dep1, tt_dep2, tt_dep3
     input_ok = 0
   endif
   if ((tt_del1 < 0.0) .or. (tt_del2 < tt_del1) .or. (tt_del3 > tt_del2)) then
     print *, 'Input error (travel-time table param.) : min_del, max_del, d_del'
     print *, tt_del1, tt_del2, tt_del3
     input_ok = 0
   endif
   
   ! check rmsmax and delmax
   if ((rmsmax <= 0.0) .or. (delmax <= 0.0) ) then
     print *, 'Input error (GrowClust param.) : rmsmax, delmax'
     print *, rmsmax, delmax
     input_ok = 0
   endif
   
   ! check iponly
   if ((iponly < 0) .or. (iponly > 2)) then 
     print *, 'Input error (GrowClust param.) : iponly'
     print *, iponly
     input_ok = 0
   endif
   
   ! check nboot
   if ((nboot < 0) .or. (nboot > maxboot)) then 
     print *, 'Input error: nboot, maxboot'
     print *, nboot, maxboot
     input_ok = 0
   endif    
   
   ! return
   if (input_ok == 0) then
     print *, 'Ending program due to above input errors!'
     print *, 'Please check input control file for errors.'
     stop
   else
     print *, 'Input parameters ok!'
     return
   endif
   
   end       
!------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------     
! **** PARAM_CHECK: Checks set of parameters from grow_params module for invalid values 
! ****** (note that not all parameters are checked, only those that may prove problematic
! ****** later in the algorithm AND are not checked elsewhere. 
! ****** E.g., it's okay to have distmax < 0 ....)
!------------------------------------------------------------------------------------
!  Inputs:  conparam  =  minimum connection fraction to join clusters
!           hshiftmax =  maximum permitted horizontal cluster shifts (km)
!           vshiftmax =  maximum permitted vertical cluster shifts (km)
!           rmedmax   =  maximum median absolute tdif residual to join clusters          
!           boxwid    =  initial "shrinking-box" width (km) in relative relocation subroutines
!           nit       =  number of iterations in relative relocation subroutines
!           samp_type =  bootstrap sampling type (1 or 2)
!           irelonorm =  relocation norm (1=L1 norm, 2=L2 norm)           
!-----------------------------------------------------------------------------------------
   subroutine PARAM_CHECK(conparam, hshiftmax, vshiftmax, rmedmax, boxwid, nit, samp_type, irelonorm)
     
   implicit none
   real :: conparam, hshiftmax, vshiftmax, rmedmax, boxwid 
   integer :: nit, samp_type, params_ok , irelonorm 
   
   print *, 'Checking grow_params.mod parameters...'
   params_ok = 1 ! input is ok unless problem is found
     
   ! check conparam
   if ((conparam <= 0.0) .or. (conparam > 1.0) ) then
     print *, 'grow_params.mod error: conparam'
     print *, conparam
     params_ok = 0
   endif
   
   ! check hshiftmax, vshiftmax
   if ((hshiftmax <= 0.0) .or. (vshiftmax <= 0.0)) then 
     print *, 'grow_params.mod error: hshiftmax, vshiftmax'
     print *, hshiftmax, vshiftmax
     params_ok = 0
   endif
   
   ! check rmedmax
   if (rmedmax < 0.0) then 
    print *, 'grow_params.mod error: rmedmax'
     print *, rmedmax
     params_ok = 0
   endif 
   
   ! check boxwid, nit
   if ((boxwid <= 0.0) .or. (nit >=200)) then 
     print *, 'grow_params.mod error: boxwid, nit'
     print *, boxwid, nit
     params_ok = 0
   endif
   
   ! check samp_type
   if ((samp_type < 1) .or. (samp_type > 2)) then 
     print *, 'grow_params.mod error: samp_type'
     print *, samp_type
     params_ok = 0
   endif
   
   ! check irelonorm
   if ((irelonorm < 1) .or. (irelonorm > 3)) then 
     print *, 'grow_params.mod error: irelonorm'
     print *, irelonorm
     params_ok = 0
   endif
   
   ! return
   if (params_ok == 0) then
     print *, 'Ending program due to above parameter errors!'
     print *, ' Please check grow_params.f90 for errors.'
     stop
   else
     print *, 'grow_params.mod parameters ok!'
     return
   endif
   
   end       
!------------------------------------------------------------------------------------

!------------------------------------------------------------------------------------     
! **** READ_EVFILE: Reads event list to get starting (catalog) locations
! -- accepts different event list formats using irxform, istform
! -- initiliazes relocated positions (lat, lon, dep, time) to catalog positions
! -- initializes the clustertree arrays index, nbranch, tlat, tlon, tdep, 
!     assigning each event to its own unique cluster 
!-----------------
!  Inputs:  ievform   =  event list format: 0 for evlist, 1 for phase
!           evfile    =  event list file name
!           nq0       =  maximum number of events
!           
!  Returns: idcusp    =  array(len=nq) of event ID #s (e.g. CUSPID or EVID)
!           qyr_cat   =  array(len=nq) of event year (from input catalog)
!           qmon_cat  =  array(len=nq) of event month (from input catalog)
!           qdy_cat   =  array(len=nq) of event day (from input catalog)
!           qhr_cat   =  array(len=nq) of event hour (from input catalog)
!           qmn_cat   =  array(len=nq) of event minute (from input catalog)
!           qsc_cat   =  array(len=nq) of event second (from input catalog)
!           qmag_cat  =  array(len=nq) of event mag (from input catalog)
!           qlat_cat  =  array(len=nq) of event latitude (from input catalog)
!           qlon_cat  =  array(len=nq) of event longitude (from input catalog)
!           qdep_cat  =  array(len=nq) of event depth, km (from input catalog)
! 
!-----------------------------------------------------------------------------------------
   subroutine READ_EVFILE(ievform, evfile, nq0,idcusp, qyr_cat, qmon_cat, &
    qdy_cat, qhr_cat, qmn_cat, qsc_cat, qmag_cat, qlat_cat, qlon_cat, qdep_cat, nq)
   
   implicit none
   integer :: i, nq0, nq, ievform
   real :: EH_cat, EZ_cat, RMS_cat
   integer, dimension(nq0) :: idcusp
   real, dimension(nq0) ::  qdep_cat
   real (kind=8), dimension(nq0) ::  qlat_cat, qlon_cat
   integer, dimension(nq0) :: qyr_cat, qmon_cat, qdy_cat, qhr_cat, qmn_cat
   real, dimension (nq0) :: qsc_cat, qmag_cat
   character (len=100) :: evfile
   character (len=2) :: evtype
   character (len=1) :: qmagtype
   
   if (ievform > 3) then
        print *, 'Error! Unknown event list format! ievform = ' , ievform
        print *, 'Ending program'
        stop
    endif
    
    evtype = '  '
    qmagtype = ' '
   
   ! open file
   open (12, file=evfile, status='old')
   print *, 'First 10 events follow: '
   
   ! read, line by line 
   do i = 1, nq0
   
    if (ievform == 0) then   ! evlist format
      read (12, 8, end=22) idcusp(i), qyr_cat(i), qmon_cat(i), qdy_cat(i), qhr_cat(i), &
             qmn_cat(i), qsc_cat(i), evtype, qmag_cat(i), qmagtype, &
             qlat_cat(i), qlon_cat(i), qdep_cat(i)
8     format (10x,i9,x,i4,x,i2,x,i2,x,i2,x,i2,x,f6.3,x,a2,x,f5.2,x,a1,x, &
            f9.5,x,f10.5,x,f7.3)
 
    elseif (ievform == 1) then ! phase format
      read (12, *, end=22) qyr_cat(i), qmon_cat(i), qdy_cat(i), qhr_cat(i), &
         qmn_cat(i), qsc_cat(i), qlat_cat(i), qlon_cat(i), qdep_cat(i), qmag_cat(i), &
         EH_cat, EZ_cat, RMS_cat, idcusp(i)      
    
    elseif (ievform == 2) then ! growclust_cat format
     read (12, *, end=22) qyr_cat(i), qmon_cat(i), qdy_cat(i), qhr_cat(i), qmn_cat(i), &
      qsc_cat(i), idcusp(i), qlat_cat(i), qlon_cat(i), qdep_cat(i), qmag_cat(i)
!860   format (i4, 4i3, f7.3, i10, f10.5, f11.5, f8.3, f6.2)

    elseif (ievform == 3) then ! hypoinverse format
    ! example: 554488 2016 08 04 10 22 50.57 40.125500 -120.206333 12.07 0.04 1.00 1.13 0.29
    
    read (12, *, end=22) idcusp(i), qyr_cat(i), qmon_cat(i), qdy_cat(i), qhr_cat(i), & 
    qmn_cat(i), qsc_cat(i), qlat_cat(i), qlon_cat(i), qdep_cat(i), &
    RMS_cat, EH_cat, qmag_cat(i)
    
    endif
    
    ! print 1st 10 events
    if (i < 11) print 8, idcusp(i), qyr_cat(i), qmon_cat(i), qdy_cat(i), qhr_cat(i), &
             qmn_cat(i), qsc_cat(i), evtype, qmag_cat(i), qmagtype, &
             qlat_cat(i), qlon_cat(i), qdep_cat(i)
              
   enddo
   i = nq0 + 1
   print *, '***Warning: more events than array dimension'
22 nq = i - 1
  close(12)

  end
  
!-------------------------------------
!-----------------------------------------------------------------------------------------
!  **** READ_XCORDATA reads cross-correlation data from an input file and organizes into
!   the arrays used by the GrowClust main program. Xcor data and station lists can be
!   input in several different formats (specified by ixcor, istat). Output arrays are of
!   two sizes: npair for event-pair arrays and nk for phase (tdif) arrays.
!
!  Inputs:  irxform   =  xcor data file format: 0 for xcortobin 
!           istform   =  station list format: 0 for stlist.filtefs 
!           it12form  =  sign convention for tdif: 12 (tt1-tt2) or 21 (tt2-tt1)
!           xcorfile  =  xcor data file name
!           stfile    =  station list file name
!           npair0    =  maximum number of event pairs
!           ndif0     =  maximum number of differential times
!           nq        =  total number of events
!           idcusp    =  event IDs
!           qlat      =  event latitude
!           qlon      =  event longitude
!           rmincut   =  minimum rxcor value to keep differential time
!           rmin_ngood=  minimum rxcor value to be considered a "good" xcor result
!           delmax    =  maximum station distance to keep differential time
!           rpsavgmin =  minimum avg. rxcor for event pair to keep differential times
!           iponly    =  0 = use P and S differential times, iponly: 1 = use P only
!           ngoodmin  =  minimum number of differential times with rxcor > rmin_ngood to keep event pair
!           
!  Returns: npair     =  actual number of event pairs kept
!           nk        =  actual number of differential times kept
!           iqq1      =  array(len=npair) of serial ID #s (starting from 1) for event i in pair i,j
!           iqq2      =  array(len=npair) of serial ID #s (starting from 1) for event j in pair i,j
!           idcusp11  =  array(len=npair) of event ID #s (e.g. CUSPID) for event i in pair i,j
!           idcusp22  =  array(len=npair) of event ID #s (e.g. CUSPID) for event j in pair i,j
!           index1    =  array(len=npair) of starting indices of tdif observations for pair i,j in the tdif arrays
!           index2    =  array(len=npair) of ending indices of tdif observations for pair i,j in the tdif arrays
!           stname    =  array(len=nk) of station names for each tdif observation
!           ipp       =  array(len=nk) of phase (1=P, 2=S) for each tdif observation
!           tdif      =  array(len=nk) of differential times (tt_j-tt_i) for each tdif observation
!           rxcor     =  array(len=nk) of cross-corr. coefficient for each tdif observation
!           dist      =  array(len=nk) of average station distance for each tdif observation
!           slat      =  array(len=nk) of station latitudes for each tdif observation
!           slon      =  array(len=nk) of station longitudes for each tdif observation
!-----------------------------------------------------------------------------------------

   subroutine READ_XCORDATA(irxform, istform, it12form, xcorfile, stfile, npair0, ndif0,&
     nq, idcusp, qlat, qlon, rmincut, rmin_ngood, delmax, rpsavgmin, iponly, ngoodmin, & 
     npair, nk, iqq1, iqq2, idcusp11, idcusp22, index1, index2, &
     stname, ipp, tdif, rxcor, dist, slat, slon)
    
    implicit none
    
! -------- Spherical geometry parameters --------    
    integer, parameter     :: dp=kind(0.d0)                   ! double precision
    real(dp), parameter    :: degrad = (180.0_dp)/(3.1415927_dp)
    real(dp), parameter    :: degkm = 111.1949266_dp
!------------------------------------------------   
    
    integer :: irxform, istform, it12form, npair0, ndif0, npair, nk, iponly, ngoodmin
    integer :: ngoodminB, iponlyB, input_ok, npairB, nkB, ix1, ix2, ip, k, nq, is, ss
    integer ::  kk, k1, k2, qnum1, qnum2, qcusp1, qcusp2, j, ngood_pr, nkeep_pr
    integer :: iq1Bmax, iq2Bmax
    real :: rmincut, rmin_ngood, delmax, rpsavgmin
    real :: rminB, rmincutB, rmin_ngoodB, delmaxB, rpsavgminB, xcorminB
    real :: rps_pr, otc12
    real(dp) :: qlat0, qlon0, reflat0, cosreflat0, dx, dy
    integer, dimension(npair0) :: idcusp11, idcusp22, index1, index2, iqq1, iqq2
    integer, dimension(npair0) :: idcusp11B, idcusp22B, index1B, index2B, iqq1B, iqq2B
    integer, dimension(ndif0) :: ipp, ippB
    integer, dimension(nq) :: idcusp
    real, dimension(ndif0) :: rxcor, tdif, dist, rxcorB, tdifB, distB
    real(dp), dimension(nq) :: qlat, qlon
    real(dp), dimension(ndif0) :: slat, slon
    character (len=100) :: xcorfile, stfile, linebuf
    character (len=12), dimension(ndif0) :: stname, stnameB
    character (len=12) :: stname_kk
    character (len=1) :: ippchar
    integer :: npair_cut=0, nk_cut=0
    
    integer, dimension (1000) :: keepkk_pr, ipp_pr
    real, dimension (1000) :: tdif_pr, rxcor_pr, dist_pr
    real(dp), dimension (1000) :: slat_pr, slon_pr 
    character (len=12), dimension(1000) :: stname_pr
    
     
    ! check for invalid file types
    if (istform > 1) then
        print *, 'Error! Unknown station format! istat = ' , istform
        print *, 'Ending program'
        stop
    endif
    if (irxform > 1) then
        print *, 'Error! Unknown xcor data format! ixcor = ' , irxform
        print *, 'Ending program'
        stop
    endif
    
   if (irxform == 0) then !--------------- precomputed binary xcor format -------------
       
       print *, 'Reading binary xcorfile: ', xcorfile   
       open (14, file=xcorfile, form='unformatted', status='old')
       
       ! first read xcortobin run parameters (denoted w/B) ------------    
       read (14) rminB, rmin_ngoodB, delmaxB, rpsavgminB, iponlyB, ngoodminB, xcorminB
       rmincutB = min(rminB, xcorminB) ! the implied cutoff in the binary file is the smaller of the two
       print *, 'Binary file parameters:'
       print *, 'rmincut, rmin_ngood, delmax, rpsavgmin, iponly, ngoodmin = '
       write (*, '(f7.3, f7.3, 1x, f10.3, 2x, f7.3, 4x, i4, 4x, i4)') &
         rmincutB, rmin_ngoodB, delmaxB, rpsavgminB, iponlyB, ngoodminB   
       
       ! check to see if input parameters are consistent with binary file
            ! some tests check if desired run params are equal, others are one-sided
            ! (it's ok if rmin > rminB, but not if rmin < rminB)
       input_ok = 1
       !if ( rmincut < rmincutB  ) then
       if ( (rmincut .ne. 0) .and. (rmincut-rmincutB < -0.001) ) then
         print *, 'Error! Desired rmin inconsistent with binary xcor file (too low).'
         print *, 'Desired rmin, binary rmin, binary rmin = ', rmincut,rmincutB
         input_ok = 0
       endif      
       !if (delmax < delmaxB ) then
       if ( (delmax .ne. 0 ) .and. (delmax-delmaxB < -0.001)) then
          print *, 'Error! Desired delmax inconsistent with binary xcor file (too low).'
          print *, 'Desired delmax, binary delmax = ', delmax, delmaxB
          input_ok = 0
       endif
       !if (rpsavgmin < rpsavgminB) then
       if ((rpsavgmin .ne. 0) .and. (rpsavgmin-rpsavgminB < -0.001)) then
          print *, 'Error! Desired rpsavgmin inconsistent with binary xcor file (too low).'
          print *, 'Desired rpsavgmin, binary rpsavgmin = ', rpsavgmin, rpsavgminB
          input_ok = 0
       endif
       !if ( rmin_ngood < rmin_ngoodB ) then
       if ( (rmin_ngood .ne. 0) .and. (rmin_ngood-rmin_ngoodB < -0.001) ) then
         print *, 'Error! Desired rmin_ngood inconsistent with binary xcor file (too low).'
         print *, 'Desired rmin_ngood, binary rmin_ngood = ', rmin_ngood, rmin_ngoodB
         input_ok = 0
       endif 
       !if (ngoodmin < ngoodminB) then
       if ((ngoodmin .ne. 0) .and. (ngoodmin < ngoodminB)) then
          print *, 'Error! Desired ngoodmin inconsistent with binary xcor file (too low).'
          print *, 'Desired ngoodmin, binary ngoodmin = ', ngoodmin, ngoodminB
          input_ok = 0
       endif      
       if (iponly .ne. iponlyB) then
          print *, 'Error! Desired iponly inconsistent with binary xcor file (not equal).'
          print *, 'Desired iponly, binary iponly = ', iponly, iponlyB
          input_ok = 0
       endif
       
       ! stop program if input inconsistent with binary file
       if (input_ok < 1) then
         print *, 'Adjust input file or modify xcor file, ', trim(xcorfile), ' accordingly'
         print *, 'Program terminated'
         close(14)
         stop
       endif
       !---------------------------------
       
       ! everything ok: continue reading file
       read (14) npairB, nkB
       print *, 'npairB, nkB = ', npairB, nkB
       if (npairB > npair0 .or. nkB > ndif0) then
          print *, '***Error, exceeds: ', npair0, ndif0
          stop
       endif
       
       ! next, read event-pair arrays (length npair)
       read (14) iqq1B(1:npairB), iqq2B(1:npairB), idcusp11B(1:npairB), idcusp22B(1:npairB)
       read (14) index1B(1:npairB), index2B(1:npairB)
       
       ! check maximum iqq1 and iqq2
       iq1Bmax = maxval(iqq1B(1:npairB))
       iq2Bmax = maxval(iqq2B(1:npairB))
       if (iq1Bmax > nq) then
         print *, 'Error! max iqq1 from binary file > nQ'
         print *, 'max(iqq1), nq: ', iq1Bmax, nq
         print *, 'Ending program. Check event list.'
         close(14)
         stop
       elseif (iq2Bmax > nq) then
         print *, 'Error! max iqq2 from binary file > nQ'
         print *, 'max(iqq2), nq: ', iq2Bmax, nq
         print *, 'Ending program. Check event list.'
         close(14) 
         stop
       else
         print *, 'Maximum quake number in binary file: ', max(iq1Bmax, iq2Bmax)
       endif
   
       ! finally, read xcor arrays (length nk)
            ! note that index1, index2 (above) map events to min/max indices in xcor arrays
       read (14) stnameB(1:nkB), ippB(1:nkB), tdifB(1:nkB), rxcorB(1:nkB), distB(1:nkB)
       
       close (14)
       print *, 'Finished reading binary file: ' , xcorfile
       
       
       ! now loop over pairs, decide which observations to keep ---------------------
        ! initialized indices, npair, nk
        k1 = 0
        k2 = 0
        npair = 0
        nk = 0
        kk = 0
      npair_cut = 0
     nk_cut = 0
       print *, 'Selecting data from binary file...'
       
       print *, 'Selection parameters:'
       print *, 'rmin, rmin_ngood, delmax, rpsavgmin, iponly, ngoodmin = '
       write (*, '(f7.3, f7.3, 1x, f10.3, 2x, f7.3, 4x, i4, 4x, i4)') &
         rmincut, rmin_ngood, delmax, rpsavgmin, iponly, ngoodmin
       
       do ip = 1, npairB
         
         ! indices in phase arrays for this pair
         ix1 = index1B(ip)
         ix2 = index2B(ip)
         
         if ((ix2 .eq. 0) .or. (ix1 > ix2)) cycle ! check for empty pair
    
        ! set counters for tdif observations for this pair  
        kk = 0         ! counts number of differential times for this pair
        rps_pr = 0.0   ! stores rpsavg for this pair
        ngood_pr = 0   ! counts number of "good" observations for this pair
        nkeep_pr = 0   ! counts number of kept observations (may not quite be "good")
        
        ! loop over all tdif observations for this pair, check if
        !    (a) the tdif observation should be kept (phase ok and rxcor >= rmin)
        !    (b) the tdif observation is "good" (dist <= delmax and rxcor >= rmin_ngood)
        do k = ix1, ix2
          
          rps_pr = rps_pr + rxcorB(k)
          kk = kk + 1
          
          if ((iponly == 1) .and. (ippB(k) .ne. 1)) then ! must be P if iponly == 1
            keepkk_pr(kk) = 0
            nk_cut = nk_cut + 1
            
          elseif (rxcorB(k) < rmincut) then ! must have high enough rxcor value 
            keepkk_pr(kk) = 0
            nk_cut = nk_cut + 1
            
          else ! yep, keep it. Now test if it's "good"
            nkeep_pr = nkeep_pr + 1
            keepkk_pr(kk) = 1 
            if ((distB(k) <= delmax) .and. (rxcorB(k) >= rmin_ngood)) then
             ngood_pr = ngood_pr + 1
            endif
            
          endif
        enddo
        rps_pr = rps_pr/real(kk)  ! compute rpsavg for pair
        
         ! check to see if we should use this pair:
         !   keep at least 1, enough good observations, and rpsavg high enough  
        if ((nkeep_pr > 0).and.(ngood_pr >= ngoodmin).and.(rps_pr >= rpsavgmin)) then
           
               ! yep, keep pair
               npair = npair + 1
               iqq1(npair) = iqq1B(ip)
               iqq2(npair) = iqq2B(ip)
               idcusp11(npair) = idcusp11B(ip)
               idcusp22(npair) = idcusp22B(ip)
           
               ! index1(npair) gives the start index for this pair in the tdif arrays
               k1 = k2 + 1 
               index1(npair) = k1
          
              ! keep the "good" phase data for this pair
              kk = 0
              do k = ix1, ix2
                kk = kk + 1
                if (keepkk_pr(kk)>0) then
                  k2 = k2 + 1
                  rxcor(k2) = rxcorB(k)
                  tdif(k2) = tdifB(k) ! note: sign convention should be ok b/c read from xcorbin file
                  ipp(k2) = ippB(k)
                  dist(k2) = distB(k)
                  stname(k2)(1:12) = stnameB(k)(1:12)
                 endif
              enddo
          
              ! index2(npair) gives the end index for this pair in the tdif arrays
              index2(npair) = k2  
        else
          npair_cut = npair_cut+1 ! count number of cut pairs
        endif  ! endif on goodpair
          
            
    enddo ! enddo on event pairs -------------
    
    ! print numbers kept
    nk = k2
    print *, 'npair, nk = ', npair, nk
    print *, 'npair_cut, nk_cut = ', npair_cut, nk_cut 
       
    ! get station locations
      do k = 1, nk
        call GETSTAT_STLIST(stfile, stname(k), slat(k), slon(k), istform)
        if (slat(k) < -99.) then
         print *,'ERROR: station not found, stname = ', stname(k)
         print *, 'Ending program: xcor data and station list inconsistent'
         stop
         endif    
      enddo
      !print *, 'Finished reading station list to get slat, slon'

!-----------------------------------------------------------------------------------------
          
   elseif (irxform == 1) then !--------------- dt.cc format --------------------------------
   
    print *, 'Reading xcor-data text file: ', xcorfile  
     open (14, file=xcorfile, status='old')
     
     ! initialized indices, npair, nk
     k1 = 0
     k2 = 0
     npair = 0
     nk = 0
     kk = 0
     npair_cut = 0
     nk_cut = 0
     
     ! read file, line by line
     do j = 1, ndif0
      read(14, '(a100)', end=52) linebuf 

      if(linebuf(1:1) .eq. '#') then ! ------------ event pair line --------------
         
        ! first, process previous pair -------
        if (kk > 0) rps_pr = rps_pr/real(kk)  ! compute rpsavg for pair
        
         ! check to see if we should use this pair (enough good observations and rpsavg high enough)  
         if ((nkeep_pr > 0) .and. (ngood_pr >= ngoodmin) .and. (rps_pr >= rpsavgmin)) then
           
           ! yep, keep pair
           npair = npair + 1
           iqq1(npair) = qnum1
           iqq2(npair) = qnum2
           idcusp11(npair) = qcusp1
           idcusp22(npair) = qcusp2
           
           ! index1(npair) gives the start index for this pair in the tdif arrays
           k1 = k2 + 1 
           index1(npair) = k1
          
          ! keep the "good" phase data for this pair
          do k = 1, kk
            if (keepkk_pr(k)>0) then
              k2 = k2 + 1
              rxcor(k2) = rxcor_pr(k)
              tdif(k2) = tdif_pr(k) ! note: possible sign convention difference is fixed on read...
              ipp(k2) = ipp_pr(k)
              dist(k2) = dist_pr(k)
              slat(k2) = slat_pr(k)
              slon(k2) = slon_pr(k)
              stname(k2)(1:12) = stname_pr(k)(1:12)
             endif
          enddo
          
          ! index2(npair) gives the end index for this pair in the tdif arrays
          index2(npair) = k2  
        
        elseif (j > 1) then ! skip first entry
           npair_cut = npair_cut + 1 ! cut this pair
           !print *, 'Cutting pair: ', npair
        
        endif  ! endif on goodpair
        
        ! ok, move on to the new pair ---------------
        read(linebuf(2:100), *), qcusp1, qcusp2, otc12 ! (note that we don't need the origin
                                                       ! time correction b/c we don't mix xcor 
                                                       ! and catalog differential times)
        
        ! look up serial ID numbers (1:nq), given cuspID
        call LOOKUP_PAIR(qcusp1, qcusp2, nq, idcusp, qnum1, qnum2) 
        
        ! check for missing events...
        if ((qnum1 < 1) .or. (qnum2 < 1)) then
          print *, 'ERROR: event pair in xcor file missing from event list!'
          print *, 'ID numbers: ', qcusp1, qcusp2
          print *, 'Ending program: xcor data and event list inconsistent'
          close(14)
          stop
        endif
        
        ! get event pair centroid (needed to compute station distance)
        qlat0 = (qlat(qnum1) + qlat(qnum2))/2.0
        qlon0 = (qlon(qnum1) + qlon(qnum2))/2.0
        
        ! set counters for tdif observations for this pair  
        kk = 0         ! counts number of differential times for this pair
        rps_pr = 0.0   ! stores rpsavg for this pair
        ngood_pr = 0   ! counts number of "good" observations for this pair
        nkeep_pr = 0   ! counts number of "kept" observations for this pair
        
      else           ! --------------------- differential time line -----------------
        
        kk = kk + 1 ! increment nobs within pair
        
        ! get station name, tdif, rxcor, and phase character (P or S)
        stname_pr(kk)(1:12) = '            ' ! initialize name with blanks 
        read(linebuf(1:7), '(a7)') stname_kk(1:7) ! assumes 7 char, right-justified stname
        
        ! skip blanks at beginning
        do is = 1,7
          if (stname_kk(is:is) .ne. ' ') exit
        enddo
        
        ! copy over to stname starting with its 4th char
        do ss = is,7
         if (stname_kk(ss:ss) .eq. ' ') exit
         stname_pr(kk)(4+ss-is:4+ss-is) = stname_kk(ss:ss)
        enddo
        
        
        ! get station name, tdif, rxcor, and phase character (P or S)
            ! note, tdif likely has opposite sign convention (t1-t2) if formatted for hypoDD
        read(linebuf(ss:100),*), tdif_pr(kk), rxcor_pr(kk), ippchar
        if (it12form == 12) tdif_pr(kk) = -tdif_pr(kk) ! correct to 21 sign convention
        
        ! convert P/S to 1/2
        if (ippchar(1:1) == 'S' .or. ippchar(1:1) == 's') then 
            ipp_pr(kk) = 2
        else ! defaults to P-phase
            ipp_pr(kk) = 1
        endif
        
        ! ---- get station distance ------
        
        call GETSTAT_STLIST(stfile, stname_pr(kk), slat_pr(kk), slon_pr(kk), istform)
        
        ! make sure station exists....
        if (slat(kk) < -99.) then
         print *,'ERROR: station not found, stname = ', stname_pr(kk)
         print *, 'Ending program: xcor data and station list inconsistent'
         close(14)
         stop
        endif
    
        reflat0 = (slat_pr(kk)+ qlat0)/2.0 ! reference latitude
        !cosqlat = cos( 0.5*(qlat(iq1)+qlat(iq2)) / degrad )
        cosreflat0 = cos(reflat0/degrad) ! cosine of this angle, for dx calculation
        dx = (slon_pr(kk)-qlon0)*cosreflat0 ! x distance (deg)
        dy = (slat_pr(kk)-qlat0) ! y distance (deg)
        dist_pr(kk) = sqrt(dx**2 + dy**2)*degkm ! total distance to centroid (km)
        
        !---------------------------------
        
        ! now check to see if we should a) keep this observation, and if so, if b) it is a "good" observation
        !   a) keep criteria: P-phase if iponly == 1 AND r >= rmin 
        !   b) good criteria: keep criteria AND r >= rmin_ngood AND del <= delmax
        ! ----------
        rps_pr = rps_pr+rxcor_pr(kk) ! increment running sum of rps_pr (even if we don't keep...)
        
        if ((iponly == 1) .and. (ipp_pr(kk) .ne. 1)) then ! must be P if iponly == 1
           keepkk_pr(kk) = 0
           nk_cut = nk_cut + 1
           
        elseif (rxcor_pr(kk) < rmincut) then ! must have high enough rxcor value to keep
            keepkk_pr(kk) = 0
            nk_cut = nk_cut + 1   
        
        else ! yep, keep observation
           nkeep_pr = nkeep_pr + 1
           keepkk_pr(kk) = 1
           
           ! now test if it's "good"
           if ((rxcor_pr(kk) >= rmin_ngood) .and. (dist_pr(kk) <= delmax)) then ! good: must have rxcor > rmin and dist < delmax
            ngood_pr = ngood_pr + 1 ! good observation
           endif
           
        endif 
        
!       ! print out every 1000 or so pairs to check
!         if (mod(npair, 1000) == 0) then
!           if (kk==1) then
!             print *, 'kept pair # ', npair
!             print *, qcusp1, qcusp2, otc12
!           endif 
!            print *, stname_pr(kk)(1:7),  tdif_pr(kk), rxcor_pr(kk), &
!             dist_pr(kk), ipp_pr(kk)
!         endif
          
      endif !------------- close if statement for event/phase line --------
       
     enddo !------------------ end loop over file lines---------------------------------------------
     
     ! error check for too many lines...
     print *, 'ERROR reading xcorfile: ', xcorfile
     print *, 'TOO MANY LINES,', j
     print *, 'Ending program...'
     close(14)
     stop
        
52     close(14) ! ------- ok, file is read, close it ----------

        !  finish by processing last pair -------
        
         if (kk > 0) rps_pr = rps_pr/real(kk)  ! compute rpsavg for pair
        
         ! check to see if we should use this pair (enough good observations and rpsavg high enough  
          if ((ngood_pr >= ngoodmin) .and. (rps_pr >= rpsavgmin)) then
           
           ! yep, keep pair
           npair = npair + 1
           iqq1(npair) = qnum1
           iqq2(npair) = qnum2
           idcusp11(npair) = qcusp1
           idcusp22(npair) = qcusp2
           
           ! index1(npair) gives the start index for this pair in the tdif arrays
           k1 = k2 + 1 
           index1(npair) = k1
          
          ! keep the "good" phase data for this pair
          do k = 1, kk
            if (keepkk_pr(k)>0) then
              k2 = k2 + 1
              rxcor(k2) = rxcor_pr(k)
              tdif(k2) = tdif_pr(k) ! note possible sign convention difference is fixed on read...
              ipp(k2) = ipp_pr(k)
              dist(k2) = dist_pr(k)
              slat(k2) = slat_pr(k)
              slon(k2) = slon_pr(k)
              stname(k2)(1:12) = stname_pr(k)(1:12)
             endif
          enddo
          
          ! index2(npair) gives the end index for this pair in the tdif arrays
          index2(npair) = k2   
        
          else
             npair_cut = npair_cut + 1
             !print *, 'Cutting pair:', npair
        
          endif  ! endif on goodpair
          
       
     nk = k2
     print *, 'Finished processing xcor text file: ', xcorfile
     print *, 'npair, nk = ', npair, nk
     print *, 'npair_cut, nk_cut = ', npair_cut, nk_cut
   
   endif ! ------------ close of loop on the xcor file type
   
   print *, 'First two pairs follow:'
   do ip = 1,2
     write (*, 407) iqq1(ip), idcusp11(ip), iqq2(ip), idcusp22(ip) 
407    format (i8, 1x, i9,'   <----> ', i8, 1x, i9 )
     do k = index1(ip), index2(ip)
      write (*, 408) stname(k)(1:12), ipp(k), tdif(k), rxcor(k), dist(k)
408    format (a12, 1x, i2, 1x, f7.3, 1x, f7.3, 1x, f10.3)
     enddo
   enddo     
   
   !stop
   return
   end
   
   
!-----------------------------------------------------------------------
!  **** GETSTAT_STLIST gets station location from stlist file generated
!  from sactogfs_stp program.  It uses only first 7 characters of
!  12 character station name. (This is ok because we are unconcerned
!  with channels at this point).

!  Inputs:  stfile  =  name of station list (unsorted)
!           snam    =  char*12 station code
!           iform   =  station format (0 = evlist)

!  Returns: flat    =  station latitude (-999. if no match)
!           flon    =  station longitude
!-----------------------------------------------------------------------
!
      subroutine GETSTAT_STLIST(stfile,snam,flat,flon,iform)
      implicit none

      integer i, iform, j, k, ss
      integer nsta

      real*8 flat
      real*8 flon
      real*8 slat(5000)
      real*8 slon(5000)

      character*12 snam
      character*100 stfile
!      character*12 stname(5000)
      character*15 stname(5000)
      character*100 linebuf

      logical firstcall
      save firstcall,stname,slat,slon,nsta
      data firstcall/.true./
!
      if (firstcall) then
         firstcall=.false.
         open (19,file=stfile,status='old')
         do i=1,5000
         
            if (iform==0) then !iform=0 for stlist format
            
               stname(i)='               '
               read (19,13,end=12) stname(i)(4:8), &
                    slat(i),slon(i)
13             format (3x,a5,10x,f10.5,f12.5)

            elseif (iform==1) then ! stations.dat format
            
               stname(i)='               '
               read(19, '(a100)', end=12) linebuf
               
               ! check to see if in NW_STNAME format
               if (linebuf(3:3)== ' ' .and. linebuf(1:1) .ne. ' ') then
                 
                 ! if so, only copy over STNAME
                 do j = 4,8
                  if (linebuf(j:j) == ' ') exit
                  stname(i)(j:j) = linebuf(j:j)
                 enddo
                 
               else ! nope: STNAME only
               
                ! skip the blanks at the beginning
                do j = 1,12
                  if (linebuf(j:j) .ne. ' ') exit
                enddo
                k = j
                
                ! now copy over STNAME (output starting a char 4, until input blank or at 5 char)
                do j = k,k+4
                  if (linebuf(j:j) == ' ') exit
                  ss = 4+(j-k)
                  stname(i)(ss:ss) = linebuf(j:j)
                enddo
               endif
               
               ! now read latitude and longitude
               read(linebuf(j:100), *) slat(i), slon(i) 

            else
               read (19,11,end=12) stname(i),slat(i),slon(i)
11             format (a15,3x,f10.5,f12.5)                   
            endif
         enddo
         print *,'***Warning:  station file has more than 5000 lines'
         i=5001
12       nsta=i-1
         close (19)
         print *,'Station locations read.  Nsta = ', nsta 
         print *, '     (may include repeats for channels...)'
         do k = 1, nsta
            write(*, '(a12, f10.4, f10.4)') stname(k), slat(k), slon(k)
         enddo
      end if

      do 30 i=1,nsta
         if (snam(4:8) .eq. stname(i)(4:8)) then     !***only check station name part, max = 5 chars
            flat=slat(i)
            flon=slon(i)
            return
         endif
30    continue
      print *,'***GETSTAT_STLIST station not found ',snam

      flat = -999.
      flon = -999.
      return
      end
   
   !-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!  **** LOOKUP_PAIR looks up serial ID numbers (1:nq) for an event pair, given
! their cuspids and an array (length nq) of cuspids, in serial order 
!
!  Inputs:  qcusp1  =  event ID (not serial) of first event in pair
!           qcusp2  =  event ID (not serial) of first event in pair
!           nq      =  total number of events in idcusp array
!           idcusp  =  array of event IDs (not serial, e.g. cuspids)

!  Returns: qnum1   = serial number of event 1 (-1 means not found)
!           qnum2   = serial number of event 2 (-1 means not found)

    subroutine LOOKUP_PAIR(qcusp1, qcusp2, nq, idcusp, qnum1, qnum2)
    
    implicit none
    integer :: k, qcusp1, qcusp2, nq, qnum1, qnum2
    integer, dimension(nq) :: idcusp
    
    ! neither are found yet
    qnum1 = -1
    qnum2 = -1
    
    ! check each idcusp in sequence (could be speeded up...)
    do k = 1, nq
    
      if (idcusp(k) == qcusp1) qnum1 = k
      if (idcusp(k) == qcusp2) qnum2 = k
      
      if (qnum1 > 0 .and. qnum2 > 0) exit ! break loop once both found

    enddo
    
    return
    end   
      

   
   
   
      