!-----------------------------------------------------------------------
! Copyright 2017 Daniel Trugman
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

! ******** GROWCLUST: Main driver program to implement the GrowClust algorithm that uses 
! a hybrid, hierarchical clustering algorithm to perform relative relocation of  
! earthquake hypocenters based on waveform cross-correlation data. ***********************
!-----------------------------------------------------------------------------------------
! Daniel Trugman, February 2017
!-----------------------------------------------------------------------------------------
! object file and module dependencies:
!       input_subs.o (subroutines for handling input file I/0)
!         vel_subs.o (subroutines for handling velocity model and travel-time table creation)
!       stats_subs.o (subroutines for handling statistical computations, including resampling)
!    grow_params.mod (module for array size and auxiliary parameters)
! --> to compile, simply type make, or (assuming a gfortran compiler):
!   rm -f growclust *.o *.mod
!   gfortran -O -c grow_params.f90
!   gfortran -O -c input_subs.f90 -o input_subs.o
!   gfortran -O -c vel_subs.f90 -o vel_subs.o
!   gfortran -O -c stats_subs.f90 -o stats_subs.o
!   gfortran -O input_subs.o vel_subs.o stats_subs.o growclust.f90 -o growclust
!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
program growclust
   
   use grow_params ! auxiliary run parameters, including array sizes (adjust as necessary)
   
   implicit none
   
! -------- Spherical geometry parameters --------    
   integer, parameter     :: dp=kind(0.d0)                   ! double precision
   real(dp), parameter    :: degrad = (180.0_dp)/(3.1415927_dp)
   real(dp), parameter    :: degkm = 111.1949266_dp
!------------------------------------------------   
   
! generic integers
   integer :: i, iponly, iqmax, ip, ii, &
              k, kk, j, jj, npair, stlist_fmt, &
              npick, nq, ntree, kraw, &
              nconnect, ib, ngoodmin, ntot, &
              nb, iq, itree, n, nk, kraw8, ip8, i8, j8, &
              npr8, npk8, k8, itreeold, nbranch_i, nbranch_j, idcusp_i, idcusp_j, &
              nbranch_min, iflag, evlist_fmt, xcordat_fmt, tdif_fmt, ierror
              

! generic single-precision reals                 
   real :: rmin, delmax, rgoodavg, qdep0, &
           qdep1, qdep2, torgdif, rms, rmed, resol, frac,  &
           qdep_off, rpsavgmin, delkm, distance, &
           qdep_off1, qdep_off2, qtim_off1, qtim_off2, tavg, &
           cdist, cdep1, cdep2, rmsmax, tsec1, tsec2, rmincut, &
           min_qdep, max_qdep

! generic double-precision reals           
   real(dp) :: qlat0, qlon0, qlat1, qlon1, qlat2, qlon2, qlat_off, qlon_off, &
                    qlat_off1, qlat_off2, qlon_off1, qlon_off2, &
                    clat1, clat2, clon1, clon2, dx, dy, dz, cosqlat 

! file names   
   character(len=100) :: infile_evlist, infile_xcordat, infile_ctl, infile_stlist, linebuf
   character(len=100) :: outfile_cat, outfile_clust, outfile_log, outfile_boot
     
! variables to store xcor file data
   integer, dimension(npair0) :: idcusp11, idcusp22, index1, index2, iqq1, iqq2
   integer, dimension(ndif0) :: ipp
   real, dimension(ndif0) :: rxcor, tdif, dist
   real(dp), dimension(ndif0) :: slat, slon
   character (len=12), dimension(ndif0) :: stname
   integer, dimension(10) :: ip810
   integer, dimension(npair0) :: ngood, indxr 
   real, dimension(npair0) :: rfactor
   real, dimension(n0) :: resid
   
!variables for clustering
   integer, dimension(nbmax) :: inbranch
   integer, dimension(nq0) :: index
   integer, dimension(ntmax) :: nbranch=0, nbranch2=0
   real, dimension(ntmax) :: tdep, torg=0., tdep00, torg00
   real(dp), dimension(ntmax) :: tlat, tlon, tlat00, tlon00 
   real, dimension(nq0) :: qdep, qtim,  qdep_cat, qtim_cat
   real(dp), dimension(nq0) :: qlat, qlon, qlat_cat, qlon_cat
   integer, dimension(nq0) :: qyr_cat, qmon_cat, qdy_cat, qhr_cat, qmn_cat, idcusp
   real, dimension (nq0) :: qsc_cat, qmag_cat
   
! variables to go into DIFCLUST
   integer, dimension(n08) :: ipp8
   real, dimension(n08) :: tdif8,  qdep18, qtim18, qdep28, qtim28, resid8
   real(dp), dimension(n08) :: slat8, slon8, qlat18, qlon18, qlat28, qlon28

! variables for to store output for first growclust estimate (no resampling/bootstrap)
   real(dp), dimension(nq0) :: qlatR, qlonR ! relocated lat/lon (no bootstrap)
   real(dp), dimension(nq0) :: qcxR, qcyR ! x,y coordinates w.r.t cluster centroid
   real(dp), dimension(nq0) :: qczR ! z coordinates w.r.t. cluster centroid
   real, dimension (nq0) :: qdepR, qtimR ! relocated dep, time, indexR
   integer, dimension(nq0) :: nbranchqR, itreeqR ! number of events/branch, cluster/tree ID
   integer, dimension(ntmax) :: nbranchtR=0 ! number of events/branches for each tree
   real, dimension(ntmax) :: tdepR, torgR=0. ! tree centroid: dep/org (no bootstrap)
   real(dp), dimension(ntmax) :: tlatR, tlonR ! tree centroid: relocated lat/lon (no bootstrap)
   integer :: ntreeR ! number of trees (no bootstrap)
   
! variables to store copies of binary file phase information for later resampling
   integer, dimension(npair0) :: idcusp1100, idcusp2200
   integer, dimension(npair0) :: index100, index200, iqq100, iqq200
   integer, dimension(ndif0) :: ipp00
   real, dimension(ndif0) :: rxcor00, tdif00, dist00
   real(dp), dimension(ndif0) :: slat00, slon00

! variables for resampling
    integer :: nboot, nitb, ix1, ix2, ix1n
    real(dp) :: SEx, SEy, SEt
    integer, dimension (ndif0) :: samp_vec
    real(dp), dimension (nq0) :: SEh_boot, SEz_boot, MADh_boot, MADz_boot
    real(dp), dimension (nq0) :: SEt_boot, MADt_boot
    real(dp), dimension (nq0) :: BlatMEAN, BlonMEAN, BdepMEAN
    real(dp), dimension (nq0) :: BtimMEAN
    real(dp), dimension (maxboot,nq0) :: Blon, Blat, Bdep
    real(dp), dimension (maxboot,nq0) :: Btim
    real(dp) :: MADx, MADy, MADz, MADt
    integer, dimension (maxboot, nq0) :: Bnb, Bitreeq
    real, dimension (nq0) :: BnbMEAN
    integer, dimension (nq0) :: BnbMAX, BnbMIN, BnbISO

! variables for RMS
    real, dimension (nq0) :: qrmsP, qrmsS, qrms
    integer, dimension (nq0) :: qndiffP, qndiffS, qnpair
    real, dimension (ndif0) :: tdif_pred
    real :: prmsP, prmsS
    integer :: iq1, iq2, kcount, npickgoodP, npickgoodS, itree1, itree2
    integer, dimension (ntmax) :: nb_indxr

! variables for log file output 
    real :: CrmsP, CrmsS, resPmean, resSmean, resPmed, resSmed
    integer :: nreloc, ntree2, ntree5, ntree10, ntree20, ntree50, ntree100
    integer :: npairC, npC, nsC, jjP, jjS
    real, dimension (ndif0) :: resP, resS
    
! variables for velocity model and ttables
    character(len=100) :: infile_vzmdl, vzfinefile
    character(len=100), dimension(2) :: TTtabfile
    real :: vpvs_factor, vz_dz, plongcutP, plongcutS, vzmin, vzmax
    real :: tt_dep0, tt_dep1, tt_ddep, tt_del0, tt_del1, tt_ddel
    integer :: vzinterp_type, rayparam_min

!---------------------------------------------------------------------------!
   
!---------------------------------------------------------------------------!
   
!    !----------- Get initial inputs ------------------------!
   print *, ' '
   !read (*, '(a)') infile ! reads from command line
   call getarg(1, infile_ctl) ! reads argument directly
   open (11, file=infile_ctl, status='old')
   print *, 'Reading input control file: ', infile_ctl
   print *, '--------------------------------------------------'
     
     
     !------- event list format, file name
     call GET_INP(11, linebuf)
     read(linebuf, *) evlist_fmt
     call GET_INP(11, linebuf)
     infile_evlist(1:100) = linebuf(1:100)
     print *, 'Input event file format, name  : '
     print *, evlist_fmt, ', ', trim(infile_evlist)
     
     ! ------ station list format, file name
     call GET_INP(11, linebuf)
     read (linebuf, *) stlist_fmt
     call GET_INP(11, linebuf)
     infile_stlist(1:100) = linebuf(1:100)
     print *, 'Input station file format, name: '
     print *,  stlist_fmt, ', ', trim(infile_stlist)
     
     ! ----- xcor data format, file name
     call GET_INP(11, linebuf)
     read (linebuf, *) xcordat_fmt, tdif_fmt
     call GET_INP(11, linebuf) 
     infile_xcordat(1:100) = linebuf(1:100)
     print *, 'Xcor data file format, name    :  '
     print *, xcordat_fmt, ', ', trim(infile_xcordat)
     if ((tdif_fmt .eq. 12) .or. (tdif_fmt .eq. 21)) then 
        print *, 'tdif sign convention: ', tdif_fmt
     else
        print *, 'Error! Unknown tdif sign convention.'
        print *, '   12: event 1 - event 2'
        print *, '   21: event 2 - event 1'
        print *, 'Program terminated. Please fix input file: ', infile_ctl
        close(11)
        stop
     endif 
     
     ! ----- velocity model and travel time files
     call GET_INP(11, linebuf)
     infile_vzmdl(1:100) = linebuf(1:100)
     print *, 'Name of input velocity model file: '
     print *, trim(infile_vzmdl)
     call GET_INP(11, linebuf)
     vzfinefile(1:100) = linebuf(1:100)
     call GET_INP(11, linebuf)
     TTtabfile(1)(1:100) = linebuf(1:100)
     call GET_INP(11, linebuf)
     TTtabfile(2)(1:100) = linebuf(1:100)
     print *, 'Name of output travel-time tables: '
     print *, trim(TTtabfile(1)), ' ,  ', trim(TTtabfile(2))
     
     
     !----  velocity model specifications
        ! vpvs_factor: controls Vp/Vs ratio if input velocity model is 2 column
        ! rayparam_min allows for control min ray param plongcut: (-1) use default values, (2) read from input file 
        !  - note if plongcut < u_moho, first arrivals will be a diffracted phase)
        !       e.g., Pn,Sn (diffracted) will arrive before Pg,Sg (direct) beyond the crossover distance
     call GET_INP(11, linebuf)
     read(linebuf,*) vpvs_factor, rayparam_min
     if (vpvs_factor < 0.01) vpvs_factor = sqrt(3.0) ! default Vp/Vs for Poisson solid
     if (rayparam_min < 0) then    ! use default values  
        plongcutP = .133    ! assumes a sub-moho Vp of 7.5 km/s (.133 = 1.0/7.5)
        plongcutS = .238    ! assumes a sub-moho Vp of 7.5 km/s (.238 = sqrt(3)/7.5)
     else                   ! read from input 
        plongcutP = rayparam_min ! p-phase
        plongcutS = rayparam_min*vpvs_factor ! s-phase
     endif
     print *, 'Assumed Vp/Vs if Vs not specified in input file: '
     write(*, '(f8.4)') vpvs_factor 
     print *,'Min ray params (p-P, p-S) at long range (.133 = no Pn, .238 = no Sn):'
     write(*, '(f8.4,1x,f8.4)'), plongcutP, plongcutS
     
     !------- travel-time table size parameters
        ! tt_dep0, tt_dep1, tt_ddep: min dep, max dep, dep spacing
        ! tt_del0, tt_del1, tt_ddel: min del, max del, del spacing
        ! vz model interpolation spacing vz_dz is set to tt_ddep
     call GET_INP(11, linebuf)
     read(linebuf,*) tt_dep0, tt_dep1, tt_ddep
     call GET_INP(11, linebuf)
     read(linebuf,*) tt_del0, tt_del1, tt_ddel
     print *, 'Travel-time table depths (km):  min, max, space' 
     write (*, '(f8.3, f8.3, f8.3)') tt_dep0, tt_dep1, tt_ddep
     print *, 'TT table distances (del):       min, max, space' 
     write (*, '(f8.3, f8.3, f8.3)') tt_del0, tt_del1, tt_ddel
     vz_dz = tt_ddep ! set interpolation distance

     ! ------- GrowClust run parameters
        ! rmin        =  minimum xcor coefficient to count in computation of GrowClust rfactor similarity coefficient
        ! delmax      =  maximum station distance to count in computation of GrowClust rfactor similarity coefficient
        ! rmsmax      =  maximum rms diff. time residual for a cluster merger to be allowed  
        ! rpsavg      =  minimum desired average xcor coefficient to keep event pair  (use 0 to keep all in input file) 
        ! rmincut     =  minimum desired xcor coefficient to keep tdif observation  (use 0 to keep all in input file)
        ! ngoodmin    =  minimum number of xcor obs. with r>rmin for a pair to be label "good"  (use 0 to keep all in input file) 
        ! iponly      =  (0) use both P and S, (1): keep only P-phase diff times (use 0 to keep all in input file)
        ! nboot       =  number of bootstrap iterations desired 
        ! nbranch_min =  minimum nbranch (number of events/cluster) to be output in cluster file 
     !------  
     call GET_INP(11, linebuf)
     read(linebuf,*) rmin, delmax, rmsmax
      print *, 'rmin, delmax, rmsmax'
     write (*, '(f6.3, f9.3, f6.3)'), rmin, delmax, rmsmax
     !------
     call GET_INP(11, linebuf)
     read(linebuf,*) rpsavgmin, rmincut, ngoodmin, iponly
     print *, 'rpsavgmin, rmincut, ngoodmin, iponly'
     write (*, '(f10.3, f10.3, i8, i8)'), rpsavgmin, rmincut, ngoodmin, iponly
     !-----
     call GET_INP(11, linebuf)
     read(linebuf,*) nboot, nbranch_min
     print *, 'nboot, nbranch_min'
     write (*, '(i7, i7)'), nboot, nbranch_min
     !-----
      
     ! ------- output file names
     call GET_INP(11, linebuf)
     outfile_cat(1:100) = linebuf(1:100)
     call GET_INP(11, linebuf)
     outfile_clust(1:100) = linebuf(1:100)
     call GET_INP(11, linebuf)
     outfile_log(1:100) = linebuf(1:100)
     call GET_INP(11, linebuf)
     outfile_boot(1:100) = linebuf(1:100)
     print *, 'Output catalog file: ', trim(outfile_cat)   
    
    close(11)

    print *, 'Input control file read.'
    
! call subroutine to check inputs for errors
    call INPUT_CHECK(plongcutP, plongcutS, vpvs_factor, tt_dep0, tt_dep1, tt_ddep, &
      tt_del0, tt_del1, tt_ddel, rmsmax, delmax, iponly, nboot, maxboot)

! call subroutine to check grow_params module for errors
    call PARAM_CHECK(conparam, hshiftmax, vshiftmax, rmedmax, boxwid, nit, samp_type, irelonorm) 
    
     print *, '--------------------------------------------------'
     print *, ' '

!-----------------------------------------------------------------------------------------    
    
! READ_infile_evlist: Read event list to get starting (catalog) locations --------------------
        ! -- accepts different event list formats specified by xcordat_fmt, stlist_fmt 
    print *, ' '
    print *, 'Reading event list: ' 
    call READ_EVFILE(evlist_fmt, infile_evlist, nq0,idcusp, qyr_cat, qmon_cat, qdy_cat,  &
    qhr_cat, qmn_cat, qsc_cat, qmag_cat, qlat_cat, qlon_cat, qdep_cat, nq, min_qdep, max_qdep)
    
    ! initialize relocated positions and clustertree arrays (1 event/cluster) -----
    do i = 1, nq   
     
     ! initialize relocated positions: qlat, qlon, qdep, qtim 
      qtim_cat(i) = 0.0
      qlat(i) = qlat_cat(i)
      qlon(i) = qlon_cat(i)
      qdep(i) = qdep_cat(i)
      qtim(i) = qtim_cat(i)
      
      ! initialize clustertree arrays
      index(i) = i             !set initial cluster number to event number
      nbranch(i) = 1           !one event per cluster to start
      tlat(i) = qlat(i)        !set cluster centroid to event location
      tlon(i) = qlon(i)
      tdep(i) = qdep(i)
      torg(i) = qtim(i)            
   enddo
   ntree = nq                 !start out with number of trees equal to number of quakes
!---------------------------

    ! READ_XCORDATA: Reads xcor data and associated station locations
    ! -- accepts different xcor, station formats using xcordat_fmt, stlist_fmt
    ! -- output event pair arrays: iqq1, iqq2, idcusp11, idcusp22, index1, index2
    ! -- output tdif/phase arrays: stname, ipp, tdif, rxcor, dist, slat, slon
    print *, ' '
    print *, 'Reading xcor data and associated station list: '
    call READ_XCORDATA(xcordat_fmt, stlist_fmt, tdif_fmt, infile_xcordat, infile_stlist, npair0, ndif0, &
     nq, idcusp, qlat, qlon, rmincut, rmin, delmax, rpsavgmin, iponly, ngoodmin, & 
     npair, nk, iqq1, iqq2, idcusp11, idcusp22, index1, index2, &
     stname, ipp, tdif, rxcor, dist, slat, slon) 
     
  ! for robustness, final check to make sure nq, npair, nk are not too large (edit 11/2016)
  ierror = 0
  if (nq > nq0) then
      print *, 'Input error: too many events! Increase nq0 in grow_params.f90...'
      print *, 'nq, nq0 = ', nq, nq0
      ierror = 1
  endif
  if (npair > npair0) then
      print *, 'Input error: too many event pairs! Increase npair0 in grow_params.f90...'
      print *, 'npair, npair0 = ', npair, npair0
      ierror = 1
  endif
  if (nk > ndif0) then
      print *, 'Input error: too many differential times! Increase ndif0 in grow_params.f90...'
      print *, 'ndif, ndif0 = ', nk, ndif0 
      ierror = 1
  endif
  if (ierror == 1) stop ! stop on this error
  
    

  ! Save copies of in event pair, tdif arrays for later resampling -----------------!
      
   ! save copy of input event pair arrays: iqq1, iqq2, idcusp11, idcusp22, index1, index2
   do ip = 1, npair
        iqq100(ip) = iqq1(ip)
        iqq200(ip) = iqq2(ip)
        idcusp1100(ip) = idcusp11(ip)
        idcusp2200(ip) = idcusp22(ip)
        index100(ip) = index1(ip)
        index200(ip) = index2(ip)
   enddo
   
   ! save copy of input tdif/phase arrays: ipp, slat, slon, dist, tdif, rxcor
   do k = 1, nk
      ipp00(k) = ipp(k)
      slat00(k) = slat(k)
      slon00(k) = slon(k)
      tdif00(k) = tdif(k)
      rxcor00(k) = rxcor(k)
      dist00(k) = dist(k)
      
      ! quality control checks on xcor data
      if (dist(k) > tt_del1) then
      	  print *, 'ERROR: STATION DISTANCE', dist(k), '> tt_del1', tt_del1
      	  print *, 'Check input files (evlist, stlist, xcordata) for errors or '
          print *, 'increase tt_del1 in GrowClust control (.inp) file.'
      	  ierror = 1
      endif
      if (tdif(k) > tdifmax) then
          print *, 'ERROR: DIFFERENTIAL TIME', tdif(k), '> tdifmax', tdifmax
          print *, 'Check differential time data for errors or '
          print *, 'increase tdifmax parameter in grow_params.f90 and recompile.'
      	  ierror = 1
      endif
        
   enddo
   if (ierror == 1) stop ! stop on this error 
   
   !-----------------------------------------------------------!
   
   ! make travel time tables
   print *, 'Making travel-time tables...'
   
   print *, 'Interpolating velocity model: ', trim(infile_vzmdl), '->', trim(vzfinefile)
   vzmin = tt_dep0
   vzmax = tt_dep1
   call vzfillin(infile_vzmdl, vzfinefile, vpvs_factor, vz_dz, vzmin, vzmax)
   print *, ' '
   
   print *, 'Creating P-phase table: ', trim(TTtabfile(1))
   call deptable(vzfinefile, 1, plongcutP, tt_dep0,tt_dep1, tt_ddep, &
     tt_del0, tt_del1, tt_ddel, TTtabfile(1)(1:100))
   print *, ' '
     
   print *, 'Creating S-phase table: ', trim(TTtabfile(2))
   call deptable(vzfinefile, 2, plongcutS, tt_dep0,tt_dep1, tt_ddep, &
     tt_del0, tt_del1, tt_ddel, TTtabfile(2)(1:100))
     
    ! Added 2/2017 to check event depths vs velocity model, travel time tables
    print '(a26, f8.2, f8.2)', 'min, max event depth:', min_qdep, max_qdep
    print '(a26, f8.2, f8.2)', 'min, max table depth:', tt_dep0, tt_dep1
    print '(a26, f8.2, f8.2)', 'min, max vzmodel depth:', vzmin, vzmax 
    if (min_qdep < tt_dep0) then
    	print *, 'WARNING: min event depth < min table depth'
    	!ierror = 1 ! allow this, but warn user (should be ok if depth is near 0)
    endif
    if (min_qdep < vzmin)  then
    	print *, 'WARNING: min event depth < min vzmodel depth'
    	!ierror = 1 ! allow this, but warn user (should be ok if depth is near 0)
    endif
    if (max_qdep > tt_dep1) then ! note tt_dep1 is >= vzmax
    	print *, 'ERROR: max event depth > max table / velocity model depth'
    	ierror = 1 ! don't allow this
    endif
    if (tt_dep0 < vzmin) then ! for robustness, check this as well
    	print *, 'ERROR: min table depth < min vzmodel depth'
    	ierror = 1
    endif
    if (ierror == 1) stop ! stop on this error
    

   print *, '--------------------------------------------------'
   print *, ' '
   print *, ' '

  !-----------------------------------------------------------!   
   
   !------- compute event pair quality (similarity) factors ------!
   iqmax = 0 
   ! loop over pairs   
   do ip = 1, npair 
      ngood(ip) = 0
      rfactor(ip) = -1.                          !***stays at -1 if 'bad' pair      
      
      if (index2(ip)-index1(ip) <= 1) cycle         !****added to make more robust
      if (index1(ip) <= 1) cycle
      if (index2(ip) <= 1) cycle     
            
      rgoodavg = 0.
      
      ! sum over all observations for a given pair that
      !     (1) have rxcor > rmin and (2) are observed at near enough station
      do i = index1(ip), index2(ip)
         if (dist(i) <= delmax .and. rxcor(i) >= rmin) then
            ngood(ip) = ngood(ip) + 1         
            rgoodavg = rgoodavg + rxcor(i)
         endif
      enddo
      if (ngood(ip) > 0) rgoodavg = rgoodavg/real(ngood(ip)) ! sum --> average
      
      
      !rfactor: measure of 'quality' of event pair:
      !     rfactor = (# of "good" rxcor values) x (avg of "good rxcor values) 
      rfactor(ip) = real(ngood(ip))*rgoodavg                   
      
      if (iqq1(ip) > iqmax) iqmax = iqq1(ip)
      if (iqq2(ip) > iqmax) iqmax = iqq2(ip)
   enddo
   print *, 'Maximum quake number = ', iqmax

   
! now sort by rfactor (indxr has the highest quality / most similar pair) last
    ! i.e. ip = indxr(npair) gives the pair number of the best pair
   call INDEXX(npair, rfactor, indxr)   

! print best pair   
   print *, 'best pair follows'
   ip = indxr(npair) ! index number of best pair
   print *, k, iqq1(ip), iqq2(ip), rfactor(ip)
   do k = index1(ip), index2(ip)
      print '(a12,f8.3, f9.3,i3,f7.3,f6.2)', stname(k), slat(k), slon(k), ipp(k), tdif(k), rxcor(k)
   enddo
  
  !--------------------------------------------------------------------------! 
! ----- open log file, print starting parameters   
   open (16, file=outfile_log)

   
   ! write log file header (NEW)
   write(16, *) '************************ Input files ************************'
   write(16, *) '     control file:   ', trim(infile_ctl)
   write(16, *) '       event list:   ', trim(infile_evlist)
   write(16, *) '     station list:   ', trim(infile_stlist)
   write(16, *) '     xcordat file:   ', trim(infile_xcordat)
   write(16, *) '     velocity mdl:   ', trim(infile_vzmdl)
   write(16, *) '******************** Travel Time Tables *********************'
   write(16, *) '          P-phase:   ', trim(TTtabfile(1))
   write(16, *) '          S-phase:   ', trim(TTtabfile(2))
   write(16, *) '*********************** Output files ************************'
   write(16, *) '     catalog file:   ', trim(outfile_cat)
   write(16, *) '     cluster file:   ', trim(outfile_clust)
   write(16, *) '         log file:   ', trim(outfile_log)
   write(16, *) '   bootstrap file:   ', trim(outfile_boot)
   write(16, *) '**************** GROWCLUST Run Parameters *******************'
   write(16, '(a56, f6.2)') ' (min rxcor value for evpair similarity coeff.): rmin =', rmin
   write(16, '(a56, f6.1)') ' (max sta. dist for evpair similarity coeff.): delmax =', delmax
   write(16, '(a56, f6.2)') ' (max rms residual to join clusters): rmsmax =', rmsmax
   write (16, '(a56, i6)' ) ' (num. bootstrap uncertainty iterations): nboot =', nboot
   write(16, *) '**************** Auxiliary Run Parameters *******************'
   write(16, '(a56, f6.4)') ' min connnection fraction to join clusters: ', conparam 
   write(16, '(a56, f6.2)') ' max catalog dist to join clusters: ', distmax
   write(16, '(a56, f6.2)') ' max relocated dist to join clusters: ', distmax2
   write(16, '(a56, i6)') ' min number in cluster to apply shift test: ', nclustshiftmin
   write(16, '(a56, f6.2)') ' max permitted horzizontal cluster shifts: ', hshiftmax
   write(16, '(a56, f6.2)') ' max permitted vertical cluster shifts: ', vshiftmax
   write(16, '(a56, f6.2)') ' max median absolute residual to join clusters: ', rmedmax  
   write(16, *) '*************************************************************'
   write(16, *) ' '
   write(16, *) '*************************************   Cluster GrowLog   ****************************************'
   write(16, *) ' '
   write(16, *) '  CID1   CID2    NB1   NB2   CID    NB NPR NDT', & 
   '   dLATC     dLONC     dDEPC     dDISTC   RMSC  RMEDC'
   
! --------- let's relocate best pair for test --------------------!
   i = iqq1(ip)
   j = iqq2(ip)
   idcusp_i = idcusp11(ip)
   idcusp_j = idcusp22(ip)
   print *, 'quake numbers for best pair = ', i, j
   print *, 'cuspids for best pair = ', idcusp_i, idcusp_j
   print *, 'cuspids from evlist = ', idcusp(i), idcusp(j)
   if (idcusp_i /= idcusp(i) .or. idcusp_j /= idcusp(j) ) then
      print *, '***Error, index problem in evlist.  Likely not same file as used for xcor calculation.'
      write (16, *) '***Error, index problem in evlist.  Likely not same file as used for xcor calculation.'
      close(16)
      stop
   endif
   
   ! DT edit: new version that doesn't call GET_LOC (doesn't sort by CUSPID if input isn't...)
   qlat1 = qlat(i)
   qlon1 = qlon(i)
   qdep1 = qdep(i)
   print *, 'evlist loc = ', qlat_cat(i), qlon_cat(i), qdep_cat(i)
   qlat2 = qlat(j)
   qlon2 = qlon(j)
   qdep2 = qdep(j)
   print *, 'evlist loc = ', qlat_cat(j), qlon_cat(j), qdep_cat(j)
   
   ! compute centroid (mean) location   
   qlat0 = (qlat1 + qlat2)/2.
   qlon0 = (qlon1 + qlon2)/2.
   qdep0 = (qdep1 + qdep2)/2.
   print *, 'mean loc = ', qlat0, qlon0, qdep0
   
   npick = index2(ip) - index1(ip) + 1
   if (npick > n0) then
      print *, '***Error, npick too big:', npick, n0
      stop
   endif
   
   ! DIFLOC performs relative event location for a pair of events
   k = index1(ip)
   call DIFLOC(qlat0, qlon0, qdep0, npick, tdif(k), ipp(k), slat(k), slon(k), TTtabfile, boxwid, nit, &
         irelonorm, qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, torgdif, resid, rms, rmed, resol)
   print *, 'DIFLOC results: '
   print *, 'quake 1 loc = ', qlat1, qlon1, qdep1
   print *, 'quake 2 loc = ', qlat2, qlon2, qdep2   
   print *, 'torgdif, rms, rmed, resol = ', torgdif, rms, rmed, resol
   print *, '--------------------------------------------------'
   print *, ' '
   print *, ' '
   

! ------------------------------------------------------------------------------------ 

! ------------------------------------------------------------------------------------
   
! ------    now start relocations   ----------------
! assumes we have read in evlist file and have nq, qlat_cat, qlat, etc., initial locations and catalog locations
!        we also have ngood, which give number of 'good' xcors for each pair
!  notes on key variables/arrays
!     - index gives the current tree (cluster) # for each event (initialized to event #)
!     - nbranch counts number of branches in each tree
!     - tlat,tlon,tdep store centroid location for each tree
!     - current locations go into qlat, qlon, qdep (originally set to catalog locations)
!     - original locations remain in qlat_cat, qlon_cat, qdep_cat
   
  
   print *, "Starting GROWCLUST algorithm"
   call TIMER
   do nitb = 0, nboot ! BOOTSTRAP ITERATION LOOP (0 = unresampled estimate) ------------------------
   
   
   if (nitb > 0) then ! bootstrap iteration -----------------
   
       ! reset event and cluster stuff for iteration
       do i = 1, nq
          qlat(i) = qlat_cat(i)     ! set initial locations to catalog locations
          qlon(i) = qlon_cat(i)
          qdep(i) = qdep_cat(i)
          qtim(i) = qtim_cat(i)
          index(i) = i             !set initial cluster number to event number
          nbranch(i) = 1           !one event per cluster to start
          tlat(i) = qlat(i)        !set cluster centroid to event location
          tlon(i) = qlon(i)
          tdep(i) = qdep(i)
          torg(i) = qtim(i)
          
       enddo
       ntree = nq                 !start out with number of trees equal to number of quakes
   
   
        ! reset  iqq1, iqq2, idcusp11, idcusp22, index1, index2 (modified on earlier output)
       do ip = 1, npair
            iqq1(ip) = iqq100(ip)
            iqq2(ip) = iqq200(ip)
            idcusp11(ip) = idcusp1100(ip)
            idcusp22(ip) = idcusp2200(ip)
            if (samp_type == 1) index1(ip) = index100(ip)
            if (samp_type == 1) index2(ip) = index200(ip)
       enddo
   
       ! get resampling vector
            ! samp_type == 1: resample with replacement from each event pair --> index1 and index2
            !    stay the same across each iteration
            ! samp_type == 2 (preferred): resample entire phase vector (dtt, dxcor, ...) --> index1 and
            !    index2 can vary for each iteration, which may be a plus, because then ngood can vary realistically
        !--------
        !  version 1: resample each pair independently
        if (samp_type == 1) then
            do ip = 1, npair
                ix1 = index1(ip)
                ix2 = index2(ip)
                if (ix2 > 0 .and. ix1 > 0) then
                  call GET_SAMPLEVEC(samp_vec(ix1:ix2), ix2-ix1+1, iseed)
                endif
            enddo
        !  version 2: resample entire nk-sized data vectors at once, then recompute index1 and index2
        elseif (samp_type == 2) then
        
            call GET_SAMPLEVEC(samp_vec(1:nk), nk, iseed) ! compute the sampling vector
        
            ! loop over event pairs, and count number of samples for each pair --> updating index1, index2
            ix1n = 1 
            do ip = 1, npair
        
                ! get index1 and index2 for this pair from unresampled data
                ix1 = index100(ip)
                ix2 = index200(ip)
            
                ! count number of data-points in resampled data vector for this pair
                kcount = 0
                if (ix2 > 0 .and. ix1 > 0) then
                  do k = ix1, ix2
                    kcount = kcount+samp_vec(k)
                  enddo
                endif
            
                ! new index1
                index1(ip) = ix1n
            
                ! new index2 
                if (kcount == 0) then
                    index2(ip) = 0 ! set to zero
                else
                   index2(ip) = index1(ip) + kcount - 1
                endif
            
                ! update ix1n for next pair
                ix1n = ix1n + kcount 
            
            enddo
        
        else
            print *, "error: undefined resampling type", samp_type
            write (16, *) "error: undefined resampling type", samp_type
            close(16)
            stop
        endif
   
       ! resample ipp, dtt, rxcor, stlon, stlat, dist
            ! note, we are not resampling station names, for speed
            ! also assumes index1 and index2 are consistent across iterations (type1), or have been fixed above (type2)
        call RESAMPLE_IDATAVEC(ipp00(1:nk), ipp(1:nk), samp_vec(1:nk), nk)
        call RESAMPLE_FDATAVEC(dist00(1:nk), dist(1:nk), samp_vec(1:nk), nk)
        call RESAMPLE_FDATAVEC(tdif00(1:nk), tdif(1:nk), samp_vec(1:nk), nk)
        call RESAMPLE_FDATAVEC(rxcor00(1:nk), rxcor(1:nk), samp_vec(1:nk), nk)
        call RESAMPLE_DFDATAVEC(slon00(1:nk), slon(1:nk), samp_vec(1:nk), nk)
        call RESAMPLE_DFDATAVEC(slat00(1:nk), slat(1:nk), samp_vec(1:nk), nk)
    
   
       !------- compute event pair quality (similarity) factors for resampled data ------!
       iqmax = 0 
       ! loop over pairs   
       do ip = 1, npair 
          ngood(ip) = 0
          rfactor(ip) = -1.                          !***stays at -1 if 'bad' pair      
      
          if (index2(ip)-index1(ip) <= 1) cycle         !****added to make more robust
          if (index1(ip) <= 1) cycle 
          if (index2(ip) <= 1) cycle      
            
          rgoodavg = 0.
      
          ! sum over all observations for a given pair that
          !     (1) have rxcor > rmin and (2) are observed at near enough to station
          do i = index1(ip), index2(ip)
             if (dist(i) <= delmax .and. rxcor(i) >= rmin) then
                ngood(ip) = ngood(ip) + 1         
                rgoodavg = rgoodavg + rxcor(i)
             endif
          enddo
          if (ngood(ip) > 0) rgoodavg = rgoodavg/real(ngood(ip)) ! sum --> average
      
      
          !rfactor: measure of 'quality' of event pair:
          !     rfactor = (# of "good" rxcor values) x (avg of "good rxcor values) 
          rfactor(ip) = real(ngood(ip))*rgoodavg                   
      
          if (iqq1(ip) > iqmax) iqmax = iqq1(ip)
          if (iqq2(ip) > iqmax) iqmax = iqq2(ip)
       enddo
       print *, 'Maximum quake number = ', iqmax
    !-------------------
   
    ! now sort by rfactor (indxr has the highest quality / most similar pair) last
        ! i.e. ip = indxr(npair) gives the pair number of the best pair
       call INDEXX(npair, rfactor, indxr)   

    ! print best pair   
       print *, 'best pair follows (resampled data)'
       ip = indxr(npair) ! index number of best pair
       k = index1(ip)
       print *, k, ip, iqq1(ip), iqq2(ip), rfactor(ip)
   
   endif ! closing if statement on resampling for bootstrap
   
  !--------------------------------------------------------------------------!          

   
   ! loop over pairs, starting with highest rfactor values
   do kraw = npair, 1, -1                 
      
      ip = indxr(kraw)               !pair number
      
      if (rfactor(ip) < 0.) cycle                 !***added to make more robust
      
      i = iqq1(ip)                   !first event number
      j = iqq2(ip)                   !second event number

        ! print out every 1000 or so
      if (mod(kraw, 1000) == 0) then 
         print *, '(sorted) pair # / total # of pairs =          ', real(kraw)/real(npair)      
         print *, '(sorted) pair #, event i, event j, rfactor_ij= ',  &
            kraw, i, j, rfactor(ip)
         print *, 'ntree, nbranch_i, nbranch_j = ', &
                   ntree, nbranch(index(i)), nbranch(index(j))
      endif
      if (i == 0 .or. j == 0) cycle
      if (i == j) cycle
      
      ! pair centroid
      qlat0 = (qlat(i) + qlat(j))/2.
      qlon0 = (qlon(i) + qlon(j))/2.
      qdep0 = (qdep(i) + qdep(j))/2. 
           
! already in same tree, no need to relocate, right?                                             
      if (index(i) == index(j)) then  
!         print *, 'already in same tree'
         cycle
                  
!possibly combine two different existing trees, note that either or both trees can be single events         
      else     

    ! first check to see how many good connections join clusters (if each tree has >1 event)
         if (nbranch(index(i)) > 1 .and. nbranch(index(j)) > 1) then      
            
            ! find total number of pairs linking cluster
            nconnect=0         
            do k = 1, npair
               if (rfactor(k) < 0.0) cycle             !**added to make robust            
               if (ngood(k) < ngoodmin) cycle
               ii = iqq1(k)     
               jj = iqq2(k)
               if ((index(ii) == index(i) .and. index(jj) == index(j)) .or. &
                  (index(jj) == index(i) .and. index(ii) == index(j)) ) then
                  nconnect = nconnect + 1
               endif
            enddo
         
            ntot = nbranch(index(i))*nbranch(index(j))    !total possible number of connections (nbranch_i*nbranch_j)
            frac = real(nconnect)/real(ntot) ! fraction of total
            
        ! if fraction too small --> not enough good connections (don't join clusters/relocate)
            if (frac < conparam) then
!               print *, 'not joining two clusters ', nconnect, ntot, frac
               cycle                 
            endif
         end if
         
! then check to see if clusters are too far apart (recall the t-arrays are tree centroids)
         cosqlat = cos(qlat0/degrad)
         dy = (tlat(index(i)) - tlat(index(j)))*degkm
         dx = (tlon(index(i)) - tlon(index(j)))*degkm*cosqlat
         dz = tdep(index(i)) - tdep(index(j))
         distance = sqrt(dx**2 + dy**2 + dz**2)
         if (distance > distmax) cycle ! don't join or relocate     
         
! ok, now find ten best pairs that link clusters, and (possibly) relocate using those
    ! (the "8" variables/arrays are working variables/arrays for this step
         
         npr8 = 0           ! number of good pairs linking clusters
         npk8 = 0           ! total number of differential times for these pairs
         
         ! loop over all observations (starting w/ best xcor results)
         do kraw8 = kraw, 1, -1
            ip8 = indxr(kraw8)
            if (rfactor(ip8) < 0.0) cycle       !***added to make robust            
            i8 = iqq1(ip8)
            j8 = iqq2(ip8)
            
            ! found a pair linking clusters
            if (index(i8) == index(i) .and. index(j8) == index(j)) then
               npr8 = npr8 + 1
               ip810(npr8) = ip8
               
               ! loop over all observations of this pair
               do k8 = index1(ip8), index2(ip8) 
                  npk8 = npk8 + 1
                  
                  ! copy tdif, ipp, slat, slon to "8" arrays
                  tdif8(npk8) = tdif(k8)
                  ipp8(npk8) = ipp(k8)
                  slat8(npk8) = slat(k8)
                  slon8(npk8) = slon(k8)
                  
                  ! copy over locations, relative to centroid 
                  qlat18(npk8) = qlat(i8) - tlat(index(i8)) 
                  qlon18(npk8) = qlon(i8) - tlon(index(i8))
                  qdep18(npk8) = qdep(i8) - tdep(index(i8))
                  qtim18(npk8) = qtim(i8) - torg(index(i8))
                  qlat28(npk8) = qlat(j8) - tlat(index(j8))
                  qlon28(npk8) = qlon(j8) - tlon(index(j8))
                  qdep28(npk8) = qdep(j8) - tdep(index(j8))
                  qtim28(npk8) = qtim(j8) - torg(index(j8))
                  
                  if (abs(qtim28(npk8)) > 50.) then
                     print *, '***TIME ERROR3: ', qtim28(npk8), qtim(j8), torg(index(j8)), j8, index(j8)
                     print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
                     print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
                     close(16)
                     stop
                  endif                  
               enddo
               
            ! found a pair linking clusters            
            else if (index(i8) == index(j) .and. index(j8) == index(i)) then
               npr8 = npr8 + 1
               ip810(npr8) = ip8
               
               ! loop over all observations of this pair               
               do k8 = index1(ip8), index2(ip8)
                  npk8 = npk8 + 1
                  
                  ! copy tdif, ipp, slat, slon to "8" arrays
                  tdif8(npk8) = -tdif(k8) !note minus tdif here
                  ipp8(npk8) = ipp(k8)
                  slat8(npk8) = slat(k8)
                  slon8(npk8) = slon(k8) 
                  
                  ! copy over locations, relative to centroid                 
                  qlat18(npk8) = qlat(j8) - tlat(index(j8))
                  qlon18(npk8) = qlon(j8) - tlon(index(j8))
                  qdep18(npk8) = qdep(j8) - tdep(index(j8))
                  qtim18(npk8) = qtim(j8) - torg(index(j8))
                  qlat28(npk8) = qlat(i8) - tlat(index(i8))
                  qlon28(npk8) = qlon(i8) - tlon(index(i8))
                  qdep28(npk8) = qdep(i8) - tdep(index(i8))
                  qtim28(npk8) = qtim(i8) - torg(index(i8)) 
                  if (abs(qtim28(npk8)) > 50.) then
                     print *, '***TIME ERROR4: ', qtim28(npk8), qtim(i8), torg(index(i8)), i8, index(i8)
                     print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
                     print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
                     close(16)
                     stop
                  endif                                    
               enddo                                           
            endif
            if (npr8 >= 10) exit                     !stop when we have enough events
         enddo

         if (npr8 == 0) then
            print *, '***ERROR: npr8 = 0'
            close(16)
            stop
         endif
         if (npk8 > n08) then
            print *, '***Error, npk8 too big: ', npk8, n08
            close(16)
            stop
         endif

! prior to (possibly) combining two clusters, relocate cluster pair (relative relocation wrt centroid)        
         qlat0 = (tlat(index(i)) + tlat(index(j)))/2.
         qlon0 = (tlon(index(i)) + tlon(index(j)))/2.
         qdep0 = (tdep(index(i)) + tdep(index(j)))/2.

        ! DIFCLUST performs relative relocation of the two clusters (relative to the centroid of the cluster pair)
        ! using the the npr8 "best" event pairs linking the clusters (the npk8 differential travel times,
        ! etc., include all observations across these linking pairs). This is analogous to DIFLOC, which does the
        ! same thing for two events, relative to their centroid...
        ! Note the outputs here are: qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, cdist, torgdif, &
        !     resid8, rms, rmed, resol, which store the updated cluster locations, origin time difference,
        !     and residual information. (Somewhat confusingly, qlat1/2,qlon1/2,qdep1/2 are new cluster centroids,
        !     not event locations...)
         call DIFCLUST(qlat0, qlon0, qdep0, npk8, tdif8, ipp8, slat8, slon8, qlat18, qlon18, &
             qdep18, qtim18, qlat28, qlon28, qdep28, qtim28, TTtabfile, boxwid, nit, irelonorm, &
             qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, cdist, torgdif, &
             resid8, rms, rmed, resol)
!         print *, 'DIFCLUST: ',qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, cdist, torgdif, rms, rmed, resol


		! added for robustness: reject cluster merger for negative cluster centroid depths
		  if (qdep1 < tt_dep0 .or. qdep2 < tt_dep0) then 
			 print *, '***rejecting cluster merger (negative depth) ', cdist, qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, torgdif
			 cycle
		  endif

         ! added for robustness
          if (abs(torgdif) > 50.) then 
             print *, '***TIME ERROR5: ', torgdif, qtim18(1), qtim28(1), qtim18(npk8), qtim28(npk8)
             print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
             print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
             print *, 'DIAGNOSTICS:'
             print '(f10.5,f11.5,f8.3,f10.5,f11.5,f8.3,f10.5,f11.5,f8.3,f8.3, f8.4,f8.4,f10.7)', &
              qlat0, qlon0, qdep0, qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, torgdif, rms, rmed, resol
             do i = 1, npk8
                print '(i4,f9.4,f10.5,f11.5,f10.5,f11.5,f8.3,f8.3,f10.5,f11.5,f8.3,f8.3)', &
                	i, tdif8(i), slat8(i), slon8(i), qlat18(i), qlon18(i), qdep18(i), qtim18(i), &
                    qlat28(i), qlon28(i), qdep28(i), qtim28(i) 
             enddo
             close(16)
             stop
          endif
          
          ! new criterion to reject cluster merger if rms or median absolute residual is too large
          if (rms > rmsmax .or. rmed > rmedmax) then
            !print *, '***rejecting cluster merger from residual misfit: ', rms, rmed, cdist, npk8
            cycle
         endif
          
          ! reject if relocated distance is too far...
         if (cdist > distmax2) then
            !print *, '***rejecting cluster merger: ', cdist, qlat1, qlon1, qdep1, qlat2, qlon2, qdep2, torgdif
            cycle
         endif
         
! ---- ok, now test how much each cluster centroid shifts -------
         
         ! number of events/cluster
         nbranch_i = nbranch(index(i))
         nbranch_j = nbranch(index(j))
         
         !original centroid (stored in tlat/tlon/tdep arrays)
         clat1 = (tlat(index(i))*nbranch_i + tlat(index(j))*nbranch_j)/real(nbranch_i + nbranch_j)  
         clon1 = (tlon(index(i))*nbranch_i + tlon(index(j))*nbranch_j)/real(nbranch_i + nbranch_j)
         cdep1 = (tdep(index(i))*nbranch_i + tdep(index(j))*nbranch_j)/real(nbranch_i + nbranch_j)
         
         !new centroid (outputs qlat1/2,qlon1/2,qdep1/2 from DIFCLUST)
         clat2 = (qlat1*nbranch_i + qlat2*nbranch_j)/real(nbranch_i + nbranch_j)                    
         clon2 = (qlon1*nbranch_i + qlon2*nbranch_j)/real(nbranch_i + nbranch_j)
         cdep2 = (qdep1*nbranch_i + qdep2*nbranch_j)/real(nbranch_i + nbranch_j)
         
    ! if number events in cluster i is > nclustshiftmin, test if centroid moves too far
         if (nbranch_i >= nclustshiftmin) then
            dy = qlat1 - tlat(index(i)) - (clat2 - clat1)        !new loc - old loc
            dx = qlon1 - tlon(index(i)) - (clon2 - clon1)
            dz = qdep1 - tdep(index(i)) - (cdep2 - cdep1)
         
            if (abs(dz) > vshiftmax) then
               !print *, '***rejecting cluster merger from vshiftmax: ', dz, nbranch_i
               cycle
            end if
            
            cosqlat=cos(qlat0/degrad)
            dx = dx*cosqlat                                                      
            delkm = sqrt(dx**2 + dy**2)*degkm
            if (delkm > hshiftmax) then
               !print *, '***rejecting cluster merger from hshiftmax: ', delkm, nbranch_i            
               cycle
            endif
            
         endif
         
    ! if number events in cluster j is > nclustshiftmin, test if centroid moves too far        
         if (nbranch_j >= nclustshiftmin) then
            dy = qlat2 - tlat(index(j)) - (clat2 - clat1)        !new loc - old loc
            dx = qlon2 - tlon(index(j)) - (clon2 - clon1)
            dz = qdep2 - tdep(index(j)) - (cdep2 - cdep1)
         
            if (abs(dz) > vshiftmax) then
               !print *, '***rejecting cluster merger from vshiftmax: ', dz, nbranch_j
               cycle
            end if

            cosqlat=cos(qlat0/degrad)                                                      
            dx = dx*cosqlat
            delkm = sqrt(dx**2 + dy**2)*degkm
            if (delkm > hshiftmax) then
               !print *, '***rejecting cluster merger from hshiftmax: ', delkm, nbranch_j            
               cycle
            endif            

         endif         
      
 !-----------------------------------------------------------------!                   

! ok, the cluster merge is now "approved". now adjust locations of all events  
! in clusters index(i) and index(j), and join to same cluster.

        ! cluster 1: new loc - old loc
         qlat_off1 = qlat1 - tlat(index(i))         
         qlon_off1 = qlon1 - tlon(index(i))
         qdep_off1 = qdep1 - tdep(index(i))
         qtim_off1 = torg(index(i)) - torgdif/2.     !***check this
         if (abs(qtim_off1) > 50.) then
            print *, '***TIME ERROR6: ', qtim_off1, torg(index(i)), index(i), i
            print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
            print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
            close(16)
            stop
         endif
         ! cluster 2: new loc - old loc
         qlat_off2 = qlat2 - tlat(index(j))         !
         qlon_off2 = qlon2 - tlon(index(j))
         qdep_off2 = qdep2 - tdep(index(j))
         qtim_off2 = torg(index(j)) + torgdif/2.     !***check this
         if (abs(qtim_off2) > 50.) then
             print *, '***TIME ERROR7: ', qtim_off2, torg(index(j)), index(j), j
             print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
             print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
            close(16)
            stop
         endif          
                          
         ntree = ntree - 1          ! tree merger
         itree = index(i)           !use first event tree number for both
         itreeold = index(j)        ! the second tree is "evacuated"
         
         ! new total number of branches
         nbranch_i = nbranch(index(i))
         nbranch_j = nbranch(index(j))
         nb = nbranch_i + nbranch_j 
         n = 0
         
         ! reassign event locations for both trees                                                                        
         do ii = 1, nq
            kk = index(ii)               
            if (kk == itree) then ! events in new tree = tree (cluster) i
               qlat(ii) = qlat(ii) + qlat_off1
               qlon(ii) = qlon(ii) + qlon_off1
               qdep(ii) = qdep(ii) + qdep_off1
               qtim(ii) = qtim(ii) + qtim_off1
               if (abs(qtim(ii)) > 50.) then
                  print *, '***TIME ERROR8: ', qtim(ii), qtim_off1, qtim_off2, torgdif
                  print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
                  print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
                  print *, i, j, index(i), index(j), torg(index(i)), torg(index(j))
                  close(16)
                  stop
               endif
               tavg = tavg + qtim(ii) ! *****DT note - not needed, we take care of this below
               index(ii) = itree
               n = n + 1
               inbranch(n) = ii                  !save quake number for each branch
            else if (kk == itreeold) then ! events in oldtree = tree (cluster) j
               qlat(ii) = qlat(ii) + qlat_off2
               qlon(ii) = qlon(ii) + qlon_off2
               qdep(ii) = qdep(ii) + qdep_off2
               qtim(ii) = qtim(ii) + qtim_off2
               if (abs(qtim(ii)) > 50.) then
                  print *, '***TIME ERROR9: ', qtim(ii), qtim_off2, qtim_off1, torgdif
                  print *, 'LARGE ORIGIN TIME CORRECTION, LIKELY XCOR DATA OR EVLIST PROBLEM.'
                  print *, 'CHECK CURRENT PAIR (ipair,qnum1,qnum2,evid1,evid2):', ip, i, j, idcusp11(ip), idcusp22(ip)
                  print *, i, j, index(i), index(j), torg(index(i)), torg(index(j))
                  close(16)                  
                  stop
               endif                    
               index(ii) = itree
               n = n + 1
               inbranch(n) = ii !save quake number for each branch
            endif
            
        enddo
         if (n /= nb) then
            print *, '***nbranch problem: ', n, nb
            close(16)
            stop
         endif
         nbranch(itree) = nb                    !first tree now has all events
         nbranch(itreeold) = 0                  !zero out other tree
         
         
! write to growlog file (only for unresampled data)
    ! itree, nbranch(itree): new tree number (same as event itree), 
    ! itree: new (merged) tree/cluster number, which has been set to the previous tree number for event 1
    ! itreeold: old tree/cluster number, which was the previous tree number for event 2
    ! nbranch_i, nbranch_j: number of events in the previous tree/cluster for events 1 and 2, respectively
    ! itree, nbranch(itree): new (merged) tree/cluster number, and new number of events in this merged tree
    ! npr8: number of event pairs linking trees/clusters (up to 10 best)
    ! npk8: total number of differential times for these even pairs
    ! qlat2, qlat1: tree/cluster centroid latitudes for original cluster 1 and 2
    ! qlon2, qlon1: tree/cluster centroid longitudes for original cluster 1 and 2
    ! qdep2, qdep1: tree/cluster centroid depths for original cluster 1 and 2
    ! distance: distance between the original two tree/cluster centroids (km)
    ! rms, rmed of travel time residuals (added to log output July 2016)
         if (nitb == 0) then
         write (16, 166) itree, itreeold, nbranch_i, nbranch_j, itree, nbranch(itree), npr8, npk8, &
            qlat2-qlat1, qlon2-qlon1, qdep2-qdep1, distance, rms, rmed
166      format (2i7, 2i6, i7, i6, i3, i4, 4f10.5, 2f7.2)
         endif       
!         print *, 'two clusters merged ', itree, nbranch(itree), npick, ntree
      
      end if      
!      print *, 'itree, nbranch(itree) = ', itree, nbranch(itree), i, j, itreeold     
      if (nbranch(itree) >= nbmax) then
         print *, '***Error: too many branches'
         close(16)
         stop
      endif
!------------------------------------------------!
            
! ------ now shift cluster to catalog centroid, set average time to zero ---------------!
      nb = nbranch(itree)
      n = nb
      if (nb < 2) then
         print *, '***ERROR, nb, itree = ', nb, itree
         close(16)
         stop
      endif
      qlat0 = 0.
      qlon0 = 0.
      qdep0 = 0.
      tlat(itree) = 0.
      tlon(itree) = 0.
      tdep(itree) = 0.
      tavg = 0.            
      do j = 1, nb 
         iq = inbranch(j)
         if (index(iq) /= itree) then
            print *,'**ERROR: ', j, iq, inbranch(iq), itree
            close(16)
            stop
         endif
         
!          if (abs(qlat(iq)-qlat_cat(iq)) > 0.1) then
!             print *,'***WARNING: catalog latitude more than 0.1 degree away'
!             print *, j,iq,qlat(iq),qlat_cat(iq),idcusp(iq)
!          endif
!          if (abs(qlon(iq)-qlon_cat(iq)) > 0.1) then
!             print *,'***WARNING: catalog longitude more than 0.1 degree away'
!             print *, itree,nb,j,iq,qlon(iq),qlon_cat(iq),idcusp(iq)
!          endif         
         
         qlat0 = qlat0 + qlat(iq)
         qlon0 = qlon0 + qlon(iq)
         qdep0 = qdep0 + qdep(iq)
                  
         tlat(itree) = tlat(itree) + qlat_cat(iq)
         tlon(itree) = tlon(itree) + qlon_cat(iq)
         tdep(itree) = tdep(itree) + qdep_cat(iq)
         tavg = tavg + qtim(iq)
      enddo

      qlat0 = qlat0/real(n)
      qlon0 = qlon0/real(n)
      qdep0 = qdep0/real(n)
      tavg = tavg/real(n)
      tlat(itree) = tlat(itree)/real(n)
      tlon(itree) = tlon(itree)/real(n)
      tdep(itree) = tdep(itree)/real(n)
      torg(itree) = 0.                        !set cluster average time to zero
      qlat_off = tlat(itree) - qlat0
      qlon_off = tlon(itree) - qlon0
      qdep_off = tdep(itree) - qdep0

      do j = 1, nb
         iq = inbranch(j)
         qlat(iq) = qlat(iq) + qlat_off
         qlon(iq) = qlon(iq) + qlon_off
         qdep(iq) = qdep(iq) + qdep_off
         qtim(iq) = qtim(iq) - tavg
      enddo
            
   enddo    !----------end loop on kraw: done with clustering and relocations for iteration 
   
   
! now compute sequential set of tree numbers (stored in index2), check that nbranch is correct
    ! these are sorted, so new tree #1 has the most events...
    !   note that index maps from quake number to tree number, while
    !   index2 (now!) maps from the (old) tree number used in the growclust computations to a
    !   to a new sequential set of tree numbers starting at 1. 
    !       i.e. index(iq) = itree1, non-sequential tree number
    !            index2(itree1) = itree2, sequential tree number
 

  ! first, get sort indices for old tree numbers (index), based on nbranch
    ! note index(iq) gives old tree number, and nb_indxr(1) gives the old tree number with the highest nbranch
        ! the -1 makes indexx sort in descending order...
  call indexx(nq, -1.0*real(nbranch(1:nq)), nb_indxr(1:nq))
   
   ! in the new version, we need to save old tree locations so we don't overwrite
   do j = 1,nq
    if (nbranch(j) > 0) then
      tlon00(j) = tlon(j)
      tlat00(j) = tlat(j)
      tdep00(j) = tdep(j)
      torg00(j) = torg(j)
    endif
   enddo
   
   ! now assign new tree numbers to index2
    ! note that for event number iq, 
    ! old tree number: j = index(iq)
    ! new tree number: i = index2(j) = index2(index(iq))       
   n = 0
   do i = 1, nq ! loop is over NEW tree # --> i
    
    j = nb_indxr(i) ! get OLD tree # --> j
    
    if (nbranch(j) > 0) then ! check for OLD empty trees
        
        index2(j) = i ! set new tree number (j = index(iq) = OLD)
        nbranch2(i) = nbranch(j) ! save old nbranch, as a later check
        
        ! set new tree coordinates
        tlat(i) = tlat00(j)
        tlon(i) = tlon00(j)
        tdep(i) = tdep00(j)
        torg(i) = torg00(j)
        n = n + 1
        
    endif     
   enddo
   if (n /= ntree) then
      print *, '***n /= ntree ', n, ntree
      close(16)
      stop
   endif
   
  ! compute final nbranch array (nevents/tree)
   nbranch = 0   
   do iq = 1, nq
      i = index2(index(iq)) ! index2: new tree number itree2 of old tree itree1, index(iq): old tree of event iq
      nbranch(i) = nbranch(i) + 1
   enddo
   do j = 1, ntree
      if (nbranch(j) /= nbranch2(j)) then
         print *, '**nbranch error ', j, ntree, nbranch(j), nbranch2(j)
         close(16)
         stop
      endif
   enddo
   
   !-----------------------------------------------------------------------------------!  
   
   ! bootstrap iteration: save bootstrap locations
   if (nitb > 0) then
   
   ! loop over events
    do iq = 1, nq

        ! lat, lon, dep
        Blon(nitb,iq) = qlon(iq)
        Blat(nitb,iq) = qlat(iq)
        Bdep(nitb,iq) = qdep(iq)
        Btim(nitb,iq) = qtim(iq)
    
        ! number of neighbors of this event that live in the tree
        itree = index2(index(iq))
        Bnb(nitb,iq) = nbranch(itree)
        Bitreeq(nitb,iq) = itree 
    
    enddo
        
  else ! first iteration (non-resampled)
  
 ! save final positions for output file, prior to bootstrapping
   do iq = 1, nq
   
        ! lon, lat, dep, timing
        qlonR(iq) = qlon(iq)
        qlatR(iq) = qlat(iq)
        qdepR(iq) = qdep(iq)
        qtimR(iq) = qtim(iq)
        
        ! cluster number and nbranch
        itree = index2(index(iq))
        itreeqR(iq) = itree
        nbranchqR(iq) = nbranch(itree)
        
        ! coordinates with respect to tree centroid (km)
        cosqlat = cos(tlat(itree)/degrad)
        qcxR(iq) = (qlonR(iq)-tlon(itree))*degkm*cosqlat
        qcyR(iq) = (qlatR(iq)-tlat(itree))*degkm
        qczR(iq) = qdepR(iq)-tdep(itree)
        
   enddo
   
   ! save tlat, tlon, tdep, torg, nbranchtR, ntreeR before bootstrapping
   ntreeR = ntree
   do ii = 1, ntree
      tlonR(ii) = tlon(ii)
      tlatR(ii) = tlat(ii)
      tdepR(ii) = tdep(ii)
      torgR(ii) = torg(ii)
      nbranchtR(ii) = nbranch(ii)
   enddo 
  
   
   ! ------- compute RMS residuals between observed and predicted travel times ------------------
    print *, 'GROWCLUST relocation complete for unresampled data'
    call TIMER 
    print *, ' '
    print *, 'Computing travel time residuals...'
   
   ! first, zero out q arrays
   do iq = 1, nq
     qnpair(iq) = 0
     qndiffP(iq) = 0
     qndiffS(iq) = 0
     qrms(iq) = 0.0
     qrmsP(iq) = 0.0
     qrmsS(iq) = 0.0
    enddo
   
   ! variables to store running totals over all events 
   CrmsP = 0.0
   CrmsS = 0.0 
   npairC = 0 
   npC = 0
   nsC = 0
   jjP = 0
   jjS = 0
   ! loop over event pairs, and compute predicted travel times for all picks
        !(WILL NEED TO MODIFY IF WE DO THIS FOR BOOTSTRAP ITERATIONS)       
   do ip = 1, npair
      
      ! get index1 and index2 for this pair from unresampled data 
      ix1 = index100(ip)
      ix2 = index200(ip)
      
      if (ix2 > 0 .and. ix1 > 0) then ! check for "empty" pair 
      
      ! get sequential quake numbers
      iq1 = iqq1(ip)
      iq2 = iqq2(ip)
      
      ! get cluster/tree numbers
      itree1 = itreeqR(iq1)
      itree2 = itreeqR(iq2)
      
      ! compute predicted travel times 
        cosqlat = cos( 0.5*(qlat(iq1)+qlat(iq2)) / degrad )
        npickgoodP = 0 ! counts number of good picks
        npickgoodS = 0
        prmsP = 0.0 ! set RMS for pair to zero, initially
        prmsS = 0.0
        do k = ix1, ix2
            
            ! predicted travel time from event1 to station
            dy = qlat(iq1) - slat(k)
            dx = (qlon(iq1) - slon(k))*cosqlat
            distance = sqrt(dx**2 + dy**2)*degkm
            !print *, k, ipp(k), TTtabfile(ipp(k)), distance, qdep(iq1)
            call GET_TTS_FAST8(TTtabfile(ipp(k)), ipp(k), distance, qdep(iq1), tsec1, iflag)
            
            ! predicted travel time from event1 to station
            dy = qlat(iq2) - slat(k)
            dx = (qlon(iq2) - slon(k))*cosqlat
            distance = sqrt(dx**2 + dy**2)*degkm
            !print *, k, TTtabfile(ipp(k)), ipp(k), distance, qdep(iq2)
            call GET_TTS_FAST8(TTtabfile(ipp(k)), ipp(k), distance, qdep(iq2), tsec2, iflag)

            ! predicted differential time      
            tdif_pred(k) = tsec2 - tsec1
            
            ! check to see if this is a "good" pick
            if (dist(k) <= delmax .and. rxcor(k) >= rmin) then
            
              if (ipp(k) == 1) then ! P-phase
                
                ! store residual in P-residual array
                jjP = jjP + 1
                resP(jjP) = tdif(k)-tdif_pred(k)
                
                npickgoodP = npickgoodP + 1  ! increment npickgood
                prmsP = prmsP + (resP(jjP))**2 ! increment residual sum for the pair
                
                
              else                  ! S-phase
                
                ! store residual in P-residual array
                jjS = jjS + 1
                resS(jjS) = tdif(k)-tdif_pred(k)
                
                npickgoodS = npickgoodS + 1  ! increment npickgood
                prmsS = prmsS + (resS(jjS))**2 ! increment residual sum for the pair
                
              endif
                
            endif
            
          enddo
        
        ! update rms, npair, ndiff
        if ((itree1 == itree2) .and. (npickgoodP + npickgoodS > 0) ) then !--------
        
            ! first, increment npair
            qnpair(iq1) = qnpair(iq1) + 1
            qnpair(iq2) = qnpair(iq2) + 1
            npairC = npairC + 1
        
            ! P-phase
            if (npickgoodP > 0) then
            
                ! increment residual sum for each event
                qrmsP(iq1) = qrmsP(iq1) + prmsP
                qrmsP(iq2) = qrmsP(iq2) + prmsP
                
                ! increment qndiff for each event
                qndiffP(iq1) = qndiffP(iq1) + npickgoodP
                qndiffP(iq2) = qndiffP(iq2) + npickgoodP
                
                ! increment running totals
                CrmsP = CrmsP + prmsP
                npC = npC + npickgoodP
                
            endif
            
            ! S-phase
            if (npickgoodS > 0) then
                
                ! increment residual sum for each event
                qrmsS(iq1) = qrmsS(iq1) + prmsS
                qrmsS(iq2) = qrmsS(iq2) + prmsS
                
                ! increment qndiff for each event
                qndiffS(iq1) = qndiffS(iq1) + npickgoodS
                qndiffS(iq2) = qndiffS(iq2) + npickgoodS
                
                ! increment running totals
                CrmsS = CrmsS + prmsS
                nsC = nsC + npickgoodS
                
            endif   
        endif !----------------------------------------------------------------
      
      endif
   
   enddo
   
    !SSE running totals: convert SSE to RMS
   CrmsP = sqrt(CrmsP/float(npC))
   CrmsS = sqrt(CrmsS/float(nsC))
   
   ! compute mean and median residuals
   call MEAN(resP(1:jjP), jjP, resPmean)
   !call MEDIAN(resP(1:jjP), jjP, resPmed) ! removed for speed (sorting really large arrays is a chore)
   call MEAN(resS(1:jjS), jjS, resSmean)
   !call MEDIAN(resS(1:jjS), jjS, resSmed) ! removed for speed (sorting really large arrays is a chore)
   
   
   ! now loop over events and convert residual sum to rms
   do iq = 1, nq
   
     if (qndiffP(iq)+qndiffS(iq) > 0) then ! at least one good data
        
        ! total SSE --> RMS
        qrms(iq) = sqrt( (qrmsP(iq) + qrmsS(iq)) / float(qndiffP(iq) + qndiffS(iq)) )
   
        ! P SSE --> RMS
        if (qndiffP(iq) > 0) then 
            qrmsP(iq) = sqrt( qrmsP(iq) / float(qndiffP(iq)) ) ! SSE --> RMS
        else
            qrmsP(iq) = rms_nan ! set to NaN flag rather than divide by 0
        endif
        
        ! S SSE --> RMS
        if (qndiffS(iq) > 0) then 
            qrmsS(iq) = sqrt( qrmsS(iq) / float(qndiffS(iq))) ! SSE --> RMS
        else
            qrmsS(iq) = rms_nan ! set to NaN flag rather than divide by 0
        endif
     
     else ! set to NaN flag rather than divide by 0
        qrms(iq) = rms_nan
        qrmsP(iq) = rms_nan
        qrmsS(iq) = rms_nan
     endif
     
   enddo
   
   !---------- final outputs for log file  ------------------

   ! first loop over events: count if relocated
   nreloc = 0
   do iq = 1, nq
    if (nbranchqR(iq) > 1) nreloc = nreloc + 1
   enddo
   
   ! now count number of clusters with 2, 5, 10, 20, 50, 100 events
   ntree2 = 0
   ntree5 = 0
   ntree10 = 0
   ntree20 = 0
   ntree50 = 0
   ntree100 = 0
   do kk = 1, ntreeR
    if (nbranch(kk) >= 2) ntree2 = ntree2 + 1
    if (nbranch(kk) >= 5) ntree5 = ntree5 + 1
    if (nbranch(kk) >= 10) ntree10 = ntree10 + 1
    if (nbranch(kk) >= 20) ntree20 = ntree20 + 1
    if (nbranch(kk) >= 50) ntree50 = ntree50 + 1
    if (nbranch(kk) >= 100) ntree100 = ntree100 + 1
   enddo
    
   
   write(16, *) '**************************************************************************************************'
   write(16, *) ' '
   write(16, *) '************************************  GROWCLUST Run Summary  *************************************'
    write(16, *) ' '
   write(16, '(a55, i10)') 'Number of catalog events: ', nq
   write(16, '(a55, i10)') 'Number of relocated events: ', nreloc
   write(16, '(a55, i10)') 'Number of input event pairs: ', npair
   write(16, '(a55, i10)') 'Number of event pairs used: ', npairC
   write(16, '(a55, i10)') 'Number of xcor data used (total, P+S): ', npC + nsC
   write(16, '(a55, i10)') 'Number of xcor data used (P-phase): ', npC
   write(16, '(a55, f10.4)') 'RMS differential time residual (P-phase): ', CrmsP
   write(16, '(a55, f10.4)') 'Mean (signed) differential time residual (P-phase): ', resPmean
   !write(16,'(a55, f10.4)') 'Median (signed) differential time residual (P-phase): ', &
   !      resPmed
   write(16, '(a55, i10)') 'Number of xcor data used (S-phase): ', nsC
   write(16, '(a55, f10.4)') 'RMS differential time residual (S-phase): ', CrmsS
   write(16, '(a55, f10.4)') 'Mean (signed) differential time residual (S-phase): ', resSmean
   !write(16, '(a55, f10.4)') 'Median differential time residual (S-phase): ', &
   !     resSmed
   write(16, *) ' '
   write(16, '(a55, i9)') 'Number of clusters with >=   2 events: ', ntree2
   write(16, '(a55, i9)') 'Number of clusters with >=   5 events: ', ntree5
   write(16, '(a55, i9)') 'Number of clusters with >=  10 events: ', ntree10
   write(16, '(a55, i9)') 'Number of clusters with >=  20 events: ', ntree20
   write(16, '(a55, i9)') 'Number of clusters with >=  50 events: ', ntree50
   write(16, '(a55, i9)') 'Number of clusters with >= 100 events: ', ntree100
   
   ! close log-file
   close(16)
   
   endif ! end if statement for zeroth iteration (pre-bootstrap)

    print *, 'Done.'
    if (nitb > 0) print *, "End of bootstrap iteration: ", nitb
    call TIMER
    print *, '--------------------------------------------------'
    print *, ' '
    print *, ' '
    
 enddo ! end loop over bootstrap iterations -------------------------------
 
 !------------------------------------------------------------!
 
 ! ----- Compute Bootstrap statistics -----------
  if (nboot > 2) then 
    
    ! compute standard errors and median absolute deviations
        
    do iq = 1,nq
        
        ! --- first, compute bootstrap standard errors and bootstrap mean ----------
            
        ! compute mean and std. dev in lon and lat
        call DPMEAN_STDDEV(Blon(1:nboot,iq), nboot, BlonMEAN(iq),SEx)
        call DPMEAN_STDDEV(Blat(1:nboot,iq), nboot, BlatMEAN(iq),SEy)
        
        ! convert from lat-lon SE to xy (km)
        cosqlat = cos(qlatR(iq)/degrad)
        SEx = SEx*cosqlat
        SEh_boot(iq) = sqrt(SEx**2 + SEy**2)*degkm
        
        ! compute mean and std. dev in dep
        call DPMEAN_STDDEV(Bdep(1:nboot,iq), nboot, BdepMEAN(iq), SEz_boot(iq))
        
         ! compute mean and std. dev in qtime
        call DPMEAN_STDDEV(Btim(1:nboot,iq), nboot, BtimMEAN(iq), SEt )
        SEt_boot(iq) = SEt
        
        ! compute mean nbranch
        call MEAN(real(Bnb(1:nboot, iq)), nboot, BnbMEAN(iq))
        
        ! compute max, min, and number of "iso" branches (in a cluster by itself)
        BnbMAX(iq) = MAXVAL(Bnb(1:nboot, iq))
        BnbMIN(iq) = MINVAL(Bnb(1:nboot, iq))
        BnbISO(iq) = 0
        do ib = 1, nboot
            if (Bnb(ib,iq) < 2) BnbISO(iq) = BnbISO(iq) + 1
        enddo
        
        ! ---- now MAD errors ----------------!
        call GETDPMAD(Blon(1:nboot,iq), nboot, MADx)
        call GETDPMAD(Blat(1:nboot,iq), nboot, MADy)
        MADx = MADx*cosqlat
        MADh_boot(iq) = sqrt(MADx**2 + MADy**2)*degkm
        
        call GETDPMAD(Bdep(1:nboot,iq), nboot, MADz)
        MADz_boot(iq) = MADz
        
        call GETDPMAD(Btim(1:nboot,iq), nboot, MADt)
        MADt_boot(iq) = MADt
        
        ! set SE to nan if not relocated
        if (qndiffP(iq)+qndiffS(iq) < 1 .or. BnbMAX(iq) < 2) then
             SEh_boot(iq) = SE_nan
             SEz_boot(iq) = SE_nan
             SEt_boot(iq) = SE_nan
             MADh_boot(iq) = MAD_nan
             MADz_boot(iq) = MAD_nan
             MADt_boot(iq) = MAD_nan
        endif
        
    enddo
    
 else ! no bootstrapping: set SE to NaN flags
    do iq = 1,nq
        SEh_boot(iq) = SE_nan
        SEz_boot(iq) = SE_nan
        SEt_boot(iq) = SE_nan
        MADh_boot(iq) = MAD_nan
        MADz_boot(iq) = MAD_nan
        MADt_boot(iq) = MAD_nan
    enddo
endif 
 
 
!-------- Make first two output files: cat and clust (whether or not we bootstrapped)-------- 
   
! first output: relocated catalog (out.growclust_cat)
   open (13, file = outfile_cat)
   
   do iq = 1, nq
                
   write (13, 860) qyr_cat(iq), qmon_cat(iq), qdy_cat(iq), qhr_cat(iq), qmn_cat(iq), &
   qsc_cat(iq)+qtimR(iq), idcusp(iq), qlatR(iq), qlonR(iq), qdepR(iq), qmag_cat(iq), &
   iq,itreeqR(iq),nbranchqR(iq),qnpair(iq),qndiffP(iq),qndiffS(iq),qrmsP(iq),qrmsS(iq), &
   MADh_boot(iq), MADz_boot(iq), MADt_boot(iq), qlat_cat(iq), qlon_cat(iq), qdep_cat(iq)           
860   format (i4, 4i3, f7.3, i10, f10.5, f11.5, f8.3, f6.2,  &                          ! 11/2016 changed lat from f9.5 to f10.5 for negative latitudes
              3i8, 3i6, 2f6.2, &
              3f8.3, 2x, f10.5, f11.5, f8.3)   
   enddo
   
   close (13)

! second output (optional): cluster location output (out.growclust_clust)
   if (outfile_clust(1:4) .ne. "none") then    
   open (13, file = outfile_clust)
      
   do ii = 1, ntreeR
      if (nbranchtR(ii) < nbranch_min) cycle
      write (13, 870) ii, nbranchtR(ii), tlatR(ii), tlonR(ii), tdepR(ii), torgR(ii)
870   format (i8, i8, f10.5, f11.5, f8.3, f8.3)
      do iq = 1, nq
      
       if (itreeqR(iq) /= ii) cycle
                                              
       write (13, 880) ii, iq, idcusp(iq), qmag_cat(iq), qyr_cat(iq), qmon_cat(iq), &
          qdy_cat(iq), qhr_cat(iq), qmn_cat(iq), qsc_cat(iq), qlatR(iq), qlonR(iq), &
          qdepR(iq), qcxR(iq), qcyR(iq), qczR(iq), MADh_boot(iq), MADz_boot(iq), &
          MADt_boot(iq), qlat_cat(iq), qlon_cat(iq), qdep_cat(iq)                
880      format (i8, i9, i10, f6.2, i5, 4i3, f7.3, f10.5, f11.5, f8.3, &                ! 11/2016 changed lat from f9.5 to f10.5 for negative latitudes
             3f10.4, 3f8.3, 2x, f10.5, f11.5, f8.3)      
      enddo
   enddo
   close (13)
   endif
 
 
 if (nboot > 2) then ! third output (optional): full bootstrap distribution (out.growclust_bootstrap)
   
   if (outfile_boot(1:4) .ne. "none") then
       
       open (13, file = outfile_boot)
       
       write (13, '(i8, i6)') nq, nboot
       
       do iq = 1, nq
   
         write(13, 864) iq, idcusp(iq), qlatR(iq), qlonR(iq), qdepR(iq), &
         itreeqR(iq), nbranchqR(iq), BnbMEAN(iq), BnbMIN(iq), BnbMAX(iq), &
         BlatMEAN(iq), BlonMEAN(iq), BdepMEAN(iq), &
         MADh_boot(iq), MADz_boot(iq), MADt_boot(iq), &
         SEh_boot(iq), SEz_boot(iq), SEt_boot(iq) 
864   format (i8, i10, f10.5, f11.5, f8.3, 2i8, f7.2, 2i8, f10.5, f11.5, f8.3, 6f8.3)   ! 11/2016 changed lat from f9.5 to f10.5 for negative latitudes
                    
        do ib = 1, nboot            
          write(13,'(f10.5, x, f11.5, x, f8.5, x, i8, x, i8)') &                        ! 11/2016 changed lat from f9.5 to f10.5 for negative latitudes
          Blat(ib,iq), Blon(ib,iq), Bdep(ib,iq), Bitreeq(ib,iq), Bnb(ib,iq)
        enddo                          
   
       enddo
       close(13)
   endif
     
 endif
   
   end program growclust
   
   
   !------------------------------------------------------------------------!
!--------------------------END PROGRAM-----------------------------------!
!------------------------------------------------------------------------!


!--------------------------SUBROUTINES----------------------------------------------------!

! DIFLOC performs differential location for an event pair with respect to its centroid. 
! The method uses an iterative ("shrinking-box") grid search approach to obtain the best (L1)
! relative locations.
!
!  Inputs:  qlat0  =  mean event latitude
!           qlon0  =  mean event longitude
!           qdep0  =  mean event depth (km)
!           npick  =  number of differential times
!           tt     =  array (len=npick) of dif times, t2-t1 (s), 
!           ip     =  array (len=npick) with phase index numbers (1 to 10) for tt data
!           slat   =  array (len=npick) with station latitudes
!           slon   =  array (len=npick) with station longitudes
!           TTtabfile  =  array with phase names (travel time tables)
!           boxwid =  starting box width (km)
!           nit    =  number of iterations to perform
!           inorm  =  relocation norm (1=L1, 2=L2, 3=robust L2)
! Returns:  qlat1  =  best-fitting latitude for first event
!           qlon1  =  best-fitting longitude
!           qdep1  =  best-fitting depth (km)
!           qlat2  =  best-fitting latitude for second event
!           qlon2  =  best-fitting longitude
!           qdep2  =  best-fitting depth (km)
!           torgdif=  origin time difference, i.e., t2-t1 median residual (>0 when 2 is later than 1)
!           resid  =  array (len=npick) of residuals (s) between observed tdif (tt) and predicted
!           rms    =  rms residual ( sqrt( sum(resid**2)/npick) )
!           rmed   =  median absolute value residual ( median( abs(resid) ) )
!           resol  =  nominal resolution (m) of final box
!
subroutine DIFLOC(qlat0,qlon0,qdep0,npick,tt,ip,slat,slon, &
           TTtabfile,boxwid,nit,inorm,qlat1,qlon1,qdep1,qlat2,qlon2,qdep2, &
           torgdif,resid,rms,rmed,resol)
   real*8 qlat0, qlon0, qlat1, qlon1, qlat2, qlon2, dlat, dlon, dlat0, dlon0, flat1, flat2, flon1, flon2, &
          flatbest1, flonbest1, flatbest2, flonbest2
   real tt(npick), resid(npick)
   real*8 slat(npick), slon(npick)
   real rabs(3000)
   integer ip(npick)
   character*100 TTtabfile(2)
   integer inorm
   degrad = 180./3.1415927
   degkm = 111.19493

    ! define initial search box
   dlat0 = 0.
   dlon0 = 0.
   ddep0 = 0.
   dlat = 0.5*boxwid/degkm
   cosqlat = cos(qlat0/degrad)
   dlon = dlat/cosqlat
   ddep = 0.5*boxwid

    ! start iteration
   do it = 1, nit

    ! grid search over box for best fit 
      fitbest = 9.e20
      
      ! get trial locations: f1, f2 (3x3 per box)
      do iy = -1, 1
         flat1 = qlat0 + dlat0 + dlat*float(iy)
         flat2 = qlat0 - dlat0 - dlat*float(iy)
         do ix = -1, 1
            flon1 = qlon0 + dlon0 + dlon*float(ix)
            flon2 = qlon0 - dlon0 - dlon*float(ix)
            do iz = -1, 1
               fdep1 = qdep0 + ddep0 + ddep*float(iz)
               fdep2 = qdep0 - ddep0 - ddep*float(iz)
!               if (fdep1 < 0.) fdep1 = 0.                     !***now permit negative depths
!               if (fdep2 < 0.) fdep2 = 0.
        
            ! compute predicted travel time and residual with observed  
               do i = 1, npick
                  dy = flat1 - slat(i)
                  dx = (flon1 - slon(i))*cosqlat
                  delkm = sqrt(dx**2 + dy**2)*degkm
                  call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, fdep1, tsec1, iflag)
                  dy = flat2 - slat(i)
                  dx = (flon2 - slon(i))*cosqlat
                  delkm = sqrt(dx**2 + dy**2)*degkm
                  call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, fdep2, tsec2, iflag)
                  tdif = tsec2 - tsec1
                  resid(i) = tt(i) - tdif
               enddo
               
                ! compute fit --> L1 norm of residual (sum of absolute residuals) or L2 norm
               fit=0.
               if (inorm .eq. 1) then ! L1 norm
                  call MEDIAN(resid,npick,residmed) ! compute median residual between observed and predicted tdif (torgdif)
                  do i=1,npick
                     fit=fit+abs(resid(i)-residmed)
                  enddo
               else if (inorm .eq. 2) then ! L2 norm
                  call MEAN(resid,npick,residmed) ! compute mean residual between observed and predicted tdif (torgdif)
                  do i=1,npick
                     fit=fit+(resid(i)-residmed)**2
                  end do
               else if (inorm .eq. 3) then ! robomean
                 !xgap=0.1
                 call ROBOMEAN2(resid,npick,0.1,10,residmed,fit)
                  ! note no need to iterate over resid-vector - done internally, stored in fit
               else
                   print *, 'ERROR in DIFLOC: undefined inorm:', inorm
                   stop 
               end if
               
               ! update fitbest, if necessary
               if (fit < fitbest) then
                  fitbest = fit
                  flatbest1 = flat1
                  flonbest1 = flon1
                  fdepbest1 = fdep1
                  flatbest2 = flat2
                  flonbest2 = flon2
                  fdepbest2 = fdep2
                  tbest = residmed  
               end if
            enddo                   !iz loop
         enddo                   !ix loop
      enddo                   !iy loop


      dlat0 = flatbest1 - qlat0
      dlon0 = flonbest1 - qlon0
      ddep0 = fdepbest1 - qdep0
      qorg = tbest                    !origin time difference = median residual between observed and predicted tdif
      torgdif = qorg

      dlat = dlat*0.67      !shrink box by 2/3 each iteration
      dlon = dlon*0.67
      ddep = ddep*0.67
   
   enddo                !end it (iteration) loop

    ! output final location
   qlat1 = flatbest1
   qlon1 = flonbest1
   qdep1 = fdepbest1
   qlat2 =  flatbest2
   qlon2 = flonbest2
   qdep2 = fdepbest2
   resol = (dlat/0.67)*degkm

    ! output final misfits
   rms = 0.
   do i = 1, npick
   
      ! event 1: for best-fit location, compute source station distance and predicted travel time   
      dy = qlat1 - slat(i)
      dx = (qlon1 - slon(i))*cosqlat
      delkm = sqrt(dx**2 + dy**2)*degkm
      call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, qdep1, tsec1, iflag) 
      
      ! event 2: for best-fit location, compute source station distance and predicted travel time   
      dy = qlat2 - slat(i)
      dx = (qlon2 - slon(i))*cosqlat
      delkm = sqrt(dx**2 + dy**2)*degkm
      call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, qdep2, tsec2, iflag)
      
      ! predicted differential time (t2-t1)
      tdif = tsec2 - tsec1
      
      ! compute residual (observed-predicted tdif), then increment sse, update rabs vector
      resid(i) = tt(i) - tdif - qorg
      rms = rms + resid(i)**2
      rabs(i) = abs(resid(i))
   enddo
   
   ! return rms residual and median absolute residual residual
   rms = sqrt(rms/float(npick))
   call MEDIAN(rabs, npick, rmed) ! always return median abs. residual for consistency (L1 and L2 norms)
!   ! median/mean/robomean residual instead
!    if (inorm .eq. 1) then         ! L1 NORM
!       call MEDIAN(rabs,npick,rmed)
!    else if (inorm .eq. 2) then    ! L2 NORM
!       call MEAN(resid,npick,rmed)
!    else if (inorm .eq. 3) then    ! robust L2 NORM
!      xgap=0.1
!      call ROBOMEAN2(resid,npick,xgap,10,rmed,fit2) 
!    end if
           
end subroutine DIFLOC
      
      
! DIFCLUST performs relative relocation of two clusters of events (relative to the centroid of the cluster pair)
! using the the npr "best" event pairs linking the clusters. (npr is not input explicitly,
! as the npick differential travel times, etc., include all observations across these linking pairs). 
! This is analogous to DIFLOC, which does the same thing for two events, 
! relative to their centroid... As in DIFLOC, the method uses an iterative
! ("shrinking-box") grid search approach to obtain the best (L1) relative locations.
!
!  Inputs:  qlat0  =  reference center point latitude
!           qlon0  =  reference center point longitude
!           qdep0  =  reference center point depth (km)
!           npick  =  number of differential times
!           tt     =  array (len=npick) of dif times, t2-t1 (s)
!           ip     =  array (len=npick) with phase index numbers (1 to 10) for tt data
!           slat   =  array (len=npick) with station latitudes
!           slon   =  array (len=npick) with station longitudes
!           qlat1  =  array (len=npick) of events in cluster1 latitude offsets from centroid
!           qlon1  =  array (len=npick) of events in cluster1 longitude offsets
!           qdep1  =  array (len=npick) of events in cluster1 depth offsets
!           qorg1  =  array (len=npick) of events in cluster1 time offsets
!           qlat2  =  array (len=npick) of events in cluster2 latitude offsets from centroid
!           qlon2  =  array (len=npick) of events in cluster2 longitude offsets
!           qdep2  =  array (len=npick) of events in cluster2 depth offsets
!           qorg2  =  array (len=npick) of events in cluster2 time offsets
!           TTtabfile  =  array with phase names (travel time tables...)
!           boxwid =  starting box width (km)
!           nit    =  number of iterations to perform
!           inorm  =  relocation norm (1=L1, 2=L2, 3=robust L2)
! Returns:  clat1  =  best-fitting latitude for first cluster centroid
!           clon1  =  best-fitting longitude for first cluster centroid
!           cdep1  =  best-fitting depth (km) of first cluster centroid
!           clat2  =  best-fitting latitude for second cluster
!           clon2  =  best-fitting longitude
!           cdep2  =  best-fitting depth (km)
!           cdist  =  cluster separation distance (km)
!           torgdif=  origin time difference, i.e., t2-t1 median residual (>0 when 2 is later than 1)
!           resid  =  array (len=npick) of residuals (s) between observed tdif (tt) and predicted
!           rms    =  rms residual ( sqrt( sum(resid**2)/npick) )
!           rmed   =  median absolute value residual ( median( abs(resid) ) )
!           resol  =  nominal resolution (m) of final box
!
subroutine DIFCLUST(qlat0, qlon0, qdep0, npick, tt, ip, slat, slon, &
           qlat1, qlon1, qdep1, qorg1, qlat2, qlon2, qdep2, qorg2, &
           TTtabfile, boxwid, nit, inorm, clat1, clon1, cdep1, clat2, clon2, cdep2, &
           cdist, torgdif, resid, rms, rmed, resol)
   implicit none
   integer, parameter :: npickmax = 30000
   
   integer :: npick, nit, it, iy, ix, iz, i, iflag, inorm 
   
   real :: degrad, degkm,  ddep0,  cosqlat,  ddep, fitbest, &
           fdep1, fdep2, dx, dy, delkm, tsec1, tsec2, &
           tdif, residmed, fit, fit2, fdepbest1,  &
           fdepbest2, tbest, qorg, torgdif,  cdep1, cdep2, &
           rms, rmed,  qdep0, boxwid, resol, cdist
           
   real (kind=8) :: dlat0, dlon0, dlat, dlon, flat1, flat2, flon1, flon2, flatbest1, flonbest1, &
                    flatbest2, flonbest2, clat1, clat2, clon1, clon2, qlat0, qlon0

   integer, dimension(npick) :: ip
   
   real, dimension(npick) :: tt,  qdep1, qorg1, qdep2, qorg2, resid
                             
   real (kind=8), dimension(npick) :: slat, slon, qlat1, qlon1, qlat2, qlon2
                             
   real, dimension(npickmax) :: rabs
   
   character (len=100), dimension(2) :: TTtabfile
   
   degrad = 180./3.1415927
   degkm = 111.19493

   
   if (npick > npickmax) then
      print *, '***Error in DIFCLUST: npick > 30000, npick = ', npick
      stop
   endif
   
 ! define initial search box
   dlat0 = 0.
   dlon0 = 0.
   ddep0 = 0.
   dlat = 0.5*boxwid/degkm
   cosqlat = cos(qlat0/degrad)
   dlon = dlat/cosqlat
   ddep = 0.5*boxwid

! start iterations
   do it = 1, nit

    ! grid search over box for best fit 
      fitbest = 9.e20
      
      ! get trial locations: f1, f2 (3x3 grid for this iteration's box)
      do iy = -1, 1
         flat1 = qlat0 + dlat0 + dlat*float(iy)
         flat2 = qlat0 - dlat0 - dlat*float(iy)
         do ix = -1, 1
            flon1 = qlon0 + dlon0 + dlon*float(ix)
            flon2 = qlon0 - dlon0 - dlon*float(ix)
            do iz = -1, 1
               fdep1 = qdep0 + ddep0 + ddep*float(iz)
               fdep2 = qdep0 - ddep0 - ddep*float(iz)
!               if (fdep1 < 0.) fdep1 = 0.                   !****now permit negative depths
!               if (fdep2 < 0.) fdep2 = 0.
        
        ! compute predicted travel time and residual with observed  
               do i = 1, npick
               
                  dy = flat1 + qlat1(i) - slat(i)
                  dx = (flon1 + qlon1(i) - slon(i))*cosqlat
                  delkm = sqrt(dx**2 + dy**2)*degkm
                  call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, fdep1 + qdep1(i), tsec1, iflag) 
                                   
                  dy = flat2 +qlat2(i)- slat(i)
                  dx = (flon2 +qlon2(i) - slon(i))*cosqlat
                  delkm = sqrt(dx**2 + dy**2)*degkm
                  call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, fdep2 + qdep2(i), tsec2, iflag)                  
                  
                  tdif = tsec2 + qorg2(i) - (tsec1 + qorg1(i))
                  resid(i) = tt(i) - tdif
               enddo
               
               ! compute fit (L1 or L2 norms)
               fit=0.
               if (inorm .eq. 1) then ! L1 norm
                  call MEDIAN(resid,npick,residmed)
                  do i=1,npick
                     fit=fit+abs(resid(i)-residmed)
                  enddo
               else if (inorm .eq. 2) then ! L2 norm
                  call MEAN(resid,npick,residmed)
                  do i=1,npick
                     fit=fit+(resid(i)-residmed)**2
                  end do
               else if (inorm .eq. 3) then ! robomean
                 !xgap=0.1
                 call ROBOMEAN2(resid,npick,0.1,10,residmed,fit)
                  ! note no need to iterate over resid-vector - done internally, stored in fit
               else
                   print *, 'ERROR in DIFCLUST: undefined inorm:', inorm
                   stop
               end if
               
               ! test for best fit
               if (fit < fitbest) then
                  fitbest = fit
                  flatbest1 = flat1
                  flonbest1 = flon1
                  fdepbest1 = fdep1
                  flatbest2 = flat2
                  flonbest2 = flon2
                  fdepbest2 = fdep2
                  tbest = residmed  
               end if
            enddo                   !iz loop
         enddo                   !ix loop
      enddo                   !iy loop


      dlat0 = flatbest1 - qlat0
      dlon0 = flonbest1 - qlon0
      ddep0 = fdepbest1 - qdep0
      qorg = tbest                    !origin time difference
      torgdif = qorg

      dlat = dlat*0.67      !shrink box by 2/3 each iteration
      dlon = dlon*0.67
      ddep = ddep*0.67
   
   enddo                !end it (iteration) loop

! output final, best locations
   clat1 = flatbest1
   clon1 = flonbest1
   cdep1 = fdepbest1
   clat2 = flatbest2
   clon2 = flonbest2
   cdep2 = fdepbest2
   resol = (dlat/0.67)*degkm
   
   dy =  clat2 - clat1
   dx = (clon2 - clon1)*cosqlat
   delkm = sqrt(dx**2 + dy**2)*degkm
   cdist = sqrt(delkm**2 + (cdep2 - cdep1)**2) ! output distance between cluster centroids

! output final, best misfit
    ! outputs rms residual between observed and predicted travel time, and median (abs-val) residual (MAD)
   rms = 0.
   do i = 1, npick
   
   ! cluster 1: for best-fit location, compute source-station distance and predicted travel time
      dy = clat1 + qlat1(i) - slat(i)
      dx = (clon1 + qlon1(i) - slon(i))*cosqlat
      delkm = sqrt(dx**2 + dy**2)*degkm
      call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, cdep1 + qdep1(i), tsec1, iflag)      
 
    ! cluster 2: for best-fit location, compute source-station distance and predicted travel time    
      dy = clat2 + qlat2(i) - slat(i)
      dx = (clon2 + qlon2(i) - slon(i))*cosqlat
      delkm = sqrt(dx**2 + dy**2)*degkm
      call GET_TTS_FAST8(TTtabfile(ip(i)), ip(i), delkm, cdep2 + qdep2(i), tsec2, iflag)
      
    ! predicted differential time (accounts for possible otime offset)
      tdif = tsec2 + qorg2(i) - (tsec1 + qorg1(i)) 
      
      ! update residual, rabs vectors
      resid(i) = tt(i) - tdif - qorg  ! get residual
      rms = rms + resid(i)**2 ! increment SSE sum
      rabs(i) = abs(resid(i)) ! absolute value
      
   enddo
   
   ! output RMS and median/mean absolute residual between observed and predicted travel time
   rms = sqrt(rms/float(npick)) ! SSE --> RMS
   call MEDIAN(rabs, npick, rmed) ! always return median abs. residual for consistency (L1 and L2 norms)
!   ! median/mean/robomean residual instead
!    if (inorm .eq. 1) then         ! L1 NORM
!       call MEDIAN(rabs,npick,rmed)
!    else if (inorm .eq. 2) then    ! L2 NORM
!       call MEAN(resid,npick,rmed)
!    else if (inorm .eq. 3) then    ! robust L2 NORM
!      xgap=0.1
!      call ROBOMEAN2(resid,npick,xgap,10,rmed,fit2) 
!    end if

end subroutine DIFCLUST


!-----------------------------------------------------------------------
! GROUPLOC locates a single event relative to group of events. 
! This is analogous to DIFLOC, which does the same thing for two events, 
! relative to their centroid... As in DIFLOC, the method uses an iterative
! ("shrinking-box") grid search approach to obtain the best (L1) relative locations.
!
!  Inputs:  qlat0  =  starting event latitude
!           qlon0  =  starting event longitude
!           qdep0  =  starting event depth (km)
!           npick  =  number of differential times
!           tt     =  array (len=npick) of dif times, t2-t1 (s)
!           ip     =  array (len=npick) with phase index numbers (1 to 10) for tt data
!           slat   =  array (len=npick) with station latitudes
!           slon   =  array (len=npick) with station longitudes
!           qlat   =  array (len=npick) with other event latitudes
!           qlon   =  array (len=npick) with other event longitudes
!           qdep   =  array (len=npick) with other event depths
!           qorg   =  array (len=npick) with other event origin time shifts   (***new to this version)
!           TTtabfile  =  array(len=npick) with phase names
!           boxwid =  starting box width (km)
!           nit    =  number of iterations to perform
!           inorm  =  norm to use when computing fit
!           fracsk =  shrinking box fraction
! Returns:  qlat1  =  best-fitting latitude for first event
!           qlon1  =  best-fitting longitude
!           qdep1  =  best-fitting depth (km)
!           qlat2  =  best-fitting latitude for second event
!           qlon2  =  best-fitting longitude 
!           qdep2  =  best-fitting depth (km)
!           torgdif=  origin time difference, i.e., t2-t1 median residual (>0 when 2 is later than 1)
!           resid  =  array (len=npick) of residuals (s) between observed tdif (tt) and predicted
!           rms    =  rms residual ( sqrt( sum(resid**2)/npick) )
!           rmed   =  L2/L1/robomean of the absolute value of residual
!           resol  =  nominal resolution (m) of final box
!
      subroutine GROUPLOC(qlat0,qlon0,qdep0,npick,tt,ip,slat,slon, &
           qlat,qlon,qdep,qorg,TTtabfile,boxwid,nit,inorm,fracsk,        &
           qlat1,qlon1,qdep1,torgdif,resid,rms,rmed,resol)
      implicit none

      integer :: i, iflag, inorm, it, ix, iy, iz, npick, nit
      
      real :: boxwid, cosqlat, degrad, degkm, delkm, ddep,  dx, dy, &
              fdep,  fdepbest,  fitbest, fit, fit2, &
              fracsk, fdep0,  qdep0,  qdep1,  &
              qorg1, residmed, resol, rmed, rms, tbest, tdif, tsec1, tsec2, xgap, &
              torgdif
      real (kind=8) :: dlat, dlon, flat, flon, flatbest, flonbest, flat0, flon0, &
                        qlat0, qlon0, qlat1, qlon1

      integer, dimension(npick) :: ip
      
      real, dimension(npick) :: qdep, qorg, rabs, resid, tt
      real (kind=8), dimension(npick) :: qlat, qlon, slat, slon

      character*100 TTtabfile(2)
!

      degrad = 180./3.1415927
      degkm = 111.19493

!      print *,'In the subroutine GROUPLOC ...'

    ! define initial search box
      dlat=0.5*boxwid/degkm
      cosqlat=cos(qlat0/degrad)
      dlon=dlat/cosqlat
      ddep=0.5*boxwid

      flat0=qlat0
      flon0=qlon0
      fdep0=qdep0

      flatbest=qlat0
      flonbest=qlon0
      fdepbest=qdep0

    ! start iterations
      do 100 it=1,nit

    ! start grid search over trial locations f1, f2 for best fit (3x3 grid for this iteration's box)
      fitbest=9.e20
      do 60 iy=-1,1
         flat=flat0+dlat*float(iy)
         do 50 ix=-1,1
            flon=flon0+dlon*float(ix)
            do 40 iz=-1,1
               fdep=fdep0+ddep*float(iz)
!               if (fdep.lt.0.) fdep=0.001                   !***now permit negative depths
        
        ! compute predicted travel time and residual with observed  
               do 20 i=1,npick
                  dy=flat-slat(i)
                  dx=(flon-slon(i))*cosqlat
                  delkm=sqrt(dx**2+dy**2)*degkm
                  call GET_TTS_FAST8(TTtabfile(ip(i)),ip(i),delkm,fdep,tsec1,iflag)
                  if (iflag.eq.-1) then
                     print *,'***DIFLOC warn1', ip(i),delkm, fdep
                     resid(i) = 9.
                     go to 20
                  endif
                  
                  
                  dy=qlat(i)-slat(i)                  !for speed, this block could be done once before loops
                  dx=(qlon(i)-slon(i))*cosqlat
                  delkm=sqrt(dx**2+dy**2)*degkm
                  call GET_TTS_FAST8(TTtabfile(ip(i)),ip(i),delkm, qdep(i),tsec2,iflag)
                  tsec2 = tsec2 + qorg(i)                    !***new to this version
     
     
                  if (iflag.eq.-1) then
                     print *,'***DIFLOC warn2 ', ip(i),delkm,qdep(i) 
                     resid(i) = 9.                  
                     go to 20
                  endif
                  tdif=tsec2-tsec1

                  resid(i)=tt(i)-tdif
!                  print *,'DIFLOC=',tt(i),tdif,resid(i)
20             continue

               if (inorm .eq. 1) then         ! L1 NORM
               call MEDIAN(resid,npick,residmed)
               else if (inorm .eq. 2) then    ! L2 NORM
               call MEAN(resid,npick,residmed)
               else if (inorm .eq. 3) then    ! robust L2 NORM
               xgap=0.1
               call ROBOMEAN2(resid,npick,xgap,10,residmed,fit2)
               end if

               fit=0.
               if (inorm .eq. 1) then
                  do 30 i=1,npick
                     fit=fit+abs(resid(i)-residmed)
30                continue
               else if (inorm .eq. 2) then
                  do i=1,npick
                     fit=fit+(resid(i)-residmed)**2
                  end do
               else if (inorm .eq. 3) then
                  fit=fit2
               end if

               if (fit.lt.fitbest.and.abs(residmed).lt.1.0) then
                  fitbest=fit
                  flatbest=flat
                  flonbest=flon
                  fdepbest=fdep
                  tbest=residmed  
!                  print *,'tbest=',tbest
               end if

40          continue
50       continue
60    continue

      flat0=flatbest
      flon0=flonbest
      fdep0=fdepbest 
      dlat=dlat*fracsk      !shrink box by fracsk each iteration
      dlon=dlon*fracsk
      ddep=ddep*fracsk


100   continue ! end loop over iteration

    ! output final, best locations
      qlat1=flat0
      qlon1=flon0
      qdep1=fdep0
      qorg1=tbest
      resol=(dlat/fracsk)*degkm

    ! output final, best misfit
      rms=0.
      do 120 i=1,npick
      
        ! for event 1 (alone), compute predicted travel time using best location
         cosqlat=cos(qlat1/degrad)
         dy=qlat1-slat(i)
         dx=(qlon1-slon(i))*cosqlat
         delkm=sqrt(dx**2+dy**2)*degkm
         call GET_TTS_FAST8(TTtabfile(ip(i)),ip(i),delkm,qdep1,tsec1,iflag)
         
         ! for event 2 (group) compute predicted travel time using best location
         cosqlat=cos(qlat(i)/degrad)
         dy=qlat(i)-slat(i)
         dx=(qlon(i)-slon(i))*cosqlat
         delkm=sqrt(dx**2+dy**2)*degkm
         call GET_TTS_FAST8(TTtabfile(ip(i)),ip(i),delkm,qdep(i),tsec2,iflag)
         
         ! predicted differential travel time
         tdif=tsec2-tsec1
         
         ! compute residual between observed and predicted differential times
         resid(i)=tt(i)-tdif-qorg1
         rms=rms+resid(i)**2 ! update sum of squared errors (used in rms)
         rabs(i)=abs(resid(i))
120   continue
      rms=sqrt(rms/float(npick)) ! compute RMS from SSE
      
      ! median/mean/robomean residual
      if (inorm .eq. 1) then         ! L1 NORM
      call MEDIAN(rabs,npick,rmed)
      else if (inorm .eq. 2) then    ! L2 NORM
      call MEAN(resid,npick,rmed)
      else if (inorm .eq. 3) then    ! robust L2 NORM
      xgap=0.1
      call ROBOMEAN2(resid,npick,xgap,10,rmed,fit2)
      end if
      
      torgdif = qorg1 ! origin time difference

      return
      end


!-----------------------------------------------------------------------
! ******** GET_TTS_FAST8 obtains a travel time for a seismic phase
! at a specified range and earthquake depth by interpolating
! from a file containing a table of travel times.  It differs
! from GET_TT in that a number of different phase tables can
! be called without reading them again from the files. This version modified to run faster
! by assuming x and d are evenly spaced. Negative depths are permitted.
!    Inputs:    phase  =  name of file containing travel time table
!               ip     =  index number for phase (up to 10)
!               del    =  range
!               qdep   =  earthquake depth
!    Returns:   tt     =  travel time (minutes)
!               iflag  = -1 if outside depth range
!                      =  0 for interpolation
!                      =  1 for extrapolation in range
!
subroutine GET_TTS_FAST8(phase,ip,del,qdep8,tt,iflag)
   implicit none
   
   integer, parameter :: nd0=280, nx0=401
     
   integer :: id, id1, id2, iflag, ip, ix, ix1, ix2, ixbest1, ixbest2
      
   real :: del, dfrac, qdep, t1, t2, tt, tt1, tt2, xfrac, xfrac1, xfrac2, xoff, &
           xoffmin1, xoffmin2, qdep8, tdep, velsurf
           
   integer, dimension(2) :: nd, nx
   
   real, dimension(2) :: dd, dx 
           
   character (len=100) :: phase, linebuf
   character (len=100), dimension(2) :: phaseold
      
      real d(nd0,2)
      real x(nx0,2)
      real t(nx0,nd0,2)

      save t,x,d,phaseold,nx,nd,dx,dd
!
! read file if new phase file is specified
!
      if (phase.ne.phaseold(ip)) then
         print *,'reading phase file name: ',phase(1:40)
         open (3,file=phase,status='old',err=990)
         read (3,'(a40)') linebuf   !ignore first line
         read (3,*) nx(ip),nd(ip)
         if (nx(ip).gt.nx0) then
             print *,'***GET_TTS nx truncated ',nx(ip),nx0
             nx(ip)=nx0
         end if
         if (nd(ip).gt.nd0) then
             print *,'***GET_TTS nd truncated ',nd(ip),nd0
             nd(ip)=nd0
         end if
         read (3,*) (d(id,ip),id=1,nd(ip))
         do 20 ix=1,nx(ip)
         read (3,*) x(ix,ip),(t(ix,id,ip),id=1,nd(ip))
20       continue
         close (3)
         dx(ip)=x(2,ip)-x(1,ip)
         dd(ip)=d(2,ip)-d(1,ip)
      end if
      phaseold(ip)=phase
      
! negative depth check
      if (qdep8 < 0.) then
         qdep = 0.
         velsurf = (d(2,ip) - d(1,ip))/(t(1,2,ip) - t(1,1,ip))
         tdep = abs(qdep8)/velsurf
      else
         qdep = qdep8
         tdep = 0.
      end if      
      
!
! check if outside depth range
      if (qdep.lt.d(1,ip).or.qdep.gt.d(nd(ip),ip)) then
         iflag=-1
         tt=999
         return
      end if
!
! first check to see if interpolation alone will work
      id1=1+int((qdep-d(1,ip))/dd(ip))
      id2=id1+1
      ix1=1+int((del-x(1,ip))/dx(ip))
      ix2=ix1+1

37    if (t(ix1,id1,ip).eq.0.) go to 50
      if (t(ix1,id2,ip).eq.0.) go to 50
      if (t(ix2,id1,ip).eq.0.) go to 50
      if (t(ix2,id2,ip).eq.0.) go to 50
      if (x(ix2,ip).lt.del) go to 50
      iflag=0
      xfrac=(del-x(ix1,ip))/(x(ix2,ip)-x(ix1,ip))
      t1=t(ix1,id1,ip)+xfrac*(t(ix2,id1,ip)-t(ix1,id1,ip))
      t2=t(ix1,id2,ip)+xfrac*(t(ix2,id2,ip)-t(ix1,id2,ip))
      dfrac=(qdep-d(id1,ip))/(d(id2,ip)-d(id1,ip))
      tt=t1+dfrac*(t2-t1)
      tt = tt + tdep                                  !***negative depth kludge
      return
!
! extrapolate to get tt
50    iflag=1
      xoffmin1=999.
      xoffmin2=999.
      ixbest1=999
      ixbest2=999
      do 60 ix=2,nx(ip)
         if (t(ix-1,id1,ip).eq.0) go to 55
         if (t(ix,id1,ip).eq.0) go to 55
         xoff=abs((x(ix-1,ip)+x(ix,ip))/2.-del)
         if (xoff.lt.xoffmin1) then
            xoffmin1=xoff
            ixbest1=ix
         end if
55       if (t(ix-1,id2,ip).eq.0) go to 60
         if (t(ix,id2,ip).eq.0) go to 60
         xoff=abs((x(ix-1,ip)+x(ix,ip))/2.-del)
         if (xoff.lt.xoffmin2) then
            xoffmin2=xoff
            ixbest2=ix
         end if
60    continue
      if (ixbest1.eq.999.or.ixbest2.eq.999) then
         iflag=-1
         tt=999
         return
      end if

      xfrac1=(del-x(ixbest1-1,ip))/(x(ixbest1,ip)-x(ixbest1-1,ip))
      t1=t(ixbest1-1,id1,ip)
      t2=t(ixbest1,id1,ip)
      tt1=t1+xfrac1*(t2-t1)

      xfrac2=(del-x(ixbest2-1,ip))/(x(ixbest2,ip)-x(ixbest2-1,ip))
      t1=t(ixbest2-1,id2,ip)
      t2=t(ixbest2,id2,ip)
      tt2=t1+xfrac2*(t2-t1)

      dfrac=(qdep-d(id1,ip))/(d(id2,ip)-d(id1,ip))
      tt=tt1+dfrac*(tt2-tt1)
      
      tt = tt + tdep                                  !***negative depth kludge      

      go to 999
!      
990   print *,'*** phase file not found: ',phase
      print *, 'ip:', ip
      stop
999   return
      end
   
   
