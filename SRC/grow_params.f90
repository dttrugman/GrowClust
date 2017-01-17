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



! ------ This module module initializes fixed parameters used by GrowClust
! ----------- (these are not read from input file, so modify and recompile the module as necessary) -----!  

 
   MODULE grow_params
   
   ! ----------- Array size parameters ----------------------------------
   integer, parameter :: nsta0=500            !max number of stations
   integer, parameter :: nq0=80000           !max number of quakes
   integer, parameter :: npair0=2000000        !max total number of event pairs
   integer, parameter :: ntmax=80000           !maximum number of trees (clusters), normally should be same as nq0
   integer, parameter :: nbmax=10000           !max number of events (branches) per tree
   integer, parameter :: n0=1000              !max number of differential times for each event pair
   integer, parameter :: ndif0=9999999        !max total number of differential times
   integer, parameter :: n08 = 10000          !max number of diff. times for 10 event pairs
   integer, parameter :: maxboot = 300        ! max number of bootstrap resamples

    ! ------- GrowClust algorithm control parameters -------------------------------
   real, parameter    :: conparam = 0.01        ! minimum connection fraction to join clusters
   real, parameter    :: distmax = 5.0          ! maximum catalog(input) distance to join clusters (km)
   real, parameter    :: distmax2 = 3.0         ! maximum relocated distance to join clusters (km)
   integer, parameter :: nclustshiftmin = 10    ! Minimum number of events in cluster to apply cluster shift test
   real, parameter    :: hshiftmax = 2.0        ! maximum permitted horizontal cluster shifts (km)
   real, parameter    :: vshiftmax = 2.0        ! maximum permitted vertical cluster shifts (km)
   real, parameter    :: rmedmax = 0.05         ! maximum median absolute tdif residual to join clusters
   
   ! ------- Relative Relocation subroutine parameters -------------
   real, parameter    :: boxwid = 3. ! initial "shrinking-box" width (km)
   integer, parameter :: nit = 20 ! number of iterations
   integer, parameter :: irelonorm = 1 ! relocation norm (L1 norm=1, L2 norm=2, 3=robust L2)
   
   ! -------- Bootstrap resampling parameters -------------------
   integer, parameter  :: iseed = 0 ! random number seed
   real, parameter     :: rms_nan = 0.0 ! rms flag in output file for unrelocated event
   real, parameter     :: MAD_nan = -1.0 ! bootstrap MAD error flag in output file for unrelocated events
   real, parameter     :: SE_nan = MAD_nan ! bootstrap SE error flag in output file for unrelocated events
   integer, parameter  :: samp_type = 2 ! bootstrap sampling type
        !  1: Resample each event pair independently (each pair always has the same # of picks in each resample)'
        !  2: Resample the entire data vectors at once (event pairs may have different # of picks in each resample)
   
   end MODULE grow_params