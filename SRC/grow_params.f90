!-----------------------------------------------------------------------
! Copyright 2022 Daniel Trugman
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
   integer, parameter :: nsta0=1000            !max number of stations
   integer, parameter :: nq0=100000           !max number of quakes
   integer, parameter :: npair0=2000000        !max total number of event pairs
   integer, parameter :: ntmax=100000           !maximum number of trees (clusters), normally should be same as nq0
   integer, parameter :: nbmax=100000           !max number of events (branches) per tree, should be of order nq0
   integer, parameter :: n0=1000              !max number of differential times for each event pair
   integer, parameter :: ndif0=15000000        !max total number of differential times
   integer, parameter :: n08 = 10000          !max number of diff. times for 10 event pairs
   integer, parameter :: maxboot = 100        ! max number of bootstrap resamples
   integer, parameter :: maxevid = 100000000 ! maximum event id number

    ! ------- GrowClust algorithm control parameters -------------------------------
   real, parameter    :: distmax = 5.0          ! maximum catalog(input) distance to join clusters (km)
   real, parameter    :: distmax2 = 3.0         ! maximum relocated distance to join clusters (km)
   integer, parameter :: nclustshiftmin = 1    ! minimum number of events in cluster to apply cluster shift test
   real, parameter    :: hshiftmax = 2.0        ! maximum permitted horizontal cluster shifts (km)
   real, parameter    :: vshiftmax = 2.0        ! maximum permitted vertical cluster shifts (km)
   real, parameter    :: rmedmax = 0.05         ! maximum median absolute tdif residual to join clusters
   
   ! ------- Relative Relocation subroutine parameters -------------
   real, parameter    :: boxwid = 3. ! initial "shrinking-box" width (km)
   integer, parameter :: nit = 15 ! number of iterations
   integer, parameter :: irelonorm = 1 ! relocation norm (L1 norm=1, L2 norm=2, 3=robust L2)
   real, parameter    :: tdifmax = 30. ! maximum differential time value allowed (for error-checking on xcor data input)
   
   ! -------- Bootstrap resampling parameters -------------------
   integer, parameter  :: iseed = 0 ! random number seed
   real, parameter     :: rms_nan = 0.0 ! rms flag in output file for unrelocated event
   real, parameter     :: MAD_nan = -1.0 ! bootstrap MAD error flag in output file for unrelocated events
   real, parameter     :: SE_nan = MAD_nan ! bootstrap SE error flag in output file for unrelocated events
   integer, parameter  :: samp_type = 2 ! bootstrap sampling type
        !  1: Resample each event pair independently (each pair always has the same # of picks in each resample)'
        !  2: Resample the entire data vectors at once (event pairs may have different # of picks in each resample)
        
   ! ------- Velocity model parameters (added 04/2018) -------------------------
   integer, parameter  :: vzmodel_type = 1 ! velocity model type: 1 = flat earth, (Z,Vp,Vs)
                                           !                  or  2 = radial, (R,Vp,Vs): 
                                           !                           note: option 2 has not been extensively tested!!!
   
   end MODULE grow_params
