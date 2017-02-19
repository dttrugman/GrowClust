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




!--------------------------SUBROUTINES used in GROWCLUST---------------------------------!

!*************** GET_SAMPLEVEC *********************************
! Subroutine to get sampv, a random resampling vector (in "bootstrap notation", e.g. Tibshirani, 1990),
! for a data vector of length nsamp.
    ! The initial data set (sample) corresponds to a sampling vector of [1 1 1 ..... 1],
    ! whereas each resampling vector has the same sum (resampling with replacement), but
    ! with different integer weights on each entry, e.g., [0 2 1 1 0 1 1 .... ]
    ! In this notation, resampling a data vector is analogous to weighting and input data
    ! vector by the sampling vector.
!
! Inputs:
!       • nsamp: number of samples in input data vector  
!       • iseed: random number generator seed (usually 0)
! Returns:
!       • sampv: output sampling vector (indexed from 1, length nsamp)
!
! Daniel Trugman, June 2016
   
subroutine GET_SAMPLEVEC(sampv, nsamp, iseed)
      
      implicit none
      integer :: nsamp, iseed, irand, i
      integer, dimension(nsamp) :: sampv
      real :: urand
      
      ! initialize sampling vector
      do i = 1,nsamp
        sampv(i) = 0
      enddo
      
      ! get nsamp samples with replacement
      do i = 1, nsamp
        
        urand = rand(iseed) ! random number between 0 and 1
        irand = 1 + nint((nsamp-1)*urand) ! convert to int between 1 and nsamp
        
        sampv(irand) = sampv(irand)+1 ! increment sampling vector at this index
        
      enddo
      
      return  
end subroutine GET_SAMPLEVEC

!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------    
!*************** RESAMPLE_DATA *********************
! Subroutine to apply sampv, a random resampling vector (in "bootstrap notation", e.g. Tibshirani, 1990),
! to an input data vector, data_in, to obtain an output resampled data vector, data_out.
    ! The initial data set (sample) corresponds to a sampling vector of [1 1 1 ..... 1],
    ! whereas each resampling vector has the same sum (resampling with replacement), but
    ! with different integer weights on each entry, e.g., [0 2 1 1 0 1 1 .... ]
    ! In this notation, resampling a data vector is analogous to weighting the input data
    ! vector by the sampling vector.
! Note, there are separate subroutines for INT*4, REAL*4, and REAL*8 data types
!
! Inputs:
!       • data_in: input data vector, length nsamp, indexed from 1
!       • sampv: input sampling vector (this is reset on output)
!       • nsamp: number of samples in input data vector  

! Returns:
!       • data_out: resampled data vector, length nsamp, indexed from 1
!
! Daniel Trugman, June 2016
!  
subroutine RESAMPLE_IDATAVEC(data_in, data_out, sampv, nsamp) ! integer data
      
      implicit none
      integer :: i,j,k, nsamp
      integer, dimension(nsamp) :: sampv
      integer, dimension(nsamp) :: data_in, data_out
      
      k = 1 ! index in output vector
      do i = 1,nsamp ! loop over input data
      
        do j = 1, sampv(i) ! copy sampv(i) values to data_out
            data_out(k) = data_in(i)
            k = k + 1
        enddo
        
      enddo
      
      return  
end subroutine RESAMPLE_IDATAVEC

subroutine RESAMPLE_FDATAVEC(data_in, data_out, sampv, nsamp) ! real*4 data
      
      implicit none
      integer :: i,j,k, nsamp
      integer, dimension(nsamp) :: sampv
      real, dimension(nsamp) :: data_in, data_out
      
      k = 1 ! index in output vector
      do i = 1,nsamp ! loop over input data
      
        do j = 1, sampv(i) ! copy sampv(i) values to data_out
            data_out(k) = data_in(i)
            k = k + 1
        enddo
        
      enddo
      
      return  
end subroutine RESAMPLE_FDATAVEC

subroutine RESAMPLE_DFDATAVEC(data_in, data_out, sampv, nsamp) ! real*8 data
      
      implicit none
      integer :: i,j,k, nsamp
      integer, dimension(nsamp) :: sampv
      real*8, dimension(nsamp) :: data_in, data_out
      
      k = 1 ! index in output vector
      do i = 1,nsamp ! loop over input data
      
        do j = 1, sampv(i) ! copy sampv(i) values to data_out
            data_out(k) = data_in(i)
            k = k + 1
        enddo
        
      enddo
      
      return  
end subroutine RESAMPLE_DFDATAVEC
    
!---------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------    
! **** DPMEAN_STDDEV returns both the mean xmean and std. dev. of input array x (length n). ****
!                   Input arrays assumed to be double precision.
!----------------------------------------------------------------------------------------------

subroutine DPMEAN_STDDEV(x,n,xmean, xstd)
      implicit none

      integer i
      integer n
      real*8 sum
      real*8 xmean
      real*8 xstd
      real*8 x(n)

     ! check for length zero vector
      if (n.eq.0) then
         xmean=0.
         xstd=0.
         return
      endif   
      ! compute mean   
      sum=0.
      do i=1,n
         sum=sum+x(i)
      enddo
      xmean=sum/real(n)
      
      ! check for length one vector
      if (n.eq.1) then
        xstd=0.
        return
      endif
      
      ! compute standard deviation (sample, N-1 normalization)
      sum = 0.
      do i=1,n
        sum = sum + (x(i) - xmean)**2
      enddo
      xstd = dsqrt(sum/real(n-1))
      
      return
end subroutine DPMEAN_STDDEV

!-----------------------------------------------------------------------
!********* TIMER: subroutine to keep track of program run time *********
!-----------------------------------------------------------------------
      subroutine TIMER
      implicit none

      integer if1
  
      real diftim
      real etime
      real oldtim
      real tnow
      real tarray(2)

      save oldtim
      data if1/1/
!

      if(if1.eq.0) goto 5
      if1=0
      oldtim=etime(tarray)
      return
    5 continue
      tnow=etime(tarray)
      diftim=tnow-oldtim
      oldtim=tnow
      print *,'elapsed seconds =',diftim
      return
      end subroutine TIMER
!------------------------------------------------------------------------!


!-------------------------------------------------------------------------
!******* MEAN returns the mean xmean of input array x (length n). ********
!-------------------------------------------------------------------------
      subroutine MEAN(x,n,xmean)
      implicit none

      integer i
      integer n

      real sum
      real xmean
      real x(n)
!

      if (n.eq.0) then
         xmean=0.
         return
      end if
      sum=0.
      do 10 i=1,n
         sum=sum+x(i)
10    continue
      xmean=sum/float(n)
      return
      end
      
! ------ double precision version
      subroutine DPMEAN(x,n,xmean)
      implicit none

      integer i
      integer n

      real*8 sum
      real*8 xmean
      real*8 x(n)
!

      if (n.eq.0) then
         xmean=0.
         return
      end if
      sum=0.
      do 10 i=1,n
         sum=sum+x(i)
10    continue
      xmean=sum/float(n)
      return
      end 

!---------------------------------------------------------------------------------------------------
!**** MEDIAN returns the median xmed of input array x (length n). Modified from on Numerical Recipes
!---------------------------------------------------------------------------------------------------
      subroutine MEDIAN(x,n,xmed)
      
      implicit none
      integer i
      integer n
      integer n2
      real xmed
      real x(n)
      real x2(n)

      do i=1,n
        x2(i)=x(i)
      enddo
      call PIKSRT(n,x2)
      n2=n/2
      if (2*n2.eq.n) then
         xmed=0.5*(x2(n2)+x2(n2+1))
      else
         xmed=x2(n2+1)
      end if
      return
      end
      
      ! double precision I/O (single precision copy passed in PIKSRT, which is probably OK)
      subroutine DPMEDIAN(x,n,xmed)
      
      implicit none
      integer i
      integer n
      integer n2
      
      real*8 x(n)
      real*8 xmed
      real x2(n)

      do i=1,n
        x2(i)=x(i)
      enddo
      call PIKSRT(n,x2)
      n2=n/2
      if (2*n2.eq.n) then
         xmed=0.5*dble((x2(n2)+x2(n2+1)))
      else
         xmed=dble(x2(n2+1))
      end if
      return
      end
      
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
!*****  PIKSRT sorts array x, length n. Modified from on Numerical Recipes ******
!--------------------------------------------------------------------------------
      subroutine PIKSRT(n,x)
      implicit none
      integer i
      integer j
      integer n
      real x(n)
      real a

      do 12 j=2,n
         a=x(j)
         do 11 i=j-1,1,-1
            if (x(i).le.a) go to 10
            x(i+1)=x(i)
11       continue
         i=0
10       x(i+1)=a
12    continue
      return
      end
!----------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------      
!**** ROBOMEAN2 returns the robust mean (L1/L2 hybrid) xmean of input array x (length n), 
!      (estimated using an IRLS-type iterative weighting scheme that transitions smoothly
!       between L2/mean weights for normal data pts and L1/median weights for outliers)
!
! Inputs:
!       • x,n: input data vector
!       • n:   length of data vector
!       • xgap: distance/weighting parameter (max weight is 1/xgap when x(i)=xmean)
!                         d=xgap+abs(x(i)-xmean)
!                         xw(i)=1./d
!       • nit: number of iterations

! Returns:
!       • xmean: weighted mean
!       • fit2 :  robust L2 norm 
!----------------------------------------------------------------------------------------
      subroutine ROBOMEAN2(x,n,xgap,nit,xmean,fit2)
      implicit none

      integer :: i, it, n, nit
      real :: d, fit2, sum, sumw, xgap, xmean
      real :: x(n), xw(n)
      
      ! for robustness
      if (n.eq.0) then
         xmean=0.
         return
      end if
      
      ! set initial weights to zero
      do i=1,n
         xw(i)=1.
      enddo
      
      do it=1,nit ! ---- iterations loop
         
         ! compute weighted data sum and sum of weights
         sum=0.
         sumw=0.
         do i=1,n
            sum=sum+x(i)*xw(i) ! weighted data sum
            sumw=sumw+xw(i) ! sum of weights
         enddo
         xmean=sum/sumw ! weighted mean (note that in limit that xw=1, sumW=n)
         
         ! compute weights for next iteration
         if (it.ne.nit) then  
            do i=1,n
               d=xgap+abs(x(i)-xmean) ! "distance"
               xw(i)=1./d             ! weights
            enddo
         endif
      enddo

      ! compute weighted misfit (weighted L2 norm)
      fit2=0.0
      do i=1,n
         fit2=fit2+xw(i)*(x(i)-xmean)**2 
      end do
!      fit2=sqrt(fit2/float(n))

      return
      end
!------------------------------------------------------------------

!------------------------------------------------------------------
! **** INDEXX returns an array of indices INDX, sorted low to high
! ****** based on values in ARRIN 
!------------------------------------------------------------------
      SUBROUTINE INDEXX(N,ARRIN,INDX)
      DIMENSION ARRIN(N),INDX(N)
      DO 11 J=1,N
        INDX(J)=J
11    CONTINUE
      L=N/2+1
      IR=N
10    CONTINUE
        IF(L.GT.1)THEN
          L=L-1
          INDXT=INDX(L)
          Q=ARRIN(INDXT)
        ELSE
          INDXT=INDX(IR)
          Q=ARRIN(INDXT)
          INDX(IR)=INDX(1)
          IR=IR-1
          IF(IR.EQ.1)THEN
            INDX(1)=INDXT
            RETURN
          ENDIF
        ENDIF
        I=L
        J=L+L
20      IF(J.LE.IR)THEN
          IF(J.LT.IR)THEN
            IF(ARRIN(INDX(J)).LT.ARRIN(INDX(J+1)))J=J+1
          ENDIF
          IF(Q.LT.ARRIN(INDX(J)))THEN
            INDX(I)=INDX(J)
            I=J
            J=J+J
          ELSE
            J=IR+1
          ENDIF
        GO TO 20
        ENDIF
        INDX(I)=INDXT
      GO TO 10
      END
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------
! **** GETSTD: program to calculate the standard deviation of an array ****
!--------------------------------------------------------------------------
      subroutine GETSTD(x,ndata,xstd)
      implicit none

      integer nmax
      parameter (nmax=100)
      integer ndata
      integer ii

      real x(nmax)
      real x0(nmax)
      real xstd
      real xmn

      if(ndata.gt.nmax) then
         print *,'Increase nmax !',nmax,ndata
         stop
      endif

      do ii=1,ndata
         x0(ii)=x(ii)
      enddo
      call MEAN(x0,ndata,xmn)

      xstd=0.
      do ii=1,ndata
         xstd=xstd+(x0(ii)-xmn)**2
      enddo
      xstd=sqrt(xstd/real(ndata))

      return
      end
!-----------------------------------------------------------

!-----------------------------------------------------------
! **** GETMAD: program to calculate the MAD of an array ****
!-----------------------------------------------------------
      subroutine GETMAD(x,ndata,xmad)
      implicit none
    
      integer ndata
      integer ii

      real x(ndata)
      real x0(ndata)
      real x2(ndata)
      real xmad
      real xmed

      do ii=1,ndata
         x0(ii)=x(ii)
      enddo
      call MEDIAN(x0,ndata,xmed)

      do ii=1,ndata
         x2(ii)=abs(x0(ii)-xmed)
      enddo
      call MEDIAN(x2,ndata,xmad)

      return
      end
      
      ! double precision version
      subroutine GETDPMAD(x,ndata,xmad)
      implicit none
    
      integer ndata
      integer ii

      real*8 x(ndata)
      real*8 x0(ndata)
      real*8 x2(ndata)
      real*8 xmad
      real*8 xmed

      do ii=1,ndata
         x0(ii)=x(ii)
      enddo
      call DPMEDIAN(x0,ndata,xmed)

      do ii=1,ndata
         x2(ii)=abs(x0(ii)-xmed)
      enddo
      call DPMEDIAN(x2,ndata,xmad)

      return
      end

!-----------------------------------------------------------

