!-----------------------------------------------------------------------
! Copyright 2020 Daniel Trugman
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

! VZFILLIN is designed to read (z, Vp, Vs) models and resamples them at any
! desired finer depth interval.  This is done with linear interpolation.
!
!  Inputs:  infile  = name of input velocity model (z, Vp, Vs) 
!           outfile = name of output, intepolated model (z, Vp, Vs)
!           sfact   = Vp/Vs ratio if Vs is unlisted (set to 0.0)
!           dz      = interpolation distance
!  Outputs: zmin    = minimum depth in input model
!           zmax    = maximum depth in input model (also used to specify minimum bottom depth)
!
      subroutine vzfillin(infile, outfile, sfact, dz, zmin, zmax)
      
      implicit none
      character*100 infile,outfile
      integer npts0
      parameter (npts0=5000)
      real :: d(npts0),v(npts0,2),v0(2),buf(npts0,3)
      real :: sfact,dz,z,z1,z2,fact,zmin,zmax
      integer :: i,j,k,kk,npts,nk,iz,nz

      open (11,file=infile,status='old')
      open (12,file=outfile)

! read input file (depth, Vp, Vs). If Vs .eq. 0.0, Vs = Vp/sfact 
      do i=1,100
         read (11,*,end=21) d(i),(v(i,j),j=1,2)
         if (v(i,2).eq.0.) v(i,2)=v(i,1)/sfact
      enddo
21    npts=i-1
      close (11)

! edited (2/18/17) to put an extra depth point at zmax km, if vzmodel isn't this deep
!   this has the effect of making the bottom-most point specify the top of a constant layer from
!    d(npts) to zmax 
      if (d(npts) < zmax) then
        npts = npts + 1
        d(npts) = zmax
        v(npts,1) = v(npts-1,1)
        v(npts,2) = v(npts-1,2)
      endif
! edited (2/18/17) to return zmin and zmax
        zmin = d(1)
        zmax = d(npts)
        
! interpolate the input model, while preserving layer interfaces
     kk=0
     nz = floor(d(npts)/dz)+1 ! number of interpolation pts
     do iz=1,nz
         z = (iz-1)*dz ! interpolation point
         do i=1,npts-1 ! loop over velocity model layers
            z1=d(i) ! top layer
            z2=d(i+1) ! bottom layer
            if (z1.le.z.and.z2.ge.z) then ! match with interpolation pt 
               if (z1.eq.z2) then ! layer interface
                  fact=0.
               else               ! not an interface: interpolate
                  fact=(z-z1)/(z2-z1)
               end if
               do j=1,2
                  v0(j)=v(i,j)+fact*(v(i+1,j)-v(i,j)) ! compute v at interpolation point
               enddo
               ! format output buffer
               kk=kk+1
               buf(kk,1)=z
               buf(kk,2)=v0(1)
               buf(kk,3)=v0(2)               
33             format (f10.3,2f8.4)
               do k=i+1,npts
                  if (d(k)-z.lt.dz) then
                     kk=kk+1
                     buf(kk,1)=d(k)
                     buf(kk,2)=v(k,1)
                     buf(kk,3)=v(k,2)
                  end if
               enddo
               exit ! exit inner loop
            end if
         enddo ! end loop on velocity model layers
     enddo ! end loop on interpolation files

! write to output file
      nk=kk
      write (12,33) (buf(1,j),j=1,3)
      do k=2,nk
         if (buf(k,1).ne.buf(k-1,1).or. &
            buf(k,2).ne.buf(k-1,2).or. &
            buf(k,3).ne.buf(k-1,3)) then
           write (12,33) (buf(k,j),j=1,3)
         end if
      enddo
      close (12)
      return
      end
      
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------
! DEPTABLE creates tables of seismic travel times as a function of distance
! and source depth.  These tables are evenly spaced and designed to be read
! with the GET_TTS subroutines.
!
!    Inputs:    vmodel   = filename of velocity model (interpolated to spacing dep3 or finer)
!               iw       = phase index: 1 = P, 2 = S
!               plongcut = minimum ray parameter
!               dep1     = min dep for table
!               dep2     = max dep for table
!               dep3     = dep spacing for table
!               de11     = min del for table
!               de12     = max del for table
!               de13     = del spacing for table
!               ideprad  = velocity model fmt, first column of input:  1=depth, 2=radius (added 04/2018)
!               ttfile   = name of output travel time table file
!
!-----------------------------------------------------------------------
      subroutine deptable(vmodel,iw,plongcut,dep1,dep2,dep3,del1,del2,del3,ideprad,ttfile)
      
      implicit none

! parameters
! npts0 = maximum number of lines in input velocity model
! nz0 = maximum number of depths in desired output tables
! nx0 = maximum number of distances is desired output tables
! nray0 = maximum number of rays during ray tracing
! ncount0 = 2*nray0
      integer npts0,nz0,nx0,nray0,ncount0
      parameter (npts0=1000,nz0=201,nx0=501,nray0=40002,ncount0=80002)

      real erad,pi
      parameter (erad=6371., pi=3.1415927)

      integer npts,ndep,nray,ncount,ndel
      integer ideptype,idep,itype,icount,idel,iw,imth,irtr
      integer i,j,i2,ip,ideprad

      real ecircum,kmdeg,degrad
      real p,pmin,pmax,pstep,plongcut,zmax,frac,h
      real dep,dep1,dep2,dep3,del,del1,del2,del3,deldel
      real scr1,scr2,scr3,scr4,xold,x,x1,x2,dx,t,t2,dt
      real tbest
      real xsave(ncount0),tsave(ncount0),psave(ncount0)
      real deptab(nz0),ptab(nray0),delttab(nray0)
      real z(npts0),alpha(npts0),beta(npts0)
      real z_s(npts0),alpha_s(npts0),beta_s(npts0)
      real slow(npts0,2),deltab(nray0),tttab(nray0)
      real tt(nx0,nz0)
      real depxcor(nray0,nz0),deptcor(nray0,nz0)

      character*100 vmodel
      character*100 ttfile

!-----------------------------------------------------------------------
!     define fixed variables
      ecircum=2.*pi*erad
      kmdeg=ecircum/360.
      degrad=180./pi
      
      zmax = 9999 ! maximum depth
      nray = 10000 ! number of rays
      ideptype = 1 ! Source depths:  (1) Range, (2) Exact
      itype = 1 ! set TT table units: (1) del in km, t in sec  or  (2) del in deg., t in min

!---------------- Get initial velocity model input ---------------------------!

! read velocity model and transform to flat earth
!          value at center of earth is removed to avoid singularity in transformation
      open (7,file=vmodel,status='old')
      do i=1,npts0
         read (7,*,end=30) z_s(i),alpha_s(i),beta_s(i)
         if (ideprad.eq.2) z_s(i)=erad-z_s(i) ! note erad is earth radius
         if (z_s(i).eq.erad) go to 30
         call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i)) ! flat-earth transform for Vp
         call FLATTEN(z_s(i),beta_s(i),z(i),beta(i)) ! flat-earth transform for Vs
      enddo   
      print *,'***',npts0,' point maximum exceeded in model'
      stop
30    close (7)
      print *,'finished reading model'

!set up dummy interface at bottom
      z_s(i)=z_s(i-1)                  
      alpha_s(i)=alpha_s(i-1)
      beta_s(i)=beta_s(i-1)
      call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i))
      call FLATTEN(z_s(i),beta_s(i),z(i),beta(i))
      npts=i
      print *,'Depth points in model= ',npts

! compute slowness
      do i=1,npts
         slow(i,1)=1./alpha(i)
         if (beta(i).ne.0.) then
            slow(i,2)=1./beta(i)
         else
            slow(i,2)=1./alpha(i)              !fluid legs are always P!
         end if       
      enddo

! print out velocity table
      print *,'************************* Table of Model Interfaces *****', &
     '*****************'
      print *,' Depth  Top Velocities  Bot Velocities    -----Flat Earth ', &
     ' Slownesses-----'
      print *,'             vp1  vs1        vp2  vs2       p1      p2  ', &
          '      s1      s2'

      do i=2,npts
         if (i.eq.2.or.z(i).eq.z(i-1)) then
            scr1=1./alpha(i-1)
            scr2=1./alpha(i)
            scr3=999.
            if (beta(i-1).ne.0.) scr3=1./beta(i-1)
            scr4=999.
            if (beta(i).ne.0.) scr4=1./beta(i)
            print 42,z_s(i),i-1,alpha_s(i-1),beta_s(i-1), &
              i,alpha_s(i),beta_s(i), scr1,scr2,scr3,scr4
42          format (f6.1,2(i5,f6.2,f5.2),2x,2f8.5,2x,2f8.5)
         end if
      enddo

!-----------------------------------------------------------------------

! ------------------- More input -------------
      
!  get range of source depths
      ndep = floor((dep2+dep3/10.-dep1)/dep3)+1 ! (note, this value has been checked on input)
      do idep = 1,ndep
        dep = dep1 + (idep-1)*dep3
        deptab(idep)=dep
      enddo

! get number of rays to compute     
      pmin=0.
      pmax=slow(1,iw)
      print *,'pmin, pmax = ', pmin, pmax
      print *,'Number of rays to compute:'
      print *, nray
      pstep=(pmax-pmin)/float(nray)
       

!-----------------------------------------------------------------------      
      

! --------------------------------- ray tracing -------------------------
!    shoot rays with different ray parameters (nray total) from the surface,
!    and keep track of the offset distance (x) and travel time (t) to different depths
200   do ip = 1, nray ! ------- loop over ray parameters (p) ---------------
      
         ! current ray parameter
         p = pmin + (ip-1)*pstep
         ptab(ip)=p

         ! rays start at x,t = 0,0
         x=0.
         t=0.
         
         !preferred v(z) interpolation method, optimal for flat-earth transform
         imth=3
         
         ! initialize arrays: depxcor, deptcor, depucor (size nray by ndep), which track
         ! the offset (x) and travel time (t) for different rays to different depths 
         do idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(ip,idep)=0.
               deptcor(ip,idep)=0.
            else
               depxcor(ip,idep)=-999.
               deptcor(ip,idep)=-999.
            end if
        enddo

         do i=1,npts-1 ! ------ loop over layers (i) ----------------------------
         
             !check to see if z exceeds zmax
             if (z_s(i).ge.zmax) then                          
                deltab(ip)=-999.
                tttab(ip)=-999.
                go to 200
             end if

             ! layer thickness
             h=z(i+1)-z(i)                           
             if (h.eq.0.) cycle    !skip if interface
             
            ! LAYERTRACE calculates the travel time and range offset for ray tracing through a single layer.
            !  Input:   p     =  horizontal slowness
            !           h     =  layer thickness
            !           utop  =  slowness at top of layer
            !           ubot  =  slowness at bottom of layer
            !           imth  =  interpolation method
            !                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)  fastest to compute
            !                         = 2,  v(z) = a - b*z            linear gradient
            !                         = 3,  v(z) = a*exp(-b*z)        referred when Earth Flattening is applied
            !  Returns: dx    =  range offset
            !           dt    =  travel time
            !           irtr  =  return code (-1: zero thickness layer, 0: ray turned above layer, 
            !                 =    1: ray passed through layer, 2: ray turned within layer, 1 segment counted)
             call LAYERTRACE(p,h,slow(i,iw),slow(i+1,iw),imth,dx,dt,irtr) ! compute dx, dt for layer
             
             ! update x,t after tracing through layer
             x=x+dx
             t=t+dt
             
             ! exit when ray has turned
             if (irtr.eq.0.or.irtr.eq.2) exit  
         
            ! save current x,t,u for ray sampling source depths (stored in deptab)
             do idep=1,ndep
                if (abs(z_s(i+1)-deptab(idep)).lt.dep3/10.) then
                   depxcor(ip,idep)=x
                   deptcor(ip,idep)=t
                end if
             enddo

!   
         enddo !------------ end loop on layers----------------------------------
         
         ! compute final (surface-to-surface) two-way offset and travel times for ray
         deltab(ip)=2*x                   !stored in km
         tttab(ip)=2*t                    !stored in seconds

      enddo             !---------------- end loop on ray parameter p --------------------
!      print *,'Completed ray tracing loop'

!----------------------------------------------
      
! special section outside loop for p = pmax for (0,0 point)
         nray=nray+1
         ip = nray
         ptab(ip)=slow(1,iw)
         deltab(ip)=0.
         tttab(ip)=0.
         do idep=1,ndep
            if (deptab(idep).eq.0.) then ! surface ray hitting (0,0)
               depxcor(ip,idep)=0.
               deptcor(ip,idep)=0.
            else                         ! surface ray does not reach this depth
               depxcor(ip,idep)=-999.
               deptcor(ip,idep)=-999.
            end if         
         enddo

!-------------------------------------------------------------------------


!-----  Now compute T(X,Z) for first arriving rays
      
      do idep=1,ndep ! loop over source depth's (Z)
         
         
         icount=0 ! save array index
         xold=-999. ! current x
         
         ! if source is at 0 depth, go skip the upgoing ray loop
         if (deptab(idep).eq.0.) then
            i2=nray ! start next loop here
            go to 223
         end if
         
         ! loop for upgoing rays from the source
         do i=1,nray                           ! loop from steep to shallow angles
            x2=depxcor(i,idep)               ! offset at this source depth
            if (x2.eq.-999.) exit            ! stop when run out of x's
            if (x2.le.xold) exit             ! stop when ray heads inward
            t2=deptcor(i,idep)               ! travel time to this depth
            icount=icount+1                  ! increment save index
            xsave(icount)=x2                 ! save offset from this depth to surface
            tsave(icount)=t2                 ! save travel time
            psave(icount)=-ptab(i)           ! save p as negative for upgoing from source
            xold=x2                          ! update x value
         enddo
         i2=i-1 ! start next loop here
         
         ! loop for downgoing rays from the source
223      continue         
         do i=i2,1,-1                       ! loop from shallow to steep angles
            if (depxcor(i,idep).eq.-999.) cycle ! skip
            if (deltab(i).eq.-999.) cycle       ! skip
            x2=deltab(i)-depxcor(i,idep)    ! source-surface offset is total offset minus offset from downgoing leg 
            t2=tttab(i)-deptcor(i,idep)     ! same for source-surface travel time
            icount=icount+1                 ! increment save index
            xsave(icount)=x2                ! save offset from this depth to surface
            tsave(icount)=t2                ! save p as postive for downgoing from source
            psave(icount)=ptab(i)           ! save p as positve for downgoing from source
         enddo

         ! total count of data points saved
         ncount=icount
         
         
         ! interpolate offsets to the desired spacing and find the first-arriving ray
         ndel = floor((del2+del3/10.-del1)/del3) + 1 ! (note, this value has been checked on input)
         
         do idel = 1, ndel !---------- loop over offsets
            
            ! get current offset
            deldel = del1 + (idel-1)*del3 ! current offset in km
            delttab(idel)=deldel
            del=deldel
            if (itype.eq.2) del=deldel*kmdeg ! convert from km to degree, if desired

            ! initialize travel time to a high value
            tt(idel,idep)=99999.
            
            ! search for first arriving ray at this offset
            do i=2,ncount

               ! get adjacent x values
               x1=xsave(i-1)
               x2=xsave(i)

               ! check to see if they span the desired del
               if (x1.gt.del.or.x2.lt.del) cycle

               ! # skip downgoing rays (p > 0) with this incidence: from refracted rays
               if (psave(i).gt.0..and.psave(i).lt.plongcut) cycle
               
               ! compute tbest
               frac=(del-x1)/(x2-x1)
               tbest=tsave(i-1)+frac*(tsave(i)-tsave(i-1))
              
               ! check to see if this is the best arrival so far
               if (tbest.lt.tt(idel,idep)) then
                  tt(idel,idep)=tbest
               end if
            enddo
            
            ! no ray arrivals
            if (tt(idel,idep).eq.99999.) tt(idel,idep)=0.    
            if (itype.eq.2) tt(idel,idep)=tt(idel,idep)/60.            

         enddo                                    !end loop on offsets
                  
      enddo                                        !end loop on depth

! fix potential edge case
if ( (delttab(1).eq.0.) .and. (deptab(1).eq.0.) ) then
    tt(1,1)=0.                 !set tt to zero at (0,0)
end if

!-------------- Make output files -----------------

! get file name from input
      print *,'Output file name for travel-time table:'
      print *, ttfile
      open (11,file=ttfile)     
      
! write headers for all files
     
! first line: model name, iw (1/2 for P/S), pmin, pmax, nray
      write (11,407) vmodel(1:20),iw,pmin,pmax,nray
407    format ('From deptable, file= ',a20,' iw =',i2,' pmin=',f8.5, &
              ' pmax=',f8.5,' nray=',i6)
     
! second line: table size ndel,ndep:
!       ndel rows (different X/DEL offsets)
!       ndep columns (different source depths)
      write (11,408) ndel,ndep
408      format (2i5)
      
!     third line: row of source depths
      write (11,409) (deptab(j),j=1,ndep)
409      format (8x,100f9.1)
      
!   fill in table:
!       first column is X/DEL of each row
!       then TT_ij = TT(x_i, z_j), where: TT_ij is the travel time of first arriving ray 
!       from a source at horizontal distance x_i and depth z_j
      do i=1,ndel
         if (itype.eq.1) then
            write (11,410) delttab(i),(tt(i,j),j=1,ndep)
         else
            write (11,413) delttab(i),(tt(i,j),j=1,ndep)
         end if

410      format (101f9.4)
413      format (f9.3,100f9.4)

      enddo
      close (11)

      end

!-----------------------------------------------------------------------
! INTERP finds the y3 value between y1 and y2, using the
! position of x3 relative to x1 and x2.
      subroutine INTERP(x1,x2,x3,y1,y2,y3)
      fract=(x3-x1)/(x2-x1)
      y3=y1+fract*(y2-y1)
      return
      end

!-----------------------------------------------------------------------
! FLATTEN calculates flat earth tranformation.
      subroutine FLATTEN(z_s,vel_s,z_f,vel_f)
      erad=6371.
      r=erad-z_s
      z_f=-erad*alog(r/erad)
      vel_f=vel_s*(erad/r)
      return
      end

!-----------------------------------------------------------------------
! UNFLATTEN is inverse of FLATTEN.
      subroutine UNFLATTEN(z_f,vel_f,z_s,vel_s)
      erad=6371.
      r=erad*exp(-z_f/erad)
      z_s=erad-r
      vel_s=vel_f*(r/erad)
      return
      end

!
!-----------------------------------------------------------------------
! LAYERTRACE calculates the travel time and range offset
! for ray tracing through a single layer.
!
! Background:
!  - Total slowness u = 1 / v
!  - Horizontal slowness = p = u(z)*sin(theta) = dT/dX = u_tp = constant for a given ray
!  - Vertical slowness = eta = sqrt(u**2-p**2)
!  - Theta = ray incidence from vertical (0 = straight down, pi/2=90deg = horizontal)
!
! For each layer, integrate ray-tracing equations:
!   dT = integral(u^2/eta,dz)
!   dX = integral(p/eta,dz)
!   dTau = integral(eta,dz)
!     where Tau(p) = T(p) - p*X(p)
!
! Input:    p     =  horizontal slowness
!           h     =  layer thickness
!           utop  =  slowness at top of layer
!           ubot  =  slowness at bottom of layer
!           imth  =  interpolation method
!                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)     fastest to compute
!                         = 2,  v(z) = a - b*z               linear gradient
!                         = 3,  v(z) = a*exp(-b*z)           preferred when Earth Flattening is applied
!
! Returns:  dx    =  range offset
!           dt    =  travel time
!           irtr  =  return code
!                 = -1, zero thickness layer
!                 =  0,  ray turned above layer
!                 =  1,  ray passed through layer
!                 =  2,  ray turned within layer, 1 segment counted
!
! Note:  This version does calculation in double precision,
!        but all i/o is still single precision
!
      subroutine LAYERTRACE(p1,h1,utop1,ubot1,imth,dx1,dt1,irtr)
      implicit real*8 (a-h,o-z)
      real*4 p1,h1,utop1,ubot1,dx1,dt1
      
      ! double precision
      p=dble(p1)
      h=dble(h1)
      utop=dble(utop1)
      ubot=dble(ubot1)

      !check for zero thickness layer
      if (h.eq.0.) then                  
         dx1=0.
         dt1=0.
         irtr=-1 ! return code: zero thickness
         return         
      end if

      ! slowness of top layer
      u=utop
      
      !check for complex vertical slowness: ray turned above layer
      if (u-p.le.0.) then
         dx1=0.
         dt1=0.
         irtr=0 ! return code: ray turned above layer
         return
      end if

      ! eta = vertical slowness: sqrt(u^2-p^2)
      eta2=(u-p)*(u+p)
      eta=dsqrt(eta2)

      ! special function needed for integral at top of layer
      if (imth.eq.2) then
         y=u+eta
         if (p.ne.0.) y=y/p
         qr=dlog(y)
      else if (imth.eq.3) then ! flat-earth
         qr=atan2(eta,p)  ! atan2 has a (y,x) formulation in fortran
      end if      


      ! b factor (ray tracing integral constant in the denominator)
      if (imth.eq.1) then
          b=-(utop**2-ubot**2)/(2.*h)
      else if (imth.eq.2) then
          vtop=1./utop
          vbot=1./ubot
          b=-(vtop-vbot)/h
      else                     
          b=-dlog(ubot/utop)/h   ! flat earth
      end if  
!
      !constant velocity layer: no need to integrate
      if (b.eq.0.) then                         
         b=1./h    ! b = 1/dz
         etau=eta  ! "integrand" for Tau equation
         ex=p/eta  ! "integrand" for X equation
         irtr=1    ! return code: ray passed through layer
         go to 160
      end if

    ! ray tracing integral at upper limit, 1/b factor omitted until end
      if (imth.eq.1) then
         etau=-eta2*eta/3.
         ex=-eta*p
      else if (imth.eq.2) then
         ex=eta/u                       !*** - in some versions (wrongly)
         etau=qr-ex
         if (p.ne.0.) ex=ex/p
      else
         etau=eta-p*qr                  ! flat-earth
         ex=qr
      end if

     ! check lower limit to see if we have turning point
      u=ubot
      if (u.le.p) then  ! if turning point, there is no contribution from bottom point
         irtr=2         ! return code: ray turned within layer
         go to 160
      end if
      
      ! no turning point: ray passes through 
      irtr=1        ! return code: ray passed through layer
      eta2=(u-p)*(u+p) ! update vertical slowness for ubot
      eta=dsqrt(eta2)   ! update vertical slowness for ubot
!
      if (imth.eq.1) then
         etau=etau+eta2*eta/3.
         ex=ex+eta*p
      else if (imth.eq.2) then
         y=u+eta
         z=eta/u
         etau=etau+z
         if (p.ne.0.) then
            y=y/p
            z=z/p
         end if
         qr=dlog(y)
         etau=etau-qr
         ex=ex-z
      else
         qr=atan2(eta,p) ! atan2 has a (y,x) formulation in fortran
         etau=etau-eta+p*qr
         ex=ex-qr
      end if      


      ! ray tracing equations to get dt, dx
160   dx=ex/b         ! horizontal offset in km
      dtau=etau/b     ! delay time in seconds
      dt=dtau+p*dx    ! convert delay time to travel time

      ! back to single precision at the end
      dx1=sngl(dx)
      dt1=sngl(dt)
      return
      end

!-----------------------------------------------------------------------
