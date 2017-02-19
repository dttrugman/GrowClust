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
        
! interpolate the input model
     kk=0
     nz = floor(d(npts)/dz)+1 ! number of interpolation pts
     do iz=1,nz
         z = (iz-1)*dz
         do i=1,npts-1
            z1=d(i)
            z2=d(i+1)
            if (z1.le.z.and.z2.ge.z) then
               if (z1.eq.z2) then
                  fact=0.
               else
                  fact=(z-z1)/(z2-z1)
               end if
               do j=1,2
                  v0(j)=v(i,j)+fact*(v(i+1,j)-v(i,j))
               enddo
               kk=kk+1
               buf(kk,1)=z
               buf(kk,2)=v0(1)
               buf(kk,3)=v0(2)               
!               write (12,33) z,(v0(j),j=1,2)
33             format (f10.3,2f8.4)
               do k=i+1,npts
                  if (d(k)-z.lt.dz) then
                     kk=kk+1
                     buf(kk,1)=d(k)
                     buf(kk,2)=v(k,1)
                     buf(kk,3)=v(k,2)
!                     write (12,33) d(k),(v(k,j),j=1,2)
                  end if
               enddo
               exit ! exit inner loop
            end if
         enddo
     enddo

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
!               ttfile   = name of output travel time table file
!
!-----------------------------------------------------------------------
      subroutine deptable(vmodel,iw,plongcut,dep1,dep2,dep3,del1,del2,del3,ttfile)
      
      implicit none

! parameters
! npts0 = maximum number of lines in input velocity model
! nz0 = maximum number of depths in desired output tables
! nx0 = maximum number of distances is desired output tables
! nray0 = maximum number of rays during ray tracing
! ncount0 = 2*nray0
      integer npts0,nz0,nx0,nray0,ncount0
      parameter (npts0=1000,nz0=100,nx0=500,nray0=40002,ncount0=80002)

      real erad,pi
      parameter (erad=6371., pi=3.1415927)

      integer npts,ndep,nump,np,ncount,ndel
      integer ideptype,idep,itype,icount,idel,iw,imth,irtr
      integer i,j,i2,ideprad

      real ecircum,kmdeg,degrad,angle
      real p,pmin,pmax,pstep,plongcut,zmax,frac,xcore,tcore,h
      real dep,dep1,dep2,dep3,del,del1,del2,del3,deldel
      real scr1,scr2,scr3,scr4,xold,x,x1,x2,dx,t,t1,t2,dt
      real tbest,pbest,ubest
      real xsave(ncount0),tsave(ncount0),psave(ncount0),usave(ncount0)
      real deptab(nz0),ptab(nray0),delttab(nray0)
      real z(npts0),alpha(npts0),beta(npts0)
      real z_s(npts0),r_s(npts0),alpha_s(npts0),beta_s(npts0)
      real slow(npts0,2),deltab(nray0),tttab(nray0)
      real angang(nx0,nz0),tt(nx0,nz0),rayray(nx0,nz0),etaeta(nx0,nz0)
      real depxcor(nray0,nz0),deptcor(nray0,nz0)
      real depucor(nray0,nz0)

      character*100 vmodel
      character*100 ttfile

!-----------------------------------------------------------------------
!     define fixed variables
      ecircum=2.*pi*erad
      kmdeg=ecircum/360.
      degrad=180./pi
      
      ideprad = 1! velocity model fmt, first column of input:  1=depth, 2=radius
      zmax = 9999 ! maximum depth
      nump = 5000 ! number of rays
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
         call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i))
         call FLATTEN(z_s(i),beta_s(i),z(i),beta(i))
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
         dep2=dep2+dep3/20.
         idep=0
         ndep = floor((dep2-dep1)/dep3)+1
         do idep = 1,ndep
            dep = dep1 + (idep-1)*dep3
            deptab(idep)=dep
         enddo   
         !ndep=idep


! get number of rays to compute     
      pmin=0.
      pmax=slow(1,iw)
      print *,'pmin, pmax = ', pmin, pmax
      print *,'Number of rays to compute:'
      print *, nump   
      pstep=(pmax-pmin)/float(nump)
       

!-----------------------------------------------------------------------      
      

! --------------------------------- ray tracing -------------------------
      np=0
200   do np = 1, nump ! ------- loop over ray parameters (p) ---------------
      p = pmin + (np-1)*pstep

         ptab(np)=p

         x=0.
         t=0.
         xcore=0.
         tcore=0.
         imth=3              !preferred v(z) interpolation method
         do idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(np,idep)=0.
               deptcor(np,idep)=0.
               depucor(np,idep)=slow(1,iw)
            else
               depxcor(np,idep)=-999.
               deptcor(np,idep)=-999.
               depucor(np,idep)=-999.
            end if
        enddo

         do i=1,npts-1 ! ------ loop over layers (i) ----------------------------
         
             !check to see if z exceeds zmax
             if (z_s(i).ge.zmax) then                          
                deltab(np)=-999.
                tttab(np)=-999.
                go to 200
             end if

    ! LAYERTRACE calculates the travel time and range offset for ray tracing through a single layer.
    !  Input:   p     =  horizontal slowness
    !           h     =  layer thickness
    !           utop  =  slowness at top of layer
    !           ubot  =  slowness at bottom of layer
    !           imth  =  interpolation method
    !                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)     fastest to compute
    !                         = 2,  v(z) = a - b*z               linear gradient
    !                         = 3,  v(z) = a*exp(-b*z)           preferred when Earth Flattening is applied
    !  Returns: dx    =  range offset
    !           dt    =  travel time
    !           irtr  =  return code (-1: zero thickness layer, 0: ray turned above layer, 
    !                 =    1: ray passed through layer, 2: ray turned within layer, 1 segment counted)
             h=z(i+1)-z(i)
             if (h.eq.0.) cycle                       !skip if interface
             call LAYERTRACE(p,h,slow(i,iw),slow(i+1,iw),imth,dx,dt,irtr)
             x=x+dx
             t=t+dt

             if (irtr.eq.0.or.irtr.eq.2) exit           !ray has turned
         
             do idep=1,ndep
                if (abs(z_s(i+1)-deptab(idep)).lt.0.1) then
                   depxcor(np,idep)=x
                   deptcor(np,idep)=t
                   depucor(np,idep)=slow(i+1,iw)            
                end if
             enddo

!   
         enddo !------------ end loop on layers----------------------------------
         
         x=2.*x
         t=2.*t

         deltab(np)=x                   !stored in km
         tttab(np)=t                    !stored in seconds

      enddo             !---------------- end loop on ray parameter p --------------------
!      print *,'Completed ray tracing loop'

!----------------------------------------------
      
! special section to get (0,0) point
         np=np+1
         ptab(np)=slow(1,iw)
         deltab(np)=0.
         tttab(np)=0.
         do idep=1,ndep
            if (deptab(idep).eq.0.) then
               depxcor(np,idep)=0.
               deptcor(np,idep)=0.
               depucor(np,idep)=slow(1,iw)
            else
               depxcor(np,idep)=-999.
               deptcor(np,idep)=-999.
               depucor(np,idep)=-999.
            end if         
         enddo

!-------------------------------------------------------------------------

!------------- Format for output ----------


! now compute T(X,Z) for first arriving rays
      do idep=1,ndep
         icount=0
         xold=-999.
         if (deptab(idep).eq.0.) then
            i2=np
            go to 223
         end if
         do i=1,np                         !upgoing rays from source
            x2=depxcor(i,idep)
            if (x2.eq.-999.) exit
            if (x2.le.xold) exit            !stop when heads inward
            t2=deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=-ptab(i)             !sav as negative for upgoing from source
            usave(icount)=depucor(i,idep)
            xold=x2
         enddo
         i2=i-1
         
223      continue         
         do i=i2,1,-1                    !downgoing rays from source
            if (depxcor(i,idep).eq.-999.) cycle
            if (deltab(i).eq.-999.) cycle
            x2=deltab(i)-depxcor(i,idep)
            t2=tttab(i)-deptcor(i,idep)
            icount=icount+1
            xsave(icount)=x2
            tsave(icount)=t2
            psave(icount)=ptab(i)
            usave(icount)=depucor(i,idep)
            xold=x2
         enddo   
         ncount=icount
         
         ndel = floor((del2-del1)/del3) + 1 ! number of interpolation pts
         do idel = 1, ndel
            deldel = del1 + (idel-1)*del3
            
            del=deldel
            if (itype.eq.2) del=deldel*kmdeg
            delttab(idel)=deldel
            
            tt(idel,idep)=99999.
            do i=2,ncount
               x1=xsave(i-1)
               x2=xsave(i)
               if (x1.gt.del.or.x2.lt.del) cycle
               if (psave(i).gt.0..and.psave(i).lt.plongcut) cycle
               frac=(del-x1)/(x2-x1)
               tbest=tsave(i-1)+frac*(tsave(i)-tsave(i-1))
               if (psave(i-1).le.0..and.psave(i).le.0. .or. &
                         psave(i-1).ge.0..and.psave(i).ge.0.) then
                  pbest=psave(i-1)+frac*(psave(i)-psave(i-1))
                  ubest=usave(i-1)+frac*(usave(i)-usave(i-1)) 
               else
                  if (frac.lt.0.5) then
                     pbest=psave(i-1)
                     ubest=usave(i-1)
                  else
                     pbest=psave(i)
                     ubest=usave(i)
                  end if
               end if
              
               if (tbest.lt.tt(idel,idep)) then
                  tt(idel,idep)=tbest
                  scr1=pbest/ubest
                  if (scr1.gt.1.) then
!                     print *,'***Warning: p>u in angle calculation'
!                     print *,'   Ray assumed horizontal'
!                     print *,deptab(idep),del,tbest,pbest,ubest
                     scr1=1.
                  end if
                  angle=asin(scr1)*degrad
                  if (angle.lt.0.) then
                     angle=-angle
                  else
                     angle=180.-angle
                  end if
                  angang(idel,idep)=angle
                  rayray(idel,idep)=pbest
                  etaeta(idel,idep)=ubest*sqrt(1.-scr1**2)                  
                  if (angang(idel,idep).lt.90.) then
                     etaeta(idel,idep)=-etaeta(idel,idep)
                  endif
               end if
            enddo
            if (tt(idel,idep).eq.99999.) tt(idel,idep)=0.
            if (itype.eq.2) tt(idel,idep)=tt(idel,idep)/60.            

         enddo                                    !end loop on range
         !ndel=idel
                  
      enddo                                        !end loop on depth

      if (delttab(1).eq.0.) then
         if (deptab(1).eq.0.) tt(1,1)=0.                 !set tt to zero at (0,0)
         do idep=1,ndep
            angang(1,idep)=0.                 !straight up at zero range
            etaeta(1,idep)=-abs(etaeta(1,idep))
         enddo

      end if

!-------------- Make output files -----------------

! get file name from input
      print *,'Output file name for travel-time table:'
      print *, ttfile
      open (11,file=ttfile)     
      
! write headers for all files
     
! first line: model name, iw (1/2 for P/S), pmin, pmax, np   
      write (11,407) vmodel(1:20),iw,pmin,pmax,np  
407    format ('From deptable, file= ',a20,' iw =',i2,' pmin=',f8.5, &
              ' pmax=',f8.5,' np=',i6)
     
! second line: table size ndel,ndep:
!       ndel rows (different X/DEL offsets)
!       ndep columns (different source depths)
      write (11,408) ndel,ndep
408      format (2i5)
      
!     third line: row of source depths
      write (11,409) (deptab(j),j=1,ndep)
409      format (8x,100f8.1)
      
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

410      format (101f8.4)
413      format (f8.3,100f8.4)

      enddo
      close (11)

!999   stop
999   return

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
      p=dble(p1)
      h=dble(h1)
      utop=dble(utop1)
      ubot=dble(ubot1)
!
      if (h.eq.0.) then                  !check for zero thickness layer
         dx1=0.
         dt1=0.
         irtr=-1
         return         
      end if
!
      u=utop
      y=u-p
      if (y.le.0.) then                       !complex vertical slowness
         dx1=0.
         dt1=0.
         irtr=0
         return
      end if
!
      q=y*(u+p)
      qs=dsqrt(q)
!
! special function needed for integral at top of layer
      if (imth.eq.2) then
         y=u+qs
         if (p.ne.0.) y=y/p
         qr=dlog(y)
      else if (imth.eq.3) then
         qr=atan2(qs,p)
      end if      
!
      if (imth.eq.1) then
          b=-(utop**2-ubot**2)/(2.*h)
      else if (imth.eq.2) then
          vtop=1./utop
          vbot=1./ubot
          b=-(vtop-vbot)/h
      else
          b=-dlog(ubot/utop)/h
      end if  
!
      if (b.eq.0.) then                         !constant velocity layer
         b=1./h
         etau=qs
         ex=p/qs
         irtr=1
         go to 160
      end if
!
! integral at upper limit, 1/b factor omitted until end
      if (imth.eq.1) then
         etau=-q*qs/3.
         ex=-qs*p
      else if (imth.eq.2) then
         ex=qs/u                       !*** - in some versions (wrongly)
         etau=qr-ex
         if (p.ne.0.) ex=ex/p
      else
         etau=qs-p*qr
         ex=qr
      end if
!
! check lower limit to see if we have turning point
      u=ubot
      if (u.le.p) then                                !if turning point,
         irtr=2                                    !then no contribution
         go to 160                                    !from bottom point
      end if 
      irtr=1
      q=(u-p)*(u+p)
      qs=dsqrt(q)
!
      if (imth.eq.1) then
         etau=etau+q*qs/3.
         ex=ex+qs*p
      else if (imth.eq.2) then
         y=u+qs
         z=qs/u
         etau=etau+z
         if (p.ne.0.) then
            y=y/p
            z=z/p
         end if
         qr=dlog(y)
         etau=etau-qr
         ex=ex-z
      else
         qr=atan2(qs,p)
         etau=etau-qs+p*qr
         ex=ex-qr
      end if      
!
160   dx=ex/b
      dtau=etau/b
      dt=dtau+p*dx                                     !convert tau to t
!
      dx1=sngl(dx)
      dt1=sngl(dt)
      return
      end

!-----------------------------------------------------------------------