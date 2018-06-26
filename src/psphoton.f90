! Program to model seismic phases as energy packets
! This version handles P and S of arbitrary polarization (SH and SV)
!
!  Originally written 2003-2004 by Peter Shearer (pshearer@ucsd.edu).
!    This version managed by Nick Mancinelli (n.j.mancinelli@gmail.com) 
!    to whom questions should be addressed.
!
!    Shearer, P.M and P.S. Earle, The global short-period wavefield
!      modelled with a Monte Carlo seismic phonon method, Geophys.
!      J. Int., 158, 1103-1117, 2004.
!
!    Converted to Mac G77 Fortran from Sun Fortran in 2006
!    Implicit none, variable initialization and debug options added
!
!    Converted to f90 in 2012 by Nick Mancinelli (n.j.mancinelli@gmail.com)
!
!    Updated 2012-2017 by Nick Mancinelli
!      -Max array sizes increased (gives smoother results for PKP precursors)
!      -EXPSATO generalized to use Von Karman ACF (kappa currently hardwired to 0.5)
!      -Velocity ratio now calculated for each scattering layer (previously approximated as sqrt(3))
!      -Email Nick for versions of this code that model interface topography.
!    
!    This file includes all the necessary subroutines; no linking
!    to external routines is required.  Compiling with optimization
!    is recommended for maximum speed.
!
!    The following files are automatically generated:
!
!    out.photon                  all energy
!    out.photon_z                vertical component energy
!    out.photon_zcore            vertical for p.le.5 s/deg.
!    out.photon_rad              radial component energy
!    out.photon_p                P wave energy
!    out.photon_sv               SV wave energy
!    out.photon_sh               SH wave energy
!
!    Note:  The program will run forever and must be manually killed.
!           It writes and overwrites the output files at regular
!           intervals.  Try to kill the program between these writes.
!
      program psphoton
	  
	  implicit none
	
      integer, parameter :: nlay0=655,nray0=50000
	  
!      parameter (nlay0=1000,nray0=15000)      !nray0 also in UPDATE_P and SCATRAYPOL
      
	  real :: z_s(nlay0),alpha_s(nlay0),beta_s(nlay0),rho(nlay0)
      real :: z(nlay0),alpha(nlay0),beta(nlay0),q(nlay0,2),slow(nlay0,2)
      real :: p(nray0,2),dx(nlay0,nray0,2),dt(nlay0,nray0,2)
      real :: dtstar(nlay0,nray0,2)
      
      real    :: amp0(nray0,2),depsave(1000,2)
      integer :: iwsave(1000)
!      real amp0(nray0,2),depsave(10000,2)
!      integer iwsave(10000)
      
      integer, parameter :: nface0=6          !also is maximum number of scattering zones, q zones
!      parameter (nface0=10)

      real :: depface(nface0),vp(nface0,2),vs(nface0,2),den(nface0,2)
      real :: scatprob(nface0,2,0:2)
      real :: zminscat(nface0),zmaxscat(nface0),xmaxscat(nface0)
      real :: el(nface0),nu(nface0),gam0(nface0),eps(nface0),alen(nface0) 
      real :: zminq(nface0),zmaxq(nface0),qalpha(nface0)
      real :: rt(2,2,3,3,nface0,nray0,2)
      
      integer             :: iddir(nlay0,nray0,2),iflag(nlay0),iflagvol(nlay0)
      character (len=100) :: vmodel

      integer, parameter  :: ntdim=3001,nxdim=360   
      real                :: rbin(ntdim,nxdim),dummy(10)
      real                :: rbin_p(ntdim,nxdim),rbin_sv(ntdim,nxdim),rbin_sh(ntdim,nxdim)
      real                :: rbin_z(ntdim,nxdim),rbin_rad(ntdim,nxdim)
      real                :: rbin_zcore(ntdim,nxdim)
      character (len=100) :: outfile,debugfile
!
      integer, parameter  :: nslowdim=49,nslowout=5
      real                :: slowstack(nslowdim,nslowdim,nslowout)
      real                :: xslow1(nslowout),xslow2(nslowout)
      real                :: tslow1(nslowout),tslow2(nslowout)
      real                :: psecdeg,azidum1,azidum2,slowang,xdum,psecdegtran,psecdegrad      

      real :: kmdeg,erad,pi,ecircum,degrad,freq,om,zsource
      real :: spenergy_ratio,pmax_takeoff,xwind1,xwind2,twind1,twind2
      real :: pvelref,svelref,dep,vpmin,vsmin,vpsource,vssource
      real :: gpp0,gps0,gsp0,gss0,gppm,gpsm,gspm,gssm,expdev,spol
      real :: pmin,vel0,pmax,dp,p1,p2,vel0minus,dang0,dang,angcor,sinthe
      real :: h,r,flatfact,pcor,tmax,ampstart,svsh,ran1,fran,x,t,plat,plon
      real :: azi,tstar,xdeg,plat2,plon2,xkm,afp,freepath,distsum,dist,del
      real :: azidum
      real :: psi2,zeta2,volang,azi2,vp0,vs0,rt1,rt2,rt3,rt4,test,tmin
      real :: test1,test2,test3,svfrac,shfrac,svamp,shamp,amp,amp2
      real :: amp_vert,amp_rad,energy_rad,energy_vert,rsum,rmax
      real :: stnmin,stnmax,pcoremax
                           
      integer :: ix,it,iwstart,np,iforcerefl,nface,nscatvol,i,j,k,irtr,idum
      integer :: nscatmin,nscatmax,nqlay,nraydump,izsource,npts,iscat
      integer :: iscatold,iface,iwave,iw,imth,ilay,ip,nlay,nray,nsurf,idir
      integer :: iww1,iww2,idirbeg,idirend,nsweep,iw_init,idir_init,ipp
      integer :: ixdir,isave,nscat,nscatold,iwave1,iwoff,iwave2,idir2
      integer :: nsave,iw2,ip2,n,idebug,idebug0,iwrap
      
      do ilay=1,nlay0
         z_s(ilay)=0.
         alpha_s(ilay)=0.
         beta_s(ilay)=0.
         rho(ilay)=0.
         z(ilay)=0.
         alpha(ilay)=0.
         beta(ilay)=0.
         iflag(ilay)=0
         iflagvol(ilay)=0
         do i=1,2
            q(ilay,i)=0.
            slow(ilay,i)=0.
            do nray=1,nray0
               p(nray,i)=0.
               dx(ilay,nray,i)=0.
               dt(ilay,nray,i)=0.
               dtstar(ilay,nray,i)=0.
               amp0(nray,i)=0.
               iddir(ilay,nray,i)=0
            enddo
         enddo
      enddo
      
      do i=1,1000
         depsave(i,1)=0.
         depsave(i,2)=0.
         iwsave(i)=0
      enddo
      
      do iface=1,nface0
         depface(iface)=0.
         do i=1,2
            vp(iface,i)=0.
            vs(iface,i)=0.
            den(iface,i)=0.
            do k=0,2
               scatprob(iface,i,k)=0.
            enddo
         enddo
         zminscat(iface)=0.
         zmaxscat(iface)=0.
         xmaxscat(iface)=0.
         el(iface)=0.
         nu(iface)=0.
         gam0(iface)=0.
         eps(iface)=0.
         alen(iface)=0.
         zminq(iface)=0.
         zmaxq(iface)=0.
         qalpha(iface)=0.
         do ix=1,2
         do it=1,2
         do iw=1,3
         do ip=1,3
         do nray=1,nray0
         do i=1,2         
            rt(ix,it,iw,ip,iface,nray,i)=0.
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      enddo
      
      do i=1,10
         dummy(i)=0.
      enddo
      do it=1,ntdim
         do ix=1,nxdim
            rbin(it,ix)=0.
            rbin_p(it,ix)=0.
            rbin_sv(it,ix)=0.
            rbin_sh(it,ix)=0.
            rbin_z(it,ix)=0.
            rbin_rad(it,ix)=0.
            rbin_zcore(it,ix)=0.
         enddo
      enddo
      do i=1,nslowdim
         do j=1,nslowdim
            do k=1,nslowout
               slowstack(i,j,k)=0.
            enddo
         enddo
      enddo
      vmodel=' '
      outfile=' '
      debugfile=' '
         
                              
      kmdeg=0.
      erad=0.
      pi=0.
      ecircum=0.
      degrad=0.
      freq=0.
      om=0.
      zsource=0.
      spenergy_ratio=0.
      pmax_takeoff=0.
      xwind1=0.
      xwind2=0.
      twind1=0.
      twind2=0.
      pvelref=0.
      svelref=0.
      dep=0.
      vpmin=0.
      vsmin=0.
      vpsource=0.
      vssource=0.
      gpp0=0.
      gps0=0.
      gsp0=0.
      gss0=0.
      gppm=0.
      gpsm=0.
      gspm=0.
      gssm=0.
      spol=0.
      pmin=0.
      vel0=0.
      pmax=0.
      dp=0.
      p1=0.
      p2=0.
      vel0minus=0.
      dang0=0.
      dang=0.
      angcor=0.
      sinthe=0.
      h=0.
      r=0.
      flatfact=0.
      pcor=0.
      tmax=0.
      ampstart=0.
      svsh=0.
      fran=0.
      x=0.
      t=0.
      plat=0.
      plon=0.
      azi=0.
      tstar=0.
      xdeg=0.
      plat2=0.
      plon2=0.
      xkm=0.
      afp=0.
      freepath=0.
      distsum=0.
      dist=0.
      del=0.
      azidum=0.
      psi2=0.
      zeta2=0.
      volang=0.
      azi2=0.
      vp0=0.
      vs0=0.
      rt1=0.
      rt2=0.
      rt3=0.
      rt4=0.
      test=0.
      tmin=0.
      test1=0.
      test2=0.
      test3=0.
      svfrac=0.
      shfrac=0.
      svamp=0.
      shamp=0.
      amp=0.
      amp2=0.
      amp_vert=0.
      amp_rad=0.
      energy_rad=0.
      energy_vert=0.
      rsum=0.
      rmax=0.
      stnmin=0.
      stnmax=0.     
                  
      ix=0
      it=0
      iwstart=0
      np=0
      iforcerefl=0
      nface=0
      nscatvol=0
      i=0
      k=0
      irtr=0
      idum=0
      nscatmin=0
      nscatmax=0
      nqlay=0
      nraydump=0
      izsource=0
      npts=0
      iscat=0
      iscatold=0
      iface=0
      iwave=0
      iw=0
      imth=0
      ilay=0
      ip=0
      nlay=0
      nray=0
      nsurf=0
      idir=0
      iww1=0
      iww2=0
      idirbeg=0
      idirend=0
      nsweep=0
      iw_init=0
      idir_init=0
      ipp=0
      ixdir=0
      isave=0
      nscat=0
      nscatold=0
      iwave1=0
      iwoff=0
      iwave2=0
      idir2=0
      nsave=0
      iw2=0
      ip2=0
      n=0
      idebug=0
      idebug0=0                    

      idum=0
      erad=6371.
      pi=3.1415927
      ecircum=2.*pi*erad
      kmdeg=ecircum/360.
      degrad=180./pi
      stnmin=0.               !used only for output files
      stnmax=0.
      pcoremax=5.             !5 s/deg is max slowness for rbin_zcore

      do ix=1,nxdim            !some compilers don't initialize to zero
      do it=1,ntdim
         rbin(it,ix)=0.
         rbin_p(it,ix)=0.
         rbin_sv(it,ix)=0.
         rbin_sh(it,ix)=0.
         rbin_z(it,ix)=0.
         rbin_rad(it,ix)=0.
      enddo
      enddo

      print *,'Enter velocity model file name'
      read (*,'(a)') vmodel

      print *,'Enter frequency (Hz)'
      read *,freq
      om=2.*pi*freq

      print *,'Enter source depth (km)'
      read *,zsource

      print *,'Radiate:  (1) P,  (2) SH, (3) SV, (4) SH/SV, (5) P+S'
      print *,'          or (6) P+S with variable ratio'
      print *,'          or (7) P with restricted p takeoff range'
      read *,iwstart
      if (iwstart.eq.6) then
         print *,'Enter S/P energy ratio (e.g., 10)'
         read *,spenergy_ratio
      else if (iwstart.eq.7) then
         print *,'Enter maximum p for P radiation'
         read *,pmax_takeoff
      end if

      print *,'Enter number of ray parameters (max=50000)'
      read *,np

      print *,'Enter coordinates of box to dump ray info'
      print *,'(Useful for identifying mystery phases)'
      print *,'Enter zeros if no dump wanted'
      print *,'Enter xwind1,xwind2,twind1,twind2 (deg and minutes)'
      read *,xwind1,xwind2,twind1,twind2

      print *,'Enter:  (0) normal  or  (1) force equal refl/trans'
      print *,'(  enter (-1) to input slowness stack numbers)'
      read *,iforcerefl
      
      if (iforcerefl.eq.-1) then
         do k=1,nslowout
            print *,'Enter xwind1,xwind2,twind1,twind2 (deg and sec)'
            read *,xslow1(k),xslow2(k),tslow1(k),tslow2(k)
         enddo
      endif

      print *,'Enter number of interfaces for refl/trans'
      read *,nface
      do i=1,nface
         print *,'Enter depth for interface ',i
         read *,depface(i)
      enddo

      print *,'Enter number of scattering volumes'
      read *,nscatvol
      do i =1,nscatvol
         print *,'ENTER PARAMETERS FOR SCATTERING LAYER: ',i
         print *,'Enter min scattering depth (km)'
         read *,zminscat(i)
         print *,'Enter max scattering depth (km)'
         read *,zmaxscat(i)
         print *,'Enter max range from source (km)'
         read *,xmaxscat(i)
! If outside xmaxscat, then iscat=iscat-1 will revert to previous volume
! Thus overlapping depth ranges can be used

         print *,'Enter reference P & S velocities for layer'
         print *,'These are use to compute Poisson ratio'
         print *,'  and S wavenumber)'
         read *, pvelref,svelref
         el(i)=om/svelref
         gam0(i)=pvelref/svelref
	 print *,'Velocity ratio = ', gam0(i)
         print *,'Enter relative size of density pert (e.g., 0.8)'
         read *,nu(i)
         print *,'Enter rms perturbation (e.g., 0.01 for 1%)'
         read *,eps(i)
         print *,'Enter scale length (km)'
         read *,alen(i)
      enddo

      print *,'Enter min,max number of scattering events for output'
      read *,nscatmin,nscatmax

      print *,'Enter number of Q layers'
      read *,nqlay
      do i=1,nqlay
         print *,'Enter zmin,max for Q layer: ',i
         read *,zminq(i),zmaxq(i)
         print *,'Enter Qalpha for layer: ',i
         read *,qalpha(i)
      enddo

      print *,'Write output at ray multiples of: (e.g., 20000)'
      read *,nraydump
                  
      print *,'debug file: (0) none, (1) full dump, (2) scat events'
      read *,idebug
      if (idebug.ne.0) then
         debugfile='out.debug'      
         open (18,file=debugfile)
         open (19,file='out.surface')
      endif
      idebug0=idebug
                     
! read velocity model and transform to flat earth
! value at center of earth is removed to
! avoid singularity in transformation
      open (11,file=vmodel,status='old')
      do 10 i=1,nlay0
         read (11,*,end=12) z_s(i),dep,alpha_s(i),beta_s(i),rho(i)
         if (z_s(i).eq.zsource) izsource=i
         if (z_s(i).eq.erad) go to 12
         if (beta_s(i).eq.0.) beta_s(i)=0.001       !****KLUGE for r/t calc.
         call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i))
         if (beta_s(i).ne.0.) then
            call FLATTEN(z_s(i),beta_s(i),z(i),beta(i))
         else
            beta(i)=0.
         end if
         iflag(i)=0
         iflagvol(i)=0

         q(i,1)=999999.

         do k=1,nqlay
            if (z_s(i).ge.zminq(k).and.z_s(i).le.zmaxq(k)) then
               q(i,1)=qalpha(k)
            end if
         enddo

         q(i,2)=(4./9.)*q(i,1)    !set Q for S using approx. for Poisson solid

10    continue
      print *,'***nlay0 maximum exceeded in model'
      stop
12    close (11)
      vpmin=alpha(1)
      vsmin=beta(1)
      vpsource=alpha(izsource)
      vssource=beta(izsource)
      print *,'finished reading model'
      print *,'vpmin,vsmin = ',vpmin,vsmin
      print *,'izsource = ',izsource
      print *,'vpsource,vssource = ',vpsource,vssource
      if (izsource.eq.0) then
         print *,'***Problem.  izsource = 0.  zsource = ',zsource
         stop
      end if
!
      z_s(i)=z_s(i-1)              !set up dummy interface at bottom
      alpha_s(i)=alpha_s(i-1)
      beta_s(i)=beta_s(i-1)
      q(i,1)=q(i-1,1)
      q(i,2)=q(i-1,2)

      call FLATTEN(z_s(i),alpha_s(i),z(i),alpha(i))
      call FLATTEN(z_s(i),beta_s(i),z(i),beta(i))
      npts=i
      print *,'Depth points in model= ',npts

      do 20 i=1,npts
         slow(i,1)=1./alpha(i)
         if (beta(i).ne.0.) then
            slow(i,2)=1./beta(i)
         else
            slow(i,2)=1./alpha(i)    !fluid legs are always P!
         end if    
20    continue

! set up scattering volume flags and other stuff
! Note that if scattering volumes overlap, iflagvol is set to
! last of scattering volumes
      do 30 iscat=1,nscatvol

!          call GAVGSATO(el(iscat),nu(iscat),gam0(iscat),eps(iscat),
!     &      alen(iscat), gpp0,gps0,gsp0,gss0)

          call GAVGSATO2(el(iscat),nu(iscat),gam0(iscat),eps(iscat), &
           alen(iscat), gpp0,gps0,gsp0,gss0, &
           gppm,gpsm,gspm,gssm)

            print *,'SCATTERING STRENGTHS for LAYER: ',iscat

            print *,'  gpp0,gps0,gsp0,gss0= ',gpp0,gps0,gsp0,gss0
            print *,'  afp = ', 1./(gpp0+gps0), 1./(gsp0+gss0)
            print *,'  gppm,gpsm,gspm,gssm= ',gppm,gpsm,gspm,gssm
            print *,'  mom. afp = ', 1./(gppm+gpsm), 1./(gspm+gssm)

            scatprob(iscat,1,0)=gpp0+gps0
            scatprob(iscat,2,0)=gsp0+gss0
            scatprob(iscat,1,1)=gpp0
            scatprob(iscat,1,2)=gps0
            scatprob(iscat,2,1)=gsp0
            scatprob(iscat,2,2)=gss0

         do 25 i=2,npts
            if (z_s(i-1).ge.zminscat(iscat).and. &
                  z_s(i).le.zmaxscat(iscat)) then
               iflagvol(i)=iscat
            end if
25       continue
30    continue            

! ---------------------------------- now define interface values
      i=1
      depface(i)=0.
      vp(i,1) = 0.1            !**kluge to get surface reflections
      vs(i,1) = 0.01
      den(i,1)= 0.01
      vp(i,2) = alpha_s(1)
      vs(i,2) = beta_s(1)
      den(i,2)= rho(1)
      iflag(1)=1

      do 45 iface=2,nface

      do 44 i=2,npts
!         print *,depface(iface),z_s(i),z_s(i-1)
         if (z_s(i).eq.z_s(i-1).and.z_s(i).eq.depface(iface)) then
            vp(iface,1) = alpha_s(i-1)
            vs(iface,1) = beta_s(i-1)
            den(iface,1)= rho(i-1)
            vp(iface,2) = alpha_s(i)
            vs(iface,2) = beta_s(i)
            den(iface,2)= rho(i)
!            print *,depface(iface),(vp(iface,k),vs(iface,k),
!     &              den(iface,k),k=1,2)
            iflag(i-1)=iface
            go to 45
         end if
44    continue
      print *,'***Interface depth not found: ',depface(iface)
      stop

45    continue
!--------------------------------------------------------------------- 

      do 110 iwave=1,2          !=1 for P, =2 for S
         iw=iwave 

      pmin=1.25E-05    !avoid problems with flat earth near p=0
      if (iwave.eq.1) then
         vel0=vpsource
         pmax=1./(vpmin+0.0001)
      else
         vel0=vssource
         pmax=1./(vsmin+0.0001)
      end if
      dp=(pmax-pmin)/float(np-1)
      p1=pmin
      p2=pmin+dp
      vel0minus=vel0-0.01       !used in takeoff angle calculation

!
      dang0=asin(p2*vel0minus)-asin(p1*vel0minus)
      do 48 i=1,np
         amp0(i,iw)=0.
         p(i,iw)=pmin+float(i-1)*dp
         if (p(i,iw)*vel0.gt.1.) go to 48   !amp0=0 for rays turning above source
         if (i.eq.1) then
            angcor=1.
         else
            dang=asin(p(i,iw)*vel0minus)-asin(p(i-1,iw)*vel0minus)
            angcor=dang/dang0
         end if
         sinthe=p(i,iw)*vel0minus
         amp0(i,iw)=sqrt(angcor*sinthe)      !amp not energy here
!         print *,i,p(i,iw),angcor,sinthe
!         print *,iw,i,p(i,iw),angcor,sinthe,amp0(i,iw)
48    continue
!      stop


      imth=3                  !best for spherical models
      do 100 i=1,npts-1

         h=z(i+1)-z(i)
         ilay=i

         do 50 ip=1,np
            call LAYERTRACE(p(ip,iw),h,slow(i,iw),slow(i+1,iw),imth, &
                            dx(ilay,ip,iw),dt(ilay,ip,iw),irtr)
            if (irtr.eq.-1) then      !zero thickness layer
               iddir(ilay,ip,iw)=1
            else if (irtr.eq.1) then  !passed through
               iddir(ilay,ip,iw)=1
            else if (irtr.eq.0) then  !turned above
               iddir(ilay,ip,iw)=-1
            else if (irtr.eq.2) then  !turned within, 1 segment counted
               iddir(ilay,ip,iw)=-1
               dx(ilay,ip,iw)=2.*dx(ilay,ip,iw)
               dt(ilay,ip,iw)=2.*dt(ilay,ip,iw)
            else
               print *,'***ERROR: irtr = ',irtr,i,ilay,ip,p(ip,iw)
               stop
            end if

            dtstar(ilay,ip,iw)=dt(ilay,ip,iw)/q(i+1,iw)

            if (iflag(i).ne.0) then   !we have interface
               iface=iflag(i)

               r=erad-depface(iface)
               flatfact=r/erad
               pcor=p(ip,iw)/flatfact

               call RTCOEF_POW(vp(iface,1),vs(iface,1),den(iface,1),   &
                           vp(iface,2),vs(iface,2),den(iface,2),  &
                           pcor,rt(1,1,1,1,iface,ip,iw) )

               if (iforcerefl.eq.1) then
                  call FORCE_RT(rt(1,1,1,1,iface,ip,iw),iforcerefl)
               end if
               
               if (idebug.eq.3.and.iw.eq.1.and.iface.eq.2) then
                  print *,iw,ip,p(ip,iw),rt(1 ,1,1,1,iface,ip,iw), &
                       rt(2,1,1,1,iface,ip,iw)   
               endif

            end if
   
50       continue

100   continue
      nlay=ilay
      print *,'finished ray tracing, iwave,nlay,np = ',iwave,nlay,np

110   continue    !end loop on iwave for P and S

      tmax=3000.     !3000 s = 50 minutes
      nray=0
      nsurf=0
!
!--------------------------------------------------------------------
!--------------------------------------- start spraying rays
      if (iwstart.eq.1.or.iwstart.eq.7) then        !P only
         iww1=1
         iww2=1
      else if (iwstart.le.4) then   !S only
         iww1=2
         iww2=2
      else                          !P and S
         iww1=1
         iww2=2
      end if
      if (izsource.eq.1) then    !surface source, takeoff down only
         idirbeg=1
         idirend=1
      else                       !source at depth, takeoff both up & down
         idirbeg=-1
         idirend=1
      end if

      nsweep=0
150   do 800 iw_init=iww1,iww2
         do 780 idir_init=idirbeg,idirend,2 
            nsweep=nsweep+1
            do 760 ipp=1,np
               ip=ipp
               iw=iw_init
               if (iwstart.eq.7.and.p(ip,iw).gt.pmax_takeoff) go to 760     !*****
               idir=idir_init

      ampstart=amp0(ip,iw)
      if (ampstart.eq.0.) go to 760   !ray parameters not radiated by source
      svsh=0.                        !pure SV polarization
      if (iw.eq.2) then
         if (iwstart.ne.6) then
            ampstart=ampstart*sqrt(23.4)
         else
            ampstart=ampstart*sqrt(spenergy_ratio)
         end if
         svsh=0.          !pure SV polarization
         if (iwstart.eq.2) then
            svsh=pi/2.
         else if (iwstart.ge.4.and.iwstart.lt.7) then
            fran=RAN1(idum)
            svsh=fran*pi     !start with random S polarization
         end if
      end if
      nray=nray+1
!      idum=-1                   !****reset random numbers (remove this when debugged!!!!)
            
      if (nray.le.1000.or.mod(nray,500).eq.0) then      !****limit debug output
         idebug=idebug0
      else
         idebug=0
      endif      

!      print *,'150: nray,ipp,p,iw,idir = ',nray,ipp,p(ip,iw),iw,idir

      if (idebug.eq.1.or.idebug.eq.3) then
         write (18,*) 'SPRAY:',nray,ipp,p(ip,iw),iw,idir
      endif

      ixdir=1
      if (idir.eq.1) then
         ilay=izsource
      else
        ilay=izsource-1
      end if
      x=0.
      t=0.
      plat=0.    !assume source on equator
      plon=0.
      azi=90.
      tstar=0.
      isave=0
      nscat=0
      iscatold=-99

200   x=x+dx(ilay,ip,iw)*float(ixdir)       !x in km
      t=t+dt(ilay,ip,iw)                    !t in seconds
      if (t.gt.tmax) go to 350
      tstar=tstar+dtstar(ilay,ip,iw)
      idir=idir*iddir(ilay,ip,iw)

!      print *,'200 = ',iw,ip,ilay,x,t,idir
!      print *,'ampstart = ',ampstart
!      if (nray.eq.467.or.nray.eq.648) print *,'2: ',ip,iw,idir,ilay

      if (idebug.eq.1) write (18,*) nray,ip,iw,idir,ilay,t,z_s(ilay)

! first check for volume scattering
         if (iflagvol(ilay).ne.0.and.dt(ilay,ip,iw).ne.0.) then  
            iscat=iflagvol(ilay)

! the following little section is kluge to allow stronger near-source scat.
            if (xmaxscat(iscat).lt.180.) then
               xdeg=x/kmdeg
               call SPH_LOC(plat,plon,xdeg,azi,plat2,plon2)
               xkm=(90.-plat2)*kmdeg
               if (x.gt.xmaxscat(iscat)) iscat=iscat-1
               if (iscat.lt.1) then
                  print *,'**iscat sync problem: ',iscat,ilay,x
                  stop
               end if
            end if     !------- end klugey section

            if (iscat.ne.iscatold) then     !entered new scattering volume
               afp=1./scatprob(iscat,iw,0)
               freepath=EXPDEV(idum)*afp
               distsum=0.
            end if
            iscatold=iscat

            if (iw.eq.1) then
               dist=dt(ilay,ip,iw)*alpha(ilay)
            else
               dist=dt(ilay,ip,iw)*beta(ilay)
            end if
            distsum=distsum+dist

!--------------------------------------------------------------------------             
!--------------------SCATTERING SECTION-------------------------------------
            if (distsum.ge.freepath) then           !SCATTERED
                                    
               nscat=nscat+1
               iscatold=-99           !***need to reset this (this bug still in Sun version)
               
               if (idebug.ne.0) then
                 write (18,*) 'SC: ',nray,ip,iw,p(ip,iw),idir,ilay,iscat
                 write (18,*) '  ',nscat,distsum,freepath,afp,x,t
               endif
                                             
               distsum=0.

               iwave1=iw
               iwoff=(iw-1)*2     !for RANDANG2 input               

               fran=RAN1(idum)*scatprob(iscat,iw,0)
                              
               if (idebug.ne.0) then
                 write (18,*) '  ',fran,scatprob(iscat,iw,1)
               endif                               
               
               if (fran.le.scatprob(iscat,iw,1)) then
                  call RANDANG2(iscat,iwoff+1,el(iscat),nu(iscat), &
                    gam0(iscat),eps(iscat),alen(iscat),psi2,zeta2,spol)
                  iwave2=1
               else
                  call RANDANG2(iscat,iwoff+2,el(iscat),nu(iscat), &
                   gam0(iscat),eps(iscat),alen(iscat),psi2,zeta2,spol)
                  iwave2=2
               end if
            
               volang=psi2*degrad      !input angle for subroutines
!               print *,'Scat: iw,iwave2,ang = ',iw,iwave2,volang
               if (idebug.eq.1) write (18,*) 'Scat: ',iw,iwave2,volang
               
               if (idebug.ne.0) then
                  write (18,*) 'A: ',volang,psi2,zeta2,spol,iwave2
               endif                              

               xdeg=x/kmdeg
               call SPH_LOC(plat,plon,xdeg,azi,plat2,plon2)

! *** the following fixes numerical problem near podes
               if (xdeg.lt.0.5) then
                  azi2=azi
               else if (xdeg.gt.179.5.and.xdeg.lt.180.5) then
                  azi2=180.-azi
                  if (azi2.lt.0.) azi2=azi2+360.
               else
                  call SPH_AZI(plat2,plon2,plat,plon,del,azi2)
                  azi2=azi2+180.
                  if (azi2.gt.360.) azi2=azi2-360.
               end if
               
               if (idebug.eq.2) then
                  write (18,*) 'L: ',plat,plon,azi,xdeg
                  write (18,*) '   ',plat2,plon2,azi2
               endif                                            

! azimuth calculation has problems near N and S pole
! the following check should exit to next ray if NaN problem occurs
	!print *,'azi2 = ',azi2
      !if (nint(azi2/1000.).ne.0.) then
	  if (isnan(azi2)) then
         print *,'**problem with azi2 updating, skipping to next ray'
         print *,'azi2 = ',azi2
         print *,plat,plon,xdeg,azi
         print *,plat2,plon2,del,azi2
         go to 760
      end if

               plat=plat2
               plon=plon2
               azi=azi2
               x=0.
               distsum=0.

210            if (idir.eq.1) then    !scat point at end of ray path
                  vp0=alpha(ilay+1)
                  vs0=beta(ilay+1)  
               else
                  vp0=alpha(ilay)
                  vs0=beta(ilay)
               end if
               idir2=idir
               azi2=azi
               
               if (idebug.ne.0) then               
                  write (18,*) 'InSCATRAYPOL: ',np,iwave1,svsh,iwave2,ip
                  write (18,*) psi2,zeta2,spol,vp0,vs0
                  write (18,*) idir2,azi2                  
               endif

               call SCATRAYPOL(p,np,iwave1,svsh,iwave2,ip,psi2,zeta2, &
                               spol,vp0,vs0,idir2,azi2)
     
               if (idebug.ne.0) then               
                  write (18,*) 'SCATRAYPOL: ',np,iwave1,svsh,iwave2,ip
                  write (18,*) psi2,zeta2,spol,vp0,vs0
                  write (18,*) idir2,azi2                                               
                  write (18,*) 'S: ',iw,ipp,ip,azi,azi2,plat,plon,idir
                  write (18,*) ' '
               endif
               

               if (idir2.eq.0) then
                  print *,'DEBUG SCATRAYPOL'
                  print *,iwave1,iwave2,ip,vp0,vs0,idir
                  print *,ilay,alpha(ilay),alpha(ilay+1),z_s(ilay)
                  print *,p(ip,iwave1),p(ipp,iw_init)
                  print *,iw_init,idir_init,ipp,nscat
                  idir=-1          !****kluge
                  go to 210        !****kluge to fix rare problem     
                  stop
               end if
                
               azi=azi2               

               iwave=iwave2
               iw=iwave
               if (idir.ne.idir2) then  !switched dir, go through layer again
! I thought next line would fix DEBUG SCATRAYP error but it did not
                  if (idir2.eq.-1.and.iddir(ilay,ip,iw).eq.-1) idir2=1
                  idir=idir2
                  go to 200
               end if

            end if     !end section wheRAN1(idum)om scattering occurred

         end if      !end section on scattering volume
! ------------------------------------ end section on volume scattering
300      continue


!--------------------------------- INTERFACE SECTION
! then deal with interfaces
      if (idir.eq.1) then              !downgoing

         if (iflag(ilay+1).eq.0) then       !no interface to worry about
            ilay=ilay+1
         else
            iface=iflag(ilay+1)

            isave=isave+1                      !****for dumping ray info
            depsave(isave,1)=depface(iface)
            depsave(isave,2)=0.
            iwsave(isave)=iw

            if (iw.eq.1) then  ! incident P-wave

              fran=RAN1(idum)
              rt1=rt(1,1,iw,iw,iface,ip,iw)    !down trans., no conversion
              rt2=rt(1,2,iw,iw,iface,ip,iw)    !reflected, no conversion
              rt3=rt(1,1,iw,3-iw,iface,ip,iw)  !down trans., phase conv.
              rt4=rt(1,2,iw,3-iw,iface,ip,iw)    !reflected, phase conv.
              test=rt1+rt2+rt3+rt4
              if (abs(test-1.).gt.0.01) then
                 print *,'***Warning in downgoing ref/trans coef'
                 print *,iw,iface,ip,depface(iface)
                 print *,rt1,rt2,rt3,rt4,test
                 r=erad-depface(iface)
                 flatfact=r/erad
                 pcor=p(ip,iw)/flatfact
                 print *,p(ip,iw),pcor,1./vp(iface,1),1./vp(iface,2)
                 xdeg=x/kmdeg
                 tmin=t/60.
                 print *,'xdeg,tmin = ',xdeg,tmin
                 print *,iw_init,idir_init,ipp,nsweep,nscat
                 rt2=1.
              end if
              test1=rt1          !down trans., no conversion
              test2=test1+rt2    !reflected, no conversion
              test3=test2+rt3    !down trans., phase conv.

!            if (ip.eq.50) then
!               print *,'**down dump, ip,iw,dep = ',ip,iw,depface(iface)
!               print *,fran,test1,test2,test3
!            end if
              if (idebug.eq.1) then
                 write (18,*) 'Pdown face: ',ip,iw,depface(iface)
              endif

              if (fran.le.test1) then          !transmitted
                 ilay=ilay+2          
              else if (fran.le.test2) then     !reflected
                 idir=-1              
              else if (fran.le.test3) then     !transmitted, converted
                 ilay=ilay+2          
                 iw2=3-iw
                 call UPDATE_P(p,np,iw,iw2,ip,ip2)
                 if (ip2.eq.0) then
                    print *,'DEBUG1: ',ilay,fran,test1,test2,test3
                    stop
                 end if
                 iw=iw2
                 ip=ip2
                 svsh=0.
              else                             !reflected, converted
                 idir=-1
                 iw2=3-iw
                 call UPDATE_P(p,np,iw,iw2,ip,ip2)
                 if (ip2.eq.0) then
                    print *,'DEBUG2: ',ilay,fran,test1,test2,test3
                    stop
                 end if
                 iw=iw2
                 ip=ip2
                 svsh=0.
!                 print *,'P reflected to SV at depth: ',depface(iface)
!                 print *,'iw,ip,svsh = ',iw,ip,svsh

              end if

            else           !incident S-wave
              svfrac=cos(svsh)**2
              shfrac=1.-svfrac
              test1=svfrac*rt(1,1,2,1,iface,ip,iw)         !SV to P trans
              test2=test1+svfrac*rt(1,2,2,1,iface,ip,iw)   !SV to P refl
              test3=test2+svfrac*rt(1,1,2,2,iface,ip,iw) &   !SV to SV trans
                         +shfrac*rt(1,1,3,3,iface,ip,iw)     !SH to SH trans
              fran=RAN1(idum)

!            if (ip.eq.50) then
!               print *,'**down dump, ip,iw,dep = ',ip,iw,depface(iface)
!               print *,fran,test1,test2,test3
!            end if
              if (idebug.eq.1) then
                 write (18,*) 'Sdown face: ',ip,iw,depface(iface)
              endif

              if (fran.le.test1) then          !SV to P trans
                 ilay=ilay+2
                 iw2=3-iw
                 call UPDATE_P(p,np,iw,iw2,ip,ip2)
                 if (ip2.eq.0) then
                    print *,'DEBUG3: ',ilay,fran,test1,test2,test3
                    stop
                 end if
                 iw=iw2
                 ip=ip2
                 svsh=0.                           
              else if (fran.le.test2) then     !SV to P refl
                 idir=-1
                 iw2=3-iw
                 call UPDATE_P(p,np,iw,iw2,ip,ip2)
                 if (ip2.eq.0) then
                    print *,'DEBUG4: ',ilay,fran,test1,test2,test3
                    stop
                 end if
                 iw=iw2
                 ip=ip2
                 svsh=0.
!                 print *,'S reflected to P at depth: ',depface(iface)
!                 print *,'iw,ip,svsh = ',iw,ip,svsh
              else if (fran.le.test3) then     !S to S transmitted
                 ilay=ilay+2
                 svamp=sqrt(svfrac*rt(1,1,2,2,iface,ip,iw))   !SV to SV trans
                 shamp=sqrt(shfrac*rt(1,1,3,3,iface,ip,iw))   !SH to SH trans
                 svsh=atan2(shamp,svamp)
              else
                 idir=-1                       !S to S reflected
                 svamp=sqrt(svfrac*rt(1,2,2,2,iface,ip,iw))   !SV to SV refl
                 shamp=sqrt(shfrac*rt(1,2,3,3,iface,ip,iw))   !SH to SH refl
                 svsh=atan2(shamp,svamp)
              end if

            end if     !P vs. S incident wave

         end if         !interface vs. nointerface

      else                         !---------------upgoing (idir = -1)

         if (ilay.eq.1) then    !at surface, need to output t,x
         
            isave=isave+1                  !****for dumping ray info
            depsave(isave,1)=0.
            depsave(isave,2)=0.
            iwsave(isave)=iw

            tmin=t/60.
            xdeg=x/kmdeg
            iwrap=mod(int(xdeg/180.),2)            
            call SPH_LOC(plat,plon,xdeg,azi,plat2,plon2)
            call SPH_AZI(0.,0.,plat2,plon2,xdeg,azidum)
            
            psecdeg=p(ip,iw)*111.19
            call SPH_AZI(plat2,plon2,0.,0.,xdum,azidum1)      !to source
            call SPH_AZI(plat2,plon2,plat,plon,xdum,azidum2)  !to last scatterer (or origin if none)

            if (iwrap.eq.1) azidum2=azidum2+180.         !wrap-around correction
            slowang=azidum2-azidum1
            psecdegtran=sin(slowang/degrad)*psecdeg
            psecdegrad=cos(slowang/degrad)*psecdeg
            
            if (idebug.eq.1) then
               write (18,*) 'at surface: ',ip,iw,xdeg,tmin
               write (19,*) nray,ip,iw,xdeg,tmin
            endif
                        
!            print *,'Surface p,xdeg,tmin = ',p(ip,iw),xdeg,tmin
!            xdeg=amod(xdeg+3600.,360.)                !wraparound fix
!            if (xdeg.gt.180.) xdeg=360.-xdeg
            ix=nint(xdeg*2.+0.5)
            if (ix.lt.1) ix=1
            if (ix.gt.nxdim) ix=nxdim   
            it=nint(t+0.5)
            if (it.lt.1) it=1
            if (it.gt.ntdim) it=ntdim

            amp=25.*ampstart*exp(-freq*pi*tstar)          !attenuation
            if (amp.lt.1.e-12) amp=0.        !avoid underflow errors
            amp2=amp**2
!            if (xdeg.gt.20..and.amp2.gt.10.) then
!               print *,' '
!               print *,'Big amplitude: ',xdeg,amp2
!               print *,'tmin = ',tmin
!               print *,'ampstart, iw0, iw = ',ampstart,iwstart,iw
!               print *,'ipp,nscat = ',ipp,nscat
!               print *,'starting p = ',p(ipp,iwstart)
!            end if

            if (nscat.ge.nscatmin.and.nscat.le.nscatmax) then
               rbin(it,ix)=rbin(it,ix)+amp2
               if (iw.eq.1) then
                  rbin_p(it,ix)=rbin_p(it,ix)+amp2
                  sinthe=p(ip,iw)*vpmin
                  amp_rad=sinthe*sqrt(amp2)
                  energy_rad=amp_rad**2
                  energy_vert=amp2-energy_rad
                  rbin_z(it,ix)=rbin_z(it,ix)+energy_vert
                  if (psecdeg.le.5.) then
                     rbin_zcore(it,ix)=rbin_zcore(it,ix)+energy_vert
                  endif
                  rbin_rad(it,ix)=rbin_rad(it,ix)+energy_rad

                  if (iforcerefl.eq.-1) then                  
                   do k=1,nslowout
                     if (xdeg.ge.xslow1(k).and.xdeg.le.xslow2(k).and. &
                         t.ge.tslow1(k).and.t.le.tslow2(k)) then
                        i=(nslowdim+1)/2+nint(psecdegrad*2.)
                        j=(nslowdim+1)/2+nint(psecdegtran*2.)
                        if (i.ge.1.and.i.le.nslowdim.and. &
                            j.ge.1.and.j.le.nslowdim) then
                           slowstack(i,j,k)=slowstack(i,j,k)+energy_vert
                        endif
                     endif
                   enddo
                  endif
                  
               else
                  svfrac=cos(svsh)**2
                  shfrac=1.-svfrac
                  rbin_sv(it,ix)=rbin_sv(it,ix)+amp2*svfrac
                  rbin_sh(it,ix)=rbin_sh(it,ix)+amp2*shfrac
                  sinthe=p(ip,iw)*vsmin
                  amp_vert=sinthe*sqrt(amp2)
                  energy_vert=amp_vert**2
                  energy_rad=amp2-energy_vert
                  rbin_z(it,ix)=rbin_z(it,ix)+energy_vert
                  if (psecdeg.le.5.) then
                     rbin_zcore(it,ix)=rbin_zcore(it,ix)+energy_vert
                  endif
                  rbin_rad(it,ix)=rbin_rad(it,ix)+energy_rad
               end if
            end if

            if (xdeg.ge.xwind1.and.xdeg.le.xwind2.and. &
                tmin.ge.twind1.and.tmin.le.twind2) then
               print *,'Match to window, xdeg,tmin,p:'
               print *,xdeg,tmin,p(ip,iw),tstar,amp,nscat
               nsave=isave
               do 8888 isave=1,nsave
                  print *,depsave(isave,1),depsave(isave,2), &
                          iwsave(isave)
8888           continue
            end if

            nsurf=nsurf+1

            if (iw.eq.1) then         !upgoing P-wave at surface

              fran=RAN1(idum)
              iface=1
              if (iforcerefl.ne.2) then
                 test1=rt(2,1,iw,iw,iface,ip,iw)    !reflected, no conversion
              else 
                 test1=1.1
              end if

              if (fran.le.test1) then
                 idir=1
              else                            !phase conversion
                 idir=1
                 iw2=3-iw

                 call UPDATE_P(p,np,iw,iw2,ip,ip2)
                 if (ip2.eq.0) then
                    print *,'DEBUG3: ',ilay,fran,test1
                    stop
                 end if
                 iw=iw2
                 ip=ip2
                 svsh=0.
              end if

            else                     !upgoing S-wave at surface

              svfrac=cos(svsh)**2
              shfrac=1.-svfrac
              iface=1
              test1=svfrac*rt(2,1,2,1,iface,ip,iw)   !SV to P refl
              fran=RAN1(idum)
              if (fran.le.test1) then              !SV to P refl
                 idir=1
                 iw2=3-iw
                 call UPDATE_P(p,np,iw,iw2,ip,ip2)
                 if (ip2.eq.0) then
                    print *,'DEBUG3b: ',ilay,fran,test1,test2,test3
                    stop
                 end if
                 iw=iw2
                 ip=ip2
                 svsh=0.
              else
                 idir=1                       !S to S reflected
                 svamp=sqrt(svfrac*rt(2,1,2,2,iface,ip,iw))   !SV to SV refl
                 shamp=sqrt(shfrac*rt(2,1,3,3,iface,ip,iw))   !SH to SH refl
                 svsh=atan2(shamp,svamp)
              end if

            end if               !P or S at surface
             
         else      !upgoing, not at surface

            if (ilay.eq.2.or.iflag(ilay-1).eq.0) then   !no interface 
               ilay=ilay-1
            else
               iface=iflag(ilay-1)

               isave=isave+1                      !****for dumping ray info
               depsave(isave,1)=depface(iface)
               depsave(isave,2)=0.
               iwsave(isave)=iw

               if (iw.eq.1) then  ! incident P-wave

                 fran=RAN1(idum)
                 test1=rt(2,2,iw,iw,iface,ip,iw)      !up trans., no conversion
                 test2=test1+rt(2,1,iw,iw,iface,ip,iw)  !refl., no conversion
                 test3=test2+rt(2,2,iw,3-iw,iface,ip,iw)  !up trans., phase conv.
                 
              if (idebug.eq.1) then
                 write (18,*) 'Pup face: ',ip,iw,depface(iface)
              else if (idebug.eq.3) then
                  write (18,*) 'Pup face: ',ip,iw,depface(iface)
                  write (18,*) '  ',p(ip,iw),rt(2,1,iw,iw,iface,ip,iw)
              endif


                 if (fran.le.test1) then          !transmitted
                    ilay=ilay-2          
                 else if (fran.le.test2) then     !reflected
                    idir=1              
                 else if (fran.le.test3) then     !transmitted, converted
                    ilay=ilay-2
                    iw2=3-iw
                    call UPDATE_P(p,np,iw,iw2,ip,ip2)
                    if (ip2.eq.0) then
                       print *,'DEBUG6: ',ilay,fran,test1,test2,test3
                       stop
                    end if
                    iw=iw2
                    ip=ip2
                    svsh=0.          
                 else                             !reflected, converted
                    idir=1
                    iw2=3-iw
                    call UPDATE_P(p,np,iw,iw2,ip,ip2)
                    if (ip2.eq.0) then
                       print *,'DEBUG7: ',ilay,fran,test1,test2,test3
                       stop
                    end if
                    iw=iw2
                    ip=ip2
                    svsh=0.
                 end if

               else           !incident S-wave
                 svfrac=cos(svsh)**2
                 shfrac=1.-svfrac
                 test1=svfrac*rt(2,2,2,1,iface,ip,iw)         !SV to P trans
                 test2=test1+svfrac*rt(2,1,2,1,iface,ip,iw)   !SV to P refl
                 test3=test2+svfrac*rt(2,2,2,2,iface,ip,iw) &     !SV to SV trans
                            +shfrac*rt(2,2,3,3,iface,ip,iw)     !SH to SH trans
                 fran=RAN1(idum)
                 
                 if (idebug.eq.1) then
                    write (18,*) 'Sup face: ',ip,iw,depface(iface)
                 endif
                 
                 if (fran.le.test1) then          !SV to P trans
                    ilay=ilay-2
                    iw2=3-iw
                    call UPDATE_P(p,np,iw,iw2,ip,ip2)
                    if (ip2.eq.0) then
                       print *,'DEBUG8: ',ilay,fran,test1,test2,test3
                       stop
                    end if
                    iw=iw2
                    ip=ip2
                    svsh=0.                           
                 else if (fran.le.test2) then     !SV to P refl
                    idir=1
                    iw2=3-iw
                    call UPDATE_P(p,np,iw,iw2,ip,ip2)
                    if (ip2.eq.0) then
                       print *,'DEBUG9: ',ilay,fran,test1,test2,test3
                       stop
                    end if
                    iw=iw2
                    ip=ip2
                    svsh=0.
                 else if (fran.le.test3) then     !S to S transmitted
                    ilay=ilay-2
                    svamp=sqrt(svfrac*rt(2,2,2,2,iface,ip,iw))
                    shamp=sqrt(shfrac*rt(2,2,3,3,iface,ip,iw))
                    svsh=atan2(shamp,svamp)
                 else
                    idir=1                       !S to S reflected
                    svamp=sqrt(svfrac*rt(2,1,2,2,iface,ip,iw))   !SV to SV refl
                    shamp=sqrt(shfrac*rt(2,1,3,3,iface,ip,iw))   !SH to SH refl
                    svsh=atan2(shamp,svamp)
                 end if

               end if   !incident P vs. incident S

            end if    !no interface vs. interface

         end if     !surface vs. not surface check

      end if        !upgoing vs downgoing if      

      if (t.lt.tmax) go to 200

350   continue                    !not end of do loop

      if (mod(nray,nraydump).eq.0) then
         print *,'writing rbin files...   nray,nsurf = ',nray,nsurf

         outfile='out.photon'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)

         outfile='out.photon_p'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin_p(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)

         outfile='out.photon_sv'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin_sv(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)

         outfile='out.photon_sh'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin_sh(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)

         outfile='out.photon_z'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin_z(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)
         
         outfile='out.photon_zcore'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin_zcore(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)         

         outfile='out.photon_rad'
         open (12,file=outfile,form='unformatted')
         write (12) 0.,180.,nxdim,0.,3000.,ntdim, &
            iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
         write (12) dummy
         write (12) ((rbin_rad(it,ix),it=1,ntdim),ix=1,nxdim)
         close (12)
         
         if (iforcerefl.eq.-1) then        !output slowness stacks
            do k=1,nslowout
               write (outfile,'(a8,i1)') 'out.slow',k
               open (12,file=outfile)
               write (12,*) xslow1(k),xslow2(k),tslow1(k),tslow2(k)
               do i=1,nslowdim
                  write (12,397) (slowstack(i,j,k),j=1,nslowdim)
397               format (50e12.4)
               enddo
            enddo
         endif
         
         print *,'Finished writing.'

         n=0
         rsum=0.
         rmax=0.
         do 500 it=1,ntdim
            do 490 ix=1,nxdim
               if (rbin(it,ix).eq.0.) go to 490
               if (rbin(it,ix).gt.rmax) rmax=rbin(it,ix)
               n=n+1
               rsum=rsum+rbin(it,ix)
490         continue
500      continue
         print *,'n,rmax,rsum = ',n,rmax,rsum

      end if

760         continue                !end loop on starting ray param
780      continue                !end loop on starting idir
800   continue                !end loop on starting iw

      go to 150

      stop
	  
      end program psphoton
	  

! UPDATE_P gets new ray index number for P/S or S/P conversion
!  Inputs:    p  =  table of ray parameters for both P and S
!            np  =  number of p values in p array
!            iw1 =  incident wave type (1=P, 2=S)
!            iw2 =  scattered wave type (1=P, 2=S)
!            ip1  =  ray index for incident wave
!  Returns:  ip2  =  ray index for scattered wave
!
!  Note:  Selects closest ray parameter SMALLER than target p to avoid
!         possibility of non-existent ray (for p close to material slowness)
!
      subroutine UPDATE_P(p,np,iw1,iw2,ip1,ip2)
      implicit none
      integer, parameter :: nray0=50000
!      parameter (nray0=15000)
      real    :: p(nray0,2)       !ray table for both P and S
      integer :: np,iw1,iw2,ip1,ip2
      real    :: ptarg,frac
      ptarg=p(ip1,iw1)
      frac=(ptarg-p(1,iw2))/(p(np,iw2)-p(1,iw2))
!      ip2=nint(1.+frac*float(np-1))
      ip2=int(1.+frac*float(np-1))      !err on side of smaller ray param.
      if (p(ip2,iw2).gt.ptarg) ip2=ip2-1
      if (ip2.lt.1) ip2=ip2+1
      if (ip2.gt.np) ip2=ip2-1
      if (ip2.lt.1.or.ip2.gt.np.or.p(ip2,iw2).gt.ptarg) then
         print *,'***ERROR in UPDATE_P'
         print *,iw1,iw2,ip1,ip2,np,frac
         print *,p(ip2,iw2),ptarg,p(1,iw2),p(np,iw2)
         ip2=0
         return
      end if
      return
      end
      


! SCATRAYPOL computes change in ray parameter, ray azimuth,
! and vertical direction for 3-D scattering for both P, S 
! and P/S conversions.  This version includes S polarization.
!    Inputs:  p  =  table of ray parameters for both P and S
!            np  =  number of p values in p array
!            iw1 =  incident wave type (1=P, 2=S)
!            svsh = incident S polarization (radians, measured from SV pol.)
!            iw2 =  scattered wave type (1=P, 2=S)
!            ip  =  ray index for incident wave
!            psi =  scattering angle from ray direction (radians)
!            zeta = scattering angle from x1 axis (initial S polarization dir.)
!            spol = scattered S polarization (radians, 0=pure psi)
!            vp0  =  P velocity at scattering point
!            vs0  =  S velocity at scattering point
!            idir =  incident ray up/down direction (1=downgoing, 2=upgoing)
!            azi  =  azimuth of incident ray (degrees)
!   Returns: ip   =  ray index for scattered wave
!            idir =  scattered wave up/down direction
!            azi  =  azimuth of scattered wave (degrees)
!            svsh =  scattered S polarization (radians)
!
      subroutine SCATRAYPOL(p,np,iw1,svsh,iw2,ip,psi,zeta,spol, &
                            vp0,vs0,idir,azi)
      implicit none
      integer, parameter :: nray0=50000
!      parameter (nray0=15000)      
      real    :: svsh,psi,zeta,spol,vp0,vs0,azi,pi,degrad
      real    :: p(nray0,2)                   !ray table for both P and S
      real    :: vel0(2),slow(2),a1(3),a2(3),a3(3),b1(3),b2(3),b3(3)
      real    :: azi0,svsh0,p1,the,phi,the2,phi2,r,p2,dot,ptarg,frac
      integer :: i,np,iw1,iw2,ip,idir,iret
      
      do i=1,2
         vel0(i)=0.
         slow(i)=0.
      enddo
      do i=1,3
         a1(i)=0.
         a2(i)=0.
         a2(i)=0.
         b1(i)=0.
         b2(i)=0.
         b3(i)=0.
      enddo
      azi0=0.
      svsh0=0.
      p1=0.
      the=0.
      phi=0.
      the2=0.
      phi2=0.
      r=0.
      p2=0.
      dot=0.
      ptarg=0.
      frac=0.
      i=0
      iret=0
                              
      pi=3.1415927
      degrad=180./pi

      vel0(1)=vp0
      vel0(2)=vs0
      slow(1)=1./vel0(1)
      slow(2)=1./vel0(2)
      azi0=azi
      svsh0=svsh
      p1=p(ip,iw1)
      if (slow(iw1).lt.p1) then
         print *,'***ERROR in SCATRAYP: ',slow(iw1),p1,ip,iw1,vp0,vs0
         idir=0
         return
      end if

      the=asin(p1*vel0(iw1))*degrad  !incident ray angle from vert (deg.)
      if (idir.eq.1) the=180.-the
      phi=azi            !incident ray angle from North (degrees)

      call TO_CAR(the,phi,1.,a3(1),a3(2),a3(3))    !a3 is ray direction vector
      the2=abs(the-90.)
      phi2=phi
      if (the.lt.90.) phi2=phi+180.

      call TO_CAR(the2,phi2,1.,a1(1),a1(2),a1(3))  !a1 is SV polarization
      call CROSS(a3,a1,a2)                         !a2 is SH polarization

      call ORTHOCHECK(a1,a2,a3,iret)
      if (iret.ne.0) then
         print *,'***Problem in SCATRAYPOL'
         print *,'   a vectors not unit orthogonal'
         print *,'a1 = ',a1
         print *,'a2 = ',a2
         print *,'a3 = ',a3
         print *,the,phi,the2,phi2
         print *,np,iw1,svsh0,iw2,ip
         print *,psi,zeta,spol,vp0,vs0
         print *,idir,azi0
         stop
      end if

      do 10 i=1,3
         b3(i)=a3(i)                               !b3 is ray direction
         b1(i)= cos(svsh)*a1(i)+sin(svsh)*a2(i)    !b1 is S pol. dir.
         b2(i)=-sin(svsh)*a1(i)+cos(svsh)*a2(i)
10    continue

      call BENDRAY(b1,b2,b3,psi,zeta,spol)        !overwrites b arrays

      call TO_POL(b3(1),b3(2),b3(3),the,phi,r)

      p2=slow(iw2)*sin(the/degrad)
      if (the.lt.90.) then
         idir=2
      else
         idir=1
      end if
      azi=phi

      if (nint(azi/1000.).ne.0.) then
         print *,'**SCATRAYPOL problem with azi'
         print *,'azi = ',azi
         print *,'b1 = ',b1
         print *,'a3 = ',a3
         print *,'a1 = ',a1
         print *,'the,phi = ',the,phi
         print *,np,iw1,svsh0,iw2,ip
         print *,psi,zeta,spol,vp0,vs0
         print *,idir,azi0
         stop
      end if


      the2=the-90.
      phi2=phi
      if (the.lt.90.) phi2=phi+180.
      call TO_CAR(the2,phi2,1.,a1(1),a1(2),a1(3))  !a1 is SV polarization

      call VDOT(a1,b1,dot)
      if (dot.gt.1) dot=1.
      if (dot.lt.-1) dot=-1.
      svsh=acos(dot)               !angular difference between a1 and b1

      if (nint(svsh/10.).ne.0.) then
         print *,'**SCATRAYPOL problem1'
         print *,'svsh,dot = ',svsh,dot
         print *,'b1 = ',b1
         print *,'a3 = ',a3
         print *,'a1 = ',a1
         print *,'the,phi = ',the,phi
         print *,np,iw1,svsh0,iw2,ip
         print *,psi,zeta,spol,vp0,vs0
         print *,idir,azi0
         stop
      end if

      ptarg=p2
      frac=(ptarg-p(1,iw2))/(p(np,iw2)-p(1,iw2))
      if (frac.lt.0.) frac=0.

      ip=nint(1.+frac*float(np-1))

      if (ip.le.0.or.ip.gt.np+20) then
         print *,'**SCATRAYPOL problema'
         print *,ip,np,frac,iw1,iw2
         print *,ptarg,p(1,iw2),p(np,iw2),p(1,iw2)
         print *,the,slow(iw2),b3
         print *,'a3 = ',a3
         print *,'a1 = ',a1
         print *,np,iw1,svsh0,iw2,ip
         print *,psi,zeta,spol,vp0,vs0
         print *,idir,azi0
         stop
      end if

      if (p(ip,iw2).gt.slow(iw2)) ip=ip-1
      if (p(ip,iw2).gt.slow(iw2)) ip=ip-1
      if (ip.lt.1.or.ip.gt.np) then
         print *,'*WARNING SCATRAYPOL: ',ip,iw2,slow,ptarg, &
                    p(1,iw2),p(np,iw2)
      end if

      return
      end



      subroutine TO_POL(x,y,z,the,phi,r)
      implicit none
      real :: x,y,z,the,phi,r,degrad
      degrad=3.1415927/180.
      r=sqrt(x**2+y**2+z**2)
      the=acos(z/r)/degrad
      if (x.ne.0..or.y.ne.0.) then
         phi=atan2(y,x)/degrad
      else
         phi=0.
      end if
      return
      end

      subroutine TO_CAR(the,phi,r,x,y,z)
      implicit none
      real :: the,phi,r,x,y,z,degrad
      degrad=3.1415927/180.
      z=r*cos(the*degrad)
      x=r*sin(the*degrad)*cos(phi*degrad)
      y=r*sin(the*degrad)*sin(phi*degrad)
      return
      end

      subroutine ORTHOCHECK(v1,v2,v3,iret)
      implicit none
      real    :: v1(3),v2(3),v3(3),dot12,dot13,dot23,d1,d2,d3
      integer :: iret
      call VDOT(v1,v2,dot12)
      call VDOT(v1,v3,dot13)
      call VDOT(v2,v3,dot23)
      d1=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
      d2=sqrt(v2(1)**2+v2(2)**2+v2(3)**2)
      d3=sqrt(v3(1)**2+v3(2)**2+v3(3)**2)
      iret=0                                !orthogonal
      if (abs(dot12).gt.0.0001) iret=1      !not orthogonal
      if (abs(dot13).gt.0.0001) iret=1
      if (abs(dot23).gt.0.0001) iret=1
      if (abs(d1-1.).gt.0.0001) iret=1
      if (abs(d2-1.).gt.0.0001) iret=1
      if (abs(d3-1.).gt.0.0001) iret=1
      return
      end

!
!
!
! FLATTEN calculates flat earth tranformation.
      subroutine FLATTEN(z_s,vel_s,z_f,vel_f)
      implicit none
      real :: z_s,vel_s,z_f,vel_f,erad,r
      erad=6371.
      r=erad-z_s
      z_f=-erad*alog(r/erad)
      vel_f=vel_s*(erad/r)
      return
      end


! LAYERTRACE calculates the travel time and range offset
! for ray tracing through a single layer.
!
! Input:    p     =  horizontal slowness
!           h     =  layer thickness
!           utop  =  slowness at top of layer
!           ubot  =  slowness at bottom of layer
!           imth  =  interpolation method
!                    imth = 1,  v(z) = 1/sqrt(a - 2*b*z)
!                         = 2,  v(z) = a - b*z
!                         = 3,  v(z) = a*exp(-b*z)
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
      implicit none
      real*4  :: p1,h1,utop1,ubot1,dx1,dt1
      real*8  :: p,h,utop,ubot,u,y,q,qs,qr,b,vtop,vbot,etau,ex,z,dx,dtau,dt
      integer :: irtr,imth
      p=dble(p1)
      h=dble(h1)
      utop=dble(utop1)
      ubot=dble(ubot1)

      if (h.eq.0.) then      !check for zero thickness layer
         dx1=0.
         dt1=0.
         irtr=-1
         return         
      end if

      u=utop
      y=u-p
      if (y.le.0.) then   !complex vertical slowness
         dx1=0.
         dt1=0.
         irtr=0
         return
      end if

      q=y*(u+p)
      qs=dsqrt(q)

! special function needed for integral at top of layer
      if (imth.eq.2) then
         y=u+qs
         if (p.ne.0.) y=y/p
         qr=dlog(y)
      else if (imth.eq.3) then
         qr=atan2(qs,p)
      end if      

      if (imth.eq.1) then
          b=-(utop**2-ubot**2)/(2.*h)
      else if (imth.eq.2) then
          vtop=1./utop
          vbot=1./ubot
          b=-(vtop-vbot)/h
      else
          b=-dlog(ubot/utop)/h
      end if  

      if (b.eq.0.) then     !constant velocity layer
         b=1./h
         etau=qs
         ex=p/qs
         irtr=1
         go to 160
      end if

! integral at upper limit, 1/b factor omitted until end
      if (imth.eq.1) then
         etau=-q*qs/3.
         ex=-qs*p
      else if (imth.eq.2) then
         ex=qs/u                    !*** - in some versions (wrongly)
         etau=qr-ex
         if (p.ne.0.) ex=ex/p
      else
         etau=qs-p*qr
         ex=qr
      end if

! check lower limit to see if we have turning point
      u=ubot
      if (u.le.p) then   !if turning point,
         irtr=2          !then no contribution
         go to 160       !from bottom point
      end if 
      irtr=1
      q=(u-p)*(u+p)
      qs=dsqrt(q)

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

160   dx=ex/b
      dtau=etau/b
      dt=dtau+p*dx     !convert tau to t

      dx1=sngl(dx)
      dt1=sngl(dt)
      return
      end


! FORCE_RT forces RT coefficients to identical values
!   iforce=1 equalizes all scattered waves (0.25 for P/SV if all present)
!   iforce=2 for no conversions, equal scat. waves of same type
!
!  Note:  Any input rt that are zero are kept at zero
!         Remaining rt are scaled to still sum to 1
!
      subroutine FORCE_RT(rt,iforce)
      implicit none
      real    :: rt(2,2,3,3),rtfix
      integer :: iforce,i,j,k,l,n

      do 100 i=1,2               !do P/SV
      do 90 k=1,2
         if (iforce.eq.2) then
            rt(i,1,k,3-k)=0.
            rt(i,2,k,3-k)=0.
         end if
         n=0
         do 30 j=1,2
         do 20 l=1,2
            if (rt(i,j,k,l).ne.0) n=n+1
20       continue
30       continue
         if (n.ne.0) then
            rtfix=1./float(n)
         else
            rtfix=0.
         end if
         do 50 j=1,2
         do 40 l=1,2
            if (rt(i,j,k,l).ne.0) rt(i,j,k,l)=rtfix
40       continue
50       continue
90    continue
100   continue

      do 200 i=1,2               !do SH
         k=3
         n=0
         do 130 j=1,2
            l=3
            if (rt(i,j,k,l).ne.0) n=n+1
130      continue
         if (n.ne.0) then
            rtfix=1./float(n)
         else
            rtfix=0.
         end if
         do 150 j=1,2
            l=3
            if (rt(i,j,k,l).ne.0) rt(i,j,k,l)=rtfix
150      continue
200   continue         

      return
      end


! subroutine RTCOEF_POW calculates reflection/transmission coefficients
! for interface between two solid layers, based on the equations on 
! p. 150 of Aki and Richards.  This version returns energy-normalized
! power as a real number.
!
!  Inputs:    vp1     =  P-wave velocity of layer 1 (top layer)
!  (real)     vs1     =  S-wave velocity of layer 1
!             den1    =  density of layer 1
!             vp2     =  P-wave velocity of layer 2 (bottom layer)
!             vs2     =  S-wave velocity of layer 2
!             den2    =  density of layer 2
!             p       =  horizontal slowness (ray parameter)
!  Returns:   refcoef =  2x2x3x3 array with coefficients
!                        i = incident wave direction (1=down, 2=up)
!                        j = scattered wave direction (1=down, 2=up)
!                        k = incident wave type (1=P, 2=SV, 3=SH)
!                        l = scattered wave type (1=P, 2=SV, 3=SH)
!
!  Requires:  RTCOEF, RTCOEF_SH
!
      subroutine RTCOEF_POW(vp1,vs1,den1,vp2,vs2,den2,p,refcoef)
      implicit none
      real :: vp1,vs1,den1,vp2,vs2,den2,p
      real :: v(2,3),den(2),refcoef(2,2,3,3)
      complex :: rt(16),rtsh(4)
      integer :: idir1,idir2,ilay1,ilay2,iw1,iw2,i,j,k,l
      real :: slow1,slow2,cos1,cos2,f1,f2,rtnorm,sum
      v(1,1)=vp1
      v(1,2)=vs1
      v(1,3)=vs1
      den(1)=den1
      v(2,1)=vp2
      v(2,2)=vs2
      v(2,3)=vs2
      den(2)=den2

      call RTCOEF(vp1,vs1,den1,vp2,vs2,den2,p,rt)

      call RTCOEF_SH(vs1,den1,vs2,den2,p,rtsh)

      do i=1,2
        do j=1,2
          do k=1,3
            do l=1,3
              refcoef(i,j,k,l)=0.
            enddo
          enddo
        enddo
      enddo

      refcoef(1,2,1,1)=cabs(rt(1))   !down P to P up
      refcoef(1,2,1,2)=cabs(rt(2))   !down P to SV up
      refcoef(1,1,1,1)=cabs(rt(3))   !down P to P down
      refcoef(1,1,1,2)=cabs(rt(4))   !down P to SV down
      refcoef(1,2,2,1)=cabs(rt(5))    !down SV to P up
      refcoef(1,2,2,2)=cabs(rt(6))    !down SV to SV up
      refcoef(1,1,2,1)=cabs(rt(7))    !down SV to P down
      refcoef(1,1,2,2)=cabs(rt(8))    !down SV to SV down
      refcoef(2,2,1,1)=cabs(rt(9))   !up P to P up
      refcoef(2,2,1,2)=cabs(rt(10))  !up P to SV up
      refcoef(2,1,1,1)=cabs(rt(11))  !up P to P down
      refcoef(2,1,1,2)=cabs(rt(12))  !up P to SV down
      refcoef(2,2,2,1)=cabs(rt(13))   !up SV to P up
      refcoef(2,2,2,2)=cabs(rt(14))   !up SV to SV up
      refcoef(2,1,2,1)=cabs(rt(15))   !up SV to P down
      refcoef(2,1,2,2)=cabs(rt(16))   !up SV to SV down

      refcoef(1,2,3,3)=cabs(rtsh(1))     !down SH to SH up 
      refcoef(1,1,3,3)=cabs(rtsh(2))     !down SH to SH down
      refcoef(2,2,3,3)=cabs(rtsh(3))     !up SH to SH up
      refcoef(2,1,3,3)=cabs(rtsh(4))     !up SH to SH down

      do 40 idir1=1,2
         ilay1=idir1
         do 30 idir2=1,2
            ilay2=3-idir2
            do 20 iw1=1,3
               do 10 iw2=1,3
                  if (refcoef(idir1,idir2,iw1,iw2).eq.0.) go to 10
                  slow1=1./v(ilay1,iw1)
                  slow2=1./v(ilay2,iw2)
                  if (p.ge.slow1.or.p.ge.slow2) then
                     refcoef(idir1,idir2,iw1,iw2)=0.
                     go to 10
                  end if

          cos1=sqrt(1.-p**2*v(ilay1,iw1)**2)
          cos2=sqrt(1.-p**2*v(ilay2,iw2)**2)
          f1=v(ilay1,iw1)*den(ilay1)*cos1
          f2=v(ilay2,iw2)*den(ilay2)*cos2

          rtnorm=refcoef(idir1,idir2,iw1,iw2)
          rtnorm=(f2/f1)*rtnorm**2
          refcoef(idir1,idir2,iw1,iw2)=rtnorm

10             continue
20          continue
30       continue
40    continue

! the following section is a kluge to eliminate transmitted
! energy at surfaces for the case where we have used very small
! velocities to simulate air or fluid.
      do 60 idir1=1,2
         ilay1=idir1
         do 50 idir2=1,2
            ilay2=3-idir2
            if (v(ilay2,1).lt.0.2) then     !small P velocity
               refcoef(idir1,idir2,1,1)=0.
               refcoef(idir1,idir2,2,1)=0.
            end if
            if (v(ilay2,2).lt.0.1) then     !small SV velocity
               refcoef(idir1,idir2,1,2)=0.
               refcoef(idir1,idir2,2,2)=0.
            end if
            if (v(ilay2,3).lt.0.1) then     !small SH velocity
               refcoef(idir1,idir2,3,3)=0.
            end if
50       continue
60    continue

      do 150 i=1,2    !now renormalize P/SV so sum of energy is one
      do 140 k=1,2
         sum=0.
         do 80 j=1,2
         do 70 l=1,2
            sum=sum+refcoef(i,j,k,l)
70       continue
80       continue
         if (sum.eq.0.) go to 140
         do 100 j=1,2
         do 90 l=1,2
            refcoef(i,j,k,l)=refcoef(i,j,k,l)/sum
90       continue
100      continue
140   continue
150   continue

      do 250 i=1,2    !now renormalize SH so sum of energy is one
         k=3
         sum=0.
         do 180 j=1,2
            l=3
            sum=sum+refcoef(i,j,k,l)
180      continue
         if (sum.eq.0.) go to 250
         do 200 j=1,2
            l=3
            refcoef(i,j,k,l)=refcoef(i,j,k,l)/sum
200      continue
250   continue

      return
      end


! subroutine RTCOEF calculates reflection/transmission coefficients
! for interface between two solid layers, based on the equations on 
! p. 150 of Aki and Richards.
!
!  Inputs:    vp1     =  P-wave velocity of layer 1 (top layer)
!  (real)     vs1     =  S-wave velocity of layer 1
!             den1    =  density of layer 1
!             vp2     =  P-wave velocity of layer 2 (bottom layer)
!             vs2     =  S-wave velocity of layer 2
!             den2    =  density of layer 2
!             hslow   =  horizontal slowness (ray parameter)
!  Returns:   rt(1)   =  down P to P up     (refl)
!  (complex)  rt(2)   =  down P to S up     (refl)
!             rt(3)   =  down P to P down   (tran)
!             rt(4)   =  down P to S down   (tran)
!             rt(5)   =  down S to P up     (refl)
!             rt(6)   =  down S to S up     (refl)
!             rt(7)   =  down S to P down   (tran)
!             rt(8)   =  down S to S down   (tran)
!             rt(9)   =    up P to P up     (tran)
!             rt(10)  =    up P to S up     (tran)
!             rt(11)  =    up P to P down   (refl)
!             rt(12)  =    up P to S down   (refl)
!             rt(13)  =    up S to P up     (tran)
!             rt(14)  =    up S to S up     (tran)
!             rt(15)  =    up S to P down   (refl)
!             rt(16)  =    up S to S down   (refl)
!
! NOTE:  All input variables are real.  
!        All output variables are complex!
!        Coefficients are not energy normalized.
!
      SUBROUTINE RTCOEF(vp1,vs1,den1,vp2,vs2,den2,hslow,rt)
      implicit none
      real    :: vp1,vs1,den1,vp2,vs2,den2,hslow
      complex :: rt(16),alpha1,beta1,rho1,alpha2,beta2,rho2,p,cone,ctwo
      complex :: si1,si2,sj1,sj2,ci1,ci2,cj1,cj2,term1,term2,a,b,c,d
      complex :: e,f,g,h,den,trm1,trm2

      alpha1=cmplx(vp1,0.)
      beta1=cmplx(vs1,0.)
      rho1=cmplx(den1,0.)
      alpha2=cmplx(vp2,0.)
      beta2=cmplx(vs2,0.)
      rho2=cmplx(den2,0.)
      p=cmplx(hslow,0.)

      cone=cmplx(1.,0.)
      ctwo=cmplx(2.,0.)

      si1=alpha1*p
      si2=alpha2*p          
      sj1=beta1*p
      sj2=beta2*p          

      ci1=csqrt(cone-si1**2)
      ci2=csqrt(cone-si2**2)
      cj1=csqrt(cone-sj1**2)
      cj2=csqrt(cone-sj2**2)         

      term1=(cone-ctwo*beta2*beta2*p*p)
      term2=(cone-ctwo*beta1*beta1*p*p)
      a=rho2*term1-rho1*term2
      b=rho2*term1+ctwo*rho1*beta1*beta1*p*p
      c=rho1*term2+ctwo*rho2*beta2*beta2*p*p
      d=ctwo*(rho2*beta2*beta2-rho1*beta1*beta1)
      E=b*ci1/alpha1+c*ci2/alpha2
      F=b*cj1/beta1+c*cj2/beta2
      G=a-d*ci1*cj2/(alpha1*beta2)
      H=a-d*ci2*cj1/(alpha2*beta1)
      DEN=E*F+G*H*p*p

      trm1=b*ci1/alpha1-c*ci2/alpha2          
      trm2=a+d*ci1*cj2/(alpha1*beta2)
      rt(1)=(trm1*F-trm2*H*p*p)/DEN              !refl down P to P up

      trm1=a*b+c*d*ci2*cj2/(alpha2*beta2)       
      rt(2)=(-ctwo*ci1*trm1*p)/(beta1*DEN)       !refl down P to S up
 
      rt(3)=ctwo*rho1*ci1*F/(alpha2*DEN)         !trans down P to P down

      rt(4)=ctwo*rho1*ci1*H*p/(beta2*DEN)        !trans down P to S down

      trm1=a*b+c*d*ci2*cj2/(alpha2*beta2)       
      rt(5)=(-ctwo*cj1*trm1*p)/(alpha1*DEN)      !refl down S to P up

      trm1=b*cj1/beta1-c*cj2/beta2               
      trm2=a+d*ci2*cj1/(alpha2*beta1)
      rt(6)=-(trm1*E-trm2*G*p*p)/DEN             !refl down S to S up

      rt(7)=-ctwo*rho1*cj1*G*p/(alpha2*DEN)      !trans down S to P down 

      rt(8)=ctwo*rho1*cj1*E/(beta2*DEN)          !trans down S to S down


      trm1=b*ci1/alpha1-c*ci2/alpha2          
      trm2=a+d*ci2*cj1/(alpha2*beta1)
      rt(11)=-(trm1*F+trm2*G*p*p)/DEN            !refl up P to P down

      trm1=a*c+b*d*ci1*cj1/(alpha1*beta1)       
      rt(12)=(ctwo*ci2*trm1*p)/(beta2*DEN)       !refl up P to S down

      rt(9)=ctwo*rho2*ci2*F/(alpha1*DEN)         !trans up P to P up

      rt(10)=-ctwo*rho2*ci2*G*p/(beta1*DEN)      !trans up P to S up

      trm1=a*c+b*d*ci1*cj1/(alpha1*beta1)       
      rt(15)=(ctwo*cj2*trm1*p)/(alpha2*DEN)      !refl up S to P down

      trm1=b*cj1/beta1-c*cj2/beta2               
      trm2=a+d*ci1*cj2/(alpha1*beta2)
      rt(16)=(trm1*E+trm2*H*p*p)/DEN             !refl up S to S down
                                               
      rt(13)=ctwo*rho2*cj2*H*p/(alpha1*DEN)      !trans up S to P up

      rt(14)=ctwo*rho2*cj2*E/(beta1*DEN)         !trans up S to S up


      return
      end



! subroutine RTCOEF_SH calculates SH reflection/transmission coefficients
! for interface between two solid layers, based on the equations on 
! p. 144 of Aki and Richards.
!
!  Inputs:    vs1     =  S-wave velocity of layer 1 (top layer)
!  (real)     den1    =  density of layer 1
!             vs2     =  S-wave velocity of layer 2
!             den2    =  density of layer 2
!             hslow   =  horizontal slowness (ray parameter)
!  Returns:   rt(1)   =  down S to S up     (refl)
!  (complex)  rt(2)   =  down S to S down   (tran)
!             rt(3)  =   up S to S up       (tran)
!             rt(4)  =   up S to S down     (refl)
!
! NOTE:  All input variables are real.  
!        All output variables are complex!
!        Coefficients are not energy normalized.
!
      SUBROUTINE RTCOEF_SH(vs1,den1,vs2,den2,hslow,rt)
      implicit none
      real    :: vs1,den1,vs2,den2,hslow
      complex :: rt(4),beta1,rho1,beta2,rho2,p,cone,sj1,sj2,cj1,cj2,d

      beta1=cmplx(vs1,0.)
      rho1=cmplx(den1,0.)
      beta2=cmplx(vs2,0.)
      rho2=cmplx(den2,0.)
      p=cmplx(hslow,0.)

      cone=cmplx(1.,0.)

      sj1=beta1*p
      sj2=beta2*p          
      cj1=csqrt(cone-sj1**2)
      cj2=csqrt(cone-sj2**2) 

      D=rho1*beta1*cj1+rho2*beta2*cj2 
      rt(1)=(rho1*beta1*cj1-rho2*beta2*cj2)/D
      rt(4)=-rt(1)
      rt(2)=2.*rho1*beta1*cj1/D
      rt(3)=2.*rho2*beta2*cj2/D

      return
      end       




! SPH_LOC finds location of second point on sphere, given range 
! and azimuth at first point.
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          del    =  angular separation between points (degrees)
!          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
! Returns: flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
!
      subroutine SPH_LOC(flat1,flon1,del,azi,flat2,flon2)
      implicit none
      real :: flat1,flon1,del,azi,flat2,flon2,pi,raddeg,delr,azr
      real :: theta1,phi1,ctheta2,theta2,phi2,sphi2,cphi2               
      if (del.eq.0.) then
         flat2=flat1
         flon2=flon1
         return
      end if
      pi=3.141592654
      raddeg=pi/180.
      delr=del*raddeg
      azr=azi*raddeg
      theta1=(90.-flat1)*raddeg
      phi1=flon1*raddeg      
      ctheta2=sin(delr)*sin(theta1)*cos(azr)+cos(theta1)*cos(delr)
      theta2=acos(ctheta2)
      if (theta1.eq.0.) then
         phi2=azr
      else if (theta2.eq.0.) then
         phi2=0.0
      else
         sphi2=sin(delr)*sin(azr)/sin(theta2)
         cphi2=(cos(delr)-cos(theta1)*ctheta2)/(sin(theta1)*sin(theta2))
         phi2=phi1+atan2(sphi2,cphi2)
      end if
      flat2=90.-theta2/raddeg
      flon2=phi2/raddeg
      if (flon2.gt.360.) flon2=flon2-360.
      if (flon2.lt.0.) flon2=flon2+360.
      return
      end


! SPH_AZI computes distance and azimuth between two points on sphere
!
! Inputs:  flat1  =  latitude of first point (degrees) 
!          flon2  =  longitude of first point (degrees)
!          flat2  =  latitude of second point (degrees)
!          flon2  =  longitude of second point (degrees)
! Returns: del    =  angular separation between points (degrees)
!          azi    =  azimuth at 1st point to 2nd point, from N (deg.)
!
! Note:  This routine is inaccurate for del less than about 0.5 degrees. 
!        For greater accuracy, use SPH_AZIDP or perform a separate
!        calculation for close ranges using Cartesian geometry.
!
      subroutine SPH_AZI(flat1,flon1,flat2,flon2,del,azi)
      implicit none
      real :: flat1,flon1,flat2,flon2,del,azi,pi,raddeg
      real :: theta1,theta2,phi1,phi2,stheta1,stheta2,ctheta1,ctheta2
      real :: cang,ang,sang,caz,saz,az 
      if ((flat1.eq.flat2.and.flon1.eq.flon2).or.  &
          (flat1.eq.90..and.flat2.eq.90.).or.  &
          (flat1.eq.-90..and.flat2.eq.-90.))  then
         del=0.
         azi=0.
         return
      end if
      pi=3.141592654
      raddeg=pi/180.
      theta1=(90.-flat1)*raddeg
      theta2=(90.-flat2)*raddeg
      phi1=flon1*raddeg
      phi2=flon2*raddeg
      stheta1=sin(theta1)
      stheta2=sin(theta2)
      ctheta1=cos(theta1)
      ctheta2=cos(theta2)
      cang=stheta1*stheta2*cos(phi2-phi1)+ctheta1*ctheta2
      ang=acos(cang)
      del=ang/raddeg
      sang=sqrt(1.-cang*cang)
      if (sang.eq.0.or.stheta1.eq.0.) then
!         print *,'***SPH_AZI zero divide error'
!         print *,sang,stheta1,flat1,flon1,flat2,flon2
!         stop
          del=0.
          azi=0.
          return
      endif      
      caz=(ctheta2-ctheta1*cang)/(sang*stheta1)
      saz=-stheta2*sin(phi1-phi2)/sang
      az=atan2(saz,caz)
      azi=az/raddeg
      if (azi.lt.0.) azi=azi+360.
      return
      end


! ANGTRAN computes new direction (theta2,phi2), given initial
! direction (theta1,phi1) and change in direction (theta,phi).
! All angles are in radians.
!
      subroutine ANGTRAN(theta1,phi1,theta,phi,theta2,phi2)
      implicit none
      real :: theta1,phi1,theta,phi,theta2,phi2,delr,azr
      real :: ctheta2,sphi2,cphi2
      if (theta.eq.0.) then
         theta2=theta1
         phi2=phi1
         return
      end if
      delr=theta
      azr=phi
      ctheta2=sin(delr)*sin(theta1)*cos(azr)+cos(theta1)*cos(delr)
      theta2=acos(ctheta2)
      if (theta1.eq.0.) then
         phi2=azr
      else if (theta2.eq.0.) then
         phi2=0.0
      else
         sphi2=sin(delr)*sin(azr)/sin(theta2)
         cphi2=(cos(delr)-cos(theta1)*ctheta2)/(sin(theta1)*sin(theta2))
         phi2=phi1+atan2(sphi2,cphi2)
      end if
      if (phi2.lt.0.) phi2=phi2+6.283185307
      if (phi2.gt.6.283185307) phi2=phi2-6.283185307
      return
      end



! RANDANG2 computes random scattering angle using Sato and Fehler
! equations.  This version works for multiple scattering layers
!   Input: iscat  = scattering layer number (max=6)
!          itype  = 1 for gpp
!                 = 2 for gps
!                 = 3 for gsp
!                 = 4 for gss
!            el   =  S-wave wavenumber (=om/beta0)
!            nu   =  density vs. velocity pert. scaling (see 4.48)
!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
!            eps  =  RMS velocity perturbation
!            a    =  correlation distance
!   Returns: psi  =  spherical coor. angle (radians)
!            zeta =  sph. coor angle from x3 axis (radians)
!            spol =  scattered S polarization (radians)
!                 =  0 for pure psi direction
!                    (only meaningful for gps and gss)
!
! Note:  Program computes g arrays upon first call.  Subsequent
!        calls use this array and will not be correct if 
!        el, nu, etc., are changed
!
      subroutine RANDANG2(iscat,itype,el,nu,gam0,eps,a,psi,zeta,spol)
      implicit none
      integer, parameter :: n=100,npts0=66000  !n=100 for 1.8 degree spacing
      integer, parameter :: nscat0=6
      logical :: firstcall(nscat0)
      integer :: iscat,itype,i,j,k,nk,it,k1,k2,idum
      real :: el,nu,gam0,eps,a,sumg1,sumg2,sumg3,sumg4,area,fran,ran1
      real :: psi,zeta,spol,pi,psi2(npts0,nscat0),zeta2(npts0,nscat0)
      real :: g(npts0,4,nscat0),sumg(npts0,4,nscat0),spol2(npts0,nscat0)
      save firstcall,psi2,zeta2,sumg,nk
      data firstcall/.true.,.true.,.true.,.true.,.true.,.true./
      data pi/3.141592754/
      
      idum=0            
      
      if (firstcall(iscat)) then
         print *,'first call  RANDANG2: setting up iscat arrays:', iscat
         firstcall(iscat)=.false.
         k=0
         sumg1=0.
         sumg2=0.
         sumg3=0.
         sumg4=0.     
         do 20 i=0,n
            do 10 j=-n,n
               k=k+1
               psi2(k,iscat) = float(i)*pi/float(n)
               zeta2(k,iscat) = float(j)*pi/float(n)

               call GSATO(psi2(k,iscat),zeta2(k,iscat),el,nu,gam0,eps,a, &
                g(k,1,iscat),g(k,2,iscat),g(k,3,iscat),g(k,4,iscat), &
                spol2(k,iscat))

               area=sin(psi2(k,iscat))*pi**2/float(n)**2
               sumg1=sumg1+area*g(k,1,iscat)
               sumg2=sumg2+area*g(k,2,iscat)
               sumg3=sumg3+area*g(k,3,iscat)
               sumg4=sumg4+area*g(k,4,iscat)
               sumg(k,1,iscat)=sumg1
               sumg(k,2,iscat)=sumg2
               sumg(k,3,iscat)=sumg3
               sumg(k,4,iscat)=sumg4
10          continue
20       continue
         nk=k
         do 40 k=1,nk
            sumg(k,1,iscat)=sumg(k,1,iscat)/sumg1
            sumg(k,2,iscat)=sumg(k,2,iscat)/sumg2
            sumg(k,3,iscat)=sumg(k,3,iscat)/sumg3
            sumg(k,4,iscat)=sumg(k,4,iscat)/sumg4
40       continue
         print *,'Finished setting up arrays'
      end if

      fran=RAN1(idum)

      k1=1
      k2=nk
      do 50 it=1,16
         k=(k1+k2)/2
         if (fran.lt.sumg(k,itype,iscat)) then
            k2=k
         else if (fran.gt.sumg(k,itype,iscat)) then
            k1=k
         else
            k2=k
         end if
50    continue
!      print *,'k,k1,k2 = ',k,k1,k2
!      print *,'  ',fran,sumg(k1,itype,iscat),sumg(k2,itype,iscat)
      zeta=zeta2(k2,iscat)
      psi=psi2(k2,iscat)
      spol=0.
      if (itype.eq.4) spol=spol2(k2,iscat)

      return
      end



! GAVGSATO computes average over solid angle of g functions
! in (4.52) of Sato and Fehler (exponential autocor)
! Inputs:    el   =  S-wave wavenumber (=om/beta0)
!            nu   =  density vs. velocity pert. scaling (see 4.48)
!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
!            eps  =  RMS velocity perturbation
!            a    =  correlation distance
!   Returns: gpp0,gps0,gsp0,gss0  =  averaged over solid angle
!
       subroutine GAVGSATO(el,nu,gam0,eps,a, &
             gpp0,gps0,gsp0,gss0)
      implicit none
      real :: psi,zeta,el,nu,gam0,eps,a,gpp,gps,gsp,gss,spol
      real :: pi,pi4,el4,gam2
      real :: area,sum1,sumgpp,sumgps,sumgsp,sumgss
      real :: gpp0,gps0,gsp0,gss0
      integer n,i,j
      n=60
      pi=3.141592754
      pi4=4.*pi
      el4=el**4
      gam2=gam0**2

      sum1=0.
      sumgpp=0.
      sumgps=0.
      sumgsp=0.
      sumgss=0.      
      do 100 i=0,n
         psi = float(i)*pi/float(n)
         do 90 j=-n,n
            zeta = float(j)*pi/float(n)

            call GSATO(psi,zeta,el,nu,gam0,eps,a, &
             gpp,gps,gsp,gss,spol)

            area=sin(psi)*pi**2/float(n)**2
            sum1=sum1+area
            sumgpp=sumgpp+area*gpp
            sumgps=sumgps+area*gps
            sumgsp=sumgsp+area*gsp
            sumgss=sumgss+area*gss
90       continue
100   continue
      print *,'test one = ',sum1/pi4
      gpp0=sumgpp/sum1
      gps0=sumgps/sum1
      gsp0=sumgsp/sum1
      gss0=sumgss/sum1

      return
      end


! GAVGSATO2 computes average over solid angle of g functions
! in (4.52) of Sato and Fehler (exponential autocor)
!
! This version also does momentum scattering coef.
!
! Inputs:    el   =  S-wave wavenumber (=om/beta0)
!            nu   =  density vs. velocity pert. scaling (see 4.48)
!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
!            eps  =  RMS velocity perturbation
!            a    =  correlation distance
!   Returns: gpp0,gps0,gsp0,gss0  =  averaged over solid angle
!            gppm,gpsm,gspm,gssm  =  momentum scattering coef.
!
       subroutine GAVGSATO2(el,nu,gam0,eps,a, &
             gpp0,gps0,gsp0,gss0, &
             gppm,gpsm,gspm,gssm)
      implicit none
      real :: psi,zeta,el,nu,gam0,eps,a,gpp,gps,gsp,gss,spol
      real :: pi,pi4,el4,gam2
      real :: area,sum1,sumgpp,sumgps,sumgsp,sumgss
      real :: sum2,sumgpp2,sumgps2,sumgsp2,sumgss2,fact
      real :: gpp0,gps0,gsp0,gss0,gppm,gpsm,gspm,gssm
      integer :: n,i,j
      n=60
      pi=3.141592754
      pi4=4.*pi
      el4=el**4
      gam2=gam0**2

      sum1=0.
      sum2=0.
      sumgpp=0.
      sumgps=0.
      sumgsp=0.
      sumgss=0.
      sumgpp2=0.
      sumgps2=0.
      sumgsp2=0.
      sumgss2=0.      
      do 100 i=0,n
         psi = float(i)*pi/float(n)
         do 90 j=-n,n
            zeta = float(j)*pi/float(n)

            call GSATO(psi,zeta,el,nu,gam0,eps,a, &
             gpp,gps,gsp,gss,spol)

            area=sin(psi)*pi**2/float(n)**2
            sum1=sum1+area
            sumgpp=sumgpp+area*gpp
            sumgps=sumgps+area*gps
            sumgsp=sumgsp+area*gsp
            sumgss=sumgss+area*gss
            fact=1.-cos(psi)
            sum2=sum2+area*fact
            sumgpp2=sumgpp2+area*gpp*fact
            sumgps2=sumgps2+area*gps*fact
            sumgsp2=sumgsp2+area*gsp*fact
            sumgss2=sumgss2+area*gss*fact
90       continue
100   continue
      print *,'test one = ',sum1/pi4
      print *,'test one (2) = ',sum2/pi4
      gpp0=sumgpp/sum1
      gps0=sumgps/sum1
      gsp0=sumgsp/sum1
      gss0=sumgss/sum1
      gppm=sumgpp2/sum2
      gpsm=sumgps2/sum2
      gspm=sumgsp2/sum2
      gssm=sumgss2/sum2

      return
      end




! GSATO computes (4.52) from Sato and Fehler (exponential autocor)
!   Inputs:  psi  =  spherical coor. angle (radians)
!            zeta =  sph. coor angle from x3 (radians)
!            el   =  S-wave wavenumber (=om/beta0)
!            nu   =  density vs. velocity pert. scaling (see 4.48)
!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
!            eps  =  RMS velocity perturbation
!            a    =  correlation distance
!   Returns: gpp,gps,gsp,gss  =  from eqn. (4.52)
!            spol =  S-to-S scattered S polarization (radians)
!                 =  0 for pure psi direction
!
       subroutine GSATO(psi,zeta,el,nu,gam0,eps,a, &
             gpp,gps,gsp,gss,spol)
      implicit none
      real :: psi,zeta,el,nu,gam0,eps,a,gpp,gps,gsp,gss,spol
      real :: pi4,el4,gam2,arg
      real :: xpp,xps,xsp,xss_psi,xss_zeta
      real :: EXPSATO
      pi4=4.*3.141592754
      el4=el**4
      gam2=gam0**2

      call XSATO(psi,zeta,nu,gam0, &
                       xpp,xps,xsp,xss_psi,xss_zeta)

      arg=(2.*el/gam0)*sin(psi/2.)
      gpp = (el4/pi4)*xpp**2*EXPSATO(eps,a,arg)
      if (gpp.lt.1.e-30) gpp=0.

      arg=(el/gam0)*sqrt(1.+gam2-2.*gam0*cos(psi))
      gps = (1./gam0)*(el4/pi4) * xps**2 * EXPSATO(eps,a,arg)
      if (gps.lt.1.e-30) gps=0.

      gsp = gam0*(el4/pi4) * xsp**2 * EXPSATO(eps,a,arg)
      if (gsp.lt.1.e-30) gsp=0.

      arg=2.*el*sin(psi/2.)
      gss = (el4/pi4)*(xss_psi**2+xss_zeta**2)*EXPSATO(eps,a,arg)
      if (gss.lt.1.e-30) gss=0.
      spol=atan2(xss_zeta,xss_psi)

      return
      end



! XSATO computes (4.50) from Sato and Fehler
!   Inputs:  psi  =  spherical coor. angle (radians)
!            zeta =  sph. coor angle from x3 (radians)
!            nu   =  density vs. velocity pert. scaling (see 4.48)
!            gam0 =  Pvel/Svel (often assumed = sqrt(3) )
!   Returns: xpp, xps, xsp, xss_psi, xss_zeta  =  from eqn. (4.50)
!
      subroutine XSATO(psi,zeta,nu,gam0, &
                       xpp,xps,xsp,xss_psi,xss_zeta)
      implicit none
      real :: psi,zeta,nu,gam0,xpp,xps,xsp,xss_psi,xss_zeta
      real :: gam2,cpsi,c2psi,spsi,czeta,szeta

      gam2=gam0**2

      cpsi=cos(psi)
      c2psi=cos(2.*psi)
      spsi=sin(psi)
      czeta=cos(zeta)
      szeta=sin(zeta)

      xpp = (1./gam2) * (nu * (-1. + cpsi   &
           + (2./gam2)*spsi**2) -2. + (4./gam2)*spsi**2)

      xps = -spsi*(nu*(1.-(2./gam0)*cpsi)-(4./gam0)*cpsi)

      xsp = (1./gam2) * spsi*czeta*(nu*(1.-(2./gam0)*cpsi)  &
                  -(4./gam0)*cpsi)

      xss_psi = czeta*(nu*(cpsi-c2psi)-2.*c2psi)

      xss_zeta = szeta*(nu*(cpsi-1.)+2.*cpsi)

      return
      end

! function EXPSATO computes (2.10) from Sato and Fehler
!    Inputs:  eps  =  RMS velocity perturbation
!             a    =  correlation distance
!             m    =  wavenumber
!    Returns: P(m) =  PSDF (Power Spectral Density Function)
!
      real function EXPSATO(eps,a,m)
      implicit none
      real :: eps,a,m
      real :: kappa,gamma_k1,gamma_k2,expon
      real, parameter :: pi=3.141592654
      real, parameter :: e =2.718281828
          kappa    = 0.5 !Hardwired to be equivalent to exponential ACF 
          gamma_k1 = GAMMA(kappa)
          gamma_k2 = GAMMA(kappa+0.5)
          expon = kappa + 1.5
      expsato=(8.*pi**1.5 * eps**2 * a**3 * gamma_k2) / (gamma_k1*(1.+a**2*m**2)**expon)
      return
      end


!      subroutine ISOSCAT(theta,phi)
!      implicit none
!      real pi,phi,theta
!      pi=3.141592754
!      phi=RAN1(idum)*2.*pi           !returns radians
!      theta=acos(RAN1(idum)*2.-1.)
!      return
!      end
        

      FUNCTION EXPDEV(IDUM)
      implicit none
      integer :: idum
      real    :: expdev,ran1
      EXPDEV=-LOG(RAN1(IDUM))
      RETURN
      END



! BENDRAY updates the ray and polarization direction vectors for
! the scattered ray.  The a1,a2,a3 vectors are overwritten with
! their new values.
! Inputs:
!     a1 = x1 axis in Sato and Fehler (S polarization direction)
!     a2 = x2 axis in S&F
!     a3 = x3 axis in S&F (ray direction)
!     psi = scattering angle from a3 (radians)
!     zeta = scattering angle from a1 (radians)
!     spol = scattered S polarization (radians) 
!            (=0 for pure psi polarization)
!            Note that this polarization is w.r.t the plane defined
!            by the incident and scattered ray vectors.  It is not
!            w.r.t the initial S polarization direction.
!
      subroutine BENDRAY(a1,a2,a3,psi,zeta,spol)
      implicit none
      real    :: a1(3),a2(3),a3(3),b1(3),b2(3),b3(3),c1(3),c2(3),c3(3)
      real    :: psi,zeta,spol,s_psi,s_zeta
      integer :: i
      do i=1,3
         b1(i)=0.
         b2(i)=0.
         b3(i)=0.
         c1(i)=0.
         c2(i)=0.
         c3(i)=0.
      enddo
      s_psi=0.
      s_zeta=0.
      i=0

! define b vectors to be cartesian coor. of e_r, e_psi, e_zeta (S&F Fig. 4.1)
      b3(3)=cos(psi)             !new ray direction
      b3(1)=sin(psi)*cos(zeta)
      b3(2)=sin(psi)*sin(zeta)
      b1(3)=-sin(psi)            !sets psi direction to b1
      b1(1)= cos(psi)*cos(zeta)
      b1(2)= cos(psi)*sin(zeta)
      call CROSS(b3,b1,b2)       !sets zeta direction to b2

! define c vectors in terms of S polarization coordinates
      s_psi=cos(spol)
      s_zeta=sin(spol)
      c1(1) = b1(1)*s_psi + b2(1)*s_zeta
      c1(2) = b1(2)*s_psi + b2(2)*s_zeta
      c1(3) = b1(3)*s_psi + b2(3)*s_zeta
      call VECVEC(b3,c3)
      call CROSS(c3,c1,c2)

! now set b vectors to new ray based vectors using a and c vectors
      do i=1,3
         b1(i) = c1(1)*a1(i) + c1(2)*a2(i) + c1(3)*a3(i)
         b2(i) = c2(1)*a1(i) + c2(2)*a2(i) + c2(3)*a3(i)
         b3(i) = c3(1)*a1(i) + c3(2)*a2(i) + c3(3)*a3(i)
      enddo

! now set a vectors to b vectors
      call VECVEC(b1,a1)
      call VECVEC(b2,a2)
      call VECVEC(b3,a3)

      return
      end



      subroutine VDOT(v1,v2,dot)
      implicit none
      real :: v1(3),v2(3),dot
      dot=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
      return
      end

      subroutine CROSS(v1,v2,v3)
      implicit none
      real :: v1(3),v2(3),v3(3)
      v3(1)=v1(2)*v2(3)-v1(3)*v2(2)   !set v3 to v1 cross v2
      v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
      v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
      return
      end 
     
      subroutine VECVEC(v1,v2)
      implicit none
      real :: v1(3),v2(3)
      v2(1)=v1(1)          !set v2 to v1
      v2(2)=v1(2)
      v2(3)=v1(3)
      return
      end 



      FUNCTION RAN1(IDUM)
      implicit none
      save
      integer :: idum,iff,ix1,ix2,ix3,j
      real    :: r(97),ran1
      integer, PARAMETER :: M1=259200,IA1=7141,IC1=54773
      integer, PARAMETER :: M2=134456,IA2=8121,IC2=28411
      integer, PARAMETER :: M3=243000,IA3=4561,IC3=51349
	  real, parameter    :: RM1=3.8580247E-6, RM2=7.4373773E-6
      DATA IFF /0/
      IF (IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF=1
        IX1=MOD(IC1-IDUM,M1)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX2=MOD(IX1,M2)
        IX1=MOD(IA1*IX1+IC1,M1)
        IX3=MOD(IX1,M3)
        DO 11 J=1,97
          IX1=MOD(IA1*IX1+IC1,M1)
          IX2=MOD(IA2*IX2+IC2,M2)
          R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
11      CONTINUE
        IDUM=1
      ENDIF
      IX1=MOD(IA1*IX1+IC1,M1)
      IX2=MOD(IA2*IX2+IC2,M2)
      IX3=MOD(IA3*IX3+IC3,M3)
      J=1+(97*IX3)/M3
!      IF(J.GT.97.OR.J.LT.1)PAUSE
      IF(J.GT.97.OR.J.LT.1) THEN
	    write(*,'(" 1 ")') ! this line and next to replace
	    read(*,*)          ! above PAUSE statement
	  ENDIF
      RAN1=R(J)
      R(J)=(FLOAT(IX1)+FLOAT(IX2)*RM2)*RM1
      RETURN
      END
