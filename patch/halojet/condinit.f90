subroutine condinit(x,u,dx,nn)
  use amr_commons
  use random
  use hydro_parameters
  use pm_commons
  use poisson_parameters, ONLY: gravity_params
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar,i,it
  real(dp),dimension(1:nvector,1:nvar+4),save::q   ! Primitive variables
  real(dp)::rx,ry,rz,rr,x0,y0,z0,eps,rcut,rmax,lr,lrmax,lrcut
  real(dp)::pi,v200,spin,c,fgas,dmax,Bnfw,Bz,dxmin,rcore,rc,rs
  real(dp)::fcut,vcirc,vrot,pnfw,dnfw,mnfw,fc,cosphi
  real(dp)::Axll,Axlr,Axrl,Axrr
  real(dp)::Ayll,Aylr,Ayrl,Ayrr
  real(dp)::xl,xr,yl,yr,zl,zr,xc,yc,zc
  real(dp)::xxl,xxr,yyl,yyr,zzl,zzr,xx,yy
  logical, save:: init_nml=.false.
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::RandNum,phi,Rrand,SS,CC,UU,csound,turb
  real(dp)::dens_min
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Add here, if you wish, some user-defined initial conditions

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  dens_min=1d-6  !Density outside rcut

  eps=2.*boxlen*0.5**nlevelmax
  x0=boxlen/2.0
  y0=boxlen/2.0
  z0=boxlen/2.0

  v200=gravity_params(1)
  c   =gravity_params(2)
  spin=gravity_params(3)
  fgas=gravity_params(4)
  rcore=gravity_params(5)
  turb =gravity_params(6)
  rcut =gravity_params(7)*3.08d21/scale_l  !Input in kpc
  dmax=1d0/eps/(1d0+eps)**2
  pi  =acos(-1d0)
  rc=rcore*3.08d21/scale_l
  rs=1.0d0

  vcirc=sqrt(4.*pi*(log(1d0+c)-c/(1d0+c))/c)
  fc=c*(0.5-0.5/(1d0+c)**2-log(1d0+c)/(1d0+c))/(log(1d0+c)-c/(1d0+c))**2
  Bnfw=B_ave*vcirc ! Scaled magnetic field

  dxmin=boxlen*0.5d0**nlevelmax
  rmax=10*boxlen
  do i=1,nn

     xl=x(i,1)-0.5*dx-x0
     xr=x(i,1)+0.5*dx-x0
     xc=x(i,1)-x0
     yl=x(i,2)-0.5*dx-y0
     yr=x(i,2)+0.5*dx-y0
     yc=x(i,2)-y0
     zl=x(i,3)-0.5*dx-z0
     zr=x(i,3)+0.5*dx-z0
     zc=x(i,3)-z0

#if NDIM==2
     rr=sqrt(xc**2+yc**2)
#endif
#if NDIM==3
     rr=sqrt(xc**2+yc**2+zc**2)
#endif
     fcut=exp(-(rr/rcut)**50)
     if(rcore .ne. 0)then
        dnfw=1d0/(rc+rr)/(rs+rr)**2
     else
        dnfw=1d0/rr/(1d0+rr)**2
     endif
     dnfw=MIN(dmax,dnfw)
     dnfw=dnfw*fcut+dens_min    !Any changes to the density need to be incorporated into fffy below!
     mnfw=(log(1d0+rr)-rr/(1d0+rr))/(log(1d0+c)-c/(1d0+c))
     vrot=2d0*sqrt(2d0/fc)*spin*vcirc*c*mnfw/rr*fcut
     cosphi=sqrt(xc**2+yc**2)/rr
     lr=log(rr)
     lrmax=log(rmax)
     pnfw=romberg(lr,lrmax,1d-5)  !integrate from infinity to get the pressure

     ! Euler variables
     q(i,1)=fgas*dnfw
     q(i,2)=-vrot/rr*yc
     q(i,3)=+vrot/rr*xc
     q(i,4)=0.0
     q(i,5)=pnfw
#if NENER>0
     !Radiation pressure, if used
     q(i,6)=pnfw*1d-5
     if(metal)q(i,7)=z_ave*0.02
#else
     if(metal)q(i,6)=z_ave*0.02
#endif

     if(turb.gt.0.0d0)then
        csound=sqrt(gamma*q(i,5)/q(i,1))
        ! Random velocities
        call ranf(localseed,RandNum)
        SS =(RandNum-0.5)*2.
        call ranf(localseed,RandNum)
        phi=(RandNum-0.5)*2.*pi
        call ranf(localseed,RandNum)
        UU =RandNum
        Rrand=UU**(1./3.) * turb * csound
        CC=Rrand*sqrt(1.-SS**2.)
        q(i,2)=q(i,2)+CC*cos(phi)
        q(i,3)=q(i,3)+CC*sin(phi)
        q(i,4)=q(i,4)+Rrand*SS
     endif

  enddo

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do ivar=1,nener
     u(1:nn,5+ivar)=q(1:nn,5+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+q(1:nn,5+ivar)/(gamma_rad(ivar)-1.0d0)
  enddo
#endif
#if NVAR>5+NENER
  ! passive scalars
  do ivar=6+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

contains
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
  function fffy(rint)
    implicit none
    !      Computes the integrand
    real(dp)::fffy
    real(dp)::rint,rrr,M,rho,vphi,rmass
    rrr=exp(rint)

    !Find maximum radius for mass calculation, if rcut is used
    !if(rrr>rcut)then
    !   rmass=rcut
    !else
    !   rmass=rrr
    !endif
    rmass=rrr

    if(rcore .ne. 0)then
       rho=1d0/(rc+rrr)/(rs+rrr)**2
       rho=rho*exp(-(rrr/rcut)**50)+dens_min
       M=4d0*pi/((rc-rs)**2*(rmass+rs)) * &
            &(rmass*(rc-rs)*rs + (rmass+rs)*(-rc**2*log(rc)+rc**2*log(rmass+rc)+rs*(-2d0*rc+rs)*(-log(rs)+log(rmass+rs))))

    else
       rho=1d0/rrr/(1d0+rrr)**2
       rho=MIN(rho,dmax)
       rho=rho*exp(-(rrr/rcut)**50)+dens_min
       M=4d0*pi*(log(1d0+rmass)-rmass/(1d0+rmass))
    endif
    rho=rho*fgas   !Convert to gas density from DM density

    vphi=2d0*sqrt(2d0/fc)*spin*vcirc*c/rrr &
         & *exp(-(rrr/rcut)**50) &
         & *(log(1d0+rrr)-rrr/(1d0+rrr))/(log(1d0+c)-c/(1d0+c))
    vphi=vphi*cosphi

    fffy=M/rrr*rho-vphi**2*rho
    return
  end function fffy
  !cccccccccccccccccccccccccccccccccccccccccccccccccccc
  function romberg(a,b,tol)
    implicit none
    real(dp)::romberg
    !
    !     Romberg returns the integral from a to b of f(x)dx using Romberg
    !     integration. The method converges provided that f(x) is continuous
    !     in (a,b). The function f must be double precision and must be
    !     declared external in the calling routine.
    !     tol indicates the desired relative accuracy in the integral.
    !
    integer::maxiter=16,maxj=5
    real(dp),dimension(100):: g
    real(dp)::a,b,tol,fourj
    real(dp)::h,error,gmax,g0,g1
    integer::nint,i,j,k,jmax

    h=0.5d0*(b-a)
    gmax=h*(fffy(a)+fffy(b))
    g(1)=gmax
    nint=1
    error=1.0d20
    i=0
10  i=i+1
    if(.not.  (i>maxiter.or.(i>5.and.abs(error)<tol)))then
       !	Calculate next trapezoidal rule approximation to integral.

       g0=0.0d0
       do k=1,nint
          g0=g0+fffy(a+(k+k-1)*h)
       end do
       g0=0.5d0*g(1)+h*g0
       h=0.5d0*h
       nint=nint+nint
       jmax=min(i,maxj)
       fourj=1.0d0

       do j=1,jmax
          ! Use Richardson extrapolation.
          fourj=4.0d0*fourj
          g1=g0+(g0-g(j))/(fourj-1.0d0)
          g(j)=g0
          g0=g1
       enddo
       if (abs(g0).gt.tol) then
          error=1.0d0-gmax/g0
       else
          error=gmax
       end if
       gmax=g0
       g(jmax+1)=g0
       go to 10
    end if
    romberg=g0

    if (i>maxiter.and.abs(error)>tol) &
         &    write(*,*) 'Romberg failed to converge; integral, error=', &
         &    romberg,error



    return
  end function romberg

end subroutine condinit


subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)

     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana
