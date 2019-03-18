subroutine condinit(x,u,dx,nn)
  use amr_commons
  use random
  use hydro_parameters
  use pm_commons
  use poisson_parameters, ONLY: gravity_params
  use cooling_module, ONLY: kB
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
  integer::ivar,i,it,nticks
  real(dp),dimension(1:nvector,1:nvar+4),save::q   ! Primitive variables
  real(dp)::xl,xr,yl,yr,zl,zr,xc,yc,zc
  real(dp)::rx,ry,rz,rr,x0,y0,z0,eps
  real(dp)::mgas,dens,press,temp,sigma,dens_outside
  real(dp)::fcut
  real(dp)::T_sphere,R0,rho_R0,T_outside,rho_out
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !============================================================
  ! Ricarda Beckmann, April 2018
  ! The initial conditions in this file are based on the isothermal
  ! sphere in Federrath2010, Section 3.4. 
  !
  ! INPUT PARAMETERS: added via gravity_params, in the order
  ! T_sphere: temperature of the sphere in K
  ! R0: Truncation radius of the sphere in units of 5E16 cm
  ! rho_R0: Density at the truncation radius, in units of 3.83E-18 g/cm**3
  ! rho_out: Density outside the sphere in units of rho_R0
  !
  ! To reproduce Federrath, set:
  ! gravity_params= 10,1,1,1E-2
  !============================================================

  !constants  
  eps=2.*boxlen*0.5**nlevelmax
  x0=boxlen*0.5
  y0=x0
  z0=x0


  !Load the variables from the gravity params
  T_sphere  = gravity_params(1)         !temperature of the sphere in K
  R0        = gravity_params(2)        !truncation radius in 5E16 cm
  rho_R0    = gravity_params(3)        !density at the truncation radius in units of 3.83E-18 g/cm**3
  rho_out   = gravity_params(4)        !density outside, in units of rho_R0

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
     fcut=exp(-(rr/1.05)**20)

     dens=rho_R0*fcut*(1/rr)**2
     dens=MAX(dens,rho_out)

     t_outside=T_sphere*rho_R0/rho_out
     fcut=MAX(fcut,1E-10)  !To avoid floating point exception in division
     temp=T_sphere/fcut
     temp=MIN(temp,t_outside)

     press=(temp/(2.6*scale_T2))*dens

     ! Euler variables
     q(i,1)=dens  !density
     q(i,2)=0.0        !vx
     q(i,3)=0.0        !vy
     q(i,4)=0.0        !vz
     q(i,5)=press  !pressure
#if NENER>0
     q(i,6)=0.0
#endif
  
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
end subroutine condinit
