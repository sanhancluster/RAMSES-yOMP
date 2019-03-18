!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::emass,xmass,ymass,zmass,rr,rx,ry,rz,rcut
  real(dp)::v200,c,spin,fgas,dnfw,dmax,fcut,rcore,rc,rs
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  emass=2.*boxlen*0.5d0**nlevelmax
  xmass=boxlen/2.0
  ymass=boxlen/2.0
  zmass=boxlen/2.0

  v200=gravity_params(1)
  c   =gravity_params(2)
  spin=gravity_params(3)
  fgas=gravity_params(4)
  rcore=gravity_params(5)
  rcut =gravity_params(7)*3.08d21/scale_l  !Input in kpc  
  rc=rcore*3.08d21/scale_l
  rs=1.0d0
  dmax=1.0/emass/(1.0+emass)**2

  do i=1,ncell
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
     rz=x(i,3)-zmass
#if NDIM==2
     rr=sqrt(rx**2+ry**2)
#endif
#if NDIM==3
     rr=sqrt(rx**2+ry**2+rz**2)
#endif
     fcut=exp(-(rr/rcut)**50)
     if(rcore.ne.0.)then
        dnfw=1d0/(rc+rr)/(rs+rr)**2
     else
        dnfw=1d0/rr/(1d0+rr)**2
     endif
     dnfw=MIN(dmax,dnfw)*fcut
     dnfw=MAX(dnfw,1d-15)
     d(i)=(1d0-fgas)*dnfw
  end do

end subroutine rho_ana
