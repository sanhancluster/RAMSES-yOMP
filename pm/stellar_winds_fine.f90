!################################################################
!################################################################
!################################################################
!################################################################
module stellar_commons
   use amr_commons, ONLY:dp
   integer(kind=4):: nt_SW, nz_SW  ! number of grids for time and metallicities
   real(dp),allocatable,dimension(:):: log_tSW ! log10 yr
   real(dp),allocatable,dimension(:):: log_zSW ! log10 z
   ! Notice that the values below should be normalised to 1Msun
   real(dp),allocatable,dimension(:,:):: log_cmSW  ! cumulative mass fraction 
   real(dp),allocatable,dimension(:,:):: log_ceSW  ! cumulative energy per 1Msun SSP     
   real(dp),allocatable,dimension(:,:):: log_cmzSW  ! cumulative metal mass fraction
   real(dp),allocatable,dimension(:,:):: log_cmdSW ! cumulative dust mass fraction
   real(dp),allocatable:: log_cmSW_spec(:,:,:)      ! cumulative mass fraction for several species
   real(dp),allocatable:: log_cmdSW_spec(:,:,:)     ! cumulative dust mass fraction for several species
end module
!################################################################
!################################################################
!################################################################
!################################################################
subroutine init_stellar_winds
   use stellar_commons
   use amr_commons
   implicit none
   integer :: iz, ich, i
   real(dp),allocatable,dimension(:,:):: log_cmHSW  ! cumulative H mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmCSW  ! cumulative C mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmNSW  ! cumulative N mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmOSW  ! cumulative O mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmMgSW ! cumulative Mg mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmSiSW ! cumulative Si mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmSSW  ! cumulative S mass fraction 
   real(dp),allocatable,dimension(:,:):: log_cmFeSW ! cumulative Fe mass fraction
   real(dp),allocatable,dimension(:,:):: log_cmdCSW  ! cumulative C  dust mass fraction
   real(dp),allocatable,dimension(:,:):: log_cmdSiSW ! cumulative Si dust mass fraction
   real(kind=8),allocatable:: dum1d(:)
   logical::ok
   if(.not.use_initial_mass)then
      write(*,*) ' ERROR occured in stellar_winds: use_initial_mass=.false.'
      write(*,*) ' set use_initial_mass=.true. or sn2_real_delay=.true.'
      call clean_stop
   endif

   ! Read stellar winds table
   inquire(FILE=TRIM(stellar_winds_file),exist=ok)
   if(.not.ok)then
      if(myid.eq.1)then
         write(*,*)'Cannot access the file ', TRIM(stellar_winds_file)
      endif
      call clean_stop
   endif

   open(unit=10,file=TRIM(stellar_winds_file),status='old',form='unformatted')
   read(10) nt_SW, nz_SW

   allocate(log_tSW  (1:nt_SW))          ! log Gyr
   allocate(log_zSW  (1:nz_SW))          ! log absolute Z
   allocate(log_cmSW (1:nt_SW,1:nz_SW))  ! log cumulative mass fraction per Msun
   allocate(log_ceSW (1:nt_SW,1:nz_SW))  ! log cumulative energy in erg per Msun
   allocate(log_cmzSW(1:nt_SW,1:nz_SW))  ! log cumulative mass fraction per Msun

   if(nchem>0)then
      allocate(log_cmHSW(1:nt_SW,1:nz_SW))
      allocate(log_cmCSW(1:nt_SW,1:nz_SW))
      allocate(log_cmNSW(1:nt_SW,1:nz_SW))
      allocate(log_cmOSW(1:nt_SW,1:nz_SW))
      allocate(log_cmMgSW(1:nt_SW,1:nz_SW))
      allocate(log_cmSiSW(1:nt_SW,1:nz_SW))
      allocate(log_cmSSW(1:nt_SW,1:nz_SW))
      allocate(log_cmFeSW(1:nt_SW,1:nz_SW))
     
      allocate(log_cmSW_spec(1:nchem,1:nt_SW,1:nz_SW))
   endif
   if(dust)then
      allocate(log_cmdSW  (1:nt_SW,1:nz_SW))
      allocate(log_cmdCSW (1:nt_SW,1:nz_SW))
      allocate(log_cmdSiSW(1:nt_SW,1:nz_SW))
      allocate(log_cmdSW_spec(1:2,1:nt_SW,1:nz_SW))
   endif

   allocate(dum1d (1:nt_SW))
   read(10) dum1d
   log_tSW(:) = dum1d(:)
   deallocate(dum1d)
   allocate(dum1d (1:nz_SW))
   read(10) dum1d
   log_zSW(:) = dum1d(:)
   deallocate(dum1d)

   allocate(dum1d (1:nt_SW))
   !  cumulative stellar mass loss
   do iz=1,nz_SW
      read(10) dum1d
      log_cmSW(:,iz) = dum1d(:)
   enddo
   ! cumulative mechanical energy from winds
   do iz=1,nz_SW
      read(10) dum1d
      log_ceSW(:,iz) = dum1d(:)
   enddo
   ! cumulative metal mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      log_cmzSW(:,iz) = dum1d(:)
   enddo

   ! cumulative H mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmHSW(:,iz) = dum1d(:)
   enddo
   ! cumulative He mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      !log_cmHeSW(:,iz) = dum1d(:)
   enddo
   ! cumulative C mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmCSW(:,iz) = dum1d(:)
   enddo
   ! cumulative N mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmNSW(:,iz) = dum1d(:)
   enddo
   ! cumulative O mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmOSW(:,iz) = dum1d(:)
   enddo
   ! cumulative Mg mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmMgSW(:,iz) = dum1d(:)
   enddo
   ! cumulative Si mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmSiSW(:,iz) = dum1d(:)
   enddo
   ! cumulative S mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmSSW(:,iz) = dum1d(:)
   enddo
   ! cumulative Fe mass from winds
   do iz=1,nz_SW
      read(10) dum1d
      if(nchem>0)log_cmFeSW(:,iz) = dum1d(:)
   enddo

   if(dust)then
      ! cumulative Dust mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmdSW(:,iz) = dum1d(:)
      enddo
      ! cumulative C Dust mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmdCSW(:,iz) = dum1d(:)
      enddo
      ! cumulative Si Dust mass from winds
      do iz=1,nz_SW
         read(10) dum1d
         log_cmdSiSW(:,iz) = dum1d(:)
      enddo
   endif


   if(nchem>0)then
      ich=0
      do i=1,nchem
          ich=ich+1
          if(myid==1) print *, 'SW CHEM(',int(i,kind=2),') = '//TRIM(chem_list(i))
          if(TRIM(chem_list(i))=='H')  log_cmSW_spec(ich,:,:)=log_cmHSW 
          if(TRIM(chem_list(i))=='C')  log_cmSW_spec(ich,:,:)=log_cmCSW 
          if(TRIM(chem_list(i))=='N')  log_cmSW_spec(ich,:,:)=log_cmNSW 
          if(TRIM(chem_list(i))=='O')  log_cmSW_spec(ich,:,:)=log_cmOSW 
          if(TRIM(chem_list(i))=='Mg') log_cmSW_spec(ich,:,:)=log_cmMgSW 
          if(TRIM(chem_list(i))=='Si') log_cmSW_spec(ich,:,:)=log_cmSiSW 
          if(TRIM(chem_list(i))=='S')  log_cmSW_spec(ich,:,:)=log_cmSSW 
          if(TRIM(chem_list(i))=='Fe') log_cmSW_spec(ich,:,:)=log_cmFeSW 
          if(TRIM(chem_list(i))=='D') log_cmSW_spec(ich,:,:)=0d0
      end do
      
      deallocate(log_cmHSW,log_cmCSW,log_cmNSW,log_cmOSW)
      deallocate(log_cmMgSW,log_cmSiSW,log_cmSSW,log_cmFeSW)
   endif
   if(nchem>2.and.dust)then
      log_cmdSW_spec(1,:,:)=log_cmdCSW
      log_cmdSW_spec(2,:,:)=log_cmdSiSW
   endif
   if(dust)then
      deallocate(log_cmdCSW,log_cmdSiSW)
   endif

   deallocate(dum1d)

end subroutine init_stellar_winds
!################################################################
!################################################################
!################################################################
!################################################################
subroutine stellar_winds_fine(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
#ifdef _OPENMP
   use omp_lib
#endif
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the thermal energy, the kinetic energy and 
  ! the metal mass dumped in the gas by stars (winds).
  ! This routine is called every fine time step.
  ! NOTICE:
  ! 1) double-check that only stars have tp.ne.0
  ! 2) set use_initial_mass=.true.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,i,ngrid,iskip,ind,ivar,ipvar,ncache
  integer::ig,ip,npart1,npart2,icpu
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part,ind_cell
  logical::ok_star
  real(dp)::zz,dd1,dd2,zzg ! YD Debug WARNING
  real(dp)::rand
  integer,dimension(1:IRandNumSize),save :: ompseed
!$omp threadprivate(ompseed)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  if(ndim.ne.3) return

#ifdef _OPENMP
!$omp parallel
    ! Give slight offsets for each OMP threads
    ompseed=MOD(tracer_seed+omp_get_thread_num()+1,4096)
!$omp end parallel
#else
    ompseed=MOD(tracer_seed+1,4096)
#endif
    call ranf(tracer_seed,rand)

  ! MC Tracer =================================================
  ! Reset tmpp array that contains the probability to be detached from the particle
  if (MC_tracer) then
     ! Loop over grids
!$omp parallel do private(igrid,npart1,ipart) schedule(dynamic,nchunk)
     do jgrid=1,active(ilevel)%ngrid
        igrid=active(ilevel)%igrid(jgrid)
        npart1=numbp(igrid)  ! Number of particles in the grid
        ! Count star particles
        if(npart1>0)then
           ipart=headp(igrid)
           do jpart=1,npart1
              if(is_star(typep(ipart))) tmpp(ipart) = 0d0
              ipart=nextp(ipart)
           end do
        end if
     end do
  end if
  ! End MC Tracer

  ! Update the rest of the passive scalars so that the fractional quantities are not changed
!!$  ipvar=ichem+nchem
  ipvar=ichem+nchem+ndust
  ncache=active(ilevel)%ngrid
!$omp parallel do private(ngrid,ind_grid,iskip,ind_cell)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           do ivar=ipvar,nvar
              unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar)/max(unew(ind_cell(i),1),smallr)
           end do
        end do
     end do
  end do

!$omp parallel private(igrid,ig,ip,npart1,npart2,ipart,next_part,ok_star,ind_grid,ind_part,ind_grid_part) &
!$omp & default(none) shared(active,ilevel,numbp,headp,nextp,typep,nchunk) reduction(+:dM_prod_SW)
  ig=0
  ip=0
!$omp do schedule(dynamic,nchunk)
  do jgrid=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(jgrid)
     npart1=numbp(igrid)  ! Number of particles in the grid
     npart2=0

     ! Count star particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           ! Select star particles
           ok_star = is_star(typep(ipart))
           if(ok_star)then
              npart2=npart2+1
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif

     ! Gather star particles
     if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           ! Select only star particles
           ok_star = is_star(typep(ipart))
           if(ok_star)then
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
           endif
           if(ip==nvector)then
              call stellar_winds_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,dM_prod_SW)
              ip=0
              ig=0
           end if
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
  end do
  if(ip>0)call stellar_winds_dump(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel,dM_prod_SW)
!$omp end parallel


  if (MC_tracer) then
     ! MC Tracer =================================================
     ! Detach tracer particles from stars
!$omp parallel do private(igrid,npart1,ipart,next_part,rand)
     do jgrid=1,active(ilevel)%ngrid
        igrid=active(ilevel)%igrid(jgrid)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ipart=headp(igrid)
           do jpart=1,npart1
              next_part=nextp(ipart)
              if(is_star_tracer(typep(ipart))) then
                 call ranf(ompseed, rand)
                 if (rand < tmpp(partp(ipart))) then
                    typep(ipart)%family = FAM_TRACER_GAS
                    typep(ipart)%tag = typep(ipart)%tag + 1
                    move_flag(ipart) = 1
                 end if
              end if
              ipart = next_part
           end do
        end if
     end do
  end if

  ! Update the rest of the passive scales so that the fractional quantities are not changed
!$omp parallel do private(ngrid,ind_grid,iskip,ind_cell)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        do i=1,ngrid
           do ivar=ipvar,nvar
              unew(ind_cell(i),ivar) = unew(ind_cell(i),ivar)*unew(ind_cell(i),1)
           end do
        end do
     end do
  end do

111 format('   Entering stlelar winds for level ',I2)

end subroutine stellar_winds_fine
!################################################################
!################################################################
!################################################################
!################################################################
subroutine stellar_winds_dump(ind_grid,ind_part,ind_grid_part,ng,np,ilevel,dM_prod_SW_local)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine stellar_winds_fine. Each stellar particle
  ! dumps mass, momentum and energy in the nearest grid cell using array
  ! unew.
  !-----------------------------------------------------------------------
  integer::i,ii,j,idim,nx_loc,ich,ivar
  real(dp)::dx_min,vol_min
  real(dp)::dx,dx_loc,scale,birth_time
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim)::x0
  integer ,dimension(1:nvector)::ind_cell
  integer ,dimension(1:nvector,1:threetondim)::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim)::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector)::ok
  real(dp),dimension(1:nvector)::mloss,mzloss,mdloss,ethermal,ekinetic,dteff
  real(dp),dimension(1:nvector)::vol_loc
  real(dp),dimension(1:nvector,1:ndim)::x
  integer ,dimension(1:nvector,1:ndim)::id,igd,icd
  integer ,dimension(1:nvector)::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::msun2g
  ! stellar library
  real(dp)::mejecta,dfmloss,dfmzloss,dfmdloss,log_deloss_erg
  real(dp)::mstar_ini,mstar_ini_msun,zstar
  real(dp)::unit_e_code
  real(dp)::dfmloss_spec(1:nchem)
  real(dp)::dfmdloss_spec(1:2)
  real(dp),dimension(1:nchem,1:nvector)::mloss_spec
  real(dp),dimension(1:2,1:nvector)::mdloss_spec
  real(dp),dimension(1:nvar)::uadd
  integer::indp_now
  real(dp)::zz,dd1,dd2,zzg ! YD Debug WARNING
  real(dp),dimension(1:ndust)::dM_prod_SW_local

  msun2g=2d33
  ! starting index for passive variables except for imetal and chem

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc(1:nvector)=dx_loc**ndim
  dx_min=(0.5D0**nlevelmax)*scale
  vol_min=dx_min**ndim

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather 27 neighboring father cells (should be present anytime !)
  do i=1,ng
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

  ! Rescale position at level ilevel
  do idim=1,ndim
     do j=1,np
        x(j,idim)=xp(ind_part(j),idim)/scale+skip_loc(idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)-x0(ind_grid_part(j),idim)
     end do
  end do
  do idim=1,ndim
     do j=1,np
        x(j,idim)=x(j,idim)/dx
     end do
  end do

  ! NGP at level ilevel
  do idim=1,ndim
     do j=1,np
        id(j,idim)=int(x(j,idim))
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igd(j,idim)=id(j,idim)/2
     end do
  end do
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
  do j=1,np
     igrid(j)=son(nbors_father_cells(ind_grid_part(j),kg(j)))
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do j=1,np
     ok(j)=ok(j).and.igrid(j)>0
  end do

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        end if
     end do
  end do
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do

  ! Compute parent cell adresses
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     else
        indp(j) = nbors_father_cells(ind_grid_part(j),kg(j))
        vol_loc(j)=vol_loc(j)*2**ndim ! ilevel-1 cell volume
     end if
  end do

  ! Compute individual time steps
  do j=1,np
     dteff(j)=dtnew(levelp(ind_part(j)))
  end do

  if(use_proper_time)then
     do j=1,np
        dteff(j)=dteff(j)*aexp**2
     end do
  endif

  ! Reset ejected mass, metallicity, thermal energy
  do j=1,np
     mloss(j)=0d0
     mzloss(j)=0d0
     mdloss(j)=0d0
     ethermal(j)=0d0
     mloss_spec(:,j)=0d0
     mdloss_spec(:,j)=0d0
  end do

  ! Compute stellar mass loss and thermal feedback due to stellar winds
  do j=1,np

     ! properties of a star particle
     birth_time=tp(ind_part(j))
     zstar = zp(ind_part(j))
     mstar_ini = mp0(ind_part(j)) ! use_initial_mass=.true. 
     mstar_ini_msun = mstar_ini*(scale_l**3*scale_d/msun2g)

!!$     call cmp_stellar_wind_props (birth_time,dteff(j),zstar,dfmloss,log_deloss_erg,dfmzloss,dfmdloss,dfmloss_spec,dfmdloss_spec)
!!$     call cmp_stellar_wind_props (birth_time,dteff(j),zstar,dfmloss,log_deloss_erg,dfmzloss,dfmdloss,dfmloss_spec,dfmdloss_spec,j)
     call cmp_stellar_wind_props_linearinterpol (birth_time,dteff(j),zstar,dfmloss,log_deloss_erg,dfmzloss,dfmdloss,dfmloss_spec,dfmdloss_spec,j)

     ! Stellar mass loss
     mejecta= mstar_ini*MAX(dfmloss,0.0d0)  ! dfmloss: mass loss fraction during dteff(j)
     mloss(j)=mloss(j)+mejecta/vol_loc(j)
     if(.not. no_wind_energy)then
        ! Thermal energy
        unit_e_code = mstar_ini_msun*(10d0**dble(log_deloss_erg)/dble(msun2g)/scale_v**2)
        ethermal(j)=ethermal(j)+MAX(unit_e_code*(mejecta/vol_loc(j)),0.0d0)
     end if
     ! Metallicity
     if(metal)then
        mzloss(j)=mzloss(j)+mstar_ini*MAX(dfmzloss,0.0d0)/vol_loc(j)
     endif
     ! Dust
     if(dust)then
        mdloss(j)=mdloss(j)+mstar_ini*MAX(dfmdloss,0.0d0)/vol_loc(j)
        ! Dust chemical species
        if(dust_chem)then
           do ich=1,2
              mdloss_spec(ich,j)=mdloss_spec(ich,j)+mstar_ini*MAX(dfmdloss_spec(ich),0.0d0)/vol_loc(j)
           end do
        endif
     endif
     ! Chemical species
     if(nchem>0)then
        do ich=1,nchem
           mloss_spec(ich,j)=mloss_spec(ich,j)+mstar_ini*MAX(dfmloss_spec(ich),0.0d0)/vol_loc(j)
        end do
     endif

     ! Reduce star particle mass
     if(MC_tracer) tmpp(ind_part(j)) = mejecta / mp(ind_part(j))
     mp(ind_part(j))=mp(ind_part(j))-mejecta
  end do

  ! Update hydro variables due to feedback
  uadd=0d0
  indp_now=0
  do j=1,np
     if(indp(j)/=indp_now) then
        if(indp_now > 0) then
           do ivar=1,ichem+nchem+ndust-1
!$omp atomic update
              unew(indp_now,ivar)=unew(indp_now,ivar)+uadd(ivar)
           end do
        end if
        uadd=0d0
        indp_now=indp(j)
     end if

     ! Specific kinetic energy of the star
     ekinetic(j)=0.5*(vp(ind_part(j),1)**2 &
          &          +vp(ind_part(j),2)**2 &
          &          +vp(ind_part(j),3)**2)
     uadd(1)=uadd(1)+mloss(j)
     uadd(2)=uadd(2)+mloss(j)*vp(ind_part(j),1)
     uadd(3)=uadd(3)+mloss(j)*vp(ind_part(j),2)
     uadd(4)=uadd(4)+mloss(j)*vp(ind_part(j),3)
     uadd(5)=uadd(5)+mloss(j)*ekinetic(j)+ ethermal(j)

     ! Add metals
     if(metal)then
        uadd(imetal)=uadd(imetal)+mzloss(j)
     endif

     ! Add Dust
     if(dust)then
        if(dust_chem)then
#if NDUST==2
           uadd(idust  )=uadd(idust  )+mdloss_spec(1,j) ! carbon   (one size)
           uadd(idust+1)=uadd(idust+1)+mdloss_spec(2,j) ! silicate (one size)
           dM_prod_SW_local(1)=dM_prod_SW_local(1)+mdloss_spec(1,j)*vol_loc(j)
           dM_prod_SW_local(2)=dM_prod_SW_local(2)+mdloss_spec(2,j)*vol_loc(j)
#endif
#if NDUST==4
           uadd(idust  )=uadd(idust  )+fsmall_ej*mdloss_spec(1,j) ! carbon   (small size)
           uadd(idust+1)=uadd(idust+1)+flarge_ej*mdloss_spec(1,j) ! carbon   (large size)
           uadd(idust+2)=uadd(idust+2)+fsmall_ej*mdloss_spec(2,j) ! silicate (small size)
           uadd(idust+3)=uadd(idust+3)+flarge_ej*mdloss_spec(2,j) ! silicate (large size)
           dM_prod_SW_local(1)=dM_prod_SW_local(1)+fsmall_ej*mdloss_spec(1,j)*vol_loc(j)
           dM_prod_SW_local(2)=dM_prod_SW_local(2)+flarge_ej*mdloss_spec(1,j)*vol_loc(j)
           dM_prod_SW_local(3)=dM_prod_SW_local(3)+fsmall_ej*mdloss_spec(2,j)*vol_loc(j)
           dM_prod_SW_local(4)=dM_prod_SW_local(4)+flarge_ej*mdloss_spec(2,j)*vol_loc(j)
#endif
        else
#if NDUST==1
           uadd(idust  )=uadd(idust  )+mdloss(j) ! one size
           dM_prod_SW_local(1)=dM_prod_SW_local(1)+mdloss(j)*vol_loc(j)
#endif
#if NDUST==2
           uadd(idust  )=uadd(idust  )+fsmall_ej*mdloss(j) ! small size
           uadd(idust+1)=uadd(idust+1)+flarge_ej*mdloss(j) ! large size
           dM_prod_SW_local(1)=dM_prod_SW_local(1)+fsmall_ej*mdloss(j)*vol_loc(j)
           dM_prod_SW_local(2)=dM_prod_SW_local(2)+flarge_ej*mdloss(j)*vol_loc(j)
#endif
!!$        write(*,*)j,mdloss(j),mdloss(j)*vol_loc(j)*scale_d*scale_l**3/2d33
        endif
     endif

     ! Add individual species
     if(nchem>0)then
        do ich=1,nchem
           uadd(ichem+ich-1)=uadd(ichem+ich-1)+mloss_spec(ich,j)
        end do
     endif
  end do
  if(indp_now > 0) then
     do ivar=1,ichem+nchem+ndust-1
!$omp atomic update
        unew(indp_now,ivar)=unew(indp_now,ivar)+uadd(ivar)
     end do
  end if

end subroutine stellar_winds_dump
!################################################################
!################################################################
!################################################################
!################################################################
!!$subroutine cmp_stellar_wind_props (birth_time,dteff, zstar,dfmloss, log_deloss_erg, dfmzloss, dfmdloss, dfmloss_spec, dfmdloss_spec)
subroutine cmp_stellar_wind_props (birth_time,dteff, zstar,dfmloss, log_deloss_erg, dfmzloss, dfmdloss, dfmloss_spec, dfmdloss_spec,j)
   use amr_commons
   use stellar_commons
   implicit none
   real(dp),intent(in)::birth_time, dteff, zstar
   real(dp)::dfmloss, dfmzloss, dfmdloss, log_deloss_erg
   real(dp),dimension(1:nchem)::dfmloss_spec
   real(dp),dimension(1:2)::dfmdloss_spec
   real(dp)::age1, age2, log_age1,log_age2,log_met
   real(dp)::ft1, ft2, fz
   real(dp)::cm1,cm2,cmz1,cmz2,cmd1,cmd2,ce1,ce2,dum1,dum2
   integer:: itg1, itg2, izg, ich
   integer:: j ! YD Debug

   ! initialise
   dfmloss = 0d0
   dfmzloss = 0d0
   dfmloss_spec = 0d0
   log_deloss_erg = -99.

   ! convert the time to physical units
   call getStarAgeGyr(birth_time+dteff, age1)
   call getStarAgeGyr(birth_time      , age2) ! double-checked.

   log_age1    = log10(max(age1*1d9,1.d0))
   log_age2    = log10(max(age2*1d9,1.d0))
   log_met     = log10(max(zstar,z_ave*0.02))
  
   ! search for the time index from stellar winds library
   call binary_search(log_tSW, log_age1, nt_SW, itg1)
   call binary_search(log_tSW, log_age2, nt_SW, itg2)

   ! search for the metallicity index from stellar winds library
   call binary_search(log_zSW, log_met , nz_SW, izg )

   ! find where we are
   ft1 = (log_tSW(itg1+1) - log_age1)/(log_tSW(itg1+1)-log_tSW(itg1))
   ft2 = (log_tSW(itg2+1) - log_age2)/(log_tSW(itg2+1)-log_tSW(itg2))
   fz  = (log_zSW(izg +1) - log_met )/(log_zSW(izg +1)-log_zSW(izg ))

   ! no extrapolation
   if (ft1 < 0.0) ft1 = 0.0 
   if (ft1 > 1.0) ft1 = 1.0 
   if (ft2 < 0.0) ft2 = 0.0 
   if (ft2 > 1.0) ft2 = 1.0 
   if (fz  < 0.0) fz  = 0.0 
   if (fz  > 1.0) fz  = 1.0 

   ! if a star particle is younger than log_tSW(1), no mass loss 
   if(itg2.eq.1.and.ft2>0.999) return

   ! mass loss fraction during [birth_time, birth_time+dteff]
   dum1 = log_cmSW(itg1,izg  )*ft1 + log_cmSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = log_cmSW(itg1,izg+1)*ft1 + log_cmSW(itg1+1,izg+1)*(1d0-ft1)
   cm1  = dum1*fz + dum2*(1d0-fz)
   dum1 = log_cmSW(itg2,izg  )*ft2 + log_cmSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = log_cmSW(itg2,izg+1)*ft2 + log_cmSW(itg2+1,izg+1)*(1d0-ft2)
   cm2  = dum1*fz + dum2*(1d0-fz)
   dfmloss  = 10d0**cm2  - 10d0**cm1

   ! metal mass loss fraction during [birth_time, birth_time+dteff]
   dum1 = log_cmzSW(itg1,izg  )*ft1 + log_cmzSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = log_cmzSW(itg1,izg+1)*ft1 + log_cmzSW(itg1+1,izg+1)*(1d0-ft1)
   cmz1 = dum1*fz + dum2*(1d0-fz)
   dum1 = log_cmzSW(itg2,izg  )*ft2 + log_cmzSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = log_cmzSW(itg2,izg+1)*ft2 + log_cmzSW(itg2+1,izg+1)*(1d0-ft2)
   cmz2 = dum1*fz + dum2*(1d0-fz)
   dfmzloss = 10d0**cmz2 - 10d0**cmz1

   ! energy during [birth_time, birth_time+dteff]
   dum1 = log_ceSW(itg1,izg  )*ft1 + log_ceSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = log_ceSW(itg1,izg+1)*ft1 + log_ceSW(itg1+1,izg+1)*(1d0-ft1)
   ce1  = dum1*fz + dum2*(1d0-fz)
   dum1 = log_ceSW(itg2,izg  )*ft2 + log_ceSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = log_ceSW(itg2,izg+1)*ft2 + log_ceSW(itg2+1,izg+1)*(1d0-ft2)
   ce2  = dum1*fz + dum2*(1d0-fz)
   if(ce1 < ce2) then
      log_deloss_erg = log10(10d0**dble(ce2) - 10d0**dble(ce1) + 1d-50)
   end if

   if(dust)then
      ! Dust mass loss fraction during [birth_time, birth_time+dteff]
      dum1 = log_cmdSW(itg1,izg  )*ft1 + log_cmdSW(itg1+1,izg  )*(1d0-ft1)
      dum2 = log_cmdSW(itg1,izg+1)*ft1 + log_cmdSW(itg1+1,izg+1)*(1d0-ft1)
      cmd1 = dum1*fz + dum2*(1d0-fz)
      dum1 = log_cmdSW(itg2,izg  )*ft2 + log_cmdSW(itg2+1,izg  )*(1d0-ft2)
      dum2 = log_cmdSW(itg2,izg+1)*ft2 + log_cmdSW(itg2+1,izg+1)*(1d0-ft2)
      cmd2 = dum1*fz + dum2*(1d0-fz)
      dfmdloss = 10d0**cmd2 - 10d0**cmd1
      if(dust_chem)then
         ! carbon
         dum1 = log_cmdSW_spec(1,itg1,izg  )*ft1 + log_cmdSW_spec(1,itg1+1,izg  )*(1d0-ft1)
         dum2 = log_cmdSW_spec(1,itg1,izg+1)*ft1 + log_cmdSW_spec(1,itg1+1,izg+1)*(1d0-ft1)
         cmd1 = dum1*fz + dum2*(1d0-fz)
         dum1 = log_cmdSW_spec(1,itg2,izg  )*ft2 + log_cmdSW_spec(1,itg2+1,izg  )*(1d0-ft2)
         dum2 = log_cmdSW_spec(1,itg2,izg+1)*ft2 + log_cmdSW_spec(1,itg2+1,izg+1)*(1d0-ft2)
         cmd2 = dum1*fz + dum2*(1d0-fz)
         dfmdloss_spec(1) = 10d0**cmd2 - 10d0**cmd1
         ! silicate
         dum1 = log_cmdSW_spec(2,itg1,izg  )*ft1 + log_cmdSW_spec(2,itg1+1,izg  )*(1d0-ft1)
         dum2 = log_cmdSW_spec(2,itg1,izg+1)*ft1 + log_cmdSW_spec(2,itg1+1,izg+1)*(1d0-ft1)
         cmd1 = dum1*fz + dum2*(1d0-fz)
         dum1 = log_cmdSW_spec(2,itg2,izg  )*ft2 + log_cmdSW_spec(2,itg2+1,izg  )*(1d0-ft2)
         dum2 = log_cmdSW_spec(2,itg2,izg+1)*ft2 + log_cmdSW_spec(2,itg2+1,izg+1)*(1d0-ft2)
         cmd2 = dum1*fz + dum2*(1d0-fz)
         dfmdloss_spec(2) = 10d0**cmd2 - 10d0**cmd1
!!$         if(myid==71.and.j==90)then
!!$            write(*,'(A,4es15.7,I9)')'age ',age1,age2,zstar,10d0**log_met,izg
!!$            write(*,'(A,5es15.7)')'SWD ',dfmdloss_spec(2)*FeoverSil/SioverSil&
!!$                 &,log10(10d0**cmd1*FeoverSil/SioverSil) &
!!$                 &,log10(10d0**cmd2*FeoverSil/SioverSil) &
!!$                 &,cmd1,cmd2
!!$         endif
      endif
   endif

   do ich=1,nchem
      ! mass loss fraction during [birth_time, birth_time+dteff]
      dum1 = log_cmSW_spec(ich,itg1,izg  )*ft1 + log_cmSW_spec(ich,itg1+1,izg  )*(1d0-ft1)
      dum2 = log_cmSW_spec(ich,itg1,izg+1)*ft1 + log_cmSW_spec(ich,itg1+1,izg+1)*(1d0-ft1)
      cm1  = dum1*fz + dum2*(1d0-fz)
      dum1 = log_cmSW_spec(ich,itg2,izg  )*ft2 + log_cmSW_spec(ich,itg2+1,izg  )*(1d0-ft2)
      dum2 = log_cmSW_spec(ich,itg2,izg+1)*ft2 + log_cmSW_spec(ich,itg2+1,izg+1)*(1d0-ft2)
      cm2  = dum1*fz + dum2*(1d0-fz)
      dfmloss_spec(ich) = 10d0**cm2  - 10d0**cm1
   end do

!777 format(f5.2,1x,f5.2,1x,5(e15.7,1x))
!   write(*,777) log_age1, log_age2, dfmloss, dfmzloss, log_deloss_erg,ce1,ce2

end subroutine cmp_stellar_wind_props
!################################################################
!################################################################
!################################################################
!################################################################
subroutine cmp_stellar_wind_props_linearinterpol (birth_time,dteff, zstar,dfmloss, log_deloss_erg, dfmzloss, dfmdloss, dfmloss_spec, dfmdloss_spec,j)
   use amr_commons
   use stellar_commons
   implicit none
   real(dp),intent(in)::birth_time, dteff, zstar
   real(dp)::dfmloss, dfmzloss, dfmdloss, log_deloss_erg, dc
   real(dp),dimension(1:nchem)::dfmloss_spec
   real(dp),dimension(1:2)::dfmdloss_spec
   real(dp)::age1, age2, log_age1,log_age2,log_met
   real(dp)::ft1, ft2, fz
   real(dp)::cm1,cm2,cmz1,cmz2,cmd1,cmd2,ce1,ce2,dum1,dum2
   integer:: itg1, itg2, izg, ich
   integer:: j ! YD Debug

   ! initialise
   dfmloss = 0d0
   dfmzloss = 0d0
   dfmdloss = 0d0
   dfmloss_spec = 0d0
   dfmdloss_spec = 0d0
   log_deloss_erg = -99.

   ! convert the time to physical units
   call getStarAgeGyr(birth_time+dteff, age1)
   call getStarAgeGyr(birth_time      , age2) ! double-checked.

   log_age1    = log10(max(age1*1d9,1.d0))
   log_age2    = log10(max(age2*1d9,1.d0))
   log_met     = log10(max(zstar,z_ave*0.02))

   ! search for the time index from stellar winds library
   call binary_search(log_tSW, log_age1, nt_SW, itg1)
   call binary_search(log_tSW, log_age2, nt_SW, itg2)

   ! search for the metallicity index from stellar winds library
   call binary_search(log_zSW, log_met , nz_SW, izg )

   ! find where we are
   ft1 = (10d0**log_tSW(itg1+1) - 10d0**log_age1)/(10d0**log_tSW(itg1+1)-10d0**log_tSW(itg1))
   ft2 = (10d0**log_tSW(itg2+1) - 10d0**log_age2)/(10d0**log_tSW(itg2+1)-10d0**log_tSW(itg2))
   fz  = (10d0**log_zSW(izg +1) - 10d0**log_met )/(10d0**log_zSW(izg +1)-10d0**log_zSW(izg ))

   ! no extrapolation
   if (ft1 < 0.0) ft1 = 0.0
   if (ft1 > 1.0) ft1 = 1.0
   if (ft2 < 0.0) ft2 = 0.0
   if (ft2 > 1.0) ft2 = 1.0
   if (fz  < 0.0) fz  = 0.0
   if (fz  > 1.0) fz  = 1.0

   ! if a star particle is younger than log_tSW(1), no mass loss
   if(itg2.eq.1.and.ft2>0.999) return

   ! mass loss fraction during [birth_time, birth_time+dteff]
   dum1 = 10d0**log_cmSW(itg1,izg  )*ft1 + 10d0**log_cmSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = 10d0**log_cmSW(itg1,izg+1)*ft1 + 10d0**log_cmSW(itg1+1,izg+1)*(1d0-ft1)
   cm1  = dum1*fz + dum2*(1d0-fz)
   dum1 = 10d0**log_cmSW(itg2,izg  )*ft2 + 10d0**log_cmSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = 10d0**log_cmSW(itg2,izg+1)*ft2 + 10d0**log_cmSW(itg2+1,izg+1)*(1d0-ft2)
   cm2  = dum1*fz + dum2*(1d0-fz)
   dfmloss  = cm2 - cm1

   ! metal mass loss fraction during [birth_time, birth_time+dteff]
   dum1 = 10d0**log_cmzSW(itg1,izg  )*ft1 + 10d0**log_cmzSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = 10d0**log_cmzSW(itg1,izg+1)*ft1 + 10d0**log_cmzSW(itg1+1,izg+1)*(1d0-ft1)
   cmz1 = dum1*fz + dum2*(1d0-fz)
   dum1 = 10d0**log_cmzSW(itg2,izg  )*ft2 + 10d0**log_cmzSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = 10d0**log_cmzSW(itg2,izg+1)*ft2 + 10d0**log_cmzSW(itg2+1,izg+1)*(1d0-ft2)
   cmz2 = dum1*fz + dum2*(1d0-fz)
   dfmzloss = cmz2 - cmz1

   ! energy during [birth_time, birth_time+dteff]
   dum1 = 10d0**log_ceSW(itg1,izg  )*ft1 + 10d0**log_ceSW(itg1+1,izg  )*(1d0-ft1)
   dum2 = 10d0**log_ceSW(itg1,izg+1)*ft1 + 10d0**log_ceSW(itg1+1,izg+1)*(1d0-ft1)
   ce1  = dum1*fz + dum2*(1d0-fz)
   dum1 = 10d0**log_ceSW(itg2,izg  )*ft2 + 10d0**log_ceSW(itg2+1,izg  )*(1d0-ft2)
   dum2 = 10d0**log_ceSW(itg2,izg+1)*ft2 + 10d0**log_ceSW(itg2+1,izg+1)*(1d0-ft2)
   ce2  = dum1*fz + dum2*(1d0-fz)
   if(ce1 < ce2) then
      log_deloss_erg = log10(dble(ce2) - dble(ce1) + 1d-50)
   end if

   if(dust)then
      ! Dust mass loss fraction during [birth_time, birth_time+dteff]
      dum1 = 10d0**log_cmdSW(itg1,izg  )*ft1 + 10d0**log_cmdSW(itg1+1,izg  )*(1d0-ft1)
      dum2 = 10d0**log_cmdSW(itg1,izg+1)*ft1 + 10d0**log_cmdSW(itg1+1,izg+1)*(1d0-ft1)
      cmd1 = dum1*fz + dum2*(1d0-fz)
      dum1 = 10d0**log_cmdSW(itg2,izg  )*ft2 + 10d0**log_cmdSW(itg2+1,izg  )*(1d0-ft2)
      dum2 = 10d0**log_cmdSW(itg2,izg+1)*ft2 + 10d0**log_cmdSW(itg2+1,izg+1)*(1d0-ft2)
      cmd2 = dum1*fz + dum2*(1d0-fz)
      dfmdloss = cmd2 - cmd1
      if(dust_chem)then
         ! carbon
         dum1 = 10d0**log_cmdSW_spec(1,itg1,izg  )*ft1 + 10d0**log_cmdSW_spec(1,itg1+1,izg  )*(1d0-ft1)
         dum2 = 10d0**log_cmdSW_spec(1,itg1,izg+1)*ft1 + 10d0**log_cmdSW_spec(1,itg1+1,izg+1)*(1d0-ft1)
         cmd1 = dum1*fz + dum2*(1d0-fz)
         dum1 = 10d0**log_cmdSW_spec(1,itg2,izg  )*ft2 + 10d0**log_cmdSW_spec(1,itg2+1,izg  )*(1d0-ft2)
         dum2 = 10d0**log_cmdSW_spec(1,itg2,izg+1)*ft2 + 10d0**log_cmdSW_spec(1,itg2+1,izg+1)*(1d0-ft2)
         cmd2 = dum1*fz + dum2*(1d0-fz)
         dfmdloss_spec(1) = cmd2 - cmd1
         if(cmd1 < 1d-9 .and. cmd2 < 1d-9) dfmdloss_spec(1) = 0d0
         ! silicate
         dum1 = 10d0**log_cmdSW_spec(2,itg1,izg  )*ft1 + 10d0**log_cmdSW_spec(2,itg1+1,izg  )*(1d0-ft1)
         dum2 = 10d0**log_cmdSW_spec(2,itg1,izg+1)*ft1 + 10d0**log_cmdSW_spec(2,itg1+1,izg+1)*(1d0-ft1)
         cmd1 = dum1*fz + dum2*(1d0-fz)
         dum1 = 10d0**log_cmdSW_spec(2,itg2,izg  )*ft2 + 10d0**log_cmdSW_spec(2,itg2+1,izg  )*(1d0-ft2)
         dum2 = 10d0**log_cmdSW_spec(2,itg2,izg+1)*ft2 + 10d0**log_cmdSW_spec(2,itg2+1,izg+1)*(1d0-ft2)
         cmd2 = dum1*fz + dum2*(1d0-fz)
         dfmdloss_spec(2) = cmd2 - cmd1
         if(cmd1 < 1d-9 .and. cmd2 < 1d-9) dfmdloss_spec(2) = 0d0
!!$         if(myid==71.and.j==90)then
!!$            write(*,'(A,4es15.7,I9)')'age ',age1,age2,zstar,10d0**log_met,izg
!!$            write(*,'(A,5es15.7)')'SWD ',dfmdloss_spec(2)*FeoverSil/SioverSil&
!!$                 &,cmd1*FeoverSil/SioverSil &
!!$                 &,cmd2*FeoverSil/SioverSil &
!!$                 &,cmd1,cmd2
!!$         endif
      endif
   endif

   do ich=1,nchem
      ! mass loss fraction during [birth_time, birth_time+dteff]
      dum1 = 10d0**log_cmSW_spec(ich,itg1,izg  )*ft1 + 10d0**log_cmSW_spec(ich,itg1+1,izg  )*(1d0-ft1)
      dum2 = 10d0**log_cmSW_spec(ich,itg1,izg+1)*ft1 + 10d0**log_cmSW_spec(ich,itg1+1,izg+1)*(1d0-ft1)
      cm1  = dum1*fz + dum2*(1d0-fz)
      dum1 = 10d0**log_cmSW_spec(ich,itg2,izg  )*ft2 + 10d0**log_cmSW_spec(ich,itg2+1,izg  )*(1d0-ft2)
      dum2 = 10d0**log_cmSW_spec(ich,itg2,izg+1)*ft2 + 10d0**log_cmSW_spec(ich,itg2+1,izg+1)*(1d0-ft2)
      cm2  = dum1*fz + dum2*(1d0-fz)
      dfmloss_spec(ich) = cm2 - cm1
      if(cm1 < 1d-9 .and. cm2 < 1d-9) dfmloss_spec(ich) = 0d0
!!$      if(myid==71.and.j==90)then
!!$         if(ich==3)write(*,'(A,3es15.7)')'SWM ',dfmloss_spec(ich),cm1,cm2
!!$      endif
   end do

!777 format(f5.2,1x,f5.2,1x,5(e15.7,1x))
!   write(*,777) log_age1, log_age2, dfmloss, dfmzloss, log_deloss_erg,ce1,ce2

end subroutine cmp_stellar_wind_props_linearinterpol
