subroutine move_fine(ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,local_counter,npart1,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
#ifndef WITHOUTMPI
  integer::info
#endif
  real(dp), dimension(1:ndim) :: xtmp, xsink_tmp
  real(dp) :: dx, twodx, scale
  real(dp), dimension(1:ndim) :: skip_loc, xbound
  integer :: nx_loc
  type(part_t) :: part_type

  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set new sink variables to old ones
  if(sink)then
     vsink_new=0d0
     oksink_new=0d0
     sink_stat(:,ilevel,:)=0d0
  endif

  ! Update particles position and velocity
  ig=0
  ip=0

  ! Loop over particles that are not tracers
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     if(npart1>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        local_counter=0
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle  <---- Very important !!!
           next_part=nextp(ipart)
           ! Classical particles
           if(ig==0)then
              ig=1
              ind_grid(ig)=igrid
           end if
           ! Skip tracers (except "classic" tracers)
           if (.not. (MC_tracer .and. is_tracer(typep(ipart)))) then
              local_counter=local_counter+1
              ip=ip+1
              ind_part(ip)=ipart
              ind_grid_part(ip)=ig
              if(ip==nvector)then
                 call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 local_counter=0
                 ip=0
                 ig=0
              end if
           end if

           ! End MC Tracer patch
           ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles

        ! If there was no particle in the grid, remove the grid from the buffer
        if(local_counter==0 .and. ig>0)then
           ig=ig-1
        end if
     end if
     igrid=next(igrid)   ! Go to next grid
  end do

  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)

  !--------------------------------------------------------------------------------
  ! Moving sinks
  !--------------------------------------------------------------------------------
  if(sink)then
     if(nsink>0)then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
        oksink_all=oksink_new
        vsink_all=vsink_new
#endif
     endif
     do isink=1,nsink
        if(oksink_all(isink)==1d0.and.(.not.fix_smbh_position))then
           vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
           xsink(isink,1:ndim)=xsink(isink,1:ndim)+vsink(isink,1:ndim)*dtnew(ilevel)
        endif
     end do
  endif

  !--------------------------------------------------------------------------------
  ! Moving tracers
  !--------------------------------------------------------------------------------
  if (MC_tracer) then  ! Loop over grids for MC tracers
     ig=0
     ip=0
     ind_grid=0
     ind_part=0
     ind_grid_part=0

     igrid=headl(myid,ilevel)
     do jgrid=1,numbl(myid,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1

           ind_grid(ig)=igrid
           ipart=headp(igrid)
           local_counter=0
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle  <---- Very important !!!
              next_part=nextp(ipart)
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if

              ! call debug_part(ipart, '@move_fine')
              part_type = typep(ipart)
              if (is_gas_tracer(part_type) .and. move_flag(ipart) == 0) then
                 local_counter=local_counter+1
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
                 if(ip==nvector)then
                    call move_gas_tracer(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                    local_counter=0
                    ip=0
                    ig=0
                 end if
              else if (is_star_tracer(part_type)) then
                 xp(ipart, :) = xp(partp(ipart), :)
                 vp(ipart, :) = vp(partp(ipart), :)
              else if (is_cloud_tracer(part_type)) then
                 call move_sink_tracer(ipart, ilevel)
              end if

              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles

           ! If there was no particle in the grid, remove the grid from the buffer
           if (local_counter == 0 .and. ig>0) then
              ig=ig-1
           end if
        end if
        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids
     if(ip>0) call move_gas_tracer(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel) ! MC Tracer

  end if

111 format('   Entering move_fine for level ',I2)

contains

  subroutine move_sink_tracer(ipart, ilevel)
    ! Move the sinks
    integer, intent(in) :: ipart, ilevel
    real(dp) :: d2
    integer :: idim

    dx = 0.5**ilevel
    twodx = 2*dx

    xtmp(:) = xp(ipart, :)
    xsink_tmp(:) = xsink(partp(ipart), :)

    d2 = 0
    do idim = 1, ndim
       d2 = d2 + (xtmp(idim) - xsink_tmp(idim))**2
    end do

    ! Closeby, moving directly
    if (d2 < dx**2) then
       xp(ipart, :) = xsink_tmp(:)
    else
       do idim = 1, ndim
          if (xtmp(idim)+dx < xsink_tmp(idim)) then
             xtmp(idim) = xtmp(idim) + dx
          else if (xtmp(idim)-dx > xsink_tmp(idim)) then
             xtmp(idim) = xtmp(idim) - dx
          else
             xtmp(idim) = xsink_tmp(idim)
          end if
       end do
       xp(ipart, :) = xtmp(:)
    end if

    vp(ipart, :) = vsink(partp(ipart), :)

    ! TODO: there is an issue there when a sink gets close to the
    ! boundaries, as the code does not take into account periodicity
  end subroutine move_sink_tracer

end subroutine move_fine
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move_fine_static(ilevel)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! Update particle position and time-centred velocity at level ilevel.
  ! If particle sits entirely in level ilevel, then use fine grid force
  ! for CIC interpolation. Otherwise, use coarse grid (ilevel-1) force.
  !----------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,ig,ip,local_counter,npart1,npart2,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
#ifndef WITHOUTMPI
  integer::info
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set new sink variables to old ones
  if(sink)then
     vsink_new=0d0
     oksink_new=0d0
     sink_stat(:,ilevel,:)=0d0
  endif

  ! Update particles position and velocity
  ig=0
  ip=0
  ! Loop over grids
  igrid=headl(myid,ilevel)
  do jgrid=1,numbl(myid,ilevel)
     npart1=numbp(igrid)  ! Number of particles in the grid
     npart2=0

     ! Count particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           if(star) then
              if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                   & (.not. static_stars .and. is_not_DM(typep(ipart)) )  ) then
                 ! FIXME: there should be a static_sink as well
                 ! FIXME: what about debris?
                 npart2=npart2+1
              endif
           else
              if(.not.static_DM) then
                 npart2=npart2+1
              endif
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif

     ! Gather DM and star particles
     if(npart2>0)then
        ig=ig+1
        ind_grid(ig)=igrid
        ipart=headp(igrid)
        local_counter=0
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           ! Select particles
           if(star) then
              if ( (.not. static_DM .and. is_DM(typep(ipart))) .or. &
                   & (.not. static_stars .and. is_not_DM(typep(ipart)) )  ) then
                 ! Note: is_not_DM only returns stars and clouds, but not tracers
                 ! FIXME: there should be a static_sink as well
                 ! FIXME: what about debris?
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 local_counter=local_counter+1
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
                 if(ip==nvector) then
                    call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                    local_counter=0
                    ip=0
                    ig=0
                 end if
              endif
           else
              if(.not.static_dm) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 local_counter=local_counter+1
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
                 if(ip==nvector) then
                    call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                    local_counter=0
                    ip=0
                    ig=0
                 end if
              endif
           endif
           ipart=next_part  ! Go to next particle
        end do

        ! If there was no particle in the grid, remove the grid from the buffer
        if (local_counter==0 .and. ig > 0)then
           ig=ig-1
        end if

        ! End loop over particles
     end if
     igrid=next(igrid)   ! Go to next grid
  end do
  ! End loop over grids
  if(ip>0)call move1(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)

  if(sink)then
     if(nsink>0)then
#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
        oksink_all=oksink_new
        vsink_all=vsink_new
#endif
     endif
     do isink=1,nsink
        if(oksink_all(isink)==1d0.and.(.not.fix_smbh_position))then
           vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
           xsink(isink,1:ndim)=xsink(isink,1:ndim)+vsink(isink,1:ndim)*dtnew(ilevel)
        endif
     end do
  endif

111 format('   Entering move_fine_static for level ',I2)

end subroutine move_fine_static
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move1(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, ONLY: uold,smallr
  use cooling_module, ONLY: XH=>X, rhoc, mH
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine computes the force on each particle by
  ! inverse CIC and computes new positions for all particles.
  ! If particle sits entirely in fine level, then CIC is performed
  ! at level ilevel. Otherwise, it is performed at level ilevel-1.
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector),save::father_cell
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x,ff,new_xp,new_vp,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::xbound,skip_loc
  integer ::nlevelmax_loc
  real(dp)::dx_min,vol_min,dx_temp,dx_min_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::d0,mstar,nISM,nCOM
  real(dp)::xx

  ! Family
  logical,dimension(1:nvector),save :: classical_tracer

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**3

  ! Finest cell size
  dx_min=(0.5D0**(nlevelmax-nlevelsheld))*scale
  vol_min=dx_min**ndim
  ! Typical ISM mass density from H/cc to code units
  nISM = n_star
  nCOM = del_star*omega_b*rhoc/aexp**3*XH/mH
  nISM = MAX(nCOM,nISM)
  d0   = nISM/scale_nH
  ! Star particle mass
  mstar=MAX(del_star*omega_b*rhoc*XH/mH,n_star)/(scale_nH*aexp**3)*vol_min
  do i=1,nlevelmax
     dx_temp=scale*0.5D0**i
     ! Test is designed so that nlevelmax is activated at aexp \simeq 0.8
     if(d0*(dx_temp/2.0)**ndim.ge.mstar/2d0)nlevelmax_loc=i+1
  enddo
  dx_min_loc=scale*0.5d0**nlevelmax_loc

  ! Lower left corner of 3x3x3 grid-cube
  do idim=1,ndim
     do i=1,ng
        x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
     end do
  end do

  ! Gather neighboring father cells (should be present anytime !)
  do i=1,ng
     father_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(father_cell,nbors_father_cells,nbors_father_grids,&
       & ng,ilevel)

  ! Rescale particle position at level ilevel
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in move'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
              write(*,*)x(j,1:ndim)
           endif
        end do
     end do
     stop
  end if

  ! CIC at level ilevel (dd: right cloud boundary; dg: left cloud boundary)
  do idim=1,ndim
     do j=1,np
        dd(j,idim)=x(j,idim)+0.5D0
        id(j,idim)=int(dd(j,idim))
        dd(j,idim)=dd(j,idim)-id(j,idim)
        dg(j,idim)=1.0D0-dd(j,idim)
        ig(j,idim)=id(j,idim)-1
     end do
  end do

   ! Compute parent grids
  do idim=1,ndim
     do j=1,np
        igg(j,idim)=ig(j,idim)/2
        igd(j,idim)=id(j,idim)/2
     end do
  end do
#if NDIM==1
  do j=1,np
     kg(j,1)=1+igg(j,1)
     kg(j,2)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j,1)=1+igg(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,2)=1+igd(j,1)+3*igg(j,2)+9*igg(j,3)
     kg(j,3)=1+igg(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,4)=1+igd(j,1)+3*igd(j,2)+9*igg(j,3)
     kg(j,5)=1+igg(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,6)=1+igd(j,1)+3*igg(j,2)+9*igd(j,3)
     kg(j,7)=1+igg(j,1)+3*igd(j,2)+9*igd(j,3)
     kg(j,8)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
  do ind=1,twotondim
     do j=1,np
        igrid(j,ind)=son(nbors_father_cells(ind_grid_part(j),kg(j,ind)))
     end do
  end do

  ! Check if particles are entirely in level ilevel
  ok(1:np)=.true.
  do ind=1,twotondim
     do j=1,np
        ok(j)=ok(j).and.igrid(j,ind)>0
     end do
  end do

  ! If not, rescale position at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           x(j,idim)=x(j,idim)/2.0D0
        end if
     end do
  end do
  ! If not, redo CIC at level ilevel-1
  do idim=1,ndim
     do j=1,np
        if(.not.ok(j))then
           dd(j,idim)=x(j,idim)+0.5D0
           id(j,idim)=int(dd(j,idim))
           dd(j,idim)=dd(j,idim)-id(j,idim)
           dg(j,idim)=1.0D0-dd(j,idim)
           ig(j,idim)=id(j,idim)-1
        end if
     end do
  end do

 ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           icg(j,idim)=ig(j,idim)-2*igg(j,idim)
           icd(j,idim)=id(j,idim)-2*igd(j,idim)
        else
           icg(j,idim)=ig(j,idim)
           icd(j,idim)=id(j,idim)
        end if
     end do
  end do
#if NDIM==1
  do j=1,np
     icell(j,1)=1+icg(j,1)
     icell(j,2)=1+icd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
        icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
        icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
        icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
        icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     else
        icell(j,1)=1+icg(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,2)=1+icd(j,1)+3*icg(j,2)+9*icg(j,3)
        icell(j,3)=1+icg(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,4)=1+icd(j,1)+3*icd(j,2)+9*icg(j,3)
        icell(j,5)=1+icg(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,6)=1+icd(j,1)+3*icg(j,2)+9*icd(j,3)
        icell(j,7)=1+icg(j,1)+3*icd(j,2)+9*icd(j,3)
        icell(j,8)=1+icd(j,1)+3*icd(j,2)+9*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adresses
  do ind=1,twotondim
     do j=1,np
        if(ok(j))then
           indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
        else
           indp(j,ind)=nbors_father_cells(ind_grid_part(j),icell(j,ind))
        end if
     end do
  end do

  ! Compute cloud volumes
#if NDIM==1
  do j=1,np
     vol(j,1)=dg(j,1)
     vol(j,2)=dd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)
     vol(j,2)=dd(j,1)*dg(j,2)
     vol(j,3)=dg(j,1)*dd(j,2)
     vol(j,4)=dd(j,1)*dd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     vol(j,1)=dg(j,1)*dg(j,2)*dg(j,3)
     vol(j,2)=dd(j,1)*dg(j,2)*dg(j,3)
     vol(j,3)=dg(j,1)*dd(j,2)*dg(j,3)
     vol(j,4)=dd(j,1)*dd(j,2)*dg(j,3)
     vol(j,5)=dg(j,1)*dg(j,2)*dd(j,3)
     vol(j,6)=dd(j,1)*dg(j,2)*dd(j,3)
     vol(j,7)=dg(j,1)*dd(j,2)*dd(j,3)
     vol(j,8)=dd(j,1)*dd(j,2)*dd(j,3)
  end do
#endif

  ! Boolean flag for classical tracers (velocity advected ones)
  ! Note: this is enough as the subroutine is only called if `MC_tracer` is false
  do j = 1, np
     classical_tracer(j) = is_gas_tracer(typep(ind_part(j)))
  end do

  ! Gather 3-force
  ff(1:np,1:ndim)=0.0D0
  if(tracer.and.hydro)then
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              if (classical_tracer(j)) then
                 ff(j,idim)=ff(j,idim) + &
                      uold(indp(j,ind),idim+1)/max(uold(indp(j,ind),1),smallr)*vol(j,ind)
              end if
           end do
        end do
     end do
  endif
  if(poisson)then
     do ind=1,twotondim
        do idim=1,ndim
           do j=1,np
              if (.not. classical_tracer(j)) then
                 ff(j,idim)=ff(j,idim)+f(indp(j,ind),idim)*vol(j,ind)
              end if
           end do
        end do
#ifdef OUTPUT_PARTICLE_POTENTIAL
        do j=1,np
           ptcl_phi(ind_part(j)) = phi(indp(j,ind))
        end do
#endif
     end do
  endif

  ! Update velocity
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_vp(j,idim)=ff(j,idim)
        end do
     else
        do j=1,np
           if (classical_tracer(j)) then
              new_vp(j,idim)=ff(j,idim)
           else
              new_vp(j,idim)=vp(ind_part(j),idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
           end if
        end do
     endif
  end do

  ! For sink cloud particle, overwrite velocity with sink velocity
  if(sink)then
     do j=1,np
        if( is_cloud(typep(ind_part(j))) ) then
           isink=-idp(ind_part(j))
           do idim=1,ndim
              new_vp(j,idim)=vsink(isink,idim)+ff(j,idim)*0.5D0*dtnew(ilevel)
           end do
        end if
     end do
  end if

  ! Store velocity
  do idim=1,ndim
     do j=1,np
        vp(ind_part(j),idim)=new_vp(j,idim)
     end do
  end do

  ! Update sink particle velocity using closest cloud particle
  if(sink)then
     do j=1,np
        if( is_cloud(typep(ind_part(j)))) then
           isink=-idp(ind_part(j))
           !Central cloud particle has tag 0
           if(is_central_cloud(typep(ind_part(j))))then
              vsink_new(isink,1:ndim)=vp(ind_part(j),1:ndim)
              oksink_new(isink)=1.0
           endif
           sink_stat(isink,ilevel,ndim*2+1)=sink_stat(isink,ilevel,ndim*2+1)+1d0
           do idim=1,ndim
              xx=xp(ind_part(j),idim)+vp(ind_part(j),idim)*dtnew(ilevel)-xsink(isink,idim)
              if(xx>scale*xbound(idim)/2.0)then
                 xx=xx-scale*xbound(idim)
              endif
              if(xx<-scale*xbound(idim)/2.0)then
                 xx=xx+scale*xbound(idim)
              endif
              sink_stat(isink,ilevel,idim     )=sink_stat(isink,ilevel,idim     )+xsink(isink,idim)+xx
              sink_stat(isink,ilevel,idim+ndim)=sink_stat(isink,ilevel,idim+ndim)+vp(ind_part(j),idim)
           enddo
       endif
     end do
  end if

  ! Update position
  do idim=1,ndim
     if(static)then
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)
        end do
     else
        do j=1,np
           new_xp(j,idim)=xp(ind_part(j),idim)+new_vp(j,idim)*dtnew(ilevel)
        end do
     endif
  end do
  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=new_xp(j,idim)
     end do
  end do

end subroutine move1
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine move_gas_tracer(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  use hydro_commons, only: fluxes
  use tracer_utils, only: safe_move, relative_level, get_cells_on_face
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !------------------------------------------------------------
  ! This routine moves the tracer following the fluxes
  ! This routine is called by move_fine.
  !------------------------------------------------------------
  integer::i,j,idim,nx_loc
  real(dp)::dx,scale
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  ! Particle-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::new_xp,old_xp
  real(dp),dimension(1:nvector,1:3),save::x
  real(dp),dimension(1:3)::skip_loc

  ! MC tracer
  real(dp),dimension(1:nvector,1:twondim,1:twotondim),save::flux
  real(dp),dimension(1:nvector,1:twondim),save::flux_part
  real(dp),dimension(1:nvector),save::rand1, rand2, rand3
  real(dp),dimension(1:nvector,1:ndim),save::proba_correction
  real(dp),dimension(1:twotondim/2)::neighborflux
  integer,dimension(1:twotondim/2)::ncell
  real(dp)::proba1, proba2, proba3
  integer::ison,ipart,iskip,dir,ndir,itmp
  integer,dimension(1:nvector, 0:twondim), save :: ind_ngrid
  integer,dimension(1:nvector, 1:twotondim), save :: ind_cell

  integer,dimension(0:twondim), save :: tmp_ind_ngrid2
  integer,dimension(0:twondim), save :: tmp_ind_ncell2
  real(dp), dimension(1:ndim) :: tmp_old_xp, tmp_new_xp

  integer,dimension(1:nvector) :: ind_parent_part, ison_part, ind_ggrid_part
  integer,dimension(1:nvector, 1:twotondim, 0:twondim) :: ind_ncell
  integer,dimension(1:nvector, 0:twondim) :: tmp_ncell
  integer, dimension(1:nvector) :: tmp_ind_cell
  integer, dimension(1:twotondim/2) :: tmp_ind_ncell
  integer, dimension(1:nvector, 1:twotondim, 1:twondim) :: rel_lvl
  integer, dimension(1:nvector) :: new_partp
  real(dp) :: rand, rtmp ! temporary real
  real(dp), dimension(1:nvector) :: outflux, factor
  logical, dimension(1:nvector) :: move

  logical :: ok

  integer :: ix, iy, iz

  ! Mesh spacing in that level
  dx = 0.5D0**ilevel
  nx_loc = (icoarse_max - icoarse_min + 1)
  skip_loc = (/0.0d0, 0.0d0, 0.0d0/)
  if (ndim > 0) skip_loc(1) = dble(icoarse_min)
  if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
  if (ndim > 2) skip_loc(3) = dble(kcoarse_min)
  scale = boxlen / dble(nx_loc)

  !=======================================!
  ! Grid specific code                    !
  !=======================================!
  ! Get neighbor grids
  call getnborgrids(ind_grid, ind_ngrid, ng)

  ! Compute the index of the cells of each grid
  do ison = 1, twotondim
     iskip = ncoarse + (ison - 1)*ngridmax
     do j = 1, ng
        ind_cell(j, ison) = ind_grid(j) + iskip
     end do
  end do

  ! Get each cell's neighbour in all 6 directions
  do ison = 1, twotondim
     do j = 1, ng
        tmp_ind_cell(j) = ind_cell(j, ison)
     end do
     call getnborfather(tmp_ind_cell, tmp_ncell, ng, ilevel)
     do dir = 0, twondim
        do j = 1, ng
           ind_ncell(j, ison, dir) = tmp_ncell(j, dir)
        end do
     end do
  end do

  ! Loop over dimension
  do ison = 1, twotondim
     iskip = ncoarse + (ison-1)*ngridmax

     ! Compute level of cell 'ison' in direction 'dir'
     do dir = 1, twondim
        do j = 1, ng
           tmp_ind_ngrid2(0:twondim) = ind_ngrid(j, 0:twondim)
           tmp_ind_ncell2(0:twondim) = ind_ncell(j, ison, :)
           call relative_level( &
                tmp_ind_ngrid2, tmp_ind_ncell2, dir, &
                rel_lvl(j, ison, dir))
        end do
     end do

     ! Get the flux for the cells in the grid
     do dir = 1, twondim
        do j = 1, ng
           flux(j, dir, ison) = fluxes(ind_cell(j, ison), dir)
        end do
     end do
  end do

  !=======================================!
  ! Particle specific code                !
  !=======================================!
  ! Create global index of grid
  do ipart = 1, np
     ind_ggrid_part(ipart) = ind_grid(ind_grid_part(ipart))
  end do

  ! Generate random numbers for each particle
  do ipart = 1, np
     call ranf(tracer_seed, rand1(ipart))
     call ranf(tracer_seed, rand2(ipart))
     call ranf(tracer_seed, rand3(ipart))
  end do

  ! Movable particles have a flag == 0
  do ipart = 1, np
     move(ipart) = .true.
  end do

  !=======================================!
  ! Force particles to be attached        !
  !=======================================!
  do idim = 1, ndim
     do j = 1, ng
        x0(j, idim) = xg(ind_grid(j), idim)
     end do
  end do

  ! Compute the location of the particle relative to its grid
  do idim = 1, ndim
     do ipart = 1, np
        ! Get the location in the grid in dx unit from center
        x(ipart, idim) = xp(ind_part(ipart), idim) / scale + skip_loc(idim)
        x(ipart, idim) = x(ipart, idim) - x0(ind_grid_part(ipart), idim)
        x(ipart, idim) = x(ipart, idim) / dx
     end do
  end do

  ! Reset ison_part to 1
  do ipart = 1, np
     ison_part(ipart) = 1
  end do

  ! Move particles to cell centers
  do idim = 1, ndim
     do ipart = 1, np
        ! Particle in the center of the grid
        if (x(ipart, idim) == 0.0D0) then
           call ranf(tracer_seed, rand)

           ! Project the particle either to the right or the left
           if (rand < 0.5) then
              ! do nothing
              x(ipart, idim) = -0.5D0
           else
              ! ison_part += 2^idim / 2
              x(ipart, idim) = +0.5D0
              ison_part(ipart) = ison_part(ipart) + 2**(idim-1)
           end if

        else if (x(ipart, idim) < 0.0D0) then
           x(ipart, idim) = -0.5D0
        else if (x(ipart, idim) > 0.0D0) then
           x(ipart, idim) = +0.5D0
           ison_part(ipart) = ison_part(ipart) + 2**(idim-1)
        end if
     end do
  end do

  ! Recompute the index of the parent cell
  do ipart = 1, np
     iskip = ncoarse + (ison_part(ipart) - 1)*ngridmax
     new_partp(ipart) = ind_ggrid_part(ipart) + iskip
     ind_parent_part(ipart) = new_partp(ipart)
  end do

  ! Get the flux for each gas particle
  do ipart = 1, np
     do dir = 1, twondim
        flux_part(ipart, dir) = flux(ind_grid_part(ipart), dir, ison_part(ipart))
     end do
  end do

  do ipart = 1, np
     factor(ipart) = 1.0D0
  end do

  !=======================================!
  ! Move tracers attached to grid         !
  !=======================================!
  ! Compute the outgoing fluxes for each particle in a cell + mass of cell
  outflux(1:np) = 0
  do dir = 1, twondim
     do ipart = 1, np
        rtmp = flux_part(ipart, dir)
        if (rtmp < 0) outflux(ipart) = outflux(ipart) + rtmp
     end do
  end do

  ! Compute correction factor for probability of moving
  proba_correction = 1.0d0

  ! for each direction ...
  do dir = 1, twondim
     ! 'Reverse' the direction
     if (mod(dir, 2) == 1) then ! 1<->2, 3<->4, 5<->6
        ndir = dir + 1
     else
        ndir = dir - 1
     end if
     ! Get the index in grid of cells on the common face of neighbor grid
     call get_cells_on_face(ndir, tmp_ind_ncell(1:twotondim/2))

     ! ... for each particle ...
     do ipart = 1, np

        ! if (xp(ind_part(ipart), 3) >= 0.59375 .and. xp(ind_part(ipart), 3) <= 0.625 .and. &
        !     (xp(ind_part(ipart), 1) <= 2*0.015625)) then
        !    print*, 'foo!', idp(ind_part(ipart)), ilevel, move(ipart)
        ! end if

        proba1 = -outflux(ipart)


        ! Store the relative level of the neighbor cell
        itmp = rel_lvl(ind_grid_part(ipart), ison_part(ipart), dir)

        ! ... decide whether it'll move ...
        if (itmp == -1) then
           ! Correct bug with subcycling
           ok = move(ipart) .and. rand1(ipart) < proba1 * nsubcycle(ilevel)
        else
           ok = move(ipart) .and. rand1(ipart) < proba1
        end if


        ! if (idp(ind_part(ipart)) == 24114) then
        !    print*, 'debugging', ilevel, xp(ind_part(ipart), :), move_flag(ind_part(ipart)), move(ipart), rand1(ipart), proba1
        ! end if


        if (ok) then
           proba2 = flux_part(ipart, dir)/outflux(ipart) * proba_correction(ipart, 1+(dir-1)/2)

           ! ... pick a direction ...
           if (rand2(ipart) < proba2 .and. move(ipart)) then
              ! if (idp(ind_part(ipart)) == 24114) then
              !    print*, 'moving in direction', dir
              ! end if
              ! Tag the particle as moved
              move(ipart) = .false.

              ! === Move to coarser === !
              if (itmp == -1) then
                 factor(ipart) = 2.0D0
                 new_partp(ipart) = ind_ncell(ind_grid_part(ipart), ison_part(ipart), dir)

                 ! === Move to same level === !
              else if (itmp == 0) then
                 factor(ipart) = 1.0D0
                 new_partp(ipart) = ind_ncell(ind_grid_part(ipart), ison_part(ipart), dir)


                 ! === Move to finer === !
              else
                 factor(ipart) = 0.5D0
                 ! Get the fine-to-coarse flux (and invert the sign)
                 do i = 1, twotondim/2
                    ! Compute the neighbor cell index
                    iskip = ncoarse + (tmp_ind_ncell(i)-1)*ngridmax
                    ncell(i) = son(ind_ncell(ind_grid_part(ipart), ison_part(ipart), dir)) &
                         + iskip
                 end do

                 ! Compute the flux from neighbor cell
                 do i = 1, twotondim/2
                    neighborflux(i) = -fluxes(ncell(i), ndir) / twotondim
                 end do

                 ! Recompute the flux in the direction
                 flux_part(ipart, dir) = sum(neighborflux)

                 ! TODO: fix this!
                 if (flux_part(ipart, dir) == 0) then
                    flux_part(ipart, dir) = 1
                    do i = 1, twotondim/2
                       neighborflux(i) = 2d0/twotondim
                    end do
                 end if

                 ! Chose randomly the target cell
                 do i = 1, twotondim/2
                    proba3 = neighborflux(i) / flux_part(ipart, dir)

                    if (rand3(ipart) == 0) cycle
                    if (rand3(ipart) < proba3) then
                       new_partp(ipart) = ncell(i)

                       ! Prevent other directions
                       rand3(ipart) = 0
                    else
                       rand3(ipart) = rand3(ipart) - max(0d0, proba3)
                    end if
                 end do

              end if
           else
              rand2(ipart) = rand2(ipart) - max(0d0, proba2)
           end if
        end if
     end do
  end do

  ! Actually move the particles to their new location
  do ipart = 1, np
     ! Compute new location in grid + father grid position
     ison_part(ipart) = (new_partp(ipart)-ncoarse-1)/ngridmax + 1
     ind_ggrid_part(ipart) = new_partp(ipart) - ncoarse - (ison_part(ipart)-1)*ngridmax

     iz = (ison_part(ipart)-1)/4
     iy = (ison_part(ipart)-1-4*iz)/2
     ix = (ison_part(ipart)-1-4*iz-2*iy)

     x(ipart, 1) = (dble(ix)-0.5d0)
     x(ipart, 2) = (dble(iy)-0.5d0)
     x(ipart, 3) = (dble(iz)-0.5d0)

     do idim = 1, ndim
        new_xp(ipart, idim) = (xg(ind_ggrid_part(ipart), idim) - skip_loc(idim) &
             + x(ipart, idim)*dx*factor(ipart)) * scale
     end do
  end do

  ! Save old positions
  do idim = 1, ndim
     do ipart = 1, np
        old_xp(ipart, idim) = xp(ind_part(ipart), idim)
     end do
  end do

  ! Safely move particles -- taking care of boundaries
  do ipart = 1, np
     tmp_new_xp(:) = new_xp(ipart, :)
     tmp_old_xp(:) = old_xp(ipart, :)
     call safe_move(tmp_new_xp, tmp_old_xp, scale)
     new_xp(ipart, :) = tmp_new_xp(:)
  end do

  ! Velocity field -- mean of previous speed plus current speed
  do idim = 1, ndim
     do ipart = 1, np
        ! Update speed
        vp(ind_part(ipart), idim) = &
             (new_xp(ipart, idim) - old_xp(ipart, idim)) / dtnew(ilevel)
     end do
  end do

  !--------------------------------------------------------------------
  ! Actually move particles
  !--------------------------------------------------------------------
  do idim = 1, ndim
     do ipart = 1, np
        xp(ind_part(ipart), idim) = new_xp(ipart, idim)
     end do
  end do

  ! Store the new parent (here a cell) of the particle
  do ipart = 1, np
     partp(ind_part(ipart)) = new_partp(ipart)
  end do

end subroutine move_gas_tracer
