!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: rhoc, mH
  use mpi_mod
  implicit none
  !-------------------------------------------------------------------------------
  ! Description: This subroutine create sink particle in cells where some density
  ! threshold has been crossed. It also removes from the gas the corresponding
  ! particle mass. On exit, all fluid variables in the cell are modified.
  ! This routine is called only once per coarse step by routine amr_step.
  ! Romain Teyssier, October 7th, 2007
  !-------------------------------------------------------------------------------
  ! local constants
  integer::ilevel,ivar

  if(verbose)write(*,*)' Entering create_sink'
  ! Move particles to finer levels
  do ilevel=levelmin,nlevelmax
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! Get the star density value in each cell
  do ilevel=levelmin,nlevelmax
     call get_rho_star(ilevel)
  enddo

  ! Create new sink particles
  ! and attach all the particles to the grid at level 1
  call make_sink(nlevelmax)
  do ilevel=nlevelmax-1,1,-1
     if(ilevel>=levelmin)call make_sink(ilevel)
     call merge_tree_fine(ilevel)
  end do

  ! Remove particle clouds around old sinks
  call kill_entire_cloud(1)

  ! update sink position before merging sinks and creating clouds
  call update_sink_position_velocity

  ! Merge sink using FOF
  call merge_sink(1)

  ! Create new particle clouds
  call create_cloud_from_sink

  ! Scatter particle back onto to the grid
  do ilevel=1,nlevelmax
     call make_tree_fine(ilevel)
     call kill_tree_fine(ilevel)
     call virtual_tree_fine(ilevel)
  end do

  ! Update hydro quantities for split cells
  if(hydro)then
     do ilevel=nlevelmax,levelmin,-1
        call upload_fine(ilevel)
#ifdef SOLVERmhd
        do ivar=1,nvar+3
#else
        do ivar=1,nvar
#endif
           call make_virtual_fine_dp(uold(1,ivar),ilevel)
        end do
        ! Update boundaries
        if(simple_boundary)call make_boundary_hydro(ilevel)
     end do
  end if

  jsink=0d0
  ! Compute Bondi parameters and gather particle
  do ilevel=nlevelmax,levelmin,-1
     if(bondi)call bondi_hoyle(ilevel)
     call merge_tree_fine(ilevel)
  end do

end subroutine create_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine make_sink(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module, ONLY: rhoc, mH, twopi
  use mpi_mod
  implicit none
  integer::ilevel
  !----------------------------------------------------------------------
  ! Description: This subroutine create sink particle in cells where some
  ! density threshold is crossed. It also removes from the gas the
  ! corresponding particle mass. On exit, all fluid variables in the cell
  ! are modified. This is done only in leaf cells.
  ! Romain Teyssier, October 7th, 2007
  !----------------------------------------------------------------------
  ! local constants
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ::ncache,nnew,ngrid,icpu,index_sink,index_sink_tot
  integer ::igrid,ix,iy,iz,ind,i,j,iskip,isink,nx_loc
  integer ::ind_cloud,itype,idim
  integer ::ntot,ntot_all,info,ntot_tmp,izero_myid,ninc,ntot_myid
  integer:: ivar
  logical ::ok_free
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::d,x,y,z,u,v,w,e,temp
  real(dp)::d_jeans,d_thres,dd_sink
  real(dp)::birth_epoch,rmax_sink2,x_half,y_half,z_half,x_box,y_box,z_box
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp)::dx,dx_loc,scale,dxx,dyy,dzz,dr_sink,rmax_sink
  real(dp)::factG,pi,d_star,star_ratio
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector),save::ind_grid_new,ind_cell_new
  integer ,dimension(1:nvector),save::ind_part_cloud,ind_grid_cloud
  logical ,dimension(1:nvector),save::ok,ok_true=.true.
  integer ,dimension(1:ncpu)::ntot_sink_cpu,ntot_sink_all
  real(dp),dimension(:,:),allocatable::x_tmp,x_tmp_all
  real(dp),dimension(:)  ,allocatable::dens_tmp,dens_tmp_all
  integer ,dimension(:)  ,allocatable::flag_tmp,flag_tmp_all,point2flag2
  ! MC tracer
  integer :: ip, itracer, np
  real(dp) :: proba
  real(dp), dimension(1:nvector), save :: proba_part
  real(dp), dimension(1:nvector, 1:ndim), save :: xsink_loc
  integer, dimension(1:nvector), save :: ind_sink, ind_tracer

  if(numbtot(1,ilevel)==0) return
  if(.not. hydro)return
  if(ndim.ne.3)return

  if(verbose)write(*,*)' Entering make_sink for level ',ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  ! Minimum radius to create a new sink from another
  ! (physical kpc -> code units)
  rmax_sink=r_gal*3.08d21/scale_l*aexp

  pi=twopi/2d0
  factG=1d0
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Density threshold for sink particle creation
  dd_sink=n_sink/scale_nH
  d_star=0d0
  if (star)d_star=n_star/scale_nH

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

  x_half=scale*xbound(1)/2.0; y_half=scale*xbound(2)/2.0; z_half=scale*xbound(3)/2.0
  x_box =scale*xbound(1); y_box =scale*xbound(2); z_box =scale*xbound(3)
  rmax_sink2=rmax_sink**2

  ! Birth epoch
  birth_epoch=t

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Set new sink variables to old ones
  msink_new=0d0; xsink_new=0d0; dMsmbh_new=0d0; Esave_new=0d0; vsink_new=0d0; oksink_new=0d0; tsink_new=0d0; idsink_new=0
  bhspin_new=0d0; spinmag_new=0d0; Efeed_new=0d0
  ! Particles dynamical friction (HP)
  if(drag_part) then
     v_DFnew(1:nsink, 1:ndim, 1:2) = 0d0
     mass_DFnew(1:nsink, 1:2) = tiny(0d0)
     mass_lowspeednew(1:nsink, 1:2) = 0d0
     fact_fastnew(1:nsink, 1:2) = 0d0
     n_partnew(1:nsink, 1:2) = 0
  end if
  !/Particles dynamical friction (HP)

#if NDIM==3

  !------------------------------------------------
  ! Convert hydro variables to primitive variables
  !------------------------------------------------
  ncache=active(ilevel)%ngrid
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
           d=max(uold(ind_cell(i),1), smallr)
           u=uold(ind_cell(i),2)/d
           v=uold(ind_cell(i),3)/d
           w=uold(ind_cell(i),4)/d
           e=uold(ind_cell(i),5)/d
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e-0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=u
           uold(ind_cell(i),3)=v
           uold(ind_cell(i),4)=w
           uold(ind_cell(i),5)=e
           ! No AGN formation site by default
           flag2(ind_cell(i))=1
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=max(uold(ind_cell(i),1), smallr)
              w=uold(ind_cell(i),ivar)/d
              uold(ind_cell(i),ivar)=w
           end do
        end do

     end do
  end do

  !------------------------------
  ! Determine AGN formation sites
  !------------------------------
  ! call quenching(ilevel)

  !----------------------------
  ! Compute number of new sinks
  !----------------------------
  ntot=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Density threshold crossed ---> logical array ok(i)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Flag leaf cells
        do i=1,ngrid
           ok(i)=son(ind_cell(i))==0
        end do
        ! Create new sink if the gas density exceed some threshold
        do i=1,ngrid
           d=max(uold(ind_cell(i),1), smallr)
           star_ratio=rho_star(ind_cell(i))/(d+rho_star(ind_cell(i)))

           ! Jeans length related density threshold
           temp=max(uold(ind_cell(i),5)*(gamma-1.0),smallc**2)
           d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
           d_thres=d_jeans

           ! User defined density threshold
           !d_thres=dd_sink

           ! Check if the stellar density is higher than a star/sink formation threshold
           if(star.and.rho_star(ind_cell(i))<max(d_star, dd_sink))ok(i)=.false.
!!$        ! Check if the star ratio is higher than a certain fraction
!!$        if(star.and.star_ratio<star_ratio_floor)ok(i)=.false.
           ! Check if the density is higher than star/sink formation threshold
           if(d    <max(d_star, dd_sink) )ok(i)=.false.
           ! Check if gas is Jeans unstable
           if(d    <d_thres)ok(i)=.false.
           ! Quenching criterion
           !!!!!!!if(flag2(ind_cell(i))==1)ok(i)=.false.

           ! ENFORCE MASS CRITERION (MT)
           ! Check if the mass criterion will be met
           if (force_exact_mseed) then
              if ((d-d_thres/4.0)*dx_loc**3 < Mseed*2d33/scale_m) ok(i)=.false.
           endif
           ! /ENFORCE MASS CRITERION

           if(ok(i))then
              x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
              y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
              z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
              do isink=1,nsink
                 dxx=x-xsink(isink,1)
                 if(dxx> x_half)then
                    dxx=dxx-x_box
                 endif
                 if(dxx<-x_half)then
                    dxx=dxx+x_box
                 endif
                 dr_sink=dxx*dxx
                 dyy=y-xsink(isink,2)
                 if(dyy> y_half)then
                    dyy=dyy-y_box
                 endif
                 if(dyy<-y_half)then
                    dyy=dyy+y_box
                 endif
                 dr_sink=dyy*dyy+dr_sink
                 dzz=z-xsink(isink,3)
                 if(dzz> z_half)then
                    dzz=dzz-z_box
                 endif
                 if(dzz<-z_half)then
                    dzz=dzz+z_box
                 endif
                 dr_sink=dzz*dzz+dr_sink
                 if(dr_sink .le. rmax_sink2)ok(i)=.false.
              enddo
           endif

        end do
        ! Calculate number of new sinks in each cell plus cloud particles
        do i=1,ngrid
           flag2(ind_cell(i))=0
           if(ok(i))then
              ntot=ntot+1
              flag2(ind_cell(i))=1
           endif
        enddo
     end do
  end do

  if(verbose)write(*,*)'Check multiple sink creation',ntot
  !--------------------------------------------------------------------------------------
  !------NEW: This part avoids multiple sink creation at same coarse time step ----------
  !--------------------------------------------------------------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_tmp,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  ntot_tmp=ntot
#endif
  allocate(point2flag2(1:ntot))
  allocate(x_tmp(1:ntot_tmp,1:3),dens_tmp(1:ntot_tmp),flag_tmp(1:ntot_tmp))
  x_tmp=0d0;dens_tmp=0d0;flag_tmp=0
#ifndef WITHOUTMPI
  ntot_sink_cpu=0; ntot_sink_all=0
  ntot_sink_cpu(myid)=ntot
  call MPI_ALLREDUCE(ntot_sink_cpu,ntot_sink_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_sink_cpu(1)=ntot_sink_all(1)
  do icpu=2,ncpu
     ntot_sink_cpu(icpu)=ntot_sink_cpu(icpu-1)+ntot_sink_all(icpu)
  end do
#endif
  if(myid.gt.1)then
     izero_myid=ntot_sink_cpu(myid-1)
  else
     izero_myid=0
  endif
  ninc=0
  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Gather cells with a new sink
        nnew=0
        do i=1,ngrid
           if (flag2(ind_cell(i))>0)then
              ninc=ninc+1
              x_tmp(izero_myid+ninc,1)=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
              x_tmp(izero_myid+ninc,2)=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
              x_tmp(izero_myid+ninc,3)=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale
              dens_tmp(izero_myid+ninc)=uold(ind_cell(i),1)
              flag_tmp(izero_myid+ninc)=1
              point2flag2(ninc)=ind_cell(i) ! This is a local pointer that is not shared with other cpus
           end if
        end do
     enddo
  enddo
#ifndef WITHOUTMPI
  allocate(x_tmp_all(1:ntot_tmp,1:3),dens_tmp_all(1:ntot_tmp),flag_tmp_all(1:ntot_tmp))
  call MPI_ALLREDUCE(x_tmp   ,x_tmp_all   ,ntot_tmp*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dens_tmp,dens_tmp_all,ntot_tmp  ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(flag_tmp,flag_tmp_all,ntot_tmp  ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  x_tmp   =x_tmp_all
  dens_tmp=dens_tmp_all
  flag_tmp=flag_tmp_all
  deallocate(x_tmp_all,dens_tmp_all,flag_tmp_all)
#endif
  do i=1,ntot_tmp
     if(flag_tmp(i).eq.1)then
        x=x_tmp(i,1);y=x_tmp(i,2);z=x_tmp(i,3)
        do j=i+1,ntot_tmp
           if(flag_tmp(j).eq.1)then
              dxx=x-x_tmp(j,1)
              if(dxx> x_half)then
                 dxx=dxx-x_box
              endif
              if(dxx<-x_half)then
                 dxx=dxx+x_box
              endif
              dr_sink=dxx*dxx
              dyy=y-x_tmp(j,2)
              if(dyy> y_half)then
                 dyy=dyy-y_box
              endif
              if(dyy<-y_half)then
                 dyy=dyy+y_box
              endif
              dr_sink=dyy*dyy+dr_sink
              dzz=z-x_tmp(j,3)
              if(dzz> z_half)then
                 dzz=dzz-z_box
              endif
              if(dzz<-z_half)then
                 dzz=dzz+z_box
              endif
              dr_sink=dzz*dzz+dr_sink
              if(dr_sink .le. rmax_sink2)then
                 ! Keep the largest gas density region
                 if(dens_tmp(i).ge.dens_tmp(j))then
                    flag_tmp(j)=0
                 else
                    flag_tmp(i)=0
                 endif
              endif
           endif
        enddo
     endif
  enddo
  ntot_myid=ntot
  ntot=0
  do i=1,ntot_myid
     if(flag_tmp(izero_myid+i).eq.1)then
        ntot=ntot+1
     else
        flag2(point2flag2(i))=0
     endif
  enddo
  deallocate(x_tmp,dens_tmp,flag_tmp,point2flag2)

  !--------------------------------------------------------------------------------------
  !------NEW: This part avoids multiple sink creation at same coarse time step ----------
  !--------------------------------------------------------------------------------------
  if(verbose)write(*,*)'multiple sink creation checked',ntot

  !---------------------------------
  ! Check for free particle memory
  !---------------------------------
  ok_free=(numbp_free-ntot*ncloud_sink)>=0
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(numbp_free,numbp_free_tot,1,MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  numbp_free_tot=numbp_free
#endif
  if(.not. ok_free)then
     write(*,*)'No more free memory for particles'
     write(*,*)ncloud_sink,ntot
     write(*,*)'Increase npartmax'
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
  end if

  !---------------------------------
  ! Compute global sink statistics
  !---------------------------------
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ntot,ntot_all,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  ntot_all=ntot
#endif
#ifndef WITHOUTMPI
  ntot_sink_cpu=0; ntot_sink_all=0
  ntot_sink_cpu(myid)=ntot
  call MPI_ALLREDUCE(ntot_sink_cpu,ntot_sink_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  ntot_sink_cpu(1)=ntot_sink_all(1)
  do icpu=2,ncpu
     ntot_sink_cpu(icpu)=ntot_sink_cpu(icpu-1)+ntot_sink_all(icpu)
  end do
#endif

  nsink=nsink+ntot_all
  nindsink=nindsink+ntot_all
  if(nsink>nsinkmax)then
     write(*,'(" Increase nsinkmax ",I7,I7)')nsink,nsinkmax
#ifndef WITHOUTMPI
    call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
    stop
#endif
 endif
  if(myid==1)then
     if(ntot_all.gt.0)then
        write(*,'(" Level = ",I6," New sink = ",I6," Tot =",I8)')ilevel,ntot_all,nsink
     endif
  end if

  !------------------------------
  ! Create new sink particles
  !------------------------------
  ! Starting identity number
  if(myid==1)then
     index_sink=nsink-ntot_all
     index_sink_tot=nindsink-ntot_all
  else
     index_sink=nsink-ntot_all+ntot_sink_cpu(myid-1)
     index_sink_tot=nindsink-ntot_all+ntot_sink_cpu(myid-1)
  end if

  ! Loop over grids
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do

        ! Gather cells with a new sink
        nnew=0
        do i=1,ngrid
           if (flag2(ind_cell(i))>0)then
              nnew=nnew+1
              ind_grid_new(nnew)=ind_grid(i)
              ind_cell_new(nnew)=ind_cell(i)
           end if
        end do

        ! Create new sink particles
        do i=1,nnew
           index_sink=index_sink+1
           index_sink_tot=index_sink_tot+1

           ! Get gas variables
           d=max(uold(ind_cell_new(i),1), smallr)
           u=uold(ind_cell_new(i),2)
           v=uold(ind_cell_new(i),3)
           w=uold(ind_cell_new(i),4)
           e=uold(ind_cell_new(i),5)

           ! Get gas cell position
           x=(xg(ind_grid_new(i),1)+xc(ind,1)-skip_loc(1))*scale
           y=(xg(ind_grid_new(i),2)+xc(ind,2)-skip_loc(2))*scale
           z=(xg(ind_grid_new(i),3)+xc(ind,3)-skip_loc(3))*scale

           ! Mass of the new sink

           ! Jeans length related density threshold
           temp=max(e*(gamma-1.0),smallc**2)
           d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
           d_thres=d_jeans
           !d_thres=0.25d0*d

           ! User defined density threshold
           !d_thres=dd_sink

           msink_new (index_sink)=min((d-d_thres/4.0)*dx_loc**3,Mseed*2d33/scale_m)
           dMsmbh_new(index_sink)=0d0
           Esave_new (index_sink)=0d0
           Efeed_new (index_sink)=0d0
           oksink_new(index_sink)=1d0
           bhspin_new(index_sink,1:ndim)=0d0
           spinmag_new(index_sink)=0d0
           ! Particles dynamical friction (HP)
           if (drag_part) then
              do itype = 1, 2
                 v_DFnew(index_sink, 1, itype) = u
                 v_DFnew(index_sink, 2, itype) = v
                 v_DFnew(index_sink, 3, itype) = w
                 mass_DFnew(isink, itype) = 0d0
                 mass_lowspeednew(isink, itype) = 0d0
                 fact_fastnew(isink, itype) = 0d0
                 n_partnew(isink, itype) = 0
              end do
           end if
           !/Particles dynamical friction (HP)

           ! Update linked list
           ind_grid_cloud(1)=ind_grid_new(i)
           call remove_free(ind_part_cloud,1)
           call add_list(ind_part_cloud,ind_grid_cloud,ok_true,1)
           ind_cloud=ind_part_cloud(1)

           ! Set new sink particle variables
           tp(ind_cloud)=0d0             ! Birth epoch
           mp(ind_cloud)=msink_new(index_sink) ! Mass
           levelp(ind_cloud)=ilevel      ! Level
           idp(ind_cloud)=-index_sink    ! Identity
           typep(ind_cloud)%family = FAM_CLOUD
           typep(ind_cloud)%tag = 0
           xp(ind_cloud,1)=x
           xp(ind_cloud,2)=y
           xp(ind_cloud,3)=z
           vp(ind_cloud,1)=u
           vp(ind_cloud,2)=v
           vp(ind_cloud,3)=w
           idsink_new(index_sink) =index_sink_tot
           tsink_new(index_sink)  =birth_epoch
           xsink_new(index_sink,1)=x
           xsink_new(index_sink,2)=y
           xsink_new(index_sink,3)=z
           vsink_new(index_sink,1)=u
           vsink_new(index_sink,2)=v
           vsink_new(index_sink,3)=w

           ! Created a sink: update cell density and put some tracers onto the sink
           ! Note that the mass accreted is not the mass of the black
           ! hole as a (small) part is radiated away
           proba = msink_new(index_sink) / (d*dx_loc**3)
           uold(ind_cell_new(i),1)=d-msink_new(index_sink)/dx_loc**3

           ! Loop over tracer particles in cell
           if (MC_tracer) then
              itracer = headp(ind_grid_new(i))
              np = 0
              do ip = 1, numbp(ind_grid_new(i))
                 ! Only move MC (gas) tracer particles pointing to the current cell
                 if (is_gas_tracer(typep(itracer)) .and. partp(itracer) == ind_cell_new(i) .and. &
                      move_flag(itracer) == 0) then
                    np = np + 1
                    ind_tracer(np)   = itracer
                    proba_part(np)   = proba
                    ind_sink(np)     = index_sink
                    xsink_loc(np, :) = xsink(index_sink, :)
                    if (np == nvector) then
                       call tracer2sink(ind_tracer, proba_part, xsink_loc, ind_sink, np, dx_loc)
                       np = 0
                    end if
                 end if
                 itracer = nextp(itracer)
              end do
              if (np > 0) call tracer2sink(ind_tracer, proba_part, xsink_loc, ind_sink, np, dx_loc)
           end if
        end do
        ! End loop over new sink particle cells

     end do
     ! End loop over cells
  end do
  ! End loop over grids

  !---------------------------------------------------------
  ! Convert hydro variables back to conservative variables
  !---------------------------------------------------------
  ncache=active(ilevel)%ngrid
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
           d=max(uold(ind_cell(i),1), smallr)
           u=uold(ind_cell(i),2)
           v=uold(ind_cell(i),3)
           w=uold(ind_cell(i),4)
           e=uold(ind_cell(i),5)
#ifdef SOLVERmhd
           bx1=uold(ind_cell(i),6)
           by1=uold(ind_cell(i),7)
           bz1=uold(ind_cell(i),8)
           bx2=uold(ind_cell(i),nvar+1)
           by2=uold(ind_cell(i),nvar+2)
           bz2=uold(ind_cell(i),nvar+3)
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(ind_cell(i),1)=d
           uold(ind_cell(i),2)=d*u
           uold(ind_cell(i),3)=d*v
           uold(ind_cell(i),4)=d*w
           uold(ind_cell(i),5)=d*e
        end do
        do ivar=imetal,nvar
           do i=1,ngrid
              d=uold(ind_cell(i),1)
              w=uold(ind_cell(i),ivar)
              uold(ind_cell(i),ivar)=d*w
           end do
        end do
     end do
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(msink_new ,msink_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(xsink_new ,xsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vsink_new ,vsink_all ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(tsink_new ,tsink_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(idsink_new,idsink_all,nsinkmax     ,MPI_INTEGER         ,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dMsmbh_new,dMsmbh_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Esave_new ,Esave_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Efeed_new ,Efeed_all ,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(bhspin_new,bhspin_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(spinmag_new,spinmag_all,nsinkmax   ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  ! Particles dynamical friction (HP)
  if (drag_part) then
     call MPI_ALLREDUCE(v_DFnew,v_DFall,nsinkmax*ndim*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mass_DFnew, mass_DFall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(n_partnew, n_partall,nsinkmax*2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mass_lowspeednew, mass_lowspeedall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(fact_fastnew, fact_fastall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  end if
  !/Particles dynamical friction (HP)

#else
  oksink_all=oksink_new
  msink_all =msink_new
  xsink_all =xsink_new
  vsink_all =vsink_new
  tsink_all =tsink_new
  idsink_all=idsink_new
  dMsmbh_all=dMsmbh_new
  Esave_all =Esave_new
  Efeed_all =Efeed_new
  bhspin_all=bhspin_new
  spinmag_all=spinmag_new
  ! Particles dynamical friction (HP)
  if (drag_part) then
     v_DFall(1:nsink, 1:ndim, 1:2) = v_DFnew(1:nsink, 1:ndim, 1:2)
     mass_DFall(1:nsink, 1:2) = mass_DFnew(1:nsink, 1:2)
     n_partall(1:nsink, 1:2) = n_partnew(1:nsink, 1:2)
     mass_lowspeedall(1:nsink, 1:2) = mass_lowspeednew(1:nsink, 1:2)
     fact_fastall(1:nsink, 1:2) = fact_fastnew(1:nsink, 1:2)
  end if
  !/Particles dynamical friction (HP)
#endif
  do isink=1,nsink
     if(oksink_all(isink)==1)then
        tsink(isink) =tsink_all(isink)
        idsink(isink)=idsink_all(isink)
        msink(isink) =msink_all(isink)
        dMsmbh(isink)=dMsmbh_all(isink)
        Esave (isink)=Esave_all (isink)
        Efeed(isink)=Efeed_all(isink)
        xsink(isink,1:ndim)=xsink_all(isink,1:ndim)
        vsink(isink,1:ndim)=vsink_all(isink,1:ndim)
        bhspin(isink,1:ndim)=bhspin_all(isink,1:ndim)
        spinmag(isink)=spinmag_all(isink)
        ! Particles dynamical friction (HP)
        if (drag_part) then
           do itype = 1,2
              mass_DF(isink, ilevel, itype) = mass_DFall(isink, itype)
              n_part(isink, ilevel, itype) = n_partall(isink, itype)
              mass_lowspeed(isink, ilevel, itype) = mass_lowspeedall(isink, itype)
              fact_fast(isink, ilevel, itype) = fact_fastall(isink, itype)
              do idim = 1, ndim
                 v_DF(isink, ilevel, idim, itype) = v_DFall(isink, idim, itype)
              end do
           end do
        end if
        !/Particles dynamical friction (HP)
     endif
  end do

#endif

end subroutine make_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine merge_sink(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: rhoc, mH, twopi
  use mpi_mod
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine merges sink usink the FOF algorithm.
  ! It keeps only the group centre of mass and remove other sinks.
  !------------------------------------------------------------------------
  integer::isink,ii,idim,new_sink
  real(dp)::dx_loc,scale,dx_min,xx,yy,zz,rr
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,nx_loc
  integer::igrp,icomp,gndx,ifirst,ilast,indx,itype
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  integer,dimension(:),allocatable::psink,gsink,psink_inv
  real(dp),dimension(1:3)::xbound,skip_loc

  integer,dimension(:),allocatable::rank_old,idsink_old
  real(dp),dimension(:),allocatable::tsink_old
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::factG
  real(dp)::egrav,ekin,uxcom,uycom,uzcom,v2rel1,v2rel2,dx_min2
  real(dp)::xc,yc,zc,vxc,vyc,vzc,xx1,yy1,zz1,xx2,yy2,zz2,vvx1,vvy1,vvz1,vvx2,vvy2,vvz2,Lx1,Ly1,Lz1,Lx2,Ly2,Lz2
  real(dp)::q,M1,M2,a1,a2,a1a2,a1L,a2L,Lx,Ly,Lz,Lmod,a1mod,a2mod,mu,af,Lmodana
  real(dp)::ax1,ay1,az1,ax2,ay2,az2,x1,y1,z1,x2,y2,z2,vx1,vy1,vz1,vx2,vy2,vz2
  real(dp)::s4=-0.129d0,s5=-0.384d0,t0=-2.686d0,t2=-3.454d0,t3=+2.353d0,pi

  ! MC Tracer
  integer :: isink_new, isink_old

  pi=twopi/2.0d0
  factG=1d0
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  if(nsink==0)return
  if(verbose)write(*,111)ilevel

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp
  dx_min2=dx_min*dx_min

  allocate(psink(1:nsink),gsink(1:nsink))
  if (MC_tracer) then
     allocate(psink_inv(1:nsink))
  end if

  allocate(rank_old(1:nsink),idsink_old(1:nsink),tsink_old(1:nsink))

  !-------------------------------
  ! Merge sinks using FOF
  !-------------------------------
  do isink=1,nsink
     psink(isink)=isink
     gsink(isink)=0
  end do

  igrp=0
  icomp=1
  ifirst=2
  do while(icomp.le.nsink)
     gndx=psink(icomp)
     if(gsink(gndx)==0)then
        igrp=igrp+1
        gsink(gndx)=igrp
     endif
     ilast=nsink
     do while((ilast-ifirst+1)>0)
        indx=psink(ifirst)
        xx=xsink(indx,1)-xsink(gndx,1)
        if(xx>scale*xbound(1)/2.0)then
           xx=xx-scale*xbound(1)
        endif
        if(xx<-scale*xbound(1)/2.0)then
           xx=xx+scale*xbound(1)
        endif
        rr=xx**2
#if NDIM>1
        yy=xsink(indx,2)-xsink(gndx,2)
        if(yy>scale*xbound(2)/2.0)then
           yy=yy-scale*xbound(2)
        endif
        if(yy<-scale*xbound(2)/2.0)then
           yy=yy+scale*xbound(2)
        endif
        rr=yy**2+rr
#endif
#if NDIM>2
        zz=xsink(indx,3)-xsink(gndx,3)
        if(zz>scale*xbound(3)/2.0)then
           zz=zz-scale*xbound(3)
        endif
        if(zz<-scale*xbound(3)/2.0)then
           zz=zz+scale*xbound(3)
        endif
        rr=zz**2+rr
#endif
 ! Uncomment this line if you want a criterion based on escape velocity
        if(rr.le.rmerge**2.*dx_min2)then
           if(vrel_merge)then
              egrav=msink(indx)*msink(gndx)/(rr+tiny(0.d0))*factG
              uxcom=(msink(indx)*vsink(indx,1)+msink(gndx)*vsink(gndx,1))/(msink(indx)+msink(gndx))
              uycom=(msink(indx)*vsink(indx,2)+msink(gndx)*vsink(gndx,2))/(msink(indx)+msink(gndx))
              uzcom=(msink(indx)*vsink(indx,3)+msink(gndx)*vsink(gndx,3))/(msink(indx)+msink(gndx))
              v2rel1=(vsink(indx,1)-uxcom)**2+(vsink(indx,2)-uycom)**2+(vsink(indx,3)-uzcom)**2
              v2rel2=(vsink(gndx,1)-uxcom)**2+(vsink(gndx,2)-uycom)**2+(vsink(gndx,3)-uzcom)**2
              ekin=0.5d0*(msink(indx)*v2rel1+msink(gndx)*v2rel2)
              if(ekin.lt.egrav)then
                 ifirst=ifirst+1
                 gsink(indx)=igrp
              else
                 psink(ifirst)=psink(ilast)
                 psink(ilast)=indx
                 ilast=ilast-1
              endif
           else
              ifirst=ifirst+1
              gsink(indx)=igrp
           endif
        else
           psink(ifirst)=psink(ilast)
           psink(ilast)=indx
           ilast=ilast-1
        endif
     end do
     icomp=icomp+1
  end do
  new_sink=igrp
  if(myid==1)then
     write(*,*)'Found ',new_sink,' groups'
     !do isink=1,nsink
     !   write(*,'(3(I4,1x),3(1PE10.3))')isink,psink(isink),gsink(isink),xsink(isink,1:ndim)
     !end do
  endif

  ! Mapping between the old position and the new ones: psink_inv(iold) = inew
  ! Note: this is the opposite of psink(inew) = iold
  if (MC_tracer) then
     do isink = 1, nsink
        psink_inv(psink(isink)) = isink
     end do
  end if


  !----------------------------------------------------
  ! Compute group centre of mass and average velocity
  !----------------------------------------------------
  xsink_new=0d0; vsink_new=0d0; msink_new=0d0; dMsmbh_new=0d0; Esave_new=0d0; idsink_new=0
  oksink_all=0d0; oksink_new=0d0; tsink_new=0d0
  rank_old=0; idsink_old=0; tsink_old=0
  bhspin_new=0d0; spinmag_new=0d0
  ! Particles dynamical friction (HP)
  if (drag_part) then
     v_DFnew_all = 0d0; mass_DFnew_all = 0d0; n_partnew_all=0
     mass_lowspeednew_all=0d0; fact_fastnew_all=0d0
     most_massive_sink = 0d0
  end if
  !/Particles dynamical friction (HP)
  do isink=1,nsink
     igrp=gsink(isink)

     !----------------------------------------------------
     ! This is done to keep track of the most massive sink
     ! after a merger with a companion
     !----------------------------------------------------
     if ( rank_old(igrp) .eq. 0)then
        rank_old(igrp)=isink
        idsink_old(igrp)=idsink(isink)
        tsink_old(igrp) =tsink(isink)
     endif
     if ( msink(isink) .gt. msink(rank_old(igrp)) )then
        rank_old(igrp)=isink
        idsink_new(igrp)=idsink(isink)
        idsink_old(igrp)=idsink(isink)
        tsink_new(igrp) =tsink(isink)
        tsink_old(igrp) =tsink(isink)
     else
        idsink_new(igrp)=idsink_old(igrp)
        tsink_new(igrp) =tsink_old(igrp)
     endif
     !----------------------------------------------------
     !----------------------------------------------------

     !idsink_new(igrp)=idsink(isink)
     if(oksink_new(igrp)==0d0)then
        oksink_all(isink)=igrp
        oksink_new(igrp)=isink
     endif
     !YDspin----------------------------------------------------
     if(msink_new(igrp).ge.msink(isink))then
        if(msink(isink).eq.0d0)then
           write(*,*)'Problem in merge_sink for spins, msink=0'
           stop
        endif
        M1=msink_new(igrp)
        M2=msink(isink)
        a1=ABS(spinmag_new(igrp))
        a2=ABS(spinmag(isink))
        ax1=bhspin_new(igrp,1)
        ay1=bhspin_new(igrp,2)
        az1=bhspin_new(igrp,3)
        ax2=bhspin(isink,1)
        ay2=bhspin(isink,2)
        az2=bhspin(isink,3)
        xx1=xsink_new(igrp,1)
        yy1=xsink_new(igrp,2)
        zz1=xsink_new(igrp,3)
        xx2=xsink(isink,1)
        yy2=xsink(isink,2)
        zz2=xsink(isink,3)
        vvx1=vsink_new(igrp,1)
        vvy1=vsink_new(igrp,2)
        vvz1=vsink_new(igrp,3)
        vvx2=vsink(isink,1)
        vvy2=vsink(isink,2)
        vvz2=vsink(isink,3)
     else
        if(msink_new(igrp).gt.0d0)then
           M1=msink(isink)
           M2=msink_new(igrp)
           a1=ABS(spinmag(isink))
           a2=ABS(spinmag_new(igrp))
           ax1=bhspin(isink,1)
           ay1=bhspin(isink,2)
           az1=bhspin(isink,3)
           ax2=bhspin_new(igrp,1)
           ay2=bhspin_new(igrp,2)
           az2=bhspin_new(igrp,3)
           xx1=xsink(isink,1)
           yy1=xsink(isink,2)
           zz1=xsink(isink,3)
           xx2=xsink_new(igrp,1)
           yy2=xsink_new(igrp,2)
           zz2=xsink_new(igrp,3)
           vvx1=vsink(isink,1)
           vvy1=vsink(isink,2)
           vvz1=vsink(isink,3)
           vvx2=vsink_new(igrp,1)
           vvy2=vsink_new(igrp,2)
           vvz2=vsink_new(igrp,3)
        else
           M2=0d0
        endif
     endif
     if(M2.ne.0d0.and.bhspinmerge)then
        q=M2/M1
        mu=q/(1d0+q)**2
        a1mod=SQRT(ax1**2+ay1**2+az1**2)
        a2mod=SQRT(ax2**2+ay2**2+az2**2)
        if(a1mod>0.and.a2mod>0)then
           a1a2=(ax1*ax2+ay1*ay2+az1*az2)/(a1mod*a2mod)
        else
           a1a2=0d0
        endif
        ! C.O.M.
        ! Check for periodicity
        xx=xx2-xx1
        if(xx>scale*xbound(1)/2.0)then
           xx2=xx2-scale*xbound(1)
        endif
        if(xx<-scale*xbound(1)/2.0)then
           xx1=xx1+scale*xbound(1)
        endif
        yy=yy2-yy1
        if(yy>scale*xbound(2)/2.0)then
           yy2=yy2-scale*xbound(2)
        endif
        if(yy<-scale*xbound(2)/2.0)then
           yy1=yy1+scale*xbound(2)
        endif
        zz=zz2-zz1
        if(zz>scale*xbound(3)/2.0)then
           zz2=zz2-scale*xbound(3)
        endif
        if(zz<-scale*xbound(3)/2.0)then
           zz1=zz1+scale*xbound(3)
        endif
        xc =(M1*xx1+M2*xx2)/(M1+M2)
        yc =(M1*yy1+M2*yy2)/(M1+M2)
        zc =(M1*zz1+M2*zz2)/(M1+M2)
        vxc=(M1*vvx1+M2*vvx2)/(M1+M2)
        vyc=(M1*vvy1+M2*vvy2)/(M1+M2)
        vzc=(M1*vvz1+M2*vvz2)/(M1+M2)
        x1 =xx1-xc
        y1 =yy1-yc
        z1 =zz1-zc
        x2 =xx2-xc
        y2 =yy2-yc
        z2 =zz2-zc
        vx1=vvx1-vxc
        vy1=vvy1-vyc
        vz1=vvz1-vzc
        vx2=vvx2-vxc
        vy2=vvy2-vyc
        vz2=vvz2-vzc
        Lx1=M1*(y1*vz1-z1*vy1)
        Ly1=M1*(z1*vx1-x1*vz1)
        Lz1=M1*(x1*vy1-y1*vx1)
        Lx2=M2*(y2*vz2-z2*vy2)
        Ly2=M2*(z2*vx2-x2*vz2)
        Lz2=M2*(x2*vy2-y2*vx2)
        Lx=Lx1+Lx2
        Ly=Ly1+Ly2
        Lz=Lz1+Lz2
        Lmod =SQRT(Lx**2+Ly**2+Lz**2)
        if(a1mod>0.and.a2mod>0.and.Lmod>0.)then
           a1L=(ax1*Lx+ay1*Ly+az1*Lz)/(a1mod*Lmod)
           a2L=(ax2*Lx+ay2*Ly+az2*Lz)/(a2mod*Lmod)
        else
           a1L=0d0
           a2L=0d0
        endif
        Lmodana = s4 / ( 1d0 + q**2 )**2 * ( a1**2 + a2**2*q**4 + 2d0*a1*a2*q**2*a1a2 ) &
             & + ( s5*mu + t0 + 2d0 ) / ( 1d0 + q**2 ) * ( a1*a1L + a2*q**2*a2L ) &
             & + 2d0*SQRT(3d0) + t2*mu + t3*mu**2
        af = 1d0 / ( 1d0 + q )**2 * SQRT( a1**2 + a2**2*q**4 + 2d0*a1*a2*q**2*a1a2 &
             & + 2d0 * ( a1*a1L + a2*q**2*a2L ) * Lmodana*q + Lmodana**2*q**2 )
        bhspin_new(igrp,1)=1d0/(1d0+q)**2 * ( a1*ax1+a2*ax2*q**2+(Lx/Lmod)*Lmodana*q )
        bhspin_new(igrp,2)=1d0/(1d0+q)**2 * ( a1*ay1+a2*ay2*q**2+(Ly/Lmod)*Lmodana*q )
        bhspin_new(igrp,3)=1d0/(1d0+q)**2 * ( a1*az1+a2*az2*q**2+(Lz/Lmod)*Lmodana*q )
        spinmag_new(igrp)=SQRT(bhspin_new(igrp,1)**2 + bhspin_new(igrp,2)**2 + bhspin_new(igrp,3)**2 )
        if(spinmag_new(igrp).gt.+maxspin) spinmag_new(igrp)=+maxspin
        if(spinmag_new(igrp).lt.-maxspin) spinmag_new(igrp)=-maxspin
        ax1=bhspin_new(igrp,1)
        ay1=bhspin_new(igrp,2)
        az1=bhspin_new(igrp,3)
        a1mod=SQRT(ax1**2+ay1**2+az1**2)
        bhspin_new(igrp,1)=ax1/a1mod
        bhspin_new(igrp,2)=ay1/a1mod
        bhspin_new(igrp,3)=az1/a1mod
     else
        if ( msink(isink) .gt. msink_new(igrp) )then
           ! Case it's not a merger
           bhspin_new(igrp,1:ndim)=bhspin(isink,1:ndim)
           spinmag_new(igrp)=spinmag(isink)
        endif
     endif
     !YDspin----------------------------------------------------

     ! Particles dynamical friction (HP)
     if (drag_part) then
        if (msink(isink).gt.most_massive_sink(igrp)) then
           most_massive_sink(igrp) = msink(isink)
           do ii = levelmin, nlevelmax
              do itype = 1, 2
                 do idim = 1, ndim
                    v_DFnew_all(igrp, ii, idim, itype) = v_DF(isink, ii, idim, itype)
                 end do
                 mass_DFnew_all(igrp, ii, itype) = mass_DF(isink, ii, itype)
                 n_partnew_all(igrp, ii, itype) = n_part(isink, ii, itype)
                 mass_lowspeednew_all(igrp, ii, itype) = mass_lowspeed(isink, ii, itype)
                 fact_fastnew_all(igrp, ii, itype) = fact_fast(isink, ii, itype)
              end do
           end do
        end if
     end if
     !/Particles dynamical friction (HP)

     msink_new (igrp)=msink_new (igrp)+msink (isink)
     dMsmbh_new(igrp)=dMsmbh_new(igrp)+dMsmbh(isink)
     Esave_new (igrp)=Esave_new (igrp)+Esave (isink)
     xx=xsink(isink,1)-xsink(int(oksink_new(igrp)),1)
     if(xx>scale*xbound(1)/2.0)then
        xx=xx-scale*xbound(1)
     endif
     if(xx<-scale*xbound(1)/2.0)then
        xx=xx+scale*xbound(1)
     endif
     xsink_new(igrp,1)=xsink_new(igrp,1)+msink(isink)*xx
     vsink_new(igrp,1)=vsink_new(igrp,1)+msink(isink)*vsink(isink,1)
#if NDIM>1
     yy=xsink(isink,2)-xsink(int(oksink_new(igrp)),2)
     if(yy>scale*xbound(2)/2.0)then
        yy=yy-scale*xbound(2)
     endif
     if(yy<-scale*xbound(2)/2.0)then
        yy=yy+scale*xbound(2)
     endif
     xsink_new(igrp,2)=xsink_new(igrp,2)+msink(isink)*yy
     vsink_new(igrp,2)=vsink_new(igrp,2)+msink(isink)*vsink(isink,2)
#endif
#if NDIM>2
     zz=xsink(isink,3)-xsink(int(oksink_new(igrp)),3)
     if(zz>scale*xbound(3)/2.0)then
        zz=zz-scale*xbound(3)
     endif
     if(zz<-scale*xbound(3)/2.0)then
        zz=zz+scale*xbound(3)
     endif
     xsink_new(igrp,3)=xsink_new(igrp,3)+msink(isink)*zz
     vsink_new(igrp,3)=vsink_new(igrp,3)+msink(isink)*vsink(isink,3)
#endif
  end do
  do isink=1,new_sink
     xsink_new(isink,1)=xsink_new(isink,1)/msink_new(isink)+xsink(int(oksink_new(isink)),1)
     vsink_new(isink,1)=vsink_new(isink,1)/msink_new(isink)
#if NDIM>1
     xsink_new(isink,2)=xsink_new(isink,2)/msink_new(isink)+xsink(int(oksink_new(isink)),2)
     vsink_new(isink,2)=vsink_new(isink,2)/msink_new(isink)
#endif
#if NDIM>2
     xsink_new(isink,3)=xsink_new(isink,3)/msink_new(isink)+xsink(int(oksink_new(isink)),3)
     vsink_new(isink,3)=vsink_new(isink,3)/msink_new(isink)
#endif
  end do
  nsink=new_sink
  msink (1:nsink)=msink_new (1:nsink)
  dMsmbh(1:nsink)=dMsmbh_new(1:nsink)
  Esave (1:nsink)=Esave_new (1:nsink)
  idsink(1:nsink)=idsink_new(1:nsink)
  tsink (1:nsink)=tsink_new (1:nsink)
  xsink(1:nsink,1:ndim)=xsink_new(1:nsink,1:ndim)
  vsink(1:nsink,1:ndim)=vsink_new(1:nsink,1:ndim)
  bhspin(1:nsink,1:ndim)=bhspin_new(1:nsink,1:ndim)
  spinmag(1:nsink)=spinmag_new(1:nsink)
  ! Particles dynamical friction (HP)
  if (drag_part) then
     v_DF(1:nsink, levelmin:nlevelmax, 1:ndim, 1:2) = v_DFnew_all(1:nsink, levelmin:nlevelmax, 1:ndim, 1:2)
     mass_DF(1:nsink, levelmin:nlevelmax, 1:2) = mass_DFnew_all(1:nsink, levelmin:nlevelmax, 1:2)
     n_part(1:nsink, levelmin:nlevelmax, 1:2) = n_partnew_all(1:nsink, levelmin:nlevelmax, 1:2)
     mass_lowspeed(1:nsink, levelmin:nlevelmax, 1:2) = mass_lowspeednew_all(1:nsink, levelmin:nlevelmax, 1:2)
     fact_fast(1:nsink, levelmin:nlevelmax, 1:2) = fact_fastnew_all(1:nsink, levelmin:nlevelmax, 1:2)
  end if
  !/Particles dynamical friction (HP)

  ! Periodic boundary conditions
  do isink=1,nsink
     xx=xsink(isink,1)
     if(xx<-scale*skip_loc(1))then
        xx=xx+scale*(xbound(1)-skip_loc(1))
     endif
     if(xx>scale*(xbound(1)-skip_loc(1)))then
        xx=xx-scale*(xbound(1)-skip_loc(1))
     endif
     xsink(isink,1)=xx
#if NDIM>1
     yy=xsink(isink,2)
     if(yy<-scale*skip_loc(2))then
        yy=yy+scale*(xbound(2)-skip_loc(2))
     endif
     if(yy>scale*(xbound(2)-skip_loc(2)))then
        yy=yy-scale*(xbound(2)-skip_loc(2))
     endif
     xsink(isink,2)=yy
#endif
#if NDIM>2
     zz=xsink(isink,3)
     if(zz<-scale*skip_loc(3))then
        zz=zz+scale*(xbound(3)-skip_loc(3))
     endif
     if(zz>scale*(xbound(3)-skip_loc(3)))then
        zz=zz-scale*(xbound(3)-skip_loc(3))
     endif
     xsink(isink,3)=zz
#endif
  enddo

  deallocate(rank_old,idsink_old,tsink_old)

  !-----------------------------------------------------
  ! Remove sink particles that are part of a FOF group.
  !-----------------------------------------------------
  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather sink particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call kill_sink(ind_part,ind_grid_part,ip)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call kill_sink(ind_part,ind_grid_part,ip)
  end do
  ! End loop over cpus

  !-----------------------------------------------------
  ! Take care of the tracer particles
  !-----------------------------------------------------
  if (MC_tracer) then
     ! Loop over cpus
     do icpu=1,ncpu
        igrid=headl(icpu,ilevel)
        ig=0
        ip=0
        ! Loop over grids
        do jgrid=1,numbl(icpu,ilevel)
           npart1=numbp(igrid)  ! Number of particles in the grid
           npart2=0

           ! Count sink particles
           if(npart1>0)then
              ipart=headp(igrid)
              ! Loop over particles
              do jpart=1,npart1
                 ! Save next particle   <--- Very important !!!
                 next_part=nextp(ipart)
                 if (is_cloud_tracer(typep(ipart))) then
                    npart2=npart2+1
                 endif
                 ipart=next_part  ! Go to next particle
              end do
           end if

           ! Gather sink particles
           if (npart2 > 0) then
              ipart = headp(igrid)
              ! Loop over particles
              do jpart = 1, npart1
                 ! Save next particle   <--- Very important !!!
                 next_part = nextp(ipart)
                 ! Select only tracer particles
                 if (is_cloud_tracer(typep(ipart))) then
                    ! Get old sink index
                    isink_old = partp(ipart)
                    ! Compute new sink index
                    isink_new = gsink(psink_inv(isink_old))

                    partp(ipart) = isink_new

                 end if

                 ipart = next_part  ! Go to next particle
              end do
              ! End loop over particles
           end if
           igrid = next(igrid)   ! Go to next grid
        end do

     end do
     ! End loop over cpus
  end if

  deallocate(psink,gsink)

  if (MC_tracer) then
     deallocate(psink_inv)
  end if


111 format('   Entering merge_sink for level ',I2)

end subroutine merge_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine kill_sink(ind_part,ind_grid_part,np)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine merge_sink
  ! It removes sink particles that are part of a FOF group.
  !-----------------------------------------------------------------------
  integer::j,isink,isink_new
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok

  do j=1,np
     isink=-idp(ind_part(j))
     ok(j)=(oksink_all(isink)==0)
     if(.not. ok(j))then
        isink_new=int(oksink_all(isink))
        idp(ind_part(j))=-isink_new
        mp(ind_part(j))=msink(isink_new)
        xp(ind_part(j),1)=xsink(isink_new,1)
        vp(ind_part(j),1)=vsink(isink_new,1)
#if NDIM>1
        xp(ind_part(j),2)=xsink(isink_new,2)
        vp(ind_part(j),2)=vsink(isink_new,2)
#endif
#if NDIM>2
        xp(ind_part(j),3)=xsink(isink_new,3)
        vp(ind_part(j),3)=vsink(isink_new,3)
#endif
     endif
  end do

  ! Remove particles from parent linked list
  call remove_list(ind_part,ind_grid_part,ok,np)
  call add_free_cond(ind_part,ok,np)
end subroutine kill_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine create_cloud_from_sink
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none


  !----------------------------------------------------------------------
  ! This routine creates the whole cloud of particles for each sink,
  ! Particles are produced in the right MPI domain and inserted in the
  ! linked list at level 1.
  ! The cloud radius is dble(ir_cloud)*dx_min, where dx_min is
  ! the cell size at levelmax. For cosmo runs, the cloud radius is
  ! dx_min/aexp (therefore it is constant in *physical* units).
  !----------------------------------------------------------------------

  real(dp)::scale,dx_min,rr,rmax,rmass
  integer ::icpu,isink,indp,ii,jj,kk,nx_loc,idim
  real(dp),dimension(1:ndim)::xrel
  real(dp),dimension(1:nvector,1:ndim)::xtest
  integer ,dimension(1:nvector)::ind_grid,cc,ind_cloud
  logical ,dimension(1:nvector)::ok_true
  logical,dimension(1:ndim)::period
  logical::in_box
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::cloud_counter

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ok_true=.true.

  if(numbtot(1,1)==0) return
  if(verbose)write(*,*)' Entering create_cloud_from_sink'

#if NDIM==3

  ! Level 1 linked list
  do icpu=1,ncpu
     if(numbl(icpu,1)>0)then
        ind_grid(1)=headl(icpu,1)
     endif
  end do

  period(1)=(nx==1)
  period(2)=(ny==1)
  period(3)=(nz==1)

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp

  rmax=dble(ir_cloud)*dx_min
  rmass=dble(ir_cloud_massive)*dx_min

  cloud_counter=0

  do kk=-2*ir_cloud,2*ir_cloud
     xrel(3)=dble(kk)*dx_min/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        xrel(2)=dble(jj)*dx_min/2.0
        do ii=-2*ir_cloud,2*ir_cloud
           xrel(1)=dble(ii)*dx_min/2.0
           rr=sqrt(sum(xrel**2))
           if(rr<=rmax)then
              do isink=1,nsink
                 xtest(1,1:3)=xsink(isink,1:3)+xrel(1:3)
                 in_box=.true.
                 do idim=1,ndim
                    if (period(idim) .and. xtest(1,idim)>boxlen)xtest(1,idim)=xtest(1,idim)-boxlen
                    if (period(idim) .and. xtest(1,idim)<0.)xtest(1,idim)=xtest(1,idim)+boxlen
                    if (xtest(1,idim)<0.0 .or. xtest(1,idim)>boxlen)in_box=.false.
                 end do
                 cc(1)=0
                 if(in_box)call cmp_cpumap(xtest,cc,1)
                 if(cc(1).eq.myid)then
                    call remove_free(ind_cloud,1)
                    call add_list(ind_cloud,ind_grid,ok_true,1)
                    cloud_counter=cloud_counter+1
                    indp=ind_cloud(1)
                    idp(indp)=-isink
                    levelp(indp)=nlevelmax_current
                    mp(indp)=msink(isink)/dble(ncloud_sink)
                    xp(indp,1:3)=xtest(1,1:3)
                    vp(indp,1:3)=vsink(isink,1:3)
                    tp(indp)=tsink(isink)     ! Birth epoch
                    typep(indp)%family = FAM_CLOUD
                    if((ii.eq.0).and.(jj.eq.0).and.(kk.eq.0))then
                       typep(indp)%tag = TAG_CLOUD_CENTRAL  !Central cloud particle
                    else
                       typep(indp)%tag = 0
                    endif
                 end if
              end do
           end if
        end do
     end do
  end do

#endif
end subroutine create_cloud_from_sink
!################################################################
!################################################################
!################################################################
subroutine kill_entire_cloud(ilevel)
  use pm_commons
  use amr_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine removes cloud particles (including the central one).
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::ig,ip,npart1,npart2,icpu,ncache,istart
  integer,dimension(1:nvector)::ind_grid,ind_part,ind_grid_part
  logical,dimension(1:nvector)::ok=.true.

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Gather sink and cloud particles.
  ! Loop over cpus
  do icpu=1,ncpu+nboundary
     if(icpu<=ncpu)then
        ncache=numbl(icpu,ilevel)
        istart=headl(icpu,ilevel)
     else
        ncache=numbb(icpu-ncpu,ilevel)
        istart=headb(icpu-ncpu,ilevel)
     end if
     igrid=istart
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,ncache
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0
        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if(idp(ipart).lt.0)then
              if (is_cloud(typep(ipart))) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif
       ! Gather sink and cloud particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              !if(idp(ipart).lt.0)then
              if (is_cloud(typep(ipart))) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call remove_list(ind_part,ind_grid_part,ok,ip)
                 call add_free_cond(ind_part,ok,ip)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
          ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)then
        call remove_list(ind_part,ind_grid_part,ok,ip)
        call add_free_cond(ind_part,ok,ip)
     end if
  end do
111 format('   Entering kill_cloud for level ',I2)
end subroutine kill_entire_cloud
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_hoyle(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: rhoc, mH , twopi
  use mpi_mod
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine computes the parameters of Bondi-Hoyle
  ! accretion for sink particles.
  ! It calls routine bondi_veocity and average_density.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,idim,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::r2,dx_loc,dx_min,scale
  real(dp)::factG,pi
  real(dp)::RandNum,phi,Rrand,SS,CC,UU
  integer ,dimension(1:ncpu,1:IRandNumSize)::allseed

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! If necessary, initialize random number generator
  if(localseed(1)==-1)then
     call rans(ncpu,iseed,allseed)
     localseed=allseed(myid,1:IRandNumSize)
  end if

  ! Mesh spacing in that level
  dx_loc=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp

  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Reset new sink variables
  v2sink_new=0d0; c2sink_new=0d0; oksink_new=0d0

  ! Gather sink particles only.

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count only sink particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 isink=-idp(ipart)

                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    npart2=npart2+1
                 end if
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather only sink particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 isink=-idp(ipart)
                 r2=0.0
                 do idim=1,ndim
                    r2=r2+(xp(ipart,idim)-xsink(isink,idim))**2
                 end do
                 if(r2==0.0)then
                    if(ig==0)then
                       ig=1
                       ind_grid(ig)=igrid
                    end if
                    ip=ip+1
                    ind_part(ip)=ipart
                    ind_grid_part(ip)=ig
                 endif
              endif
              if(ip==nvector)then
                 call bondi_velocity(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call bondi_velocity(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(oksink_new,oksink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(c2sink_new,c2sink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(v2sink_new,v2sink_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     oksink_all=oksink_new
     c2sink_all=c2sink_new
     v2sink_all=v2sink_new
#endif
  endif

  do isink=1,nsink
     if(oksink_all(isink)==1d0)then
        c2sink(isink)=c2sink_all(isink)
        v2sink(isink)=v2sink_all(isink)
        ! Compute sink radius
        r2sink(isink)=(factG*msink(isink)/(v2sink(isink)+c2sink(isink)))**2
        !r2sink(isink)=(factG*msink(isink)/(c2sink(isink)))**2

        ! If radius is far smaller than the resolution, spurious velocity
        ! must not be taken into account
        !if (sqrt(r2sink(isink)) .lt. dx_min/4d0) then
        !   r2sink(isink)=(factG*msink(isink)/c2sink(isink))**2
        !endif

        r2k(isink)   =min(max(r2sink(isink),(dx_min/4.0)**2),(2.*dx_min)**2)
     endif
  end do

  ! Gather sink and cloud particles.
  wdens=0d0; wvol =0d0; wc2=0d0; wmom=0d0; jsink_new=0d0

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather sink and cloud particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call average_density(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 if((.not.random_jet)) call jet_AGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call average_density(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
     !if(ip>0 .and. sink_AGN .and.(.not.random_jet)) call jet_AGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
     if(ip>0 .and.(.not.random_jet)) call jet_AGN(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(random_jet)then
     if(myid==1)then
        do isink=1,nsink
           ! Random directions
           call ranf(localseed,RandNum)
           SS =(RandNum-0.5)*2.
           call ranf(localseed,RandNum)
           phi=(RandNum-0.5)*2.*pi
           call ranf(localseed,RandNum)
           UU =RandNum
           Rrand=UU**(1./3.)
           CC=Rrand*sqrt(1.-SS**2.)
           jsink_new(isink,1)=CC*cos(phi)
           jsink_new(isink,2)=CC*sin(phi)
           jsink_new(isink,3)=Rrand*SS
        enddo
     endif
  endif

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(wdens,wdens_new,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wvol ,wvol_new ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wc2  ,wc2_new  ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(wmom ,wmom_new ,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(jsink_new,jsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     wdens_new=wdens
     wvol_new=wvol
     wc2_new =wc2
     wmom_new=wmom
     jsink_all=jsink_new
#endif
  endif

  do isink=1,nsink
     weighted_density(isink,ilevel)=wdens_new(isink)
     weighted_volume (isink,ilevel)=wvol_new (isink)
     weighted_momentum(isink,ilevel,1:ndim)=wmom_new(isink,1:ndim)
     weighted_c2     (isink,ilevel)=wc2_new  (isink)
     do i=1,ndim
        jsink(isink,i)=jsink(isink,i)+jsink_all(isink,i)
     enddo
  end do

111 format('   Entering bondi_hoyle for level ',I2)

end subroutine bondi_hoyle
!################################################################
!################################################################
!################################################################
!################################################################
subroutine bondi_velocity(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: rhoc, mH
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  ! It computes the gas velocity and soud speed in the cell
  ! each sink particle sits in.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::v2,c2,d,u,v,w,e
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp)::dx,dx_loc,scale,vol_loc
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

#if NDIM==3
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in bondi_velocity'
     write(*,*)ilevel,ng,np
     write(*,*)-idp(ind_part(j))
     write(*,*)x(j,1:3)
     write(*,*)vp(ind_part(j),1:3)
     stop
  end if

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

  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Gather hydro variables
  do j=1,np
     if(ok(j))then
        d=max(uold(indp(j),1), smallr)
        u=uold(indp(j),2)/d
        v=uold(indp(j),3)/d
        w=uold(indp(j),4)/d
        e=uold(indp(j),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        isink=-idp(ind_part(j))
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2
        c2=MAX(gamma*(gamma-1.0)*e,smallc**2)
        ! Relative velocity of the gas in regards of the sink
        u=u-vsink(isink,1)
        v=v-vsink(isink,2)
        w=w-vsink(isink,3)
        v2=(u**2+v**2+w**2)
        v2sink_new(isink)=v2
        c2sink_new(isink)=c2
        oksink_new(isink)=1d0
     endif
  end do

#endif

end subroutine bondi_velocity
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_density(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: rhoc, mH
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle. Each cloud particle
  ! reads up the value of density, sound speed and velocity from its
  ! position in the grid.
  !-----------------------------------------------------------------------
  logical::error
  integer::i,j,ind,idim,nx_loc,isink
  real(dp)::dx,scale,weight,r2
  ! Grid-based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::dgas,ugas,vgas,wgas,c2gas
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::u,v,w,d,e,v2,c2
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        !!! MT
	!!!if(x(j,idim)<1.5D0.or.x(j,idim)>4.5D0)error=.true.
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in average_density'
     do idim=1,ndim
        do j=1,np
           if(x(j,idim)<2D0.or.x(j,idim)>4D0)then
          write(*,*) 'tp =', tp(ind_part(j))
          write(*,*) 'id =', idp(ind_part(j))
!!!!!           if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)then
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

  ! Gather gas density
  dgas(1:np)=0.0D0
  ugas(1:np)=0.0D0
  vgas(1:np)=0.0D0
  wgas(1:np)=0.0D0
  c2gas(1:np)=0.0D0
  do ind=1,twotondim
     do j=1,np
        dgas(j)=dgas(j)+uold(indp(j,ind),1)*vol(j,ind)
        d=max(uold(indp(j,ind),1), smallr)
        u=uold(indp(j,ind),2)/d
        v=uold(indp(j,ind),3)/d
        w=uold(indp(j,ind),4)/d
        e=uold(indp(j,ind),5)/d
#ifdef SOLVERmhd
        bx1=uold(indp(j,ind),6)
        by1=uold(indp(j,ind),7)
        bz1=uold(indp(j,ind),8)
        bx2=uold(indp(j,ind),nvar+1)
        by2=uold(indp(j,ind),nvar+2)
        bz2=uold(indp(j,ind),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2
        c2=MAX(gamma*(gamma-1.0)*e,smallc**2)
        ! --------------------------
        ! Volume-weighted quantities
        ! (if you change this, be sure to renormalize
        ! correctly these quantities in grow_bondi)
        ! --------------------------
        !ugas(j)=ugas(j)+u*vol(j,ind)
        !vgas(j)=vgas(j)+v*vol(j,ind)
        !wgas(j)=wgas(j)+w*vol(j,ind)
        !c2gas(j)=c2gas(j)+c2*vol(j,ind)
        ! --------------------------
        ! Mass-weighted quantities
        ! --------------------------
        ugas(j)=ugas(j)+d*u*vol(j,ind)
        vgas(j)=vgas(j)+d*v*vol(j,ind)
        wgas(j)=wgas(j)+d*w*vol(j,ind)
        c2gas(j)=c2gas(j)+d*c2*vol(j,ind)
     end do
  end do

  do j=1,np
     isink=-idp(ind_part(j))
     r2=0d0
     do idim=1,ndim
        r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
     end do
     weight=exp(-r2/r2k(isink))
     wdens(isink)=wdens(isink)+weight*dgas(j)
     wmom(isink,1)=wmom(isink,1)+weight*ugas(j)
     wmom(isink,2)=wmom(isink,2)+weight*vgas(j)
     wmom(isink,3)=wmom(isink,3)+weight*wgas(j)
     wc2  (isink)=wc2  (isink)+weight*c2gas(j)
     wvol (isink)=wvol (isink)+weight
  end do

end subroutine average_density
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_bondi(ilevel)
  use pm_commons
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH, twopi
  ! AGNRT
#ifdef RT
  use rt_parameters,only: rt_AGN
#endif
  !/AGNRT
  use mpi_mod

  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine performs Bondi-Hoyle accretion of the gas onto
  ! sink particles. On exit, sink mass and velocity are modified.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::i,ig,ip,npart1,npart2,icpu,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  real(dp)::density,volume
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  integer::ind,ivar,iskip
  real(dp)::alpha,d_star,pi,factG,c2mean,nfloor,sigmav2,v2mean
  real(dp)::ZZ1,ZZ2,r_lso,epsilon_r,onethird
  real(dp),dimension(1:3)::velocity
  real(dp)::fourpi,prefact
  real(dp)::r_bondi,r_acc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nfloor = omega_b*rhoc/aexp**3*XH/mH/scale_nH

  onethird=1d0/3d0
  fourpi=2d0*twopi
  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp
  prefact=fourpi*6.67d-8*scale_m*1.66d-24/(6.652d-25*3d10)/(scale_m/scale_t)

  ! Reset new sink variables
  msink_new=0d0; vsink_new=0d0; dMBH_coarse_new=0d0; dMEd_coarse_new=0d0; dMsmbh_new=0d0
  ! AGNRT
#ifdef RT
  if (rt_AGN) dMeff_coarse_new = 0.d0
#endif

  ! Set unew to uold for myid cells
  ! Need unew to get initial density before Bondi accreting mass
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do ivar=1,nvar
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        enddo
     enddo
  enddo

  d_star=0d0
  if (star)d_star=n_star/scale_nH
  ! From km/s to user units: assume a 10 km/s vel disp in the ISM
  sigmav2=(sigmav_max*1d5/scale_v)**2d0

  c_avgptr=0d0;v_avgptr=0d0;d_avgptr=0d0

  ! Compute Bondi-Hoyle accretion rate
  do isink=1,nsink
     density=0d0
     c2mean=0d0
     velocity=0d0
     volume=0d0
     do i=levelmin,nlevelmax
        density=density+weighted_density(isink,i)
        c2mean =c2mean +weighted_c2     (isink,i)
        velocity(1)=velocity(1)+weighted_momentum(isink,i,1)
        velocity(2)=velocity(2)+weighted_momentum(isink,i,2)
        velocity(3)=velocity(3)+weighted_momentum(isink,i,3)
        volume =volume +weighted_volume (isink,i)
     end do
     density=max(density/volume, smallr)
     ! --------------------------
     ! If volume-weighted
     ! --------------------------
     !velocity(1:3)=velocity(1:3)/volume
     !c2mean =c2mean /volume
     ! --------------------------
     ! If mass-weighted
     ! --------------------------
     velocity(1:3)=velocity(1:3)/volume/density
     c2mean =c2mean /volume/density
     ! --------------------------
     v_avgptr(isink)=dsqrt(SUM((velocity(1:3)-vsink(isink,1:3))**2))
     v2mean =min(SUM((velocity(1:3)-vsink(isink,1:3))**2),sigmav2)
     total_volume(isink)=volume
     alpha=max((density/(d_boost/scale_nH))**boost_acc,1d0)
     if(Esave(isink).eq.0d0)then
        dMBHoverdt(isink)=alpha * fourpi *density* (factG*msink(isink))**2 &
             & / (c2mean+v2mean)**1.5d0
     else
        ! Prevent the accretion of new material onto the BH if
        ! energy has not been released in the previous coarse time step
        dMBHoverdt(isink)=0d0
     endif

     if(spin_bh)then
        ZZ1=1d0+(1d0-spinmag(isink)**2)**onethird*((1d0+spinmag(isink))**onethird &
             & +(1d0-spinmag(isink))**onethird)
        ZZ2=SQRT(3d0*spinmag(isink)**2+ZZ1**2)
        if(spinmag(isink).ge.0d0)then
           r_lso=3d0+ZZ2-SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
        else
           r_lso=3d0+ZZ2+SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
        endif
        epsilon_r=1d0-sqrt(1d0-2d0/(3d0*r_lso))
     else
        epsilon_r=0.1d0
     endif
     dMEdoverdt(isink)=prefact*msink(isink)/epsilon_r
     ! Force accretion at Eddington (MT)
     if (force_accretion) then
        dMBHoverdt(isink) = dMEdoverdt(isink)  !!! SHOULD THIS BE 0 IF ESAVE .NEQ. 0d0?
     endif
     ! Force accretion at Eddington (MT)
     if(dMBHoverdt(isink)/dMEdoverdt(isink) .lt. X_floor)dMBHoverdt(isink)=f_bondi*dMBHoverdt(isink)
     c_avgptr(isink)=dsqrt(c2mean)
     d_avgptr(isink)=density

     ! Calculate whether to switch drag force on or off
     r_bondi=factG*msink(isink)/c2mean  ! Bondi radius
     if(v2mean>0)then
        r_acc=2*factG*msink(isink)/v2mean ! Accretion radius
        rg_scale(isink)=MIN(r_bondi,r_acc)           ! gravitational scale radius
     else
        rg_scale(isink)=r_bondi
     end if
  end do

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather sink and cloud particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call accrete_bondi(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call accrete_bondi(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dMBH_coarse_new,dMBH_coarse_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dMEd_coarse_new,dMEd_coarse_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(dMsmbh_new     ,dMsmbh_all     ,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
! AGNRT
#ifdef RT
     if(rt_AGN) call MPI_ALLREDUCE(dMeff_coarse_new,dMeff_coarse_all,nsinkmax,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#endif
!/AGNRT
#else
     msink_all=msink_new
     vsink_all=vsink_new
     dMBH_coarse_all=dMBH_coarse_new
     dMEd_coarse_all=dMEd_coarse_new
     dMsmbh_all     =dMBsmbh_new
#ifdef RT
     if(rt_AGN) dMeff_coarse_all = dMeff_coarse_new
#endif
#endif
  endif

  do isink=1,nsink
     if(.not.fix_smbh_position)then
        vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)+vsink_all(isink,1:ndim)
        msink(isink)       =msink(isink)       +msink_all(isink)
        vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)
     else
        msink(isink)       =msink(isink)       +msink_all(isink)
     endif
     dMBH_coarse(isink)=dMBH_coarse(isink)+dMBH_coarse_all(isink)
     dMEd_coarse(isink)=dMEd_coarse(isink)+dMEd_coarse_all(isink)
     dMsmbh     (isink)=dMsmbh     (isink)+dMsmbh_all     (isink)
#ifdef RT
     if(rt_AGN) dMeff_coarse(isink) = dMeff_coarse(isink) + dMeff_coarse_all(isink)
#endif
  end do

  ! YDspin
  if(spin_bh)call growspin
  ! YDspin

111 format('   Entering grow_bondi for level ',I2)

end subroutine grow_bondi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_bondi(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: rhoc, mH, twopi
  ! AGNRT
#ifdef RT
  use rt_parameters,only: nGroups, iGroups, group_egy, rt_AGN
  use rt_hydro_commons,only: rtuold
#endif
  !/AGNRT
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::r2,d,u,v,w,e
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp)::dx,dx_loc,scale,vol_loc,weight,acc_mass
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::floorB,d_ini,dmsink
  real(dp)::dx_min
  integer ::counter
  real(dp)::dmom,ekk,dvdrag,vnorm_rel,factor,mach,alpha,cgas,fudge,factG
  real(dp),dimension(1:ndim),save::dpdrag,vrel,fdrag
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,pi,d_star
  real(dp)::epsilon_r, Deltat, ddt, dvdrag_norm
#ifdef RT
  ! AGNRT
  real(dp)::Np_inj
  integer ::mygroup
  !/AGNRT
#endif
  real(dp), dimension(imetal:nvar):: utmp_passive
  integer::ivar

  ! Tracers
  integer, dimension(1:nvector), save :: ind_tracer, ind_sink
  real(dp), dimension(1:nvector), save :: proba_tracer
  real(dp), dimension(1:nvector, 1:ndim), save :: xsink_loc
  real(dp) :: proba
  integer :: itracer, ii

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  pi=twopi/2d0
  factG=1d0
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp
  fudge=4d0*twopi*factG**2
  d_star=1d100
  if (star)then
     d_star=n_star/scale_nH
  else
     d_star=5d1/scale_nH
  endif


  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp

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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in accrete_bondi'
     write(*,*)ilevel,ng,np
     stop
  end if

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
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
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
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! MC tracer
  ! Init local index
  ii = 0

  ! Remove mass from hydro cells
  do j=1,np
     if(ok(j))then
        isink=-idp(ind_part(j))
        r2=0d0
        do idim=1,ndim
           r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
        end do
        weight=exp(-r2/r2k(isink))

        d=max(uold(indp(j),1), smallr)
        u=uold(indp(j),2)/d
#if NDIM>1
        v=uold(indp(j),3)/d
#endif
#if NDIM>2
        w=uold(indp(j),4)/d
#endif
        e=uold(indp(j),ndim+2)/d
        do ivar=imetal,nvar
           utmp_passive(ivar) = uold(indp(j), ivar)/d
        end do

#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        d_ini=max(unew(indp(j),1), smallr)

        !The hard density floor is needed in case the sink particle
        !remains in the same cell throughout, e.g. via fix_smbh_position=.true.
        floorB=max(0.75*d_ini,1d-10)

        ! Allow for zero accretion (HP)
        if (no_accretion) then
           acc_mass = 0d0
        else if (maximum_accretion) then
           ! Accrete up to 99 % of the cell mass, kernel weighted
           ! for smoother transition at accretion region edge
           floorB=max(0.01*d_ini,1d-10)
           acc_mass=(d-floorB)*vol_loc
        else if (eddington_limit) then
           acc_mass=min(dMBHoverdt(isink),dMEdoverdt(isink)) &
                & *weight/total_volume(isink)*dtnew(ilevel)
        else
           acc_mass=dMBHoverdt(isink)* weight/total_volume(isink)*dtnew(ilevel)
        end if
        acc_mass=max( min(acc_mass, (d-floorB)*vol_loc), 0d0)

        if(spin_bh)then
           epsilon_r=eps_sink(isink)
        else
           epsilon_r=0.1d0
        endif
        if(mad_jet)then
           if(dMBHoverdt(isink)/dMEdoverdt(isink).lt.X_floor)then
              ! from Benson & Babul 2009, for an ADAF
              epsilon_r=epsilon_r*(dMBHoverdt(isink)/(X_floor*dMEdoverdt(isink)))
           endif
        endif

        ! If there is no feedback, accrete 100% of the removed mass
        if(.not.sink_AGN)then
           epsilon_r=0d0
        endif

        dmsink=acc_mass*(1d0-epsilon_r)

        ! Add the accreted mass to the total accreted mass over
        ! a coarse time step
        dMBH_coarse_new(isink)=dMBH_coarse_new(isink) + &
             & dMBHoverdt(isink)*weight/total_volume(isink)*dtnew(ilevel)
        dMEd_coarse_new(isink)=dMEd_coarse_new(isink) + &
             & dMEdoverdt(isink)*weight/total_volume(isink)*dtnew(ilevel)
        dMsmbh_new     (isink)=dMsmbh_new     (isink) + dmsink

! AGNRT
#ifdef RT
        if (rt_AGN) dMeff_coarse_new(isink) = dMeff_coarse_new     (isink) + dmsink
#endif
!/AGNRT

        msink_new(isink  )=msink_new(isink  )+dmsink
        vsink_new(isink,1)=vsink_new(isink,1)+dmsink*u
#if NDIM>1
        vsink_new(isink,2)=vsink_new(isink,2)+dmsink*v
#endif
#if NDIM>2
        vsink_new(isink,3)=vsink_new(isink,3)+dmsink*w
#endif

        vp(ind_part(j),1)=mp(ind_part(j))*vp(ind_part(j),1)+dmsink*u
        vp(ind_part(j),2)=mp(ind_part(j))*vp(ind_part(j),2)+dmsink*v
        vp(ind_part(j),3)=mp(ind_part(j))*vp(ind_part(j),3)+dmsink*w
        mp(ind_part(j))=mp(ind_part(j))+dmsink
        vp(ind_part(j),1)=vp(ind_part(j),1)/mp(ind_part(j))
        vp(ind_part(j),2)=vp(ind_part(j),2)/mp(ind_part(j))
        vp(ind_part(j),3)=vp(ind_part(j),3)/mp(ind_part(j))

        if (MC_tracer) then
           ! MC Tracer
           ! The gas is coming from the central cell (and the tracers)
           proba = acc_mass / vol_loc / d
           itracer = headp(igrid(j))

           do i = 1, numbp(igrid(j))
              ! Select gas tracer particles within current cell that
              ! havn't been moved.
              if (is_gas_tracer(typep(itracer)) .and. move_flag(itracer) == 0 .and. &
                   partp(itracer) == indp(j)) then
                 ii = ii + 1
                 ind_tracer(ii) = itracer
                 proba_tracer(ii) = proba
                 ! TODO: check whether to use xsink/xsink_new
                 xsink_loc(ii, :) = xsink_new(isink, :)
                 ind_sink(ii) = isink

                 if (ii == nvector) then
                    call tracer2sink(ind_tracer, proba_tracer, &
                         xsink_loc, ind_sink, ii, dx_loc)
                    ii = 0
                 end if
              end if
              itracer = nextp(itracer)
           end do
        end if

        d=d-acc_mass/vol_loc

#ifdef SOLVERmhd
        e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        uold(indp(j),1)=max(d, smallr)
        uold(indp(j),2)=d*u
#if NDIM>1
        uold(indp(j),3)=d*v
#endif
#if NDIM>2
        uold(indp(j),4)=d*w
#endif
        uold(indp(j),ndim+2)=d*e
        do ivar=imetal,nvar
           uold(indp(j), ivar) = d*utmp_passive(ivar)
        end do

        ! Only use gas drag when scale radius is unresolved
        if(drag .and. rg_scale(isink) <= dx_min*0.2)then
           ! HP: updates on the gas DF
           Deltat=0d0; counter=0
           do while ((Deltat .lt. dtnew(ilevel)) .and. (counter .le. 10))
              ddt=dtnew(ilevel) - Deltat
              d=max(uold(indp(j),1), smallr)
              u=uold(indp(j),2)/d
              v=uold(indp(j),3)/d
              w=uold(indp(j),4)/d
              ekk=0.5d0*d*(u*u+v*v+w*w)
              e=uold(indp(j),ndim+2)-ekk
#ifdef SOLVERmhd
           e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif
              cgas=max(sqrt(gamma*(gamma-1d0)*e/d), smallc)
              ! Compute the drag force exerted by the gas on the sink particle
              vrel(1)=vp(ind_part(j),1)-u
              vrel(2)=vp(ind_part(j),2)-v
              vrel(3)=vp(ind_part(j),3)-w
              vnorm_rel=max(sqrt( vrel(1)**2 + vrel(2)**2 + vrel(3)**2 ), smallc)
              mach=vnorm_rel/cgas
              alpha=max((d/d_boost)**boost_drag,1d0)
              factor=alpha * fudge * d*msink(isink)/cgas**2 / vnorm_rel
              ! TODO: deal with very low mach number
              if(mach.le.0.950d0)factor=factor/mach**2*(0.5d0*log((1d0+mach)/(1d0-mach))-mach)
              if(mach.ge.1.007d0)factor=factor/mach**2*(0.5d0*log(mach**2-1d0)+3.2d0)
              ! HP: updates on the gas DF
              dvdrag_norm = factor * ddt * vnorm_rel
              if ((dvdrag_norm/vnorm_rel .ge. 0.1) .and. (counter .le. 9)) then
                 ddt=0.1/(dvdrag_norm/vnorm_rel)*ddt
              end if
              factor = min(1/ddt, factor)
              !/HP

              do idim=1,ndim
                 fdrag(idim)= - factor * vrel(idim)
                 dvdrag = fdrag(idim)*ddt  ! HP: replaced dtnew(ilevel) by ddt
                 dpdrag(idim)=mp(ind_part(j))*dvdrag
                 vp(ind_part(j),idim)=vp(ind_part(j),idim)+dvdrag
                 vsink_new(isink,idim)=vsink_new(isink,idim)+dvdrag*mp(ind_part(j))
              enddo

              ! HP: updates on the gas DF
              Deltat = Deltat + ddt
              counter = counter +1
              !/HP
              ekk=0.0d0
              do idim=1,ndim
                 dmom=dpdrag(idim)/vol_loc
                 uold(indp(j),idim+1)=uold(indp(j),idim+1)-dmom
                 ekk=ekk+0.5d0*uold(indp(j),idim+1)**2/d
              enddo
              ! Update the total energy with new corresponding kinetic energy
              e=e+ekk
#ifdef SOLVERmhd
              e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)
#endif

              uold(indp(j),ndim+2)=e
           end do
        endif

        ! AGNRT
#ifdef RT
        if(rt_AGN) then
           do mygroup = 1, nGroups
              Np_inj = LAGN_coarse(isink)*dtnew(ilevel) * weight/total_volume(isink) * lumfrac_AGN(isink, mygroup)/group_egy(mygroup)
              rtuold(indp(j),iGroups(mygroup))=rtuold(indp(j),iGroups(mygroup)) + (Np_inj/vol_loc)
           end do
        endif
#endif
        !/AGNRT

     endif

  end do

  if (MC_tracer .and. ii > 0) then
     call tracer2sink(ind_tracer, proba_tracer, &
          xsink_loc, ind_sink, ii, dx_loc)
  end if


end subroutine accrete_bondi
!################################################################
!################################################################
!################################################################
!################################################################
subroutine grow_jeans(ilevel)
  use pm_commons
  use amr_commons
  use cooling_module, ONLY: XH=>X, rhoc, mH
  use mpi_mod
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine determines if a cell covered by an old sink particle
  ! cross the density threshold for sink particle formation.If so, the
  ! density is reduced and the corresponding mass is given to the sink.
  ! On exit, sink mass and velocity are modified.
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part,info
  integer::i,ig,ip,npart1,npart2,icpu,nx_loc,isink
  integer,dimension(1:nvector),save::ind_grid,ind_part,ind_grid_part
  integer::nlevelmax_loc
  real(dp)::dx_min,vol_min,dx_temp
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale
  real(dp)::d0,mstar,nISM,nCOM

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Reset new sink variables
  msink_new=0d0; vsink_new=0d0

  ! Finest cell size
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)
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

  ! Ensure that sinks accrete according to jeans only from the most refined cells
  if(ilevel .lt. nlevelmax_loc)return

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count sink and cloud particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 npart2=npart2+1
              endif
              ipart=next_part  ! Go to next particle
           end do
        endif

        ! Gather sink and cloud particles
        if(npart2>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              ! Select only sink particles
              !if(idp(ipart).lt.0 .and. tp(ipart).eq.0.d0)then
              if (is_cloud(typep(ipart))) then
                 if(ig==0)then
                    ig=1
                    ind_grid(ig)=igrid
                 end if
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
              endif
              if(ip==nvector)then
                 call accrete_jeans(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
                 ip=0
                 ig=0
              end if
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
        end if

        igrid=next(igrid)   ! Go to next grid
     end do

     ! End loop over grids
     if(ip>0)call accrete_jeans(ind_grid,ind_part,ind_grid_part,ig,ip,ilevel)
  end do
  ! End loop over cpus

  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(msink_new,msink_all,nsinkmax     ,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(vsink_new,vsink_all,nsinkmax*ndim,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     msink_all=msink_new
     vsink_all=vsink_new
#endif
  endif
  do isink=1,nsink
     vsink(isink,1:ndim)=vsink(isink,1:ndim)*msink(isink)+vsink_all(isink,1:ndim)
     msink(isink)       =msink(isink)       +msink_all(isink)
     vsink(isink,1:ndim)=vsink(isink,1:ndim)/msink(isink)
  end do

111 format('   Entering grow_jeans for level ',I2)

end subroutine grow_jeans
!################################################################
!################################################################
!################################################################
!################################################################
subroutine accrete_jeans(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use cooling_module, ONLY: rhoc, mH, twopi
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine bondi_hoyle.
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::v2,d,u,v=0d0,w=0d0,e
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif
  real(dp)::dx,dx_loc,scale,vol_loc,temp,d_jeans,acc_mass,dd_sink,d_thres
  real(dp)::factG,pi
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp), dimension(imetal:nvar):: utmp_passive
  integer::ivar


  ! Tracers
  integer, dimension(1:nvector), save :: ind_tracer, ind_sink
  real(dp), dimension(1:nvector), save :: proba_tracer
  real(dp), dimension(1:nvector, 1:ndim), save :: xsink_loc
  integer :: ii, itracer

  if (MC_tracer) then
     if (myid == 1) then
        print*, '__________ WARNING ! __________'
        print*, 'jeans accretion not tested with MC tracers'
     end if
  end if
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0

  pi=twopi/2d0
  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Density threshold for sink particle formation
  dd_sink=n_sink/scale_nH

#if NDIM==3
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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in accrete_jeans'
     write(*,*)ilevel,ng,np
     stop
  end if

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
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
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
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! MC tracer
  ! Init local index
  ii = 0

  ! Gather hydro variables
  do j=1,np
     if(ok(j))then

        ! Convert to primitive variables
        d=max(uold(indp(j),1), smallr)
        u=uold(indp(j),2)/d
#if NDIM>1
        v=uold(indp(j),3)/d
#endif
#if NDIM>2
        w=uold(indp(j),4)/d
#endif
        e=uold(indp(j),ndim+2)/d
        do ivar=imetal,nvar
           utmp_passive(ivar) = uold(indp(j), ivar)/d
        end do
#ifdef SOLVERmhd
        bx1=uold(indp(j),6)
        by1=uold(indp(j),7)
        bz1=uold(indp(j),8)
        bx2=uold(indp(j),nvar+1)
        by2=uold(indp(j),nvar+2)
        bz2=uold(indp(j),nvar+3)
        e=e-0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
        v2=(u**2+v**2+w**2)
        e=e-0.5d0*v2

        ! Jeans length related density threshold
        temp=max(e*(gamma-1.0),smallc**2)
        d_jeans=temp*3.1415926/(4.0*dx_loc)**2/factG
        d_thres=d_jeans

        ! User defined density threshold
        !d_thres=dd_sink

        if(d .ge. d_thres)then
           isink=-idp(ind_part(j))
           acc_mass=min((d-d_thres/4.0)*vol_loc,1d5*2d33/scale_m)
           !write(*,*)'Masse accretee par accretion de Jeans=',acc_mass*scale_m/2d33,' Msun'
           msink_new(isink  )=msink_new(isink  )+acc_mass
           vsink_new(isink,1)=vsink_new(isink,1)+acc_mass*u
#if NDIM>1
           vsink_new(isink,2)=vsink_new(isink,2)+acc_mass*v
#endif
#if NDIM>2
           vsink_new(isink,3)=vsink_new(isink,3)+acc_mass*w
#endif
           !d=d_thres/4.0
           ! MC Tracer
           if (MC_tracer) then
              itracer = headp(igrid(j))
              do i = 1, numbp(igrid(j))
                 ! Select tracer particles within current cell that
                 ! haven't been moved.
                 if (is_gas_tracer(typep(itracer)) .and. move_flag(itracer) == 0 .and. &
                      partp(itracer) == icell(j)) then
                    ii = ii + 1
                    ind_tracer(ii) = itracer
                    proba_tracer(ii) = acc_mass / d
                    ! TODO: check it's xsink and not xsink_new
                    xsink_loc(ii, :) = xsink_new(isink, :)
                    ind_sink(ii) = isink

                    if (ii == nvector) then
                       call tracer2sink(ind_tracer, proba_tracer, &
                            xsink_loc, ind_sink, ii, dx_loc)
                       ii = 0
                    end if
                 end if
                 itracer = nextp(itracer)
              end do
           end if

           ! Accrete mass (Jeans accretion)
           d=d-acc_mass

           mp(ind_part(j))=mp(ind_part(j))+acc_mass

           ! Convert back to conservative variable
#ifdef SOLVERmhd
           e=e+0.125d0*((bx1+bx2)**2+(by1+by2)**2+(bz1+bz2)**2)/d
#endif
           e=e+0.5d0*(u**2+v**2+w**2)
           uold(indp(j),1)=d
           uold(indp(j),2)=d*u
#if NDIM>1
           uold(indp(j),3)=d*v
#endif
#if NDIM>2
           uold(indp(j),4)=d*w
#endif
           uold(indp(j),ndim+2)=d*e
           do ivar=imetal,nvar
              uold(indp(j), ivar) = d*utmp_passive(ivar)
           end do
        endif

     endif
  end do
  ! End loop over particles

  ! MC tracer
  ! Empty buffer
  if (MC_tracer .and. ii > 0) then
     call tracer2sink(ind_tracer, proba_tracer, &
          xsink_loc, ind_sink, ii, dx_loc)
  end if

#endif

end subroutine accrete_jeans
!################################################################
!################################################################
!################################################################
!################################################################
subroutine jet_AGN(ind_grid,ind_part,ind_grid_part,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ng,np,ilevel
  integer,dimension(1:nvector)::ind_grid
  integer,dimension(1:nvector)::ind_grid_part,ind_part
  !-----------------------------------------------------------------------
  ! This routine computes the angular momentum vector of the gas around
  ! the sink
  !-----------------------------------------------------------------------
  integer::i,j,idim,nx_loc,isink
  real(dp)::r2,d,u,v,w
#ifdef SOLVERmhd
  real(dp)::bx1,bx2,by1,by2,bz1,bz2
#endif

  real(dp)::dx,dx_loc,scale,vol_loc,weight
  logical::error
  ! Grid based arrays
  real(dp),dimension(1:nvector,1:ndim),save::x0
  integer ,dimension(1:nvector),save::ind_cell
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle based arrays
  logical,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector,1:ndim),save::x
  integer ,dimension(1:nvector,1:ndim),save::id,igd,icd
  integer ,dimension(1:nvector),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc
  real(dp)::d_x,d_y,d_z
  real(dp)::uc,vc,wc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

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

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<=0.0D0.or.x(j,idim)>=6.0D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in jet_AGN'
     write(*,*)ilevel,ng,np
     stop
  end if

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
#if NDIM==1
  do j=1,np
     kg(j)=1+igd(j,1)
  end do
#endif
#if NDIM==2
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     kg(j)=1+igd(j,1)+3*igd(j,2)+9*igd(j,3)
  end do
#endif
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
#if NDIM==1
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)
     end if
  end do
#endif
#if NDIM==2
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)
     end if
  end do
#endif
#if NDIM==3
  do j=1,np
     if(ok(j))then
        icell(j)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
     end if
  end do
#endif

  ! Compute parent cell adress
  do j=1,np
     if(ok(j))then
        indp(j)=ncoarse+(icell(j)-1)*ngridmax+igrid(j)
     end if
  end do

  ! Remove mass from hydro cells
  do j=1,np
     if(ok(j))then
        isink=-idp(ind_part(j))
        r2=0d0
        do idim=1,ndim
           r2=r2+(xp(ind_part(j),idim)-xsink(isink,idim))**2
        end do
        weight=exp(-r2/r2k(isink))

        d=max(uold(indp(j),1), smallr)
        u=uold(indp(j),2)/d
#if NDIM>1
        v=uold(indp(j),3)/d
#endif
#if NDIM>2
        w=uold(indp(j),4)/d
#endif

        d_x=xp(ind_part(j),1)-xsink(isink,1)
        d_y=xp(ind_part(j),2)-xsink(isink,2)
        d_z=xp(ind_part(j),3)-xsink(isink,3)
        !dr =sqrt(d_x*d_x+d_y*d_y+d_z*d_z)
        uc =u-vsink(isink,1)
        vc =v-vsink(isink,2)
        wc =w-vsink(isink,3)

        !j_sp= d_x*uc + d_y*vc + d_z*wc
        jsink_new(isink,1)=jsink_new(isink,1) + (d_y*wc-d_z*vc)*d*vol_loc
        jsink_new(isink,2)=jsink_new(isink,2) + (d_z*uc-d_x*wc)*d*vol_loc
        jsink_new(isink,3)=jsink_new(isink,3) + (d_x*vc-d_y*uc)*d*vol_loc

     endif
  enddo

  do j=1,np
     isink=-idp(ind_part(j))
  enddo

end subroutine jet_AGN
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine get_rho_star(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  use cooling_module
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Particles that are not entirely in
  ! level ilevel contribute also to the level density field
  ! (boundary particles) using buffer grids.
  ! Array rho_star is stored with:
  ! - rho_star containing the poisson source term fo stars
  !------------------------------------------------------------------
  integer::iskip,icpu,ind,i,nx_loc,ibound
  real(dp)::dx,scale,dx_loc

  if(.not. poisson)return
  if(numbtot(1,ilevel)==0)return

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  !--------------------------
  ! Initialize rho_star to zero
  !--------------------------
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
        rho_star(active(ilevel)%igrid(i)+iskip)=0.0D0
     end do
  end do

  !-------------------------------------------------------
  ! Initialize rho_star to zero in virtual boundaries
  !-------------------------------------------------------
  do icpu=1,ncpu
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,reception(icpu,ilevel)%ngrid
           rho_star(reception(icpu,ilevel)%igrid(i)+iskip)=0.0D0
        end do
     end do
  end do

  !---------------------------------------------------------
  ! Compute star particle contribution to density field
  !---------------------------------------------------------
  ! Compute density due to current level particles
  call rhostar_from_current_level(ilevel)
  ! Update boudaries
  call make_virtual_reverse_dp(rho_star(1),ilevel)
  call make_virtual_fine_dp   (rho_star(1),ilevel)

  !----------------------------------------------------
  ! Reset rho_star in physical boundaries
  !----------------------------------------------------
  do ibound=1,nboundary
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,boundary(ibound,ilevel)%ngrid
           rho_star(boundary(ibound,ilevel)%igrid(i)+iskip)=0.0
        end do
     end do
  end do

end subroutine get_rho_star
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine rhostar_from_current_level(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme from particles that are not entirely in
  ! level ilevel (boundary particles).
  !------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,idim,icpu,local_count
  integer::i,ig,ip,npart1
  real(dp)::dx

  integer,dimension(1:nvector),save::ind_grid,ind_cell
  integer,dimension(1:nvector),save::ind_part,ind_grid_part
  real(dp),dimension(1:nvector,1:ndim),save::x0

  ! Mesh spacing in that level
  dx=0.5D0**ilevel

  ! Loop over cpus
  do icpu=1,ncpu
     ! Loop over grids
     igrid=headl(icpu,ilevel)
     ig=0
     ip=0
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        if(npart1>0)then
           ig=ig+1
           ind_grid(ig)=igrid
           ipart=headp(igrid)
           
           local_count = 0
           ! Loop over particles
           do jpart=1,npart1
              if(ig==0)then
                 ig=1
                 ind_grid(ig)=igrid
              end if

              ! Select only stars
              if (is_star(typep(ipart))) then
                 local_count = local_count + 1
                 ip=ip+1
                 ind_part(ip)=ipart
                 ind_grid_part(ip)=ig
                 if(ip==nvector)then
                    ! Lower left corner of 3x3x3 grid-cube
                    do idim=1,ndim
                       do i=1,ig
                          x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
                       end do
                    end do
                    do i=1,ig
                       ind_cell(i)=father(ind_grid(i))
                    end do
                    call cic_star(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
                    ip=0
                    ig=0
                    local_count = 0
                 end if
              end if

              ipart=nextp(ipart)  ! Go to next particle
           end do
           ! End loop over particles

           ! Decrease grid counter for grids with no stars
           if (local_count == 0) then
              ind_grid(ig) = 0
              ig = ig - 1
           end if
        end if

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

     if(ip>0)then
        ! Lower left corner of 3x3x3 grid-cube
        do idim=1,ndim
           do i=1,ig
              x0(i,idim)=xg(ind_grid(i),idim)-3.0D0*dx
           end do
        end do
        do i=1,ig
           ind_cell(i)=father(ind_grid(i))
        end do
        call cic_star(ind_cell,ind_part,ind_grid_part,x0,ig,ip,ilevel)
     end if

  end do
  ! End loop over cpus

end subroutine rhostar_from_current_level
!##############################################################################
!##############################################################################
!##############################################################################
!##############################################################################
subroutine cic_star(ind_cell,ind_part,ind_grid_part,x0,ng,np,ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer, intent(in)::ng,np,ilevel
  integer ,dimension(1:nvector), intent(in)::ind_cell,ind_grid_part,ind_part
  real(dp),dimension(1:nvector,1:ndim), intent(in)::x0
  !------------------------------------------------------------------
  ! This routine computes the density field at level ilevel using
  ! the CIC scheme. Only cells that are in level ilevel
  ! are updated by the input particle list.
  !------------------------------------------------------------------
  logical::error
  integer::j,ind,idim,nx_loc
  real(dp)::dx,dx_loc,scale,vol_loc
  ! Grid-based arrays
  integer ,dimension(1:nvector,1:threetondim),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim),save::nbors_father_grids
  ! Particle-based arrays
  logical ,dimension(1:nvector),save::ok
  real(dp),dimension(1:nvector),save::mmm
  real(dp),dimension(1:nvector),save::vol2
  real(dp),dimension(1:nvector,1:ndim),save::x,dd,dg
  integer ,dimension(1:nvector,1:ndim),save::ig,id,igg,igd,icg,icd
  real(dp),dimension(1:nvector,1:twotondim),save::vol
  integer ,dimension(1:nvector,1:twotondim),save::igrid,icell,indp,kg
  real(dp),dimension(1:3)::skip_loc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Gather neighboring father cells (should be present anytime !)
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ng,ilevel)

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

  ! Gather particle mass
  do j=1,np
     mmm(j)=mp(ind_part(j))
  end do

  ! Check for illegal moves
  error=.false.
  do idim=1,ndim
     do j=1,np
        if(x(j,idim)<0.5D0.or.x(j,idim)>5.5D0)error=.true.
     end do
  end do
  if(error)then
     write(*,*)'problem in cic_star'
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

  ! Compute parent cell position
  do idim=1,ndim
     do j=1,np
        icg(j,idim)=ig(j,idim)-2*igg(j,idim)
        icd(j,idim)=id(j,idim)-2*igd(j,idim)
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
     icell(j,1)=1+icg(j,1)+2*icg(j,2)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)
  end do
#endif
#if NDIM==3
  do j=1,np
     icell(j,1)=1+icg(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,2)=1+icd(j,1)+2*icg(j,2)+4*icg(j,3)
     icell(j,3)=1+icg(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,4)=1+icd(j,1)+2*icd(j,2)+4*icg(j,3)
     icell(j,5)=1+icg(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,6)=1+icd(j,1)+2*icg(j,2)+4*icd(j,3)
     icell(j,7)=1+icg(j,1)+2*icd(j,2)+4*icd(j,3)
     icell(j,8)=1+icd(j,1)+2*icd(j,2)+4*icd(j,3)
  end do
#endif

  ! Compute parent cell adress
  do ind=1,twotondim
     do j=1,np
        indp(j,ind)=ncoarse+(icell(j,ind)-1)*ngridmax+igrid(j,ind)
     end do
  end do

  ! Update mass density field
  do ind=1,twotondim
     do j=1,np
        ok(j)=igrid(j,ind)>0
     end do

     do j=1,np
        vol2(j)=mmm(j)*vol(j,ind)/vol_loc
     end do

     do j=1,np
        if (ok(j)) then
           rho_star(indp(j,ind))=rho_star(indp(j,ind))+vol2(j)
        end if
     end do
  end do

end subroutine cic_star
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine quenching(ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !------------------------------------------------------------------------
  ! This routine selects regions which are eligible for SMBH formation.
  ! It is based on a stellar density threshold and on a stellar velocity
  ! dispersion threshold.
  ! On exit, flag2 array is set to 0 for AGN sites and to 1 otherwise.
  !------------------------------------------------------------------------
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::str_d,tot_m,ave_u,ave_v,ave_w,sig_u,sig_v,sig_w
  integer::igrid,ipart,jpart,next_part,ind_cell,iskip,ind
  integer::i,npart1,npart2,nx_loc
  real(dp),dimension(1:3)::skip_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

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
  vol_loc=dx_loc**ndim

#if NDIM==3
  ! Gather star particles only.

  ! Loop over grids
  do i=1,active(ilevel)%ngrid
     igrid=active(ilevel)%igrid(i)
     ! Number of particles in the grid
     npart1=numbp(igrid)
     npart2=0

     ! Reset velocity moments
     str_d=0.0
     tot_m=0.0
     ave_u=0.0
     ave_v=0.0
     ave_w=0.0
     sig_u=0.0
     sig_v=0.0
     sig_w=0.0

     ! Count star particles
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
           ! Save next particle   <--- Very important !!!
           next_part=nextp(ipart)
           !if(idp(ipart).gt.0.and.tp(ipart).ne.0)then
           !if(tp(ipart).ne.0)then
           if(is_star(typep(ipart))) then
              npart2=npart2+1
              tot_m=tot_m+mp(ipart)
              ave_u=ave_u+mp(ipart)*vp(ipart,1)
              ave_v=ave_v+mp(ipart)*vp(ipart,2)
              ave_w=ave_w+mp(ipart)*vp(ipart,3)
              sig_u=sig_u+mp(ipart)*vp(ipart,1)**2
              sig_v=sig_v+mp(ipart)*vp(ipart,2)**2
              sig_w=sig_w+mp(ipart)*vp(ipart,3)**2
           endif
           ipart=next_part  ! Go to next particle
        end do
     endif

     ! Normalize velocity moments
     if(npart2.gt.0)then
        ave_u=ave_u/tot_m
        ave_v=ave_v/tot_m
        ave_w=ave_w/tot_m
        sig_u=sqrt(max(sig_u/tot_m-ave_u**2, 0._dp))*scale_v/1d5
        sig_v=sqrt(max(sig_v/tot_m-ave_v**2, 0._dp))*scale_v/1d5
        sig_w=sqrt(max(sig_w/tot_m-ave_w**2, 0._dp))*scale_v/1d5
        str_d=tot_m/(2**ndim*vol_loc)*scale_nH
     endif

     ! Loop over cells
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        ! AGN formation sites
        ! if n_star>0.1 H/cc and v_disp>75 km/s
        if(str_d>0.1.and.MAX(sig_u,sig_v,sig_w)>20.)then
           flag2(ind_cell)=0
        else
           flag2(ind_cell)=1
        end if
     end do
  end do
  ! End loop over grids

#endif

111 format('   Entering quenching for level ',I2)

end subroutine quenching
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine growspin
  use amr_commons
  use pm_commons
  use hydro_commons
  use mpi_mod
  implicit none
  !----------------------------------------------------------------------
  ! Description: This subroutine computes spin values
  ! Yohan Dubois, May 29th, 2013
  !----------------------------------------------------------------------
  ! local constants
  integer::isink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::onethird,epsilon_r,nu_ratio,fourpi,prefact
  real(dp)::mbh,dmcoarse,mdisc,dm,chi,ax2,ay2,az2,jsinknorm,amod,amod2,mass_8,t_nu1
  real(dp)::rwarp,rsg,msg,dmacc,ax1,ay1,az1,bhspinnorm,theta
  real(dp)::ZZ1,ZZ2,r_lso,aac,jbh,jd,alpha,a01,chieps,yeartosec,dmini
  real(dp)::mbh2,dm2,dmdeplete,spinup

  if(.not. hydro)return
  if(ndim.ne.3)return
  if(nsink.eq.0)return

  if(verbose)write(*,*)'Entering growspin'

  yeartosec=3600d0*24d0*365d0
  onethird=1d0/3d0
  alpha=0.1d0 ! King Pringle Livio 2007 (alpha=0.1-0.4)
  a01=alpha/0.1d0
  nu_ratio=4d0*(1d0+7d0*alpha**2)/(4d0+alpha**2)/(2d0*alpha**2) ! Ogilvie (1999)
  fourpi=4d0*acos(-1d0)
  prefact=fourpi*6.67d-8*1.66d-24/(6.652d-25*3d10)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3

  do isink=1,nsink

     mdisc=0d0
     dm   =0d0
     chi=MIN(dMBHoverdt(isink)/dMEdoverdt(isink),1d0)
     !if(mad_jet)reduction_factor=MIN(chi/X_floor,1d0)
     !YD: accreted mass over dt(ilevel). Coarse is actually not coarse but correspond to dt(ilevel)
     dmcoarse=dMsmbh_all(isink)*scale_m
     dmini=dmcoarse
     !YD: BH mass before gas accretion over dt(ilevel)
     mbh=msink(isink)*scale_m-dmcoarse
     !YD: I assume epsilon_r constant over dt(ilevel) (can be updated in the end of the routine for Eddington rate)

     ax2=jsink(isink,1)
     ay2=jsink(isink,2)
     az2=jsink(isink,3)
     jsinknorm=SQRT(ax2**2+ay2**2+az2**2)
     ax2=ax2/jsinknorm
     ay2=ay2/jsinknorm
     az2=az2/jsinknorm

     do while(dmcoarse.gt.0d0)

        !------------ Compute the anti-alignement criterion ------------------------------------------------
        amod=ABS(spinmag(isink))
        if(amod.ne.0d0)then

           ZZ1=1d0+(1d0-spinmag(isink)**2)**onethird*((1d0+spinmag(isink))**onethird &
                & +(1d0-spinmag(isink))**onethird)
           ZZ2=SQRT(3d0*spinmag(isink)**2+ZZ1**2)
           if(spinmag(isink).ge.0d0)then
              r_lso=3d0+ZZ2-SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
           else
              r_lso=3d0+ZZ2+SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
           endif
           epsilon_r=1d0-sqrt(1d0-2d0/(3d0*r_lso))
           !if(mad_jet)epsilon_r=epsilon_r*reduction_factor
           chieps=chi/(epsilon_r/0.1d0)

           mass_8=mbh/(1d8*2d33)
           t_nu1=5.3d5*amod**0.875*chieps**(-0.75)*nu_ratio**(-0.875) & ! Fanidakis et al. (2011), equation (19)
                & *a01**(-1.5)*mass_8**1.375                            ! expressed in years
           rwarp=6.4d3*amod**0.625d0*mass_8**0.125*chieps**(-0.25)*nu_ratio**(-0.625)/sqrt(a01) ! equation (15)

           if(selfgrav)then
              rsg=5d2*a01**(28d0/45d0)*mass_8**(-52d0/45d0)*chieps**(-22d0/45d0) ! Dotti et al (2013)
              if(rsg.lt.rwarp)then
                 msg=6d5*a01**(-1d0/45d0)*mass_8**(34d0/45d0)*chieps**(4d0/45d0)*2d33 ! in g Dotti et al (2013)
                 mdisc=msg
              else
                 !dmacc=fourpi*6.67d-8*mbh*1.66d-24/(0.1d0*6.652d-25*3d10)*chieps ! in g/s
                 dmacc=prefact*mbh/epsilon_r*chi ! in g/s
                 mdisc=dmacc*(t_nu1*yeartosec) ! in g
              endif
           else
              dmacc=prefact*mbh/epsilon_r*chi ! in g/s
              mdisc=dmacc*(t_nu1*yeartosec) ! in g
           endif

           ax1=bhspin(isink,1)
           ay1=bhspin(isink,2)
           az1=bhspin(isink,3)
           bhspinnorm=SQRT(ax1**2+ay1**2+az1**2)
           ax1=ax1/bhspinnorm
           ay1=ay1/bhspinnorm
           az1=az1/bhspinnorm
           theta=ax1*ax2+ay1*ay2+az1*az2

           if(spinmag(isink).ge.maxspin.and.theta.ge.0d0)then
              ! if alignement and max spin: Add up disc (dmcoarse) and BH AM to get the new BH AM
              dm=dmcoarse
           else
              ! if not: Add up disc (dm) and BH AM to get the new BH AM
              ! Repeat the step until dmcoarse is consumed
              !dm=MIN(mdisc,dmcoarse)
              !dm=MAX(MIN(mdisc,dmcoarse),dmcoarse*1d-3)
              dm=MIN(MAX(mdisc,dmini*1d-3),dmcoarse) ! no more than 1000 iterations to update one spin
           endif
           if(selfgrav.and.(rsg.lt.rwarp))then
              aac=dm/mbh*SQRT(rsg  )/amod
           else
              aac=dm/mbh*SQRT(rwarp)/amod
           endif

           ! J_BH is expressed in units of jbh=amod*G*mbh**2/clight
           jbh=1d0
           jd=2d0*jbh*aac

           if(-aac>theta)then
              ! Anti-alignement occurs: BH AM unit vector is the opposite of the total AM (BH+D)
              ! and BH spin becomes negative
              bhspin(isink,1) = ax1*jbh + ax2*jd
              bhspin(isink,2) = ay1*jbh + ay2*jd
              bhspin(isink,3) = az1*jbh + az2*jd
              ax1=bhspin(isink,1)
              ay1=bhspin(isink,2)
              az1=bhspin(isink,3)
              bhspinnorm=SQRT(ax1**2+ay1**2+az1**2)
              ax1=ax1/bhspinnorm
              ay1=ay1/bhspinnorm
              az1=az1/bhspinnorm
              if(spinmag(isink).gt.0d0) spinmag(isink)  = - spinmag(isink)
           else
              ! Alignement occurs: BH AM unit vector is the total AM (BH+D)
              ! and BH spin becomes positive
              bhspin(isink,1) = ax1*jbh + ax2*jd
              bhspin(isink,2) = ay1*jbh + ay2*jd
              bhspin(isink,3) = az1*jbh + az2*jd
              if(spinmag(isink).lt.0d0) spinmag(isink)  = - spinmag(isink)
           endif

           amod2=SQRT(bhspin(isink,1)**2+bhspin(isink,2)**2+bhspin(isink,3)**2)
           if(amod2.ne.0d0)then
              bhspin(isink,1:3)=bhspin(isink,1:3)/amod2
           else
              bhspin(isink,1:3)=0d0
           endif

        else
           ! if BH has no spin, it gets the gas AM direction (no warp by definition)
           dm=dmcoarse
           bhspin(isink,1:3)=jsink(isink,1:3)
        endif
        !------------ Compute the anti-alignement criterion ------------------------------------------------

        if(mad_jet.and.(chi.lt.X_floor))then
           ! Fourth-order polynomial fit to the spinup parameters of McKinney et al, 2012
           dmdeplete=dm
           mbh2=mbh
           do while(dmdeplete.gt.0d0)
              spinup=0.97166d0-12.0026d0*spinmag(isink)-4.04337d0*spinmag(isink)**2 &
                   & +5.81317*spinmag(isink)**3+2.50482*spinmag(isink)**4
              dm2=MIN(1d-2*mbh2/(abs(spinup)+1d-2),dmdeplete)
              spinmag(isink)=spinmag(isink)+spinup*dm2/mbh2
              mbh2=mbh2+dm2
              dmdeplete=dmdeplete-dm2
           enddo
        else
           !------------ Compute the radius of last stable orbit --------
           ZZ1=1d0+(1d0-spinmag(isink)**2)**onethird*((1d0+spinmag(isink))**onethird &
                & +(1d0-spinmag(isink))**onethird)
           ZZ2=SQRT(3d0*spinmag(isink)**2+ZZ1**2)
           ! spinmag is a signed amplitude
           ! depending on co(counter)-rotation with the surrounding gas
           if(spinmag(isink).ge.0d0)then
              r_lso=3d0+ZZ2-SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
           else
              r_lso=3d0+ZZ2+SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
           endif
           !------------ Compute the radius of last stable orbit --------

           !------------ Update spin value ------------------------------
           if((mbh+dm)/mbh.le.SQRT(r_lso))then
              spinmag(isink)=MIN(onethird*SQRT(r_lso)*mbh/(mbh+dm) &
                   & *(4d0-SQRT(3d0*r_lso*(mbh/(mbh+dm))**2-2d0)),maxspin)
           else
              spinmag(isink)=maxspin
           endif
           !------------ Update spin value ------------------------------
        endif

        mbh=mbh+dm
        dmcoarse=dmcoarse-dm
     enddo
     ! End of loop over dm

     ZZ1=1d0+(1d0-spinmag(isink)**2)**onethird*((1d0+spinmag(isink))**onethird &
          & +(1d0-spinmag(isink))**onethird)
     ZZ2=SQRT(3d0*spinmag(isink)**2+ZZ1**2)
     if(spinmag(isink).ge.0d0)then
        r_lso=3d0+ZZ2-SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
     else
        r_lso=3d0+ZZ2+SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
     endif
     eps_sink(isink)=1d0-sqrt(1d0-2d0/(3d0*r_lso))
  enddo

end subroutine growspin
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine AGN_feedback
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons, ONLY: cic_levelmax
  use cooling_module, only: twopi
  use sink_particle_tracer, only : prepare_MC_tracer_to_jet
  ! AGNRT
#ifdef RT
  use rt_parameters,only: rt_AGN, nGroups
  use sed_module,only: getAGNfraclum
#endif
  !/AGNRT
  use mpi_mod

  implicit none
  !----------------------------------------------------------------------
  ! Description: This subroutine checks AGN events in cells where a
  ! sink particle lies.
  ! Yohan Dubois, December 15th, 2010
  !----------------------------------------------------------------------
  ! local constants
  integer::nAGN,iAGN,ilevel,ivar,info
  logical ,dimension(:),allocatable::ok_blast_agn
  integer ,dimension(:),allocatable::ind_blast,iAGN_myid
  real(dp),dimension(:),allocatable::mAGN,dAGNcell,vol_gas,mass_gas,vol_blast,mass_blast,psy_norm
  real(dp),dimension(:),allocatable::dMBH_AGN,dMEd_AGN,dMsmbh_AGN,EsaveAGN,Msmbh,spinmagAGN,eps_AGN,X_radio
  real(dp),dimension(:,:),allocatable::xAGN,vAGN,jAGN
  real(dp),dimension(:,:),allocatable::passiveAGN
  real(dp)::temp_blast,volume
#ifdef RT
  ! AGNRT
  real(dp)::Erad_AGN, epsilon_r
  real(dp)::scale_vol,scale_evtocode
  real(dp)::scale_Np,scale_Fp
  real(dp)::scale_m
  !/AGNRT
#endif
  integer,dimension(1:nsink)::itemp
  integer::isort,isink,idim,ilun,ii,itype
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,Mfrac,ttsta,ttend
  character(LEN=5)::nchar
  character(LEN=80)::filename
  real(dp),allocatable,dimension(:)::xdp
  integer  :: nx_loc
  real(dp) :: dx_min, rmax, r2_cloud, vol_cloud
  real(dp) :: dx_dm, rmax_dm, r2_dm, vol_dm, onepi, scale

  if(.not. hydro)return
  if(ndim.ne.3)return
  if(nsink.eq.0)return

!!$  call MPI_BARRIER(MPI_COMM_WORLD,info)
  if(myid.eq.1)ttsta=MPI_WTIME()

  if(verbose)write(*,*)'Entering make_AGN'

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! HP: dynamical friction from particles
  if (drag_part) then
     onepi = twopi /2.d0
     nx_loc=(icoarse_max-icoarse_min+1)
     scale=boxlen/dble(nx_loc)
     dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp
     rmax = DF_ncells*dx_min
     r2_cloud = rmax**2
     vol_cloud = 4d0/3d0*onepi*(rmax)**3

     if (cic_levelmax.gt.0) then
        dx_dm = scale*0.5d0**cic_levelmax/aexp
     else
        dx_dm = dx_min
     end if
     rmax_dm = DF_ncells*dx_dm
     r2_dm = rmax_dm**2
     vol_dm = 4d0/3d0*onepi*rmax_dm**3
  end if
  !/HP: dynamical friction from particles

  ! AGNRT
#ifdef RT
  if(rt_AGN) then
     scale_m=scale_d*scale_l**3d0
     call rt_units(scale_Np, scale_Fp)
     scale_vol=scale_l**ndim
     ! 1.60217646d-12 erg = 1 eV
     scale_evtocode=1.60217646d-12/(scale_d*scale_l**5/scale_t**2)
     !if (myid==1)write(*,*)'UNITS scalevol, scaleevtocode, np', scale_vol, scale_evtocode, scale_Np
  endif
#endif
  !/AGNRT


  if(myid==1.and.nsink>0.and.sinkprops.and. (.not.finestep_AGN))then
     if(.not.(nstep_coarse==nstep_coarse_old.and.nstep_coarse>0))then
     call title(nstep_coarse,nchar)
     filename='sink_'//TRIM(nchar)//'.dat'
     ilun=ncpu*4+11
     open(unit=ilun,file=TRIM(filename),form='unformatted')
     write(ilun)nsink     ! Number of sink
     write(ilun)ndim      ! Number of dimensions
     if (cosmo) then
        write(ilun)aexp      ! expansion factor
     else
        write(ilun)t
     end if
     write(ilun)scale_l   ! length scale
     write(ilun)scale_d   ! density scale
     write(ilun)scale_t   ! time scale
     allocate(xdp(1:nsink))
     write(ilun)idsink(1:nsink) ! Identities
     write(ilun)msink (1:nsink) ! Masses
     do idim=1,ndim
        xdp(1:nsink)=xsink(1:nsink,idim)
        write(ilun)xdp    ! Positions
     enddo
     do idim=1,ndim
        xdp(1:nsink)=vsink(1:nsink,idim)
        write(ilun)xdp    ! Velocities
     enddo
     do idim=1,ndim
        xdp(1:nsink)=jsink(1:nsink,idim)
        write(ilun)xdp    ! gas AM
     enddo
     write(ilun)dMBHoverdt(1:nsink) ! Bondi accretion rate
     write(ilun)dMEdoverdt(1:nsink) ! Eddington accretion rate
     write(ilun)dMsmbh    (1:nsink) ! Total accreted mass
     write(ilun)d_avgptr  (1:nsink) ! Mean gas density
     write(ilun)c_avgptr  (1:nsink) ! Mean sound speed
     write(ilun)v_avgptr  (1:nsink) ! Relative BH-gas velocity
     write(ilun)Esave     (1:nsink) ! Energy saved from last coarse step
     do idim=1,ndim
        xdp(1:nsink)=bhspin(1:nsink,idim)
        write(ilun)xdp    ! BH AM
     enddo
     write(ilun)spinmag   (1:nsink) ! BH spin parameter
     write(ilun)eps_sink  (1:nsink) ! BH radiative efficiency
#ifdef RT
     if (rt_AGN) then
        write(ilun)LAGN_coarse(1:nsink) ! Luminosity
     end if
#endif
     ! HP: dynamical friction from particles
     if (drag_part) then
        v_background(1:nsink, 1:ndim, 1:2) = 0d0
        m_background(1:nsink, 1:2) = tiny(0d0)
        n_background(1:nsink, 1:2) = 0
        mass_lowspeed_background(1:nsink, 1:2) = 0d0
        fact_fast_background(1:nsink, 1:2) = 0d0

        do ii = levelmin, nlevelmax
           do isink = 1,nsink
              do itype = 1, 2 !1 for stars - 2 for DM
                 do idim = 1, ndim
                    v_background(isink, idim, itype) = &
                         v_background(isink, idim, itype) * m_background(isink, itype) + &
                         v_DF(isink, ii, idim, itype) * mass_DF(isink, ii, itype)
                 end do
                 m_background(isink, itype) = m_background(isink, itype) + mass_DF(isink, ii, itype)
                 do idim = 1, ndim
                    v_background(isink, idim, itype) = v_background(isink, idim, itype) / m_background(isink, itype)
                 end do
                 n_background(isink, itype) = n_background(isink, itype) + n_part(isink, ii, itype)
                 mass_lowspeed_background(isink, itype) = mass_lowspeed_background(isink, itype) + mass_lowspeed(isink, ii, itype)
                 fact_fast_background(isink, itype) = fact_fast_background(isink, itype) + fact_fast(isink, ii, itype)
              end do
           end do
        end do

        do itype = 1, 2
           if (itype.eq.1) then
              volume = vol_cloud
           else
              volume = vol_dm
           end if
           write(ilun) m_background(1:nsink, itype) / volume ! star/DM density around BH
        end do

        do itype = 1, 2
           do idim= 1, ndim
              write(ilun) v_background(1:nsink, idim, itype) ! star/DM velocity
           end do
        end do

        do itype = 1, 2
           write(ilun) n_background(1:nsink, itype) ! number of stars/DM
        end do

        do itype = 1, 2
           if (itype.eq.1) then
              volume = vol_cloud
           else
              volume = vol_dm
           end if
           !density of slow moving stars/DM
           write(ilun) mass_lowspeed_background(1:nsink, itype) / volume
        end do

        do itype = 1, 2
           if (itype.eq.1) then
              volume = vol_cloud
           else
              volume = vol_dm
           end if
           !factor for fast moving stars/DM (see Chandrasekar+43)
           write(ilun) fact_fast_background(1:nsink, itype) / volume

        end do
     end if
     !/HP: dynamical friction from particles
     write(ilun)t         ! Simulation time
     close(ilun)
     deallocate(xdp)
     endif
  endif

  ! Get only AGN that are in my CPU or close to the border of my cpu (huge speed-up!)
  call getAGNonmyid(itemp,nAGN)

  ! Allocate the arrays for the position and the mass of the AGN in myid
  allocate(xAGN(1:nAGN,1:3),vAGN(1:nAGN,1:3),mAGN(1:nAGN),dAGNcell(1:nAGN),iAGN_myid(1:nAGN) &
       & ,dMBH_AGN(1:nAGN),dMEd_AGN(1:nAGN),dMsmbh_AGN(1:nAGN),Msmbh(1:nAGN),EsaveAGN(1:nAGN),jAGN(1:nAGN,1:3) &
       & ,spinmagAGN(1:nAGN),eps_AGN(1:nAGN), X_radio(1:nAGN))
  allocate(passiveAGN(1:nAGN,imetal:nvar))
  xAGN=0; vAGN=0; mAGN=0d0; dAGNcell=0; iAGN_myid=0; dMBH_AGN=0d0; dMEd_AGN=0d0; EsaveAGN=0d0; jAGN=0d0; spinmagAGN=0; eps_AGN=0; X_radio(1:nAGN)=0; passiveAGN = 0

  ! AGNRT
  ! Compute the radiation to release over the next coarse timestep (for all sinks?)
#ifdef RT
  if (rt_AGN) then
     ! Total radiation from accretion over the last coarse timestep
     !   -> Etot = epsilon_r*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
     ! Technically, dMsmbh_AGN is the accreted mass, so (1-epsilon_r) * Mdot * dt
     !   -> Etot ~ Etot/(1-epsilon_r)
     ! Normalization for later
     !   -> Etot_code = Etot / (scale_evtocode) / (scale_vol * scale_np)
     ! Energy -> luminosity
     !   -> Lcoarse = Etot_code / dtnew(levelmin)
     ! Need to change rtuold by !!!(the last vol_loc is to get the number *density* of photons
     !   -> Np_inj = Lcoarse*dtnew(ilevel) * group_agnegy_frac(grp)/group_egy(mygroup) * weight/volume !!/vol_loc
     do isink=1, nsink
        if(dMBH_coarse(isink)/dMEd_coarse(isink) .ge. X_floor) then
           if(spin_bh)then
              epsilon_r=eps_sink(isink)
           else
              epsilon_r=0.1d0
           endif

           Erad_agn = epsilon_r*dMeff_coarse(isink) * (3d10/scale_v)**2d0 / (1-epsilon_r)
           Erad_agn = Erad_agn  / (scale_evtocode) / (scale_vol * scale_np)
           LAGN_coarse(isink) = Erad_agn / dtnew(levelmin)

           ! Compute the group_agnegy_frac for each sink
           call getAGNfraclum(msink(isink)*scale_m/2d33, &
                dMBH_coarse(isink)/dMEd_coarse(isink), &
                lumfrac_AGN(isink, :))
           if(verbose .and. myid==1) write(*,*)'Luminosity fraction for sink #', isink, lumfrac_AGN(isink,:)

        else
           LAGN_coarse(isink) = 0.d0
        end if
     end do
  end if
#endif
  !/AGNRT



  do iAGN=1,nAGN
     isort=itemp(iAGN)
     iAGN_myid(iAGN)=isort
     xAGN(iAGN,1)=xsink(isort,1)
     xAGN(iAGN,2)=xsink(isort,2)
     xAGN(iAGN,3)=xsink(isort,3)
     vAGN(iAGN,1)=vsink(isort,1)
     vAGN(iAGN,2)=vsink(isort,2)
     vAGN(iAGN,3)=vsink(isort,3)

     if(spin_bh)then
        jAGN(iAGN,1)=bhspin(isort,1)
        jAGN(iAGN,2)=bhspin(isort,2)
        jAGN(iAGN,3)=bhspin(isort,3)
     else
        jAGN(iAGN,1)=jsink(isort,1)
        jAGN(iAGN,2)=jsink(isort,2)
        jAGN(iAGN,3)=jsink(isort,3)
     endif
     Msmbh     (iAGN)=msink(isort)
     dMBH_AGN  (iAGN)=dMBH_coarse(isort)
     dMEd_AGN  (iAGN)=dMEd_coarse(isort)
     dMsmbh_AGN(iAGN)=dMsmbh(isort)
     EsaveAGN  (iAGN)=Esave (isort)
     spinmagAGN(iAGN)=spinmag(isort)
     eps_AGN   (iAGN)=eps_sink(isort)
  enddo

  ! Allocate arrays that are outputs of average_AGN (initialised in average_AGN)
  allocate(vol_gas(1:nAGN),mass_gas(1:nAGN),psy_norm(1:nAGN),vol_blast(1:nAGN),mass_blast(1:nAGN))
  allocate(ind_blast(1:nAGN),ok_blast_agn(1:nAGN))

  ! Check if AGN goes into jet mode
  ok_blast_agn(1:nAGN)=.false.
  do iAGN=1,nAGN
     if (dMEd_AGN(iAGN) == 0) then
        X_radio(iAGN) = huge(X_radio(iAGN))
     else
        X_radio(iAGN) = dMBH_AGN(iAGN)/dMEd_AGN(iAGN)
     end if

     if(X_radio(iAGN).lt.X_floor .and. EsaveAGN(iAGN).eq.0d0)then
        Mfrac=dMsmbh_AGN(iAGN)/(Msmbh(iAGN)-dMsmbh_AGN(iAGN))
        if(Mfrac.ge.jetfrac)ok_blast_agn(iAGN)=.true.
     endif
  enddo

  ! Compute some averaged quantities before doing the AGN energy input
  call average_AGN(xAGN,dMBH_AGN,dMEd_AGN,mAGN,dAGNcell,passiveAGN,jAGN,vol_gas,mass_gas,psy_norm,vol_blast &
       & ,mass_blast,ind_blast,nAGN,iAGN_myid,ok_blast_agn,EsaveAGN,X_radio)

  ! Check if AGN goes into thermal blast wave mode
  do iAGN=1,nAGN
     if(X_radio(iAGN) .ge. X_floor .and. EsaveAGN(iAGN).eq.0d0)then
        ! Compute estimated average temperature in the blast
        temp_blast=0.0
        if(vol_gas(iAGN)>0.0)then
           temp_blast=eAGN_T*1d12*dMsmbh_AGN(iAGN)/mass_gas(iAGN)
        else
           if(ind_blast(iAGN)>0)then
              temp_blast=eAGN_T*1d12*dMsmbh_AGN(iAGN)/mass_blast(iAGN)
           endif
        endif
        if(temp_blast>TAGN)then
           ok_blast_agn(iAGN)=.true.
        endif
     endif
  end do


  ! Modify hydro quantities to account for a Sedov blast wave
  call AGN_blast(xAGN, vAGN, dMsmbh_AGN, dMBH_AGN, dMEd_AGN, mAGN,    &
       dAGNcell, passiveAGN, jAGN, ind_blast, vol_gas, psy_norm, vol_blast, &
       nAGN, iAGN_myid, ok_blast_agn, EsaveAGN, spinmagAGN, eps_AGN, X_radio)

  ! Store AGN quantities for later use for Monte Carlo tracers
  if (MC_tracer) then
     call prepare_MC_tracer_to_jet(xAGN, jAGN, mAGN, dAGNcell, &
          X_radio, ind_blast, nAGN)
  end if

  ! Reset total accreted mass if AGN input has been done
  do iAGN=1,nAGN
     if(ok_blast_agn(iAGN))then
        isort=iAGN_myid(iAGN)
        dMsmbh(isort)=0d0
     endif
  end do
  ! Important: initialise coarse Bondi and Eddington mass for the next coarse step
#ifndef WITHOUTMPI
  dMsmbh_new=dMsmbh
  call MPI_ALLREDUCE(dMsmbh_new,dMsmbh_all,nsink,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
  dMsmbh=dMsmbh_all
#endif
  dMBH_coarse=0d0; dMEd_coarse=0d0
  ! AGNRT
#ifdef RT
  if(rt_AGN) dMeff_coarse = 0.d0
#endif
  !/AGNRT

  ! Important: save the energy for the next time step that has not been put in the current time step
#ifndef WITHOUTMPI
  Esave_new=0d0
  do iAGN=1,nAGN
     isink=iAGN_myid(iAGN)
     Esave_new(isink)=EsaveAGN(iAGN)
  enddo
  call MPI_ALLREDUCE(Esave_new,Esave_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(Efeed_new,Efeed_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  Esave= Esave_all
  Efeed= Efeed_all
#else
  Esave=EsaveAGN
#endif
  Efeed_new=0d0 !Reset for next coarse step

  ! Deallocate everything
  deallocate(vol_gas,mass_gas,psy_norm,vol_blast,mass_blast)
  deallocate(ind_blast)
  deallocate(xAGN,vAGN,mAGN,jAGN)
  deallocate(passiveAGN)
  deallocate(iAGN_myid,ok_blast_agn,dMBH_AGN,dMEd_AGN,dMsmbh_AGN,Msmbh,EsaveAGN)
  deallocate(spinmagAGN,eps_AGN)

  ! Update hydro quantities for split cells
  do ilevel=nlevelmax,levelmin,-1
     call upload_fine(ilevel)
#ifdef SOLVERmhd
     do ivar=1,nvar+3
#else
     do ivar=1,nvar
#endif
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     enddo

  enddo

  if (myid.eq.1 .and. (.not.finestep_AGN)) then
     ttend=MPI_WTIME()
     write(*,*) ' Time elapsed in AGN_feedback [sec]:', sngl(ttend-ttsta)
  endif

end subroutine AGN_feedback
!################################################################
!################################################################
!################################################################
!################################################################
subroutine average_AGN(xAGN,dMBH_AGN,dMEd_AGN,mAGN,dAGNcell,passiveAGN,jAGN,vol_gas,mass_gas,psy_norm,vol_blast &
     & ,mass_blast,ind_blast,nAGN,iAGN_myid,ok_blast_agn,EsaveAGN,X_radio)
  use pm_commons
  use amr_commons
  use hydro_commons

  use utils, only : distance3d, AGN_integration, AGN_integration_depth, AGN_VOLUME_INTEGRATION_bisect, AGN_VOLUME_INTEGRATION_MC

  use geom_types, only : Capsule_t, Box_t
  use geom_distances
  use geom_volumes

  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine average the hydro quantities inside the AGN bubble
  ! Jet case    : remove some gas from the BH's cell
  ! Thermal case: get the mass of gas in the AGN bubble
  ! Yohan Dubois, December 15th, 2010
  !------------------------------------------------------------------------
  integer, intent(in) :: nAGN
  integer, intent(in), dimension(1:nAGN) :: iAGN_myid
  real(dp), intent(in), dimension(1:nAGN) :: X_radio, EsaveAGN
  real(dp), dimension(1:nAGN, 1:3), intent(in) :: jAGN, xAGN, dMBH_AGN, dMEd_AGN
  ! real(dp), dimension(1:nAGN, 1:3), intent(out) :: xAGN
  real(dp), dimension(1:nAGN), intent(out) :: mAGN, dAGNcell
  real(dp), dimension(1:nAGN, imetal:nvar), intent(out) :: passiveAGN
  real(dp), dimension(1:nAGN), intent(out) :: vol_gas, mass_gas, psy_norm, vol_blast, mass_blast

  integer::ilevel,ncache,iAGN,ind,ix,iy,iz,ngrid,iskip,ivar
  integer::i,isink,nx_loc,igrid,info
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dr_AGN,d,u,v,w,dr_cell
  real(dp)::scale,dx,dxx,dyy,dzz,dx_min,dx_loc,vol_loc,rmax2,rmax
  real(dp)::x_half,y_half,z_half,x_box,y_box,z_box
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nAGN)::ind_blast

#ifndef WITHOUTMPI
  real(dp),dimension(1:nsink)::vol_gas_mpi,mass_gas_mpi,mAGN_mpi,psy_norm_mpi
  real(dp),dimension(1:nsink)::vol_gas_all,mass_gas_all,mAGN_all,psy_norm_all
  real(dp),dimension(1:nsink, imetal:nvar)::passiveAGN_mpi, passiveAGN_all
#endif
  logical ,dimension(1:nAGN)::ok_blast_agn
  logical ,dimension(1:nvector),save::ok
  real(dp)::jtot,j_x,j_y,j_z,psy, tmp_volume
  real(dp)::eint,ekk, d2, drjet, dzjet

  type(Capsule_t) :: capsule
  type(Box_t) :: box

  ! Configure geometry calculator
  ! 3 level of refinement (2**3 cells in each direction)
  call set_depth(AGN_integration_depth)
  call set_ndraw(100)

  ! Set box direction (we'll work in box's frame)
  box%u = [1, 0, 0]
  box%v = [0, 1, 0]
  box%w = [0, 0, 1]

  if(verbose)write(*,*)'Entering average_AGN'

  ! Mesh spacing in that level
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)
  x_half=scale*xbound(1)/2.0; y_half=scale*xbound(2)/2.0; z_half=scale*xbound(3)/2.0
  x_box =scale*xbound(1); y_box =scale*xbound(2); z_box =scale*xbound(3)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  rmax=MAX(1d0*dx_min*scale_l/aexp,rAGN*3.08d21)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  ! Initialize the averaged variables
  vol_gas=0d0;mass_gas=0d0;vol_blast=0d0;mass_blast=0d0;ind_blast=-1;psy_norm=0d0
  passiveAGN = 0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale

                 box%origin = [x, y, z]
                 box%extents = [dx_loc, dx_loc, dx_loc] / 2._dp ! This is half the width of the box

                 do iAGN=1,nAGN
                    ! Here we ignore periodicity.
                    ! This leads to a huge speedup, since we only need
                    ! to get the AGN within the CPU or near the CPU
                    ! boundary. This can be done by a direct lookup at
                    ! the Hilbert curve. Close to the boundary, it's
                    ! more complicated.  It is also worth noting that
                    ! this *conserves* the deposited mass as well as
                    ! the energy, etc... the deposition is just
                    ! asymetric
                    call distance3d(x, y, z, xAGN(iAGN, 1), xAGN(iAGN, 2), xAGN(iAGN, 3), &
                         dxx, dyy, dzz, ignore_periodicity=.true.)
                    ! dxx=x-xAGN(iAGN,1)
                    ! dyy=y-xAGN(iAGN,2)
                    ! dzz=z-xAGN(iAGN,3)
                    dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz

                    ! ------------------------------------------
                    ! case 0: Some energy has not been released
                    ! ------------------------------------------
                    if(EsaveAGN(iAGN).gt.0d0)then
                       if(dr_AGN.le.rmax2)then
                          vol_gas (iAGN)=vol_gas (iAGN)+vol_loc*uold(ind_cell(i),1)
                       endif
                       dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                       if(dr_cell.le.dx_loc/2.0)then
                          ind_blast (iAGN)=ind_cell(i)
                          vol_blast (iAGN)=vol_loc
                       endif

                    ! ------------------------------------------
                    ! case 1: All energy has been released in
                    ! previous time step -> choose between kinetic
                    ! or thermal feedback
                    ! ------------------------------------------
                    else

                       ! ------------------------------------------
                       ! Jet feedback
                       ! ------------------------------------------
                       if(X_radio(iAGN).lt.X_floor)then
                          if(ok_blast_agn(iAGN))then
                             jtot=norm2(jAGN(iAGN, 1:ndim))
                             if (jtot > 0d0) then
                                j_x=jAGN(iAGN,1)/jtot
                                j_y=jAGN(iAGN,2)/jtot
                                j_z=jAGN(iAGN,3)/jtot
                                dr_cell=max(abs(dxx),abs(dyy),abs(dzz))

                                dzjet= dxx*j_x + dyy*j_y + dzz*j_z
                                
                                ! Work in simulation's frame
                                capsule%segment%start = xAGN(iAGN, :) - rmax * jAGN(iAGN, :) / jtot
                                capsule%segment%end = xAGN(iAGN, :) + rmax * jAGN(iAGN, :) / jtot
                                capsule%r = rmax

                                ! Compute distance between cell and capsule containing jet
                                call LineSegmentBoxDistanceSquared(box, capsule%segment, d2)
                                ! call CapsuleBoxDistanceSquared(box, capsule, d2)

                                ! Check if the cell lies within the AGN jet cylindre
                                if (d2 <= rmax2)then
                                   if (AGN_integration == AGN_VOLUME_INTEGRATION_MC) then
                                      call CapsuleBoxIntegrateMC(psy_function_loc, &
                                           box, capsule, tmp_volume, psy)
                                      if (tmp_volume > 0) then
                                         psy = psy / tmp_volume
                                      else
                                         psy = 1d-10
                                      end if
                                      vol_gas(iAGN) = vol_gas(iAGN) + tmp_volume
                                   else if (AGN_integration == AGN_VOLUME_INTEGRATION_bisect) then
                                      call CapsuleBoxIntegrate(psy_function_loc, &
                                           box, capsule, tmp_volume, psy)
                                      if (tmp_volume > 0) then
                                         psy = psy / tmp_volume
                                      else
                                         psy = 1d-10
                                      end if
                                      vol_gas(iAGN) = vol_gas(iAGN) + tmp_volume
                                   else
                                      vol_gas(iAGN)=vol_gas(iAGN)+vol_loc
                                      drjet=sqrt(max(dr_AGN-dzjet*dzjet, 0._dp))
                                      psy = exp(-drjet**2/2d0/rmax2)
                                   end if
                                   psy_norm(iAGN)=psy_norm(iAGN)+psy*vol_loc
                                endif

                                if(dr_cell <= dx_loc/2.0)then
                                   ind_blast(iAGN)=ind_cell(i)
                                   d=max(uold(ind_blast(iAGN),1), smallr)
                                   u=uold(ind_blast(iAGN),2)/d
                                   v=uold(ind_blast(iAGN),3)/d
                                   w=uold(ind_blast(iAGN),4)/d
                                   ekk=0.5d0*d*(u*u+v*v+w*w)
                                   eint=uold(ind_blast(iAGN),5)-ekk
                                   vol_blast  (iAGN)=vol_loc
                                   mAGN(iAGN)=min(mloadAGN*dMsmbh(iAGN),0.25d0*d*vol_loc)
                                   dAGNcell(iAGN) = d * vol_loc

                                   do ivar = imetal, nvar
                                      passiveAGN(iAGN, ivar) = uold(ind_blast(iAGN), ivar)/d
                                      uold(ind_blast(iAGN), ivar) = uold(ind_blast(iAGN), ivar) &
                                           & - passiveAGN(iAGN, ivar)*mAGN(iAGN)/vol_loc
                                   end do

                                   d=max(uold(ind_blast(iAGN),1)-mAGN(iAGN)/vol_loc, smallr)
                                   uold(ind_blast(iAGN),1)=d
                                   uold(ind_blast(iAGN),2)=d*u
                                   uold(ind_blast(iAGN),3)=d*v
                                   uold(ind_blast(iAGN),4)=d*w
                                   uold(ind_blast(iAGN),5)=eint+0.5d0*d*(u*u+v*v+w*w)
                                endif
                                ! If no spin for the jet then put all thermal
                             else
                                if(dr_AGN <= rmax2) then
                                   vol_gas(iAGN)=vol_gas(iAGN)+vol_loc
                                endif
                                dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                             endif
                          endif

                       ! ------------------------------------------
                       ! Thermal feedback
                       ! ------------------------------------------
                       else
                          if(dr_AGN.le.rmax2)then
                             vol_gas (iAGN)=vol_gas (iAGN)+vol_loc*uold(ind_cell(i),1)
                             mass_gas(iAGN)=mass_gas(iAGN)+vol_loc*uold(ind_cell(i),1)
                          endif
                          dr_cell=MAX(ABS(dxx),ABS(dyy),ABS(dzz))
                          if(dr_cell.le.dx_loc/2.0)then
                             ind_blast (iAGN)=ind_cell(i)
                             vol_blast (iAGN)=vol_loc
                             mass_blast(iAGN)=vol_loc*uold(ind_cell(i),1)
                          endif
                       endif
                    endif

                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  !################################################################
#ifndef WITHOUTMPI
  vol_gas_mpi=0d0; mass_gas_mpi=0d0; mAGN_mpi=0d0; psy_norm_mpi=0d0
  passiveAGN_mpi = 0._dp
  ! Put the nAGN size arrays into nsink size arrays to synchronize processors
  do iAGN=1,nAGN
     isink=iAGN_myid(iAGN)
     vol_gas_mpi (isink)=vol_gas (iAGN)
     mass_gas_mpi(isink)=mass_gas(iAGN)
     mAGN_mpi    (isink)=mAGN    (iAGN)
     do ivar = imetal, nvar
        passiveAGN_mpi(isink, ivar) = passiveAGN(iAGN, ivar)
     end do
     psy_norm_mpi(isink)=psy_norm(iAGN)
  enddo
  call MPI_ALLREDUCE(vol_gas_mpi ,vol_gas_all ,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mass_gas_mpi,mass_gas_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mAGN_mpi    ,mAGN_all    ,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(passiveAGN_mpi, passiveAGN_all,nsink*(nvar-imetal+1), &
       & MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(psy_norm_mpi,psy_norm_all,nsink,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  vol_gas_mpi =vol_gas_all
  mass_gas_mpi=mass_gas_all
  mAGN_mpi    =mAGN_all
  passiveAGN_mpi = passiveAGN_all
  psy_norm_mpi=psy_norm_all
  ! Put the nsink size arrays into nAGN size arrays
  do iAGN=1,nAGN
     isink=iAGN_myid(iAGN)
     vol_gas (iAGN)=vol_gas_mpi (isink)
     mass_gas(iAGN)=mass_gas_mpi(isink)
     mAGN    (iAGN)=mAGN_mpi    (isink)
     do ivar = imetal, nvar
        passiveAGN(iAGN, ivar) = passiveAGN_mpi(isink, ivar)
     end do
     psy_norm(iAGN)=psy_norm_mpi(isink)
  enddo
#endif
  !################################################################

  if(verbose)write(*,*)'Exiting average_AGN'
contains
  real(dp) function psy_function_loc(X) result(psy)
    ! A note on this function: it is called in the frame of the
    ! capsule, in which the last dimension is // to
    real(dp), intent(in) :: X(3)

    call psy_function(X, rmax2, psy)

  end function psy_function_loc
end subroutine average_AGN

!################################################################
!################################################################
!################################################################
!################################################################
subroutine AGN_blast(xAGN,vAGN,dMsmbh_AGN,dMBH_AGN,dMEd_AGN,mAGN,dAGNcell,passiveAGN,jAGN,ind_blast,vol_gas &
     & ,psy_norm,vol_blast,nAGN,iAGN_myid,ok_blast_AGN,EsaveAGN,spinmagAGN,eps_AGN, X_radio)
  use pm_commons
  use hydro_commons
  use cooling_module, only: XH=>X, rhoc, mH

  use utils, only : AGN_integration, AGN_integration_depth, &
       AGN_VOLUME_INTEGRATION_MC, AGN_VOLUME_INTEGRATION_bisect, &
       draw_normal, pi
  use geom_types, only : Capsule_t, Box_t
  use geom_distances
  use geom_volumes
  use amr_commons

  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine do the AGN energy inupt
  ! Jet case: do a kinetic jet solution as in Dubois et al. (2010)
  ! by depositing mass, momentum and energy into a small cylindre around the BH
  ! Thermal case: do a thermal energy input as in Teyssier et al. (2010)
  ! by depositing internal energy into small bubble
  ! Yohan Dubois, December 15th, 2010
  !------------------------------------------------------------------------
  integer, intent(in) :: nAGN
  integer, dimension(1:nAGN), intent(in) :: iAGN_myid, ind_blast
  logical, dimension(1:nAGN), intent(in) :: ok_blast_AGN
  real(dp),dimension(1:nAGN),intent(in)::mAGN,dAGNcell,vol_blast
  real(dp),dimension(1:nAGN, imetal:nvar), intent(in)::passiveAGN
  real(dp),dimension(1:nAGN,1:3), intent(in)::xAGN,vAGN,jAGN
  real(dp),dimension(1:nAGN), intent(in)::psy_norm,dMsmbh_AGN,dMBH_AGN,dMEd_AGN,spinmagAGN,eps_AGN,X_radio

  real(dp), dimension(1:nAGN), intent(inout) :: EsaveAGN, vol_gas

  integer::ilevel,iAGN,ind,ix,iy,iz,ngrid,iskip
  integer::i,nx_loc,igrid,ncache
  integer::ivar
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::x,y,z,dx,dxx,dyy,dzz,dr_AGN,d,u,v,w,d_gas
  real(dp)::scale,dx_min,dx_loc,vol_loc,rmax2,rmax,x_half,y_half,z_half,x_box,y_box,z_box
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp),dimension(1:3)::xbound,skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  logical ,dimension(1:nvector),save::ok
  logical, dimension(1:nAGN) :: ok_save
  real(dp), dimension(1:nAGN) :: p_gas, EAGN, uBlast
  real(dp)::jtot,j_x,j_y,j_z,psy,nCOM,T2_1,T2_2,ekk,eint,etot
  real(dp)::ekkold,T2maxAGNz,epsilon_r,onethird,eff_mad
  integer::idim

  type(Capsule_t) :: capsule
  type(Box_t) :: box
  real(dp) :: d2, dzjet, drjet, tmp_volume

  ! Configure geometry calculator
  ! 3 level of refinement (2**3 cells in each direction)
  call set_depth(AGN_integration_depth)
  call set_ndraw(100)

  ! Set box direction (we'll work in box's frame)
  box%u = [1, 0, 0]
  box%v = [0, 1, 0]
  box%w = [0, 0, 1]

  if(verbose)write(*,*)'Entering AGN_blast'

  onethird=1d0/3d0
  ! Mesh spacing in that level
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)
  x_half=scale*xbound(1)/2.0; y_half=scale*xbound(2)/2.0; z_half=scale*xbound(3)/2.0
  x_box =scale*xbound(1); y_box =scale*xbound(2); z_box =scale*xbound(3)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nCOM = 1d-1*omega_b*rhoc/aexp**3*XH/mH/scale_nH
  T2maxAGNz=T2maxAGN*aexp


  ! Maximum radius of the ejecta
  rmax=MAX(1d0*dx_min*scale_l/aexp,rAGN*3.08d21)
  rmax=rmax/scale_l
  rmax2=rmax*rmax

  uBlast=0d0
  ok_save=.false.
  do iAGN=1,nAGN
     if(EsaveAGN(iAGN).gt.0d0)then
        ok_save(iAGN)=.true.
        EAGN (iAGN)=EsaveAGN(iAGN)
        if (vol_gas(iAGN) .gt. 0d0) then !!!
           p_gas(iAGN)=EAGN    (iAGN) / vol_gas(iAGN)
        else
           p_gas(iAGN) = 0.0d0
        end if
     else if(ok_blast_agn(iAGN))then
        if(spin_bh)then
           epsilon_r=eps_AGN(iAGN)
        else
           epsilon_r=0.1d0
        endif
!!$        ! WARNING: Take care of spinmag on my CPU
!!$        ZZ1=1d0+(1d0-spinmagAGN(iAGN)**2)**onethird*((1d0+spinmagAGN(iAGN))**onethird &
!!$             & +(1d0-spinmagAGN(iAGN))**onethird)
!!$        ZZ2=SQRT(3d0*spinmagAGN(iAGN)**2+ZZ1**2)
!!$        if(spinmagAGN(iAGN).ge.0d0)then
!!$           r_lso=3d0+ZZ2-SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
!!$        else
!!$           r_lso=3d0+ZZ2+SQRT((3d0-ZZ1)*(3d0+ZZ1+2d0*ZZ2))
!!$        endif
!!$        epsilon_r=1d0-sqrt(1d0-2d0/(3d0*r_lso))
        if(X_radio(iAGN).lt.X_floor)then
           if(mad_jet)then
              ! Fourth-order polynomial fit to McKinney et al, 2012 (jet+wind)
              eff_mad=4.10507+0.328712*spinmagAGN(iAGN)+76.0849*spinmagAGN(iAGN)**2d0 &
                   & +47.9235*spinmagAGN(iAGN)**3d0+3.86634*spinmagAGN(iAGN)**4d0
              eff_mad=eff_mad/100d0
              EAGN(iAGN)=eff_mad*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
              p_gas(iAGN)=(1d0-f_ekAGN)*EAGN(iAGN) / vol_gas(iAGN)
              if(mAGN(iAGN).gt.0d0)uBlast(iAGN)=sqrt(2d0*f_ekAGN*EAGN(iAGN)/mAGN(iAGN))
           else
              EAGN(iAGN)=eAGN_K*epsilon_r*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
              p_gas(iAGN)=(1d0-f_ekAGN)*EAGN(iAGN) / vol_gas(iAGN)
              if(mAGN(iAGN).gt.0d0)uBlast(iAGN)=sqrt(2d0*f_ekAGN*EAGN(iAGN)/mAGN(iAGN))
           endif
        else
           EAGN  (iAGN)=eAGN_T*epsilon_r*dMsmbh_AGN(iAGN)*(3d10/scale_v)**2d0
           if (vol_gas(iAGN) .gt. 0d0) then
              p_gas(iAGN)=EAGN    (iAGN) / vol_gas(iAGN)
           else
              ! If the volume is 0, then p_gas is not used, but we
              ! don't want to trigger the line bellow (because of
              ! floating point exception dark magic)
              p_gas(iAGN) = 0.0d0
           end if
           ! Note that p_gas is \propto v**2, since vol_gas is a mass
        endif
     endif
     Efeed_new(iAGN)=Efeed_new(iAGN)+EAGN(iAGN)
  end do
  EsaveAGN=0d0

  ! Loop over levels
  do ilevel=levelmin,nlevelmax
     ! Computing local volume (important for averaging hydro quantities)
     dx=0.5D0**ilevel
     dx_loc=dx*scale
     vol_loc=dx_loc**ndim
     ! Cells center position relative to grid center position
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        xc(ind,1)=(dble(ix)-0.5D0)*dx
        xc(ind,2)=(dble(iy)-0.5D0)*dx
        xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Loop over grids
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do

        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Flag leaf cells
           do i=1,ngrid
              ok(i)=son(ind_cell(i))==0
           end do

           do i=1,ngrid
              if(ok(i))then
                 ! Get gas cell position
                 x=(xg(ind_grid(i),1)+xc(ind,1)-skip_loc(1))*scale
                 y=(xg(ind_grid(i),2)+xc(ind,2)-skip_loc(2))*scale
                 z=(xg(ind_grid(i),3)+xc(ind,3)-skip_loc(3))*scale

                 box%origin = [x, y, z]
                 box%extents = [dx_loc, dx_loc, dx_loc] / 2._dp ! This is half the width of the box

                 do iAGN=1,nAGN

                    ! ------------------------------------------
                    ! case 0: Some energy has not been released
                    ! in previous time step -> do thermal input
                    ! ------------------------------------------
                    if(ok_save(iAGN))then
                       dxx=x-xAGN(iAGN,1)
                       dyy=y-xAGN(iAGN,2)
                       dzz=z-xAGN(iAGN,3)
                       dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz
                       if(dr_AGN.le.rmax2)then
                          ekk=0d0
                          d=max(uold(ind_cell(i),1), smallr)
                          do idim=1,ndim
                             ekk=ekk+0.5d0*uold(ind_cell(i),idim+1)**2/d
                          end do
                          etot=uold(ind_cell(i),5)
                          eint=etot - ekk
                          T2_1=(gamma-1d0)*eint/d*scale_T2
                          if(T2_1 .lt. T2maxAGNz)then
                             eint=eint + p_gas(iAGN)*d
                             T2_2=(gamma-1d0)*eint/d*scale_T2
                             if(T2_2 .le. T2maxAGNz)then
                                uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_gas(iAGN)*d
                             else
                                uold(ind_cell(i),5)=uold(ind_cell(i),5)+T2maxAGNz/scale_T2/(gamma-1d0)*d
                                EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGNz)/scale_T2/(gamma-1d0)*d*vol_loc
                             endif
                          else
                             EsaveAGN(iAGN)=EsaveAGN(iAGN)+p_gas(iAGN)*d*vol_loc
                          endif
                       endif

                    ! ------------------------------------------
                    ! case 1: All energy has been released in
                    ! previous time step -> choose between kinetic
                    ! or thermal feedback
                    ! ------------------------------------------
                    else if(ok_blast_agn(iAGN))then
                       ! ------------------------------------------
                       ! Jet case
                       ! ------------------------------------------
                       if(X_radio(iAGN).lt.X_floor)then

                          dxx=x-xAGN(iAGN,1)
                          dyy=y-xAGN(iAGN,2)
                          dzz=z-xAGN(iAGN,3)
                          dr_AGN=norm2([dxx, dyy, dzz])
                          jtot=norm2(jAGN(iAGN, :))
                          if(jtot > 0._dp)then
                             j_x=jAGN(iAGN,1)/jtot
                             j_y=jAGN(iAGN,2)/jtot
                             j_z=jAGN(iAGN,3)/jtot

                             dzjet= dxx*j_x + dyy*j_y + dzz*j_z

                             capsule%segment%start = xAGN(iAGN, :) - rmax * jAGN(iAGN, :) / jtot
                             capsule%segment%end = xAGN(iAGN, :) + rmax * jAGN(iAGN, :) / jtot
                             capsule%r = rmax
                             call LineSegmentBoxDistanceSquared(box, capsule%segment, d2)
                             ! call CapsuleBoxDistanceSquared(box, capsule, d2)

                             ! Check if the cell lies within the AGN jet cylindre
                             if (d2 <= rmax2) then
                                if (AGN_integration == AGN_VOLUME_INTEGRATION_MC) then
                                   call CapsuleBoxIntegrateMC(psy_function_loc, &
                                        box, capsule, tmp_volume, psy)
                                   if (tmp_volume > 0) then
                                      psy = psy / tmp_volume
                                   else
                                      psy = 1d-10
                                   end if
                                   vol_gas(iAGN) = vol_gas(iAGN) + tmp_volume
                                else if (AGN_integration == AGN_VOLUME_INTEGRATION_bisect) then
                                   call CapsuleBoxIntegrate(psy_function_loc, &
                                        box, capsule, tmp_volume, psy)
                                   if (tmp_volume > 0) then
                                      psy = psy / tmp_volume
                                   else
                                      psy = 1d-10
                                   end if
                                   vol_gas(iAGN) = vol_gas(iAGN) + tmp_volume
                                else
                                   vol_gas(iAGN) = vol_gas(iAGN) + vol_loc
                                   drjet=sqrt(max(dr_AGN-dzjet*dzjet, 0._dp))
                                   psy = exp(-drjet**2/2d0/rmax2)
                                end if

                                if (dzjet < 0._dp) then
                                   u=-j_x*uBlast(iAGN)
                                   v=-j_y*uBlast(iAGN)
                                   w=-j_z*uBlast(iAGN)
                                else if (dzjet > 0._dp) then
                                   u= j_x*uBlast(iAGN)
                                   v= j_y*uBlast(iAGN)
                                   w= j_z*uBlast(iAGN)
                                endif

                                ekk=0d0
                                ! Compute kinetic energy before velocity input
                                do idim=1,ndim
                                   ekk=ekk+0.5d0*uold(ind_cell(i),idim+1)**2/max(uold(ind_cell(i),1), smallr)
                                end do
                                ekkold=ekk
                                ! Compute total energy before energy input
                                etot=uold(ind_cell(i),5)
                                ! Compute internal energy before energy input
                                eint=etot - ekk
                                ! Compute temperature T2=T/mu in Kelvin before energy input
                                T2_1=(gamma-1d0)*eint/max(uold(ind_cell(i),1), smallr)*scale_T2

                                d_gas=mAGN(iAGN)/psy_norm(iAGN)
                                ! Compute the density and the metal density of the cell
                                uold(ind_cell(i),1)=uold(ind_cell(i),1) + d_gas * psy
                                d=max(uold(ind_cell(i),1), smallr)
                                do ivar= imetal, nvar
                                   uold(ind_cell(i), ivar) = uold(ind_cell(i),ivar) + passiveAGN(iAGN, ivar)*d_gas*psy
                                end do
                                ! Velocity at a given dr_AGN linearly interpolated between zero and uBlast
                                u= u + vAGN(iAGN,1)
                                v= v + vAGN(iAGN,2)
                                w= w + vAGN(iAGN,3)

                                ! Add each momentum component of the jet to the gas
                                uold(ind_cell(i),2)=uold(ind_cell(i),2)+d_gas*u * psy
                                uold(ind_cell(i),3)=uold(ind_cell(i),3)+d_gas*v * psy
                                uold(ind_cell(i),4)=uold(ind_cell(i),4)+d_gas*w * psy

                                ekk=0d0
                                ! Compute kinetic energy after velocity input
                                do idim=1,ndim
                                   ekk=ekk+0.5*uold(ind_cell(i),idim+1)**2/d
                                end do
                                if(T2_1 .ge. T2maxAGNz)then
                                   etot=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)*psy &
                                        & + p_gas(iAGN)
                                   ! Update total energy with new kinetic energy and
                                   ! old internal energy (Temperature does not increase!)
                                   uold(ind_cell(i),5)=ekk+eint
                                   T2_2=(gamma-1d0)*(etot-uold(ind_cell(i),5))/d*scale_T2
                                   EsaveAGN(iAGN)=EsaveAGN(iAGN)+T2_2/scale_T2/(gamma-1d0)*d*vol_loc
                                else
                                   ! Compute new total energy
                                   etot=uold(ind_cell(i),5)+0.5*d_gas*(u*u+v*v+w*w)*psy &
                                        & + p_gas(iAGN)
                                   ! Compute new internal energy
                                   eint=etot - ekk
                                   ! Compute T2=T/mu in Kelvin
                                   T2_2=(gamma-1d0)*eint/d*scale_T2
                                   if(T2_2 .le. T2maxAGNz)then
                                      uold(ind_cell(i),5)=ekk+T2_2/scale_T2/(gamma-1d0)*d
                                   else
                                      uold(ind_cell(i),5)=ekk+T2maxAGNz/scale_T2/(gamma-1d0)*d
                                      EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGNz)/scale_T2/(gamma-1d0)*d*vol_loc
                                   endif
                                endif

                             endif
                             ! Jet case with jsink=0
                          else
                             if(dr_AGN.le.rmax2)then
                                uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_gas(iAGN)
                             endif
                          endif

                          ! ------------------------------------------
                          ! Thermal case
                          ! ------------------------------------------
                       else
                          dxx=x-xAGN(iAGN,1)
                          dyy=y-xAGN(iAGN,2)
                          dzz=z-xAGN(iAGN,3)
                          dr_AGN=dxx*dxx+dyy*dyy+dzz*dzz
                          if(dr_AGN.le.rmax2)then
                             ekk=0d0
                             d=max(uold(ind_cell(i),1), smallr)
                             do idim=1,ndim
                                ekk=ekk+0.5d0*uold(ind_cell(i),idim+1)**2/d
                             end do
                             etot=uold(ind_cell(i),5)
                             eint=etot - ekk
                             T2_1=(gamma-1d0)*eint/d*scale_T2
                             if(T2_1 .lt. T2maxAGNz)then
                                eint=eint + p_gas(iAGN)*d
                                T2_2=(gamma-1d0)*eint/d*scale_T2
                                if(T2_2 .le. T2maxAGNz)then
                                   uold(ind_cell(i),5)=uold(ind_cell(i),5)+p_gas(iAGN)*d
                                else
                                   uold(ind_cell(i),5)=uold(ind_cell(i),5)+T2maxAGNz/scale_T2/(gamma-1d0)*d
                                   EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGNz)/scale_T2/(gamma-1d0)*d*vol_loc
                                endif
                             else
                                EsaveAGN(iAGN)=EsaveAGN(iAGN)+p_gas(iAGN)*d*vol_loc
                             endif

                          endif
                       endif

                    endif
                    !End of ok_blast_agn

                 end do
              endif
           end do

        end do
        ! End loop over cells
     end do
     ! End loop over grids
  end do
  ! End loop over levels

  do iAGN=1,nAGN

     if(ind_blast(iAGN)>0)then

     ! ------------------------------------------
     ! case 0: Some energy has not been released
     ! in previous time step -> do thermal input
     ! ------------------------------------------
     if(ok_save(iAGN))then
        if(vol_gas(iAGN)==0d0)then
           ekk=0d0
           d=max(uold(ind_blast(iAGN),1), smallr)
           do idim=1,ndim
              ekk=ekk+0.5d0*uold(ind_blast(iAGN),idim+1)**2/d
           end do
           etot=uold(ind_blast(iAGN),5)
           eint=etot - ekk
           T2_1=(gamma-1d0)*eint/d*scale_T2
           if(T2_1 .lt. T2maxAGNz)then
              eint=eint + EAGN(iAGN)/vol_blast(iAGN)
              T2_2=(gamma-1d0)*eint/d*scale_T2
              if(T2_2 .le. T2maxAGNz)then
                 uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+EAGN(iAGN)/vol_blast(iAGN)
              else
                 uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+T2maxAGNz/scale_T2/(gamma-1d0)*d
                 EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGNz)/scale_T2/(gamma-1d0)*d*vol_loc
              endif
           else
              EsaveAGN(iAGN)=EsaveAGN(iAGN)+EAGN(iAGN)
           endif

        endif
     ! ------------------------------------------
     ! case 1: All energy has been released in
     ! previous time step -> choose between kinetic
     ! or thermal feedback
     ! ------------------------------------------
     else
        if(vol_gas(iAGN)==0d0.and.ok_blast_agn(iAGN))then
           ! ------------------------------------------
           ! Jet case
           ! ------------------------------------------
           if(X_radio(iAGN).lt.X_floor)then
              ! Here vol_blast lies for the cell volume where the AGN sits in
              d_gas=mAGN(iAGN)/vol_blast(iAGN)
              u=vAGN(iAGN,1)
              v=vAGN(iAGN,2)
              w=vAGN(iAGN,3)

              uold(ind_blast(iAGN),1)=uold(ind_blast(iAGN),1)+d_gas
              do ivar= imetal, nvar
                 uold(ind_blast(iAGN), ivar) = uold(ind_blast(iAGN),ivar) + passiveAGN(iAGN, ivar)*d_gas
              end do
              ekk=0d0
              d=max(uold(ind_blast(iAGN),1), smallr)
              do idim=1,ndim
                 ekk=ekk+0.5d0*uold(ind_blast(iAGN),idim+1)**2/d
              end do
              etot=uold(ind_blast(iAGN),5)
              eint=etot - ekk
              T2_1=(gamma-1d0)*eint/d*scale_T2
              if(T2_1 .lt. T2maxAGNz)then
                 eint=eint + EAGN(iAGN)/vol_blast(iAGN)
                 T2_2=(gamma-1d0)*eint/d*scale_T2
                 if(T2_2 .le. T2maxAGNz)then
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+EAGN(iAGN)/vol_blast(iAGN)
                 else
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+T2maxAGNz/scale_T2/(gamma-1d0)*d
                    EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGNz)/scale_T2/(gamma-1d0)*d*vol_loc
                 endif
              else
                 EsaveAGN(iAGN)=EsaveAGN(iAGN)+EAGN(iAGN)
              endif

           ! ------------------------------------------
           ! Thermal case
           ! ------------------------------------------
           else
              ekk=0d0
              d=max(uold(ind_blast(iAGN),1), smallr)
              do idim=1,ndim
                 ekk=ekk+0.5d0*uold(ind_blast(iAGN),idim+1)**2/d
              end do
              etot=uold(ind_blast(iAGN),5)
              eint=etot - ekk
              T2_1=(gamma-1d0)*eint/d*scale_T2
              if(T2_1 .lt. T2maxAGNz)then
                 eint=eint + EAGN(iAGN)/vol_blast(iAGN)
                 T2_2=(gamma-1d0)*eint/d*scale_T2
                 if(T2_2 .le. T2maxAGNz)then
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+EAGN(iAGN)/vol_blast(iAGN)
                 else
                    uold(ind_blast(iAGN),5)=uold(ind_blast(iAGN),5)+T2maxAGNz/scale_T2/(gamma-1d0)*d
                    EsaveAGN(iAGN)=EsaveAGN(iAGN)+(T2_2-T2maxAGNz)/scale_T2/(gamma-1d0)*d*vol_loc
                 endif
              else
                 EsaveAGN(iAGN)=EsaveAGN(iAGN)+EAGN(iAGN)
              endif

           endif
        endif
     endif

     endif
  end do

  if(verbose)write(*,*)'Exiting AGN_blast'
contains
  real(dp) function psy_function_loc(X) result(psy)
    ! A note on this function: it is called in the frame of the
    ! capsule, in which the last dimension is // to
    real(dp), intent(in) :: X(3)

    call psy_function(X, rmax2, psy)

  end function psy_function_loc

end subroutine AGN_blast
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getAGNonmyid(isink_myid,nsink_myid)
  use amr_commons
  use pm_commons
  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine check which BHs stand in the cpu myid
  ! Yohan Dubois, December 15th, 2010
  !------------------------------------------------------------------------
  integer,dimension(1:nsink), intent(out)::isink_myid
  integer, intent(out)::nsink_myid

  integer ::ii
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  integer::lmin,isink,nx_loc,ilevel,bit_length,maxdom,icpu
  integer::imin,jmin,kmin,imax,jmax,kmax,ndom,impi,i,j,k,ncpu_read
  integer,dimension(1:ncpu)::cpu_list
  logical,dimension(1:ncpu)::cpu_read
  real(dp)::scale,dx,dx_min,drsink
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp),dimension(1:3)::skip_loc
  real(dp)::xxmin,yymin,zzmin,xxmax,yymax,zzmax,dmax
  real(qdp),dimension(1:8)::bounding_min,bounding_max
  real(qdp)::dkey,order_min(1),oneqdp=1.0

  ! Mesh spacing in that level
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  ! Maximum radius of the ejecta
  drsink=2d0*MAX(1d0*dx_min*scale_l/aexp,rAGN*3.08d21)
  drsink=drsink/scale_l
  !-----------------------
  ! Map parameters
  !-----------------------
  isink_myid=0
  ii=0
  do isink=1,nsink

     cpu_read=.false.
     ! Compute boundaries for the sink cube of influence
     xxmin=(xsink(isink,1)-drsink)/scale ; xxmax=(xsink(isink,1)+drsink)/scale
     yymin=(xsink(isink,2)-drsink)/scale ; yymax=(xsink(isink,2)+drsink)/scale
     zzmin=(xsink(isink,3)-drsink)/scale ; zzmax=(xsink(isink,3)+drsink)/scale

     if(TRIM(ordering).eq.'hilbert')then

        dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
        do ilevel=1,nlevelmax
           dx=0.5d0**ilevel
           if(dx.lt.dmax)exit
        end do
        lmin=ilevel
        bit_length=lmin-1
        maxdom=2**bit_length
        imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
        if(bit_length>0)then
           imin=int(xxmin*dble(maxdom))
           imax=imin+1
           jmin=int(yymin*dble(maxdom))
           jmax=jmin+1
           kmin=int(zzmin*dble(maxdom))
           kmax=kmin+1
        endif

        !dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
        dkey=(real(2**(nlevelmax+1),kind=qdp)/real(maxdom,kind=qdp))**ndim
        ndom=1
        if(bit_length>0)ndom=8
        idom(1)=imin; idom(2)=imax
        idom(3)=imin; idom(4)=imax
        idom(5)=imin; idom(6)=imax
        idom(7)=imin; idom(8)=imax
        jdom(1)=jmin; jdom(2)=jmin
        jdom(3)=jmax; jdom(4)=jmax
        jdom(5)=jmin; jdom(6)=jmin
        jdom(7)=jmax; jdom(8)=jmax
        kdom(1)=kmin; kdom(2)=kmin
        kdom(3)=kmin; kdom(4)=kmin
        kdom(5)=kmax; kdom(6)=kmax
        kdom(7)=kmax; kdom(8)=kmax

        do i=1,ndom
           if(bit_length>0)then
              call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
           else
              order_min(1)=0
           endif
           bounding_min(i)=(order_min(1))*dkey
           bounding_max(i)=(order_min(1)+oneqdp)*dkey
        end do

        cpu_min=0; cpu_max=0
        do impi=1,ncpu
           do i=1,ndom
              if (   bound_key(impi-1).le.bounding_min(i).and.&
                   & bound_key(impi  ).gt.bounding_min(i))then
                 cpu_min(i)=impi
              endif
              if (   bound_key(impi-1).lt.bounding_max(i).and.&
                   & bound_key(impi  ).ge.bounding_max(i))then
                 cpu_max(i)=impi
              endif
           end do
        end do

        ncpu_read=0
        do i=1,ndom
           do j=cpu_min(i),cpu_max(i)
              if(.not. cpu_read(j))then
                 ncpu_read=ncpu_read+1
                 cpu_list(ncpu_read)=j
                 cpu_read(j)=.true.
              endif
           enddo
        enddo
     else
        ncpu_read=ncpu
        do j=1,ncpu
           cpu_list(j)=j
        end do
     end  if

     ! Create the index array for sinks in processor myid
     do k=1,ncpu_read
        icpu=cpu_list(k)
        if(icpu==myid)then
           ii=ii+1
           isink_myid(ii)=isink
        endif
     enddo

  enddo

  ! Number of sinks in processor myid
  nsink_myid=ii

end subroutine getAGNonmyid
!################################################################
!################################################################
!################################################################
!################################################################
subroutine update_sink_position_velocity
  use pm_commons
  use amr_commons
  use mpi_mod
  implicit none
  !------------------------------------------------------------------------
  ! This routine updates position and velocity of sink particles.
  !------------------------------------------------------------------------
  integer::isink,idim,size_mpi,info,nx_loc
  real(dp)::vdum,ncloud,scale
  real(dp),dimension(1:3)::xbound

  ! Mesh spacing in that level
  xbound(1:3)=(/dble(nx),dble(ny),dble(nz)/)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)

  if(fix_smbh_position)return

  ! update the sink particle position based on total sink particles
#ifndef WITHOUTMPI
  size_mpi=nsinkmax*(nlevelmax-levelmin+1)*(ndim*2+1)
  call MPI_ALLREDUCE(sink_stat,sink_stat_all,size_mpi,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
  sink_stat_all=sink_stat
#endif

  if(myid==1.and.debug)write(*,*)"NCLOUD first sink in update_sink_position_velocity:",sum(sink_stat_all(1,levelmin:nlevelmax,7))
  do isink=1,nsink
     ncloud = sum(sink_stat_all(isink,levelmin:nlevelmax,ndim*2+1))
     if(ncloud>1)then
        do idim=1,ndim
           vdum = sum(sink_stat_all(isink,levelmin:nlevelmax,idim))
           xsink(isink,idim)=vdum/ncloud
           if(xsink(isink,idim)>scale*xbound(idim))then
              xsink(isink,idim)=xsink(isink,idim)-scale*xbound(idim)
           endif
           if(xsink(isink,idim)<0.0)then
              xsink(isink,idim)=xsink(isink,idim)+scale*xbound(idim)
           endif
        enddo
        do idim=1,ndim
           vdum = sum(sink_stat_all(isink,levelmin:nlevelmax,idim+ndim))
           vsink(isink,idim)=vdum/ncloud
        enddo
     endif
  enddo

end subroutine update_sink_position_velocity

! This routine is necessary for clump_merger
!################################################################
!################################################################
!################################################################
!################################################################
subroutine true_max(x,y,z,ilevel)
  use amr_commons
  use pm_commons
  use hydro_commons
  use clfind_commons, only:ivar_clump
  use poisson_commons, only:rho
  implicit none
  real(dp)::x,y,z
  integer::ilevel

  !----------------------------------------------------------------------------
  ! Description: This subroutine takes the cell of maximum density and computes
  ! the true maximum by expanding the density around the cell center to second order.
  !----------------------------------------------------------------------------

  integer::k,j,i,nx_loc,counter
  integer,dimension(1:nvector)::cell_index,cell_lev
  real(dp)::det,dx,dx_loc,scale,disp_max
  real(dp),dimension(-1:1,-1:1,-1:1)::cube3
  real(dp),dimension(1:nvector,1:ndim)::xtest
  real(dp),dimension(1:ndim)::gradient,displacement
  real(dp),dimension(1:ndim,1:ndim)::hess,minor

#if NDIM==3

  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale

  counter=0
  do i=-1,1
     do j=-1,1
        do k=-1,1
           counter=counter+1
           xtest(counter,1)=x+i*dx_loc
           xtest(counter,2)=y+j*dx_loc
           xtest(counter,3)=z+k*dx_loc
        end do
     end do
  end do

  call get_cell_index(cell_index,cell_lev,xtest,ilevel,counter)

  counter=0
  if(ivar_clump==0)then
     do i=-1,1
        do j=-1,1
           do k=-1,1
              counter=counter+1
              cube3(i,j,k)=rho(cell_index(counter))
           end do
        end do
     end do
  else if(hydro)then
     do i=-1,1
        do j=-1,1
           do k=-1,1
              counter=counter+1
              cube3(i,j,k)=uold(cell_index(counter),1)
           end do
        end do
     end do
  else
     return
  end if

! compute gradient
  gradient(1)=0.5*(cube3(1,0,0)-cube3(-1,0,0))/dx_loc
  gradient(2)=0.5*(cube3(0,1,0)-cube3(0,-1,0))/dx_loc
  gradient(3)=0.5*(cube3(0,0,1)-cube3(0,0,-1))/dx_loc

  if (maxval(abs(gradient(1:3)))==0.)return

  ! compute hessian
  hess(1,1)=(cube3(1,0,0)+cube3(-1,0,0)-2*cube3(0,0,0))/dx_loc**2.
  hess(2,2)=(cube3(0,1,0)+cube3(0,-1,0)-2*cube3(0,0,0))/dx_loc**2.
  hess(3,3)=(cube3(0,0,1)+cube3(0,0,-1)-2*cube3(0,0,0))/dx_loc**2.

  hess(1,2)=0.25*(cube3(1,1,0)+cube3(-1,-1,0)-cube3(1,-1,0)-cube3(-1,1,0))/dx_loc**2.
  hess(2,1)=hess(1,2)
  hess(1,3)=0.25*(cube3(1,0,1)+cube3(-1,0,-1)-cube3(1,0,-1)-cube3(-1,0,1))/dx_loc**2.
  hess(3,1)=hess(1,3)
  hess(2,3)=0.25*(cube3(0,1,1)+cube3(0,-1,-1)-cube3(0,1,-1)-cube3(0,-1,1))/dx_loc**2.
  hess(3,2)=hess(2,3)

  !determinant
  det=hess(1,1)*hess(2,2)*hess(3,3)+hess(1,2)*hess(2,3)*hess(3,1)+hess(1,3)*hess(2,1)*hess(3,2) &
       -hess(1,1)*hess(2,3)*hess(3,2)-hess(1,2)*hess(2,1)*hess(3,3)-hess(1,3)*hess(2,2)*hess(3,1)

  !matrix of minors
  minor(1,1)=hess(2,2)*hess(3,3)-hess(2,3)*hess(3,2)
  minor(2,2)=hess(1,1)*hess(3,3)-hess(1,3)*hess(3,1)
  minor(3,3)=hess(1,1)*hess(2,2)-hess(1,2)*hess(2,1)

  minor(1,2)=-1.*(hess(2,1)*hess(3,3)-hess(2,3)*hess(3,1))
  minor(2,1)=minor(1,2)
  minor(1,3)=hess(2,1)*hess(3,2)-hess(2,2)*hess(3,1)
  minor(3,1)=minor(1,3)
  minor(2,3)=-1.*(hess(1,1)*hess(3,2)-hess(1,2)*hess(3,1))
  minor(3,2)=minor(2,3)


  !displacement of the true max from the cell center
  displacement=0.
  do i=1,3
     do j=1,3
        displacement(i)=displacement(i)-minor(i,j)/(det+1.d10*tiny(0.d0))*gradient(j)
     end do
  end do

  !clipping the displacement in order to keep max in the cell
  disp_max=maxval(abs(displacement(1:3)))
  if (disp_max > dx_loc*0.499999)then
     displacement(1)=displacement(1)/disp_max*dx_loc*0.499999
     displacement(2)=displacement(2)/disp_max*dx_loc*0.499999
     displacement(3)=displacement(3)/disp_max*dx_loc*0.499999
  end if

  x=x+displacement(1)
  y=y+displacement(2)
  z=z+displacement(3)

#endif
end subroutine true_max
!################################################################
!################################################################
!################################################################
subroutine get_drag_part(ilevel)
  !------------------------------------------------------------------------
  ! Hugo Pfister
  ! This routine computes the drag from collisionless particles
  ! onto sink particles.
  ! On exit velocity of clouds is modified
  ! We distinguish stars and DM: index 2 is DM - 1 is stars
  !------------------------------------------------------------------------
  use pm_commons
  use amr_commons
  use hydro_commons
  use cooling_module, ONLY: rhoc, mH, twopi
  use poisson_commons, only: cic_levelmax
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif

  integer::ilevel

  integer :: icpu, igrid, next_part, nx_loc, jgrid
  integer :: npart1, isink, ipart, jpart, idim, itype
  integer :: ii,nsink_cell

  real(dp) :: pi, factG, c, dx_min, scale, r2_cloud, vol_cloud
  real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp) :: b90, Rsh, bmin, CoulombLog, factor
  real(dp) :: dx_dm, r2_dm, vol_dm, rmax_dm, d2_sink_part
  real(dp) :: vrel_part_norm, volume, xc,yc,zc
  real(dp) :: rmax, dx, dx2cell_sink
  real(dp),dimension(1:3)::skip_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,*) 'Entering get_drag_part for level', ilevel

  pi=twopi/2d0
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx = scale*0.5d0**ilevel/aexp
  dx_min=scale*0.5d0**(nlevelmax-nlevelsheld)/aexp
  rmax = DF_ncells*dx_min
  r2_cloud = rmax**2
  vol_cloud = 4d0/3d0*pi*(rmax)**3

  if (cic_levelmax.gt.0) then
    dx_dm = scale*0.5d0**cic_levelmax/aexp
  else
    dx_dm = dx_min
  end if
  rmax_dm = DF_ncells*dx_dm
  r2_dm = rmax_dm**2
  vol_dm = 4d0/3d0*pi*rmax_dm**3
  dx2cell_sink = (2*max(dx, rmax_dm)+rmax_dm)**2

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  factG=1
  if(cosmo)factG=3d0/8d0/pi*omega_m*aexp
  c=3d10/scale_v

  !We compute the relative velocity of sinks wrt the background
  !using the value from previous timestep to avoid doubling
  !loops on cells.
  v_background(1:nsink, 1:ndim, 1:2) = 0d0
  m_background(1:nsink, 1:2) = tiny(0d0)
  do ii = levelmin, nlevelmax
    do isink = 1,nsink
        do itype = 1, 2
            do idim = 1, ndim
                v_background(isink, idim, itype) =&
                 v_background(isink, idim, itype) * m_background(isink, itype)+&
                 v_DF(isink, ii, idim, itype) * mass_DF(isink, ii, itype)
            end do
            m_background(isink, itype) = m_background(isink, itype) + mass_DF(isink, ii, itype)
            do idim = 1, ndim
                v_background(isink, idim, itype) =&
                 v_background(isink, idim, itype) / m_background(isink, itype)
            end do
        end do
    end do
  end do

  vrel_sink (1:nsink, 1:ndim, 1:2) = 0d0
  vrel_sink_norm(1:nsink, 1:2) = tiny(0d0)
  do isink = 1, nsink
    do itype = 1, 2
        do idim = 1, ndim
            vrel_sink(isink, idim, itype) = vsink(isink, idim) - v_background(isink, idim, itype)
            vrel_sink_norm(isink, itype) = vrel_sink_norm(isink, itype) + vrel_sink(isink, idim, itype)**2
        end do
        vrel_sink_norm(isink, itype) = max(sqrt(vrel_sink_norm(isink, itype)), tiny(0d0))
    end do
  end do

  !We compute the different quantities needed or DF: fraction of stars moving
  !slower, factor for fast moving stars. We also update the value of v_DF and
  !mass_DF which will be used to compute DF at next timestep and higher levels
  !CARE: Here what we call v_DF is in fact the momentum to correctly
  !mass-weight between CPUS
  v_DFnew(1:nsink, 1:ndim, 1:2) = 0d0
  mass_DFnew(1:nsink, 1:2) = tiny(0d0)
  mass_lowspeednew(1:nsink, 1:2) = 0d0
  fact_fastnew(1:nsink, 1:2) = 0d0
  n_partnew(1:nsink, 1:2) = 0
  ! Loop over cpus
  do icpu=1, ncpu
     igrid=headl(icpu,ilevel)

     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        nsink_cell=0
      ! Loop over particles
      if(npart1>0)then
         !We check which sink might be close (speeds up calculation)
         xc = (xg(igrid, 1) - skip_loc(1))*scale
         yc = (xg(igrid, 2) - skip_loc(2))*scale
         zc = (xg(igrid, 3) - skip_loc(3))*scale
         do isink = 1, nsink
             dx = abs(xc-xsink(isink, 1))
             if (dx.gt.0.5*scale) dx = scale - dx
             d2_sink_part = dx**2
             if (d2_sink_part.lt.dx2cell_sink) then
                 dx = abs(yc-xsink(isink, 2))
                 if (dx.gt.0.5*scale) dx = scale - dx
                 d2_sink_part = d2_sink_part+ dx**2
                 if (d2_sink_part.lt.dx2cell_sink) then
                     dx = abs(zc-xsink(isink, 3))
                     if (dx.gt.0.5*scale) dx = scale - dx
                     d2_sink_part = d2_sink_part+ dx**2
                     if (d2_sink_part.lt.dx2cell_sink) then
                         nsink_cell = nsink_cell+1
                         sink_cell(nsink_cell) = isink
                     end if
                 end if
             end if
         end do

         ipart=headp(igrid)
         ! Loop over particles
         do jpart=1, npart1
             ! Save next particle   <--- Very important !!!
             next_part=nextp(ipart)
             !Stellar case (1)
             !if (tp(ipart).ne.0.d0) then
             if (is_star(typep(ipart))) then
                do ii = 1, nsink_cell
                   isink = sink_cell(ii)
                   dx = abs(xp(ipart, 1)-xsink(isink, 1))
                   if (dx.gt.0.5*scale) dx = scale - dx
                   d2_sink_part = dx**2

                   if (d2_sink_part.lt.r2_cloud) then
                      dx = abs(xp(ipart, 2)-xsink(isink, 2))
                      if (dx.gt.0.5*scale) dx = scale - dx
                      d2_sink_part = d2_sink_part+ dx**2

                      if (d2_sink_part.lt.r2_cloud) then
                         dx = abs(xp(ipart, 3)-xsink(isink, 3))
                         if (dx.gt.0.5*scale) dx = scale - dx
                         d2_sink_part = d2_sink_part+ dx**2

                         if (d2_sink_part.lt.r2_cloud) then
                            v_DFnew(isink, 1:ndim, 1) = v_DFnew(isink, 1:ndim, 1) + mp(ipart) * vp(ipart, 1:ndim)
                            mass_DFnew(isink, 1)      = mass_DFnew(isink, 1)      + mp(ipart)

                            vrel_part_norm = 0d0
                            do idim = 1, ndim
                               vrel_part_norm = vrel_part_norm + (vp(ipart, idim) - v_background(isink, idim, 1))**2
                            end do
                            vrel_part_norm = max(tiny(0d0), sqrt(vrel_part_norm))

                            n_partnew(isink, 1) = n_partnew(isink, 1) +1
                            if (vrel_part_norm.le.vrel_sink_norm(isink, 1)) then
                               mass_lowspeednew(isink, 1) = mass_lowspeednew(isink, 1) + mp(ipart)
                            else
                               fact_fastnew(isink, 1) = fact_fastnew(isink, 1) + mp(ipart) *&
                                    (log((vrel_part_norm + vrel_sink_norm(isink, 1))/(vrel_part_norm-vrel_sink_norm(isink,1))) - 2*vrel_sink_norm(isink, 1) / vrel_part_norm)
                            end if
                         end if
                      end if
                   end if
                end do

             !DM case (2)
             !else if (idp(ipart).gt.0 .and. tp(ipart).eq.0.d0) then
             else if (is_dm(typep(ipart))) then
                do ii = 1, nsink_cell
                   isink = sink_cell(ii)
                   dx = abs(xp(ipart, 1)-xsink(isink, 1))
                   if (dx.gt.0.5*scale) dx = scale - dx
                   d2_sink_part = dx**2

                   if (d2_sink_part.lt.r2_dm) then
                      dx = abs(xp(ipart, 2)-xsink(isink, 2))
                      if (dx.gt.0.5*scale) dx = scale - dx
                      d2_sink_part = d2_sink_part+ dx**2

                      if (d2_sink_part.lt.r2_dm) then
                         dx = abs(xp(ipart, 3)-xsink(isink, 3))
                         if (dx.gt.0.5*scale) dx = scale - dx
                         d2_sink_part = d2_sink_part+ dx**2

                         if (d2_sink_part.lt.r2_dm) then
                            v_DFnew(isink, 1:ndim, 2) = v_DFnew(isink, 1:ndim, 2) + mp(ipart) * vp(ipart, 1:ndim)
                            mass_DFnew(isink, 2)      = mass_DFnew(isink, 2)      + mp(ipart)

                            vrel_part_norm = 0d0
                            do idim = 1, ndim
                               vrel_part_norm = vrel_part_norm + (vp(ipart, idim) - v_background(isink, idim, 2))**2
                            end do
                            vrel_part_norm = max(tiny(0d0), sqrt(vrel_part_norm))

                            n_partnew(isink, 2) = n_partnew(isink, 2) +1
                            if (vrel_part_norm.le.vrel_sink_norm(isink, 2)) then
                               mass_lowspeednew(isink, 2) = mass_lowspeednew(isink, 2) + mp(ipart)
                            else
                               fact_fastnew(isink, 2) = fact_fastnew(isink, 2) + mp(ipart) *&
                                    (log((vrel_part_norm + vrel_sink_norm(isink, 2))/(vrel_part_norm-vrel_sink_norm(isink,2))) - 2*vrel_sink_norm(isink, 2) / vrel_part_norm)
                            end if
                         end if
                      end if
                   end if
                end do
             end if
             ipart=next_part  ! Go to next particle
          end do
       endif

       igrid=next(igrid)   ! Go to next grid
    end do
    ! End loop over grids

  end do
  ! End loop over cpus


  !We gather all the quantities: mass, mass moving slow etc... from all CPUs
  if(nsink>0)then
#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(v_DFnew,v_DFall,nsinkmax*ndim*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mass_DFnew, mass_DFall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(n_partnew, n_partall,nsinkmax*2,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(mass_lowspeednew, mass_lowspeedall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     call MPI_ALLREDUCE(fact_fastnew, fact_fastall,nsinkmax*2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
#else
     v_DFall(1:nsink, 1:ndim, 1:2)  = v_DFnew(1:nsink, 1:ndim, 1:2)
     mass_DFall(1:nsink, 1:2)       = mass_DFnew(1:nsink, 1:2)
     n_partall(1:nsink, 1:2)        = n_partnew(1:nsink, 1:2)
     mass_lowspeedall(1:nsink, 1:2) = mass_lowspeednew(1:nsink, 1:2)
     fact_fastall(1:nsink, 1:2)     = fact_fastnew(1:nsink, 1:2)
#endif
  endif

  !Three things are done in this loop:
  !   1) we fill the table that contains all the quantities of all levels with the
  !   ones we just computed
  !   2) we switch v_DF which was momentum to realy velocity)
  !   3) we remove the momentum due to matter which is at THIS level
  !   using eq. (6) and (7) from Antonini+11
  do isink = 1 ,nsink
     do itype = 1,2
        mass_DF(isink, ilevel, itype)       = mass_DFall(isink, itype)
        n_part(isink, ilevel, itype)        = n_partall(isink, itype)
        mass_lowspeed(isink, ilevel, itype) = mass_lowspeedall(isink, itype)
        fact_fast(isink, ilevel, itype)     = fact_fastall(isink, itype)

        b90 = factG * msink(isink) / vrel_sink_norm(isink, itype)**2
        Rsh = factG * msink(isink) / c**2
        bmin = max(b90, Rsh)
        if (itype.eq.1) then
           CoulombLog = DF_ncells*dx_min / bmin
           volume = vol_cloud
        else
           CoulombLog = DF_ncells*dx_dm / bmin
           volume = vol_dm
        end if
        if (CoulombLog .gt. 1) then
           factor = 4*pi*factG**2*msink(isink)/vrel_sink_norm(isink, itype)**3 * &
                & (mass_lowspeedall(isink, itype)*log(CoulombLog)+&
                &  fact_fastall(isink, itype)) / volume
        else
           factor = 0
        end if
        ! At the first timestep, the dt hasn't been computed yet
        if (dtnew(ilevel) > 0) then
           factor = min(1/dtnew(ilevel), factor)
        else
           factor = 0
        end if

        do idim = 1, ndim
           vsink(isink, idim) = vsink(isink, idim) - factor*dtnew(ilevel)*vrel_sink(isink, idim, itype)
           v_DF(isink, ilevel, idim, itype) = v_DFall(isink, idim, itype) / mass_DFall(isink, itype)
        end do
     end do
  end do

  ! We remove momentum from cloud particles from this level
  ! Loop over cpus
  do icpu=1, ncpu
     igrid=headl(icpu,ilevel)

     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid

        ! Loop over particles
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1, npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
              !if (tp(ipart).eq.0d0  .and. idp(ipart).lt.0) then
              if (is_cloud(typep(ipart))) then
                 isink = -idp(ipart)
                 do itype = 1, 2
                    b90 = factG * msink(isink) / vrel_sink_norm(isink, itype)**2
                    Rsh = factG * msink(isink) / c**2
                    bmin = max(b90, Rsh)
                    if (itype.eq.1) then
                       CoulombLog = DF_ncells*dx_min / bmin
                       volume = vol_cloud
                    else
                       CoulombLog = DF_ncells*dx_dm / bmin
                       volume = vol_dm
                    end if
                    if (CoulombLog .gt. 1) then
                       factor = 0d0
                       do ii = levelmin, nlevelmax
                          factor = factor + (&
                               &mass_lowspeed(isink, ii, itype)*log(CoulombLog)+&
                               &fact_fast(isink, ii, itype)) / volume
                       end do
                       factor = factor * 4*pi*factG**2*msink(isink)/vrel_sink_norm(isink, itype)**3
                    else
                       factor = 0
                    end if
                    if (dtnew(ilevel) > 0) then
                       factor = min(1/dtnew(ilevel), factor)
                    else
                       factor = 0
                    end if

                    do idim = 1, ndim
                       vp(ipart, idim) = vp(ipart,idim) - factor*dtnew(ilevel) * vrel_sink(isink, idim, itype)
                    end do
                 end do
              end if

              ipart=next_part  ! Go to next particle
           end do
        endif

        igrid=next(igrid)   ! Go to next grid
     end do
     ! End loop over grids

  end do
  ! End loop over cpus

end subroutine get_drag_part

subroutine psy_function(X, rmax2, psy)
  ! This function is used to compute how much mass is deposited w.r.t
  ! to AGN-jet radial distance. It expects X to be in frame of jet with last axis being aligned
  use amr_commons, only : dp
  real(dp), intent(in) :: X(3), rmax2
  real(dp), intent(out) :: psy

  psy = exp(-(X(1)**2 + X(2)**2) / rmax2 / 2._dp)

end subroutine psy_function
