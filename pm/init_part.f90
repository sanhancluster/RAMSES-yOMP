 subroutine init_part
  use amr_commons
  use pm_commons
  use clfind_commons
  use random
  use hydro_commons, only : uold
#ifdef DICE
  ! DICE patch
  use dice_commons
  use cooling_module
  use gadgetreadfilemod
  ! DICE patch
#endif
#ifdef RT
  use rt_parameters,only: convert_birth_times
#endif
  use mpi_mod
  implicit none
  !------------------------------------------------------------
  ! Allocate particle-based arrays.
  ! Read particles positions and velocities from grafic files
  !------------------------------------------------------------
  integer::npart2,ndim2,ncpu2
  integer::ipart,jpart,ipart_old,ilevel,idim
  integer::i,j,igrid,ncache,ngrid,iskip
  integer::ind,ix,iy,iz,ilun,icpu
  integer::i1,i2,i3
  integer::i1_min=0,i1_max=0,i2_min=0,i2_max=0,i3_min=0,i3_max=0
  integer::buf_count,indglob
  real(dp)::dx,xx1,xx2,xx3,vv1,vv2,vv3,mm1
  real(dp)::min_mdm_cpu,min_mdm_all
  real(dp),dimension(1:twotondim,1:3)::xc
  integer ,dimension(1:nvector)::ind_grid,ind_cell,ii
  integer,dimension(1:nvector)::pp
  integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  integer(i8b),allocatable,dimension(:)::isp8
  integer(1),allocatable,dimension(:)::ii1
  real(kind=4),allocatable,dimension(:,:)::init_plane,init_plane_x
  integer(i8b),allocatable,dimension(:,:)::init_plane_id
  real(dp),allocatable,dimension(:,:,:)::init_array,init_array_x
  integer(i8b),allocatable,dimension(:,:,:)::init_array_id
  real(kind=8),dimension(1:nvector,1:3)::xx,vv
  real(kind=8),dimension(1:nvector)::mm
  type(part_t)::tmppart
  real(kind=8)::dispmax=0.0
#ifndef WITHOUTMPI
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer,dimension(1:nvector)::cc
  integer,dimension(MPI_STATUS_SIZE,2*ncpu)::statuses
  integer,dimension(2*ncpu)::reqsend,reqrecv
  integer,dimension(ncpu)::sendbuf,recvbuf
  integer::dummy_io,info,info2,npart_new
  integer::countsend,countrecv
  integer::ibuf,tagu=102
  integer,parameter::tagg=1109,tagg2=1110,tagg3=1111
#endif
  logical::error,keep_part,eof,read_pos=.false.,ok,read_ids=.false.
  character(LEN=1024)::filename,filename_x, filename_id
  character(LEN=80)::fileloc
  character(LEN=20)::filetype_loc
  character(LEN=5)::nchar,ncharcpu
  integer, dimension(1:ncpu,1:IRandNumSize)::allseed

#ifdef DICE
  ! DICE patch
  integer::type_index
  integer::dummy_int,blck_size,jump_blck,blck_cnt,stat,ifile
  integer::head_blck,pos_blck,vel_blck,id_blck,mass_blck,u_blck,metal_blck,age_blck
  integer::head_size,pos_size,vel_size,id_size,mass_size,u_size,metal_size,age_size
  integer::kpart,lpart,mpart,opart,gpart,ngas,nhalo
  !integer, dimension(nvector)::ids
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp),dimension(1:nvector)::tt,zz,uu
  real,dimension(1:nvector,1:3)::xx_sp,vv_sp
  real,dimension(1:nvector)::mm_sp,tt_sp,zz_sp,uu_sp
  real(dp)::mgas_tot
  real::dummy_real,ipbar
  character(LEN=12)::ifile_str
  character(LEN=4)::blck_name
  logical::eob,file_exists,skip
  TYPE(gadgetheadertype)::header
  real(dp)::tyoung
  real(dp),parameter::myr2s=3.1536000d+13
  ! DICE patch
#endif

  ! MC_tracer
  integer :: maxidp, minidp
  if(verbose)write(*,*)'Entering init_part'

  if(myid.eq.1)write(*,*)'WARNING: NEVER USE FAMILY CODES / TAGS > 127.'
  if(myid.eq.1)write(*,*)'See https://bitbucket.org/rteyssie/ramses/wiki/Particle%20Families'

  if(allocated(xp))then
     if(verbose)write(*,*)'Initial conditions already set'
     return
  end if

  ! Allocate particle variables
  allocate(xp    (npartmax,ndim))
  allocate(vp    (npartmax,ndim))
  allocate(mp    (npartmax))
  if (MC_tracer) then
     allocate(tmpp  (npartmax))
     allocate(itmpp (npartmax))
     allocate(partp (npartmax))
     allocate(move_flag(npartmax))
     move_flag = 0
  end if
  allocate(nextp (npartmax))
  allocate(prevp (npartmax))
  allocate(levelp(npartmax))
  allocate(idp   (npartmax))
  allocate(typep (npartmax))
#ifdef OUTPUT_PARTICLE_POTENTIAL
  allocate(ptcl_phi(npartmax))
#endif
#ifdef DICE
  ! DICE patch
  allocate(up(npartmax))
  if(ic_mask_ptype.gt.-1)then
     allocate(maskp(npartmax))
  endif
  ! DICE patch
#endif
  xp=0.0; vp=0.0; mp=0.0; levelp=0; idp=0;
  typep(1:npartmax)%family=FAM_UNDEF; typep(1:npartmax)%tag=0
  if(star.or.sink)then
     allocate(tp(npartmax))
     tp=0.0
     if(write_stellar_densities) then
        allocate(st_n_tp(npartmax))
        st_n_tp=0.0
        allocate(st_n_SN(npartmax))
        st_n_SN=0.0
        allocate(st_e_SN(npartmax))
        st_e_SN=0.0
     endif
     if(metal)then
        allocate(zp(npartmax))
        zp=0.0
     end if
     if(use_initial_mass)then
        allocate(mp0(npartmax))
        mp0=0.0
     endif
  end if

  !--------------------
  ! Read part.tmp file
  !--------------------

  if(nrestart>0)then

     ilun=2*ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/part_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/part_'//TRIM(nchar)//'.out'
     endif

     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)
     ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     open(unit=ilun,file=fileloc,form='unformatted')
     rewind(ilun)
     read(ilun)ncpu2
     read(ilun)ndim2
     read(ilun)npart2
     if (MC_tracer) then
        read(ilun)localseed, tracer_seed
     else
        read(ilun)localseed
     end if

     read(ilun)nstar_tot
     read(ilun)mstar_tot
     read(ilun)mstar_lost
     read(ilun)nsink
     if(ncpu2.ne.ncpu.or.ndim2.ne.ndim.or.npart2.gt.npartmax)then
        write(*,*)'File part.tmp not compatible'
        write(*,*)'Found   =',ncpu2,ndim2,npart2
        write(*,*)'Expected=',ncpu,ndim,npartmax
        call clean_stop
     end if
     ! Read position
     allocate(xdp(1:npart2))
     do idim=1,ndim
        read(ilun)xdp
        xp(1:npart2,idim)=xdp
     end do
     ! Read velocity
     do idim=1,ndim
        read(ilun)xdp
        vp(1:npart2,idim)=xdp
     end do
     ! Read mass
     read(ilun)xdp
     mp(1:npart2)=xdp
     deallocate(xdp)
     ! Read identity
     allocate(isp8(1:npart2))
     read(ilun)isp8
     idp(1:npart2)=isp8
     deallocate(isp8)

     ! Read level
     allocate(isp(1:npart2))
     read(ilun)isp
     levelp(1:npart2)=isp
     deallocate(isp)

     ! Read family
     allocate(ii1(1:npart2))
     read(ilun)ii1
     typep(1:npart2)%family = ii1
     ! Read tag
     read(ilun)ii1
     typep(1:npart2)%tag = ii1
     deallocate(ii1)

#ifdef OUTPUT_PARTICLE_POTENTIAL
     ! We don't need the potential, but read it anyway (to get the records correctly for tp/zp)
     read(ilun)
#endif
     if(star.or.sink)then
        ! Read birth epoch
        allocate(xdp(1:npart2))
        read(ilun)xdp
        tp(1:npart2)=xdp
        if(convert_birth_times) then
           do i = 1, npart2 ! Convert birth time to proper for RT postpr.
              call getProperTime(tp(i),tp(i))
           enddo
        endif
        ! BEGIN SD PATCH----------------------------------------------------!SD
        if(write_stellar_densities) then
           ! Read gas density at birth
           read(ilun)xdp
           st_n_tp(1:npart2) = xdp
           ! Read gas density at SN
           read(ilun)xdp
           st_n_SN(1:npart2) = xdp
           ! Read SN energy injected
           read(ilun)xdp
           st_e_SN(1:npart2) = xdp
        endif
        ! END SD PATCH------------------------------------------------------!SD

        if(metal)then
           ! Read metallicity
           read(ilun)xdp
           zp(1:npart2)=xdp
        end if
        if(use_initial_mass)then
           read(ilun)xdp
           mp0(1:npart2)=xdp
        endif
        deallocate(xdp)
     end if

     if (MC_tracer) then
        ! We need to convert the ids of the star to their local index

        ! Read partp
        allocate(isp(1:npart2))

        ! Now read partp
        read(ilun)isp
        partp(1:npart2) = isp

        minidp = npart2
        maxidp = 0
        do i = 1, npart2
           if (is_star(typep(i))) then
              minidp = min(minidp, idp(i))
              maxidp = max(maxidp, idp(i))
           end if
        end do

        ! We now need to convert partp for star tracers (idp -> local index)
        ! Either the difference between the smallest star id and the largest is
        ! small enough (here, less than 50,000,000 -- that's 50M in memory)...
        if (maxidp - minidp < 50000000) then
           allocate(isp8(minidp:maxidp))
           isp8(:) = -1
           do i = 1, npart2
              if (is_star(typep(i))) then
                 isp8(idp(i)) = i
              end if
           end do

           ok = .true.
           do i = 1, npart2
              if (is_star_tracer(typep(i))) then
                 if (isp8(partp(i)) == -1) then
                    write(*, *) 'An error occured while loading star tracers. Aborting.'
                    stop 1
                 end if
                 partp(i) = isp8(partp(i))
              end if
           end do

           deallocate(isp8)
        else
        ! ... or for each tracers we loop on *all* the particles but
        ! it costs 0 memory and is only runned once.)
           do i = 1, npart2
              ! Get star tracers
              if (is_star_tracer(typep(i))) then
                 star_loop: do j = 1, npart2
                    if (is_star(typep(j))) then
                       ! Check that star's id == tracer partp
                       if (partp(i) == idp(j)) then
                          partp(i) = j
                          exit star_loop
                       end if
                    end if
                 end do star_loop
                 if (.not. is_star(typep(partp(i)))) then
                    write(*, *) 'An error occured while loading star tracers. Aborting.'
                    stop 1
                 end if
              end if
           end do
        end if

        deallocate(isp)
        ! End MC Patch
     end if
     close(ilun)

     ! Send the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
           dummy_io=1
           call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg, &
                & MPI_COMM_WORLD,info2)
        end if
     endif
#endif
     ! Get nlevelmax_part from cosmological inital conditions
     if(cosmo)then
        min_mdm_cpu = 1.0
        do ipart=1,npart2
           ! Get dark matter only
           if (is_DM(typep(ipart))) then
              ! note: using two nested if so that the second one is only evaluated for DM particles
              if (mp(ipart) .lt. min_mdm_cpu) min_mdm_cpu = mp(ipart)
           end if
        end do

#ifndef WITHOUTMPI
        call MPI_ALLREDUCE(min_mdm_cpu,min_mdm_all,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
#else
        min_mdm_all = min_mdm_cpu
#endif
        ilevel = 1
        do while(.true.)
           mm1 = 0.5d0**(3*ilevel)*(1.0d0-omega_b/omega_m)
           if((mm1.GT.0.90*min_mdm_all).AND.(mm1.LT.1.10*min_mdm_all))then
              nlevelmax_part = ilevel
              exit
           endif
           ilevel = ilevel+1
        enddo
        if(myid==1) write(*,*) 'nlevelmax_part=',nlevelmax_part
     endif

     if(debug)write(*,*)'part.tmp read for processor ',myid
     npart=npart2

     if (tracer .and. MC_tracer) then
        ! Attempt to read mass from binary file
        if (myid == 1) then
           if (trim(tracer_feed_fmt) == 'binary' .and. tracer_mass < 0) then
              open(unit=10, file=trim(tracer_feed), form='unformatted', status='old')
              read(10) ! ntot
              read(10) tracer_mass
              close(10)
           end if
        end if
        ! Broadcast to all CPUs the value of the tracer mass
#ifndef WITHOUTMPI
        call MPI_BCAST(tracer_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif
        if (myid == 1) write(*, *) 'Using a tracer mass of ', tracer_mass
     end if

  else

     filetype_loc=filetype
#ifdef DICE
     if(.not. cosmo)filetype_loc='dice'
#else
     if(.not. cosmo)filetype_loc='ascii'
#endif

     select case (filetype_loc)

     case ('grafic')
        call load_grafic
     case ('ascii')
        call load_ascii
#ifdef DICE
     case ('dice')
        call load_dice
#endif
     case ('gadget')
        call load_gadget

     case DEFAULT
        write(*,*) 'Unsupported format file ' // filetype
        call clean_stop

     end select

     ! Initialize tracer particles
     if(tracer) then
        if(tracer_seed(1)==-1)then
           call rans(ncpu, tseed, allseed)
           tracer_seed = allseed(myid, 1:IRandNumSize)
        end if

        if (trim(tracer_feed_fmt) == 'binary') then
           call load_tracers_bin(1)
        else if (trim(tracer_feed_fmt) == 'binary2') then
           call load_tracers_bin(2)
        else if (trim(tracer_feed_fmt) == 'inplace') then
           call load_tracers_inplace
        else if (trim(tracer_feed_fmt) == 'ascii') then
           call load_tracers
        else
           write(*, '(a,a,a)')'Data input format not understood: "', (tracer_feed_fmt), '"'
           stop
        end if

        ! Reset first balance flags
        tracer_first_balance_levelmin = nlevelmax + 1
     end if

  end if

  if(sink)call init_sink


contains

  subroutine load_grafic
    ! Read data in the grafic format. The particle type is derived
    ! following conversion rules (see pm_commons:props2type)
    ! Grafic format for Ramses assumes the following unit for particles:
    ! - Header lengths (boxszie, pixel size...) in comoving Mpc
    ! - Velocities in proper km s**-1 (file ic_velc*)
    ! - Displacements from cell centers in comoving Mpc/h
    !    (file ic_posc* if present, if not generated through Zeldovich approximation)
    ! - Ids in int or long int
    !    (file ic_particle_ids if present, if not generated internally)

    !----------------------------------------------------
    ! Reading initial conditions GRAFIC2 multigrid arrays
    !----------------------------------------------------
    ipart=0
    ! Loop over initial condition levels
    do ilevel=levelmin,nlevelmax

       if(initfile(ilevel)==' ')cycle

       ! Mesh size at level ilevel in coarse cell units
       dx=0.5D0**ilevel

       ! Set position of cell centers relative to grid center
       do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
          if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
          if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
       end do

       !--------------------------------------------------------------
       ! First step: compute level boundaries and particle positions
       !--------------------------------------------------------------
       i1_min=n1(ilevel)+1; i1_max=0
       i2_min=n2(ilevel)+1; i2_max=0
       i3_min=n3(ilevel)+1; i3_max=0
       ipart_old=ipart

       ! Loop over grids by vector sweeps
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
             do i=1,ngrid
                xx1=xg(ind_grid(i),1)+xc(ind,1)
                xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                xx2=xg(ind_grid(i),2)+xc(ind,2)
                xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                xx3=xg(ind_grid(i),3)+xc(ind,3)
                xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                i1_min=MIN(i1_min,int(xx1)+1)
                i1_max=MAX(i1_max,int(xx1)+1)
                i2_min=MIN(i2_min,int(xx2)+1)
                i2_max=MAX(i2_max,int(xx2)+1)
                i3_min=MIN(i3_min,int(xx3)+1)
                i3_max=MAX(i3_max,int(xx3)+1)
                keep_part=son(ind_cell(i))==0
                if(keep_part)then
                   ipart=ipart+1
                   if(ipart>npartmax)then
                      write(*,*)'Maximum number of particles incorrect'
                      write(*,*)'npartmax should be greater than',ipart
                      stop
                   endif
                   if(ndim>0)xp(ipart,1)=xg(ind_grid(i),1)+xc(ind,1)
                   if(ndim>1)xp(ipart,2)=xg(ind_grid(i),2)+xc(ind,2)
                   if(ndim>2)xp(ipart,3)=xg(ind_grid(i),3)+xc(ind,3)
                   mp(ipart)=0.5d0**(3*ilevel)*(1.0d0-omega_b/omega_m)
                end if
             end do
          end do
          ! End loop over cells
       end do
       ! End loop over grids

       ! Check that all grids are within initial condition region
       error=.false.
       if(active(ilevel)%ngrid>0)then
          if(i1_min<1.or.i1_max>n1(ilevel))error=.true.
          if(i2_min<1.or.i2_max>n2(ilevel))error=.true.
          if(i3_min<1.or.i3_max>n3(ilevel))error=.true.
       end if
       if(error) then
          write(*,*)'Some grid are outside initial conditions sub-volume'
          write(*,*)'for ilevel=',ilevel
          write(*,*)i1_min,i1_max
          write(*,*)i2_min,i2_max
          write(*,*)i3_min,i3_max
          write(*,*)n1(ilevel),n2(ilevel),n3(ilevel)
          call clean_stop
       end if
       if(debug)then
          write(*,*)myid,i1_min,i1_max,i2_min,i2_max,i3_min,i3_max
       endif

       !---------------------------------------------------------------------
       ! Second step: read initial condition file and set particle velocities
       !---------------------------------------------------------------------
       ! Allocate initial conditions array
       if(active(ilevel)%ngrid>0)then
          allocate(init_array(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
          allocate(init_array_x(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
          init_array=0d0
          init_array_x=0d0
       end if
       allocate(init_plane(1:n1(ilevel),1:n2(ilevel)))
       allocate(init_plane_x(1:n1(ilevel),1:n2(ilevel)))

       filename_id=TRIM(initfile(ilevel))//'/ic_particle_ids'
       INQUIRE(file=filename_id,exist=read_ids)
       if(read_ids) then
         if(myid==1)write(*,*)'Reading particle ids from file '//TRIM(filename_id)
         allocate(init_plane_id(1:n1(ilevel),1:n2(ilevel)))
         allocate(init_array_id(i1_min:i1_max,i2_min:i2_max,i3_min:i3_max))
       end if

       ! Loop over input variables
       do idim=1,ndim

          ! Read dark matter initial displacement field
          if(multiple)then
             call title(myid,nchar)
             if(idim==1)filename=TRIM(initfile(ilevel))//'/dir_velcx/ic_velcx.'//TRIM(nchar)
             if(idim==2)filename=TRIM(initfile(ilevel))//'/dir_velcy/ic_velcy.'//TRIM(nchar)
             if(idim==3)filename=TRIM(initfile(ilevel))//'/dir_velcz/ic_velcz.'//TRIM(nchar)
          else
             if(idim==1)filename=TRIM(initfile(ilevel))//'/ic_velcx'
             if(idim==2)filename=TRIM(initfile(ilevel))//'/ic_velcy'
             if(idim==3)filename=TRIM(initfile(ilevel))//'/ic_velcz'

             if(idim==1)filename_x=TRIM(initfile(ilevel))//'/ic_poscx'
             if(idim==2)filename_x=TRIM(initfile(ilevel))//'/ic_poscy'
             if(idim==3)filename_x=TRIM(initfile(ilevel))//'/ic_poscz'

             INQUIRE(file=filename_x,exist=ok)
             if(.not.ok)then
                read_pos = .false.
             else
                read_pos = .true.
                if(myid==1)write(*,*)'Reading file '//TRIM(filename_x)
             end if

          endif

          if(myid==1)write(*,*)'Reading file '//TRIM(filename)

          if(multiple)then
             ilun=myid+10
             ! Wait for the token
#ifndef WITHOUTMPI
             if(IOGROUPSIZE>0) then
                if (mod(myid-1,IOGROUPSIZE)/=0) then
                   call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tagg2,&
                        & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
                end if
             endif
#endif
             open(ilun,file=filename,form='unformatted')
             rewind ilun
             read(ilun) ! skip first line
             do i3=1,n3(ilevel)
                read(ilun)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                if(active(ilevel)%ngrid>0)then
                   if(i3.ge.i3_min.and.i3.le.i3_max)then
                      init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                           & init_plane(i1_min:i1_max,i2_min:i2_max)
                   end if
                endif
             end do
             close(ilun)
             ! Send the token
#ifndef WITHOUTMPI
             if(IOGROUPSIZE>0) then
                if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
                   dummy_io=1
                   call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tagg2, &
                        & MPI_COMM_WORLD,info2)
                end if
             endif
#endif

          else
             if(myid==1)then
                open(10,file=filename,form='unformatted')
                rewind 10
                read(10) ! skip first line
             end if
             do i3=1,n3(ilevel)
                if(myid==1)then
                   if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                   read(10)((init_plane(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                else
                   init_plane=0.0
                endif
                buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                call MPI_BCAST(init_plane,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif

                if(active(ilevel)%ngrid>0)then
                   if(i3.ge.i3_min.and.i3.le.i3_max)then
                      init_array(i1_min:i1_max,i2_min:i2_max,i3) = &
                           & init_plane(i1_min:i1_max,i2_min:i2_max)
                   end if
                endif
             end do
             if(myid==1)close(10)

             if(read_pos) then
                if(myid==1)then
                   open(10,file=filename_x,form='unformatted')
                   rewind 10
                   read(10) ! skip first line
                end if
                do i3=1,n3(ilevel)
                   if(myid==1)then
                      if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                      read(10)((init_plane_x(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                   else
                      init_plane_x=0.0
                   endif
                   buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
                   call MPI_BCAST(init_plane_x,buf_count,MPI_REAL,0,MPI_COMM_WORLD,info)
#endif
                   if(active(ilevel)%ngrid>0)then
                      if(i3.ge.i3_min.and.i3.le.i3_max)then
                         init_array_x(i1_min:i1_max,i2_min:i2_max,i3) = &
                              & init_plane_x(i1_min:i1_max,i2_min:i2_max)
                      end if
                   endif
                end do
                if(myid==1)close(10)
             end if

             if(read_ids) then
                if(myid==1)then
                   open(10,file=filename_id,form='unformatted')
                   rewind 10
                   read(10) ! skip first line
                end if
                do i3=1,n3(ilevel)
                   if(myid==1)then
                      if(debug.and.mod(i3,10)==0)write(*,*)'Reading plane ',i3
                      read(10)((init_plane_id(i1,i2),i1=1,n1(ilevel)),i2=1,n2(ilevel))
                   else
                      init_plane_id=0
                   endif
                   buf_count=n1(ilevel)*n2(ilevel)
#ifndef WITHOUTMPI
#ifndef LONGINT
                   call MPI_BCAST(init_plane_id,buf_count,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#else
                   call MPI_BCAST(init_plane_id,buf_count,MPI_INTEGER8,0,MPI_COMM_WORLD,info)
#endif
#endif
                   if(active(ilevel)%ngrid>0)then
                      if(i3.ge.i3_min.and.i3.le.i3_max)then
                         init_array_id(i1_min:i1_max,i2_min:i2_max,i3) = &
                              & init_plane_id(i1_min:i1_max,i2_min:i2_max)
                      end if
                   endif
                end do
                if(myid==1)close(10)
              end if
          endif

          if(active(ilevel)%ngrid>0)then
             ! Rescale initial displacement field to code units
             init_array=dfact(ilevel)*dx/dxini(ilevel)*init_array/vfact(ilevel)
             if(read_pos)then
                init_array_x = init_array_x/boxlen_ini
             endif
             ! Loop over grids by vector sweeps
             ipart=ipart_old
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
                   do i=1,ngrid
                      xx1=xg(ind_grid(i),1)+xc(ind,1)
                      xx1=(xx1*(dxini(ilevel)/dx)-xoff1(ilevel))/dxini(ilevel)
                      xx2=xg(ind_grid(i),2)+xc(ind,2)
                      xx2=(xx2*(dxini(ilevel)/dx)-xoff2(ilevel))/dxini(ilevel)
                      xx3=xg(ind_grid(i),3)+xc(ind,3)
                      xx3=(xx3*(dxini(ilevel)/dx)-xoff3(ilevel))/dxini(ilevel)
                      i1=int(xx1)+1
                      i1=int(xx1)+1
                      i2=int(xx2)+1
                      i2=int(xx2)+1
                      i3=int(xx3)+1
                      i3=int(xx3)+1
                      keep_part=son(ind_cell(i))==0
                      if(keep_part)then
                         ipart=ipart+1
                         vp(ipart,idim)=init_array(i1,i2,i3)
                         if(.not. read_pos)then
                            dispmax=max(dispmax,abs(init_array(i1,i2,i3)/dx))
                         else
                            xp(ipart,idim)=xg(ind_grid(i),idim)+xc(ind,idim)+init_array_x(i1,i2,i3)
                            dispmax=max(dispmax,abs(init_array_x(i1,i2,i3)/dx))
                          if (read_ids) then
                            idp(ipart) = init_array_id(i1,i2,i3)
                          end if
                         endif
                      end if
                   end do
                end do
                ! End loop over cells
             end do
             ! End loop over grids
          endif

       end do
       ! End loop over input variables

       ! Deallocate initial conditions array
       if(active(ilevel)%ngrid>0)then
          deallocate(init_array,init_array_x)
       end if
       deallocate(init_plane,init_plane_x)

       if(read_ids) then
         deallocate(init_plane_id)
         deallocate(init_array_id)
       end if


       if(debug)write(*,*)'npart=',ipart,'/',npartmax,' for PE=',myid

    end do
    ! End loop over levels

    ! Initial particle number
    npart=ipart

    ! Move particle according to Zeldovich approximation
    if(.not. read_pos)then
       xp(1:npart,1:ndim)=xp(1:npart,1:ndim)+vp(1:npart,1:ndim)
    endif

    ! Scale displacement to velocity
    vp(1:npart,1:ndim)=vfact(1)*vp(1:npart,1:ndim)

    ! Periodic box
    do ipart=1,npart
#if NDIM>0
       if(xp(ipart,1)<  0.0d0  )xp(ipart,1)=xp(ipart,1)+dble(nx)
       if(xp(ipart,1)>=dble(nx))xp(ipart,1)=xp(ipart,1)-dble(nx)
#endif
#if NDIM>1
       if(xp(ipart,2)<  0.0d0  )xp(ipart,2)=xp(ipart,2)+dble(ny)
       if(xp(ipart,2)>=dble(ny))xp(ipart,2)=xp(ipart,2)-dble(ny)
#endif
#if NDIM>2
       if(xp(ipart,3)<  0.0d0  )xp(ipart,3)=xp(ipart,3)+dble(nz)
       if(xp(ipart,3)>=dble(nz))xp(ipart,3)=xp(ipart,3)-dble(nz)
#endif
    end do

#ifndef WITHOUTMPI
    ! Compute particle Hilbert ordering
    sendbuf=0
    do ipart=1,npart
       xx(1,1:3)=xp(ipart,1:3)
       xx_dp(1,1:3)=xx(1,1:3)
       call cmp_cpumap(xx_dp,cc,1)
       if(cc(1).ne.myid)sendbuf(cc(1))=sendbuf(cc(1))+1
    end do

    ! Allocate communication buffer in emission
    do icpu=1,ncpu
       ncache=sendbuf(icpu)
       if(ncache>0)then
          allocate(emission(icpu,1)%up(1:ncache,1:twondim+1))
          allocate(emission(icpu,1)%fp(1:ncache,1:2))
       end if
    end do

    ! Fill communicators
    jpart=0
    sendbuf=0
    do ipart=1,npart
       xx(1,1:3)=xp(ipart,1:3)
       xx_dp(1,1:3)=xx(1,1:3)
       call cmp_cpumap(xx_dp,cc,1)
       if(cc(1).ne.myid)then
          icpu=cc(1)
          sendbuf(icpu)=sendbuf(icpu)+1
          ibuf=sendbuf(icpu)
          emission(icpu,1)%up(ibuf,1)=xp(ipart,1)
          emission(icpu,1)%up(ibuf,2)=xp(ipart,2)
          emission(icpu,1)%up(ibuf,3)=xp(ipart,3)
          emission(icpu,1)%up(ibuf,4)=vp(ipart,1)
          emission(icpu,1)%up(ibuf,5)=vp(ipart,2)
          emission(icpu,1)%up(ibuf,6)=vp(ipart,3)
          emission(icpu,1)%up(ibuf,7)=mp(ipart)
          emission(icpu,1)%fp(ibuf,1)=part2int(typep(ipart))
          emission(icpu,1)%fp(ibuf,2)=idp(ipart)
       else
          jpart=jpart+1
          xp(jpart,1:3)=xp(ipart,1:3)
          vp(jpart,1:3)=vp(ipart,1:3)
          mp(jpart)    =mp(ipart)
          idp(jpart)   =idp(ipart)
       endif
    end do

    ! Communicate virtual particle number to parent cpu
    call MPI_ALLTOALL(sendbuf,1,MPI_INTEGER,recvbuf,1,MPI_INTEGER,MPI_COMM_WORLD,info)

    ! Compute total number of newly created particles
    npart_new=0
    do icpu=1,ncpu
       npart_new=npart_new+recvbuf(icpu)
    end do

    if(jpart+npart_new.gt.npartmax)then
       write(*,*)'No more free memory for particles'
       write(*,*)'Increase npartmax'
       write(*,*)myid
       write(*,*)jpart,npart_new
       write(*,*)bound_key
       call MPI_ABORT(MPI_COMM_WORLD,1,info)
    end if

    ! Allocate communication buffer in reception
    do icpu=1,ncpu
       ncache=recvbuf(icpu)
       if(ncache>0)then
          allocate(reception(icpu,1)%up(1:ncache,1:twondim+1))
          allocate(reception(icpu,1)%fp(1:ncache,1:2))
       end if
    end do

    ! Taking care of real values
    ! Receive particles
    countrecv=0
    do icpu=1,ncpu
       ncache=recvbuf(icpu)
       if(ncache>0)then
          buf_count=ncache*(twondim+1)
          countrecv=countrecv+1
          call MPI_IRECV(reception(icpu,1)%up,buf_count, &
               & MPI_DOUBLE_PRECISION,icpu-1,&
               & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
       end if
    end do

    ! Send particles
    countsend=0
    do icpu=1,ncpu
       ncache=sendbuf(icpu)
       if(ncache>0)then
          buf_count=ncache*(twondim+1)
          countsend=countsend+1
          call MPI_ISEND(emission(icpu,1)%up,buf_count, &
               & MPI_DOUBLE_PRECISION,icpu-1,&
               & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
       end if
    end do

    ! Wait for full completion of receives
    call MPI_WAITALL(countrecv,reqrecv,statuses,info)

    ! Wait for full completion of sends
    call MPI_WAITALL(countsend,reqsend,statuses,info)

    ! Taking care of int values
    ! Receive particles
    countrecv=0
    do icpu=1,ncpu
       ncache=recvbuf(icpu)
       if(ncache>0)then
          buf_count=ncache * 2
          countrecv=countrecv+1
#ifndef LONGINT
          call MPI_IRECV(reception(icpu,1)%fp,buf_count, &
                & MPI_INTEGER,icpu-1,&
                & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
#else
          call MPI_IRECV(reception(icpu,1)%fp,buf_count, &
                & MPI_INTEGER8,icpu-1,&
                & tagu,MPI_COMM_WORLD,reqrecv(countrecv),info)
#endif

       end if
    end do

    ! Send particles
    countsend=0
    do icpu=1,ncpu
       ncache=sendbuf(icpu)
       if(ncache>0)then
          buf_count=ncache * 2
          countsend=countsend+1
#ifndef LONGINT
                    call MPI_ISEND(emission(icpu,1)%fp,buf_count, &
                          & MPI_INTEGER,icpu-1,&
                          & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
#else
                    call MPI_ISEND(emission(icpu,1)%fp,buf_count, &
                          & MPI_INTEGER8,icpu-1,&
                          & tagu,MPI_COMM_WORLD,reqsend(countsend),info)
#endif
       end if
    end do


    ! Wait for full completion of receives
    call MPI_WAITALL(countrecv,reqrecv,statuses,info)

    ! Wait for full completion of sends
    call MPI_WAITALL(countsend,reqsend,statuses,info)

    ! Create new particles
    do icpu=1,ncpu
       do ibuf=1,recvbuf(icpu)
          jpart=jpart+1
          xp(jpart,1)=reception(icpu,1)%up(ibuf,1)
          xp(jpart,2)=reception(icpu,1)%up(ibuf,2)
          xp(jpart,3)=reception(icpu,1)%up(ibuf,3)
          vp(jpart,1)=reception(icpu,1)%up(ibuf,4)
          vp(jpart,2)=reception(icpu,1)%up(ibuf,5)
          vp(jpart,3)=reception(icpu,1)%up(ibuf,6)
          mp(jpart)  =reception(icpu,1)%up(ibuf,7)
          idp(jpart) =reception(icpu,1)%fp(ibuf,2)
       end do
    end do

    ! Erase old particles
    do ipart=jpart+1,npart
       xp(ipart,1)=0d0
       xp(ipart,2)=0d0
       xp(ipart,3)=0d0
       vp(ipart,1)=0d0
       vp(ipart,2)=0d0
       vp(ipart,3)=0d0
       mp(ipart)=0d0
       idp(ipart)=0
    end do

    npart=jpart

    ! Deallocate communicators
    do icpu=1,ncpu
       if(sendbuf(icpu)>0) then
        deallocate(emission(icpu,1)%up)
        deallocate(emission(icpu,1)%fp)
       end if

       if(recvbuf(icpu)>0)then
         deallocate(reception(icpu,1)%up)
         deallocate(reception(icpu,1)%fp)
       end if
    end do

    write(*,*)'npart=',ipart,'/',npartmax,' for PE=',myid
#endif

    ! Compute particle initial level
    do ipart=1,npart
       levelp(ipart)=levelmin
    end do

    ! Setup DM for all particles
    do ipart=1, npart
       typep(ipart)%family = FAM_DM
       typep(ipart)%tag = 0
    end do

    ! Compute particle initial age and metallicity
    if(star.or.sink)then
       do ipart=1,npart
          tp(ipart)=0d0
          if(metal)then
             zp(ipart)=0d0
          end if
       end do
    end if

    ! Compute particle initial identity
    if(.not.read_ids) then
      npart_cpu=0; npart_all=0
      npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
      call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
      call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
      npart_cpu(1)=npart_all(1)
#endif
      do icpu=2,ncpu
         npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
      end do
      if(myid==1)then
         do ipart=1,npart
            idp(ipart)=ipart
         end do
      else
         do ipart=1,npart
            idp(ipart)=npart_cpu(myid-1)+ipart
         end do
      end if
    end if

  end subroutine load_grafic

  subroutine load_ascii
    ! This function load from ASCII file. As is, you can only load dark matter particles
    ! Local particle count
    ipart=0

    if(trim(initfile(levelmin)).ne.' ')then

       filename=TRIM(initfile(levelmin))//'/ic_part'
       if(myid==1)then
          open(10,file=filename,form='formatted')
          indglob=0
       end if
       eof=.false.

       do while (.not.eof)
          xx=0.0
          if(myid==1)then
             jpart=0
             do i=1,nvector
                read(10,*,end=100)xx1,xx2,xx3,vv1,vv2,vv3,mm1
                jpart=jpart+1
                indglob=indglob+1
                xx(i,1)=xx1+boxlen/2.0
                xx(i,2)=xx2+boxlen/2.0
                xx(i,3)=xx3+boxlen/2.0
                vv(i,1)=vv1
                vv(i,2)=vv2
                vv(i,3)=vv3
                mm(i  )=mm1
                ii(i  )=indglob
                tmppart%family = FAM_DM
                tmppart%tag    = 0
                pp(i  )=part2int(tmppart)
             end do
100          continue
             if(jpart<nvector)eof=.true.
          endif
          buf_count=nvector*3
#ifndef WITHOUTMPI
          call MPI_BCAST(xx,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(vv,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(mm,nvector  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(ii,nvector  ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(pp,nvector  ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(eof,1       ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(jpart,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call cmp_cpumap(xx,cc,jpart)
#endif

          do i=1,jpart
#ifndef WITHOUTMPI
             if(cc(i)==myid)then
#endif
                ipart=ipart+1
                if(ipart>npartmax)then
                   write(*,*)'Maximum number of particles incorrect'
                   write(*,*)'npartmax should be greater than',ipart
                   stop
                endif
                xp(ipart,1:3)= xx(i,1:3)
                vp(ipart,1:3)= vv(i,1:3)
                mp(ipart)    = mm(i)
                levelp(ipart)= levelmin
                idp(ipart)   = ii(i)
                ! Get back the particle type from the communicated
                ! shortened integer
                typep(ipart) = int2part(pp(i))
                if(use_initial_mass) mp0(ipart) = mm(i)
#ifndef WITHOUTMPI
             endif
#endif
          enddo

       end do
       if(myid==1)close(10)

    end if
    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    npart_cpu(1)=npart_all(1)
#endif
    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do
    if(debug)write(*,*)'npart=',npart,'/',npart_cpu(ncpu)
  end subroutine load_ascii

#ifdef DICE
  subroutine load_dice

    dice_init=.true.
    ! Conversion factor from user units to cgs units
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_m = scale_d*scale_l**3
    ! Reading header of the Gadget file
    error=.false.
    ipart    = 0

    !
    if (use_proper_time)then
       tyoung = t_sne*myr2s/(scale_t/aexp**2)
    else
       tyoung = t_sne*myr2s/scale_t
    endif
    tyoung = ic_t_restart - tyoung


    do ifile=1,ic_nfile
       write(ifile_str,*) ifile
       if(ic_nfile.eq.1) then
          filename=TRIM(initfile(levelmin))//'/'//TRIM(ic_file)
       else
          filename=TRIM(initfile(levelmin))//'/'//TRIM(ic_file)//'.'//ADJUSTL(ifile_str)
       endif
       INQUIRE(FILE=filename,EXIST=file_exists)
       if(.not.file_exists) then
          if(myid==1) write(*,*) TRIM(filename)," not found"
          call clean_stop
       endif
       if(myid==1)then
          write(*,'(A12,A)') " Opening -> ",trim(filename)
          if((ic_format.ne.'Gadget1').and.(ic_format.ne.'Gadget2')) then
             if(myid==1) write(*,*) 'Specify a valid IC file format [ic_format=Gadget1/Gadget2]'
             error=.true.
          endif
          OPEN(unit=1,file=filename,status='old',action='read',form='unformatted',access="stream")
          ! Init block address
          head_blck  = -1
          pos_blck   = -1
          vel_blck   = -1
          id_blck    = -1
          u_blck     = -1
          mass_blck  = -1
          metal_blck = -1
          age_blck   = -1

          if(ic_format .eq. 'Gadget1') then
             ! Init block counter
             jump_blck = 1
             blck_cnt = 1
             do while(.true.)
                ! Reading data block header
                read(1,POS=jump_blck,iostat=stat) blck_size
                if(stat /= 0) exit
                ! Saving data block positions
                if(blck_cnt .eq. 1) then
                   head_blck  = int(jump_blck+sizeof(blck_size))
                   head_size  = blck_size
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 2) then
                   pos_blck   = int(jump_blck+sizeof(blck_size))
                   pos_size   = int(blck_size/(3*sizeof(dummy_real)))
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 3) then
                   vel_blck   = int(jump_blck+sizeof(blck_size))
                   vel_size   = int(blck_size/(3*sizeof(dummy_real)))
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 4) then
                   id_blck    = int(jump_blck+sizeof(blck_size))
                   id_size    = int(blck_size/sizeof(dummy_int))
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 5) then
                   u_blck     = int(jump_blck+sizeof(blck_size))
                   u_size     = int(blck_size/sizeof(dummy_real))
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 6) then
                   mass_blck  = int(jump_blck+sizeof(blck_size))
                   mass_size  = int(blck_size/sizeof(dummy_real))
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 7) then
                   metal_blck = int(jump_blck+sizeof(blck_size))
                   metal_size = int(blck_size/sizeof(dummy_real))
                   write(*,*)blck_cnt,blck_size
                endif
                if(blck_cnt .eq. 8) then
                   age_blck   = int(jump_blck+sizeof(blck_size))
                   age_size   = int(blck_size/sizeof(dummy_real))
                   write(*,*)blck_cnt,blck_size
                endif
                jump_blck = int(jump_blck+blck_size+2*sizeof(dummy_int))
                blck_cnt = blck_cnt+1
             enddo
          endif

          if(ic_format .eq. 'Gadget2') then
             ! Init block counter
             jump_blck = 1
             write(*,'(A50)')"__________________________________________________"
             do while(.true.)
                ! Reading data block header
                read(1,POS=jump_blck,iostat=stat) dummy_int
                if(stat /= 0) exit
                read(1,POS=jump_blck+sizeof(dummy_int),iostat=stat) blck_name
                if(stat /= 0) exit
                read(1,POS=jump_blck+sizeof(dummy_int)+sizeof(blck_name),iostat=stat) dummy_int
                if(stat /= 0) exit
                read(1,POS=jump_blck+2*sizeof(dummy_int)+sizeof(blck_name),iostat=stat) dummy_int
                if(stat /= 0) exit
                read(1,POS=jump_blck+3*sizeof(dummy_int)+sizeof(blck_name),iostat=stat) blck_size
                if(stat /= 0) exit
                ! Saving data block positions
                if(blck_name .eq. ic_head_name) then
                   head_blck  = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   head_size  = blck_size
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_pos_name) then
                   pos_blck   = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   pos_size   = int(blck_size/(3*sizeof(dummy_real)))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_vel_name) then
                   vel_blck  = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   vel_size  = int(blck_size/(3*sizeof(dummy_real)))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_id_name) then
                   id_blck    = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   id_size    = int(blck_size/sizeof(dummy_int))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_mass_name) then
                   mass_blck  = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   mass_size  = int(blck_size/sizeof(dummy_real))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_u_name) then
                   u_blck     = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   u_size     = int(blck_size/sizeof(dummy_real))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_metal_name) then
                   metal_blck = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   metal_size = int(blck_size/sizeof(dummy_real))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                if(blck_name .eq. ic_age_name) then
                   age_blck   = int(jump_blck+sizeof(blck_name)+4*sizeof(dummy_int))
                   age_size   = int(blck_size/sizeof(dummy_real))
                   write(*,*) '-> Found ',blck_name,' block'
                endif
                jump_blck = int(jump_blck+blck_size+sizeof(blck_name)+5*sizeof(dummy_int))
             enddo
          endif

          if((head_blck.eq.-1).or.(pos_blck.eq.-1).or.(vel_blck.eq.-1)) then
             write(*,*) 'Gadget file does not contain handful data'
             error=.true.
          endif
          if(head_size.ne.256) then
             write(*,*) 'Gadget header is not 256 bytes'
             error=.true.
          endif

          ! Byte swapping doesn't appear to work if you just do READ(1)header
          READ(1,POS=head_blck) header%npart,header%mass,header%time,header%redshift, &
               header%flag_sfr,header%flag_feedback,header%nparttotal, &
               header%flag_cooling,header%numfiles,header%boxsize, &
               header%omega0,header%omegalambda,header%hubbleparam, &
               header%flag_stellarage,header%flag_metals,header%totalhighword, &
               header%flag_entropy_instead_u, header%flag_doubleprecision, &
               header%flag_ic_info, header%lpt_scalingfactor

          nstar_tot = sum(header%npart(3:5))
          npart     = sum(header%npart)
          ngas      = header%npart(1)
          nhalo     = header%npart(2)
          if(cosmo) T2_start = 1.356d-2/aexp**2

          write(*,'(A50)')"__________________________________________________"
          write(*,*)"Found ",npart," particles"
          skip=.false.
          do j=1,6
             if(ic_skip_type(j).eq.0) skip=.true.
          enddo
          if(.not.skip) write(*,*)"----> ",header%npart(1)," type 0 particles with header mass ",header%mass(1)
          skip=.false.
          do j=1,6
             if(ic_skip_type(j).eq.1) skip=.true.
          enddo
          if(.not.skip) write(*,*)"----> ",header%npart(2)," type 1 particles with header mass ",header%mass(2)
          skip=.false.
          do j=1,6
             if(ic_skip_type(j).eq.2) skip=.true.
          enddo
          if(.not.skip) write(*,*)"----> ",header%npart(3)," type 2 particles with header mass ",header%mass(3)
          skip=.false.
          do j=1,6
             if(ic_skip_type(j).eq.3) skip=.true.
          enddo
          if(.not.skip) write(*,*)"----> ",header%npart(4)," type 3 particles with header mass ",header%mass(4)
          skip=.false.
          do j=1,6
             if(ic_skip_type(j).eq.4) skip=.true.
          enddo
          if(.not.skip) write(*,*)"----> ",header%npart(5)," type 4 particles with header mass ",header%mass(5)
          skip=.false.
          do j=1,6
             if(ic_skip_type(j).eq.5) skip=.true.
          enddo
          if(.not.skip) write(*,*)"----> ",header%npart(6)," type 5 particles with header mass ",header%mass(6)

          write(*,'(A50)')"_____________________progress_____________________"
          if((pos_size.ne.npart).or.(vel_size.ne.npart)) then
             write(*,*) 'POS =',pos_size
             write(*,*) 'VEL =',vel_size
             write(*,*) 'Number of particles does not correspond to block sizes'
             error=.true.
          endif

       endif
       if(error) call clean_stop
#ifndef WITHOUTMPI
       call MPI_BCAST(nstar_tot,1,MPI_INTEGER,0,MPI_COMM_WORLD,info)
#endif
       eob      = .false.
       kpart    = 0
       lpart    = 0
       mpart    = 0
       gpart    = 0
       opart    = 0
       mgas_tot = 0.
       ipbar    = 0.
       do while(.not.eob)
          xx=0.
          vv=0.
          ii=0
          mm=0.
          tt=0.
          zz=0.
          uu=0.
          if(myid==1)then
             jpart=0
             do i=1,nvector
                jpart=jpart+1

                ! All particles counter
                kpart=kpart+1
                if(kpart.le.header%npart(1)) type_index = 1
                do j=1,5
                   if(kpart.gt.sum(header%npart(1:j)).and.kpart.le.sum(header%npart(1:j+1))) type_index = j+1
                enddo
                if((sum(header%npart(3:5)).gt.0).and.(kpart.gt.(header%npart(1)+header%npart(2)))) mpart=mpart+1
                if(type_index.ne.2) gpart=gpart+1

                ! Reading Gadget1 or Gadget2 file line-by-line
                ! Mandatory data
                read(1,POS=pos_blck+3*sizeof(dummy_real)*(kpart-1)) xx_sp(i,1:3)
                read(1,POS=vel_blck+3*sizeof(dummy_real)*(kpart-1)) vv_sp(i,1:3)
                if(header%mass(type_index).gt.0) then
                   mm_sp(i) = real(header%mass(type_index))
                else
                   opart=opart+1
                   read(1,POS=mass_blck+sizeof(dummy_real)*(opart-1)) mm_sp(i)
                endif
                ! Optional data
                if(id_blck.ne.-1) then
                   read(1,POS=id_blck+sizeof(dummy_int)*(kpart-1)) ii(i)
                else
                   ii(i) = kpart
                endif
                if(kpart.le.header%npart(1)) then
                   if((u_blck.ne.-1).and.(u_size.eq.header%npart(1))) then
                      read(1,POS=u_blck+sizeof(dummy_real)*(kpart-1)) uu_sp(i)
                   endif
                endif
                if(metal) then
                   if((metal_blck.ne.-1).and.(metal_size.eq.npart)) then
                      read(1,POS=metal_blck+sizeof(dummy_real)*(kpart-1)) zz_sp(i)
                   endif
                   if((metal_blck.ne.-1).and.(metal_size.eq.ngas+nstar_tot)) then
                      read(1,POS=metal_blck+sizeof(dummy_real)*(gpart-1)) zz_sp(i)
                   endif
                endif
                if(star) then
                   if((age_blck.ne.-1).and.(age_size.eq.sum(header%npart(3:5)))) then
                      if((sum(header%npart(3:5)).gt.0).and.(kpart.gt.(header%npart(1)+header%npart(2)))) then
                         read(1,POS=age_blck+sizeof(dummy_real)*(mpart-1)) tt_sp(i)
                      endif
                   endif
                endif
                ! Scaling to ramses code units
                if(cosmo) then
                   gadget_scale_l = scale_l/header%boxsize
                   gadget_scale_v = 1e3*SQRT(aexp)/header%boxsize*aexp/100.
                endif
                xx(i,:)   = xx_sp(i,:)*(gadget_scale_l/scale_l)*ic_scale_pos
                vv(i,:)   = vv_sp(i,:)*(gadget_scale_v/scale_v)*ic_scale_vel
                mm(i)     = mm_sp(i)*(gadget_scale_m/scale_m)*ic_scale_mass
                if(cosmo) then
                   if(type_index .eq. 1) mass_sph = mm(i)
                   if(xx(i,1)<  0.0d0  )xx(i,1)=xx(i,1)+dble(nx)
                   if(xx(i,1)>=dble(nx))xx(i,1)=xx(i,1)-dble(nx)
                   if(xx(i,2)<  0.0d0  )xx(i,2)=xx(i,2)+dble(ny)
                   if(xx(i,2)>=dble(ny))xx(i,2)=xx(i,2)-dble(ny)
                   if(xx(i,3)<  0.0d0  )xx(i,3)=xx(i,3)+dble(nz)
                   if(xx(i,3)>=dble(nz))xx(i,3)=xx(i,3)-dble(nz)
                endif

                if(metal) then
                   if(metal_blck.ne.-1) then
                      zz(i) = zz_sp(i)*ic_scale_metal
                   else
                      zz(i) = 0.02*z_ave
                   endif
                endif
                if(kpart.gt.header%npart(1)+header%npart(2)) then
                   if(age_blck.ne.-1) then
                      if(cosmo) then
                         tt(i) = tt_sp(i)
                      else
                         tt(i) = tt_sp(i)*(gadget_scale_t/(scale_t/aexp**2))*ic_scale_age
                      endif
                   else
                      tt(i) = -13.8*1d9*3.15360d7/scale_t ! Age of the universe
                   endif
                endif
                if(kpart.le.header%npart(1)) then
                   if(cosmo) then
                      uu(i) = T2_start/scale_T2
                   else
                      ! Temperature stored in units of K/mu
                      uu(i) = uu_sp(i)*mu_mol*(gadget_scale_v/scale_v)**2*ic_scale_u
                   endif

                endif
                if(kpart.le.header%npart(1)) mgas_tot = mgas_tot+mm(i)
                ! Check the End Of Block
                if(kpart.ge.ipbar*(npart/49.0))then
                   write(*,'(A1)',advance='no') "_"
                   ipbar = ipbar+1.0
                endif
                if(kpart.ge.npart) then
                   write(*,'(A1)') " "
                   write(*,'(A,A7,A)') ' ',TRIM(ic_format),' file successfully loaded'
                   write(*,'(A50)')"__________________________________________________"
                   eob=.true.
                   exit
                endif
             enddo
          endif
#ifndef WITHOUTMPI
          call MPI_BCAST(eob,1         ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(xx,nvector*3  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(vv,nvector*3  ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(ii,nvector    ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(mm,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(zz,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(tt,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(uu,nvector    ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(jpart,1       ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call MPI_BCAST(header%npart,6,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
          call cmp_cpumap(xx,cc,jpart)
#endif
          do i=1,jpart
#ifndef WITHOUTMPI
             ! Check the CPU map
             if(cc(i)==myid)then
#endif
                ! Determine current particle type
                if((lpart+i).le.header%npart(1)) type_index = 1
                do j=1,5
                   if((lpart+i).gt.sum(header%npart(1:j)).and.(lpart+i).le.sum(header%npart(1:j+1))) type_index = j+1
                enddo
                skip           = .false.
                do j=1,6
                   if(ic_skip_type(j).eq.type_index-1) skip=.true.
                enddo
                if(.not.skip) then
                   if(abs(xx(i,1)-ic_center(1)).ge.boxlen/2d0) cycle
                   if(abs(xx(i,2)-ic_center(2)).ge.boxlen/2d0) cycle
                   if(abs(xx(i,3)-ic_center(3)).ge.boxlen/2d0) cycle
                   ipart          = ipart+1
                   if(ipart.gt.npartmax) then
                      write(*,*) "Increase npartmax"
#ifndef WITHOUTMPI
                      call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
                      stop
#endif
                   endif
                   xp(ipart,1:3)  = xx(i,1:3)+boxlen/2.0D0-ic_center(1:3)
                   vp(ipart,1:3)  = vv(i,1:3)
                   ! Flag gas particles with idp=1
                   if(type_index.gt.1)then
                      idp(ipart)   = ii(i)+1
                   else
                      idp(ipart)   = 1
                   endif
                   mp(ipart)      = mm(i)
                   if(use_initial_mass) mp0(ipart) = mm(i)
                   levelp(ipart)  = levelmin
                   if(star) then
                      tp(ipart)    = tt(i)
                      ! Particle metallicity
                      if(metal) then
                         zp(ipart)  = zz(i)
                      endif
                   endif
                   if(type_index.gt.2)then
                      if(star)then
                         typep(ipart)%family = FAM_STAR
                         typep(ipart)%tag    = 0
                         if(sn_ic) then
                            ! Only young enough stars should be active
                            if (tp(ipart).gt.tyoung) typep(ipart)%tag = TAG_STAR_ACTIVE
                         end if
                      end if
                   else if(type_index.eq.2)then
                      typep(ipart)%family = FAM_DM
                      typep(ipart)%tag    = 0
                   end if
                   up(ipart)      = uu(i)
                   if(ic_mask_ptype.gt.-1)then
                      if(ic_mask_ptype.eq.type_index-1)then
                         maskp(ipart) = 1.0
                      else
                         maskp(ipart) = 0.0
                      endif
                   endif
                   ! Add a gas particle outside the zoom region
                   if(cosmo) then
                      do j=1,6
                         if(type_index.eq.cosmo_add_gas_index(j)) then
                            ! Add a gas particle
                            xp(ipart+1,1:3) = xp(ipart,1:3)
                            vp(ipart+1,1:3) = vp(ipart,1:3)
                            idp(ipart+1)    = -1
                            mp(ipart+1)     = mp(ipart)*(omega_b/omega_m)
                            levelp(ipart+1) = levelmin
                            up(ipart+1)     = T2_start/scale_T2
                            if(metal) then
                               zp(ipart+1)  = z_ave*0.02
                            endif
                            ! Remove mass from the DM particle
                            mp(ipart) = mp(ipart)-mp(ipart+1)
                            ! Update index
                            ipart           = ipart+1
                         endif
                      end do
                   endif
                endif
#ifndef WITHOUTMPI
             endif
#endif
          enddo
          lpart = lpart+jpart
       enddo
       if(myid==1)then
          write(*,'(A,E10.3,A)') ' Gas mass in AMR grid -> ',mgas_tot,' unit_m'
          write(*,'(A50)')"__________________________________________________"
          close(1)
       endif
    enddo
    npart = ipart
    ! Compute total number of particle
    npart_cpu       = 0
    npart_all       = 0
    npart_cpu(myid) = npart
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
    npart_cpu(1) = npart_all(1)
#else
    npart_all       = npart
#endif
    if(myid==1)then
       write(*,*) ' npart_tot -> ',sum(npart_all)
       write(*,'(A50)')"__________________________________________________"
       close(1)
    endif
    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do
    if(debug)write(*,*)'npart=',npart,'/',npart_cpu(ncpu)
    ifout = ic_ifout
    t = ic_t_restart
    ! DICE patch
  end subroutine load_dice

#endif
  subroutine load_tracers
    implicit none

    real(dp):: xx1, xx2, xx3
    integer(1), dimension(1:nvector) :: ixx

    ! The tracers are loaded after all the other particles, so the first tracer
    ! is the particle number 'npart'
    ipart=npart

    if(myid==1)then
       open(10, file=trim(tracer_feed), form='formatted', status='old')
       write(*, *) 'Reading initial tracers from ', trim(tracer_feed)
    end if
    eof=.false.

    ! Binary mode
    do while (.not.eof)
       xx=0.0
       if(myid==1)then
          jpart=0
          do i=1,nvector
             read(10,*,end=100)xx1,xx2,xx3
             jpart=jpart+1
             indglob=indglob+1
             xx(i,1)=xx1*boxlen !+boxlen/2.0
             xx(i,2)=xx2*boxlen !+boxlen/2.0
             xx(i,3)=xx3*boxlen !+boxlen/2.0
             ii(i  )=indglob
             ixx(i )=FAM_TRACER_GAS
          end do
100       continue
          if(jpart<nvector)eof=.true.
       endif
       buf_count=nvector*3
#ifndef WITHOUTMPI
       call MPI_BCAST(xx,buf_count,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(ii,nvector  ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(ixx,nvector ,MPI_INTEGER1        ,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(eof,1       ,MPI_LOGICAL         ,0,MPI_COMM_WORLD,info)
       call MPI_BCAST(jpart,1     ,MPI_INTEGER         ,0,MPI_COMM_WORLD,info)
       call cmp_cpumap(xx,cc,jpart)
#endif

       do i=1,jpart
#ifndef WITHOUTMPI
          if(cc(i)==myid)then
#endif
             ipart=ipart+1
             if(ipart>npartmax)then
                write(*,*)'Maximum number of particles incorrect'
                write(*,*)'npartmax should be greater than',ipart, 'got', npartmax
                stop
             endif
             xp(ipart,:)  = xx(i,:)
             vp(ipart,:)  = vv(i,:)
             mp(ipart)    = tracer_mass
             levelp(ipart)= levelmin
             idp(ipart)   = ii(i)
             typep(ipart)%family  = ixx(i)
#ifndef WITHOUTMPI
          endif
#endif
       enddo

    end do
    if(myid==1)close(10)
    ! end if
    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=count(is_tracer(typep(:)) .and. (levelp(:) > 0))
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    npart_cpu(1)=npart_all(1)
#endif
    write(*,*)'npart=',npart_cpu(myid),'/',sum(npart_cpu), '(tracers)'

    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do

  end subroutine load_tracers

  !=================================================================
  ! Loads the tracer mass and decide whether we're using v1 or v2 of
  ! initial tracers
  !=================================================================
  subroutine load_tracers_bin(iversion)
    integer, intent(in) :: iversion
    integer :: unit_record, ntot, ipos
    real(dp) :: tmp_tracer_mass

    if (myid == 1) then
       open(newunit=unit_record, file=trim(tracer_feed), &
            form='unformatted', status='old')
       read(unit_record) ntot
       read(unit_record) tmp_tracer_mass

       if (tracer_mass > 0) then
          if (tmp_tracer_mass /= tracer_mass) then
             write(*, *) 'WARNING: the tracer mass from file differs from the one from namelist. Keeping latter.'
          end if
          write(*, *) 'Using a tracer mass of ', tracer_mass
       else
          if (tmp_tracer_mass > 0) then
             tracer_mass = tmp_tracer_mass
             write(*, *) 'Using a tracer mass of ', tracer_mass
          end if
       end if

       call ftell(unit_record, ipos)
       close(unit_record)
    end if

#ifndef WITHOUTMPI
    call MPI_BCAST(ntot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(tracer_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif

    if (iversion == 2) then
       call load_tracers_bin_v2(ntot)
    else
       if (myid == 1)  write(*, *) 'Reading initial tracers (binary v1) from ', trim(tracer_feed)

       call load_tracers_bin_v1(ntot)
    end if

  end subroutine load_tracers_bin

  !------------------------------------------------------------
  ! Create the tracer inplace by looping on the amr grid
  subroutine load_tracers_inplace
    integer :: nx_loc, icpu, jgrid, igrid, j, icell, iskip
    integer :: ix, iy, iz, npart_tot
    real(dp) :: scale, dx, dx_loc, vol_loc, d

    real(dp) :: xcell(ndim), skip_loc(ndim)
    integer(i8b) :: itracer_start
    integer(i8b) :: ntracer_loc, ntracer_cpu(ncpu), idp_start

    real(dp) :: dx_cell(twotondim, ndim)

    real(dp) :: npart_loc_real, rand
    integer :: npart_loc

    ! Broadcast the number of particles for the id of the tracers
#ifndef WITHOUTMPI
    call MPI_ALLREDUCE(npart, npart_tot, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, info)
#endif

    ! Build the positions of the cells w.r.t. their grid in dx unit
    do ind = 1, twotondim
       iz = (ind-1)/4
       iy = (ind-1-4*iz)/2
       ix = (ind-1-4*iz-2*iy)

       dx_cell(ind, 1) = (real(ix, dp)-0.5_dp)
       dx_cell(ind, 2) = (real(iy, dp)-0.5_dp)
       dx_cell(ind, 3) = (real(iz, dp)-0.5_dp)
    end do

    nx_loc = (icoarse_max - icoarse_min + 1)
    scale = boxlen / dble(nx_loc)

    skip_loc = [0.0d0, 0.0d0, 0.0d0]

    if (ndim > 0) skip_loc(1) = dble(icoarse_min)
    if (ndim > 1) skip_loc(2) = dble(jcoarse_min)
    if (ndim > 2) skip_loc(3) = dble(kcoarse_min)

    ! Store index of first tracer
    itracer_start = npart
    ipart = npart

    npart_loc_real = 0

    ! Loop over levels
    do ilevel = levelmin, nlevelmax
       dx = 0.5_dp**(ilevel)
       dx_loc = dx * scale
       vol_loc = dx_loc**3

       do jgrid = 1, active(ilevel)%ngrid
          igrid = active(ilevel)%igrid(jgrid)
          ! Loop on cells
          do ind = 1, twotondim
             iskip = ncoarse + (ind-1) * ngridmax
             icell = iskip + igrid

             ! Select leaf cells
             if (son(icell) == 0) then
                ! In zoomed region (if any)
                if (ivar_refine > 0) then
                   if (uold(icell, ivar_refine) / uold(icell, 1) < var_cut_refine) then
                      cycle
                   end if
                end if

                ! Compute number of tracers to create
                d = uold(icell, 1) * vol_loc
                npart_loc_real = d / tracer_mass
                npart_loc = int(npart_loc_real)

                ! The number of tracer is real, so we have to decide
                ! whether the number is the floor or ceiling of the
                ! real number.
                call ranf(tracer_seed, rand)

                if (rand < npart_loc_real-npart_loc) then
                   npart_loc = npart_loc + 1
                end if

                ! Get cell position
                xcell(:) = (xg(igrid, :) - skip_loc(:) + dx_cell(ind, :) * dx) * scale

                ! Now create the right number of tracers
                !
                ! Note: we don't create the idp of the tracers here. See below.
                do j = 1, npart_loc
                   ipart = ipart+1
                   if (ipart > npartmax) then
                      write(*,*) 'Maximum number of particles incorrect'
                      write(*,*) 'npartmax should be greater than', ipart, 'got', npartmax
                      stop
                   end if
                   xp(ipart, 1) = xcell(1)
                   xp(ipart, 2) = xcell(2)
                   xp(ipart, 3) = xcell(3)

                   vp(ipart, :) = 0._dp
                   mp(ipart) = tracer_mass
                   levelp(ipart) = ilevel
                   typep(ipart)%family = FAM_TRACER_GAS
                end do
             end if
          end do
          ! Get next grid
       end do
       ! End loop over active grids
    end do ! End loop over levels

    ! Store total number of particules
    npart = ipart

    ! Count tracers and scatter to other CPUs
    ntracer_loc = npart - itracer_start
    ntracer_cpu(myid) = ntracer_loc

#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLGATHER(ntracer_loc, 1, MPI_INTEGER, ntracer_cpu, 1, MPI_INTEGER, MPI_COMM_WORLD, info)
#else
    call MPI_ALLGATHER(ntracer_loc, 1, MPI_INTEGER8, ntracer_cpu, 1, MPI_INTEGER8, MPI_COMM_WORLD, info)
#endif
#endif

    ! Compute number of tracer in CPUs of lesser rank
    do icpu = 2, ncpu
       ntracer_cpu(icpu) = ntracer_cpu(icpu-1) + ntracer_cpu(icpu)
    end do

    ! Get first available index: this is the total number of
    ! particules + the number of tracers in CPUs with smaller ranks
    if (myid == 1) then
       idp_start = npart_tot
    else
       idp_start = npart_tot + ntracer_cpu(myid-1)
    end if

    ! Now loop on the created particles and give them an id
    do ipart = itracer_start+1, npart+1
       idp_start = idp_start + 1
       idp(ipart) = idp_start
    end do

    ! Update the global counter (useless if nothing is loaded after the tracers)
    if (myid == 1) then
       indglob = indglob + ntracer_cpu(ncpu)
    end if

    if (myid == 1 .and. ntracer_cpu(ncpu) == 0) then
       write(*,*) '______ NO TRACER CREATED! ______'
    end if
    if (ntracer_loc > 0) &
         write(*,'(a,i15,a,i15,a,i7)') 'ntracer=', ntracer_loc, '/', ntracer_cpu(ncpu), &
         '(tracers) for PE=', myid

  end subroutine load_tracers_inplace

  !------------------------------------------------------------
  ! Read the tracer in version 1 of the format
  !
  ! The version 1 is a record based format with the following structure
  ! * ntracer [integer]
  ! * mtracer [float64]
  ! * x[ntracer] [float64]
  ! * y[ntracer] [float64]
  ! * z[ntracer] [float64]
  !
  subroutine load_tracers_bin_v1(ntot)
    integer, intent(in) :: ntot
    integer :: unit_in

    real(dp), dimension(nvector) :: xx1, xx2, xx3
    integer(1), dimension(1:nvector) :: ixx

    integer, dimension(nvector) :: ii, cmap
    real(dp), dimension(1:nvector, 1:ndim) :: tmpxx
    integer :: icpu

    integer :: j, jj, nbuffer
    integer :: ix, iy, iz

    ! The tracers are loaded after all the other particles, so the first tracer
    ! is the particle number 'npart'
    ipart=npart

    ! Because we read the file in stream mode, we need to know where
    ! each element is using this map:
    !
    !  Length | Starting position | Comment
    ! --------+-------------------+-----------------------
    !       4 | 1                 | record length
    !       4 | 5                 | number of particles N
    !       4 | 9                 | record end
    ! --------+-------------------+-----------------------
    !       4 | 14                | record length
    !       4 | 17                | tracer mass
    !       4 | 21                | record end
    ! --------+-------------------+-----------------------
    !       4 | 25                | record length
    !      8N | ix=29             | x
    !       4 | ix+8N             | record end
    ! --------+-------------------+-----------------------
    !       4 | ix+8N+4           | record length
    !      8N | iy=ix+8N+8        | y
    !       4 | iy+8N             | record end
    ! --------+-------------------+-----------------------
    !       4 | iy+8N+4           | record length
    !      8N | iz=iy+8N+8        | z
    !       4 | iz+8N             | record end

    if (myid == 1) then
       ix = 1 + 28
       iy = ix + ntot * 8 + 8
       iz = iy + ntot * 8 + 8

       open(newunit=unit_in, file=trim(tracer_feed), access='stream', status='old')
    end if


    do jj = 1, ntot, nvector
       if (jj + nvector > ntot) then
          nbuffer = ntot - jj + 1
       else
          nbuffer = nvector
       end if

       ! Read nbuffer elements from file
       if(myid == 1) then
          do j = 1, nbuffer
             read(unit_in, pos=ix + 8*(j+jj-1)) xx1(j)
          end do

          do j = 1, nbuffer
             read(unit_in, pos=iy + 8*(j+jj-1)) xx2(j)
          end do

          do j = 1, nbuffer
             read(unit_in, pos=iz + 8*(j+jj-1)) xx3(j)
          end do

          do j = 1, nbuffer
             indglob = indglob + 1
             ii(j) = indglob
             ixx(j) = FAM_TRACER_GAS
          end do
       end if

#ifndef WITHOUTMPI
       call MPI_BCAST(xx1, nbuffer, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(xx2, nbuffer, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(xx3, nbuffer, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(ii,  nbuffer, MPI_INTEGER         , 0, MPI_COMM_WORLD, info)
       call MPI_BCAST(ixx, nbuffer, MPI_INTEGER         , 0, MPI_COMM_WORLD, info)
#endif
       tmpxx(1:nbuffer, 1) = xx1(1:nbuffer) * boxlen
       tmpxx(1:nbuffer, 2) = xx2(1:nbuffer) * boxlen
       tmpxx(1:nbuffer, 3) = xx3(1:nbuffer) * boxlen

       call cmp_cpumap(tmpxx, cmap, nbuffer)

       do j = 1, nbuffer
#ifndef WITHOUTMPI
          if (cmap(j)==myid) then
#endif
             ipart=ipart+1
             if(ipart>npartmax)then
                write(*,*)'Maximum number of particles incorrect'
                write(*,*)'npartmax should be greater than',ipart, 'got', npartmax
                stop
             end if
             xp(ipart, 1)  = xx1(j)
             xp(ipart, 2)  = xx2(j)
             xp(ipart, 3)  = xx3(j)

             vp(ipart,:)  = 0._dp
             mp(ipart)    = tracer_mass
             levelp(ipart)= levelmin
             idp(ipart)   = ii(j)
             typep(ipart)%family = int(ixx(j), 1)
#ifndef WITHOUTMPI
          endif
#endif
       end do ! End loop on buffer
    end do ! End loop on particle number

    if (myid == 1) close(unit_in)

    ! end if
    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
    call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
    npart_cpu(1)=npart_all(1)
#endif
    do icpu=2,ncpu
       npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
    end do
    write(*,*)'npart=',npart,'/',npart_cpu(ncpu), '(tracers)'

  end subroutine load_tracers_bin_v1

  !------------------------------------------------------------
  ! Read the tracer in version 2 of the format
  !
  ! This version starts by 2 records containing:
  ! * ntracer [integer]
  ! * mtracer [float64]
  !
  ! The data is then written directly in binary, encoded as float64 in
  ! C order. If you read 3 contiguous parts, you'll get x, y and z of
  ! a given particle.
  subroutine load_tracers_bin_v2(ntot)
    integer, intent(in) :: ntot

    integer :: unit_in
    real(dp), dimension(:, :), allocatable :: allpos
    integer, dimension(:), allocatable :: allcmap
    integer, dimension(nvector) :: cmap
    real(dp), dimension(1:nvector, 1:ndim) :: tmpxx
    real(dp):: xx1, xx2, xx3

    type(communicator), dimension(:), allocatable :: sender  ! To send data
    type(communicator) :: receiver                           ! To receive data

#ifndef WITHOUTMPI
    integer, dimension(MPI_STATUS_SIZE) :: status
    integer, dimension(MPI_STATUS_SIZE, 2:ncpu) :: statuses
    integer :: ierror
    integer, dimension(2:ncpu) :: reqsend1, reqsend2
#endif

    integer :: icpu
    integer :: iwrite

    integer(i8b) :: pos
    integer :: j, jj, nbuffer, ibuff
    integer, dimension(ncpu) :: cpu_count
    integer :: cpu_count_loc

    ! The tracers are loaded after all the other particles, so the first tracer
    ! is the particle number 'npart'
    ipart=npart

    if (myid == 1) then
       allocate(allpos(ntot, 3), allcmap(ntot))
    end if

#ifndef WITHOUTMPI
    call MPI_BCAST(ntot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, info)
    call MPI_BCAST(tracer_mass, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, info)
#endif

    cpu_count(:) = 0
    if (myid == 1) then
       write(*, 321, advance='no') trim(tracer_feed)
       unit_in = 11
       open(unit=unit_in, file=trim(tracer_feed), access='stream', status='old')

       iwrite = 0
       do jj = 1, ntot, nvector
          if (iwrite < int((100. * jj) / ntot)) then
             iwrite = int((100. * jj) / ntot)
             write(*, '(".")', advance='no')
          end if

          nbuffer = min(nvector, ntot-jj+1)

          ! Read data from file
          do j = 1, nbuffer
             pos = 33 + int(jj+j-2, i8b)*24

             read(unit_in, pos=pos) xx1, xx2, xx3

             tmpxx(j, 1) = xx1 * boxlen
             tmpxx(j, 2) = xx2 * boxlen
             tmpxx(j, 3) = xx3 * boxlen
          end do

          ! Compute cpu_map
          call cmp_cpumap(tmpxx, cmap, nbuffer)

          ! Fill the sender
          do j = 1, nbuffer
             cpu_count(cmap(j)) = cpu_count(cmap(j)) + 1
             allcmap(jj+j-1) = cmap(j)
             allpos(jj+j-1, :) = tmpxx(j, :)
          end do
       end do

       close(unit_in)

       ! Force termination of line
       write (*,*) ''
    end if

#ifndef WITHOUTMPI
    ! Send count of particles to each CPU
    call MPI_SCATTER(&
         cpu_count, 1, MPI_INTEGER,     &
         cpu_count_loc, 1, MPI_INTEGER, &
         0, MPI_COMM_WORLD, ierror)
#endif
    write(*,'(a,i10,a,i10,a,i6)') 'ntracer=',cpu_count_loc,' /',ntot, ' for PE=', myid

    ! Allocate sender/receiver
    if (myid == 1) then
       allocate(sender(1:ncpu))
       do icpu = 1, ncpu
          allocate(&
               sender(icpu)%f(cpu_count(icpu), 1:1), &
               sender(icpu)%up(cpu_count(icpu), 1:3))
       end do
    end if
    allocate(&
         receiver%f(cpu_count_loc, 1:1), &
         receiver%up(cpu_count_loc, 1:3))

    ! Fill the send buffer
    if (myid == 1) then
       cpu_count(:) = 0
       do j = 1, ntot
          ibuff = cpu_count(allcmap(j)) + 1

          cpu_count(allcmap(j)) = ibuff
          indglob = indglob + 1

          sender(allcmap(j))%up(ibuff, :) = allpos(j, :)
          sender(allcmap(j))%f(ibuff, :) = indglob
       end do
       deallocate(allpos, allcmap)

#ifndef WITHOUTMPI
       do icpu = 2, ncpu
          call MPI_ISEND(sender(icpu)%up, 3*cpu_count(icpu), MPI_DOUBLE_PRECISION, &
               icpu-1, 0, MPI_COMM_WORLD, reqsend1(icpu), ierror)
          call MPI_ISEND(sender(icpu)%f,    cpu_count(icpu), MPI_INTEGER, &
               icpu-1, 0, MPI_COMM_WORLD, reqsend2(icpu), ierror)
       end do
       receiver%up = sender(myid)%up
       receiver%f = sender(myid)%f
#endif
    else
#ifndef WITHOUTMPI
       call MPI_RECV(receiver%up, 3*cpu_count_loc, MPI_DOUBLE_PRECISION, &
            0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
       call MPI_RECV(receiver%f,    cpu_count_loc, MPI_INTEGER, &
            0, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierror)
#endif
    end if

    ! Save the particles
    do j = 1, cpu_count_loc
       ipart=ipart+1
       if(ipart>npartmax)then
          write(*,*)'Maximum number of particles incorrect'
          write(*,*)'npartmax should be greater than', ipart, 'got', npartmax, 'for PE=', myid
          stop
       end if
       xp(ipart, 1)  = receiver%up(j, 1)
       xp(ipart, 2)  = receiver%up(j, 2)
       xp(ipart, 3)  = receiver%up(j, 3)

       vp(ipart,:)  = 0._dp
       mp(ipart)    = tracer_mass
       levelp(ipart)= levelmin
       idp(ipart)   = receiver%f(j, 1)
       typep(ipart)%family = FAM_TRACER_GAS
    end do

    npart=ipart

    ! Compute total number of particle
    npart_cpu=0; npart_all=0
    npart_cpu(myid)=npart

#ifndef WITHOUTMPI
    ! Wait for transmission end
    if (myid == 1) then
       call MPI_WAITALL(ncpu-1, reqsend1, statuses, ierror)
       call MPI_WAITALL(ncpu-1, reqsend2, statuses, ierror)
       do icpu = 2, ncpu
          deallocate(sender(icpu)%f, sender(icpu)%up)
       end do
       deallocate(sender)
    end if
    deallocate(receiver%f, receiver%up)
#endif

321 format('Reading initial tracers (binary v2) from ', A)
  end subroutine load_tracers_bin_v2

end subroutine init_part
#define TIME_START(cs) call SYSTEM_CLOCK(COUNT=cs)
#define TIME_END(ce) call SYSTEM_CLOCK(COUNT=ce)
#define TIME_SPENT(cs,ce,cr) REAL((ce-cs)/cr)

subroutine load_gadget
  ! This routine only creates DM particles
  use amr_commons
  use pm_commons
  use gadgetreadfilemod
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
  integer,dimension(1:nvector)::cc
#endif

  logical::ok
  TYPE(gadgetheadertype)::gadgetheader
  integer::numfiles
  integer::ifile
  real,dimension(:,:),allocatable:: pos, vel
  real(dp)::massparticles
  integer(kind=8)::allparticles
  integer(i8b),dimension(:),allocatable:: ids
  integer::nparticles
  integer::i,icpu,ipart,start
  integer(i8b),dimension(1:ncpu)::npart_cpu,npart_all
  character(LEN=256)::filename
  real(dp),dimension(1:nvector,1:3)::xx_dp
  integer::clock_start,clock_end,clock_rate
  real(dp)::gadgetvfact

  ! Local particle count
  ipart=0
  call SYSTEM_CLOCK(COUNT_RATE=clock_rate)

  if(TRIM(initfile(levelmin)).NE.' ')then
     filename=TRIM(initfile(levelmin))
     ! read first header to get information
     call gadgetreadheader(filename, 0, gadgetheader, ok)
     if(.not.ok) call clean_stop
     numfiles = gadgetheader%numfiles
     gadgetvfact = sqrt(aexp) / gadgetheader%boxsize * aexp / 100.
#ifndef LONGINT
     allparticles=int(gadgetheader%nparttotal(2),kind=8)
#else
     allparticles=int(gadgetheader%nparttotal(2),kind=8) &
          & +int(gadgetheader%totalhighword(2),kind=8)*4294967296_i8b !2^32
#endif
     massparticles=1d0/dble(allparticles)
     do ifile=0,numfiles-1
        call gadgetreadheader(filename, ifile, gadgetheader, ok)
        nparticles = gadgetheader%npart(2)
        allocate(pos(3,nparticles))
        allocate(vel(3,nparticles))
        allocate(ids(nparticles))
        TIME_START(clock_start)
        call gadgetreadfile(filename,ifile,gadgetheader, pos, vel, ids)
        TIME_END(clock_end)
        if(debug) write(*,*) myid, ':Read ', nparticles, ' from gadget file ', ifile, ' in ', &
             TIME_SPENT(clock_start, clock_end, clock_rate)
        start = 1
        TIME_START(clock_start)
        do i=1,nparticles
           xx_dp(1,1) = pos(1,i)/gadgetheader%boxsize
           xx_dp(1,2) = pos(2,i)/gadgetheader%boxsize
           xx_dp(1,3) = pos(3,i)/gadgetheader%boxsize
#ifndef WITHOUTMPI
           call cmp_cpumap(xx_dp,cc,1)
           if(cc(1)==myid)then
#endif
              ipart=ipart+1
#ifndef WITHOUTMPI
              if (ipart .ge. size(mp)) then
                 write(*,*) "For ", myid, ipart, " exceeds ", size(mp)
                 call clean_stop
              end if
#endif
              xp(ipart,1:3)=xx_dp(1,1:3)
              vp(ipart,1)  =vel(1, i) * gadgetvfact
              vp(ipart,2)  =vel(2, i) * gadgetvfact
              vp(ipart,3)  =vel(3, i) * gadgetvfact
              mp(ipart)    = massparticles
              levelp(ipart)=levelmin
              idp(ipart)   =ids(i)

              ! Get the particle type
              typep(ipart)%family = FAM_DM
              typep(ipart)%tag    = 0
#ifndef WITHOUTMPI
           endif
#endif
        enddo
#ifndef WITHOUTMPI
        TIME_END(clock_end)
        if(debug) write(*,*) myid, ':Processed ', nparticles, ' in ',&
             &  TIME_SPENT(clock_start, clock_end, clock_rate), " ipart now ", ipart
#endif
        deallocate(pos,vel,ids)
     end do

  end if
  npart=ipart
  ! Compute total number of particleclock_rate
  npart_cpu=0; npart_all=0
  npart_cpu(myid)=npart
#ifndef WITHOUTMPI
#ifndef LONGINT
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
#else
  call MPI_ALLREDUCE(npart_cpu,npart_all,ncpu,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,info)
#endif
  npart_cpu(1)=npart_all(1)
#endif
  do icpu=2,ncpu
     npart_cpu(icpu)=npart_cpu(icpu-1)+npart_all(icpu)
  end do
  write(*,*)'npart=',npart,'/',npartmax

end subroutine load_gadget
