subroutine read_params
  use amr_commons
  use pm_parameters
  use poisson_parameters
  use hydro_parameters
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef DICE
  use dice_commons
#endif
  use mpi_mod
  implicit none
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,narg,levelmax,levelhold=-1
  character(LEN=80)::infile, info_file
  character(LEN=80)::cmdarg
  character(LEN=5)::nchar
  integer(kind=8)::ngridtot=0
  integer(kind=8)::nparttot=0
  real(kind=8)::delta_tout=0,tend=0
  real(kind=8)::delta_aout=0,aend=0
  logical::nml_ok, info_ok, restart_file_ok
  integer,parameter::tag=1134
  integer::mythr
#ifndef WITHOUTMPI
  integer::dummy_io,ierr,info2
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/run_params/clumpfind,cosmo,pic,sink,sinkprops,lightcone,poisson,hydro,rt,verbose,debug &
       & ,nrestart,nrestart_seek,ncontrol,nstepmax,nsubcycle,load_weights,part_univ_cost,exact_timer,nremap,ordering &
       & ,bisec_tol,static,overload,cost_weighting,aton,nrestart_quad,restart_remap &
       & ,static_dm,static_gas,static_stars,convert_birth_times,use_proper_time,remap_pscalar &
       & ,dtstop,magic_number,nchunk,dtmax,sinkprops_dir,remove_invalid_particle
  namelist/output_params/output,noutput,foutput,aout,tout &
       & ,tend,delta_tout,aend,delta_aout,gadget_output,walltime_hrs,minutes_dump &
       & ,early_stop_hrs,dump_stop,foutput_timer,wallstep,output_now
  namelist/amr_params/levelmin,levelmax,ngridmax,ngridtot &
       & ,npartmax,nparttot,nexpand,boxlen,nlevel_collapse &
       & ,nsinkmax,levelhold,holdback
  namelist/poisson_params/epsilon,gravity_type,gravity_params &
       & ,cg_levelmin,cic_levelmax,npartmax_rho
  namelist/lightcone_params/thetay_cone,thetaz_cone,zmax_cone
  namelist/movie_params/levelmax_frame,nw_frame,nh_frame,ivar_frame &
       & ,xcentre_frame,ycentre_frame,zcentre_frame,movie_vars &
       & ,deltax_frame,deltay_frame,deltaz_frame,movie,zoom_only_frame &
       & ,imovout,imov,tstartmov,astartmov,tendmov,aendmov,proj_axis,movie_vars_txt &
       & ,theta_camera,phi_camera,dtheta_camera,dphi_camera,focal_camera,dist_camera,ddist_camera &
       & ,perspective_camera,smooth_frame,shader_frame,tstart_theta_camera,tstart_phi_camera &
       & ,tend_theta_camera,tend_phi_camera,method_frame,varmin_frame,varmax_frame
#ifdef DICE
  namelist/dice_params/ ic_file,ic_nfile,ic_format,IG_rho,IG_T2,IG_metal &
       & ,ic_head_name,ic_pos_name,ic_vel_name,ic_id_name,ic_mass_name &
       & ,ic_u_name,ic_metal_name,ic_age_name &
       & ,gadget_scale_l, gadget_scale_v, gadget_scale_m ,gadget_scale_t &
       & ,ic_scale_pos,ic_scale_vel,ic_scale_mass,ic_scale_u,ic_scale_age &
       & ,ic_scale_metal,ic_center,ic_ifout,amr_struct,ic_t_restart,ic_mag_const &
       & ,ic_mag_center_x,ic_mag_center_y,ic_mag_center_z &
       & ,ic_mag_axis_x,ic_mag_axis_y,ic_mag_axis_z &
       & ,ic_mag_scale_R,ic_mag_scale_H,ic_mag_scale_B,cosmo_add_gas_index,ic_skip_type &
       & ,ic_mask_ivar,ic_mask_min,ic_mask_max,ic_mask_ptype,analytic_gas_profile
#endif
  namelist/tracer_params/ MC_tracer,tracer,tracer_feed,tracer_feed_fmt,tracer_mass, &
       tracer_first_balance_part_per_cell,tracer_first_balance_levelmin,tracer_per_cell, &
       tracer_level,no_init_gas_tracer,tracer_to_jet
  ! MPI initialization
#ifndef WITHOUTMPI
#ifdef _OPENMP
  call MPI_INIT_THREAD(MPI_THREAD_SERIALIZED,info2,ierr)
  if(info2<MPI_THREAD_SERIALIZED) then
      write(*,*) 'Error: MPI_THREAD_SERIALIZED is not supported in current MPI library'
      stop
  end if
#else
  call MPI_INIT(ierr)
#endif
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)
  myid=myid+1 ! Careful with this...
#endif
#ifdef WITHOUTMPI
  ncpu=1
  myid=1
#endif
#ifdef _OPENMP
!$omp parallel private(mythr)
  mythr=omp_get_thread_num()+1
  if(mythr==1) then
     nthr=omp_get_num_threads()
  end if
!$omp end parallel
#endif
  !--------------------------------------------------
  ! Advertise RAMSES
  !--------------------------------------------------
  if(myid==1)then
  write(*,*)'_/_/_/       _/_/     _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'_/    _/    _/  _/    _/_/_/_/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/ _/ _/   _/        _/         _/       '
  write(*,*)'_/_/_/     _/_/_/_/   _/    _/     _/_/    _/_/_/       _/_/   '
  write(*,*)'_/    _/   _/    _/   _/    _/         _/  _/               _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/   _/    _/  _/         _/    _/ '
  write(*,*)'_/    _/   _/    _/   _/    _/    _/_/_/   _/_/_/_/    _/_/_/  '
  write(*,*)'                        Version 3.0                            '
  write(*,*)'       written by Romain Teyssier (University of Zurich)       '
  write(*,*)'               (c) CEA 1999-2007, UZH 2008-2014                '
  write(*,*)' '
#ifdef _OPENMP
  write(*,'(" Working with nproc = ",I4," and nthr = ",I3," for ndim = ",I1)')ncpu,nthr,ndim
  write(*,'(" With nvector = ",I3)')nvector
#else
  write(*,'(" Working with nproc = ",I4," for ndim = ",I1)')ncpu,ndim
  write(*,'(" With nvector = ",I3)')nvector
#endif
#ifdef SOLVERhydro
  write(*,'(" Using solver = hydro with nvar = ",I2)')nvar
  ! Check nvar is not too small
  if(nvar<ndim+2)then
     write(*,*)'You should have: nvar>=ndim+2'
     write(*,'(" Please recompile with -DNVAR=",I2)')ndim+2
     call clean_stop
  endif
#endif

#ifdef SOLVERmhd
  write(*,'(" Using solver = mhd with nvar = ",I2)')nvar
  if(nvar<8)then
     write(*,*)'You should have: nvar>=8'
     write(*,'(" Please recompile with -DNVAR=8")')
     call clean_stop
  endif
#endif

  !Write I/O group size information
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '
  if(IOGROUPSIZE>0) write(*,*)'IOGROUPSIZE=',IOGROUPSIZE
  if(IOGROUPSIZECONE>0) write(*,*)'IOGROUPSIZECONE=',IOGROUPSIZECONE
  if(IOGROUPSIZEREP>0) write(*,*)'IOGROUPSIZEREP=',IOGROUPSIZEREP
  if(IOGROUPSIZE>0.or.IOGROUPSIZECONE>0.or.IOGROUPSIZEREP>0)write(*,*)' '

  ! Write information about git version
  call write_gitinfo

  ! Read namelist filename from command line argument
  narg = command_argument_count()
  IF(narg .LT. 1)THEN
     write(*,*)'You should type: ramses3d input.nml [nrestart]'
     write(*,*)'File input.nml should contain a parameter namelist'
     write(*,*)'nrestart is optional'
     call clean_stop
  END IF
  CALL GET_COMMAND_ARGUMENT(1,infile)
  endif
#ifndef WITHOUTMPI
  call MPI_BCAST(infile,80,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Read the namelist
  !-------------------------------------------------

  ! Wait for the token
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif


  namelist_file=TRIM(infile)
  INQUIRE(file=infile,exist=nml_ok)
  if(.not. nml_ok)then
     if(myid==1)then
        write(*,*)'File '//TRIM(infile)//' does not exist'
     endif
     call clean_stop
  end if

  !-------------------------------------------------
  ! Default passive scalar map
  !-------------------------------------------------
#if NVAR>NDIM+2
  allocate(remap_pscalar(1:nvar-(ndim+2)))
  do i=1,nvar-(ndim+2)
     remap_pscalar(i) = i+ndim+2
  enddo
#endif

  open(1,file=infile)
  rewind(1)
  read(1,NML=run_params)
  rewind(1)
  read(1,NML=output_params)
  rewind(1)
  read(1,NML=amr_params)
  rewind(1)
  read(1,NML=tracer_params,END=84)
84 continue
  if (tracer_first_balance_levelmin <= 0) tracer_first_balance_levelmin = levelmax + 1
  rewind(1)
  read(1,NML=lightcone_params,END=83)
83 continue
  rewind(1)
  read(1,NML=movie_params,END=82)
82 continue
  rewind(1)
  read(1,NML=poisson_params,END=81)
81 continue
#ifdef DICE
  rewind(1)
  read(1,NML=dice_params,END=106)
106 continue
#endif

  !-------------------------------------------------
  ! Read optional nrestart command-line argument
  !-------------------------------------------------
  if (myid==1 .and. narg == 2) then
     CALL GET_COMMAND_ARGUMENT(2,cmdarg)
     read(cmdarg,*) nrestart
  endif
  
  ! check for the most recent nout and restart from it. 
  if (myid==1 .and. nrestart == -1) then
     do while(nrestart == -1)
        call title(nrestart_seek,nchar)
        info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
        INQUIRE(FILE=info_file, EXIST=restart_file_ok)
        if (restart_file_ok) then
           nrestart = nrestart_seek
        else
           nrestart_seek = nrestart_seek - 1
        endif
        if(nrestart_seek <= 0) then
           nrestart = 0
        end if
     enddo
  endif
    
  if (myid==1 .and. nrestart .gt. 0) then
     call title(nrestart,nchar)
     info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
     inquire(file=info_file, exist=info_ok)
     do while(.not. info_ok .and. nrestart .gt. 1)
        nrestart = nrestart - 1
        call title(nrestart,nchar)
        info_file='output_'//TRIM(nchar)//'/info_'//TRIM(nchar)//'.txt'
        inquire(file=info_file, exist=info_ok)
     enddo
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(info_ok,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
#endif

  if (nrestart .gt. 0 .and. .not. info_ok) then
     if (myid==1) then
         write(*,*) "Error: Could not find restart file"
     endif
     call clean_stop
  endif

#ifndef WITHOUTMPI
  call MPI_BCAST(nrestart,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

  !-------------------------------------------------
  ! Compute time step for outputs
  !-------------------------------------------------
  if(tend>0)then
     if(delta_tout==0)delta_tout=tend
     noutput=MIN(int(tend/delta_tout),MAXOUT)
     do i=1,noutput
        tout(i)=dble(i)*delta_tout
     end do
  else if(aend>0)then
     if(delta_aout==0)delta_aout=aend
     noutput=MIN(int(aend/delta_aout),MAXOUT)
     do i=1,noutput
        aout(i)=dble(i)*delta_aout
     end do
  endif
  noutput=MIN(noutput,MAXOUT)
  if(imovout>0) then
     allocate(tmovout(0:imovout))
     allocate(amovout(0:imovout))
     tmovout=1d100
     amovout=1d100
     if(tendmov>0)then
        do i=0,imovout
           tmovout(i)=(tendmov-tstartmov)*dble(i)/dble(imovout)+tstartmov
        enddo
     endif
     if(aendmov>0)then
        do i=0,imovout
           amovout(i)=(aendmov-astartmov)*dble(i)/dble(imovout)+astartmov
        enddo
     endif
     if(tendmov==0.and.aendmov==0)movie=.false.
  endif
  !--------------------------------------------------
  ! Check for errors in the namelist so far
  !--------------------------------------------------
  levelmin=MAX(levelmin,1)
  nlevelmax=levelmax
  nlevelmax_current=levelmin
  if(levelhold>0)then
     nlevelsheld=nlevelmax-levelhold
     if(myid==1)write(*,*)'Using an effictive maximum level of',levelhold
  endif


  nml_ok=.true.
  if(levelmin<1)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmin should not be lower than 1 !!!'
     nml_ok=.false.
  end if
  if(nlevelmax<levelmin)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelmax should not be lower than levelmin'
     nml_ok=.false.
  end if
  if(levelhold> nlevelmax)then
     if(myid==1)write(*,*)'Error in the namelist:'
     if(myid==1)write(*,*)'levelhold has to be lower than levelmax'
     nml_ok=.false.
  endif

  if(ngridmax==0)then
     if(ngridtot==0)then
        if(myid==1)write(*,*)'Error in the namelist:'
        if(myid==1)write(*,*)'Allocate some space for refinements !!!'
        nml_ok=.false.
     else
        ngridmax=int(ngridtot/int(ncpu,kind=8),kind=4)
     endif
  end if
  if(npartmax==0)then
     npartmax=int(nparttot/int(ncpu,kind=8),kind=4)
  endif
  if(myid>1)verbose=.false.
  if(sink.and.(.not.pic))then
     pic=.true.
  endif
  !if(clumpfind.and.(.not.pic))then
  !   pic=.true.
  !endif
  !if(pic.and.(.not.poisson))then
  !   poisson=.true.
  !endif
  ! Check that we allow for sinks if needed
  if (sink .and. nsinkmax <= 0) then
     if (myid==1) write(*,*) 'Error in the namelist:'
     if (myid==1) write(*,*) 'With sinks, nsinkmax should be > 0'
     nml_ok = .false.
  end if

  call read_hydro_params(nml_ok)
#ifdef RT
  call rt_read_hydro_params(nml_ok)
#endif
#if NDIM==3
  if (clumpfind)call read_clumpfind_params
#endif
  if (movie)call set_movie_vars

  close(1)

  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif

  !-----------------
  ! Max size checks
  !-----------------
  if(nlevelmax>MAXLEVEL)then
     write(*,*) 'Error: nlevelmax>MAXLEVEL'
     call clean_stop
  end if
  if(nregion>MAXREGION)then
     write(*,*) 'Error: nregion>MAXREGION'
     call clean_stop
  end if
#ifndef QUADHILBERT
  if(nlevelmax>21) then
     write(*,*)'Error: nlevelmax>21 requires compiling with -DQUADHILBERT'
     call clean_stop
  end if
#endif

  !-----------------
  ! MC tracer
  !-----------------
  if(MC_tracer .and. (.not. tracer))then
     write(*,*)'You have activate the MC tracer but not the tracers.'
     call clean_stop
  end if

  if(MC_tracer .and. (.not. pic)) then
     write(*,*)'You have activate the MC tracer PIC is false.'
     call clean_stop
  end if

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     nexpand     (i)=nexpand     (i-levelmin+1)
     nsubcycle   (i)=nsubcycle   (i-levelmin+1)
     load_weights(i)=load_weights(i-levelmin+1)
     r_refine    (i)=r_refine    (i-levelmin+1)
     a_refine    (i)=a_refine    (i-levelmin+1)
     b_refine    (i)=b_refine    (i-levelmin+1)
     x_refine    (i)=x_refine    (i-levelmin+1)
     y_refine    (i)=y_refine    (i-levelmin+1)
     z_refine    (i)=z_refine    (i-levelmin+1)
     m_refine    (i)=m_refine    (i-levelmin+1)
     exp_refine  (i)=exp_refine  (i-levelmin+1)
     initfile    (i)=initfile    (i-levelmin+1)
  end do
  do i=1,levelmin-1
     nexpand     (i)= 1
     nsubcycle   (i)= 1
     load_weights(i)= 1
     r_refine    (i)=-1.0
     a_refine    (i)= 1.0
     b_refine    (i)= 1.0
     x_refine    (i)= 0.0
     y_refine    (i)= 0.0
     z_refine    (i)= 0.0
     m_refine    (i)=-1.0
     exp_refine  (i)= 2.0
     initfile    (i)= ' '
  end do

  if(.not.cosmo)then
     use_proper_time=.false.
     convert_birth_times=.false.
  endif

  if(.not. nml_ok)then
     if(myid==1)write(*,*)'Too many errors in the namelist'
     if(myid==1)write(*,*)'Aborting...'
     call clean_stop
  end if

#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#endif

end subroutine read_params
