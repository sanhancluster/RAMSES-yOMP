module amr_parameters

  ! Define real types
  integer,parameter::sp=kind(1.0E0)
#ifndef NPRE
  integer,parameter::dp=kind(1.0E0) ! default
#else
#if NPRE==4
  integer,parameter::dp=kind(1.0E0) ! real*4
#else
  integer,parameter::dp=kind(1.0D0) ! real*8
#endif
#endif
#ifdef QUADHILBERT
  integer,parameter::qdp=kind(1.0_16) ! real*16
#else
  integer,parameter::qdp=kind(1.0_8) ! real*8
#endif
  integer,parameter::MAXOUT=1000
  integer,parameter::MAXLEVEL=100

  ! Define integer types (for particle IDs mostly)
  ! Warning: compiler needs to accept fortran standard 2003.
  ! Specific numbers for fortran kinds are, in principle, implementation
  ! dependent, so "i8b=8" with "integer(i8b)" is implementation-dependent.
  ! See portability note for non-gcc compilers: (boud 2016-11-29) -
  ! https://gcc.gnu.org/onlinedocs/gfortran/ISO_005fFORTRAN_005fENV.html
  ! The selected_int_kind approach below is more likely to be portable:
  integer,parameter::i4b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
#ifndef LONGINT
  !integer,parameter::i8b=4  ! default long int are 4-byte int
  integer,parameter::i8b=selected_int_kind(9) ! since log(2*10^9)/log(2)=30.9
#else
  integer, parameter :: i8b=selected_int_kind(18) ! since log(2*10^18)/log(2)=60.8
  !integer,parameter::i8b=8  ! long int are 8-byte int
#endif /* LONGINT */

  ! Number of dimensions
#ifndef NDIM
  integer,parameter::ndim=1
#else
  integer,parameter::ndim=NDIM
#endif
  integer,parameter::twotondim=2**ndim
  integer,parameter::threetondim=3**ndim
  integer,parameter::twondim=2*ndim

  ! Vectorization parameter
#ifndef NVECTOR
  integer,parameter::nvector=500  ! Size of vector sweeps
#else
  integer,parameter::nvector=NVECTOR
#endif
  integer::nchunk=2  ! Size of OpenMP chunk size

  ! Precompute powers of 2 in the right prec for hilbert curve (MT)
  ! Careful, starts at zero, so that powof2(0) == 1
  real(qdp),parameter,dimension(0:ndim*32-1) :: powof2 = (/ (2.0_qdp**ipower, ipower=0,ndim*32-1) /)
  ! real(qdp),parameter:: zero_qdp = 0._qdp, one_qdp = 1._qdp
  !/MT

  integer, parameter :: nstride = 65536

  ! Run control
  logical::verbose =.false.   ! Write everything
  logical::hydro   =.false.   ! Hydro activated
  logical::pic     =.false.   ! Particle In Cell activated
  logical::poisson =.false.   ! Poisson solver activated
  logical::cosmo   =.false.   ! Cosmology activated
  logical::star    =.false.   ! Star formation activated
  logical::sink    =.false.   ! Sink particles activated
  logical::sinkprops=.false.  ! Write sink properties at each coarse step
  logical::rt      =.false.   ! Radiative transfer activated
  logical::debug   =.false.   ! Debug mode activated
  logical::static  =.false.   ! Static mode activated
  logical::static_dm=.false.  ! Static mode for dm only activated
  logical::static_gas=.false. ! Static mode for gas only activated
  logical::static_stars=.false.! Static mode for stars only activated
  logical::tracer  =.false.   ! Tracer particles activated
  ! PATCH MC Tracer !
  logical::MC_tracer = .false. ! Tracer using Monte Carlo
  ! =============== !
  logical::lightcone=.false.  ! Enable lightcone generation
  logical::clumpfind=.false.  ! Enable clump finder
  logical::aton=.false.       ! Enable ATON coarse grid radiation transfer
  logical::exact_timer=.false.! Measure exact timer for each routine using barriers
  real(dp)::part_univ_cost=0d0  ! Give this amount of extra cost weights to particle (regardless of level)

  ! Mesh parameters
  integer::nx=1,ny=1,nz=1     ! Number of coarse cells in each dimension
  integer::levelmin=1         ! Full refinement up to levelmin
  integer::nlevelmax=1        ! Maximum number of level
  integer::nlevelsheld=0      ! Intermediate level at which to cut refinement
  integer::ngridmax=0         ! Maximum number of grids
  integer,dimension(1:MAXLEVEL)::nexpand=1 ! Number of mesh expansion
  integer::nexpand_bound=1    ! Number of mesh expansion for virtual boundaries
  real(dp)::boxlen=1.0D0      ! Box length along x direction
  character(len=128)::ordering='hilbert'
  logical::cost_weighting=.true. ! Activate load balancing according to cpu time
  ! Recursive bisection tree parameters
  integer::nbilevelmax=1      ! Max steps of bisection partitioning
  integer::nbinodes=3         ! Max number of internal nodes
  integer::nbileafnodes=2     ! Max number of leaf (terminal) nodes
  real(dp)::bisec_tol=0.05d0  ! Tolerance for bisection load balancing

  ! Step parameters
  integer::nrestart=0         ! New run or backup file number
  integer::nrestart_seek=5000 ! Search for the final snapshot in descending order 
  integer::nrestart_quad=0    ! Restart with double precision Hilbert keys
  real(dp)::trestart=0.0      ! Restart time
  logical::restart_remap=.false. ! Force load balance on restart
  integer::nstepmax=1000000   ! Maximum number of time steps
  integer::ncontrol=1         ! Write control variables
  integer::nremap=0           ! Load balancing frequency (0: never)
  integer,allocatable,dimension(:)::remap_pscalar
  real(dp)::dtstop=0          ! Halt the simulation when dt is smaller than this value

  ! Output parameters
  logical::output=.true.      ! Write output files
  integer::iout=1             ! Increment for output times
  integer::ifout=1            ! Increment for output files
  integer::iback=1            ! Increment for backup files
  integer::noutput=1          ! Total number of outputs
  integer::foutput=1000000    ! Frequency of outputs
  logical::gadget_output=.false. ! Output in gadget format
  logical::output_now=.false. ! write output next step
  real(dp)::walltime_hrs=-1.  ! Wallclock time for submitted job
  real(dp)::minutes_dump=1.   ! Dump an output minutes before walltime ends
  logical::dump_stop=.true.   ! Stop when non-scheduled dump is finished
  integer::foutput_timer=-1   ! Output frequency of timer
  integer::wallstep=-1       ! Number of step foward to finish

  ! Lightcone parameters
  real(dp)::thetay_cone=12.5
  real(dp)::thetaz_cone=12.5
  real(dp)::zmax_cone=2.0

  ! Cosmology and physical parameters
  real(dp)::boxlen_ini        ! Box size in h-1 Mpc
  real(dp)::omega_b=0.045     ! Omega Baryon
  real(dp)::omega_m=1.0D0     ! Omega Matter
  real(dp)::omega_l=0.0D0     ! Omega Lambda
  real(dp)::omega_k=0.0D0     ! Omega Curvature
  real(dp)::h0     =1.0D0     ! Hubble constant in km/s/Mpc
  real(dp)::aexp   =1.0D0     ! Current expansion factor
  real(dp)::hexp   =0.0D0     ! Current Hubble parameter
  real(dp)::texp   =0.0D0     ! Current proper time
  real(dp)::n_sink = -1.d0    ! Sink particle density threshold in H/cc
  real(dp)::ns_sink = -1.d0   ! Sink particle stellar density threshold in H/cc
  real(dp)::rho_sink = -1.D0  ! Sink particle density threshold in g/cc
  real(dp)::d_sink = -1.D0    ! Sink particle density threshold in user units
  real(dp)::m_star =-1.0      ! Star particle mass in units of mass_sph
  real(dp)::n_star =0.1D0     ! Star formation density threshold in H/cc
  real(dp)::t_star =0.0D0     ! Star formation time scale in Gyr
  real(dp)::eps_star=0.0D0    ! Star formation efficiency (0.02 at n_star=0.1 gives t_star=8 Gyr)
  real(dp)::T2_star=0.0D0     ! Typical ISM polytropic temperature
  real(dp)::g_star =1.6D0     ! Typical ISM polytropic index
  real(dp)::jeans_ncells=-1   ! Jeans polytropic EOS
  real(dp)::del_star=2.D2     ! Minimum overdensity to define ISM
  real(dp)::eta_sn =0.0D0     ! Supernova mass fraction
  real(dp)::eta_ssn=0.95      ! Single supernova ejected mass fraction (sf_imf=.true. only)
  real(dp)::yield  =0.0D0     ! Supernova yield
  real(dp)::f_ek   =1.0D0     ! Supernovae kinetic energy fraction (only between 0 and 1)
  real(dp)::rbubble=0.0D0     ! Supernovae superbubble radius in pc
  real(dp)::f_w    =0.0D0     ! Supernovae mass loading factor
  integer ::ndebris=1         ! Supernovae debris particle number
  real(dp)::mass_gmc=-1.0     ! Stochastic exploding GMC mass
  real(dp)::z_ave  =0.0D0     ! Average metal abundance
  real(dp)::B_ave  =0.0D0     ! Average magnetic field
  real(dp)::z_reion=8.5D0     ! Reionization redshift
  real(dp)::T2_start          ! Starting gas temperature
  real(dp)::T2max=huge(1._dp) ! Temperature ceiling for cooling_fine
  real(dp)::t_diss =20.0D0    ! Dissipation timescale for feedback
  real(dp)::t_sne =10.0D0     ! Supernova blast time
  real(dp)::J21    =0.0D0     ! UV flux at threshold in 10^21 units
  real(dp)::a_spec =1.0D0     ! Slope of the UV spectrum
  real(dp)::beta_fix=0.0D0    ! Pressure fix parameter
  real(dp)::kappa_IR=0d0      ! IR dust opacity
  real(dp)::ind_rsink=4.0d0   ! Number of cells defining the radius of the sphere where AGN feedback is active
  real(dp)::ir_eff=0.75       ! efficiency of the IR feedback (only when ir_feedback=.true.)
  real(dp)::sf_trelax=0.0D0   ! Relaxation time for star formation (cosmo=.false. only)
  real(dp)::sf_tdiss=0.0D0    ! Dissipation timescale for subgrid turbulence in units of turbulent crossing time
  integer::sf_model=3         ! Virial star formation model
  integer::nlevel_collapse=3  ! Number of levels to follow initial dark matter collapse (cosmo=.true. only)
  real(dp)::mass_star_max=120.0D0 ! Maximum mass of a star in solar mass
  real(dp)::mass_sne_min=10.0D0   ! Minimum mass of a single supernova in solar mass
  real(dp)::sn_trelax=0.D0    ! Blocks SN until t>sn_trelax
  real(dp)::f_ekAGN=1.0D0     ! AGN kinetic energy fraction (only between 0 and 1)
  real(dp)::rAGN   =0.0D0     ! AGN superbubble radius in kpc
  real(dp)::eAGN_K =1d0       ! AGN energy efficiency in Radio mode
  real(dp)::eAGN_T =0.15d0    ! AGN energy efficiency in Quasar mode
  real(dp)::TAGN   =0.0d0     ! AGN temperature floor for energy release in Quasar mode
  real(dp)::T2maxAGN=1d10     ! Maximum temperature allowed after AGN energy input
  real(dp)::jetfrac=0.0d0     ! Fraction of accreted mass before releasing jet AGN energy
  real(dp)::Mseed  =1.0d5     ! Mass of the initial sink particle in solar mass
  real(dp)::boost_acc=2.0d0   ! Boost power factor for the accretion rate
  real(dp)::boost_drag=2.0d0  ! Boost power factor for drag force
  real(dp)::r_gal  =1.0d2     ! SMBH sphere radius of influence in kpc
  real(dp)::sigmav_max=10d15  ! Maximum relative velocity in the Bondi accretion rate in kpc
  real(dp)::mloadAGN=1d2      ! Mass loading factor of the jet
  real(dp)::f_bondi=1d0       ! Fraction of the Bondi rate
  real(dp)::maxspin=0.998_dp  ! Maximum spin amplitude
  real(dp)::X_floor=1.0d-2    ! radio/quasar floor
  real(dp)::rmerge=1.0d0      ! Number of dx_min to allow for BH coalescence
  real(dp)::star_ratio_floor=0.25d0
  real(dp)::del_jeans=0.d0    ! Gas density contrast to trigger Jeans-based refinement criterion
  real(dp)::dens_jeans=0.d0   ! Gas density threshold trigger Jeans-based refinement criterion
  real(dp) :: d_boost=1       ! Density for boost Booth and Schaye (HP)
  integer  :: df_ncells=4     ! number of cells for DF
  logical::point_mass_sink=.false.          ! deposit the mass on central sink for gravity 

  real(dp)::zdmax=-1d0        ! Maximum allowed dust ratio
  
  logical ::self_shielding=.false.
  logical ::pressure_fix=.false.
  logical ::nordlund_fix=.true.
  logical ::cooling=.false.
  logical ::neq_chem=.false.  ! Non-equilbrium chemistry activated
  logical ::isothermal=.false.
  logical ::metal=.false.
  logical ::dust=.false.
  logical ::dust_cooling=.false.   ! Activate high-T cooling of the gas from dust
  logical ::metal_gasonly=.false.  ! Only count gas-phase metals for cooling
  character(LEN=20)::sticking_coef='chaabouni'
  logical ::haardt_madau=.false.
  logical ::delayed_cooling=.false.
  logical ::smbh=.false.
  logical ::agn=.false.
  logical ::use_proper_time=.false.
  logical ::convert_birth_times=.false. ! Convert stellar birthtimes: conformal -> proper
  logical ::ir_feedback=.false. ! Activate ir feedback from accreting sinks
  logical ::sf_virial=.false.   ! Activate SF Virial criterion
  logical ::sf_log_properties=.false. ! Log in ascii files birth properties of stars and supernovae
  logical ::sf_imf=.false.      ! Activate IMF sampling for SN feedback when resolution allows it
  logical ::sf_compressive=.false. ! Advect compressive and solenoidal turbulence terms separately
  logical ::sn_IC=.false.     ! Allow SN from initial stars (MT: only with dice)
  logical ::bondi=.true.      ! Activate Bondi accretion onto sink particle
  logical ::eddington_limit=.true.           ! Switch for Eddington limit for the smbh case
  logical ::sink_AGN=.false.  ! Activate AGN feedback on sink particles
  logical ::finestep_AGN=.false.  !Activate AGN feedback on the fine timestep. needs sink_AGN=.true.
  logical ::drag=.false.      ! Activate drag force
  logical ::random_jet=.false.
  logical ::spin_bh=.true.
  logical ::vrel_merge=.true.
  logical ::bhspinmerge=.true.
  logical ::selfgrav=.true.
  logical ::mad_jet=.false.
  logical ::force_exact_mseed=.false. ! Enforce Mass Criterion (MT)
  logical ::force_accretion=.false. ! Enforce Eddington-capped accretion (MT)
  logical ::no_accretion=.false. ! with grow_jeans=.false. remove BH accretion (HP)
  logical ::maximum_accretion=.false.  !Accrete all gas within the BH's accretion region. A type of mass flux accretion using the sink particle
  logical ::fix_smbh_position=.false.         ! Allows the sink particle to be fixed relative to the box
  logical ::momentum_feedback=.false.
  logical :: drag_part = .false. !activate the friction from stars/DM (HP)

  ! Output times
  real(dp),dimension(1:MAXOUT)::aout=1.1       ! Output expansion factors
  real(dp),dimension(1:MAXOUT)::tout=0.0       ! Output times

  ! Movie
  integer,parameter::NMOV=5
  integer::imovout=0             ! Increment for output times
  integer::imov=1                ! Initialize
  real(kind=8)::tstartmov=0.,astartmov=0.
  real(kind=8)::tendmov=0.,aendmov=0.
  real(kind=8),allocatable,dimension(:)::amovout,tmovout
  logical::movie=.false.
  integer::nw_frame=512 ! prev: nx_frame, width of frame in pixels
  integer::nh_frame=512 ! prev: ny_frame, height of frame in pixels
  integer::levelmax_frame=0
  real(kind=8),dimension(1:4*NMOV)::xcentre_frame=0d0
  real(kind=8),dimension(1:4*NMOV)::ycentre_frame=0d0
  real(kind=8),dimension(1:4*NMOV)::zcentre_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltax_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltay_frame=0d0
  real(kind=8),dimension(1:2*NMOV)::deltaz_frame=0d0
  real(kind=8),dimension(1:NMOV)::dtheta_camera=0d0
  real(kind=8),dimension(1:NMOV)::dphi_camera=0d0
  real(kind=8),dimension(1:NMOV)::theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::tstart_theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::tstart_phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::tend_theta_camera=0d0
  real(kind=8),dimension(1:NMOV)::tend_phi_camera=0d0
  real(kind=8),dimension(1:NMOV)::focal_camera=0d0
  real(kind=8),dimension(1:NMOV)::dist_camera=0d0
  real(kind=8),dimension(1:NMOV)::ddist_camera=0d0
  real(kind=8),dimension(1:NMOV)::smooth_frame=1d0
  real(kind=8),dimension(1:NMOV)::varmin_frame=-1d60
  real(kind=8),dimension(1:NMOV)::varmax_frame=1d60
  integer,dimension(1:NMOV)::ivar_frame=0
  logical,dimension(1:NMOV)::perspective_camera=.false.
  logical,dimension(1:NMOV)::zoom_only_frame=.false.
  character(LEN=NMOV)::proj_axis='z' ! x->x, y->y, projection along z
  character(LEN=6),dimension(1:NMOV)::shader_frame='square'
  character(LEN=10),dimension(1:NMOV)::method_frame='mean_mass'
#ifdef SOLVERmhd
  integer,dimension(0:NVAR+7)::movie_vars=0
  character(len=5),dimension(0:NVAR+7)::movie_vars_txt=''
#else
  integer,dimension(0:NVAR+3)::movie_vars=0
  character(len=5),dimension(0:NVAR+3)::movie_vars_txt=''
#endif

  ! Refinement parameters for each level
  real(dp),dimension(1:MAXLEVEL)::m_refine =-1.0 ! Lagrangian threshold
  real(dp),dimension(1:MAXLEVEL)::r_refine =-1.0 ! Radius of refinement region
  real(dp),dimension(1:MAXLEVEL)::x_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::y_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::z_refine = 0.0 ! Center of refinement region
  real(dp),dimension(1:MAXLEVEL)::exp_refine = 2.0 ! Exponent for distance
  real(dp),dimension(1:MAXLEVEL)::a_refine = 1.0 ! Ellipticity (Y/X)
  real(dp),dimension(1:MAXLEVEL)::b_refine = 1.0 ! Ellipticity (Z/X)
  real(dp)::var_cut_refine=-1.0 ! Threshold for variable-based refinement
  real(dp)::mass_cut_refine=-1.0 ! Mass threshold for particle-based refinement
  integer::ivar_refine=-1 ! Variable index for refinement
  logical::sink_refine=.false. ! Fully refine on sink particles

  ! Initial condition files for each level
  logical::multiple=.false.
  character(LEN=80),dimension(1:MAXLEVEL)::initfile=' '
  character(LEN=20)::filetype='ascii'

  ! Initial condition regions parameters
  integer,parameter::MAXREGION=100
  integer                           ::nregion=0
  character(LEN=10),dimension(1:MAXREGION)::region_type='square'
  real(dp),dimension(1:MAXREGION)   ::x_center=0.
  real(dp),dimension(1:MAXREGION)   ::y_center=0.
  real(dp),dimension(1:MAXREGION)   ::z_center=0.
  real(dp),dimension(1:MAXREGION)   ::length_x=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_y=1.E10
  real(dp),dimension(1:MAXREGION)   ::length_z=1.E10
  real(dp),dimension(1:MAXREGION)   ::exp_region=2.0

  ! Boundary conditions parameters
  integer,parameter::MAXBOUND=100
  logical                           ::simple_boundary=.false.
  integer                           ::nboundary=0
  integer                           ::icoarse_min=0
  integer                           ::icoarse_max=0
  integer                           ::jcoarse_min=0
  integer                           ::jcoarse_max=0
  integer                           ::kcoarse_min=0
  integer                           ::kcoarse_max=0
  integer ,dimension(1:MAXBOUND)    ::boundary_type=0
  integer ,dimension(1:MAXBOUND)    ::ibound_min=0
  integer ,dimension(1:MAXBOUND)    ::ibound_max=0
  integer ,dimension(1:MAXBOUND)    ::jbound_min=0
  integer ,dimension(1:MAXBOUND)    ::jbound_max=0
  integer ,dimension(1:MAXBOUND)    ::kbound_min=0
  integer ,dimension(1:MAXBOUND)    ::kbound_max=0
  logical                           ::no_inflow=.false.

  !Number of processes sharing one token
  !Only one process can write at a time in an I/O group
  integer::IOGROUPSIZE=0          ! Main snapshot
  integer::IOGROUPSIZECONE=0       ! Lightcone
  integer::IOGROUPSIZEREP=0        ! Subfolder size
  logical::withoutmkdir=.false.    !If true mkdir should be done before the run
  logical::print_when_io=.false.   !If true print when IO
  logical::synchro_when_io=.false. !If true synchronize when IO

  ! Stochastic feedback:
  real(dp)::SN_dT2_min = 0d0        ! Temperature increase for stochastic SNe
  ! Velocity kick maximum for new stellar particles (in random direction)
  real(dp)::SF_kick_kms = -1d0
  ! Keep track of densities at which stellar particles are born and go SN:
  logical::write_stellar_densities=.false.
  ! Kimm feedback stuff:
  ! Efficiency of stellar feedback energy used to heat up/blow out the gas:
  !  real(dp)::eff_sfbk=1.0D0
  real(dp)::E_SNII=1d51        ! different from ESN used in feedback.f90
  integer ::loading_type = 0 ! 0: uniform, 1: turbulence-based
  ! Realistic time delay for individual star particle. t_sne is oldest age to consider if set:
  logical ::sn2_real_delay=.false.
  logical ::use_initial_mass=.false. ! read/write initial mass of particles
  ! Activate mechanical feedback (cannot use with star_particle_winds):
  integer :: mechanical_feedback=0
  logical :: mechanical_geen = .false.
  logical :: log_mfb=.false.
  logical :: log_mfb_mega=.false.
  real(dp) :: A_SN=3d5
  real(dp) :: expN_SN=-2d0/17d0
  real(dp) :: A_SN_geen = 5d5
  real(dp) :: expN_SN_boost = -0.15    ! parameter needed to make sure
  real(dp) :: expE_SN_boost = 16./17.  ! that pSN_geen = A_SN_Geen during an adiabatic phase
  real(dp) :: porosity = 1.0d0

  ! Trick to fill a spherical region with passive scalar
  integer  :: level_zoom=0       ! Zoom level
  real(dp) :: xzoom=0.5       ! x position
  real(dp) :: yzoom=0.5       ! y position
  real(dp) :: zzoom=0.5       ! z position
  real(dp) :: rzoom=0.0       ! radius of the zoom

  ! Star formation stuff:
  ! character(len=10)::star_maker='density' ! density,hopkins, cen, padoan
  ! real(dp)::T2thres_SF=1d10  ! Temperature threshold
  real(dp)::fstar_min=1d0     ! Mstar,min = nH*dx_min^3*fstar_min
  real(dp)::M_SNII=10d0       ! Mean progenitor mass of the TypeII SNe
  real(dp)::sf_lamjt=1d0      ! Number of cells below which the Jeans length trigger SF
  character(len=8 )::star_imf=''          ! salpeter,kroupa,(chabrier)
  real(dp)::sf_mach_threshold=2.0d0       ! Minimal Mach number for turbulent SF
!!$  ! Density above which SF will occur regardless of the kinematic condition (dark cloud):
!!$  real(dp)::n_dc=1d10
!!$  ! Density above which Hopkins and Cen SF routines are evaluated:
!!$  real(dp)::n_gmc=1D-1
  ! Star particle mass in units of the number of SN:
  integer ::nsn2mass=-1


end module amr_parameters
