module pm_commons

  use amr_parameters
  use pm_parameters
  use random

  implicit none

  ! Sink particle related arrays
  real(dp), allocatable, dimension(:) :: msink, r2sink, v2sink, c2sink, oksink_new, oksink_all, tsink
  real(dp), allocatable, dimension(:) :: msink_new, msink_all, r2k, v2sink_new, c2sink_new, tsink_new, tsink_all
  real(dp), allocatable, dimension(:) :: v2sink_all, c2sink_all
  real(dp), allocatable, dimension(:) :: dMBHoverdt, dMEdoverdt, wdens, wvol, wc2
  real(dp), allocatable, dimension(:) :: wdens_new, wvol_new, wc2_new, total_volume
  real(dp), allocatable, dimension(:,:) :: wmom, wmom_new
  real(dp), allocatable, dimension(:,:) :: vsink, vsink_new, vsink_all
  real(dp), allocatable, dimension(:,:) :: xsink, xsink_new, xsink_all
  real(dp), allocatable, dimension(:,:) :: weighted_density, weighted_volume, weighted_c2
  real(dp), allocatable, dimension(:,:) :: jsink, jsink_new, jsink_all
  real(dp), allocatable, dimension(:) :: dMBH_coarse, dMEd_coarse, dMsmbh, dMBH_coarse_new
  real(dp), allocatable, dimension(:) :: dMEd_coarse_new, dMsmbh_new, dMBH_coarse_all, dMEd_coarse_all, dMsmbh_all
  real(dp), allocatable, dimension(:) :: Esave, Esave_new, Esave_all, Efeed, Efeed_new, Efeed_all
  ! AGNRT
  real(dp),allocatable,dimension(:)::LAGN_coarse
  real(dp),allocatable,dimension(:)::dMeff_coarse, dMeff_coarse_new, dMeff_coarse_all
  real(dp),allocatable,dimension(:,:)::lumfrac_AGN
  !/AGNRT
  real(dp), allocatable, dimension(:,:,:) :: weighted_momentum
  real(dp), allocatable, dimension(:,:,:) :: sink_stat, sink_stat_all
  real(dp), allocatable, dimension(:) :: c_avgptr, v_avgptr, d_avgptr
  real(dp), allocatable, dimension(:) :: spinmag, spinmag_new, spinmag_all
  real(dp), allocatable, dimension(:,:) :: bhspin, bhspin_new, bhspin_all
  real(dp), allocatable, dimension(:) :: eps_sink
  real(dp), allocatable, dimension(:) :: rg_scale    ! Gravitational scale radius of the black hole

  integer , allocatable, dimension(:) :: idsink, idsink_new, idsink_all
  ! Particles dynamical friction (HP)
  real(dp), allocatable, dimension(:,:,:) :: v_background, vrel_sink
  integer , allocatable, dimension(:,:) :: n_background
  real(dp), allocatable, dimension(:,:) :: m_background, vrel_sink_norm
  real(dp), allocatable, dimension(:,:) :: mass_lowspeed_background, fact_fast_background
  real(dp), allocatable, dimension(:,:,:) :: v_DFnew, v_DFall
  real(dp), allocatable, dimension(:,:,:,:) :: v_DF, v_DFnew_all
  real(dp), allocatable, dimension(:,:) :: mass_DFnew, mass_DFall
  real(dp), allocatable, dimension(:,:) :: fact_fastnew, fact_fastall
  real(dp), allocatable, dimension(:,:) :: mass_lowspeednew, mass_lowspeedall
  integer , allocatable, dimension(:,:) :: n_partnew, n_partall
  integer , allocatable, dimension(:, :, :) :: n_part, n_partnew_all
  real(dp), allocatable, dimension(:,:,:) :: mass_DF, fact_fast, mass_lowspeed
  real(dp), allocatable, dimension(:,:,:) :: mass_DFnew_all
  real(dp), allocatable, dimension(:,:,:) :: mass_lowspeednew_all, fact_fastnew_all
  real(dp), allocatable, dimension(:) :: most_massive_sink
  integer , allocatable, dimension(:) :: sink_cell
  integer::ncloud_sink                       !Number of cloud particles
  !/Particles dynamical friction (HP)
  integer :: nindsink = 0


  ! Particles related arrays
  real(dp), allocatable, dimension(:,:) :: xp        ! Positions
  real(dp), allocatable, dimension(:,:) :: vp        ! Velocities
  real(dp), allocatable, dimension(:)   :: mp        ! Masses
  integer,  allocatable, dimension(:)   :: move_flag ! Move flag (for particles), >0 means don't move!
  real(dp), allocatable, dimension(:)   :: mp0       ! Initial masses (for TK multiple SN)
#ifdef OUTPUT_PARTICLE_POTENTIAL
  real(dp), allocatable, dimension(:)   :: ptcl_phi  ! Potential of particle added by AP for output purposes
#endif
  real(dp), allocatable, dimension(:)   :: tp       ! Birth epoch
  real(dp), allocatable, dimension(:,:) :: weightp  ! weight of cloud parts for sink accretion only
  real(dp), allocatable, dimension(:)   :: zp       ! Birth metallicity
  real(dp), allocatable, dimension(:)   :: tmpp     ! Working array
  integer,  allocatable, dimension(:)   :: itmpp    ! Working array
  integer,  allocatable, dimension(:)   :: partp    ! Particle parent (for tracers only)
  integer, allocatable, dimension(:)  :: nextp     ! Next particle in list
  integer, allocatable, dimension(:)  :: prevp     ! Previous particle in list
  integer, allocatable, dimension(:)  :: levelp    ! Current level of particle
  integer(i8b), allocatable, dimension(:) :: idp   ! Identity of particle
  real(dp),allocatable,dimension(:)  ::st_n_tp  ! Gas density at birth epoch         !SD
  real(dp),allocatable,dimension(:)  ::st_n_sn  ! Gas density at SN epoch            !SD
  real(dp),allocatable,dimension(:)  ::st_e_sn  ! SN energy injected                 !SD

  ! Tree related arrays
  integer, allocatable, dimension(:)   :: headp    ! Head particle in grid
  integer, allocatable, dimension(:)   :: tailp    ! Tail particle in grid
  integer, allocatable, dimension(:)   :: numbp    ! Number of particles in grid
  ! Global particle linked lists
  integer :: headp_free, tailp_free, numbp_free = 0, numbp_free_tot = 0
  ! Local and current seed for random number generator
  integer, dimension(IRandNumSize) :: localseed = -1
  integer, dimension(IRandNumSize) :: tracer_seed = -1

  ! Particle types
  integer, parameter :: NFAMILIES=5
  integer(1),parameter :: FAM_DM=1, FAM_STAR=2, FAM_CLOUD=3, FAM_DEBRIS=4, FAM_OTHER=5, FAM_UNDEF=127
  integer(1),parameter :: FAM_TRACER_GAS=0
  integer(1),parameter :: FAM_TRACER_DM=-1, FAM_TRACER_STAR=-2, FAM_TRACER_CLOUD=-3, FAM_TRACER_DEBRIS=-4, FAM_TRACER_OTHER=-5

  ! Customize here for particle tags within particle types (e.g. different kind of stars).
  ! Note that the type should be integer(1) (1 byte integers) for memory concerns.
  ! Also don't forget to create a function is_<type>_<tag>. See the wiki for a more complete example.
  ! By default, the tag is always 0.
  integer(1),parameter :: TAG_STAR_ACTIVE=1
  integer(1),parameter :: TAG_CLOUD_CENTRAL=1

  ! Particle keys for outputing. They should match the above particle
  ! types, except for 'under' family
  character(len=13), dimension(-NFAMILIES:NFAMILIES), parameter :: particle_family_keys = (/ &
       ' other_tracer', 'debris_tracer', ' cloud_tracer', '  star_tracer', ' other_tracer', &
       '   gas_tracer', &
       '           DM', '         star', '        cloud', '       debris', '        other'/)

  type(part_t), allocatable, dimension(:) :: typep  ! Particle type array

contains
  function cross(a, b)
    use amr_parameters, only:dp
    real(dp), dimension(1:3) :: a, b
    real(dp), dimension(1:3) :: cross
    ! computes the cross product c= a x b
    cross(1) = a(2)*b(3)-a(3)*b(2)
    cross(2) = a(3)*b(1)-a(1)*b(3)
    cross(3) = a(1)*b(2)-a(2)*b(1)
  end function cross

  elemental logical pure function is_DM(typep)
    type(part_t), intent(in) :: typep
    is_DM = typep%family == FAM_DM
  end function is_DM

  elemental logical pure function is_not_DM(typep)
    ! Check that the particle is not DM and not a tracer
    type(part_t), intent(in) :: typep
    is_not_DM = typep%family /= FAM_DM .and. is_not_tracer(typep)
  end function is_not_DM

  elemental logical pure function is_star(typep)
    type(part_t), intent(in) :: typep
    is_star = typep%family == FAM_STAR
  end function is_star

  elemental logical pure function is_cloud(typep)
    type(part_t), intent(in) :: typep
    is_cloud = typep%family == FAM_CLOUD
  end function is_cloud

  elemental logical pure function is_debris(typep)
    type(part_t), intent(in) :: typep
    is_debris = typep%family == FAM_DEBRIS
  end function is_debris

  elemental logical pure function is_tracer(typep)
    type(part_t), intent(in) :: typep
    is_tracer = typep%family <= 0
  end function is_tracer

  elemental logical pure function is_not_tracer(typep)
    type(part_t), intent(in) :: typep
    is_not_tracer = typep%family > 0
  end function is_not_tracer

  elemental logical pure function is_gas_tracer(typep)
    type(part_t), intent(in) :: typep
    is_gas_tracer = typep%family == FAM_TRACER_GAS
  end function is_gas_tracer

  elemental logical pure function is_star_tracer(typep)
    type(part_t), intent(in) :: typep
    is_star_tracer = typep%family == FAM_TRACER_STAR
  end function is_star_tracer

  elemental logical pure function is_cloud_tracer(typep)
    type(part_t), intent(in) :: typep
    is_cloud_tracer = typep%family == FAM_TRACER_CLOUD
  end function is_cloud_tracer

  elemental logical pure function is_star_active(typep)
    type(part_t), intent(in) :: typep
    is_star_active = (typep%family == FAM_STAR) .and. (typep%tag == TAG_STAR_ACTIVE)
  end function is_star_active

  elemental logical pure function is_central_cloud(typep)
    type(part_t), intent(in) :: typep
    is_central_cloud = (typep%family == FAM_CLOUD) .and. (typep%tag == TAG_CLOUD_CENTRAL)
  end function is_central_cloud

  elemental function part2int (part)
    ! Convert a particle into an integer
    ! This saves some space e.g. when communicating
    integer :: part2int
    type(part_t), intent(in) :: part

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    part2int = (int(part%family) + a) * b + (int(part%tag) + a)
  end function part2int

  elemental function int2part(index)
    ! Convert from an index to particle type
    type(part_t) :: int2part
    integer, intent(in) :: index

    ! This is the largest value for integer(1)
    integer, parameter :: a = 128, b = 2*a

    int2part%family = int(index / b - a, 1)
    int2part%tag = int(mod(index, b) - a, 1)
  end function int2part

  function props2type(idpii, tpii, mpii)
    use amr_commons
    use pm_parameters, only : part_t

    ! Converts from "old" ramses to "new" ramses
    !
    ! Here's the match, add yours here for backward compatibility purposes
    ! DM     tpii == 0
    ! stars  tpii != 0 and idpii > 0
    ! sinks  tpii != 0 and idpii < 0
    !
    ! This is mostly for support of GRAFFIC I/O.
    ! The reason we use idpii instead of idp is to prevent name clashes
    real(dp), intent(in) :: tpii, mpii
    integer, intent(in)  :: idpii

    type(part_t) :: props2type

    props2type%tag = 0
    if (tpii == 0) then
       if (idpii > 0) then
          props2type%family = FAM_DM
       else if (idpii < 0) then
          props2type%family = FAM_CLOUD
       end if
    else if (idpii > 0) then
       props2type%family = FAM_STAR
       props2type%tag = 0
    else if (idpii < 0) then
       props2type%family = FAM_STAR
       props2type%tag = TAG_STAR_ACTIVE
    else if (mpii == 0) then
       props2type%family = FAM_TRACER_GAS
    end if

  end function props2type

end module pm_commons
