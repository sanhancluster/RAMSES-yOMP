module utils
  ! A module containing a few reusable functions
  use amr_parameters, only : icoarse_min, icoarse_max, jcoarse_min, jcoarse_max, kcoarse_min, kcoarse_max, &
       boxlen, nx, ny, nz, dp, ndim, nvector
  use random, only : IRandNumSize, ranf

  implicit none

  private

  ! Constants
  type Constants_t
     real(dp) :: c = 29979245800._dp  ! cm / s
  end type Constants_t

  type Units_t
     real(dp) :: pc  = 3.0856775814913673d18 ! cm
     real(dp) :: kpc = 3.0856775814913673d21 ! cm
     real(dp) :: Mpc = 3.0856775814913673d24 ! cm
  end type Units_T

  type(Constants_t) :: constants
  type(Units_t)     :: cgs

  real(dp) :: x_box(3), x_half(3), scale, xbound(3)
  integer  :: nx_loc, skip_loc(3)

  ! Public parameters
  real(dp), parameter, public :: pi = atan(1._dp) * 4._dp
  integer, parameter, public :: AGN_integration_depth = 6
  integer, parameter, public :: AGN_VOLUME_INTEGRATION_MC = 1, AGN_VOLUME_INTEGRATION_bisect = 2
  integer, parameter, public :: AGN_integration = AGN_VOLUME_INTEGRATION_bisect

  public :: distance3d, draw_normal, find_grid_containing, normalize_position, cell2grid, grid_level
  public :: debug_grid, find_grid_containing_max
  public :: constants, cgs


contains
  subroutine init_utils()
    logical, save :: called = .false.
    if (called) return

    called = .true.
    xbound = [real(nx, dp), real(ny, dp), real(nz, dp)]

    nx_loc = (icoarse_max - icoarse_min + 1)
    scale = boxlen/dble(nx_loc)

    x_half = scale*xbound/2.0
    x_box  = scale*xbound

    skip_loc(1) = int(icoarse_min)
    skip_loc(2) = int(jcoarse_min)
    skip_loc(3) = int(kcoarse_min)
  end subroutine init_utils

  subroutine distance3d(x1, y1, z1, x2, y2, z2, dx, dy, dz, ignore_periodicity)
    ! Compute the distance between two positions, taking care of periodicity
    real(dp), intent(in) :: x1, y1, z1
    real(dp), intent(in) :: x2, y2, z2
    logical, intent(in), optional :: ignore_periodicity

    real(dp), intent(out) :: dx, dy, dz

    logical :: ignore

    call init_utils()

    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1

    if (present(ignore_periodicity)) then
       ignore = ignore_periodicity
    else
       ignore = .false.
    end if

    if (ignore) return

    if (dx > x_half(1)) then
       dx = dx - x_box(1)
    else if (dx < -x_half(1)) then
       dx = dx + x_box(1)
    end if

    if (dy > x_half(2)) then
       dy = dy - x_box(2)
    else if (dy < -x_half(2)) then
       dy = dy + x_box(2)
    end if

    if (dz > x_half(3)) then
       dz = dz - x_box(3)
    else if (dz < -x_half(3)) then
       dz = dz + x_box(3)
    end if

  end subroutine distance3d

  subroutine normalize_position(pos, np)
    ! Normalize positions assuming periodical boundaries
    real(dp), dimension(1:nvector, 1:ndim), intent(inout) :: pos
    integer, intent(in) :: np

    integer :: idim, ipos

    call init_utils()

    do idim = 1, ndim
       do ipos = 1, np
          if (pos(ipos, idim) < 0._dp) then
             pos(ipos, idim) = pos(ipos, idim) + x_box(idim)
          else if (pos(ipos, idim) > x_box(idim)) then
             pos(ipos, idim) = pos(ipos, idim) - x_box(idim)
          end if
       end do
    end do

  end subroutine normalize_position

  ! elemental real(dp) function erfinv(x) result(p)
  !   ! Compute the inverse error function
  !   !
  !   ! From Giles, Mike. (2012). Approximating the Erfinv Function. GPU
  !   ! Computing                        Gems                       Jade
  !   ! Edition. . 10.1016/B978-0-12-385963-1.00010-1.
  !   real(dp), intent(in) :: x

  !   real(dp) :: w

  !   w = - log((1.0_dp - x) * (1.0_dp + x))
  !   if ( w < 6.250000_dp ) then
  !      w = w - 3.125000_dp
  !      p = -3.6444120640178196996d-21
  !      p =  -1.685059138182016589d-19    + p*w;
  !      p =  1.2858480715256400167d-18    + p*w;
  !      p =   1.115787767802518096d-17    + p*w;
  !      p =  -1.333171662854620906d-16    + p*w;
  !      p =  2.0972767875968561637d-17    + p*w;
  !      p =  6.6376381343583238325d-15    + p*w;
  !      p = -4.0545662729752068639d-14    + p*w;
  !      p = -8.1519341976054721522d-14    + p*w;
  !      p =  2.6335093153082322977d-12    + p*w;
  !      p = -1.2975133253453532498d-11    + p*w;
  !      p = -5.4154120542946279317d-11    + p*w;
  !      p =   1.051212273321532285d-09    + p*w;
  !      p = -4.1126339803469836976d-09    + p*w;
  !      p = -2.9070369957882005086d-08    + p*w;
  !      p =  4.2347877827932403518d-07    + p*w;
  !      p = -1.3654692000834678645d-06    + p*w;
  !      p = -1.3882523362786468719d-05    + p*w;
  !      p =   0.0001867342080340571352_dp + p*w;
  !      p = -0.00074070253416626697512_dp + p*w;
  !      p =  -0.0060336708714301490533_dp + p*w;
  !      p =     0.24015818242558961693_dp + p*w;
  !      p =      1.6536545626831027356_dp + p*w;
  !   else if ( w < 16._dp ) then
  !      w = sqrt(w) - 3.25_dp;
  !      p =   2.2137376921775787049d-09;
  !      p =   9.0756561938885390979d-08    + p*w;
  !      p =  -2.7517406297064545428d-07    + p*w;
  !      p =   1.8239629214389227755d-08    + p*w;
  !      p =   1.5027403968909827627d-06    + p*w;
  !      p =   -4.013867526981545969d-06    + p*w;
  !      p =   2.9234449089955446044d-06    + p*w;
  !      p =   1.2475304481671778723d-05    + p*w;
  !      p =  -4.7318229009055733981d-05    + p*w;
  !      p =   6.8284851459573175448d-05    + p*w;
  !      p =   2.4031110387097893999d-05    + p*w;
  !      p =   -0.0003550375203628474796_dp + p*w;
  !      p =   0.00095328937973738049703_dp + p*w;
  !      p =   -0.0016882755560235047313_dp + p*w;
  !      p =    0.0024914420961078508066_dp + p*w;
  !      p =   -0.0037512085075692412107_dp + p*w;
  !      p =     0.005370914553590063617_dp + p*w;
  !      p =       1.0052589676941592334_dp + p*w;
  !      p =       3.0838856104922207635_dp + p*w;
  !   else
  !      w = sqrt(w) - 5._dp
  !      p = -2.7109920616438573243d-11;
  !      p = -2.5556418169965252055d-10    + p*w;
  !      p =  1.5076572693500548083d-09    + p*w;
  !      p = -3.7894654401267369937d-09    + p*w;
  !      p =  7.6157012080783393804d-09    + p*w;
  !      p = -1.4960026627149240478d-08    + p*w;
  !      p =  2.9147953450901080826d-08    + p*w;
  !      p = -6.7711997758452339498d-08    + p*w;
  !      p =  2.2900482228026654717d-07    + p*w;
  !      p = -9.9298272942317002539d-07    + p*w;
  !      p =  4.5260625972231537039d-06    + p*w;
  !      p = -1.9681778105531670567d-05    + p*w;
  !      p =  7.5995277030017761139d-05    + p*w;
  !      p = -0.00021503011930044477347    + p*w;
  !      p = -0.00013871931833623122026_dp + p*w;
  !      p =      1.0103004648645343977_dp + p*w;
  !      p =      4.8499064014085844221_dp + p*w;
  !   end if
  ! end function erfinv

  subroutine draw_normal(seed, rand, vmin, vmax)
    use random, only: gaussdev
    ! Draw a number from a normal distribution with mean 0 and
    ! variance 1, with optional boundaries
    integer, intent(in), dimension(IRandNumSize) :: seed

    real(dp), intent(out) :: rand

    real(dp), intent(in), optional :: vmin, vmax

    call gaussdev(seed, rand)

    ! Redraw until in the boundary
    if (present(vmin) .and. present(vmax)) then
       do while (rand < vmin .or. rand > vmax)
          call gaussdev(seed, rand)
       end do
    else if (present(vmin)) then
       do while (rand < vmin)
          call gaussdev(seed, rand)
       end do
    else if (present(vmax)) then
       do while (rand > vmax)
          call gaussdev(seed, rand)
       end do
    end if

  end subroutine draw_normal

  subroutine find_grid_containing_max(pos, ind_grid, ind_cell, ind_level, npart, maxlevel)
    ! Find the finest grid and cell containing the positions.
    use amr_commons, only : levelmin, ncoarse, son, xg
    use amr_parameters, only : ndim, ngridmax, nvector
    real(dp), intent(in), dimension(1:nvector, 1:ndim) :: pos
    integer, intent(in) :: npart, maxlevel

    integer, intent(out), dimension(1:nvector) :: ind_grid, ind_cell, ind_level

    integer :: igrid, ipart, idim, ilevel, istart, ind

    real(dp) :: xgrid(3)

    call init_utils()

    ind_grid(1:npart) = 0; ind_cell(1:npart) = 0; ind_level(1:npart) = 0

    ! Find coarse grid containing position
    istart = 1
    igrid = 1

    ! Initialize at level = 1
    xgrid = (xg(igrid, :) - skip_loc) * scale

    if (any(pos(:, :) < 0._dp)) then
       print*, 'Got a position < 0'
       stop
    else if (any(pos(:, :) > boxlen)) then ! assuming x_box == [1, 1, 1]
       print*, 'Got a position > 1'
       stop
    end if

    do ipart = 1, npart
       ind = 1
       ! Loop over cells to compute cell containing object
       do idim = 1, ndim
          if (pos(ipart, idim) > xgrid(idim)) then
             ind = ind + 2**(idim-1)
          end if
       end do

       ! Compute cell and grid indexes
       ind_grid(ipart)  = igrid
       ind_cell(ipart)  = igrid + ncoarse + (ind - 1) * ngridmax
       ind_level(ipart) = levelmin
    end do

    ! Loop over levels
    do ipart = 1, npart
       level_loop: do ilevel = istart+1, maxlevel+1
          igrid = son(ind_cell(ipart))

          ! True if cell is refined → go down the tree
          if (igrid > 0) then
             xgrid = (xg(igrid, :) - skip_loc) * scale
             ind = 1

             ! Loop over cells to compute cell containing object
             do idim = 1, ndim
                if (pos(ipart, idim) > xgrid(idim))then
                   ind = ind + 2**(idim-1)
                end if
             end do

             ! Compute cell and grid indexes
             ind_grid(ipart)  = igrid
             ind_cell(ipart)  = igrid + ncoarse + (ind - 1) * ngridmax
             ind_level(ipart) = ilevel
          else
             exit level_loop
          end if
       end do level_loop
    end do ! particle loop

  end subroutine find_grid_containing_max

  subroutine find_grid_containing(pos, ind_grid, ind_cell, ind_level, npart)
    use amr_commons, only: nlevelmax
    ! Find the finest grid and cell containing the positions.
    real(dp), intent(in), dimension(1:nvector, 1:ndim) :: pos
    integer, intent(in) :: npart

    integer, intent(out), dimension(1:nvector) :: ind_grid, ind_cell, ind_level
    call find_grid_containing_max(pos, ind_grid, ind_cell, ind_level, npart, nlevelmax)
  end subroutine find_grid_containing

  elemental function cell2grid(icell) result(igrid)
    use amr_commons, only : ngridmax, ncoarse
    ! Return the index of the grid containing the cell

    integer, intent(in) :: icell
    integer :: igrid

    integer :: ipos

    ipos = (icell - ncoarse - 1) / ngridmax + 1
    igrid = icell - ncoarse - (ipos-1)*ngridmax

  end function cell2grid

  function grid_level(igrid) result(ilevel)
    use amr_commons, only : father
    integer, intent(in) :: igrid
    integer              :: ilevel

    integer :: icell_father, igrid_tmp

    ilevel = 0
    if (igrid <= 0) return

    igrid_tmp = igrid

    ! Walk up the tree until reaching the root grid
    do while (igrid_tmp > 0)
       icell_father = father(igrid_tmp)
       igrid_tmp = cell2grid(icell_father)
       ilevel = ilevel + 1
    end do

  end function grid_level

  function debug_grid(igrid) result(res)
    use amr_commons
    implicit none
    integer, intent(in) :: igrid
    logical :: res

    if ( &
         (myid == 1 .and. igrid == 2210) .or. (myid == 1 .and. igrid == 3) &
         .or. (myid == 4 .and. igrid == 5648)) then
       write(*, '(i3,"-",i5)', advance='no') myid, igrid
       res = .true.
    else
       res = .false.
    end if
  end function debug_grid

end module utils
