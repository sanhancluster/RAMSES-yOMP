module geom_volumes
  use geom_types
  use geom_distances
  use geom_utils
  implicit none

  private
  logical, parameter :: debug = .false.
  integer :: MAX_DEPTH = 9, NTODRAW = 1000
  public :: BoxContains, CapsuleContains
  public :: BoxVolume, CapsuleVolume
  public :: CapsuleBoxVolume, CapsuleBoxIntegrate, CapsuleBoxIntegrateMC

  public :: set_depth, get_depth
  public :: set_ndraw, get_ndraw

  interface
     real(dp) function callback_fun(x)
       use geom_types, only: dp
       real(dp), intent(in), dimension(3) :: x
     end function callback_fun
  end interface

contains
  subroutine CapsuleBoxVolume(box, capsule, volume)
    ! Compute the volume of the intersection of a box with a capsule
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: volume

    real(dp) :: useless

    volume = 0._dp
    ! Some interaction (capsule in box or the opposite)
    call divide(box, capsule, 1, 0, volume, useless, dummy, force_refine=.false.)
  contains
    real(dp) function dummy(X)
      real(dp), intent(in) :: X(3)
      dummy = 1._dp
    end function dummy
  end subroutine CapsuleBoxVolume

  subroutine CapsuleBoxIntegrate(integrand, box, capsule, volume, integral)
    ! Compute the volume of the intersection of a box with a capsule
    !
    ! Note
    ! ----
    ! The integrand function is called with arguments in the frame of
    ! the capsule. The starting point is at [0, 0, 0], the ending
    ! point at [0, 0, L]
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: integral, volume

    procedure(callback_fun) :: integrand

    volume   = 0._dp
    integral = 0._dp
    ! Some interaction (capsule in box or the opposite)
    call divide(box, capsule, 1, 0, volume, integral, integrand, force_refine=.true.)
  end subroutine CapsuleBoxIntegrate

  subroutine CapsuleBoxIntegrateMC(integrand, box, capsule, volume, integral)
    ! Compute the volume of the intersection of a box with a capsule (using Monte Carlo sampling)
    !
    ! Note
    ! ----
    ! The integrand function is called with arguments in the frame of
    ! the capsule. The starting point is at [0, 0, 0], the ending
    ! point at [0, 0, L]
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: integral, volume

    procedure(callback_fun) :: integrand
    integer :: ndrawn
    real(dp) :: vec(3), tmp

    volume   = 0._dp
    integral = 0._dp

    do ndrawn = 1, NTODRAW
       ! Draw random number within the box
       call random_number(vec)
       vec = (vec - 0.5_dp) * box%extents                   ! Rescale by extents
       vec = vec(1)*box%u + vec(2)*box%v + vec(3)*box%w     ! Apply on axis
       vec = vec + box%origin                               ! Add origin
       if (CapsuleContains(capsule, vec)) then
          integral = integral + integrand(vec) / NTODRAW
          volume = volume + 1._dp / NTODRAW
       end if
    end do

    ! Rescale by box volume
    tmp = BoxVolume(box)
    volume = volume * tmp
    integral = integral * tmp

  end subroutine CapsuleBoxIntegrateMC


  recursive subroutine divide(box, capsule, axis, depth, volume, integral, callback, force_refine)
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    integer, intent(in) :: depth, axis
    procedure(callback_fun) :: callback
    logical, intent(in) :: force_refine

    real(dp), intent(inout) :: volume, integral

    type(Box_t) :: box1, box2
    real(dp) :: dir(3), d2
    integer :: naxis, nin, i, j, k

    real(dp) :: point(3), dV, u1(3), u2(3), u3(3)
    character(len=1) :: ikey(3)

    ! Compute distance from box to capsule
    call CapsuleBoxDistanceSquared(box, capsule, d2)
    if (d2 > 0) then
       if (debug) call write_padding(depth, 'disjoint')
    else
       ! Compute which point falls in capsule
       nin = 0
       do k = -1, 1, 2
          do j = -1, 1, 2
             do i = -1, 1, 2
                point = box%origin &
                     + i * box%u * box%extents(1) &
                     + j * box%v * box%extents(2) &
                     + k * box%w * box%extents(3)
                if (CapsuleContains(capsule, point)) then
                   nin = nin + 1
                end if
             end do
          end do
       end do

       ! Convert position in capsule's frame
       point = box%origin - capsule%segment%start

       ! Build direct basis (z-axis along capsule segment)
       u3 = capsule%segment%end - capsule%segment%start
       u3 = u3 / norm2(u3)
       u1 = [&
            -u3(2) - u3(3), &
             u3(1) - u3(3), &
             u3(1) + u3(2) &
             ]
       u1 = u1 / norm2(u1)
       u2 = vector_product(u3, u1)
       u2 = u2 / norm2(u2)
       ! Project point on basis
       point = [dot_product(u1, point), dot_product(u2, point), dot_product(u3, point)]

       if (nin == 8 .and. .not. force_refine) then
          dV = BoxVolume(box)
          if (debug) call write_padding(depth, 'all inside')
          volume = volume + dV
          integral = integral + callback(point) * dV
       else if (depth == MAX_DEPTH) then
          if (nin == 0) return
          ! Estimate volume by number of points in volume
          dV = BoxVolume(box) * nin / 8._dp

          volume = volume + dV
          integral = integral + callback(point) * dV
          if (debug) call write_padding(depth, 'max depth')
       else
          ! Split box in 2
          box1 = box
          box2 = box

          if      (axis == 1) then
             dir = box%u
          else if (axis == 2) then
             dir = box%v
          else  ! (axis == 3)
             dir = box%w
          end if

          box1%origin = box%origin - dir * box%extents(axis) / 2
          box1%extents(axis) = box%extents(axis) / 2

          box2%origin = box%origin + dir * box%extents(axis) / 2
          box2%extents(axis) = box%extents(axis) / 2

          naxis = next_axis(axis)
          ikey = (/'x', 'y', 'z'/)
          if (debug) call write_padding(depth, '> LEFT ' // ikey(axis))
          call divide(box1, capsule, naxis, depth+1, volume, integral, callback, force_refine)
          if (debug) call write_padding(depth, '> RIGHT ' // ikey(axis))
          call divide(box2, capsule, naxis, depth+1, volume, integral, callback, force_refine)
       end if
    end if

  end subroutine divide

  !----------------------------------------
  ! Containing functions
  !----------------------------------------
  logical function CapsuleContains(capsule, point) result(inside)
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(in) :: point(3)

    real(dp) :: t, d2

    call PointLineSegDistanceSquared(capsule%segment, point, t, d2)
    ! print*, t, sqrt(d2)
    inside = d2 <= capsule%r**2

  end function CapsuleContains

  logical function BoxContains(box, point) result(inside)
    type(Box_t), intent(in) :: box
    real(dp), intent(in) :: point(3)

    real(dp) :: vec(3)
    integer :: idim

    ! Project on 3 axis
    vec = (/ &
         dot_product(box%u, point) / box%extents(1), &
         dot_product(box%v, point) / box%extents(2), &
         dot_product(box%w, point) / box%extents(3)  &
         /)

    inside = .true.
    do idim = 1, 3
       if (vec(idim) < -1._dp .or. vec(idim) > 1._dp) then
          inside = .false.
          return
       end if
    end do

  end function BoxContains

  !----------------------------------------
  ! Volume functions
  !----------------------------------------
  real(dp) function CapsuleVolume(capsule) result(volume)
    type(Capsule_t), intent(in) :: capsule

    real(dp) :: vec(3)

    vec = capsule%segment%end - capsule%segment%start

    volume = 4*pi/3*capsule%r**3 + pi * capsule%r**2 * norm2(vec)
  end function CapsuleVolume

  real(dp) function BoxVolume(box) result(volume)
    type(Box_t), intent(in) :: box

    ! Each extent is the half extent
    volume = box%extents(1) * box%extents(2) * box%extents(3) * 8
  end function BoxVolume

  !----------------------------------------
  ! Setter/getter
  !----------------------------------------
  subroutine set_depth(depth)
    integer, intent(in) :: depth

    MAX_DEPTH = depth
  end subroutine set_depth

  subroutine get_depth(depth)
    integer, intent(out) :: depth
    depth = MAX_DEPTH
  end subroutine get_depth

  subroutine set_ndraw(val)
    integer, intent(in) :: val

    NTODRAW = val
  end subroutine set_ndraw

  subroutine get_ndraw(val)
    integer, intent(out) :: val

    val = NTODRAW
  end subroutine get_ndraw

end module geom_volumes
