program test_volume
  use geom_distances
  use geom_types
  use geom_volumes
  use geom_utils
  use test_utils
  implicit none

  type(Segment_t) :: segment
  type(Capsule_t) :: capsule
  type(Box_t) :: box
  type(Line_t) :: line
  real(dp) :: t, d2, V, integral

  !----------------------------------------
  ! Point to line
  !----------------------------------------
  call test_group('Point to line')
  call reset(box, segment, capsule)
  line%origin = (/0, 0, 0/)
  line%direction = (/1, 0, 0/)
  call PointLineDistanceSquared(line, (/0d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point at origin')

  call PointLineDistanceSquared(line, (/1d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point at tip')

  line%origin = (/0, 0, 0/)
  line%direction = (/1, 1, 1/)
  call PointLineDistanceSquared(line, (/1d0, 1d0, 1d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point at tip (2)')

  line%origin = (/0, 0, 0/)
  line%direction = (/1, 0, 0/)
  call PointLineDistanceSquared(line, (/0d0, 1d0, 0d0/), t, d2)
  call assert(is_close(d2, 1._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point on y-axis')

  call PointLineDistanceSquared(line, (/0d0, 0d0, 1d0/), t, d2)
  call assert(is_close(d2, 1._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point on z-axis')

  call PointLineDistanceSquared(line, (/0d0, 0d0, 1d0/), t, d2)
  call assert(is_close(d2, 1._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point on z-axis')

  call PointLineDistanceSquared(line, (/100d0, 10d0, 20d0/), t, d2)
  call assert(is_close(d2, 500._dp, 1d-14) .and. &
       is_close(t, 100._dp, 1d-14), 'Point far')

  line%origin = (/0, 0, 0/)
  line%direction = (/2, 0, 0/)
  call PointLineDistanceSquared(line, (/100d0, 10d0, 20d0/), t, d2)
  call assert(is_close(d2, 500._dp, 1d-14) .and. &
       is_close(t, 50._dp, 1d-14), 'Point far')
  call test_group_close()

  !----------------------------------------
  ! Point to segment
  !----------------------------------------
  call test_group('Point to segment')
  call reset(box, segment, capsule)
  call PointLineSegDistanceSquared(segment, (/0d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 0.5_dp, 1d-14), 'Point at origin')

  call PointLineSegDistanceSquared(segment, (/-1d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 0._dp, 1d-14), 'Point at start')

  call PointLineSegDistanceSquared(segment, (/1d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 0._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point at end')

  call PointLineSegDistanceSquared(segment, (/3d0, 0d0, 0d0/), t, d2)
  call assert(is_close(d2, 4._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point along axis')

  call PointLineSegDistanceSquared(segment, (/0d0, 1d0, 2d0/), t, d2)
  call assert(is_close(d2, 5._dp, 1d-14) .and. &
       is_close(t, 0.5_dp, 1d-14), 'Point perp. to axis')

  call PointLineSegDistanceSquared(segment, (/1d0, 1d0, 2d0/), t, d2)
  call assert(is_close(d2, 5._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point perp. to axis (at tip)')

  call PointLineSegDistanceSquared(segment, (/2d0, 1d0, 2d0/), t, d2)
  call assert(is_close(d2, 6._dp, 1d-14) .and. &
       is_close(t, 1._dp, 1d-14), 'Point perp. to axis (beyond tip)')

  call test_group_close()

  !----------------------------------------
  ! CpasuleContains
  !----------------------------------------
  call test_group('Capsule contains')
  call reset(box, segment, capsule)
  call assert(CapsuleContains(capsule, (/0._dp, 0._dp, 0._dp/)), 'center of capsule')
  call assert(CapsuleContains(capsule, (/1._dp, 0._dp, 0._dp/)), 'in capsule (x-axis)')
  call assert(CapsuleContains(capsule, (/0._dp, 0.5_dp, 0._dp/)), 'in capsule (y-axis)')
  call assert(CapsuleContains(capsule, (/0._dp, 0._dp, 0.5_dp/)), 'in capsule (z-axis)')
  call assert(CapsuleContains(capsule, (/2._dp, 0._dp, 0._dp/)), 'side of capsule')
  call assert(.not. CapsuleContains(capsule, (/2.00001_dp, 0._dp, 0._dp/)), 'outside of capsule')

  call reset(box, segment, capsule)
  capsule%segment%start = (/-10, 0, 0/)
  capsule%segment%end = (/10, 0, 0/)
  capsule%r = 1.d0
  call assert(CapsuleContains(capsule, (/0._dp, 0._dp, 0._dp/)), 'center of capsule')
  call assert(CapsuleContains(capsule, (/1._dp, 0._dp, 0._dp/)), 'in capsule')
  call assert(CapsuleContains(capsule, (/2._dp, 0._dp, 0._dp/)), 'side of capsule')
  call assert(.not. CapsuleContains(capsule, (/-11.00001_dp, 0._dp, 0._dp/)), 'outside of capsule')

  call test_group_close()

  !----------------------------------------
  ! CpasuleContains
  !----------------------------------------
  call test_group('Box contains')
  call reset(box, segment, capsule)
  box%extents = (/2, 3, 4/)
  call assert(BoxContains(box, (/0._dp, 0._dp, 0._dp/)), 'center of box')
  call assert(BoxContains(box, (/2._dp, 0._dp, 0._dp/)), 'face-x')
  call assert(BoxContains(box, (/0._dp, 3._dp, 0._dp/)), 'face-y')
  call assert(BoxContains(box, (/0._dp, 0._dp, 4._dp/)), 'face-z')
  call assert(BoxContains(box, (/2._dp, 3._dp, 0._dp/)), 'edge-xy')
  call assert(BoxContains(box, (/2._dp, 0._dp, 4._dp/)), 'edge-xz')
  call assert(BoxContains(box, (/0._dp, 3._dp, 4._dp/)), 'edge-yz')
  call assert(.not. BoxContains(box, (/10._dp, 3._dp, 4._dp/)), 'outside')


  call test_group_close()

  !----------------------------------------
  ! Volumes
  !----------------------------------------
  call test_group('Volumes')

  call reset(box, segment, capsule)
  call assert(is_close(BoxVolume(box), 8._dp, 1d-14), 'box volume')
  capsule%segment%end = capsule%segment%start
  call assert(is_close(CapsuleVolume(capsule), 4*pi/3*capsule%r**3, 1d-14), &
       'capsule volume (sphere)')

  call reset(box, segment, capsule)
  call assert(is_close(CapsuleVolume(capsule), 4*pi/3*capsule%r**3 &
       + pi * capsule%r**2 * norm2(capsule%segment%end - capsule%segment%start), 1d-14), &
       'capsule volume')

  ! Compute intersection volume
  call reset(box, segment, capsule)
  capsule%r = 3
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, BoxVolume(box), 1d-14), 'box in capsule')

  call reset(box, segment, capsule)
  capsule%segment%start = (/-10, -2, 0/)
  capsule%segment%end = (/10, -2, 0/)
  capsule%r = 2._dp

  box%extents = (/2, 2, 2/)
  box%origin = (/0, 0, 0/)
  ! TODO: better constrain on precision
  call set_depth(6)
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, pi/2*capsule%r**2 * 4, 5d-2*V), 'Capsule on box side (depth 6)')

  call set_depth(9)
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, pi/2*capsule%r**2 * 4, 5d-2*V), 'Capsule on box side (depth 9)')

  call set_depth(12)
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, pi/2*capsule%r**2 * 4, 3d-2*V), 'Capsule on box side (depth 12)')

  call set_depth(15)
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, pi/2*capsule%r**2 * 4, 1.5d-2*V), 'Capsule on box side (depth 15)')

  call set_depth(18)
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, pi/2*capsule%r**2 * 4, 3d-3*V), 'Capsule on box side (depth 18)')

  call set_depth(21)
  call CapsuleBoxVolume(box, capsule, V)
  call assert(is_close(V, pi/2*capsule%r**2 * 4, 1.5d-3*V), 'Capsule on box side (depth 21)')

  call test_group_close()
  !----------------------------------------
  ! Volumes
  !----------------------------------------
  call test_group('Integrate')
  call reset(box, segment, capsule)
  call set_depth(9)
  capsule%r = 3
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, BoxVolume(box) * 2, 1d-14), 'box in capsule')

  call reset(box, segment, capsule)
  capsule%segment%start = (/-10, -2, 0/)
  capsule%segment%end = (/10, -2, 0/)
  capsule%r = 2._dp

  box%extents = (/2, 2, 2/)
  box%origin = (/0, 0, 0/)
  ! TODO: better constrain on precision
  call set_depth(6)
  call CapsuleBoxVolume(box, capsule, V)
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, pi/2*capsule%r**2 * 4 * 2, 5d-2*integral), 'Capsule on box side (depth 6)')

  call set_depth(9)
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, pi/2*capsule%r**2 * 4 * 2, 5d-2*integral), 'Capsule on box side (depth 9)')

  call set_depth(12)
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, pi/2*capsule%r**2 * 4* 2, 3d-2*integral), 'Capsule on box side (depth 12)')

  call set_depth(15)
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, pi/2*capsule%r**2 * 4 * 2, 1.5d-2*integral), 'Capsule on box side (depth 15)')

  call set_depth(18)
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, pi/2*capsule%r**2 * 4 * 2, 3d-3*integral), 'Capsule on box side (depth 18)')

  call set_depth(21)
  call CapsuleBoxIntegrate(dummyx2, box, capsule, V, integral)
  call assert(is_close(integral, pi/2*capsule%r**2 * 4 * 2, 1.5d-3*integral), 'Capsule on box side (depth 21)')

  call set_depth(3)
  call reset(box, segment, capsule)
  capsule%segment%start = (/-1, -1, -1/)
  capsule%segment%end = (/1, 1, 1/)
  capsule%r = 10._dp

  box%origin = (/0, 0, 0/)
  box%extents = (/1, 1, 1/)
  call CapsuleBoxIntegrate(tophalf, box, capsule, V, integral)
  call assert(is_close(integral, 4._dp, 4d-14), 'Non constant function')

  call test_summary()
contains

  subroutine reset(box, segment, capsule)
    type(Box_t), intent(inout) :: box
    type(Segment_t), intent(inout) :: segment
    type(Capsule_t), intent(inout) :: capsule

    box%origin = (/0, 0, 0/)
    box%u = (/1, 0, 0/)
    box%v = (/0, 1, 0/)
    box%w = (/0, 0, 1/)
    box%extents = (/1, 1, 1/)

    segment%start = (/-1, 0, 0/)
    segment%end = (/1, 0, 0/)

    capsule%segment = segment
    capsule%r = 1._dp

    call set_depth(3)
  end subroutine reset

  subroutine print_details
    print*, 'Box:'
    print*, '----'
    print*, box%origin
    print*, box%u
    print*, box%v
    print*, box%w
    print*, box%extents
    print*, ''
    print*, 'Capsule:'
    print*, '--------'
    print*, capsule%segment%start
    print*, capsule%segment%end
    print*, capsule%r
    print*, ''
    print*, 'Intersection'
    print*, '------------'
    print*, 'd=', sqrt(d2)
  end subroutine print_details

  real(dp) function dummyx2(X)
    real(dp), intent(in) :: X(3)
    dummyx2 = 2._dp
  end function dummyx2

  real(dp) function tophalf(X)
    ! This is only relevant to the test using it.
    real(dp), intent(in) :: X(3)

    if (X(3) > sqrt(3._dp)) then
       tophalf = 1._dp
    else
       tophalf = 0._dp
    end if
  end function tophalf
end program
