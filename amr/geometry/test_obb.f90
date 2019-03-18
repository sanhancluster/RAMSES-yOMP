program testBoxAndCapsule
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit
  use geom_distances
  use geom_types
  use test_utils
  implicit none

  type(Box_t) :: box
  ! type(Line_t) :: line
  type(Segment_t) :: segment

  real(dp) :: d2

  call test_group('Line and Box')
  !----------------------------------------
  ! Lines parallel to edge
  !----------------------------------------
  call test_group('// Box face & Box edge')
  call reset(box, segment)
  segment%start = (/1, -1, 0/)
  segment%end = (/1, 1, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Passing through edge")

  call reset(box, segment)
  segment%start = (/2, -10, 0/)
  segment%end = (/2, 10, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 1._dp, 1d-14), "Parallel to edge, above face")

  call reset(box, segment)
  segment%start = (/2, 1, -10/)
  segment%end = (/2, 1, 10/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 1._dp, 1d-14), "Parallel to other edge, above face")

  call reset(box, segment)
  segment%start = (/2, 2, -10/)
  segment%end = (/2, 2, 10/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 2._dp, 1d-14), "Parallel to edge")

  call reset(box, segment)
  segment%start = (/-10, 2, 2/)
  segment%end = (/10, 2, 2/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 2._dp, 1d-14), "Line along segment")

  call reset(box, segment)
  segment%start = (/-1, 0, 0/) / 2._dp
  segment%end = (/1, 0, 0/) / 2._dp
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Within box")

  call test_group_close()

  !----------------------------------------
  ! Lines parallel to face only
  !----------------------------------------
  call test_group('// Box face only')
  call reset(box, segment)
  segment%start = (/1, 1, -10/)
  segment%end = (/1, 2, 10/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Passing through node")

  call reset(box, segment)
  segment%start = (/0, -1, -10/)
  segment%end = (/0, 1, 10/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Passing through face")

  call reset(box, segment)
  segment%start = (/1, 0, 0/)
  segment%end = (/10, 0, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Touching face")

  call reset(box, segment)
  segment%start = (/2, 0, 0/)
  segment%end = (/10, 0, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 1._dp, 1d-14), "Pointing out")

  call reset(box, segment)
  segment%start = (/0, 0, 0/)
  segment%end = (/10, 0, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Crossing face")

  call reset(box, segment)
  segment%start = (/-1, 0, 0/) / 2._dp
  segment%end = (/1, 1, 0/) / 2._dp
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Within box")
  
  call test_group_close()

  !----------------------------------------
  ! Point like segments
  !----------------------------------------
  call test_group('Point-like')
  call reset(box, segment)
  segment%start = (/1, 0, 0/)
  segment%end = (/1, 0, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "On face")

  call reset(box, segment)
  segment%start = (/2, 1, 0/) / 2._dp
  segment%end = (/2, 1, 0/) / 2._dp
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "On edge")

  call reset(box, segment)
  segment%start = (/1, 1, 1/)
  segment%end = (/1, 1, 1/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "On node")

  call reset(box, segment)
  segment%start = (/0, 0, 0/)
  segment%end = (/0, 0, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Within box")

  call reset(box, segment)
  segment%start = (/2, 0, 0/)
  segment%end = (/2, 0, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 1._dp, 1d-14), "Outside box (x direction)")

  call reset(box, segment)
  segment%start = (/0, 2, 0/)
  segment%end = (/0, 2, 0/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 1._dp, 1d-14), "Outside box (y direction)")

  call reset(box, segment)
  segment%start = (/0, 0, 2/)
  segment%end = (/0, 0, 2/)
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 1._dp, 1d-14), "Outside box (z direction)")

  call test_group_close()

  !----------------------------------------
  ! Any line
  !----------------------------------------
  call test_group('Any segment')
  call reset(box, segment)
  segment%start = (/-1, -2, -3/) / 3._dp
  segment%end = (/1, 4, 8/) / 8._dp
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Inside box")

  call reset(box, segment)
  segment%start = (/-1, -2, -3/) / 3._dp
  segment%end = (/1, 4, 8/) / 6._dp
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0._dp, 1d-14), "Trough box")

  call reset(box, segment)
  segment%start = (/0, 0, 11/) / 10._dp
  segment%end = (/1, 2, 20/) / 10._dp
  call LineSegmentBoxDistanceSquared(box, segment, d2)
  call assert(is_close(d2, 0.01_dp, 1d-14), "Outside box")

  call test_group_close()

  call test_group_close()
contains

  subroutine reset(box, segment)
    type(Box_t), intent(inout) :: box
    type(Segment_t), intent(inout) :: segment
    box%origin = (/0, 0, 0/)
    box%u = (/1, 0, 0/)
    box%v = (/0, 1, 0/)
    box%w = (/0, 0, 1/)
    box%extents = (/1, 1, 1/)

    segment%start = (/0, 0, 0/)
    segment%end = (/0, 0, 0/)
  end subroutine reset

  subroutine print_details
    print*, 'Box:'
    print*, '----'
    print*, box%origin
    print*, box%u
    print*, box%v
    print*, box%w
    print*, ''
    print*, 'Segment:'
    print*, '--------'
    print*, segment%start
    print*, segment%end
    print*, ''
    print*, 'Intersection'
    print*, '------------'
    print*, 'd=', sqrt(d2)
  end subroutine print_details

end program testBoxAndCapsule

