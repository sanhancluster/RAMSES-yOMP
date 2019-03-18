module geom_types
  implicit none
  integer, parameter :: dp = selected_real_kind(15)

  type Point_t
     real(dp), dimension(3) :: pos
  end type Point_t

  type Box_t
     real(dp), dimension(3) :: origin ! Center
     real(dp), dimension(3) :: u, v, w
     real(dp), dimension(3) :: extents
  end type Box_t

  type Segment_t
     real(dp), dimension(3) :: start, end
  end type Segment_t

  type Line_t
     real(dp), dimension(3) :: origin ! Starting point
     real(dp), dimension(3) :: direction ! Direction
  end type Line_T

  type Capsule_t
     type(Segment_t) :: segment
     real(dp) :: r
  end type Capsule_t

end module geom_types
