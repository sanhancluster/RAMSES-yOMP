module geom_distances
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit
  use geom_utils
  use geom_types
  implicit none
  private
  public :: LineBoxDistanceSquared, LineSegmentBoxDistanceSquared, PointBoxDistanceSquared, CapsuleBoxDistanceSquared
  public :: PointLineDistanceSquared, PointLineSegDistanceSquared

contains
  subroutine CapsuleBoxDistanceSquared(box, capsule, d2)
    type(Box_t), intent(in) :: box
    type(Capsule_t), intent(in) :: capsule
    real(dp), intent(out) :: d2

    call LineSegmentBoxDistanceSquared(box, capsule%segment, d2)

    ! Get min distance - capsule's sphere radius
    d2 = max(d2-capsule%r**2, 0._dp)

  end subroutine CapsuleBoxDistanceSquared

  subroutine LineSegmentBoxDistanceSquared(box, seg, d2)
    type(Box_t), intent(in) :: box
    type(Segment_t), intent(in) :: seg
    real(dp), intent(out) :: d2

    real(dp) :: t
    type(Line_t) :: line

    line%origin = seg%start
    line%direction = seg%end - seg%start

    call LineBoxDistanceSquared(box, line, d2, t)

    if (t < 0) then
       call PointBoxDistanceSquared(box, seg%start, d2)
    else if (t > 1) then
       call PointBoxDistanceSquared(box, seg%end, d2)
    end if
  end subroutine LineSegmentBoxDistanceSquared

  subroutine LineBoxDistanceSquared(box, line, d2, t)
    type(Box_t), intent(in) :: box
    type(Line_t), intent(in) :: line
    real(dp), intent(out) :: d2, t

    real(dp) :: Pprime(3), tmp(3), dprime(3)
    integer :: n0, axis
    logical :: ind(3)

    type(Box_t) :: newBox
    type(Line_t) :: newLine

    ! Project segment onto box own frame
    tmp = (line%origin - box%origin)
    Pprime = (/ &
         dot_product(tmp, box%u), &
         dot_product(tmp, box%v), &
         dot_product(tmp, box%w)  &
         /)
    dprime = (/ &
         dot_product(line%direction, box%u), &
         dot_product(line%direction, box%v), &
         dot_product(line%direction, box%w)  &
         /)

    if (dprime(1) < 0) then
       Pprime(1) = -Pprime(1)
       dprime(1) = -dprime(1)
    end if
    if (dprime(2) < 0) then
       Pprime(2) = -Pprime(2)
       dprime(2) = -dprime(2)
    end if
    if (dprime(3) < 0) then
       Pprime(3) = -Pprime(3)
       dprime(3) = -dprime(3)
    end if

    ! Store new Box and line
    newBox = box
    newBox%u = (/1, 0, 0/)
    newBox%v = (/0, 1, 0/)
    newBox%w = (/0, 0, 1/)

    newLine = line
    newLine%origin = Pprime
    newLine%direction = dprime

    ! Count number of null component
    ind = dprime > 0._dp
    n0 = count(.not. ind)
    if (n0 == 3) then
       call PointBoxDistanceSquared(newBox, newLine%origin, d2)
       t = 0
    else if (n0 == 2) then
       ! Select non null axis
       if (ind(1)) axis = 1
       if (ind(2)) axis = 2
       if (ind(3)) axis = 3
       call CaseTwoZeroComponents(newLine, newBox, axis, d2, t)
    else if (n0 == 1) then
       ! Select non null axis
       if (.not. ind(1)) axis = 1
       if (.not. ind(2)) axis = 2
       if (.not. ind(3)) axis = 3
       call CaseOneZeroComponents(newLine, newBox, axis, d2, t)
    else if (n0 == 0) then
       call CaseNoZeroComponents(newLine, newBox, d2, t)
    end if

  contains
    subroutine CaseTwoZeroComponents(line, box, X, d2, t)
      type(Box_t), intent(in) :: box
      type(Line_t), intent(in) :: line
      integer, intent(in) :: X

      real(dp), intent(out) :: d2, t

      real(dp) :: delta
      integer :: Y, Z

      ! In the comments, we assume axis to be x-axis
      ! Maps 2->3, 3->1 and 1->2
      Y = next_axis(X)
      Z = next_axis(Y)

      d2 = 0

      t = (box%extents(X) - line%origin(X)) / line%direction(X)
      ! Distance in y direction
      if (line%origin(Y) < -box%extents(Y)) then
         delta = line%origin(Y) + box%extents(Y)
         d2 = d2 + delta**2
      else if (line%origin(Y) > box%extents(Y)) then
         delta = line%origin(Y) - box%extents(Y)
         d2 = d2 + delta**2
      end if

      ! Distance in z direction
      if (line%origin(Z) < -box%extents(Z)) then
         delta = line%origin(Z) + box%extents(Z)
         d2 = d2 + delta**2
      else if (line%origin(Z) > box%extents(Z)) then
         delta = line%origin(Z) - box%extents(Z)
         d2 = d2 + delta**2
      end if

    end subroutine CaseTwoZeroComponents

    subroutine CaseOneZeroComponents(line, box, Z, d2, t)
      type(Box_t), intent(in) :: box
      type(Line_t), intent(in) :: line
      integer, intent(in) :: Z

      real(dp), intent(out) :: d2, t

      real(dp) :: prod0, prod1, ptMinusExtents(3), tmp, delta, invLSquared, inv
      integer :: X, Y

      d2 = 0

      ! In the comments, we assume axis to be the z axis
      ! Get next axis in order
      X = next_axis(Z)
      Y = next_axis(Y)

      ptMinusExtents = line%origin - box%extents

      prod0 = line%direction(Y) * ptMinusExtents(X)
      prod1 = line%direction(X) * ptMinusExtents(Y)

      if (prod0 >= prod1) then
         ! Line intersects x-axis of box
         ! Closest point is along bottom edge of right face of box
         tmp = line%origin(Y) + box%extents(Y)
         delta = prod0 - line%direction(X) * tmp

         if (delta >= 0) then
            ! No intersection, so compute distance
            invLSquared = 1._dp / (&
                 line%direction(X) * line%direction(X) + &
                 line%direction(Y) * line%direction(Y))
            d2 = d2 + delta*delta * invLSquared
            t = - (line%direction(X) * ptMinusExtents(X) + line%direction(Y) * tmp) * invLSquared
         else
            ! Line intersects box, distance is 0
            ! d2 = d2 + 0
            inv = 1 / line%direction(X)
            t = -ptMinusExtents(X) * inv
         end if
      else
         ! Line intersects y-axis of box
         ! Closest point is along top edge of left face of box
         tmp = line%origin(X) + box%extents(X)
         delta = prod1 - line%direction(Y) * tmp

         if (delta >= 0) then
            ! No intersection, so compute distance
            invLSquared = 1._dp / (&
                 line%direction(X) * line%direction(X) + &
                 line%direction(Y) * line%direction(Y))
            d2 = d2 + delta*delta * invLSquared
            t = - (line%direction(Y) * ptMinusExtents(Y) + line%direction(X) * tmp) * invLSquared
         else
            ! Line intersects box, distance is 0
            inv = 1 / line%direction(Y)
            t = -ptMinusExtents(Y) * inv
         end if
      end if

      ! Now consider the z-direction
      if (line%origin(Z) < -box%extents(Z)) then
         delta = line%origin(Z) + box%extents(Z)
         d2 = d2 + delta*delta
      else if (line%origin(Z) > box%extents(Z)) then
         delta = line%origin(Z) - box%extents(Z)
         d2 = d2 + delta*delta
      end if

    end subroutine CaseOneZeroComponents

    subroutine CaseNoZeroComponents(line, box, d2, t)
      type(Box_t), intent(in) :: box
      type(Line_t), intent(in) :: line

      real(dp), intent(out) :: d2, t

      real(dp) :: ptMinusExtents(3), dyEx, dxEy, dzEx, dxEz, dzEy, dyEz
      integer :: X=1, Y=2, Z=3

      ptMinusExtents = line%origin - box%extents
      dyEx = line%direction(Y) * ptMinusExtents(X)
      dxEy = line%direction(X) * ptMinusExtents(Y)

      if (dyEx >= dxEy) then
         dzEx = line%direction(Z) * ptMinusExtents(X)
         dxEz = line%direction(X) * ptMinusExtents(Z)

         if (dzEx>= dxEz) then
            ! Line intersects x = box.extent.x plane
            call Face(line, box, ptMinusExtents, X, d2, t)
         else
            ! Line intersects z = box.extent.z plane
            call Face(line, box, ptMinusExtents, Z, d2, t)
         end if
      else
         dzEy = line%direction(Z) * ptMinusExtents(Y)
         dyEz = line%direction(Y) * ptMinusExtents(Z)

         if (dzEy >= dyEz) then
            ! Line intersects y = box.extent.y plane
            call Face(line, box, ptMinusExtents, Y, d2, t)
         else
            ! Line intersects z = box.extent.z plane
            call Face(line, box, ptMinusExtents, Z, d2, t)
         end if
      end if

    end subroutine CaseNoZeroComponents

    subroutine Face(line, box, ptMinusExtents, axis, d2, t)
      type(Box_t), intent(in) :: box
      type(Line_t), intent(in) :: line
      real(dp), intent(in) :: ptMinusExtents(3)
      integer, intent(in) :: axis

      real(dp), intent(out) :: d2, t

      real(dp) :: ptPlusExtents(3), lSqr, tmp, delta, inverse
      integer :: X, Y, Z

      X = axis
      Y = next_axis(X)
      Z = next_axis(Y)

      ptPlusExtents = line%origin + box%extents

      d2 = 0

      if (line%direction(X) * ptPlusExtents(Y) >= line%direction(Y) * ptMinusExtents(X)) then
         ! Region 0, 5 or 4
         if (line%direction(X) * ptPlusExtents(Z) >= line%direction(Z) * ptMinusExtents(X)) then
            ! Region 0 - line intersects face
            ! Distance is 0
            inverse = 1.0_dp / line%direction(X)
            t = -ptMinusExtents(X) * inverse
         else
            ! Region 4 or 5
            lSqr = line%direction(X)**2 + line%direction(Z)**2
            tmp = lSqr * ptPlusExtents(Y) - line%direction(Y) * &
                 (line%direction(X) * ptMinusExtents(X) + line%direction(Z) * ptPlusExtents(Z))

            if (tmp <= 2 * lSqr * box%extents(Y)) then
               ! Region 4
               call region4code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    tmp, lSqr, delta, t, d2)
            else
               ! Region 5
               call region5code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    lSqr, delta, t, d2)
            end if
         end if
      else
         if (line%direction(X) * ptPlusExtents(Z) >= line%direction(Z) * ptMinusExtents(X)) then
            ! Region 1 or 2
            lSqr = line%direction(X)**2 + line%direction(Y)**2
            tmp = lSqr * ptPlusExtents(Z) - line%direction(Z) * &
                 (line%direction(X) * ptMinusExtents(X) + line%direction(Y) * ptPlusExtents(Y))
            if (tmp <= 2 * lSqr * box%extents(Z)) then
               ! Region 2
               call region2code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    tmp, lSqr, delta, t, d2)
            else
               ! Region 1
               call region1code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    lSqr, delta, t, d2)
            end if
         else
            lSqr = line%direction(X)**2 &
                 + line%direction(Z)**2
            tmp = lSqr * ptPlusExtents(Y) - line%direction(Y) * &
                 (line%direction(X) * ptMinusExtents(X) + line%direction(Z) * ptPlusExtents(Z))

            if (tmp >= 0) then
               ! Region 4 or 5
               if (tmp <= 2 * lSqr * box%extents(Y)) then
                  ! Region 4. Same as before (TODO)
                  call region4code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    tmp, lSqr, delta, t, d2)
               else
                  ! Region 5. Same as before (TODO)
                  call region5code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    lSqr, delta, t, d2)
               end if
            end if

            lSqr = line%direction(X)**2 + line%direction(Y)**2
            tmp = lSqr * ptPlusExtents(Z) - line%direction(Z) * &
                 (line%direction(X) * ptMinusExtents(X) + line%direction(Y) * ptPlusExtents(Y))

            if (tmp >= 0) then
               ! Region 1 or 2
               if (tmp <= 2 * lSqr * box%extents(Z)) then
                  ! Region 2. Same as before
                  call region2code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                    tmp, lSqr, delta, t, d2)
               else
                  ! Region 1. Same as before
                  call region1code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
                       lSqr, delta, t, d2)
               end if
               return ! d2
            end if

            ! Region 3
            lSqr = lSqr + line%direction(Y)**2
            delta = line%direction(X) * ptMinusExtents(X) &
                 + line%direction(Y) * ptPlusExtents(Y) &
                 + line%direction(Z) * ptPlusExtents(Z)
            t = - delta / lSqr
            d2 = d2 + ptMinusExtents(X)**2 &
                 + ptPlusExtents(Y)**2 &
                 + ptPlusExtents(Z)**2 &
                 + delta * t
         end if
      end if
    end subroutine Face

    subroutine region1code(line, ptMinusExtents, ptPlusExtents, X, Y, Z,&
         lSqr, delta, t, d2)
      type(Line_t), intent(in) :: line
      real(dp), intent(in) :: ptMinusExtents(3), ptPlusExtents(3)
      integer, intent(in) :: X, Y, Z

      real(dp), intent(inout) :: lSqr, delta, t, d2

      lSqr = lSqr + line%direction(Z)**2
      delta = line%direction(X) * ptMinusExtents(X) &
           + line%direction(Y) * ptMinusExtents(Y) &
           + line%direction(Z) * ptPlusExtents(Z)

      t = -delta / lSqr
      d2 = d2 + ptMinusExtents(X)**2 &
           + ptPlusExtents(Y)**2 &
           + ptMinusExtents(Z)**2 &
           + delta * t

    end subroutine region1code

    subroutine region2code(line, ptMinusExtents, ptPlusExtents, X, Y, Z, &
         tmp, lSqr, delta, t, d2)
      type(Line_t), intent(in) :: line
      real(dp), intent(in) :: ptMinusExtents(3), ptPlusExtents(3)
      integer, intent(in) :: X, Y, Z

      real(dp), intent(inout) :: tmp, lSqr, delta, t, d2

      tmp = ptPlusExtents(Z) - tmp / lSqr
      lSqr = lSqr + line%direction(Z)**2
      delta = line%direction(X) * ptMinusExtents(X) &
           + line%direction(Y) * ptPlusExtents(Y) &
           + line%direction(Z) * tmp
      t = -delta/ lSqr
      d2 = d2 + ptMinusExtents(X)**2 &
           + ptMinusExtents(Y)**2 &
           + tmp**2 &
           + delta * t
    end subroutine region2code

    subroutine region4code(line, ptMinusExtents, ptPlusExtents, X, Y, Z, &
         tmp, lSqr, delta, t, d2)
      type(Line_t), intent(in) :: line
      real(dp), intent(in) :: ptMinusExtents(3), ptPlusExtents(3)
      integer, intent(in) :: X, Y, Z

      real(dp), intent(inout) :: tmp, lSqr, delta, t, d2

      tmp = ptPlusExtents(Y) - tmp / lSqr
      lSqr = lSqr + line%direction(Y)**2
      delta = line%direction(X) * ptMinusExtents(X) &
           + line%direction(Y) * tmp &
           + line%direction(Z) * ptPlusExtents(Z)
      t = -delta / lSqr
      d2 = d2 + ptMinusExtents(X)**2 + tmp**2 + ptPlusExtents(Z)**2 + delta * t
    end subroutine region4code

    subroutine region5code(line, ptMinusExtents, ptPlusExtents, X, Y, Z, &
         lSqr, delta, t, d2)
      type(Line_t), intent(in) :: line
      real(dp), intent(in) :: ptMinusExtents(3), ptPlusExtents(3)
      integer, intent(in) :: X, Y, Z

      real(dp), intent(inout) :: lSqr, delta, t, d2

      lSqr = lSqr + line%direction(Y)**2
      delta = line%direction(X) * ptMinusExtents(X) &
           + line%direction(Y) * ptMinusExtents(Y) &
           + line%direction(Z) * ptPlusExtents(Z)
      t = - delta/ lSqr

      d2 = d2 + ptMinusExtents(X)**2 &
           + ptMinusExtents(Y) **2 &
           + ptPlusExtents(Z)**2 &
           + delta * t

    end subroutine region5code

  end subroutine LineBoxDistanceSquared

  subroutine PointBoxDistanceSquared(box, point, d2)
    type(Box_t), intent(in) :: box
    real(dp), intent(in) :: point(3)
    real(dp), intent(out) :: d2

    real(dp) :: offset(3), pPrime(3), d

    integer :: dim

    d2 = 0
    offset = point - box%origin
    pPrime = (/&
         dot_product(offset, box%u), &
         dot_product(offset, box%v), &
         dot_product(offset, box%w)  &
         /)

    do dim = 1, 3
       if (pPrime(dim) < -box%extents(dim)) then
          d = pPrime(dim) + box%extents(dim)
          d2 = d2 + d**2
       else if (pPrime(dim) > box%extents(dim)) then
          d = pPrime(dim) - box%extents(dim)
          d2 = d2 + d**2
       else
          ! nothing to do
       end if
    end do

  end subroutine PointBoxDistanceSquared

  subroutine PointLineDistanceSquared(line, point, t, d2)
    type(Line_t), intent(in) :: line
    real(dp), intent(in), dimension(3) :: point

    real(dp), intent(out) :: d2, t

    real(dp) :: qPrime(3), vec(3)

    t = dot_product(line%direction, point - line%origin) / dot_product(line%direction, line%direction)
    qPrime = line%origin + t * line%direction
    vec = point - qPrime
    d2 = dot_product(vec, vec)
  end subroutine PointLineDistanceSquared

  subroutine PointLineSegDistanceSquared(segment, point, t, d2)
    type(Segment_t), intent(in) :: segment
    real(dp), intent(in), dimension(3) :: point

    real(dp), intent(out) :: d2, t

    type(Line_t) :: line
    real(dp) :: vec(3)

    line%origin = segment%start
    line%direction = segment%end - segment%start

    call PointLineDistanceSquared(line, point, t, d2)
    if (t < 0) then
       t = 0
       vec = point - segment%start
       d2 = dot_product(vec, vec)
    else if (t > 1) then
       t = 1
       vec = point - segment%end
       d2 = dot_product(vec, vec)
    end if

  end subroutine PointLineSegDistanceSquared

end module geom_distances

