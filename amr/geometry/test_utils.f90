module test_utils
  use, intrinsic :: iso_fortran_env, only : stdin=>input_unit, &
       stdout=>output_unit, &
       stderr=>error_unit
  use geom_types
  implicit none
  integer :: test_depth = 0, NTEST_OK, NTEST
contains
  subroutine test_group(name)
    character(len=*), intent(in) :: name
    integer :: i
    do i = 1, test_depth
       write(stdout, '(4x)', advance='no')
    end do
    write(stdout, *) name
    test_depth = test_depth + 1
  end subroutine test_group

  subroutine test_group_close
    test_depth = test_depth - 1
  end subroutine test_group_close


  subroutine assert(bool, msg, silent)
    logical, intent(in) :: bool
    character(len=*), intent(in) :: msg
    logical, intent(in), optional :: silent

    integer, save :: itest = 1

    integer :: i
    logical :: print_msg

    NTEST = NTEST + 1
    if (bool) NTEST_OK = NTEST_OK + 1
    if (present(silent)) then
       print_msg = .not. silent
    else
       print_msg = .true.
    end if

    do i = 1, test_depth
       write(stdout, '(4x)', advance='no')
    end do

    if (.not. bool) then
       write(stdout, '(4x,i3,1x,a)') itest, msg // '...failed!'
    else
       if (print_msg) &
            write(stdout, '(4x,i3,1x,a)') itest, msg // '...ok!'
    end if
    itest = itest + 1
  end subroutine assert

  function is_close(A, B, aerr)
    real(dp), intent(in) :: A, B, aerr
    logical :: is_close

    if (abs(A - B) < aerr) then
       is_close = .true.
    else
       is_close = .false.
    end if
  end function is_close

  subroutine yolo
    write(stderr, *) 'An error occured. Crashing.'
    stop 1
  end subroutine yolo

  subroutine test_summary
    write(*, '(i3,"/",i3," tests")') NTEST_OK, NTEST
    if (NTEST /= NTEST_OK) then
       ! Stop with an error code
       stop 1
    end if
  end subroutine test_summary

end module test_utils
