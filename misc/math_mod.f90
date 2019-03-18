module math_mod
  use amr_parameters, only : dp
  implicit none
  
contains
  elemental function safe_exp(x)
    real(dp), intent(in) :: x
    real(dp) :: safe_exp
    if (x < -300._dp) then
       safe_exp = 0._dp
    else
       safe_exp = exp(x)
    end if
  end function safe_exp
  
end module math_mod
