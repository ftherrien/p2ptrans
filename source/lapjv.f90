! Fortran interface for lapjv c++ code

module jv
 interface
  subroutine lapjv(sumSol, jSol, cost, n) bind(C, name="lapjv_fort")
    use iso_c_binding
    implicit none

    real(c_double), intent (out) :: sumSol
    integer(c_int), intent (out), dimension (*) :: jSol
    real(c_double), intent (in), dimension (*) :: cost
    integer(c_int), intent (in) :: n
  end subroutine
 end interface
end module jv