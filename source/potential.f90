! README: This file contains two body potentials to minimize. If you wish to make your own potential,
! add a case to the `distance` and `derivative` functions.

module potential

  use utils

  implicit none

contains

  function distance(Apos, Bpos, mat, vec, pot, param)

    double precision, intent(in), dimension(:,:) :: &
         Apos, &  ! position matrix
         Bpos

    double precision, dimension(size(Apos,2)) :: &
         d

    double precision, intent(in), dimension(3,3) :: &
         mat

    double precision, intent(in), dimension(3) :: &
         vec

    double precision :: &
         distance

    double precision, intent(in) :: &
         param ! Parameter of the potential

    character*(*), intent(in) :: &
         pot

    select case (pot)
    case ("LJ")

       d = sum((Apos - free_trans(Bpos,mat,vec))**2,1)

       distance = sum(param**12/d**6 - 2*param**6/d**3)
       
    case ("Euclidean")
   
       distance = sum(sqrt(sum((Apos - free_trans(Bpos,mat,vec))**2,1)))
    end select

  end function distance

  function derivative(Apos, Bpos, mat, vec, pot, param) result(E)

    double precision, intent(in), dimension(:,:) :: &
         Apos, &  ! position matrix
         Bpos

    double precision, dimension(size(Apos,2)) :: &
         P

    double precision, dimension(3, size(Apos,2)) :: &
         E

    double precision, intent(in), dimension(3,3) :: &
         mat

    double precision, intent(in), dimension(3) :: &
         vec

    character*(*), intent(in) :: &
         pot

    double precision, intent(in) :: &
         param ! Parameter of the potential

    select case (pot)
    case ("LJ")
       E = Apos - free_trans(Bpos,mat,vec)
       P = param**6/(sum(E**2,1))**7 - 1/(sum(E**2,1))**4
       E = -E*spread(P,1,3)
    case ("Euclidean")
       E = Apos - free_trans(Bpos,mat,vec)
    end select

  end function derivative


end module potential
