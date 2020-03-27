module utils

  implicit none

contains

  function free_trans(pos, mat, vec)

    double precision, intent(in), dimension(:,:) :: &
         pos  ! position matrix

    double precision, intent(in), dimension(3,3) :: &
         mat

    double precision, intent(in), dimension(3) :: &
         vec

    double precision, dimension(size(pos,1), size(pos,2)) :: &
         free_trans

    free_trans = matmul(mat, pos) + spread(vec, 2, size(pos,2))

  end function free_trans

  function rot_mat(angles) result(R)

    ! Creates the rotation matrix from the 3 angles

    double precision, intent(in), dimension(3) :: &
         angles     ! Rotation axis (unitary vector)

    double precision, dimension(3,1) :: &
         u     ! Rotation axis (unitary vector)

    double precision, dimension(3,3) :: &
         Q, & ! cross product matrix
         P    ! u.u**T

    double precision, dimension(3,3) :: &
         R    ! Transformation matrix

    u(1,1) = sin(angles(2)) * cos(angles(3))
    u(2,1) = sin(angles(2)) * sin(angles(3))
    u(3,1) = cos(angles(2))

    P = matmul(u,transpose(u))
    Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

    R = P + (eye() - P)*cos(angles(1)) + Q*sin(angles(1))

  end function rot_mat

  subroutine center(pos,n)

    ! Center the structures on their center of mass

    integer, intent(in) :: &
         n ! Number of atoms

    double precision, intent(inout), dimension(3,n) :: &
         pos  ! position matrix

    pos = pos - spread(sum(pos,2) / size(pos,2),2,size(pos,2))

  end subroutine center

  function eye() result(a)
    ! Creates the identity matrix
    ! Copied from Rosetta Code at: https://rosettacode.org/wiki/Identity_matrix#Fortran
    ! Checked and modified for double

    double precision :: a(3,3)
    integer :: i,j

    forall(i = 1:3, j = 1:3) a(i,j) = (i/j)*(j/i)

  end function eye

  function norm(a)

    ! Calculates the norm of a vector

    double precision :: norm
    double precision, dimension(3), intent(in) :: a

    norm = sqrt(sum(a**2))

  end function norm

  function split(a,char)

    ! Only keep what is before "char" 

    character*200 :: split
    character*1, intent(in) :: char
    character*200, intent(in) :: a
    integer :: i

    i=0
    do while (a(len(a) - i:len(a) - i) /= char .and. i < len(a))
       i = i + 1
    enddo

    split = trim(a(1:len(a)-i))

  end function split

  recursive function det(a,n) result(accumulation)
    ! Finds the determinant of an n by n matrix
    ! Copied from Rosetta Code at: https://rosettacode.org/wiki/Matrix_arithmetic#Fortran
    ! Checked and modified for determinant only, and double precision
    double precision, dimension(n,n), intent(in) :: a
    integer, intent(in) :: n
    double precision, dimension(n-1, n-1) :: b
    double precision :: accumulation
    integer :: i, sgn
    if (n == 1) then
       accumulation = a(1,1)
    else
       accumulation = 0
       sgn = 1
       do i=1, n
          b(:, :(i-1)) = a(2:, :i-1)
          b(:, i:) = a(2:, i+1:)
          accumulation = accumulation + sgn * a(1, i) * det(b, n-1)
          sgn = -sgn
       enddo
    endif
  end function det

  function sort(array) result(idx)

    double precision, intent(in), dimension(:) :: &
         array
    ! Main subroutine to exported in python. Minimizes the distance and resturns the optimal
    double precision, dimension(size(array)) :: &
         order
    ! result
    integer, dimension(size(array)) :: &
         idx

    integer :: &
         k,i,n

    n = size(array)

    order(1) = array(1)
    idx(1) = 1

    do i=2, n 
       k=0
       do while (array(i) < order(i-1-k) .and. k < i-1)
          k=k+1
       enddo
       order(i-k+1:n) = order(i-k:n-1)
       order(i-k) = array(i)
       idx(i-k+1:n) = idx(i-k:n-1)
       idx(i-k) = i
    enddo

  end function sort

end module utils

