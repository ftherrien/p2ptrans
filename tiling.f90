module tiling

  implicit none

  public :: &
       parallelepiped, &
       sphere
       
  private :: &
       calc_angles, &
       cross, &
       det, &
       norm

  contains

  subroutine parallelepiped(Acell,nx,ny,nz,ASC)

    double precision, intent(in), dimension(3,3) :: &
         Acell ! Matrix containing the 3 cell vectors

    integer, intent(in) :: &
         nx,ny,nz ! Number of repetition in each direction

    double precision, intent(out), dimension(3,nx*ny*nz) :: &
         ASC ! Position of each cell in the SC
    
    integer :: &
         i,j,k,l
       
    l = 0
    do i=1,nx
       do j=1,ny
          do k=1,nz
             l=l+1
             ASC(:,l) = matmul(Acell,(/i,j,k/))
          enddo
       enddo
    enddo
    
  end subroutine parallelepiped


  subroutine sphere(Acell,ncell,ASC)

    double precision, intent(in), dimension(3,3) :: &
         Acell ! Matrix containing the 3 cell vectors

    integer, intent(in) :: &
         ncell ! Number of repetition in each direction

    double precision, intent(out), dimension(3,ncell) :: &
         ASC ! Position of each cell in the SC

    double precision, dimension(ncell) :: &
         dists
    
    integer :: &
         i,j,k,l,nnx,nny,nnz, &
         pos

    double precision :: &
         rad, dist, &
         diag

    double precision, parameter :: &
         pi = 3.141592653589793d+0

    ! Diagonals (This is probably overkill)
    ! Finds the longest diagonal to add it to the radius of the sphere such that there is
    ! at least n intger number of cells in the sphere.
    diag = 0
    diag = max(norm(Acell(:,1) + Acell(:,2) + Acell(:,3)),diag)
    diag = max(norm(-Acell(:,1) + Acell(:,2) + Acell(:,3)),diag)
    diag = max(norm(Acell(:,1) - Acell(:,2) + Acell(:,3)),diag)
    diag = max(norm(-Acell(:,1) - Acell(:,2) + Acell(:,3)),diag)
    
    
    rad = (det(Acell,3)*ncell*3/(4*pi))**(1.0d0/3.0d0) + diag

    nnx = int(2*rad*norm(cross(Acell(:,1),Acell(:,2)))/(2*dot_product(Acell(:,3),cross(Acell(:,1),Acell(:,2)))))+1
    nny = int(2*rad*norm(cross(Acell(:,3),Acell(:,1)))/(2*dot_product(Acell(:,2),cross(Acell(:,3),Acell(:,1)))))+1
    nnz = int(2*rad*norm(cross(Acell(:,2),Acell(:,3)))/(2*dot_product(Acell(:,1),cross(Acell(:,2),Acell(:,3)))))+1

    print*, 'Test',nnx, nny, nnz, rad, ncell

    l = 0
    dists = 0
    do i=-nnx,nnx
       do j=-nny,nny
          do k=-nnz,nnz
             dist = norm(matmul(Acell,(/i+0.5d0,j+0.5d0,k+0.5d0/)))
             if (dist <= rad) then
                l=l+1
                if (l<=ncell) then
                   dists(l) = dist
                   ASC(:,l) = matmul(Acell,(/i,j,k/))
                elseif (maxval(dists) > dist) then
                   pos = maxloc(dists,1)
                   dists(pos) = dist
                   ASC(:,pos) = matmul(Acell,(/i,j,k/))
                endif
             endif
          enddo
       enddo
    enddo

  end subroutine sphere  

  subroutine calc_angles(A,lengths,angles)

    double precision, parameter :: &
         pi = 3.141592653589793d+0

    double precision, intent(in), dimension(3,3) :: &
         A         ! Cell structure

    double precision, intent(out), dimension(3) :: &
         lengths, &     ! Length vector 
         angles      ! Angle vector

    double precision, dimension(3,3) :: &
         m

    integer :: &
         i,j,k 

    m = transpose(A)
    lengths = sqrt(sum(m**2,2))
    angles = 0
    do i=1,3
       j = modulo(i, 3)+1
       k = modulo((i + 1), 3)+1
       angles(i) = max(-1.0d+0,min(1.0d+0,(dot_product(m(j,:), m(k,:)) / (lengths(j) * lengths(k)))))
    enddo
    
    angles = acos(angles) * 180.d+0 / pi

  end subroutine calc_angles

    function cross(a, b)
    !Copied from Rosetta Code at: https://rosettacode.org/wiki/Vector_products#Fortran
    ! Checked and modified for double precision
    double precision, dimension(3) :: cross
    double precision, dimension(3), intent(in) :: a, b

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - b(1)*a(2)
  end function cross

  recursive function det(a,n) result(accumulation)
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

  function norm(a)

    double precision :: norm
    double precision, dimension(3), intent(in) :: a
    
    norm = sqrt(sum(a**2))
    
  end function norm

end module tiling

!Made by Louis Popovic. Special thanks to Felix Therrien. 
