module tiling

  implicit none

  public :: &
       parallelepiped, &
       sphere, &
       circle
       
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


  subroutine sphere(Acell, ncell, center, ASC)

    double precision, intent(in), dimension(3,3) :: &
         Acell ! Matrix containing the 3 cell vectors

    integer, intent(in) :: &
         ncell ! Number of repetition in each direction
    
    double precision, intent(in), dimension(3) :: &
         center

    double precision, intent(out), dimension(3,ncell) :: &
         ASC ! Position of each cell in the SC

    double precision, dimension(ncell) :: &
         dists
    
    integer :: &
         i,j,k,l,m,nnx,nny,nnz, &
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
    
    rad = (abs(det(Acell,3))*ncell*3/(4*pi))**(1.0d0/3.0d0) + diag

    nnz = abs(int(2*rad*norm(cross(Acell(:,1),Acell(:,2)))/(2*dot_product(Acell(:,3),cross(Acell(:,1),Acell(:,2))))))+1
    nny = abs(int(2*rad*norm(cross(Acell(:,3),Acell(:,1)))/(2*dot_product(Acell(:,2),cross(Acell(:,3),Acell(:,1))))))+1
    nnx = abs(int(2*rad*norm(cross(Acell(:,2),Acell(:,3)))/(2*dot_product(Acell(:,1),cross(Acell(:,2),Acell(:,3))))))+1

    ! write(*,*) 'Square cell dimensions:',nnx, nny, nnz

    l = 0
    dists = -1.0d0
    do i=-nnx,nnx
       do j=-nny,nny
          do k=-nnz,nnz
             dist = norm(matmul(Acell,(/i,j,k/)) + center)
             if (dist <= rad) then
                if (l < ncell) l=l+1
                if (l == 1) then
                   dists(1) = dist
                else if (dist < dists(l) .or. dists(l) == -1.0d0) then
                   m=0
                   do while (dist < dists(l-1-m) .and. m < l-1)
                      m=m+1
                   enddo
                   dists(l-m+1:ncell) = dists(l-m:ncell-1)
                   dists(l-m) = dist
                   ASC(:,l-m+1:ncell) = ASC(:,l-m:ncell-1)
                   ASC(:,l-m) = matmul(Acell,(/i,j,k/)) + center
                endif
             endif
          enddo
       enddo
    enddo

  end subroutine sphere  

  subroutine circle(Acell,ncell,ASC)

    double precision, intent(in), dimension(2,2) :: &
         Acell ! Matrix containing the 2 cell vectors

    integer, intent(in) :: &
         ncell ! Number of repetition in each direction

    double precision, intent(out), dimension(2,ncell) :: &
         ASC ! Position of each cell in the SC

    double precision, dimension(ncell) :: &
         dists
    
    integer :: &
         i,j,k,l,nnx,nny, &
         pos

    double precision :: &
         rad, dist, &
         diag, area

    double precision, parameter :: &
         pi = 3.141592653589793d+0

    ! Diagonals (This is probably overkill)
    ! Finds the longest diagonal to add it to the radius of the circle such that there is
    ! at least n intger number of cells in the sphere.
    diag = 0
    diag = max(norm(Acell(:,1) + Acell(:,2)),diag)
    diag = max(norm(-Acell(:,1) + Acell(:,2)),diag)

    area = sqrt((Acell(1,1)*Acell(2,2))**2 + (Acell(2,1)*Acell(1,2))**2)
    
    rad = sqrt(area*ncell/pi) + diag

    nnx = abs(int(rad*norm(Acell(:,2))/area))+1
    nny = abs(int(rad*norm(Acell(:,1))/area))+1

    ! print*, 'Square cell dimensions:',nnx ,nny

    l = 0
    dists = -1
    do i=-nnx,nnx
       do j=-nny,nny
          dist = norm(matmul(Acell,(/i+0.5d0,j+0.5d0/)))
          if (dist <= rad) then
             
             if (l < ncell) l=l+1
             if (l == 1) then
                dists(1) = dist
             else if (dist < dists(l) .or. dists(l) == -1) then
                k=0
                do while (dist < dists(l-1-k) .and. k < l-1)
                   k=k+1
                enddo
                dists(l-k+1:ncell) = dists(l-k:ncell-1)
                dists(l-k) = dist
                ASC(:,l-k+1:ncell) = ASC(:,l-k:ncell-1)
                ASC(:,l-k) = matmul(Acell,(/i,j/))
             endif
          endif
       enddo
    enddo

  end subroutine circle
  
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
    double precision, dimension(:), intent(in) :: a
    
    norm = sqrt(sum(a**2))
    
  end function norm

end module tiling

!Made by Louis Popovic. Special thanks to Felix Therrien. 
