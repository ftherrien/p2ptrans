module tiling

  use utils
  
  implicit none

  public :: &
       sphere, &
       circle
       
  private :: &
       cross
  
  contains
    
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
         i,j,k,l,m,nnx,nny,nnz

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
                else if (dist < dists(l) .or. dists(l) < 0.0d0) then
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

  subroutine circle(Acell, ncell, center, ASC)

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
         i,j,k,l,nnx,nny

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
          dist = norm(matmul(Acell,(/i+0.5d0,j+0.5d0,0.0d0/)) + center)
          if (dist <= rad) then
             
             if (l < ncell) l=l+1
             if (l == 1) then
                dists(1) = dist
             else if (dist < dists(l) .or. dists(l) < 0.0d0) then
                k=0
                do while (dist < dists(l-1-k) .and. k < l-1)
                   k=k+1
                enddo
                dists(l-k+1:ncell) = dists(l-k:ncell-1)
                dists(l-k) = dist
                ASC(:,l-k+1:ncell) = ASC(:,l-k:ncell-1)
                ASC(:,l-k) = matmul(Acell,(/dble(i),dble(j),0.0d0/)) + center
             endif
          endif
       enddo
    enddo

  end subroutine circle
  
  function cross(a, b)
    !Copied from Rosetta Code at: https://rosettacode.org/wiki/Vector_products#Fortran
    ! Checked and modified for double precision
    double precision, dimension(3) :: cross
    double precision, dimension(3), intent(in) :: a, b

    cross(1) = a(2)*b(3) - a(3)*b(2)
    cross(2) = a(3)*b(1) - a(1)*b(3)
    cross(3) = a(1)*b(2) - b(1)*a(2)
  end function cross
  
end module tiling

!Made by Louis Popovic. Special thanks to Felix Therrien. 
