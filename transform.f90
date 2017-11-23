module transform

  implicit none

  public :: &
       init_random_seed, &
       norm, &
       trans, &
       center, &
       dist

  private :: &
       eye

  contains

  subroutine init_random_seed()
  ! Copied form the GCC docs: https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
            INTEGER :: i, n, clock
            INTEGER, DIMENSION(:), ALLOCATABLE :: seed
          
            CALL RANDOM_SEED(size = n)
            ALLOCATE(seed(n))
          
            CALL SYSTEM_CLOCK(COUNT=clock)
          
            seed = clock + 37 * (/ (i - 1, i = 1, n) /)
            CALL RANDOM_SEED(PUT = seed)
          
            DEALLOCATE(seed)

  end subroutine init_random_seed

  subroutine trans(pos,tetha,u,vec)

    double precision, intent(inout), dimension(:,:), allocatable :: &
         pos  ! position matrix

    double precision, intent(in), dimension(3,1) :: &
         u, &    ! Rotation axis (unitary vector)
         vec     ! Translation vector

    double precision, dimension(3) :: &
         tvec    ! displacement from origin
    
    double precision, intent(in) :: &
         tetha    ! angle of rotation

    double precision, dimension(3,3) :: &
         Q, & ! cross product matrix
         P, & ! u.u**T
         R    ! Transformation matrix

    P = matmul(u,transpose(u))
    Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

    R = P + (eye() - P)*cos(tetha) + Q*sin(tetha)

    tvec = vec(:,1) + sum(pos,2) / size(pos,2)

    call center(pos)
    
    pos = matmul(R,pos) + spread(tvec,2,size(pos,2))
    
  end subroutine trans

  subroutine center(pos)

    double precision, intent(inout), dimension(:,:), allocatable :: &
         pos  ! position matrix

    pos = pos - spread(sum(pos,2) / size(pos,2),2,size(pos,2))
    
  end subroutine center

  function dist(posA,posB)

    double precision, intent(inout), dimension(:,:), allocatable :: &
         posA, posB  ! Rotation axis
    double precision, dimension(:,:), allocatable :: &
         dmat
    double precision :: &
         dist
    integer :: &
         i,j ! Iterators

    allocate(dmat(size(posA,2),size(posA,2)))
    
    dmat = 0
    
    do i=1, size(posA,2)
       do j=1, size(posA,2)
          dmat(i,j) = norm(posA(:,i)-posB(:,j))
       enddo
    enddo

    ! write(*,'(200(F10.5,F10.5,F10.5,F10.5,F10.5,/))') dmat

    dist = max(maxval(minval(dmat,2)), maxval(minval(dmat,1))) 
    
  end function dist

  function dist2(posA,posB)

    double precision, intent(inout), dimension(:,:), allocatable :: &
         posA, posB  ! Rotation axis
    double precision, dimension(:,:), allocatable :: &
         dmat
    double precision :: &
         dist2
    integer :: &
         i,j ! Iterators
    
    allocate(dmat(size(posA,2),size(posA,2)))
    
    dmat = 0
    
    do i=1,size(posA,2)
       do j=1,i-1
          dmat(i,j) = norm(posA(:,i)-posB(:,j))
       enddo
    enddo

    dist2 = sum(dmat) 
    
  end function dist2

  function eye() result(a)
    !Copied from Rosetta Code at: https://rosettacode.org/wiki/Identity_matrix#Fortran
    ! Checked and modified for double
    double precision :: a(3,3)
    integer :: i,j

    forall(i = 1:3, j = 1:3) a(i,j) = (i/j)*(j/i)

  end function eye

  function norm(a)

    double precision :: norm
    double precision, dimension(3), intent(in) :: a

    norm = sqrt(sum(a**2))

  end function norm

end module transform
