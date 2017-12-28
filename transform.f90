module transform

  implicit none

  public :: &
       trans, &
       center, &
       mapping
       
  private :: &
       eye, &
       init_random_seed, &
       norm, &
       dist, &
       dist2

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

  subroutine trans(pos,n,tetha,u,vec)

    integer, intent(in) :: &
         n ! Number of atoms
    
    double precision, intent(inout), dimension(3,n) :: &
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

    call center(pos,n)
    
    pos = matmul(R,pos) + spread(tvec,2,size(pos,2))
    
  end subroutine trans

  subroutine center(pos,n)

    integer, intent(in) :: &
         n ! Number of atoms
    
    double precision, intent(inout), dimension(3,n) :: &
         pos  ! position matrix

    pos = pos - spread(sum(pos,2) / size(pos,2),2,size(pos,2))
    
  end subroutine center

  function cost(posA,posB,n)

    integer, intent(in) :: &
         n ! Number of atoms
    double precision, intent(inout), dimension(3,n) :: &
         posA, posB  ! Rotation axis
    double precision, dimension(n,n) :: &
         cost
    integer :: &
         i,j ! Iterators
    
    cost = 0
    
    do i=1,size(posA,2)
       do j=1,size(posA,2)
          cost(i,j) = norm(posA(:,i)-posB(:,j))
       enddo
    enddo
    
  end function cost
  
  function dist(posA,posB,n)

    integer, intent(in) :: &
         n ! Number of atoms
    double precision, intent(inout), dimension(3,n) :: &
         posA, posB  ! Rotation axis
    double precision, dimension(n,n) :: &
         dmat
    double precision :: &
         dist
    
    dmat = cost(posA, posB, n)
    
    dist = max(maxval(minval(dmat,2)), maxval(minval(dmat,1))) 
    
  end function dist
  
  function dist2(posA,posB,n)

    integer, intent(in) :: &
         n ! Number of atoms
    double precision, intent(inout), dimension(3,n) :: &
         posA, posB  ! Rotation axis
    double precision, dimension(n,n) :: &
         dmat
    double precision :: &
         dist2
    
    dmat = cost(posA, posB, n)
    
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

  subroutine mapping(Apos, Bpos, n, map, dmin)

    use hungarian
    
    integer, intent(in) :: &
         n ! Total number of atoms 

    double precision, intent(inout), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms

    integer, intent(out), dimension(n) :: &
         map ! List of index

    double precision, intent(out) :: &
         dmin
    
    double precision, dimension(3,n) :: &
         newBpos ! new positions of the atoms 

    double precision, dimension(3,1) :: &
         vec, & ! Translation vecto
         u      ! Rotation axis

    double precision :: &
         d, &   ! Distance between the two structures
         tetha   ! Rotation angle

    integer :: &
         i   ! Iterator

    integer, parameter :: &
         n_iter = 100000 ! Number of Monte-Carlo iterations

    double precision, parameter :: &
         pi = 3.141592653589793d0

    ! Center both cells at the geometric center
    call center(Bpos,n)
    call center(Apos,n)

    ! Initital distance as the reference
    dmin = dist(Apos,Bpos,n)
    
    ! Optimization iterations
    do i=1,n_iter

       ! Random translation, rotation axis and angle
       call random_number(tetha)
       call random_number(vec)
       call random_number(u)
       u = (u - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))) ! Recasts vector in the 8 octans
       vec = vec - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))
       u = u / norm(u) ! Normalizes
       vec = dmin * vec / sqrt(3.0d0) ! Moves less when close
       tetha = 2*pi*tetha

       newBpos = Bpos

       call trans(newBpos,n,tetha,u,vec) ! MC step

       d = dist(Apos,newBpos,n)

       ! Keep only better results no need for a temp for now
       if (d <= dmin) then
          dmin = d
          Bpos = newBpos
       endif

    enddo

    ! Finds the best mapping
    call munkres(dmin,map,cost(Apos,Bpos,n),n)
    
  end subroutine mapping
   

end module transform
