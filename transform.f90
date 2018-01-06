module transform

  implicit none

  public :: &
       trans, &
       center, &
       mapping, &
       fastmapping
       
  private :: &
       eye, &
       init_random_seed, &
       norm, &
       dist, &
       dist2, &
       gradient_descent, &
       u2angles, &
       angles2u

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
    double precision, intent(in), dimension(3,n) :: &
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
    double precision, intent(in), dimension(3,n) :: &
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
    double precision, intent(in), dimension(3,n) :: &
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

    call init_random_seed()

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

    dmin = dist(Apos,Bpos,n)
    
  end subroutine mapping

  ! subroutine old_fastmapping(Apos, Bpos, n, map, dmin, rate, n_iter)

  !   use hungarian

  !   integer, intent(in) :: &
  !        n, & ! Total number of atoms
  !        n_iter ! Number of iterations 

  !   double precision, intent(inout), dimension(3,n) :: &
  !        Apos, Bpos ! Position of the atoms

  !   double precision, intent(in) :: &
  !        rate

  !   integer, intent(out), dimension(n) :: &
  !        map ! List of index

  !   double precision, intent(out) :: &
  !        dmin

  !   double precision, dimension(3,n) :: &
  !        Y ! new positions of the atoms

  !   double precision, dimension(3,3) :: &
  !        W, & ! new positions of the atoms
  !        grad

  !   double precision :: &
  !        mse

  !   double precision, dimension(3) :: &
  !        tmp

  !   integer :: &
  !        i,j   ! Iterator

  !   double precision, parameter :: &
  !        pi = 3.141592653589793d0

  !   ! Center both cells at the geometric center
  !   call center(Bpos,n)
  !   call center(Apos,n)

  !   ! Initital distance as the reference
  !   dmin = dist(Apos,Bpos,n)

  !   call init_random_seed()
  !   call random_number(W)

  !   ! Optimization iterations
  !   do i=1,n_iter
  !      Y = matmul(W,Bpos)
  !      grad=0
  !      mse=0
  !      do j=1,n
  !         tmp = Y(:,n)
  !         Y(:,2:n) = Y(:,1:n-1)
  !         Y(:,1)= tmp
  !         mse = mse + sum(sum((Apos - Y)**2,1))
  !         grad = grad + matmul(Apos - Y, transpose(Bpos))
  !      enddo
  !      W = W + rate*grad
  !      print*, mse
  !   enddo

  !   Bpos = Y

    
  !   ! Finds the best mapping
  !   call munkres(dmin,map,cost(Apos,Bpos,n),n)

  ! end subroutine old_fastmapping

  subroutine fastmapping(Apos, Bpos, n, map, dmin, n_iter, rate1, rate2)

    use hungarian
    
    integer, intent(in) :: &
         n ! Total number of atoms

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2 ! Rate of the gradient descent for displacement

    integer, intent(in) :: &
         n_iter   ! Number of iteration of the gradient descent

    double precision, intent(inout), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms

    integer, intent(out), dimension(n) :: &
         map ! List of index

    double precision, intent(out) :: &
         dmin
    
    double precision, dimension(3,1) :: &
         vec, & ! Translation vecto
         u      ! Rotation axis

    double precision :: &
         tetha   ! Rotation angle

    double precision, parameter :: &
         pi = 3.141592653589793d0

    call init_random_seed()

    ! Center both cells at the geometric center
    call center(Bpos,n)
    call center(Apos,n)

    ! Random initial step
    call random_number(tetha)
    vec = 0 
    call random_number(u)
    u = (u - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))) ! Recasts vector in the 8 octans
    u = u / norm(u) ! Normalizes
    tetha = 2*pi*tetha
    
    call gradient_descent(tetha, u, vec, Apos, Bpos, n, n_iter, rate1, rate2)

    call trans(Bpos,n,tetha,u,vec)
    
    ! Finds the best mapping
    call munkres(dmin,map,cost(Apos,Bpos,n),n)

    dmin = dist(Apos,Bpos,n)
    
  end subroutine fastmapping

  subroutine gradient_descent(tetha, u, vec, Apos, Bpos, n, n_iter, rate1, rate2)
    
    integer, intent(in) :: &
         n, n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2 ! Rate for disp

    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms
    
    double precision, dimension(3,n) :: &
         pos  ! position matrix
    
    double precision, intent(inout), dimension(3,1) :: &
         u, &    ! Rotation axis (unitary vector)
         vec     ! Translation vector

    double precision, dimension(3,1) :: &
         vec_tmp, & ! Temporary vector
         vec_out  ! Output vec
    
    double precision, intent(inout) :: &
         tetha    ! angle of rotation

    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus, & ! distance when substracting dx
         tetha_out, & ! Output angle
         a1, & ! Initial value of phi in u
         a2, & ! Initial value of tetha in u
         a1_out, & ! First output angle in u
         a2_out ! Second output angle in u

    integer :: &
         i,j ! Iterator

    double precision, parameter :: &
         dx = 1d-5  ! Derivative step
    
    do j=1, n_iter

       ! Tetha
       pos = Bpos
       call trans(pos,n,tetha + dx,u,vec)
       dist_plus = dist(Apos,pos,n)

       pos = Bpos
       call trans(pos,n,tetha - dx,u,vec)
       dist_minus = dist(Apos,pos,n)

       tetha_out = tetha - rate1*(dist_plus - dist_minus) / ( 2 * dx )

       ! u vector
       call u2angles(a1,a2,u)

       !a1
       pos = Bpos
       call trans(pos,n,tetha,angles2u(a1+dx,a2),vec)
       dist_plus = dist(Apos,pos,n)

       pos = Bpos
       call trans(pos,n,tetha,angles2u(a1-dx,a2),vec)
       dist_minus = dist(Apos,pos,n)

       a1_out = a1 - rate1*(dist_plus - dist_minus) / ( 2 * dx )

       !a2
       pos = Bpos
       call trans(pos,n,tetha,angles2u(a1,a2+dx),vec)
       dist_plus = dist(Apos,pos,n)

       pos = Bpos
       call trans(pos,n,tetha,angles2u(a1,a2-dx),vec)
       dist_minus = dist(Apos,pos,n)

       a2_out = a2 - rate1*(dist_plus - dist_minus) / ( 2 * dx )
       
       ! vec
       do i=1,3

          vec_tmp = vec

          pos = Bpos
          vec_tmp(i,1) = vec(i,1) + dx
          call trans(pos,n,tetha,u,vec_tmp)
          dist_plus = dist(Apos,pos,n)

          pos = Bpos
          vec_tmp(i,1) = vec(i,1) - dx
          call trans(pos,n,tetha,u,vec_tmp)
          dist_minus = dist(Apos,pos,n)

          vec_out(i,1) = vec(i,1) - rate2*(dist_plus - dist_minus) / ( 2 * dx )
          
       enddo
       
       tetha = tetha_out
       u = angles2u(a1_out, a2_out)
       vec = vec_out

    enddo
    
  end subroutine gradient_descent

  subroutine u2angles(a1,a2,u)

    double precision, intent(in), dimension(3) :: &
         u

    double precision, intent(out) :: &
         a1, a2

    a1 = atan2(sqrt(u(1)**2+u(2)**2),u(3))
    a2 = atan2(u(2),u(1))    
    
  end subroutine u2angles

  function angles2u(a1,a2) result(u)

    double precision, dimension(3,1) :: &
         u

    double precision, intent(in) :: &
         a1, a2

    u(1,1) = sin(a1)*cos(a2)
    u(2,1) = sin(a1)*sin(a2)
    u(3,1) = cos(a1)
    
  end function angles2u
  
end module transform
