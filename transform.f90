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
       cost, &
       cost_map, &
       gradient_descent_rand, &
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

  function cost(Apos,Bpos,n)

    integer, intent(in) :: &
         n ! Number of atoms
    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos  ! Rotation axis
    double precision, dimension(n,n) :: &
         cost
    integer :: &
         i,j ! Iterators
    
    cost = 0
    
    do i=1,size(Apos,2)
       do j=1,size(Apos,2)
          cost(i,j) = exp(-(norm(Apos(:,i))/5)**2)*norm(Apos(:,i)-Bpos(:,j)) ! PARAM
       enddo
    enddo
    
  end function cost

   function cost_map(Apos,Bpos) result(cost)

    double precision, intent(in), dimension(:,:) :: &
         Apos  ! Rotation axis
    double precision, intent(in), dimension(:,:) :: &
         Bpos  ! Rotation axis
    double precision, dimension(size(Apos,2),size(Apos,2)) :: &
         cost
    integer :: &
         i,j ! Iterators
    
    cost = 0
    
    do i=1,size(Apos,2)
       do j=1,size(Bpos,2)
          cost(i,j) = norm(Apos(:,i)-Bpos(:,j))
       enddo
    enddo
    
  end function cost_map
  
  function dist(Apos, Bpos, n, atoms, n_atoms)

    integer, intent(in) :: &
         n, & ! Number of atoms
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms
    double precision, allocatable, dimension(:,:) :: &
         dmat
    double precision :: &
         dist
    integer :: &
         n_cell, & ! Number of cells
         i
    integer, parameter :: &
         type = 1 ! Type of distance. 1: Hausdorff, 2. Sum of distances 
         

    if (type == 1) then
    
       dist = 0

       n_cell = n/sum(atoms)
    
       do i=0,n_atoms-1

          allocate(dmat(n_cell*atoms(i+1),n_cell*atoms(i+1)))
          dmat = cost(Apos(:,n_cell*i+1:n_cell*(i+atoms(i+1))), &
               Bpos(:,n_cell*i+1:n_cell*(i+atoms(i+1))), &
               n_cell*atoms(i+1))
          dist = max(maxval(minval(dmat,2)), maxval(minval(dmat,1)),dist) 
          deallocate(dmat)
       
       enddo

    else if (type == 2) then

       dist = 0

       n_cell = n/sum(atoms)
    
       do i=0,n_atoms-1

          allocate(dmat(n_cell*atoms(i+1),n_cell*atoms(i+1)))
          dmat = cost(Apos(:,n_cell*i+1:n_cell*(i+atoms(i+1))), &
               Bpos(:,n_cell*i+1:n_cell*(i+atoms(i+1))), &
               n_cell*atoms(i+1))
          dist = sum(dmat) + dist
          deallocate(dmat)
       
       enddo

    endif
    
  end function dist
  
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

  subroutine mapping(map, dmin, Apos, Bpos, n, atoms, n_atoms)

    use hungarian
    
    integer, intent(in) :: &
         n, & ! Total number of atoms
         n_atoms

    integer, dimension(n_atoms), intent(in) :: &
         atoms

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
    dmin = dist(Apos,Bpos,n, atoms, n_atoms)
    
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

       d = dist(Apos,newBpos,n,atoms, n_atoms)

       ! Keep only better results no need for a temp for now
       if (d <= dmin) then
          dmin = d
          Bpos = newBpos
       endif

    enddo

    ! Finds the best mapping
    call munkres(dmin,map,cost(Apos,Bpos,n),n)

    dmin = dist(Apos,Bpos,n,atoms, n_atoms)
    
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

    subroutine gradient_descent(tetha, u, vec, Apos, Bpos, n,atoms,n_atoms, n_iter, rate1, rate2)
    
    integer, intent(in) :: &
         n, n_iter, n_atoms ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2 ! Rate for disp
             
    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms
    
    double precision, dimension(3,n) :: &
         pos, postmp ! position matrix
    
    double precision, intent(inout), dimension(3,1) :: &
         u, &    ! Rotation axis (unitary vector)
         vec     ! Translation vector

    double precision, dimension(3,1) :: &
         vec_tmp, & ! Temporary vector
         vec_out, & ! Output vec
         u_min

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
       dist_plus = dist(Apos,pos,n,atoms, n_atoms)
       
       pos = Bpos
       call trans(pos,n,tetha - dx,u,vec)
       dist_minus = dist(Apos,pos,n,atoms, n_atoms)
       
       tetha_out = tetha - rate1*(dist_plus - dist_minus) / ( 2 * dx )

       
       ! vec
       do i=1,2

          vec_tmp = vec

          pos = Bpos
          vec_tmp(i,1) = vec(i,1) + dx
          call trans(pos,n,tetha,u,vec_tmp)
          dist_plus = dist(Apos,pos,n,atoms, n_atoms)

          pos = Bpos
          vec_tmp(i,1) = vec(i,1) - dx
          call trans(pos,n,tetha,u,vec_tmp)
          dist_minus = dist(Apos,pos,n,atoms, n_atoms)

          vec_out(i,1) = vec(i,1) - rate2*(dist_plus - dist_minus) / ( 2 * dx )
                    
       enddo
       
       tetha = tetha_out
       vec = vec_out

    enddo
    
  end subroutine gradient_descent
  
  subroutine gradient_descent_rand(tetha, u, vec, Apos, Bpos, n,atoms,n_atoms, n_iter, rate1, rate2, T)
    
    integer, intent(in) :: &
         n, n_iter, n_atoms ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         T
         
    double precision :: &
         rand_rate1, & ! Rate for angles
         rand_rate2 ! Rate for disp
    
    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms
    
    double precision, dimension(3,n) :: &
         pos, postmp ! position matrix
    
    double precision, intent(inout), dimension(3,1) :: &
         u, &    ! Rotation axis (unitary vector)
         vec     ! Translation vector

    double precision, dimension(3,1) :: &
         vec_tmp, & ! Temporary vector
         vec_out, & ! Output vec
         vec_min, &
         u_min
    
    double precision, intent(inout) :: &
         tetha    ! angle of rotation

    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus, & ! distance when substracting dx
         tetha_out, & ! Output angle
         a1, & ! Initial value of phi in u
         a2, & ! Initial value of tetha in u
         a1_out, & ! First output angle in u
         a2_out, & ! Second output angle in u
         accept, & ! Accept step
         dist_cur, &
         dist_prev, &
         dist_min, &
         tetha_min
    
    integer :: &
         i,j ! Iterator

    double precision, parameter :: &
         dx = 1d-5  ! Derivative step

    call init_random_seed()
    
    do j=1, n_iter

       call random_number(rand_rate1)

       rand_rate2 = 10**rand_rate1*rate2
       rand_rate1 = 10**rand_rate1*rate1

       ! Tetha
       pos = Bpos
       call trans(pos,n,tetha + dx,u,vec)
       dist_plus = dist(Apos,pos,n,atoms, n_atoms)
       
       pos = Bpos
       call trans(pos,n,tetha - dx,u,vec)
       dist_minus = dist(Apos,pos,n,atoms, n_atoms)
       
       tetha_out = tetha - rand_rate1*(dist_plus - dist_minus) / ( 2 * dx )

       
       ! vec
       do i=1,2

          vec_tmp = vec

          pos = Bpos
          vec_tmp(i,1) = vec(i,1) + dx
          call trans(pos,n,tetha,u,vec_tmp)
          dist_plus = dist(Apos,pos,n,atoms, n_atoms)

          pos = Bpos
          vec_tmp(i,1) = vec(i,1) - dx
          call trans(pos,n,tetha,u,vec_tmp)
          dist_minus = dist(Apos,pos,n,atoms, n_atoms)

          vec_out(i,1) = vec(i,1) - rand_rate2*(dist_plus - dist_minus) / ( 2 * dx )

          
       enddo

       pos = Bpos
       call trans(pos,n,tetha_out,u,vec_out)
       dist_cur = dist(Apos,pos,n,atoms, n_atoms)

       if (j==1) then
          dist_prev = dist_cur
          dist_min = dist_cur
          tetha_min = tetha
          vec_min = vec
       endif
       
       call random_number(accept)

       if (dist_cur <= dist_prev .or. accept <= exp(-(dist_cur-dist_prev)/T)) then

          dist_prev = dist_cur
          
          tetha = tetha_out
          vec = vec_out

          if (dist_cur < dist_min) then
             dist_min = dist_cur
             tetha_min = tetha
             vec_min = vec
          endif
       else
          !print*,"SKIP"
          
       endif

    enddo

    tetha = tetha_min
    vec = vec_min
    
  end subroutine gradient_descent_rand

  subroutine fastmapping(map, nmap, dmin, Apos, Bpos, n, atoms, n_atoms, n_iter, rate1, rate2, T)

    use hungarian
    
    integer, intent(in) :: &
         n, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2, & ! Rate of the gradient descent for displacement
         T
         
    integer, intent(in) :: &
         n_iter   ! Number of iteration of the gradient descent
         
    double precision, intent(inout), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms

    double precision, allocatable, dimension(:,:) :: &
         inBpos ! Position of the atoms

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type
    
    integer, intent(out), dimension(n) :: &
         map ! List of index
    
    double precision, intent(out) :: &
         dmin

    integer, intent(out), dimension(n_atoms) :: &
         nmap
    
    double precision, dimension(3,1) :: &
         vec, & ! Translation vecto
         u      ! Rotation axis

    double precision :: &
         tetha, & ! Rotation angle
         d

    double precision, allocatable, dimension(:,:) :: &
         dmat

    double precision, dimension(n,n) :: &
         mat
    
    integer :: &
         i, &
         n_cell
    
    double precision, parameter :: &
         pi = 3.141592653589793d0

    call init_random_seed()

    ! Center both cells at the geometric center
    call center(Bpos,n)
    call center(Apos,n)

    ! Random initial step
    call random_number(tetha)
    vec = 0 
    u = reshape((/0,0,1/),(/3,1/))
    tetha = 2*pi*tetha
    
    call gradient_descent_rand(tetha, u, vec, Apos, Bpos, n,atoms,n_atoms,n_iter, rate1, rate2, T)
    
    call gradient_descent(tetha, u, vec, Apos, Bpos, n,atoms,n_atoms,n_iter, rate1/100.0d0, rate2/100.0d0)

    call trans(Bpos,n,tetha,u,vec)
    
    ! Finds the best mapping
    n_cell = n/sum(atoms)

    dmin = 0
    
    do i=0,n_atoms-1

       nmap(i+1) = n_cell*atoms(i+1)*0.5d0
       
       allocate(dmat(n_cell*atoms(i+1),n_cell*atoms(i+1)))
       
       
       dmat = cost_map(Apos( : , n_cell*i+1 : n_cell*(i+atoms(i+1)) ), &
                       Bpos( : , n_cell*i+1 : n_cell*i+nmap(i+1)))
       
       call munkres(d,map(n_cell*i+1:n_cell*(i+atoms(i+1))),dmat,n_cell*atoms(i+1))

       dmin = dmin + d

       deallocate(dmat)

       map(n_cell*i+1:n_cell*(i+atoms(i+1))) = map(n_cell*i+1:n_cell*(i+atoms(i+1))) + n_cell*i

    enddo

    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(*,"(10(F5.3,X))") mat   
    ! call munkres(dmin,map,cost(Apos,Bpos,n),n)

    dmin = dist(Apos, Bpos, n, atoms, n_atoms)
    
  end subroutine fastmapping
  
end module transform
