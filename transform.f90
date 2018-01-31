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
       free_trans, &
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
  
  subroutine equivalent_center(pos,cell,icell)
    
    double precision, intent(inout), dimension(:,:) :: &
         pos  ! position matrix

    double precision, intent(in), dimension(3,3) :: &
         cell, &  ! cell
         icell ! inverse of cell
    
    pos = pos - spread(matmul(cell,nint(matmul(icell,sum(pos,2)/size(pos,2)))),2,size(pos,2))
    
  end subroutine equivalent_center

  function cost(Apos,Bpos)

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! Rotation axis
    double precision, dimension(size(Apos,2),size(Bpos,2)) :: &
         cost
    integer :: &
         i,j ! Iterators
    
    cost = 0
    
    do i=1,size(Apos,2)
       do j=1,size(Bpos,2)
          cost(i,j) = norm(Apos(:,i)-Bpos(:,j))
       enddo
    enddo
    
  end function cost

   function cost_map(Apos,Bpos) result(cost)

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! Rotation axis
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
  
  function dist(Apos, Bpos, frac, atoms, n_atoms)

    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms
    double precision, intent(in) :: &
         frac
    double precision, allocatable, dimension(:,:) :: &
         dmat
    double precision :: &
         dist
    integer :: &
         An_cell, Bn_cell, & ! Number of cells
         i
    integer, parameter :: &
         type = 4 ! Type of distance. 1: Hausdorff, 2: Sum of distances, 3: semi-Hausdorff 4: Semi-Semi-Hausdorff(ssh) 

    if (type == 1) then
    
       dist = 0

       An_cell = size(Apos,2)/sum(atoms)
       Bn_cell = size(Bpos,2)/sum(atoms)
    
       do i=0,n_atoms-1

          allocate(dmat(An_cell*atoms(i+1),Bn_cell*atoms(i+1)))
          dmat = cost(Apos(:,An_cell*i+1:An_cell*(i+atoms(i+1))), &
               Bpos(:,Bn_cell*i+1:Bn_cell*(i+atoms(i+1))))
          dist = max(maxval(minval(dmat,2)), maxval(minval(dmat,1)),dist) 
          deallocate(dmat)
       
       enddo

    else if (type == 2) then

       dist = 0

       An_cell = size(Apos,2)/sum(atoms)
       Bn_cell = size(Bpos,2)/sum(atoms)
    
       do i=0,n_atoms-1

          allocate(dmat(An_cell*atoms(i+1),Bn_cell*atoms(i+1)))
          dmat = cost(Apos(:,An_cell*i+1:An_cell*(i+atoms(i+1))), &
               Bpos(:,Bn_cell*i+1:Bn_cell*(i+atoms(i+1))))
          dist = sum(dmat) + dist
          deallocate(dmat)
       
       enddo

       else if (type == 3) then
    
       dist = 0

       An_cell = size(Apos,2)/sum(atoms)
       Bn_cell = size(Bpos,2)/sum(atoms)
    
       do i=0,n_atoms-1

          allocate(dmat(An_cell*atoms(i+1),Bn_cell*atoms(i+1)))
          dmat = cost(Apos(:,An_cell*i+1:An_cell*(i+atoms(i+1))), &
               Bpos(:,Bn_cell*i+1:Bn_cell*(i+atoms(i+1))))
          dist = max(maxval(minval(dmat,1)),dist) 
          deallocate(dmat)
       enddo

       else if (type == 4) then

          dist = 0

          An_cell = size(Apos,2)/sum(atoms)
          Bn_cell = size(Bpos,2)/sum(atoms)

          do i=0,n_atoms-1

             allocate(dmat(An_cell*atoms(i+1),Bn_cell*atoms(i+1)))
             dmat = cost(Apos(:,An_cell*i+1:An_cell*(i+atoms(i+1))), &
                  Bpos(:,Bn_cell*i+1:Bn_cell*(i+atoms(i+1))))
             dist = max(maxval(minval(dmat(1:int(size(Apos,2)*frac),:),2)), maxval(minval(dmat(:,1:int(size(Bpos,2)*frac)),1)),dist)
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
         pi = 3.141592653589793d0, &
         frac = 1
    

    call init_random_seed()

    ! Center both cells at the geometric center
    call center(Bpos,n)
    call center(Apos,n)

    ! Initital distance as the reference
    dmin = dist(Apos,Bpos,frac, atoms, n_atoms)
    
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

       d = dist(Apos,newBpos,frac,atoms, n_atoms)

       ! Keep only better results no need for a temp for now
       if (d <= dmin) then
          dmin = d
          Bpos = newBpos
       endif

    enddo

    ! Finds the best mapping
    call munkres(dmin,map,cost(Apos,Bpos),n)

    dmin = dist(Apos,Bpos,frac,atoms, n_atoms)
    
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

    subroutine gradient_descent(mat, vec, Apos, Bpos, frac, atoms,n_atoms, n_iter, rate1, rate2)
    
    integer, intent(in) :: &
         n_iter, n_atoms ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         frac ! fraction of A and B to use for optimisation
             
    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Position of the atoms
    
    double precision, dimension(3,size(Bpos,2)) :: &
         pos, postmp ! position matrix
    
    double precision, intent(inout), dimension(3,1) :: &         
         vec     ! Translation vector

    double precision, intent(inout), dimension(3,3) :: &
         mat    ! Transformation matrix

    double precision, dimension(3,3) :: &
         mat_tmp, & ! Temporary transformation matrix
         mat_out
    
    double precision, dimension(3,1) :: &
         vec_tmp, & ! Temporary vector
         vec_out    ! Output vec

    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus   ! distance when substracting dx
    
    integer :: &
         i,j, & ! Iterator
         n

    double precision, parameter :: &
         dx = 1d-5  ! Derivative step

    n = size(Bpos,2)

    mat_out = 0
    vec_out = 0
    
    do j=1, n_iter

       ! mat
       do i=1,4

          mat_tmp = mat

          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) + dx
          dist_plus = dist(Apos, free_trans(Bpos, mat_tmp, vec),frac,atoms,n_atoms)

          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - dx
          dist_minus = dist(Apos, free_trans(Bpos, mat_tmp, vec),frac,atoms,n_atoms)

          mat_out((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - rate1*(dist_plus - dist_minus) / ( 2 * dx )

       enddo
       
       ! vec
       do i=1,2

          vec_tmp = vec

          vec_tmp(i,1) = vec(i,1) + dx
          dist_plus = dist(Apos,free_trans(Bpos, mat, vec_tmp),frac,atoms, n_atoms)

          vec_tmp(i,1) = vec(i,1) - dx
          dist_minus = dist(Apos,free_trans(Bpos, mat, vec_tmp),frac,atoms, n_atoms)

          vec_out(i,1) = vec(i,1) - rate2*(dist_plus - dist_minus) / ( 2 * dx )
          
       enddo

       vec = vec_out
       mat = mat_out
       
    enddo
    
  end subroutine gradient_descent
  
  subroutine gradient_descent_rand(mat, vec, Apos, Bpos, frac, atoms, n_atoms, n_iter, rate1, rate2, T)
    ! New Gradient Descent Random  
    
    integer, intent(in) :: &
         n_iter, n_atoms ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         T, & ! Monte-Carlo temperature
         frac ! Fraction of A and B to use in optimisation
         
    double precision :: &
         rand_rate1, & ! Rate for angles
         rand_rate2 ! Rate for disp
    
    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Centered position of the atoms
    
    double precision, dimension(3,size(Bpos,2)) :: &
         pos, postmp ! position matrix
    
    double precision, intent(inout), dimension(3,1) :: &
         vec     ! Translation vector

    double precision, intent(inout), dimension(3,3) :: &
         mat    ! Transformation matrix

    double precision, dimension(3,3) :: &
         mat_tmp, & ! Temporary transformation matrix
         mat_out, & ! Output transformation matrix
         mat_min
    
    double precision, dimension(3,1) :: &
         vec_tmp, & ! Temporary vector
         vec_out, & ! Output vec
         vec_min
    
    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus, & ! distance when substracting dx
         accept, & ! Accept step
         dist_cur, &
         dist_prev, &
         dist_min
    
    integer :: &
         i,j, & ! Iterator
         n ! Size of Bpos

    double precision, parameter :: &
         dx = 1d-5  ! Derivative step

    n = size(Bpos,2)
    
    call init_random_seed()

    dist_cur = dist(Apos,free_trans(Bpos, mat, vec),frac,atoms, n_atoms)
    
    dist_prev = dist_cur
    dist_min = dist_cur
    mat_min = mat
    vec_min = vec

    vec_out = 0
    mat_out = 0
    
    do j=1, n_iter

       call random_number(rand_rate1)

       rand_rate2 = 10**rand_rate1*rate2
       rand_rate1 = 10**rand_rate1*rate1

       ! mat
       do i=1,4

          mat_tmp = mat

          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) + dx
          dist_plus = dist(Apos,free_trans(Bpos, mat_tmp, vec),frac,atoms,n_atoms)
          
          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - dx
          dist_minus = dist(Apos,free_trans(Bpos, mat_tmp, vec) ,frac,atoms,n_atoms)

          mat_out((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - rand_rate1*(dist_plus - dist_minus) / ( 2 * dx )

       enddo

       
       ! vec
       do i=1,2

          vec_tmp = vec

          vec_tmp(i,1) = vec(i,1) + dx
          dist_plus = dist(Apos,free_trans(Bpos, mat, vec_tmp),frac,atoms, n_atoms)

          vec_tmp(i,1) = vec(i,1) - dx
          dist_minus = dist(Apos,free_trans(Bpos, mat, vec_tmp),frac,atoms, n_atoms)

          vec_out(i,1) = vec(i,1) - rand_rate2*(dist_plus - dist_minus) / ( 2 * dx )
          
       enddo
       
       dist_cur = dist(Apos,free_trans(Bpos, mat_out, vec_out),frac,atoms, n_atoms)

       
       call random_number(accept)

       if (dist_cur <= dist_prev .or. accept <= exp(-(dist_cur-dist_prev)/T)) then

          dist_prev = dist_cur
          
          mat = mat_out
          vec = vec_out

          if (dist_cur < dist_min) then
             dist_min = dist_cur
             mat_min = mat
             vec_min = vec
          endif
       else
          
       endif

    enddo

    mat = mat_min
    vec = vec_min
    
  end subroutine gradient_descent_rand


  subroutine fastmapping(map, dmin, Apos,na, Bpos, nb, frac, Acell, iAcell, atoms, n_atoms, n_iter, rate1, rate2, T)

    use hungarian
    
    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2, & ! Rate of the gradient descent for displacement
         T
         
    integer, intent(in) :: &
         n_iter   ! Number of iteration of the gradient descent
         
    double precision, intent(inout), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(inout), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, intent(in), dimension(3,3) :: &
         Acell, & ! Unit cell of A
         iAcell
    
    double precision, intent(in) :: &
         frac ! Fraction of A and B to use in optimisation
    
    double precision, allocatable, dimension(:,:) :: &
         inBpos ! Position of the atoms

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type
    
    integer, intent(out), dimension(na) :: &
         map ! List of index
    
    double precision, intent(out) :: &
         dmin
    
    double precision, dimension(3,1) :: &
         vec, & ! Translation vecto
         u      ! Rotation axis

    double precision :: &
         tetha, & ! Rotation angle
         d

    double precision, allocatable, dimension(:,:) :: &
         dmat

    double precision, dimension(na,nb) :: &
         mat

    double precision, dimension(3,3) :: &
         tmat ! Transformation matrix
    
    integer :: &
         i, &
         An_cell, Bn_cell

    real :: &
         start, &
         finish
    
    double precision, parameter :: &
         pi = 3.141592653589793d0

    call init_random_seed()

    ! Center both cells at the geometric center
    call center(Bpos,nb)
    call center(Apos,na)

    ! Random initial step
    call random_number(tmat)
    tmat(3,:) = 0 ! For 2D only
    tmat(:,3) = 0 ! For 2D only
        
    call gradient_descent_rand(tmat, vec, Apos, Bpos, &
         frac, atoms,n_atoms,n_iter, rate1, rate2, T)
        
    call gradient_descent(tmat, vec, Apos, Bpos, &
         frac, atoms,n_atoms,n_iter, rate1/100.0d0, rate2/100.0d0)

    Bpos = free_trans(Bpos, tmat, vec)
    
    call equivalent_center(Bpos,Acell,iAcell)
    
    ! Finds the best mapping
    An_cell = size(Apos,2)/sum(atoms)
    Bn_cell = size(Bpos,2)/sum(atoms)

    dmin = 0
    
    do i=0,n_atoms-1
       
       allocate(dmat(An_cell*atoms(i+1),An_cell*atoms(i+1)))

       call cpu_time(start)
       dmat = cost_map(Apos( : , An_cell*i+1 : An_cell*(i+atoms(i+1)) ), &
                       Bpos( : , Bn_cell*i+1 : Bn_cell*i + Bn_cell*atoms(i+1)*frac))
       call cpu_time(finish)

       print*, "dist_time", finish - start 
       
       call cpu_time(start)
       call munkres(d,map(An_cell*i+1:An_cell*(i+atoms(i+1))),dmat,An_cell*atoms(i+1))
       call cpu_time(finish)

       print*, "munkres_time", finish - start
       
       dmin = dmin + d

       deallocate(dmat)

       map(An_cell*i+1:An_cell*(i+atoms(i+1))) = map(An_cell*i+1:An_cell*(i+atoms(i+1))) + An_cell*i

    enddo

    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(*,"(10(F5.3,X))") mat   
    ! call munkres(dmin,map,cost(Apos,Bpos,n),n)


    dmin = dist(Apos, Bpos,frac,atoms, n_atoms)
    
  end subroutine fastmapping
  
end module transform
