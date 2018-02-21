module transform

  use hungarian
  
  implicit none

  public :: &
       trans, &
       center, &
       mapping, &
       fastmapping, &
       test_ds
       
  private :: &
       eye, &
       init_random_seed, &
       norm, &
       dist, &
       dist_separate, &
       det, &
       free_trans, &
       cost, &
       cost_map, &
       analytical_gd_rot, &
       gradient_descent_rand, &
       gradient_descent_explore, &
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

  subroutine trans(pos,n,theta,u,vec)

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
         theta    ! angle of rotation

    double precision, dimension(3,3) :: &
         Q, & ! cross product matrix
         P, & ! u.u**T
         R    ! Transformation matrix

    P = matmul(u,transpose(u))
    Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

    R = P + (eye() - P)*cos(theta) + Q*sin(theta)

    tvec = vec(:,1) !+ sum(pos,2) / size(pos,2) TMP
    
    !call center(pos,n) TMP
    
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

    print*, "Center:",  matmul(cell,nint(matmul(icell,sum(pos,2)/size(pos,2))))
    
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

   function cost_map(Apos,Bpos,frac_in) result(cost)

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! Rotation axis
    double precision, optional, intent(in) :: &
         frac_in
    double precision :: &
         frac
    double precision, dimension(size(Apos,2),size(Apos,2)) :: &
         cost
    integer :: &
         i,j ! Iterators
    
    cost = 0

    frac = 1
    
    if (present(frac_in)) frac = frac_in
    
    do i=1,size(Apos,2)
       do j=1,size(Bpos,2)
          if (i <= int(frac*size(Apos,2)) .or. j <= int(frac*size(Bpos,2))) then 
             cost(i,j) = norm(Apos(:,i)-Bpos(:,j))
          endif
       enddo
    enddo
    
  end function cost_map
  
  function dist(Apos, Bpos, frac, atoms, n_atoms, tmat, vec)

    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! position matrix
    double precision, dimension(size(Bpos,1),size(Bpos,2)) :: &
         tBpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms
    double precision, optional, intent(in), dimension(3,3) :: &
         tmat
    double precision, optional, intent(in), dimension(3) :: &
         vec
    double precision, dimension(3,3) :: &
         rot_mat
    double precision, intent(in) :: &
         frac
    double precision, allocatable, dimension(:,:) :: &
         dmat
    double precision, dimension(3) :: &
         u_per
    integer, dimension(size(Apos,2)) :: &
         map ! List of index
    double precision :: &
         dist, &
         d
    integer :: &
         An_cell, Bn_cell, & ! Number of cells
         i
    integer :: &
         type_dist ! Type of distance. 1: Hausdorff, 2: Sum of distances, 3: semi-Hausdorff 4: Semi-Semi-Hausdorff(ssh) 5: Mapping distance 6. Mapping + stretching

    if (present(tmat) .and. present(vec)) then
       type_dist = 6
    else
       type_dist = 5
    endif
    
    if (type_dist == 1) then
    
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

    else if (type_dist == 2) then

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

       else if (type_dist == 3) then
    
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

       else if (type_dist == 4) then

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

       else if (type_dist == 5) then
          
          ! Finds the best mapping
          An_cell = size(Apos,2)/sum(atoms)
          Bn_cell = size(Bpos,2)/sum(atoms)

          dist = 0

          do i=0,n_atoms-1

             allocate(dmat(An_cell*atoms(i+1),An_cell*atoms(i+1)))

             dmat = cost_map(Apos( : , An_cell*i+1 : An_cell*(i+atoms(i+1)) ), &
                  Bpos( : , Bn_cell*i+1 : Bn_cell*i + int(Bn_cell*atoms(i+1)*frac)))

             call munkres(d,map(An_cell*i+1:An_cell*(i+atoms(i+1))),dmat,An_cell*atoms(i+1))

             dist = dist + d

             deallocate(dmat)

          enddo

       else if (type_dist == 6) then

          ! This distance takes in the unrotated Bpos
          
          ! Finds the best mapping
          An_cell = size(Apos,2)/sum(atoms)
          Bn_cell = size(Bpos,2)/sum(atoms)

          tBpos = free_trans(Bpos, tmat, vec) 
          
          dist = 0

          do i=0,n_atoms-1

             allocate(dmat(An_cell*atoms(i+1),An_cell*atoms(i+1)))

             dmat = cost_map(Apos( : , An_cell*i+1 : An_cell*(i+atoms(i+1)) ), &
                  tBpos( : , Bn_cell*i+1 : Bn_cell*(i+atoms(i+1))),frac)

             call munkres(d,map(An_cell*i+1:An_cell*(i+atoms(i+1))),dmat,An_cell*atoms(i+1))
             rot_mat(:,1) = norm(tmat(:,1))*(/1,0,0/)
             u_per = tmat(:,2) -  tmat(:,1)*dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1))
             rot_mat(:,2) = (/ dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1)), &
                  norm(u_per), 0.0d0/)
             rot_mat(:,3) = (/ dot_product(tmat(:,1), tmat(:,3))/norm(tmat(:,1)), &
                  dot_product(u_per, tmat(:,3))/norm(u_per), &
                  det(tmat,3)/(norm(tmat(:,1))*norm(u_per))/)

             ! This is to avoid seeing an inversion in 2D as a rotation in a
             ! prohibited axis. This is not necessary in 3D (but not wrong)
             rot_mat(2,:) = rot_mat(3,3)/abs(rot_mat(3,3))*rot_mat(2,:)
             rot_mat(3,3) = abs(rot_mat(3,3))
             
             ! Adding the mapping distance d and the stretching 
             
             dist = dist + d !+ sum(sqrt(sum((matmul(rot_mat,Bpos)-Bpos)**2,1))) TMP             
             deallocate(dmat)

          enddo
          
       endif
    
     end function dist

    subroutine dist_separate(dist, stretch, Apos, Bpos, frac, atoms, n_atoms, rot_mat, tmat, vec)

    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! position matrix
    double precision, dimension(size(Bpos,1),size(Bpos,2)) :: &
         tBpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms
    double precision, intent(in), dimension(3,3) :: &
         tmat
    double precision, intent(in), dimension(3) :: &
         vec
    double precision, intent(out), dimension(3,3) :: &
         rot_mat
    double precision, intent(in) :: &
         frac
    double precision, allocatable, dimension(:,:) :: &
         dmat
    double precision, dimension(3) :: &
         u_per
    integer, dimension(size(Apos,2)) :: &
         map ! List of index
    double precision, intent(out) :: &
         dist, &
         stretch
    double precision :: &
         d
    integer :: &
         An_cell, Bn_cell, & ! Number of cells
         i

          ! This distance takes in the unrotated Bpos
          
          ! Finds the best mapping
          An_cell = size(Apos,2)/sum(atoms)
          Bn_cell = size(Bpos,2)/sum(atoms)

          tBpos = free_trans(Bpos, tmat, vec) 
          
          dist = 0
          stretch = 0

          do i=0,n_atoms-1

             allocate(dmat(An_cell*atoms(i+1),An_cell*atoms(i+1)))

             dmat = cost_map(Apos( : , An_cell*i+1 : An_cell*(i+atoms(i+1)) ), &
                  tBpos( : , Bn_cell*i+1 : Bn_cell*(i+atoms(i+1))), frac)

             call munkres(d,map(An_cell*i+1:An_cell*(i+atoms(i+1))),dmat,An_cell*atoms(i+1))
             rot_mat(:,1) = norm(tmat(:,1))*(/1,0,0/)
             u_per = tmat(:,2) -  tmat(:,1)*dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1))
             rot_mat(:,2) = (/ dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1)), &
                  norm(u_per), 0.0d0/)
             rot_mat(:,3) = (/ dot_product(tmat(:,1), tmat(:,3))/norm(tmat(:,1)), &
                  dot_product(u_per, tmat(:,3))/norm(u_per), &
                  det(tmat,3)/(norm(tmat(:,1))*norm(u_per))/)

             ! This is to avoid seeing an inversion in 2D as a rotation in a
             ! prohibited axis. This is not necessary in 3D (but not wrong)
             rot_mat(2,:) = rot_mat(3,3)/abs(rot_mat(3,3))*rot_mat(2,:)
             rot_mat(3,3) = abs(rot_mat(3,3))

             ! Adding the mapping distance d and the stretching 
             
             stretch = stretch + sum(sqrt(sum((matmul(rot_mat,Bpos)-Bpos)**2,1)))
             dist = dist + d
             
             deallocate(dmat)

          enddo
          
        end subroutine dist_separate

        subroutine test_ds(dist, stretch, Apos, Bpos,n, frac, atoms, n_atoms, rot_mat, tmat, vec)

    integer, intent(in) :: &
         n_atoms, & ! Number of types of atoms
         n
    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos  ! position matrix
    double precision, dimension(size(Bpos,1),size(Bpos,2)) :: &
         tBpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms
    double precision, intent(in), dimension(3,3) :: &
         tmat
    double precision, intent(in), dimension(3) :: &
         vec
    double precision, intent(out), dimension(3,3) :: &
         rot_mat
    double precision, intent(in) :: &
         frac
    double precision, allocatable, dimension(:,:) :: &
         dmat
    double precision, dimension(3) :: &
         u_per
    integer, dimension(size(Apos,2)) :: &
         map ! List of index
    double precision, intent(out) :: &
         dist, &
         stretch
    double precision :: &
         d
    integer :: &
         An_cell, Bn_cell, & ! Number of cells
         i

          ! This distance takes in the unrotated Bpos
          
          ! Finds the best mapping
          An_cell = size(Apos,2)/sum(atoms)
          Bn_cell = size(Bpos,2)/sum(atoms)

          tBpos = free_trans(Bpos, tmat, vec) 
          
          dist = 0
          stretch = 0

          do i=0,n_atoms-1

             allocate(dmat(An_cell*atoms(i+1),An_cell*atoms(i+1)))

             dmat = cost_map(Apos( : , An_cell*i+1 : An_cell*(i+atoms(i+1)) ), &
                  tBpos( : , Bn_cell*i+1 : Bn_cell*(i+atoms(i+1))), frac)

             call munkres(d,map(An_cell*i+1:An_cell*(i+atoms(i+1))),dmat,An_cell*atoms(i+1))
             rot_mat(:,1) = norm(tmat(:,1))*(/1,0,0/)
             u_per = tmat(:,2) -  tmat(:,1)*dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1))
             rot_mat(:,2) = (/ dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1)), &
                  norm(u_per), 0.0d0/)
             rot_mat(:,3) = (/ dot_product(tmat(:,1), tmat(:,3))/norm(tmat(:,1)), &
                  dot_product(u_per, tmat(:,3))/norm(u_per), &
                  det(tmat,3)/(norm(tmat(:,1))*norm(u_per))/)

             ! This is to avoid seeing an inversion in 2D as a rotation in a
             ! prohibited axis. This is not necessary in 3D (but not wrong)
             rot_mat(2,:) = rot_mat(3,3)/abs(rot_mat(3,3))*rot_mat(2,:)
             rot_mat(3,3) = abs(rot_mat(3,3))
             
             ! Adding the mapping distance d and the stretching 
             
             stretch = stretch + sum(sqrt(sum((matmul(rot_mat,Bpos)-Bpos)**2,1)))
             dist = dist + d
             
             deallocate(dmat)

          enddo
          
  end subroutine test_ds
  
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
  
  subroutine mapping(map, dmin, Apos, Bpos, n, atoms, n_atoms)
    
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
         theta   ! Rotation angle

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
       call random_number(theta)
       call random_number(vec)
       call random_number(u)
       u = (u - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))) ! Recasts vector in the 8 octans
       vec = vec - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))
       u = u / norm(u) ! Normalizes
       vec = dmin * vec / sqrt(3.0d0) ! Moves less when close
       theta = 2*pi*theta

       newBpos = Bpos

       call trans(newBpos,n,theta,u,vec) ! MC step

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
         mat_out, &
         rmat
    
    double precision, dimension(3,1) :: &
         vec_tmp, & ! Temporary vector
         vec_out    ! Output vec

    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus, &   ! distance when substracting dx
         dist_map, &
         dist_stretch
    
    integer :: &
         i,j, & ! Iterator
         n

    double precision, parameter :: &
         dx = 1d-5  ! Derivative step

    n = size(Bpos,2)

    mat_out = mat
    vec_out = vec
    
    do j=1, n_iter

       ! mat
       do i=1,4

          mat_tmp = mat

          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) + dx
          dist_plus = dist( Apos, Bpos,frac,atoms,n_atoms, mat_tmp, vec)

          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - dx
          dist_minus = dist( Apos, Bpos,frac,atoms,n_atoms, mat_tmp, vec)

          mat_out((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - &
               rate1*(dist_plus**2 - dist_minus**2) / ( 2 * dx )

       enddo
       
       ! vec
       do i=1,2

          vec_tmp = vec

          vec_tmp(i,1) = vec(i,1) + dx
          dist_plus = dist( Apos, Bpos,frac,atoms, n_atoms, mat, vec_tmp)

          vec_tmp(i,1) = vec(i,1) - dx
          dist_minus = dist( Apos, Bpos,frac,atoms, n_atoms, mat, vec_tmp)

          vec_out(i,1) = vec(i,1) - &
               rate2*(dist_plus**2 - dist_minus**2) / ( 2 * dx )
          
       enddo

       ! call dist_separate(dist_map, dist_stretch, Apos, Bpos, frac,atoms, n_atoms,rmat, mat_out, vec_out)

       ! write(*,"(3(F5.3,X))") rmat
       
       ! print*, dist_map, dist_stretch 
       
       vec = vec_out
       mat = mat_out
       
    enddo
    
  end subroutine gradient_descent


    subroutine analytical_gd_rot(theta, u, vec, Apos, Bpos, n_iter, rate1, rate2)
    
    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2    ! Rate for disp
    
    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping
    
    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix
    
    double precision, intent(inout), dimension(3,1) :: &         
         vec, &     ! Translation vector
         u
         
    double precision, intent(inout) :: &
         theta    ! Transformation matrix

    double precision, dimension(3,3) :: &
         Px, Py, Pt, P, & ! Temporary transformation matrix
         Qx, Qy, Qt, Q, &
         Mx, My, Mt, M

    double precision, dimension(size(Bpos,2),1) :: &
         ones
    
    double precision :: &
         dist
    
    integer :: &
         i,j ! Iterator

    double precision, parameter :: &
         dx = 1d-5  ! Derivative step

    ones = 1.0d0
    
    do j=1, n_iter

       P = matmul(u,transpose(u))
       Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

       M = P + (eye() - P)*cos(theta) + Q*sin(theta)
       
       dist = sum(sqrt(sum((Apos - free_trans(Bpos,M,vec))**2,1)))

       print*, dist
       
       E = Apos - free_trans(Bpos,M,vec)
       E = E / spread(sqrt(sum(E**2,1)),1,3)

       ! ! ux and uy is for 3D only
       ! Px = transpose(reshape((/2*u(1,1), &
       !      u(2,1) , &
       !      (u(1,1)*(1-u(2,1)**2) - 2*u(1,1)**3)/sqrt(u(1,1)**2*(1-u(2,1)**2)-u(1,1)**4), &
       !      u(2,1)  , &
       !      u(2,1)**2, &
       !      -u(1,1)*u(2,1)**2/sqrt(u(2,1)**2*(1-u(1,1)**2)-u(2,1)**4), &
       !      (u(1,1)*(1-u(2,1)**2) - 2*u(1,1)**3)/sqrt(u(1,1)**2*(1-u(2,1)**2)-u(1,1)**4) , &
       !      -u(1,1)*u(2,1)**2/sqrt(u(2,1)**2*(1-u(1,1)**2)-u(2,1)**4), &
       !      -2*u(1,1)/), &
       !      (/3,3/)))
       ! Qx = transpose(reshape((/0.0d0, &
       !      u(1,1)/sqrt(1-u(1,1)**2-u(2,1)**2), &
       !      u(2,1), &
       !      -u(1,1)/sqrt(1-u(1,1)**2-u(2,1)**2), &
       !      0.0d0, &
       !      -1.0d0, &
       !      -u(2,1), &
       !      1.0d0, &
       !      0.0d0/), &
       !      (/3,3/)))

       ! Mx = Px + (eye() - Px)*cos(theta) + Qx*sin(theta)

       ! Py = transpose(reshape((/u(1,1)**2, &
       !      u(1,1) , &
       !      -u(2,1)*u(1,1)**2/sqrt(u(1,1)**2*(1-u(2,1)**2)-u(1,1)**4), &
       !      u(1,1)  , &
       !      2*u(2,1), &
       !      (u(2,1)*(1-u(1,1)**2) - 2*u(2,1)**3)/sqrt(u(2,1)**2*(1-u(1,1)**2)-u(2,1)**4), &
       !      -u(2,1)*u(1,1)**2/sqrt(u(1,1)**2*(1-u(2,1)**2)-u(1,1)**4) , &
       !      (u(2,1)*(1-u(1,1)**2) - 2*u(2,1)**3)/sqrt(u(2,1)**2*(1-u(1,1)**2)-u(2,1)**4), &
       !      -2*u(2,1)/), &
       !      (/3,3/)))
       ! Qy = transpose(reshape((/0.0d0, &
       !      u(2,1)/sqrt(1-u(1,1)**2-u(2,1)**2), &
       !      1.0d0, &
       !      -u(2,1)/sqrt(1-u(1,1)**2-u(2,1)**2), &
       !      0.0d0, &
       !      u(1,1), &
       !      -1.0d0, &
       !      u(1,1), &
       !      0.0d0/), &
       !      (/3,3/)))

       ! My = Py + (eye() - Py)*cos(theta) + Qy*sin(theta)

       Pt = matmul(u,transpose(u))
       Qt = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

       Mt = Pt - (eye() - Pt)*sin(theta) + Qt*cos(theta)
       
       ! u(1,1) = u(1,1) + rate1*dist*sum(matmul(E,transpose(Bpos)) * Mx)
       ! u(2,1) = u(2,1) + rate1*dist*sum(matmul(E,transpose(Bpos)) * My)
       ! u(3,1) = sqrt(1-u(1,1)**2-u(2,1)**2)
       theta = theta + rate1*dist*sum(matmul(E,transpose(Bpos)) * Mt)
       vec = vec + rate2*dist*matmul(E,ones)

       
    enddo
    
  end subroutine analytical_gd_rot
  
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

    dist_cur = dist( Apos, Bpos, frac,atoms, n_atoms, mat, vec)
    
    dist_prev = dist_cur
    dist_min = dist_cur
    mat_min = mat
    vec_min = vec

    vec_out = 0
    mat_out = 0
    
    do j=1, n_iter

       call random_number(rand_rate1)

       rand_rate2 = 100**rand_rate1*rate2
       rand_rate1 = 100**rand_rate1*rate1

       ! mat
       do i=1,4

          mat_tmp = mat

          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) + dx
          dist_plus = dist( Apos, Bpos,frac,atoms,n_atoms, mat_tmp, vec)
          
          mat_tmp((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - dx
          dist_minus = dist( Apos, Bpos,frac,atoms,n_atoms, mat_tmp, vec)

          mat_out((i-1)/2+1,modulo(i-1,2)+1) = mat((i-1)/2+1,modulo(i-1,2)+1) - rand_rate1*(dist_plus - dist_minus) / ( 2 * dx )

       enddo

       
       ! vec
       do i=1,2

          vec_tmp = vec

          vec_tmp(i,1) = vec(i,1) + dx
          dist_plus = dist( Apos, Bpos, frac,atoms, n_atoms, mat, vec_tmp)

          vec_tmp(i,1) = vec(i,1) - dx
          dist_minus = dist( Apos, Bpos, frac,atoms, n_atoms, mat, vec_tmp)

          vec_out(i,1) = vec(i,1) - rand_rate2*(dist_plus - dist_minus) / ( 2 * dx )
          
       enddo
       
       dist_cur = dist( Apos, Bpos, frac,atoms, n_atoms, mat, vec_tmp)

       
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

    subroutine gradient_descent_explore(mat, vec, Apos, Bpos, frac, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2)
    ! New Gradient Descent Random  
    
    integer, intent(in) :: &
         n_iter, n_atoms, & ! Number of atoms
         n_ana, &
         n_conv

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         frac ! Fraction of A and B to use in optimisation
         
    double precision :: &
         rand_rate1, & ! Rate for angles
         rand_rate2 ! Rate for disp
    
    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Centered position of the atoms
    
    double precision, dimension(3,size(Bpos,2)) :: &
         pos, postmp, & ! position matrix
         tBpos

    double precision, dimension(3,sum(int(frac*atoms*size(Apos,2)/sum(atoms)))) :: &
         Apos_mapped, & ! position matrix
         Bpos_opt, &
         Bpos_test, &
         Bpos_test2
    
    double precision, intent(out), dimension(3,1) :: &
         vec     ! Translation vector

    double precision, intent(out), dimension(3,3) :: &
         mat    ! Transformation matrix

    double precision, dimension(3,3) :: &
         mat_tmp, & ! Temporary transformation matrix
         mat_out, & ! Output transformation matrix
         mat_min, &
         P,Q
    
    double precision, dimension(3,1) :: &
         vec_out, & ! Output vec
         vec_min, &
         u, &
         u_min

    double precision, allocatable, dimension(:,:) :: &
         dmat

    integer, allocatable, dimension(:) :: &
         map
         
    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus, & ! distance when substracting dx
         accept, & ! Accept step
         dist_cur, &
         dist_map, &
         dist_stretch, &
         dist_min, &
         mul_vec, &
         theta, &
         theta_min
    
    integer :: &
         i,j,k,l, & ! Iterator
         id, idx, &
         n, & ! Size of Bpos
         An_cell, Bn_cell

    double precision, parameter :: &
         dx = 1d-5, &  ! Derivative step
         max_vec = 0.3d0, &
         pi = 3.141592653589793d0

    An_cell = size(Apos,2)/sum(atoms)
    Bn_cell = size(Bpos,2)/sum(atoms)
    
    call init_random_seed()

    mul_vec = max_vec*2*maxval(sqrt(sum(Bpos**2,1)))/sqrt(2.0d0) ! 2D only sqrt(3) in 3D
    
    dist_min = sum(sqrt(sum((Apos - Bpos)**2,1)))
    
    do j=1, n_iter

       call random_number(theta)
       call random_number(u)
       call random_number(vec)
       
       theta = theta*2*pi

       vec = vec - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))
       vec = vec * mul_vec
       vec(3,1) = 0.0d0 ! 2D only

       ! ! 3D only
       ! u = u - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))
       ! u = u / norm(u)
       ! u(3,1) = abs(u(3,1))

       u =  reshape((/0.0d0,0.0d0,1.0d0/),(/3,1/)) ! 2D only

       print*, "New u"
       
       do k=1, n_conv
          
          P = matmul(u,transpose(u))
          Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1), &
               u(3,1),0.0d0,-u(1,1), &
               -u(2,1),u(1,1),0.0d0/),(/3,3/)))

          mat = P + (eye() - P)*cos(theta) + Q*sin(theta)
          
          tBpos = free_trans(Bpos,mat,vec)

          id = 0
          idx = 0
          n = 0
          
          do i=0,n_atoms-1

             id = id + n

             n = An_cell*atoms(i+1)

             allocate(dmat(n,n), map(n))

             dmat = cost_map(Apos( : , id + 1 : id + n ), &
                  tBpos( : , id + 1 : id + n ),frac)

             call munkres(dist_map, map, dmat, n)

             map = map + id

             do l=1, int(n*frac)
                idx = idx + 1 
                Apos_mapped(:,idx) = Apos(:,map(l))
                Bpos_opt(:,idx) = Bpos(:,id+l)
                Bpos_test2(:,idx) = Bpos(:,id+l)
                Bpos_test(:,idx:idx) = free_trans(Bpos(:,id+l:id+l),mat,vec)
             enddo

             deallocate(dmat,map)
          
          enddo

          print*, "Remapped", dist_map
          
          call analytical_gd_rot(theta, u, vec, Apos_mapped, Bpos_opt, n_ana, rate1, rate2)
          
       enddo
          
       
       P = matmul(u,transpose(u))
       Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1), &
            u(3,1),0.0d0,-u(1,1), &
            -u(2,1),u(1,1),0.0d0/),(/3,3/)))
       
       mat = P + (eye() - P)*cos(theta) + Q*sin(theta)

       dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,mat,vec))**2,1)))
       
       if (dist_cur < dist_min) then
          dist_min = dist_cur
          mat_min = mat
          vec_min = vec
       endif

    enddo
    
    mat = mat_min
    
    vec = vec_min
    
  end subroutine gradient_descent_explore

  
  subroutine fastmapping(tmat, rmat, map, dmin, &
       Apos, na, Bpos, nb, &
       frac, Acell, iAcell, atoms, n_atoms, &
       n_iter, n_ana, n_conv, &
       rate1, rate2)
    
    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2     ! Rate of the gradient descent for displacement
         
    integer, intent(in) :: &
         n_iter, &   ! Number of iteration of the gradient descent
         n_ana, &
         n_conv
         
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
         theta, & ! Rotation angle
         d, &
         dist_map, &
         dist_stretch

    double precision, allocatable, dimension(:,:) :: &
         dmat

    double precision, dimension(na,nb) :: &
         mat

    double precision, intent(out), dimension(3,3) :: &
         tmat, & ! Transformation matrix
         rmat    ! Rotation matrix
    
    integer :: &
         i, &
         n, id, &
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
    tmat(3,:) = 0.0d0 ! For 2D only
    tmat(:,3) = 0.0d0 ! For 2D only
    tmat(3,3) = 1.0d0 ! For 2D only

    vec = 0

    call gradient_descent_explore(tmat, vec, Apos, Bpos, &
         frac, atoms,n_atoms,n_iter, n_ana, n_conv, rate1, rate2)
    
    ! call gradient_descent(tmat, vec, Apos, Bpos, &
    !      frac, atoms,n_atoms,n_iter, rate1/10, rate2/10)
    
    !call dist_separate(dist_map, dist_stretch, Apos, Bpos,frac,atoms, n_atoms,rmat, tmat, vec)


    !print*, "Final distances:", dist_map, dist_stretch
    
    !dmin = dist_map + dist_stretch
    
    Bpos = free_trans(Bpos, tmat, vec)
    
    call equivalent_center(Bpos,Acell,iAcell)
    
    ! Finds the best mapping
    An_cell = size(Apos,2)/sum(atoms)
    Bn_cell = size(Bpos,2)/sum(atoms)

    !dmin = 0

    id = 0
    n = 0
    
    do i=0,n_atoms-1

       id = id + n
       
       n = An_cell*atoms(i+1)
       
       allocate(dmat(n,n))

       call cpu_time(start)
       dmat = cost_map(Apos( : , id + 1 : id + n ), &
                       Bpos( : , id + 1 : id + n ), frac)
       call cpu_time(finish)

       print*, "dist_time", finish - start 
       
       call cpu_time(start)
       call munkres(d,map(id + 1 : id + n),dmat,An_cell*atoms(i+1))
       call cpu_time(finish)

       print*, "munkres_time", finish - start
       
       !dmin = dmin + d

       deallocate(dmat)

       map(id + 1 : id + n) = map(id + 1 : id + n) + id

    enddo

    dmin = d

    rmat = 0

    print*, "tmat",tmat 
    print*, "rmat",rmat
    print*, "map ",map 
    print*, "dmin",dmin
    print*, "Apos",Apos
    print*, "Bpos",Bpos
    
    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(*,"(10(F5.3,X))") mat   
    ! call munkres(dmin,map,cost(Apos,Bpos,n),n)
    
  end subroutine fastmapping
  
end module transform
