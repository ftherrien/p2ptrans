module transform
  
  use utils
  !use hungarian
  use jv 
  use potential

  implicit none

  public :: &
       fastoptimization, &
       intoptimization, &
       closest, &
       fixed_tmat, &
       fixed_tmat_int, &
       optimize_vec

  private :: &
       init_random_seed, &
       classify,&
       center, &
       cost_map, &
       mapping, &
       no_mapping, &
       analytical_gd_rot, &
       analytical_gd_vec, &
       analytical_gd_free, &
       analytical_gd_slant, &
       analytical_gd_std, &
       gradient_descent_explore, &
       gradient_descent_explore_free, &
       find_peaks, &
       adjust, &
       remove_slanting, &
       classification, &
       pi, rotmin

  double precision,  parameter :: &
       pi = 3.141592653589793d0  

  ! Rotmin can easily be stuck, it should be replaced
  ! with a proper RU decomposition
  logical, parameter ::  rotmin = .true.

contains
  
  subroutine init_random_seed()
    ! Copied form the GCC docs: https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED

    use omp_lib

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (OMP_get_thread_num()+1) * (/ (i - 1, i = 1, n) /)
    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)

  end subroutine init_random_seed

  subroutine classify(n, n_classes, classes_list, disps, tol)

    ! Classifies the displacements into n classes and gives them a label classes_list

    double precision, intent(in), dimension(:,:) :: &
         disps ! Displacements (aka connections)

    double precision, dimension(size(disps,1),size(disps,2)) :: &
         classes ! Vectors representing each class

    integer, intent(out), dimension(size(disps,2)) :: &
         classes_list ! List of classes

    double precision, intent(in) :: &
         tol ! Tolerence factor for classification

    integer, intent(out), dimension(size(disps,2)) :: &
         n_classes ! Number of vectors in each class

    integer, intent(out) :: &
         n ! Number of classes

    logical :: &
         classified

    integer :: &
         i,j

    classes(:,1) = disps(:,1)
    n_classes = 1
    n = 1
    classes_list = 1
    do i=2,size(disps,2)
       classified = .false.
       do j=1,n

          if (norm(classes(:,j) - disps(:,i)) < tol) then

             classes(:,j) = (n_classes(j)*classes(:,j) + disps(:,i)) / (n_classes(j) + 1)
             n_classes(j) = n_classes(j) + 1
             classes_list(i) = j
             classified = .true.
             exit
          endif
       enddo

       if (.not. classified) then
          n = n + 1    
          classes(:,n) = disps(:,i)
          classes_list(i) = n
       endif
    enddo

  end subroutine classify

  function cost_map(Apos,Bpos, n_A, n_B, pot, param) result(cost)

    ! Creates the cost matrix (Euclidean distance)

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! Rotation axis
    integer, intent(in) :: &
         n_A, n_B
    double precision, dimension(size(Apos,2),size(Apos,2)) :: &
         cost
    integer :: &
         i,j ! Iterators
    character*(*), intent(in) :: &
         pot
    double precision, intent(in) :: &
         param

    cost = 0

    do i=1,size(Apos,2)
       do j=1,size(Bpos,2)
          if (j <= n_B) then
             cost(i,j) = distance(Apos(:,i:i),Bpos(:,j:j), eye(), (/0.0d0,0.0d0,0.0d0/), pot, param)
          elseif (i <= n_A) then  
             cost(i,j) = 1000 ! Very high cost
          endif
       enddo
    enddo

  end function cost_map


  subroutine closest(cost, Apos, Bpos, n)

    integer, intent(in) :: &
         n ! Total number of atoms

    double precision, intent(in), dimension(3,n) :: &
         Apos, Bpos ! Position of the atoms

    ! integer, intent(out), dimension(n) :: &
    !      closest_list ! Position of the atoms

    double precision, intent(out), dimension(n,n) :: &
         cost

    cost = cost_map(Apos, Bpos, n, n, "Euclidean", 0.0d0)

    ! closest_list =  minloc(cost,2)

  end subroutine closest
  
  subroutine no_mapping(Apos_mapped, Bpos_opt, Apos, Bpos, fracB, atoms, n_atoms)

    ! Consider the order in which the atoms are given to be the mapping, creates
    ! the equivalent outputs to mapping

    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms ! Number of sections per atom (same as in p2ptrans.py)
    double precision, intent(in) :: &
         fracB ! Fraction of the spheres to be used
    double precision, intent(out), &
         dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
         Apos_mapped, Bpos_opt  ! position matrix

    integer :: &
         An_cell, & ! Number of cells
         i,l, &
         id, idx, &
         n, n_B

    An_cell = size(Apos,2)/sum(atoms)

    id = 0
    idx = 0
    n = 0

    do i=0,n_atoms-1

       id = id + n

       n = An_cell*atoms(i+1)
       n_B = int(An_cell*fracB)*atoms(i+1)

       do l=1, n_B
          idx = idx +  1
          Apos_mapped(:,idx) = Apos(:,id+l)
          Bpos_opt(:,idx) = Bpos(:,id+l)
       enddo

    enddo

  end subroutine no_mapping

  subroutine mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, fracA, fracB, atoms, n_atoms, &
       pot, param)

    ! Mapps the atom according the the cost function using the lap package

    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos, tBpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms ! Number of sections per atom (same as in p2ptrans.py)
    double precision, intent(in) :: &
         fracA, fracB ! Fraction of the spheres to be used
    double precision, intent(out), &
         dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
         Apos_mapped, Bpos_opt  ! position matrix
    double precision, allocatable, dimension(:,:) :: &
         dmat
    integer, allocatable, dimension(:) :: &
         map
    double precision :: &
         dist_map
    integer :: &
         An_cell, & ! Number of cells
         i,l, &
         id, idx, &
         n, n_A, n_B
    character*(*), intent(in) :: &
         pot
    double precision, intent(in) :: &
         param

    An_cell = size(Apos,2)/sum(atoms)

    id = 0
    idx = 0
    n = 0

    do i=0,n_atoms-1

       id = id + n

       n = An_cell*atoms(i+1)
       n_B = int(An_cell*fracB)*atoms(i+1)
       n_A = int(An_cell*fracA)*atoms(i+1)

       allocate(dmat(n,n), map(n))

       dmat = cost_map(Apos( : , id + 1 : id + n ), &
            tBpos( : , id + 1 : id + n ),n_A, n_B, &
            pot, param)

       !call munkres(dist_map, map, dmat, n)
       call lapjv(dist_map, map, dmat, n)

       map = map + id

       do l=1, n_B
          idx = idx + 1
          Apos_mapped(:,idx) = Apos(:,map(l))
          Bpos_opt(:,idx) = Bpos(:,id+l)
       enddo

       deallocate(dmat,map)

    enddo

  end subroutine mapping

  subroutine analytical_gd_rot(twodim, angles, vec, Apos, Bpos, n_iter, rate1, rate2, tol, pot, param)

    ! Gradient descent with respect to the 3 angles of rotation
    ! ADJSUT has been implemented
    ! TODO: Could be replaced with direct SVD minimization

    logical, intent(in) :: &
         twodim
    
    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         tol

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec       ! Translation vector

    double precision, dimension(3,1) :: &         
         u         ! Axis of rotation

    double precision, intent(inout), dimension(3) :: &
         angles    ! 3 vectors of the rotation: 1. angle of rotation 2&3. Spherical angles of the axis of rotation

    double precision, dimension(3,3) :: &
         P1, P2, P3, & ! Temporary transformation matrix
         Q1, Q2, Q3, &
         M1, M2, M3, M

    double precision, dimension(3) :: &
         angles_prev

    double precision, dimension(3,1) :: &
         vec_prev

    double precision, dimension(size(Bpos,2),1) :: &
         ones

    double precision :: &
         dist, &
         dist_prev, &
         dist_init, &
         s2, s3, & ! Sine of angle 1 and 2
         c2, c3, & ! Cosine of angle 1 and 2
         nat

    integer :: &
         j      ! Iterator

    character*(*) :: &
         pot

    double precision :: &
         param


    nat = dble(size(Apos,2))

    ones = 1.0d0

    M = rot_mat(angles)
    dist = 0.0d0
    dist_init = 0.0d0
    dist_prev = tol + 1.0d0
    angles_prev = angles
    vec_prev = vec

    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       M = rot_mat(angles)

       dist_prev = dist
       dist = distance(Apos, Bpos, M, vec, pot, param)

       if (j==1) then
          dist_prev = dist + tol + 1.0d0
          dist_init = dist
       endif

       E = derivative(Apos, Bpos, M, vec, pot, param)

       s2 = sin(angles(2))
       s3 = sin(angles(3))
       c2 = cos(angles(2))
       c3 = cos(angles(3))

       P2 = transpose(reshape((/2*s2*c2*c3**2, &
            2*s2*c2*s3*c3, &
            (1 - 2*s2**2) * c3, &
            2*s2*c2*s3*c3, &
            2*s2*c2*s3**2, &
            (1 - 2*s2**2) * s3, &
            (1 - 2*s2**2) * c3, &
            (1 - 2*s2**2) * s3, &
            -2*c2*s2/), &
            (/3,3/)))
       Q2 = transpose(reshape((/0.0d0, &
            s2, &
            c2 * s3, &
            -s2, &
            0.0d0, &
            -c2*c3, &
            -c2*s3, &
            c2*c3, &
            0.0d0/), &
            (/3,3/)))

       M2 = P2 - P2*cos(angles(1)) + Q2*sin(angles(1))

       P3 = transpose(reshape((/-2*s2**2*c3*s3, &
            (1 - 2*s3**2) * s2**2, &
            -s2*c2*s3, &
            (1 - 2*s3**2) * s2**2, &
            2*s2**2*s3*c3, &
            s2*c2*c3, &
            -s2*c2*s3, &
            s2*c2*c3, &
            0.0d0/), &
            (/3,3/)))
       Q3 = transpose(reshape((/0.0d0, &
            0.0d0, &
            s2*c3, &
            0.0d0, &
            0.0d0, &
            s2*s3, &
            -s2*c3, &
            -s2*s3, &
            0.0d0/), &
            (/3,3/)))

       M3 = P3 - P3*cos(angles(1)) + Q3*sin(angles(1))

       u(1,1) = s2 * c3
       u(2,1) = s2 * s3
       u(3,1) = c2

       P1 = matmul(u,transpose(u))
       Q1 = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

       M1 = - (eye() - P1)*sin(angles(1)) + Q1*cos(angles(1))

       ! print*, dist, angles(1), dist_prev - dist, dist_init, "ROT"

       if(dist > dist_prev .or. abs(dist_init) < 1.0d-8) then
          ! print*, "ADJUST", dist, dist_init
          dist_init = 10 * dist_init
          angles = angles_prev
          vec = vec_prev
       else

          angles_prev = angles
          vec_prev = vec

          angles(1) = angles(1) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M1)
          if (.not. twodim) then
             angles(2) = angles(2) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M2)
             angles(3) = angles(3) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M3)
          endif
          vec       = vec       + rate2 * dist / dist_init * matmul(E,ones)

       endif

    enddo

    ! print*, "Fini", j, dist_prev - dist, dist

  end subroutine analytical_gd_rot

  subroutine analytical_gd_free(twodim, tmat, vec, Apos, Bpos, n_iter, rate1, rate2, tol, pot, param)

    ! Gradient descent with respect to a linear transfromation (3x3 matrix)
    ! ADJUST has been implemented

    logical, intent(in) :: &
         twodim

    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         tol

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(inout), dimension(3,3) :: &         
         tmat     ! Transformation matrix

    double precision, dimension(3,3) :: &         
         tmat_prev     ! Transformation matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec     ! Translation vector

    double precision, dimension(3,1) :: &         
         vec_prev     ! Translation vector

    double precision, dimension(size(Bpos,2),1) :: &
         ones

    character*(*), intent(in) :: &
         pot

    double precision, intent(in) :: &
         param

    double precision :: &
         dist, &
         dist_prev, &
         dist_init, &
         nat

    integer :: &
         j ! Iterator

    nat = dble(size(Apos,2))

    ones = 1.0d0

    dist = 0.0d0
    dist_init = 0.0d0
    dist_prev = tol+1
    tmat_prev = tmat
    vec_prev = vec

    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1
       
       dist_prev = dist
       dist = distance(Apos, Bpos, tmat, vec, pot, param)

       if (j==1) then
          dist_prev = dist + tol + 1.0d0
          dist_init = dist
       endif

       ! print*, dist, dist_prev - dist, "FREE", j

       E = derivative(Apos, Bpos, tmat, vec, pot, param)

       if(dist > dist_prev .or. abs(dist_init) < 1.0d-8) then
          dist_init = 10 * dist_init
          tmat = tmat_prev
          vec = vec_prev
       else

          tmat_prev = tmat
          vec_prev = vec

          if (twodim) then
             tmat(1:2,1:2) = tmat(1:2,1:2) + rate1 * dist / dist_init * matmul(E(1:2,:),transpose(Bpos(1:2,:)))
          else
             tmat = tmat + rate1 * dist / dist_init * matmul(E,transpose(Bpos))
          endif

          vec = vec   + rate2 * dist / dist_init * matmul(E,ones)
          
       endif

    enddo

    ! print*, "free", j, dist  

  end subroutine analytical_gd_free
  
  subroutine analytical_gd_vec(straighten, tmat, vec, Apos, Bpos, n_iter, rate2, tol, pot, param)

    ! Gradient descent with respect to the vector only (3x3 matrix)
    ! ADJUST has been implemented
    logical, intent(in) :: &
         straighten

    integer, intent(in) :: &
         n_iter ! Number of atoms

    
    double precision, intent(in) :: &
         rate2, & ! Rate for disp
         tol

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(in), dimension(3,3) :: &         
         tmat     ! Transformation matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec     ! Translation vector

    double precision, dimension(3,1) :: &         
         vec_prev, &    ! Translation vector
         vec_in

    double precision, dimension(size(Bpos,2),1) :: &
         ones

    character*(*), intent(in) :: &
         pot

    double precision, intent(in) :: &
         param

    double precision :: &
         dist, &
         dist_prev, &
         dist_init, &
         nat, &
         sum_pot, &
         sum_grad

    integer :: &
         j, i, k ! Iterator

    double precision, parameter :: &
         c = 0.001

    nat = dble(size(Apos,2))

    ones = 1.0d0

    dist = 0.0d0
    dist_init = 0.0d0
    dist_prev = tol+1
    vec_prev = vec
    vec_in = vec

    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1
       
       dist_prev = dist
       sum_pot = 0.0d0
       if (straighten) then
       do k=1,size(Apos,2)
          if (abs(Apos(3,k) - dot_product(tmat(3,:),Bpos(:,k)) - vec(3,1)) < param) then
             sum_pot = sum_pot + c/(2**6 - 1)*(param**6/(Apos(3,k) - dot_product(tmat(3,:),Bpos(:,k)) - vec(3,1))**6 - 1.0)
          endif
       enddo
       endif
       
       dist = distance(Apos, Bpos, tmat, vec, pot, param) + sum_pot

       if (j==1) then
          dist_prev = dist + tol + 1.0d0
          dist_init = dist
          ! print*, "free", j, dist, sum_pot  
       endif

       ! print*, dist, sum_pot, vec(3,1), dist_prev - dist, "FREE", j

       E = derivative(Apos, Bpos, tmat, vec, pot, param)

       if(dist > dist_prev .or. abs(dist_init) < 1.0d-8 .or. vec(3,1)*vec_in(3,1) < 0.0d0) then
          dist = dist_prev + tol + 1.0d0
          dist_init = 10 * dist_init
          vec = vec_prev
       else

          vec_prev = vec

          do i=1,3
            if (i==3 .and. straighten) then
          
               sum_grad = 0.0d0
               do k=1,size(Apos,2)
                  if (abs(Apos(3,k) - dot_product(tmat(3,:),Bpos(:,k)) - vec(3,1)) < param) then
                     sum_grad = sum_grad + c/(2**6 -1) / &
                          (Apos(3,k) - dot_product(tmat(3,:), Bpos(:,k)) - vec(3,1))**7 ! Divided by 6*param**6 because derivative of LJ also gets rid of that factor
                  endif
               enddo
               
               vec(i,1) = vec(i,1) + rate2 * dist / dist_init * ( sum(E(i,:)) - sum_grad )
            else
               vec(i,1) = vec(i,1) + rate2 * dist / dist_init * sum(E(i,:))
            endif
          enddo
       endif

    enddo

    ! print*, "free", j, dist, sum_pot  

  end subroutine analytical_gd_vec

  subroutine analytical_gd_std(twodim, std, tmat, vec, Apos, Bpos, n_classes, classes, n_iter, rate1, rate2, tol)

    ! Gradient descent with respect to a linear transformation to minimize the class specific
    ! standard deviation
    ! ADJUST has been implemented

    logical, intent(in) :: &
         twodim

    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         tol

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, allocatable, dimension(:,:) :: &
         E ! position matrix

    integer, intent(in), dimension(:) :: &
         n_classes

    integer, intent(in), dimension(:) :: &
         classes

    double precision, intent(inout), dimension(3,3) :: &         
         tmat     ! Transformation matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec     ! Translation vector

    double precision, allocatable, dimension(:,:) :: &
         Apos_class, &
         Bpos_class

    double precision, dimension(3,3) :: &         
         tmat_grad, &     ! Transformation matrix
         tmat_prev

    double precision, dimension(3,1) :: &         
         vec_grad, &      ! Translation vector
         vec_prev

    double precision, allocatable, dimension(:,:) :: &
         ones

    double precision :: &
         std, &
         std_prev, &
         std_2prev, &
         std_init, &
         nat

    integer :: &
         j, k, l, i ! Iterator

    nat = dble(size(Apos,2))

    std  = 1.0d0
    std_prev = 2.0d0 + tol
    std_init = std
    tmat_prev = tmat
    vec_prev = vec

    j=0
    do while (j < n_iter .and. abs(std - std_prev) > tol)
       j=j+1

       std_2prev = std_prev
       std_prev = std

       tmat_grad = 0.0d0
       vec_grad = 0.0d0
       std = 0.0d0
       do i=1,size(n_classes)
          allocate(Apos_class(3,n_classes(i)), Bpos_class(3,n_classes(i)), &
               ones(1,n_classes(i)), E(3,n_classes(i)))

          ones = 1.0d0

          l=0
          do k=1,size(classes)
             if (classes(k) == i) then
                l=l+1
                Apos_class(:,l) = Apos(:,k)
                Bpos_class(:,l) = Bpos(:,k)
             endif
          enddo

          std = std + sum((Apos_class - free_trans(Bpos_class,tmat,vec))**2) - &
               norm(sum(Apos_class - free_trans(Bpos_class,tmat,vec),2))**2/size(Apos_class,2)


          E = Apos_class - free_trans(Bpos_class,tmat,vec)
          E = E - matmul(reshape(sum(E,2),(/3,1/)), ones)/size(E,2)

          tmat_grad = tmat_grad + rate1 * std_prev / std_init * matmul(E,transpose(Bpos_class))
          vec_grad = vec_grad + rate2 * std_prev / std_init * matmul(E,transpose(ones))

          deallocate(Apos_class, Bpos_class, ones, E)

       enddo

       if (j==1) then
          std_prev = std + tol + 1.0d0
          std_init = std
       endif

       if(std > std_prev) then
          std_init = 10 * std_init
          !print*, "ADJUST", std, std_prev, std_init
          tmat = tmat_prev
          vec = vec_prev
          std = std_2prev
          std_prev = std + 1.0d0 + tol
       elseif(std_init < tol .or. std < tol) then
          std_prev = std
       else
          tmat_prev = tmat
          vec_prev = vec

          !print*, std, std_prev, std_prev - std, "STD"
          if (twodim) then
             tmat(1:2,1:2) = tmat(1:2,1:2) + tmat_grad(1:2,1:2)
          else
             tmat = tmat + tmat_grad
             vec = vec + vec_grad
          endif
          
       endif

    enddo

    !print*, "end", j, std
    
  end subroutine analytical_gd_std

  subroutine analytical_gd_slant(slant, vec, Apos, Bpos, n_iter, rate1, rate2, tol)

    ! Gradient descent that minimizes the distance with respect to a slanting matrix:
    ! RKRt where K = [[1k0],[010],[001]]
    ! This is NOT used by default
    ! ADJUST has not been implemented

    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         tol

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, intent(out), dimension(3,3) :: &         
         slant     ! Transformation matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec       ! Translation vector

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, dimension(3,1) :: &         
         v1, v2, & ! Vectors of the plane of slanting
         dv1d1, &  ! Derivative of the plane vector
         dv1d2, &
         dv2d1, &
         dv2d2, &
         dv2d3

    double precision, dimension(3) :: &
         angles    ! 3 vectors of the rotation: 1. angle of rotation 2&3. Spherical angles of the axis of rotation

    double precision, dimension(3,3) :: &
         M1, M2, M3, Mk

    double precision, dimension(size(Bpos,2),1) :: &
         ones

    double precision :: &
         dist, &
         dist_prev, &
         dist_init, &
         s1, s2, s3, & ! Sine of angles
         c1, c2, c3, & ! Cosine of angles
         k, & ! Slanting factor
         nat

    integer :: &
         j ! Iterator

    nat = dble(size(Apos,2))

    ones = 1.0d0

    call init_random_seed()

    ! call random_number(angles)

    ! angles = (/0.0d0, pi/2.0d0, pi/2.0d0/) 

    angles = 0.0d0

    k = 1.0d0

    slant = eye()
    dist = 0.0d0
    dist_init = 0.0d0
    dist_prev = tol + 1.0d0

    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       s1 = sin(angles(1))
       s2 = sin(angles(2))
       s3 = sin(angles(3))
       c1 = cos(angles(1))
       c2 = cos(angles(2))
       c3 = cos(angles(3))

       v1 = reshape((/s2*c1, s2*s1, c2/),(/3,1/))

       v2 = reshape((/s3*c2*c1-c3*s1, &
            s3*c2*s1+c3*c1, &
            -s3*s2/),(/3,1/))       

       dv1d1 = reshape((/-s2*s1, s2*c1, 0.0d0/),(/3,1/))
       dv1d2 = reshape((/c2*c1, c2*s1, -s2/),(/3,1/))

       dv2d1 = reshape((/-s3*c2*s1 - c3*c1, &
            s3*c2*c1 - c3*s1, &
            0.0d0/),(/3,1/))
       dv2d2 = reshape((/-s3*s2*c1, &
            -s3*s2*s1, &
            -s3*c2/),(/3,1/))
       dv2d3 = reshape((/c3*c2*c1 + s3*s1, &
            c3*c2*s1 - s3*c1, &
            -c3*s2/),(/3,1/))

       slant = k*matmul(v1,transpose(v2)) + eye()

       dist_prev = dist
       dist = distance(Apos, Bpos, slant, vec, "Euclidean", 0.0d0)

       if (j==1) then
          dist_init = dist
       endif

       ! print*, dist, dist_prev - dist, "SLANT"

       E = Apos - free_trans(Bpos,slant,vec)
       E = E / spread(sqrt(sum(E**2,1)),1,3)

       M1 = k * (matmul(dv1d1,transpose(v2)) + matmul(v1, transpose(dv2d1))) 
       M2 = k * (matmul(dv1d2,transpose(v2)) + matmul(v1, transpose(dv2d2)))
       M3 = k * (matmul(v1, transpose(dv2d3)))
       Mk = matmul(v1,transpose(v2))

       angles(1) = angles(1) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M1)
       angles(2) = angles(2) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M2)
       angles(3) = angles(3) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M3)
       k = k + rate1 * 0.01 * dist / dist_init  * sum(matmul(E,transpose(Bpos)) * Mk)
       vec = vec + rate2 * dist / dist_init  * matmul(E,ones)

    enddo

    ! print*, dist

  end subroutine analytical_gd_slant

  subroutine gradient_descent_explore(twodim, angles, vec, Apos, Bpos, cell, icell, &
       fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2, tol, zdist, &
       pot, param, outdir)
    ! New Gradient Descent Random

    ! Randomly prepare the structures and then minimize the distance with respect to a rotation

    use omp_lib

    integer, intent(in) :: &
         n_iter, n_atoms, & ! Number of atoms
         n_ana, &
         n_conv

    logical, intent(in) :: &
         twodim

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         fracA, fracB, & ! Fraction of A and B to use in optimisation
         tol, &
         zdist ! Initial guess of separation distance

    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Centered position of the atoms

    double precision, intent(in), dimension(3,3) :: &
         cell, &  ! cell
         icell ! inverse of cell

    double precision, dimension(3,size(Bpos,2)) :: &
         tBpos

    double precision, dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
         Apos_mapped, & ! position matrix
         Apos_mapped_prev, &
         Bpos_opt

    double precision, intent(out), dimension(3) :: &
         vec, &     ! Translation vector
         angles     ! Angles of rotation

    double precision, dimension(3) :: &
         vec_local, &     ! Translation vector
         angles_local

    double precision :: &
         dist_cur, &
         mul_vec, &
         diag

    double precision, allocatable, dimension(:) :: &
         dist_min

    double precision, allocatable, dimension(:,:) :: &
         angles_min, &
         vec_min

    integer :: &
         j,k, & ! Iterator
         n_threads, thread, &
         pos

    character*(*), intent(in) :: &
         pot, &
         outdir

    character*200 :: &
         cont_filename

    logical :: &
         file_exists
    
    double precision, intent(in) :: &
         param

    if(twodim) then
       diag = 0
       diag = max(norm(cell(:,1) + cell(:,2)),diag)
       diag = max(norm(-cell(:,1) + cell(:,2)),diag)
       
       mul_vec = diag*2/sqrt(2.0d0) !Only in 2D sqrt(3) in 3D
    else
       diag = 0
       diag = max(norm(cell(:,1) + cell(:,2) + cell(:,3)),diag)
       diag = max(norm(-cell(:,1) + cell(:,2) + cell(:,3)),diag)
       diag = max(norm(cell(:,1) - cell(:,2) + cell(:,3)),diag)
       diag = max(norm(-cell(:,1) - cell(:,2) + cell(:,3)),diag)
       
       mul_vec = diag*2/sqrt(3.0d0)
    endif

    !$omp parallel default(private) shared(dist_min, &
    !$omp angles_min, vec_min, angles, vec) &
    !$omp firstprivate(zdist, pot, param, twodim, n_iter, mul_vec, cell, fracA, fracB, tol, &
    !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms, outdir)

    call init_random_seed()

    n_threads = OMP_get_num_threads()

    !$omp single
    allocate(dist_min(n_threads), &
         angles_min(3,n_threads), &
         vec_min(3,n_threads))
    !$omp end single

    thread = OMP_get_thread_num() + 1

    write(cont_filename, "(A, A, I0, A)") outdir, "/continue_", thread, ".dat"
    
    inquire(file=trim(cont_filename), exist=file_exists)
    
    if (file_exists) then
    
       open(20+thread, file = trim(cont_filename), access="stream")

       read(20+thread) angles_min(:,thread), vec_min(:,thread), dist_min(thread)

       close(20+thread)
  
    else

       angles_min(:,thread) = 0.0d0
       vec_min(:,thread) = 0.0d0
       dist_min(thread) = sum(sqrt(sum((Apos - Bpos)**2,1)))

    endif
    
    open(20+thread, file=trim(cont_filename), status='replace', access="stream")

    write(20+thread) angles_min(:,thread), vec_min(:,thread), dist_min(thread)
    flush(20+thread)
    
    !$omp do
    do j=1, n_iter

       call random_number(angles_local)
       call random_number(vec_local)

       angles_local = angles_local*2*pi

       vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)

       if(twodim) then
          vec_local(3) = 0.0d0 ! 2D only
          angles_local(2:3) = 0.0d0
       endif

       vec_local = vec_local * mul_vec

       vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

       if (twodim) then
          vec_local(3) = zdist
       endif
       
       write(13,"(A, I4, A, I6)") "New initial step on thread", thread, ", iteration", j
       flush(13)

       Apos_mapped = 1.0d0
       Apos_mapped_prev = 0.0d0
       Bpos_opt = 0.0d0
       k=1

       do while ( k <= n_conv .and. any(abs(Apos_mapped - Apos_mapped_prev) > 1.0d-12))
          k=k+1

          Apos_mapped_prev = Apos_mapped

          tBpos = free_trans(Bpos,rot_mat(angles_local),vec_local)

          call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms, pot, param)

          call analytical_gd_rot(twodim,angles_local, vec_local, Apos_mapped, Bpos_opt, &
               n_ana, rate1, rate2, tol, pot, param)

       enddo
       
       dist_cur = distance(Apos_mapped, Bpos_opt, rot_mat(angles_local), vec_local, pot, param) 
       
       write(13,"(A, I4, A, I6, A, F8.3)") &
            "Opt dist found for thread", thread,", iteration", j,":", dist_cur
       flush(13)

       if (dist_cur < dist_min(thread)) then
          dist_min(thread) = dist_cur
          angles_min(:,thread) = angles_local
          vec_min(:,thread) = vec_local
          write(20+thread, pos=1) angles_min(:,thread), vec_min(:,thread), dist_min(thread)
          flush(20+thread)
       endif

    enddo
    !$omp end do

    close(20+thread, status='delete')
    
    !$omp barrier
    !$omp single

    pos = minloc(dist_min, 1)

    angles = angles_min(:,pos)
    vec = vec_min(:,pos)

    deallocate(dist_min, &
         angles_min,&
         vec_min)
    !$omp end single 

    !$omp end parallel

  end subroutine gradient_descent_explore

  subroutine gradient_descent_explore_free(twodim, tmat,vec,stats,tmats,vecs,Apos, Bpos, cell, icell, &
       fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, vecrep_in, rate1, rate2, tol, max_vol, &
       zdist, pot, param)
    ! New Gradient Descent Random

    ! Randomly prepare the structures and then minimize the distance with respect to a transformation matrix

    use omp_lib

    logical, intent(in) :: &
         twodim

    integer, intent(in) :: &
         n_iter, n_atoms, & ! Number of atoms
         n_ana, &
         n_conv

    double precision, intent(in) :: &
         rate1, & ! Rate for tmat
         rate2, & ! Rate for disp
         fracA, fracB, & ! Fraction of A and B to use in optimisation
         tol, &
         max_vol, &
         zdist

    integer, intent(in) :: &
         vecrep_in ! Number of repetition of the random vector for 2D per rep of the rotation

    integer :: &
         vecrep ! Number of repetition of the random vector per rep. of the rotation

    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Centered position of the atoms

    double precision, intent(in), dimension(3,3) :: &
         cell, &  ! cell
         icell ! inverse of cell

    double precision, dimension(3,size(Bpos,2)) :: &
         tBpos

    double precision, dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
         Apos_mapped, & ! position matrix
         Apos_mapped_prev, &
         Bpos_opt
    
    double precision, intent(out), dimension(3) :: &
         vec             ! Translation vector

    double precision, intent(out), dimension(n_iter,3) :: &
         stats, &
         vecs

    double precision, intent(out), dimension(n_iter,3,3) :: &
         tmats
    double precision, dimension(3) :: &
         vec_local, &       ! Translation vector
         vec_rot_local, &
         angles_local

    double precision, intent(out), dimension(3,3) :: &
         tmat    ! Transformation matrix

    double precision, dimension(3,3) :: &
         tmat_local, &
         tmat_local_reset

    double precision, dimension(2,2) :: &
         tt
    double precision :: &
         dist_cur, &
         dist_cur_rot, &
         mul_vec, &
         diag, &
         dt_tmat, &
         eig1, &
         eig2

    double precision, allocatable, dimension(:) :: &
         dist_min

    double precision, allocatable, dimension(:,:,:) :: &
         tmat_min

    double precision, allocatable, dimension(:,:) :: &
         vec_min

    integer :: &
         j,k,vl, & ! Iterator
         n_threads, thread, &
         pos

    logical :: &
         vol

    character*(*), intent(in) :: &
         pot

    double precision, intent(in) :: &
         param

    if(twodim) then
       diag = 0
       diag = max(norm(cell(:,1) + cell(:,2)),diag)
       diag = max(norm(-cell(:,1) + cell(:,2)),diag)
       
       mul_vec = diag*2/sqrt(2.0d0)
       
       vecrep = vecrep_in
    else
       diag = 0
       diag = max(norm(cell(:,1) + cell(:,2) + cell(:,3)),diag)
       diag = max(norm(-cell(:,1) + cell(:,2) + cell(:,3)),diag)
       diag = max(norm(cell(:,1) - cell(:,2) + cell(:,3)),diag)
       diag = max(norm(-cell(:,1) - cell(:,2) + cell(:,3)),diag)
       
       mul_vec = diag*2/sqrt(3.0d0)

       vecrep = 1
    endif

    !$omp parallel default(private) shared(dist_min, &
    !$omp tmat_min, vec_min, vec, tmat, stats, vecs, tmats) &
    !$omp firstprivate(zdist, param, pot, twodim, vecrep, n_iter, mul_vec, max_vol, cell, &
    !$omp fracA, fracB, tol, &
    !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms)

    call init_random_seed()

    n_threads = OMP_get_num_threads()

    !$omp single
    allocate(dist_min(n_threads), &
         tmat_min(3,3,n_threads), &
         vec_min(3,n_threads))
    !$omp end single

    thread = OMP_get_thread_num() + 1

    dist_min(thread) = 2*sum(sqrt(sum((Apos - Bpos)**2,1))) ! A large number

    !$omp do
    do j=1, n_iter
       
       vol = .false.
       do while (.not. vol)

          call random_number(tmat_local)
          call random_number(vec_local)
          
          if (twodim) then

             angles_local(1) = tmat_local(1,2)*2*pi ! Using tmat local as a random number
             angles_local(2:3) = 0.0d0

             tmat_local = (1.0d0 + 0.1d0 * max_vol * (2 * tmat_local(1,1) - 1)) ** (1.0d0/2.0d0) * eye()
             tmat_local(:,3) = 0.0d0
             tmat_local(3,:) = 0.0d0
             tmat_local(3,3) = 1.0d0

             tmat_local_reset = matmul(rot_mat(angles_local), tmat_local)
             
          else
             vol = .true. 
             
             tmat_local = (1.0d0 + max_vol) * (tmat_local - 0.5d0) ! TDOO: Is this optimal? 

             dt_tmat = det(tmat_local,3)

             ! The volume of the cell is made smaller then the ratio or max volume
             ! It seems like smaller volume lead to better mapping when one to one is imposed

             tmat_local = dt_tmat / abs(dt_tmat) * tmat_local * (abs(modulo(abs(dt_tmat), max_vol) / dt_tmat))**(1.0d0/3.0d0)

          endif

          do vl=1, vecrep
             call random_number(vec_local)

             vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)
             if (twodim) then
                tmat_local = tmat_local_reset
                vec_local(3) = 0.0d0
             endif

             vec_local = vec_local * mul_vec

             vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

             if (twodim) then
                vec_local(3) = zdist
             endif

             write(13,"(A, I4, A, I6, A, I6)") "New initial step on thread", thread, &
                  ", iteration", j, ".", vl
             flush(13)

             Apos_mapped = 1.0d0
             Apos_mapped_prev = 0.0d0
             k=1
             
             do while ( k <= n_conv .and. any(abs(Apos_mapped - Apos_mapped_prev) > 1.0d-12))
                k=k+1

                Apos_mapped_prev = Apos_mapped

                tBpos = free_trans(Bpos,tmat_local,vec_local)
                
                call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
                     fracA, fracB, atoms, n_atoms, pot, param)
                
                call analytical_gd_free(twodim, tmat_local, vec_local, Apos_mapped, Bpos_opt, &
                     n_ana, rate1, rate2, tol, pot, param)
                
             enddo

             angles_local = 0
             vec_rot_local = vec_local

             if (rotmin) then
                call analytical_gd_rot(twodim, angles_local, vec_rot_local, tmat_local, eye(), &
                     100000, 1.0d0, 1.0d0, tol, "Euclidean", 0.0d0)
             else
                call analytical_gd_rot(twodim, angles_local, vec_rot_local, Apos_mapped, Bpos_opt, &
                     n_ana*1000, rate1, rate2, tol, pot, param)
             endif

             vec_rot_local = vec_local
             
             dist_cur = distance(Apos_mapped, Bpos_opt, tmat_local, vec_local, pot, param)
             dist_cur_rot = distance(Apos_mapped, Bpos_opt, rot_mat(angles_local), vec_rot_local, pot, param)
             
             if (twodim) then
                tt = matmul(transpose(tmat_local(1:2,1:2)), tmat_local(1:2,1:2))
                eig1 = (tt(1,1)+tt(2,2))/2 + sqrt(((tt(1,1)-tt(2,2))/2)**2 +&
                     tt(1,2)*tt(2,1))
                eig2 = (tt(1,1)+tt(2,2))/2 - sqrt(((tt(1,1)-tt(2,2))/2)**2 +&
                     tt(1,2)*tt(2,1))
                if (abs(det(tmat_local,3) - 1.0d0) < max_vol .and. &
                     abs(sqrt(eig1) - 1) < max_vol .and. &
                     abs(sqrt(eig2) - 1) < max_vol) then
             write(13,"(A, I4, A, I6, A, I6, A, F8.3, A, F8.3)") &
                  "Opt dist found for thread", thread,", iteration", j,".",vl,":", &
                  dist_cur, " ", dist_cur_rot
             flush(13)

                   if (vl == 1 .or. dist_cur < stats(j,3)) then !TMP vec
                      stats(j,1) = angles_local(1)
                      stats(j,2) = det(tmat_local,3)
                      stats(j,3) = dist_cur
                      tmats(j,:,:) = tmat_local
                      vecs(j,:) = vec_local
                      vol = .true.
                   endif !TMP vec

                   if (dist_cur < dist_min(thread)) then
                      dist_min(thread) = dist_cur
                      tmat_min(:,:,thread) = tmat_local
                      vec_min(:,thread) = vec_local
                   endif
                endif
             else
                write(13,"(A, I4, A, I6, A, I6, A, F8.3, A, F8.3)") &
                     "Opt dist found for thread", thread,", iteration", j,".",vl,":", &
                     dist_cur, " ", dist_cur_rot
                flush(13)
                if (dist_cur_rot + dist_cur < dist_min(thread)) then ! TMP TESTING
                   ! if (dist_cur < dist_min(thread)) then ! TMP TESTING
                   dist_min(thread) = dist_cur_rot + dist_cur
                   tmat_min(:,:,thread) = tmat_local
                   vec_min(:,thread) = vec_local
                endif
             endif

          enddo
       enddo
    enddo
    !$omp end do
    !$omp barrier
    !$omp single

    pos = minloc(dist_min, 1)

    tmat = tmat_min(:,:,pos)
    vec = vec_min(:,pos)
    
    deallocate(dist_min, &
         tmat_min,&
         vec_min)
    !$omp end single 

    !$omp end parallel

  end subroutine gradient_descent_explore_free

  subroutine find_peaks(array, n_array, tol_prom, n, peaks_out, idx_out, prom_out) 

    double precision, intent(in), dimension(:) :: &
         array

    integer, intent(in) :: &
         n_array

    double precision, intent(in) :: &
         tol_prom

    double precision, intent(out), dimension(size(array)) :: &
         peaks_out, prom_out

    double precision, dimension(size(array)) :: &
         peaks, peaks_copy, prom

    integer, intent(out), dimension(size(array)) :: &
         idx_out

    integer, dimension(size(array)) :: &
         idx, idx_copy, idxx

    integer :: &
         i,j,n,id, &
         start, &
         finish

    integer, parameter :: &
         peak_size = 10

    n=0
    do i = 1, n_array
       start = 1 + modulo(i-peak_size-1, n_array)
       finish = 1 + modulo(i+peak_size-1, n_array)
       if (start > finish) then
           if (all(array(start:n_array) >= array(i)) .and. all(array(1:finish) >= array(i))) then
              n=n+1
              idx(n) = i
              peaks(n) = array(i)
           endif
        else
           if (all(array(start:finish) >= array(i))) then
              n=n+1
              idx(n) = i
              peaks(n) = array(i)
           endif
        endif

    enddo

    idxx(1:n) = sort(peaks(1:n))

    peaks_copy = peaks
    idx_copy = idx


    do i=1,n
       peaks(i) = peaks_copy(idxx(i))
       idx(i) = idx_copy(idxx(i))
    enddo

    prom(1) = 1

    do i=2,n
       id = idx(minloc(abs(idx(1:i-1)-idx(i)), 1))
       if (id < idx(i)) then
          prom(i) = (maxval(array(id:idx(i))) - peaks(i)) / (maxval(array) - peaks(1))

       else
          prom(i) = (maxval(array(idx(i):id)) - peaks(i)) / (maxval(array) - peaks(1))
       endif
    enddo

    j = 0
    do i=1,n
       if (prom(i) > tol_prom) then
          j = j + 1
          peaks_out(j) = peaks(i)
          idx_out(j) = idx(i)
          prom_out(j) = prom(i)
       endif
    enddo

    n = j

  end subroutine find_peaks

  subroutine adjust(twodim, Apos, Bpos, Apos_mapped, Bpos_opt, & ! Output
       tmat, vec, & ! Output
       na, nb, check, &
       fracA, fracB, &
       atoms, n_atoms, rate1, rate2, &
       n_ana, n_conv, n_adjust, &
       tol, zdist, &
       pot, param)

    logical, intent(in) :: &
         twodim

    character*(*), intent(in) :: &
         pot

    double precision, intent(in) :: &
         param, &
         zdist

    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2, &     ! Rate of the gradient descent for displacement
         tol

    integer, intent(in):: &
         n_ana, &
         n_conv, &
         n_adjust

    double precision, intent(in), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(in), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, dimension(3,nb) :: &
         tBpos ! Position of the atoms

    double precision, intent(in) :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, intent(inout), dimension(3) :: &
         vec        ! Translation vecto

    double precision, dimension(3) :: &
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         angles  ! Rotation angles

    double precision, intent(inout), dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, intent(inout), dimension(:,:) :: &
         Apos_mapped, &
         Bpos_opt

    double precision, dimension(size(Apos_mapped,1),size(Apos_mapped,2)) :: &
         Apos_mapped_prev, Apos_mapped_2prev, &
         tBpos_opt

    integer :: &
         i, j

    logical, intent(in) :: &
         check

    vec_rot = vec
    Apos_mapped = 1.0d0
    Apos_mapped_prev = 0.0d0
    i=0

    if (n_adjust == 0) then
       tBpos = free_trans(Bpos,tmat,vec)
       call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
            fracA, fracB, atoms, n_atoms, pot, param)
    endif
    
    do while (any(abs(Apos_mapped - Apos_mapped_2prev) > 1.0d-12) .and. i < n_adjust)
       i = i + 1
       
       Apos_mapped_2prev = Apos_mapped
       Apos_mapped_prev = 0.0d0
       j=0
       do while (any(abs(Apos_mapped - Apos_mapped_prev) > 1.0d-12) .and. j <= n_conv)
          j = j + 1

          Apos_mapped_prev = Apos_mapped

          tBpos = free_trans(Bpos,tmat,vec)

          if (check) then

             call no_mapping(Apos_mapped, Bpos_opt, Apos, Bpos, &
                  fracB, atoms, n_atoms)

          else

             call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
                  fracA, fracB, atoms, n_atoms, pot, param)

          endif

          angles = 0
          tBpos_opt = free_trans(Bpos_opt,tmat,(/0.0d0,0.0d0,0.0d0/))



          call analytical_gd_rot(twodim, angles, vec, Apos_mapped, tBpos_opt, &
               n_ana*100, rate1, rate2, tol, pot, param)

          tmat = matmul(rot_mat(angles),tmat)

       enddo

       angles = 0.0d0

       ! This step is just to get the "unstretched distance"
       if (twodim .and. rotmin) then
          vec_rot = 0.0d0
          call analytical_gd_rot(twodim, angles, vec_rot, tmat, eye(), &
               100000, 1.0d0, 1.0d0, &
               tol, "Euclidean", 0.0d0)
          vec_rot = 0.0d0
          vec_rot(3) = zdist
       else
          vec_rot = vec
       endif

       
       call analytical_gd_rot(twodim, angles, vec_rot, Apos_mapped, Bpos_opt, &
               n_ana*1000, rate1, rate2, tol, pot, param)

       write(13,*) "-->", i

       write(13,*) "Stretched distance:", &
            distance(Apos_mapped, Bpos_opt, tmat, vec, pot, param)

       write(13,*) "Unstretched distance:", &
            distance(Apos_mapped, Bpos_opt, rot_mat(angles), vec_rot, pot, param)
       flush(13)

       call analytical_gd_free(twodim, tmat, vec, Apos_mapped, Bpos_opt, n_ana*1000, &
            rate1, rate2, tol, pot, param)

    enddo

  end subroutine adjust

  subroutine remove_slanting(Apos_mapped, Bpos_opt, & ! Output
       tmat, vec, & ! Output
       rate1, rate2, &
       n_ana, n_adjust, &
       tol)

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2, &     ! Rate of the gradient descent for displacement
         tol

    integer, intent(in):: &
         n_ana, &
         n_adjust

    double precision, intent(inout), dimension(3) :: &
         vec        ! Translation vecto

    double precision, dimension(3) :: &
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         angles  ! Rotation angles

    double precision, intent(inout), dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, dimension(3,3) :: &
         slant ! Transformation matrix

    double precision, intent(inout), dimension(:,:) :: &
         Apos_mapped, &
         Bpos_opt

    double precision, dimension(size(Apos_mapped,1),size(Apos_mapped,2)) :: &
         tBpos_opt, rBpos_opt

    integer :: &
         i

    rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)

    do i=1, n_adjust

       write(13,*) "-->", i

       tBpos_opt = free_trans(Bpos_opt, tmat, (/0.0d0,0.0d0,0.0d0/))

       angles = 0.0d0

       call analytical_gd_rot(.false., angles, vec, rBpos_opt, tBpos_opt, &
            n_ana*1000, rate1, rate2, tol, "Euclidean", 0.0d0)

       tBpos_opt = free_trans(tBpos_opt, rot_mat(angles), (/0.0d0,0.0d0,0.0d0/))

       tmat = matmul(rot_mat(angles), tmat) 

       call analytical_gd_slant(slant, vec, rBpos_opt, tBpos_opt, n_ana*1000, rate1, rate2, &
            tol)

       tmat = matmul(slant, tmat)

       ! print*, "Slant vol.", det(slant,3)
       ! print*, "Volume", det(tmat,3)

    enddo

    write(13,*) "Stretched distance after slanting:", &
         distance(Apos_mapped, Bpos_opt, tmat, vec, "Euclidean", 0.0d0)
    flush(13)

  end subroutine remove_slanting

  subroutine classification(minimize_in, twodim, Apos, Bpos, Apos_mapped, Bpos_opt, & ! Output
       tmat, vec, & ! Output
       classes_list, &
       fracA, fracB, n_frac, &
       na, nb, n_class_in, &
       remap, check, &
       atoms, n_atoms, rate1, rate2, &
       n_ana, n_out, &
       tol, tol_class, tol_std, &
       init_class, &
       pot, param)

    logical, intent(in) :: &
         twodim
    
    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms, & ! Number of types of atoms per cell
         n_frac

    double precision, intent(in) :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2, &     ! Rate of the gradient descent for displacement
         tol, tol_class, &
         init_class, & ! Initial class separation criteria
         tol_std

    integer, intent(in):: &
         n_ana, &
         n_out, &
         n_class_in

    double precision, intent(in), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(in), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, dimension(3,nb) :: &
         tBpos ! Position of the atoms

    double precision, intent(in) :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, intent(inout), dimension(3) :: &
         vec        ! Translation vecto


    double precision, intent(inout), dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, intent(inout), dimension(3,n_out) :: &
         Apos_mapped, &
         Bpos_opt

    double precision, dimension(3,n_out) :: &
         Apos_mapped_prev, &
         disps

    integer :: &
         i, j, &
         id, &
         n_B, n_tot

    logical, intent(in) :: &
         remap, &
         check, &
         minimize_in

    logical :: &
         minimize

    integer :: &
         n, n_class

    double precision :: &
         tol_adjust, &
         std, std_prev

    integer, dimension(n_out) :: &
         classes_list, &
         classes_list_prev, &
         n_classes_trail

    integer, allocatable, dimension(:) :: &
         n_classes

    character*(*), intent(in) :: &
         pot
    double precision, intent(in) :: &
         param

    n_class = n_class_in
    minimize = minimize_in

    tol_adjust = init_class
    std = 1.0d0
    classes_list = 0
    classes_list_prev = 1
    j=0
    if (n_class == 0) then
       minimize = .false.
       n_class = 1
       tol_adjust = tol_class
    endif
    
    do while ( std > tol_class .and. j < n_class)
       j = j + 1

       write(13,*) "-->", j

       if ((remap .and. j/=1) .or. .not. minimize) then

          
          Apos_mapped_prev = Apos_mapped

          tBpos = free_trans(Bpos,tmat,vec)

          if (check) then

             call no_mapping(Apos_mapped, Bpos_opt, Apos, Bpos, &
                  fracB, atoms, n_atoms)

          else

             call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
                  fracA, fracB, atoms, n_atoms, pot, param)

          endif

          if (any(abs(Apos_mapped_prev - Apos_mapped) < 1.0d-12)) then
             write(13,*) "REMAP!"
          endif

       endif

       classes_list_prev = classes_list
       std_prev = std

       disps = Apos_mapped-free_trans(Bpos_opt,tmat,vec)

       n_B = 0 
       n_tot = 0
       id = 0
       do i=0,n_atoms-1

          n_B = n_frac*atoms(i+1)

          call classify(n, n_classes_trail(n_tot+1:n_tot + n_B), classes_list(id+1:id + n_B), &
               disps(:, id+1:id + n_B), tol_adjust)

          classes_list(id+1:id + n_B) =  n_tot + classes_list(id+1:id + n_B)

          id = id + n_B
          n_tot = n_tot + n

       enddo

       allocate(n_classes(n_tot))

       do i=1,n_tot
          n_classes(i) = n_classes_trail(i)
       enddo

       if (minimize) then
       
          call analytical_gd_std(twodim, std, tmat, vec, Apos_mapped, Bpos_opt, &
               n_classes, classes_list, n_ana*1000, rate1, rate2, tol_std)

       else

          std = 0.0d0
          
       endif
       
       deallocate(n_classes)

       write(13,*) "Tolerance for classification:", tol_adjust
       write(13,*) "Final standard deviation:", sqrt(std/dble(size(Apos_mapped,2)-1))
       write(13,*) "Number of classes:", n_tot
       write(13,*) "Stretched distance:", distance(Apos_mapped, Bpos_opt, tmat, vec, pot, param)
       flush(13)

       ! tol_adjust = min(3.0d0*sqrt(std/dble(size(Apos_mapped,2)-1)), tol_adjust/2.0d0) !3 sigma 99%


       !       if (all(classes_list == classes_list_prev) .and. std_prev - std < tol) then
       if (std_prev - std < tol) then
          tol_adjust = tol_adjust/2.0d0
       endif


    enddo

  end subroutine classification

  subroutine optimize_vec(vec_out, Apos, Bpos, n, vec, filename)

    integer, intent(in) :: &
         n           ! Number of atoms
    
    double precision, intent(in), dimension(3,n) :: &
         Apos, &
         Bpos        ! Bpos ordered according to the mapping
    
    double precision, intent(in), dimension(3) :: &         
         vec         ! Translation vector

    double precision, intent(out), dimension(3) :: &         
         vec_out     ! Translation vector

    character*200, intent(in) :: &
         filename   ! Number of atoms

    double precision :: &
         tol, rate2

    integer :: &
         n_ana

    logical :: &
         exist

    namelist /input/ &
         tol, &
         rate2, &
         n_ana

    ! Namelist default values
    tol = 1.0d-4
    rate2 = 1.0d-3
    n_ana = 300

    inquire(file = filename, exist=exist)

    if (exist) then
       open (unit = 11, file = filename, status = 'OLD')
       read (11, input)
       close (11)
    endif

    vec_out = vec

    call analytical_gd_vec(.false., eye(), vec_out, &
         Apos, Bpos, n_ana*1000, rate2, &
         tol, "Euclidean", 0.0d0)
  
  end subroutine optimize_vec
  
  subroutine fastoptimization(Apos_out, Bpos_out, Bpos_out_stretch, & ! Output
       n_out, n_A, classes_list_out, & ! Output
       tmat, dmin, vec, & ! Output
       Apos, na, Bpos, nb, &
       Acell, iAcell, &
       atoms, n_atoms, &
       filename, outdir)

    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2     ! Rate of the gradient descent for displacement

    integer :: &
         n_iter, &   ! Number of iteration of the gradient descent
         n_ana, &
         n_conv, &
         n_adjust, &
         n_class

    character*200, intent(in) :: &
         filename, &
         outdir

    double precision, intent(inout), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(inout), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, intent(out), dimension(3,na) :: &
         Apos_out ! Position of the atoms

    double precision, intent(out), dimension(3,nb) :: &
         Bpos_out, & ! Position of the atoms
         Bpos_out_stretch

    integer, intent(out) :: &
         n_out, &
         n_A

    double precision, intent(in), dimension(3,3) :: &
         Acell, & ! Unit cell of A
         iAcell

    double precision :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, dimension(3), intent(out) :: &
         dmin

    double precision, intent(out), dimension(3) :: &
         vec        ! Translation vecto

    double precision, dimension(3) :: &
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         angles  ! Rotation angles

    double precision :: &
         tol, &
         tol_class, &
         tol_std, &
         init_class, &
         max_vol, &
         dmin_half

    double precision, intent(out), dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, allocatable, dimension(:,:) :: &
         Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, &
         tBpos_opt, &
         rBpos_opt, &
         Bpos_opt_stretch  ! position matrix

    integer, intent(out), dimension(size(Apos,2)) :: &         
         classes_list_out

    integer, allocatable, dimension(:) :: &
         classes_list, &
         classes_list_prev, &
         n_classes_trail

    integer :: &
         n_frac

    double precision, allocatable, dimension(:,:) :: &
         stats, &
         vecs
    
    double precision, allocatable, dimension(:,:,:) :: &
         tmats

    logical :: &
         slanting, &
         remap, &
         free, &
         check, &
         usebest, &
         exist

    character*200 :: &
         savebest, &
         progressfile

    double precision, parameter, dimension(3) :: &
         zeros = (/0.0d0, 0.0d0, 0.0d0/)


    namelist /input/ &
         fracA, fracB, &
         tol, tol_std, &
         tol_class, &
         rate1, rate2, &
         slanting, &
         n_iter, n_ana, &
         n_conv, n_adjust, &
         n_class, &
         max_vol, free, &
         check, &
         savebest, &
         usebest, &
         remap, &
         init_class

    ! Namelist default values
    init_class = 1.0d0
    tol = 1.0d-4
    tol_std = tol*1.0d-3
    tol_class = 1.0d-3
    rate1 = 1.0d-3
    rate2 = 1.0d-3
    slanting = .false.
    fracA = 0.09d0
    fracB = 0.3d0
    n_iter = 1000
    n_ana = 300
    n_conv = 5
    n_class = 30
    n_adjust = 10
    max_vol = 4.0d0
    free = .false.
    check = .false.
    savebest = trim(outdir)//"/best.dat"
    usebest = .false.
    remap = .true.


    progressfile = trim(outdir)//"/progress.txt"
    open(13, file = trim(progressfile), status='replace')

    inquire(file = filename, exist=exist)

    if (exist) then
       open (unit = 11, file = filename, status = 'OLD')
       read (11, input)
       close (11)
    endif


    n_frac = int(fracB*size(Apos,2)/sum(atoms))
    n_A = int(fracA*size(Apos,2)/sum(atoms))
    n_out = n_frac * sum(atoms)
    Apos_out = 0.0d0
    Bpos_out = 0.0d0
    angles = 0.0d0

    allocate(Apos_mapped(3,n_out), Bpos_opt(3,n_out), &
         Apos_mapped_prev(3,n_out), &
         tBpos_opt(3,n_out), &
         rBpos_opt(3,n_out), &
         Bpos_opt_stretch(3,n_out), &
         classes_list(n_out), &
         classes_list_prev(n_out), &
         n_classes_trail(n_out), &
         stats(n_iter,3), &
         vecs(n_iter,3), &
         tmats(n_iter,3,3))


    ! Center both cells at the geometric center
    call center(Bpos,nb)
    call center(Apos,na)


    if (check) then

       call init_random_seed()
       call random_number(angles)

       angles = angles*2*pi

       tmat = rot_mat(angles)

       !tmat = eye()

       vec = 0

    else

       if (usebest) then

          open(12, file = trim(savebest))
          read(12,*) tmat, vec
          close(12)
          
       else
          if (free) then

             call gradient_descent_explore_free(.false., tmat, vec, stats, tmats, vecs, Apos, Bpos, Acell, iAcell, &
                  fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, 1, rate1, rate2, tol, &
                  max_vol, 0.0d0, "Euclidean", 0.0d0)

          else

             call gradient_descent_explore(.false., angles, vec, Apos, Bpos, Acell, iAcell, &
                  fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, rate1, rate2, tol, &
                  0.0d0, "Euclidean", 0.0d0, trim(outdir))

             tmat = rot_mat(angles)
             
          endif
          
          if (trim(savebest)/="") then
             open(12, file = trim(savebest))
             write(12,*) tmat, vec
             close(12)

          endif

       endif

    endif

    write(13,*) "/======== Stretching Adjustment ========\\"

    flush(13)

    call adjust(.false., Apos, Bpos, Apos_mapped, Bpos_opt, &
         tmat, vec, &
         na, nb, check, &
         fracA, fracB, &
         atoms, n_atoms, rate1, rate2, &
         n_ana, n_conv, n_adjust, &
         tol, 0.0d0, "Euclidean", 0.0d0)

    if (slanting) then

       write(13,*) "/======== Deslanting ========\\"
       flush(13)

       call remove_slanting(Apos_mapped, Bpos_opt, & ! Output
            tmat, vec, & ! Output
            rate1, rate2, &
            n_ana, n_adjust, &
            tol)
    endif


    write(13,*) "/======== Classification ========\\"

    call classification(.true., .false., Apos, Bpos, Apos_mapped, Bpos_opt, & ! Output
         tmat, vec, & ! Output
         classes_list, &
         fracA, fracB, n_frac, &
         na, nb, n_class, &
         remap, check, &
         atoms, n_atoms, rate1, rate2, &
         n_ana, n_out, &
         tol, tol_class, tol_std, &
         init_class, &
         "Euclidean", 0.0d0)
    write(13,*) "Final stretched distance:", distance(Apos_mapped, Bpos_opt,tmat,vec,"Euclidean",0.0d0)

    ! Reshift after classification
    call analytical_gd_vec(.false.,tmat, vec, &
         Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
         tol, "Euclidean", 0.0d0)

    vec_rot = 0.0d0

    ! This is only to get closer since angles start from scratch.
    ! It needs more strict convergence since this is only three points
    call analytical_gd_rot(.false., angles, vec_rot, tmat, eye(), &
         100000, 1.0d0, 1.0d0, tol/100, "Euclidean", 0.0d0)
    
    vec_rot = vec
    
    call analytical_gd_rot(.false., angles, vec_rot, Apos_mapped, Bpos_opt, &
         n_ana*1000, rate1, rate2, tol, "Euclidean", 0.0d0)
    
    rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)

    write(13,*) "Final rmat"
    write(13,"(3(F7.3,X))") rot_mat(angles)

    classes_list_out = 0
    classes_list_out(1:n_out) = classes_list

    write(13,*) "Final tmat"
    write(13,"(3(F7.3,X))") tmat

    write(13,*) "Final stretched distance:", distance(Apos_mapped, Bpos_opt,tmat,vec,"Euclidean",0.0d0)
    
    Bpos_opt_stretch = free_trans(Bpos_opt,tmat,vec)

    Bpos_out(:,1:n_out) = rBpos_opt
    Bpos_out_stretch(:,1:n_out) = Bpos_opt_stretch
    Apos_out(:,1:n_out) = Apos_mapped

    dmin(1) = distance(Apos_mapped, rBpos_opt,eye(),zeros,"Euclidean",0.0d0)
    
    dmin_half = distance(Apos_mapped(:, 1:int(n_out/2)), rBpos_opt(:, 1:int(n_out/2)), eye(), zeros, "Euclidean", 0.0d0)
       
    dmin(2) = (dmin(1)/n_out - dmin_half/int(n_out/2))/(n_out**(1.0/3.0) - int(n_out/2)**(1.0/3.0))

    dmin(3) = dmin(1)/n_out - dmin(2)*n_out**(1.0/3.0)

    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(13,"(10(F5.3,X))") mat

    close(13)

    deallocate(Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, &
         tBpos_opt, &
         rBpos_opt, &
         Bpos_opt_stretch, &
         classes_list, &
         classes_list_prev, &
         n_classes_trail, &
         stats, &
         vecs, &
         tmats)

  end subroutine fastoptimization

    subroutine fixed_tmat(Apos_out, Bpos_out, Bpos_out_stretch, & ! Output
       n_out, n_A, classes_list_out, & ! Output
       tmat, dmin, vec, & ! Output
       Apos, na, Bpos, nb, &
       tmat_in, vec_in, &
       atoms, n_atoms, &
       filename, outdir)

    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2     ! Rate of the gradient descent for displacement

    integer :: &
         n_iter, &   ! Number of iteration of the gradient descent
         n_ana, &
         n_conv, &
         n_adjust, &
         n_class

    character*200, intent(in) :: &
         filename, &
         outdir

    double precision, intent(inout), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(inout), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, intent(out), dimension(3,na) :: &
         Apos_out ! Position of the atoms

    double precision, intent(out), dimension(3,nb) :: &
         Bpos_out, & ! Position of the atoms
         Bpos_out_stretch

    integer, intent(out) :: &
         n_out, &
         n_A

    double precision :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, dimension(3), intent(out) :: &
         dmin

    double precision, intent(out), dimension(3) :: &
         vec        ! Translation vector

    double precision, intent(in), dimension(3) :: &
         vec_in        ! Translation vector

    double precision, dimension(3) :: &
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         angles  ! Rotation angles

    double precision :: &
         tol, &
         tol_class, &
         tol_std, &
         init_class, &
         max_vol, &
         dmin_half

    double precision, intent(in), dimension(3,3) :: &
         tmat_in ! Transformation matrix
    
    double precision, intent(out), dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, allocatable, dimension(:,:) :: &
         Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, &
         tBpos_opt, &
         rBpos_opt, &
         Bpos_opt_stretch  ! position matrix

    integer, intent(out), dimension(size(Apos,2)) :: &         
         classes_list_out

    integer, allocatable, dimension(:) :: &
         classes_list, &
         classes_list_prev, &
         n_classes_trail

    integer :: &
         n_frac

    double precision, allocatable, dimension(:,:) :: &
         stats, &
         vecs
    
    double precision, allocatable, dimension(:,:,:) :: &
         tmats

    logical :: &
         slanting, &
         remap, &
         free, &
         check, &
         usebest, &
         exist

    character*200 :: &
         savebest, &
         progressfile

    double precision, parameter, dimension(3) :: &
         zeros = (/0.0d0, 0.0d0, 0.0d0/)


    namelist /input/ &
         fracA, fracB, &
         tol, tol_std, &
         tol_class, &
         rate1, rate2, &
         slanting, &
         n_iter, n_ana, &
         n_conv, n_adjust, &
         n_class, &
         max_vol, free, &
         check, &
         savebest, &
         usebest, &
         remap, &
         init_class

    ! Namelist default values
    init_class = 1.0d0
    tol = 1.0d-4
    tol_std = tol*1.0d-3
    tol_class = 1.0d-3
    rate1 = 1.0d-3
    rate2 = 1.0d-3
    slanting = .false.
    fracA = 0.09d0
    fracB = 0.3d0
    n_iter = 1000
    n_ana = 300
    n_conv = 5
    n_class = 30
    n_adjust = 10
    max_vol = 4.0d0
    free = .false.
    check = .false.
    savebest = trim(outdir)//"/best.dat"
    usebest = .false.
    remap = .true.


    progressfile = trim(outdir)//"/progress_find_cell.txt"
    open(13, file = trim(progressfile), status='replace')

    inquire(file = filename, exist=exist)

    if (exist) then
       open (unit = 11, file = filename, status = 'OLD')
       read (11, input)
       close (11)
    endif


    n_frac = int(fracB*size(Apos,2)/sum(atoms))
    n_A = int(fracA*size(Apos,2)/sum(atoms))
    n_out = n_frac * sum(atoms)
    Apos_out = 0.0d0
    Bpos_out = 0.0d0
    angles = 0.0d0

    allocate(Apos_mapped(3,n_out), Bpos_opt(3,n_out), &
         Apos_mapped_prev(3,n_out), &
         tBpos_opt(3,n_out), &
         rBpos_opt(3,n_out), &
         Bpos_opt_stretch(3,n_out), &
         classes_list(n_out), &
         classes_list_prev(n_out), &
         n_classes_trail(n_out), &
         stats(n_iter,3), &
         vecs(n_iter,3), &
         tmats(n_iter,3,3))


    ! Center both cells at the geometric center
    call center(Bpos,nb)
    call center(Apos,na)


    tmat = tmat_in
    
    vec = vec_in

    write(13,*) "/======== Classification ========\\"

    Apos_mapped = Apos(:,1:n_out)
    Bpos_opt = Bpos(:,1:n_out)
    
    call classification(.false., .false., Apos, Bpos, Apos_mapped, Bpos_opt, & ! Output
         tmat, vec, & ! Output
         classes_list, &
         fracA, fracB, n_frac, &
         na, nb, n_class, &
         remap, check, &
         atoms, n_atoms, rate1, rate2, &
         n_ana, n_out, &
         tol, tol_class, tol_std, &
         tol_class, &
         "Euclidean", 0.0d0)

    vec_rot = 0.0d0

    call init_random_seed()

    call random_number(angles) ! So that the rotation min doesn't get stuck                                             
    
    ! This is only to get closer since angles start from scratch.
    ! It needs more strict convergence since this is only three points
    call analytical_gd_rot(.false., angles, vec_rot, tmat, eye(), &
         100000, 1.0d0, 1.0d0, tol/100, "Euclidean", 0.0d0)
    
    vec_rot = vec
    
    call analytical_gd_rot(.false., angles, vec_rot, Apos_mapped, Bpos_opt, &
         n_ana*1000, rate1, rate2, tol, "Euclidean", 0.0d0)
    
    rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)

    write(13,*) "Final rmat"
    write(13,"(3(F7.3,X))") rot_mat(angles)

    classes_list_out = 0
    classes_list_out(1:n_out) = classes_list

    write(13,*) "Final tmat"
    write(13,"(3(F7.3,X))") tmat
    
    Bpos_opt_stretch = free_trans(Bpos_opt,tmat,vec)

    Bpos_out(:,1:n_out) = rBpos_opt
    Bpos_out_stretch(:,1:n_out) = Bpos_opt_stretch
    Apos_out(:,1:n_out) = Apos_mapped

    dmin(1) = distance(Apos_mapped, rBpos_opt,eye(),zeros,"Euclidean",0.0d0)
    
    dmin_half = distance(Apos_mapped(:, 1:int(n_out/2)), rBpos_opt(:, 1:int(n_out/2)), eye(), zeros, "Euclidean", 0.0d0)
       
    dmin(2) = (dmin(1)/n_out - dmin_half/int(n_out/2))/(n_out**(1.0/3.0) - int(n_out/2)**(1.0/3.0))

    dmin(3) = dmin(1)/n_out - dmin(2)*n_out**(1.0/3.0)

    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(13,"(10(F5.3,X))") mat

    close(13)

    deallocate(Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, &
         tBpos_opt, &
         rBpos_opt, &
         Bpos_opt_stretch, &
         classes_list, &
         classes_list_prev, &
         n_classes_trail, &
         stats, &
         vecs, &
         tmats)

  end subroutine fixed_tmat
  
  subroutine intoptimization(Apos_out, Bpos_out, Bpos_out_stretch, &
       n_out, n_A, classes_list_out, ttrans, rtrans, dmin, stats, n_peaks, &
       peak_thetas, sym, n_iter, &
       Apos, na, Bpos, nb, &
       Acell, iAcell, &
       atoms, n_atoms, &
       switched, &
       filename, outdir)

    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms, & ! Number of types of atoms per cell
         sym, & ! Max symmetry of the two structures
         n_iter

    double precision :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2     ! Rate of the gradient descent for displacement

    integer :: &
         n_ana, &
         n_conv, &
         n_adjust, &
         n_class

    character*200, intent(in) :: &
         filename, &
         outdir

    double precision, intent(inout), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(inout), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, intent(out), dimension(n_iter,3,na) :: &
         Apos_out ! Position of the atoms

    double precision, intent(out), dimension(n_iter,3,nb) :: &
         Bpos_out, & ! Position of the atoms
         Bpos_out_stretch

    integer, intent(out) :: &
         n_out, &
         n_A, &
         n_peaks

    logical, intent(in) :: &
         switched
    
    double precision, intent(in), dimension(3,3) :: &
         Acell, & ! Unit cell of A
         iAcell

    double precision :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, intent(out), dimension(n_iter,3) :: &
         dmin

    double precision, intent(out), dimension(n_iter,3) :: &
         stats

    double precision, dimension(n_iter,3,3) :: &
         tmats, &
         tmats_ordered

    double precision, dimension(3) :: &
         vec, &     ! Translation vecto
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         angles       ! Rotation axis

    double precision :: &
         tol, &
         tol_class, &
         tol_std, &
         dmin_half, &
         init_class

    double precision, dimension(n_iter) :: &
         dists, &
         peaks, prom, &
         thetas

    double precision, intent(out), dimension(n_iter) :: &
         peak_thetas

    double precision, dimension(n_iter,3,3) :: &
         peak_tmats

    double precision, dimension(n_iter,3) :: &
         vecs, &
         vecs_ordered, &
         peak_vecs

    integer, dimension(n_iter) :: &
         idx

    double precision, intent(out), dimension(n_iter,3,4) :: &
         ttrans, &
         rtrans

    double precision, dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, allocatable, dimension(:,:) :: &
         Apos_mapped, Apos_mapped_prev, &
         Bpos_opt, &
         tBpos_opt, &
         Bpos_opt_stretch, & ! position matrix
         rBpos_opt

    integer :: &
         i, j, k, &
         j_in, &
         id, &
         n_bins, &
         n_frac, &
         vecrep, &
         n_empty

    double precision, parameter, dimension(3) :: &
         zeros = (/0.0d0, 0.0d0, 0.0d0/)

    double precision, allocatable, dimension(:) :: &
         bins
    
    integer, allocatable, dimension(:) :: &
         classes_list, &
         bin_index

    integer, intent(out), dimension(n_iter, size(Apos,2)) :: &
         classes_list_out

    logical :: &
         remap, &
         free, &
         usebest, &
         exist, &
         findpeaks, &
         check

    character*200 :: &
         savebest, &
         progressfile

    character*20 :: &
         pot

    double precision :: &
         param, &
         zdist, &
         max_vol, &
         size_bin, &
         min_prom

    namelist /input2d/ &
         fracA, fracB, &
         tol, tol_std, &
         tol_class, &
         rate1, rate2, &
         n_ana, &
         n_conv, n_adjust, &
         n_class, &
         max_vol, &
         free, &
         savebest, &
         usebest, &
         remap, &
         pot, &
         param, &
         findpeaks, &
         check, &
         vecrep, &
         zdist, &
         min_prom, &
         init_class

    tol = 1.0d-6
    tol_std = tol*1.0d-3
    tol_class = 1.0d-3
    init_class = 1.0d0
    rate1 = 1.0d2
    rate2 = 1.0d2
    fracA = 0.0d0
    fracB = 0.25d0
    n_ana = 300
    n_conv = 5
    n_class = 30
    n_adjust = 10
    findpeaks = .false.
    free = .true.
    max_vol = 0.08d0
    savebest = trim(outdir)//"/best2d.dat"
    usebest = .false.
    remap = .true.
    pot = "LJ"
    param = 2.5d0
    check = .false.
    vecrep = 10
    min_prom = 0.6d0

    progressfile = trim(outdir)//"/progress.txt"
    open(13, file = trim(progressfile), status='replace')

    inquire(file = filename, exist=exist)

    if (exist) then
       open (unit = 11, file = filename, status = 'OLD')
       read (11, input2d)
       close (11)
    endif

    if (switched) then
       zdist = - param - maxval(Bpos(3,:))
    else
       zdist = param + maxval(Apos(3,:))
    endif
    
    n_frac = int(fracB*size(Apos,2)/sum(atoms))
    n_A = int(fracA*size(Apos,2)/sum(atoms))
    n_out = n_frac*sum(atoms)
    Apos_out = 0
    Bpos_out = 0

    allocate(Apos_mapped(3,n_out), Bpos_opt(3,n_out), &
         Apos_mapped_prev(3,n_out), &
         tBpos_opt(3,n_out), &
         rBpos_opt(3,n_out), &
         Bpos_opt_stretch(3,n_out), &
         classes_list(n_out))

    ! Center both cells at the geometric center
    call center(Bpos,nb)
    call center(Apos,na)

    if (check) then

       call init_random_seed()
       call random_number(angles)

       angles = angles*2*pi

       tmat = rot_mat(angles)

       vec = 0

    else

       if (usebest) then

          open(12, file = trim(savebest))
          read(12,*) tmat, vec, stats, tmats, vecs
          close(12)

       else
          if (free) then

             call gradient_descent_explore_free(.true., tmat, vec, stats, tmats,vecs, Apos, Bpos, Acell, iAcell, &
                  fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, vecrep, rate1, rate2, tol, &
                  max_vol, zdist, trim(pot), param)
             
          else

             call gradient_descent_explore(.true., angles, vec, Apos, Bpos, Acell, iAcell, &
                  fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, rate1, rate2, tol, &
                  zdist, trim(pot), param, trim(outdir))

             tmat = rot_mat(angles)
             
          endif

          if (trim(savebest)/="") then
             open(12, file = trim(savebest))
             write(12,*) tmat, vec, stats, tmats, vecs
             close(12)

          endif

       endif

    endif

    if (findpeaks) then

       write(13,*) "/======== Finding Peaks ========\\"

       flush(13)

       idx = sort(modulo(stats(:,1), 2.0d0*pi/dble(sym)))

       do i=1,n_iter
          thetas(i) = modulo(stats(idx(i),1), 2.0d0*pi/dble(sym))
          dists(i) = stats(idx(i),3)
          vecs_ordered(i,:) = vecs(idx(i),:)
          tmats_ordered(i,:,:) = tmats(idx(i),:,:)
       enddo

       vecs = vecs_ordered
       tmats = tmats_ordered

       n_bins = int(n_iter/10.0d0) + 1

       size_bin = 2.0d0*pi/dble(sym)/n_bins

       allocate(bins(n_bins), bin_index(n_bins+1))

       j=1
       n_empty=0
       bin_index(1) = 0
       do i=1, n_bins
          bins(i) = 0
          j_in = j
          do while (thetas(j) < i*size_bin .and. j < n_iter + 1)
             bins(i) = bins(i) + dists(j)
             j=j+1
          enddo
          bin_index(i+1) = j - 1
          if (j == j_in) then
             n_empty = n_empty + 1
             idx(n_empty) = i
          else
             bins(i) = bins(i) / (j - j_in)
          endif
       enddo

       do k=1, n_empty
          bins(idx(k)) = (bins(1+modulo(idx(k)-2, n_bins)) + bins(1 + modulo(idx(k), n_bins)))/2
       enddo
       
       call find_peaks(bins, n_bins, min_prom, n_peaks, peaks, idx, prom)
       
       do i=1, n_peaks

          id = bin_index(idx(i)) + minloc(dists(bin_index(idx(i))+1:bin_index(idx(i)+1)),1)
          
          peak_tmats(i,:,:) = tmats(id, :, :)
          peak_thetas(i) = thetas(id)
          peak_vecs(i,:) = vecs(id, :)

       enddo

       write(13,*) "Peak Angles:", 180 * modulo(peak_thetas(1:n_peaks), 2.0d0*pi/dble(sym)) / pi
       write(13,*) "Peak Distances:", peaks(1:n_peaks)
       write(13,*) "Peak Prominences:", prom(1:n_peaks)

       deallocate(bins, bin_index)

    else

       n_peaks = 1
       peak_vecs(1,:) = vec
       peak_tmats(1,:,:) = tmat
      
    endif

    classes_list_out = 0
    angles = 0.0d0

    do k=1,n_peaks

       if (findpeaks) then
          write(13,*) "New Peak #", k, "++++++++++++++++++++++++++++++++++++"
          write(13,*) "Angle:", modulo(peak_thetas(k)*180/pi, 360.0/sym)
       endif
       
       tmat = peak_tmats(k,:,:)
       vec = peak_vecs(k,:)

       write(13,*) "/======== Stretching Adjustment ========\\"

       flush(13)

       call adjust(.true., Apos, Bpos, Apos_mapped, Bpos_opt, &
            tmat, vec, &
            na, nb, check, &
            fracA, fracB, &
            atoms, n_atoms, rate1, rate2, &
            n_ana, n_conv, n_adjust, &
            tol, zdist, trim(pot), param)

       write(13,*) "/======== Classification ========\\"

       call classification(.true., .true., Apos, Bpos, Apos_mapped, Bpos_opt, & ! Output
            tmat, vec, & ! Output
            classes_list, &
            fracA, fracB, n_frac, &
            na, nb, n_class, &
            remap, check, &
            atoms, n_atoms, rate1, rate2, &
            n_ana, n_out, &
            tol, tol_class, tol_std, &
            init_class, &
            trim(pot), param)

       write(13,*) "Stretched distance:", distance(Apos_mapped, Bpos_opt, tmat, vec, pot, param)
       
       if (trim(pot)=="LJ") then
          
          vec(3) = zdist

          call analytical_gd_vec(.true.,tmat, vec, &
               Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
               tol, pot, param)

          if (remap) then
             i=0
             ! This is necessary because .true. in analytical_gd_vec that might affect the mapping
             do while (any(abs(Apos_mapped - Apos_mapped_prev) > 1.0d-12) .and. i < n_adjust)
                i = i + 1
                

                Apos_mapped_prev = Apos_mapped
                call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, free_trans(Bpos, tmat, vec), &
                     fracA, fracB, atoms, n_atoms, pot, param)

                call analytical_gd_vec(.true.,tmat, vec, &
                     Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
                     tol, pot, param)
             enddo
          endif


       else

          call analytical_gd_vec(.false.,tmat, vec, &
               Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
               tol, pot, param)
          
       endif

       write(13,*) "Stretched distance:", distance(Apos_mapped, Bpos_opt, tmat, vec, pot, param)
       
       call init_random_seed()

       call random_number(angles(1)) ! So that the rotation min doesn't get stuck 
       
       ! This step is just to get the "unstretched distance"
       if (rotmin) then
          vec_rot = 0.0d0
          call analytical_gd_rot(.true., angles, vec_rot, tmat, eye(), &
               100000, 1.0d0, 1.0d0, &
               tol, "Euclidean", 0.0d0)
       endif
        
       vec_rot = 0.0d0
       vec_rot(3) = zdist
       
       call analytical_gd_rot(.true., angles, vec_rot, Apos_mapped, Bpos_opt, &
            n_ana*1000, rate1, rate2, tol, pot, param)

       if (.not. findpeaks) then

          peak_thetas(k) = angles(1)

       endif

       write(13,*) "Final angle (degrees):", angles(1)*180/pi
       
       rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)

       rtrans(k,:,1:3) = rot_mat(angles)
       rtrans(k,:,4) = vec_rot

       write(13,*) "Final rmat"
       write(13,"(3(F10.6,X))") rot_mat(angles)

       classes_list_out(k,1:n_out) = classes_list

       write(13,*) "Final tmat"
       write(13,"(3(F10.6,X))") tmat
       
       Bpos_opt_stretch = free_trans(Bpos_opt, tmat, vec)

       ttrans(k,:,1:3) = tmat
       ttrans(k,:,4) = vec

       Bpos_out(k,:,1:n_out) = rBpos_opt
       Bpos_out_stretch(k,:,1:n_out) = Bpos_opt_stretch
       Apos_out(k,:,1:n_out) = Apos_mapped

       write(13,*) "Stretched distance:", distance(Apos_mapped, Bpos_opt, tmat, vec, pot, param)
       
       dmin(k,1) = distance(Apos_mapped, rBpos_opt, eye(), zeros, trim(pot), param)

       select case (trim(pot))
       case ("LJ")

       dmin(k,2) = distance(Apos_mapped, Bpos_opt_stretch, eye(), zeros, trim(pot), param)

       dmin(k,3) = 0.0d0
       
       case("Euclidean")

       call analytical_gd_rot(.true., angles, vec_rot, Apos_mapped(:, 1:int(n_out/2)), &
            rBpos_opt(:, 1:int(n_out/2)), &
            n_ana*1000, rate1, rate2, tol, pot, param)
          
       dmin_half = distance(Apos_mapped(:, 1:int(n_out/2)), rBpos_opt(:, 1:int(n_out/2)), &
            rot_mat(angles), vec_rot, trim(pot), param)
       
       dmin(k,2) = (dmin(k,1)/n_out - dmin_half/int(n_out/2))/(n_out**(1.0/4.0) - int(n_out/2)**(1.0/4.0))

       dmin(k,3) = dmin(k,1)/n_out - dmin(k,1)*n_out**(1.0/4.0)

       end select
       
       ! ! Print the cost matrix
       ! mat = cost(Apos,Bpos,n)   
       ! write(13,"(10(F5.3,X))") mat

    enddo

    close(13)

    deallocate(Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, &
         tBpos_opt, &
         Bpos_opt_stretch, &
         classes_list)

  end subroutine intoptimization

  subroutine fixed_tmat_int(Apos_out, Bpos_out, Bpos_out_stretch, &
       n_out, n_A, classes_list_out, vec, dmin, &
       Apos, na, Bpos, nb, &
       ttrans, &
       atoms, n_atoms, &
       switched, &
       filename, outdir)

    integer, intent(in) :: &
         na, nb, & ! Total number of atoms
         n_atoms ! Number of types of atoms per cell

    double precision :: &
         rate1, &  ! Rate of the gradient descent for angles
         rate2     ! Rate of the gradient descent for displacement

    character*200, intent(in) :: &
         filename, &
         outdir
    
    integer :: &
         n_ana, &
         n_conv, &
         n_adjust, &
         n_class

    double precision, intent(out), dimension(3) :: &
         dmin

    double precision, intent(inout), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(inout), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, intent(out), dimension(3,na) :: &
         Apos_out ! Position of the atoms

    double precision, intent(out), dimension(3,nb) :: &
         Bpos_out, & ! Position of the atoms
         Bpos_out_stretch

    integer, intent(out) :: &
         n_out, &
         n_A

    logical, intent(in) :: &
         switched
    
    double precision :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    integer, intent(in), dimension(n_atoms) :: &
         atoms ! Number of atoms of each type

    double precision, dimension(3) :: &
         angles  ! Rotation axis

    double precision :: &
         tol, &
         tol_class, &
         tol_std, &
         dmin_half, &
         init_class

    double precision, intent(out), dimension(3) :: &
         vec

    double precision, dimension(3) :: &
         vec_rot 
    
    double precision, intent(in), dimension(3,4) :: &
         ttrans

    double precision, dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, allocatable, dimension(:,:) :: &
         Apos_mapped, Apos_mapped_prev, &
         Bpos_opt, &
         tBpos_opt, &
         Bpos_opt_stretch, & ! position matrix
         rBpos_opt

    integer :: &
         n_frac, &
         vecrep, &
         i

    double precision, parameter, dimension(3) :: &
         zeros = (/0.0d0, 0.0d0, 0.0d0/)
    
    integer, allocatable, dimension(:) :: &
         classes_list

    integer, intent(out), dimension(size(Apos,2)) :: &
         classes_list_out

    logical :: &
         remap, &
         free, &
         usebest, &
         exist, &
         findpeaks, &
         check

    character*200 :: &
         savebest, &
         progressfile

    character*20 :: &
         pot

    double precision :: &
         param, &
         zdist, &
         max_vol, &
         min_prom

    namelist /input2d/ &
         fracA, fracB, &
         tol, tol_std, &
         tol_class, &
         rate1, rate2, &
         n_ana, &
         n_conv, n_adjust, &
         n_class, &
         max_vol, &
         free, &
         savebest, &
         usebest, &
         remap, &
         pot, &
         param, &
         findpeaks, &
         check, &
         vecrep, &
         zdist, &
         min_prom, &
         init_class

    tol = 1.0d-6
    tol_std = tol*1.0d-3
    tol_class = 1.0d-3
    init_class = 1.0d0
    rate1 = 1.0d2
    rate2 = 1.0d2
    fracA = 0.0d0
    fracB = 0.25d0
    n_ana = 300
    n_conv = 5
    n_class = 30
    n_adjust = 10
    findpeaks = .false.
    free = .true.
    max_vol = 0.08d0
    savebest = trim(outdir)//"/best2d.dat"
    usebest = .false.
    remap = .true.
    pot = "LJ"
    param = 2.5d0
    check = .false.
    vecrep = 10
    min_prom = 0.6d0

    progressfile = trim(outdir)//"/progress_find_cell.txt"
    open(13, file = trim(progressfile), status='replace')

    inquire(file = filename, exist=exist)

    if (exist) then
       open (unit = 11, file = filename, status = 'OLD')
       read (11, input2d)
       close (11)
    endif

    if (switched) then
       zdist = - param - maxval(Bpos(3,:))
    else
       zdist = param + maxval(Apos(3,:))
    endif
    
    n_frac = int(fracB*size(Apos,2)/sum(atoms))
    n_A = int(fracA*size(Apos,2)/sum(atoms))
    n_out = n_frac*sum(atoms)
    Apos_out = 0
    Bpos_out = 0

    allocate(Apos_mapped(3,n_out), Bpos_opt(3,n_out), &
         Apos_mapped_prev(3,n_out), &
         tBpos_opt(3,n_out), &
         rBpos_opt(3,n_out), &
         Bpos_opt_stretch(3,n_out), &
         classes_list(n_out))

    ! Center both cells at the geometric center
    call center(Bpos,nb)
    call center(Apos,na)
    
    tmat = ttrans(:,1:3)

    vec = ttrans(:,4)

    classes_list_out = 0
    angles = 0.0d0

    write(13,*) "/======== Classification ========\\"


    Apos_mapped = Apos(:,1:n_out)
    Bpos_opt = Bpos(:,1:n_out)
    
    call classification(.false., .true., Apos, Bpos, Apos_mapped, Bpos_opt, & ! Output
         tmat, vec, & ! Output
         classes_list, &
         fracA, fracB, n_frac, &
         na, nb, n_class, &
         remap, check, &
         atoms, n_atoms, rate1, rate2, &
         n_ana, n_out, &
         tol, tol_class, tol_std, &
         tol_class, &
         trim(pot), param)

       
       if (trim(pot)=="LJ") then
          
          vec(3) = zdist

          call analytical_gd_vec(.true.,tmat, vec, &
               Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
               tol, pot, param)

          if (remap) then
             i=0
             ! This is necessary because .true. in analytical_gd_vec that might affect the mapping
             do while (any(abs(Apos_mapped - Apos_mapped_prev) > 1.0d-12) .and. i < n_adjust)
                i = i + 1

                Apos_mapped_prev = Apos_mapped
                call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, free_trans(Bpos, tmat, vec), &
                     fracA, fracB, atoms, n_atoms, pot, param)

                call analytical_gd_vec(.true.,tmat, vec, &
                     Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
                     tol, pot, param)
             enddo
          endif
          
       else

          call analytical_gd_vec(.false.,tmat, vec, &
               Apos_mapped, Bpos_opt, n_ana*1000, rate2,&
               tol, pot, param)
          
       endif

       call init_random_seed()

       call random_number(angles(1)) ! So that the rotation min doesn't get stuck 
       
       ! This step is just to get the "unstretched distance"
       if (rotmin) then
          vec_rot = 0.0d0
          call analytical_gd_rot(.true., angles, vec_rot, tmat, eye(), &
               100000, 1.0d0, 1.0d0, &
               tol, "Euclidean", 0.0d0)
       endif

       write(13,*) "After rotmin"
       write(13,"(3(F10.6,X))") rot_mat(angles)
        
       vec_rot = 0.0d0
       vec_rot(3) = zdist
       
       call analytical_gd_rot(.true., angles, vec_rot, Apos_mapped, Bpos_opt, &
            n_ana*1000, rate1, rate2, tol, pot, param)

       write(13,*) "Final angle (degrees):", angles(1)*180/pi
       
       rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)

       write(13,*) "Final rmat"
       write(13,"(3(F10.6,X))") rot_mat(angles)

       classes_list_out(1:n_out) = classes_list

       write(13,*) "Final tmat"
       write(13,"(3(F10.6,X))") tmat
       
       Bpos_opt_stretch = free_trans(Bpos_opt, tmat, vec)

       Bpos_out(:,1:n_out) = rBpos_opt
       Bpos_out_stretch(:,1:n_out) = Bpos_opt_stretch
       Apos_out(:,1:n_out) = Apos_mapped

       dmin(1) = distance(Apos_mapped, rBpos_opt, eye(), zeros, trim(pot), param)

       select case (trim(pot))
       case ("LJ")

       dmin(2) = distance(Apos_mapped, Bpos_opt_stretch, eye(), zeros, trim(pot), param)

       dmin(3) = 0.0d0
       
       case("Euclidean")

       call analytical_gd_rot(.true., angles, vec_rot, Apos_mapped(:, 1:int(n_out/2)), &
            rBpos_opt(:, 1:int(n_out/2)), &
            n_ana*1000, rate1, rate2, tol, pot, param)
          
       dmin_half = distance(Apos_mapped(:, 1:int(n_out/2)), rBpos_opt(:, 1:int(n_out/2)), &
            rot_mat(angles), vec_rot, trim(pot), param)
       
       dmin(2) = (dmin(1)/n_out - dmin_half/int(n_out/2))/(n_out**(1.0/4.0) - int(n_out/2)**(1.0/4.0))

       dmin(3) = dmin(1)/n_out - dmin(1)*n_out**(1.0/4.0)

       end select
       
       ! ! Print the cost matrix
       ! mat = cost(Apos,Bpos,n)   
       ! write(13,"(10(F5.3,X))") mat

    close(13)

    deallocate(Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, &
         tBpos_opt, &
         Bpos_opt_stretch, &
         classes_list)

  end subroutine fixed_tmat_int
  
end module transform
