module transform

  use hungarian

  implicit none

  public :: &
       center, &
       fastoptimization

  private :: &
       eye, &
       init_random_seed, &
       norm, &
       det, &
       free_trans, &
       cost_map, &
       analytical_gd_rot, &
       analytical_gd_free, &
       analytical_gd_slant, &
       analytical_gd_std, &
       analytical_gd_vol, &
       gradient_descent_explore, &
       classify, &
       pi

  double precision, parameter :: &
       pi = 3.141592653589793d0  

contains

  subroutine init_random_seed()
    ! Copied form the GCC docs: https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED

    use omp_lib

    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + OMP_get_thread_num() * 37 * (/ (i - 1, i = 1, n) /)

    CALL RANDOM_SEED(PUT = seed)

    DEALLOCATE(seed)

  end subroutine init_random_seed

  function rot_mat(angles) result(R)

    ! Creates the rotation matrix from the 3 angles
    
    double precision, intent(in), dimension(3) :: &
         angles     ! Rotation axis (unitary vector)

    double precision, dimension(3,1) :: &
         u     ! Rotation axis (unitary vector)

    double precision, dimension(3,3) :: &
         Q, & ! cross product matrix
         P    ! u.u**T

    double precision, dimension(3,3) :: &
         R    ! Transformation matrix

    u(1,1) = sin(angles(2)) * cos(angles(3))
    u(2,1) = sin(angles(2)) * sin(angles(3))
    u(3,1) = cos(angles(2))

    P = matmul(u,transpose(u))
    Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

    R = P + (eye() - P)*cos(angles(1)) + Q*sin(angles(1))

  end function rot_mat

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

  subroutine center(pos,n)

    ! Center the structures on their center of mass
    
    integer, intent(in) :: &
         n ! Number of atoms

    double precision, intent(inout), dimension(3,n) :: &
         pos  ! position matrix

    pos = pos - spread(sum(pos,2) / size(pos,2),2,size(pos,2))

  end subroutine center

  function free_trans(pos, mat, vec)

    ! Applies the transformation and transpose the structure by vec

    double precision, intent(in), dimension(:,:) :: &
         pos  ! position matrix

    double precision, intent(in), dimension(3,3) :: &
         mat ! Transfromation matrix

    double precision, intent(in), dimension(3) :: &
         vec ! Translation vector

    double precision, dimension(size(pos,1), size(pos,2)) :: &
         free_trans

    free_trans = matmul(mat, pos) + spread(vec, 2, size(pos,2))

  end function free_trans

  function cost_map(Apos,Bpos, n_A, n_B) result(cost)

    ! Creates the cost matrix (Euclidean distance)
    
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! Rotation axis
    integer, intent(in) :: &
         n_A, n_B
    double precision, dimension(size(Apos,2),size(Apos,2)) :: &
         cost
    integer :: &
         i,j,n ! Iterators

    cost = 0

    do i=1,size(Apos,2)
       do j=1,size(Bpos,2)
          if (j <= n_B) then
             cost(i,j) = norm(Apos(:,i)-Bpos(:,j))
          elseif (i <= n_A) then  
             cost(i,j) = 1000 ! Very high cost
          endif
       enddo
    enddo

  end function cost_map

  subroutine no_mapping(Apos_mapped, Bpos_opt, Apos, Bpos, fracA, fracB, atoms, n_atoms)

    ! Consider the order in which the atoms are given to be the mapping, creates
    ! the equivalent outputs to mapping
    
    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms ! Number of sections per atom (same as in p2ptrans.py)
    double precision, intent(in) :: &
         fracA, fracB ! Fraction of the spheres to be used
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

  subroutine mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, fracA, fracB, atoms, n_atoms)

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
         dist_map, d_tot
    integer :: &
         An_cell, & ! Number of cells
         i,l, &
         id, idx, &
         n, n_A, n_B

    An_cell = size(Apos,2)/sum(atoms)

    id = 0
    idx = 0
    n = 0
    d_tot = 0.0d0

    do i=0,n_atoms-1

       id = id + n

       n = An_cell*atoms(i+1)
       n_B = int(An_cell*fracB)*atoms(i+1)
       n_A = int(An_cell*fracA)*atoms(i+1)

       allocate(dmat(n,n), map(n))

       dmat = cost_map(Apos( : , id + 1 : id + n ), &
            tBpos( : , id + 1 : id + n ),n_A, n_B)

       call munkres(dist_map, map, dmat, n)

       d_tot = d_tot + dist_map

       map = map + id

       do l=1, n_B
          idx = idx + 1
          Apos_mapped(:,idx) = Apos(:,map(l))
          Bpos_opt(:,idx) = Bpos(:,id+l)
       enddo

       deallocate(dmat,map)

    enddo

  end subroutine mapping
  
  function eye() result(a)
    ! Creates the identity matrix
    ! Copied from Rosetta Code at: https://rosettacode.org/wiki/Identity_matrix#Fortran
    ! Checked and modified for double
    
    double precision :: a(3,3)
    integer :: i,j

    forall(i = 1:3, j = 1:3) a(i,j) = (i/j)*(j/i)

  end function eye

  function norm(a)

    ! Calculates the norm of a vector
    
    double precision :: norm
    double precision, dimension(3), intent(in) :: a

    norm = sqrt(sum(a**2))

  end function norm

  function split(a,char)

    ! Only keep what is before "char" 
    
    character*200 :: split
    character*1, intent(in) :: char
    character*200, intent(in) :: a
    integer :: i
    
    i=0
    do while (a(len(a) - i:len(a) - i) /= char .and. i < len(a))
       i = i + 1
    enddo
    
    split = trim(a(1:len(a)-i))
    
  end function split

  recursive function det(a,n) result(accumulation)
    ! Finds the determinant of an n by n matrix
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
         i, j ! Iterator

    nat = dble(size(Apos,2))

    ones = 1.0d0

    call init_random_seed()

    ! call random_number(angles)

    ! angles = (/0.0d0, pi/2.0d0, pi/2.0d0/) 

    angles = 0.0d0

    k = 1.0d0

    slant = eye()
    dist = 0.0d0
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
       dist = sum(sqrt(sum((Apos - free_trans(Bpos,slant,vec))**2,1)))

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

  subroutine analytical_gd_rot(angles, vec, Apos, Bpos, n_iter, rate1, rate2, tol)

    ! Gradient descent with respect to the 3 angles of rotation
    ! ADJSUT has been implemented
    ! TODO: Could be replaced with direct SVD minimization
    
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
         dist_2prev, &
         dist_init, &
         s2, s3, & ! Sine of angle 1 and 2
         c2, c3, & ! Cosine of angle 1 and 2
         nat

    integer :: &
         i, j      ! Iterator


    nat = dble(size(Apos,2))

    ones = 1.0d0

    M = rot_mat(angles)
    dist = 0.0d0
    dist_prev = tol + 1.0d0
    dist_2prev = dist_prev
    angles_prev = angles
    vec_prev = vec
    
    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       M = rot_mat(angles)

       dist_prev = dist
       dist = sum(sqrt(sum((Apos - free_trans(Bpos,M,vec))**2,1)))

       if (j==1) then
          dist_prev = dist + tol + 1.0d0
          dist_init = dist
       endif

       E = Apos - free_trans(Bpos,M,vec)
       E = E / spread(sqrt(sum(E**2,1)),1,3)

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
              
       ! print*, dist, dist_prev - dist, "ROT"

       if(dist > dist_prev .or. dist_init < 1.0d-8) then
          ! print*, "ADJUST", dist, dist_init
          dist_init = 10 * dist_init
          angles = angles_prev
          vec = vec_prev
       else

          angles_prev = angles
          vec_prev = vec
          
          angles(1) = angles(1) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M1)
          angles(2) = angles(2) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M2)
          angles(3) = angles(3) + rate1 * dist / dist_init * sum(matmul(E,transpose(Bpos)) * M3)
          vec       = vec       + rate2 * dist / dist_init * matmul(E,ones)

       endif
       
    enddo

  end subroutine analytical_gd_rot

  subroutine analytical_gd_free(tmat, vec, Apos, Bpos, sq, n_iter, rate1, rate2, tol)

    ! Gradient descent with respect to a linear transfromation (3x3 matrix)
    ! ADJUST has been implemented
    
    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         tol
    
    logical, intent(in) :: &
         sq ! Square distance mode

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(inout), dimension(3,3) :: &         
         tmat     ! Transformation matrix

    double precision, dimension(3,3) :: &         
         tmat_prev     ! Transformation matrix

    double precision, dimension(3,3) :: &         
         ddet ! Determinant matrix (cofactor matrix)
    
    double precision, dimension(2,2) :: &         
         ttmat ! Determinant matrix (similar to cofactor matrix)

    double precision, intent(inout), dimension(3,1) :: &         
         vec     ! Translation vector

    double precision, dimension(3,1) :: &         
         vec_prev     ! Translation vector
    
    double precision, dimension(size(Bpos,2),1) :: &
         ones

    double precision :: &
         dist, &
         dist_prev, &
         dist_init, &
         nat

    integer :: &
         j, k, l ! Iterator

    nat = dble(size(Apos,2))

    ones = 1.0d0

    dist = 0.0d0
    dist_prev = tol+1
    tmat_prev = tmat
    vec_prev = vec

    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       dist_prev = dist
       dist = sum(sqrt(sum((Apos - free_trans(Bpos,tmat,vec))**2,1)))
       
       if (j==1) then
          dist_prev = dist + tol + 1.0d0
          dist_init = dist
       endif

       ! print*, dist, dist_prev - dist, "FREE", j

       E = Apos - free_trans(Bpos,tmat,vec)

       if (.not. sq) then
          E = E / spread(sqrt(sum(E**2,1)),1,3)
       endif

       if(dist > dist_prev .or. dist_init < 1.0d-8) then
          dist_init = 10 * dist_init
          tmat = tmat_prev
          vec = vec_prev
       else

          tmat_prev = tmat
          vec_prev = vec
       
          tmat = tmat + rate1 * dist / dist_init * matmul(E,transpose(Bpos))
          vec = vec   + rate2 * dist / dist_init * matmul(E,ones)

       endif
       
    enddo

    ! print*, "free", j, dist  

  end subroutine analytical_gd_free
  
subroutine analytical_gd_vol(tmat, vec, Apos, Bpos, sq, n_iter, rate1, rate2, tol, ratio)

  ! Gradient descent with respect to a linear transformation where the determinant must
  ! be "ratio"
  ! ADJUST has not been implemented
  
  integer, intent(in) :: &
       n_iter ! Number of atoms

  double precision, intent(in) :: &
       rate1, & ! Rate for angles
       rate2, & ! Rate for disp
       tol, &
       ratio
  
  logical, intent(in) :: &
       sq ! Square distance mode

  double precision, intent(in), dimension(:,:) :: &
       Apos, &
       Bpos ! Bpos ordered according to the mapping

  double precision, dimension(3,size(Bpos,2)) :: &
       E ! position matrix

  double precision, intent(inout), dimension(3,3) :: &         
       tmat     ! Transformation matrix

  double precision, dimension(3,3) :: &         
       ddet ! Determinant matrix (cofactor matrix)
  
  double precision, dimension(2,2) :: &         
       ttmat ! Determinant matrix (similar to cofactor matrix)

  double precision, intent(inout), dimension(3,1) :: &         
       vec     ! Translation vector

  double precision, dimension(size(Bpos,2),1) :: &
       ones

  double precision :: &
       dist, &
       dist_prev, &
       dist_init, &
       nat

  integer :: &
       j, k, l ! Iterator

  nat = dble(size(Apos,2))

  ones = 1.0d0

  dist = 0.0d0
  dist_prev = tol+1

  j=0
  do while (j < n_iter .and. abs(dist - dist_prev) > tol)
     j=j+1

     tmat(1,1) = (ratio + tmat(1,2)*(tmat(2,1)*tmat(3,3)-tmat(3,1)*tmat(2,3)) + &
     - tmat(1,3)*(tmat(2,1)*tmat(3,2)-tmat(2,2)*tmat(3,1))) &
     / (tmat(2,2)*tmat(3,3) - tmat(3,2)*tmat(2,3))

     dist_prev = dist
     dist = sum(sqrt(sum((Apos - free_trans(Bpos,tmat,vec))**2,1)))
     
     if (j==1) then
        dist_init = dist
     endif

     ! print*, dist, dist_prev - dist, "VOL", j, det(tmat,3)

     E = Apos - free_trans(Bpos,tmat,vec)

     if (.not. sq) then
        E = E / spread(sqrt(sum(E**2,1)),1,3)
     endif

     ! Cofactor matrix
     do k=1,3
        do l=1,3
           ttmat(:k-1,:l-1) = tmat(:k-1,:l-1)
           ttmat(:k-1,l:) = tmat(:k-1,l+1:)
           ttmat(k:,:l-1) = tmat(k+1:,:l-1)
           ttmat(k:,l:) = tmat(k+1:,l+1:)
           
           ddet(k,l) = (1-2*modulo(k+l,2))*det(ttmat,2)
        enddo
     enddo

     tmat = tmat + rate1 * dist / dist_init * (matmul(E,transpose(Bpos)) - &
          dot_product(E(1,:),Bpos(1,:))*ddet/ddet(1,1))

     ! tmat = tmat + rate1 * dist / dist_init * (matmul(E,transpose(Bpos)))

     vec = vec   + rate2 * dist / dist_init * matmul(E,ones)

  enddo

  end subroutine analytical_gd_vol

  subroutine analytical_gd_std(std, tmat, vec, Apos, Bpos, n_classes, classes, n_iter, rate1, rate2, tol)

    ! Gradient descent with respect to a linear transformation to minimize the class specific
    ! standard deviation
    ! ADJUST has been implemented

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

    double precision, dimension(3,3) :: &         
         ddet ! Determinant matrix (cofactor matrix)
    
    double precision, dimension(2,2) :: &         
         ttmat ! Determinant matrix (similar to cofactor matrix)

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
          E = E - matmul(reshape(sum(E,2),(/3,1/))/size(E,2), ones)

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
          ! print*, "ADJUST", std, std_prev, std_init
          tmat = tmat_prev
          vec = vec_prev
          std = std_2prev
          std_prev = std + 1.0d0 + tol
       elseif(std_init < tol .or. std < tol) then
          std_prev = std
       else
          tmat_prev = tmat
          vec_prev = vec
          
          ! print*, std, std_prev, std_prev - std, "STD"
          tmat = tmat + tmat_grad
          vec = vec + vec_grad

       endif
       
    enddo

  end subroutine analytical_gd_std
  
  subroutine gradient_descent_explore(angles, vec, Apos, Bpos, cell, icell, &
       fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2, tol)
    ! New Gradient Descent Random

    ! Randomly prepare the structures and then minimize the distance with respect to a rotation
    
    use omp_lib

    integer, intent(in) :: &
         n_iter, n_atoms, & ! Number of atoms
         n_ana, &
         n_conv

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         fracA, fracB, & ! Fraction of A and B to use in optimisation
         tol

    double precision :: &
         rand_rate1, & ! Rate for angles
         rand_rate2 ! Rate for disp

    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Centered position of the atoms

    double precision, intent(in), dimension(3,3) :: &
         cell, &  ! cell
         icell ! inverse of cell
    
    double precision, dimension(3,size(Bpos,2)) :: &
         postmp, & ! position matrix
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

    double precision, dimension(3,3) :: &
         mat    ! Transformation matrix

    double precision, dimension(3,3) :: &
         mat_tmp, & ! Temporary transformation matrix
         mat_out, & ! Output transformation matrix
         mat_min, &
         P,Q

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
         mul_vec, &
         diag

    double precision, allocatable, dimension(:) :: &
         dist_min

    double precision, allocatable, dimension(:,:) :: &
         angles_min, &
         vec_min
    
    integer :: &
         i,j,k,l, & ! Iterator
         id, idx, &
         n, & ! Size of Bpos
         An_cell, Bn_cell, &
         n_threads, thread, &
         pos

    double precision, parameter :: &
         pi = 3.141592653589793d0
    
    diag = 0
    diag = max(norm(cell(:,1) + cell(:,2) + cell(:,3)),diag)
    diag = max(norm(-cell(:,1) + cell(:,2) + cell(:,3)),diag)
    diag = max(norm(cell(:,1) - cell(:,2) + cell(:,3)),diag)
    diag = max(norm(-cell(:,1) - cell(:,2) + cell(:,3)),diag)

    mul_vec = diag*2/sqrt(3.0d0)
    
    !$omp parallel default(private) shared(dist_min, &
    !$omp angles_min, vec_min, angles, vec) &
    !$omp firstprivate(n_iter, mul_vec, cell, fracA, fracB, tol, &
    !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms)

    call init_random_seed()
    
    n_threads = OMP_get_num_threads()

    !$omp single
    allocate(dist_min(n_threads), &
         angles_min(3,n_threads), &
         vec_min(3,n_threads))
    !$omp end single

    thread = OMP_get_thread_num() + 1
    
    dist_min(thread) = sum(sqrt(sum((Apos - Bpos)**2,1)))
    
    !$omp do
    do j=1, n_iter
       
       call random_number(angles_local)
       call random_number(vec_local)

       angles_local = angles_local*2*pi

       vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)

       vec_local = vec_local * mul_vec
       
       vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

       write(13,"(A, I4, A, I6)") "New initial step on thread", thread, ", iteration", j
       flush(13)
       
       Apos_mapped = 1.0d0
       Apos_mapped_prev = 0.0d0
       k=1

       do while ( k <= n_conv .and. any(Apos_mapped /= Apos_mapped_prev))
          k=k+1

          Apos_mapped_prev = Apos_mapped
          
          tBpos = free_trans(Bpos,rot_mat(angles_local),vec_local)

          call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms)

          call analytical_gd_rot(angles_local, vec_local, Apos_mapped, Bpos_opt, &
               n_ana, rate1, rate2, tol)

       enddo

       dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,rot_mat(angles_local),vec_local))**2,1)))

       write(13,"(A, I4, A, I6, A, F8.3)") &
            "Opt dist found for thread", thread,", iteration", j,":", dist_cur
       flush(13)
       
       if (dist_cur < dist_min(thread)) then
          dist_min(thread) = dist_cur
          angles_min(:,thread) = angles_local
          vec_min(:,thread) = vec_local
       endif

    enddo
    !$omp end do
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

  subroutine gradient_descent_explore_free(tmat, vec, Apos, Bpos, cell, icell, &
       fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2, tol, max_vol, &
       fixed_vol, ratio)
    ! New Gradient Descent Random

    ! Randomly prepare the structures and then minimize the distance with respect to a transformation matrix

    
    use omp_lib

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
         ratio

    logical, intent(in) :: &
         fixed_vol

    double precision :: &
         rand_rate1, & ! Rate for tmat
         rand_rate2 ! Rate for disp

    integer, intent(in), dimension(n_atoms) :: &
         atoms

    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos ! Centered position of the atoms

    double precision, intent(in), dimension(3,3) :: &
         cell, &  ! cell
         icell ! inverse of cell
    
    double precision, dimension(3,size(Bpos,2)) :: &
         postmp, & ! position matrix
         tBpos

    double precision, dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
         Apos_mapped, & ! position matrix
         Apos_mapped_prev, &
         Bpos_opt
    
    double precision, intent(out), dimension(3) :: &
         vec             ! Translation vector

    double precision, dimension(3) :: &
         vec_local, &       ! Translation vector
         vec_rot_local, &
         angles_local
         
    double precision, intent(out), dimension(3,3) :: &
         tmat    ! Transformation matrix

    double precision, dimension(3,3) :: &
         mat_tmp, & ! Temporary transformation matrix
         mat_out, & ! Output transformation matrix
         mat_min, &
         P,Q, &
         tmat_local

    double precision, allocatable, dimension(:,:) :: &
         dmat

    integer, allocatable, dimension(:) :: &
         map

    double precision :: &
         dist_plus, & ! distance when adding dx
         dist_minus, & ! distance when substracting dx
         accept, & ! Accept step
         dist_cur, &
         dist_cur_rot, &
         dist_map, &
         dist_stretch, &
         mul_vec, &
         diag, &
         dt_tmat

    double precision, allocatable, dimension(:) :: &
         dist_min

    double precision, allocatable, dimension(:,:,:) :: &
         tmat_min

    double precision, allocatable, dimension(:,:) :: &
         vec_min
    
    integer :: &
         i,j,k,l, & ! Iterator
         id, idx, &
         n, & ! Size of Bpos
         An_cell, Bn_cell, &
         n_threads, thread, &
         pos
    
    double precision, parameter :: &
         pi = 3.141592653589793d0

    diag = 0
    diag = max(norm(cell(:,1) + cell(:,2) + cell(:,3)),diag)
    diag = max(norm(-cell(:,1) + cell(:,2) + cell(:,3)),diag)
    diag = max(norm(cell(:,1) - cell(:,2) + cell(:,3)),diag)
    diag = max(norm(-cell(:,1) - cell(:,2) + cell(:,3)),diag)

    mul_vec = diag*2/sqrt(3.0d0)
    
    !$omp parallel default(private) shared(dist_min, &
    !$omp tmat_min, vec_min, vec, tmat) &
    !$omp firstprivate(n_iter, mul_vec, fixed_vol, ratio, max_vol, cell, fracA, fracB, tol, &
    !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms)

    call init_random_seed()
    
    n_threads = OMP_get_num_threads()

    !$omp single
    allocate(dist_min(n_threads), &
         tmat_min(3,3,n_threads), &
         vec_min(3,n_threads))
    !$omp end single

    thread = OMP_get_thread_num() + 1
    
    ! dist_min(thread) = sum(sqrt(sum((Apos - Bpos)**2,1)))
    dist_min(thread) = 2*sum(sqrt(sum((Apos - Bpos)**2,1))) ! TMP TESTING
    
    !$omp do
    do j=1, n_iter
       
       call random_number(tmat_local)
       call random_number(vec_local)
       
       tmat_local = max_vol * (tmat_local - 0.5d0)

       dt_tmat = det(tmat_local,3)

       if (fixed_vol) then
          ! tmat_local = dt_tmat / abs(dt_tmat) * tmat_local * abs((ratio / dt_tmat))**(1.0d0/3.0d0)
          tmat_local = dt_tmat / abs(dt_tmat) * tmat_local * (abs(modulo(abs(dt_tmat), ratio) / dt_tmat))**(1.0d0/3.0d0)
       else
          tmat_local = dt_tmat / abs(dt_tmat) * tmat_local * (abs(modulo(abs(dt_tmat), max_vol) / dt_tmat))**(1.0d0/3.0d0)
       endif

       vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)

       vec_local = vec_local * mul_vec
       
       vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

       write(13,"(A, I4, A, I6)") "New initial step on thread", thread, ", iteration", j
       flush(13)
       
       Apos_mapped = 1.0d0
       Apos_mapped_prev = 0.0d0
       k=1

       do while ( k <= n_conv .and. any(Apos_mapped /= Apos_mapped_prev))
          k=k+1

          Apos_mapped_prev = Apos_mapped
          
          tBpos = free_trans(Bpos,tmat_local,vec_local)

          call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms)

          call analytical_gd_free(tmat_local, vec_local, Apos_mapped, Bpos_opt, &
               .false., n_ana, rate1, rate2, tol)

       enddo

       angles_local = 0 !TMP TESTING
       vec_rot_local = vec_local !TMP TESTING
       
       call analytical_gd_rot(angles_local, vec_rot_local, Apos_mapped, Bpos_opt, &
            n_ana*1000, rate1, rate2, tol) ! TMP TESTING
       
       dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,tmat_local,vec_local))**2,1)))
       dist_cur_rot = sum(sqrt(sum((Apos_mapped - &
            free_trans(Bpos_opt,rot_mat(angles_local),vec_rot_local))**2,1))) !TMP TESTING

       write(13,"(A, I4, A, I6, A, F8.3)") &
            "Opt dist found for thread", thread,", iteration", j,":", dist_cur, dist_cur_rot 
       flush(13)
       
       if (dist_cur_rot + dist_cur < dist_min(thread)) then ! TMP TESTING
       ! if (dist_cur < dist_min(thread)) then ! TMP TESTING
          dist_min(thread) = dist_cur_rot + dist_cur
          tmat_min(:,:,thread) = tmat_local
          vec_min(:,thread) = vec_local
       endif

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

  subroutine fastoptimization(Apos_out, Bpos_out, Bpos_out_stretch, & ! Output
       n_out, n_A, classes_list_out, & ! Output
       tmat, dmin, vec, & ! Output
       Apos, na, Bpos, nb, &
       Acell, iAcell, ratio, &
       atoms, n_atoms, &
       filename, outdir)

    ! Main subroutine to exported in python. Minimizes the distance and resturns the optimal
    ! result
    
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

    double precision, dimension(3,nb) :: &
         tBpos ! Position of the atoms

    integer, intent(out) :: &
         n_out, &
         n_A

    double precision, intent(in), dimension(3,3) :: &
         Acell, & ! Unit cell of A
         iAcell

    double precision :: &
         fracA, fracB ! Fraction of A and B to use in optimisation
         
    double precision, intent(in) :: &
         ratio

    double precision, allocatable, dimension(:,:) :: &
         inBpos ! Position of the atoms

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, intent(out) :: &
         dmin

    double precision, intent(out), dimension(3) :: &
         vec        ! Translation vecto
    
    double precision, dimension(3) :: &
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         angles, &  ! Rotation angles
         center_vec ! Readjustment vector

    double precision :: &
         d, &
         dist_map, &
         dist_stretch, &
         std, std_prev, &
         tol_adjust, &
         tol, &
         tol_class, &
         tol_std, &
         max_vol

    double precision, allocatable, dimension(:,:) :: &
         dmat

    double precision, dimension(na,nb) :: &
         mat

    double precision, intent(out), dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, dimension(3,3) :: &
         slant, & ! Transformation matrix
         testmat

    double precision, allocatable, dimension(:,:) :: &
         Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, Apos_mapped_2prev, &
         disps, &
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
         n, &
         n_tot, &
         n_frac, &
         n_B

    integer, allocatable, dimension(:) :: &
         n_classes

    integer :: &
         i, j, &
         id, &
         An_cell, Bn_cell

    double precision, parameter :: &
         pi = 3.141592653589793d0

    logical :: &
         slanting, &
         fixed_vol, &
         remap, &
         free, &
         check, &
         usebest, &
         exist

    character*200 :: &
         savebest, &
         progressfile

    namelist /input/ &
         fracA, fracB, &
         tol, tol_std, &
         tol_class, &
         rate1, rate2, &
         fixed_vol, slanting, &
         n_iter, n_ana, &
         n_conv, n_adjust, &
         n_class, &
         max_vol, free, &
         check, &
         savebest, &
         usebest, &
         remap

    ! Namelist default values
    tol = 1.0d-4
    tol_std = tol*1.0d-3
    tol_class = 1.0d-3
    rate1 = 1.0d-5
    rate2 = 1.0d-5
    slanting = .false.
    fixed_vol = .false.
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
    Apos_out = 0
    Bpos_out = 0
    
    allocate(Apos_mapped(3,n_out), Bpos_opt(3,n_out), &
         Apos_mapped_prev(3,n_out), Apos_mapped_2prev(3,n_out), &
         disps(3,n_out), &
         tBpos_opt(3,n_out), &
         rBpos_opt(3,n_out), &
         Bpos_opt_stretch(3,n_out), &
         classes_list(n_out), &
         classes_list_prev(n_out), &
         n_classes_trail(n_out))
         
    
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

       
       if (free) then
          if (usebest) then
             
             open(12, file = trim(savebest))
             read(12,*) tmat, vec
             close(12)

          else
             
             call gradient_descent_explore_free(tmat, vec, Apos, Bpos, Acell, iAcell, &
                  fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, rate1, rate2, tol, &
                  max_vol, fixed_vol,  ratio)

             if (trim(savebest)/="") then
                open(12, file = trim(savebest))
                write(12,*) tmat, vec
                close(12)
             endif
          endif
       else
          if (usebest) then
             open(12, file = trim(savebest))
             read(12,*) angles, vec
             close(12)

          else

             call gradient_descent_explore(angles, vec, Apos, Bpos, Acell, iAcell, &
                  fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, rate1, rate2, tol)

             if (trim(savebest)/="") then
                open(12, file = trim(savebest))
                write(12,*) angles, vec
                close(12)
             endif
          endif
          
          tmat = rot_mat(angles)
       
       endif

    endif
    
    write(13,*) "/======== Stretching Adjustment ========\\"

    flush(13)
    
    vec_rot = vec
    Apos_mapped = 1.0d0
    Apos_mapped_prev = 0.0d0
    i=0

    if (n_adjust == 0) then
       tBpos = free_trans(Bpos,tmat,vec)
       call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
            fracA, fracB, atoms, n_atoms)
    endif
    
    do while (any(Apos_mapped /= Apos_mapped_2prev) .and. i < n_adjust)
    i = i + 1

       Apos_mapped_2prev = Apos_mapped
       Apos_mapped_prev = 0.0d0
       j=0
       do while (any(Apos_mapped /= Apos_mapped_prev) .and. j <= n_conv)
       j = j + 1
       
          Apos_mapped_prev = Apos_mapped
          
          tBpos = free_trans(Bpos,tmat,vec)

          if (check) then

             call no_mapping(Apos_mapped, Bpos_opt, Apos, Bpos, &
               fracA, fracB, atoms, n_atoms)

          else
             
             call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms)

          endif
             
          angles = 0
          tBpos_opt = free_trans(Bpos_opt,tmat,(/0.0d0,0.0d0,0.0d0/))
          
          call analytical_gd_rot(angles, vec, Apos_mapped, tBpos_opt, &
               n_ana*100, rate1, rate2, tol)
          
          tmat = matmul(rot_mat(angles),tmat)

       enddo

       vec_rot = vec
       
       ! This step is just to get the "unstretched distance"
       call analytical_gd_rot(angles, vec_rot, Apos_mapped, Bpos_opt, &
               n_ana*1000, rate1, rate2, tol)

       write(13,*) "-->", i

       write(13,*) "Stretched distance:", &
            sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt, tmat, vec))**2,1)))

       write(13,*) "Unstretched distance:", &
            sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt, rot_mat(angles), vec_rot))**2,1)))
       flush(13)
       
       call analytical_gd_free(tmat, vec, Apos_mapped, Bpos_opt,.false., n_ana*1000, rate1, rate2, tol)
       
    enddo

    if (slanting) then

        write(13,*) "/======== Deslanting ========\\"
        flush(13)
        
        rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)
        
        do i=1, n_adjust
        
           write(13,*) "-->", i
        
           tBpos_opt = free_trans(Bpos_opt, tmat, (/0.0d0,0.0d0,0.0d0/))
        
           angles = 0.0d0
        
           call analytical_gd_rot(angles, vec, rBpos_opt, tBpos_opt, &
                n_ana*1000, rate1, rate2, tol)
        
           tBpos_opt = free_trans(tBpos_opt, rot_mat(angles), (/0.0d0,0.0d0,0.0d0/))
        
           tmat = matmul(rot_mat(angles), tmat) 
        
           call analytical_gd_slant(slant, vec, rBpos_opt, tBpos_opt, n_ana*1000, rate1, rate2, tol)
           
           tmat = matmul(slant, tmat)

           ! print*, "Slant vol.", det(slant,3)
           ! print*, "Volume", det(tmat,3)
           
        enddo
        
        write(13,*) "Stretched distance after slanting:", &
             sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt, tmat, vec))**2,1)))
        flush(13)
        
    endif

    write(13,*) "/======== Classification ========\\"
    
    tol_adjust = 1.0d0
    std = 1.0d0
    classes_list = 0
    classes_list_prev = 1
    j=0
    do while ( std > tol_class .and. j < n_class)
       j = j + 1

       write(13,*) "-->", j
       
       if (remap .and. j/=1) then

          center_vec = sum(free_trans(Bpos_opt,tmat,vec) - Apos_mapped,2) / n_out

          Apos_mapped_prev = Apos_mapped
          
          tBpos = free_trans(Bpos,tmat,vec - center_vec)

          if (check) then
             
             call no_mapping(Apos_mapped, Bpos_opt, Apos, Bpos, &
                  fracA, fracB, atoms, n_atoms)

          else
             
             call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
                  fracA, fracB, atoms, n_atoms)

          endif
             
          if (any(Apos_mapped_prev /= Apos_mapped)) then
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

       call analytical_gd_std(std, tmat, vec, Apos_mapped, Bpos_opt, &
            n_classes, classes_list, n_ana*1000, rate1, rate2, tol_std)

       deallocate(n_classes)

       write(13,*) "Tolerance for classification:", tol_adjust
       write(13,*) "Final standard deviation:", sqrt(std/dble(size(Apos_mapped,2)-1))
       write(13,*) "Number of classes:", n_tot
       write(13,*) "Stretched distance", sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt, tmat, vec))**2,1)))
       flush(13)
       
       ! tol_adjust = min(3.0d0*sqrt(std/dble(size(Apos_mapped,2)-1)), tol_adjust/2.0d0) !3 sigma 99%

       
!       if (all(classes_list == classes_list_prev) .and. std_prev - std < tol) then
       if (std_prev - std < tol) then
          tol_adjust = tol_adjust/2.0d0
       endif
       

    enddo

    call analytical_gd_rot(angles, vec_rot, Apos_mapped, Bpos_opt, &
         n_ana*1000, rate1, rate2, tol)

    rBpos_opt = free_trans(Bpos_opt, rot_mat(angles), vec_rot)

    write(13,*) "Final rmat"
    write(13,"(3(F7.3,X))") rot_mat(angles)
    
    classes_list_out = 0
    classes_list_out(1:n_out) = classes_list
    
    write(13,*) "Final tmat"
    write(13,"(3(F7.3,X))") tmat
    
    ! Reshift after classification
    center_vec = sum(free_trans(Bpos_opt,tmat,vec) - Apos_mapped,2) / n_out

    Bpos_opt_stretch = free_trans(Bpos_opt,tmat,vec - center_vec)

    center_vec = sum(Bpos_opt_stretch - Apos_mapped,2) / n_out

    Bpos_out(:,1:n_out) = rBpos_opt
    Bpos_out_stretch(:,1:n_out) = Bpos_opt_stretch
    Apos_out(:,1:n_out) = Apos_mapped

    dmin = sum(sqrt(sum((Apos_mapped - rBpos_opt)**2,1)))

    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(13,"(10(F5.3,X))") mat

    close(13)
    
    deallocate(Apos_mapped, Bpos_opt, &
         Apos_mapped_prev, Apos_mapped_2prev, &
         disps, &
         tBpos_opt, &
         rBpos_opt, &
         Bpos_opt_stretch, &
         classes_list, &
         classes_list_prev, &
         n_classes_trail)

  end subroutine fastoptimization

end module transform