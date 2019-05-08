module transform

  use hungarian
  use omp_lib

  implicit none

  public :: &
       trans, &
       center, &
       fastoptimization, &
       canonicalize

  private :: &
       eye, &
       init_random_seed, &
       norm, &
       det, &
       free_trans, &
       cost_map, &
       classify, &
       analytical_gd_rot, &
       analytical_gd_free, &
       analytical_gd_std, &
       gradient_descent_explore

  double precision, parameter ::  bond = 1.0d0
  
contains

  subroutine init_random_seed()
    ! Copied form the GCC docs: https://gcc.gnu.org/onlinedocs/gcc-4.6.4/gfortran/RANDOM_005fSEED.html#RANDOM_005fSEED
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed

    CALL RANDOM_SEED(size = n)
    ALLOCATE(seed(n))

    CALL SYSTEM_CLOCK(COUNT=clock)

    seed = clock + 37 * (OMP_get_thread_num()+1) * (/ (i - 1, i = 1, n) /)
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

    tvec = vec(:,1) + sum(pos,2) / size(pos,2)

    call center(pos,n)

    pos = free_trans(pos, rot_mat(theta, u), tvec)

  end subroutine trans

  function rot_mat(theta,u) result(R)

    double precision, intent(in), dimension(3,1) :: &
         u     ! Rotation axis (unitary vector)

    double precision, intent(in) :: &
         theta    ! angle of rotation

    double precision, dimension(3,3) :: &
         Q, & ! cross product matrix
         P    ! u.u**T

    double precision, dimension(3,3) :: &
         R    ! Transformation matrix

    P = matmul(u,transpose(u))
    Q = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))

    R = P + (eye() - P)*cos(theta) + Q*sin(theta)

  end function rot_mat

  subroutine canonicalize(can_mat,tmat)

    double precision, dimension(3) :: &
         u_per     ! Vector perpendicular to the 2 first vectors

    double precision, intent(out), dimension(3,3) :: &
         can_mat ! Rotated (canonical) matrix

    double precision, intent(in), dimension(3,3) :: &
         tmat    ! Transformed matrix to rotate

    can_mat(:,1) = norm(tmat(:,1))*(/1,0,0/)
    
    u_per = tmat(:,2) -  tmat(:,1)*dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1))**2
    
    can_mat(:,2) = (/ dot_product(tmat(:,1), tmat(:,2))/norm(tmat(:,1)), &
         norm(u_per), 0.0d0/)
    
    can_mat(:,3) = (/ dot_product(tmat(:,1), tmat(:,3))/norm(tmat(:,1)), &
         dot_product(u_per, tmat(:,3))/norm(u_per), &
         det(tmat,3)/(norm(tmat(:,1))*norm(u_per))/)

    ! This is to avoid seeing an inversion in 2D as a rotation in a           
    ! prohibited axis. This is not necessary in 3D (but not wrong)            
    can_mat(2,:) = can_mat(3,3)/abs(can_mat(3,3))*can_mat(2,:)
    can_mat(3,3) = abs(can_mat(3,3))
    
  end subroutine canonicalize  

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

  function cost_map(Apos,Bpos, n_A, n_B) result(cost)

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
             !cost(i,j) = norm(Apos(:,i)-Bpos(:,j))
             cost(i,j) = distance(Apos(:,i:i),Bpos(:,j:j), eye(), (/0.0d0,0.0d0,0.0d0/))
          elseif (i <= n_A) then  
             cost(i,j) = 1000 ! TODO: Should that be a param?
          endif
       enddo
    enddo

  end function cost_map

  subroutine mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, fracA, fracB, atoms, n_atoms)
    
    integer, intent(in) :: &
         n_atoms ! Number of types of atoms
    double precision, intent(in), dimension(:,:) :: &
         Apos, Bpos, tBpos  ! position matrix
    integer, intent(in), dimension(n_atoms) :: &
         atoms
    double precision, intent(in) :: &
         fracA, fracB
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
            tBpos( : , id + 1 : id + n ), n_A, n_B)
       
       call munkres(dist_map, map, dmat, n)
       
       map = map + id

       do l=1, int(n*fracB)
          idx = idx + 1 
          Apos_mapped(:,idx) = Apos(:,map(l))
          Bpos_opt(:,idx) = Bpos(:,id+l)
       enddo

       deallocate(dmat,map)

    enddo

  end subroutine mapping

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

    function distance(Apos, Bpos, mat, vec)

    double precision, intent(in), dimension(:,:) :: &
         Apos, &  ! position matrix
         Bpos

    double precision, intent(in), dimension(3,3) :: &
         mat

    double precision, intent(in), dimension(3) :: &
         vec
    
    double precision :: &
         distance

    distance = sum(-(sum((Apos - free_trans(Bpos,mat,vec))**2,1)/bond**2+1.0d0)**(-3))
    
  end function distance

  function distance_b(Apos, Bpos, mat, vec)

    double precision, intent(in), dimension(:,:) :: &
         Apos, &  ! position matrix
         Bpos

    double precision, intent(in), dimension(3,3) :: &
         mat

    double precision, intent(in), dimension(3) :: &
         vec

    double precision :: &
         distance_b
    
    distance_b = sum(-(sum((Apos - free_trans(Bpos,mat,vec))**2,1)/bond**2+1.0d0)**(-3)) - &
         minval(-(sum((Apos - free_trans(Bpos,mat,vec))**2,1)/bond**2+1.0d0)**(-3))
    
  end function distance_b
    
  subroutine classify(n, n_classes, classes_list, disps, tol)

    double precision, intent(in), dimension(:,:) :: &
         disps

    double precision, dimension(size(disps,1),size(disps,2)) :: &
         classes

    integer, intent(out), dimension(size(disps,2)) :: &
         classes_list

    double precision, intent(in) :: &
         tol
         
    integer, intent(out), dimension(size(disps,2)) :: &
         n_classes
    
    integer, intent(out) :: &
         n

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
          ! if (abs(norm(classes(:,j)) - &
          !      dot_product(classes(:,j), disps(:,i)) / norm(classes(:,j))) < tol & 
          !      .and. norm(classes(:,j))**2*norm(disps(:,i))**2 - &
          !      dot_product(classes(:,j), disps(:,i))/norm(classes(:,j)) < tol) then

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

  subroutine analytical_gd_std(std, tmat, vec, Apos, Bpos, n_classes, classes, n_iter, rate1, rate2, tol)

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
         std_init, &
         nat

    integer :: &
         j, k, l, i ! Iterator

    nat = dble(size(Apos,2))

    std  = 1.0d0
    std_prev = 2.0d0 + tol
    std_init = std
           
    j=0
    do while (j < n_iter .and. abs(std - std_prev) > tol)
       j=j+1
       
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

          ! dist = sum(-(sum((Apos - free_trans(Bpos,tmat,vec))**2,1)+9.0d0)**(-3))

          E = Apos_class - free_trans(Bpos_class,tmat,vec)
          E = E - matmul(reshape(sum(E,2),(/3,1/))/size(E,2), ones)
          ! E = -E / spread((sum(E**2,1) + 9.0d0)**4,1,3)

          tmat_grad = tmat_grad + rate1 * std_prev / std_init * matmul(E,transpose(Bpos_class))
          vec_grad = vec_grad + rate2 * std_prev / std_init * matmul(E,transpose(ones))

          deallocate(Apos_class, Bpos_class, ones, E)

       enddo
       
       if (j==1) then
          std_init = std
          std_prev = std + 1.0d0 + tol
       endif

       ! print*, std, std_prev - std, "STD"

       if(std > std_prev) then
          std_init = 10 * std_init
          ! print*, "ADJUST", std, std_prev, std_init
          tmat = tmat_prev
          vec = vec_prev
          std = std_prev + 1.0d0 + tol
       elseif(std_init < tol .or. std < tol) then
          std_prev = std
       else
          tmat_prev = tmat
          vec_prev = vec
          
          ! print*, std, std_prev - std, "STD"
          tmat = tmat + tmat_grad
          vec = vec + vec_grad

       endif
       
    enddo

    print*, "END of STD", j, std, std_prev - std, std_init

  end subroutine analytical_gd_std

  subroutine analytical_gd_rot(theta, u, vec, Apos, Bpos, n_iter, rate1, rate2, scale_in, verbose)

    logical, intent(in) :: &
         verbose

    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         scale_in ! Scaling factor

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec, &     ! Translation vector
         u

    double precision, dimension(3,1) :: &         
         vec_prev     ! Translation vector
    
    double precision, intent(inout) :: &
         theta    ! Transformation matrix

    double precision, dimension(3,3) :: &
         Px, Py, Pt, & ! Temporary transformation matrix
         Qx, Qy, Qt, &
         Mx, My, Mt, M
    
    double precision, dimension(size(Bpos,2),1) :: &
         ones

    double precision :: &
         dist, &
         dist_prev, &
         scale, &
         theta_prev, &
         loc

    integer :: &
         i,j ! Iterator

    double precision, parameter :: &
         tol = 1d-8


    scale = scale_in

    ones = 1.0d0

    dist = 0
    dist_prev = tol+1

    theta_prev = theta
    vec_prev = vec
    
    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       M = rot_mat(theta,u)

       dist_prev = dist
       ! dist = sum(sqrt(sum((Apos - free_trans(Bpos,M,vec))**2,1)))
       dist = distance(Apos, Bpos, M, vec)
!       sum(-(sum((Apos - free_trans(Bpos,M,vec))**2,1)+9.0d0)**(-3))
       
       if (j==1) then
          dist_prev = dist + 1.0d0 + tol
       endif
       
       ! if (verbose) then
       ! print*, "ROT", j, dist, dist_prev - dist, scale
       ! endif

       if (verbose) then
          print*, "Theta in", theta, vec
          print*, Apos
          print*, Bpos
       endif
    
       E = Apos - free_trans(Bpos,M,vec)
       ! E = E / spread(sqrt(sum(E**2,1)),1,3)
       E = bond**6*E / spread((sum(E**2,1) + bond**2)**4,1,3)

       if (verbose) then
          print*, "E", E
       endif
       
       Pt = matmul(u,transpose(u))
       Qt = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))
       Mt = - (eye() - Pt)*sin(theta) + Qt*cos(theta)

       if(dist > dist_prev .or. dist > -1.0d-8) then
          ! print*, "ADJUST", dist, dist_prev
          scale = 1.0d0 / 10.0d0 * scale
          theta = theta_prev
          vec = vec_prev
       else

          theta_prev = theta
          vec_prev = vec
       
          theta = theta + rate1 * scale / dist * sum(matmul(E,transpose(Bpos)) * Mt)
          vec   = vec   + rate2 * scale / dist * matmul(E,ones)

       endif

       if (verbose) then
          print*,"-->", j , scale, dist, abs(dist - dist_prev), modulo(theta*180.0d0 &
               /3.14159265358d0,360.0d0), norm(vec)
       endif
       
    enddo

    if (verbose) then
       print*,"-->", j , scale, dist, abs(dist - dist_prev), theta, norm(vec)
    endif
       
  end subroutine analytical_gd_rot

  subroutine analytical_gd_rot_old(theta, u, vec, Apos, Bpos, n_iter, rate1, rate2, scale_in, verbose)

    logical, intent(in) :: &
         verbose

    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         scale_in ! Scaling factor

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(inout), dimension(3,1) :: &         
         vec, &     ! Translation vector
         u

    double precision, dimension(3,1) :: &         
         vec_prev     ! Translation vector
    
    double precision, intent(inout) :: &
         theta    ! Transformation matrix

    double precision, dimension(3,3) :: &
         Px, Py, Pt, & ! Temporary transformation matrix
         Qx, Qy, Qt, &
         Mx, My, Mt, M
    
    double precision, dimension(size(Bpos,2),1) :: &
         ones

    double precision :: &
         dist, &
         dist_prev, &
         scale, &
         theta_prev, &
         loc

    integer :: &
         i,j ! Iterator

    double precision, parameter :: &
         tol = 1d-6

    scale = scale_in

    ones = 1.0d0

    dist = 0
    dist_prev = tol+1

    theta_prev = theta
    vec_prev = vec
    
    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       M = rot_mat(theta,u)

       dist_prev = dist
       dist = sum(sqrt(sum((Apos - free_trans(Bpos,M,vec))**2,1)))
       ! dist = distance(Apos, Bpos, M, vec)
!       sum(-(sum((Apos - free_trans(Bpos,M,vec))**2,1)+9.0d0)**(-3))
       
       if (j==1) then
          dist_prev = dist + 1.0d0 + tol
       endif
       
       ! if (verbose) then
       ! print*, "ROT", j, dist, dist_prev - dist, scale
       ! endif

       E = Apos - free_trans(Bpos,M,vec)
       
       E = E / spread(sqrt(sum(E**2,1)),1,3)
       E(3,:) = 0.0d0
       E(:,3) = 0.0d0
       ! E = bond**6*E / spread((sum(E**2,1) + bond**2)**4,1,3)

       Pt = matmul(u,transpose(u))
       Qt = transpose(reshape((/0.0d0,-u(3,1),u(2,1),u(3,1),0.0d0,-u(1,1),-u(2,1),u(1,1),0.0d0/),(/3,3/)))
       
       Mt = - (eye() - Pt)*sin(theta) + Qt*cos(theta)
       
       if(dist > dist_prev) then
          ! print*, "ADJUST", dist, dist_prev
          scale = 10.0d0 * scale
          theta = theta_prev
          vec = vec_prev
       else

          theta_prev = theta
          vec_prev = vec

          theta = theta + rate1 / scale * dist * sum(matmul(E,transpose(Bpos)) * Mt)
          vec   = vec   + rate2 / scale * dist * matmul(E,ones)

       endif

       if (.false.) then
          print*,"-->", j , scale, dist, abs(dist - dist_prev), modulo(theta*180.0d0 &
               /3.14159265358d0,360.0d0), norm(vec)
       endif
       
    enddo

    if (verbose) then
       print*,"-->", j , scale, dist, abs(dist - dist_prev), theta, norm(vec)
    endif

  end subroutine analytical_gd_rot_old

  subroutine analytical_gd_free(tmat, vec, Apos, Bpos, sq, n_iter, rate1, rate2, scale_in, max_vol)

    integer, intent(in) :: &
         n_iter ! Number of atoms

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         scale_in, &
         max_vol
    
    logical, intent(in) :: &
         sq ! Square distance mode

    double precision, intent(in), dimension(:,:) :: &
         Apos, &
         Bpos ! Bpos ordered according to the mapping

    double precision, dimension(3,size(Bpos,2)) :: &
         E ! position matrix

    double precision, intent(inout), dimension(3,3) :: &         
         tmat     ! Translation vector

    double precision, dimension(3,3) :: &         
         ddet, &     ! Translation vector
         tmat_prev
    
    double precision, intent(inout), dimension(3,1) :: &         
         vec     ! Translation vector

    double precision, dimension(3,1) :: &         
         vec_prev     ! Translation vector
    
    double precision, dimension(size(Bpos,2),1) :: &
         ones

    double precision :: &
         dist, &
         D, &
         dist_prev, &
         dettmat, &
         scale

    integer :: &
         j ! Iterator

    double precision, parameter :: &
         tol = 1d-5, &
         importance = 1.0d6

    scale = scale_in
    
    ones = 1.0d0

    dist = 0
    dist_prev = tol+1

    tmat_prev = tmat
    vec_prev = vec

    ddet = 0

    j=0
    do while (j < n_iter .and. abs(dist - dist_prev) > tol)
       j=j+1

       dettmat = det(tmat,3)
       
       dist_prev = dist
       ! if (abs(dettmat - 1.0d0) >= max_vol) then
       !    dist = sum(sqrt(sum((Apos - free_trans(Bpos,tmat,vec))**2,1))) + &
       !    importance*(abs(dettmat-1.0d0) - max_vol)**2
       ! else
       !D = sum(sqrt(sum((Apos - free_trans(Bpos,tmat,vec))**2,1)))
       !dist = D * det(tmat,3)
       !dist = sum(sqrt(sum((Apos - free_trans(Bpos,tmat,vec))**2,1)))
       !dist = sum(-(sqrt(sum((Apos - free_trans(Bpos,tmat,vec))**2,1))+1.0d0)**(-6))
       ! dist = sum(-(sum((Apos - free_trans(Bpos,tmat,vec))**2,1)+9.0d0)**(-3))
       dist = distance(Apos, Bpos, tmat, vec)
       ! endif
          
       if (j==1) then
          dist_prev = dist + tol + 1.0d0
       endif
       
       E = Apos - free_trans(Bpos,tmat,vec)
       if (.not. sq) then
!          E = E / spread(sqrt(sum(E**2,1)),1,3)
          ! E = -E / spread(sqrt(sum(E**2,1))*(sqrt(sum(E**2,1)) + 1.0d0)**7,1,3)
          E = bond**6*E / spread((sum(E**2,1) + bond)**4,1,3)
       endif

       ddet(1,1) = tmat(2,2)
       ddet(1,2) = -tmat(2,1)
       ddet(2,1) = -tmat(1,2)
       ddet(2,2) = tmat(1,1)

       ! print*, dist, dist-dist_prev, dettmat, importance*(det(tmat,3) - max_vol)**2 
       
       if(dist > dist_prev .or. dist > -1.0d-8) then
          ! print*, "ADJUST", dist, dist_prev, scale
          scale = 1.0d0 / 10.0d0 * scale
          tmat = tmat_prev
          vec = vec_prev
       else

          tmat_prev = tmat
          vec_prev = vec
       
          ! if (abs(dettmat-1.0d0) >= max_vol) then
          !    tmat = tmat + rate1 * dist / scale * ( matmul(E,transpose(Bpos)) &
          !         - abs(dettmat-1.0d0) / (dettmat-1.0d0) * &
          !         importance*(abs(dettmat-1.0d0) - max_vol)*ddet)             
          ! else
          tmat = tmat + rate1 * scale / dist * matmul(E,transpose(Bpos))
          ! tmat = tmat + rate1 * dist / scale * (matmul(E,transpose(Bpos))*det(tmat,3) - D*ddet) 
          ! endif
          
          vec  = vec  + rate2 * scale / dist * matmul(E,ones)
          
       endif       
    enddo

    print*, j, dist, scale_in, dist-dist_prev 

  end subroutine analytical_gd_free
  
  subroutine gradient_descent_explore(theta,u,vec,stats, Apos, Bpos, cell, icell, &
       fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2)
    ! New Gradient Descent Random

    integer, intent(in) :: &
         n_iter, n_atoms, & ! Number of atoms
         n_ana, &
         n_conv

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         fracA, fracB ! Fraction of A and B to use in optimisation

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
         Bpos_opt
    
    double precision, intent(out), dimension(3) :: &
         vec, &     ! Translation vector
         u

    double precision, intent(out), dimension(n_iter,5) :: &
         stats
    
    double precision, dimension(3) :: &
         vec_local, &     ! Translation vector
         u_local
    
    double precision, intent(out) :: &
         theta

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
         diag, &
         theta_local

    double precision, allocatable, dimension(:) :: &
         dist_min, &
         theta_min

    double precision, allocatable, dimension(:,:) :: &
         u_min, &
         vec_min
    
    integer :: &
         i,j,k,l,vl, & ! Iterator
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

    mul_vec = diag*2/sqrt(2.0d0) !Only in 2D sqrt(3) in 3D
    
    !$omp parallel default(private) shared(dist_min, &
    !$omp theta_min, u_min, vec_min, u, theta, vec, stats) &
    !$omp firstprivate(n_iter, mul_vec, cell, fracA, fracB, &
    !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms)

    call init_random_seed()
    
    n_threads = OMP_get_num_threads()

    !$omp single
    allocate(dist_min(n_threads), &
         theta_min(n_threads), &
         u_min(3,n_threads),&
         vec_min(3,n_threads))
    !$omp end single

    thread = OMP_get_thread_num() + 1
    
    dist_min(thread) = sum(sqrt(sum((Apos - Bpos)**2,1)))
    
    !$omp do
    do j=1, n_iter
       
       call random_number(theta_local)
       call random_number(u_local)

       theta_local = theta_local*2*pi

       ! ! ! TMP testing
       ! theta_local = floor(theta_local*12)-6.0d0
       ! if (theta_local >= 0.0d0) then
       !    theta_local = 0.0526 + abs(theta_local)*pi/3.0d0
       ! else
       !    theta_local = -0.0526 + abs(theta_local)*pi/3.0d0
       ! endif
       ! ! ! TMP testing

       stats(j,5) = theta_local

       do vl=1, 10 ! TMP vec
       call random_number(vec_local)
          
       vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)
       vec_local(3) = 0.0d0 ! 2D only

       vec_local = vec_local * mul_vec
       
       vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

       ! ! 3D only
       ! u_local = u_local - (/0.5d0,0.5d0,0.5d0/)
       ! u_local = u_local / norm(u_local)
       ! u_local(3) = abs(u_local(3))

       u_local =  (/0.0d0,0.0d0,1.0d0/) ! 2D only

       write(*,*) "New initial step", thread, j
       
       do k=1, n_conv

          tBpos = free_trans(Bpos,rot_mat(theta_local,u_local),vec_local)

          call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms)

          if (k==1) then
             dist_cur = distance(Apos_mapped, Bpos_opt, rot_mat(theta_local,u_local), &
            vec_local) / size(Apos_mapped,2)
          endif

          call analytical_gd_rot(theta_local, u_local, vec_local, Apos_mapped, Bpos_opt, &
               n_ana, rate1, rate2, dist_cur, .false.)

       enddo

       dist_cur = distance(Apos_mapped, Bpos_opt, rot_mat(theta_local,u_local), &
            vec_local) / size(Apos_mapped,2) 

       write(*,*) "Opt dist", thread, dist_cur

       if (vl == 1 .or. dist_cur < stats(j,4)) then !TMP vec
       stats(j,1) = theta_local
       stats(j,2:3) = vec_local(1:2)
       stats(j,4) = dist_cur
       endif !TMP vec
       
       if (dist_cur < dist_min(thread)) then
          dist_min(thread) = dist_cur
          theta_min(thread) = theta_local
          u_min(:,thread) = u_local
          vec_min(:,thread) = vec_local
       endif

    enddo !TMP vec
    enddo
    !$omp end do

    !$omp barrier
    !$omp single
    pos = minloc(dist_min, 1)
    u = u_min(:,pos)
    theta = theta_min(pos)
    vec = vec_min(:,pos)

    print*, "Shortest Distance", minval(dist_min)
    print*, "u", u
    print*, "theta", theta*180/pi

    deallocate(dist_min, &
         theta_min, &
         u_min,&
         vec_min)
    !$omp end single 
    
    !$omp end parallel

  end subroutine gradient_descent_explore

    subroutine gradient_descent_explore_free(tmat,vec,stats,tmats, Apos, Bpos, cell, icell, &
       fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2, max_vol)
    ! New Gradient Descent Random

    integer, intent(in) :: &
         n_iter, n_atoms, & ! Number of atoms
         n_ana, &
         n_conv

    double precision, intent(in) :: &
         rate1, & ! Rate for angles
         rate2, & ! Rate for disp
         fracA, fracB, & ! Fraction of A and B to use in optimisation
         max_vol

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
         Bpos_opt
    
    double precision, intent(out), dimension(3) :: &
         vec     ! Translation vector

    double precision, intent(out), dimension(n_iter,5) :: &
         stats

    double precision, intent(out), dimension(n_iter,3,3) :: &
         tmats
    
    double precision, dimension(3) :: &
         u, &
         vec_local, &     ! Translation vector
         vecrot_local, &
         vecrot2_local, &
         u_local
    
    double precision, intent(out), dimension(3,3) :: &
         tmat

    double precision, dimension(3,3) :: &
         tmat_local, &    ! Transformation matrix
         tmat_local_reset

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
         shift, &
         diag, &
         theta_local, &
         theta2_local,&
         theta

    double precision, allocatable, dimension(:) :: &
         dist_min

    double precision, allocatable, dimension(:,:) :: &
         vec_min

    double precision, allocatable, dimension(:,:,:) :: &
         tmat_min
    
    integer :: &
         i,j,k,l,vl, & ! Iterator
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

    mul_vec = diag*2/sqrt(2.0d0) !Only in 2D sqrt(3) in 3D

    call init_random_seed()
    call random_number(shift)
    shift = 1.0d0 + 0.03d0*shift
    
    !$omp parallel default(private) shared(dist_min, &
    !$omp tmat_min, vec_min, tmat, u, theta, vec, stats, tmats) &
    !$omp firstprivate(n_iter, max_vol, shift, mul_vec, cell, fracA, fracB, &
    !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms)

    call init_random_seed()
    
    n_threads = OMP_get_num_threads()

    !$omp single
    allocate(dist_min(n_threads), &
         tmat_min(3,3,n_threads),&
         vec_min(3,n_threads))
    !$omp end single

    thread = OMP_get_thread_num() + 1
    
    dist_min(thread) = sum(sqrt(sum((Apos - Bpos)**2,1)))
    
    !$omp do
    do j=1, n_iter

       stats(j,5) = 0
       do while (stats(j,5) == 0)
          call random_number(theta_local)
          call random_number(tmat_local)

          theta_local = theta_local*2*pi

          tmat_local = 2*tmat_local - 1
          tmat_local(:,3) = 0
          tmat_local(3,:) = 0

          u_local =  (/0.0d0,0.0d0,1.0d0/) ! 2D only

          tmat_local = (1.0d0 + 0.05d0 * tmat_local(1,1)) ** (1.0d0/2.0d0) * eye()
          tmat_local(:,3) = 0
          tmat_local(3,:) = 0
          tmat_local(3,3) = 1.0d0
          
          tmat_local_reset = matmul(rot_mat(theta_local, u_local), tmat_local)
          ! tmat_local_reset = rot_mat(theta_local, u_local)
          ! tmat_local_reset = matmul(rot_mat(theta_local,u_local), max_vol/2.0d0 * tmat_local + eye())
          print*, "NEW VOL!!!!!", det(tmat_local,3)
          
          ! ! ! TMP testing
          ! theta_local = floor(theta_local*12)-6.0d0
          ! if (theta_local >= 0.0d0) then
          !    theta_local = 0.0526 + abs(theta_local)*pi/3.0d0
          ! else
          !    theta_local = -0.0526 + abs(theta_local)*pi/3.0d0
          ! endif
          ! ! ! TMP testing

          ! stats(j,5) = theta_local

          do vl=1, 10 ! TMP vec
             tmat_local = tmat_local_reset
             call random_number(vec_local)

             vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)
             vec_local(3) = 0.0d0 ! 2D only

             vec_local = vec_local*mul_vec

             vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

             ! ! 3D only
             ! u_local = u_local - (/0.5d0,0.5d0,0.5d0/)
             ! u_local = u_local / norm(u_local)
             ! u_local(3) = abs(u_local(3))

             write(*,*) "New initial step", thread, j, det(tmat_local,3)

             ! tmat_local = rot_mat(theta_local,u_local)

             do k=1, n_conv

                tBpos = free_trans(Bpos,tmat_local,vec_local)

                call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
                     fracA, fracB, atoms, n_atoms)

                if (k==1) then
                   dist_cur = distance(Apos_mapped, Bpos_opt, tmat_local, vec_local) / size(Apos_mapped,2) 
                endif

                call analytical_gd_free(tmat_local, vec_local, Apos_mapped, Bpos_opt, &
                     .false.,n_ana*100, rate1, rate2, dist_cur*10, max_vol)
                
             enddo

             ! call analytical_gd_rot(theta_local, u_local, vecrot_local, Apos_mapped, &
             !      matmul((1.0d0 + max_vol * tmat_local(1,1)) * eye(), Bpos_opt), &
             !      n_ana*100, rate1, rate2, dist_cur, .false.)

             vecrot_local = vec_local

             call analytical_gd_rot(theta_local, u_local, vecrot_local, Apos_mapped, &
                  Bpos_opt, &
                  n_ana*100, rate1, rate2, dist_cur, .false.)

             theta2_local = 0.0d0
             vecrot2_local = 0.0d0
             
             call analytical_gd_rot_old(theta2_local, u_local, vecrot2_local, tmat_local, eye(), &
               100000, 1.0d0, 1.0d0, &
               1.0d0, .true.)

             print*, "THETA VS THETA", modulo(theta_local*180/pi, 60.0d0), &
                  modulo(theta2_local*180/pi, 60.0d0)

             ! dist_cur = sum(sqrt(sum((Apos_mapped(:,1:int(0.5*size(Apos_mapped,2))) - &
             !      free_trans(Bpos_opt(:,1:int(0.5*size(Apos_mapped,2))), &
             !      tmat_local, vec_local))**2,1))) / (0.5*size(Apos_mapped,2))

             ! dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,tmat_local,vec_local))**2,1))) &
             !      / size(Apos_mapped,2)

             dist_cur = distance(Apos_mapped, Bpos_opt, tmat_local, vec_local) / size(Apos_mapped,2)
             
             ! dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,tmat_local,vec_local))**2,1))) &
             !      * det(tmat_local,3) / size(Apos_mapped,2)

             write(*,*) "Opt dist", thread, dist_cur, det(tmat_local,3)

             if (abs(det(tmat_local,3) - 1.0d0) < max_vol) then
                if (vl == 1 .or. dist_cur < stats(j,4)) then !TMP vec
                   stats(j,1) = theta2_local
                   stats(j,2:3) = vec_local(1:2)
                   stats(j,4) = dist_cur
                   stats(j,5) = det(tmat_local, 3)
                   tmats(j,:,:) = tmat_local
                   if (det(tmat_local, 3)==0) then
                      print*,"BIZZZZZZZARE"
                      print*, "TMAT", tmat
                   endif
                endif !TMP vec

                if (dist_cur < dist_min(thread)) then
                   dist_min(thread) = dist_cur
                   tmat_min(:,:,thread) = tmat_local
                   vec_min(:,thread) = vec_local
                endif
             endif

          enddo !TMP vec
       enddo
    enddo
    !$omp end do

    !$omp barrier
    !$omp single
    pos = minloc(dist_min, 1)
    tmat = tmat_min(:,:,pos)
    vec = vec_min(:,pos)

    print*, "Shortest Distance", minval(dist_min)
    print*, "tmat", tmat
    print*, "vec" ,vec

    deallocate(dist_min, &
         tmat_min, &
         vec_min)
    !$omp end single 
    
    !$omp end parallel

  end subroutine gradient_descent_explore_free
  
  ! subroutine gradient_descent_explore_free(tmat, vec, Apos, Bpos, cell, icell, &
  !      fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2)
  !   ! New Gradient Descent Random

  !   use omp_lib

  !   integer, intent(in) :: &
  !        n_iter, n_atoms, & ! Number of atoms
  !        n_ana, &
  !        n_conv

  !   double precision, dimension(3,3), intent(out) :: &
  !        tmat

  !   double precision, intent(in) :: &
  !        rate1, & ! Rate for angles
  !        rate2, & ! Rate for disp
  !        fracA, fracB ! Fraction of A and B to use in optimisation

  !   double precision :: &
  !        rand_rate1, & ! Rate for angles
  !        rand_rate2 ! Rate for disp

  !   integer, intent(in), dimension(n_atoms) :: &
  !        atoms

  !   double precision, intent(in), dimension(:,:) :: &
  !        Apos, Bpos ! Centered position of the atoms

  !   double precision, intent(in), dimension(3,3) :: &
  !        cell, &  ! cell
  !        icell ! inverse of cell
    
  !   double precision, dimension(3,size(Bpos,2)) :: &
  !        postmp, & ! position matrix
  !        tBpos

  !   double precision, dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
  !        Apos_mapped, & ! position matrix
  !        Bpos_opt
    
  !   double precision, intent(out), dimension(3) :: &
  !        vec      ! Translation vector

  !   double precision, dimension(3) :: &
  !        vec_local, &     ! Translation vector
  !        u_local

  !   double precision, dimension(3,3) :: &
  !        tmat_local, &
  !        P,Q

  !   double precision, allocatable, dimension(:,:) :: &
  !        dmat

  !   integer, allocatable, dimension(:) :: &
  !        map
    

  !   double precision :: &
  !        dist_plus, & ! distance when adding dx
  !        dist_minus, & ! distance when substracting dx
  !        accept, & ! Accept step
  !        dist_cur, &
  !        dist_map, &
  !        dist_stretch, &
  !        mul_vec, &
  !        diag, &
  !        theta_local

  !   double precision, allocatable, dimension(:) :: &
  !        dist_min

  !   double precision, allocatable, dimension(:,:) :: &
  !        vec_min, &
  !        tmat_min
    
  !   integer :: &
  !        i,j,k,l, & ! Iterator
  !        id, idx, &
  !        n, & ! Size of Bpos
  !        An_cell, Bn_cell, &
  !        n_threads, thread, &
  !        pos

  !   double precision, parameter :: &
  !        pi = 3.141592653589793d0

  !   diag = 0
  !   diag = max(norm(cell(:,1) + cell(:,2) + cell(:,3)),diag)
  !   diag = max(norm(-cell(:,1) + cell(:,2) + cell(:,3)),diag)
  !   diag = max(norm(cell(:,1) - cell(:,2) + cell(:,3)),diag)
  !   diag = max(norm(-cell(:,1) - cell(:,2) + cell(:,3)),diag)

  !   mul_vec = diag*2/sqrt(2.0d0) !Only in 2D sqrt(3) in 3D
    
  !   !$omp parallel default(private) shared(dist_min, &
  !   !$omp tmat_min, vec_min, tmat, vec) &
  !   !$omp firstprivate(n_iter, mul_vec, cell, fracA, fracB, &
  !   !$omp icell, n_conv, n_ana, Apos, Bpos, rate1, rate2, atoms, n_atoms)

  !   call init_random_seed()
    
  !   n_threads = OMP_get_num_threads()

  !   !$omp single
  !   allocate(dist_min(n_threads), &
  !        tmat_min(9,n_threads), &
  !        vec_min(3,n_threads))
  !   !$omp end single

  !   thread = OMP_get_thread_num() + 1
    
  !   dist_min(thread) = sum(sqrt(sum((Apos - Bpos)**2,1)))
    
  !   !$omp do
  !   do j=1, n_iter
       
  !      call random_number(theta_local)
  !      call random_number(u_local)
  !      call random_number(vec_local)

  !      theta_local = theta_local*2*pi

  !      vec_local = vec_local - (/0.5d0,0.5d0,0.5d0/)
  !      vec_local(3) = 0.0d0 ! 2D only

  !      vec_local = vec_local*mul_vec
       
  !      vec_local = vec_local - matmul(cell,nint(matmul(icell,vec_local))) 

  !      ! ! 3D only
  !      ! u_local = u_local - (/0.5d0,0.5d0,0.5d0/)
  !      ! u_local = u_local / norm(u_local)
  !      ! u_local(3) = abs(u_local(3))

  !      u_local =  (/0.0d0,0.0d0,1.0d0/) ! 2D only

  !      !tmat_local = rot_mat(theta_local,u_local)

  !      call random_number(tmat_local)

  !      tmat_local = 3.0d0*tmat_local-1.5d0
       
  !      write(*,*) "New initial step", thread, j
       
  !      do k=1, n_conv

  !         tBpos = free_trans(Bpos,tmat_local,vec_local)

  !         call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
  !              fracA, fracB, atoms, n_atoms)

  !         if (k==1) then
  !            dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,tmat_local,vec_local))**2,1)))
  !         endif

  !         call analytical_gd_free(tmat_local, vec_local, Apos_mapped, Bpos_opt, &
  !              .false.,n_ana*100, rate1, rate2, dist_cur*10)

  !      enddo

  !      dist_cur = sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt,tmat_local,vec_local))**2,1)))

  !      write(*,*) "Opt dist", thread, dist_cur

  !      if (dist_cur < dist_min(thread)) then
  !         dist_min(thread) = dist_cur
  !         tmat_min(:,thread) = reshape(tmat_local,(/9/))
  !         vec_min(:,thread) = vec_local
  !      endif

  !   enddo
  !   !$omp end do

  !   !$omp barrier
  !   !$omp single
  !   pos = minloc(dist_min, 1)
  !   tmat = reshape(tmat_min(:,pos),(/3,3/))
  !   vec = vec_min(:,pos)

  !   print*, "Shortest Distance", minval(dist_min)
    
  !   deallocate(dist_min, &
  !        tmat_min, &
  !        vec_min)
  !   !$omp end single 
    
  !   !$omp end parallel

  ! end subroutine gradient_descent_explore_free

  function sort(array) result(idx)
    
    double precision, intent(in), dimension(:) :: &
         array

    double precision, dimension(size(array)) :: &
         order

    integer, dimension(size(array)) :: &
         idx

    integer :: &
         k,i,n
    
    n = size(array)

    order(1) = array(1)
    idx(1) = 1

    do i=2, n 
       k=0
       do while (array(i) < order(i-1-k) .and. k < i-1)
          k=k+1
       enddo
       order(i-k+1:n) = order(i-k:n-1)
       order(i-k) = array(i)
       idx(i-k+1:n) = idx(i-k:n-1)
       idx(i-k) = i
    enddo

  end function sort

  function angle(mat)

    double precision, intent(in), dimension(3,3) :: mat
    double precision :: angle
    double precision, dimension(3,1) :: u, vec

    angle = 0.0d0
    u = reshape((/0.0d0, 0.0d0, 1.0d0/),(/3,1/))
    vec = 0.0d0
    
    call analytical_gd_rot_old(angle, u, vec, mat, eye(), &
               300, 1.0d-3, 1.0d-3, &
               1.0d0, .false.)
    
  end function angle

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
         i,j,n,id
    
    integer, parameter :: &
         peak_size = 10

    n=0
    do i = 1, n_array
       if (all(array(max(1, i-peak_size):max(1,i-1)) >= array(i)) .and. &
            all(array(min(i+1, n_array):min(i+peak_size, n_array)) >= array(i))) then
          n=n+1
          idx(n) = i
          peaks(n) = array(i)
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
    do i=1,n-1
       if (prom(i) > tol_prom) then
          j = j + 1
          peaks_out(j) = peaks(i)
          idx_out(j) = idx(i)
          prom_out(j) = prom(i)
       endif
    enddo

    n = j

  end subroutine find_peaks

  subroutine fastoptimization(Apos_out, Bpos_out, Bpos_out_stretch, &
       n_out, classes_list_out, ttrans, rtrans, dmin, stats, n_peaks, &
       peak_thetas, &
       Apos, na, Bpos, nb, &
       fracA, fracB, Acell, iAcell, atoms, n_atoms, &
       n_iter, n_ana, n_conv, n_adjust, sym, &
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
         n_conv, &
         n_adjust, &
         sym         ! Rotational symmetry

    double precision, intent(inout), dimension(3,na) :: &
         Apos ! Position of the atoms

    double precision, intent(inout), dimension(3,nb) :: &
         Bpos ! Position of the atoms

    double precision, intent(out), dimension(n_iter,3,na) :: &
         Apos_out ! Position of the atoms

    double precision, intent(out), dimension(n_iter,3,nb) :: &
         Bpos_out, & ! Position of the atoms
         Bpos_out_stretch

    double precision, dimension(3,nb) :: &
         tBpos ! Position of the atoms

    integer, intent(out) :: &
         n_out, &
         n_peaks

    double precision, intent(in), dimension(3,3) :: &
         Acell, & ! Unit cell of A
         iAcell

    double precision, intent(in) :: &
         fracA, fracB ! Fraction of A and B to use in optimisation

    double precision, allocatable, dimension(:,:) :: &
         inBpos ! Position of the atoms

    integer, intent(in), dimension(n_atoms) :: &
         atoms !Number of atoms of each type

    double precision, intent(out) :: &
         dmin

    double precision, intent(out), dimension(n_iter,5) :: &
         stats

    double precision, dimension(n_iter,3,3) :: &
         tmats, &
         tmats_ordered

    double precision, dimension(3) :: &
         vec, &     ! Translation vecto
         vec_rot, & ! Translation vector for the rotated unstretched matrix
         vec_rot2, &
         u, &       ! Rotation axis
         center_vec

    double precision :: &
         theta, & ! Rotation angle
         theta2, &
         d, &
         dist_map, &
         dist_stretch, &
         tol_adjust, &
         tol, &
         tol_class, &
         std, std_prev, &
         theta_prev 

    double precision, dimension(n_iter) :: &
         dists, meanrun,  &
         peaks, prom, &
         thetas

    double precision, intent(out), dimension(n_iter) :: &
         peak_thetas

    double precision, dimension(n_iter,3,3) :: &
         peak_tmats

    double precision, dimension(n_iter,2) :: &
         vecs, &
         peak_vecs

    integer, dimension(n_iter) :: &
         idx

    double precision, allocatable, dimension(:,:) :: &
         dmat

    double precision, dimension(na,nb) :: &
         mat

    double precision, intent(out), dimension(n_iter,3,4) :: &
         ttrans, &
         rtrans
    
    double precision, dimension(3,3) :: &
         tmat ! Transformation matrix

    double precision, &
         dimension(3,int(fracB*size(Apos,2)/sum(atoms))*sum(atoms)) :: &
         Apos_mapped, Apos_mapped_prev, &
         Bpos_opt, &
         tBpos_opt, &
         Bpos_opt_stretch, & ! position matrix
         disps

    integer :: &
         i, j, k, &
         n, id, &
         An_cell, Bn_cell, &
         n_B, n_tot, &
         n_frac, &
         size_mean

    double precision, parameter :: &
         pi = 3.141592653589793d0

    double precision, parameter, dimension(3) :: &
         zeros = (/0.0d0, 0.0d0, 0.0d0/)
    
    integer, allocatable, dimension(:) :: &
         classes_list, &
         classes_list_prev, &
         n_classes_trail, &
         n_classes

    integer, intent(out), dimension(n_iter, size(Apos,2)) :: &
         classes_list_out

    double precision :: &
         max_vol

    tol_class = 1.0d-3
    tol = 1.0d-10

    max_vol = 0.08d0
    
    n_frac = int(fracB*size(Apos,2)/sum(atoms))
    n_out = n_frac*sum(atoms)
    Apos_out = 0
    Bpos_out = 0

    allocate(classes_list(n_out), &
         classes_list_prev(n_out), &
         n_classes_trail(n_out))

    ! Center both cells at the geometric center
    call center(Bpos,nb)
    call center(Apos,na)

    ! call gradient_descent_explore(theta, u, vec, stats, Apos, Bpos, Acell, iAcell, &
    !      fracA, fracB, atoms,n_atoms,n_iter, n_ana, n_conv, rate1, rate2)

    ! call gradient_descent_explore_free(tmat,vec,stats,tmats, Apos, Bpos, Acell, iAcell, &
    !      fracA, fracB, atoms, n_atoms, n_iter, n_ana, n_conv, rate1, rate2, max_vol)


    tmat = reshape((/0.51600960391793982, -0.89233497231869008, 0.0000000000000000, &
         0.89019796996613731, 0.51507810840580193, 0.0000000000000000, &
         0.0000000000000000, 0.0000000000000000, 1.0000000000000000/),(/3,3/))

    vec = (/-0.30628812662975724, -0.47756853264300353, 0.0000000000000000/)
    
    ! print*, tmat, vec

    if (.false.) then
    
       idx = sort(modulo(stats(:,1), 2.0d0*pi/dble(sym)))
       
       do i=1,n_iter
          thetas(i) = stats(idx(i),1)
          dists(i) = stats(idx(i),4)
          vecs(i,:) = stats(idx(i),2:3)
          tmats_ordered(i,:,:) = tmats(idx(i),:,:)
       enddo

       tmats = tmats_ordered
       
!       size_mean = int(1d-2*n_iter)
       size_mean = 1
       
       do i=1, n_iter-size_mean+1
          meanrun(i) = sum(dists(i:i+size_mean-1))/dble(size_mean)
       enddo
       
       call find_peaks(meanrun, n_iter-size_mean+1, 0.6d0, n_peaks, peaks, idx, prom)
       
       do i=1,n_peaks

          peak_tmats(i,:,:) = tmats((idx(i)+minloc(dists(idx(i):idx(i)+size_mean-1),1)-1),:,:)
          peak_thetas(i) = thetas(idx(i)+minloc(dists(idx(i):idx(i)+size_mean-1),1)-1)
          peak_vecs(i,:) = vecs(idx(i)+minloc(dists(idx(i):idx(i)+size_mean-1),1)-1,:)
       
       enddo
       
       print*, "peak angles", 180 * modulo(peak_thetas(1:n_peaks), 2.0d0*pi/dble(sym)) / pi
       print*, "peaks", peaks(1:n_peaks)
       print*, "prom", prom(1:n_peaks)

    else

       n_peaks = 1
       ! peak_thetas(1) = theta
       peak_vecs(1,:) = vec(1:2)
       peak_tmats(1,:,:) = tmat

    endif
    
    u = (/0.0d0,0.0d0,1.0d0/)

    
    do k=1,n_peaks

       print*, "NEW PEAK +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
       print*, "Peak angle:", modulo(peak_thetas(k)*180/pi, 60.0d0)

       theta = 0.0d0
       ! tmat = rot_mat(peak_thetas(k),u)
       tmat = peak_tmats(k,:,:)

       vec_rot2 = 0.0d0
       theta2 = 0.0d0
       
       call analytical_gd_rot_old(theta2, u, vec_rot2, tmat, eye(), &
            100000, 1.0d0, 1.0d0, &
            1.0d0, .false.)

       print*, "Unperturbed tmat", modulo(theta2*180/pi, 60.0d0)
       
    vec = 0.0d0
    vec(1:2) = peak_vecs(k,:)

    do i=1,n_adjust

       theta_prev = theta

       do j=1,n_conv
       
          tBpos = free_trans(Bpos,tmat,vec)
          
          call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms)

          if (i==1 .and. j==1) then
             dmin = sum(sqrt(sum((Apos_mapped - Bpos_opt)**2,1)))
          endif
          
          theta = 0.0d0
          tBpos_opt = free_trans(Bpos_opt,tmat,zeros)

          
          call analytical_gd_rot(theta, u, vec, Apos_mapped, tBpos_opt, &
               n_ana*100, rate1, rate2, &
               distance(Apos_mapped, tBpos_opt,eye(),zeros), .false.)

          tmat = matmul(rot_mat(theta,u),tmat)
          
          theta = peak_thetas(k)
          vec_rot = vec

          
          call analytical_gd_rot(theta, u, vec_rot, Apos_mapped, Bpos_opt, &
               n_ana*10, rate1, rate2, &
               distance(Apos_mapped, Bpos_opt,eye(),zeros) , .false.)

          theta2 = 0.0d0
          vec_rot2 = 0

          call analytical_gd_rot_old(theta2, u, vec_rot2, tmat, eye(), &
               100000, 1.0d0, 1.0d0, &
               1.0d0, .false.)

          print*, "THETA vs THETA",  modulo(theta*180/pi, 60.0d0),  modulo(theta2*180/pi, 60.0d0)
          
       enddo

       write(*,*) "--------------> Adjustment step:", i

       write(*,*) "Stretched distance:", distance(Apos_mapped, Bpos_opt, tmat, vec) &
            / size(Apos_mapped,2)

       write(*,*) "Unstretched distance:", distance(Apos_mapped, Bpos_opt, rot_mat(theta,u), vec_rot)

       write(*,*) "Unstretched angles:", modulo(theta*180/pi, 60.0d0)
       
       call analytical_gd_free(tmat, vec, Apos_mapped, Bpos_opt,.false., &
            n_ana*100, rate1, rate2, distance(Apos_mapped, Bpos_opt,eye(),zeros), max_vol)
       
       
    enddo

    ! call analytical_gd_free(tmat, vec, Apos_mapped, Bpos_opt,.true., n_ana*100, rate1, rate2)
    
    ! ! Recenter
    ! vec = vec + sum(free_trans(Bpos_opt,rot_mat(theta,u),vec) - Apos_mapped,2) / n_out
    
    if (.true.) then
       write(*,*) "/======== Classification ========\\"

       tol_adjust = 1.0d0
       std = 1.0d0
       classes_list = 0
       classes_list_prev = 1
       j=0
       do while ( std > tol_class .and. j < 10)
          j = j + 1

          write(*,*) "-->", j

          ! REMAP

          center_vec = sum(free_trans(Bpos_opt,tmat,vec) - Apos_mapped,2) / n_out

          center_vec(1) = center_vec(1) - 0.1 
          
          Apos_mapped_prev = Apos_mapped
          
          tBpos = free_trans(Bpos,tmat,vec - center_vec)
          
          call mapping(Apos_mapped, Bpos_opt, Apos, Bpos, tBpos, &
               fracA, fracB, atoms, n_atoms)

          if (any(Apos_mapped_prev /= Apos_mapped)) then
             print*, "REMAP!"
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
               n_classes, classes_list, n_ana*1000, 1.0d-3, 1.0d-3, tol)

          deallocate(n_classes)

          write(*,*) "Number of classes", n_tot
          write(*,*) "Tolerance for classification:", tol_adjust
          write(*,*) "Final standard deviation:", sqrt(std/dble(size(Apos_mapped,2)-1))
          write(*,*) "Stretched distance", sum(sqrt(sum((Apos_mapped - free_trans(Bpos_opt, tmat, vec))**2,1)))
          
          if (std_prev - std < tol) then
             tol_adjust = tol_adjust/2.0d0
          endif
          
       enddo

       classes_list_out = 0
       classes_list_out(k,1:n_out) = classes_list

       print*,"Class list:", classes_list
       
       write(*,*) "Final tmat"
       write(*,"(3(F7.3,X))") tmat

       ! Reshift after classification
       ! center_vec = sum(free_trans(Bpos_opt,tmat,vec) - Apos_mapped,2) / n_out

       center_vec = 0.0d0
       
       ! End of calssification -----------------------------------
    else
       center_vec = 0.0d0
    endif

    print*, "Theta before", theta*180 / pi
    
    call analytical_gd_rot(theta, u, vec_rot, Apos_mapped, Bpos_opt, &
         n_ana*10, rate1, rate2, &
         distance(Apos_mapped, Bpos_opt,eye(),zeros) , .false.)

    print*, "Final angle:", theta*180 / pi
    
    Bpos_opt_stretch = free_trans(Bpos_opt,tmat,vec - center_vec)
    
    ttrans(k,:,1:3) = tmat
    ttrans(k,:,4) = vec + center_vec

    Bpos_opt = free_trans(Bpos_opt,rot_mat(theta,u),vec_rot)
    
    rtrans(k,:,1:3) = rot_mat(theta,u)
    rtrans(k,:,4) = vec_rot
    
    Bpos_out(k,:,1:n_out) = Bpos_opt
    Bpos_out_stretch(k,:,1:n_out) = Bpos_opt_stretch
    Apos_out(k,:,1:n_out) = Apos_mapped

    dmin = distance(Apos_mapped, Bpos_opt,eye(),zeros)

    enddo

    deallocate(classes_list, &
            classes_list_prev, &
            n_classes_trail)

    ! ! Print the cost matrix
    ! mat = cost(Apos,Bpos,n)   
    ! write(*,"(10(F5.3,X))") mat   
    ! call munkres(dmin,map,cost(Apos,Bpos,n),n)

  end subroutine fastoptimization

end module transform
