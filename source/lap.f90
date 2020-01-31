! Linear assignment problem (LAP) solver
! -----------------------------------------------
! The first method (method = 1) is a brute force
! solution that computes all the permutations (n!)
! using Heap's algorithm - time complexity O(n!)
! The second method (method = 2) uses the Munkres
! algorithm (also called the Hungarian algo) to
! solve the problem - time complexity O(n^3)
! -----------------------------------------------
! F. Brieuc - March 2017

! Modified by Felix Therrien to keep only Khun-Munkres and to work
! with double precision
! January 2018

module hungarian

  implicit none

  public :: &
       munkres
  private :: &
       step1, step2, step3, step4, step5, step6

contains

  subroutine munkres(sumSol,jSol,C,n)
    ! Implementation of the Munkres algorithm (also referred to as the Hungarian
    ! algorithm). J. Munkres, Journal of the SIAM 5, 1 (1957)
    ! The following implementation is based on
    ! http://csclab.murraystate.edu/%7Ebob.pilgrim/445/munkres.html

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(in),  dimension(n,n) :: C   ! cost matrix (min sum)
    double precision, dimension(n,n) :: CC   ! cost matrix (min sum)
    
    double precision, intent(out) :: sumSol ! maximal sum
    integer, dimension(n), intent(out) :: jSol ! solution indices
 
    integer :: step, i, j, tmp
    logical :: done

    integer, dimension(:,:), allocatable :: M    ! mask matrix
    integer, dimension(:), allocatable :: rowCover, colCover !cover row and cols

    integer :: pathRow0, pathCol0  ! starting point for path finding part

    integer :: test

    pathRow0 = 0
    pathCol0 = 0
    
    done = .false.
    step = 1
    tmp = 0

    CC = C

    allocate(M(n,n))      ! mask matrix - contains starred zeros
    allocate(rowCover(n)) ! to keep track of covered rows
    allocate(colCover(n)) ! to keep track of covered columns

    do i = 1, n
       M(:,i) = 0
       rowCover(i) = 0
       colCover(i) = 0
    enddo

    do while(.not. done)
       select case(step)
       case(1)
          call step1(step,n,CC)
       case(2)           
          call step2(step,n,CC,M,rowCover, colCover)
       case(3)           
          call step3(step,n,CC,M, colCover)
       case(4)           
          call step4(step,n,CC,M,rowCover, colCover, pathRow0, pathCol0)
       case(5)           
          call step5(step,n,CC,M,rowCover, colCover, pathRow0, pathCol0)
       case(6)           
          call step6(step,n,CC,rowCover, colCover)
       case default ! done
          do i = 1, n
             do j = 1, n
                if (M(j,i) == 1) jSol(i) = j
             enddo
          enddo
          done = .true.
       end select
    enddo

    sumSol = 0
    do i = 1, n
       sumSol = sumSol + C(jSol(i),i)
    enddo

    
    deallocate(M)
    deallocate(rowCover)
    deallocate(colCover)

  end subroutine munkres

  subroutine step1(step, n, CC)
    ! row reduction : for each row find the smallest value and substract it from
    ! all elements of that row. Go to step 2.

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(inout),  dimension(n,n) :: CC   ! cost matrix (min sum)

    integer, intent(out) :: step

    integer :: i, j
    
    double precision :: minVal

    do i = 1, n
       minVal = CC(1,i)
       do j = 1, n
          if (CC(j,i) < minVal) minVal = CC(j,i)
       enddo
       CC(:,i) = CC(:,i) - minVal
    enddo

    step = 2

  end subroutine step1

  subroutine step2(step, n, CC, M, rowCover, colCover)
    ! Search for zeros.
    ! Find a zero (Z) in the matrix. If no zeros has been previously starred in
    ! its row and column then star Z. Go to step 3.

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(inout),  dimension(n,n) :: CC   ! cost matrix (min sum)

    integer, dimension(:), intent(inout) :: rowCover, colCover !cover row and cols
    
    integer, dimension(:,:), intent(inout) :: M    ! mask matrix
    
    integer, intent(out) :: step

    integer :: i, j

    do i = 1, n
       do j = 1, n
          if (CC(j,i) == 0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
             M(j,i) = 1
             rowCover(i) = 1
             colCover(j) = 1
          endif
       enddo
    enddo
    ! uncovers
    do i = 1, n
       rowCover(i) = 0
       colCover(i) = 0
    enddo

    step = 3

  end subroutine step2

  subroutine step3(step, n, CC, M, colCover)
    ! cover each column containing a starred zero. If n column are covered
    ! the starred zero describe an optimal assignment and we are done otherwise
    ! go to step 4.

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(inout),  dimension(n,n) :: CC   ! cost matrix (min sum)

    integer, dimension(:), intent(inout) :: colCover !cover row and cols
    
    integer, dimension(:,:), intent(inout) :: M    ! mask matrix

    integer, intent(out) :: step

    integer :: colCount, i, j

    colCount = 0
    do i = 1, n
       do j = 1, n
          ! if starred and column is uncovered
          if (M(j,i) == 1 .and. colCover(j) == 0) then
             colCover(j) = 1
             colCount = colCount + 1
          endif
       enddo
    enddo

    if (colCount == n) then
       step = 0
    else
       step = 4
    endif

  end subroutine step3

  subroutine step4(step, n, CC, M, rowCover, colCover,pathRow0, pathCol0)
    ! Find a uncovered zero and prime it. If there is no starred zero in the row
    ! go to step 5. Otherwise, cover the row and uncover the column containing
    ! the starred zero. Continue until no uncovered zeros is left. Go to step 6.

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(inout),  dimension(n,n) :: CC   ! cost matrix (min sum)

    integer, dimension(:), intent(inout) :: rowCover, colCover !cover row and cols
    
    integer, dimension(:,:), intent(inout) :: M    ! mask matrix

    integer, intent(inout) :: pathRow0, pathCol0  ! starting point for path finding part
    
    integer, intent(out) :: step

    logical :: done, starInRow
    integer :: i, j, row, col

    done = .false.

    do while (.not. done)
       ! find an uncovered zero
       row = 0; col = 0
       starInRow = .false.
       loop1: do i = 1, n
          loop2: do j = 1, n
             if (CC(j,i) == 0 .and. rowCover(i) == 0 .and. colCover(j) == 0) then
                row = i
                col = j
                exit loop1
             endif
          enddo loop2
       enddo loop1
       if (row == 0) then !no zero uncoverred left
          done = .true.
          step = 6
       else
          M(col,row) = 2 !primed zero
          ! search if there is a starred zero in the same row
          do j = 1, n
             if (M(j,row) == 1) then
                starInRow = .true.
                col = j
             endif
          enddo
          if (starInRow) then ! if there is a starred zero in line
             rowCover(row) = 1
             colCover(col) = 0
          else ! if no starred zero in line
             done = .true.
             step = 5
             pathRow0 = row
             pathCol0 = col
          endif
       endif
    enddo

  end subroutine step4

  subroutine step5(step, n, CC, M, rowCover, colCover, pathRow0, pathCol0)
    ! Augmenting path algorithm: construct a serie of alternating primed and
    ! starred zeros as follows. Let Z0 be the uncoverd primed zero found in
    ! step 4. Let Z1 be the starred zero in the column of Z0 (if any).
    ! Let Z2 be the primed zero in the row of Z1 (there will always be one).
    ! Continue until the series terminates at a primed zero that has no starred
    ! zero in its column. Then unstar each starred zeros of the series, star
    ! each primed zeros of the series, erase all primes and uncover every line
    ! and columns. Return to step 3.

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(inout),  dimension(n,n) :: CC   ! cost matrix (min sum)

    integer, intent(out) :: step

    integer, dimension(:,:), intent(inout) :: M    ! mask matrix

    integer, dimension(:), intent(inout) :: rowCover, colCover !cover row and cols

    integer, intent(inout) :: pathRow0, pathCol0  ! starting point for path finding part
    
    logical :: done
    integer :: i, j
    integer :: row, col
    integer :: pathCount
    integer, dimension(2*n+1,2) :: path

    path = 0
    
    pathCount = 1

    path(pathCount,1) = pathRow0
    path(pathCount,2) = pathCol0
    
    done = .false.

    do while (.not. done)

       if (path(1,1) /= pathRow0 .or. path(1,2) /= pathCol0) call abort
       
       ! search for a starred zero in column
       row = 0
       col = path(pathCount,2)
       
       do i = 1, n
          if (M(col,i) == 1) row = i
       enddo
       
       if (row /= 0) then ! update path
          pathCount = pathCount + 1
          path(pathCount,1) = row
          path(pathCount,2) = path(pathCount-1,2)
       else
          done = .true.
       endif
       
       if (.not. done) then
          ! search for a prime zero in row
          do j = 1, n
             if (M(j,row) == 2) col = j
          enddo
          ! update path
          pathCount = pathCount + 1
          path(pathCount,1) = path(pathCount-1,1)
          path(pathCount,2) = col
       endif
       
    enddo

    ! augment path
    do i = 1, pathCount
       if(M(path(i,2),path(i,1)) == 1) then
          M(path(i,2),path(i,1)) = 0
       else
          M(path(i,2),path(i,1)) = 1
       endif
    enddo

    ! clear covers and erase primes
    do i = 1, n
       rowCover(i) = 0
       colCover(i) = 0
       do j = 1, n
          if (M(j,i) == 2) M(j,i) = 0
       enddo
    enddo

    step = 3

  end subroutine step5

  subroutine step6(step, n, CC, rowCover, colCover)
    ! Search for the smallest uncovered value and add it to covered rows
    ! and substract it from uncovered columns. Return to step 4.

    integer, intent(in) :: n      ! dimension of C - assumed to be a (nxn) square matrix
    double precision, intent(inout),  dimension(n,n) :: CC   ! cost matrix (min sum)

    integer, dimension(:), intent(inout) :: rowCover, colCover !cover row and cols
    
    integer, intent(out) :: step

    integer :: i, j

    double precision :: minVal

    minVal = huge(i)

    do i = 1, n
       do j = 1, n
          if (rowCover(i) == 0 .and. colCover(j) == 0 .and. CC(j,i) < minVal) then
             minVal = CC(j,i)
          endif
       enddo
    enddo
    do i = 1, n
       do j = 1, n
          if (rowCover(i) == 1) CC(j,i) = CC(j,i) + minVal
          if (colCover(j) == 0) CC(j,i) = CC(j,i) - minVal
       enddo
    enddo

    step = 4

  end subroutine step6

end module hungarian
