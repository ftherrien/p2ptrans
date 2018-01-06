program main
  
! Save every accepted state, and distance
! Random small displacement
! Hungarian algo at the end


! 0. Small displacement, calc dist, initial transformation, run code
! 1. Tesselation: In python, library
! 2. Different types
! 3. Hunagrian 
 
  use transform

  implicit none

  integer, parameter :: &
       n = 100 !Number of atoms

  double precision, dimension(3,n) :: &
       Apos, Bpos ! Position of the atoms

  double precision, dimension(3,n) :: &
       newBpos ! new position of the atoms

  double precision, dimension(3,1) :: &
       vec, & ! Translation vector
       u      ! Rotation axis
       
  double precision :: &
       d, &   ! Distance between the two structures
       dmin, &   ! Minimal dstance between structures
       tetha   ! Rotation angle
       
  integer :: &
       i   ! Iterator

  integer, parameter :: &
       n_iter = 100000 ! Number of Monte-Carlo iterations

  double precision, parameter :: &
       pi = 3.141592653589793d0

  !u = reshape((/0.0d0,0.0d0,1.0d0/),(/3,1/)) ! in 2D 
  
  call init_random_seed()
  call random_number(Apos)
  call random_number(vec)
  call random_number(tetha)
  call random_number(u) ! Only in 3D
  u = (u - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))) ! Only in 3D
  u = u / norm(u) ! Only in 3D

  !Apos(3,:) = 0 ! in 2D for testing
  !vec(3,:) = 0 ! in 2D for testing
  
  Bpos = Apos

  print*, 'Angle (deg)',tetha*260
  print*, 'Disp', vec

  tetha = 2*pi*tetha
  call trans(Bpos,n,tetha,u,vec)

  ! Print initial structures
  open(unit=11,file='Apos.txt', status='replace', action='write')
  open(unit=12,file='Bpos.txt', status='replace', action='write')

  write(11,'(200(F10.5,F10.5,F10.5,/))') Apos
  write(12,'(200(F10.5,F10.5,F10.5,/))') Bpos

  close(11)
  close(12)
  
  call map(Apos, Bpos, n)
  
  ! Print initial structures
  open(unit=12,file='Apos_final.txt', status='replace', action='write')
  open(unit=13,file='Bpos_final.txt', status='replace', action='write')

  write(12,'(200(F10.5,F10.5,F10.5,/))') Apos  
  write(13,'(200(F10.5,F10.5,F10.5,/))') Bpos

  close(11)
  close(12)
  close(13)
  
end program main
