program main_Q1
  
  use transform

  implicit none

  double precision, allocatable, dimension(:,:) :: &
       Apos, Bpos ! Position of the atoms

  double precision, allocatable, dimension(:,:) :: &
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

  allocate(Apos(3,100), Bpos(3,100), newBpos(3,100))

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
  call trans(Bpos,tetha,u,vec)
  call center(Apos)

  ! Print initial structures
  open(unit=11,file='Apos.txt', status='replace', action='write')
  open(unit=12,file='Bpos.txt', status='replace', action='write')

  write(11,'(200(F10.5,F10.5,F10.5,/))') Apos
  write(12,'(200(F10.5,F10.5,F10.5,/))') Bpos

  close(11)
  close(12)
  
  call center(Bpos)
  call center(Apos)

  d = dist(Apos,Bpos)

  dmin = d
  
  do i=1,n_iter

     call random_number(tetha)
     call random_number(vec)
     call random_number(u) ! Only in 3D
     u = (u - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))) ! Only in 3D
     u = u / norm(u) ! Only in 3D
     vec = vec - reshape((/0.5d0,0.5d0,0.5d0/),(/3,1/))
     vec(3,:) = 0 ! in 2D for testing

     vec = dmin * vec / sqrt(3.0d0) ! in 2D sqrt(3) in 3D
     tetha = 2*pi*tetha

     newBpos = Bpos
     
     call trans(newBpos,tetha,u,vec)

     d = dist(Apos,newBpos)

     if (d <= dmin) then
        dmin = d
        Bpos = newBpos
     endif
     
  enddo

  print*,'Final dmin', dmin

  
  ! Print initial structures
  open(unit=12,file='Bpos_final.txt', status='replace', action='write')

  write(12,'(200(F10.5,F10.5,F10.5,/))') Bpos

  close(11)
  close(12)
  
end program main_Q1
