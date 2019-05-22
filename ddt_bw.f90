subroutine pingpong_usempi(buf, rank)
use mpi
implicit none
integer :: buf(0:1048577)
integer, intent(in) :: rank
double precision :: t0, t1
integer i, ierr

if (rank.eq.0) then
  do i=0,1048577
    buf(i) = i
  enddo
  call mpi_send(buf(1:1048576:2), 524288, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
  call mpi_recv(buf(1:1048576:2), 524288, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  t0 = MPI_Wtime()
  do i=1,1024
    call mpi_send(buf(1:1048576:2), 524288, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, ierr)
    call mpi_recv(buf(1:1048576:2), 524288, MPI_INTEGER, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  end do
  t1 = MPI_Wtime()
  write (*,*) 'BW fortran = ', 1/(t1-t0), ' GW/s'
elseif (rank.eq.1) then
  buf = -1
  call mpi_recv(buf(1:1048576:2), 524288, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  call mpi_send(buf(1:1048576:2), 524288, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
  do i=1,1024
    call mpi_recv(buf(1:1048576:2), 524288, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_send(buf(1:1048576:2), 524288, MPI_INTEGER, 0, 0, MPI_COMM_WORLD, ierr)
  end do
  do i=0,1048577,2
    if (buf(i).ne.-1) write (*,*) 'buf(', i, ') = ', i, ' != -1'
  enddo
  do i=1,1048576,2
    if (buf(i).ne.i) write (*,*) 'buf(', i, ') = ', i, ' != ', i
  enddo
endif
end subroutine

subroutine pingpong_ddt(buf, rank)
use mpi
implicit none
integer :: buf(0:1048577)
integer, intent(in) :: rank
double precision :: t0, t1
integer :: datatype
integer i, ierr

call mpi_type_vector(524288, 1, 2, MPI_INT, datatype, ierr)
call mpi_type_commit(datatype, ierr)

if (rank.eq.0) then
  do i=0,1048577
    buf(i) = i
  enddo
  call mpi_send(buf(1), 1, datatype, 1, 0, MPI_COMM_WORLD, ierr)
  call mpi_recv(buf(1), 1, datatype, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  t0 = MPI_Wtime()
  do i=1,1024
    call mpi_send(buf(1), 1, datatype, 1, 0, MPI_COMM_WORLD, ierr)
    call mpi_recv(buf(1), 1, datatype, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  end do
  t1 = MPI_Wtime()
  write (*,*) 'BW ddt = ', 1/(t1-t0), ' GW/s'
elseif (rank.eq.1) then
  buf = -1
  call mpi_recv(buf(1), 1, datatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
  call mpi_send(buf(1), 1, datatype, 0, 0, MPI_COMM_WORLD, ierr)
  do i=1,1024
    call mpi_recv(buf(1), 1, datatype, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
    call mpi_send(buf(1), 1, datatype, 0, 0, MPI_COMM_WORLD, ierr)
  end do
  do i=0,1048577,2
    if (buf(i).ne.-1) write (*,*) 'buf(', i, ') = ', i, ' != -1'
  enddo
  do i=1,1048576,2
    if (buf(i).ne.i) write (*,*) 'buf(', i, ') = ', i, ' != ', i
  enddo
endif
end subroutine

program pingpong
use mpi
implicit none

integer :: buf(0:1048577)
integer :: datatype, ierr
integer type_size
integer rank

call mpi_init(ierr)
call mpi_comm_rank(mpi_comm_world, rank, ierr)

call pingpong_usempi(buf, rank)
call mpi_barrier(MPI_COMM_WORLD, ierr)

call pingpong_ddt(buf, rank)
call mpi_barrier(MPI_COMM_WORLD, ierr)

call mpi_finalize(ierr)
end program
