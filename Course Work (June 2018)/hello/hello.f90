program hello
include 'mpif.h'
integer myid, ierror, num_of_procs, tag, source, destination, count
integer status(MPI_STATUS_SIZE)
integer init


call MPI_INIT(ierror)
call MPI_COMM_SIZE(MPI_COMM_WORLD, num_of_procs, ierror)
call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierror)
count=1
destination=1
source=0
tag=1234
if (myid.eq.source) then
	init=5
	call MPI_SEND(init, count, MPI_INTEGER, 1, tag, MPI_COMM_WORLD, ierror )
	call MPI_SEND(init, count, MPI_INTEGER, 2, tag, MPI_COMM_WORLD, ierror )
	call MPI_SEND(init, count, MPI_INTEGER, 3, tag, MPI_COMM_WORLD, ierror )
endif
if (myid.eq.1) then
	call MPI_RECV(init, count, MPI_INTEGER, source, tag, MPI_COMM_WORLD, status, ierror)
	print*, 'node', myid, ':Hello, world! ', 'received number is ', init
endif

if (myid.eq.2) then
	call MPI_RECV(init, count, MPI_INTEGER, source, tag, MPI_COMM_WORLD, status, ierror)
	print*, 'node', myid, ':Hello, world! ', 'received number is ', init
endif

if (myid.eq.3) then
	call MPI_RECV(init, count, MPI_INTEGER, source, tag, MPI_COMM_WORLD, status, ierror)
	print*, 'node', myid, ':Hello, world! ', 'received number is ', init
endif
call MPI_FINALIZE(ierror)
stop
end
