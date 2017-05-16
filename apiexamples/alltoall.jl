import MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
elseif rank == 1
	buffer = Int32[2,4,6]
elseif rank == 2
	buffer = Int32[7,9,11]
else
	buffer = Int32[8,10,12]
end
tbf = MPI.Alltoall(buffer,1,comm)
#tbf = Array(Int32,2)
#MPI.Alltoall!(buffer,1,tbf,1,comm)
if rank == 0
	println(tbf)	# Int32[1,2,7,8]
end
MPI.Finalize()

