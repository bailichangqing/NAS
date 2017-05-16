import MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
elseif rank == 1
	buffer = Int32[2,4,6]
end
MPI.Bcast!(buffer,length(buffer),1,comm)
if rank == 1
	sleep(2)
end
@printf("%02d: buffer = ",rank)
println(buffer)
MPI.Finalize()

