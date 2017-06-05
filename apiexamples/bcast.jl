import MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
elseif rank == 1
	buffer = Int32[2,4,6]
end
ccall((:ctimer_start,"libcclock"),Void,())
MPI.Bcast!(buffer,length(buffer),1,comm)
cputime = ccall((:ctimer_stop,"libcclock"),Float64,())
if rank == 1
	sleep(2)
end
@printf("%02d: buffer = ",rank)
println(buffer)
MPI.Finalize()
if rank == 0
	println("cputime elapsed:	",cputime)
end

