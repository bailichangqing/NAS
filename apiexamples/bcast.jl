import MPI
function tt(buffer)
	MPI.Bcast!(buffer,length(buffer),1,comm)
end
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
elseif rank == 1
	buffer = Int32[2,4,6]
end
#ccall((:ctimer_start,"libcclock"),Void,())
for i = 1:10
	tt(buffer)
end
if rank == 0
	@time tt(buffer)
else
	tt(buffer)
end
#cputime = ccall((:ctimer_stop,"libcclock"),Float64,())
MPI.Finalize()

