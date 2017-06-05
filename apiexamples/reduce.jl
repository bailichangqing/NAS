import MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
else
	buffer = Int32[2,4,6]
end
tbf = MPI.Reduce(buffer,MPI.SUM,0,comm)
if rank == 1
	println(tbf)
end
MPI.Finalize()

