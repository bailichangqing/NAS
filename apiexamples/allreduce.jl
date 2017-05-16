import MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
else
	buffer = Int32[2,4,6]
end
tbf = Int32[0,0,0]
MPI.Allreduce!(buffer,tbf,MPI.SUM,comm)
if rank == 0
	println(tbf)
end
MPI.Finalize()

