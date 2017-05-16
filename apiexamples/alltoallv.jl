import MPI
MPI.Init()
comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
if rank == 0
	buffer = Int32[1,3,5]
	scount = Int32[0,3]
	rcount = Int32[0,3]
else
	buffer = Int32[2,4,6]
	scount = Int32[3,0]
	rcount = Int32[3,0]
end
tbf = Int32[0,0,0]
tbf = MPI.Alltoallv(buffer,scount,rcount,comm)
if rank == 0
	println(tbf)
end
MPI.Finalize()

