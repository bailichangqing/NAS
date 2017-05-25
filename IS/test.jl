import MPI
function main()
 MPI.Init()
 comm = MPI.COMM_WORLD
 rank = MPI.Comm_rank(comm)
 if rank == 0

 else
   println("this is proc 1")
 end
 MPI.Finalize()
end

main()
