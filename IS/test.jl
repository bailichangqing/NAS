import MPI
function main()
 MPI.Init()
 comm = MPI.COMM_WORLD
 rank = MPI.Comm_rank(comm)
 if rank == 0
   # sender
   sendmsg = [1]
   sreq = MPI.Isend(sendmsg,1,0,comm)
   MPI.Wait!(sreq)
 else
   recvmsg = [0]
   rreq = MPI.Irecv!(recvmsg,0,0,comm)
   status = MPI.Wait!(rreq)
   println("I got msg ",recvmsg," from p 0")
   println(typeof(status))
 end
 MPI.Finalize()
end

main()
