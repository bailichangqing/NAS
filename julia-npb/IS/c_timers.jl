import MPI

start = Array(Float64,64)
elapsed = Array(Float64,64)

function timer_clear(n)
  elapsed[n + 1] = 0.0
end

function timer_start(n)
  start[n + 1] = ccall((:MPI_Wtime,"libmpi"),Float64,())
end

function timer_stop(n)
  now = ccall((:MPI_Wtime,"libmpi"),Float64,())
  t = now - start[n + 1]
  elapsed[n + 1] += t
end

function timer_read(n)
  return elapsed[n + 1]
end
