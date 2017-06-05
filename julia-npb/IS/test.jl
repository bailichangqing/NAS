import MPI
function randlc(X::Float64,A::Float64)
  global randlc_KS,randlc_R23,randlc_R46,randlc_T23,randlc_T46
  randlc_KS = 0
  randlc_R23 = 0.0
  randlc_R46 = 0.0
  randlc_T23 = 0.0
  randlc_T46 = 0.0
  T1 = 0.0
  T2 = 0.0
  T3 = 0.0
  T4 = 0.0
  A1 = 0.0
  A2 = 0.0
  X1 = 0.0
  X2 = 0.0
  Z = 0.0
  i::Int32 = 0
  j::Int32 = 0

  if randlc_KS == 0
    randlc_R23 = 1.0
    randlc_R46 = 1.0
    randlc_T23 = 1.0
    randlc_T46 = 1.0

    for i = 1:23
      randlc_R23 = 0.50 * randlc_R23
      randlc_T23 = 2.0 * randlc_T23
    end
    for i = 1:46
      randlc_R46 = 0.50 * randlc_R46
      randlc_T46 = 2.0 * randlc_T46
    end
    randlc_KS = 1
  end
#/*  Break A into two parts such that A = 2^23 * A1 + A2 and set X = N.  */

  T1 = randlc_R23 * A
  j = trunc(Int32,T1)
  A1 = j
  A2 = A - randlc_T23 * A1

  #= Break X into two parts such that X = 2^23 * X1 + X2, compute
     Z = A1 * X2 + A2 * X1  (mod 2^23), and then
     X = 2^23 * Z + A2 * X2  (mod 2^46).                          =#

  T1 = randlc_R23 * X
  j = trunc(Int64,T1)
  X1 = j
  X2 = X - randlc_T23 * X1
  T1 = A1 * X2 + A2 * X1

  j = trunc(Int32,randlc_R23 * T1)
  T2 = j
  Z = T1 - randlc_T23 * T2
  T3 = randlc_T23 * Z + A2 * X2
  j = trunc(Int32,randlc_R46 * T3)
  T4 = j
  X = T3 - randlc_T46 * T4
  return randlc_R46 * X,X
end

function find_my_seed(kn::Int32,   # my processor rank, 0<=kn<=num procs
                      np::Int32,   # np = num procs
                      nn::Int64,   # total num of ran numbers, all procs
                      s::Float64,  # Ran num seed, for ex.: 314159265.00
                      a::Float64)  # Ran num gen mult, try 1220703125.00

  i::Int64 = 0

  t1::Float64 = 0.0
  t2::Float64 = 0.0
  t3::Float64 = 0.0
  an::Float64 = 0.0

  mq::Int64 = 0
  nq::Int64 = 0
  kk::Int64 = 0
  ik::Int64 = 0

  nq = div(nn,np)
  mq = 0
  while nq > 1
    mq += 1
    nq = div(nq,2)
  end
  t1 = a
  i = 1
  while i <= mq
    t2,t1 = randlc(t1,t1)
    i += 1
  end

  an = t1
  kk = kn
  t1 = s
  t2 = an

  for i = 1:100
    ik = div(kk,2)
    if 2 * ik != kk
      t3,t1 = randlc(t1,t2)
    end
    if ik == 0
      break
    end
    t3,t2 = randlc(t2,t2)
    kk = ik
  end
  return t1
end

function fmain()
 MPI.Init()
 comm = MPI.COMM_WORLD
 rank = MPI.Comm_rank(comm)
 if rank == 0
println(find_my_seed(Int32(1),
                     Int32(2),
                     262144,
                     314159265.00,      # Random number gen seed
                     1220703125.00 ))
 else

 end
 MPI.Finalize()
end

fmain()
