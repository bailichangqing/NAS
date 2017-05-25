import MPI
include("c_timers.jl")
include("c_print_results.jl")


#/******************/
#/* default values */
#/******************/
if isdefined(:CLASS) == false
  CLASS = 'S'
  NUM_PROCS = 2
end
MIN_PROCS = 1

#/*************/
#/*  CLASS S  */
#/*************/
if CLASS == 'S'
  TOTAL_KEYS_LOG_2 = 16
  MAX_KEY_LOG_2 = 11
  NUM_BUCKETS_LOG_2 = 9
end

#/*************/
#/*  CLASS W  */
#/*************/
if CLASS == 'W'
  TOTAL_KEYS_LOG_2 = 20
  MAX_KEY_LOG_2 = 16
  NUM_BUCKETS_LOG_2 = 10
end

#/*************/
#/*  CLASS A  */
#/*************/
if CLASS == 'A'
  TOTAL_KEYS_LOG_2 = 23
  MAX_KEY_LOG_2 = 19
  NUM_BUCKETS_LOG_2 = 10
end

#/*************/
#/*  CLASS B  */
#/*************/
if CLASS == 'B'
  TOTAL_KEYS_LOG_2 = 25
  MAX_KEY_LOG_2 = 21
  NUM_BUCKETS_LOG_2 = 10
end

#/*************/
#/*  CLASS C  */
#/*************/
if CLASS == 'C'
  TOTAL_KEYS_LOG_2 = 27
  MAX_KEY_LOG_2 = 23
  NUM_BUCKETS_LOG_2 = 10
end

#/*************/
#/*  CLASS D  */
#/*************/
if CLASS == 'D'
  TOTAL_KEYS_LOG_2 = 29
  MAX_KEY_LOG_2 = 27
  NUM_BUCKETS = 10
  MIN_PROCS = 4
end

TOTAL_KEYS = 1 << TOTAL_KEYS_LOG_2
MAX_KEY = 1 << MAX_KEY_LOG_2
NUM_BUCKETS = 1 << NUM_BUCKETS_LOG_2
NUM_KEYS = div(TOTAL_KEYS , NUM_PROCS) * MIN_PROCS

#/*****************************************************************/
#/* On larger number of processors, since the keys are (roughly)  */
#/* gaussian distributed, the first and last processor sort keys  */
#/* in a large interval, requiring array sizes to be larger. Note */
#/* that for large NUM_PROCS, NUM_KEYS is, however, a small number*/
#/* The required array size also depends on the bucket size used. */
#/* The following values are validated for the 1024-bucket setup. */
#/*****************************************************************/

if NUM_PROCS < 256
  SIZE_OF_BUFFERS = 3 * div(NUM_KEYS , 2)
elseif NUM_PROCS < 512
  SIZE_OF_BUFFERS = 5 * div(NUM_KEYS , 2)
elseif NUM_PROCS < 1024
  SIZE_OF_BUFFERS = 4 * div(NUM_KEYS , 2)
else
  SIZE_OF_BUFFERS = 13 * div(NUM_KEYS , 2)
end

#/*****************************************************************/
#/* NOTE: THIS CODE CANNOT BE RUN ON ARBITRARILY LARGE NUMBERS OF */
#/* PROCESSORS. THE LARGEST VERIFIED NUMBER IS 1024. INCREASE     */
#/* MAX_PROCS AT YOUR PERIL                                       */
#/*****************************************************************/
if CLASS == 'S'
  MAX_PROCS = 128
else
  MAX_PROCS = 1024
end

MAX_ITERATIONS = 10
TEST_ARRAY_SIZE = 5

#/***********************************/
#/* Enable separate communication,  */
#/* computation timing and printout */
#/***********************************/
TIMING_ENABLED = 0
if isdefined(:NO_MTIMERS) == true
  TIMING_ENABLED = -1
end
timeron = 0
function TIMER_START(x)
  global timeron,NO_MTIMERS
  if TIMING_ENABLED == -1
  else
    if timeron == 1
      timer_start(x)
    end
  end
end

function TIMER_STOP(x)
  global timeron,NO_MTIMERS
  if(TIMING_ENABLED == -1)
  else
    if timeron == 1
      timer_stop(x)
    end
  end
end
if TIMING_ENABLED == 0
  T_TOTAL = 0
  T_RANK = 1
  T_RCOMM = 2
  T_VERIFY = 3
  T_LAST = 3
end


#/*************************************/
#/* Typedef: if necessary, change the */
#/* size of int here by changing the  */
#/* int type to, say, long            */
#/*************************************/



#/********************/
#/* MPI properties:  */
#/********************/
my_rank = 0
comm_size = 0

#/********************/
#/* Some global info */
#/********************/
total_local_keys = 0
total_lesser_keys = 0
key_buff_ptr_global = 0


key_buff_ptr_global = 0       # used by full_verify to get
total_lesser_keys = 0         # copies of rank info
passed_verification = 0

#/************************************/
#/* These are the three main arrays. */
#/* See SIZE_OF_BUFFERS def above    */
#/************************************/
key_array = Array{Int32}(SIZE_OF_BUFFERS)
fill!(key_array,0)
key_buff1 = Array{Int32}(SIZE_OF_BUFFERS)
fill!(key_buff1,0)
key_buff2 = Array{Int32}(SIZE_OF_BUFFERS)
fill!(key_buff2,0)
bucket_size = Array{Int32}(NUM_BUCKETS+TEST_ARRAY_SIZE)
fill!(bucket_size,0)
bucket_size_totals = Array{Int32}(NUM_BUCKETS+TEST_ARRAY_SIZE)
fill!(bucket_size_totals,0)
bucket_ptrs = Array{Int32}(NUM_BUCKETS)
fill!(bucket_ptrs,0)
process_bucket_distrib_ptr1 = Array{Int32}(NUM_BUCKETS+TEST_ARRAY_SIZE)
fill!(process_bucket_distrib_ptr1,0)
process_bucket_distrib_ptr2 = Array{Int32}(NUM_BUCKETS+TEST_ARRAY_SIZE)
fill!(process_bucket_distrib_ptr2,0)
send_count = Array{Int32}(MAX_PROCS)
fill!(send_count,0)
recv_count = Array{Int32}(MAX_PROCS)
fill!(recv_count,0)
send_displ = Array{Int32}(MAX_PROCS)
fill!(send_displ,0)
recv_displ = Array{Int32}(MAX_PROCS)
fill!(recv_displ,0)

#/**********************/
#/* Partial verif info */
#/**********************/
test_index_array = Array{Int64}(TEST_ARRAY_SIZE)
fill!(test_index_array,0)
test_rank_array = Array{Int64}(TEST_ARRAY_SIZE)
fill!(test_rank_array,0)

S_test_index_array = Int64[48427,17148,23627,62548,4431]
S_test_rank_array = Int64[0,18,346,64917,65463]

W_test_index_array = Int64[357773,934767,875723,898999,404505]
W_test_rank_array = Int64[1249,11698,1039987,1043896,1048018]

A_test_index_array = Int64[2112377,662041,5336171,3642833,4250760]
A_test_rank_array = Int64[104,17523,123928,8288932,8388264]

B_test_index_array = Int64[41869,812306,5102857,18232239,26860214]
B_test_rank_array = Int64[33422937,10244,59149,33135281,99]

C_test_index_array = Int64[44172927,72999161,74326391,129606274,21736814]
C_test_rank_array = Int64[61147,882988,266290,133997595,133525895]

D_test_index_array = Int64[1317351170,995930646,1157283250,1503301535,1453734525]
D_test_rank_array = Int64[1,36538729,1978098519,2145192618,2147425337]

#=
/*
 *    FUNCTION RANDLC (X, A)
 *
 *  This routine returns a uniform pseudorandom double precision number in the
 *  range (0, 1) by using the linear congruential generator
 *
 *  x_{k+1} = a x_k  (mod 2^46)
 *
 *  where 0 < x_k < 2^46 and 0 < a < 2^46.  This scheme generates 2^44 numbers
 *  before repeating.  The argument A is the same as 'a' in the above formula,
 *  and X is the same as x_0.  A and X must be odd double precision integers
 *  in the range (1, 2^46).  The returned value RANDLC is normalized to be
 *  between 0 and 1, i.e. RANDLC = 2^(-46) * x_1.  X is updated to contain
 *  the new seed x_1, so that subsequent calls to RANDLC using the same
 *  arguments will generate a continuous sequence.
 *
 *  This routine should produce the same results on any computer with at least
 *  48 mantissa bits in double precision floating point data.  On Cray systems,
 *  double precision should be disabled.
 *
 *  David H. Bailey     October 26, 1990
 *
 *     IMPLICIT DOUBLE PRECISION (A-H, O-Z)
 *     SAVE KS, R23, R46, T23, T46
 *     DATA KS/0/
 *
 *  If this is the first call to RANDLC, compute R23 = 2 ^ -23, R46 = 2 ^ -46,
 *  T23 = 2 ^ 23, and T46 = 2 ^ 46.  These are computed in loops, rather than
 *  by merely using the ** operator, in order to insure that the results are
 *  exact on all systems.  This code assumes that 0.5D0 is represented exactly.
 */


/*****************************************************************/
/*************           R  A  N  D  L  C             ************/
/*************                                        ************/
/*************    portable random number generator    ************/
/*****************************************************************/
=#
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
  i::Int64 = 0
  j::Int64 = 0

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
  j = trunc(Int64,T1)
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

  j = trunc(Int64,randlc_R23 * T1)
  T2 = j
  Z = T1 - randlc_T23 * T2
  T3 = randlc_T23 * Z + A2 * X2
  j = trunc(Int64,randlc_R46 * T3)
  T4 = j
  X = T3 - randlc_T46 * T4
  return randlc_R46 * X,X,A
end



#=*****************************************************************/
/************   F  I  N  D  _  M  Y  _  S  E  E  D    ************/
/************                                         ************/
/************ returns parallel random number seq seed ************/
/*****************************************************************/

/*
 * Create a random number sequence of total length nn residing
 * on np number of processors.  Each processor will therefore have a
 * subsequence of length nn/np.  This routine returns that random
 * number which is the first random number for the subsequence belonging
 * to processor rank kn, and which is used as seed for proc kn ran # gen.
 */=#

function find_my_seed(kn::Int64,   # my processor rank, 0<=kn<=num procs
                      np::Int64,   # np = num procs
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
    t2,t1,t1 = randlc(t1,t1)
    i += 1
  end

  an = t1
  kk = kn
  t1 = s
  t2 = an

  for i = 1:100
    ik = div(kk,2)
    if 2 * ik != kk
      t3,t1,t2 = randlc(t1,t2)
    end
    if ik == 0
      break
    end
    t3,t2,t2 = randlc(t2,t2)
    kk = ik
  end
  return t1
end

#=/*****************************************************************/
/*************      C  R  E  A  T  E  _  S  E  Q      ************/
/*****************************************************************/=#

function create_seq(seed::Float64,a::Float64)
  x::Float64 = 0.0
  i::Int64 = 0
  k::Int64 = 0

  k = div(MAX_KEY,4)

  for i = 1:NUM_KEYS
    x,seed,a = randlc(seed,a)
    xtmp,seed,a = randlc(seed,a)
    x += xtmp
    xtmp,seed,a = randlc(seed,a)
    x += xtmp
    xtmp,seed,a = randlc(seed,a)
    x += xtmp

    key_array[i] = trunc(Int32, k * x)
  end
end



#/*****************************************************************/
#/*************    F  U  L  L  _  V  E  R  I  F  Y     ************/
#/*****************************************************************/


function full_verify()
  #MYTEST
  #=
  for i = 1:MAX_PROCS
    @printf("%d\n",send_displ[i])
  end
  exit(0)
  =#

  global passed_verification,total_local_keys,total_lesser_keys
  i::Int32 = 0
  j::Int32 = 0
  k = Array(Int32,1)
  last_local_key::Int32 = 0

  TIMER_START(T_VERIFY)

  #/*  Now, finally, sort the keys:  */
  for i = 1:total_local_keys
    key_buff_ptr_global[key_buff2[i] + 1] -= 1
    key_array[key_buff_ptr_global[key_buff2[i] + 1] - total_lesser_keys] = key_buff2[i]
  end
  last_local_key = (total_local_keys<1)? 1 : (total_local_keys);  #problem


  #/*  Send largest key value to next processor  */
  if my_rank > 0
    request = MPI_Irecv(k,
                        my_rank - 1,
                        1000,
                        MPI.COMM_WORLD)
  end
  if my_rank < comm_size - 1
    MPI.Send([key_array[last_local_key]],
              my_rank + 1,
              1000,
              MPI.COMM_WORLD)
  end
  if(my_rank > 0)
    status = MPI.Wait!(request)
  end

#    /*  Confirm that neighbor's greatest key value
#        is not greater than my least key value       */
  j = 0
  if(my_rank > 0 && total_local_keys > 0)
    if(k > key_array[1])
      j += 1
    end
  end


#    /*  Confirm keys correctly sorted: count incorrectly sorted keys, if any */
  for i = 2:total_local_keys
    if(key_array[i - 1] > key_array[i])
      j += 1
    end
  end

  if(j != 0)
    @printf( "Processor %d:  Full_verify: number of keys out of sort: %d\n",
            my_rank, j );
  else
    passed_verification += 1
  end

  TIMER_STOP( T_VERIFY );
end


#/*****************************************************************/
#/*************             R  A  N  K             ****************/
#/*****************************************************************/


function rank(iteration)
  global passed_verification
  i = 0
  k = 0
  shift = MAX_KEY_LOG_2 - NUM_BUCKETS_LOG_2
  key = 0
  bucket_sum_accumulator = 0
  j = 0
  m = 0
  local_bucket_sum_accumulator = 0
  min_key_val = 0
  max_key_val = 0



  TIMER_START( T_RANK )

  #  Iteration alteration of keys
  if(my_rank == 0 )
    key_array[iteration + 1] = iteration;
    key_array[iteration + MAX_ITERATIONS + 1] = MAX_KEY - iteration;
  end
  #MYTEST
  #=
  for i = 1:SIZE_OF_BUFFERS
    @printf("%d\n",key_array[i])
  end
  exit(0)
  =#
  #  Initialize
  for i = 1:NUM_BUCKETS+TEST_ARRAY_SIZE
    bucket_size[i] = 0;
    bucket_size_totals[i] = 0;
    process_bucket_distrib_ptr1[i] = 0;
    process_bucket_distrib_ptr2[i] = 0;
  end


  # Determine where the partial verify test keys are, load into
  # top of array bucket_size
  for i = 1:TEST_ARRAY_SIZE
    if div(test_index_array[i],NUM_KEYS) == my_rank
      bucket_size[NUM_BUCKETS+i] = key_array[test_index_array[i] % NUM_KEYS + 1]
    end
  end



  # Determine the number of keys in each bucket
  for i = 1:NUM_KEYS
    bucket_size[key_array[i] >> shift + 1] += 1;
  end

  #MYTEST
  #=
  for i = 1:NUM_BUCKETS+TEST_ARRAY_SIZE
    @printf("%d\n",bucket_size[i])
  end
  exit(0)
  =#

  # Accumulative bucket sizes are the bucket pointers
  bucket_ptrs[1] = 0;
  for i = 2:NUM_BUCKETS
    bucket_ptrs[i] = bucket_ptrs[i-1] + bucket_size[i-1]
  end
  #MYTEST
  #=
  for i = 1:NUM_BUCKETS
    @printf("%d\n",bucket_ptrs[i])
  end
  exit(0)
  =#

  # Sort into appropriate bucket
  for i = 1:NUM_KEYS
    key = key_array[i]
    key_buff1[(bucket_ptrs[key >> shift + 1] += 1)] = key
  end
  #MYTEST
  #=
  for i = 1:SIZE_OF_BUFFERS
    @printf("%d\n",key_buff1[i])
  end
  exit(0)
  =#

  TIMER_STOP( T_RANK )
  TIMER_START( T_RCOMM )

  # Get the bucket size totals for the entire problem. These
  # will be used to determine the redistribution of keys
  MPI.Allreduce!(bucket_size,
                 bucket_size_totals,
                 MPI.SUM,
                 MPI.COMM_WORLD)

  #MYTEST
  #=
  for i = 1:NUM_BUCKETS+TEST_ARRAY_SIZE
    @printf("%d\n",bucket_size_totals[i])
  end
  exit(0)
  =#
  TIMER_STOP( T_RCOMM )
  TIMER_START( T_RANK )

  #=   Determine Redistibution of keys: accumulate the bucket size totals
      till this number surpasses NUM_KEYS (which the average number of keys
      per processor).  Then all keys in these buckets go to processor 0.
      Continue accumulating again until supassing 2*NUM_KEYS. All keys
      in these buckets go to processor 1, etc.  This algorithm guarantees
      that all processors have work ranking; no processors are left idle.
      The optimum number of buckets, however, does not result in as high
      a degree of load balancing (as even a distribution of keys as is
      possible) as is obtained from increasing the number of buckets, but
      more buckets results in more computation per processor so that the
      optimum number of buckets turns out to be 1024 for machines tested.
      Note that process_bucket_distrib_ptr1 and ..._ptr2 hold the bucket
      number of first and last bucket which each processor will have after
      the redistribution is done.                                          =#

  bucket_sum_accumulator = 0
  local_bucket_sum_accumulator = 0
  send_displ[1] = 0
  process_bucket_distrib_ptr1[1] = 0
  j = 1
  for i = 1:NUM_BUCKETS
    bucket_sum_accumulator       += bucket_size_totals[i]
    local_bucket_sum_accumulator += bucket_size[i]
    if bucket_sum_accumulator >= j * NUM_KEYS
      send_count[j] = local_bucket_sum_accumulator
      if j != 1
        send_displ[j] = send_displ[j-1] + send_count[j-1]
        process_bucket_distrib_ptr1[j] = process_bucket_distrib_ptr2[j-1]+1
      end
      process_bucket_distrib_ptr2[j] = i - 1
      j += 1
      local_bucket_sum_accumulator = 0
    end
  end
  #MYTEST
  #=
  for i = 1:NUM_BUCKETS+TEST_ARRAY_SIZE
    @printf("%d\n",process_bucket_distrib_ptr1[i])
  end
  exit(0)
  =#

  # When NUM_PROCS approaching NUM_BUCKETS, it is highly possible
  # that the last few processors don't get any buckets.  So, we
  #  need to set counts properly in this case to avoid any fallouts.
  while j <= comm_size
    send_count[j] = 0
    process_bucket_distrib_ptr1[j] = 1
    j += 1
  end
  #MYTEST
  #=
  for i = 1:MAX_PROCS
    @printf("%d\n",send_count[i])
  end
  exit(0)
  =#

  TIMER_STOP( T_RANK )
  TIMER_START( T_RCOMM )

  # This is the redistribution section:  first find out how many keys
  #    each processor will send to every other processor:
  recv_countbuf = MPI.Alltoall(send_count,1,MPI.COMM_WORLD)
  for i = 1:length(recv_countbuf)
    recv_count[i] = recv_countbuf[i]
  end
  #MYTEST
  #=
  for i = 1:length(recv_count)
    @printf("%d\n",recv_count[i])
  end
  exit(0)
  =#

  # Determine the receive array displacements for the buckets
  recv_displ[1] = 0
  for i = 2:comm_size
    recv_displ[i] = recv_displ[i-1] + recv_count[i-1]
  end


  # Now send the keys to respective processors
  key_buff2buff = MPI.Alltoallv(key_buff1,
                                send_count,
                                recv_count,
                                MPI.COMM_WORLD)
  for i = 1:length(key_buff2buff)     #maybe should be deleted later
    key_buff2[i] = key_buff2buff[i]
  end
  #MYTEST
  #=
  for i = 1:length(key_buff2)
    @printf("%d\n",key_buff2[i])
  end
  exit(0)
  =#

  TIMER_STOP( T_RCOMM )
  TIMER_START( T_RANK )

  # The starting and ending bucket numbers on each processor are
  # multiplied by the interval size of the buckets to obtain the
  # smallest possible min and greatest possible max value of any
  # key on each processor
  min_key_val = process_bucket_distrib_ptr1[my_rank + 1] << shift;
  max_key_val = ((process_bucket_distrib_ptr2[my_rank + 1] + 1) << shift)-1
  #MYTEST
  #=
  #@printf("%d\n",process_bucket_distrib_ptr2[my_rank + 1] + 1)
  @printf("%d\n",min_key_val)
  @printf("%d\n",max_key_val)
  exit(0)
  =#

  # Clear the work array
  for i = 1:max_key_val-min_key_val+1
    key_buff1[i] = 0
  end
  #MYTEST
  #=
  for i = 1:length(key_buff1)
    @printf("%d\n",key_buff1[i])
  end
  exit(0)
  =#

  # Determine the total number of keys on all other
  # processors holding keys of lesser value
  m = 0
  for k = 1:my_rank
    for i = process_bucket_distrib_ptr1[k] : process_bucket_distrib_ptr2[k]
      m += bucket_size_totals[i + 1]  # m has total # of lesser keys
    end
  end
  #MYTEST
  #=
  @printf("%d\n",m)
  exit(0)
  =#

  # Determine total number of keys on this processor
  j = 0
  for i = process_bucket_distrib_ptr1[my_rank + 1] : process_bucket_distrib_ptr2[my_rank + 1]
    j += bucket_size_totals[i + 1];     # j has total number of local keys
  end
  #MYTEST
  #=
  @printf("%d\n",j)
  exit(0)
  =#

  # Ranking of all keys occurs in this section:
  # shift it backwards so no subtractions are necessary in loop
  #key_buff_ptr = key_buff1 - min_key_val
  key_buff_ptr = key_buff1  #TODO

  # In this section, the keys themselves are used as their
  # own indexes to determine how many of each there are: their
  # individual population
  for i = 1:j
    key_buff_ptr[key_buff2[i] + 1 - min_key_val] += 1 # Now they have individual key population
  end
  #MYTEST
  #=
  if my_rank == 0
    @printf("%d\n",j)
    for i = 1:SIZE_OF_BUFFERS
      @printf("%d\n",key_buff1[i])
    end
  end
  MPI.Barrier(MPI.COMM_WORLD)
  exit(0)
  =#

  # To obtain ranks of each key, successively add the individual key
  # population, not forgetting the total of lesser keys, m.
  # NOTE: Since the total of lesser keys would be subtracted later
  # in verification, it is no longer added to the first key population
  # here, but still needed during the partial verify test.  This is to
  # ensure that 32-bit key_buff can still be used for class D.
  # key_buff_ptr[min_key_val] += m
  for i = min_key_val + 1:max_key_val
    key_buff_ptr[i+1] += key_buff_ptr[i]
  end
  #MYTEST
  #=
  if my_rank == 0
    for i = 1:TEST_ARRAY_SIZE + NUM_BUCKETS
      @printf("%d\n",bucket_size_totals[i])
    end
  end
  MPI.Barrier(MPI.COMM_WORLD)
  exit(0)
  =#

  # This is the partial verify test section
  # Observe that test_rank_array vals are
  # shifted differently for different cases
  for i = 1:TEST_ARRAY_SIZE
    k = bucket_size_totals[i+NUM_BUCKETS]     # Keys were hidden here
    if min_key_val <= k  &&  k <= max_key_val
      # Add the total of lesser keys, m, here
      #MYTEST
      #=
      if my_rank == 0
        println("k = ",k)
      end
      =#
      key_rank = key_buff_ptr[k] + m
      failed = 0

      if CLASS == 'S'
        if i <= 3
          if key_rank != test_rank_array[i]+iteration
            failed = 1
          else
            passed_verification += 1
          end
        else
          if key_rank != test_rank_array[i]-iteration
            failed = 1
          else
            passed_verification += 1
          end
        end
      elseif CLASS == 'W'
        if i < 3
          if key_rank != test_rank_array[i]+(iteration-2)
            failed = 1
          else
            passed_verification += 1
          end
        else
          if key_rank != test_rank_array[i]-iteration
            failed = 1
          else
            passed_verification += 1
          end
        end
      elseif CLASS == 'A'
        if i <= 3
          if key_rank != test_rank_array[i]+(iteration-1)
            failed = 1
          else
            passed_verification += 1
          end
        else
          if key_rank != test_rank_array[i]-(iteration-1)
            failed = 1
          else
            passed_verification += 1
          end
        end
      elseif CLASS == 'B'
        if i == 1 || i == 2 || i == 4
          if key_rank != test_rank_array[i]+iteration
            failed = 1
          else
            passed_verification += 1
          end
        else
          if key_rank != test_rank_array[i]-iteration
            failed = 1
          else
            passed_verification += 1
          end
        end
      elseif CLASS == 'C'
        if i <= 3
          if key_rank != test_rank_array[i]+iteration
            failed = 1
          else
            passed_verification += 1
          end
        else
          if key_rank != test_rank_array[i]-iteration
            failed = 1
          else
            passed_verification += 1
          end
        end
      elseif CLASS == 'D'
        if i < 3
          if key_rank != test_rank_array[i]+iteration
            failed = 1
          else
            passed_verification += 1
          end
        else
          if key_rank != test_rank_array[i]-iteration
            failed = 1
          else
            passed_verification += 1
          end
        end
      end
      if failed == 1
        @printf("Failed partial verification: iteration %d, processor %d, test key %d\n",
                 iteration, my_rank, i)
      end
    end
  end


  TIMER_STOP( T_RANK )


# Make copies of rank info for use by full_verify: these variables
# in rank are local; making them global slows down the code, probably
# since they cannot be made register by compiler

  if iteration == MAX_ITERATIONS
    key_buff_ptr_global = key_buff_ptr
    total_local_keys    = j
    total_lesser_keys   = 0;  # no longer set to 'm', see note above
  end

end

function debuginfo(rank,ptr)
  if my_rank == rank
    fopt = open("/home/cq/jopt","w+")
    for i = 1:length(ptr)
      atp = ptr[i]
      write(fopt,"$atp\n")
    end
    close(fopt)
  end
  MPI.Barrier(MPI.COMM_WORLD)
  MPI.Finalize()
  exit(0)
end
#*****************************************************************
#*************             M  A  I  N             ****************
#*****************************************************************

function main()
  global timeron,passed_verification
  i = 0
  iteration = 0
  itemp = 0

  timecounter = 0.0
  maxtime = 0.0

# Initialize MPI
  MPI.Init()
  my_rank = MPI.Comm_rank(MPI.COMM_WORLD)
  comm_size = MPI.Comm_size(MPI.COMM_WORLD)


# Initialize the verification arrays if a valid class
  for i = 1:TEST_ARRAY_SIZE
    if CLASS == 'S'
      test_index_array[i] = S_test_index_array[i]
      test_rank_array[i]  = S_test_rank_array[i]
    elseif CLASS == 'A'
      test_index_array[i] = A_test_index_array[i]
      test_rank_array[i]  = A_test_rank_array[i]
    elseif CLASS == 'W'
      test_index_array[i] = W_test_index_array[i]
      test_rank_array[i]  = W_test_rank_array[i]
    elseif CLASS == 'B'
      test_index_array[i] = B_test_index_array[i]
      test_rank_array[i]  = B_test_rank_array[i]
    elseif CLASS == 'C'
      test_index_array[i] = C_test_index_array[i]
      test_rank_array[i]  = C_test_rank_array[i]
    elseif CLASS == 'D'
      test_index_array[i] = D_test_index_array[i]
      test_rank_array[i]  = D_test_rank_array[i]
    end
  end



# Printout initial NPB info
  if my_rank == 0
    @printf( "\n\n NAS Parallel Benchmarks 3.3 -- IS Benchmark\n\n" )
    @printf( " Size:  %ld  (class %c)\n", TOTAL_KEYS*MIN_PROCS, CLASS )
    @printf( " Iterations:   %d\n", MAX_ITERATIONS )
    @printf( " Number of processes:     %d\n", comm_size )

    timeron = 0
    if isfile("timer.flag") == true
      timeron = 1
    end
  end

# Check that actual and compiled number of processors agree
  if comm_size != NUM_PROCS
    if my_rank == 0
      @printf("\n ERROR: compiled for %d processes\n Number of active processes: %d\n Exiting program!\n\n",
              NUM_PROCS,
              comm_size)
    end
    MPI.Finalize()
    exit(1)
  end

# Check to see whether total number of processes is within bounds.
# This could in principle be checked in setparams.c, but it is more
# convenient to do it here
  if comm_size < MIN_PROCS || comm_size > MAX_PROCS
    if my_rank == 0
      @printf( "\n ERROR: number of processes %d not within range %d-%d\n Exiting program!\n\n",
               comm_size,
               MIN_PROCS,
               MAX_PROCS);
    end
    MPI.Finalize()
    exit(1)
  end

  timeron = MPI.bcast(timeron, 0, MPI.COMM_WORLD)

  if TIMING_ENABLED == 0
    for i = 1:T_LAST
      timer_clear(i)
    end
  end

# Generate random number sequence and subsequent keys on all procs
  create_seq( find_my_seed(my_rank,
                           comm_size,
                           4*TOTAL_KEYS*MIN_PROCS,
                           314159265.00,      # Random number gen seed
                           1220703125.00 ),   # Random number gen mult
              1220703125.00 )                 # Random number gen mult
  debuginfo(1,key_array)

# Do one interation for free (i.e., untimed) to guarantee initialization of
# all data and code pages and respective tables
  rank(1)

# Start verification counter
  passed_verification = 0;
  if my_rank == 0 && CLASS != 'S'
    @printf( "\n   iteration\n" );
  end

# Initialize timer
  timer_clear(0)

# Initialize separate communication, computation timing
  if TIMING_ENABLED == 0
    for i = 1:T_LAST
      timer_clear(i)
    end
  end

# Start timer
  timer_start(0)


# This is the main iteration
  for iteration = 1:MAX_ITERATIONS
    if my_rank == 0 && CLASS != 'S'
      @printf( "        %d\n", iteration )
    end
    rank( iteration );
  end
# Stop timer, obtain time for processors
  timer_stop( 0 )

  timecounter = timer_read( 0 )

# End of timing, obtain maximum time of all processors
  maxtime = MPI.Reduce(timecounter,
                       MPI.MAX,
                       0,
                       MPI.COMM_WORLD)

# This tests that keys are in sequence: sorting of last ranked key seq
# occurs here, but is an untimed operation
  full_verify()


# Obtain verification counter sum
  itemp = passed_verification;
  passed_verification = MPI.Reduce(itemp,
                                   MPI.SUM,
                                   0,
                                   MPI.COMM_WORLD)
  #MYTEST
  #=
  @printf("%d\n",itemp)
  exit(0)
  =#

# The final printout
  if my_rank == 0
    if passed_verification != 5*MAX_ITERATIONS + comm_size
      passed_verification = 0
    end
    #=
    c_print_results( "IS",
                     CLASS,
                     TOTAL_KEYS,
                     MIN_PROCS,
                     0,
                     MAX_ITERATIONS,
                     NUM_PROCS,
                     comm_size,
                     maxtime,
                     ( (MAX_ITERATIONS)*TOTAL_KEYS*MIN_PROCS)
                                                  /maxtime/1000000.,
                     "keys ranked",
                     passed_verification,
                     NPBVERSION,
                     COMPILETIME,
                     MPICC,
                     CLINK,
                     CMPI_LIB,
                     CMPI_INC,
                     CFLAGS,
                     CLINKFLAGS );  #TODO
    =#
    c_print_results( "IS",
                     CLASS,
                     TOTAL_KEYS,
                     MIN_PROCS,
                     0,
                     MAX_ITERATIONS,
                     NUM_PROCS,
                     comm_size,
                     maxtime,
                     ( (MAX_ITERATIONS)*TOTAL_KEYS*MIN_PROCS)
                                                  /maxtime/1000000.,
                     "keys ranked",
                     passed_verification,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0,
                     0)
  end


  if TIMING_ENABLED == 0
    if timeron == 1
      t1 = Array(Float64,T_LAST + 1)
      tmin = Array(Float64,T_LAST + 1)
      tsum = Array(Float64,T_LAST + 1)
      tmax = Array(Float64,T_LAST + 1)
      t_recs = Array(Char,T_LAST + 1,9)

      for i = 1:T_LAST + 1
        t1[i] = timer_read( i )
      end

      MPI.Reduce!(t1,
                  tmin,
                  MPI.MIN,
                  0,
                  MPI.COMM_WORLD)
      MPI.Reduce!(t1,
                  tsum,
                  MPI.SUM,
                  0,
                  MPI.COMM_WORLD)
      MPI.Reduce!(t1,
                  tmax,
                  MPI.MAX,
                  0,
                  MPI.COMM_WORLD)

      if my_rank == 0
        #TODO
      end
    end
  end
  MPI.Finalize()

  return 0

end




TIMING_ENABLED = -1
main()











#=
#  test of find_my_seed

kn = 0
np = 10
nn = 10
s = 314159265.00
a = 1220703125.00
b = find_my_seed(kn,np,nn,s,a)
println("b:",b)
=#

#=
#  test of create_seq
create_seq(123.0,345.0)
for i = 1: NUM_KEYS
  println(key_array[i])
end
=#
