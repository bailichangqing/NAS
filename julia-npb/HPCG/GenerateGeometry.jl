module GenerateProblem
using Primes
importall Base
type Counter
  length::Int # number of prime factor counts(cannot exceed 32 or a 32-bit integer)
  max_counts::Array{Int,1} #maximum value for prime factor counts
  cur_counts::Array{Int,1} #current prime factor counts
end
function counternew(Counter::this, Array{Int64,1}::counts, Int64::length )
  for i = 1:32
    this.length = length
    push!(this.max_counts, counts[i])
    push!(this.cur_counts,0)
  end
  push!(this.max_counts, 0)
  push!(this.cur_counts, 0)
  this.max_counts[length] = this.cur_counts[length] = 0
end

function counternext(Counter::this)
   i = Int64
for i :this.length
    this.cur_counts[i]+= 1
    if (this.cur_counts[i] > this.max_counts[i])
      this.cur_counts[i] = 0
      end
    end
end

function counteriszero(Counter::this)
  i = Int64
  for i : this.length
    if (this.cur_counts[i])
      return 0
    else
      return 1
  end
end

function counterproduct(Counter::this , Array{Int64,1}::multipliers)
  i, j = Int64
  k=0, x=1

  for (i : this.length)
    for (j : this.cur_counts[i])
      k = 1
      x *= multipliers[i]
      end
  end

  return x * k
end

function countermaxcursub(Counter::this, Counter::that, Counter::res) {
  i = Int64;

  res.length = length;
  for (i = 0; i < res.length; ++i) {
    res.max_counts[i] = this.max_counts[i] - that.cur_counts[i];
    res.cur_counts[i] = 0;
  }
}

function
primefactor_i(x ::Int64, factors::Array{Int64,1})
  i = 0
  d = 3
  r = Int64
  sq=sqrt(x)
  #div_t r; julia has div(x,y) for a truncated result

  #remove 2 as a factor with shifts
  for ( x > 1 && (x & 1) == 0: x >>>= 1)
    factors[i+=1] = 2
  end

  #keep removing subsequent odd numbers
  for ( d <= sq : d += 2)
    while (1)
      r = div(x, d)
      if (x%d == 0)
        push!(factors, d)
        x = r
        continue
      end
      break
    end
  end
  if (x > 1 || i == 0)  # left with a prime or x==1 */
    push!(factors,x)
  end
  push!(factors,0) #/* terminate with 0 */
end
function
gen_min_area3(n::Int64, f1::Array{Int64,1}, f2::Array{Int64,1}, f3::Array{Int64,1})
  i = j = df_cnt = tf1 = tf2  = tf3 = Int64
  factors = distinct_factors = count_factors = Array{Int64, 33}
  c_main = c1 = c2 = Counter

  #at the beginning, minimum area is the maximum area */
  area = min_area = 2.0 * n + 1.0;

  primefactor_i( n, factors ) # /* factors are sorted: ascending order */

  if (1 == n || factors[1] == 0) #/* prime number */
    f1 = n
    f2 = 1
    f3 = 1
    return
  elseif (factors[2] == 0) #/* two prime factors */
    f1 = factors[0]
    f2 = factors[1]
    f3 = 1
    return
  elseif (factors[3] == 0) # /* three prime factors */
    f1 = factors[0]
    f2 = factors[1]
    f3 = factors[2]
    return
      end

  #/* we have more than 3 prime factors so we need to try all possible combinations */
j= 0
  for (i = 0: factors[i])
    push!(distinct_factors,factors[i])
    count_factors[j-1] = 0
    while (distinct_factors[j-1] == factors[++i])
    count_factors[j-1]+=1
    end
 # USING FACTOR FUNCTION FROM PRIMES.JL
  # factors =Int64[]
       #k =  collect(keys(factor(n)))
       #v = collect(values(factor(n)))
       #for i = 1:length(k)
        #            for j = 1:v[i]
        #            push!(factors,k[i])
        #            end
        #            end
    #push!(factors,0)
  df_cnt = length(count_factors)

  Counter_new( c_main, count_factors, df_cnt )

  Counter_new( c1, count_factors, df_cnt )

  for ( Counter_next( c1 ), ! Counter_is_zero( c1 ))

    Counter_max_cur_sub( c_main, &c1, &c2 );
    for ( Counter_next( c2 ), ! Counter_is_zero( c2 ))
      tf1 = Counter_product( c1, distinct_factors )
      tf2 = Counter_product( c2, distinct_factors )
      tf3 = n / tf1/ tf2

      area = tf1 * tf2 + tf2 * tf3 + tf1 * tf3
      if (area < min_area)
        min_area = area;
        *f1 = tf1;
        *f2 = tf2;
        *f3 = tf3;
              end
            end
          end
      end

#=
  Computes the factorization of the total number of processes into a
  3-dimensional process grid that is as close as possible to a cube. The
  quality of the factorization depends on the prime number structure of the
  total number of processes. It then stores this decompostion together with the
  parallel parameters of the run in the geometry data structure.

  @param[in]  size total number of MPI processes
  @param[in]  rank this process' rank among other MPI processes
  @param[in]  numThreads number of OpenMP threads in this process
  @param[in]  nx, ny, nz number of grid points for each local block in the x, y, and z dimensions, respectively
  @param[out] geom data structure that will store the above parameters and the factoring of total number of processes into three dimensions
=#
function GenerateGeometry(size::Int64,rank::Int64, numThreads::Int64, nx::Int64, ny::Int64, int::Int64, geom::Geoemtry)

  npx::Int64
  npy::Int64
  npz::Int64

  gen_min_area3( size, npx, npy, npz )

  #Now compute this process's indices in the 3D cube
  ipz::Int64 = div(rank, npx*npy)
  ipy::Int64 = div(rank-ipz*npx*npy, npx)
  ipx::Int64 = mod(rank, npx)

#if isdefined(HPCG_DEBUG)  == true
  if (rank==0)
        #had HPCG_fout commands
    println("size= %d", size)
    println( "nx  = %d \n ny  = %d \n nz  = %d \n npx = %d \n npy = %d \n npz = %d", nx, ny, nz, npx, npy, npz)
    println("For rank = %d \n ", rank)
    println("ipx = %d \n ipy = %d \n ipy = %d \n ipz = %d \n", ipx, ipy, ipz)
    #if (size == npx*npy*npz = 0)
         #assert(size==npx*npy*npz);
        #how do i say get out of here

  geom.size = size
  geom.rank = rank
  geom.numThreads = numThreads
  geom.nx = nx
  geom.ny = ny
  geom.nz = nz
  geom.npx = npx
  geom.npy = npy
  geom.npz = npz
  geom.ipx = ipx
  geom.ipy = ipy
  geom.ipz = ipz
      end


