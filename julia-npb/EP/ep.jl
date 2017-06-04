#-------------------------------------------------------------------------#
#                                                                         #
#        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3         #
#                                                                         #
#                                   E P                                   #
#                                                                         #
#-------------------------------------------------------------------------#
#                                                                         #
#    This benchmark is part of the NAS Parallel Benchmark 3.3 suite.      #
#    It is described in NAS Technical Reports 95-020 and 02-007           #
#                                                                         #
#    Permission to use, copy, distribute and modify this software         #
#    for any purpose with or without fee is hereby granted.  We           #
#    request, however, that all derived work reference the NAS            #
#    Parallel Benchmarks 3.3. This software is provided "as is"           #
#    without express or implied warranty.                                 #
#                                                                         #
#    Information on NPB 3.3, including the technical report, the          #
#    original specifications, source code, results and information        #
#    on how to submit new results, is available at:                       #
#                                                                         #
#           http://www.nas.nasa.gov/Software/NPB/                         #
#                                                                         #
#    Send comments or suggestions to  npb@nas.nasa.gov                    #
#                                                                         #
#          NAS Parallel Benchmarks Group                                  #
#          NASA Ames Research Center                                      #
#          Mail Stop: T27A-1                                              #
#          Moffett Field, CA   94035-1000                                 #
#                                                                         #
#          E-mail:  npb@nas.nasa.gov                                      #
#          Fax:     (650) 604-3957                                        #
#                                                                         #
#-------------------------------------------------------------------------#

include("../common/timers.jl")
include("../common/print_results.jl")
include("./npbparams.jl")

import MPI
import Timers
import PrintResults

const mk = 16
const mm = m - mk
const nn = 2^mm
const nk = 2^mk
const nq = 10
const epsilon = 1.0e-8
const a = 1220703125.0
const s = 271828183.0

const t_total = 1
const t_gpairs = 2
const t_randn = 3
const t_rcomm = 4
const t_last = 4

tsum  = Array{Float64}(t_last+2)
t1m   = Array{Float64}(t_last+2)
tming = Array{Float64}(t_last+2)
tmaxg = Array{Float64}(t_last+2)

# 
#  Since we're using gotos here, we need to 
#  wrap everything in a big function
# 
function main()

	MPI.Init()

	node     = MPI.Comm_rank(MPI.COMM_WORLD)
	no_nodes = MPI.Comm_size(MPI.COMM_WORLD)

	root           = 0
	timers_enabled = false

	dum  = [1.0, 1.0, 1.0]
	x = Array{Float64}(2*nk)
	q = Array{Float64}(nq)

	t_recs = ["total", "gpairs", "randn", "rcomm", "totcomp",  " totcomm"]

	if node == root
		@printf("\n NAS Parallel Benchmarks 3.3 (Julia) -- EP Benchmark\n\n")
		@printf(" Number of random numbers generated: %15.0f\n", 2.0^(m+1))
		@printf(" Number of active processes:         %13d\n\n",no_nodes)

		if isfile("timer.flag")
			timers_enabled = true
		end
	end

	verified = false

	dum[1] = NPBRand.vranlc!(0, dum[1], dum[2], dum)
	for i in 1:20
		res, newx = NPBRand.randlc(1.0, 1.0)
		@printf("rand[%d] = %15.0f\n", i, res)
	end


	MPI.bcast(timers_enabled, root, MPI.COMM_WORLD)

	np             = nn / no_nodes
	no_large_nodes = nn % no_nodes

	if node == root
		@printf("np = %d, no_large_nodes = %d\n", np, no_large_nodes)
	end
	
	if node < no_large_nodes
		np_add = 1
	else
		np_add = 0
	end

	np = np + np_add

	if np == 0 
		@printf("Too many nodes: %6d %6d", no_nodes, nn)
		MPI.Abort(MPI.COMM_WORLD)
		exit(1)
	end

	#   Call the random number generator functions and initialize
	#   the x-array to reduce the effects of paging on the timings.
	#   Also, call all mathematical functions that are used. Make
	#   sure these initializations cannot be eliminated as dead code.

	dum[1] = NPBRand.vranlc!(0, dum[1], dum[2], dum)
	dum[1], dum[2] = NPBRand.randlc(dum[2], dum[3])

	if node == root
		@printf("dum[1] is %f\n", dum[1])
	end

	for i in 1:2*nk
		x[i] = -1.0e99
	end

		if node == root
			for i in 1:nq
				@printf("x[%d] = %10.10f\n", i, x[i])
			end
		end

	Mops = log(sqrt(abs(max(1.0, 1.0))))

	if node == root
		@printf("Mops is %f\n", Mops)
	end

	#---------------------------------------------------------------------
	#      Synchronize before placing time stamp
	#---------------------------------------------------------------------
	for i in 1:t_last
		Timers.timer_clear(i)
	end

	MPI.Barrier(MPI.COMM_WORLD)
	Timers.timer_start(1)

	t1 = a

	t1 = NPBRand.vranlc!(0, t1, a, x)

	#   Compute AN = A ^ (2 * NK) (mod 2^46).

	t1 = a

	for i in 1:(mk+1)
		t2, t1 = NPBRand.randlc(t1, t1)
	end

	an = t1
	tt = s
	gc = 0.0
	sx = 0.0
	sy = 0.0

	@printf("an is %f\n", an)

	for i in 1:nq
		q[i] = 0.0
	end

	#   Each instance of this loop may be performed independently. We compute
	#   the k offsets separately to take into account the fact that some nodes
	#   have more numbers to generate than others

	if np_add == 1
		k_offset = node * np - 1
	else
		k_offset = no_large_nodes*(np+1) + (node-no_large_nodes)*np -1
	end

	for k in 1:np
		kk = k_offset + k
		t1 = s
		t2 = an

		# Find starting seed t1 for this kk.
		
		for i in 1:100
			ik = kk / 2
			if 2*ik != kk
				t3, t1 = NPBRand.randlc(t1, t2)
			end
		
			if (ik == 0) 
				@goto compute_uniform
			end

			@printf("calling randlc with %f %f\n", t2, t2)
			t3, t2 = NPBRand.randlc(t2, t2)

			kk = ik
		end

		# Compute uniform pseudorandom numbers.
		@label compute_uniform

		if timers_enabled == true
			Timers.timer_start(t_randn)
		end

		if node == root
			@printf("t1 is %f, t2 is %f, s is %f\n", t1, t2, s)
		end
		t1 = NPBRand.vranlc!(2 * nk, t1, a, x)
		if node == root
			@printf("t1 returned is %f\n", t1)
		end

		if timers_enabled == true
			Timers.timer_stop(t_randn)
		end

		if node == root
			for i in 1:nq
				@printf("x[%d] = %10.10f\n", i, x[i])
			end
		end

		# Compute Gaussian deviates by acceptance-rejection method and 
		# tally counts in concentric square annuli.  This loop is not 
		# vectorizable. 

		if timers_enabled == true
			Timers.timer_start(t_gpairs)
		end

		for i in 1:nk
			x1 = 2.0 * x[2*i-1] - 1.0
			x2 = 2.0 * x[2*i] - 1.0
			t1 = x1^2 + x2^2
			if t1 <= 1.0
				t2 = sqrt(-2.0 * log(t1) / t1)
				t3 = (x1 * t2)
				t4 = (x2 * t2)
				l = Int(round(max(abs(t3), abs(t4)))) + 1  # KCH added line here, we can't have arbitrary indexing in julia (yet)
				q[l] += 1.0
				sx += t3
				sy += t4
			end
		end

		if timers_enabled == true
			Timers.timer_stop(t_gpairs)
		end

	end # ! for k in 1,np

	if timers_enabled == true
		Timers.timer_start(t_rcomm)
	end

	sxa = [sx]
	sya = [sy]
	MPI.Allreduce!(sxa, x, 1, MPI.SUM, MPI.COMM_WORLD)
	sx = x[1]
	MPI.Allreduce!(sya, x, 1, MPI.SUM, MPI.COMM_WORLD)
	sy = x[1]
	MPI.Allreduce!(q, x, nq, MPI.SUM, MPI.COMM_WORLD)

	if timers_enabled == true
		Timers.timer_stop(t_rcomm)
	end

	if node == root
		@printf("after allreduce, sx is %f, sy is %f\n", sx, sy)
	end

	for i in 1:nq
		q[i] = x[i]
	end

	for i in 1:nq
		gc += q[i]
	end

	Timers.timer_stop(1)
	tm = Timers.timer_read(1)
	tma = [tm]

	MPI.Allreduce!(tma, x, 1, MPI.MAX, MPI.COMM_WORLD)

	tm = x[1]

	if node == root
		nit = 0
		verified = true
		 if m == 24 
			sx_verify_value = -3.247834652034740e3
			sy_verify_value = -6.958407078382297e3
		 elseif m == 25
			sx_verify_value = -2.863319731645753e3
			sy_verify_value = -6.320053679109499e3
		 elseif m == 28
			sx_verify_value = -4.295875165629892e3
			sy_verify_value = -1.580732573678431e4
		 elseif m == 30
			sx_verify_value =  4.033815542441498e4
			sy_verify_value = -2.660669192809235e4
		 elseif m == 32
			sx_verify_value =  4.764367927995374e4
			sy_verify_value = -8.084072988043731e4
		 elseif m == 36
			sx_verify_value =  1.982481200946593e5
			sy_verify_value = -1.020596636361769e5
		 elseif m == 40
			sx_verify_value = -5.319717441530e5
			sy_verify_value = -3.688834557731e5
		 else
			verified = false
		 end

		if verified == true 
			sx_err = abs((sx - sx_verify_value)/sx_verify_value)
			sy_err = abs((sy - sy_verify_value)/sy_verify_value)
			verified = ((sx_err < epsilon) && (sy_err < epsilon))
		end

		Mops = 2.0^(m+1) / tm / 1000000.0

		@printf("""EP Benchmark Results:\n\nCPU Time=%10.4f
N = 2^%5d
No. Gaussian Pairs=%15.0f
Sums = %5.15f %5.15f
Counts:
""",
			tm, m, gc, sx*10, sy*10)

		for i in 1:nq
			@printf("%3d %15.0f\n", i-1, q[i])
		end

		PrintResults.print_results("EP", cclass, m+1, 0, 0, nit, npm, no_nodes, tm, Mops, "Random numbers generated", verified, npbversion, rand)

	end # !node == root

	if timers_enabled == false
		@goto do_finalize
	end


println("tlast is $t_last")
println("reading timers")
	for i in 1:t_last
		t1m[i] = Timers.timer_read(i)
	end

println("timer read done")
	t1m[t_last+2] = t1m[t_rcomm]
	t1m[t_last+1] = t1m[t_total] - t1m[t_last+2]

println("reducing")
	tsum  = MPI.Reduce(t1m, t_last+2, MPI.SUM, 0, MPI.COMM_WORLD)
println("reducing")
	tming = MPI.Reduce(t1m, t_last+2, MPI.MIN, 0, MPI.COMM_WORLD)
println("reducing")
	tmaxg = MPI.Reduce(t1m, t_last+2, MPI.MAX, 0, MPI.COMM_WORLD)
println("reductions done")

	if node == 0
		@printf(" nprocs =%6d           minimum     maximumm     average\n", no_nodes)
		for i in 1:(t_last+2)
 			tsum[i] = tsum[i] / no_nodes
			@printf(" timer %2d(%8s) :   %10.4f  %10.4f  %10.4f\n", i, t_recs[i], tming[i], tmaxg[i], tsum[i])
		end
	end

	@label do_finalize
	MPI.Finalize()
end

main()


