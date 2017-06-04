#
#  This utility configures a NPB to be built for a specific number
#  of nodes and a specific class. It creates a file "npbparams.jl" 
#  in the source directory. This file keeps state information about 
#  which size of benchmark is currently being built (so that nothing
#  if unnecessarily rebuilt) and defines 
#  the number of nodes and class for which a 
#  benchmark is being built. 
#
#  The utility takes 3 arguments: 
#        ./setparams.jl benchmark-name nprocs class
#     benchmark-name is "sp", "bt", etc
#     nprocs is the number of processors to run on
#     class is the size of the benchmark
#  These parameters are checked for the current benchmark. If they
#  are invalid, this program prints a message and aborts. 
#  If the parameters are ok, the current npbparams.jl (actually just
#  the first line) is read in. If the new parameters are the same as 
#  the old, nothing is done, but an exit code is returned to force the
#  user to specify (otherwise the make procedure succeeds but builds a
#  binary of the wrong name).  Otherwise the file is rewritten. 
#  Errors write a message (to stdout) and abort. 
#  
#  This program makes use of two extra benchmark "classes"
#  class "X" means an invalid specification. It is returned if
#  there is an error parsing the config file. 
#  class "U" is an external specification meaning "unknown class"
#  
#  Unfortunately everything has to be case sensitive. This is
#  because we can always convert lower to upper or v.v. but
#  can't feed this information back to the makefile, so typing
#  make CLASS=a and make CLASS=A will produce different results.
#
# 
#!/usr/bin/env julia

@enum benchmark_types SP BT LU MG FT IS DT EP CG
@enum iotypes NONE FULL SIMPLE EPIO FORTRAN

FILENAME = "npbparams.jl"
VERBOSE  = true
VERSION  = "3.3.1"
DEFFILE  = "../config/make.def"

function get_info(args)

	if length(args) < 3
		@printf("Usage: %s (%d)  benchmark-name nprocs class\n", Base.PROGRAM_FILE, length(args))
		exit(1)
	end

	subtypep = NONE
	nprocs = parse(Int, args[2])
	classp = args[3]

	if args[1] == "sp" || args[1] == "SP"
		typep = SP
	elseif args[1] == "ft" || args[1] == "FT"
		typep = FT
	elseif args[1] == "lu" || args[1] == "LU"
		typep = LU
	elseif args[1] == "mg" || args[1] == "MG"
		typep = MG
	elseif args[1] == "is" || args[1] == "IS"
		typep = IS
	elseif args[1] == "dt" || args[1] == "DT"
		typep = DT
	elseif args[1] == "ep" || args[1] == "EP"
		typep = EP
	elseif args[1] == "cg" || args[1] == "CG"
		typep = CG
	elseif args[1] == "bt" || args[1] == "BT"
		typep = BT
		if length(args) != 4
			subtypep = NONE
		else
			if args[4] == "full" || args[4] == "FULL"
				subtypep = FULL
			elseif args[4] == "simple" || args[4] == "SIMPLE"
				subtypep = SIMPLE
			elseif args[4] == "epio" || args[4] == "EPIO"
				subtypep = EPIO
			elseif args[4] == "fortran" || args[4] == "FORTRAN"
				subtypep = FORTRAN
			elseif args[4] == "none" || args[4] == "NONE"
				subtypep = NONE
			else 
				@printf("setparams: Error: unknown btio type %s\n", args[4])
				exit(1)
			end
		end
	else
		@printf("setparams: Error: unknown benchmark type %s\n", args[1])
		exit(1)
	end

	return (nprocs, classp, typep, subtypep)

end

function npb_isqrt(i)
	if i <= 0
		return -1
	end

	sq = 0
	for root in 1:i
		sq = root*root
		if sq == i
			return root
		end
	end
	return -1
end

function isqrt2(i)
	if i <= 0
		return -1
	end
	sq = 0
	xdim = 1
	while sq <= i
		sq = xdim*xdim
		if sq == i
			return xdim
		end
	end

	xdim -= 1
	ydim = i / xdim
	
	while (xdim*ydim != i) && (2*ydim >= xdim)
		xdim += 1
		ydim = i / xdim
	end

	if (xdim*ydim == i) && (2*ydim >= xdim)
		return xdim
	end
	
	return -1
end

function ilog2(i)
	exp2 = 1
	
	if i <= 0
		return -1
	end

	for log2 in 0:29
		if exp2 == i
			return log2
		end
	
		if exp2 > i
			break
		end

		exp2 *= 2
	end

	return -1
end

function ipow2(i)
	pow2 = 1
	if i < 0
		return -1
	elseif i == 0
		return 1
	else
		i -= 1
		while i != 0
			pow2 *= 2
			i -= 1
		end
	end

	return pow2
end

	

function check_info(tp, np, cp)

	if np <= 0
		println("setparams: Number of processors must be greater than zero")
		exit(1)
	end
	
	if tp == SP || tp == BT
		rootprocs = npb_isqrt(np)
		if rootprocs < 0
			@printf("setparams: Number of processors %d must be a square (1,4,9,...) for this benchmark\n", np)
			exit(1)
		end
		
		if cp == "S" && np > 16
			@printf("setparams: BT and SP sample sizes cannot be run on more\n")
			@printf("           than 16 processors because the cell size would be too small.\n")
			exit(1)
		end
	elseif tp == LU
		rootprocs = isqrt2(nprocs)
		if rootprocs < 0
			@printf("setparams: Failed to determine proc_grid for nprocs=%d\n", np)
			exit(1)
		end
	elseif tp == CG || tp == FT || tp == MG || tp == IS
		logprocs = ilog2(np)
		if logprocs < 0
			@printf("setparams: Number of processors must be a power of two (1,2,4,...) for this benchmark\n")
			exit(1)
		end
	elseif tp == EP || tp == DT
		# do nothing
	else
		# never should have gotten this far with a bad name
		@printf("setparams: (Internal Error) Benchmark type %d unknown to this program\n", tp)
		exit(1)
	end

	if cp != "S" && cp != "W" && cp != "A" && cp != "B" && cp != "C" && cp != "D" && cp != "E"
		@printf("setparams: Unknown benchmark class %s\n", cp)
		@printf("setparams: Allowed classes are \"S\", \"W\", and \"A\" through \"E\"\n")
		exit(1)
	end

	if cp == "E" && (tp == IS || tp == DT)
		@printf("setparams: Benchmark class %s not defined for IS or DT\n", cp)
		exit(1)
	end
	
	if cp == "D" && tp == IS && np < 4
		@printf("setparams: IS class D size cannot be run on less than 4 processors\n")
		exit(1)
	end

end

#
# read_info(): Read previous information from file.
#   Not an error if file doesn't exist, because this
#   may be the first time we're running.
#   Assumes the first line of the file is in a special
#   format that we understand (since we wrote it).
# 
# Input:
#   @tp = type
# Return values:
#   @nprocsp = number of processors
#   @classp  = problem class
#   @subtypep = subtype 
#
function read_info(tp)
	
	subtypep = -1
	classp = "X"
	subtype = -1
	nprocsp = -1

	if !isfile(FILENAME)
		if VERBOSE == true
			@printf("setparams: INFO: configuration file %s does not exist (yet)\n", FILENAME)
		end
		@goto abrt
	end

	fp = open(FILENAME, "r")

	line = readline(fp)
	
	if tp == BT

		m = match(r"# NPROCS = (\d+) CLASS = (\w) SUBTYPE = (\w+)", line)

		if !m
			@printf("setparams: (Internal error): could not parse existing config file\n")
			exit(1)
		end

		matched_procs = m.captures[1]
		matched_class = m.captures[2]
		matched_subtp = m.captures[3]

		if length(m.captures) != 3
			if length(m.captures) != 2
				@printf("setparams: Error parsing config file %s. Ignoring previous settings\n", FILENAME)
				@goto abrt_fp
			end
			subtypep = 0
			@goto clean_finish
		end

		if matched_subtp == "full" || matched_subtp == "FULL"
			subtypep = FULL
		elseif matched_subtp == "simple" || matched_subtp == "SIMPLE"
			subtypep = SIMPLE
		elseif matched_subtp == "epio" || matched_subtp == "EPIO"
			subtypep = EPIO
		elseif matched_subtp == "fortran" || matched_subtp == "FORTRAN"
			subtypep = FORTRAN
		else
			subtypep = -1
		end
	elseif tp == IS || tp == DT || tp == SP || tp == FT || 
		   tp == MG || tp == LU || tp == EP || tp == CG

		m = match(r"# NPROCS = (\d+) CLASS = (\w)", line)
		matched_procs = m.captures[1]
		matched_class = m.captures[2]
		if length(m.captures) != 2
			@printf("setparams: Error parsing config file %s. Ignoring previous settings\n", FILENAME)
			@goto abrt_fp
		end
	else
		@printf("setparams: (Internal Error) Benchmark type %d unknown to this program\n", tp)
		exit(1)
	end
		
@label clean_finish
		close(fp)
		return (nprocsp, classp, subtypep)
				
@label abrt_fp
	close(fp)
@label abrt
	nprocsp = -1
	classp = "X"
	subtypep = -1

	return (nprocsp, classp, subtypep)
end


# KCH TODO
function write_sp_info(fp, np, clp)
	println("ALERT: write_sp_info Not implemented")
end

# KCH TODO
function write_lu_info(fp, np, clp)
	println("ALERT: write_lu_info Not implemented")
end

# KCH TODO
function write_mg_info(fp, np, clp)
	println("ALER: write_mg_info Not implemented")
end

function write_is_info(fp, np, clp)
	if  clp != "S" &&
		clp != "W" &&
		clp != "A" &&
		clp != "B" &&
		clp != "C" &&
		clp != "D"
		@printf("setparams: (Internal error): invalid class type %s\n", clp)
		exit(1)
	end
end

function write_dt_info(fp, np, clp)

	num_samps = 0
	dev       = 0
	num_src   = 0

	if clp == "S"
		num_samps = 1728
		dev       = 128
		num_src   = 4
	elseif clp == "W"
		num_samps = 1728*8
		dev       = 128*2
		num_src   = 4*2
	elseif clp == "A"
		num_samps = 1728*64
		dev       = 128*4
		num_src   = 4*4
	elseif clp == "B"
		num_samps = 1728*512
		dev       = 128*8
		num_src   = 4*8
	elseif clp == "D"
		num_samps = 1728*4096*8
		dev       = 128*32
		num_src   = 4*32
	else
		@printf("setparams: (Internal error): invalid class type %s\n", clp)
		exit(1)
	end

	@printf(fp, "const NUM_SAMPLES   = %d\n", num_samps)
	@printf(fp, "const STD_DEVIATION = %d\n", dev)
	@printf(fp, "const NUM_SOURCES   = %d\n", num_src)

end

function write_ft_info(fp, np, clp)
	# easiest way (given the way the benchmark is written)
	# is to specify log of number of grid points in each
	# direction m1, m2, m3. nt is the number of iterations
	nx = 0
	ny = 0
	nz = 0
	niter = 0

	if class == "S" 
		nx = 64
		ny = 64
		nz = 64
		niter = 6
	elseif class == "W"
		nx = 128
		ny = 128
		nz = 32
		niter = 6
	elseif class == "A"
		nx = 256
		ny = 256
		nz = 128
		niter = 6
	elseif class == "B"
		nx = 512
		ny = 256
		nz = 256
		niter = 20
	elseif class == "C"
		nx = 512
		ny = 512
		nz = 512
		niter = 20
	elseif class == "D"
		nx = 2048
		ny = 1024
		nz = 1024
		niter = 25
	elseif class == "E"
		nx = 4096
		ny = 2048
		nz = 2048
		niter = 25
	else
		@printf("setparams: (Internal Error): invalid class type %s\n", clp)
		exit(1)
	end

	maxdim = nx
	
	if ny > maxdim
		maxdim = ny
	end

	if nz > maxdim
		maxdim = nz
	end
	
	@printf(fp, "const nx = %d\n", nx)
	@printf(fp, "const ny = %d\n", ny)
	@printf(fp, "const nz = %d\n", nz)
	@printf(fp, "const maxdim = %d\n", maxdim)
	@printf(fp, "const niter_default = %d\n", niter)
	@printf(fp, "const np_min = %d\n", np)
	@printf(fp, "const ntdivnp = ((nx*ny)/np_min)*nz\n")
	@printf(fp, "const ntotal_f = 1.0*nx*ny*nz\n")

end

function write_ep_info(fp, np, clp)
	m = 0
	if clp == "S" 
		m = 24
	elseif clp == "W"
		m = 26
	elseif clp == "A" 
		m = 28
	elseif clp == "B"
		m = 30
	elseif clp == "C"
		m = 32
	elseif clp == "D"
		m = 36
	elseif clp == "E"
		m = 40
	else
		@printf("setparams: Internal error: invalid class type %s\n", clp)
		exit(1)
	end

	@printf(fp, "const cclass = \"%s\"\n", clp)
	@printf(fp, "const m      = %d\n", m)
	@printf(fp, "const npm    = %d\n", np)
end

# KCH TODO
function write_cg_info(fp, np, clp)
	println("ALERT: write_cg_info Not implemented")
end

# KCH TODO
function write_bt_info(fp, np, clp, stp)
	println("ALERT: write_bt_info Not implemented")
end

function write_misc_info(fp)

	@printf(fp, "const compiletime = \"%s\"\n", Dates.format(now(), "d u yyyy"))
	
	rand = "(none)"

	try 
		open(DEFFILE, "r") do ff
			for line in eachline(ff)
				if ismatch(r"\s*RAND\s*=\s*(\w+)", line)
					m = match(r"\s*RAND\s*=\s*(\w+)", line)
					rand = m.captures[1]
					break
				end
			end
		end
	catch 
		@printf("\nsetparams: File %s doesn't exist. To build the NAS benchmarks\n
				             you need to create it according to the instructions\n
				             in the README in the main directory and comments in\n
				             the file config/make.def.template\n", DEFFILE)
	end

	@printf(fp, "const npbversion = \"%s\"\n", VERSION)
	@printf(fp, "const rand = \"%s\"\n", rand)

	if rand == "randi8"
		@printf(fp, "include(\"../common/randi8.jl\")\n")
	elseif rand == "randdp"
		@printf(fp, "include(\"../common/randdp.jl\")\n")
	else
		@printf("setparams: (Internal error): Unknown randlc implementation %s\n", rand)
		exit(1)
	end

	@printf(fp, "import NPBRand\n")

end


# write_info(): Write new information to config file. 
#               First line is in a special format so we can read
#               it in again. Then comes a warning. The rest is all
#               specific to a particular benchmark. 
#
function write_info(tp, np, clp, stp)

	BT_TYPES = ["NONE", "FULL", "SIMPLE", "EPIO", "FORTRAN"]

	fp = 0

	try 
		fp = open(FILENAME, "w")
	catch err
		println("error opening file")
	end

	if tp == BT
		if (stp == -1 || stp == 0)
			@printf(fp, "# NPROCS = %d CLASS = %s\n", np, clp)
		else
			@printf(fp, "# NPROCS = %d CLASS = %s SUBTYPE = %s\n", np, clp, stp)
		end
	elseif tp == SP || tp == LU || tp == MG ||
		   tp == FT || tp == EP || tp == CG 
		@printf(fp, "# NPROCS = %d CLASS = %s\n", np, clp)
	elseif tp == IS || tp == DT
		@printf(fp, "const CLASS = \"%s\"\n", clp)
		@printf(fp, "const NUM_PROCS = %d\n", np)
	else
		@printf("setparams: (Internal error): Unknown benchmark type %d\n", tp)
		exit(1)	
	end

	@printf(fp, "#\n#\n")
	@printf(fp, "#  This file is generated automatically by the setparams utility\n")
	@printf(fp, "#  It sets the number of processors and the class of the NPB\n")
	@printf(fp, "#  in this directory. Do not modify it by hand.\n")
	@printf(fp, "#\n")

	if tp == SP
		write_sp_info(fp, np, clp)
	elseif tp == LU
		write_lu_info(fp, np, clp)
	elseif tp == MG
		write_mg_info(fp, np, clp)
	elseif tp == IS
		write_is_info(fp, np, clp)
	elseif tp == DT
		write_dt_info(fp, np, clp)
	elseif tp == FT
		write_ft_info(fp, np, clp)
	elseif tp == EP
		write_ep_info(fp, np, clp)
	elseif tp == CG
		write_cg_info(fp, np, clp)
	elseif tp == BT
		write_bt_info(fp, np, clp, stp)
	else
		@printf("setparams: (Internal error): Unknown benchmark type %d\n", tp)
	end

	write_misc_info(fp)
	
	close(fp)

end


###########################################
############ main program #################
###########################################

(np, clp, tp, stp) = get_info(ARGS)

if clp != "U"
	if VERBOSE == true
		@printf("setparams: For benchmark %s: number of processors = %d class = %s\n", tp, np, clp)
	end
	check_info(tp, np, clp)
end

(nprocs_old, class_old, old_subtype) = read_info(tp)

if clp != "U"
	if class_old != "X"
		if VERBOSE == true
			@printf("setparams:     old settings: number of processors = %d class = %s\n", nprocs_old, class_old)
		end
	end
else 
	println("found class $clp")
	ms = 
"""
  *********************************************************************
  * You must specify NPROCS and CLASS to build this benchmark         *
  * For example, to build a class A benchmark for 4 processors, type  *
  *       make {benchmark-name} NPROCS=4 CLASS=A                      *
  *********************************************************************
"""
	@printf("setparams:\n%s", ms)

	if class_old != "X"
		if VERBOSE == 1
			@printf("setparams: Previous settings were CLASS=%s NPROCS=%d\n", class_old, nprocs_old)
		end
	end
	exit(1)
end

if np != nprocs_old || clp != class_old || stp != old_subtype
	if VERBOSE == true
		@printf("setparams: Writing %s\n", FILENAME)
	end
	write_info(tp, np, clp, stp)
else
	if VERBOSE == true
		@printf("setparams: Settings unchanged. %s unmodified\n", FILENAME)
	end
end

