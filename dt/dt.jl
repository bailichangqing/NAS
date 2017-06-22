#=*************************************************************************
 *                                                                       *
 *        N  A  S     P A R A L L E L     B E N C H M A R K S  3.3       *
 *                                                                       *
 *                                  D T					 *
 *                                                                       *
 *************************************************************************
 *                                                                       *
 *   This benchmark is part of the NAS Parallel Benchmark 3.3 suite.     *
 *                                                                       *
 *   Permission to use, copy, distribute and modify this software        *
 *   for any purpose with or without fee is hereby granted.  We          *
 *   request, however, that all derived work reference the NAS           *
 *   Parallel Benchmarks 3.3. This software is provided "as is"          *
 *   without express or implied warranty.                                *
 *                                                                       *
 *   Information on NPB 3.3, including the technical report, the         *
 *   original specifications, source code, results and information       *
 *   on how to submit new results, is available at:                      *
 *                                                                       *
 *          http:  www.nas.nasa.gov/Software/NPB                         *
 *                                                                       *
 *   Send comments or suggestions to  npb@nas.nasa.gov                   *
 *   Send bug reports to              npb-bugs@nas.nasa.gov              *
 *                                                                       *
 *         NAS Parallel Benchmarks Group                                 *
 *         NASA Ames Research Center                                     *
 *         Mail Stop: T27A-1                                             *
 *         Moffett Field, CA   94035-1000                                *
 *                                                                       *
 *         E-mail:  npb@nas.nasa.gov                                     *
 *         Fax:     (650) 604-3957                                       *
 *                                                                       *
 *************************************************************************
 *                                                                       *
 *   Author: M. Frumkin							 *						 *
 *                                                                       *
 ************************************************************************=#
include("./npbparams.jl")
include("./DGraph.jl")
include("../common/timers.jl")
include("../common/print_results.jl")

using PrintResults.print_results
using NPBRand.randlc
import MPI

if isdefined(:CLASS) == false
	CLASS = "S"
	NUM_PROCS = 1
end

passed_verification = 0
timer_on = 0
timers_tot = 64

function verify(bmname,rnm2)
	verify_value = 0.0
	epsilon = 1.0E-8
	cls = CLASS
	verified = -1
	if cls != "U"
		if cls == "S"
			if searchindex(bmname,"BH") > 0
				verify_value = 30892725.0
			elseif searchindex(bmname,"WH") > 0
				verify_value = 67349758.0
			elseif searchindex(bmname,"SH") > 0
				verify_value = 58875767.0
			else
				@printf(STDERR,"No such benchmark as %s.\n",bmname)
			end
			verified = 0
		elseif cls == "W"
			if searchindex(bmname,"BH") > 0
				verify_value = 4102461.0
			elseif searchindex(bmname,"WH") > 0
				verify_value = 204280762.0
			elseif searchindex(bmname,"SH") > 0
				verify_value = 186944764.0
			else
				@printf(STDERR,"No such benchmark as %s.\n",bmname)
			end
			verified = 0
		elseif cls == "A"
			if searchindex(bmname,"BH") > 0
				verify_value = 17809491.0
			elseif searchindex(bmname,"WH") > 0
				verify_value = 1289925229.0
			elseif searchindex(bmname,"SH") > 0
				verify_value = 610856482.0
			else
				@printf(STDERR,"No such benchmark as %s.\n",bmname)
			end
			verified = 0
		elseif cls == "B"
			if searchindex(bmname,"BH") > 0
				verify_value = 4317114.0
			elseif searchindex(bmname,"WH") > 0
				verify_value = 7877279917.0
			elseif searchindex(bmname,"SH") > 0
				verify_value = 1836863082.0
			else
				@printf(STDERR,"No such benchmark as %s.\n",bmname)
			end
			verified = 0
		elseif cls == "C"
			if searchindex(bmname,"BH") > 0
				verify_value = 0.0
			elseif searchindex(bmname,"WH") > 0
				verify_value = 0.0
			elseif searchindex(bmname,"SH") > 0
				verify_value = 0.0
			else
				@printf(STDERR,"No such benchmark as %s.\n",bmname)
			end
			verified = -1
		elseif cls == "D"
			if searchindex(bmname,"BH") > 0
				verify_value = 0.0
			elseif searchindex(bmname,"WH") > 0
				verify_value = 0.0
			elseif searchindex(bmname,"SH") > 0
				verify_value = 0.0
			else
				@printf(STDERR,"No such benchmark as %s.\n",bmname)
			end
			verified = -1
		else
			@printf(STDERR,"No such class as %c.\n",cls)
		end
		@printf(STDERR," %s L2 Norm = %f\n",bmname,rnm2)
		if verified == -1
			@printf(STDERR," No verification was performed.\n")
		elseif rnm2 - verify_value < epsilon &&
					 rnm2 - verify_value > -epsilon	#abs here does not work on ALTIX
			verified = 1
			@printf(STDERR," Deviation = %f\n",(rnm2 - verify_value))
		else
			verified = 0
			@printf(STDERR," The correct verification value = %f\n",verify_value)
			@printf(STDERR," Got value = %f\n",rnm2)
		end
	else
		verified = -1
	end
	return verified
end

function ipowMod(a::Int32,n::Int64,md::Int32)
	seed = 1
	q = a
	r = 1
	if n < 0
		@printf(STDERR,"ipowMod: exponent must be nonnegative exp=%lld\n",n)
		n = -n	#temp fix
		#return 1
	end
	if md <= 0
		@printf(STDERR,"ipowMod: module must be positive mod=%d",md)
		return 1
	end
	if n == 0
		return 1
	end
	while n > 1
		n2 = div(n,2)
		if n2 * 2 == n
			seed = (q * q) % md
			q = seed
			n = n2
		else
			seed = (r * q) % md
			r = seed
			n = n - 1
		end
	end
	seed = (r * q) % md
	return seed
end

function buildSH(cls::String)
	#=
		Nodes of the graph must be topologically sorted
		to avoid MPI deadlock.
	=#
	numSources = NUM_SOURCES	#must be power of 2
	numOfLayers = 0
	tmpS = numSources >> 1
	firstLayerNode = 0
	mask = 0x0
	ndid = 0
	ndoff = 0
	i = 0
	j = 0
	nm = "DT_SH.$cls"
	dg = newDGraph(nm)

	while tmpS > 1
		numOfLayers += 1
		tmpS >>= 1
	end
	for i = 0:numSources - 1
		nm = "Source.$i"
		nd = newNode(nm)
		AttachNode(dg,nd)
	end
	for j = 0:numOfLayers - 1
		mask = 0x00000001 << j
		for i = 0:numSources - 1
			nm = "Comparator.$(i+j*firstLayerNode)"
			nd = newNode(nm)
			AttachNode(dg,nd)
			ndoff = i&(~mask)
			ndid = firstLayerNode + ndoff
			ar = newArc(dg.node[ndid + 1],nd)
      AttachArc(dg,ar)
      ndoff += mask
      ndid = firstLayerNode+ndoff
      ar = newArc(dg.node[ndid],nd)
      AttachArc(dg,ar)
		end
		firstLayerNode += numSources
	end
	mask = 0x00000001 << numOfLayers
	for i = 0:numSources - 1
		nm = "Sink.$i"
		nd = newNode(nm)
		AttachNode(dg,nd)
		ndoff = i&(~mask)
    ndid = firstLayerNode + ndoff
    ar = newArc(dg.node[ndid],nd)
    AttachArc(dg,ar)
    ndoff += mask
    ndid = firstLayerNode+ndoff
    ar=newArc(dg.node[ndid],nd)
    AttachArc(dg,ar)
	end
	return dg
end

function buildWH(cls::String)
#=
	Nodes of the graph must be topologically sorted
	to avoid MPI deadlock.
=#
	i = 0
	j = 0
	numSources = NUM_SOURCES
	maxInDeg = 4
	numLayerNodes = numSources
	firstLayerNode = 0
	totComparators = 0
	numPrevLayerNodes = numLayerNodes
	id = 0
	sid = 0
	nm = "DT_WH.$cls"
	dg = newDGraph(nm)
	for i = 0:numSources - 1
		nm = "Sink.$i"
		nd = newNode(nm)
		AttachNode(dg,nd)
	end
	totComparators = 0
	numPrevLayerNodes = numLayerNodes
	while numLayerNodes > maxInDeg
		numLayerNodes = div(numLayerNodes,maxInDeg)
		if numLayerNodes * maxInDeg < numPrevLayerNodes
			numLayerNodes += 1
		end
		for i = 0:numLayerNodes - 1
			nm = "Comparator.$totComparators"
			totComparators += 1
			nd = newNode(nm)
			id = AttachNode(dg,nd)
			for j = 0:maxInDeg - 1
				sid = i * maxInDeg + j
				if sid >= numPrevLayerNodes
					break
				end
				snd = dg.node[firstLayerNode + sid + 1]
				ar = newArc(dg.node[id + 1],snd)
				AttachArc(dg,ar)
			end
		end
		firstLayerNode += numPrevLayerNodes
		numPrevLayerNodes = numLayerNodes
	end
	source = newNode("source")
	AttachNode(dg,source)
	for i = 0:numPrevLayerNodes - 1
		nd = dg.node[firstLayerNode + i + 1]
		ar = newArc(source,nd)
		AttachArc(dg,ar)
	end
	for i = 0:div(dg.numNodes,2) - 1	#Topological sorting
		tmp = dg.node[i + 1]
		dg.node[i + 1] = dg.node[dg.numNodes - i]
		dg.node[i + 1].id = i
		dg.node[dg.numNodes - i] = tmp
		dg.node[dg.numNodes - i].id = dg.numNodes - 1 - i
	end
	return dg
end

function buildBH(cls::String)
#=
	Nodes of the graph must be topologically sorted
	to avoid MPI deadlock.
=#
	i = 0
	j = 0
	numSources = NUM_SOURCES
	maxInDeg = 4
	numLayerNodes = numSources
	firstLayerNode = 0
	totComparators = 0
	numPrevLayerNodes = numLayerNodes
	id = 0
	sid = 0
	nm = "DT_BH.$cls"
	dg = newDGraph(nm)

	for i = 0:numSources - 1
		nm = "Source.$i"
		nd = newNode(nm)
		AttachNode(dg,nd)
	end
	while numLayerNodes > maxInDeg
		numLayerNodes = div(numLayerNodes,maxInDeg)
		if numLayerNodes * maxInDeg < numPrevLayerNodes
			numLayerNodes += 1
		end
		for i = 0:numLayerNodes - 1
			nm = "Comparator.$totComparators"
			totComparators += 1
			nd = newNode(nm)
			id = AttachNode(dg,nd)
			for j = 0:maxInDeg - 1
				sid = i * maxInDeg + j
				if sid >= numPrevLayerNodes
					break
				end
				snd = dg.node[firstLayerNode + sid + 1]
				ar = newArc(snd,dg.node[id + 1])
				AttachArc(dg,ar)
			end
		end
		firstLayerNode += numPrevLayerNodes
		numPrevLayerNodes = numLayerNodes
	end
	sink = newNode("Sink")
	AttachNode(dg,sink)
	for i = 0:numPrevLayerNodes - 1
		nd = dg.node[firstLayerNode + i + 1]
		ar = newArc(nd,sink)
		AttachArc(dg,ar)
	end
	return dg
end

type Arr
	len::Int64
	val
end

function newArr(len::Int64)
	arr = Arr(len,Array(Float64,len))
	return arr
end

function arrShow(a)
	#RISK
	if typeof(a) == Void || typeof(a) != Arr
		@printf(STDERR,"-- NULL array\n")
	else
		@printf(STDERR,"-- length=%d\n",a->len)
	end
end

function CheckVal(feat::Arr)
	csum = 0.0
	for i = 0:feat.len - 1
		csum += (feat.val[i + 1] ^ 2) / feat.len	#= The truncation does not work since
																															 result will be 0 for large len	=#
	end
	return csum
end

function GetFNumDPar()
	return NUM_SAMPLES,STD_DEVIATION
end

function GetFeatureNum(mbname::String,id::Int32)
	tran = 314159265.0
	A = 2 * id + 1
	denom,tran = randlc(tran,Float64(A))
	cval = 'S'
	mean = NUM_SAMPLES
	stdev = 128
	rtfs = 0
	len = 0
	mean,stdev = GetFNumDPar()
	rtfs = ipowMod(Int32(round(1/denom,RoundToZero)) * Int32(cval),Int64(2 * id + 1),Int32(2 * stdev))
	if rtfs < 0
		rtfs = -rtfs
	end
	len = mean - stdev + rtfs
	return len
end

function RandomFeatures(bmname::String,fdim::Int32,id::Int32)
	len = GetFeatureNum(bmname,id) * fdim
	feat = newArr(len)
	nxg = 2
	nyg = 2
	nzg = 2
	nfg = 5
	nx = 421
	ny = 419
	nz = 1427
	nf = 3527
	expon = (len * (id + 1)) % 3141592
	seedx = ipowMod(Int32(nxg),expon,Int32(nx))
	seedy = ipowMod(Int32(nyg),expon,Int32(ny))
	seedz = ipowMod(Int32(nzg),expon,Int32(nz))
	seedf = ipowMod(Int32(nfg),expon,Int32(nf))
	#=
	if MPI.Comm_rank(MPI.COMM_WORLD) == 0
		opt = open("/Users/conghao/garbage/jopt","w")
		@printf(opt,"seedx = %d\nseedy = %d\nseedz = %d\nseedf = %d\nexpon = %lld\n",seedx,seedy,seedz,seedf,expon)
		close(opt)
	end
	=#
	i = 0
	if timer_on == 1
		Timers.timer_clear(id + 2)
		Timers.timer_start(id + 2)
	end
	while i < len - 1
		seedx = (seedx * nxg) %nx
		seedy = (seedy * nyg) %ny
		seedz = (seedz * nzg) %nz
		seedf = (seedf * nfg) %nf
		feat.val[i + 1] = seedx
		feat.val[i + 2] = seedy
		feat.val[i + 3] = seedz
		feat.val[i + 4] = seedf
		i += fdim
	end
	if timer_on == 1
		Timers.timer_stop(id + 2)
		@printf(STDERR,"** RandomFeatures time in node %d = %f\n",id,Timers.timer_read(id+2))
	end
	return feat
end

function Resample(a::Arr,blen::Int64)
	i = 0
	j = 0
	jlo = 0
	jhi = 0
	avval = 0.0
	nval = Array(Float64,blen)
	tmp = newArr(10)
	for i = 0:blen - 1
		nval[i + 1] = 0.0
	end
	for i = 1:a.len - 2
		jlo = Int64(round(0.5*(2*i-1)*(div(blen,a.len)),RoundToZero))
		jhi = Int64(round(0.5*(2*i+1)*(div(blen,a.len)),RoundToZero))

		avval = div(a.val[i + 1],jhi - jlo + 1)
		for j = jlo:jhi
			nval[j + 1] += avval
		end
	end
	nval[1] = a.val[1]
	nval[blen] = a.val[a.len]
	a.val = nval
	a.len = blen
end

fielddim = 4
function WindowFilter(a::Arr,b::Arr,w::Int32)
	i = 0
	j = 0
	k = 0
	rms0 = 0.0
	rms1 = 0.0
	rmsm1 = 0.0
	weight = Float64(div(w + 1,w + 2))

	w += 1
	if timer_on == 1
		Timers.timer_clear(w + 1)
		Timers.timer_start(w + 1)
	end
	if a.len < b.len
		Resample(a,b.len)
	end
	if a.len > b.len
		Resample(b,a.len)
	end
	i = fielddim
	while i < a.len - fielddim
		rms0 = (a.val[i + 1] - b.val[i + 1]) ^ 2
				 + (a.val[i + 2] - b.val[i + 2]) ^ 2
				 + (a.val[i + 3] - b.val[i + 3]) ^ 2
				 + (a.val[i + 4] - b.val[i + 4]) ^ 2
		j = i + fielddim
		rms1 = (a.val[j + 1] - b.val[j + 1]) ^ 2
				 + (a.val[j + 2] - b.val[j + 2]) ^ 2
				 + (a.val[j + 3] - b.val[j + 3]) ^ 2
				 + (a.val[j + 4] - b.val[j + 4]) ^ 2
		j = i - fielddim
		rmsm1 = (a.val[j + 1] - b.val[j + 1]) ^ 2
				  + (a.val[j + 2] - b.val[j + 2]) ^ 2
				  + (a.val[j + 3] - b.val[j + 3]) ^ 2
				  + (a.val[j + 4] - b.val[j + 4]) ^ 2
		k = 0
		if rms1 < rms0
			k = 1
			rms0 = rms1
		end
		if rmsm1 < rms0
			k = -1
		end
		if k == 0
			j = i + fielddim
			a.val[i + 1] = weight * b.val[i + 1]
			a.val[i + 2] = weight * b.val[i + 2]
			a.val[i + 3] = weight * b.val[i + 3]
			a.val[i + 4] = weight * b.val[i + 4]
		elseif k == 1
			j = i + fielddim
			a.val[i + 1] = weight * b.val[j + 1]
			a.val[i + 2] = weight * b.val[j + 2]
			a.val[i + 3] = weight * b.val[j + 3]
			a.val[i + 4] = weight * b.val[j + 4]
		else	#if k == -1
			j = i - fielddim
			a.val[i + 1] = weight * b.val[j + 1]
			a.val[i + 2] = weight * b.val[j + 2]
			a.val[i + 3] = weight * b.val[j + 3]
			a.val[i + 4] = weight * b.val[j + 4]
		end
		i += fielddim
	end
	if timer_on == 1
		Timers.timer_stop(w + 1)
		@printf(STDERR,"** WindowFilter time in node %d = %f\n",(w-1),Timers.timer_read(w + 1))
	end
	return a
end

function SendResults(dg::DGraph,nd::DGNode,feat)
	i = 0
	tag = 0
	if typeof(feat) == Void || typeof(feat) != Arr
		return 0
	end
	for i = 0:nd.outDegree - 1
		ar = nd.outArc[i + 1]
		if ar.tail != nd
			continue
		end
		head = ar.head
		tag = ar.id
		if head.address != nd.address
			MPI.Send(feat.len,head.address,tag,MPI.COMM_WORLD)
			MPI.Send(feat.val,head.address,tag,MPI.COMM_WORLD)
		end
	end
	return 1
end

function CombineStreams(dg::DGraph,nd::DGNode)
	resfeat = newArr(NUM_SAMPLES * fielddim)
	i = 0
	len = 0
	tag = 0
	if nd.inDegree == 0
		return 0
	end
	for i = 0:nd.inDegree - 1
		ar = nd.inArc[i + 1]
		if ar.head != nd
			continue
		end
		tail = ar.tail
		if tail.address != nd.address
			len = 0
			tag = ar.id
			len,status = MPI.Recv(Int,tail.address,tag,MPI.COMM_WORLD)
			feat = newArr(len)
			MPI.Recv!(feat.val,tail.address,tag,MPI.COMM_WORLD)
			resfeat = WindowFilter(resfeat,feat,nd.id)
		else
			featp = tail.feat
			feat = newArr(featp.len)
			for ii = 1:length(featp.val)
				feat.val[ii] = featp.val[ii]
			end
			resfeat = WindowFilter(resfeat,feat,nd.id)
		end
	end
	for i = 0:resfeat.len - 1
		tmppp = round(resfeat.val[i + 1],RoundToZero)
		if tmppp > typemax(Int32)
			println("ERROR!!!!!!!!!!!!!!!!!!					try to convert $(tmppp) into Int32")
		end
		tmpp = Int32(tmppp)
		tmp = div(tmpp,nd.inDegree)
		resfeat.val[i + 1] = tmp
	end
	nd.feat = resfeat
	return nd.feat
end

function Reduce(a::Arr,w::Int32)
	retv = 0.0
	if timer_on == 1
		Timers.timer_clear(w + 1)
		Timers.timer_start(w + 1)
	end
	retv = Int32(round(w * CheckVal(a),RoundToZero))	#= The casting needed for node
                               											 and array dependent verifcation	=#
	if timer_on == 1
		Timers.timer_stop(w + 1)
		@printf(STDERR,"** Reduce time in node %d = %f\n",(w-1),Timers.timer_read(w + 1))
	end
	return retv
end

function ReduceStreams(dg::DGraph,nd::DGNode)
	csum = 0.0
	for i = 0:nd.inDegree - 1
		ar = nd.inArc[i + 1]
		if ar.head != nd
			continue
		end
		tail = ar.tail
		if tail.address != nd.address
			tag = ar.id
			len,status = MPI.Recv(Int,tail.address,tag,MPI.COMM_WORLD)
			feat = newArr(len)
			MPI.Recv!(feat.val,tail.address,tag,MPI.COMM_WORLD)
			#=
			if nd.address == 4
				opt = open("/Users/conghao/garbage/jfeat","a+")
				for ii = 0:feat.len - 1
					@printf(opt,"%f\n",feat.val[ii + 1])
				end
				close(opt)
			end
			=#
			csum += Reduce(feat,Int32(nd.id + 1))
		else
			csum += Reduce(tail.feat,Int32(nd.id + 1))
		end
	end
	if nd.inDegree > 0
		csum = div(Int64(round(csum,RoundToZero)),nd.inDegree)
	end
	retv = (nd.id + 1) * csum
	return Float64(retv)
end

function ProcessNodes(dg::DGraph,me::Int64)
	verified::Int32 = 0
	for i = 0:dg.numNodes - 1
		nd = dg.node[i + 1]
		if nd.address != me
			continue
		end
		if contains(nd.name,"Source")
			#=
			if me == 3
				tfile = open("opt","w+")
				@printf(tfile,"dg->name = %s,fielddim = %d,nd.id = %d\n",dg.name,Int32(fielddim),Int32(nd.id))
				close(tfile)
			end
			=#
			nd.feat = RandomFeatures(dg.name,Int32(fielddim),Int32(nd.id))
			SendResults(dg,nd,nd.feat)
			#=
			if me == 3
				tfile = open("/Users/conghao/garbage/jopt","w+")
				@printf(tfile,"me:%d\n",me)
				for iterator = 0:length(nd.feat.val) - 1
					@printf(tfile,"%f\n",nd.feat.val[iterator + 1])
				end
				close(tfile)
			end
			=#
		elseif contains(nd.name,"Sink")
			chksum = ReduceStreams(dg,nd)
			tag = dg.numArcs + nd.id	#	make these to avoid clash with arc tags
			MPI.Send(chksum,0,tag,MPI.COMM_WORLD)
		else
			feat = CombineStreams(dg,nd)
			SendResults(dg,nd,feat)
		end
	end
	if me == 0
		rchksum = 0.0
		chksum = 0.0
		for i = 0:dg.numNodes - 1
			nd = dg.node[i + 1]
			if contains(nd.name,"Sink") != true
				continue
			end
			tag = dg.numArcs + nd.id	#	make these to avoid clash with arc tags
			rchksum,status = MPI.Recv(Float64,nd.address,tag,MPI.COMM_WORLD)
			chksum += rchksum
		end
		verified = verify(dg.name,chksum)
	end
	return verified
end

function main()
	MPI.Init()
	my_rank = MPI.Comm_rank(MPI.COMM_WORLD)
	comm_size = MPI.Comm_size(MPI.COMM_WORLD)
	if length(ARGS) != 1 ||
		 (	ARGS[1] != "BH"
		 	&&ARGS[1] != "WH"
			&&ARGS[1] != "SH")
		if my_rank == 0
			@printf(STDERR,"** Usage: mpirun -np N ../bin/dt.S GraphName\n")
			@printf(STDERR,"** Where \n   - N is integer number of MPI processes\n")
			@printf(STDERR,"   - S is the class S, W, or A \n")
			@printf(STDERR,"   - GraphName is the communication graph name BH, WH, or SH.\n")
			@printf(STDERR,"   - the number of MPI processes N should not be be less than \n")
			@printf(STDERR,"     the number of nodes in the graph\n")
		end
		MPI.Finalize()
		exit(0)
	end
	if ARGS[1] == "BH"
		dg = buildBH(CLASS)
	elseif ARGS[1] == "WH"
		dg = buildWH(CLASS)
	elseif ARGS[1] == "SH"
		dg = buildSH(CLASS)
	end

	if timer_on == 1 && dg.numNodes + 1 > timers_tot
		timer_on == 0
		if my_rank == 0
			@printf(STDERR,"Not enough timers. Node timeing is off. \n")
		end
	end
	if dg.numNodes > comm_size
		if my_rank == 0
			@printf(STDERR,"**  The number of MPI processes should not be less than \n")
			@printf(STDERR,"**  the number of nodes in the graph\n")
			@printf(STDERR,"**  Number of MPI processes = %d\n",comm_size)
			@printf(STDERR,"**  Number nodes in the graph = %d\n",dg.numNodes)
		end
		MPI.Finalize()
		exit(0)
	end
	for i = 0:dg.numNodes - 1
		dg.node[i + 1].address = i
	end
	if my_rank == 0
		@printf( "\n\n NAS Parallel Benchmarks 3.3 -- DT Benchmark\n\n" )
		graphShow(dg,0)
		Timers.timer_clear(1)
		Timers.timer_start(1)
	end
	verified = ProcessNodes(dg,my_rank)

	featnum = NUM_SAMPLES * fielddim
	bytes_sent = featnum * dg.numArcs
	bytes_sent = div(bytes_sent,1048576)
	if my_rank == 0
		Timers.timer_stop(1)
		tot_time = Timers.timer_read(1)
		print_results(dg.name,
										CLASS,
										featnum,
										0,
										0,
										dg.numNodes,
										0,
										comm_size,
										tot_time,
										div(bytes_sent,tot_time),
										"bytes transmitted",
										verified > 0,
										npbversion,
										rand)
	end
	MPI.Finalize()
	return 1
end

main()
#@printf("a = %d\n",ipowMod(Int32(2),7176,Int32(421)))
