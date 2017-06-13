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
include("../common/randdp.jl")
include("../common/print_results.jl")

using PrintResults.print_result
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
		if n2 * 2 == 2
			seed = (q * q) % md
			q = seed
			n = n2
		end
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
		nm = newNode(nm)
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
			for j = 0:maxInDeg
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
	for i = 0:numPrevLayerNodes
		nd = dg.node[firstLayerNode + i + 1]
		ar = newArc(nd,sink)
		AttachArc(dg,ar)
	end
	return dg
end

type Arr
	len::Int32
	val
end

function newArr(len::Int32)
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
	i = 0
	for i = 0:feat.len - 1
		csum += div(feat.val[i + 1] * feat.val[i + 1],feat.len)	#= The truncation does not work since
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
	denom,tran = randlc(tran,A)
	cval = 'S'
	mean = NUM_SAMPLES
	stdev = 128
	rtfs = 0
	len = 0
	mean,stdev = GetFNumDPar()
	rtfs = ipowMod(Int32(round(1/denom,RoundToZero)) * Int(cval),Int64(2 * id + 1),2 * stdev)
	if rtfs < 0
		rtfs = -rtfs
	end
	len = mean - stdev + rtfs
	return len
end
