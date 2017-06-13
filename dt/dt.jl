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
	CLASS = 'S'
	NUM_PROCS = 1
end

passed_verification = 0
timer_on = 0
timers_tot = 64

function verify(bmname,rnm2)
	verify_value = 0.0
	epsilon = 
