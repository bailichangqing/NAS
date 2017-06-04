module PrintResults


# KCH TODO: get rid of compilation options
function print_results(name,
			  cclass,
			  n1,
			  n2,
			  n3,
			  niter,
			  nprocs_compiled,
			  nprocs_total,
		      tt,
			  mops,
			  optype,
			  passed_verification,
			  npbversion,
			  rand)

    @printf("\n\n %s Benchmark Completed\n", name)

    @printf(" Class           =                        %s\n", cclass)

    if (n2 == 0) && (n3 == 0) 
		if name[1:2] == "EP"
			@printf(" Size            =          %15.0f\n", 2.0^n1)
		else
			@printf(" Size            = %12d\n", n1)
		end
	else
		@printf(" Size            = %4dx%4dx%4d\n", n1, n2, n3)
	end
	

    @printf(" Iterations      =             %12d\n", niter)
 
    @printf(" Time in seconds =             %12.2f\n", tt)

    @printf(" Total processes =             %12d\n", nprocs_total)

    if nprocs_compiled != 0
        @printf(" Compiled procs  =             %12d\n", nprocs_compiled)
	end

    @printf(" Mop/s total     =             %12.2f\n", mops)

    @printf(" Mop/s/process   =             %12.2f\n", mops/nprocs_total)

    @printf(" Operation type  = %24s\n", optype)

    if passed_verification == true
        @printf(" Verification    =               SUCCESSFUL\n")
    else
        @printf(" Verification    =             UNSUCCESSFUL\n")
	end

    @printf(" Version         =             %12s\n", npbversion)
	@printf(" RAND            =             %12s\n", rand)

    @printf("\n\n" )
#=
    @printf(" Please send feedbacks and/or the results of this run to:\n\n" )
    @printf(" NPB Development Team\n" )
    @printf(" npb@nas.nasa.gov\n\n\n" )
=#
end



end # !module


				
