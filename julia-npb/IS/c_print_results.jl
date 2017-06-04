function c_print_results(name,
                         class,
                         n1,
                         n2,
                         n3,
                         niter,
                         nprocs_compiled,
                         nprocs_total,
                         t,
                         mops,
                         optype,
                         passed_verification,
                         npbversion,
                         compiletime,
                         mpicc,
                         clink,
                         cmpi_lib,
                         cmpi_inc,
                         cflags,
                         clinkflags)
  #SMP part is not included in current version
  evalue = "1000"

  @printf( "\n\n %s Benchmark Completed\n", name );

  #@printf( " Class           =                        %c\n", class );

  if n3 == 0
      nn = n1;
      if  n2 != 0
        nn *= n2;
      end
      @printf( " Size            =             %12ld\n", nn );   # as in IS
  else
    @printf( " Size            =              %3dx %3dx %3d\n", n1,n2,n3 );
  end

  @printf( " Iterations      =             %12d\n", niter );

  @printf( " Time in seconds =             %12.2f\n", t );

  @printf( " Total processes =             %12d\n", nprocs_total );

  if nprocs_compiled != 0
    @printf( " Compiled procs  =             %12d\n", nprocs_compiled );
  end

  #@printf( " Mop/s total     =             %12.2f\n", mops );

  #@printf( " Mop/s/process   =             %12.2f\n", mops/((float) nprocs_total) );

  @printf( " Operation type  = %24s\n", optype);

  if passed_verification != 0
    @printf( " Verification    =               SUCCESSFUL\n" );
  else
    @printf( " Verification    =             UNSUCCESSFUL\n" );
  end
  #@printf( " Version         =             %12s\n", npbversion );

  #@printf( " Compile date    =             %12s\n", compiletime );

  #@printf( "\n Compile options:\n" );

  #@printf( "    MPICC        = %s\n", mpicc );

  #@printf( "    CLINK        = %s\n", clink );

  #@printf( "    CMPI_LIB     = %s\n", cmpi_lib );

  #@printf( "    CMPI_INC     = %s\n", cmpi_inc );

  #@printf( "    CFLAGS       = %s\n", cflags );

  #@printf( "    CLINKFLAGS   = %s\n", clinkflags );

  @printf( "\n\n" );
  @printf( " Please send feedbacks and/or the results of this run to:\n\n" );
  @printf( " NPB Development Team\n" );
  @printf( " npb@nas.nasa.gov\n\n\n" );

end
