#include<stdio.h>
#include<mpi.h>
#include<time.h>

int main()
{
	//clock_t start,end;
	time_t start,end;
	MPI_Init(NULL,NULL);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int size;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	
	int buffer[3];
	
	if(rank == 0)
	{
		buffer[0] = 1;
		buffer[1] = 3;
		buffer[2] = 5;
	}
	else if(rank == 1)
	{
		buffer[0] = 2;
		buffer[1] = 4;
		buffer[2] = 6;
	}
	//start = clock();
	start = time(NULL);
	MPI_Bcast(buffer,3,MPI_INT,1,MPI_COMM_WORLD);
	//end = clock();
	end = time(NULL);
	MPI_Finalize();
	if(rank == 0)
		//printf("cputime elapsed:	%f\n",((double)(end - start))/CLOCKS_PER_SEC);
		printf("cputime elapsed:	%ld\n",(end - start));
	return 0;
}
