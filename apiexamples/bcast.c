#include<stdio.h>
#include<mpi.h>
#include<time.h>
#include<sys/time.h>

int main()
{
	//clock_t start,end;
	struct timespec start,end;
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
	clock_gettime(CLOCK_REALTIME,&start);
	MPI_Bcast(buffer,3,MPI_INT,1,MPI_COMM_WORLD);
	//end = clock();
	clock_gettime(CLOCK_REALTIME,&end);
	MPI_Finalize();
	if(rank == 0)
		//printf("cputime elapsed:	%f\n",((double)(end - start))/CLOCKS_PER_SEC);
		printf("cputime elapsed:	%ld\n",1000000000L * (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec));
	return 0;
}
