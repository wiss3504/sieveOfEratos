#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>
#include <stdbool.h>

#define TERMINATOR -1

void strikeMult(bool *arr, int start, int n, int prime){
	for(int i=start;i<=n;i+=prime){
		arr[i]=0;
	}
}

int main(int argc, char *argv[]){
	int rank,size,n;
	double s_time,e_time;
	
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	if(argc<2){
		if(rank == 0) printf("Usage: %s <n>\n", argv[0]);
		MPI_Finalize();
		return 1;
	}
	n=atoi(argv[1]);
	int sqrt_n = (int)sqrt(n);
	if(rank==0) s_time = MPI_Wtime();
	if(size==1){
		bool *prime = (bool *)malloc((n+1) * sizeof(bool));
		for(int i=2;i<=n;i++){
			prime[i] = 1;
		}
		for(int i=2;i<=sqrt_n;i++){
			if(prime[i]){
				strikeMult(prime,i*i,n,i);
			}
		}
		printf("Primes up to %d:\n", n);
        	for (int i = 2; i <= n; i++) {
            		if (prime[i]) {
                		printf("%d ", i);
            		}
        	}
        	printf("\n");
        	free(prime);
        	e_time = MPI_Wtime();
        	printf("Execution time (sequential): %f seconds\n", e_time - s_time);
	} else{
		bool *prime = (bool *)malloc((n+1) * sizeof(bool));
		for(int i=2;i<=n;i++){
			prime[i]=1;
		}
		if(rank==0){
			for(int i=2;i<=sqrt_n;i++){
				if(prime[i]){
					strikeMult(prime,i*i,n,i);
					MPI_Send(&i,1,MPI_INT,1,0,MPI_COMM_WORLD);
				}
			}
			int end_signal = TERMINATOR;
			MPI_Send(&end_signal,1,MPI_INT,1,0,MPI_COMM_WORLD);
		}else{
			int x;
			while(true){
				MPI_Recv(&x,				1,MPI_INT,rank-1,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
				if(x==TERMINATOR) break;
				strikeMult(prime,x*x,n,x);
				if(rank<size-1){
					MPI_Send(&x,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
				}
			}
			if(rank<size-1){
				int end_signal = TERMINATOR;
				MPI_Send(&end_signal,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
			}
		}
	
		if(rank==0){
			printf("Primes up to %d:\n",n);
			/*for(int i=2;i<=n;i++){
				if(prime[i]){
					printf("%d ",i);
				}
			}
			printf("\n");*/
			e_time = MPI_Wtime();
			printf("Execution time (parallel): %f seconds\n",e_time - s_time);
		}
		free(prime);
	}
	MPI_Finalize();
	return 0;
}