#include <iostream>
#include <math.h>
#include <mpi.h>

using namespace std;

#define MASTER	0		//main thread
#define NONE	-1		//for no neighbour
#define BEGIN	1
#define STEPS	100		//total number of steps
#define CRITERIA	0.004		//convergence criteria

int update(float *chunk, int chunkRow, int chunkCol, float *aboveRow, float *belowRow, float *leftCol, float *rightCol, int above, 
int below, int left, int right, int taskid, int print);

int main(int argc, char *argv[]) {
	int taskid, numProcs, dest, source, chunkCol, chunkRow, left, right, above, below, msgtype, rc;
	MPI_Status status;

	int convergence = 0;

	MPI_Init(&argc,&argv);
	MPI_Comm_size(MPI_COMM_WORLD,&numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD,&taskid);

	int M = 1024; //rows
	int N = 2048; //columns

	//determining the partitioning of chunks along the row and column
	if(numProcs==1){
		chunkRow = M;
		chunkCol = N;
	}
	else if(numProcs==2){
		chunkRow = M;
		chunkCol = N/2;
	}
	else if(numProcs==4 || numProcs==16){
		int sqroot = sqrt(numProcs);
		chunkRow = M/sqroot;
		chunkCol = N/sqroot;
	}
	else if(numProcs==8){
		int div = numProcs/2;
		chunkRow = M/2;
		chunkCol = N/div;
	}
	else if(numProcs==32){
		int div = numProcs/4;
		chunkRow = M/4;
		chunkCol = N/div;
	}
	else{
		cout<<"ERROR: the number of tasks must be 1, 2, 4, 8, 16 or 32"<<endl;
        cout<<"Quitting..."<<endl;
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
	}

	//one single chunk for each process
	float* myChunk = new float[chunkRow * chunkCol];

	//for receiving ghost rows and cols
	float aboveRow[chunkCol], belowRow[chunkCol], leftCol[chunkRow], rightCol[chunkRow];

	//main process code
	if(taskid==MASTER){
		float ** V = new float*[M];
    	for(int i = 0; i < M; i++){
        	V[i] = new float[N];
		}
		
		//initializing the 2D array
		for(int i=0; i<M; i++){
			V[i][0] = 1.0;
			V[i][N-1] = 1.0;
		}
		for(int i=0; i<N; i++){
			V[0][i] = 1.0;
			V[M-1][i] = 1.0;
		}
		for(int i=1; i<M-1; i++){
			for(int j=1; j<N-1; j++){
				V[i][j] = 0.0;
			}
		}
		
		//array of chunks/sub-regions
		float* chunks = new float[numProcs * chunkRow * chunkCol];
		
		//making chunks/sub-regions
		int i=0, j=0, numChunkRow=0, flag = 0;
		int colOffset=0, rowOffset=0;
		for(int chunk=0; chunk<numProcs; chunk++){
			for(int r=0; r<chunkRow; r++){
				j=chunkCol*colOffset;
				for(int c=0; c<chunkCol; c++){
					*(chunks + chunk * chunkRow * chunkCol + r * chunkCol + c) = V[i][j];
					j++;
				}
				i++;
			}
			if(j==N){
				rowOffset++;
				colOffset=0;
				j=0;
				//finding num of chunks in one row
				if(!flag){
					numChunkRow = chunk+1;
					flag=1;
				}
			}
			else{
				colOffset++;
			}
			i=rowOffset*chunkRow;
		}

		//deleting the main array after making chunks to free memory
		for(int i=0; i<M; i++){
      		delete [] V[i];  
		} 
      	delete [] V;

		//figuring out neighbors for each chunk
		for(int i=1; i<numProcs && numProcs>1; i++){
			//top row
			if(i>0 && i<numChunkRow){
				//top right
				if(i==(numChunkRow-1)){
					above=right = NONE;
					left = i-1;
					if(numProcs==2){
						below=-1;
					}
					else{
						below = i+numChunkRow;
					}
				}
				//top row with no above
				else{
					above=NONE;
					right=i+1;
					left=i-1;
					below=i+numChunkRow;
				}
			}
			//bottom row
			else if(i>=(numProcs-numChunkRow)){
				//bottom left
				if(i==(numProcs-numChunkRow)){
					left=below=NONE;
					right=i+1;
					above=i-numChunkRow;
				}
				//bottom right
				else if(i==(numProcs-1)){
					right=below=NONE;
					left=i-1;
					above=i-numChunkRow;
				}
				//bottom row with no below
				else{
					below=NONE;
					left=i-1;
					right=i+1;
					above=i-numChunkRow;
				}
			}
			//left column with no left neighbour
			else if(i%numChunkRow==0){
				left=NONE;
				right=i+1;
				above=i-numChunkRow;
				below=i+numChunkRow;
			}
			//right column with no right neighbor
			else if((i+1)%numChunkRow==0){
				right=NONE;
				left=i-1;
				above=i-numChunkRow;
				below=i+numChunkRow;
				//cout<<"right column"<<endl;
			}
			//all in between
			else{
				above=i-numChunkRow;
				right=i+1;
				left=i-1;
				below=i+numChunkRow;
			}

			dest = i;
			
			MPI_Send(&above, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&below, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&left, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(&right, 1, MPI_INT, dest, BEGIN, MPI_COMM_WORLD);
			MPI_Send(chunks+i*chunkRow*chunkCol, chunkRow*chunkCol, MPI_FLOAT, dest, BEGIN, MPI_COMM_WORLD);
		}

		//figuring out neighbours of master's chunk i.e. top left
		if(numProcs==1){
			above=below=left=right = NONE;
		}
		else if(numProcs==2){
			above=left=below=NONE;
			right=0+1;
		}
		else{
			above=left=NONE;
			right=0+1;
			below=0+numChunkRow;
		}

		//making chunk for master and one chunk for updating
		for(int j=0; j<chunkRow; j++){
			for(int k=0; k<chunkCol; k++){
				*((myChunk+j*chunkCol)+k) = *(chunks + 0 * chunkRow * chunkCol + j * chunkCol + k);
			}
		}

		int stepCount = 0;
		while(stepCount<STEPS){
			//for receiving rows cols
			if(numProcs > 1){
				//for sending rows cols
				float tempRow[chunkCol], tempCol[chunkRow];

				if(above!=NONE){
					for(int i=0; i<chunkCol; i++){
						tempRow[i] = *((myChunk+0*chunkCol)+i);
					}
					//non-blocking send
					if(numProcs==2 || numProcs==4){
						MPI_Request req;
						MPI_Isend(&tempRow, chunkCol, MPI_FLOAT, above, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Irecv(&aboveRow, chunkCol, MPI_FLOAT, above, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Wait(&req, &status);
					}
					//blocking send
					else{
						MPI_Send(&tempRow, chunkCol, MPI_FLOAT, above, BEGIN, MPI_COMM_WORLD);
						MPI_Recv(&aboveRow, chunkCol, MPI_FLOAT, above, BEGIN, MPI_COMM_WORLD, &status);
					}
				}
				if(below!=NONE){
					for(int i=0; i<chunkCol; i++){
						tempRow[i] = *((myChunk+(chunkRow-1)*chunkCol)+i);
					}
					if(numProcs==2 || numProcs==4){
						MPI_Request req;
						MPI_Isend(&tempRow, chunkCol, MPI_FLOAT, below, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Irecv(&belowRow, chunkCol, MPI_FLOAT, below, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Wait(&req, &status);
					}
					else{
						MPI_Send(&tempRow, chunkCol, MPI_FLOAT, below, BEGIN, MPI_COMM_WORLD);
						MPI_Recv(&belowRow, chunkCol, MPI_FLOAT, below, BEGIN, MPI_COMM_WORLD, &status);
					}
				}
				if(left!=NONE){
					for(int i=0; i<chunkRow; i++){
						tempCol[i] = *((myChunk+i*chunkRow)+0);
					}
					if(numProcs==2 || numProcs==4){
						MPI_Request req;
						MPI_Isend(&tempCol, chunkRow, MPI_FLOAT, left, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Irecv(&leftCol, chunkRow, MPI_FLOAT, left, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Wait(&req, &status);
					}
					else{
						MPI_Send(&tempCol, chunkRow, MPI_FLOAT, left, BEGIN, MPI_COMM_WORLD);
						MPI_Recv(&leftCol, chunkRow, MPI_FLOAT, left, BEGIN, MPI_COMM_WORLD, &status);
					}
				}
				if(right!=NONE){
					for(int i=0; i<chunkRow; i++){
						tempCol[i] = *((myChunk+i*chunkRow)+(chunkCol-1));
					}
					if(numProcs==2 || numProcs==4){
						MPI_Request req;
						MPI_Isend(&tempCol, chunkRow, MPI_FLOAT, right, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Irecv(&rightCol, chunkRow, MPI_FLOAT, right, BEGIN, MPI_COMM_WORLD, &req);
						MPI_Wait(&req, &status);
					}
					else{
						MPI_Send(&tempCol, chunkRow, MPI_FLOAT, right, BEGIN, MPI_COMM_WORLD);
						MPI_Recv(&rightCol, chunkRow, MPI_FLOAT, right, BEGIN, MPI_COMM_WORLD, &status);
					}
				}
			}

			if(convergence!=1){
				int print = 0;
				if((stepCount+1)%5 == 0){
					print = 1;
				}
				convergence = update(myChunk, chunkRow, chunkCol, (float *)aboveRow, (float *)belowRow, (float *)leftCol, (float *)rightCol, above, below, left, right, taskid, print);
			}
			stepCount++;
		}

	}

	//worker process code
	if(taskid!=MASTER){
		source = MASTER;
		msgtype = BEGIN;
		//for sending rows cols
		float tempRow[chunkCol], tempCol[chunkRow];

		MPI_Recv(&above, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&below, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&left, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(&right, 1, MPI_INT, source, msgtype, MPI_COMM_WORLD, &status);
		MPI_Recv(myChunk, chunkRow*chunkCol, MPI_FLOAT, source, msgtype, MPI_COMM_WORLD, &status);

		int stepCount = 0;
		while(stepCount<STEPS){
			if(above!=NONE){
				for(int i=0; i<chunkCol; i++){
					tempRow[i] = *((myChunk+0*chunkCol)+i);
				}
				if(numProcs==2 || numProcs==4){
					MPI_Request req1;
					MPI_Isend(&tempRow, chunkCol, MPI_FLOAT, above, BEGIN, MPI_COMM_WORLD, &req1);
					MPI_Irecv(&aboveRow, chunkCol, MPI_FLOAT, above, msgtype, MPI_COMM_WORLD, &req1);
					MPI_Wait(&req1, &status);
				}
				else{
					MPI_Send(&tempRow, chunkCol, MPI_FLOAT, above, BEGIN, MPI_COMM_WORLD);
					MPI_Recv(&aboveRow, chunkCol, MPI_FLOAT, above, msgtype, MPI_COMM_WORLD, &status);
				}
			}
			if(below!=NONE){
				for(int i=0; i<chunkCol; i++){
					tempRow[i] = *((myChunk+(chunkRow-1)*chunkCol)+i);
				}
				if(numProcs==2 || numProcs==4){
					MPI_Request req1;
					MPI_Isend(&tempRow, chunkCol, MPI_FLOAT, below, BEGIN, MPI_COMM_WORLD, &req1);
					MPI_Irecv(&belowRow, chunkCol, MPI_FLOAT, below, msgtype, MPI_COMM_WORLD, &req1);
					MPI_Wait(&req1, &status);
				}
				else{
					MPI_Send(&tempRow, chunkCol, MPI_FLOAT, below, BEGIN, MPI_COMM_WORLD);
					MPI_Recv(&belowRow, chunkCol, MPI_FLOAT, below, msgtype, MPI_COMM_WORLD, &status);
				}
			}
			if(left!=NONE){
				for(int i=0; i<chunkRow; i++){
					tempCol[i] = *((myChunk+i*chunkRow)+0);
				}
				if(numProcs==2 || numProcs==4){
					MPI_Request req1;
					MPI_Isend(&tempCol, chunkRow, MPI_FLOAT, left, BEGIN, MPI_COMM_WORLD, &req1);
					MPI_Irecv(&leftCol, chunkRow, MPI_FLOAT, left, msgtype, MPI_COMM_WORLD, &req1);
					MPI_Wait(&req1, &status);
				}
				else{
					MPI_Send(&tempCol, chunkRow, MPI_FLOAT, left, BEGIN, MPI_COMM_WORLD);
					MPI_Recv(&leftCol, chunkRow, MPI_FLOAT, left, msgtype, MPI_COMM_WORLD, &status);
				}
			}
			if(right!=NONE){
				for(int i=0; i<chunkRow; i++){
					tempCol[i] = *((myChunk+i*chunkRow)+(chunkCol-1));
				}
				if(numProcs==2 || numProcs==4){
					MPI_Request req1;
					MPI_Isend(&tempCol, chunkRow, MPI_FLOAT, right, BEGIN, MPI_COMM_WORLD, &req1);
					MPI_Irecv(&rightCol, chunkRow, MPI_FLOAT, right, msgtype, MPI_COMM_WORLD, &req1);
					MPI_Wait(&req1, &status);
				}
				else{
					MPI_Send(&tempCol, chunkRow, MPI_FLOAT, right, BEGIN, MPI_COMM_WORLD);
					MPI_Recv(&rightCol, chunkRow, MPI_FLOAT, right, msgtype, MPI_COMM_WORLD, &status);
				}
			}

			if(convergence!=1){
				int print = 0;
				if((stepCount+1)%5 == 0){
					print = 1;
				}
				convergence = update(myChunk, chunkRow, chunkCol, (float *)aboveRow, (float *)belowRow, (float *)leftCol, (float *)rightCol, above, below, left, right, taskid, print);
			}
			stepCount++;
		}

	}
	
	MPI_Finalize();
}

//updates and checks convergence
int update(float *chunk, int chunkRow, int chunkCol, float *aboveRow, float *belowRow, float *leftCol, float *rightCol, int above, 
int below, int left, int right, int taskid, int print){
	float* copyChunk = new float[chunkRow * chunkCol];
	for(int i=0; i<chunkRow; i++){
		for(int j=0; j<chunkCol; j++){
			*((copyChunk+i*chunkCol) + j) = *((chunk+i*chunkCol) + j);
		}
	}

	//setting offsets for chunk
	int i=0, j=0, rows=chunkRow, cols=chunkCol;
	if(above == -1){
		i=i+1;
	}
	if(left == -1){
		j=j+1;
	}
	if(right == -1){
		cols=cols-1;
	}
	if(below == -1){
		rows=rows-1;
	}
	int it=j;

	float sum=0.0;
	while(i<rows){
		j=it;
		while(j<cols){
			//above
			if(i==0 && above!=-1){
				sum=sum+(*(aboveRow+j));
			}
			else{
				sum=sum+(*((copyChunk+(i-1)*chunkCol) + j));
			}
			//below
			if(i==rows-1 && below!=-1){
				sum=sum+(*(belowRow+j));
			}
			else{
				sum=sum+(*((copyChunk+(i+1)*chunkCol) + j));
			}
			//left
			if(j==0 && left!=-1){
				sum=sum+(*(leftCol+i));
			}
			else{
				sum=sum+(*((copyChunk+i*chunkCol) + (j-1)));
			}
			//right
			if(j==cols-1 && right!=-1){
				sum=sum+(*(rightCol+i));
			}
			else{
				sum=sum+(*((copyChunk+i*chunkCol) + (j+1)));
			}
			*((chunk+i*chunkCol) + j) = sum/4;
			sum = 0.0;
			j++;
		}
		i++;
	}

	//printing the max difference after every 5 timesteps
	if(print){
		float max=0.0, diff;
		for(int i=0; i<chunkRow; i++){
			for(int j=0; j<chunkCol; j++){
				diff = abs((*((copyChunk+i*chunkCol) + j))-(*((chunk+i*chunkCol) + j)));
				if(diff > max){
					max = diff;
				}
			}
		}
		cout<<"Current max difference for task "<<taskid<<" is: "<<max<<endl;
	}

	//convergence check
	int converged=0;
	int notConverged=1;
	for(int i=0; i<chunkRow; i++){
		for(int j=0; j<chunkCol; j++){
			if(abs((*((copyChunk+i*chunkCol) + j))-(*((chunk+i*chunkCol) + j))) > CRITERIA){
				converged = 0;
				notConverged = 0;
				break;
			}
			else{
				converged = 1;
			}
		}
		if(notConverged==0){
			break;
		}
	}

	if(converged){
		cout<<"Task "<<taskid<<" has converged!"<<endl;
	}

	return converged;
}
