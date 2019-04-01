	/******************************************************************************************************
	 * 																									  *
	 * 						Functions developed to perform 2D FFT in MPI Space							  *
	 * 																									  *
	 ******************************************************************************************************
																				 Author: Dr. Mirco Meazzo */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

typedef double FFT_SCALAR;

/*============================================= Functions Def =============================================*/
void z_aliasing(int nx, int ny, int nz, int nzd, FFT_SCALAR *U, FFT_SCALAR *U_read){
		int nz_left = 1+ (nz-1)/2; 	int reader=0;
	
		for( int stride_x = 0; stride_x < 2*nzd*ny*nx; stride_x = stride_x + 2*nzd*ny) {
			for( int stride_y = 0; stride_y < 2*nzd*ny; stride_y = stride_y + 2*nzd) {
				for (int k= (nzd-nz_left+1)*2; k < nzd*2; k++){
					U[stride_x + stride_y+k] = U_read[reader];
					reader++;
					//printf("U[%d] = %g\n", stride_x + stride_y+k, U[stride_x + stride_y+k]);
				}
				for (int k= (nz_left)*2; k < (nzd-nz_left+1)*2; k++){
					U[stride_x + stride_y+k]=0;
				}
				for (int k= 0; k < (nz_left)*2; k++){
					U[stride_x + stride_y+k] = U_read[reader];
					reader++;
					//printf("U[%d] = %g\n", stride_x + stride_y+k, U[stride_x + stride_y+k]);
				}
		}
	}
}

void x_aliasing(int nx, int ny, int nzd, int nxd, FFT_SCALAR *U, FFT_SCALAR *U_read){
	int reader = nx*ny*nzd*2-1;		//printf("reader %d", reader);
	for( int stride_z = 2*nxd*ny*nzd-1; stride_z > 0 ; stride_z = stride_z - 2*nxd*ny) {			//Backward x evitare sovrascritture
		//printf("stride_z %d\n", stride_z)
		for( int stride_y = 0; stride_y < 2*nxd*ny ; stride_y = stride_y + 2*nxd) {
			//printf("\tstride_y %d\n", stride_y);
			for( int i = 0; i < 2*(nxd-nx); i++) {
				U[stride_z - stride_y - i] = 0;		
				//printf("U0[%d] = %g\n", stride_z - stride_y - i, U[stride_z - stride_y - i]);	
			}
			for( int i = 2*(nxd-nx); i < 2*nxd; i++) {
				U[stride_z - stride_y - i] = U_read[reader];
				//printf("reader %d\n", reader);
				reader--;
				//printf("U[%d] = %g\n", stride_z - stride_y - i, U[stride_z - stride_y - i]);			
			}
		}
	}
}

void x_dealiasing(int scounts, int nx, int nxd, FFT_SCALAR *u, int des_rank, int rank) {
	/* scounts = initial number of ny*nz modes */
	int stride_x, placeholder=0;
	for (int mode =0; mode < scounts; mode++) {
		stride_x = mode*nxd*2;
		for (int i = 0; i < 2*nx; i++) {
			u[placeholder] = u[stride_x+i];
			placeholder++;
			//if (rank==des_rank) printf("u[%d]= %g\n", placeholder, u[placeholder]);
		}
	}
}

void z_dealiasing(int nx, int ny, int nz, int nzd, FFT_SCALAR *U) {
	int writer=0;	int reader=0;	int nz_left = 1+(nz-1)/2;
	FFT_SCALAR *temp = (FFT_SCALAR *) malloc(2*nz_left*sizeof(FFT_SCALAR));
	if (temp == NULL) {
		perror(".:Error while allocating temporary vector in z_dealiasing routine:.\n");
		abort();
	}
	for (int stride_x = 0; stride_x < 2*nx*ny*nzd; stride_x = stride_x + 2*ny*nzd ) {
		for (int stride_y = 0; stride_y < 2*ny*nzd; stride_y = stride_y + 2*nzd) {
			reader=0;
			for (int k= 0; k < 2*nz_left; k++) {
				temp[reader] = U[stride_x + stride_y + k];
				reader++;
			}	
			for (int k=2*(nzd-(nz_left-1)); k < 2*nzd; k++ ) {
				U[writer] = U[stride_x + stride_y + k];
				writer++;	
			}	
			reader=0;
			for (int k = 2*(nz_left-1); k < 2*nz; k++) {
				U[writer] = temp[reader];
				writer++;	reader++;	
			}
		}		
	}
}

void cores_handler( int modes, int size, int *modes_per_proc) {
	int rank =0;
	int check=0;

	for (int i = 0; i < modes; i++) {
		modes_per_proc[rank] = modes_per_proc[rank]+1;
		rank = rank+1;
		if (rank == size ) rank = 0;
	}
	
	for (int i = 0; i < size; i++){
		//printf("%d modes on rank %d\n", modes_per_proc[i], i);
		check = check+modes_per_proc[i];
	}
	if ( (int)(check - modes) != 0 ) {
			printf("[ERROR] check - modes = %d!!\nUnable to scatter modes properly\nAbort... \n", check - modes);
	}
}

void read_data(int nx, int ny, int nz, FFT_SCALAR *U_read, char file_to_read[4]) {
	//On rank 0 read the dataset
	FILE *U_dat;	U_dat = fopen( file_to_read, "r");
	printf("Reading initial field from %s...\n",file_to_read);
	for ( int i = 0; i < (nx)*(ny)*(nz)*2; i++) {
		fscanf( U_dat, "%lf", &U_read[i]);
		//printf("I've read %lf\n", U_read[i]);
	}
	printf("\tReading completed\n");
}

/* APPLY GLOBAL AA Z PENCIL VERSION */
void print_z_pencil(int nz, int in_ilo, int in_ihi, int in_jlo,
		FFT_SCALAR *u, int rank, int scounts, int desidered_rank);
void apply_AA(int nx, int ny, int nz, int nxd, int nzd, FFT_SCALAR *U, FFT_SCALAR *U_read){
	int nz_left = 1+ (nz-1)/2; 	int reader=0;
	// z-AA
	for( int stride_y = 0; stride_y < 2*nzd*nxd*ny; stride_y = stride_y + 2*nzd*nxd) {
		for( int stride_x = 0; stride_x < 2*nzd*nx; stride_x = stride_x + 2*nzd) {
			for (int k= (nzd-nz_left+1)*2; k < nzd*2; k++){
				U[stride_x + stride_y+k] = U_read[reader];
				reader++;
			}
			for (int k= (nz_left-1)*2; k < (nzd-nz_left)*2; k++){
				U[stride_x + stride_y+k]=0;
			}
			for (int k= 0; k < (nz_left)*2; k++){
				U[stride_x + stride_y+k] = U_read[reader];
				reader++;
			}
		}
		for (int stride_x = 2*nzd*nx; stride_x < 2*nzd*nxd; stride_x = stride_x + 2*nzd*nxd) {
			for (int k= 0; k < nzd*2; k++){		
				U[stride_x + stride_y+k]=0;
			}
		}
	}
	//print_z_pencil(nzd, 0, nxd-1, 0, U, 0, 2*nxd*ny*nzd, 0);
}

/*	APPLY GLOBAL AA X-PENCIL VERSION 
void apply_AA(int nx, int ny, int nz, int nxd, int nzd, FFT_SCALAR *U, FFT_SCALAR *U_read) {
	int nz_left = 1+ (nz-1)/2, reader= 2*nx*ny*(nz_left-1) ;
	  	  //Fill the array with read values and zeros for AA
	int i, stride_y, stride_z, last_index;
	for ( stride_z = 0; stride_z < nz_left*ny*nxd*2; stride_z = stride_z + ny*nxd*2) {
		//printf("\n\nstride z %d\n", stride_z );
		for ( stride_y = 0; stride_y < ny*nxd*2; stride_y = stride_y + nxd*2) {
			//printf("\nstride y %d\n", stride_y );
			for ( i = 0; i < (nx)*2; i++) {
				U[stride_z + stride_y+i] = U_read[reader];
				//printf("U[%d] =  %g\n", (stride_z + stride_y+i), U[stride_z + stride_y+i]);
				reader++;
			}
			for ( i = (nx)*2; i < nxd*2; i++) {
				U[stride_z + stride_y+i] = 0;
				// printf("U[%d] =  %g\n", (stride_z + stride_y+i), U[stride_z + stride_y+i]);
			}
		}
		last_index = stride_z + stride_y;
		//printf("last %d\n", (nzd - nz_left+1)*nxd*ny*2);
	}
	//Fill with zeros from nz to nzd
	for ( int i = last_index; i < (nzd - nz_left+1)*nxd*ny*2; i++) {
		U[i] = 0;
	}
	reader= 0;
	for ( stride_z = (nzd - nz_left+1)*nxd*ny*2; stride_z < nzd*ny*nxd*2; stride_z = stride_z + ny*nxd*2) {
		//printf("\n\nstride z %d\n", stride_z );
		for ( stride_y = 0 ; stride_y < ny*nxd*2; stride_y = stride_y + nxd*2) {
			//printf("\nstride y %d\n", stride_y );
			for ( i = 0; i < (nx)*2; i++) {
				U[stride_z + stride_y+i] = U_read[reader];
				//printf("U[%d] =  %g\n", (stride_z + stride_y+i), U[stride_z + stride_y+i]);
				reader++;
			}
			for ( i = (nx)*2; i < nxd*2; i++) {
				U[stride_z + stride_y+i] = 0;
				//printf("U[%d] =  %g\n", (stride_z + stride_y+i), U[stride_z + stride_y+i]);
			}
		}
		last_index = stride_z + stride_y;
		//printf("last %d\n", last_index);
	}
	 for (int i =0; i < nxd*nzd*ny*2; i++) {
	  		  printf("u[%d] = %g\n", i, U[i]);
	  	  } 
}*/

void print_x_pencil(int nx, int in_jlo, int in_jhi, int in_klo,
		FFT_SCALAR *u, int rank, int scounts, int desidered_rank) {
if (rank == desidered_rank) {
	  int stride_nz = in_klo;
	  int stride_ny = in_jlo;

	  for (int i = 0; i < scounts; i++) {
   	  if ( i % (nx*2) == 0) {
   		  printf("========(ny= %d, nz= %d)=======\n", stride_ny , stride_nz);
   		  if ( (stride_ny ) == in_jhi) {
   			  stride_ny = in_jlo;
   			  stride_nz ++;
   		  }
   		  else stride_ny ++;
   	  }
   	  printf("u[%d]= %.10f\n", (i), u[i]);
	  }
}
}

void print_sx_pencil(int nx, int in_jlo, int in_jhi, int in_klo,
		FFT_SCALAR *u, int rank, int scounts, int desidered_rank) {
if (rank == desidered_rank) {
	  int stride_nz = in_klo;
	  int stride_ny = in_jlo;

	  for (int i = 0; i < scounts; i++) {
   	  if ( i % (nx) == 0) {
   		  printf("========(ny= %d, nz= %d)=======\n", stride_ny , stride_nz);
   		  if ( (stride_ny ) == in_jhi) {
   			  stride_ny = in_jlo;
   			  stride_nz ++;
   		  }
   		  else stride_ny ++;
   	  }
   	  printf("u[%d]= %.10f\n", (i), u[i]);
	  }
}
}

void print_y_pencil(int ny, int y_ilo, int y_khi, int y_klo, 
		FFT_SCALAR *u, int rank, int scounts, int desidered_rank) {
if (rank == desidered_rank) {
	int stride_nx = y_ilo;
	int stride_nz = y_klo;

	for (int i = 0; i < scounts; i++) {
   		if ( i % (ny*2) == 0) {
   			printf("========(nx= %d, nz= %d)=======\n", stride_nx , stride_nz);
   			if ( (stride_nz ) == y_khi) {
   				stride_nz = y_klo;
   				stride_nx++;
   			}
   			else stride_nz++;
   			}
   		printf("u[%d]= %.10f\n", (i), u[i]);
    	}
 	}
}

void print_z_pencil(int nz, int in_ilo, int in_ihi, int in_jlo,
		FFT_SCALAR *u, int rank, int scounts, int desidered_rank) {
if (rank == desidered_rank) {
	int stride_nx = in_ilo;
	int stride_ny = in_jlo;

	for (int i = 0; i < scounts; i++) {
		if ( i % (nz*2) == 0) {
			printf("========(nx= %d, ny= %d)=======\n", stride_nx , stride_ny);
			if ( (stride_nx) == in_ihi) {
				stride_nx = in_ilo;
				stride_ny++;
			}
			else stride_nx ++;
		}
		printf("u[%d]= %.10f\n", (i), u[i]);
	}
}
}

void remove_undesidered_xmodes(FFT_SCALAR *u, int nx_old, int nx_new, int jhi, int khi) {
	int stride_nz, stride_ny, writer=0;
	for(int k = 0; k< khi; k++) {
		stride_nz = k*(jhi)*nx_old;
		for(int j = 0; j< jhi; j++) {
			stride_ny = j*nx_old;
			for(int i = 0; i < nx_new; i++) {
				u[writer] = u[i+stride_ny+stride_nz];
				//if (rank==1) printf("u[%d]=%f\n", writer, u[writer]);
				writer++;
			}
		}
	}
}


void Alltoall(int rank, int size, int in_jlo, int in_jhi, int in_ilo,
					 int in_ihi, int nz, int nx, FFT_SCALAR *arr, FFT_SCALAR *arr_recv, int flag){
	/* Flag = 1 	=> 	Scatterw
	 * Flag = -1 	=>	Gatherw */

	int *contiguous_x = (int *) malloc(sizeof(int)*size);
	int *contiguous_y = (int *) malloc(sizeof(int)*size);
	int *sendcounts = (int *) malloc(sizeof(int)*size);
	int *senddispls = (int *) malloc(sizeof(int)*size);
	int *recvdispls = (int *) malloc(sizeof(int)*size);
	int *recvcounts = (int *) malloc(sizeof(int)*size);
	if (( contiguous_x||contiguous_y||senddispls||sendcounts||recvdispls||recvcounts ) == NULL) {
		perror(".:Error while allocating memory for Alltoallw parameters:.\n");
		abort();
	}
	MPI_Datatype recvtype[size];
	contiguous_y[rank] = (in_jhi-in_jlo+1);
	contiguous_x[rank] = (in_ihi-in_ilo+1);
	for (int i = 0; i < size; i++){
		sendcounts[i] = 0;	recvdispls[i] = 0;		recvcounts[i] = 0;		recvtype[i] = MPI_DOUBLE;
	}
	// Broadcaster is the only one that send something
	if (rank == 0) {
		for (int i  = 0; i < size; i++){
			sendcounts[i] = 1;
		}
	}
	senddispls[rank] = (2*nz*in_ilo + 2*nz*nx*in_jlo )*sizeof(double);
	recvcounts[0] = 2*nz*(in_jhi-in_jlo+1)*(in_ihi-in_ilo+1);
	//printf("RECV COUNTS %d\n", recvcounts[0]);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Allgather(&contiguous_y[rank],1,MPI_INT,contiguous_y,1,MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(&contiguous_x[rank],1,MPI_INT,contiguous_x,1,MPI_INT, MPI_COMM_WORLD);
	MPI_Allgather(&senddispls[rank],1,MPI_INT,senddispls,1,MPI_INT, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Datatype vector[size], contiguous[size];
	int bytes_stride = sizeof(double)*2*nz*nx;

	for (int i = 0; i < size; i++) {
		MPI_Type_contiguous(2*nz*contiguous_x[i], MPI_DOUBLE, &contiguous[i]);
		MPI_Type_create_hvector(contiguous_y[i], 1, bytes_stride, contiguous[i], &vector[i]);
		MPI_Type_commit(&vector[i]);
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (flag == 1) {
		MPI_Alltoallw(&arr[0], sendcounts, senddispls, vector, &arr_recv[0], recvcounts, recvdispls, recvtype, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	}
	else if (flag == -1) {
		MPI_Alltoallw(&arr_recv[0], recvcounts, recvdispls, recvtype, &arr[0], sendcounts, senddispls, vector, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	} 
	else {
		perror(".:Invalid FLAG for Alltoall call:.\n\n");
		abort();
	}
	/*/Check function
	if (rank == 3){
		  for(int i = 0; i < recvcounts[0]; i++){
			  printf("arr_recv[%d]= %f\n", i, arr_recv[i]);
		  }
	  } */
	MPI_Type_free(vector);
}

