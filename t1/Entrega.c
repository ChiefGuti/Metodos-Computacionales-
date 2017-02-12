#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

int L = 5, l = 2, d = 1, V0 = 100, m, N;
double h = 0.02;

int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);

int transformer(int i, int j)
{	
	
	return i*m + j;
}

double *init(int x0, int x1, int y0, int y1, double *array)
{	
	int a;
	for(a = x0; a <= x1; a++)
	{
		array[transformer(y0, a)] = V0/2;
		array[transformer(y1, a)] = -V0/2;
	}
	return array;
}
int main(void){

MPI_Init(NULL, NULL);

int world_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
int world_size;
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
int up, down,a,left, right, x0, x1, y0, y1, i=1, j=1, n=0;
double average;

a= round(m*m/world_size);
m = L/h;
N = 2*m*m;

x0 = m/2 - l/(h*2) - 1;
x1 = m/2 + l/(h*2) - 1;
y0 = m/2 - d/(h*2) - 1;
y1 = m/2 + d/(h*2) - 1;
	
	
double *V_new = malloc(m*m*sizeof(double));
double *V_sub = malloc(a*sizeof(double));
double *V = malloc(m*m*sizeof(double));
double *V1 = malloc(m*m*sizeof(double));

while(n<N)
{
	V= init(x0,x1,y0,y1,V);
	
	if (world_rank == 0)
	{
		
		for(i=0 ; i <m+a+1 ; i++)
		{
			V_sub[p]=V[i];
			
		}

		for(i=0 ; i<a+1 ; i++)
		{	
			if(!(i>m+1))
				{
				up=i-m;
				down=i+m;
				left=i-1;
				right=i+1;
				average = (V_sub[up]+ V_sub[down] + V_sub[left] + V_sub[right])/4;
				V_new[q]=average;
				}
			else
				{	
				V_new[i]= V_sub[i];	
				}					
		} 
	}

	if (!(world_rank==0 && world_rank == world_size))
	{
		int c=0;
		for(i= a*(world_rank)-m ; i <m+(1+world_rank)*a ; i++)
		{
		V_sub[c]=V[i];
		c+=1;
		}


		int l=0;
		for(i=m+1 ; i<a+1 ; i++)
		{
			up=i-m;
			down=i+m;
			left=i-1;
			right=i+1;
			average = (V_sub[up]+ V_sub[down] + V_sub[left] + V_sub[right])/4;
			V_new[l]=average;
			l+=1;
		} 
	
	}

	if(world_rank==world_size)
	{
		int h=0;
		for(i=a*world_rank-m ; i <m*m ; i++)
		{
			V_sub[h]=V[i];
			h+=1;
		}
	
		int k=0;
		for(i=m ; i<a+1 ; i++)
		{
			up=i-m;
			down=i+m;
			left=i-1;
			right=i+1;
			average = (V_sub[up]+ V_sub[down] + V_sub[left] + V_sub[right])/4;
			V_new[k]=average;
			k+=1;
		} 	
	}
	
	MPI_Allgather(V_new,a,MPI_DOUBLE,V1,a,MPI_DOUBLE,MPI_COMM_WORLD);
	for(i=0;i<m*m;i++)
	{
	V[i]=V1[i]
	}
	n=+1;	
}	

for(i=0 ; i <m*m ;i++)
{
	printf("%f\n", V[i]);
}
MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
return 0;
}
