#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <math.h>

/*definiio las variables dads en el enunciado*/

int L = 5, l = 2, d = 1, V0 = 100, m=512, N;
double h = 5/512;

int transformer(int i, int j);
double *init(int x0, int x1, int y0, int y1, double *array);

/*Esta funcion me transforma los indices de una arreglo 2x2 a uno unidimensional*/

int transformer(int i, int j)
{	
	
	return i*m + j;
}
/*defino la funcion que me impone las consdicion inicial de las varillas */

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

/*solicito  la cantidad d procesadores y el numero correpondinte a cada procesador*/

MPI_Init(NULL, NULL);

int world_rank;
MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
int world_size;
MPI_Comm_size(MPI_COMM_WORLD, &world_size);
int up, down,a,left, right, x0, x1, y0, y1, i=1, j=1, n=0;
double average;

/*definio la cantidad de elementos que tendra cada procesador y la cantidad de iteraciones*/

a= round(m*m/world_size);
N = 2*m*m;

/* caluclui las cordenadas de las varillas en un arreglo 2x2*/

x0 = m/2 - l/(h*2) - 1;
x1 = m/2 + l/(h*2) - 1;
y0 = m/2 - d/(h*2) - 1;
y1 = m/2 + d/(h*2) - 1;
	
/* cada procesador va a tenre en memoria estos tres vectores*/	

double *V_new = malloc(m*m*sizeof(double));
double *V_sub = malloc(a*sizeof(double));
double *V = malloc(m*m*sizeof(double));
double *V1 = malloc(m*m*sizeof(double));

while(n<N)
{
	/*cada procesador restables las condiciones de las varillas despues de cada iteracion*/
	V= init(x0,x1,y0,y1,V);
	
	/*Cada procesador va a buscar los elementos que le correspondes del arreglo V que todos tiene gusrdado*/
	if (world_rank == 0)
	{
		/*el procesador 0 tendra el comenzo del arreglo V mas una fila mas de elementos para poder calcular la relajacion
		* selecciono parte del arreglo V y lo guardo en V_sub*/
		for(i=0 ; i <m+a+1 ; i++)
		{
			V_sub[i]=V[i];
			
		}
		/*Calcullo la relajacion para los puntos*/
		for(i=0 ; i<a+1 ; i++)
		{	
			/* NO calculo la relajacion para la primera fila para que sepueda calcular el promedio
			*el valor del promedio lo guardo en V_new*/
			if(!(i<m+1))
				{
				up=i-m;
				down=i+m;
				left=i-1;
				right=i+1;
				average = (V_sub[up]+ V_sub[down] + V_sub[left] + V_sub[right])/4;
				V_new[i]=average;
				}
			else /* evito primera fila y gurado los valores originales de V*/
				{	
				V_new[i]= V_sub[i];	
				}					
		} 
	}
	/* Para los porcesadores que no tengan problemas de tenere el principio o el final del arreglo*/
	if (!(world_rank==0 && world_rank == world_size))
	{
		/* Seleciono la parte de V que le corresponde a cada procesador y la guardo enV_sub que va a ser de longitud 
		* 2m+a para poder calcular la relajacion de todos su elementos*/ 
		int c=0;
		for(i= a*(world_rank)-m ; i <m+(1+world_rank)*a ; i++)
		{
		V_sub[c]=V[i];
		c+=1;
		}

		/* calculo el promedio y lo guardo en V_new y hago el for para los valores que realmente le pertencen al procesadro y me quito las dos filas mas que agrego 
		* para calcular el promedio , asi que V_new es un array de longitud a */

		int r=0;
		for(i=m+1 ; i<a+1 ; i++)
		{
			up=i-m;
			down=i+m;
			left=i-1;
			right=i+1;
			average = (V_sub[up]+ V_sub[down] + V_sub[left] + V_sub[right])/4;
			V_new[r]=average;
			r+=1;
		} 
	
	}
	/* ahora hago el caso del ultimo procesador que tiene el final del arreglo V*/

	if(world_rank==world_size)
	{	
		/* seleciono los elementos con un fila mas del arreglo anteriro hasta la cantidad total de elementos de V que es mxm*/
		int p=0;
		for(i=a*world_rank-m ; i <m*m ; i++)
		{
			V_sub[p]=V[i];
			p+=1;
		}
		/*comienzo el for desde una fila despues para poder hacer la relajacion y hasta una fila antes del tamno total del arreglo*/
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
	/*Cada procesador tiene un Vectro V_new que tiene la relajacion de una sub parte del arreglo V entonces solicito que cada procesador envie y reciva el segmento correspondinte
	* de cada procesadro y los guarde en orden del procesador 0 - world_size en el array V1 y asi cada procesadro tendra la relajacion completa de todo los procesadores*/ 
	MPI_Allgather(V_new,a,MPI_DOUBLE,V1,a,MPI_DOUBLE,MPI_COMM_WORLD);
	/* en cada procesador paso todo V1 a V para poder hacer el proceso de nuevo y hacer la relajacion*/
	for(i=0;i<m*m;i++)
	{
	V[i]=V1[i]
	}
	n=+1;	
}	
/* imprimo el array V despues de hacer todas las iteraciones*/

for(i=0 ; i <m*m ;i++)
{
	printf("%f\n", V[i]);
}
MPI_Barrier(MPI_COMM_WORLD);
MPI_Finalize();
return 0;
}
