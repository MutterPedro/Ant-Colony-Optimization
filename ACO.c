#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "mpi.h"
#include <math.h>
#define TEOR_EVAPORACAO 0.3
#define Q 1

int n;
double** qualidade;
double** feromonio;

typedef struct {
	double x;
	double y;
}Node;

typedef struct {
	int* visitadas;
	double caminho;
}Ant;

Node* cidades;
Ant formiga;

double euclidiana(double x1, double y1,double x2, double y2){
	return (sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
}

double somatorioFeromonioQualidade();

double somatorioDistanciaPercorrida();

double probabilidade(int i,int j);

double atualizarFeromonio(int i, int j, int somar);

int jaVisitou(int formigaIdx, int cidadeIdx);

int main(int argc, char** argv){
	int size, rank, rc;
	MPI_Status status;

	MPI_Init(&argc,&argv); 
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	fflush(stdout);
	
	if(rank == 0){
		scanf("%d",&n);
		cidades = malloc(n*sizeof(Node));		

		int i=0,j=0;
		//inicializando matriz de distancias
		for(i=0;i<n;i++){
			scanf("%lf",&cidades[i].x);
			scanf("%lf",&cidades[i].y);			
		}

		//inicializando matriz qualidade e feromonio
		feromonio = malloc(n*sizeof(double*));
		qualidade = malloc(n*sizeof(double*));
		
		for(i=0;i<n;i++){
			feromonio[i] = malloc(n*sizeof(double));
			qualidade[i] = malloc(n*sizeof(double));
		}
		for(i=0;i<n;i++){
			for(j=0;j<n;j++){
				if(i == j){
					qualidade[i][j] = -1;
					feromonio[i][j] = -1;
				} else {
					qualidade[i][j] = 1/euclidiana(cidades[i].x,cidades[j].x,cidades[i].y,cidades[j].y);
					feromonio[i][j] = 1;			
				}
			}
		}		
	}
	int proxima,i,j;
	MPI_Bcast(&n,1, MPI_INT,0, MPI_COMM_WORLD);

	if(rank > 0){
		feromonio = malloc(n*sizeof(double*));
		qualidade = malloc(n*sizeof(double*));
		
		for(i=0;i<n;i++){
			feromonio[i] = malloc(n*sizeof(double));
			qualidade[i] = malloc(n*sizeof(double));
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	for(i =0;i<n;i++){
		for(j=0;j<n;j++){
			MPI_Bcast(&qualidade[i][j],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
			MPI_Bcast(&feromonio[i][j],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
		}
	}	

	formiga.visitadas = malloc(n*sizeof(int));
	formiga.visitadas[0] = rank;
	double prob = 0;
	for(i =1;i<n;i++){
		for(j=0;j<n;j++){			
			double p = probabilidade(rank,j);
			if(p > prob && rank != j && jaVisitou(j) == 0){		
				prob = p;
				proxima = j;
			}
		}
		formiga.visitadas[i] = proxima;
	}
	
	printf("A formiga %d visitou: [",rank);
	for(i =1;i<n;i++){
		printf("%d ",formiga.visitadas[i]);
	}
	printf("]\n");	
	
	MPI_Finalize();
	
	
	return 0;
}

double somatorioFeromonioQualidade(){
	int i=0,j=0;
	double soma = 0;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			soma += feromonio[i][j]*qualidade[i][j];
		}
	}

	return soma;
}

double somatorioDistanciaPercorrida(){
	int i=0;
	double soma = 0;
	for(i=0;i<n;i++){
		//soma += Q/formigas[i].caminho;
	}

	return soma;
}

double probabilidade(int i,int j){
	return (
		(feromonio[i][j]*qualidade[i][j])/(somatorioFeromonioQualidade())
	);
}

double atualizarFeromonio(int i, int j, int somar){
	return(
		((1-TEOR_EVAPORACAO)*feromonio[i][j])+(somatorioDistanciaPercorrida()*somar)
	);
}

int jaVisitou(int cidadeIdx){
	int i=0;
	for(i=0;i<n;i++){
		if(formiga.visitadas[i] == cidadeIdx){
			return 1;
		}
	}	
	
	return 0;
}
