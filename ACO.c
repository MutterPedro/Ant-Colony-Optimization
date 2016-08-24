#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "mpi.h"
#include <math.h>
#include <float.h>
#define TEOR_EVAPORACAO 0.3
#define FEROMONIO_INICIAL 1
#define Q 1
#define CONT_MAX 5

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

int* melhorCaminho;
int menorDistancia=0, cont=0;

double euclidiana(double x1, double y1,double x2, double y2){
	return (sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
}

void inicializarValores();
void inicializarMatrizes();
double somatorioFeromonioQualidade();
double probabilidade(int i,int j);
void atualizarFeromonio(int i, int j);
void evaporarFeromonio();
int jaVisitou(int cidadeIdx);
void liberarMemoria();
void enviarMatrizFeromonios(int tag);
void receberMatrizFeromonios(int rank);

MPI_Status status;
int main(int argc, char** argv){
	int size, rank, rc;

	MPI_Init(&argc,&argv); 
	MPI_Comm_size (MPI_COMM_WORLD, &size); 
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	fflush(stdout);
	
	if(rank == 0){		
		inicializarValores();
	}

	//enviando N para os outros processos
	MPI_Bcast(&n,1, MPI_INT,0, MPI_COMM_WORLD);

	if(rank > 0){
		cidades = malloc(n*sizeof(Node));
		inicializarMatrizes();
	}	

	formiga.visitadas = malloc((n+1)*sizeof(int));
	int atual = rank;
	//while(cont < CONT_MAX){
	for(int k=0;k<5;k++){

		int proxima,i,j;
		//enviando as matrizes de qualidade e feromonio para os outros processos
		MPI_Barrier(MPI_COMM_WORLD);
		for(i =0;i<n;i++){
			MPI_Bcast(&cidades[i],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
			for(j=0;j<n;j++){
				MPI_Bcast(&qualidade[i][j],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
				MPI_Bcast(&feromonio[i][j],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
			}
		}
		
		//inicializando variaveis
		formiga.visitadas[0] = atual;
		for(i = 1;i<n;i++){
			formiga.visitadas[i] = -1;
		}
		formiga.caminho = 0;

		//calculando probabilidades e atualizando a tabela de feromonios
		double prob = 0;		
		for(i =1;i<n;i++){
			prob = 0;			
			for(j=0;j<n;j++){		
				double p = probabilidade(atual,j);
				if((p > prob) && (atual != j) && (jaVisitou(j) == 0)){				
					prob = p;
					proxima = j;
				}
			}
			formiga.visitadas[i] = proxima;		
			formiga.caminho += 1/qualidade[atual][proxima];
			atualizarFeromonio(atual,proxima);
			atual = proxima;		
		}

		//parte em que a formiga volta para a cidade inicial
		formiga.visitadas[n] = formiga.visitadas[0];
		formiga.caminho += 1/qualidade[atual][formiga.visitadas[0]];

		//enviando os caminhos percorridos e vendo qual o menor
		MPI_Send(&formiga.caminho, 1, MPI_DOUBLE, 0, rank, MPI_COMM_WORLD);		
		//MPI_Send(formiga.visitadas, n, MPI_INT, 0, rank, MPI_COMM_WORLD);	
		if(rank == 0){
			double distancia,menor = DBL_MAX;
			//int caminhos[n];		
			for(i=0;i<n;i++){
				MPI_Recv(&distancia, 1, MPI_DOUBLE, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				//MPI_Recv(caminhos, n, MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				if(distancia < menor){
					menor = distancia;
					//melhorCaminho = caminhos;
				}
			}
			printf("%lf\n", menor);
			if(menorDistancia == menor){
				cont++;
			} else {
				menorDistancia = menor;
				cont = 0;
			}		
		}

		MPI_Barrier(MPI_COMM_WORLD);
		enviarMatrizFeromonios(rank);
		receberMatrizFeromonios(rank);
		
		printf("A formiga %d visitou: [",rank);
		for(i =0;i<=n;i++){
			printf("%d ",formiga.visitadas[i]);
		}
		printf("] distancia = %lf\n",formiga.caminho);
	}

	liberarMemoria();
	
	MPI_Finalize();	
	
	return 0;
}

void inicializarValores(){
	scanf("%d",&n);
	cidades = malloc(n*sizeof(Node));
	melhorCaminho = malloc(n*sizeof(int));	

	int i=0,j=0;
	//inicializando matriz de distancias
	for(i=0;i<n;i++){
		scanf("%lf",&cidades[i].x);
		scanf("%lf",&cidades[i].y);			
	}
	inicializarMatrizes();
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i == j){
				qualidade[i][j] = -1;
				feromonio[i][j] = -1;
			} else {
				qualidade[i][j] = 1/euclidiana(cidades[i].x,cidades[j].x,cidades[i].y,cidades[j].y);
				feromonio[i][j] = FEROMONIO_INICIAL;			
			}
		}
	}
}

void inicializarMatrizes(){
	int i=0,j=0;
	//inicializando matriz qualidade e feromonio
	feromonio = malloc(n*sizeof(double*));
	qualidade = malloc(n*sizeof(double*));	
	for(i=0;i<n;i++){
		feromonio[i] = malloc(n*sizeof(double));
		qualidade[i] = malloc(n*sizeof(double));
	}
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

double probabilidade(int i,int j){
	return (
		(feromonio[i][j]*qualidade[i][j])/(somatorioFeromonioQualidade())
	);
}

void atualizarFeromonio(int i, int j){
	feromonio[i][j]+=Q/formiga.caminho;
}

void evaporarFeromonio(){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			feromonio[i][j] *= (1-TEOR_EVAPORACAO);
			if(feromonio[i][j] < 0)
				feromonio[i][j] = 0;
		}
	}
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

void enviarMatrizFeromonios(int tag){
	int i,j;
	for(i =0;i<n;i++){
		for(j=0;j<n;j++){
			MPI_Send(&feromonio[i][j], 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
		}
	}
}

void receberMatrizFeromonios(int rank){
	int i,j;
	if(rank == 0){
		for(i =0;i<n;i++){
			for(j=0;j<n;j++){
				double feromonioRecv;
				MPI_Recv(&feromonioRecv, 1, MPI_DOUBLE, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				if(feromonioRecv != feromonio[i][j]){
					feromonio[i][j] += feromonioRecv-FEROMONIO_INICIAL;
				}
			}
		}
		evaporarFeromonio();
	}
}

void liberarMemoria(){
	free(melhorCaminho);
	free(cidades);
	free(formiga.visitadas);
	int i;
	for(i=0;i<n;i++){
		free(qualidade[i]);
		free(feromonio[i]);
	}
	free(qualidade);
	free(feromonio);
}