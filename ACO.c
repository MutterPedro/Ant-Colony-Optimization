#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include "mpi.h"
#include <math.h>
#include <float.h>
#define TEOR_EVAPORACAO 0.2
#define FEROMONIO_INICIAL 1
#define Q 1
#define CONT_MAX 10
#define FEROMONIO_MINIMO 0.1

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
double menorDistancia=0;
int cont=0;

double euclidiana(double x1, double y1,double x2, double y2){
	return (sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
}

void inicializarValores();
void inicializarMatrizes();
double somatorioFeromonioQualidade();
double probabilidade(int i,int j);
void atualizarFeromonio(int i, int j, double valor);
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

	int proxima,i,j;
	formiga.visitadas = malloc((n+1)*sizeof(int));
	for(i =0;i<n;i++){
		MPI_Bcast(&cidades[i],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
		for(j=0;j<n;j++){
			MPI_Bcast(&qualidade[i][j],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
		}
	}
	int atual = rank;
	while(cont < CONT_MAX){

		//enviando as matrizes de qualidade e feromonio para os outros processos
		MPI_Barrier(MPI_COMM_WORLD);
		for(i =0;i<n;i++){
			for(j=0;j<n;j++){
				MPI_Bcast(&feromonio[i][j],1, MPI_DOUBLE,0, MPI_COMM_WORLD);
			}
		}
		
		//inicializando variaveis
		formiga.visitadas[0] = atual;
		for(i = 1;i<=n;i++){
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
			atualizarFeromonio(atual,proxima,formiga.caminho);
			atual = proxima;		
		}

		//parte em que a formiga volta para a cidade inicial
		formiga.visitadas[n] = formiga.visitadas[0];
		formiga.caminho += 1/qualidade[atual][formiga.visitadas[0]];
		atualizarFeromonio(atual,formiga.visitadas[0],formiga.caminho);
		atual = formiga.visitadas[0];

		//enviando os caminhos percorridos e vendo qual o menor	
		MPI_Send(formiga.visitadas, n+1, MPI_INT, 0, rank, MPI_COMM_WORLD);	
		if(rank == 0){
			double distancia,menor = DBL_MAX;
			int caminhos[n+1];
			Ant f;			
			for(i=0;i<n;i++){
				MPI_Recv(caminhos, n+1, MPI_INT, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, &status);
				f.visitadas = caminhos;
				f.caminho = 0;
				for(j=0;j<n;j++){
					f.caminho += 1/qualidade[caminhos[j]][caminhos[j+1]];
					//printf("atualizando feromonio em %d %d\n", caminhos[j],caminhos[j+1]);
					atualizarFeromonio(caminhos[j],caminhos[j+1],f.caminho);
				}
				distancia = f.caminho;
				if(distancia < menor){
					menor = distancia;
					for(j=0;j<n+1;j++){
						melhorCaminho[j] = caminhos[j];
					}
				}
			}
			if(menorDistancia == menor){
				cont++;
			} else {
				cont = 0;
			}
			menorDistancia = menor;
			/*for(i=0;i<n;i++){
				for(j=0;j<n;j++){
					printf("%lf ",feromonio[i][j]);
				}
				printf("\n");
			}*/
			printf("\n");
			evaporarFeromonio();			
		}
		MPI_Bcast(&menorDistancia,1, MPI_DOUBLE,0, MPI_COMM_WORLD);
		MPI_Bcast(&cont,1, MPI_DOUBLE,0, MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);		
		printf("A formiga %d visitou: [",rank);
		for(i =0;i<=n;i++){
			printf("%d ",formiga.visitadas[i]);
		}
		printf("] distancia = %lf\n",formiga.caminho);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank == 0){
		int i;
		printf("Melhor caminho obtido foi: ");
		for(i =0;i<n+1;i++){
			printf("%d ",melhorCaminho[i]);
		}
		printf(" -> distancia:%lf\n", menorDistancia);
	}

	liberarMemoria();
	
	MPI_Finalize();	
	
	return 0;
}

void inicializarValores(){
	scanf("%d",&n);
	cidades = malloc(n*sizeof(Node));
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
				qualidade[i][j] = 1/euclidiana(cidades[i].x,cidades[i].y,cidades[j].x,cidades[j].y);
				feromonio[i][j] = FEROMONIO_INICIAL;			
			}
		}
	}
}

void inicializarMatrizes(){
	melhorCaminho = malloc((n+1)*sizeof(int));
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

void atualizarFeromonio(int i, int j,double valor){
	feromonio[i][j]+=Q/valor;
}

void evaporarFeromonio(){
	int i,j;
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			feromonio[i][j] -= FEROMONIO_INICIAL*(1-TEOR_EVAPORACAO);
			if(feromonio[i][j] < FEROMONIO_MINIMO)
				feromonio[i][j] = FEROMONIO_MINIMO;
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