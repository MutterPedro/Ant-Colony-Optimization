#include <stdio.h>
#include <stdlib.h>
#include "omp.h"
#include <math.h>
#define TEOR_EVAPORACAO 0.3

int n;
double qualidade[n][n],feromonio[n][n];

typedef struct {
	double x;
	double y;
}Node;

typedef struct {
	Node* visitadas;
	double caminho;
}Ant;

Ant formigas[n];

double euclidiana(double x1, double y1,double x2, double y2){
	return (sqrt(pow(x1-x2,2)+pow(y1-y2,2)));
}

double somatorioFeromonioQualidade(){
	int i=0,j=0;
	double soma = 0;
	for(;i<n;){
		for(;j<n;j++){
			soma += feromonio[i][j]*qualidade[i][j];
		}
	}

	return soma;
}

double somatorioDistanciaPercorrida(){
	int i=0;
	double soma = 0;
	for(;i<n;){
		soma += 1/formigas[i].caminho;
	}

	return soma;
}

double probabilidade(int i,int j){
	return (
		(feromonio[i][j]*qualidade[i][j])/(somatorioFeromonioQualidade())
	);
}

double atualizarFeromonio(int i, int j){
	return(
		((1-TEOR_EVAPORACAO)*feromonio[i][j])+somatorioDistanciaPercorrida()
	);
}

int main(int argc, char** args){
	scanf("%d",&n);
	Node cidades[n];

	int i=0,j=0,k=0,l=0;
	for(i=0;i<n;i++){
		scanf("%lf",&cidades[i].x);
		scanf("%lf",&cidades[i].y);
		formigas[i].visitadas = malloc(n*sizeof(Node));
		formigas[i].visitadas[0] = cidades[i];
		formigas[i].caminho = 0;
	}

	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			if(i == j){
				qualidade[i][j] = -1;
				feromonio[i][j] = -1;
			} else {
				qualidade[i][j] = 1/euclidiana(cidades[i].x,cidades[j].x,cidades[i].y,cidades[j].y);
				feromonio[i][j] = 0;			
			}
		}
	}


	
	return 0;
}
