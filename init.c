/**
 * AstyanaxPopSim
 *
 * \file init.c
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief Initialization of population
 * \version 2.0
 * \copyright (C) 2014
 *
*/

#include "init.h"

/**
 * \fn Locus* initializeSimulation(Param* P)
 * \brief Create the linked list of locus to initialize simulation
 *
 * \Param P Simulation parameters
 * \returns pointer to the first element of the linked list of locus
*/
Locus* initializeSimulation(Param* P, gsl_rng* r){
	double* tableFrequence = freqTable(P);

	Locus* pop;

	if(P->generatePopInit == 1){
		if(P->mutaRate == 0 && P->nbLocus == 0){
			fprintf(stderr, "Error : You have to specify number of locus to generate or mutation rate\n");
			exit(0);
		}
		if(P->nbLocus == 0){
			P->nbLocus = nbLocusInit(P->mutaRate, P->texas_popsize);
		}
		printf("%ld\n", P->nbLocus);
		pop	= initPop(P, tableFrequence, P->nbLocus, r);
	}else{
		pop = calloc(1, sizeof(Locus));
		if(pop == NULL){
			fprintf(stderr, "ERROR : Cannot create initial population\n");
			exit(0);
		}
	}
	
	free(tableFrequence);

	return(pop);
}

/**
 * \fn int nbLocusInit(int mutaRate, int popsize)int nbLocusInit(int mutaRate, int popsize)
 * \brief calculate number of locus to generate at equilibrium
*/
int nbLocusInit(float mutaRate, int popsize){
	return (int) 4 * popsize * mutaRate * watterson(popsize);
}

/**
 * \fn Locus* initPop(Param* P, float* freq, long nbLocus)
 * \brief Initialize population
 *
 * \param P parameters
 * \param freq 
 * \param nbLocus Number of locus to create
 * \returns Pointer to the first element of the linked list of locus
*/
Locus* initPop(Param* P, double* freq, long nbLocus, gsl_rng* r){
	Locus* locus = calloc(1, sizeof(Locus));
	if(locus == NULL){
		fprintf(stderr, "Cannot allocate initial population\n");
		exit(0);
	}

	float freqLocus = freqInit(freq, P->texas_popsize, r);

	locus->cf = (P->elabra_gen == 0)?freqLocus:0;
	locus->sf = freqLocus;
	locus->texas = freqLocus;
	locus->gene = 0;
		
	if(nbLocus > 1){
		locus->next = initPop(P, freq, --nbLocus, r);
	}
	
	return locus;
}

/**
 * \fn inline float freqInit(float* table, long popSize)
 * \brief Give a initial frequency for a locus
 *
 * \param table distribution of DAF
 * \param popSize Size of population
 * \return Initial frequency of the locus
*/
inline double freqInit(double* table, long popSize, gsl_rng* r){
	//printf("%.38f\n", (float)random()/RAND_MAX); 
	//float f = (float)random()/RAND_MAX;
	double f = gsl_rng_uniform(r);
	int i;
	for(i = 1; i < 2*popSize; i++){
		if(f <= table[i]){
			return((double)i/(2*popSize));
		}
	}
	return 0;
}

/**
 * \fn float* freqTable(Param* P)
 * \brief Calculate the frequency of the 2N-1 allele
 *
 * \param P simulation parameters
 * \returns distribution of DAF
*/
double* freqTable(Param* P){
	double* table = calloc((2*P->texas_popsize),sizeof(double));
	if(table == NULL){
		fprintf(stderr, "Error at line %d of file %s\n", __LINE__, __FILE__);
		exit(0);
	}

	table[0] = 0.;
	
	double nbWatterson = watterson(P->texas_popsize);
	
	long i;
	for(i = 1; i < 2*P->texas_popsize; i++){
		double freqSite = (double) 1/(i*nbWatterson);
		table[i] = table[(i-1)] + freqSite;
	}
	
	return table;
}

/**
 * \fn float watterson(long popSize)
 * \brief Calculate Watterson coefficient
 *
 * \param popSize Size of population
 * \returns Watterson coefficient
*/
double watterson(long popSize){
	double nb = 0;
	int i;
	
	for(i = 1; i < 2*popSize; i++){
		nb += (double) 1/i;
	
	}	
	
	return nb;
}

