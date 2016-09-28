/**
 * AstyanaxPopSim
 *
 * \file popgenet.c
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief Functions of population genetics
 * \copyright (C) 2014
 * \version 2.0
 *
*/

#include "popgenet.h"
#include "fileIO.h"

/**
 * \fn int max(int a, int b)
 * \brief returns the max between a and b
 * \Param a
 * \Param b
 * \returns max
*/
int max(int a, int b){
	return a > b ? a : b;
}

/**
 * \fn Locus* evolve(Locus* pop, Param* P)
 * \brief Function that make evolve the population
 * \Param pop Starting population
 * \Param P Simulation parameters
 * \returns Final population after simulation
*/
Locus* evolve(Locus* pop, Param* P, gsl_rng* r){
	long i;
	FILE* outfile = fopen(P->filename, "w");
	FILE* outfileLab = fopen(P->filenameLab, "w");

	addHeaderOutFile(outfile);
	addHeaderOutFile(outfileLab);
	
	//ProgressBar
	loadBar(0, P->tot_gen, 100, 80); 
	
	FILE* f = fopen("freq_gen_0", "w");
	exportPop(pop, f);
	fclose(f);

	for(i = 1; i <= P->tot_gen; i++){
		printf("%ld\n", i);
		newFreq(pop, P, i, 0, r);
		pop = mutate(pop, P, i);
		if(i > P->init_gen && P->migra == 1){
		//printf("%ld: %ld\n", i, P->nbLocus);
			migrate(pop, P, i, r);
		}

		pop = deleteNonSnp(pop, P);
	
		loadBar(i, P->tot_gen, 100, 80); 
		if((i-1)%10 == 0){
			exportStats(outfile, i, pop, P, r);
			if(i >= (P->init_gen + P->elabra_gen)){
				evolveLab(outfileLab, pop, P, i, r);
			}
		}	


		if(i == 500 || i == 1500 || i == 3000 || i == 5000 || (i > 5000 && i%10000 == 0)){
			char filename[50];
			sprintf(filename, "freq_gen_%ld", i);
			FILE* f = fopen(filename, "w");
			exportPop(pop, f);
			fclose(f);
		}

		
	}
	
	fclose(outfile);
	fclose(outfileLab);
	return pop;
}

/**
 * \fn void evolveLab(FILE* f, Locus* pop, Param* P, long i)
 * \brief Genetic drift for lab fishes
 * \Param f output file handle
 * \Param pop Linked list of locus
 * \Param P Simulation parameters
 * \Param i generation
*/
void evolveLab(FILE* f, Locus* pop, Param* P, long i, gsl_rng* r){
	int nbGeneration = max(P->labCF_gen, P->labSF_gen);
	Locus* pop2 = copyList(pop);

	int j;
	for(j = 1; j <= nbGeneration; j++){
		newFreq(pop2, P, i, 1, r);	
	}

	exportStats(f, i, pop2, P, r);
	
	if(i == P->exportLocusPopGene && P->exportLocusPop == 1){
		exportLocusPop(pop2, r);
	}
	freeLocus(pop2);
}

/**
 * \fn void newFreq(Locus* pop, Param* P, long gen)
 * \brief Call freqAllele for each population of each locus
 * \Param pop Pointer of the linked list of locus
 * \Param P Simulation parameters
 * \Param gen Simulation generation
 * \returns void
*/
void newFreq(Locus* pop, Param* P, long gen, int lab, gsl_rng* r){
	long cfPopSize, sfPopSize, texasPopSize;
	int i;
	if(lab == 0){
		cfPopSize = popSize(P, gen);
		sfPopSize = P->elabra_popsize;
		texasPopSize = P->texas_popsize;

		if(gen == (P->elabra_gen + P->init_gen)){
			if(P->elabra_gen == 0){
				pop->cf = pop->texas;
			}else{
				pop->cf = pop->sf;
			}
		}
		
		for(i = 0; i < P->subGenCF; i++){
			pop->cf = freqAllele(pop->cf, cfPopSize, r);
		}

		if(gen == P->init_gen && P->elabra_gen != 0){
			pop->sf = pop->texas;
		}
		for(i = 0; i < P->subGenSF; i++){
			pop->sf = freqAllele(pop->sf, sfPopSize, r);
		}

		for(i = 0; i < P->subGenTexas; i++){
			pop->texas = freqAllele(pop->texas, texasPopSize, r);
		}
	
		if(pop->next != NULL){
			newFreq(pop->next, P, gen, lab, r);
		}
	}else if(lab == 1){ // lab == 1 when we perform genetic drift for fishes in laboratory.
		cfPopSize = P->labCF_popsize;
		sfPopSize = 0;
		texasPopSize = P->labSF_popsize;
		
		pop->texas = freqAllele(pop->texas, texasPopSize, r);
		pop->sf = freqAllele(pop->sf, sfPopSize, r);
		pop->cf = freqAllele(pop->cf, cfPopSize, r);

		if(pop->next != NULL){
			newFreq(pop->next, P, gen, lab, r);
		}
	}
}

/**
 * \fn inline double freqAllele(double freq, long popSize)
 * \brief calculate new frequency based on generation n-1 frequency and pop size
 * \Param freq Frequency at generation n-1
 * \Param popSize Population size
 * \returns new frequency
*/
/*inline double freqAllele(double freq, long popSize, gsl_rng* r){
	if(freq == 0 || freq == 1){
		return freq;
	}
	
	int nbZeros = 0;
	int i;
	
	for(i = 0; i < 2*popSize; i++){
		//double k = (double)random()/(double)RAND_MAX;
		double k = gsl_rng_uniform(r);
		nbZeros += (k<freq)?1:0;
	}

	return((double)nbZeros/(2*popSize));
}*/

inline double freqAllele(double freq, long popSize, gsl_rng* r){
	return((double) gsl_ran_binomial(r, freq, 2*popSize)/(2*popSize));
}


/**
 * \fn inline void loadBar(int x, int n, int r, int w)
 * \brief Print and update the progress bar
 * \Param x step
 * \Param n total number of step
 * \Param r precision
 * \Param w width of the bar
 * \remark function shamefully stolen from http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/
*/
inline void loadBar(int x, int n, int r, int w){
    // Only update r times.
    if ( x % (n/r) != 0 ) return;
 
    // Calculuate the ratio of complete-to-incomplete.
    double ratio = x/(double)n;
    int   c     = ratio * w;
 
    // Show the percentage complete.
    printf("%3d%% [", (int)(ratio*100) );
 
    // Show the load bar.
    int i;
	for (i=0; i<c; i++)
       printf("=");
 
    for (i=c; i<w; i++)
       printf(" ");
 
    // ANSI Control codes to go back to the
    // previous line and clear it.
    printf("]\n\033[F\033[J");
}

/**
 * \fn Locus* mutate(Locus* pop, Param* P, long gen
 * \brief Create new locus
 * \Param pop pointer to locus linked list
 * \Param P Simulation parameters
 * \Param gen generation
 * \returns Pointer to the new locus linked list
 *
 * The function adds 2N * µ new locus at the beginning
 * of the linked list for each population that exists at generation gen.
*/
Locus* mutate(Locus* pop, Param* P, long gen){
	int i;

	for(i = 0; i < (P->mutaRate * 2*P->texas_popsize * P->subGenTexas); i++){
		Locus* newLocus = calloc(1, sizeof(Locus));
		if(newLocus == NULL){	
			fprintf(stderr, "error on %s at line %d", __FILE__, __LINE__);
			exit(0);
		}
		newLocus->texas = (double)1/(double)(2*P->texas_popsize);
		newLocus->gene = gen;
		newLocus->next = pop;
		P->nbLocus++;
		pop = newLocus;
	}

	if(gen >= P->init_gen){
		for(i = 0; i < (P->mutaRate * 2*P->elabra_popsize * P->subGenSF); i++){
			Locus* newLocus = calloc(1, sizeof(Locus));
			
			if(newLocus == NULL){	
				fprintf(stderr, "error on %s at line %d", __FILE__, __LINE__);
				exit(0);
			}

			newLocus->sf = (double)1/(double)(2*P->elabra_popsize);
			newLocus->gene = gen;
			newLocus->next = pop;
			P->nbLocus++;
			pop = newLocus;
		}
	}
	
	if(gen >= (P->init_gen + P->elabra_gen)){
		// Mutation for CF ?
		long cfPopSize = popSize(P, gen);
		for(i = 0; i < (P->mutaRate * 2*cfPopSize * P->subGenCF); i++){
			Locus* newLocus = calloc(1, sizeof(Locus));
			
			if(newLocus == NULL){	
				fprintf(stderr, "error on %s at line %d", __FILE__, __LINE__);
				exit(0);
			}

			newLocus->cf = (double)1/(double)(2*cfPopSize);
			newLocus->gene = gen;
			newLocus->next = pop;
			P->nbLocus++;
			pop = newLocus;
		}	
	}
	
	return pop;	
}

/**
 * \fn long popSize(Param* P, long gen)
 * \brief Give the size of the CF population at generation gen
 * \Param P Simulation parameters
 * \Param gen the generation to be checked
*/
long popSize(Param* P, long gen){
	long size = 0;
	if(gen >= (P->init_gen + P->elabra_gen + P->bn_gen)){
		size = P->cf_popsize;
	}else if(gen >= (P->init_gen + P->elabra_gen)){
		size = P->bn_popsize;
	}

	return size;
}

/**
 * \fn Locus* deleteNonSnp(Locus* pop, Param* P)
 * \brief Delete locus that have lost derived alleles in all populations
 * \Param pop Locus linked list
 * \Param P Simulation parameters
 * \returns a pointer to the new locus linked list
*/
Locus* deleteNonSnp(Locus* pop, Param* P){
	if(pop->next != NULL){
		pop->next = deleteNonSnp(pop->next, P);
	}

	while((pop->cf == 0 && pop->sf == 0 && pop->texas == 0)){
		Locus* sup = pop;
		P->nb_snp0++;
		if(pop->next != NULL){
			pop = pop->next;
			free(sup);
			sup = NULL;
		}else{
			free(sup);
			sup = NULL;
			return NULL;
		}
	}
	return pop;
}

/**
*/
void migrate(Locus* pop, Param* P, int gen, gsl_rng* r){
	double alea[6];
	int i = 0;
	for(i = 0; i < 6; i++){
		//alea[i] = (double)rand() / RAND_MAX;
		alea[i] = gsl_rng_uniform(r);
	}
	
	int cfSize = popSize(P, gen);
	
	if(P->elabra_popsize == 0){ // Migration between Texas and Pachon only
		int migTtoC = (int)floor(gaussrand(P->perc_poissonmigraTC*P->texas_popsize, 0, r));
		int migCtoT = (int)floor(gaussrand(P->perc_poissonmigraCT*cfSize, 0, r));
		migrateFreq(pop, cfSize, 0, 0, migTtoC, migCtoT, 0, 0, alea, P, gen);
	}else{
		int migTtoS = (int)floor(gaussrand(P->perc_poissonmigraTS*P->texas_popsize, 0, r));	
		int migStoT = (int)floor(gaussrand(P->perc_poissonmigraST*P->elabra_popsize, 0, r));
		int migStoC = (gen > (P->init_gen + P->elabra_gen)) ? (int)floor(gaussrand(P->perc_poissonmigraSC*P->elabra_popsize, 0, r)) : 0; //if no cf pop, no migration from and to CF	
		int migCtoS = (gen > (P->init_gen + P->elabra_gen)) ? (int)floor(gaussrand(P->perc_poissonmigraCS*cfSize, 0, r)) : 0;
		migrateFreq(pop, cfSize, migTtoS, migStoT, 0, 0, migCtoS, migStoC, alea, P, gen);
	}
}

void migrateFreq(Locus* pop, int cfPopSize, int migTtoS, int migStoT, int migTtoC, int migCtoT, int migCtoS, int migStoC, double* alea, Param* P, int gen){	
	if(P->elabra_popsize == 0){
		double texasNewFreq = (double)((2 * migCtoT * pop->cf)+(2 * P->texas_popsize * pop->texas))/(2 * migCtoT + 2 * P->texas_popsize);
		double cfNewFreq = (double)((2 * migTtoC * pop->texas)+(2 * cfPopSize * pop->cf))/(2 * migTtoC + 2 * cfPopSize);
	
		if(alea[0] < P->CTmr){
			//printf("MIGRATION CT\n");
			pop->texas = texasNewFreq;
		}
		if(alea[1] < P->TCmr){
			//printf("MIGRATION TC\n");
			pop->cf = cfNewFreq;
		}
		//printf("toto 3\n");
		if(pop->next != NULL){
			migrateFreq(pop->next, cfPopSize, 0, 0, migTtoC, migCtoT, 0, 0, alea, P, gen);
		}
	}else{
		double texasNewFreq = (double)((2 * migStoT * pop->sf)+(2 * P->texas_popsize * pop->texas))/(2 * migStoT + 2 * P->texas_popsize);
		double cfNewFreq = (double)((2 * migStoC * pop->sf)+(2 * cfPopSize * pop->cf))/(2 * migStoC + 2 * cfPopSize);
		//double sfNewFreq;

		//Migration to surface
		if(alea[2] < P->TSmr && alea[5] > P->CSmr){
			//printf("MIGRATION\n");
			pop->sf = (double)((2 * P->elabra_popsize * pop->sf) + (2 * migTtoS * pop->texas) + (2 * migCtoS * pop->cf)) / (2 * migTtoS + 2 * migCtoS + 2 * P->elabra_popsize);
		}else if(alea[2] < P->TSmr){
			//printf("MIGRATION\n");
			pop->sf = (double)((2 * P->elabra_popsize * pop->sf) + (2 * migTtoS * pop->texas)) / (2 * migTtoS + 2 * P->elabra_popsize);
		}else if(alea[5] > P->CSmr){
			//printf("MIGRATION\n");
			pop->sf = (double)((2 * P->elabra_popsize * pop->sf) + (2 * migCtoS * pop->cf)) / (2 * migCtoS + 2 * P->elabra_popsize);
		}
		
		//Migration to texas
		if(alea[3] < P->STmr){
			//printf("MIGRATION ST\n");
			pop->texas = texasNewFreq;
		}

		//Migration to cave
		if(alea[4] < P->SCmr){
			//printf("MIGRATION ST\n");
			pop->cf = cfNewFreq;
		}
		
		if(pop->next != NULL){
			migrateFreq(pop->next, cfPopSize, migTtoS, migStoT, 0, 0, migCtoS, migStoC, alea, P, gen);
		}
	}
}
	/*
	else{ //No migration between Texas and Pachon

		if(alea < P->STmr){
			nbFish = gaussrand
		}	
*/


/**
 //
double gaussrand(int m, int sd){
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	
	X *= sd;
	X += m;

	return X;
}*/

double gaussrand(int m, int sd, gsl_rng* r){
	return(gsl_ran_gaussian(r, sd) + m);
}
