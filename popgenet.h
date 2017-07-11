/**
 * AstyanaxPopSim
 *
 * \file popgenet.h
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief header of popgenet.c. Contains functions prototypes.
 * \copyright (C) 2014
 *
*/

#ifndef DEF_POPGENET
#define DEF_POPGENET

#include "structs.h"

//Prototype
Locus* evolve(Locus*, Param*, gsl_rng*);
void newFreq(Locus*, Param*, long, int, gsl_rng*);
inline double freqAllele(double, long, gsl_rng*);
void evolveLab(FILE*, Locus*, Param*, long, gsl_rng*);
int max(int a, int b);
inline void loadBar(int, int, int, int);
Locus* mutate(Locus*, Param*, long);
long popSize(Param*, long);
Locus* deleteNonSnp(Locus*, Param*);
void migrate(Locus*, Param*, int, gsl_rng*);
double gaussrand(int m, int sd, gsl_rng*);
void migrateFreq(Locus*, int, int, int, int, int, int, int, double*, Param*, int);
#endif

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
