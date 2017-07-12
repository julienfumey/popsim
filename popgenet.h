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
extern inline double freqAllele(double, long, gsl_rng*);
void evolveLab(FILE*, Locus*, Param*, long, gsl_rng*);
int max(int a, int b);
extern inline void loadBar(int, int, int, int);
Locus* mutate(Locus*, Param*, long);
long popSize(Param*, long);
Locus* deleteNonSnp(Locus*, Param*);
void migrate(Locus*, Param*, int, gsl_rng*);
double gaussrand(int m, int sd, gsl_rng*);
void migrateFreq(Locus*, int, int, int, int, int, int, int, double*, Param*, int);
#endif
