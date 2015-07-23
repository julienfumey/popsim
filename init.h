/*
 * AstyanaxPopSim
 *
 * \file init.h
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief header of init.c
 * \copyright (C) 2014
 *
*/

#ifndef DEF_INIT
#define DEF_INIT

#include "structs.h"

//Prototypes
Locus* initPop(Param*, double*, long, gsl_rng*);
double freqInit(double*, long, gsl_rng*);
double* freqTable(Param*);
double watterson(long);
Locus* initializeSimulation(Param*, gsl_rng*);
int nbLocusInit(float, int);

#endif
