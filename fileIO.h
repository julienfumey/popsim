/*
 * AstyanaxPopSim
 *
 * \file fileIO.h
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief header of fileIO.c
 * \copyright (C) 2014
 *
*/

#ifndef DEF_FILEIO
#define DEF_FILEIO

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "structs.h"

// Macros
#define MAX_SIZE_PARAMETERS 600

void exportPop(Locus*, FILE*);
void printFreq(Locus*, FILE*, int);
Param* parseParametersFile(char*);
void addHeaderOutFile(FILE*);
void exportStats(FILE*, long, Locus*, Param*, gsl_rng*);
void statsPop(Locus*, Snp*);
void statsPop2(Locus*, Snp*);
void scoreSimu(Snp*, double*);
double calculScore(Snp*, int, int, int, int, int, int, int);
void nbSitePolyByPop(Locus*, Snp*);
void calcCorr(Snp* stats);
void exportLocusPop(Locus*, gsl_rng*);
void printFileExportLocus(FILE*, Locus*, gsl_rng*);

#endif
