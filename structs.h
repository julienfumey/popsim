/**
 * AstyanaxPopSim
 *
 * \file structs.h
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief header of structs.c. Contains functions prototypes.
 * \copyright (C) 2014
 *
*/
#ifndef DEF_STRUCTS
#define DEF_STRUCTS

#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Prototypes
typedef struct Locus Locus;
typedef struct Param Param;
typedef struct Snp Snp;

void freeLocus(Locus*);
Locus* copyList(Locus*);

// Structures

/**
 * \struct Locus
 * \brief Element of linked list representing one locus
*/
struct Locus{
	double sf; /**< Frequency of derived alleles frequency for the locus in SF */
	double cf; /**< Frequency of derived alleles frequency for the locus in CF */
	double texas; /**< Frequency of derived alleles frequency for the locus in Texas Fish */
	long gene; /**< Generation of apparition of the locus */
	Locus *next; /**< Pointer to next locus */
};

/**
 * \struct Param
 * \brief Simulation parameters
*/
struct Param{
	unsigned long nbLocus; /**< Number of locus to generate in init pop */
	int tot_gen; /**< Number of simulations generations */
	int init_gen; /**< Number of initial generations (1 pop that accumulates mutations) */
	int elabra_gen; /**< Number of generation of generations before separation between CF and El Abra SF */
	int bn_gen; /**< Number of generation of bottleneck (CF) */
	int ancestral_popsize; /**< Size of the ancestral population */
	int texas_popsize; /**< Size of Texas population */
	int elabra_popsize; /**< Size of Sierra de El Abra population */
	int bn_popsize; /**< Size of CF population during bottleneck */
	int cf_popsize; /**< Size of CF population after bottleneck */
	double mutaRate; /**< Mutation rate */
	double STmr; /**< Migration rate from SF (El Abra) to SF (Texas) */
	double TSmr; /**< Migration rate from SF (Texas) to SF (El Abra) */
	double SCmr; /**< Migration rate from SF (El Abra) to CF */
	double CSmr; /**< Migration rate from CF to SF (El Abra) */
	double TCmr; /**< Migration rate from SF (Texas) to CF */
	double CTmr; /**< Migration rate from CF to SF (Texas) */
	double migra;
	double perc_poissonmigraTC; /**< Percentage of migrant fishes from Texas to cave */
	double perc_poissonmigraCT; /**< Percentage of migrant fishes from cave to Texas */
	double perc_poissonmigraST; /**< Percentage of migrant fishes from El Abra to Texas */
	double perc_poissonmigraTS; /**< Percentage of migrant fishes from Texas to El Abra */
	double perc_poissonmigraCS; /**< Percentage of migrant fishes from cave to El Abra */
	double perc_poissonmigraSC; /**< Percentage of migrant fishes from El Abra to cave */
	char* filename; /**< output filename */
	int generatePopInit; /**< Generate a initial population of DAF ? (1 = yes, 0 = no) */
	int labCF_gen; /**< Number of generations in lab for CF */
	int labSF_gen; /**< Number of generations in lab for SF */
	int labCF_popsize; /**< Population size of cavefish in lab */
	int labSF_popsize; /**< Population size of SF in lab*/
	char* filenameLab; /**< output filename for extra lab generations */
	unsigned long nb_snp0; /**< Number of SNPs that have lost polymorphysm for all populations*/
	int exportLocusPopGene;
	int exportLocusPop;
	int subGenTexas;
	int subGenSF;
	int subGenCF;

};

/**
 * \struct Snp
 * \brief Number of SNPs by categories
*/
struct Snp{
	unsigned int snp; /**< Number of SNPs position */
	int fixedCF; /**< Number of SNPs ancestral fixed SF, derived fixed CF */
	int fixedSF; /**< Number of SNPs ancestral fixed CF, derived fixed SF */
	int shared; /**< Number of SNPs shared in SF and CF */
	int derivedFixedCF; /**< Number of SNPs polymorphic SF, derived fixed CF */
	int derivedFixedSF; /**< Number of SNPs polymorphic CF, derived fixed SF */
	int ancestralFixedCF; /**< Number of SNPs polymorphic SF, ancestral fixed CF */
	int ancestralFixedSF; /**< Number of SNPs polymorphic CF, ancestral fixed SF */
	int ancestralFixedBothPop; /**< Number of SNPs ancestral fixed in SF and CF (derived allele is lost) */
	int derivedFixedBothPop; /**< Number of SNPs derived fixed in SF and CF (ancestral allele is lost) */
	int polyCF;
	int polySF;
	int polyTexas;
	int fixedInitCF;
	int fixedNeoCF;
	int fixedInitSF;
	int fixedNeoSF;
	double corr; /**<Correlation coefficient for shared SNP between CF and SF */
	double x;
	double y;
	double xy;
   	double xsq;
   	double ysq;
};

#endif
