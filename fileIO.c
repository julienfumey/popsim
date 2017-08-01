/*
 * AstyanaxPopSim
 *
 * \file fileIO.h
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \brief header of fileIO.c
 * \copyright (C) 2014
 *
*/
#include "fileIO.h"

/**
 * \fn void exportPop(Locus* pop, FILE* f)
 * \brief export the frequency of derived alleles for each locus of texas fish in file f
 * \Param pop the population to export
 * \Param f file handle
 */
void exportPop(Locus* pop, FILE* f){
	fprintf(f, "%f\t%f\n", pop->texas, pop->cf);
	if(pop->next != NULL){
		exportPop(pop->next, f);
	}
}

/*void printFreq(Locus* pop, FILE* f, int p){
	if(pop->next != NULL){
		printFreq(pop->next, f, p);
	}

	if(p == 0){
		fprintf(f, "\t%f", pop->texas);
	}else if(p == 1){
		fprintf(f, "\t%f", pop->sf);
	}else if(p == 2){
		fprintf(f, "\t%f", pop->cf);
	}
}
*/

/**
 * \fn Param *parseParametersFile(char* fileName)
 * \brief This function parses parameters file and write the paramaters into a variable
 * \Param fileName Name of the file
 * \returns a variable of type Param that contains all parameters.
*/
Param *parseParametersFile(char* fileName){
	FILE* fichier = fopen(fileName, "r");

	if(fichier == NULL){
		fprintf(stderr, "Unable to open parameters file\n");
		exit(1);
	}
	Param *parametres = calloc(1, sizeof(Param));
	char line[2000];
	char param[1000], value[1000];

	while(fgets(line, sizeof(line), fichier) != NULL){
		if(!strchr("#\n", line[0])){
			sscanf(line, "%s = %s", param, value);
			if(strcmp(param, "nb_locus") == 0){
				parametres->nbLocus = strtol(value, NULL, 10);
			}else if(strcmp(param, "nb_gen") == 0){
				parametres->tot_gen = strtol(value, NULL, 10);
			}else if(strcmp(param, "init_gen") == 0){
				parametres->init_gen = strtol(value, NULL, 10);
			}else if(strcmp(param, "elabra_gen") == 0){
				parametres->elabra_gen = strtol(value, NULL, 10);
			}else if(strcmp(param, "bn_gen") == 0){
				parametres->bn_gen = strtol(value, NULL, 10);
			}else if(strcmp(param, "ancestral_popsize") == 0){
				parametres->ancestral_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "texas_popsize") == 0){
				parametres->texas_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "elabra_popsize") == 0){
				parametres->elabra_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "bn_popsize") == 0){
				parametres->bn_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "cf_popsize") == 0){
				parametres->cf_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "labCF_gen") == 0){
				parametres->labCF_gen = strtol(value, NULL, 10);
			}else if(strcmp(param, "labSF_gen") == 0){
				parametres->labSF_gen = strtol(value, NULL, 10);
			}else if(strcmp(param, "labCF_popsize") == 0){
				parametres->labCF_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "labSF_popsize") == 0){
				parametres->labSF_popsize = strtol(value, NULL, 10);
			}else if(strcmp(param, "mutaRate") == 0){
				parametres->mutaRate = strtof(value, NULL);
			}else if(strcmp(param, "STmr") == 0){
				parametres->STmr = strtof(value, NULL);
			}else if(strcmp(param, "TSmr") == 0){
				parametres->TSmr = strtof(value, NULL);
			}else if(strcmp(param, "SCmr") == 0){
				parametres->SCmr = strtof(value, NULL);
			}else if(strcmp(param, "CSmr") == 0){
				parametres->CSmr = strtof(value, NULL);
			}else if(strcmp(param, "TCmr") == 0){
				parametres->TCmr = strtof(value, NULL);
			}else if(strcmp(param, "CTmr") == 0){
				parametres->CTmr = strtof(value, NULL);
			}else if(strcmp(param, "migra") == 0){
				parametres->migra = strtof(value, NULL);
			}else if(strcmp(param, "migrantfish_percTC") == 0){
				parametres->perc_poissonmigraTC = strtof(value, NULL);
			}else if(strcmp(param, "migrantfish_percCT") == 0){
				parametres->perc_poissonmigraCT = strtof(value, NULL);
			}else if(strcmp(param, "migrantfish_percTS") == 0){
				parametres->perc_poissonmigraTS = strtof(value, NULL);
			}else if(strcmp(param, "migrantfish_percST") == 0){
				parametres->perc_poissonmigraST = strtof(value, NULL);
			}else if(strcmp(param, "migrantfish_percSC") == 0){
				parametres->perc_poissonmigraSC = strtof(value, NULL);
			}else if(strcmp(param, "migrantfish_percCS") == 0){
				parametres->perc_poissonmigraCS = strtof(value, NULL);
			}else if(strcmp(param, "resultfile") == 0){
				parametres->filename = (char*)calloc(1, strlen(value)+1);
				strcpy(parametres->filename, value);
			}else if(strcmp(param, "resultfileLab") == 0){
				parametres->filenameLab = (char*)calloc(1, strlen(value)+1);
				strcpy(parametres->filenameLab, value);
			}else if(strcmp(param, "generatePopInit") == 0){
				parametres->generatePopInit = strtol(value, NULL, 10);
			}else if(strcmp(param, "exportLocusPopGene") == 0){
				parametres->exportLocusPopGene = strtol(value, NULL, 10);
			}else if(strcmp(param, "exportLocusPop") == 0){
				parametres->exportLocusPop = strtol(value, NULL, 10);
			}else if(strcmp(param, "subGenSF") == 0){
				parametres->subGenSF = strtol(value, NULL, 10);
			}else if(strcmp(param, "subGenCF") == 0){
				parametres->subGenCF = strtol(value, NULL, 10);
			}else if(strcmp(param, "subGenTexas") == 0){
				parametres->subGenTexas = strtol(value, NULL, 10);
			}else{
				printf("Parameter %s does not exist\n", param);
			}
		}
	}
	fclose(fichier);
	return parametres;
}

/**
 * \fn void addHeaderOutFile(FILE* f)
 * \brief add header to output file
 * \Param f file handler
*/
void addHeaderOutFile(FILE* f){
	fprintf(f, "\"generations\", \"nb locus\", \"Polymorphism lost\", \"nb snp\",  \"fixé CF\", \"fixé SF\", \"partagé\", \"ancestral fixé CF\", \"ancestral fixé SF\", \"dérivé fixé CF\", \"dérive fixé SF\", \"Ancestral CF+SF\", \"derivé CF+SF\", \"Score NS\", \"Score S\", \"Score NC\", \"fixé CF\", \"fixé SF\", \"partagé\", \"ancestral fixé CF\", \"ancestral fixé SF\", \"dérivé fixé CF\", \"dérive fixé SF\", \"Ancestral CF+SF\", \"derivé CF+SF\", \"Score NS\", \"Score S\", \"Score NC\", \"Poly CF\", \"Poly SF\", \"Poly Texas\", \"Fixed init CF\", \"Fixed neo CF\", \"Fixed init SF\", \"Fixed neo SF\",\"Correlation inf5\", \"Correlation\"\n");
}

/**
 * \fn void exportStats(FILE* f, long i, Locus* pop)
 * \brief Write stats for generation i in f
 * \Param f file handler
 * \Param i generation
 * \Param pop Linked list of locus
*/
void exportStats(FILE* f, long i, Locus* pop, Param* P, gsl_rng* r){
	Snp* stats = calloc(1, sizeof(Snp));
	Snp* stats2 = calloc(1, sizeof(Snp));
	double *score = calloc(3, sizeof(double));
	double *score2 = calloc(3, sizeof(double));


	if(score == NULL || score2 == NULL){
		fprintf(stderr, "Error at line %d of file %s\n", __LINE__, __FILE__);
		exit(0);
	}

	statsPop(pop, stats);
	statsPop2(pop, stats2);
	scoreSimu(stats, score);
	scoreSimu(stats2, score2);
	nbSitePolyByPop(pop, stats);
	calcCorr(stats);
	calcCorr(stats2);

	fprintf(f, "%ld, %ld, %ld, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %f, %f, %f, %d, %d, %d, %d, %d, %d, %d, %d, %d, %f, %f, %f, %d, %d, %d, %d, %d, %d, %d, %f, %f\n", i, P->nbLocus, P->nb_snp0, stats->snp, stats->fixedCF, stats->fixedSF, stats->shared, stats->ancestralFixedCF, stats->ancestralFixedSF, stats->derivedFixedCF, stats->derivedFixedSF, stats->ancestralFixedBothPop, stats->derivedFixedBothPop, score[1], score[0], score[2], stats2->fixedCF, stats2->fixedSF, stats2->shared, stats2->ancestralFixedCF, stats2->ancestralFixedSF, stats2->derivedFixedCF, stats2->derivedFixedSF, stats2->ancestralFixedBothPop, stats2->derivedFixedBothPop, score2[1], score2[0], score2[2], stats->polyCF, stats->polySF, stats->polyTexas, stats->fixedInitCF, stats->fixedNeoCF, stats->fixedInitSF, stats->fixedNeoSF,stats->corr,stats2->corr);
	free(stats);
	free(score);
	free(stats2);
	free(score2);
	stats = NULL;
	score = NULL;
	stats2 = NULL;
	score2 = NULL;
}

/**
 * \fn void statsPop(Locus* pop, Snp* stats)
 * \brief Counts SNPs by SNPs categories
 * \Param pop locus linked list
 * \Param stats data structure where results are saved
 * High (>.95) and low (< .05) frequency alleles
 * are considered as fixed
*/
void statsPop(Locus* pop, Snp* stats){
	if(pop->cf < 0.05 && pop->texas > 0.95){
		stats->snp++;
		stats->fixedSF++;
	}else if(pop->cf > 0.95 && pop->texas < 0.05){
		stats->snp++;
		stats->fixedCF++;
	}else if(pop->cf < 0.05 && pop->texas < 0.05){
		stats->ancestralFixedBothPop++;
	}else if(pop->cf > 0.95 && pop->texas > 0.95){
		stats->derivedFixedBothPop++;
	}else if(pop->cf >= 0.05 && pop->cf <= 0.95){
		if(pop->texas >= 0.05 && pop->texas <= 0.95){
			stats->snp++;
			stats->shared++;

			stats->x += pop->cf;
			stats->y += pop->texas;
			stats->xy += (pop->cf * pop->texas);
			stats->xsq += pow(pop->cf, 2);
			stats->ysq += pow(pop->texas, 2);
		}else if(pop->texas < 0.05){
			stats->snp++;
			stats->ancestralFixedSF++;
		}else if(pop->texas > 0.95){
			stats->snp++;
			stats->derivedFixedSF++;
		}
	}else if(pop->texas >= 0.05 && pop->texas <= 0.95){
		if(pop->cf < 0.05){
			stats->snp++;
			stats->ancestralFixedCF++;
		}else if(pop->cf > 0.95){
			stats->snp++;
			stats->derivedFixedCF++;
		}
	}

	if(pop->next != NULL){
		statsPop(pop->next, stats);
	}
}

/**
 * \fn void statsPop2(Locus* pop, Snp* stats)
 * \brief Counts SNPs by SNPs categories
 * \Param pop locus linked list
 * \Param stats data structure where results are saved
 * High (>.95) and low (< .05) frequency alleles
 * are NOT considered as fixed contrary to function statsPop
*/

void statsPop2(Locus* pop, Snp* stats){
	if(pop->cf == 0 && pop->texas == 1){
		stats->snp++;
		stats->fixedSF++;
	}else if(pop->cf == 1 && pop->texas == 0){
		stats->snp++;
		stats->fixedCF++;
	}else if(pop->cf == 0 && pop->texas == 0){
		stats->ancestralFixedBothPop++;
	}else if(pop->cf == 1 && pop->texas == 1){
		stats->derivedFixedBothPop++;
	}else if(pop->cf > 0 && pop->cf < 1){
		if(pop->texas > 0 && pop->texas < 1){
			stats->snp++;
			stats->shared++;

			stats->x += pop->cf;
			stats->y += pop->texas;
			//printf("%f ", pop->texas);
			stats->xy += (pop->cf * pop->texas);
			stats->xsq += pow(pop->cf, 2);
			stats->ysq += pow(pop->texas, 2);

		}else if(pop->texas == 0){
			stats->snp++;
			stats->ancestralFixedSF++;
		}else if(pop->texas == 1){
			stats->snp++;
			stats->derivedFixedSF++;
		}
	}else if(pop->texas > 0 && pop->texas < 1){
		if(pop->cf == 0){
			stats->snp++;
			stats->ancestralFixedCF++;
		}else if(pop->cf == 1){
			stats->snp++;
			stats->derivedFixedCF++;
		}
	}

	if(pop->next != NULL){
		statsPop2(pop->next, stats);
	}


}

void nbSitePolyByPop(Locus* pop, Snp* stats){
	if(pop->cf > 0 && pop->cf < 1){
		stats->polyCF++;
	}
	if(pop->sf > 0 && pop->sf < 1){
		stats->polySF++;
	}
	if(pop->texas > 0 && pop->texas < 1){
		stats->polyTexas++;
	}

	if((pop->texas == 0 && pop->cf == 1) || (pop->cf == 1 && pop->texas < 1 && pop->texas > 0)){
		if(pop->gene == 0){
			stats->fixedInitCF++;
		}else{
			stats->fixedNeoCF++;
		}
	}else if((pop->cf == 0 && pop->texas == 1) || (pop->texas == 1 && pop->cf < 1 && pop->cf > 0)){
		if(pop->gene == 0){
			stats->fixedInitSF++;
		}else{
			stats->fixedNeoSF++;
		}
	}else if(pop->texas == 1 && pop->cf == 1){
		if(pop->gene == 0){
			stats->fixedInitSF++;
			stats->fixedInitCF++;
		}else{
			stats->fixedNeoSF++;
			stats->fixedNeoCF++;
		}
	}

	if(pop->next != NULL){
		nbSitePolyByPop(pop->next, stats);
	}
}

/**
 * \fn double calculScore(Snp* stats, int fixedCF, int fixedSF, int shared, int derivedFixedCF, int derivedFixedSF, int ancestralFixedCF, int ancestralFixedSF)
 * \brief Calculate the score of the simulation regarding of the percentage of SNPs in each categories in Illumina assembly
 * \Param stats Number of SNPs by categories in simulation
 * \Param * Percentage of SNPs in the different categories in Illumina
*/
double calculScore(Snp* stats, int fixedCF, int fixedSF, int shared, int derivedFixedCF, int derivedFixedSF, int ancestralFixedCF, int ancestralFixedSF){
	double score = 0;

	score += (double)pow((fixedCF - ceil(((double)stats->fixedCF / (double)stats->snp) * 100)),2) / (double)fixedCF;
	score += (double)pow((fixedSF - ceil(((double)stats->fixedSF/(double)stats->snp)*100)),2) / (double)fixedSF;
	score += (double)pow((shared - ceil(((double)stats->shared/(double)stats->snp)*100)),2) / (double)shared;
	score += (double)pow((derivedFixedCF - ceil(((double)stats->derivedFixedCF/(double)stats->snp)*100)),2) / (double)derivedFixedCF;
	score += (double)pow((derivedFixedSF - ceil(((double)stats->derivedFixedSF/(double)stats->snp)*100)),2) / (double)derivedFixedSF;
	score += (double)pow((ancestralFixedCF - ceil(((double)stats->ancestralFixedCF/(double)stats->snp)*100)),2) / (double)ancestralFixedCF;
	score += (double)pow((ancestralFixedSF - ceil(((double)stats->ancestralFixedSF/(double)stats->snp)*100)),2) / (double)ancestralFixedSF;

	return score;
}

/**
 * \fn void scoreSimu(Snp* stats, double* score)
 * \brief Call calculScore with different percentage of SNPs in each categories (Synonymous, Non-synonymous and Non-coding)
 * \Param stats Number of SNPs by categories in simulation
 * \Param * Percentage of SNPs in the different categories in Illumina
*/

void scoreSimu(Snp* stats, double* score){
	// S
	int fixedCF = 10;
	int fixedSF = 6;
	int shared = 15;
	int derivedFixedCF = 9;
	int derivedFixedSF = 3;
	int ancestralFixedCF = 46;
	int ancestralFixedSF = 13;
	score[0] = calculScore(stats, fixedCF, fixedSF, shared, derivedFixedCF, derivedFixedSF, ancestralFixedCF, ancestralFixedSF);

	// NS
	fixedCF = 11;
	fixedSF = 7;
	shared = 15;
	derivedFixedCF = 8;
	derivedFixedSF = 4;
	ancestralFixedCF = 41;
	ancestralFixedSF = 14;
	score[1] = calculScore(stats, fixedCF, fixedSF, shared, derivedFixedCF, derivedFixedSF, ancestralFixedCF, ancestralFixedSF);

	// NC
	fixedCF = 10;
	fixedSF = 6;
	shared = 15;
	derivedFixedCF = 8;
	derivedFixedSF = 3;
	ancestralFixedCF = 44;
	ancestralFixedSF = 14;
	score[2] = calculScore(stats, fixedCF, fixedSF, shared, derivedFixedCF, derivedFixedSF, ancestralFixedCF, ancestralFixedSF);
}

/**
 * \fn void calcCorr(Snp* stats)
 * \brief Calculate correlation coefficient for shared SNP between CF and SF
 * \Param stats Number of SNPs by categories in simulation
*/
void calcCorr(Snp* stats){
	float num = (stats->shared * stats->xy) - (stats->x * stats->y);
	float denomx = sqrt(stats->shared * stats->xsq - pow(stats->x,2));
	float denomy = sqrt(stats->shared * stats->ysq - pow(stats->y,2));
	//printf("(%d * %f) - (%f * %f) = %f\n",stats->shared, stats->xy, stats->x, stats->y, num);
	stats->corr = num / (denomx * denomy);
}

/*
 *
 */

void exportLocusPop(Locus* pop, gsl_rng* r){
	FILE* out = fopen("statsLocus.txt", "w");
	fprintf(out, "\"%%der CF\", \"%%der SF\", \"CF1 all 1\", \"CF1 all 2\", \"CF2 all 1\", \"CF2 all 2\", \"SF1 all 1\", \"SF1 all 2\", \"SF2 all 1\", \"SF2 all 2\"\n");

	printFileExportLocus(out, pop, r);
	fclose(out);
}

void printFileExportLocus(FILE* f, Locus* pop, gsl_rng* r){
	int all[8] = {0,0,0,0,0,0,0,0}; // all 1 cf, all 2 cf, all 1 sf, all 1 sf

	int i,j;
	for(i = 0; i < 4; i++){
		all[i] =  (gsl_rng_uniform(r) <= pop->cf) ? 1 : 0;
	}

	for(j=i; j < 8; j++){
		all[j] =  (gsl_rng_uniform(r) <= pop->texas) ? 1 : 0;
	}

	fprintf(f, "%f, %f", pop->cf, pop->texas);
	for(i = 0; i < 8; i++){
		fprintf(f, ",%d", all[i]);
	}

	fprintf(f, "\n");

	if(pop->next != NULL){
		printFileExportLocus(f, pop->next, r);
	}
}
