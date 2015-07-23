/**
 * AstyanaxPopSim
 *
 * \file AstyanaxPopSim.c
 * \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
 * \version 2.0
 *
*/

#include "AstyanaxPopSim.h"

int main(int argc, char *argv[]){
	if(argc != 2){
		fprintf(stderr, "Usage : ./simu <parameters_file>\n");
		exit(1);
	}

	//srandom(time(NULL) + getpid());
	gsl_rng* random_generator = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(random_generator, time(0)*getpid());

	clock_t begin = clock();

	Param *parameters = parseParametersFile(argv[1]);
	Locus* pop = initializeSimulation(parameters, random_generator);
	

	pop = evolve(pop, parameters, random_generator);
	
	gsl_rng_free(random_generator);

	FILE* f = fopen("done", "w");
	fprintf(f, "done at %s %s\nexecuted in %f", __DATE__, __TIME__, ((double)(clock() - begin)/CLOCKS_PER_SEC));
	fclose(f);
	return 1;
}
