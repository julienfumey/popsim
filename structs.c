/**
 	* AstyanaxPopSim 
	* \file structs.c
	* \author Julien Fumey <julien.fumey@legs.cnrs-gif.fr>, LEGS - CNRS UPR 9034
	* \brief Functions for manipulating data structures
	* \copyright (C) 2014
	* \version 2.0
*/

#include "structs.h"

/**
 * \fn void freeLocus(Locus* pop)
 * \brief Delete all locus from pop
 * \Param pop the linked list to be delete
*/
void freeLocus(Locus* pop){
	Locus* temp;
	while(pop != NULL){
		temp = pop;
		pop = pop->next;
		free(temp);
		temp = NULL;
	}
}

/**
 * \fn Locus* copyList(Locus* pop)
 * \brief Make a copy of the linked list pop
 * \Param pop the linked list to be copy
*/
Locus* copyList(Locus* pop){
	Locus* locus = calloc(1, sizeof(Locus));
	if(locus == NULL){
		fprintf(stderr, "error on %s at line %d", __FILE__, __LINE__);
		exit(0);
	}

	locus->cf = pop->cf;
	locus->sf = pop->sf;
	locus->texas = pop->texas;
	locus->gene = pop->gene;
	
	if(pop->next != NULL){
		locus->next = copyList(pop->next);
	}else{
		locus->next = NULL;
	}	
	
return locus;
}
