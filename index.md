### Welcome to PopSim.

### What is PopSim?
Popsim is a Wright-Fisher population genetics simulator.

It was first designed to study *Astyanax mexicanus* cavefish evolution (see How to cite PopSim section).

### How to install PopSim
You need to have the GNU Scientific Library (GSL) installed on your computer to use PopSim.

To install PopSim, just compile it:

`gcc *.c -o PopSim -lgsl -lgcblas -lm`

### How to launch PopSim
To launch PopSim, you need to create a parameters file. All parameters must be filled :


Lines can be commented using a #
for each parameter, just write parameter name followed by an = signs and the parameter value
ie. param = value

#List of parameters:
* **nb_locus**: number of locus to be generated (if initial state generated without burning)
* **nb_gen**: number of generation
* **init_gen**: number of burning generation
* **elabra_gen**: number of generation in the "El Abra" population
* **bn_gen**: number of generations during bottleneck
* **labCF_gen**: number of lab generation for "CF" pop
* **labSF_gen**: number of lab generations for "SF" pop
* **ancestral_popsize**: Effective size of the ancestral population
* **texas_popsize**: Effective size of the "Texas" population
* **elabra_popsize**: Effective size of the "El Abra" population
* **bn_popsize**: Effective population during bottleneck
* **cf_popsize**: Effective size of the "CF" population
* **labCF_popsize**: Effective size of the "CF" population in lab
* **labSF_popsize**: Effective size of the "Texas" and "El Abra" populations in lab
* **mutaRate**: mutation rate
* **migra**: Should migration occur between population ? 0: no; 1: yes
* **STmr**: migration rate from "El Abra" to "Texas" populations
* **migrantfish_percST**: frequency of "El Abra" population migrating to "Texas"
* **TSmr**: : migration rate from "Texas" to "El Abra" populations
* **migrantfish_percTS**: frequency of "Texas" population migrating to "El Abra"
* **SCmr**: migration rate from "El Abra" to "CF" populations
* **migrantfish_percSC**: frequency of "El Abra" population migrating to "CF"
* **CSmr**: migration rate from "CF" to "El Abra" populations
* **migrantfish_percCS**: frequency of "CF" population migrating to "El Abra"
* **CTmr**: migration rate from "CF" to "Texas" populations
* **migrantfish_percCT**: frequency of "CF" population migrating to "Texas"
* **TCmr**: migration rate from "Texas" to "CF" populations
* **migrantfish_percTC**: frequency of "Texas" population migrating to "CF"
* **resultfile**: name of output file before lab evolution
* **resultfileLab**: name of output after lab evolution
* **generatePopInit**: How should the ancestral population be generated ? 0: burning 1: Sampling in theoretical distribution (recommanded)
* **exportLocusPop**: 1: will export derived alleles frequencies at each locus at generation defined by *exportLocusPopGene*; 0: won't export
* **exportLocusPopGene**: At which generation should derived alleles frequencies be exported?
* **subGenSF**: Number of subgenerations in "El Abra" population
* **subGenCF**: Number of subgenerations in "CF" population
* **subGenTexas**: Number of subgenerations in "Texas" population

# How to cite PopSim?
If you use PopSim, please cite this paper :

Evidence of Late Pleistocene origin of Astyanax mexicanus cavefish
Julien Fumey, Hélène Hinaux, Céline Noirot, Sylvie Rétaux, Didier Casane
bioRxiv 094748; doi: https://doi.org/10.1101/094748

