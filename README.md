# two-islands

"Diversity patterns and speciation processes ina two-island system with continuous migration", by
Debora Princepe, Simone Czarnobai, Thiago M. Pradella, Rodrigo A. Caetano,Flavia M. D. Marquitti, Marcus A.M. de Aguiar, and Sabrina B. L. Araujo

We provide two codes, both written in FORTRAN:

1) Program_FiniteB_2islands.f90 
	This program considers a finite genome, and the user can set the parameters in the code. 
	The outputs are TemporalB.txt, NetworkB.txt, Ind_informatonB.txt. 
	1.1) Temporal.txt has 7 columns, the first 3 refers to  the parameters: population size (pop), genome size(B) and   migration rate(mig). The remaining columns refer to the time (generation), species code (esp), its respective abundance in island 1 (abund1) and in island2(abund2).
	1.2) Network.txt can be used to plot the network of mating compatibility. It has four columns:  the time (generation), codes of a pair of compatible individuals (ind1 and ind2), and the percentage of dissimilarity between the two individuals H/B*100 (H/B).
	1.3)Ind_information.txt has five columns: it gives at each time (generation), the code of each individual (ind), its species (esp) under a global referential, the island that it is (island), and if it is a migrant (migrate_yes/no) 

	For additional information, please, contact Sabrina Araujo (araujosbl@fisica.ufpr.br). 


2) infinite_model.f90 
	This program considers the infinite genome model, and the user can set the other parameters in the code. The outputs are:
	number.dat - time, sp.total, sp.island1, sp.island2
	numbertot.dat -  3 lines for each time step:
	                 time, abundances of species in island1
	                 time, abundances of species in island2
	                 time, abundances of all species in the system
