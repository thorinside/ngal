
/*

	pRandom.h
		- header file for C++ class pRandom


	class prandom
		- provides a Genetic Algorithm oriented interface to the 
		  random number generator of GNU C++


	L. Weaver		May 16, 1995

	This version:	@(#)pRandom.h	1.5	09 Jun 1995
*/

#ifndef pRandom_H
#define pRandom_H

#include "ACG.h"
#include <time.h>

class pRandom {
public:


	// FUNCTIONS TO ACCESS RANDOM NUMBERS
		// coin toss; 0.5 = unbiased
	int True(float p=0.5);

		// returns an integer in [lower_bound, upper_bound]
	int uniform_int(int lower_bound, int upper_bound);
		// returns a float in [lower_bound, upper_bound)
	float uniform_float(float lower_bound, float upper_bound);

		// exponential distribution : an efficient implementation
	int exp_int(int base, float p);

		// generate a specified number of random bits at an address
	int bytes_at(void *addr, int num_bytes);

	// MUTATION CLOCK FUNCTIONS
		//	- an efficient method of checking for mutation
	float set_mutation_rate(float p);	// sets the rate of mutation
	int mutation();	// returns TRUE with prob. p as set in set_mutation_rate




	// SEED RELATED FUNCTIONS
		// return the last seed used (works for autoseed() & seed())
	int last_seed(){return lastseed;};

		// change generators to one seeded with new_seed
	int seed(int new_seed);

		// change to a generator seeded automatically by the system
	int autoseed(); 



	// HOUSE-KEEPING STUFF
	pRandom();
	~pRandom();
	

private:

	ACG *generator;		// the generator being used
	int lastseed;

	long unsigned int mutation_clock;
	float mutation_prob;
	
};


extern pRandom Random;

#endif
