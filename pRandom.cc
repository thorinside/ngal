
/*

	pRandom.c
		- implementation file for C++ class pRandom


	class pRandom
		- provides a Genetic Algorithm oriented interface to the 
		  random number generator of GNU C++


	L. Weaver		May 16, 1995

	This version:	@(#)pRandom.c	1.6	10 Jun 1995
*/

#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>
#include "ACG.h"
#include <iostream>

#include "pRandom.h"

#define ACG_RANDMAX	4294967295


pRandom Random;

pRandom::pRandom()
{
	generator=NULL;
	mutation_clock = 0;
	mutation_prob = 0.0;

	autoseed();
}

pRandom::~pRandom()
{
	delete generator;
}


int pRandom::True(float p)
// returns TRUE with probability p
{
	return (generator->asLong()<=(long unsigned int)
						(p*(double)ACG_RANDMAX));
}


int pRandom::uniform_int(int lower_bound, int upper_bound)
// returns each integer in [lower_bound,upper_bound] with equal probability
{
	register double range = (double) (upper_bound-lower_bound+1);
	return lower_bound +
		(int)(
			range*
				((double)(generator->asLong())/
				(1.0+(double)ACG_RANDMAX))
		);
}


float pRandom::uniform_float(float lower_bound, float upper_bound)
// returns all values in [lower_bound,upper_bound) with equal probability
{
	register double range = (double) (upper_bound-lower_bound);
	return lower_bound +
		(
			range*
				((double)(generator->asLong())
				/(1.0+(double)ACG_RANDMAX))
		);
}


int pRandom::exp_int(int base, float p)
// returns an integer x >= base, with prob (1-p)^(x-base)*p
{
	// this method is probably far from efficient, but is a working
	// placeholder
//	while(1) {
//		if (True(p))
//			return base;
//		else
//			base++;
//	}


	// a far faster implementation than the previous method
	double p1 = (double)p;
	double total = p1;
	double randval = (double) uniform_float(0.0,1.0);

	while(total<randval){
		base++;
		total += (1.0-total)*p1;
	}

	return base;

}


int pRandom::bytes_at(void *addr, int num_bytes)
// generate num_bytes random bytes, placing them at addr
{
	int i;
	long num;
	char *base = (char *)addr;

	for(i=0;i<num_bytes;i++) {
		num = generator->asLong();
		(*base) = *(char *)&num;
		base++;
	}

	return num_bytes;
}


float pRandom::set_mutation_rate(float p)
// set the probability of mutation
//	must be done BEFORE mutation() is first called
{
	assert(p>=0.0);

	mutation_clock = exp_int(1,p);
	return (mutation_prob=p);
}


int pRandom::mutation()
// returns TRUE (1) with probability mutation_rate, as set by 
// set_mutation_rate() which must be called before this function is first called
//	works by setting mutation_clock (via exp_int()) to a randomly determined
//	number of ticks(calls to mutation()) which should occur before the next
//	TRUE is returned. When this ticks expire, a new number of ticks is
//	determined and TRUE returned
//		- this is much faster than generating a random number each time
//		  to check if mutation should occur
{
	if ((--mutation_clock)==0) {
		// time for a mutation
			// reset clock for the next one
		mutation_clock = (long unsigned int) exp_int(1,mutation_prob);
		return 1;

	} else {
		// not yet time for a mutation
		return 0;
	}
}

int pRandom::seed(int new_seed)
// change generators to one seeded with new_seed
{
	lastseed = new_seed;

	if (generator!=NULL)
		delete generator;

	generator = new ACG(lastseed);

	return lastseed;
}

int pRandom::autoseed()
// change to a generator seeded automatically by the system
{
	lastseed = (unsigned int)(time(0));

	if (generator!=NULL)
		delete generator;

	generator = new ACG(lastseed);

	return lastseed;
}

