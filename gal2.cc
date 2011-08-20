/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include "ga.h"
#include "pRandom.h"

const int NUMUSERS = 6;

long nDisks = 4;
int diskspace[4] = {2000, 2100, 2100, 50};
int userspace[NUMUSERS];

double CEvaluator::EvaluateFitness(CChromosomeBase &chrbase)
{
  CChromosome<unsigned long> &chr = (CChromosome<unsigned long> &)chrbase;
  double total = 0;
  double totalused = 0;
  int used[4] = {0,0,0,0};

  for (int i = 0; i < NUMUSERS; i++) {
    unsigned long val = chr[i];
    int disk = 0;
    if (val & 1) disk += 1;
    if (val & 2) disk += 2;

    //    int disk = (val & 1) + (val & 2)*2;
    //    unsigned long partsize = ULONG_MAX / (nDisks-1);
    //    unsigned long disk = chr[i] % nDisks;
    //    cerr << "disk = " << disk << endl;

    if ((diskspace[disk] - used[disk] - userspace[i]) >= 0) {
      used[disk] += userspace[i];
    }
  }
  for (int i = 0; i < nDisks; i++) {
    totalused += used[i];
    total += diskspace[i];
  }
  double retval = totalused;
  return retval;
}

CEvaluator &EV = *new CEvaluator;
CChromosomeFactoryBase &CF = *new CChromosomeFactory<unsigned long>;

COpMutate op(0.001);
COpCrossover op2(0.60);
COpSelectRoulette selector;

int main() {

  // Initialize user array
  long total = 0;
  for (int user = 0; user < NUMUSERS; user++) {
    userspace[user] = Random.uniform_int(2,200);
    total += userspace[user];
  }

  std::cout << "Total = " << total << std::endl;

  CEngine ga(50,NUMUSERS);
  int converged = 0;
  
  ga.SetMutateOp(op);
  ga.SetCrossoverOp(op2);
  ga.SetSelectionOp(selector);
  
  CPool &pool = ga.GetPool();

  std::cout << "[" << ga.nGenerations << "]" <<
    " Stats: max=" << pool.GetMaxFitness() <<
      " min=" << pool.GetMinFitness() <<
	" avg=" << pool.GetMeanFitness() <<
	  " sum="  << pool.GetSumFitness() <<
	    " nmutation="  << op.GetUsages() <<
	      " ncross=" << op2.GetUsages() << std::endl;
  
  for (int i = 0; i < 2500; i++)
    {
      ga.Generation();
      
      if ((i / 5) && !(i % 5))
	std::cout << "[" << ga.nGenerations << "]" <<
	  " Stats: max=" << pool.GetMaxFitness() <<
	    " min=" << pool.GetMinFitness() <<
	      " avg=" << pool.GetMeanFitness() <<
		" sum="  << pool.GetSumFitness() <<
		  " nmutation="  << op.GetUsages() <<
		    " ncross=" << op2.GetUsages() << std::endl;

      if (pool.GetVariance() == 0.0) {
//	      converged++;
      }
      else 
	converged = 0; 
      if (converged == 1) {
        break;
      }
    }
	std::cout << "Population converged: " << std::endl;
	std::cout << "  pool fitness max  = " << pool.GetMaxFitness() <<
          "  variance = " << pool.GetVariance() << std::endl <<
	    "  generations = " << ga.nGenerations << std::endl <<
	      "  children = " << ga.nChildren << std::endl;

	CChromosome<unsigned long> &chr = 
	  (CChromosome<unsigned long> &)pool.GetBest();
	std::cout << "REALFITNESS = " << (double)chr.GetObjective() << std::endl;
  return 0;
}
