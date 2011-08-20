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

void printstatline(CEngine &ga, CPool &pool, COperator &opmut, 
		   COperator &opcross)
{
  std::cout << "[" << ga.nGenerations << "]" <<
    " Stats: max=" << pool.GetMaxFitness() <<
    " min=" << pool.GetMinFitness() <<
    " avg=" << pool.GetMeanFitness() <<
    " sum="  << pool.GetSumFitness() <<
    " nmutation="  << opmut.GetUsages() <<
    " ncross=" << opcross.GetUsages() << std::endl;
}

CEvaluator &EV = *new CEvaluator;
CChromosomeFactoryBase &CF = *new CChromosomeFactory<unsigned short>;

COpMutate opmut(0.01);
COpCrossover opcross(0.90);
COpSelectStochastic selector;

int main() {

  CEngine ga(30,5,1.1); // popsize = 30, chromlen = 5, scalefactor = 1.1
  int converged = 0;
  
  ga.SetMutateOp(opmut);
  ga.SetCrossoverOp(opcross);
  ga.SetSelectionOp(selector);
  
  CPool &pool = ga.GetPool();
  printstatline(ga,pool,opmut,opcross);
  
  for (int i = 0; i < 150; i++)
    {
      ga.Generation();
      
      if ((i / 5) && !(i % 5)) {
	printstatline(ga,pool,opmut,opcross);
      }

      if (pool.GetVariance() == 0.0) {
	converged++;
      }
      else 
	converged = 0; 
      if (converged == 5) {
	std::cout << "Population converged: " << std::endl;
	break;
      }
    }

  std::cout << "  pool fitness max  = " << pool.GetMaxFitness() <<
    "  variance = " << pool.GetVariance() << std::endl <<
    "  generations = " << ga.nGenerations << std::endl <<
    "  children = " << ga.nChildren << std::endl;
  
  return 0;
}

double CEvaluator::EvaluateFitness(CChromosomeBase &chrbase)
{
  CChromosome<unsigned short> &chr = (CChromosome<unsigned short> &)chrbase;
  double fitness = 0;
  
  for (int i = 0; i < chr.GetLength(); i++) {
    fitness += chr[i];
  }
  return fabs(fitness);
}







