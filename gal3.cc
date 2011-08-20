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
CChromosomeFactoryBase &CF = *new CChromosomeFactory<unsigned long>;

COpMutate opmut(0.001);
COpCrossover opcross(0.70);
COpSelectStochastic selector;

int main() {

  CEngine ga(80,1,1.6); // popsize = 80, chromlen = 1, scalefactor = 1.6
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
  
  CChromosome<unsigned long> &chr = 
    (CChromosome<unsigned long> &)pool.GetBest();
  std::cout << "X = " << (double)chr[0]/(double)ULONG_MAX << std::endl;
  return 0;
}

double CEvaluator::EvaluateFitness(CChromosomeBase &chrbase)
{
  CChromosome<unsigned long> &chr = (CChromosome<unsigned long> &)chrbase;
  double total = 0;

  unsigned long val = chr[0];

  if (val & (1 << 1)) total += 1;
  if (val & (1 << 2)) total += 2;
  if (val & (1 << 5)) total += 3;
  if (val & (1 << 8)) total += 4;
  if (val & (1 << 7)) total -= 3;
  if (val & (1 << 30)) total += 15;
  if (val & (1 << 9)) total /= 2;
  if (total < 0) total = 0;
  return total;
}


