/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "pRandom.h"
#include "operator.h"

/////////////////////////////////////////////////////////////////

COperator::COperator()
{
  nUsages = 0;
}

long COperator::GetUsages() 
{
  return nUsages;
}

void COperator::IncUsages() 
{
  nUsages++;
}

void COperator::SetProbability(double prob)
{
  dProb = prob;
}

double COperator::GetProbability()
{
  return dProb;
}

/////////////////////////////////////////////////////////////////

CChromosomeBase &COpSelectRoulette::Select(CPool &pool)
{
  double randnum, partsum;
  int j = 0;
  int popsize = pool.GetPoolSize();

  randnum = Random.uniform_float(0.0,(float)pool.GetSumFitness());
  partsum = 0.0;
  
  for (j = 0; ((j < popsize) && (partsum <= randnum)); j++) {
    partsum += pool[j].GetFitness();
  }
  if (j == popsize) return pool[j-1];
  else
    return pool[j];
}

/////////////////////////////////////////////////////////////////

COpSelectStochastic::~COpSelectStochastic()
{
  if (selected) delete [] selected;
}

void COpSelectStochastic::Preselect(CPool &pool)
{
  double expected;
  double jassign;
  double *fraction;
  int j, k;
  
  int poolsize = pool.GetPoolSize();
  double avg = pool.GetMeanFitness();

  fraction = new double[poolsize];

  if (selected) delete [] selected;
  selected = new int[poolsize];
  nremain = poolsize;

  j = 0;
  k = 0;
  do {
    expected = pool[j].GetFitness() / avg;
    fraction[j] = modf(expected, &jassign);
    while (jassign > 0.0) {
      jassign--;
      selected[k++] = j;
    }
    j++;
  } while (j != poolsize);
  j = 0;
  int infinite = 0;
  while (k < poolsize) { // assign fractional parts
    infinite++;
    if (infinite > 20000) 
      std::cerr << "Got here... k = " << k << " avg = " << avg
	<< " fraction[j] = " << fraction[j] << std::endl;
    if (j >= poolsize) j = 0;
    if (fraction[j] > 0.0) {
      if (Random.True(fraction[j])) {
	selected[k++] = j;
	fraction[j] = fraction[j] - 1.0;
      }
    }
    j++;
  }
}

CChromosomeBase &COpSelectStochastic::Select(CPool &pool)
{
  int jpick;
  int index;

  jpick = Random.uniform_int(0,nremain-1);
  index = selected[jpick];
  selected[jpick] = selected[--nremain];
  return pool[index];
}

/////////////////////////////////////////////////////////////////

COpMutate::COpMutate() : COperator()
{
  SetProbability(0.001); // Default probability
}

COpMutate::COpMutate(double prob) : COperator()
{
  SetProbability(prob);
}

void COpMutate::Mutate(CChromosomeBase &chr)
{
  for (int i = 0; i < chr.GetNumBits(); i++) {
    if (Random.True(GetProbability())) {
      IncUsages();
      if (chr.GetBit(i)) 
	chr.SetBit(i,0);
      else 
	chr.SetBit(i,1);
    }
  }    
}

/////////////////////////////////////////////////////////////////

COpCrossover::COpCrossover() : COperator()
{
  SetProbability(0.6); // Default probability
}

COpCrossover::COpCrossover(double prob) : COperator()
{
  SetProbability(prob);
}

void COpCrossover::Crossover(CChromosomeBase &chr1, CChromosomeBase &chr2)
{
  if (Random.True(GetProbability())) {
    IncUsages();
    int len = chr1.GetNumBits();
    int site = Random.uniform_int(0,len-1);

    // Do the crossover
    for (int i = site; i < len; i++) {
      int bit = chr1.GetBit(i);
      chr1.SetBit(i, chr2.GetBit(i));
      chr2.SetBit(i, bit);
    }
  }    
}
