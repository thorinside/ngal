/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/
#include <iostream>
#include <math.h>
#include "pool.h"

/////////////////////////////////////////////////////////////////
// CPool class implementation

CPool::CPool()
{
  nSize = 0;
  apchrData = 0;
  nChromlen = 0;
}

CPool::CPool(int poolsize, int chromlen, double scalefactor)
{
  nSize = poolsize;
  nChromlen = chromlen;

  fScaleFitness = 0;
  if (scalefactor > 1.0) {
    fScaleFitness = 1;
  }
  dScaleFactor = scalefactor;

  apchrData = new CChromosomeBase *[nSize];
  
  for (int i = 0; i < nSize; i++) {
    apchrData[i] = CF.Create(nChromlen);
    InitChromosome(*apchrData[i]);
  }
  CalcPoolFitness();
}

CPool::CPool(const CPool& that)
{
  nSize = 0;
  apchrData = 0;
  iBest = 0;
  nChromlen = 0;
  fScaleFitness = 0;
  dScaleFactor = 1.0;
  *this = that;
}

CPool::~CPool()
{
  if (nSize && apchrData) {
    for (int i = 0; i < nSize; i++) {
      delete apchrData[i];
    }
    delete [] apchrData;
  }
}

void CPool::InitChromosome(CChromosomeBase &chr)
{
  // This function is meant to be overridden to
  // provide a way to customize initialization of
  // chromosomes. The default is just to randomize
  // the bits.

  chr.Randomize();
}

void CPool::CalcPoolStatistics()
{
  double max, avg, min, sumfitness;

  sumfitness = apchrData[0]->GetFitness();
  min = apchrData[0]->GetFitness();
  iBest = 0;
  max = apchrData[0]->GetFitness();

  for (int i = 1; i < nSize; i++)
    {
      double fitness = apchrData[i]->GetFitness();
      sumfitness = sumfitness + fitness;
      if (fitness > max) {
	max = fitness;
	iBest = i;
      }
      if (fitness < min) min = fitness;
    }
  avg = sumfitness/(double)nSize;
  
  dMaxFitness = max;
  dMinFitness = min;
  dMeanFitness = avg;
  dSumFitness = sumfitness;

  // Next we figure out the variance
  sumfitness = 0.0;
  for (int i = 0; i < nSize; i++) {
    double fitness = apchrData[i]->GetFitness();
    sumfitness += (fitness - avg)*(fitness - avg);
  }
  dVarFitness = sumfitness/nSize;
}

CChromosomeBase &CPool::GetBest()
{
  return *apchrData[iBest];
}

double CPool::GetSumFitness() const
{
  return dSumFitness;
}

double CPool::GetMeanFitness() const
{
  return dMeanFitness;
}

double CPool::GetMaxFitness() const
{
  return dMaxFitness;
}

double CPool::GetMinFitness() const
{
  return dMinFitness;
}

double CPool::GetVariance() const
{
  return dVarFitness;
}

double CPool::GetStdDev() const
{
  return sqrt(dVarFitness);
}

int CPool::GetPoolSize() const
{
  return nSize;
}

int CPool::GetChromLen() const
{
  return nChromlen;
}

CChromosomeBase &CPool::operator[](int index) const
{
  // We'd very likely throw an exception
  // if the index is out of range. Here we
  // will check and exit if we're out of
  // range.
  
  if ((index >= 0) && (index < nSize))
    return *apchrData[index];
  else {
    std::cerr << "CChromosomeBase &CPool::operator[](int index) " << index << std::endl;
    exit(1);
  }
  return *apchrData[0];
}

CPool &CPool::operator=(const CPool &that)
{
  if ((nSize != 0) && (apchrData != 0)) {
    for (int i = 0; i < nSize; i++) 
      delete apchrData[i];
    delete [] apchrData;
  }
  
  nSize = that.nSize;
  nChromlen = that.nChromlen;
  dSumFitness = that.dSumFitness;
  dMeanFitness = that.dMeanFitness;
  dMaxFitness = that.dMaxFitness;
  dMinFitness = that.dMinFitness;
  dVarFitness = that.dVarFitness;
  iBest = that.iBest;

  fScaleFitness = that.fScaleFitness;
  dScaleFactor = that.dScaleFactor;

  // Just for completeness we will even
  // copy over the scaling information even
  // though it may not be needed
  
  dUmin = that.dUmin;
  dUmax = that.dUmax;
  dUavg = that.dUavg;
  a = that.a;
  b = that.b;

  // Create new cromosomes and copy the important
  // data over chromosome by chromosome

  apchrData = new CChromosomeBase *[nSize];
  for (int i = 0; i < nSize; i++) {
    apchrData[i] = CF.Create(nChromlen);
    *apchrData[i] = that[i];
    apchrData[i]->SetObjective(that[i].GetObjective());
    apchrData[i]->SetFitness(that[i].GetFitness());
  }
  return *this;
}

void CPool::ScaleFitness()
{
  if (fScaleFitness == 1) {
    Prescale();
    for (int i = 0; i < nSize; i++) {
      apchrData[i]->SetFitness(Scale(apchrData[i]->GetObjective()));
    }
  }
  else {
    for (int i = 0; i < nSize; i++) {
      apchrData[i]->SetFitness(apchrData[i]->GetObjective());
    }
  }
}

double CPool::Scale(double objective)
{
  double val = (a * objective + b);

  // val may be a very small negative number at
  // this point. If so, return the floating point
  // absolute value of the number. This should
  // ensure the fitness function is always positive
  // without causing too much harm to the algorithm

  if (val < 0.0) val = fabs(val);
  return val;
}

void CPool::Prescale()
{
  double delta;

  if (dUmin > ((dScaleFactor*dUavg - dUmax) / 
		     (dScaleFactor - 1.0))) { // Non-negative test
    delta = dUmax - dUavg;

    // Here delta may be zero, so to avoid
    // floating point exceptions we'll scale
    // linearly for this part of the algorithm

    if (delta == 0.0) { a = 1; b = 0; }
    else {
      a = (dScaleFactor - 1.0) * dUavg / delta;
      b = dUavg * (dUmax - dScaleFactor*dUavg) / delta;
    }
  }
  else {
    delta = dUavg - dUmin;
    if (delta == 0.0) { a = 1; b = 0; }
    else {
      a = dUavg / delta;
      b = -dUmin * dUavg / delta;
    }
  }
}

void CPool::CalcPoolFitness()
{
  // Calculate Objective Function values and
  // determine min, max, and avg objective function
  // values for scaling function.
  
  // initialize variables to initial
  // values
  double dTotal = 0.0;
  double dObj = EV.EvaluateFitness(*apchrData[0]);
  dUmin = dObj;
  dUmax = dObj;

  for (int i = 0; i < nSize; i++) {
    dObj = EV.EvaluateFitness(*apchrData[i]);
    apchrData[i]->SetObjective(dObj);
    if (dObj < dUmin) dUmin = dObj;
    if (dObj > dUmax) dUmax = dObj;
    dTotal += dObj;
  }
  dUavg = (dTotal/(double)nSize);

  // Now we run the scaling function to scale
  // things properly.
  ScaleFitness();

  // Here is the only place that pool statistics
  // should need to be calculated
  CalcPoolStatistics();
}
