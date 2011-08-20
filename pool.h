/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#ifndef POOL_H_
#define POOL_H_

#include "chrom.h"
#include "eval.h"

// There must be a Chromosome Factory that produces the correct
// chromosome type for this application. A line such as the
// following must exist somewhere in the main program code:
//
// CChromosomeFactoryBase &CF = *new CChromosomeFactory<short>;

extern CChromosomeFactoryBase &CF; // external chromosome factory
extern CEvaluator &EV;             // external fitness evaluator

class CPool
{
public:
  CPool();
  CPool(int poolsize, int chromlen, double scalefitness = 1.0);
  CPool(const CPool &);
  virtual ~CPool();

  // Virtual Methods

  virtual void InitChromosome(CChromosomeBase &);

  // Public Methods

  void CalcPoolFitness();

  CChromosomeBase &GetBest();

  double GetSumFitness() const;
  double GetMeanFitness() const;
  double GetMaxFitness() const;
  double GetMinFitness() const;
  double GetVariance() const;
  double GetStdDev() const;
  int GetPoolSize() const;
  int GetChromLen() const;

  CChromosomeBase &operator[](int index) const;
  CPool &operator=(const CPool &);

private:
  int nSize;
  int nChromlen;
  double dSumFitness;  // Pool Fitness total (sum)
  double dMeanFitness; // Pool Fitness mean
  double dMaxFitness;  // Pool Fitness maximum
  double dMinFitness;  // Pool Fitness minimum
  double dVarFitness;  // Pool Fitness variance

  CChromosomeBase **apchrData; // Array of chromosomes
  int iBest;    // index of the pool's best chromosome

  // fitness scaling flag
  int fScaleFitness;
  double dScaleFactor;

  // Objective min, max and avg for scaling
  double dUmin;
  double dUmax;
  double dUavg;

  // Scaling coefficients
  double a, b;

  // Pre-scaling function, called by Scale should
  // not be called on its own.
  double Scale(double objective);
  void Prescale();
  void ScaleFitness();
  void CalcPoolStatistics();
};

#endif // POOL_H_

