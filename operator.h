/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#ifndef OPERATOR_H_
#define OPERATOR_H_

#include "chrom.h"
#include "pool.h"

class COperator 
{	
public:		
  COperator();
  
  long GetUsages();
  void SetProbability(double prob);
  double GetProbability();

protected:
  void IncUsages();
  
private:
  long nUsages;
  double dProb;
};

class COpSelect : public COperator
{
public:
  COpSelect() {}
  virtual void Preselect(CPool &) = 0;
  virtual CChromosomeBase &Select(CPool &) = 0;
};

class COpSelectRoulette : public COpSelect
{
public:
  COpSelectRoulette() {}
  virtual void Preselect(CPool &) {}
  virtual CChromosomeBase &Select(CPool &);
};

class COpSelectStochastic : public COpSelect
{
public:
  COpSelectStochastic() : selected(0) {}
  virtual ~COpSelectStochastic();
  virtual void Preselect(CPool &);
  virtual CChromosomeBase &Select(CPool &);
private:
  int nremain;
  int *selected;
};

class COpMutate : public COperator 
{
public:
  COpMutate();
  COpMutate(double prob);

  virtual void Mutate(CChromosomeBase &);
};

class COpCrossover : public COperator 
{
public:
  COpCrossover();
  COpCrossover(double prob);

  virtual void Crossover(CChromosomeBase &, CChromosomeBase &);
};
#endif // OPERATOR_H_
