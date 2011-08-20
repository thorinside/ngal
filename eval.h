/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#ifndef EVAL_H_
#define EVAL_H_

class CEvaluator {	
public:		
  CEvaluator() {}
  virtual double EvaluateFitness(CChromosomeBase &);
};

#endif // EVAL_H_
