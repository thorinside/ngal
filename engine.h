/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#ifndef ENGINE_H_
#define ENGINE_H_

#include "pool.h"
#include "chrom.h"
#include "operator.h"

extern CEvaluator &EV;             // external fitness evaluator

class CEngine {	
public:		
  CEngine(int poolsize,int chromlen,double scalefactor = 1.0);	
  virtual ~CEngine();

  void Generation();

  CPool &GetPool();
  COperator &GetSelectionOp();
  COperator &GetMutateOp();
  COperator &GetCrossoverOp();
  void SetSelectionOp(COpSelect &);
  void SetMutateOp(COpMutate &);
  void SetCrossoverOp(COpCrossover &);

  int nGenerations;
  int nChildren;

private:
  COpSelect *opSel;
  COpMutate *opMutate;
  COpCrossover *opCross;

  CPool *pool;
};

#endif // ENGINE_H_
