/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#include <iostream>
#include "pRandom.h"
#include "engine.h"

CEngine::CEngine (int poolsize, int chromlen, double scalefactor)
{
  nGenerations = 0;
  nChildren = 0;

  pool = new CPool(poolsize, chromlen, scalefactor);
  
  opSel = 0;
  opMutate = 0;
  opCross = 0;
}
	
CEngine::~CEngine()
{
  delete pool;
}

void CEngine::Generation()
{
  // Check if we've got all our operators

  if (opSel && opMutate && opCross) {
    CPool newpool(*pool); // Create a new pool
    opSel->Preselect(*pool); // Do preselection
    for (int i = 0; i < pool->GetPoolSize(); i += 2) {
      CChromosomeBase &child1 = newpool[i];
      CChromosomeBase &child2 = newpool[i+1];

      child1 = opSel->Select(*pool);
      child2 = opSel->Select(*pool);

      // Mutate and crossover
      opCross->Crossover(child1, child2);
      opMutate->Mutate(child1);
      opMutate->Mutate(child2);

      nChildren += 2;
    }
    
    // Be an elitist
    int loc = Random.uniform_int(0,pool->GetPoolSize()-1);
    newpool[loc] = pool->GetBest();

    newpool.CalcPoolFitness(); // Do fitness calculations

    *pool = newpool; // Copy the new pool back
    nGenerations++;
  }
  else {
    std::cerr << "void CEngine::Generation() not enough operators" << std::endl;
  }
}

CPool &CEngine::GetPool()
{
  return *pool;
}

COperator &CEngine::GetSelectionOp()
{
  return *opSel;
}

COperator &CEngine::GetMutateOp()
{
  return *opMutate;
}

COperator &CEngine::GetCrossoverOp()
{
  return *opCross;
}

void CEngine::SetSelectionOp(COpSelect &select)
{
  opSel = &select;
}

void CEngine::SetMutateOp(COpMutate &mutate)
{
  opMutate = &mutate;
}

void CEngine::SetCrossoverOp(COpCrossover &cross)
{
  opCross = &cross;
}
