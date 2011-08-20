/*****************************************************************
 ** SGAL Simple Genetic Algorithm Library
 **
 ** Copyright (c) 1995, Neal A. Sanche
 *****************************************************************/

#ifndef CHROM_H_
#define CHROM_H_

#include <stdlib.h>
#include <iostream>

// Chromosome base class, for polymorphism

//////////////////////////////////////////////////////////////////
// CChromosomeBase class declaration

class CChromosomeBase {
public:
  CChromosomeBase(): dFitness(0.0), dObjective(0.0), nLength(0), nBits(0) {}
  virtual ~CChromosomeBase() {}

  // Query the length and the number of bits
  virtual int GetLength() { return nLength; }
  virtual int GetNumBits() { return nBits; }

  // Get and set the value of a specific bit
  virtual int GetBit(int) = 0;
  virtual void SetBit(int, int) = 0;

  // Get and set fitness value (possibly scaled)
  virtual const double GetFitness() const { return dFitness; }
  virtual void SetFitness(const double fitness) { dFitness = fitness; }

  // Get and set raw objective funtion values
  virtual const double GetObjective() const { return dObjective; }
  virtual void SetObjective(const double obj) { dObjective = obj; }

  // Randomize the chromosome
  virtual void Randomize() = 0;

  // Make a chromosome equal to another
  virtual CChromosomeBase& operator=(const CChromosomeBase&) = 0;

protected:
  double dFitness;
  double dObjective;
  int nLength;
  int nBits;
};

//////////////////////////////////////////////////////////////////
// CChromosome<TGENE> class declaration

template<class TGENE>
class CChromosome: public CChromosomeBase {
public:
  CChromosome();
  CChromosome(int length);
  CChromosome(const CChromosome&);
  CChromosomeBase& operator=(const CChromosomeBase&);
  virtual ~CChromosome();
  
  TGENE &operator[](int index);  // Array operator
  
  virtual int GetBit(int index);
  virtual void SetBit(int index, int value);
  virtual void Randomize();
protected:
  TGENE *geneData;
};


//////////////////////////////////////////////////////////////////
// CChromosome<TGENE> class implementation

template<class TGENE>
CChromosome<TGENE>::CChromosome(): CChromosomeBase()
{
  nLength = 0;
  nBits = 8*sizeof(TGENE);
  geneData = 0;
}

template<class TGENE>
CChromosome<TGENE>::CChromosome(int length)
{
  nLength = length;
  nBits = 8*sizeof(TGENE) * nLength;
  
  geneData = new TGENE[nLength];
  
  for(int i = 0; i < nLength; i++) {
    geneData[i] = (TGENE)0;
  }
}  

template<class TGENE>
CChromosome<TGENE>::CChromosome(const CChromosome& that)
{
  nLength = that.nLength;
  geneData = 0;
  *this = that;
}

template<class TGENE>
CChromosomeBase& CChromosome<TGENE>::operator=(const CChromosomeBase& thatbase)
{
  CChromosome<TGENE> &that = (CChromosome<TGENE> &)thatbase;

  if ((geneData != 0) && (nLength > 0)) {
    delete[] geneData;
  }
  nLength = that.nLength;
  nBits = 8*sizeof(TGENE) * nLength;
  if (nBits != that.nBits) {
    std::cerr << "Bad data type in copy." << std::endl;
    return *this;
  }
  
  geneData = new TGENE[nLength];
  for (int i = 0; i < nLength; i++) {
    geneData[i] = that.geneData[i];
  }
  return *this;
}

template<class TGENE>
CChromosome<TGENE>::~CChromosome()
{
  if ((geneData != 0) && (nLength > 0)) {
    delete[] geneData;
  }
}

template<class TGENE>
TGENE &CChromosome<TGENE>::operator[](int index)  // Array operator
{
  return geneData[index];
}

template<class TGENE>
int CChromosome<TGENE>::GetBit(int index)
{
  // Do range checking
  if (index < 0 || index >= nBits) return -1;
  
  // Figure out which gene we're mangling
  int iGene = index/8/sizeof(TGENE);
  int iBit = index % (sizeof(TGENE)*8);
  
  char *acGene = (char *)&geneData[iGene];
  if (acGene[iBit/8] & ((char)1 << (iBit%8))) return 1;
  else return 0;
}

template<class TGENE>
void CChromosome<TGENE>::SetBit(int index, int value)
{
  // Do range checking
  if ((index < 0) || (index >= nBits)) return;
  
  // Figure out which gene we're mangling
  int iGene = index/8/sizeof(TGENE);
  int iBit = index % (sizeof(TGENE)*8);
  char val = 0;
  if (value > 0) val = 1;
  
  char *acGene = (char *)&geneData[iGene];
  acGene[iBit/8] = acGene[iBit/8] | ((char)val << (iBit%8));
}

template<class TGENE>
void CChromosome<TGENE>::Randomize()
{
  for (int i = 0; i < nBits; i++) {
    float test = (float)rand() / (float)RAND_MAX;
    if (test >= 0.5) {
      SetBit(i,1);
    }
    else {
      SetBit(i,0);
    }
  }
}

/////////////////////////////////////////////////////////////////
// CChromosomeFactoryBase class definition

class CChromosomeFactoryBase
{
public:
  CChromosomeFactoryBase() {}
  virtual CChromosomeBase *Create() = 0;
  virtual CChromosomeBase *Create(int) = 0;
};

/////////////////////////////////////////////////////////////////
// CChromosomeFactory<TGENE> class definition

template<class TGENE>
class CChromosomeFactory : public CChromosomeFactoryBase {
public:
  CChromosomeFactory() {}
  virtual CChromosomeBase *Create();
  virtual CChromosomeBase *Create(int length);
};

/////////////////////////////////////////////////////////////////
// CChromosomeFactory<TGENE> class implementation

template<class TGENE>
CChromosomeBase *CChromosomeFactory<TGENE>::Create()
{
  // Create a chromosome for use in the other
  // parts of the program
  return new CChromosome<TGENE>;
}

template<class TGENE>
CChromosomeBase *CChromosomeFactory<TGENE>::Create(int length)
{
  // Create a chromosome for use in the other
  // parts of the program
  return new CChromosome<TGENE>(length);
}

#endif // CHROM_H_
