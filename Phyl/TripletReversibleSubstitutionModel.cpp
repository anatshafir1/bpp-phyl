//
// File: TripletReversibleSubstitutionModel.cpp
// Created by:  Laurent Gueguen
// Created on: Feb 2009
//

/*
Copyright or � or Copr. CNRS, (November 16, 2004)
This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "TripletReversibleSubstitutionModel.h"

// From SeqLib:
#include <Seq/SequenceContainerTools.h>
#include <Seq/WordAlphabet.h>

// From NumCalc:
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/EigenValue.h>

using namespace bpp;

// From the STL:
#include <cmath>

using namespace std;

/******************************************************************************/

TripletReversibleSubstitutionModel::TripletReversibleSubstitutionModel(const CodonAlphabet* palph,
                                                                       NucleotideSubstitutionModel* pmod)  : WordReversibleSubstitutionModel(palph,"Triplet.")
{
  _VAbsRevMod.push_back(pmod);
  _VAbsRevMod.push_back(pmod);
  _VAbsRevMod.push_back(pmod);

  int i;

  addParameters_(pmod->getParameters());

  _rate=new double[3];
  for (i=0;i< 3; i++)
    _rate[i]=1.0/3;

}

TripletReversibleSubstitutionModel::TripletReversibleSubstitutionModel(const CodonAlphabet* palph,
                                                                       NucleotideSubstitutionModel* pmod1,  
                                                                       NucleotideSubstitutionModel* pmod2,
                                                                       NucleotideSubstitutionModel* pmod3) : WordReversibleSubstitutionModel(palph,"Triplet.")
{
  string st="Triplet.";
  
  if ((pmod1==pmod2) || (pmod2==pmod3) || (pmod1==pmod3)){
    int i;
    for (i=0;i< 3; i++){
      _VAbsRevMod.push_back(pmod1);
      _VnestedPrefix.push_back(pmod1->getNamespace());
    }
  
    pmod1->setNamespace(st+"012_"+_VnestedPrefix[0]);
    addParameters_(pmod1->getParameters());
  }
  else {
    _VAbsRevMod.push_back(pmod1);
    _VnestedPrefix.push_back(pmod1->getNamespace());
    _VAbsRevMod[0]->setNamespace(st+"0_"+_VnestedPrefix[0]);
    addParameters_(pmod1->getParameters());
    
    _VAbsRevMod.push_back(pmod2);
    _VnestedPrefix.push_back(pmod2->getNamespace());
    _VAbsRevMod[1]->setNamespace(st+"1_"+_VnestedPrefix[1]);
    addParameters_(pmod2->getParameters());
    
    _VAbsRevMod.push_back(pmod3);
    _VnestedPrefix.push_back(pmod3->getNamespace());
    _VAbsRevMod[2]->setNamespace(st+"2_"+_VnestedPrefix[2]);
    addParameters_(pmod3->getParameters());
  }
  
  _rate=new double[3];
  for (unsigned int i=0;i< 3; i++)
    _rate[i]=1.0/3;

}
  
TripletReversibleSubstitutionModel::TripletReversibleSubstitutionModel(const TripletReversibleSubstitutionModel& wrsm) : WordReversibleSubstitutionModel(wrsm.alphabet_,"Triplet")
{
  int i;

  _VAbsRevMod.push_back(wrsm._VAbsRevMod[0]);
  _VnestedPrefix.push_back(wrsm._VnestedPrefix[0]);

  _VAbsRevMod.push_back(wrsm._VAbsRevMod[1]);
  _VnestedPrefix.push_back(wrsm._VnestedPrefix[1]);

  _VAbsRevMod.push_back(wrsm._VAbsRevMod[2]);
  _VnestedPrefix.push_back(wrsm._VnestedPrefix[2]);

  addParameters_(wrsm.getParameters());

  _rate=new double[3];
  for (i=0;i< 3; i++)
    _rate[i]=1.0/3;
}

string TripletReversibleSubstitutionModel::getName() const
{
  string s="TripletReversibleSubstitutionModel model:";
  for (int i=0; i< _VAbsRevMod.size(); i++)
    s+=" "+ _VAbsRevMod[i]->getName();
  
  return s;
}


