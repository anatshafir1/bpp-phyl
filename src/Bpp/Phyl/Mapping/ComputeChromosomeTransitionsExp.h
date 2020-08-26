//
// File: ComputeChromosomeTransitionsExp.h
// Created by: Anat Shafir
// Created on: Mon August 23 16:48 2020
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004, 2005, 2006)

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
#ifndef _COMPUTECHROMOSOMETRANSITIONSEXP_H_
#define _COMPUTECHROMOSOMETRANSITIONSEXP_H_


#include "Bpp/Phyl/Model/ChromosomeSubstitutionModel.h"
#include "Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h"
#include "Bpp/Phyl/Mapping/ComputeChangesExpectations.h"
#include "Bpp/Phyl/Likelihood/MarginalNonRevAncestralStateReconstruction.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <vector>
#include <map>
#include <utility>
#include <string>

using namespace std;

namespace bpp
{
    class ComputeChromosomeTransitionsExp{
        private:
            DRNonHomogeneousTreeLikelihood* lik_;
            MarginalNonRevAncestralStateReconstruction* ancr_;
            const TreeTemplate<Node>* tree_;
            const ChromosomeSubstitutionModel* model_;
            const ChromosomeAlphabet* alphabet_;
            std::map <int, std::map <int, double>> expNumOfChangesPerBranch_;
            std::map <int, double> expNumOfChanges_;
            ComputeChangesExpectations* simExpectations_;
            int jumpTypeMethod_;
        public:
            ComputeChromosomeTransitionsExp(DRNonHomogeneousTreeLikelihood* lik, MarginalNonRevAncestralStateReconstruction* ancr, ComputeChangesExpectations* simExpectations, int method = 0)
            :lik_(lik), ancr_(ancr), tree_(dynamic_cast<const TreeTemplate<Node>*>(&(lik->getTree()))), model_(dynamic_cast <const ChromosomeSubstitutionModel*>(lik->getSubstitutionModelSet()->getModel(0))), alphabet_(dynamic_cast <const ChromosomeAlphabet*>(lik->getAlphabet())),
            expNumOfChangesPerBranch_(), expNumOfChanges_(), simExpectations_(simExpectations), jumpTypeMethod_(method){}

            ComputeChromosomeTransitionsExp(const ComputeChromosomeTransitionsExp& exp):
                lik_(exp.lik_),
                ancr_(exp.ancr_),
                tree_ (exp.tree_),
                model_ (exp.model_),
                alphabet_(exp.alphabet_),
                expNumOfChangesPerBranch_ (exp.expNumOfChangesPerBranch_),
                expNumOfChanges_ (exp.expNumOfChanges_),
                simExpectations_(exp.simExpectations_),
                jumpTypeMethod_(exp.jumpTypeMethod_)
            {}
            ComputeChromosomeTransitionsExp& operator=(const ComputeChromosomeTransitionsExp& exp){
                lik_ = exp.lik_;
                ancr_ = exp.ancr_;
                tree_ = exp.tree_;
                model_ = exp.model_;
                alphabet_ = exp.alphabet_;
                expNumOfChangesPerBranch_ = exp.expNumOfChangesPerBranch_;
                expNumOfChanges_ = exp.expNumOfChanges_;
                simExpectations_ = exp.simExpectations_;
                jumpTypeMethod_ = exp.jumpTypeMethod_;
                return *this;
            }
            virtual ~ComputeChromosomeTransitionsExp(){};
            void init();
            void computeExpectationOfChangePerBranch(int nodeId, int jumpStateStart, int jumpStateEnd, VVdouble &jointProbFatherNode, int transitionType = -1);
            ChromosomeSubstitutionModel::typeOfTransition getTypeOfTransition(int startState, int endState);
            //more sophisticated function: if there is an overlap between different transition types-> the chosen state is sampled according to probabilities
            ChromosomeSubstitutionModel::typeOfTransition getTypeOfTransitionWithProb(int startState, int endState);
            void computeBaseNumExpectation(int nodeId, int jumpStateStart, VVdouble jointProbFatherNode);
            void computeDemiDuplExpectation(int nodeId, int jumpStateStart, VVdouble jointProbFatherNode);

            void computeExpectationPerType();
            void printResults();
            


    };
}
#endif // _COMPUTECHROMOSOMETRANSITIONSEXP_H_