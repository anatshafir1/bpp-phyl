//
// File: ComputeChangesExpectations.h
// Created by: Anat Shafir
// Created on: Mon August 22 15:46 2020
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


#ifndef _COMPUTECHANGESEXPECTATIONS_H_
#define _COMPUTECHANGESEXPECTATIONS_H_

#include "../Tree.h"
#include "../TreeTemplate.h"
#include "Bpp/Phyl/Model/SubstitutionModel.h"
#include "Bpp/Phyl/Model/AbstractSubstitutionModel.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <vector>
#include <map>
#include <utility>

using namespace std;

namespace bpp
{
    class ComputeChangesExpectations{
      private:
        const Alphabet* alphabet_;
        const TreeTemplate<Node>* tree_;
        const SubstitutionModel* model_;
        std::vector<double> waitingTimes_;// waiting time parameter in the exponential distribution
        std::map <int, std::map<std::pair<int, int>, int>> ancestralTerminalsCounts_; // for each node Id, map <<firstState, LastState>, occurence>
        typedef std::pair<std::pair<int,int>,std::pair<int, int>> pairOfpairs; // first= ancestral states termainals. Second = jumps between each pair of states
        std::map <int, map <pairOfpairs, double>> branchTransitionsExp_; //for each node, for each possible combination of ancestral terminals and jump states- the expected number of transitions per jumps
        VVdouble jumpProbs_;  // probability to transit from state i to state j Qij/sum{Qij}
        std::vector<Node> branchOrder_;
      public:
        ComputeChangesExpectations(const Alphabet* alpha, const TreeTemplate<Node>* tree, const SubstitutionModel* model)
        :alphabet_(alpha), tree_(tree), model_(model), waitingTimes_(), ancestralTerminalsCounts_(),
        branchTransitionsExp_(), jumpProbs_(), branchOrder_(){}
        ComputeChangesExpectations(const ComputeChangesExpectations& exp):
          alphabet_ (exp.alphabet_),
          tree_ (exp.tree_),
          model_ (exp.model_),
          waitingTimes_ (exp.waitingTimes_),
          ancestralTerminalsCounts_ (exp.ancestralTerminalsCounts_),
          branchTransitionsExp_(exp.branchTransitionsExp_),
          jumpProbs_ (exp.jumpProbs_),
          branchOrder_ (exp.branchOrder_)
        {}
        ComputeChangesExpectations& operator=(const ComputeChangesExpectations& exp){
          alphabet_ = exp.alphabet_;
          tree_ = exp.tree_;
          model_ = exp.model_;
          waitingTimes_ = exp.waitingTimes_;
          ancestralTerminalsCounts_ = exp.ancestralTerminalsCounts_;
          branchTransitionsExp_ = exp.branchTransitionsExp_;
          jumpProbs_ = exp.jumpProbs_;
          branchOrder_ = exp.branchOrder_;
          return *this;
        }
        virtual ~ComputeChangesExpectations(){};

        void init();
        void runIteration(int state);
        void computeExpectationAndPosterior();
        void runSimulations(int numOfSimulations);
        static bool compareBranches(Node& node1, Node& node2);
        int getRandomState(int currentState);
        double getExpectation(int nodeId, int startAncestral, int endAncestral, int jumpStart, int jumpEnd);
     
    };

}
#endif // _COMPUTECHANGESEXPECTATIONS_H_