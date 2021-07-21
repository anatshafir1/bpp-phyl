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
#include "Bpp/Phyl/Tree/PhyloTree.h"
#include "Bpp/Phyl/Tree/PhyloTreeTools.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/VectorTools.h>

// From Seqlib:
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <vector>
#include <map>
#include <utility>
#include <string>

#define THRESHOLD_EXP 0.5
#define THRESHOLD_HEURISTIC 0.95
#define MAX_ITER_HEURISTICS 10
#define MAX_SIM_HEURISTICS 1000
#define BRANCH_MULTIPLIER_FACTOR 2

using namespace std;

namespace bpp
{
    class ComputeChromosomeTransitionsExp{
        // a class used to calculate the expectation of transitions per each type of 
        // transtion along the tree or per branch
        typedef pair<uint, PhyloBranch> Branch;
        private:
            map<uint, map<size_t, VVdouble>> jointProbabilitiesFatherSon_;
            const PhyloTree* tree_;
            const ChromosomeSubstitutionModel* model_;
            const ChromosomeAlphabet* alphabet_;
            vector<double> waitingTimes_;
            VVdouble jumpProbs_;  // probability to transit from state i to state j Qij/sum{Qij}
            vector<Branch> branchOrder_;
            map <uint, map<pair<int, int>, int>> ancestralTerminalsCounts_;  // for each node Id, map <<firstState, LastState>, occurence>
            map <uint, map <pair<int, int>, Vdouble>> branchTransitionsExp_; //for each node, for each possible pair of terminals-> Vdouble: the index is the type of transition. The double is the expectation
            map <uint, map <int, double>> expNumOfChangesPerBranch_;   // node -> jump type (gain, loss, dupl, demi-dupl, baseNum, maxChr) -> expecation
            map <uint, double> expNumOfChanges_;    //node->expectation per branch induced by the node
            int jumpTypeMethod_;    // which function to use for type classification- 0 if deterministic, 1 if probabilistic
            map <pair<int, int>, map<int, double>> stateJumpTypeProb_; // key = jump states i->j. value = map of change type and probability
            double isNeededHeuristics(uint nodeId, map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals); // check if we need to run heuristics (in case there are not enough transitions)
            double getCumulativeProbability(uint nodeId, vector <pair<int, int>>* terminalsToAccount = 0);
            vector <int> setVectorOfInitStatesForHeuristics(map <uint, vector<pair<int,int>>>& unAccountedNodesAndTerminals) const; // get the init states for which we have to rerun the simulation
            void updateNumNonAccountedBranches(map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals, int iteration, const string FilePath);
            void updateBranchLengths(int initState, int iteration, map <int, double>* ratesPerState);
            void getPosteriorAndExpForNonAccountedFor(map <uint, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun);
            void computeExpPerTypeHeuristics(map <uint, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun);
            bool isMaxStateValid(int prevState) const;
            
        public:
            ComputeChromosomeTransitionsExp(const std::shared_ptr<ChromosomeSubstitutionModel> model,  const PhyloTree* tree, const ChromosomeAlphabet* alphabet, map<uint, map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int method = 0)
            :jointProbabilitiesFatherSon_(jointProbabilitiesFatherSon), tree_(tree), model_(model.get()), alphabet_(alphabet),
            waitingTimes_(), jumpProbs_(), branchOrder_(), ancestralTerminalsCounts_(), branchTransitionsExp_(), expNumOfChangesPerBranch_(), expNumOfChanges_(), jumpTypeMethod_(method), stateJumpTypeProb_(){}

            ComputeChromosomeTransitionsExp(const ComputeChromosomeTransitionsExp& exp):
                jointProbabilitiesFatherSon_(exp.jointProbabilitiesFatherSon_),
                tree_ (exp.tree_),
                model_ (exp.model_),
                alphabet_(exp.alphabet_),
                waitingTimes_(exp.waitingTimes_),
                jumpProbs_(exp.jumpProbs_),
                branchOrder_(exp.branchOrder_),
                ancestralTerminalsCounts_(exp.ancestralTerminalsCounts_),
                branchTransitionsExp_(exp.branchTransitionsExp_),
                expNumOfChangesPerBranch_ (exp.expNumOfChangesPerBranch_),
                expNumOfChanges_ (exp.expNumOfChanges_),
                jumpTypeMethod_(exp.jumpTypeMethod_),
                stateJumpTypeProb_(exp.stateJumpTypeProb_)

            {}
            ComputeChromosomeTransitionsExp& operator=(const ComputeChromosomeTransitionsExp& exp){
                jointProbabilitiesFatherSon_ = exp.jointProbabilitiesFatherSon_;
                tree_ = exp.tree_;
                model_ = exp.model_;
                alphabet_ = exp.alphabet_;
                waitingTimes_ = exp.waitingTimes_;
                jumpProbs_ = exp.jumpProbs_;
                branchOrder_ = exp.branchOrder_;
                ancestralTerminalsCounts_ = exp.ancestralTerminalsCounts_;
                branchTransitionsExp_ = exp.branchTransitionsExp_;
                expNumOfChangesPerBranch_ = exp.expNumOfChangesPerBranch_;
                expNumOfChanges_ = exp.expNumOfChanges_;
                jumpTypeMethod_ = exp.jumpTypeMethod_;
                stateJumpTypeProb_ = exp.stateJumpTypeProb_;
                return *this;
            }
            ComputeChromosomeTransitionsExp* clone() const { return new ComputeChromosomeTransitionsExp(*this); }
            virtual ~ComputeChromosomeTransitionsExp(){};
            void init();
            void computeExpectationOfChangePerBranch(uint nodeId, VVdouble &jointProbFatherNode, int transitionType);
            ChromosomeSubstitutionModel::typeOfTransition getTypeOfTransition(int startState, int endState);
            //more sophisticated function: if there is an overlap between different transition types-> the chosen state is sampled according to probabilities
            ChromosomeSubstitutionModel::typeOfTransition getTypeOfTransitionWithProb(int startState, int endState);


            void computeExpectationPerType();
            void printResults(const string path = "none");
            PhyloTree* getResultTree();
            // // from previous used class
            void runIteration(int state, map <uint, vector<pair<int,int>>>* unAccountedNodesAndTerminals = 0);
            void computeExpectationAndPosterior();
            void runSimulations(int numOfSimulations);
            static bool compareBranches(Branch& edge1, Branch& edge2);//sorting function to sort the branches in ascending order of length
            int getRandomState(int currentState);
            double getExpectation(uint nodeId, int startAncestral, int endAncestral, int typeOfChange);
            void updateMapOfJumps(int startState, int endState);
            void updateExpectationsPerBranch(uint nodeId, pair<int, int> ancestralTerminals, pair<int, int> jumpStates);
            void runHeuristics(const string FilePath = "none");
            


    };
}
#endif // _COMPUTECHROMOSOMETRANSITIONSEXP_H_