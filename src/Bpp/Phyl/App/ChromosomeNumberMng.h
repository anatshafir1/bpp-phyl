//
// File: ChromosomeNumberMng.h
// Created by: Anat Shafir
// Created on: Mon September 11 14:57 2020
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
#ifndef _CHROMOSOMENUMBERMNG_H_
#define _CHROMOSOMENUMBERMNG_H_


#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>



//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/chrFasta.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

//from bpp-phyl
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Tree/PhyloTreeTools.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/Io/IoTree.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/MarginalAncestralReconstruction.h>
#include <Bpp/Phyl/NewLikelihood/JointMLAncestralReconstruction.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>
//#include <Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
//#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>
//#include <Bpp/Phyl/Likelihood/MLAncestralStateReconstruction.h>
//#include <Bpp/Phyl/Likelihood/MarginalNonRevAncestralStateReconstruction.h>
#include <Bpp/Phyl/Likelihood/ChromosomeNumberOptimizer.h>
#include <Bpp/Phyl/Mapping/ComputeChromosomeTransitionsExp.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>


//standard libraries
#include <string>
#include <vector>
#include <iostream>
#include <time.h>

using namespace std;
namespace bpp{
    class ChromosomeNumberMng{
        private:
            PhyloTree* tree_;
            ChromosomeAlphabet* alphabet_;
            VectorSiteContainer* vsc_;
            std::map<uint, uint> chrRange_; //maxObserved-minObserved chromosome number
            unsigned int numberOfUniqueStates_; // number of unique states (number of chromosomes)



        public:
            //constructor
            ChromosomeNumberMng(): tree_(0), alphabet_(0), vsc_(0), chrRange_(), numberOfUniqueStates_(0){}
            ChromosomeNumberMng(const ChromosomeNumberMng& mng):
                tree_(mng.tree_->clone()), alphabet_(mng.alphabet_->clone()), vsc_(mng.vsc_->clone()), chrRange_(mng.chrRange_), numberOfUniqueStates_(mng.numberOfUniqueStates_)
            {}
            ChromosomeNumberMng& operator=(const ChromosomeNumberMng& mng){
                tree_ = tree_->clone();
                alphabet_ = mng.alphabet_->clone();
                vsc_ = mng.vsc_->clone();
                chrRange_ = mng.chrRange_;
                numberOfUniqueStates_ = mng.numberOfUniqueStates_;
                return *this;
            }
            virtual ~ChromosomeNumberMng(){
                if (vsc_){
                    delete vsc_;
                }
                delete tree_;
                delete alphabet_;
            };

            //Functions for initialization of the model
            void getCharacterData(const string &path);
            static void setMaxChrNum(unsigned int maxNumberOfChr);
            static void setMinChrNum(unsigned int minNumberOfChr);
            void getTree(const string &path, double treeLength = 0);

            // getters for testers
            const ChromosomeAlphabet* getAlphabet() const {return alphabet_;}
            const VectorSiteContainer* getSeqData() const {return vsc_;}
            const map<uint, uint> getChromosomeRange() const {return chrRange_;}
            const PhyloTree* getPhyloTree() const {return tree_;}
            


            //core functions of ChromEvol
            //void runTest();
            void runChromEvol();
            ChromosomeNumberOptimizer* optimizeLikelihoodMultiStartPoints() const;
            void getJointMLAncestralReconstruction(ChromosomeNumberOptimizer* optimizer) const;
            void getMarginalAncestralReconstruction(ChromosomeNumberOptimizer* chrOptimizer, const string &filePath);
            // map<int, map<size_t, VVdouble>> getMarginalAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik) const;
            void computeExpectations(ChromosomeNumberOptimizer* chrOptimizer, int numOfSimulations) const;
            void simulateData();
            void printSimulatedData(vector<size_t> leavesStates, vector<string> leavesNames, size_t iter);
            void printTreeWithStates(PhyloTree tree, std::map<uint, std::vector<size_t>> &ancestors, const string &filePath) const;
            void convertNodesNames(PhyloTree &tree, uint nodeId, std::map<uint, std::vector<size_t>> &ancestors, bool alphabetStates = true) const;
            void writeOutputToFile(ChromosomeNumberOptimizer* chrOptimizer) const;
            void printLikParameters(ChromosomeNumberOptimizer* chrOptimizer, SingleProcessPhyloLikelihood* lik, ofstream &outFile) const;
            uint findMinCladeSize(std::map<uint, vector<uint>> mapModelNodesIds) const;
            std::map<uint, uint> findMRCAForEachModelNodes(std::map<uint, vector<uint>> mapOfModelsAndNodes) const;
            void writeTreeWithCorrespondingModels(PhyloTree tree, std::map<uint, vector<uint>> &modelAndNodes) const;

        protected:
            void setNodeIdsForAllModels(string &path);
            //shared_ptr<PhyloNode> getMRCA(PhyloTree* tree, std::vector<shared_ptr<PhyloNode>> nodes);
            void getNodeIdsPerModelFromLine(string &content, PhyloTree* tree, std::map<uint, std::pair<uint, std::vector<uint>>> &modelAndNodeIds);
            std::shared_ptr<LikelihoodCalculationSingleProcess> setHeterogeneousLikInstance(SingleProcessPhyloLikelihood* likProcess, ParametrizablePhyloTree* parTree, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, bool forAncestral = false) const;
            std::shared_ptr<NonHomogeneousSubstitutionProcess> setHeterogeneousModel(ParametrizablePhyloTree* tree, SingleProcessPhyloLikelihood* ntl, ValueRef <Eigen::RowVectorXd> rootFreqs,  std::map<int, vector<pair<uint, int>>> sharedParams) const;
            VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* initialAlpha);
            void rescale_tree(PhyloTree* tree, double chrRange);
            void getMaxParsimonyUpperBound(double* parsimonyScore) const;
            // functions to print the tree with ancestral reconstruction
            void printSimulatedDataAndAncestors(SiteSimulationResult* simResult) const;
            void printSimulatedEvoPath(const string outPath, SiteSimulationResult* simResult) const;
            static string printTree(const PhyloTree& tree);
            static string nodeToParenthesis(const uint nodeId, const PhyloTree& tree);
            std::map<int, vector <double>> getVectorToSetModelParams(SingleProcessPhyloLikelihood* lik, size_t modelIndex = 1) const;
            double getOriginalTreeLength(string &path) const;
            // void printPosteriorProbNodes(std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, vector<double>& rootPosterior) const;




    };
}
#endif // _CHROMOSOMENUMBERMNG_H_