//
// File: ChromosomeNumberOptimizer.h
// Created by: Anat Shafir
// Created on: Wednesday September 2 15:05 2020
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
#ifndef _CHROMOSOMENUMBEROPTIMIZER_H_
#define _CHROMOSOMENUMBEROPTIMIZER_H_

//from bpp-core
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>

//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>


//from bpp-phyl
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
//#include <Bpp/Phyl/NewLikelihood/JointMLAncestralReconstruction.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/App/ChromEvolOptions.h>
#include <Bpp/Phyl/OptimizationTools.h>
// From Seqlib:
#include <vector>
#include <map>
#include <utility>
#include <string>
#include <omp.h>
using namespace std;
namespace bpp
{
    class ChromosomeNumberOptimizer{
    // a class which is used for ChromEvol to run the likelihood optimization procedures
    // with different options available in ChromEvol

        private:
            vector <SingleProcessPhyloLikelihood*> vectorOfLikelohoods_;
            //vector <Context> vectorOfContexts_;
            const PhyloTree* tree_;
            const ChromosomeAlphabet* alphabet_;
            const VectorSiteContainer* vsc_;
            bool optimizeBaseNumber_;
            vector<unsigned int> numOfPoints_;
            vector<unsigned int> numOfIterations_;
            string typeOfOptimizer_;
            string baseNumOptimizationMethod_;
            mutable std::map<uint, uint> baseNumberUpperBound_;
            double tolerance_;
            bool standardOptimization_;
            int BrentBracketing_;
            vector <double> probsForMixedOptimization_;
            std::map<uint, vector<int>> fixedParams_;
            mutable std::map<int, std::vector<std::pair<uint, int>>> sharedParams_;
            
            

        public:
            ChromosomeNumberOptimizer(
                const PhyloTree* tree,
                const ChromosomeAlphabet* alpha,
                const VectorSiteContainer* vsc,
                std::map<uint, uint> baseNumberUpperBound):
                    vectorOfLikelohoods_(),
                    //vectorOfContexts_(),
                    tree_(tree),
                    alphabet_(alpha),
                    vsc_(vsc),
                    optimizeBaseNumber_(),
                    numOfPoints_(),
                    numOfIterations_(),
                    typeOfOptimizer_(),
                    baseNumOptimizationMethod_(),
                    baseNumberUpperBound_(baseNumberUpperBound),
                    tolerance_(),
                    standardOptimization_(),
                    BrentBracketing_(),
                    probsForMixedOptimization_(),
                    fixedParams_(),
                    sharedParams_()
            {}

            ChromosomeNumberOptimizer(const ChromosomeNumberOptimizer& opt):
                vectorOfLikelohoods_(opt.vectorOfLikelohoods_),
                //vectorOfContexts_(opt.vectorOfContexts_),
                tree_ (opt.tree_),
                alphabet_(opt.alphabet_),
                vsc_(opt.vsc_),
                optimizeBaseNumber_(opt.optimizeBaseNumber_),
                numOfPoints_(opt.numOfPoints_),
                numOfIterations_(opt.numOfIterations_),
                typeOfOptimizer_(opt.typeOfOptimizer_),
                baseNumOptimizationMethod_(opt.baseNumOptimizationMethod_),
                baseNumberUpperBound_(opt.baseNumberUpperBound_),
                tolerance_(opt.tolerance_),
                standardOptimization_(opt.standardOptimization_),
                BrentBracketing_(opt.BrentBracketing_),
                probsForMixedOptimization_(opt.probsForMixedOptimization_),
                fixedParams_(opt.fixedParams_),
                sharedParams_(opt.sharedParams_)
            {}
            ChromosomeNumberOptimizer& operator=(const ChromosomeNumberOptimizer& opt){
                vectorOfLikelohoods_ = opt.vectorOfLikelohoods_;
                //vectorOfContexts_ = opt.vectorOfContexts_;
                tree_ = opt.tree_;
                alphabet_ = opt.alphabet_;
                vsc_ = opt.vsc_;
                optimizeBaseNumber_ = opt.optimizeBaseNumber_;
                numOfPoints_ = opt.numOfPoints_;
                numOfIterations_ = opt.numOfIterations_;
                typeOfOptimizer_ = opt.typeOfOptimizer_;
                baseNumOptimizationMethod_ = opt.baseNumOptimizationMethod_;
                baseNumberUpperBound_ = opt.baseNumberUpperBound_;
                tolerance_ = opt.tolerance_;
                standardOptimization_ = opt.standardOptimization_;
                BrentBracketing_ = opt.BrentBracketing_;
                probsForMixedOptimization_ = opt.probsForMixedOptimization_;
                fixedParams_ = opt.fixedParams_;
                sharedParams_ = opt.sharedParams_;
                return *this;
            }
            ChromosomeNumberOptimizer* clone() const { return new ChromosomeNumberOptimizer(*this); }
            virtual ~ChromosomeNumberOptimizer(){clearVectorOfLikelihoods(0);};
            //init models
                        // std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelComplexParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds
            //void initModels(std::map<uint, std::pair<int, map<int, std::vector<double>>>> modelComplexParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numberOfModels, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
        //     //initialize all the optimization specific members
            void initOptimizer(
                vector<unsigned int> numOfPoints,
                vector<unsigned int> numOfIterations,
                string typeOfOptimizer,
                string baseNumOptimizationMethod,
                double tolerance,
                bool standardOptimization,
                int BrentBracketing,
                vector <double>& probsForMixedOptimization)
            {
                numOfPoints_ = numOfPoints;
                numOfIterations_ = numOfIterations;
                typeOfOptimizer_ = typeOfOptimizer;
                baseNumOptimizationMethod_ = baseNumOptimizationMethod;
                tolerance_ = tolerance;
                standardOptimization_ = standardOptimization;
                BrentBracketing_ =BrentBracketing;
                probsForMixedOptimization_ = probsForMixedOptimization;
                

            }
            const std::map<int, std::vector<pair<uint, int>>> getSharedParams(){return sharedParams_;}
            const double getAICOfBestModel() const {
                std::map<uint, vector<int>> fixedParams = fixedParams_;
                size_t numOfFixedParams = getNumberOfFixedParams(vectorOfLikelohoods_[0], fixedParams);
                return calculateAICc(vectorOfLikelohoods_[0], numOfFixedParams);
            }
            void runNewBranchModel(omp_lock_t &mutex, SingleProcessPhyloLikelihood* lik, std::vector<SingleProcessPhyloLikelihood*> &newShiftLikCandidates, vector<uint> &candidateShiftNodesIds, size_t i, uint numOfShifts, double parsimonyBound, uint numOfPoints);
            void optimizeMultiProcessModel(std::map<int, std::vector<pair<uint, int>>>* sharedParams ,std::map<uint, vector<int>>* fixedParams, vector<SingleProcessPhyloLikelihood*>* perCandidateLikVec = 0);
            //void optimizeHeterogeneous();
            void optimize(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
            void optimizeInParallel(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds);
            vector<SingleProcessPhyloLikelihood*> getVectorOfLikelihoods(){return vectorOfLikelohoods_;}
            static vector <double> setFixedRootFrequencies(const std::string &path, std::shared_ptr<ChromosomeSubstitutionModel> chrModel);
            static std::map<uint, std::vector<string>> getRelatedParameterNamesForEachModel(ParameterList &params, std::string pattern, uint numOfModels, std::map<int, vector<std::pair<uint, int>>>* mapSharedParams = 0);
            //static std::map<uint, std::pair<int, std::map<int, vector<double>>>> getModelParameters(SingleProcessPhyloLikelihood* tl);
            // get the map of models and the corresponding nodes.
            static void getMutableMapOfModelAndNodeIds(std::map<uint, vector<uint>> &mapModelNodesIds, SingleProcessPhyloLikelihood* lik, uint rootId = 0);
            static std::map<uint, pair<int, std::map<int, std::vector<double>>>> getMapOfParamsForComplexModel(SingleProcessPhyloLikelihood* lik, std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames, uint numOfModels);
            static void updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, SingleProcessPhyloLikelihood* tl, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams = 0);
            //void writeOutputToFile() const;
            void printRootFrequencies(SingleProcessPhyloLikelihood* lik, ofstream &outFile) const;
            static std::string getStringParamName(int type);
            static uint getNumberOfParametersPerParamType(int paramType, vector<int> &funcTypes);
            static uint getModelFromParamName(string name);
            static int getTypeOfParamFromParamName(string name);
            static size_t getNumberOfFixedParams(SingleProcessPhyloLikelihood* lik, std::map<uint, vector<int>> &fixedParams);


        protected:
        //     // for model initiation
            //SingleProcessPhyloLikelihood* getLikelihoodFunction(const PhyloTree* tree, const VectorSiteContainer* vsc, std::shared_ptr<ChromosomeSubstitutionModel> &chrModel, DiscreteDistribution* rdist, const string& fixedRootFreqPath);
            
            
        //     // //functions of optimization
            static void createMapOfSharedParameterNames(std::map<int, std::vector<std::pair<uint, int>>> &sharedParams, std::map<string, vector<std::pair<uint, int>>> &sharedParamsNames);
            void initLikelihoods(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds, uint numOfModels, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams);
            static void setRandomPoints(SingleProcessPhyloLikelihood* lik, uint nodeToSplit, std::map<int, std::vector<uint>>* sharedParams, std::map<int, std::vector<uint>>* updatedSharedParams, int numOfPoints);
            static void updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName);
            static int getEnumOfParamName(std::string pattern);
            static void setParamsNameInForMultiProcess(std::map<uint, std::map<int, vector<string>>> &mapOfParamsNamesPerModelType, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams);
            static void aliasParametersInSubstitutionProcess(std::map<uint, std::map<int, vector<string>>> &mapOfParamsNamesPerModelType, std::map<int, vector<std::pair<uint, int>>>* updatedSharedParams, std::shared_ptr<NonHomogeneousSubstitutionProcess> process);
            
            unsigned int optimizeModelParameters(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector<unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams);//, unsigned int inwardBracketing, bool standardOptimization);
            unsigned int optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector<unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, bool mixed = false, unsigned int currentIterNum = 0);
            unsigned int optimizeMultiDimensions(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams, bool mixed = false, unsigned int currentIterNum = 0);
            unsigned int useMixedOptimizers(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<uint, vector<int>>* fixedParams);
            void optimizeBaseNum(SingleProcessPhyloLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, double upperBound, const string &paramName, ParameterList& params, uint model);

            // // function working on the likelihoods vector object
            void clearVectorOfLikelihoods(size_t new_size);
            void clearVectorOfLikelihoods(size_t new_size, std::vector<SingleProcessPhyloLikelihood*> &likelihoodsVec);
            void deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del);
            static bool compareLikValues(SingleProcessPhyloLikelihood* lik1, SingleProcessPhyloLikelihood* lik2);

            // // helper functions for optimization
            void checkLegalUseOfGradientOptimization();
            vector <string> getNonFixedParams(SingleProcessPhyloLikelihood* tl, ParameterList &allParams, map<uint, vector<int>>* fixedParams) const;
            void fillVectorOfBaseNumCandidates(vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const;
            uint getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) const;
            void getAllPossibleChrRanges(vector <unsigned int> &baseNumCandidates) const;
            //string findParameterNameInModel(string fullParameterName) const;
            //void setNewBounds(const ParameterList params, Parameter &param, map<string, pair<string, bool>> &paramPairsMap, double* lowerBound, const ChromosomeSubstitutionModel* model);

            //print functions
            void printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, const string path = "none") const;
           
            void printLikelihoodVectorValues(vector <SingleProcessPhyloLikelihood*> lik_vec) const;

            /*********************************************************
             * Functions associated with heterogenous ChromEvol models
            **********************************************************/
            double calculateAICc(SingleProcessPhyloLikelihood* lik, size_t numOfFixedParams) const;
            //getNewLikObject(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<int, uint>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound)
            void getNewLikObjectForParallelRuns(std::vector<SingleProcessPhyloLikelihood*> &newShiftLikCandidates, size_t index, std::vector<SingleProcessPhyloLikelihood*> &perCandidateLikVec, SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound);
            void getNewLikObject(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams, std::map<int, std::vector<std::pair<uint, int>>>* updatedSharedParams, uint numOfPoints, std::map<uint, vector<int>> &fixedParams, double parsimonyBound);
            static SingleProcessPhyloLikelihood* setHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, std::map<int, vector<std::pair<uint,int>>>* updatedSharedParams);
            static SingleProcessPhyloLikelihood* setRandomHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, std::map<uint, vector<int>> &fixedParams, std::map<int, std::vector<std::pair<uint, int>>>* sharedParams);
            //void optimizeSingleHeterogeneousModel(size_t index, int maxNumOfModels, std::vector<uint> &candidateShiftNodesIds, vector<uint> &baseNumCandidates);
            void getValidCandidatesForShift(std::vector<uint> &candidateShiftNodesIds, int minCladeSize);
            void updateSharedParameters(std::map<int, vector<std::pair<uint, int>>> &sharedParams, uint prevShift, uint numOfShifts) const;

    };
}
#endif // _CHROMOSOMENUMBEROPTIMIZER_H_