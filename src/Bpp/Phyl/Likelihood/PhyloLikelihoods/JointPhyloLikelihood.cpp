//
// File: JointPhyloLikelihood.cpp
// Authors:
//   Anat Shafir
// Created: 2023-04-16 11:57:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#include "JointPhyloLikelihood.h"

using namespace std;
using namespace bpp;

JointPhyloLikelihood::JointPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings,
 bool inCollection) :
  AbstractPhyloLikelihood(context),
  SetOfAbstractPhyloLikelihood(context, pC, {}, inCollection),
  likCal_(new LikelihoodCalculation(context)),
  traitOptimization_(true),
  expectedHistory_(expectedHistory),
  numOfMappings_(numOfMappings),
  weightedFrequencies_(weightedFrequencies),
  tempLik_(0),
  firstLikChange_(true)

{
  addPhyloLikelihood(1, "_1");
  addPhyloLikelihood(2, "");
  tempLik_ = dynamic_cast<SingleProcessPhyloLikelihood*>(getAbstractPhyloLikelihood(nPhylo_[1]));
  likCal_->setLikelihoodNode(makeLikelihoods());



}


JointPhyloLikelihood::JointPhyloLikelihood(const JointPhyloLikelihood& sd) :
  AbstractPhyloLikelihood(sd),
  SetOfAbstractPhyloLikelihood(sd),
  likCal_(sd.likCal_),
  traitOptimization_(sd.traitOptimization_),
  expectedHistory_(sd.expectedHistory_),
  numOfMappings_(sd.numOfMappings_),
  weightedFrequencies_(sd.weightedFrequencies_),
  tempLik_(sd.tempLik_),
  firstLikChange_(sd.firstLikChange_)
{}

void JointPhyloLikelihood::fireParameterChanged(const ParameterList& params)
  {
    //getAbstractPhyloLikelihood(nPhylo_[0])->matchParametersValues(params);
    // if we optimize trait parameters, we need to change the tree because of the mappings based one new parameters
    if (expectedHistory_){
      if (traitOptimization_){
        // getting trait an updated trait model
        getAbstractPhyloLikelihood(nPhylo_[0])->matchParametersValues(params);
        getAbstractPhyloLikelihood(nPhylo_[0])->getValue();
        StochasticMapping* stm = new StochasticMapping(std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(getAbstractPhyloLikelihood(nPhylo_[0])->getLikelihoodCalculation()), numOfMappings_);
        double seedUb = 10000000;
        RandomTools::setSeed(static_cast<long int>(seedUb));
        stm->generateStochasticMapping();
        auto mappings = stm->createMappingHistoryTrees();
        auto expectedMapping = stm->generateExpectedMapping(mappings);
        //size_t numberOfStates = std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(getAbstractPhyloLikelihood(nPhylo_[0])->getLikelihoodCalculation())->getStateMap().getNumberOfModelStates();
        auto nodes = expectedMapping->getAllNodes();
        // creating map of states with the relevant branches (for the creation of the heterogeneous model)
        std::map <uint, std::vector<uint>> nodeModels;
        for (size_t i = 0; i < nodes.size(); i ++){
          uint nodeId = expectedMapping->getNodeIndex(nodes[i]);
          int nodeState = stm->getNodeState(nodes[i]);
          if (!(expectedMapping->isLeaf(nodeId))){
            auto sons = expectedMapping->getSons(nodeId);
            for (size_t j = 0; j < sons.size(); j++){
              nodeModels[static_cast<uint>(nodeState)].push_back(sons[j]);
            }
          }
        }
        // creating the new likelihood object with the different tree
        const NonHomogeneousSubstitutionProcess* prevSubstitutionModel = dynamic_cast<const NonHomogeneousSubstitutionProcess*>(&(std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>((tempLik_)->getLikelihoodCalculation()))->getSubstitutionProcess());
        ValueRef <Eigen::RowVectorXd> rootFreqs = std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation())->getRootFreqs();
        auto rootFreqsValues =  rootFreqs->getTargetValue();
        Vdouble rootFreqsBpp;
        copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(prevSubstitutionModel->getModel(1)->getStateMap(), false)), rootFreqsBpp);
        std::shared_ptr<FrequencySet> rootFrequencies = std::shared_ptr<FrequencySet>(rootFreqsFixed->clone());

        ParametrizablePhyloTree tree =  ParametrizablePhyloTree(*expectedMapping);
        std::shared_ptr<ParametrizablePhyloTree> parTree = std::shared_ptr<ParametrizablePhyloTree>((&tree)->clone());
        std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(prevSubstitutionModel->getRateDistribution()->clone());
        std::shared_ptr<NonHomogeneousSubstitutionProcess> subPro = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);
        // now adding the models- there should be numberOfStates models
        size_t index = 1;
        auto it = nodeModels.begin();
        while (it != nodeModels.end()){
          auto model = prevSubstitutionModel->getModel(index);
          subPro->addModel(std::shared_ptr<BranchModel>(model->clone()), nodeModels[it->first]);
          index ++;
          it ++;
        }
        delete stm;
        SubstitutionProcess* nsubPro= subPro->clone();
        Context* context = new Context();
        auto data = std::dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tempLik_->getLikelihoodCalculation())->getData();
        auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *data->clone(), *nsubPro, weightedFrequencies_);
        SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
        auto lik_to_del = tempLik_;
        tempLik_ = newLik;
        if (!(firstLikChange_)){
          delete lik_to_del;
        }else{
          firstLikChange_ = false;
        }



        
      }
      tempLik_->matchParametersValues(params);
      likCal_->setLikelihoodNode(makeLikelihoods());
      tempLik_->getLikelihoodNode();
    
    }

  }

 ValueRef<DataLik> JointPhyloLikelihood::makeLikelihoods(){
    auto phylo1 = getAbstractPhyloLikelihood(nPhylo_[0]);
    auto lik1 = phylo1->getLikelihoodNode();
    auto lik2 = tempLik_->getLikelihoodNode();
    // add the likelihoods because the target values are log likelihoods
    auto jointLik = CWiseAdd<DataLik, std::tuple<DataLik, DataLik> >::create(context_, {lik1, lik2}, Dimension<DataLik> ());

    return jointLik;

 }

// // this function can be abstract here, and then redefined by the derived classes
//  void JointPhyloLikelihood::optimizeTraitModel(double tol, uint numOfIterations){
//   traitOptimization_ = true;
//   OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
//   OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));
//   OptimizationTools::optimizeNumericalParameters2(
//     this, getParameters(), 0,
//     tol, numOfIterations, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);
//   delete profiler;
//   delete messenger;

//  }

// void JointPhyloLikelihood::optimizeSecondModel(double tol, uint numOfIterations){
//   traitOptimization_ = false;
//   OutputStream* profiler  = new StlOutputStream(new ofstream("profile.txt", ios::out));
//   OutputStream* messenger = new StlOutputStream(new ofstream("messages.txt", ios::out));
//   OptimizationTools::optimizeNumericalParameters2(
//   this, getParameters(), 0,
//   tol, numOfIterations, messenger, profiler, false, false, 1, OptimizationTools::OPTIMIZATION_NEWTON);
//   delete profiler;
//   delete messenger;
  
// }

//  void JointPhyloLikelihood::optimize(double tol, uint numOfIterationsPerModel, uint numberOfIterationsPerOneOptimization){
//   for (uint i = 0; i < numberOfIterationsPerOneOptimization; i++){
//     optimizeTraitModel(tol, numOfIterationsPerModel);
//     optimizeSecondModel(tol,numOfIterationsPerModel);

//   }

//  }