//
// File: StochasticMapping.cpp
// Authors:
//

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Text/TextTools.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric> // to sum over items in a vector
#include <cmath>

#include "../Simulation/MutationProcess.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"
#include "Reward.h"
#include "RewardMappingTools.h"
#include "StochasticMapping.h"

using namespace bpp;
using namespace std;

#define STATE "state"

/******************************************************************************/

StochasticMapping::StochasticMapping(std::shared_ptr<LikelihoodCalculationSingleProcess> drl, size_t numOfMappings, size_t numOfMappingTrials) :
  likelihood_(drl),
  tree_ (make_shared<PhyloTree>(*drl->getSubstitutionProcess().getParametrizablePhyloTree())),
//  mappingParameters_(drl->getSubstitutionProcess()),
  ConditionalProbabilities_(),
  nodesCounter_(0),
  numOfMappings_(numOfMappings),
  ancetralStates_(),
  mappings_(),
  jumpsProbs_(),
  notRepresentedNodes_(),
  numOfMappingTrials_(numOfMappingTrials),
  MLAncr_(0)// ,
  // nodeIdToIndex_()
{
  //giveNamesToInternalNodes(*tree_);                     // set names for the internal nodes of the tree, in case of absence
  ComputeConditionals();
  initJumpProbs();
  //initMappings();
}

/******************************************************************************/

StochasticMapping::~StochasticMapping()
{}

/******************************************************************************/
void StochasticMapping::sampleAllAncestals(){
  // initializing the ancestral states for each node and mapping
  auto nodeIndices = tree_->getNodeIndexes(tree_->getAllNodes());
  for(size_t i = 0; i < nodeIndices.size(); i++){
    ancetralStates_[nodeIndices[i]].resize(numOfMappings_);
  }
  for (size_t i = 0; i < numOfMappings_; ++i)
  {
    /* step 1+2: simulate a set of ancestral states, based on the fractional likelihoods from step 1 */
    
    sampleAncestrals(i);
  }
  

}


/******************************************************************************/

void StochasticMapping::generateStochasticMapping()
{
  // initializing the ancestral states for each node and mapping
  auto nodeIndices = tree_->getNodeIndexes(tree_->getAllNodes());
  for(size_t i = 0; i < nodeIndices.size(); i++){
    ancetralStates_[nodeIndices[i]].resize(numOfMappings_);
  }
  for (size_t i = 0; i < numOfMappings_; ++i)
  {
    /* step 1+2: simulate a set of ancestral states, based on the fractional likelihoods from step 1 */
    
    sampleAncestrals(i);

    /* step 3: simulate mutational history of each lineage of the phylogeny, conditional on the ancestral states */
    vector<uint> failedNodeIds;
    //bool success = sampleMutationsGivenAncestrals(i);
    bool success = sampleMutationsGivenAncestrals(i, &failedNodeIds);
    // if (!success){
    //   sampleAncestrals(i); // verify that it doesn't push any elements again
    //   clearMapping(i);
    //   success = sampleMutationsGivenAncestrals(i, &failedNodeIds);
    // }
    if (!success){
      for (size_t k = 0; k < failedNodeIds.size(); k++){
        notRepresentedNodes_[failedNodeIds[k]].push_back(i);
        auto fatherNode = tree_->getFatherOfNode(tree_->getNode(failedNodeIds[k]));
        uint father = tree_->getNodeIndex(fatherNode);
        size_t fatherState = ancetralStates_[father][i];
        auto branchPtr = tree_->getIncomingEdges(tree_->getNode(failedNodeIds[k]))[0];
        auto branchLength = branchPtr->getLength();
        auto alphabet = likelihood_->getData()->getAlphabet();
        MutationPath tryMapping(alphabet, fatherState, branchLength);
        // add an empty mmutation path instance
        mappings_[failedNodeIds[k]].push_back(tryMapping);
       
      }
  
    }
  }
}
/******************************************************************************/
void StochasticMapping::clearMapping(size_t mappingIndex){
  auto nodeIndices = tree_->getNodeIndexes(tree_->getAllNodes());
  for (size_t n = 0; n < nodeIndices.size(); n++){
    if (mappings_.find(nodeIndices[n]) == mappings_.end()){
      continue;
    }
    if (mappings_[nodeIndices[n]].size() == mappingIndex + 1){
      mappings_[nodeIndices[n]].pop_back();
    }
    else if (mappings_[nodeIndices[n]].size() > mappingIndex + 1){
      throw Exception("StochasticMapping::clearMapping(): mapping index out of range!");
    }
  }
}

/******************************************************************************/
void StochasticMapping::initJumpProbs(){
  auto nbState = likelihood_-> getStateMap().getNumberOfModelStates();
  auto nbModels = likelihood_->getSubstitutionProcess().getNumberOfModels();
  for (size_t modelId = 1; modelId <= nbModels; modelId++){
    if (jumpsProbs_.find(modelId) == jumpsProbs_.end()){
      jumpsProbs_[modelId].resize(nbState);
      for (size_t i = 0; i < nbState; i++){
        jumpsProbs_[modelId][i].resize(nbState);
        for (size_t j = 0; j < nbState; j++){
          auto model = dynamic_pointer_cast<const SubstitutionModel>(likelihood_->getSubstitutionProcess().getModel(modelId));
          if (i == j){
            jumpsProbs_[modelId][i][j] = 0;
          }else{
            jumpsProbs_[modelId][i][j] = model->Qij(i, j)/(-1*model->Qij(i,i));
          }
        }
      }

    }
  }
  
}
/******************************************************************************/
size_t StochasticMapping::giveRandomState(size_t beginState, size_t modelIndex) const 
{
  auto nbState = likelihood_-> getStateMap().getNumberOfModelStates();
	for (size_t i = 0 ; i < 100000 ; i++) 
	{
		double u = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
		double cumulativeProb = 0.0;
		for (size_t state = 0; state < nbState; state ++) 
		{
			cumulativeProb += jumpsProbs_.at(modelIndex)[beginState][state];
			if (u < cumulativeProb) {
				return state;
			}
		}
	}
	throw Exception("StochasticMapping::giveRandomState: could not give random character. The reason is unknown.");
	return 1;

}
/*****************************************************************************/
double StochasticMapping::getRateToLeaveState(uint nodeId, size_t mapping){
  auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
  uint father = tree_->getNodeIndex(fatherNode);
  auto model = dynamic_pointer_cast<const SubstitutionModel>(likelihood_->getSubstitutionProcess().getModel(father, 0)); // father or son??? Should be a father, because the models start at particular nodes, and I should get the model of the preceeding branch
  size_t fatherState = ancetralStates_[father][mapping];
  auto rateToLeave = -1* model->Qij(fatherState, fatherState);
  return rateToLeave;

}
/******************************************************************************/
bool StochasticMapping::setExpectedAncestrals(shared_ptr<PhyloTree> expectedMapping, std::map<uint, std::vector<double>> &ancestralStatesFrequencies){
  bool allAncestralsConsistent = true;
  // find all the possible combination of father son states
  std::map<uint, vector<pair<size_t,size_t>>> allowedStates; // son id is the key. The value is a vector of all possible combination of states, where the son is the
  // second in the pair, and the father state is the first in the pair
  auto nodeIndices = expectedMapping->getNodeIndexes(expectedMapping->getAllNodes());
  for (size_t n = 0; n < nodeIndices.size(); n++){
    if (expectedMapping->isLeaf(expectedMapping->getNode(nodeIndices[n]))){
      continue;
    }
    auto sons = expectedMapping->getSons(nodeIndices[n]);
    uint fatherId = nodeIndices[n];
    for (size_t j = 0; j < sons.size(); j++){
      allowedStates[sons[j]];
      for (size_t m = 0; m < numOfMappings_; m++){
        std::pair<size_t,size_t> fatherSonCombStates(ancetralStates_[fatherId][m], ancetralStates_[sons[j]][m]);
        auto itComb = std::find(allowedStates[sons[j]].begin(), allowedStates[sons[j]].end(), fatherSonCombStates);
        if (itComb == allowedStates[sons[j]].end()){
          allowedStates[sons[j]].push_back(fatherSonCombStates);

        }
      }
    }
  }
  std::map<uint, size_t> expectedAncestrals;
  for (size_t n = 0; n < nodeIndices.size(); n++){
    if (expectedMapping->isLeaf(expectedMapping->getNode(nodeIndices[n]))){
      continue;
    }
    if (expectedMapping->getRootIndex() == nodeIndices[n]){
      auto d = std::distance(ancestralStatesFrequencies[nodeIndices[n]].begin(),std::max_element(ancestralStatesFrequencies[nodeIndices[n]].begin(), ancestralStatesFrequencies[nodeIndices[n]].end()));
      if (MLAncr_){
        expectedAncestrals[nodeIndices[n]] = (*MLAncr_)[nodeIndices[n]][0];

      }else{
        expectedAncestrals[nodeIndices[n]] = static_cast<size_t>(d);

      }
    }
    auto sons = expectedMapping->getSons(nodeIndices[n]);
    uint father = nodeIndices[n];
    for (size_t j = 0; j < sons.size(); j++){
      size_t state;
      uint nodeId = sons[j];
      if (MLAncr_){
        state = (*MLAncr_)[nodeId][0];
      }else{
        auto d = std::distance(ancestralStatesFrequencies[nodeId].begin(),std::max_element(ancestralStatesFrequencies[nodeId].begin(), ancestralStatesFrequencies[nodeId].end()));
        state = static_cast<size_t>(d);
      }
      std::pair<size_t,size_t> statesFatherSon(expectedAncestrals[father], state);
      auto it = std::find(allowedStates[nodeId].begin(), allowedStates[nodeId].end(), statesFatherSon);
      if (it != allowedStates[nodeId].end()){
        expectedAncestrals[nodeId] = state;

      }else{
        allAncestralsConsistent = false;
        return allAncestralsConsistent;
      }
    }
  }
  auto nodes = expectedMapping->getAllNodes();
  for (size_t i = 0; i < nodes.size(); i++){
    uint currNodeId = expectedMapping->getNodeIndex(nodes[i]);
    if (expectedMapping->isLeaf(nodes[i])){
      nodes[i]->setName(nodes[i]->getName()+"-"+ std::to_string(expectedAncestrals[currNodeId]));
    }else{
      nodes[i]->setName("N"+ std::to_string(currNodeId)+"-"+ std::to_string(expectedAncestrals[currNodeId]));
    }
  }
  return allAncestralsConsistent;

}

/******************************************************************************/

shared_ptr<PhyloTree> StochasticMapping::generateExpectedMapping()
{
  sampleAllAncestals();
  // // initialize the expected history
  Newick writer;
  Newick reader;
  std::string tree_str = writer.writeTreeToParenthesis(*tree_);
  std::shared_ptr<PhyloTree> expectedMapping = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree(tree_str));
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());


  // compute a vector of the posterior asssignment probabilities for each inner node
  std::map<uint, std::vector<double>> ancestralStatesFrequencies;
  computeStatesFrequencies(ancestralStatesFrequencies);
  // set the ancestral states according to the maximal posterior (i.e, conditional) probability
  bool allConsistent = setExpectedAncestrals(expectedMapping, ancestralStatesFrequencies);
  // also need to handle cases where the most posterior states of father and son do not come together.
  size_t counter = 0;
  while((!allConsistent) && (counter < 10)){
    ancestralStatesFrequencies.clear();
    ancetralStates_.clear();
    sampleAllAncestals();
    computeStatesFrequencies(ancestralStatesFrequencies);
    allConsistent = setExpectedAncestrals(expectedMapping, ancestralStatesFrequencies);
    counter ++;
  }

  ConditionalProbabilities_.clear();
  if (!allConsistent){
    throw Exception("StochasticMapping::generateExpectedMapping(): did not find any consistent states between father and son!");
  }
  
  // find mappings
  auto nodeIndices = tree_->getNodeIndexes(tree_->getAllNodes());
  // calculate expected duration times
  //std::map<uint, std::vector<double>> dwellingTimes;
  for (size_t i = 0; i < nodeIndices.size(); i++){
    if (tree_->isLeaf(tree_->getNode(nodeIndices[i]))){
      continue;
    }
    //auto father = nodeIndices[i];
    auto sons = tree_->getSons(nodeIndices[i]);
    for (size_t j = 0; j < sons.size(); j++){
      uint fatherId;
      auto node = expectedMapping->getNode(sons[j]);
      auto fatherNode = expectedMapping->getFatherOfNode(node);
      if (expectedMapping->getRootIndex() == expectedMapping->getNodeIndex(fatherNode)){
        fatherId = tree_->getRootIndex();
      }else{
        fatherId = expectedMapping->getNodeIndex(fatherNode);
      }
      vector<MutationPath> nodeMappings;
      vector<double> nodeDwellingTimes;
      nodeDwellingTimes.resize(likelihood_-> getStateMap().getNumberOfModelStates());
      for (size_t mappingIndex= 0; mappingIndex < numOfMappings_; mappingIndex++){
        bool success = sampleMutationsGivenAncestralsPerBranch(fatherId, sons[j], mappingIndex, nodeMappings, numOfMappingTrials_);
        if (!success){
          throw Exception("Mapping failure!!!");
        }
        assignDewellingTimesUnderEachStatePerMappingPerBranch(sons[j], ancetralStates_[fatherId][mappingIndex], nodeDwellingTimes, nodeMappings[mappingIndex]);

      }
      // get dwelling times
      for (size_t s = 0; s < nodeDwellingTimes.size(); s++){
        nodeDwellingTimes[s] /= static_cast<double>(numOfMappings_);
      }

      size_t startState = static_cast<size_t>(getNodeState(fatherNode));
      size_t endState = static_cast<size_t>(getNodeState(node));
      std::map<pair<size_t, size_t>, double> transitionOcurrences;
      std::map<pair<size_t, size_t>, double> timeDurations;
      vector<size_t> mostFreqPath;
      getExpectedNumberOfTransitionsPerBranchGivenTerminals(sons[j], fatherId, startState, endState, transitionOcurrences, timeDurations, mostFreqPath, nodeMappings);
      vector<double> timeDurationsPerState;
      size_t nbStates = likelihood_-> getStateMap().getNumberOfModelStates();
      timeDurationsPerState.resize(nbStates);
      auto it = timeDurations.begin();
      while (it != timeDurations.end()){
        auto &transition = it->first;
        size_t currState = transition.first;
        timeDurationsPerState[currState] += timeDurations[transition];
        it ++;
      }
      findExpectedPathOnBranch(expectedMapping, startState, endState, fatherId, sons[j], transitionOcurrences, timeDurations, timeDurationsPerState, nodeMappings, mostFreqPath);
      
      
    }
  }
  ancetralStates_.clear();

  return expectedMapping;
}
/******************************************************************************/
void StochasticMapping::findExpectedPathOnBranch(std::shared_ptr<PhyloTree> expectedMapping, size_t fatherState, size_t sonState, uint fatherId, uint nodeId, std::map<pair<size_t, size_t>, double> &transitionOcurrences, std::map<pair<size_t, size_t>, double> &timeDurations, vector<double> &timeDurationsPerState, vector<MutationPath> &nodeMappings, vector<size_t> &mostFreqPath){  
  auto branch = expectedMapping->getEdgeToFather(nodeId);
  double branchLength = branch->getLength();
  bool foundPath = true;
  auto mappingStates = MultiStateMappingPath::findExpectedMappingPathForEachNode(fatherState, sonState, transitionOcurrences, timeDurationsPerState, branchLength, foundPath);
  if (mappingStates.size() > 0){
    mappingStates.push_back(sonState); // this is a dummy transition, just to create the transition of the last state to itself.

  }else{
    if (foundPath){
      return;

    }else{
      std::cout << "Most frequent path:\n";
      if (mostFreqPath.size() == 0){
        return;
      }
      for (auto &state : mostFreqPath){
        mappingStates.push_back(state);
        std::cout << state << ",";
      }
      mappingStates.push_back(sonState);
      std::cout << sonState << std::endl;
    }
      
  }
  std::map<pair<size_t, size_t>, double> occurrencesOfTrnasitionsInExpectedPath;
  for (size_t j = 0; j < mappingStates.size()-1; j++){
    std::pair<size_t,size_t> transition(mappingStates[j], mappingStates[j+1]);
    if (occurrencesOfTrnasitionsInExpectedPath.find(transition) != occurrencesOfTrnasitionsInExpectedPath.end()){
      occurrencesOfTrnasitionsInExpectedPath[transition] += 1;
    }else{
      occurrencesOfTrnasitionsInExpectedPath[transition] = 1;
    }
  }
  double sumOfChosenTransitionsTimes = 0;
  std::map<std::pair<size_t, size_t>, double> newBranchLengths;
  for (size_t j = 0; j < mappingStates.size()-1; j++){
    std::pair<size_t,size_t> transition(mappingStates[j], mappingStates[j+1]);
    newBranchLengths[transition] = timeDurations[transition]/occurrencesOfTrnasitionsInExpectedPath[transition];
    sumOfChosenTransitionsTimes += newBranchLengths[transition];
  }
  // now fragmenting the edge
  double segmentBranchLength;
  //double dwellingTime;
  uint newNodeId;
  double sumOfTransitionsTime = 0;
  std::pair<size_t,size_t> transition;
  for (size_t j = 0; j < mappingStates.size()-2; j++){
    transition = pair<size_t,size_t>(mappingStates[j], mappingStates[j+1]);
    segmentBranchLength = newBranchLengths[transition]/(sumOfChosenTransitionsTimes/branchLength);
    sumOfTransitionsTime += segmentBranchLength;
    auto edge_to_fragment = expectedMapping->getEdgeToFather(nodeId);
    newNodeId = expectedMapping->createNodeOnEdge(expectedMapping->getEdgeIndex(edge_to_fragment), segmentBranchLength);
    (expectedMapping->getNode(newNodeId))->setName("N_dummy_"+ std::to_string(newNodeId)+"-"+ std::to_string(mappingStates[j+1]));


  }
  auto lastTransition = pair<size_t,size_t>(mappingStates[mappingStates.size()-2], mappingStates[mappingStates.size()-1]);
  double estimatedRemained = newBranchLengths[lastTransition]/(sumOfChosenTransitionsTimes/branchLength);
  double truelyRemained = branchLength-sumOfTransitionsTime;
  double epsilon = 1e-6;
  if (std::abs(estimatedRemained - truelyRemained) > epsilon){
    throw Exception("StochasticMapping::findExpectedHistoryTransitionsAndTimeDurationsMultiState(): sum of segments is"+std::to_string(sumOfTransitionsTime)+ " while branch length is "+ std::to_string(branchLength) + "\n");
  }

}

/******************************************************************************/
void StochasticMapping::getTimeDurationsPerStateGivenAncestrals(std::map<uint, std::map<pair<size_t, size_t>, double>> &timeDurations, std::map<uint, std::vector<double>> &timeDurationsPerState){
  auto nodes = tree_->getAllNodes();
  auto nbStates = likelihood_->getStateMap().getNumberOfModelStates();
  for (size_t i = 0; i < nodes.size(); i++){
    uint nodeId = tree_->getNodeIndex(nodes[i]);
    if (nodeId == tree_->getRootIndex()){
      continue;
    }
    timeDurationsPerState[nodeId].resize(nbStates);
    auto &timeDurationsPerNode = timeDurations[nodeId];
    auto it = timeDurationsPerNode.begin();
    while (it != timeDurationsPerNode.end()){
      auto &transition = it->first;
      size_t currState = transition.first;
      timeDurationsPerState[nodeId][currState] += timeDurationsPerNode[transition];
      it ++;
    }

  }


}

/******************************************************************************/
void StochasticMapping::findExpectedHistoryTransitionsAndTimeDurationsMultiState(std::shared_ptr<PhyloTree> expectedMapping, std::map<uint, std::vector<double>> &dwellingTimes){
  std::map<uint, std::map<pair<size_t, size_t>, double>> transitionOcurrences;
  std::map<uint, std::map<pair<size_t, size_t>, double>> timeDurations;
  std::unordered_map<uint, vector<size_t>> mostFreqPaths;

  getExpectedNumberOfTransitionsPerGivenTermianls(expectedMapping, transitionOcurrences, timeDurations, mostFreqPaths);
  std::map<uint, std::vector<double>> timeDurationsPerState;
  getTimeDurationsPerStateGivenAncestrals(timeDurations, timeDurationsPerState);
  
  auto nodes = tree_->getAllNodes();
  for (size_t i = 0; i < nodes.size(); i++){
    uint nodeId = tree_->getNodeIndex(nodes[i]);
    if (nodeId == tree_->getRootIndex()){
      continue;
    }
    bool foundPath = true;
    auto father = tree_->getFatherOfNode(tree_->getNode(nodeId));
    uint fatherId = tree_->getNodeIndex(father);
    if (fatherId == tree_->getRootIndex()){
      fatherId = expectedMapping->getRootIndex();
    }
    auto branch = expectedMapping->getEdgeToFather(nodeId);
    double branchLength = branch->getLength();
    size_t fatherState = (size_t)getNodeState(expectedMapping->getNode(fatherId));
    size_t sonState = (size_t)getNodeState(expectedMapping->getNode(nodeId));
    auto mappingStates = MultiStateMappingPath::findExpectedMappingPathForEachNode(fatherState, sonState, transitionOcurrences[nodeId], timeDurationsPerState[nodeId], branchLength, foundPath);
    if (mappingStates.size() > 0){
      mappingStates.push_back(sonState); // this is a dummy transition, just to create the transition of the last state to itself.

    }else{
      if (foundPath){
        continue;

      }else{
        std::cout << "Most frequent path:\n";
        if (mostFreqPaths[nodeId].size() == 0){
          continue;
        }
        for (auto &state : mostFreqPaths[nodeId]){
          mappingStates.push_back(state);
          std::cout << state << ",";
        }
        mappingStates.push_back(sonState);
        std::cout << sonState << std::endl;
      }
      
    }
    std::map<pair<size_t, size_t>, double> occurrencesOfTrnasitionsInExpectedPath;
    for (size_t j = 0; j < mappingStates.size()-1; j++){
      std::pair<size_t,size_t> transition(mappingStates[j], mappingStates[j+1]);
      if (occurrencesOfTrnasitionsInExpectedPath.find(transition) != occurrencesOfTrnasitionsInExpectedPath.end()){
        occurrencesOfTrnasitionsInExpectedPath[transition] += 1;
      }else{
        occurrencesOfTrnasitionsInExpectedPath[transition] = 1;
      }
    }
    auto &timeDurationPerTransitionPerNode = timeDurations[nodeId];
    double sumOfChosenTransitionsTimes = 0;
    std::map<std::pair<size_t, size_t>, double> newBranchLengths;
    for (size_t j = 0; j < mappingStates.size()-1; j++){
      std::pair<size_t,size_t> transition(mappingStates[j], mappingStates[j+1]);
      newBranchLengths[transition] = timeDurationPerTransitionPerNode[transition]/occurrencesOfTrnasitionsInExpectedPath[transition];
      sumOfChosenTransitionsTimes += newBranchLengths[transition];
    }
    // now fragmenting the edge
    double segmentBranchLength;
    //double dwellingTime;
    uint newNodeId;
    double sumOfTransitionsTime = 0;
    std::pair<size_t,size_t> transition;
    for (size_t j = 0; j < mappingStates.size()-2; j++){
      transition = pair<size_t,size_t>(mappingStates[j], mappingStates[j+1]);
      segmentBranchLength = newBranchLengths[transition]/(sumOfChosenTransitionsTimes/branchLength);
      sumOfTransitionsTime += segmentBranchLength;
      auto edge_to_fragment = expectedMapping->getEdgeToFather(nodeId);
      newNodeId = expectedMapping->createNodeOnEdge(expectedMapping->getEdgeIndex(edge_to_fragment), segmentBranchLength);
      (expectedMapping->getNode(newNodeId))->setName("N_dummy_"+ std::to_string(newNodeId)+"-"+ std::to_string(mappingStates[j+1]));


    }
    auto lastTransition = pair<size_t,size_t>(mappingStates[mappingStates.size()-2], mappingStates[mappingStates.size()-1]);
    double estimatedRemained = newBranchLengths[lastTransition]/(sumOfChosenTransitionsTimes/branchLength);
    double truelyRemained = branchLength-sumOfTransitionsTime;
    double epsilon = 1e-6;
    if (std::abs(estimatedRemained - truelyRemained) > epsilon){
      throw Exception("StochasticMapping::findExpectedHistoryTransitionsAndTimeDurationsMultiState(): sum of segments is"+std::to_string(sumOfTransitionsTime)+ " while branch length is "+ std::to_string(branchLength) + "\n");
    }

  }

}

/******************************************************************************/
void StochasticMapping::findTransitionsAndTimeDurationsForBinary(std::shared_ptr<PhyloTree> expectedMapping, std::map<uint, std::vector<double>> &dwellingTimes, std::map<uint, std::vector<double>> &ancestralStatesFrequencies){
  auto nodes = tree_->getAllNodes();
  for (size_t i = 0; i < nodes.size(); i++){
    uint nodeId = tree_->getNodeIndex(nodes[i]);
    if (nodeId == tree_->getRootIndex()){
      continue;
    }
    auto father = tree_->getFatherOfNode(tree_->getNode(nodeId));
    uint fatherId = tree_->getNodeIndex(father);
    if (fatherId == tree_->getRootIndex()){
      fatherId = expectedMapping->getRootIndex();
    }
    size_t fatherState = (size_t)getNodeState(expectedMapping->getNode(fatherId));
    size_t sonState = (size_t)getNodeState(expectedMapping->getNode(nodeId));
    if (fatherState != sonState){
      // we assume that only one change had occurred, and add one node
      auto edge_to_fragment = expectedMapping->getEdgeToFather(nodeId);
      uint newNodeId = expectedMapping->createNodeOnEdge(expectedMapping->getEdgeIndex(edge_to_fragment), dwellingTimes[nodeId][fatherState]);
      (expectedMapping->getNode(newNodeId))->setName("N_dummy_"+ std::to_string(newNodeId)+"-"+ std::to_string(sonState));

    }else{
      // in case we have a branch with same terminal states, we consider two options:
      // 1. x->y->x (two transitions)
      // 2. no transitions if the duration time of the other state was too small.
      size_t otherState = 1-fatherState;
      auto branch = expectedMapping->getEdgeToFather(nodeId);
      double branchLength = branch->getLength();
      if (dwellingTimes[nodeId][otherState] < EPSILON_THRESHOLD * branchLength){
        // most probably a noise, and no tranition had occurred
        continue;
      }
      double fatherStatePosterior = ancestralStatesFrequencies[tree_->getNodeIndex(father)][fatherState];
      double sonStatePosterior = ancestralStatesFrequencies[nodeId][fatherState];
      
      auto edge_to_fragment = expectedMapping->getEdgeToFather(nodeId);
      // dividing the proportions of time according to the posterior probabilities
      double weightFather = fatherStatePosterior/(fatherStatePosterior+sonStatePosterior);
      uint newNodeId = expectedMapping->createNodeOnEdge(expectedMapping->getEdgeIndex(edge_to_fragment), weightFather* dwellingTimes[nodeId][fatherState]);
      (expectedMapping->getNode(newNodeId))->setName("N_dummy_"+ std::to_string(newNodeId)+"-"+ std::to_string(otherState));
      auto second_edge_to_fragment = expectedMapping->getEdgeToFather(nodeId);
      newNodeId = expectedMapping->createNodeOnEdge(expectedMapping->getEdgeIndex(second_edge_to_fragment), dwellingTimes[nodeId][otherState]);
      (expectedMapping->getNode(newNodeId))->setName("N_dummy_"+ std::to_string(newNodeId)+"-"+ std::to_string(sonState));
    }
  }

}


/*******************************************************************************/

shared_ptr<PhyloTree> StochasticMapping::generateAnalyticExpectedMapping(size_t divMethod)
{
  /* Compute the posterior assignment probabilities to internal nodes, based on the fractional probablities computed earlier */
  // const vector<int> states =  tl_->getAlphabetStates();
  vector<int> states = likelihood_->getStateMap().getAlphabetStates();
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
  std::map<uint, std::vector<double>> posteriorProbabilities;
   // because the sum of partial likelihoods (i.e, the fractional probabilities) is in fact the probablity of the data, it is sufficient to standardize the vector of fractional probabilires for each node to obtain the posterior probabilities
  getPosteriorProbabilities(posteriorProbabilities);
  // /* Assign states to internal nodes based on the majority rule over the posterior probabilities */
  shared_ptr<PhyloTree> expectedMapping(make_shared<PhyloTree>(*tree_));
  setExpectedAncestrals(expectedMapping, posteriorProbabilities);
  /* Compute the reward per state per site - expect two entries per site (that is, two entries in total).
  // Let r0 be the reward of state 0 nd r1 the reward of state 1. */
  
  //std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new ConstantRateDistribution());
  //likelihood_->getSubstitutionProcess().getModel(0, 0);
  const std::shared_ptr<const TransitionModel> model = dynamic_pointer_cast<const TransitionModel>(likelihood_->getSubstitutionProcess().getModel(0, 0));

  VVDouble expectedDwellingTimes;
  expectedDwellingTimes.clear();
  expectedDwellingTimes.resize(nodeIds.size(), VDouble(states.size()));
  for (size_t s = 0; s < states.size(); ++s)
  {
    UserAlphabetIndex1 alpha = UserAlphabetIndex1(likelihood_->getData()->getAlphabet());
    alpha.setIndex(states[s], 1); // set the reward of the state as 1 and the reward for the rest of the states as 0
    for (size_t m = 0; m < states.size(); ++m)
    {
      if (m != s)
      {
        alpha.setIndex(states[m], 0); //Note@Laurent (Julien 17/06/20): can you chack my correction there and above? I changed s/m to states[s] and states[m], is that correct?
      }
    }
    DecompositionReward reward(dynamic_cast<const SubstitutionModel*>(model.get()), &alpha); // TO FIX 20.6: this line attempts to delete alpha which doesn't belong to it. cloning it didn't help - get help from Itay / Anat
    shared_ptr<LikelihoodCalculationSingleProcess> rewardLik = make_shared<LikelihoodCalculationSingleProcess>(*likelihood_);
    ProbabilisticRewardMapping mapping(RewardMappingTools::computeRewardVectors(*rewardLik, tree_->getAllNodesIndexes(), reward, false));
    for (size_t n = 0; n < nodeIds.size(); ++n)
    {
      uint nodeId = nodeIds[n];
      if (nodeId != tree_->getRootIndex()) // for any node except to the root
      {
        expectedDwellingTimes[static_cast<size_t>(nodeId)][s] = mapping.getReward(nodeId, 0); //Note@Laurent (Julien 17/06/20): what is nodeId is negative?
      }
    }
  }

  // standardize expected dwelling itmes, if needed, and update the mapping accorgingly
  double sumOfDwellingTimes;
  // bool updateBranch;
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  // for (size_t n = 0; n < nodes.size(); ++n)
  // {
  //   node = nodes[n];
  //   if (node->hasFather()) // for any node except to the root
  //   {
  //     branchLength = node->getDistanceToFather();
  //     sumOfDwellingTimes = 0;
  //     updateBranch = true;
  //     for (size_t s = 0; s < states.size(); ++s)
  //     {
  //       if (expectedDwellingTimes[static_cast<size_t>(node->getId())][s] == 0) //Note@Laurent (Julien 17/06/20): what is nodeId is negative?

  //       {
  //         updateBranch =  false;
  //       }
  //       sumOfDwellingTimes = sumOfDwellingTimes + expectedDwellingTimes[static_cast<size_t>(node->getId())][s]; //Note@Laurent (Julien 17/06/20): what is nodeId is negative?

  //     }

  //     if (branchLength < 0.00001) // branch length is 0 -> no need to update mapping on the branch
  //     {
  //       node->setDistanceToFather(branchLength);
  //       updateBranch = false;
  //     }
  //     else
  //     {
  //       if (sumOfDwellingTimes != branchLength)
  //       {
  //         for (size_t s = 0; s < states.size(); ++s)
  //         {
  //           expectedDwellingTimes[static_cast<size_t>(node->getId())][s] =  branchLength * (expectedDwellingTimes[static_cast<size_t>(node->getId())][s]) / sumOfDwellingTimes;
  //         }
  //       }
  //     }

  //     if (updateBranch)
  //     {
  //       updateBranchByDwellingTimes(node, expectedDwellingTimes[static_cast<size_t>(node->getId())], posteriorProbabilities, divMethod);
  //     }
  //   }
  // }
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;

  // /* free the resources */
  // delete alpha;
  // delete rDist;
  // delete tlModel;
  // delete drtl;

  return expectedMapping;
}

/******************************************************************************/

void StochasticMapping::giveNamesToInternalNodes(PhyloTree& tree)
{
  auto nodes = tree.getAllInnerNodes();
  for (auto& node:nodes)
  {
    if (!node->hasName())
      node->setName("_baseInternal_" + TextTools::toString(tree.getNodeIndex(node)));
  }
}

/******************************************************************************/

void StochasticMapping::setLeafsStates(std::shared_ptr<PhyloTree> mapping)
{
  // auto leafsStates = likelihood_->getData();
  auto leaves = mapping->getAllLeaves();

  for (auto& leaf: leaves)
  {
    string nodeName = leaf->getName();
//    auto lstates = likelihoods_->getNode(mapping->getNodeIndex(leaf))
//   size_t leafState = static_cast<size_t>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
//     //note@Laurent (Julien on 17/06/20): I thing the above line is incorrect, in particulat the use of the getAlphabetStateAsInt function. It is supposed to take as input a state index (size_t) and return the corresponding character state as an integer. Here you give as input to the method already a sequence character (integer). In most cases that will still work as the characters states for resolved characters are usually 0..n, and there corresponding states 0..n. But it will fail for models with gaps (character state -1) and Markov modulated models (character states 0..n, but state index 0..k*n)
//     setNodeState(node, leafState);
//   }
  }
}

/******************************************************************************/

void StochasticMapping::ComputeConditionals()
{
  // some auxiliiary variables
  size_t nbState = likelihood_->getStateMap().getNumberOfModelStates();
  auto flt = likelihood_->getForwardLikelihoodTree(0);
  auto speciesIndices = tree_->getNodeIndexes(tree_->getAllNodes());
  ConditionalProbabilities_.clear();
  ConditionalProbabilities_.resize((speciesIndices.size()));
  for (size_t i = 0; i < speciesIndices.size(); i++){
    auto speciesId = speciesIndices[i];
    ConditionalProbabilities_[speciesId].resize(nbState);
    if (tree_->getRootIndex() == speciesId){
        auto rootLik = numeric::cwise(((flt->getForwardLikelihoodArrayAtRoot())->getTargetValue()).col(0));
        Eigen::RowVectorXd rootFreqsDouble = likelihood_->getRootFreqs()->getTargetValue();
        auto rootFreqs = numeric::cwise(rootFreqsDouble.row(0).transpose());
        auto prodFreqRootLik = rootFreqs * rootLik;
        auto conditionals = prodFreqRootLik/prodFreqRootLik.sum();
        fillRootConditionals(conditionals);

    }else{
      auto& dagIndexes = flt->getDAGNodesIndexes(speciesId);
      if(dagIndexes.size() > 1){
        throw Exception("StochasticMapping::ComputeConditionals(): not implemented for mixture models!");
      }
      for (const auto& index : dagIndexes){
        auto edgeIndex =  flt->getIncomingEdges(index)[0]; // extracting incoming edge from the father
        auto processEdge = flt->getProcessTree()->getEdge(edgeIndex); //getting the edge from which the Pij(t) matrix should be extracted
        auto transitionMatrix = processEdge->getTransitionMatrix()->getTargetValue(); // Pij(t)
        // std::cout << "Transition matrix:" << std::endl;
        // std::cout << transitionMatrix << std::endl;
        // std::cout << "*** *** ***" << std::endl;
        auto LikNodeMat = flt->getForwardLikelihoodArray(index)->getTargetValue(); // getting likelihood calculations for the specific node (speciesId)
        auto sonLik = numeric::cwise(LikNodeMat.col(0)); // getting the likelihood for the first site
        for (size_t fatherState = 0; fatherState < nbState; fatherState++){
          ConditionalProbabilities_[speciesId][fatherState].resize(nbState);
          auto fatherStatePijt = numeric::cwise(transitionMatrix.row(fatherState).transpose()); // getting P(fatherState->j)(t)
          auto fatherSonJoint = fatherStatePijt * sonLik; //  vector of P(fatherState->j)(t) * Lik(j) for each son state j
          auto conditionals = fatherSonJoint/fatherSonJoint.sum(); 
          for (size_t sonState = 0; sonState < nbState; sonState ++){
            auto conditional = ExtendedFloat(conditionals.float_part()(sonState), conditionals.exponent_part());
            ConditionalProbabilities_[speciesId][fatherState][sonState] = ExtendedFloat::convert(conditional);
          }
        }
      }      
    }
  }

}


/******************************************************************************/
void StochasticMapping::fillRootConditionals(ExtendedFloatArrayXd &conditionals){
  size_t nbState = likelihood_-> getStateMap().getNumberOfModelStates();
  for (size_t fatherState = 0; fatherState < nbState; fatherState++){
    ConditionalProbabilities_[tree_->getRootIndex()][fatherState].resize(nbState);
    for (size_t sonState = 0; sonState < nbState; sonState++){
      auto conditional = ExtendedFloat(conditionals.float_part()(sonState), conditionals.exponent_part());
      ConditionalProbabilities_[tree_->getRootIndex()][fatherState][sonState] = ExtendedFloat::convert(conditional);
    }
  }
}
/******************************************************************************/
void StochasticMapping::getPosteriorProbabilities(std::map<uint, std::vector<double>> &posteriorProbs){
  // some auxiliiary variables
  size_t nbState = likelihood_->getStateMap().getNumberOfModelStates();
  auto flt = likelihood_->getForwardLikelihoodTree(0);
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
  for (size_t i = 0; i < nodeIds.size(); i++){
    auto nodeId = nodeIds[i];
    posteriorProbs[nodeId];
    posteriorProbs[nodeId].resize(nbState);
    if (tree_->getRootIndex() == nodeId){
        auto rootLik = numeric::cwise(((flt->getForwardLikelihoodArrayAtRoot())->getTargetValue()).col(0));
        auto rootPosterior = rootLik/rootLik.sum();
        for (size_t s = 0; s < nbState; s++){
          auto efPosterior = ExtendedFloat(rootPosterior.float_part()(s), rootPosterior.exponent_part());
          posteriorProbs[nodeId][s] = ExtendedFloat::convert(efPosterior);
        }

    }else{
      auto& dagIndexes = flt->getDAGNodesIndexes(nodeId);
      if(dagIndexes.size() > 1){
        throw Exception("StochasticMapping::getPosteriorProbabilities(): not implemented for mixture models!");
      }
      for (const auto& index : dagIndexes){
        auto LikNodeMat = flt->getForwardLikelihoodArray(index)->getTargetValue(); // getting likelihood calculations for the specific node (speciesId)
        auto nodeLik = numeric::cwise(LikNodeMat.col(0)); // getting the likelihood for the first site
        auto nodePosterior = nodeLik/nodeLik.sum(); 
        for (size_t s = 0; s < nbState; s ++){
          auto efPosteriorProb = ExtendedFloat(nodePosterior.float_part()(s), nodePosterior.exponent_part());
          posteriorProbs[nodeId][s] = ExtendedFloat::convert(efPosteriorProb);
        }
      }      
    }
  }
  return;

}

/******************************************************************************/

void StochasticMapping::computeStatesFrequencies(std::map<uint, std::vector<double>> &ancestralStatesFreqs)
{
  size_t nbStates = likelihood_-> getStateMap().getNumberOfModelStates();
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
  for (size_t i = 0; i < nodeIds.size(); i++){
    ancestralStatesFreqs[nodeIds[i]];
    ancestralStatesFreqs[nodeIds[i]].resize(nbStates);
    std::fill(ancestralStatesFreqs[nodeIds[i]].begin(), ancestralStatesFreqs[nodeIds[i]].end(), 0);
    if (tree_->isLeaf(tree_->getNode(nodeIds[i]))){
      // for leaves there is a frequency of 1 in one of the states. It does
      // not matter which mapping to choose. I choose arbitrary mapping 0.
      ancestralStatesFreqs[nodeIds[i]][ancetralStates_[nodeIds[i]][0]] = 1.0;
    }else{
      for (size_t j = 0; j < numOfMappings_; j++){
        ancestralStatesFreqs[nodeIds[i]][ancetralStates_[nodeIds[i]][j]]++;
      
      }
      for (size_t k = 0; k < nbStates; k++){
        ancestralStatesFreqs[nodeIds[i]][k] /= static_cast<double>(numOfMappings_);
      }

    }


  }

  

  // // some auxiliiary variables
  // size_t statesNum = tl_->getNumberOfStates();
  // const SiteContainer* leafsStates = tl_->getData();
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree_);
  // vector<Node*> nodes = ttree->getNodes();

  // // compute the node assignment probabilities based on their frequency in the mappings
  // for (size_t i = 0; i < nodes.size(); ++i)
  // {
  //   Node* node = nodes[i];
  //   int nodeId = node->getId();
  //   string nodeName = node->getName();
  //   // in leafs - don't iterate to save time, as the frequency of a state is either 0 or 1 based on the known character data
  //   if (node->isLeaf())
  //   {
  //     size_t leafState = static_cast<int>(tl_->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
  //     for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
  //     {
  //       if (nodeState != leafState)
  //       {
  //         ancestralStatesFrequencies[nodeId][nodeState] = 0;
  //       }
  //       else
  //       {
  //         ancestralStatesFrequencies[nodeId][nodeState] = 1;
  //       }
  //     }
  //   }
  //   else
  //   {
  //     // else, go over all the mappings and collect the number of states assignment per state
  //     fill(ancestralStatesFrequencies[nodeId].begin(), ancestralStatesFrequencies[nodeId].end(), 0); // reset all the values to 0
  //     for (size_t h = 0; h < mappings.size(); ++h)
  //     {
  //       Node* nodeInMapping = dynamic_cast<TreeTemplate<Node>*>(mappings[h])->getNode(nodeName);
  //       ancestralStatesFrequencies[nodeId][static_cast<size_t>(getNodeState(nodeInMapping))]++; //Note@Laurent (Julien 17/06/20): assuming node state is positive, is that so?
  //     }
  //     // now divide the vector entries by the number of mappings
  //     for (size_t nodeState = 0; nodeState < statesNum; ++nodeState)
  //     {
  //       ancestralStatesFrequencies[nodeId][nodeState] = ancestralStatesFrequencies[nodeId][nodeState] / static_cast<int>(mappings.size());
  //     }
  //   }
  // }
}


/******************************************************************************/

size_t StochasticMapping::sampleState(const VDouble& distibution)
{
  size_t state = 0;        // the default state is 0
  double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);

  for (size_t i = 0; i < distibution.size(); ++i)
  {
    prob -= distibution[i];
    if (prob < 0)  // if the the sampled probability is smaller than the probability to choose state i -> set state to be i
    {
      state = i;
      break;
    }
  }
  return state;
}

/******************************************************************************/

void StochasticMapping::sampleAncestrals(size_t mappingIndex)
{

  // iterate over the nodes in preorder, and sampling the ancestral states
  sampleAncestralsRecursively(tree_->getRootIndex(), mappingIndex);

}

/******************************************************************************/
void StochasticMapping::sampleAncestralsRecursively(uint nodeId, size_t mappingIndex, uint *fatherIndex){
  auto node = tree_->getNode(nodeId);
  if (nodeId == tree_->getRootIndex()){
    size_t rootState = sampleState(ConditionalProbabilities_[nodeId][0]);
    ancetralStates_[nodeId][mappingIndex] = rootState;
  }else{
    size_t fatherState = ancetralStates_[*fatherIndex][mappingIndex];
    size_t sonState = sampleState(ConditionalProbabilities_[nodeId][fatherState]);
    ancetralStates_[nodeId][mappingIndex] = sonState;
  }
  if (!(tree_->isLeaf(nodeId))){
    auto nodeIdsSons = tree_->getSons(nodeId);
    for (size_t i = 0; i < nodeIdsSons.size(); i++){
      sampleAncestralsRecursively(nodeIdsSons[i], mappingIndex, &nodeId);
    }  
  }
}

/******************************************************************************/


bool StochasticMapping::sampleMutationsGivenAncestrals(size_t mappingIndex, vector<uint>* failedNodes)
{
  bool allSuccess = true;
  auto nodeIndices = tree_->getNodeIndexes(tree_->getAllNodes());
  for (size_t i = 0; i < nodeIndices.size(); i++){
    if (tree_->isLeaf(tree_->getNode(nodeIndices[i]))){
      continue;
    }
    auto father = nodeIndices[i];
    auto sons = tree_->getSons(father);
    for (size_t j = 0; j < sons.size(); j++){
      // 1. change to bool return type
      // 2. If false -> the simulation has failed.
      // 3. If the simulation has failed -> resample the ancestral states.
      // 4. Once the ancestral states are resampled -> call again to sampleMutationsGivenAncestrals()
       
      bool success = sampleMutationsGivenAncestralsPerBranch(father, sons[j], mappingIndex, mappings_[sons[j]], numOfMappingTrials_);
      if (!success){
        allSuccess = false;
        if (failedNodes){
          failedNodes->push_back(sons[j]);          
          continue;
        }else{
          return allSuccess;
        }      
      }
    }    
  }
  return allSuccess;
}

/******************************************************************************/

void StochasticMapping::updateBranchMapping(PhyloNode* son, const MutationPath& branchMapping)
{
  // const vector<size_t> states = branchMapping.getStates();
  // const VDouble times = branchMapping.getTimes();
  // Node* curNode = son;
  // Node* nextNode;
  // int eventsNum = static_cast<int>(branchMapping.getNumberOfEvents());

  // if (eventsNum == 0) // if there are no events >-return nothing
  //   return;
  // else
  // {
  //   for (int i = eventsNum - 1; i > -1; --i) // add a new node to represent the transition
  //   {
  //     nodesCounter_ = nodesCounter_ + 1;
  //     const string name = "_mappingInternal" + TextTools::toString(nodesCounter_) + "_";
  //     nextNode = new Node(static_cast<int>(nodesCounter_), name);
  //     setNodeState(nextNode, states[i]);
  //     nextNode->setDistanceToFather(times[i]);

  //     // set the father to no longer be the father of curNode
  //     Node* originalFather = curNode->getFather();
  //     originalFather->removeSon(curNode); // also removes originalFather as the father of curNode
  //     // set nextNode to be the new father of curNode
  //     curNode->setFather(nextNode); // also adds curNode to the sons of nextNode
  //     // set curNode's original father ot be the father of nextNode
  //     nextNode->setFather(originalFather); // also adds nextNode to the sons of originalFather - or not? make sure this is the father at all times
  //     curNode = nextNode;
  //   }
  //   return;
  // }
}

/******************************************************************************/

bool StochasticMapping::sampleMutationsGivenAncestralsPerBranch(uint father, uint son, size_t mappingIndex, vector<MutationPath> &mappings, size_t maxIterNum)
{
  
  size_t fatherState = ancetralStates_[father][mappingIndex];
  size_t sonState = ancetralStates_[son][mappingIndex];

  auto branchPtr = tree_->getIncomingEdges(tree_->getNode(son))[0];
  auto branchLength = branchPtr->getLength();

  /* simulate mapping on a branch until you manage to finish at the son's state */
  bool success = sampleEvolutionaryPathForBranch(sonState, fatherState, father, son, branchLength, mappingIndex, mappings, maxIterNum); //TODO put the following lines (inside the for loop) into the new function
  if (!success){
    std::cout << "Mapping failure! " << "Mapping index: " << mappingIndex;
    std::cout << ", nodeId: " << son << ", fatherState: " << fatherState << ", sonState: " << sonState << ", branchLength: " << branchLength;
    //std::cout << ", probability of son given father: " << ConditionalProbabilities_[son][fatherState][sonState];
    // if (!(father == tree_->getRootIndex())){
    //   auto grandFather = tree_->getFatherOfNode (tree_->getNode(father));
    //   uint grandFatherId = tree_->getNodeIndex(grandFather);
    //   size_t grandFatherState = ancetralStates_[grandFatherId][mappingIndex];
    //   std::cout << ", father id: " << father << ", grand father id: " << grandFatherId << ", grandFather state: " << grandFatherState;
    //   std::cout << ", probability of father given grandFather: " << ConditionalProbabilities_[father][grandFatherState][fatherState] << std::endl;

    // }

  }
  return success;
}

/******************************************************************************/
void StochasticMapping::getDewellingTimesUnderEachStatePerMapping(vector<double> *dwellingTimes, size_t mappingIndex){
  auto rootId = tree_->getRootIndex();
  auto sons = tree_->getSons(rootId);
  for (size_t i = 0; i < sons.size(); i++){
    getDewellingTimesUnderEachStatePerMappingRecursively(sons[i], ancetralStates_[rootId][mappingIndex], dwellingTimes, mappingIndex, 0);
  }

}
/******************************************************************************/
void StochasticMapping::getDewellingTimesUnderEachStatePerNode(std::map<uint, vector<double>> *dwellingTimes){
  // initialize map
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
  for (size_t i = 0; i < nodeIds.size(); i++){
    (*dwellingTimes)[nodeIds[i]];
    (*dwellingTimes)[nodeIds[i]].resize(likelihood_-> getStateMap().getNumberOfModelStates());
  }
  
  for (size_t m = 0; m < numOfMappings_; m++){
    auto rootId = tree_->getRootIndex();
    auto sons = tree_->getSons(rootId);
    for (size_t i = 0; i < sons.size(); i++){
      getDewellingTimesUnderEachStatePerMappingRecursively(sons[i], ancetralStates_[rootId][m], 0, m, dwellingTimes);
    }
  }
  auto it = dwellingTimes->begin();
  while (it != dwellingTimes->end()){
    (*dwellingTimes)[it->first] /= static_cast<double>(numOfMappings_);
    it ++;
  }
}
// /*****************************************************************************/
bool StochasticMapping::sampleEvolutionaryPathForBranch(size_t sonState, size_t fatherState, uint father, uint son, double branchLength, size_t mappingIndex, vector<MutationPath>& mappings, size_t maxIterNum, bool replace){
  bool success = true;
  auto alphabet = likelihood_->getData()->getAlphabet();

  auto model = dynamic_pointer_cast<const SubstitutionModel>(likelihood_->getSubstitutionProcess().getModel(father, 0)); // father or son??? Should be a father, because the models start at particular nodes, and I should get the model of the preceeding branch
  size_t modelIndex;
  if (father == tree_->getRootIndex()){
    modelIndex = 1; // I guess the model index should be 1 for the root. Is it true?
  }else{
    modelIndex = likelihood_->getSubstitutionProcess().getModelNumberForNode(father);
  }  
  for (size_t i = 0; i < maxIterNum; i++){
    double disFromNode = 0.0;
    size_t curState = fatherState;
    MutationPath tryMapping(alphabet, fatherState, branchLength);
  
    double timeTillChange;
    // if the father's state is not the same as the son's state -> use the correction corresponding to equation (11) in the paper
    if (fatherState != sonState)
    {   // sample timeTillChange conditional on it being smaller than branchLength
      double u = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
      double waitingTimeParam = -1 * model->Qij(fatherState, fatherState); // get the parameter for the exoponential distribution to draw the waiting time from
      double tmp = u * (1.0 - exp(branchLength * -waitingTimeParam));
      timeTillChange =  -log(1.0 - tmp) / waitingTimeParam;
      assert (timeTillChange < branchLength);
    }
    else
    {
      timeTillChange = RandomTools::randExponential(-1. / model->Qij(fatherState, fatherState));// draw the time until a transition from exponential distribution with the rate of leaving fatherState
    }

    while (disFromNode + timeTillChange < branchLength)  // a jump occured but not passed the whole branch ->
    {

      curState = giveRandomState(curState, modelIndex);  // draw the state to transition to after from initial state curState based on the relative tranistion rates distribution (see MutationProcess.cpp line 50)
      tryMapping.addEvent(curState, timeTillChange);                                        // add the current state and time to branch history
      disFromNode += timeTillChange;
      timeTillChange = RandomTools::randExponential(-1. / model->Qij(curState, curState));        // draw the time until a transition from exponential distribution with the rate of leaving curState
                                          
    }
    // the last jump passed the length of the branch -> finish the simulation and check if it's sucessfull (i.e, mapping is finished at the son's state)
    if (curState != sonState) // if the simulation failed, try again
    {
      continue;
    }
    else                      // if the simulation was sucessfully, add it to the build mapping
    {
      mappings.push_back(tryMapping);
      if (mappings.size() != mappingIndex+1){
        throw Exception("StochasticMapping::sampleEvolutionaryPathForBranch(): Unsuccessful mapping was sampled!!!");
      }
      // if (replace){
      //   mappings_[son][mappingIndex] = tryMapping;

      // }else{
      //   mappings_[son].push_back(tryMapping);
      //   // *** debug ***//
      //   if (mappings_[son].size() != mappingIndex+1){
      //     throw Exception ("StochasticMapping::sampleMutationsGivenAncestralsPerBranch: Something went wrong when filling mappings_ object!");
      //   }

      // }

      return success;
    }
  }
  return false;
}
/******************************************************************************/
bool StochasticMapping::isAccounted(uint nodeId, size_t mappingIndex){
  bool accounted = true;
  if (notRepresentedNodes_.find(nodeId) != notRepresentedNodes_.end()){
    auto &failedMappings = notRepresentedNodes_[nodeId];
    if (std::find(failedMappings.begin(), failedMappings.end(), mappingIndex) != failedMappings.end()){
      accounted = false;

    }
  }
  return accounted;
}
/******************************************************************************/
void StochasticMapping::assignDewellingTimesUnderEachStatePerMappingPerBranch(uint nodeId, size_t initialState, vector<double> &dwellingTimes, MutationPath &mutationPath){
  vector<size_t> states;
  states = mutationPath.getStates();
  auto times = mutationPath.getTimes();
  auto branchPtr = tree_->getIncomingEdges(tree_->getNode(nodeId))[0];
  auto branchLength = branchPtr->getLength();
  double spentTimeOnBranch = 0;
  for (size_t i = 0; i < states.size(); i++){
    spentTimeOnBranch += times[i];
    if (i == 0){
      // dwelling time of the intial state (that lasts from the previous branch)
      dwellingTimes[initialState] += times[i];
    }else{
      dwellingTimes[states[i-1]] += times[i];
    }
  }
  double timeSpentUnderLastState = branchLength-spentTimeOnBranch;
  if (states.size() == 0){
    dwellingTimes[initialState] += timeSpentUnderLastState;
  }else{
    dwellingTimes[states[states.size()-1]] += timeSpentUnderLastState;
  }

}



/******************************************************************************/
void StochasticMapping::getDewellingTimesUnderEachStatePerMappingRecursively(uint nodeId, size_t initialState, vector<double> *dwellingTimesStates, size_t mappingIndex, std::map<uint, std::vector<double>> *dwellingTimesPerNode){
  auto mutationPath = mappings_[nodeId][mappingIndex];
  auto & dwellingTimes = (dwellingTimesPerNode) ? (*dwellingTimesPerNode)[nodeId] : *dwellingTimesStates;
  //if ((notRepresentedNodes_.find(nodeId) == notRepresentedNodes_.end()) && ()){
  bool accountedBranch = isAccounted(nodeId, mappingIndex);
  vector<size_t> states;
  if (accountedBranch){
    states = mutationPath.getStates();
    auto times = mutationPath.getTimes();
    auto branchPtr = tree_->getIncomingEdges(tree_->getNode(nodeId))[0];
    auto branchLength = branchPtr->getLength();
    double spentTimeOnBranch = 0;
    for (size_t i = 0; i < states.size(); i++){
      spentTimeOnBranch += times[i];
      if (i == 0){
        // dwelling time of the intial state (that lasts from the previous branch)
        dwellingTimes[initialState] += times[i];
      }else{
        dwellingTimes[states[i-1]] += times[i];
      }
    }
    double timeSpentUnderLastState = branchLength-spentTimeOnBranch;
    if (states.size() == 0){
      dwellingTimes[initialState] += timeSpentUnderLastState;
    }else{
      dwellingTimes[states[states.size()-1]] += timeSpentUnderLastState;
    }

  }else{
    if (tree_->isLeaf(nodeId)){
      return;

    }

  }

  if (!(tree_->isLeaf(tree_->getNode(nodeId)))){
    auto sons = tree_->getSons(nodeId);
    for (size_t n = 0; n < sons.size(); n++){
      size_t initialStateForSon;
      if (accountedBranch){
        if (states.size() > 0){
          initialStateForSon = states[states.size()-1];
        }else{
          initialStateForSon = initialState;
        }
      }else{
        initialStateForSon = ancetralStates_[nodeId][mappingIndex];
      }
      getDewellingTimesUnderEachStatePerMappingRecursively(sons[n], initialStateForSon, dwellingTimesStates, mappingIndex, dwellingTimesPerNode);
      
    }
  }
}
/******************************************************************************/
void StochasticMapping::getNumOfOcuurencesForEachTransitionPerMapping(size_t mappingIndex, std::map<uint, std::map<std::pair<size_t, size_t>, double>> &transitionOcurrences){
  //size_t nbState = likelihood_->getStateMap().getNumberOfModelStates();

  // fill the map with the corresponding occurences (starting from the root)
  uint rootId = tree_->getRootIndex();
  size_t rootState = ancetralStates_[rootId][mappingIndex];
  auto sons = tree_->getSons(rootId);
  for (size_t n = 0; n < sons.size(); n++){
    getNumOfOcuurencesForEachTransitionPerMappingRecursively(sons[n], rootState, mappingIndex, transitionOcurrences);
  }
}
/******************************************************************************/
void StochasticMapping::getExpectedNumberOfTransitionsPerGivenTermianls(std::shared_ptr<PhyloTree> expectedTree, std::map<uint, std::map<pair<size_t, size_t>, double>> &transitionOcurrences, std::map<uint, std::map<pair<size_t, size_t>, double>> &timeDurations, std::unordered_map<uint, vector<size_t>> &mostFreqPaths){
  vector<uint> nodeIndexes = tree_->getNodeIndexes(tree_->getAllNodes());
  for (size_t i = 0; i < nodeIndexes.size(); i++){
    if (nodeIndexes[i] == tree_->getRootIndex()){
      continue;
    }
    auto node = expectedTree->getNode(nodeIndexes[i]);
    auto fatherNode = expectedTree->getFatherOfNode(node);
    uint fatherId;
    if (expectedTree->getRootIndex() == expectedTree->getNodeIndex(fatherNode)){
      fatherId = tree_->getRootIndex();
    }else{
      fatherId = expectedTree->getNodeIndex(fatherNode);
    }
    size_t startState = static_cast<size_t>(getNodeState(fatherNode));
    size_t endState = static_cast<size_t>(getNodeState(node));
    mostFreqPaths[nodeIndexes[i]];
    getExpectedNumberOfTransitionsPerBranchGivenTerminals(nodeIndexes[i], fatherId, startState, endState, transitionOcurrences[nodeIndexes[i]], timeDurations[nodeIndexes[i]], mostFreqPaths[nodeIndexes[i]], mappings_[nodeIndexes[i]]);
    

  }

}
/******************************************************************************/
// std::vector<size_t> StochasticMapping::decimalToBinaryPowers(int decimalNumber) {
//     std::vector<size_t> powers;
//     int power = 0;
//     while (decimalNumber > 0) {
//         if (decimalNumber % 2 == 1) {
//             powers.push_back((size_t)(power));
//         }

//         decimalNumber /= 2;
//         ++power;
//     }

//     return powers;
// }
/******************************************************************************************************/
// std::map<size_t, vector<size_t>> StochasticMapping::createEdges(std::map<size_t, double> &vertices, std::map<std::pair<size_t, size_t>, double> &transitions){
//   std::map<size_t, vector<size_t>> edges;
//   auto it = transitions.begin();
//   while (it != transitions.end()){
//     auto transition = it->first;
//     if ((transitions[transition] <= 0) || (transition.first == transition.second)){
//       it ++;
//       continue;

//     }
//     size_t outgoing = transition.first;
//     if (vertices[outgoing] < EPSILON_THRESHOLD){
//         it ++;
//         continue;
//     }
//     size_t incoming = transition.second;
//     if (vertices[incoming] < EPSILON_THRESHOLD){
//         it ++;
//         continue;
//     }
//     if (edges.find(outgoing) == edges.end()){
//       edges[outgoing];
      
//     }
//     edges[outgoing].push_back(incoming);
//     it ++;
//   }
//   return edges;
// }
/******************************************************************************/
// void StochasticMapping::findBestPath(std::pair<size_t,size_t> &bestCandidatePathId, std::map<std::pair<size_t, size_t>, double> &paths, size_t desiredPathId, std::map<size_t, vector<size_t>> &edges, std::map<std::pair<size_t, size_t>, double> &transitions, size_t end, bool &foundPath){
//   std::pair<size_t,size_t> candidatePathId;
//   double weightBest = 0;
//   auto itPath = paths.begin();
//   vector<size_t> neighbors;
//   while (itPath != paths.end()){
//     size_t pathId = (itPath->first).first;
//     if (pathId != desiredPathId){
//       itPath++;
//       continue;
//     }
//     size_t lastMid = (itPath->first).second;
//     neighbors = edges[lastMid];
//     double weight;
    
//     if (std::find(neighbors.begin(), neighbors.end(), end) != neighbors.end()){
//       candidatePathId = itPath->first;
//       std::pair<size_t, size_t> lastEdge(lastMid, end);
//       weight = paths[candidatePathId] + transitions[lastEdge];
//       foundPath = true;
//       if (weight > weightBest){
//         bestCandidatePathId = itPath->first;
//         weightBest = weight;

//       }
//     }
//     itPath++;
//   }
// }
/******************************************************************************/
// void StochasticMapping::reconstructBestPath(std::vector<size_t> &bestPath, size_t lengthOfPath, std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> &pathReconstruction, std::pair<size_t,size_t> bestCandidatePathId, size_t start, size_t end){
//   bestPath.resize(lengthOfPath);
//   std::pair<size_t, size_t> fatherPath;
//   size_t currState;
//   for (size_t k = 1; k <= lengthOfPath; k++){
//     size_t reverseIndex = lengthOfPath-k;
//     if (k == 1){
//       bestPath[reverseIndex] = end;
//     }else if (k == lengthOfPath){
//       bestPath[reverseIndex] = start;
//     }else{
//       if (k  == 2){
//         fatherPath = bestCandidatePathId;
//       }else{
//         fatherPath = pathReconstruction[fatherPath];
//       }
//       currState = fatherPath.second;
//       bestPath[reverseIndex] = currState;

//     }
//   }

// }

/******************************************************************************/
// vector<size_t> StochasticMapping::findExpectedMappingPathForEachNode(size_t start, size_t end, std::map<std::pair<size_t, size_t>, double> &transitions, vector<double> &dwellingTimes, double totalDurationTime){
//   std::map<size_t, double> relativeTimeDuration;
//   std::vector<size_t> bestPath;
//   size_t desiredPathId = 0;
  
//   size_t numOfNotAdded = 0;
//   for (size_t i = 0; i < dwellingTimes.size(); i++){
//     relativeTimeDuration[i] = dwellingTimes[i]/totalDurationTime;
//     if (relativeTimeDuration[i] < EPSILON_THRESHOLD){
//       numOfNotAdded ++;
//     }else{
//       if (i != end){
//         desiredPathId += std::pow(2, i);
//       }else{
//         if (start == end){
//           desiredPathId += std::pow(2, i);
//         }
//       }
//     }
//   }
//   auto edges = createEdges(relativeTimeDuration, transitions);
//   if (edges.size() == 0){
//     return bestPath;
//   }
//   size_t pathLength = dwellingTimes.size()-numOfNotAdded-1; // we don't include the final end node.
  
//   std::map<std::pair<size_t, size_t>, double> paths;
//   std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> pathReconstruction;
//   auto neighbors = edges[start];
//   if (neighbors.size() == 0){  
//     //bestPath.push_back(end);
//     return bestPath;
//   }
//   if ((pathLength == 1) && (start != end)){
//     if (std::find(neighbors.begin(), neighbors.end(), end) != neighbors.end()){
//       bestPath.push_back(start);
//       bestPath.push_back(end);
//       return bestPath;
//     }
//     throw Exception("StochasticMapping::findExpectedMappingPathForEachNode: No such path!");
//   }
//   for (size_t i = 0; i < neighbors.size(); i++){
//     if (neighbors[i] == end){
//       continue;
//     }
//     size_t pathId = std::pow(2, start) + std::pow(2, neighbors[i]);
//     size_t fatherPathId = 0;
//     size_t endOfFather = start;
//     std::pair<size_t, size_t> fatherPathWithEnd(fatherPathId, endOfFather);

    
//     std::pair<size_t, size_t> pathWithEnd(pathId, neighbors[i]);
//     std::pair<size_t, size_t> edge(start, neighbors[i]);
//     paths[pathWithEnd] = transitions[edge];
//     pathReconstruction[pathWithEnd] = fatherPathWithEnd;
//   }
//   if (start == end){
//     pathLength -= 1; // path length excludes start and an additional neighbor but not the dest node, since start and end are the same.

//   }else{
//     pathLength -= 2; // path length excludes start and an additional neighbor in addition to the dest node

//   }
//   vector<std::pair<size_t, size_t>> pathsIdsOfInitialLength;
//   auto it = paths.begin();
//   while (it != paths.end()){
//     size_t id = (it->first).first;
//     if (id != 0){
//         pathsIdsOfInitialLength.push_back(it->first);
//     }

//     it ++;
//   }
//   auto pathsIdsOfGivenLength = pathsIdsOfInitialLength;
  
//   for (size_t i = 0; i < pathLength; i++){
//     vector<std::pair<size_t, size_t>> pathsIdsOfCurrentLength;
//     for (size_t j = 0; j < pathsIdsOfGivenLength.size(); j++){
//       size_t fatherIdPath = pathsIdsOfGivenLength[j].first;
//       size_t newStart = pathsIdsOfGivenLength[j].second;
//       neighbors = edges[newStart];
//       for (size_t k = 0; k< neighbors.size(); k++){
//         if (neighbors[k] == end){
//           continue;
//         }
//         auto used = decimalToBinaryPowers(fatherIdPath);
//         if (std::find(used.begin(), used.end(), neighbors[k]) != used.end()){
//           continue;
//         }
//         size_t currentPathId = fatherIdPath + std::pow(2, neighbors[k]);

//         std::pair<size_t, size_t> currentPathWithEnd(currentPathId, neighbors[k]);
//         std::pair<size_t, size_t> edge(newStart, neighbors[k]);
//         double weight = paths[pathsIdsOfGivenLength[j]] + transitions[edge];
//         if (paths.find(currentPathWithEnd) != paths.end()){
//           if (weight <= paths[currentPathWithEnd]){
//             continue;
//           }
//         }else{
//           paths[currentPathWithEnd] = weight;
//           pathReconstruction[currentPathWithEnd] = pathsIdsOfGivenLength[j];
//           pathsIdsOfCurrentLength.push_back(currentPathWithEnd);
//         }
        
//       }

//     }
//     pathsIdsOfGivenLength = pathsIdsOfCurrentLength;
//   }
//   // Find the best path
  
//   bool foundPath = false;
//   std::pair<size_t,size_t> bestCandidatePathId; 
//   findBestPath(bestCandidatePathId, paths, desiredPathId, edges, transitions, end, foundPath);

//   // reconstruct the best path
//   size_t lengthOfPath = dwellingTimes.size()-numOfNotAdded;
//   if (start == end){
//     lengthOfPath ++;
//   }
//   if (!foundPath){
//     throw Exception("Path does not exist!");
//   }
//   reconstructBestPath(bestPath, lengthOfPath, pathReconstruction, bestCandidatePathId, start, end);
//   return bestPath;


// }

/******************************************************************************/
// This fucntion is needed for the multi-state heuristic approach 
void StochasticMapping::getExpectedNumberOfTransitionsPerBranchGivenTerminals(uint nodeId, uint fatherId, size_t startState, size_t endState, std::map<pair<size_t, size_t>, double> &transitionOcurrences, std::map<pair<size_t, size_t>, double> &timeDurations, vector<size_t> &mostFreqPath, vector<MutationPath> &mappings){
  size_t counter = 0;
  auto branch = tree_->getEdgeToFather(nodeId);
  auto branchLength = branch->getLength();
  std::unordered_map<string, int> mappingsFrequencies;
  for (size_t i = 0; i < numOfMappings_; i++){
    if (!(isAccounted(nodeId, i))){
      continue;
    }
    double remainedTime = 0;
    size_t fatherState = ancetralStates_[fatherId][i];
    size_t sonState = ancetralStates_[nodeId][i];

    if ((fatherState == startState) && (sonState == endState)){
      counter ++;
      auto mutationPath = mappings[i];
      vector<size_t> states = mutationPath.getStates();
      string states_str = "";

      for (auto &state : states){
        states_str += std::to_string(state);
      }
      if (mappingsFrequencies.find(states_str) != mappingsFrequencies.end()){
        mappingsFrequencies[states_str] += 1;
      }else{
        mappingsFrequencies[states_str] = 1;

      }
      auto times = mutationPath.getTimes();
      if (states.size() == 0){
        std::pair<size_t, size_t> noTransition(fatherState, sonState);


        if (timeDurations.find(noTransition) == timeDurations.end()){
          timeDurations[noTransition] = branchLength;
        }else{
          timeDurations[noTransition] += branchLength;
        }
        continue;
      }
      std::pair<size_t, size_t> firstTransitionOnBranch(fatherState, states[0]);
      // updating number of occurrences

      if (transitionOcurrences.find(firstTransitionOnBranch) == transitionOcurrences.end()){
        transitionOcurrences[firstTransitionOnBranch] = 1;
      }else{
        transitionOcurrences[firstTransitionOnBranch] += 1;

      }
      // updating time durations
      double timeSoFar = 0;
      timeSoFar += times[0];

      if (timeDurations.find(firstTransitionOnBranch) == timeDurations.end()){
        timeDurations[firstTransitionOnBranch] = times[0];
      }else{
        timeDurations[firstTransitionOnBranch] += times[0];
      }


      for (size_t j = 0; j < states.size()-1; j++){
        std::pair<size_t, size_t> transition(states[j], states[j+1]);
        if (transitionOcurrences.find(transition) == transitionOcurrences.end()){
          transitionOcurrences[transition] = 1;
          timeDurations[transition] = times[j+1];
        }else{
          transitionOcurrences[transition] += 1;
          timeDurations[transition] += times[j+1];

        }
        timeSoFar += times[j+1];
        
      }
      std::pair<size_t, size_t> selfTransition(states[states.size()-1], states[states.size()-1]);
      remainedTime = branchLength - timeSoFar;
      if (timeDurations.find(selfTransition) == timeDurations.end()){
        timeDurations[selfTransition] = remainedTime;
      }else{
        timeDurations[selfTransition] += remainedTime;

      }

    }

  }
  auto it = transitionOcurrences.begin();
  while (it != transitionOcurrences.end()){
    transitionOcurrences[it->first] /= static_cast<double>(counter);
    it ++;
  }
  auto itTime = timeDurations.begin();
  // just for debug!!!
  double sumOfTimeDuration = 0;
  //
  while (itTime != timeDurations.end()){
    timeDurations[itTime->first] /= static_cast<double>(counter);
    sumOfTimeDuration += timeDurations[itTime->first];
    itTime ++;
  }
  if (sumOfTimeDuration <= 0){
    std::cout << "problem is here!" << std::endl;
  }
  auto itPath = mappingsFrequencies.begin();
  string mostFrequent = mappingsFrequencies.begin()->first;
  int maxFreq = 0;

  while (itPath != mappingsFrequencies.end()){
    if (itPath->second > maxFreq){
      mostFrequent = itPath->first;
      maxFreq = itPath->second;
    }
    itPath ++;
  }

  if (mostFrequent != ""){
    mostFreqPath.push_back(startState);
    stringToVector(mostFrequent, mostFreqPath);
  }
  
}
/******************************************************************************/
void StochasticMapping::stringToVector(const std::string& str, vector<size_t> &res) {
    for (char c : str) {
        if (c >= '0' && c <= '9') {
            res.push_back(c - '0');
        }
    }
}

/******************************************************************************/
void StochasticMapping::getNumOfOcuurencesForEachTransitionPerMappingRecursively(uint nodeId, size_t initialState, size_t mappingIndex, std::map<uint, std::map<pair<size_t, size_t>, double>> &transitionOcurrences){
  auto mutationPath = mappings_[nodeId][mappingIndex];
  bool accountedBranch = isAccounted(nodeId, mappingIndex);
  if (tree_->isLeaf(nodeId) && (!accountedBranch)){
    return;
  }
  vector<size_t> states;
  if (accountedBranch){
    states = mutationPath.getStates();
    if (states.size() > 0){ // the size is 0 in case no transitions have occured
      std::pair<size_t, size_t> firstTransitionOnBranch(initialState, states[0]);
      if (transitionOcurrences.find(nodeId) == transitionOcurrences.end()){
        transitionOcurrences[nodeId][firstTransitionOnBranch] = 0;
      }else{
        if (transitionOcurrences[nodeId].find(firstTransitionOnBranch) == transitionOcurrences[nodeId].end()){
          transitionOcurrences[nodeId][firstTransitionOnBranch] = 0;

        }
      }
      transitionOcurrences[nodeId][firstTransitionOnBranch] += 1;
      for (size_t i = 0; i < states.size()-1; i++){
        std::pair<size_t, size_t> transition(states[i], states[i+1]);
        if (transitionOcurrences[nodeId].find(transition) == transitionOcurrences[nodeId].end()){
          transitionOcurrences[nodeId][transition] = 0;
        }

        transitionOcurrences[nodeId][transition] += 1;
     }
    }
  }
  if (!(tree_->isLeaf(tree_->getNode(nodeId)))){
    size_t sonsInitState;
    if (accountedBranch){
      if (states.size() > 0){
        sonsInitState = states[states.size()-1];

      }else{
        sonsInitState = initialState;
      }

    }else{
      sonsInitState = ancetralStates_[nodeId][mappingIndex];
    }

    auto sons = tree_->getSons(nodeId);
    for (size_t n = 0; n < sons.size(); n++){
      getNumOfOcuurencesForEachTransitionPerMappingRecursively(sons[n], sonsInitState, mappingIndex, transitionOcurrences);
    }
  }
}
/******************************************************************************/
std::map<pair<size_t, size_t>, double> StochasticMapping::sumTotalOccurences(std::map<uint, std::map<pair<size_t, size_t>, double>>* transitionOcurrencesPerNodePtr){
  std::map<pair<size_t, size_t>, double> totalExpectedTransitions;
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
  //auto nbState = likelihood_->getStateMap().getNumberOfModelStates();
  auto itNode = transitionOcurrencesPerNodePtr->begin();
  while (itNode != transitionOcurrencesPerNodePtr->end()){
    auto &transitions = (*transitionOcurrencesPerNodePtr)[itNode->first];
    auto itTransitions = transitions.begin();
    while (itTransitions != transitions.end()){
      auto expectedNum = (*transitionOcurrencesPerNodePtr)[itNode->first][itTransitions->first];
      if (totalExpectedTransitions.find(itTransitions->first) == totalExpectedTransitions.end()){
        totalExpectedTransitions[itTransitions->first] = 0;
      }
      totalExpectedTransitions[itTransitions->first] += expectedNum;
      itTransitions ++;
    }
    itNode ++;
  }

  // for (size_t i = 0; i < nbState; i++){
  //   for (size_t j = 0; j < nbState; j++){
  //     if (i == j){
  //       continue;
  //     }
  //     pair<size_t, size_t> pairStates(i,j);
  //     totalExpectedTransitions[pairStates] = 0;
  //     for (size_t n = 0; n < nodeIds.size(); n++){
  //       if ((*transitionOcurrencesPerNodePtr)[nodeIds[n]][pairStates] > 0){
  //         auto expectedNum = (*transitionOcurrencesPerNodePtr)[nodeIds[n]][pairStates];
  //         totalExpectedTransitions[pairStates] += expectedNum;
          
  //       }

  //     }
  //   } 
  // }
  return totalExpectedTransitions;
}

/******************************************************************************/
std::map<pair<size_t, size_t>, double> StochasticMapping::getTotalNumOfOcuurencesForEachTransition(std::map<uint, std::map<pair<size_t, size_t>, double>>* transitionOcurrencesPerNodePtr){
  std::map<uint, std::map<pair<size_t, size_t>, double>> transitionsOccurencesPerNode;
  if (!transitionOcurrencesPerNodePtr){
    transitionsOccurencesPerNode = getNumOfOcuurencesForEachTransitionPerNode();
    transitionOcurrencesPerNodePtr = &transitionsOccurencesPerNode;
  }
  auto totalExpectedTransitions = sumTotalOccurences(transitionOcurrencesPerNodePtr);
  return totalExpectedTransitions;


}
/******************************************************************************/
std::map<pair<size_t, size_t>, double> StochasticMapping::getExpectedNumOfOcuurencesForEachTransition(std::map<uint, std::map<pair<size_t, size_t>, double>>* transitionOcurrencesPerNodePtr){
  std::map<uint, std::map<pair<size_t, size_t>, double>> transitionsOccurencesPerNode;
  if (!transitionOcurrencesPerNodePtr){
    transitionsOccurencesPerNode = getExpectedNumOfOcuurencesForEachTransitionPerNode();
    transitionOcurrencesPerNodePtr = &transitionsOccurencesPerNode;
  }
  auto totalExpectedTransitions = sumTotalOccurences(transitionOcurrencesPerNodePtr);
  return totalExpectedTransitions;
}
/*******************************************************************************/
// for each mapping get the number of occurences of a transition
std::map<size_t, std::map<std::pair<size_t, size_t>, double>> StochasticMapping::getNumOfOccurencesForEachTransitionForEachMapping(){
  std::map<size_t, std::map<std::pair<size_t, size_t>, double>> numOfOccurences;
  //auto nbState = likelihood_->getStateMap().getNumberOfModelStates();
  for (size_t i = 0; i < numOfMappings_; i++){
    std::map<uint, std::map<pair<size_t, size_t>, double>> transitionOcurrences;
    getNumOfOcuurencesForEachTransitionPerMapping(i, transitionOcurrences);
    std::map<pair<size_t, size_t>, double> totalNumOfTransitions = sumTotalOccurences(&transitionOcurrences);
    auto it = totalNumOfTransitions.begin();
    while(it != totalNumOfTransitions.end()){
      numOfOccurences[i][it->first] = totalNumOfTransitions[it->first];
      it ++;
    }
  }
  return numOfOccurences;
}
/*******************************************************************************/
std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> StochasticMapping::getNumOfOccurrencesFromRootToNode(std::map<uint, std::map<size_t, bool>> &presentMapping){
  std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> occurrencesFromRootToNode;
  uint rootIndex = tree_->getRootIndex();
  auto nodeIdsSons = tree_->getSons(rootIndex);
  for (size_t i = 0; i < numOfMappings_; i++){
    std::map<uint, std::map<pair<size_t, size_t>, double>> occurrencesPerMapping;
    getNumOfOcuurencesForEachTransitionPerMapping(i, occurrencesPerMapping);
    ////////////////////// DEBUG //////////////////////////////////////////
    
    // std::cout << "Mapping #" << i << std::endl;
    // auto it = occurrencesPerMapping.begin();
    // while (it != occurrencesPerMapping.end()){
    //   auto &nodeOccurrences = occurrencesPerMapping[it->first];
    //   string nodeName;
    //   if (tree_->isLeaf(it->first)){
    //     nodeName = (tree_->getNode(it->first))->getName();
    //   }else{
    //     nodeName = "N" + std::to_string(it->first);

    //   }
    //   std::cout << "\tNode: " << nodeName << std::endl;
    //   auto itTrans = nodeOccurrences.begin();
    //   while (itTrans != nodeOccurrences.end()){
    //     std::cout << "\t" << (itTrans->first).first << " -> " << (itTrans->first).second << ": " << nodeOccurrences[itTrans->first] << std::endl;
    //     itTrans ++;
    //   }

    //   it ++;
    // }

    //////////////////////////////////////////////////////////////////////
    presentMapping[rootIndex][i] = true;
    std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> occurrencesFromRootPerMapping;
    for (size_t n = 0; n < nodeIdsSons.size(); n++){
      updateFromRootToNodeRecursively(presentMapping, occurrencesPerMapping, i, nodeIdsSons[n], occurrencesFromRootPerMapping);
      auto itNodes = occurrencesFromRootPerMapping.begin();
      while (itNodes != occurrencesFromRootPerMapping.end()){
        //if (tree_->isLeaf(itNodes->first)){
        occurrencesFromRootToNode[itNodes->first][i] = occurrencesFromRootPerMapping[itNodes->first][i];
        //}
        itNodes ++;

      }

    }   
  }
  return occurrencesFromRootToNode;
}
/*******************************************************************************/
void StochasticMapping::updateFromRootToNodeRecursively(std::map<uint, std::map<size_t, bool>> &presentMapping, std::map<uint, std::map<pair<size_t, size_t>, double>> &occurrencesPerMapping, size_t mappingIndex, uint nodeId, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &occurrencesFromRootToNode){
  auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
  uint father = tree_->getNodeIndex(fatherNode);
  bool transitionsFromFather = false;
  if (!(presentMapping[father][mappingIndex])){
    presentMapping[nodeId][mappingIndex] = false;
  
  }else{
    if (notRepresentedNodes_.find(nodeId) != notRepresentedNodes_.end()){
      if (std::find(notRepresentedNodes_[nodeId].begin(), notRepresentedNodes_[nodeId].end(), mappingIndex) != notRepresentedNodes_[nodeId].end()){
        presentMapping[nodeId][mappingIndex] = false;   
      }else{
        presentMapping[nodeId][mappingIndex] = true;
      }
    }else{
      presentMapping[nodeId][mappingIndex] = true;
    }
  }
  if (presentMapping[nodeId][mappingIndex]){
    if (father != tree_->getRootIndex()){
      if (occurrencesFromRootToNode.find(father) != occurrencesFromRootToNode.end()){
        if (occurrencesFromRootToNode[father].find(mappingIndex) != occurrencesFromRootToNode[father].end()){
          auto fatherFromRoot = occurrencesFromRootToNode[father][mappingIndex];
          auto itTransitions = fatherFromRoot.begin();
          while (itTransitions != fatherFromRoot.end()){
            occurrencesFromRootToNode[nodeId][mappingIndex][itTransitions->first] = fatherFromRoot[itTransitions->first];
            transitionsFromFather = true;
            itTransitions ++;
          }
        }
      }
    }
    auto &nodeOccurrences = occurrencesPerMapping[nodeId];
    auto it = nodeOccurrences.begin();

    while (it != nodeOccurrences.end()){
      if (!transitionsFromFather){
        occurrencesFromRootToNode[nodeId][mappingIndex][it->first] = nodeOccurrences[it->first];

      }else{
        if (occurrencesFromRootToNode[nodeId][mappingIndex].find(it->first) != occurrencesFromRootToNode[nodeId][mappingIndex].end()){
          occurrencesFromRootToNode[nodeId][mappingIndex][it->first] += nodeOccurrences[it->first];
        }else{
          occurrencesFromRootToNode[nodeId][mappingIndex][it->first] = nodeOccurrences[it->first];
        }
      }
      it ++;
    }
    
  }
  
  if (!(tree_->isLeaf(nodeId))){
    auto sons =  tree_->getSons(nodeId);
    for (size_t n = 0; n < sons.size(); n++){
      updateFromRootToNodeRecursively(presentMapping, occurrencesPerMapping, mappingIndex, sons[n], occurrencesFromRootToNode);
    }
  }
}


/*******************************************************************************/
std::map<uint, std::map<pair<size_t, size_t>, double>> StochasticMapping::getNumOfOcuurencesForEachTransitionPerNode(){
  std::map<uint, std::map<pair<size_t, size_t>, double>> transitionOcurrences;
  //initMapOfNumOfOccurences(transitionOcurrences);
  for (size_t i = 0; i < numOfMappings_; i++){
    getNumOfOcuurencesForEachTransitionPerMapping(i, transitionOcurrences);
  }
  return transitionOcurrences;
}
/******************************************************************************/
void StochasticMapping::initMapOfNumOfOccurences(std::map<uint, std::map<pair<size_t, size_t>, double>> &transitionOcurrences){
  auto nodeIds = tree_->getNodeIndexes(tree_->getAllNodes());
  auto nbState = likelihood_->getStateMap().getNumberOfModelStates();
  for (size_t n = 0; n < nodeIds.size(); n++){
    if (nodeIds[n] == tree_->getRootIndex()){
      continue;
    }
    for (size_t i = 0; i < nbState; i++){
      for (size_t j = 0; j < nbState; j++){
        if (i == j){
          continue;
        }
        pair<size_t, size_t> pairStates(i,j);
        transitionOcurrences[nodeIds[n]][pairStates] = 0;
      } 
    }   
  }
}
/******************************************************************************/
std::map<uint, std::map<pair<size_t, size_t>, double>> StochasticMapping::getExpectedNumOfOcuurencesForEachTransitionPerNode(){
  // intializing the map
  std::map<uint, std::map<pair<size_t, size_t>, double>> transitionOcurrences = getNumOfOcuurencesForEachTransitionPerNode();
  auto itNode = transitionOcurrences.begin();
  while(itNode != transitionOcurrences.end()){
    auto nodeTransitions = transitionOcurrences[itNode->first];
    auto itTransition = nodeTransitions.begin();
    while(itTransition != nodeTransitions.end()){
      if (notRepresentedNodes_.find(itNode->first) != notRepresentedNodes_.end()){
        transitionOcurrences[itNode->first][itTransition->first] /= (double)(numOfMappings_ - notRepresentedNodes_[itNode->first].size());
      }else{
        transitionOcurrences[itNode->first][itTransition->first] /= (double)numOfMappings_;

      }
      
      itTransition ++;
    }    
    itNode ++;
  }
  return transitionOcurrences;
}
/******************************************************************************/
VVdouble StochasticMapping::getExpectedRateOfTransitionGivenState(Vdouble &dwellingTimesPerState, std::map<std::pair<size_t, size_t>, double> &numOfOccurencesPerTransition){
  auto nbState = likelihood_->getStateMap().getNumberOfModelStates();
  // initializing
  VVdouble expectedRate;
  expectedRate.resize(nbState);
  for (size_t i = 0; i < nbState; i++){
    expectedRate[i].resize(nbState);
    std::fill(expectedRate[i].begin(), expectedRate[i].end(), 0);
  }
  // filling with actual expected rates
  auto itTransitions = numOfOccurencesPerTransition.begin();
  while(itTransitions != numOfOccurencesPerTransition.end()){
    auto beginState = (itTransitions->first).first;
    auto endState = (itTransitions->first).second;
    expectedRate[beginState][endState] += (numOfOccurencesPerTransition[itTransitions->first]/dwellingTimesPerState[beginState]);
    itTransitions ++;
  }
  return expectedRate;

}
/******************************************************************************/
VVdouble StochasticMapping::getDwellingTimeOfStatePerEachMapping(){
  VVdouble dwellingTimesPerMapping;
  auto nbState = likelihood_->getStateMap().getNumberOfModelStates();
  dwellingTimesPerMapping.resize(numOfMappings_);
  for (size_t i = 0; i < numOfMappings_; i++){
    dwellingTimesPerMapping[i].resize(nbState);
    std::fill(dwellingTimesPerMapping[i].begin(), dwellingTimesPerMapping[i].end(), 0);
    getDewellingTimesUnderEachStatePerMapping(&dwellingTimesPerMapping[i], i);
  }
  return dwellingTimesPerMapping;
}

/******************************************************************************/
Vdouble StochasticMapping::getDwellingTimesUnderEachState(bool expectedDuration){
  Vdouble expectedDwellingTimes;
  auto nbState = likelihood_->getStateMap().getNumberOfModelStates();
  expectedDwellingTimes.resize(nbState);
  // set all values to zero initially
  std::fill(expectedDwellingTimes.begin(), expectedDwellingTimes.end(), 0);
  for (size_t i = 0; i < numOfMappings_; i++){
    getDewellingTimesUnderEachStatePerMapping(&expectedDwellingTimes, i);

  }

  if (expectedDuration){
    for (size_t i = 0; i < nbState; i++){
      expectedDwellingTimes[i] /= (double)numOfMappings_;
    }

  }
  return expectedDwellingTimes;
}
/******************************************************************************/
void StochasticMapping::printUnrepresentedLeavesWithCorrespondingMappings(ofstream &stream){
  stream << "# Unrepresennted nodes:" << std::endl;
  auto it = notRepresentedNodes_.begin();
  while (it != notRepresentedNodes_.end()){
    if (tree_->isLeaf(it->first)){
      stream << "\t" << (tree_->getNode(it->first))->getName() << " node id: " << it->first << std::endl;

    }else{
      stream << "\tN" << it->first << " node id: " << it->first << std::endl;
    }
    
    auto &mappingIndices = notRepresentedNodes_[it->first];
    for (size_t i = 0; i < mappingIndices.size(); i++){
      if (i == mappingIndices.size()-1){
        stream << mappingIndices[i] << std::endl;
      }else{
        stream << mappingIndices[i] << ", ";
      }
    }
    stream << "Unrepresented in " << (double)mappingIndices.size()/(double)numOfMappings_ << std::endl;
    stream << "****" << std::endl;
    it ++;
  }
  

}
/******************************************************************************/
bool StochasticMapping::tryToReplaceMapping(double branchLength, uint nodeId, size_t mappingIndex, vector<MutationPath> &mappings, size_t maxNumOfIterations){
  auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
  uint father = tree_->getNodeIndex(fatherNode);
  size_t fatherState = ancetralStates_[father][mappingIndex];
  size_t sonState = ancetralStates_[nodeId][mappingIndex];
  bool success = sampleEvolutionaryPathForBranch(sonState, fatherState, father, nodeId, branchLength, mappingIndex, mappings, maxNumOfIterations, true);
  return success;
}
/******************************************************************************/
std::vector<std::shared_ptr<PhyloTree>> StochasticMapping::createMappingHistoryTrees() const{
  std::vector<std::shared_ptr<PhyloTree>> trees;
  for (size_t i = 0; i < numOfMappings_; i ++){
    auto tree = createMappingHistoryTree(i);
    trees.push_back(tree);
  }
  return trees;
  
}

/******************************************************************************/
std::shared_ptr<PhyloTree> StochasticMapping::createMappingHistoryTree(size_t mappingIndex) const{
  Newick writer;
  Newick reader;
  std::string tree_str = writer.writeTreeToParenthesis(*tree_);
  std::shared_ptr<PhyloTree> tree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree(tree_str));
  //std::shared_ptr<PhyloTree> tree = std::shared_ptr<PhyloTree>(tree_->clone());
  uint rootId = tree_->getRootIndex();
  size_t initialState = ancetralStates_.at(rootId)[mappingIndex];
  assignTransitionOnHistoryTreeRec(rootId, mappingIndex, tree);
  // chnage name here for the root
  auto rootNode = tree->getNode(tree->getRootIndex());
  rootNode->setName("N"+std::to_string(tree->getRootIndex())+"-"+ std::to_string(initialState));
  return tree;

}
/******************************************************************************/
void StochasticMapping::assignTransitionOnHistoryTreeRec(uint nodeId, size_t mappingIndex, std::shared_ptr<PhyloTree> tree) const{

  // get the state of the initial node, change its name so it will include the state
  // go over the transitions, and create the relevant internal nodes with names that
  // include the states
  auto node = tree->getNode(nodeId);
  if (tree->getRootIndex() != nodeId){
    size_t initialState = ancetralStates_.at(nodeId)[mappingIndex];
    if (tree->isLeaf(node)){
       node->setName(node->getName()+"-"+ std::to_string(initialState));
       return;
    }else{
       node->setName("N-"+ std::to_string(initialState));
    }
  }

  auto sonsIds = tree->getSons(nodeId);
  for (size_t i = 0; i < sonsIds.size(); i++){
    auto son = tree->getNode(sonsIds[i]);
    auto edge_to_fragment = tree->getEdgeToFather(son);
    auto mutationPath = mappings_.at(sonsIds[i])[mappingIndex];
    std::vector<size_t> states = mutationPath.getStates();
    std::vector<double> times = mutationPath.getTimes();
    for (size_t j = 0; j < times.size(); j++){
      uint newNodeId = tree->createNodeOnEdge(tree->getEdgeIndex(edge_to_fragment), times[j]);
      (tree->getNode(newNodeId))->setName("N_dummy_"+ std::to_string(newNodeId)+"-"+ std::to_string(states[j]));

    }
    assignTransitionOnHistoryTreeRec(sonsIds[i], mappingIndex, tree);

  }

}

/******************************************************************************/
void StochasticMapping::updateBranchByDwellingTimes(PhyloNode* node, VDouble& dwellingTimes, VVDouble& ancestralStatesFrequencies, size_t divMethod)
{
  // /* first, convert the dwelling times vector to a mutation path of the branch */
  // size_t statesNum = tl_->getNumberOfStates();
  // size_t sonState = getNodeState(node);
  // size_t fatherState = getNodeState(node->getFather());
  // double branchLength = node->getDistanceToFather();
  // double Pf = 1;
  // double Ps = 1;
  // double shareOfFather = 0;
  // double shareOfSon = 0;
  // MutationPath branchMapping(mappingParameters_->getSubstitutionModel()->getAlphabet(), fatherState, branchLength);
  // // set the first event with the dwelling time that matches the state of the father
  // if (fatherState == sonState)
  // {
  //   if (divMethod == 0)
  //   {
  //     if (node->hasFather())
  //     {
  //       Pf = ancestralStatesFrequencies[node->getFather()->getId()][fatherState];
  //     }
  //     Ps = 1;
  //     if (!node->isLeaf())
  //     {
  //       Ps = ancestralStatesFrequencies[node->getId()][sonState];
  //     }
  //     shareOfFather = Pf / (Pf + Ps);
  //     branchMapping.addEvent(fatherState, dwellingTimes[fatherState] * shareOfFather);
  //   }
  //   else
  //   {
  //     branchMapping.addEvent(fatherState, 0);
  //   }
  // }
  // else
  // {
  //   branchMapping.addEvent(fatherState, dwellingTimes[fatherState]);
  // }
  // // set all events except for the one entering the son
  // for (size_t state = 0; state < statesNum; ++state)
  // {
  //   if (state != fatherState && state != sonState && dwellingTimes[state] > 0) // if the state matches an event which is not the first or the last -> add it
  //   {
  //     branchMapping.addEvent(state, dwellingTimes[state]);
  //   }
  // }
  // // change the length of the branch whose bottom node is the son according to the dwelling time of the relevant state
  // if (fatherState == sonState)
  // {
  //   if (divMethod == 0)
  //   {
  //     shareOfSon = 1 - shareOfFather;
  //     node->setDistanceToFather(dwellingTimes[sonState] * shareOfSon);
  //   }
  //   else
  //   {
  //     node->setDistanceToFather(dwellingTimes[sonState]);
  //   }
  // }
  // else
  // {
  //   node->setDistanceToFather(dwellingTimes[sonState]);
  // }

  // /* secondly, update the expected history with the dwelling times-based mutation path */
  // updateBranchMapping(node, branchMapping);
}
