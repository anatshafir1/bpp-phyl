#include "StochasticMapping.h"
#include "../Simulation/MutationProcess.h"
#include "RewardMappingTools.h"
#include "Reward.h"
#include "DecompositionReward.h"
#include "ProbabilisticRewardMapping.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>

#include <iostream>
#include <fstream>
#include <algorithm>
#include <numeric> // to sum over items in a vector

using namespace bpp;
using namespace std;

#define STATE "state"

/******************************************************************************/

StochasticMapping::StochasticMapping(std::shared_ptr<LikelihoodCalculationSingleProcess> drl, size_t numOfMappings, size_t numOfMappingTrials) :
  likelihood_(drl),
  tree_ (make_shared<PhyloTree>(drl->getSubstitutionProcess().getParametrizablePhyloTree())),
//  mappingParameters_(drl->getSubstitutionProcess()),
  ConditionalProbabilities_(),
  nodesCounter_(0),
  numOfMappings_(numOfMappings),
  ancetralStates_(),
  mappings_(),
  jumpsProbs_(),
  notRepresentedNodes_(),
  numOfMappingTrials_(numOfMappingTrials)// ,
  // nodeIdToIndex_()
{
  //giveNamesToInternalNodes(*tree_);                     // set names for the internal nodes of the tree, in case of absence
  ComputeConditionals();
  initJumpProbs();
  //initMappings();
}

/******************************************************************************/

StochasticMapping::~StochasticMapping()
{
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
          auto model = dynamic_cast<const SubstitutionModel*>(likelihood_->getSubstitutionProcess().getModel(modelId));
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
  auto model = dynamic_cast<const SubstitutionModel*>(likelihood_->getSubstitutionProcess().getModel(father, 0)); // father or son??? Should be a father, because the models start at particular nodes, and I should get the model of the preceeding branch
  size_t fatherState = ancetralStates_[father][mapping];
  auto rateToLeave = -1* model->Qij(fatherState, fatherState);
  return rateToLeave;

}
/******************************************************************************/

void StochasticMapping::setExpectedAncestrals(shared_ptr<PhyloTree> expectedMapping, VVDouble& ancestralStatesFrequencies)
{
  // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(expectedMapping);
  // vector<Node*> nodes = ttree->getNodes();
  // for (size_t i = 0; i < nodes.size(); ++i)
  // {
  //   Node* node = nodes[i];
  //   int nodeId = node->getId();
  //   size_t j = static_cast<size_t>(nodeId); //Note@Laurent (Julien 17/06/20): is this really intended, as nodeIds can be discontinuous? Should there be some can of index instead? 
  //   auto d = distance(ancestralStatesFrequencies[j].begin(), max_element(ancestralStatesFrequencies[j].begin(), ancestralStatesFrequencies[j].end()));
  //   size_t state = static_cast<size_t>(d); //Note@Laurent (Julien 17/06/20): assimuming this is always positive, is that so? 
  //   setNodeState(node, state); // in the case of a leaf, the assigned state must be sampled
  // }
}

/******************************************************************************/

shared_ptr<PhyloTree> StochasticMapping::generateExpectedMapping(vector<shared_ptr<PhyloTree>>& mappings, size_t divMethod)
{
  // // initialize the expected history
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  shared_ptr<PhyloTree> expectedMapping(make_shared<PhyloTree>(*tree_));
  // setLeafsStates(expectedMapping);

  // // compute a vector of the posterior asssignment probabilities for each inner node
  // VVDouble ancestralStatesFrequencies;
  // ancestralStatesFrequencies.clear();
  // vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(expectedMapping)->getNodes();
  // size_t statesNum = tl_->getNumberOfStates();
  // ancestralStatesFrequencies.resize(nodes.size(), VDouble(statesNum));
  // computeStatesFrequencies(ancestralStatesFrequencies, mappings);

  // // set the ancestral states accrdonig to the maximal posterior (i.e, conditional) probability
  // setExpectedAncestrals(expectedMapping, ancestralStatesFrequencies);

  // // update the expected history with the dwelling times
  // for (size_t n = 0; n < nodes.size(); ++n)
  // {
  //   Node* node = nodes[n];
  //   if (node->hasFather()) // for any node except to the root
  //   {
  //     // initialize vector of average dwelling times for the branch stemming from node
  //     VDouble AverageDwellingTimes;
  //     AverageDwellingTimes.clear();
  //     AverageDwellingTimes.resize(statesNum, 0);
  //     // compute the average dwelling times of all the states
  //     for (size_t i = 0; i < mappings.size(); ++i)
  //     {
  //       TreeTemplate<Node>* mapping =  dynamic_cast<TreeTemplate<Node>*>(mappings[i]);
  //       // get the pointers to the node and its father in the i'th mapping
  //       Node* curNode = mapping->getNode(node->getName());
  //       Node* father = mapping->getNode(node->getFather()->getName()); // the original father of the node (according to the base tree) in the mapping
  //       while (curNode != father)
  //       {
  //         AverageDwellingTimes[static_cast<size_t>(getNodeState(curNode))] += curNode->getDistanceToFather(); //Note@Laurent (Julien 17/06/20): assuming state is positive, is that so? 
  //         curNode = curNode->getFather();
  //       }
  //     }
  //     double branchLength = node->getDistanceToFather();   // this is the length of the original branch in the base tree
  //     bool updateBranch = true;
  //     for (size_t state = 0; state < statesNum; ++state)
  //     {
  //       AverageDwellingTimes[state] /= static_cast<double>(mappings.size());
  //       if (AverageDwellingTimes[state] == branchLength) // if one of the dwelling times equals the branch length, then there is only one state along te branch and there is no need to edit it
  //       {
  //         updateBranch = false;
  //       }
  //     }
  //     // break the branch according to average dwelling times
  //     if (updateBranch)
  //     {
  //       updateBranchByDwellingTimes(node, AverageDwellingTimes, ancestralStatesFrequencies, divMethod);
  //     }
  //   }
  // }
  // nodesCounter_ = dynamic_cast<TreeTemplate<Node>*>(baseTree_)->getNodes().size() - 1;
  return expectedMapping;
}

/******************************************************************************/

shared_ptr<PhyloTree> StochasticMapping::generateAnalyticExpectedMapping(size_t divMethod)
{
  // /* Compute the posterior assignment probabilities to internal nodes, based on the fractional probablities computed earlier */
  // const vector<int> states =  tl_->getAlphabetStates();
  // vector<int> nodeIds = baseTree_->getNodesId();
  // size_t nodeId;
  // VVDouble posteriorProbabilities;
  // posteriorProbabilities.clear();
  // posteriorProbabilities.resize(baseTree_->getNumberOfNodes(), VDouble(states.size()));
  // double nodeDataProb;
  // // because the sum of partial likelihoods (i.e, the fractional probabilities) is in fact the probablity of the data, it is sufficient to standardize the vector of fractional probabilires for each node to obtain the posterior probabilities
  // for (size_t n = 0; n < baseTree_->getNumberOfNodes(); ++n)
  // {
  //   nodeId = static_cast<size_t>(nodeIds[n]); //Note@Laurent (Julien 17/06/20): what is nodeId is negative? 
  //   nodeDataProb = 0;
  //   for (size_t s = 0; s < states.size(); ++s)
  //   {
  //     nodeDataProb = nodeDataProb + fractionalProbabilities_[nodeId][s];
  //   }
  //   for (size_t nodeState = 0; nodeState < states.size(); ++nodeState)
  //   {
  //     posteriorProbabilities[nodeId][nodeState] = fractionalProbabilities_[nodeId][nodeState] / nodeDataProb;
  //   }
  // }

  // /* Assign states to internal nodes based on the majority rule over the posterior probabilities */
  shared_ptr<PhyloTree> expectedMapping(make_shared<PhyloTree>(*tree_));

  // setLeafsStates(expectedMapping);
  // setExpectedAncestrals(expectedMapping, posteriorProbabilities);

  // /* Compute the reward per state per site - expect two entries per site (that is, two entries in total).
  //    Let r0 be the reward of state 0 nd r1 the reward of state 1. */
  // UserAlphabetIndex1* alpha = new UserAlphabetIndex1(tl_->getAlphabet());
  // DiscreteDistribution* rDist = new ConstantRateDistribution();
  // TransitionModel* tlModel = tl_->getModelForSite(0, 0)->clone();
  // DRTreeLikelihood* drtl = new DRHomogeneousTreeLikelihood(*baseTree_, *(tl_->getData()), tlModel, rDist, false);
  // drtl->initialize();
  // vector<int> ids = baseTree_->getNodesId();
  // const SubstitutionModel* model = dynamic_cast<const SubstitutionModel*>(tl_->getModelForSite(0, 0));

  // /* Compute the expected dwelling times per branch and state as follows:
  //    For branch b of length t, the average welling time in state 0 is r0*t (based on Minin and Suchard paper).
  //    The average dwelling time in state 1 should complement to t (make sure of it!) */
  // vector<Node*> nodes = dynamic_cast<TreeTemplate<Node>*>(expectedMapping)->getNodes();
  // double branchLength;
  // Node* node;
  // VVDouble expectedDwellingTimes;
  // expectedDwellingTimes.clear();
  // expectedDwellingTimes.resize(nodes.size(), VDouble(states.size()));
  // for (size_t s = 0; s < states.size(); ++s)
  // {
  //   alpha->setIndex(states[s], 1); // set the reward of the state as 1 and the reward for the rest of the states as 0
  //   for (size_t m = 0; m < states.size(); ++m)
  //   {
  //     if (m != s)
  //     {
  //       alpha->setIndex(states[m], 0); //Note@Laurent (Julien 17/06/20): can you chack my correction there and above? I changed s/m to states[s] and states[m], is that correct?
  //     }
  //   }
  //   unique_ptr<Reward> reward(new DecompositionReward(model, alpha));
  //   unique_ptr<ProbabilisticRewardMapping> mapping(RewardMappingTools::computeRewardVectors(*drtl, ids, *reward, false));
  //   for (size_t n = 0; n < nodes.size(); ++n)
  //   {
  //     node = nodes[n];
  //     if (node->hasFather()) // for any node except to the root
  //     {
  //       expectedDwellingTimes[static_cast<size_t>(node->getId())][s] = mapping->getReward(node->getId(), 0); //Note@Laurent (Julien 17/06/20): what is nodeId is negative? 
  //     }
  //   }
  // }

  // // standardize expected dwelling itmes, if needed, and update the mapping accorgingly
  // double sumOfDwellingTimes;
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
  //         updateBranch = false;
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

int StochasticMapping::getNodeState(const PhyloNode* node) const
{
  return (dynamic_cast<const BppInteger*>(node->getProperty(STATE)))->getValue();
}

/******************************************************************************/

void StochasticMapping::setNodeState(PhyloNode* node, size_t state)
{
  BppInteger stateProperty(static_cast<int>(state));
  node->setProperty(STATE, stateProperty);
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
  //auto leafsStates = likelihood_->getData();
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

void StochasticMapping::computeStatesFrequencies(VVDouble& ancestralStatesFrequencies, vector<shared_ptr<PhyloTree>>& mappings)
{
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
       
      bool success = sampleMutationsGivenAncestralsPerBranch(father, sons[j], mappingIndex, numOfMappingTrials_);
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

bool StochasticMapping::sampleMutationsGivenAncestralsPerBranch(uint father, uint son, size_t mappingIndex, size_t maxIterNum)
{
  
  size_t fatherState = ancetralStates_[father][mappingIndex];
  size_t sonState = ancetralStates_[son][mappingIndex];

  auto branchPtr = tree_->getIncomingEdges(tree_->getNode(son))[0];
  auto branchLength = branchPtr->getLength();

  /* simulate mapping on a branch until you manage to finish at the son's state */
  bool success = sampleEvolutionaryPathForBranch(sonState, fatherState, father, son, branchLength, mappingIndex, maxIterNum); //TODO put the following lines (inside the for loop) into the new function
  if (!success){
    std::cout << "Mapping failure! " << "Mapping index: " << mappingIndex;
    std::cout << ", nodeId: " << son << ", fatherState: " << fatherState << ", sonState: " << sonState << ", branchLength: " << branchLength;
    std::cout << ", probability of son given father: " << ConditionalProbabilities_[son][fatherState][sonState];
    if (!(father == tree_->getRootIndex())){
      auto grandFather = tree_->getFatherOfNode (tree_->getNode(father));
      uint grandFatherId = tree_->getNodeIndex(grandFather);
      size_t grandFatherState = ancetralStates_[grandFatherId][mappingIndex];
      std::cout << ", father id: " << father << ", grand father id: " << grandFatherId << ", grandFather state: " << grandFatherState;
      std::cout << ", probability of father given grandFather: " << ConditionalProbabilities_[father][grandFatherState][fatherState] << std::endl;

    }

  }
  return success;
}

/******************************************************************************/
void StochasticMapping::getDewellingTimesUnderEachStatePerMapping(vector<double> &dwellingTimes, size_t mappingIndex){
  auto rootId = tree_->getRootIndex();
  auto sons = tree_->getSons(rootId);
  for (size_t i = 0; i < sons.size(); i++){
    getDewellingTimesUnderEachStatePerMappingRecursively(sons[i], ancetralStates_[rootId][mappingIndex], dwellingTimes, mappingIndex);
  }

}
/*****************************************************************************/
bool StochasticMapping::sampleEvolutionaryPathForBranch(size_t sonState, size_t fatherState, uint father, uint son, double branchLength, size_t mappingIndex, size_t maxIterNum, bool replace){
  bool success = true;
  auto alphabet = likelihood_->getData()->getAlphabet();

  auto model = dynamic_cast<const SubstitutionModel*>(likelihood_->getSubstitutionProcess().getModel(father, 0)); // father or son??? Should be a father, because the models start at particular nodes, and I should get the model of the preceeding branch
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
      if (replace){
        mappings_[son][mappingIndex] = tryMapping;

      }else{
        mappings_[son].push_back(tryMapping);
        // *** debug ***//
        if (mappings_[son].size() != mappingIndex+1){
          throw Exception ("StochasticMapping::sampleMutationsGivenAncestralsPerBranch: Something went wrong when filling mappings_ object!");
        }

      }

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
void StochasticMapping::getDewellingTimesUnderEachStatePerMappingRecursively(uint nodeId, size_t initialState, vector<double> &dwellingTimes, size_t mappingIndex){
  auto mutationPath = mappings_[nodeId][mappingIndex];
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
      // auto &failedMappings = notRepresentedNodes_[nodeId];
      // if (std::find(failedMappings.begin(), failedMappings.end(), mappingIndex) != failedMappings.end()){
      //   return;
      // }

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


      getDewellingTimesUnderEachStatePerMappingRecursively(sons[n], initialStateForSon, dwellingTimes, mappingIndex);
      
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
std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> StochasticMapping::getNumOfOccurrencesFromRootToTip(std::map<uint, std::map<size_t, bool>> &presentMapping){
  std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> occurrencesFromRootToLeaf;
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
      updateFromRootToLeafRecursively(presentMapping, occurrencesPerMapping, i, nodeIdsSons[n], occurrencesFromRootPerMapping);
      auto itNodes = occurrencesFromRootPerMapping.begin();
      while (itNodes != occurrencesFromRootPerMapping.end()){
        if (tree_->isLeaf(itNodes->first)){
          occurrencesFromRootToLeaf[itNodes->first][i] = occurrencesFromRootPerMapping[itNodes->first][i];
        }
        itNodes ++;

      }

    }   
  }
  return occurrencesFromRootToLeaf;
}
/*******************************************************************************/
void StochasticMapping::updateFromRootToLeafRecursively(std::map<uint, std::map<size_t, bool>> &presentMapping, std::map<uint, std::map<pair<size_t, size_t>, double>> &occurrencesPerMapping, size_t mappingIndex, uint nodeId, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &occurrencesFromRootToLeaf){
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
      if (occurrencesFromRootToLeaf.find(father) != occurrencesFromRootToLeaf.end()){
        if (occurrencesFromRootToLeaf[father].find(mappingIndex) != occurrencesFromRootToLeaf[father].end()){
          auto fatherFromRoot = occurrencesFromRootToLeaf[father][mappingIndex];
          auto itTransitions = fatherFromRoot.begin();
          while (itTransitions != fatherFromRoot.end()){
            occurrencesFromRootToLeaf[nodeId][mappingIndex][itTransitions->first] = fatherFromRoot[itTransitions->first];
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
        occurrencesFromRootToLeaf[nodeId][mappingIndex][it->first] = nodeOccurrences[it->first];

      }else{
        if (occurrencesFromRootToLeaf[nodeId][mappingIndex].find(it->first) != occurrencesFromRootToLeaf[nodeId][mappingIndex].end()){
          occurrencesFromRootToLeaf[nodeId][mappingIndex][it->first] += nodeOccurrences[it->first];
        }else{
          occurrencesFromRootToLeaf[nodeId][mappingIndex][it->first] = nodeOccurrences[it->first];
        }
      }
      it ++;
    }
    
  }
  
  if (!(tree_->isLeaf(nodeId))){
    auto sons =  tree_->getSons(nodeId);
    for (size_t n = 0; n < sons.size(); n++){
      updateFromRootToLeafRecursively(presentMapping, occurrencesPerMapping, mappingIndex, sons[n], occurrencesFromRootToLeaf);
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
    getDewellingTimesUnderEachStatePerMapping(dwellingTimesPerMapping[i], i);
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
    getDewellingTimesUnderEachStatePerMapping(expectedDwellingTimes, i);

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
bool StochasticMapping::tryToReplaceMapping(double branchLength, uint nodeId, size_t mappingIndex, size_t maxNumOfIterations){
  auto fatherNode = tree_->getFatherOfNode(tree_->getNode(nodeId));
  uint father = tree_->getNodeIndex(fatherNode);
  size_t fatherState = ancetralStates_[father][mappingIndex];
  size_t sonState = ancetralStates_[nodeId][mappingIndex];
  bool success = sampleEvolutionaryPathForBranch(sonState, fatherState, father, nodeId, branchLength, mappingIndex, maxNumOfIterations, true);
  return success;
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
