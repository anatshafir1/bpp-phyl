#include "MarginalNonRevAncestralStateReconstruction.h"
#include <Bpp/Numeric/VectorTools.h>

using namespace bpp;
using namespace std;


//****************************************************************************************************************************************/

void MarginalNonRevAncestralStateReconstruction::computePosteriorProbabilitiesOfNodesForEachStatePerSite(){
  postProbNode_ = new map<int, map<size_t, std::vector<double>>>;
  jointProbabilities_ = new std::map<int, std::map<size_t, VVdouble>>;

  std::vector <int> nodesIds = tree_.getNodesId();
  //initializing the nodePostProb map
  for (size_t n = 0; n < nodesIds.size(); n++){
    for (size_t i = 0; i < nbDistinctSites_; i++){
      (*postProbNode_)[nodesIds[n]][i].reserve(nbStates_);
      (*jointProbabilities_)[nodesIds[n]][i].reserve(nbStates_);
      for (size_t s = 0; s < nbStates_; s++){
        (*postProbNode_)[nodesIds[n]][i].push_back(0);
        Vdouble fatherJointProb;
        fatherJointProb.reserve(nbStates_);
        for (size_t fatherState = 0; fatherState < nbStates_; fatherState++){
          Vdouble fatherStateProb;
          fatherJointProb.push_back(0);

        }
        (*jointProbabilities_)[nodesIds[n]][i].push_back(fatherJointProb);

      }
    }   
  }
  // Filling the posterior probabilities
  for (size_t rootState = 0; rootState < nbStates_; rootState++){
    map<int, map<size_t, std::vector<double>>>* postProbPerRootState = getPosteriorProbabilitiesOfNodesForEachRootStatePerSite(likelihood_, rootState, nodesIds);
    for (size_t n = 0; n < nodesIds.size(); n++){
      for (size_t i = 0; i < nbDistinctSites_; i++){
        for (size_t state = 0; state < nbStates_; state++){
          (*postProbNode_)[nodesIds[n]][i][state] += (*postProbPerRootState)[nodesIds[n]][i][state];
        }
      }
    }
    delete postProbPerRootState;
  }

}

//***********************************************************************************************************************/
double MarginalNonRevAncestralStateReconstruction::getJointLikelihoodFatherNode(
  DRNonHomogeneousTreeLikelihood* NHdrl, 
  int nodeId, 
  size_t nodeState, 
  size_t fatherState,
  size_t rootState,
  size_t site,
  size_t nbClasses)
{

  //std::vector<double> p = NHdrl->getRateDistribution()->getProbabilities();
  //std::vector<double> treeLikelihood = NHdrl->getLikelihoodData()->getRootRateSiteLikelihoodArray();

  double rootFreq = NHdrl->getRootFrequencies(site)[rootState];
  VVVdouble Pijt = NHdrl->getTransitionProbabilitiesPerRateClass(nodeId, site);
  

  int fatherId = tree_.getFatherId(nodeId);
  
  DRASDRTreeLikelihoodData* likelihoodData = NHdrl->getLikelihoodData();
  map<int, VVVdouble>* likelihoods_node = &likelihoodData->getLikelihoodArrays(nodeId);
  map<int, VVVdouble>* likelihoodData_father = &likelihoodData->getLikelihoodArrays(fatherId);
  VVVdouble* likelihoods_node_father = &(*likelihoods_node)[fatherId];
  VVVdouble* likelihoods_father_node = &(*likelihoodData_father)[nodeId];

  VVdouble* likelihoods_node_father_site = &(*likelihoods_node_father)[site];
  VVdouble* likelihoods_father_node_site = &(*likelihoods_father_node)[site];
  double postProb = 0;
  for (size_t c = 0; c < nbClasses; c++)
  {
    Vdouble* likelihoods_node_father_site_c = &(*likelihoods_node_father_site)[c];
    Vdouble* likelihoods_father_node_site_c = &(*likelihoods_father_node_site)[c];

    double FatherTerm_y = (*likelihoods_node_father_site_c)[fatherState];
    double node_lik_x = (*likelihoods_father_node_site_c)[nodeState];
    double Pyx = Pijt[c][fatherState][nodeState];
    postProb += ((FatherTerm_y * Pyx * node_lik_x * r_[c] * rootFreq)/l_[site]);


  }
  return postProb;


}
//***************************************************************************************************************************************/
map<int, map<size_t, std::vector<double>>>* MarginalNonRevAncestralStateReconstruction::getPosteriorProbabilitiesOfNodesForEachRootStatePerSite(
  DRNonHomogeneousTreeLikelihood* NHdrl,
  size_t rootState,
  std::vector<int> nodeIds)
{

  map<int, map<size_t, std::vector<double>>>* postProbs = new map<int, map<size_t, std::vector<double>>>;
  const TreeTemplate<Node> tree = NHdrl->getTree();
  const Node* rootNode = tree.getRootNode();
  NHdrl->computeLikelihoodPrefixConditionalOnRoot(rootNode, rootState);
  for (size_t n = 0; n < nodeIds.size(); n++){
    int nodeId = nodeIds[n];
    for (size_t i = 0; i < nbDistinctSites_; i++){
      (*postProbs)[nodeId][i].reserve(nbStates_);
      for (size_t nodeState = 0; nodeState < nbStates_; nodeState++){
        double nodePostProb = 0;
        //calculate for the root
        if (!tree.hasFather(nodeId)){
          if (nodeState == rootState){
            VVVdouble rootLikelihoods  = NHdrl->getLikelihoodData()->getRootLikelihoodArray();
            double rootFreq = NHdrl->getRootFrequencies(i)[rootState];
            for (size_t c = 0; c < nbClasses_; c++){
              nodePostProb += (rootFreq * rootLikelihoods[i][c][rootState] * r_[c]);
            }
            nodePostProb /= l_[i];
          }
          (*postProbs)[nodeId][i].push_back(nodePostProb);
          continue;
        }
        //if it is not the root
        for (size_t fatherState = 0; fatherState < nbStates_; fatherState ++){
          double jointNodeFatherRootProb = getJointLikelihoodFatherNode(NHdrl, nodeId, nodeState, fatherState, rootState, i, nbClasses_);
          (*jointProbabilities_)[nodeId][i][nodeState][fatherState] += jointNodeFatherRootProb;
          nodePostProb += jointNodeFatherRootProb;
        }
        (*postProbs)[nodeId][i].push_back(nodePostProb);
      }
    }
  }
  return postProbs;
}
//*****************************************************************************************************************************/

const std::map<int, std::vector<size_t> > MarginalNonRevAncestralStateReconstruction::getAllAncestralStates() const{
  std::map<int, std::vector<size_t> > ancestors;
  std::vector<int> nodesIds = tree_.getNodesId();
  for (size_t n = 0; n < nodesIds.size(); n++){
    int nodeId = nodesIds[n];
    ancestors[nodeId].reserve(nbSites_);
    for (size_t i = 0; i <nbSites_; i++){
      size_t bestState = VectorTools::whichMax((*postProbNode_)[nodeId][i]);
      ancestors[nodeId].push_back(bestState);

    }
  }
  return ancestors;
}