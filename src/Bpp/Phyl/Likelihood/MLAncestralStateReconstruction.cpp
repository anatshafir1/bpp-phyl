#include "MLAncestralStateReconstruction.h"

using namespace bpp;
using namespace std;


void MLAncestralStateReconstruction::computeJointLikelihood(){
    ancestors_ = new std::map<int, std::map<size_t, std::vector<size_t>>>;
    std::vector<int> nodesIds = tree_.getNodesId();
    for (size_t i = 0; i < nodesIds.size(); i++){
      for (size_t j = 0; j < nbDistinctSites_; j ++){
        std::vector <size_t> ancestorsPerSite;
        ancestorsPerSite.reserve(nbStates_);
        //(*ancestors_)[nodesIds[i]][j].reserve(nbStates_);
        for (size_t s = 0; s < nbStates_; s++){
          ancestorsPerSite.push_back(0);
          //(*ancestors_)[nodesIds[i]][j].push_back(0);
        }
        (*ancestors_)[nodesIds[i]][j] = ancestorsPerSite;
      }     
    }
    likelihoodData_ = new DRASRTreeLikelihoodData(&tree_, nbCategories_);
    likelihoodData_->initLikelihoodsForAncestralReconstruction(*data_, *model_, pxy_, ancestors_);
    computeSubtreeBestLikelihood(tree_.getRootNode());

}
/*********************************************************************************************************/
void MLAncestralStateReconstruction::computeSubtreeBestLikelihood(const Node* node)
{
  if (node->isLeaf()){
    return;
  }

  size_t nbSites  = likelihoodData_->getLikelihoodArray(node->getId()).size();
  size_t nbNodes  = node->getNumberOfSons();

  // Must reset the likelihood array first (i.e. set all of them to 1):
  VVVdouble* _likelihoods_node = &likelihoodData_->getLikelihoodArray(node->getId());
  for (size_t i = 0; i < nbSites; i++)
  {
    //For each site in the sequence,
    VVdouble* _likelihoods_node_i = &(*_likelihoods_node)[i];
    for (size_t c = 0; c < nbClasses_; c++)
    {
      //For each rate classe,
      Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
      for (size_t x = 0; x < nbStates_; x++)
      {
        //For each initial state,
        (*_likelihoods_node_i_c)[x] = 1.;
      }
    }
  }

  for (size_t l = 0; l < nbNodes; l++)
  {
    //For each son node,
    const Node* son = node->getSon(l);
    computeSubtreeBestLikelihood(son); //Recursive method:
  }

  if (node->hasFather()){
    VVVdouble* currNodeProb = &((*pxy_)[node->getId()]);
    fillLikelihoodsArrays(node, _likelihoods_node, currNodeProb, nbSites, nbClasses_, nbStates_, nbNodes);
  }else{
    fillRootLikelihoodsArrays(node, _likelihoods_node, nbSites, nbClasses_, nbStates_, nbNodes);
  }

}
/*****************************************************************************************************************/
//For conviniency I can fill only the array of the best reconstructed state of the root,
// thus in case of the root the LR(i) is the best likelihood overall given the state i of the root
void MLAncestralStateReconstruction::fillRootLikelihoodsArrays(const Node* node, VVVdouble* likNode, size_t nbSites, size_t nbClasses, size_t nbStates, size_t nbNodes){
  for (size_t i = 0; i < nbSites; i++){
  //For each site in the sequence,
    VVdouble* _likelihoods_node_i = &(*likNode)[i];
    for (size_t c = 0; c < nbClasses; c++){
      //For each rate classe,
      Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
      for (size_t x = 0; x < nbStates; x++){
        //For each state of the root
        (*_likelihoods_node_i_c)[x] = rootFrequiencies_[x];
        size_t nbSons = node->getNumberOfSons();
        for (size_t n = 0; n < nbSons; n++){
          const Node* son = node->getSon(n);
          vector<size_t> * _patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
          VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
          VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
          Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
          (*_likelihoods_node_i_c)[x] *= (*_likelihoods_son_i_c)[x];

        }
        
      }

    }

  }
}

/*****************************************************************************************************************/
void MLAncestralStateReconstruction::fillLikelihoodsArrays(const Node* node, VVVdouble* likNode,  VVVdouble* currNodeProb, size_t nbSites, size_t nbClasses, size_t nbStates, size_t nbNodes){
  for (size_t i = 0; i < nbSites; i++){
    //For each site in the sequence,
    VVdouble* _likelihoods_node_i = &(*likNode)[i];
    for (size_t c = 0; c < nbClasses; c++){
      //For each rate classe,
      Vdouble* _likelihoods_node_i_c = &(*_likelihoods_node_i)[c];
      for (size_t x = 0; x < nbStates; x++){
        //For each father state
        double maxLikelihood_y = 0;
        //int maxState = 0;
        for (size_t y = 0; y < nbStates; y++){
          //for each current node state
          double likelihood_y = 1.;
          for (size_t n = 0; n < nbNodes; n++){
            //iterating over the sons           
            const Node* son = node->getSon(n);
            vector<size_t> * _patternLinks_node_son = &likelihoodData_->getArrayPositions(node->getId(), son->getId());
            VVVdouble* _likelihoods_son = &likelihoodData_->getLikelihoodArray(son->getId());
            VVdouble* _likelihoods_son_i = &(*_likelihoods_son)[(*_patternLinks_node_son)[i]];
            Vdouble* _likelihoods_son_i_c = &(*_likelihoods_son_i)[c];
            likelihood_y *= (*_likelihoods_son_i_c)[y];
              
          }
          likelihood_y *= (*currNodeProb)[c][x][y];
          if (likelihood_y > maxLikelihood_y){
            maxLikelihood_y = likelihood_y;
            (*ancestors_)[node->getId()][i][x] = y;
          }
        }
        (*_likelihoods_node_i_c)[x] *= maxLikelihood_y;
        
      }

    }

  }
}
/****************************************************************************************************/
std::map<int, std::vector<size_t> > MLAncestralStateReconstruction::getAllAncestralStates() const{
  std::map<int, std::vector<size_t> > ancestors;
  getAllAncestralStatesRec(tree_.getRootNode(), ancestors);
  return ancestors;

}
/*****************************************************************************************************/
void MLAncestralStateReconstruction::getAllAncestralStatesRec(const Node* node, std::map<int, std::vector<size_t> > &ancestors) const{
  if (node->isLeaf()){
    int leafId = node->getId();
    ancestors[leafId].reserve(nbDistinctSites_);
    for (size_t i = 0; i < nbDistinctSites_; i++){
      //we do not care which state is the father
      ancestors[leafId].push_back((*ancestors_)[leafId][i][0]);   
    }
    return;
  }else{
    int nodeId = node->getId();
    if (tree_.isRoot(nodeId)){
      for (size_t i = 0; i < nbDistinctSites_; i++){
        VVVdouble* likelihoodNode = &likelihoodData_->getLikelihoodArray(node->getId());
        Vdouble likValues = (*likelihoodNode)[i][0];//assuming there is only one class      
        size_t bestState = 0;
        double bestLikValue = 0;
        for (size_t s = 0; s < nbStates_; s++){
          if (likValues[s] > bestLikValue){
            bestLikValue = likValues[s];
            bestState = s;
          }
        }       
        ancestors[nodeId].push_back(bestState);
      }      
    }else{
      int fatherId = node->getFatherId();
      for (size_t i = 0; i < nbDistinctSites_; i++){
        size_t bestFatherState = ancestors[fatherId][i];
        ancestors[nodeId].push_back((*ancestors_)[nodeId][i][bestFatherState]);

      }
    }
    size_t nbSons = node->getNumberOfSons();
    for (size_t n = 0; n < nbSons; n++){
      const Node* son = node->getSon(n);
      getAllAncestralStatesRec(son, ancestors);
    }
  }
}
  
