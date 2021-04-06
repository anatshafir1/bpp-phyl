#include "FwLikMLAncestralReconstruction.h"

using namespace bpp;
using namespace std;


void FwLikMLAncestralReconstruction::initialize(const AlignedValuesContainer& sites){
  nbSites_ = Eigen::Index(sites.getNumberOfSites ()); 
  likelihoodMatrixDim_ = conditionalLikelihoodDimension (nbState_, nbSites_);
  ConditionalLikelihoodForwardRef bidonRoot=ConstantZero<MatrixLik>::create(context_, MatrixDimension(1,1));
  createNode(bidonRoot);
      
  rootAt(bidonRoot); // for construction, temporary top node for new edges
  auto n = makeForwardComputationAtRoot (processTree_->getRoot(), sites);
  auto sonsMulNode = n;
  //std::cerr << "   -> N : " <<  getNodeIndex(sonsMulNode) << "; " << sonsMulNode->getTargetValue() << std::endl;
  rootAt(sonsMulNode);
  deleteNode(bidonRoot);

  //RootConditionalLikelihoodsRef rootNode;
  //auto rootFreqs = CWiseFill<MatrixLik, RowLik>::create(context_, {rFreqs_}, likelihoodMatrixDim_);
  //std::cerr << "   -> Root frequencies : "<< " " << rootFreqs->getTargetValue() << std::endl;
  //rootNode = MatrixMaxProduct<RowLik, MatrixLik, MatrixLik>::create (
                           //context_, {rootFreqs, sonsMulNode}, RowVectorDimension (nbSites_));

  // just for debug:
  //auto vec = rootFreqs->getTargetValue().array() * sonsMulNode->getTargetValue().array();
  //cout << "****" << endl;
  //cout << "The multiplication at the root: " << endl;
  //cout << vec << endl;
  //cout <<"*****" << endl;
  
  //////////////////////////////////////////
  //std::cerr << "   rootNode->likelihood: " << "; " << rootNode->getTargetValue() << std::endl;
  //std::cerr << "   rootNode->likelihood: " << "; " << rootNode->getTargetValue() << std::endl;
  //ValueRef<double> val;
  //val = SumOfLogarithms<RowLik>::create (context_, {rootNode}, RowVectorDimension (Eigen::Index (nbSites_)));
  //std::cerr << "   final loglikelihood: " << "; " << val->getTargetValue() << std::endl;
  //return val;
}


ConditionalLikelihoodForwardRef FwLikMLAncestralReconstruction::makeForwardLikelihoodAtNode (shared_ptr<ProcessNode> processNode, const AlignedValuesContainer & sites, ConditionalLikelihoodForwardRef* forwardEdge)
{
  // if (processTree_->getRootIndex() == processTree_->getNodeIndex(processNode)){
  //   return setRootFrequencies(processNode, sites);

  // }
  const auto childBranches = processTree_->getBranches (processNode);
  auto spIndex=processNode->getSpeciesIndex();
  ConditionalLikelihoodForwardRef forwardNode;
  const auto edge = processTree_->getIncomingEdges(processNode)[0];
  const auto brlen= edge->getBrLen();
  const auto model= edge->getModel();
  const auto nMod = edge->getNMod();
  auto zero=NumericConstant<size_t>::create(context_, size_t(0));

  if (childBranches.empty ())
  {
    ConditionalLikelihoodForwardRef leafNode;
    leafNode = makeInitialConditionalLikelihood (processNode->getName (), sites);

    //const auto edge = processTree_->getEdge(edgeIndex);
    if (dynamic_cast<const TransitionModel*>(model->getTargetValue()))
    {
      auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, zero, nMod}, transitionMatrixDimension (size_t(nbState_)));
      
      edge->setTransitionMatrix(transitionMatrix);
      forwardNode = ForwardTransition::create (
        context_, {transitionMatrix, leafNode}, likelihoodMatrixDim_);
      if (forwardEdge){
        *forwardEdge = MatrixArgMaxProduct<MatrixLik, MatrixLik, MatrixLik>::create (
          context_, {transitionMatrix, leafNode}, likelihoodMatrixDim_);
      }
      
    }
    if (!hasNodeIndex(forwardNode)) 
    {
      setNodeIndices(forwardNode, processNode, spIndex);
    }
  }
  else {
    // depE are edges used to link ForwardLikelihoodTree edges to this
    // node
    std::vector<ConditionalLikelihoodForwardRef> depE(childBranches.size());
    NodeRefVec depsForSonsMul(childBranches.size());
    std::vector<ConditionalLikelihoodForwardRef> sonsEdges(childBranches.size());
    
    for (size_t i = 0; i < childBranches.size (); ++i) {
        ConditionalLikelihoodForwardRef sonEdge;
        depE[i] =  makeForwardLikelihoodAtNode(processTree_->getSon(childBranches[i]), sites, &sonEdge);
        depsForSonsMul[i] = depE[i];
        sonsEdges[i] = sonEdge;
        
    }
    // L_son1(j) * L_son2(j) --> j is the state of the father
    auto sonsMulNode = SpeciationForward::create(context_, std::move(depsForSonsMul),
                                                          likelihoodMatrixDim_);
    auto transitionMatrix = ConfiguredParametrizable::createMatrix<ConfiguredModel, TransitionMatrixFromModel> (context_, {model, brlen, zero, nMod, factorNode_}, transitionMatrixDimension (size_t(nbState_)));
      
    edge->setTransitionMatrix(transitionMatrix);

    // L_node(i) = max_j{P_ij * (L_son1(j) * L_son2(j))}
    // Note: the MaxJointLik operator is very similar to MatrixProduct.
    // The only difference is that instead of summation, MaxJointLik uses the operator max().
    forwardNode = MaxJointLik::create (
        context_, {transitionMatrix, sonsMulNode}, likelihoodMatrixDim_);
    if (forwardEdge){
      *forwardEdge = MatrixArgMaxProduct<MatrixLik, MatrixLik, MatrixLik>::create (
        context_, {transitionMatrix, sonsMulNode}, likelihoodMatrixDim_);
    }


    if (!hasNodeIndex(forwardNode)) 
    {
      setNodeIndices(forwardNode, processNode, spIndex);
    }
    for (size_t n = 0; n < sonsEdges.size(); n++){
      if (!hasEdgeIndex(sonsEdges[n])){
        link (forwardNode, depE[n], sonsEdges[n]);
        setEdgeIndex(sonsEdges[n], processTree_->getEdgeIndex(childBranches[n]));
        auto spIndexEdge = childBranches[n]->getSpeciesIndex();
    
        if (brlen)
        {
          if (mapEdgesIndexes_.find(spIndexEdge)==mapEdgesIndexes_.end())
            mapEdgesIndexes_[spIndexEdge]=DAGindexes();
          mapEdgesIndexes_[spIndexEdge].push_back(getEdgeIndex(sonsEdges[n]));
        }
      }
        

    }
      //linkNodes(forwardNode, childBranches, depE);
      
    
    //std::cerr << "   -> son1 of N "<< getNodeIndex(forwardNode) << " is: " << getNodeIndex(depE[0]) << "; " << depE[0]->getTargetValue() << std::endl; 
    //std::cerr << "   -> son2 of N " << getNodeIndex(forwardNode) << " is: " << getNodeIndex(depE[1]) << "; " << depE[1]->getTargetValue() << std::endl;  
  }
  std::cerr <<  " -> N " << getNodeIndex(forwardNode) << " same as  -> N " << processTree_->getNodeIndex(processNode) << std::endl;

  //std::cerr << "   -> N " << getNodeIndex(forwardNode) << "; " << forwardNode->getTargetValue() << std::endl;
  
  return(forwardNode);
}
/********************************************************************************************************************************************/

ConditionalLikelihoodForwardRef FwLikMLAncestralReconstruction::makeInitialConditionalLikelihood (const string & sequenceName, const AlignedValuesContainer & sites)
{
  size_t nbSites=sites.getNumberOfSites();
  const auto sequenceIndex = sites.getSequencePosition (sequenceName);
  MatrixLik initCondLik (nbState_, nbSites);
  for (size_t site = 0; site < nbSites; ++site) {
    for (auto state = 0; state < nbState_; ++state) {
      initCondLik (Eigen::Index (state), Eigen::Index (site)) =
        sites (site, sequenceIndex, statemap_.getAlphabetStateAsInt(size_t(state)));
    }
  }

  return Sequence_DF::create (context_, std::move(initCondLik), sequenceName);
}

/**************************************************************************************/
ConditionalLikelihoodForwardRef FwLikMLAncestralReconstruction::makeForwardComputationAtRoot(std::shared_ptr<ProcessNode> processNode, const AlignedValuesContainer & sites){
  const auto childBranches = processTree_->getBranches (processNode);
  auto spIndex=processNode->getSpeciesIndex();
  std::vector<ConditionalLikelihoodForwardRef> depE(childBranches.size());
  std::vector<ConditionalLikelihoodForwardRef> sonsEdges(childBranches.size());
  NodeRefVec depsForSonsMul(childBranches.size());
    
  for (size_t i = 0; i < childBranches.size (); ++i) {
    ConditionalLikelihoodForwardRef sonEdge;
    depE[i] =  makeForwardLikelihoodAtNode(processTree_->getSon(childBranches[i]), sites, &sonEdge);
    depsForSonsMul[i] = depE[i];
    sonsEdges[i] = sonEdge;      
  }
  auto sonsMulNode = SpeciationForward::create(context_, std::move(depsForSonsMul),
                                                          likelihoodMatrixDim_);

  setNodeIndices(sonsMulNode, processNode, spIndex);
  for (size_t n = 0; n < sonsEdges.size(); n++){
    if (!hasEdgeIndex(sonsEdges[n])){
      link (sonsMulNode, depE[n], sonsEdges[n]);
      setEdgeIndex(sonsEdges[n], processTree_->getEdgeIndex(childBranches[n]));
      auto spIndexEdge = childBranches[n]->getSpeciesIndex();
    
      //if (brlen)
      //{
      if (mapEdgesIndexes_.find(spIndexEdge)==mapEdgesIndexes_.end())
        mapEdgesIndexes_[spIndexEdge]=DAGindexes();
      mapEdgesIndexes_[spIndexEdge].push_back(getEdgeIndex(sonsEdges[n]));
      //}
    }       

  }
  //linkNodes(sonsMulNode, childBranches, depE);
  return sonsMulNode;
  // RootConditionalLikelihoodsRef rootNode;
  // auto rootFreqs = CWiseFill<MatrixLik, RowLik>::create(context_, {rFreqs_}, likelihoodMatrixDim_);
  // rootNode = LikelihoodRootConditional::create (
  //                          context_, {rootFreqs, sonsMulNode}, RowVectorDimension (nbSites_));
  // ValueRef<double> val;
  // val = SumOfLogarithms<RowLik>::create (context_, {rootNode}, RowVectorDimension (Eigen::Index (nbSites_)));



  //return val;
}
/**************************************************************************************/
void FwLikMLAncestralReconstruction::linkNodes(ConditionalLikelihoodForwardRef forwardNode, std::vector<std::shared_ptr<ProcessEdge>> childBranches, std::vector<ConditionalLikelihoodForwardRef> depE){
  
  auto currNodeIndex = getNodeIndex(forwardNode);
  auto edgesFrom = getOutgoingEdges(currNodeIndex);
  for (size_t i = 0; i < childBranches.size (); ++i) {
  // first we need to check if the link already exists
    auto sonIndex = getNodeIndex(depE[i]);
    auto edgesTo = getIncomingEdges(sonIndex);
    for (size_t j = 0; j < edgesTo.size(); j++){
      if (std::find(edgesFrom.begin(), edgesFrom.end(), edgesTo[j]) == edgesFrom.end()){
        link(forwardNode, depE[i]);
      }
    }
    if (edgesTo.size() > 1){
      throw Exception("FwLikMLAncestralReconstruction: more than one incoming edge!!");
    }
  } 
}
/**************************************************************************************/
void FwLikMLAncestralReconstruction::setNodeIndices(ConditionalLikelihoodForwardRef forwardNode, std::shared_ptr<ProcessNode> processNode, unsigned int spIndex){
  createNode(forwardNode);
  setNodeIndex(forwardNode, processTree_->getNodeIndex(processNode));
  if (mapNodesIndexes_.find(spIndex)==mapNodesIndexes_.end()){
    mapNodesIndexes_[spIndex]=DAGindexes();
  }
  mapNodesIndexes_[spIndex].push_back(getNodeIndex(forwardNode));

}
/**************************************************************************************/
