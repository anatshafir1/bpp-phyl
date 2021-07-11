#include "JointMLAncestralReconstruction.h"

using namespace bpp;
using namespace std;

void JointMLAncestralReconstruction::init(){
    ancestors_ = new std::map<uint, std::vector<size_t>>;
    fathers_ = new std::map<uint, uint>;
    shared_ptr<AlignedValuesContainer> data(dynamic_cast<AlignedValuesContainer*>(likelihood_->getShrunkData()->clone()));
    recursiveJointMLAncestralStates(tree_->getRoot(), *data);
    
}
/**********************************************************************************/
void JointMLAncestralReconstruction::recursiveJointMLAncestralStates(
                    const std::shared_ptr<PhyloNode> node,
                    AlignedValuesContainer& data) const{

    vector<size_t> ancestorsPerNode;
    ancestorsPerNode.reserve(nbDistinctSites_);
    uint nodeIndex = tree_->getNodeIndex(node);
    if (tree_->isLeaf(node)){ 
        ancestorsPerNode = getAncestralStatesForNode(tree_->getNodeIndex(node));
        (*ancestors_)[nodeIndex] = ancestorsPerNode;
    }else
    {
        ancestorsPerNode = getAncestralStatesForNode(tree_->getNodeIndex(node));
        (*ancestors_)[nodeIndex] = ancestorsPerNode;
        vector<shared_ptr<PhyloNode> > vsons=tree_->getSons(node);
        for (size_t n = 0; n < vsons.size(); n++){
            uint sonId = tree_->getNodeIndex(vsons[n]);
            (*fathers_)[sonId] = nodeIndex;
        }
  
        for (size_t i = 0; i < vsons.size(); i++){
            recursiveJointMLAncestralStates(vsons[i], data);
        }
    }
}
/**********************************************************************************/
std::vector<size_t> JointMLAncestralReconstruction::getAncestralStatesForNode(uint nodeId) const{
    vector <size_t> ancestors;
    if (nodeId == tree_->getNodeIndex(tree_->getRoot())){
        auto rootReconstruction = likelihood_->getJointMLAncestralReconstructionFromRoot(nodeId)->getTargetValue();
        std::cerr <<  " -> Root: N " << nodeId <<": " << rootReconstruction << std::endl;
        for (size_t i = 0; i < nbDistinctSites_; i++){
            size_t rootStateForCurrSite = (size_t)rootReconstruction(i);
            ancestors.push_back(rootStateForCurrSite);
        }

    }else{
        auto nodeReconstruction = likelihood_->getJointMLAncestralReconstructionAtNode(nodeId)->getTargetValue();
        //std::cerr <<  " ->  N " << nodeId <<": " << nodeReconstruction << std::endl;
        uint fatherId = (*fathers_).at(nodeId);
        
        //shared_ptr<PhyloNode> father = tree_->getNode(fatherId);
        //uint fatherId = tree_->getNodeIndex(father);
        for (size_t i = 0; i < nbDistinctSites_; i++){
            size_t fatherState = (*ancestors_)[fatherId][i];
            //std::cerr <<  " ->      Father state is: "  << fatherState << std::endl;
            size_t currNodeState = (size_t)(nodeReconstruction(fatherState, i));
            //std::cerr <<  " ->      Current state is: "  << currNodeState << std::endl;
            ancestors.push_back(currNodeState);
        }
    }
    return ancestors;


}
/**************************************************************************************/
std::map<uint, std::vector<size_t> > JointMLAncestralReconstruction::getAllAncestralStates() const{
    if (ancestors_ == 0){
        throw Exception("getAllAncestralStates: Should use first init()");
    }
    return (*ancestors_);
}
/**************************************************************************************/
Sequence* JointMLAncestralReconstruction::getAncestralSequenceForNode(uint nodeId) const{
    string name = tree_->getNode(nodeId)->hasName() ? tree_->getNode(nodeId)->getName() : ("" + TextTools::toString(nodeId));
    vector<int> allStates(nbSites_);

    const auto& rootPatternLinks = likelihood_->getRootArrayPositions();

    const auto& statemap = likelihood_->getStateMap();

    auto states = getAncestralStatesForNode(nodeId);
    for (size_t i = 0; i < nbSites_; i++){
        allStates[i] = statemap.getAlphabetStateAsInt(states[rootPatternLinks(Eigen::Index(i))]);
    }
    return new BasicSequence(name, allStates, alphabet_);
}
/**********************************************************************************/
AlignedSequenceContainer* JointMLAncestralReconstruction::getAncestralSequences() const
{
  AlignedSequenceContainer* asc = new AlignedSequenceContainer(alphabet_);
  vector<shared_ptr<PhyloNode> > inNodes = tree_->getAllInnerNodes();
  for (size_t i = 0; i < inNodes.size(); i++)
  {
    Sequence* seq = getAncestralSequenceForNode(tree_->getNodeIndex(inNodes[i]));
    asc->addSequence(*seq);
    delete seq;
  }
  return asc;
}