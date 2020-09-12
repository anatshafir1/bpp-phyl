#include "ChromosomeNumberMng.h"
#include "ChromEvolOptions.h"

using namespace bpp;

void ChromosomeNumberMng::getCharacterData (const string& path){
    ChromosomeAlphabet* alphaInitial = new ChromosomeAlphabet(ChromEvolOptions::minAlpha_, ChromEvolOptions::maxAlpha_);
    VectorSequenceContainer* initialSetOfSequences = chrFasta::readSequencesFromFile(path, alphaInitial);
    size_t numOfSequences = initialSetOfSequences->getNumberOfSequences();
    vector <string> sequenceNames = initialSetOfSequences->getSequencesNames();

    unsigned int maxNumberOfChr = 1; //the minimal number of chromosomes cannot be zero
    unsigned int minNumOfChr = ChromEvolOptions::maxAlpha_;

    std::vector <int> UniqueCharacterStates;
    cout<<"vector size is "<< UniqueCharacterStates.size()<<endl;
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = initialSetOfSequences->getSequence(sequenceNames[i]);
        int character = seq.getValue(0);
        if (character == -1){
            continue;
        }
        if (character == static_cast<int>(ChromEvolOptions::maxAlpha_)+1){
            continue;
        }
        // if it is a composite state
        if (character > static_cast<int>(ChromEvolOptions::maxAlpha_) +1){
            const std::vector<int> compositeCharacters = alphaInitial->getSetOfStatesForAComposite(character);
            for (size_t j = 0; j < compositeCharacters.size(); j++){
                if ((unsigned int) compositeCharacters[j] > maxNumberOfChr){
                    maxNumberOfChr = compositeCharacters[j];
                }
                if ((unsigned int) compositeCharacters[j] < minNumOfChr){
                    minNumOfChr = compositeCharacters[j];
                }
                
            }
            continue;
        }

        if (!std::count(UniqueCharacterStates.begin(), UniqueCharacterStates.end(), character)){
            UniqueCharacterStates.push_back(character);
        }
        if ((unsigned int) character > maxNumberOfChr){
            maxNumberOfChr = character;
        }
        if ((unsigned int) character < minNumOfChr){
            minNumOfChr = character;
        }

    }
    numberOfUniqueStates_ = (unsigned int)UniqueCharacterStates.size() + alphaInitial->getNumberOfCompositeStates();
    chrRange_ = maxNumberOfChr - minNumOfChr;
    cout <<"Number of unique states is " << numberOfUniqueStates_ <<endl;

    setMaxChrNum(maxNumberOfChr);
    setMinChrNum(minNumOfChr);

    vsc_ = resizeAlphabetForSequenceContainer(initialSetOfSequences, alphaInitial);
    delete initialSetOfSequences;
    delete alphaInitial;
    return;
}
// /*******************************************************************************************************************/
VectorSiteContainer* ChromosomeNumberMng::resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* alphaInitial){
    size_t numOfSequences = vsc->getNumberOfSequences();
    vector <string> sequenceNames = vsc->getSequencesNames();
    alphabet_ = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
        // fill with composite values
    if (alphaInitial->getNumberOfCompositeStates() > 0){
        const std::map <int, std::map<int, double>> compositeStates = alphaInitial->getCompositeStatesMap();
        std::map <int, std::map<int, double>>::const_iterator it = compositeStates.begin();
        while (it != compositeStates.end()){
            int compositeState = it->first;
            std::string charComposite = alphaInitial->intToChar(compositeState);
            alphabet_->setCompositeState(charComposite);
            it++;
        }
    }
    VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(alphabet_);
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = vsc->getSequence(sequenceNames[i]);
        BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), alphabet_);
        resized_alphabet_site_container->addSequence(new_seq);

    }
    return resized_alphabet_site_container;
}
/*******************************************************************************************/
void ChromosomeNumberMng::setMaxChrNum(unsigned int maxNumberOfChr){
    if (ChromEvolOptions::maxChrNum_ < 0){
        ChromEvolOptions::maxChrNum_ = maxNumberOfChr + abs(ChromEvolOptions::maxChrNum_);
    }else{
        if ((int)maxNumberOfChr > ChromEvolOptions::maxChrNum_){
            ChromEvolOptions::maxChrNum_ = maxNumberOfChr;
        }
    }

}
/****************************************************************************/
void ChromosomeNumberMng::setMinChrNum(unsigned int minNumberOfChr){
    if (ChromEvolOptions::minChrNum_ < 0){
        ChromEvolOptions::minChrNum_ = minNumberOfChr - abs(ChromEvolOptions::minChrNum_);
    }else{
        if ((int)minNumberOfChr < ChromEvolOptions::minChrNum_){
            ChromEvolOptions::minChrNum_ = minNumberOfChr;
        }

    }
}
/********************************************************************************************/

void ChromosomeNumberMng::getTree(const string& path){
    Newick newick;
    tree_ = newick.readTree(path);
    rescale_tree(tree_, numberOfUniqueStates_);
    return;

}
/****************************************************************************/
void ChromosomeNumberMng::rescale_tree(TreeTemplate<Node>* tree, unsigned int chrRange){
    double scale_tree_factor = 1.0;
    //string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    //std :: cout << tree_str << endl;
    bool rooted = tree->isRooted();
    if (!rooted){
        throw UnrootedTreeException("The given input tree is unrooted. Tree must be rooted!", tree);
    }
    if (ChromEvolOptions::branchMul_ == 1.0){
        return;
    }else{
        //tree must be rescaled
        double treeLength = tree->getTotalLength();
        if (ChromEvolOptions::branchMul_ == 999){
            scale_tree_factor = (double)chrRange/treeLength;
        }
        tree->scaleTree(scale_tree_factor);

    }

}
/*****************************************************************************************/
void ChromosomeNumberMng::getMaxParsimonyUpperBound(double* parsimonyBound) const{

    TreeTemplate<Node>* treeForParsimonyBound = tree_->clone();
    DRTreeParsimonyScore maxParsimonyObject = DRTreeParsimonyScore(*treeForParsimonyBound, *vsc_);
    *parsimonyBound = (maxParsimonyObject.getScore())/(tree_->getTotalLength());
    delete treeForParsimonyBound;
    return;   

}
/*****************************************************************************************/
ChromosomeNumberOptimizer* ChromosomeNumberMng::optimizeLikelihoodMultiStartPoints() const{
    
    vector<double> modelParams;
    modelParams.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    ChromEvolOptions::initVectorOfChrNumParameters(modelParams);
    double parsimonyBound = 0;
    if (ChromEvolOptions::maxParsimonyBound_){
        getMaxParsimonyUpperBound(&parsimonyBound);
    }
    bool calculateDerivatives = true;
    if (ChromEvolOptions::optimizationMethod_ == "Brent"){
        calculateDerivatives  = false;
    }

    ChromosomeNumberOptimizer* opt = new ChromosomeNumberOptimizer(tree_, alphabet_, vsc_, chrRange_);
    opt->initModels(modelParams, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], calculateDerivatives, ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_);

    //initialize all the optimization specific parameters
    opt->initOptimizer(ChromEvolOptions::OptPointsNum_, ChromEvolOptions::OptIterNum_, ChromEvolOptions::optimizationMethod_, ChromEvolOptions::baseNumOptimizationMethod_,
        ChromEvolOptions::tolerance_, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
        ChromEvolOptions::probsForMixedOptimization_);
    //optimize models
    opt->optimize();
    // it is safe to delete the chrOptimizer, because the destructor doesn't delete nothing associated with the vector of likelihoods
    return opt;
       
}
/******************************************************************************************************/
void ChromosomeNumberMng::getJointMLAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik) const{
    std::cout << "ML Joint Ancestral Reconstruction"<<endl;
    TransitionModel* model = lik->getSubstitutionModelSet()->getModel(0);
    std::vector<double> rootFreqs = lik->getRootFrequencies(0);
    std::map<int, VVVdouble> Pijt;
    TreeTemplate<Node> tree = lik->getTree();
    vector <int> nodesIds = tree.getNodesId();
    for (size_t n = 0; n < nodesIds.size(); n++){
        Pijt[nodesIds[n]] = lik->getTransitionProbabilitiesPerRateClass(nodesIds[n], 0);
    }
    MLAncestralStateReconstruction* ancr = new MLAncestralStateReconstruction(lik, model, rootFreqs, &Pijt);
    ancr->computeJointLikelihood();
    std::map<int, std::vector<size_t> > ancestors = ancr->getAllAncestralStates();
    printTreeWithStates(lik, lik->getTree(), ancestors);

    delete ancr;

}
/**************************************************************************************************************/
std::map<int, std::map<size_t, VVdouble>> ChromosomeNumberMng::getMarginalAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik) const{
    std::cout << "Marginal Ancestral Reconstruction"<<endl;
    MarginalNonRevAncestralStateReconstruction* ancr = new MarginalNonRevAncestralStateReconstruction(lik);
    ancr->computePosteriorProbabilitiesOfNodesForEachStatePerSite();
    std::map<int, std::vector<size_t> > ancestors = ancr->getAllAncestralStates();
    std::map<int, map<size_t, std::vector<double>>>* probs = ancr->getPosteriorProbForAllNodesAndStatesPerSite();
    printTreeWithStates(lik, lik->getTree(), ancestors, probs);
    std::map<int, std::map<size_t, VVdouble>> jointProbabilitiesFatherSon = ancr->getAllJointFatherNodeProbabilities();
    delete ancr;
    return jointProbabilitiesFatherSon;

}
/**********************************************************************************************************************/
void ChromosomeNumberMng::computeExpectations(DRNonHomogeneousTreeLikelihood* lik, std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int numOfSimulations) const{

    ComputeChromosomeTransitionsExp* expCalculator = new ComputeChromosomeTransitionsExp(lik, jointProbabilitiesFatherSon, ChromEvolOptions::jumpTypeMethod_);
    expCalculator->runSimulations(numOfSimulations);
    expCalculator->computeExpectationPerType();
    expCalculator->printResults();

    //delete
    delete expCalculator;
}
/***********************************************************************************/
void ChromosomeNumberMng::runChromEvol() const{
    // optimize likelihood
    ChromosomeNumberOptimizer* chrOptimizer = optimizeLikelihoodMultiStartPoints();
    std::vector <DRNonHomogeneousTreeLikelihood> lik_vec = chrOptimizer->getVectorOfLikelihoods();
    // get joint ML ancestral reconstruction
    getJointMLAncestralReconstruction(&lik_vec[0]);
    //get Marginal ML ancestral reconstruction, and with the help of them- calculate expectations of transitions
    std::map<int, std::map<size_t, VVdouble>>  jointProbabilitiesFatherSon = getMarginalAncestralReconstruction(&lik_vec[0]);
    //compute expectations
    computeExpectations(&lik_vec[0], jointProbabilitiesFatherSon, ChromEvolOptions::NumOfSimulations_);
    //Clear the vector of likelihoods entirely
    delete chrOptimizer;


}
/**************************************************************************************/
void ChromosomeNumberMng::printTreeWithStates(DRNonHomogeneousTreeLikelihood* lik, TreeTemplate<Node> tree, std::map<int, std::vector<size_t> > ancestors, std::map<int, map<size_t, std::vector<double>>>* probs){
    map <string, double> mapOfNodeNameProb;
    vector<int> nodesIds = tree.getNodesId();
    for (size_t n= 0; n < nodesIds.size(); n++){
        for (size_t i = 0; i < ancestors[nodesIds[0]].size(); i++){
            size_t state = ancestors[nodesIds[n]][i] + (dynamic_cast<const ChromosomeAlphabet*>(lik->getAlphabet()))->getMin();
            string prevName;
            
            if (tree.isLeaf(nodesIds[n])){
                prevName = tree.getNodeName(nodesIds[n]);
                const string newName = (prevName + "-"+ std::to_string(state));
                tree.setNodeName(nodesIds[n], newName);
                if (probs){
                    double postProb = (*probs)[nodesIds[n]][i][ancestors[nodesIds[n]][i]];
                    cout <<newName<< " state index: " << state << " : " << postProb <<endl;
                    mapOfNodeNameProb[newName] = postProb;

                }else{
                    cout <<newName<< " state index: " << state <<endl;
                }
            }else{
                prevName = "N" + std::to_string(nodesIds[n]);
                const string newName = (prevName + "-"+ std::to_string(state));
                tree.setNodeName(nodesIds[n], newName);
                if (probs){
                    double postProb = (*probs)[nodesIds[n]][i][ancestors[nodesIds[n]][i]];
                    cout <<newName<< " state index: " << state << " : " << postProb <<endl;
                    mapOfNodeNameProb[newName] = postProb;

                }else{
                    cout <<newName<< " state index: " << state <<endl;
                }
            }
            
            
        }
        
    }
    string tree_str;
    if (probs){
        tree_str = printTree(tree, &mapOfNodeNameProb);
    }
    else{
        tree_str = printTree(tree);
    }
    cout << tree_str << endl;

}
/****************************************************************************************/
string ChromosomeNumberMng::printTree(const TreeTemplate<Node>& tree, map<string, double>* mapNameProb)
{
  ostringstream s;
  s << "(";
  const Node* node = tree.getRootNode();
  if (node->isLeaf() && node->hasName()) // In case we have a tree like ((A:1.0)); where the root node is an unamed leaf!
  {
    s << node->getName();
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), mapNameProb);
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0), mapNameProb);
    for (size_t i = 1; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i), mapNameProb);
    }
  }
  s << ")";
  s << tree.getRootNode()->getName();
  if (! mapNameProb){
    if (node->hasDistanceToFather()){
        s << ":" << node->getDistanceToFather();

    }
        
  }else{
      s << ":" << (*mapNameProb)[node->getName()];
  }

 
  
  s << ";" << endl;
  return s.str();
}

/******************************************************************************/
string ChromosomeNumberMng::nodeToParenthesis(const Node& node, map<string, double>* mapNameProb)
{
  ostringstream s;
  if (node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(*node[0], mapNameProb);
    for (int i = 1; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      s << "," << nodeToParenthesis(*node[i], mapNameProb);
    }
    s << ")";
  }
  if (! node.isLeaf()){
      s << node.getName();
  }
  if (node.hasDistanceToFather()){
    if (!mapNameProb){
        s << ":" << node.getDistanceToFather();
    }else{
        s << ":" << (*mapNameProb)[node.getName()];

    }
  }

  return s.str();
}
