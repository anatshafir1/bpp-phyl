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
    if (ChromEvolOptions::baseNum_ != IgnoreParam){
        if (ChromEvolOptions::baseNum_ > (int)chrRange_){
            chrRange_ = ChromEvolOptions::baseNum_ + 1;
        }
    }
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
        ChromEvolOptions::maxChrNum_ = maxNumberOfChr + std::abs(ChromEvolOptions::maxChrNum_);
    }else{
        if ((int)maxNumberOfChr > ChromEvolOptions::maxChrNum_){
            ChromEvolOptions::maxChrNum_ = maxNumberOfChr;
        }
    }

}
/****************************************************************************/
void ChromosomeNumberMng::setMinChrNum(unsigned int minNumberOfChr){
    if (ChromEvolOptions::minChrNum_ < 0){
        ChromEvolOptions::minChrNum_ = minNumberOfChr - std::abs(ChromEvolOptions::minChrNum_);
    }else{
        if ((int)minNumberOfChr < ChromEvolOptions::minChrNum_){
            ChromEvolOptions::minChrNum_ = minNumberOfChr;
        }

    }
}
/********************************************************************************************/

void ChromosomeNumberMng::getTree(const string& path, double treeLength){
    Newick reader;
    tree_ = reader.readPTree(path);
    double treeLengthToScale = (treeLength > 0) ? treeLength : (double) numberOfUniqueStates_;
    rescale_tree(tree_, treeLengthToScale);
    return;

}
/****************************************************************************/
void ChromosomeNumberMng::rescale_tree(PhyloTree* tree, double chrRange){
    double scale_tree_factor = 1.0;
    //string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    //std :: cout << tree_str << endl;
    bool rooted = tree->isRooted();
    if (!rooted){
        //throw UnrootedTreeException("The given input tree is unrooted. Tree must be rooted!", tree);
        throw Exception("The given input tree is unrooted. Tree must be rooted!\n");
    }
    if (ChromEvolOptions::branchMul_ == 1.0){
        return;
    }else{
        //tree must be rescaled
        double treeLength = tree->getTotalLength();
        if (ChromEvolOptions::branchMul_ == 999){
            scale_tree_factor = chrRange/treeLength;
        }
        tree->scaleTree(scale_tree_factor);

    }

}
/*****************************************************************************************/
void ChromosomeNumberMng::getMaxParsimonyUpperBound(double* parsimonyBound) const{
    Newick reader;
    TreeTemplate<Node>* tree = reader.readTree(ChromEvolOptions::treeFilePath_);
    double factor = tree_->getTotalLength()/tree->getTotalLength();
    tree->scaleTree(factor);
    DRTreeParsimonyScore maxParsimonyObject = DRTreeParsimonyScore(*tree, *vsc_);
    *parsimonyBound = (maxParsimonyObject.getScore())/(tree->getTotalLength());
    delete tree;
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
    //bool calculateDerivatives = true;
    //if (ChromEvolOptions::optimizationMethod_ == "Brent"){
        //calculateDerivatives  = false;
    //}
    unsigned int maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    ChromosomeNumberOptimizer* opt = new ChromosomeNumberOptimizer(tree_, alphabet_, vsc_, maxBaseNumTransition);
    opt->initModels(modelParams, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_);

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
void ChromosomeNumberMng::getJointMLAncestralReconstruction(ChromosomeNumberOptimizer* optimizer) const{
    vector<SingleProcessPhyloLikelihood*> vectorOfLikelihoods = optimizer->getVectorOfLikelihoods();
    // get the best likelihood
    SingleProcessPhyloLikelihood* lik = vectorOfLikelihoods[0];
    ValueRef <Eigen::RowVectorXd> rootFreqs = lik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    std::cout << "*** Root frequencies !!!! ****" << std::endl;
    auto rootFreqsValues =  rootFreqs->getTargetValue();
    Vdouble rootFreqsBpp;
    copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    vector <double> modelParams = getVectorToSetModelParams(lik);
    unsigned int maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelParams, maxBaseNumTransition, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    std::shared_ptr<SubstitutionModel> model(static_pointer_cast<SubstitutionModel>(chrModel)->clone());

    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqsBpp);
    std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    
    //const SubstitutionModel* modelRaw = dynamic_cast<const SubstitutionModel*>(lik->getLikelihoodCalculationSingleProcess()->getSubstitutionProcess().getModel(1));
    //std::shared_ptr<SubstitutionModel> model(modelRaw->clone());
    //ParameterList paramsUpdated = lik->getLikelihoodCalculationSingleProcess()->getParameters();
    //ParameterList params = model->getParameters();
    //NonHomogeneousSubstitutionProcess* subProSim;
    ParametrizablePhyloTree parTree(*tree_);
    auto subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone(), shared_ptr<FrequencySet>(rootFrequencies->clone()));
    //subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone());
    SubstitutionProcess* subProcess = subProSim->clone();
    Context context;
    auto likAncestralRec = std::make_shared<LikelihoodCalculationSingleProcess>(context, *vsc_->clone(), *subProcess, rootFreqs);
    ParameterList paramsUpdated = likAncestralRec->getParameters();
    likAncestralRec->makeJointMLAncestralReconstruction();
    JointMLAncestralReconstruction* ancr = new JointMLAncestralReconstruction(likAncestralRec);
    ancr->init();
    std::map<uint, std::vector<size_t>> ancestors = ancr->getAllAncestralStates();
    std::map<uint, std::vector<size_t>>::iterator it = ancestors.begin();
    std::cout <<"******* ******* ANCESTRAL RECONSTRUCTION ******* ********" << endl;
    while(it != ancestors.end()){
        uint nodeId = it->first;
        if(!(tree_->isLeaf(tree_->getNode(nodeId)))){
            cout << "   ----> N-" << nodeId <<" states are: " << endl;
            for (size_t s = 0; s < ancestors[nodeId].size(); s++){
                cout << "           state: "<< ancestors[nodeId][s] + alphabet_->getMin() << endl;
            }
        }else{
            cout << "   ----> " << (tree_->getNode(nodeId))->getName() << " states are: " << endl;
            for (size_t s = 0; s < ancestors[nodeId].size(); s++){
                cout << "           state: "<< ancestors[nodeId][s]+ alphabet_->getMin() << endl;
            }
        }
        it++;
    }
    const string outFilePath = ChromEvolOptions::resultsPathDir_ + "//" + "MLAncestralReconstruction.tree";
    PhyloTree* treeWithStates = tree_->clone();
    printTreeWithStates(*treeWithStates, ancestors, outFilePath);
    delete treeWithStates;


    delete ancr;
    //double likVal = likAncestralRec->makeJointMLAncestralReconstructionTest();
    std::cout << "********************************************\n";
    std::cout << " * * * * * * * * * * * * * * * * * * * * *\n";
    std::cout << "********************************************\n";
    //std::cout << "Ancestral reconstruction best for root is : " << likVal << endl;
    delete subProSim;
}
/***********************************************************************************/
vector <double> ChromosomeNumberMng::getVectorToSetModelParams(SingleProcessPhyloLikelihood* lik) const{
    vector <double> params;
    params.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        switch(i){
            case ChromosomeSubstitutionModel::BASENUM:
                if (ChromEvolOptions::baseNum_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::baseNum_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.baseNum_1")).getValue());
                }       
                break;
            case ChromosomeSubstitutionModel::BASENUMR:
                if (ChromEvolOptions::baseNumR_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::baseNumR_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.baseNumR_1")).getValue());
                }          
                break;
            case ChromosomeSubstitutionModel::DUPL:
                if (ChromEvolOptions::constDupl_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::constDupl_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.dupl_1")).getValue());
                }
                break;
            case ChromosomeSubstitutionModel::LOSS:
                if (ChromEvolOptions::constLoss_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::constLoss_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.loss_1")).getValue());
                }
                break;
            case ChromosomeSubstitutionModel::GAIN:
                if (ChromEvolOptions::constGain_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::constGain_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.gain_1")).getValue());
                }
                break;
            case ChromosomeSubstitutionModel::DEMIDUPL:
                if (ChromEvolOptions::constDemiDupl_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::constDemiDupl_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.demi_1")).getValue());
                }
                break;
            case ChromosomeSubstitutionModel::LOSSR:
                if (ChromEvolOptions::lossR_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::lossR_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.lossR_1")).getValue());
                }
                break;
            case ChromosomeSubstitutionModel::GAINR:
                if (ChromEvolOptions::gainR_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::gainR_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.gainR_1")).getValue());
                }
                break;
            case ChromosomeSubstitutionModel::DUPLR:
                if (ChromEvolOptions::duplR_ == IgnoreParam){
                    params.push_back(ChromEvolOptions::duplR_);
                }else{
                    params.push_back((lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.duplR_1")).getValue());
                }
                break;
            default:
                throw Exception("ChromosomeNumberMng::getVectorToSetModelParams(): Invalid rate type!");
                break;
        }

    }
    return params; 


}

/***********************************************************************************/
void ChromosomeNumberMng::runChromEvol(){
    // if (ChromEvolOptions::simulateData_){
    //     //simulate data using a tree and a set of model parameters
    //     RandomTools::setSeed(static_cast<long>(ChromEvolOptions::seed_));
    //     simulateData();
    //     if (ChromEvolOptions::numOfDataToSimulate_ > 1){
    //         return;
    //     }

    // }
    // optimize likelihood
    ChromosomeNumberOptimizer* chrOptimizer = optimizeLikelihoodMultiStartPoints();
    // get joint ML ancestral reconstruction
    getJointMLAncestralReconstruction(chrOptimizer);
    //get Marginal ML ancestral reconstruction, and with the help of them- calculate expectations of transitions
    const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+ "ancestorsProbs.txt";
    getMarginalAncestralReconstruction(chrOptimizer, outFilePath);
    
    //compute expectations

    delete chrOptimizer;


}
/**************************************************************************************/
void ChromosomeNumberMng::getMarginalAncestralReconstruction(ChromosomeNumberOptimizer* chrOptimizer, const string &filePath){
    vector<SingleProcessPhyloLikelihood*> vectorOfLikelihoods = chrOptimizer->getVectorOfLikelihoods();
    // get the best likelihood
    SingleProcessPhyloLikelihood* lik = vectorOfLikelihoods[0];
    auto singleLikProcess = lik->getLikelihoodCalculationSingleProcess();
    vector<shared_ptr<PhyloNode> > nodes = tree_->getAllNodes();
    size_t nbNodes = nodes.size();
    MarginalAncestralReconstruction *asr = new MarginalAncestralReconstruction(singleLikProcess);
    std::map<uint, VVdouble> posteriorProbs;
    std::map<uint, vector<size_t>> mapOfAncestors;
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        posteriorProbs[nodeId].reserve(1);//one site
        mapOfAncestors[nodeId] = asr->getAncestralStatesForNode(nodeId, posteriorProbs[nodeId], false); 
    }
    ofstream outFile;
    outFile.open(filePath);
    outFile << "NODE";
    for (size_t i = 0; i < alphabet_->getSize(); i ++){
        outFile << "\t" << (i + alphabet_->getMin());
    }
    outFile <<"\n";
    std::map<uint, std::vector<size_t>>::iterator it = mapOfAncestors.begin();
    while(it != mapOfAncestors.end()){
        uint nodeId = it->first;
        if(!(tree_->isLeaf(tree_->getNode(nodeId)))){
            outFile << "N-" << nodeId;
        }else{
            outFile << (tree_->getNode(nodeId))->getName();
        }
        for (size_t i = 0; i < posteriorProbs[nodeId][0].size(); i ++){
            outFile << "\t" << (posteriorProbs[nodeId][0][i]);

        }
        outFile << "\n";

        it++;
    }

    outFile.close();
    const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+"MarginalAncestralReconstruction.tree";
    printTreeWithStates(*tree_, mapOfAncestors, outFilePath);
    delete asr;
}
/**************************************************************************************/
/* void ChromosomeNumberMng::printPosteriorProbNodes(std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, vector<double>& rootPosterior) const{
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        throw Exception("Error in ChromosomeNumberMng::printPosteriorProbNodes(): No results file path!\n");
    }
    const string outPath = ChromEvolOptions::resultsPathDir_+"//"+ "ancestorsProbs.txt";
    ofstream outFile;
    outFile.open(outPath);
    vector<int> nodesIds = tree_->getNodesId();
    outFile <<"NODE\t";
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        (i < alphabet_->getSize()-1) ? (outFile << (int)i + alphabet_->getMin() << "\t") : (outFile << (int)i + alphabet_->getMin() <<"\n");
    }
    for (size_t n = 0; n < nodesIds.size(); n++){
        int nodeId = nodesIds[n];
        string nodeName;
        if (tree_->isLeaf(nodeId)){
            nodeName = tree_->getNodeName(nodeId);

        }else{
            nodeName = "N" + std::to_string(nodeId);
        }
        outFile << nodeName <<"\t";
        if (tree_->getRootId() == nodeId){
            for (size_t state = 0; state < alphabet_->getSize(); state++){
                (state < alphabet_->getSize()-1) ? (outFile << rootPosterior[state] << "\t") : (outFile << rootPosterior[state] <<"\n");
            }
            continue;
        }
        
        for (size_t son = 0; son < alphabet_->getSize(); son++){
            double posteriorProb = 0;
            for (size_t father = 0; father < alphabet_->getSize(); father++){
                posteriorProb += jointProbabilitiesFatherSon[nodeId][0][son][father];

            }
            (son == alphabet_->getSize()-1) ? (outFile << posteriorProb <<"\n") : (outFile << posteriorProb <<"\t");          
        }
        
    }
    outFile.close();

} */
/**************************************************************************************/
/* void ChromosomeNumberMng::printSimulatedEvoPath(TreeTemplate<Node> tree, const string outPath, RASiteSimulationResult* simResult) const{
    ofstream outFile;
    outFile.open(outPath);
    size_t totalNumTransitions = 0;
    vector<int> nodesIds = tree.getNodesId();
    for (size_t n = 0; n < nodesIds.size(); n++){
        if (tree.isRoot(nodesIds[n])){
            outFile << "N-" + to_string(nodesIds[n]) << endl;
            outFile <<"\tThe root state is: "<< ((int)(simResult->getRootAncestralState())+ alphabet_->getMin()) <<endl;
        }else{
            if (tree.isLeaf(nodesIds[n])){
                outFile << tree.getNodeName(nodesIds[n]) << endl;
            }else{
                outFile << "N-" + to_string(nodesIds[n]) <<endl;

            }
            MutationPath mutPath = simResult->getMutationPath(nodesIds[n]);
            vector<size_t> states = mutPath.getStates();
            vector<double> times = mutPath.getTimes();
            totalNumTransitions += static_cast<int>(times.size());
            for (size_t i = 0; i < states.size(); i++){
                outFile <<"\tt = "<<times[i] << " to state = "<< ((int)(states[i]) + alphabet_->getMin()) << endl;

            }
            outFile <<"# Number of transitions per branch: "<< times.size() <<endl;   
            
        }
        
        outFile <<"*************************************"<<endl;
        
    }
    outFile <<"Total number of transitions is: "<< totalNumTransitions << endl;
    outFile.close();

} */
void ChromosomeNumberMng::printTreeWithStates(PhyloTree tree, std::map<uint, std::vector<size_t>> &ancestors, const string &filePath) const{
    uint rootId = tree.getRootIndex();
    convertNodesNames(tree, rootId, ancestors);
    string tree_str = printTree(tree);
    cout << tree_str << endl;
    if (filePath != "none"){
       ofstream outFile;
       outFile.open(filePath);
       outFile << tree_str << endl;
       outFile.close();
    }

}
/**************************************************************************************/
void ChromosomeNumberMng::convertNodesNames(PhyloTree &tree, uint nodeId, std::map<uint, std::vector<size_t>> &ancestors) const{
    size_t state = ancestors[nodeId][0] + alphabet_->getMin();
    if (tree.isLeaf(nodeId)){
        string prevName = tree.getNode(nodeId)->getName();
        const string newName = (prevName + "-"+ std::to_string(state));
        tree.getNode(nodeId)->setName(newName);

    }else{
        // internal node -> N[nodeId]-[state]
        string prevName = "N" + std::to_string(nodeId);
        const string newName = (prevName + "-"+ std::to_string(state));
        tree.getNode(nodeId)->setName(newName);
        auto sons = tree.getSons(tree.getNode(nodeId));
        //auto sons = tree.getNode(nodeId)->getSons();
        for (size_t i = 0; i < sons.size(); i++){
            uint sonId = tree.getNodeIndex(sons[i]);
            convertNodesNames(tree, sonId, ancestors);

        }
    }
}


/****************************************************************************************/
string ChromosomeNumberMng::printTree(const PhyloTree& tree)
{
  ostringstream s;
  s << "(";
  uint rootId = tree.getRootIndex();
  auto node = tree.getNode(rootId);
  if (tree.isLeaf(rootId) && node->hasName()) // In case we have a tree like ((A:1.0)); where the root node is an unamed leaf!
  {
    s << node->getName();
    auto sons = tree.getSons(node);
    for (size_t i = 0; i < sons.size(); ++i)
    {
        uint sonId = tree.getNodeIndex(sons[i]);
        s << "," << nodeToParenthesis(sonId, tree);
    }
  }
  else
  {
    auto sons = tree.getSons(node);
    uint firstSonId = tree.getNodeIndex(sons[0]);
    s << nodeToParenthesis(firstSonId, tree);
    for (size_t i = 1; i < sons.size(); ++i)
    {
        uint sonId = tree.getNodeIndex(sons[i]);
        s << "," << nodeToParenthesis(sonId, tree);
    }
  }
  s << ")";
  s << tree.getNode(rootId)->getName();

  s << ";" << endl;
  return s.str();
}

/******************************************************************************/
string ChromosomeNumberMng::nodeToParenthesis(const uint nodeId, const PhyloTree &tree)
{
  ostringstream s;
  if (tree.isLeaf(tree.getNode(nodeId)))
  {
    s << tree.getNode(nodeId)->getName();
  }
  else
  {
    s << "(";
  
    auto sons = tree.getSons(tree.getNode(nodeId));
    uint firstSonId = tree.getNodeIndex(sons[0]);
    s << nodeToParenthesis(firstSonId, tree);
    for (size_t i = 1; i < sons.size(); i++)
    {
        uint sonId = tree.getNodeIndex(sons[i]);
        
        s << "," << nodeToParenthesis(sonId, tree);
    }
    s << ")";
  }
  if (! tree.isLeaf(tree.getNode(nodeId))){
      s << tree.getNode(nodeId)->getName();
  }
  shared_ptr<PhyloBranch> branch=tree.getEdgeToFather(nodeId);
  if (branch->hasLength()){
    s << ":" << branch->getLength();

  }

  return s.str();
}
/*********************************************************************************/
/* void ChromosomeNumberMng::simulateData(){
    if ((ChromEvolOptions::minChrNum_ <= 0) || (ChromEvolOptions::maxChrNum_ < 0)){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): minimum and maximum chromsome number should be positive!");
    }
    if (ChromEvolOptions::maxChrNum_ <= ChromEvolOptions::minChrNum_){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): maximum chromsome number should be larger than minimum chromosome number!");
    }
    alphabet_ = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
    vector<double> modelParams;
    modelParams.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    ChromEvolOptions::initVectorOfChrNumParameters(modelParams);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    SubstitutionModelSet* modelSet = new SubstitutionModelSet(alphabet_);
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(alphabet_, modelParams, (unsigned int)ChromEvolOptions::maxBaseNumTransition_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL,  ChromEvolOptions::rateChangeType_);
    vector <int> nodeIds = tree_->getNodesId();
    nodeIds.pop_back();
    modelSet->addModel(chrModel, nodeIds);
    if (ChromEvolOptions::fixedFrequenciesFilePath_ == "none"){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): You need to path the file of fixed root frequencies!!");
    }
    ChromosomeNumberOptimizer::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, modelSet);

    for (size_t i = 0; i < (size_t)ChromEvolOptions::numOfDataToSimulate_; i++){
        NonHomogeneousSequenceSimulator* sim = new NonHomogeneousSequenceSimulator(modelSet, rdist, tree_);
        RASiteSimulationResult* simResult = sim->dSimulateSite();
        vector <size_t> leavesStates = simResult->getFinalStates();
        vector<string> leavesNames = simResult->getLeaveNames();
        printSimulatedData(leavesStates, leavesNames, i);
        printSimulatedDataAndAncestors(simResult);
        if (ChromEvolOptions::resultsPathDir_ != "none"){
            printSimulatedEvoPath(*tree_, ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedEvolutionPaths.txt", simResult);
        }
        delete simResult;
        delete sim;

    }
    delete modelSet;
    delete rdist;

} */
/*******************************************************************************/
/* void ChromosomeNumberMng::printSimulatedData(vector<size_t> leavesStates, vector<string> leavesNames, size_t iter){
    cout << "Simulated data #" << iter << endl;
    for (size_t i = 0; i < leavesNames.size(); i++){
        cout << leavesNames[i] << " "<< leavesStates[i] + alphabet_->getMin() <<endl;
    }
    cout << "******************************"<<endl;
    
    if (ChromEvolOptions::resultsPathDir_ != "none"){
        //create vector site container object and save fasta file.
        VectorSiteContainer* simulatedData = new VectorSiteContainer(alphabet_);
        for (size_t i = 0; i < leavesNames.size(); i++){
            int state = (int)leavesStates[i] + alphabet_->getMin();
            BasicSequence seq = BasicSequence(leavesNames[i], alphabet_->intToChar(state), static_cast <const Alphabet*>(alphabet_));
            simulatedData->addSequence(seq);
        }
        vsc_ = simulatedData;
        string pathForSimulatedData = ChromEvolOptions::resultsPathDir_ + "//"+ "chr_counts"+ to_string(static_cast<int>(iter)) +".fasta";
        Fasta fasta;
        fasta.writeSequences(pathForSimulatedData, *simulatedData);

    }


    
} */
/****************************************************************************/
/* void ChromosomeNumberMng::printSimulatedDataAndAncestors(RASiteSimulationResult* simResult) const{
    std::map<int, std::vector<size_t> > ancestors;
    vector<int> nodesIds = tree_->getNodesId();
    for (size_t i = 0; i < nodesIds.size(); i++){
        vector<size_t> nodesStates;
        if (tree_->isRoot(nodesIds[i])){
            nodesStates.push_back(simResult->getRootAncestralState());
        }else{
            nodesStates.push_back(simResult->getAncestralState(nodesIds[i]));
        }
        ancestors[nodesIds[i]] = nodesStates;
    }
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        printTreeWithStates(*tree_, ancestors, ChromEvolOptions::resultsPathDir_);
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedDataAncestors.tree";
        printTreeWithStates(*tree_, ancestors, outFilePath);
    }
  
} */