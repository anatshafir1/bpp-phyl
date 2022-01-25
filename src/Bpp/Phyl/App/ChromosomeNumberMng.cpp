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
    uint chrRangeNum = maxNumberOfChr - minNumOfChr;
    for (uint j = 1; j <= static_cast<uint>(ChromEvolOptions::numOfModels_); j++){      
        if (ChromEvolOptions::baseNum_[j] != IgnoreParam){
            if (ChromEvolOptions::baseNum_[j] > (int)chrRangeNum){
                chrRange_[j] = ChromEvolOptions::baseNum_[j] + 1;
            }else{
                chrRange_[j] = chrRangeNum;
            }
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
/*************************************************************************************************************/
void ChromosomeNumberMng::setNodeIdsForAllModels(string &path){
    if (path == "none"){
        auto nodes = tree_->getAllNodes();
        for (size_t i = 0; i < nodes.size(); i++){
            uint nodeId = tree_->getNodeIndex(nodes[i]);
            if (nodeId == tree_->getRootIndex()){
                continue;
            }else{
                ChromEvolOptions::mapModelNodesIds_[1].push_back(nodeId);
            }
        }
        return;
    }
    ifstream stream;
    stream.open(path.c_str());
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    PhyloTree* tree = tree_->clone();
    std::map<uint, std::pair<uint, std::vector<uint>>> mapOfModelMRCAAndNodes;
    std::map<uint, uint> mapNodeModel;
    for (size_t i = 0; i < lines.size(); i ++){
        if (lines[i] == ""){
            continue;
        }
        getNodeIdsPerModelFromLine(lines[i], tree, mapOfModelMRCAAndNodes);

    }
    auto it_ModelNodes = mapOfModelMRCAAndNodes.begin();
    while (it_ModelNodes != mapOfModelMRCAAndNodes.end()){
        uint model = it_ModelNodes->first;
        uint nodeId = mapOfModelMRCAAndNodes[model].first;
        mapNodeModel[nodeId] = model;
        ChromEvolOptions::mapModelNodesIds_[model] = mapOfModelMRCAAndNodes[model].second;
        it_ModelNodes ++;
    }
    auto it = mapOfModelMRCAAndNodes.begin();
    while (it != mapOfModelMRCAAndNodes.end()){
        uint model = it->first;
        vector<uint> nodeIds = mapOfModelMRCAAndNodes[model].second;
        for (size_t i = 0; i < nodeIds.size(); i++){
            auto itNodeModel = mapNodeModel.find(nodeIds[i]);
            if (itNodeModel != mapNodeModel.end()){
                uint modelOfDescendant = mapNodeModel[nodeIds[i]];
                if (modelOfDescendant != model){
                    vector<uint> subtree = mapOfModelMRCAAndNodes[modelOfDescendant].second;
                    for (size_t j = 0; j < subtree.size(); j++){
                        auto nodeToDelIt = std::find(ChromEvolOptions::mapModelNodesIds_[model].begin(), ChromEvolOptions::mapModelNodesIds_[model].end(), subtree[j]);
                        if (nodeToDelIt != ChromEvolOptions::mapModelNodesIds_[model].end()){
                            ChromEvolOptions::mapModelNodesIds_[model].erase(std::remove(ChromEvolOptions::mapModelNodesIds_[model].begin(), ChromEvolOptions::mapModelNodesIds_[model].end(), subtree[j]),ChromEvolOptions::mapModelNodesIds_[model].end());

                        }
                    }
                }
            }
        }
        it ++;
    }

    delete tree;


}
/**************************************************************************************************************/
// shared_ptr<PhyloNode> ChromosomeNumberMng::getMRCA(PhyloTree* tree, std::vector<shared_ptr<PhyloNode>> nodes){
//     shared_ptr<PhyloNode> mrca;
//     auto nodesInIndices = tree->getNodeIndexes(nodes);
//     size_t numOfFound = 0;

//     for (size_t i = 0; i < nodes.size(); i++){
//         if (tree->isLeaf(nodes[i])){
//             continue;
//         }
//         auto nodesUnderSubtree = tree->getSubtreeNodes(nodes[i]);
//         auto subtreeNodesIndices = tree->getNodeIndexes(nodesUnderSubtree);
//         for (size_t j = 0; j < subtreeNodesIndices.size(); j++){
//             auto it = std::find(nodesInIndices.begin(), nodesInIndices.end(), subtreeNodesIndices[j]);
//             if (it != nodesInIndices.end()){
//                 numOfFound ++;
//                 if (numOfFound == nodesInIndices.size()){
//                     mrca = nodes[i];
//                     return mrca;
//                 }
//             }
//         }
//     }
//     // MRCA was not among the nodes
//     // choose for example the first node
//     uint nodeId = nodesInIndices[0];
    
//     while (numOfFound < nodesInIndices.size()){
//         numOfFound = 0;
//         if (nodeId == tree->getRootIndex()){
//             mrca = tree->getRoot();
//             break;
//         }
//         auto edgeIndex =  tree->getIncomingEdges(nodeId)[0]; 
//         nodeId = tree->getFatherOfEdge(edgeIndex);
//         auto nodesOfSubtree = tree->getSubtreeNodes(tree->getNode(nodeId));
//         for (size_t i = 0; i < nodesOfSubtree.size(); i++){
//             auto subtreeNodeId = tree->getNodeIndex(nodesOfSubtree[i]);
//             auto it = std::find(nodesInIndices.begin(), nodesInIndices.end(), subtreeNodeId);
//             if (it != nodesInIndices.end()){
//                 numOfFound ++;
//                 mrca = tree->getNode(nodeId);
//             }
//         }

//     }
//     return mrca;
// }
/**************************************************************************************************************/
void ChromosomeNumberMng::getNodeIdsPerModelFromLine(string &content, PhyloTree* tree, std::map<uint, std::pair<uint, std::vector<uint>>> &modelAndNodeIds){
    vector<string> paramValues;
    std::regex modelPattern ("([\\d]+)");
    std::regex treePattern ("\\(([\\S]+)\\)");
    StringTokenizer stoken = StringTokenizer(content, "=");
    while (stoken.hasMoreToken()){
        paramValues.push_back(stoken.nextToken());
    }
    shared_ptr<PhyloNode> mrca_node;
    uint model;
    vector<uint> nodes;
    for (size_t i = 0; i < paramValues.size(); i++){
        std::smatch sm;
        if (i == 0){
            std::regex_search(paramValues[i], sm, modelPattern);
            model = std::stoi(sm[0]);

        }else{
            std::regex_search(paramValues[i], sm, treePattern);
            string speciesNonSepWithBrackets = sm[0];
            string speciesNonSep = speciesNonSepWithBrackets.substr(1, speciesNonSepWithBrackets.length()-2);

            vector<string> speciesNames;
            StringTokenizer speciesToken = StringTokenizer(speciesNonSep, ",");
            while (speciesToken.hasMoreToken()){
                speciesNames.push_back(speciesToken.nextToken());
            }
            std::map<std::string, shared_ptr<PhyloNode>> subtreeLeavesAsNodes;
            vector<shared_ptr<PhyloNode>> allLeaves = tree->getAllLeaves();
            for (size_t k = 0; k < allLeaves.size(); k++){
                subtreeLeavesAsNodes[allLeaves[k]->getName()] = allLeaves[k];
            }
            vector<shared_ptr<PhyloNode>> leaveNodesForMrca;
            for (size_t k = 0; k < speciesNames.size(); k++){
                leaveNodesForMrca.push_back(subtreeLeavesAsNodes[speciesNames[k]]);
            }
            
            //mrca_node = tree->MRCA(leaveNodesForMrca);
            mrca_node = ChromEvolOptions::getMRCA(tree, leaveNodesForMrca);
            // DUBUG !!!! !!!! !!!!!!!//////////////////////////////

            /////////////////////////////////////////////////////////

            // just for meantime !!

            //////////// DEBUG ////////////////////////////
            // uint mrca_id = tree->getNodeIndex(mrca_node);
            // auto sons = tree->getSons(tree->getNode(tree->getRootIndex()));
            // auto mrca_should_beRoot = tree->MRCA(sons);
            // auto rootIndex = tree->getNodeIndex(mrca_should_beRoot);
            // std::cout << tree->getRootIndex() << std::endl;
            ///////////////////////////////////////////////
            auto mrca_id = tree->getNodeIndex(mrca_node);
            ChromEvolOptions::initialModelNodes_.push_back(mrca_id);
            auto subtreeNodes = tree->getSubtreeNodes(mrca_node);
            auto allNodeIds = tree->getNodeIndexes(subtreeNodes);
            for (size_t j = 0; j  < allNodeIds.size(); j++){
                if (allNodeIds[j] == tree->getRootIndex()){
                    continue;
                }
                nodes.push_back(allNodeIds[j]);
            }
            auto leavesUnderNode = tree->getLeavesUnderNode(mrca_node);
            std:: cout << "Model #" << model << std::endl;
            for (size_t j= 0; j < leavesUnderNode.size(); j++){
                std::cout << leavesUnderNode[j]->getName() << std::endl;
            }
        }    

    }
    modelAndNodeIds[model].first = tree->getNodeIndex(mrca_node);
    modelAndNodeIds[model].second = nodes;
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
        if (minNumberOfChr == 1){
            ChromEvolOptions::minChrNum_ = 0;
            std::cout << "Warning !!!! minChrNum_ should be at least 1!!" << std::endl;
            std::cout << "The mininal chromosome number was determined to be 1" << std::endl;
        }
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

        if ((ChromEvolOptions::branchMul_ == 999) || (ChromEvolOptions::treeLength_)){
            scale_tree_factor = chrRange/treeLength;
        }else{
            scale_tree_factor = ChromEvolOptions::branchMul_;
        }
        if (scale_tree_factor == 0){
            throw Exception("ChromosomeNumberMng::rescale_tree(): ERROR!!! Tree will be scaled to 0!!!!");
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
    std::map<uint, std::pair<int, std::map<int, vector<double>>>> complexParamsValues;
    ChromEvolOptions::getInitialValuesForComplexParams(complexParamsValues);
    
    double parsimonyBound = 0;
    if (ChromEvolOptions::maxParsimonyBound_){
        getMaxParsimonyUpperBound(&parsimonyBound);
    }
    //bool calculateDerivatives = true;
    //if (ChromEvolOptions::optimizationMethod_ == "Brent"){
        //calculateDerivatives  = false;
    //}
    vector<uint> numOfIterationsForBackward = ChromEvolOptions::OptIterNumNextRounds_;
    vector<uint> numOfPointsForForward = ChromEvolOptions::OptPointsNumNextRounds_;
    if (!ChromEvolOptions::forwardPhase_){
        ChromEvolOptions::OptIterNumNextRounds_ = {0};
        ChromEvolOptions::OptPointsNumNextRounds_ = {1};

    }
    std::map<uint, uint> maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    ChromosomeNumberOptimizer* opt = new ChromosomeNumberOptimizer(tree_, alphabet_, vsc_, maxBaseNumTransition);
    //opt->initModels(complexParamsValues, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_, ChromEvolOptions::mapModelNodesIds_);
    //initialize all the optimization specific parameters
    opt->initOptimizer(ChromEvolOptions::OptPointsNum_, ChromEvolOptions::OptIterNum_, ChromEvolOptions::OptPointsNumNextRounds_, ChromEvolOptions::OptIterNumNextRounds_, ChromEvolOptions::optimizationMethod_, ChromEvolOptions::baseNumOptimizationMethod_,
        ChromEvolOptions::tolerance_, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
        ChromEvolOptions::probsForMixedOptimization_);
    //optimize models
    //opt->optimize(complexParamsValues, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_ ,ChromEvolOptions::mapModelNodesIds_);
    time_t t1;
    time(&t1);
    time_t t2;
    if (ChromEvolOptions::parallelization_){
        opt->optimizeInParallel(complexParamsValues, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_ ,ChromEvolOptions::mapModelNodesIds_);

    }else{
        opt->optimize(complexParamsValues, parsimonyBound, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], ChromEvolOptions::fixedFrequenciesFilePath_, ChromEvolOptions::fixedParams_ ,ChromEvolOptions::mapModelNodesIds_);

    }
    ChromEvolOptions::OptIterNumNextRounds_ = numOfIterationsForBackward;
    ChromEvolOptions::OptPointsNumNextRounds_ = numOfPointsForForward;
    opt->setIterNumForNextRound(ChromEvolOptions::OptIterNumNextRounds_);
    opt->setPointsNumForNextRound(ChromEvolOptions::OptPointsNumNextRounds_);
    if (ChromEvolOptions::backwardPhase_){
        opt->optimizeBackwards(parsimonyBound, ChromEvolOptions::parallelization_);
    }


    time(&t2);
    std::cout <<"**** **** Total running time of the optimization procedure is: "<< (t2-t1) <<endl;
    //initialize all the optimization specific parameters

    
    // it is safe to delete the chrOptimizer, because the destructor doesn't delete nothing associated with the vector of likelihoods

    return opt;
       
}
/******************************************************************************************************/
void ChromosomeNumberMng::getJointMLAncestralReconstruction(ChromosomeNumberOptimizer* optimizer) const{
    vector<SingleProcessPhyloLikelihood*> vectorOfLikelihoods = optimizer->getVectorOfLikelihoods();
    // get the best likelihood
    SingleProcessPhyloLikelihood* lik = vectorOfLikelihoods[0];
    //ValueRef <Eigen::RowVectorXd> rootFreqs = lik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    //std::cout << "*** Root frequencies !!!! ****" << std::endl;
    //auto rootFreqsValues =  rootFreqs->getTargetValue();
    //Vdouble rootFreqsBpp;
    //copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
    //DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    std::map<int, vector<pair<uint, int>>> sharedParams = optimizer->getSharedParams();
    uint numOfModels = static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, lik, &sharedParams);
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelsParams = ChromosomeNumberOptimizer::getMapOfParamsForComplexModel(lik, typeWithParamNames, numOfModels);
    //std::map<int, vector<double>> modelCompositeParams = getVectorToSetModelParams(lik);
    // int baseNumber;
    // (ChromEvolOptions::baseNum_ == IgnoreParam) ? (baseNumber = IgnoreParam) : (baseNumber = static_cast<int>(lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.baseNum_1").getValue()));
    // std::map<uint, uint> maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;

    // std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelCompositeParams, baseNumber, maxBaseNumTransition, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    // std::shared_ptr<SubstitutionModel> model(static_pointer_cast<SubstitutionModel>(chrModel)->clone());

    // std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqsBpp);
    // std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    
    //ParametrizablePhyloTree parTree(*tree_);
    //auto subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone(), shared_ptr<FrequencySet>(rootFrequencies->clone()));
    //subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone());
    //SubstitutionProcess* subProcess = subProSim->clone();
    //auto sequenceData = vsc_->clone();
    //Context context;
    ParametrizablePhyloTree parTree = ParametrizablePhyloTree(*tree_);
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    ChromosomeNumberOptimizer::getMutableMapOfModelAndNodeIds(mapModelNodesIds, lik);
    std::map <uint, uint> baseNumberUpperBound;
    for (size_t m = 1; m <= numOfModels; m ++){
        auto branchProcess = lik->getSubstitutionProcess().getModel(m);
        baseNumberUpperBound[static_cast<uint>(m)] = dynamic_cast<const ChromosomeSubstitutionModel*>(branchProcess)->getMaxChrRange();
    }
    std::shared_ptr<LikelihoodCalculationSingleProcess> likAncestralRec = setHeterogeneousLikInstance(lik, &parTree, baseNumberUpperBound, mapModelNodesIds, modelsParams, true);
    //auto likAncestralRec = std::make_shared<LikelihoodCalculationSingleProcess>(context, *sequenceData, *subProcess, rootFreqs);
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
    auto sequenceData = likAncestralRec->getData();
    auto process = &(likAncestralRec->getSubstitutionProcess());
    auto context = &(likAncestralRec->getContext());
    delete process;
    delete sequenceData;
    delete context;
    
}
/***********************************************************************************/
std::map<int, vector<double>> ChromosomeNumberMng::getVectorToSetModelParams(SingleProcessPhyloLikelihood* lik, size_t modelIndex) const{
    
    ParameterList substitutionParams = lik->getSubstitutionModelParameters();
    std::map<int, vector <double>> compositeParams;

    for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        vector<string> paramNames;
        switch(i){
            case ChromosomeSubstitutionModel::BASENUM:   
                break;
            case ChromosomeSubstitutionModel::BASENUMR:
                paramNames = compositeParameter::getRelatedParameterNames(substitutionParams, "baseNumR");      
                break;
            case ChromosomeSubstitutionModel::DUPL:
                paramNames = compositeParameter::getRelatedParameterNames(substitutionParams, "dupl"); 
                break;
            case ChromosomeSubstitutionModel::LOSS:
                paramNames = compositeParameter::getRelatedParameterNames(substitutionParams, "loss");
                break;
            case ChromosomeSubstitutionModel::GAIN:
                paramNames = compositeParameter::getRelatedParameterNames(substitutionParams, "gain"); 
                break;
            case ChromosomeSubstitutionModel::DEMIDUPL:
                paramNames = compositeParameter::getRelatedParameterNames(substitutionParams, "demi"); 
                break;
            default:
                throw Exception("ChromosomeNumberMng::getVectorToSetModelParams(): Invalid rate type!");
                break;
        }
        if (i == ChromosomeSubstitutionModel::BASENUM){
            continue;
        }
        vector<double> paramValues;
        for (size_t j = 0; j < paramNames.size(); j++){
            paramValues.push_back(lik->getLikelihoodCalculationSingleProcess()->getParameter(paramNames[j]).getValue());
        }
        compositeParams[static_cast<int>(i)] = paramValues;


    }
    return compositeParams; 


}
/***********************************************************************************/
std::shared_ptr<NonHomogeneousSubstitutionProcess> ChromosomeNumberMng::setHeterogeneousModel(ParametrizablePhyloTree* parTree, SingleProcessPhyloLikelihood* ntl, ValueRef <Eigen::RowVectorXd> rootFreqs,  std::map<int, vector<pair<uint, int>>> sharedParams) const{
    uint numOfModels = static_cast<uint>(ntl->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, ntl, &sharedParams);
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelsParams = ChromosomeNumberOptimizer::getMapOfParamsForComplexModel(ntl, typeWithParamNames, numOfModels);
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    ChromosomeNumberOptimizer::getMutableMapOfModelAndNodeIds(mapModelNodesIds, ntl);
    std::map <uint, uint> baseNumberUpperBound;
    for (size_t m = 1; m <= numOfModels; m ++){
        auto branchProcess = ntl->getSubstitutionProcess().getModel(m);
        baseNumberUpperBound[static_cast<uint>(m)] = dynamic_cast<const ChromosomeSubstitutionModel*>(branchProcess)->getMaxChrRange();
    }
    auto rootFreqsValues =  rootFreqs->getTargetValue();
    Vdouble rootFreqsBpp;
    copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelsParams[1].second, modelsParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqsBpp);
    //std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    FrequencySet* rootFrequencies = rootFreqsFixed->clone();
    //std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies->clone());
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);

    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelsParams[i].second, modelsParams[i].first, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
        }   
        subProSim->addModel(std::shared_ptr<ChromosomeSubstitutionModel>(chrModel->clone()), mapModelNodesIds[i]);
    }
    return subProSim;


}
/***********************************************************************************/
std::shared_ptr<LikelihoodCalculationSingleProcess> ChromosomeNumberMng::setHeterogeneousLikInstance(SingleProcessPhyloLikelihood* likProcess, ParametrizablePhyloTree* tree, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, bool forAncestral) const{
    ValueRef <Eigen::RowVectorXd> rootFreqs = likProcess->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    auto rootFreqsValues =  rootFreqs->getTargetValue();
    Vdouble rootFreqsBpp;
    copyEigenToBpp(rootFreqsValues, rootFreqsBpp);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    ParametrizablePhyloTree* parTree = tree->clone();
    
    uint numOfModels = static_cast<uint>(likProcess->getSubstitutionProcess().getNumberOfModels());
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelParams[1].second, modelParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqsBpp);
    //std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    FrequencySet* rootFrequencies = rootFreqsFixed->clone();
    //std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies->clone());
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);

    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelParams[i].second, modelParams[i].first, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
        }   
        subProSim->addModel(std::shared_ptr<ChromosomeSubstitutionModel>(chrModel->clone()), mapModelNodesIds[i]);
    }


    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    std::shared_ptr<LikelihoodCalculationSingleProcess> lik;
    if (forAncestral){
        lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc_->clone(), *nsubPro, rootFreqs);

    }else{
        lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc_->clone(), *nsubPro, false);
    }
    

    //delete subProSim;
    return lik;

}

/***********************************************************************************/
void ChromosomeNumberMng::runChromEvol(){
    setNodeIdsForAllModels(ChromEvolOptions::nodeIdsFilePath_);
    if (ChromEvolOptions::simulateData_){
        //simulate data using a tree and a set of model parameters
        RandomTools::setSeed(static_cast<long>(ChromEvolOptions::seed_));
        simulateData();
        return;

    }
    
    // optimize likelihood
    ChromosomeNumberOptimizer* chrOptimizer = optimizeLikelihoodMultiStartPoints();
    ////////////////////////////////////////////////////////////////
    // !!!!! Note !!!!! The first model should be the root model!!!
    ////////////////////////////////////////////////////////////////
    writeOutputToFile(chrOptimizer);
    // get joint ML ancestral reconstruction
    getJointMLAncestralReconstruction(chrOptimizer);
    //get Marginal ML ancestral reconstruction, and with the help of them- calculate expectations of transitions
    const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+ "ancestorsProbs.txt";
    getMarginalAncestralReconstruction(chrOptimizer, outFilePath);
    // Only temporary ///////////////////////////////////////////////
    //delete chrOptimizer;

    /////////////////////////////////////////////////////////////////
    //compute expectations
    computeExpectations(chrOptimizer, ChromEvolOptions::NumOfSimulations_);
    //The optimizer is deleted inside the computeExpectations object!

    //delete chrOptimizer;


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
void ChromosomeNumberMng::printSimulatedEvoPath(const string outPath, SiteSimulationResult* simResult) const{
    ofstream outFile;
    outFile.open(outPath);
    size_t totalNumTransitions = 0;
    vector<shared_ptr<PhyloNode> > nodes = tree_->getAllNodes();
    size_t nbNodes = nodes.size();
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        if (tree_->getRootIndex() == nodeId){
            outFile << "N-" + std::to_string(nodeId) << endl;
            
            outFile <<"\tThe root state is: "<< ((int)(simResult->getRootAncestralState()+ alphabet_->getMin())) <<endl;
        }else{
            if (tree_->isLeaf(nodeId)){
                outFile << tree_->getNode(nodeId)->getName() << endl;
            }else{
                outFile << "N-" + std::to_string(nodeId) <<endl;

            }
            MutationPath mutPath = simResult->getMutationPath(nodeId);
            vector<size_t> states = mutPath.getStates();
            vector<double> times = mutPath.getTimes();
            totalNumTransitions += static_cast<int>(times.size());

            auto edgeIndex =  tree_->getIncomingEdges(nodeId)[0]; 
            auto fatherIndex = tree_->getFatherOfEdge(edgeIndex);
            outFile << "Father is: " << "N-" << fatherIndex << std::endl;
            size_t fatherState;
            if (fatherIndex == tree_->getRootIndex()){
                fatherState = simResult->getRootAncestralState() + alphabet_->getMin();    
            }else{
                fatherState = simResult->getAncestralState(fatherIndex) + alphabet_->getMin(); 
            }
            for (size_t i = 0; i < states.size(); i++){
                outFile << "from state: "<< fatherState  <<"\tt = "<<times[i] << " to state = "<< ((int)(states[i]) + alphabet_->getMin()) << endl;
                fatherState = ((int)(states[i]) + alphabet_->getMin());
            }
            outFile <<"# Number of transitions per branch: "<< times.size() <<endl;   
            
        }
        
        outFile <<"*************************************"<<endl;
        
    }
    outFile <<"Total number of transitions is: "<< totalNumTransitions << endl;
    outFile.close();

}

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
void ChromosomeNumberMng::convertNodesNames(PhyloTree &tree, uint nodeId, std::map<uint, std::vector<size_t>> &ancestors, bool alphabetStates) const{
    size_t state = ancestors[nodeId][0];
    if (alphabetStates){
        state += alphabet_->getMin();
    }
    
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
            convertNodesNames(tree, sonId, ancestors, alphabetStates);

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
    auto sons = tree.getSons(node) ;
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
  if (!tree.isLeaf(tree.getNode(nodeId))){
      s << tree.getNode(nodeId)->getName();
  }
  shared_ptr<PhyloBranch> branch=tree.getEdgeToFather(nodeId);
  if (branch->hasLength()){
    s << ":" << branch->getLength();

  }

  return s.str();
}
/*********************************************************************************/
void ChromosomeNumberMng::simulateData(){
    if ((ChromEvolOptions::minChrNum_ <= 0) || (ChromEvolOptions::maxChrNum_ < 0)){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): minimum and maximum chromsome number should be positive!");
    }
    if (ChromEvolOptions::maxChrNum_ <= ChromEvolOptions::minChrNum_){
        throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): maximum chromsome number should be larger than minimum chromosome number!");
    }

    alphabet_ = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    ParametrizablePhyloTree* parTree =  new ParametrizablePhyloTree(*tree_);
    std::map<uint, std::pair<int, std::map<int, vector<double>>>> complexParamsValues;
    ChromEvolOptions::getInitialValuesForComplexParams(complexParamsValues);
    std::map<uint, uint> maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    //1. ChromEvolOptions::mapModelNodesIds_: already calculated
    
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, complexParamsValues[1].second, complexParamsValues[1].first, maxBaseNumTransition[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    if (chrModel->getBaseNumber() != IgnoreParam){
        chrModel->correctBaseNumForSimulation(ChromEvolOptions::maxChrInferred_);

    }
    
    if (ChromEvolOptions::fixedFrequenciesFilePath_ == "none"){
        throw Exception("ChromosomeNumberMng::simulateData(): ERROR! The file of fixed root frequencies is missing!!!");

        

    }
    vector <double> rootFreqs = ChromosomeNumberOptimizer::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
    std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies->clone());

    // adding models
    for (uint i = 1; i <= (uint)(ChromEvolOptions::numOfModels_); i++){
        if (i > 1){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, complexParamsValues[i].second, complexParamsValues[i].first, maxBaseNumTransition[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
            if (chrModel->getBaseNumber() != IgnoreParam){
                chrModel->correctBaseNumForSimulation(ChromEvolOptions::maxChrInferred_);

            }
            
        }   
        subProSim->addModel(chrModel, ChromEvolOptions::mapModelNodesIds_[i]);
    }
    SimpleSubstitutionProcessSiteSimulator* simulator = new SimpleSubstitutionProcessSiteSimulator(*subProSim);
    SiteSimulationResult* simResult = simulator->dSimulateSite();
    vector <size_t> leavesStates = simResult->getFinalStates();
    vector<string> leavesNames = simResult->getLeaveNames();
    printSimulatedData(leavesStates, leavesNames, 0);
    printSimulatedDataAndAncestors(simResult);
    if (ChromEvolOptions::resultsPathDir_ != "none"){
        printSimulatedEvoPath(ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedEvolutionPaths.txt", simResult);
    }
    delete simResult;
    delete simulator;
    
    





    
    
    // std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, ChromEvolOptions::gain_, ChromEvolOptions::loss_, 
    //                                                             ChromEvolOptions::dupl_, ChromEvolOptions::demiDupl_, ChromEvolOptions::baseNum_, 
    //                                                             ChromEvolOptions::baseNumR_, ChromEvolOptions::maxBaseNumTransition_, 
    //                                                             ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, 
    //                                                             ChromEvolOptions::rateChangeType_);
    // std::shared_ptr<SubstitutionModel> model(static_pointer_cast<SubstitutionModel>(chrModel)->clone());

    // if (ChromEvolOptions::fixedFrequenciesFilePath_ == "none"){
    //     throw Exception("ERROR!!! ChromosomeNumberMng::simulateData(): You need to supply the path for the file of fixed root frequencies!!");
    // }

    // vector <double> rootFreqs = ChromosomeNumberOptimizer::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
    // std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
    // std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    
    // ParametrizablePhyloTree parTree(*tree_);
    // auto process= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone(), shared_ptr<FrequencySet>(rootFrequencies->clone()));

    // for (size_t i = 0; i < (size_t)ChromEvolOptions::numOfDataToSimulate_; i++){
    //     SimpleSubstitutionProcessSiteSimulator* simulator = new SimpleSubstitutionProcessSiteSimulator(*process);
    //     SiteSimulationResult* simResult = simulator->dSimulateSite();
    //     vector <size_t> leavesStates = simResult->getFinalStates();
    //     vector<string> leavesNames = simResult->getLeaveNames();
    //     printSimulatedData(leavesStates, leavesNames, i);
    //     printSimulatedDataAndAncestors(simResult);
    //     if (ChromEvolOptions::resultsPathDir_ != "none"){
    //         printSimulatedEvoPath(ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedEvolutionPaths.txt", simResult);
    //     }
    //     delete simResult;
    //     delete simulator;

    // }
    // delete rdist;

}
/*******************************************************************************/
void ChromosomeNumberMng::printSimulatedData(vector<size_t> leavesStates, vector<string> leavesNames, size_t iter){
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
        string pathForSimulatedData;
        if (ChromEvolOptions::characterFilePath_ == "none"){
            pathForSimulatedData = ChromEvolOptions::resultsPathDir_ + "//"+ "chr_counts"+ std::to_string(iter) +".fasta";

        }else{
            pathForSimulatedData = ChromEvolOptions::characterFilePath_;
        }
        
        Fasta fasta;
        fasta.writeSequences(pathForSimulatedData, *simulatedData);

    }


    
}
/****************************************************************************/
void ChromosomeNumberMng::printSimulatedDataAndAncestors(SiteSimulationResult* simResult) const{
    std::map<uint, std::vector<size_t> > ancestors;
    vector<shared_ptr<PhyloNode> > nodes = tree_->getAllNodes();
    size_t nbNodes = nodes.size();
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        vector<size_t> nodesStates;
        if (nodeId == tree_->getRootIndex()){
            nodesStates.push_back(simResult->getRootAncestralState()); 
        }else{
            nodesStates.push_back(simResult->getAncestralState(nodeId));
        }       
        ancestors[nodeId] = nodesStates;
    }
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        printTreeWithStates(*tree_, ancestors, ChromEvolOptions::resultsPathDir_);
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_ +"//"+ "simulatedDataAncestors.tree";
        printTreeWithStates(*tree_, ancestors, outFilePath);
    }
  
}

/*****************************************************************************************************/

void ChromosomeNumberMng::computeExpectations(ChromosomeNumberOptimizer* chrOptimizer, int numOfSimulations) const{
    std::cout << "Strating the computation of expectations ...." << std::endl;
    vector<SingleProcessPhyloLikelihood*> vectorOfLikelihoods = chrOptimizer->getVectorOfLikelihoods();
    // get the best likelihood
    SingleProcessPhyloLikelihood* ntl = vectorOfLikelihoods[0];
    auto lik = ntl->getLikelihoodCalculationSingleProcess();
    
    //////////////////////////////////////////////////
    std::map<int, vector<pair<uint, int>>> sharedParams = chrOptimizer->getSharedParams();
    ParametrizablePhyloTree tree =  ParametrizablePhyloTree(*tree_);
    ParametrizablePhyloTree* parTree = (&tree)->clone();
    //ParametrizablePhyloTree parTree = tree;
    ValueRef <Eigen::RowVectorXd> rootFreqs = ntl->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    std::shared_ptr<NonHomogeneousSubstitutionProcess> multiModelProcess =  setHeterogeneousModel(parTree, ntl, rootFreqs, sharedParams);


    size_t nbStates = alphabet_->getSize();
    std::map <uint, std::map<size_t, VVdouble>> jointProbabilitiesFatherSon;
    uint rootId = tree_->getRootIndex();
    vector<shared_ptr<PhyloNode> > nodes = tree_->getAllNodes();
    size_t nbNodes = nodes.size();
    for (size_t n = 0; n < nbNodes; n++){
        uint nodeId = tree_->getNodeIndex(nodes[n]);
        if (nodeId == rootId){
            continue;
        }
        jointProbabilitiesFatherSon[nodeId][0].reserve(nbStates);
        lik->makeJointLikelihoodFatherNode_(nodeId, jointProbabilitiesFatherSon[nodeId][0], 0, 0);
      
    }
    std::cout << "Finished with the calculation of joint likelihoods of father and son..."<< std::endl;
    std::cout << "Starting running simulations ... " << std::endl;
    //creating the model with MLE parameters
    
    // map<int, vector <double>> modelCompositeParams = getVectorToSetModelParams(lik);
    // int baseNumber;
    // (ChromEvolOptions::baseNum_ == IgnoreParam) ? (baseNumber = IgnoreParam) : (baseNumber = static_cast<int>(lik->getLikelihoodCalculationSingleProcess()->getParameter("Chromosome.baseNum_1").getValue()));
    // std::map<uint, uint> maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange_;
    // std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet_, modelCompositeParams, baseNumber, maxBaseNumTransition, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    //the optimizer object in no longer needed
    //auto substitutionProcess = &(lik->getSubstitutionProcess());
    //const NonHomogeneousSubstitutionProcess* substitutionProcessPtr = dynamic_cast<const NonHomogeneousSubstitutionProcess*>(substitutionProcess);


    //const NonHomogeneousSubstitutionProcess* multiModelProcess = &(dynamic_cast<const NonHomogeneousSubstitutionProcess>(lik->getSubstitutionProcess()));
    //std::shared_ptr<NonHomogeneousSubstitutionProcess> multiModelProcess = std::shared_ptr<NonHomogeneousSubstitutionProcess>(substitutionProcessPtr->clone());
    delete chrOptimizer;
    //initializing the expectation instance
    ComputeChromosomeTransitionsExp* expCalculator = new ComputeChromosomeTransitionsExp(multiModelProcess, tree_, alphabet_, jointProbabilitiesFatherSon, ChromEvolOptions::jumpTypeMethod_);
    expCalculator->runSimulations(numOfSimulations);
    std::cout << "Simulations are done. Now strating with the conputation of expectation per type..." << std::endl;
    expCalculator->computeExpectationPerType();
    std::cout << "Computation is done ... Printing results ... " << std::endl;
    if (ChromEvolOptions::resultsPathDir_ == "none"){
        expCalculator->printResults();
    }else{
        const string outFilePath = ChromEvolOptions::resultsPathDir_+"//"+ "expectations.txt";
        //const string outFilePathForNonAccounted = ChromEvolOptions::resultsPathDir_+"//"+ "exp_nonAccounted_branches.txt";
        const string outFilePathHeuristics = ChromEvolOptions::resultsPathDir_+"//"+ "expectations_second_round.txt";
        const string outTreePath = ChromEvolOptions::resultsPathDir_+"//"+ "exp.tree";
        expCalculator->printResults(outFilePath);
        std::cout << "Run heuristics if needed ..." << std::endl;
        expCalculator->runHeuristics();
        std::cout << "Printing heuristics results ... " << std::endl;
        expCalculator->printResults(outFilePathHeuristics);
        std::cout << "Printing the tree of expectations ... " << std::endl;
        PhyloTree* expTree = expCalculator->getResultTree();
        string tree_str = printTree(*expTree);
        std::cout << "Done! Deleting all unneeded objects!" << std::endl;
        delete expTree;
        ofstream treeFile;
        treeFile.open(outTreePath);
        treeFile << tree_str;
        treeFile.close();
    }
    

    delete expCalculator;
}
/******************************************************************************/
void ChromosomeNumberMng::writeOutputToFile(ChromosomeNumberOptimizer* chrOptimizer) const{
    double AICc = chrOptimizer->getAICOfBestModel();
    auto bestLik = chrOptimizer->getVectorOfLikelihoods()[0];
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    ChromosomeNumberOptimizer::getMutableMapOfModelAndNodeIds(mapModelNodesIds, bestLik, tree_->getRootIndex());

    const string outPath = (ChromEvolOptions::resultsPathDir_ == "none") ? (ChromEvolOptions::resultsPathDir_) : (ChromEvolOptions::resultsPathDir_ + "//" + "chromEvol.res");
    ofstream outFile;
    if (outPath != "none"){
        outFile.open(outPath);
    }
    outFile << "Min allowed chromosome number  = " << alphabet_->getMin() << std::endl;
    outFile << "Max allowed chromosome number = " << alphabet_->getMax() << std::endl;
    auto originalTreeLength = getOriginalTreeLength(ChromEvolOptions::treeFilePath_);
    outFile << "Original tree length was: " << originalTreeLength <<std::endl;
    outFile << "Tree scaling factor is: " << tree_->getTotalLength()/originalTreeLength << std::endl;
    outFile << "tree Length was scaled to: " << tree_->getTotalLength() << std::endl;
    auto numOfModels = bestLik->getSubstitutionProcess().getNumberOfModels();
    outFile << "Number of models in the best model = " << numOfModels << std::endl;
    outFile << "Min clade size specified in the parameter file = " << ChromEvolOptions::minCladeSize_ << std::endl;
    // not all the assignements of the nodes induce clades, therefore the min clade size will represent the 
    // number of species under a specific model (better ask Itay)
    // ChromosomeNumberOptimizer::getMutableMapOfModelAndNodeIds(mapModelNodesIds, bestLik);
    uint minSizeOfClade = findMinCladeSize(mapModelNodesIds);
    outFile << "Min clade size in the best model = " << minSizeOfClade << std::endl;
    auto modelAndRepresentitives = findMRCAForEachModelNodes(mapModelNodesIds);
    outFile << "Shifting nodes are: " << std::endl;
    for (uint i = 1; i <= numOfModels; i++){
        for (size_t j = 0; j < modelAndRepresentitives[i].size(); j++){
            outFile << "# Model $" << i << " = " << "N" << modelAndRepresentitives[i][j] << std::endl;
        }
        
    }
    writeTreeWithCorrespondingModels(*tree_, mapModelNodesIds);

    chrOptimizer->printRootFrequencies(bestLik, outFile);
    printLikParameters(chrOptimizer, bestLik, outFile);
    outFile << "AICc of the best model = "<< AICc << std::endl;
    outFile.close();

}
void ChromosomeNumberMng::writeTreeWithCorrespondingModels(PhyloTree tree, std::map<uint, vector<uint>> &modelAndNodes) const{
    std::map<uint, std::vector<size_t>> mapOfNodeAndModel;
    auto it = modelAndNodes.begin();
    while (it != modelAndNodes.end()){
        size_t model = (size_t)(it->first);
        for (size_t i = 0; i < modelAndNodes[it->first].size(); i++){
            uint nodeId = modelAndNodes[it->first][i];
            mapOfNodeAndModel[nodeId].push_back(model);
        }
        it ++;
    }
    uint rootId = tree.getRootIndex();
    convertNodesNames(tree, rootId, mapOfNodeAndModel, false);
    string tree_str = printTree(tree);
    //outFile << tree_str << std::endl;
    string pathForTree = ChromEvolOptions::resultsPathDir_ +"//"+ "treeWithShifts.tree";
    ofstream outFileTree;
    outFileTree.open(pathForTree);
    outFileTree << tree_str << std::endl;
    outFileTree.close();

}
/******************************************************************************/
std::map<uint, std::vector<uint>> ChromosomeNumberMng::findMRCAForEachModelNodes(std::map<uint, vector<uint>> mapOfModelsAndNodes) const{
    std::map <uint, std::vector<uint>> modelWithRepresentitives;
    auto it = mapOfModelsAndNodes.begin();
    while (it != mapOfModelsAndNodes.end()){
        auto nodes = mapOfModelsAndNodes[it->first];
        for (size_t i = 0; i < nodes.size(); i++){
            uint nodeId = nodes[i];
            if (nodeId == tree_->getRootIndex()){
                modelWithRepresentitives[it->first].push_back(nodeId);
                break;
            }
            auto edgeIndex =  tree_->getIncomingEdges(nodeId)[0]; 
            auto fatherIndex = tree_->getFatherOfEdge(edgeIndex);
            if (std::find(nodes.begin(), nodes.end(), fatherIndex) == nodes.end()){
                // if the node has no father in the list, this is the representitive in the current model
                modelWithRepresentitives[it->first].push_back(nodeId);
            }

        }
        it ++;
    }
    return modelWithRepresentitives;

}
/******************************************************************************/
uint ChromosomeNumberMng::findMinCladeSize(std::map<uint, vector<uint>> mapModelNodesIds) const{
    auto it = mapModelNodesIds.begin();
    uint minNumSpecies = (uint)(tree_->getAllLeavesNames().size());
    if (mapModelNodesIds.size() == 1){
        return minNumSpecies;
    }
    while (it != mapModelNodesIds.end()){
        auto nodes = mapModelNodesIds[it->first];
        uint numOfSpecies = 0;
        for (size_t i = 0; i < nodes.size(); i++){
            uint nodeId = nodes[i];
            if (tree_->isLeaf(nodeId)){
                numOfSpecies ++;
            }
        }
        if (minNumSpecies > numOfSpecies){
            minNumSpecies = numOfSpecies;
        }

        it ++;
    }
    return minNumSpecies;
}
/******************************************************************************/
void ChromosomeNumberMng::printLikParameters(ChromosomeNumberOptimizer* chrOptimizer, SingleProcessPhyloLikelihood* lik, ofstream &outFile) const{

    outFile << "Final optimized likelihood is: "<< lik->getValue() << endl;
    outFile << "Final model parameters are:"<<endl;
    ParameterList substitutionModelParams = lik->getSubstitutionModelParameters();
    size_t numOfModels = lik->getSubstitutionProcess().getNumberOfModels();
    std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    std::map<uint, std::map<int, std::vector<std::string>>> mapOfTypeAndName;
    std::map<pair<uint, int>, vector<string>> mapOfAliasedTypeModelAndParam;

    auto sharedParams = chrOptimizer->getSharedParams();
    auto it  = sharedParams.begin();
    while(it != sharedParams.end()){
        // get the string full names of the first parameter to which the other ones are aliased
        auto sharedParametersBlock = sharedParams[it->first];
        uint firstParamModel = sharedParametersBlock[0].first;
        int firstType = sharedParametersBlock[0].second;
        uint numOfSubParams = ChromosomeNumberOptimizer::getNumberOfParametersPerParamType(firstType, ChromEvolOptions::rateChangeType_);      
        string basicName = ChromosomeNumberOptimizer::getStringParamName(firstType);
        std::vector<string> firstParamNames;
        if (firstType == ChromosomeSubstitutionModel::BASENUM){
            firstParamNames.push_back("Chromosome." + basicName + "_"+ std::to_string(firstParamModel));
        }else{
            for (size_t i = 0; i < numOfSubParams; i++){
                firstParamNames.push_back("Chromosome." + basicName +std::to_string(i)+"_"+ std::to_string(firstParamModel));
            }

        }
        // add the aliased parameters to the map, such that they will correspond to the name of the first parameter to which they are aliased
        for (size_t i = 0; i < sharedParametersBlock.size(); i++){
            pair<uint, int> paramModelAndType;
            paramModelAndType.first = sharedParametersBlock[i].first;
            paramModelAndType.second = sharedParametersBlock[i].second;
            string shortName = ChromosomeNumberOptimizer::getStringParamName(sharedParametersBlock[i].second);
            if (paramModelAndType.second == ChromosomeSubstitutionModel::BASENUM){
                mapOfAliasedTypeModelAndParam[paramModelAndType].push_back("Chromosome." + shortName + "_"+ std::to_string(paramModelAndType.first));
                mapOfTypeAndName[sharedParametersBlock[i].first][sharedParametersBlock[i].second].push_back(firstParamNames[0]);
            }else{
                for (size_t j = 0; j < numOfSubParams; j++){
                    mapOfAliasedTypeModelAndParam[paramModelAndType].push_back("Chromosome." + shortName +std::to_string(j)+ "_"+ std::to_string(paramModelAndType.first));
                    mapOfTypeAndName[sharedParametersBlock[i].first][sharedParametersBlock[i].second].push_back(firstParamNames[j]);
                }

            }
            
        }
        // remove the names of the first parameter from the list of names
        for (size_t i = 0; i < numOfSubParams; i++){
            paramsNames.erase(std::remove(paramsNames.begin(), paramsNames.end(), firstParamNames[i]), paramsNames.end());

        }
        
        it ++;
        
     }
    // add the independent parameters
    for (int i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        uint numOfSubParams = ChromosomeNumberOptimizer::getNumberOfParametersPerParamType(i, ChromEvolOptions::rateChangeType_);
        if (numOfSubParams == 0){ // parameter is ignored
            continue;
        }
        string paramBaseName = ChromosomeNumberOptimizer::getStringParamName(i);
        for (uint j = 1; j <= numOfModels; j++){
            for (size_t k = 0; k < numOfSubParams; k++){
                string paramName;
                if (i == ChromosomeSubstitutionModel::BASENUM){
                    paramName = "Chromosome."+  paramBaseName + "_"+ std::to_string(j);
                }else{
                    paramName = "Chromosome."+  paramBaseName + std::to_string(k)+ "_"+ std::to_string(j);
                }
                if (std::find(paramsNames.begin(), paramsNames.end(), paramName) !=  paramsNames.end()){
                    mapOfTypeAndName[j][i].push_back(paramName);
                    pair<uint, int> paramModelAndType;
                    paramModelAndType.first = j;
                    paramModelAndType.second = i;
                    mapOfAliasedTypeModelAndParam[paramModelAndType].push_back(paramName);
                }
            }

        }
        
    }

    // print the values:
    for (uint m = 1; m <= numOfModels; m++){
        auto types = mapOfTypeAndName[m];
        for (int i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
            auto foundType = types.find(i);
            if (foundType == types.end()){
                continue;
            }
            for (size_t k = 0; k < mapOfTypeAndName[m][i].size(); k++){
                std::pair<uint, int> paramModelAndType;
                paramModelAndType.first = m;
                paramModelAndType.second = i;
                if (i != ChromosomeSubstitutionModel::BASENUM){
                
                    outFile << mapOfAliasedTypeModelAndParam[paramModelAndType][k] << " = "<< lik->getLikelihoodCalculation()->getParameter(mapOfTypeAndName[m][i][k]).getValue() << endl;

                }else{
                    outFile << mapOfAliasedTypeModelAndParam[paramModelAndType][k] << " = "<< (int)(lik->getLikelihoodCalculation()->getParameter(mapOfTypeAndName[m][i][k]).getValue()) << endl;
                }
            }        
        }
    }
}
/******************************************************************************/
double ChromosomeNumberMng::getOriginalTreeLength(string &path) const{
    Newick reader;
    auto originalTree = reader.readPTree(path);
    double treeLength = originalTree->getTotalLength();
    return treeLength;
}
