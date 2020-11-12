#include "ComputeChromosomeTransitionsExp.h"
using namespace bpp;

/****************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationOfChangePerBranch(int nodeId, VVdouble &jointProbFatherNode, int jumpType){

    double expectation = 0;
    for (size_t x = 0; x < alphabet_->getSize(); x ++){
        for (size_t y = 0; y < alphabet_->getSize(); y++){
            expectation += jointProbFatherNode[x][y] * getExpectation(nodeId, (int)y, (int)x, jumpType);

        }
    }
    expNumOfChangesPerBranch_[nodeId][jumpType] += expectation;
    expNumOfChanges_[jumpType] += expectation;
    return;

}
/**********************************************************************************************/
ChromosomeSubstitutionModel::typeOfTransition ComputeChromosomeTransitionsExp::getTypeOfTransition(int startState, int endState){
    // convert from state index to real chromsome number
    int chrStart = startState + alphabet_->getMin();
    int chrEnd = endState + alphabet_->getMin();
    //baseNumber transitions
    int baseNumber = model_->getBaseNumber();
    if (baseNumber != IgnoreParam){
        if (chrEnd > chrStart){
            if ((chrEnd - chrStart) % baseNumber == 0){
                if (chrEnd != 2 * chrStart){
                    return ChromosomeSubstitutionModel::BASENUM_T;
                }
            }
        }        
    }
    //gain
    if (chrStart + 1 == chrEnd){
        if ((model_->getConstGain() != IgnoreParam) || (model_->getChangeRateGain() != IgnoreParam)){
            return ChromosomeSubstitutionModel::GAIN_T;
        }
        
    //loss
    }else if (chrStart == chrEnd + 1){
        if ((model_->getConstLoss() != IgnoreParam) || (model_->getChangeRateLoss() != IgnoreParam)){
            return ChromosomeSubstitutionModel::LOSS_T;
        }
        
    //dupl
    }else if (chrStart * 2 == chrEnd){
        if ((model_->getConstDupl() != IgnoreParam) || (model_->getChangeRateDupl() != IgnoreParam)){
            return ChromosomeSubstitutionModel::DUPL_T;
        }
        
    //demi dupl
    }else if (model_->getDemiDupl() != IgnoreParam){
        if ((chrStart % 2 == 0) && (chrEnd == chrStart * 1.5)){
            return ChromosomeSubstitutionModel::DEMIDUPL_T;
        }else if ((chrStart % 2 != 0) && (chrEnd == (int)ceil(chrStart * 1.5))){
            return ChromosomeSubstitutionModel::DEMIDUPL_T;
        }else if ((chrStart % 2 != 0) && (chrEnd == (int)floor(chrStart * 1.5))){
            return ChromosomeSubstitutionModel::DEMIDUPL_T;
        }

    }

    else if (chrEnd == alphabet_->getMax()){
        return ChromosomeSubstitutionModel::MAXCHR_T;
    }
    return ChromosomeSubstitutionModel::ILLEGAL;

}

/**************************************************************************************/
ChromosomeSubstitutionModel::typeOfTransition ComputeChromosomeTransitionsExp::getTypeOfTransitionWithProb(int startState, int endState){
    std::vector<ChromosomeSubstitutionModel::typeOfTransition> jumpType;
    bool legalMove = false;
    // convert from state index to real chromsome number
    int chrStart = startState + alphabet_->getMin();
    int chrEnd = endState + alphabet_->getMin();
    pair <int, int> jumpStates;
    jumpStates.first = chrStart;
    jumpStates.second = chrEnd;

    if (stateJumpTypeProb_.find(jumpStates) != stateJumpTypeProb_.end()){

    }
    //gain
    if (chrStart + 1 == chrEnd){
        if ((model_->getConstGain() != IgnoreParam) || (model_->getChangeRateGain() != IgnoreParam)){
            return ChromosomeSubstitutionModel::GAIN_T;
        }
    }
    //loss
    if (chrStart - 1 == chrEnd){
        if ((model_->getConstLoss() != IgnoreParam) || (model_->getChangeRateLoss() != IgnoreParam)){
            return ChromosomeSubstitutionModel::LOSS_T;
        }
    }
    //baseNumber transitions
    int baseNumber = model_->getBaseNumber();
    if (baseNumber != IgnoreParam){
        if (chrEnd > chrStart){
            if ((chrEnd - chrStart) % baseNumber == 0){
                legalMove = true;
                jumpType.push_back(ChromosomeSubstitutionModel::BASENUM_T);                
            }
        }        
    }
    // duplication
    if (chrEnd == 2 * chrStart){
        if ((model_->getConstDupl() != IgnoreParam) || (model_->getChangeRateDupl() != IgnoreParam)){
            legalMove = true;
            jumpType.push_back(ChromosomeSubstitutionModel::DUPL_T);
        }

    }
    //Demi-duplication
    if (model_->getDemiDupl() != IgnoreParam){
        if (chrStart % 2 == 0){
            if (chrEnd == chrStart * 1.5){
                legalMove = true;
                jumpType.push_back(ChromosomeSubstitutionModel::DEMIDUPL_T);
            }
        }else{
            if ((chrEnd == (int)ceil(chrStart * 1.5)) || (chrEnd == (int)floor(chrStart * 1.5))){
                legalMove = true;
                jumpType.push_back(ChromosomeSubstitutionModel::DEMIDUPL_T);
            }
        }
            
    }

    // maxChr not assigned to any of the possible transitions
    if ((chrEnd  == alphabet_->getMax()) && (!legalMove)){
        legalMove = true;
        return ChromosomeSubstitutionModel::MAXCHR_T;
    }if(!legalMove){
        return ChromosomeSubstitutionModel::ILLEGAL;
    }
    // choose the most probable transition
     if (jumpType.size() == 1)
        return jumpType[0];
    
    std::vector<double> weights;
    std::vector<size_t> indices;
    std::vector<size_t> result;
    result.push_back(0);
    double sumOfRates = 0;
   
    
    //sample the jump type according randomely according to probabilities
    
    for (size_t i = 0; i < jumpType.size(); i++){
        if (jumpType[i] == ChromosomeSubstitutionModel::BASENUM_T){
            sumOfRates += model_->getBaseNumR();
            weights.push_back(model_->getBaseNumR());
            indices.push_back(i);
        }else if(jumpType[i] == ChromosomeSubstitutionModel::DUPL_T){
            double rate =  model_->getRate(chrStart, model_->getConstDupl(), model_->getChangeRateDupl());
            sumOfRates += rate;
            weights.push_back(rate);
            indices.push_back(i);
        }else if (jumpType[i] == ChromosomeSubstitutionModel::DEMIDUPL_T){
            sumOfRates += model_->getDemiDupl();
            weights.push_back(model_->getDemiDupl());
            indices.push_back(i);

        }
        // add also gain 1->2

    }
    for (size_t j = 0; j < weights.size(); j++){
        weights[j] /= sumOfRates;
    }
    
    RandomTools::getSample(indices, weights, result);
    
        
    
    return jumpType[result[0]];    
}
/**********************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpPerTypeHeuristics(map <int, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun){
    map <int, vector<pair<int, int>>>::iterator it = nonAccountedForBranchesFromFirstRun.begin();
    while (it != nonAccountedForBranchesFromFirstRun.end()){
        int nodeId = it->first;
        vector <pair<int,int>> terminals = nonAccountedForBranchesFromFirstRun[nodeId];
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
            for (size_t j = 0; j < terminals.size(); j++){
                pair <int, int> terminal = terminals[j];
                size_t father = terminal.first;
                size_t son = terminal.second;          
                double expectation = jointProbabilitiesFatherSon_[nodeId][0][son][father] * getExpectation(nodeId, (int)father, (int)son, i);
                expNumOfChangesPerBranch_[nodeId][i] += expectation;
                expNumOfChanges_[i] += expectation;


            }

        }
        it++;
    }

}

/**********************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationPerType(){
    //init();
    vector<int> nodeIds = tree_->getNodesId();
    for (size_t n = 0; n < nodeIds.size(); n++){
        if (tree_->getRootId() == nodeIds[n]){
            continue;
        }
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
            computeExpectationOfChangePerBranch(nodeIds[n], jointProbabilitiesFatherSon_[nodeIds[n]][0], i);

        }
    }

}

/*********************************************************************************/
void ComputeChromosomeTransitionsExp::printResults(const string path) {
    ofstream outFile;
    if (path == "none"){
        throw Exception("ERROR!!! ComputeChromosomeTransitionsExp::printResults(): not provided file path!\n");
    }
    outFile.open(path);
    vector<string> typeNames;
    std::vector<int> nodesIds = tree_->getNodesId();
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
        vector<string> statementes;
        string typeOfRate;
        if (i == ChromosomeSubstitutionModel::GAIN_T){
            typeOfRate = "GAIN";
            statementes.push_back("#Nodes with GAIN events with expectation above ");
            statementes.push_back("Total number of gain events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::LOSS_T){
            typeOfRate = "LOSS";
            statementes.push_back("#Nodes with LOSS events with expectation above ");
            statementes.push_back("Total number of loss events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::DUPL_T){
            typeOfRate = "DUPLICATION";
            statementes.push_back("#Nodes with duplication events with expectation above ");
            statementes.push_back("Total number of duplication events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::DEMIDUPL_T){
            typeOfRate = "DEMI-DUPLICATION";
            statementes.push_back("#Nodes with demi-duplication events with expectation above ");
            statementes.push_back("Total number of demi-duplications events: ");
            typeNames.push_back(typeOfRate);
        }else if (i == ChromosomeSubstitutionModel::BASENUM_T){
            typeOfRate = "BASE-NUMBER";
            statementes.push_back("#Nodes with transitions in base number events with expectation above ");
            statementes.push_back("Total number of transitions in base number: ");
            typeNames.push_back(typeOfRate);
        }else{
            typeOfRate = "TOMAX";
            typeNames.push_back(typeOfRate);
            continue;
        }
        outFile << statementes[0] <<THRESHOLD_EXP<<endl;
        
        for (size_t n = 0; n < nodesIds.size(); n++){
            if (expNumOfChangesPerBranch_[nodesIds[n]][i] > THRESHOLD_EXP){
                string nodeName;
                if (tree_->getRootId() == nodesIds[n]){
                    continue;
                }
                if (tree_->isLeaf(nodesIds[n])){
                    nodeName = tree_->getNodeName(nodesIds[n]);

                }else{
                    nodeName = "N" + std::to_string(nodesIds[n]);
                }
                outFile << nodeName <<": "<<expNumOfChangesPerBranch_[nodesIds[n]][i]<<endl; 
            }
        }
        outFile << statementes[1] <<  expNumOfChanges_[i] <<endl;
        outFile <<"#+++++++++++++++++++++++++++++\n\n";
    }
    outFile << "EVENTS NOT ACCOUNTED FOR IN THE SIMULATIONS: "<< endl;
    for (size_t n = 0; n < nodesIds.size(); n++){
        if (nodesIds[n] == tree_->getRootId()){
            continue;
        }
        double cumulativeProb = getCumulativeProbability(nodesIds[n], 0);
        if (cumulativeProb < THRESHOLD_HEURISTIC){
            string nodeName;
            if (tree_->isLeaf(nodesIds[n])){
                nodeName = tree_->getNodeName(nodesIds[n]);

            }else{
                nodeName = "N" + std::to_string(nodesIds[n]);
            }
            outFile << nodeName <<": "<< cumulativeProb <<endl;

        }

    }
    outFile <<"#+++++++++++++++++++++++++++++\n\n";
    outFile << "#ALL EVENTS EXPECTATIONS PER NODE"<<endl;
    outFile <<"NODE\t";
    for (size_t i = 0; i < typeNames.size(); i++){
        if (i == typeNames.size()-1){
           outFile <<typeNames[i] << endl;
           continue; 
        }
        outFile <<typeNames[i] <<"\t";
    }
    for (size_t n = 0; n < nodesIds.size(); n++){
        string nodeName;
        if (tree_->getRootId() == nodesIds[n]){
            continue;
        }
        if (tree_->isLeaf(nodesIds[n])){
            nodeName = tree_->getNodeName(nodesIds[n]);
        }else{
            nodeName = "N" + std::to_string(nodesIds[n]);
        }
        outFile << nodeName <<"\t";
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
            (i == ChromosomeSubstitutionModel::NUMTYPES-1) ? (outFile << expNumOfChangesPerBranch_[nodesIds[n]][i] <<endl) : (outFile <<expNumOfChangesPerBranch_[nodesIds[n]][i] <<"\t");
        }
    }
    outFile << "#+++++++++++++++++++++++++++++\n\n"<<endl;
    //get expected number of changes from root to tip
    outFile << "#EXPECTED NUMBER OF EVENTS FROM ROOT TO LEAF"<<endl;
    outFile <<"NODE\t";
    for (size_t i = 0; i < typeNames.size(); i++){
        if (i == typeNames.size()-1){
           outFile <<typeNames[i] << endl;
           continue; 
        }
        outFile <<typeNames[i] <<"\t";
    }
    vector<int> leavesIds = tree_->getLeavesId();
    for (size_t j = 0; j < leavesIds.size(); j++){
        string leafName = tree_->getNodeName(leavesIds[j]);
        outFile << leafName <<"\t";
        for (int k = 0; k < ChromosomeSubstitutionModel::NUMTYPES; k++){
            double expectedRootToTip = 0;
            int currNodeId = leavesIds[j];
            while (tree_->getRootId() != currNodeId){
                expectedRootToTip += expNumOfChangesPerBranch_[currNodeId][k];
                currNodeId = tree_->getFatherId(currNodeId);
            }
            if (k == ChromosomeSubstitutionModel::NUMTYPES-1){
                outFile << expectedRootToTip <<endl;

            }else{
                outFile << expectedRootToTip <<"\t";
            }
                
        }

    }

    outFile <<"#+++++++++++++++++++++++++++++\n\n";
    outFile <<"#TOTAL EXPECTATIONS:\n";
    // print total expectations
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
        outFile << typeNames[i] <<": "<< expNumOfChanges_[i] <<endl;

    }
    outFile.close();
    

}
/*****************************************************************************************/
TreeTemplate<Node>* ComputeChromosomeTransitionsExp::getResultTree(){

    TreeTemplate<Node>* printTree = tree_->clone();
    //string branchProp = "expectation";
    vector <int> nodeIds = printTree->getNodesId();
    for (size_t n = 0; n < nodeIds.size(); n++){
        string nodeName;
        if (printTree->getRootId() == nodeIds[n]){
            nodeName = "N" + std::to_string(nodeIds[n]);
            printTree->setNodeName(nodeIds[n], nodeName);
            continue;
        }else if (printTree->isLeaf(nodeIds[n])){
            nodeName = printTree->getNodeName(nodeIds[n]);
            
        }else{
            nodeName = "N" + std::to_string(nodeIds[n]);

        }
        string expected = "[";
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
            if (i == ChromosomeSubstitutionModel::NUMTYPES - 1){
                expected = expected + to_string(expNumOfChangesPerBranch_[nodeIds[n]][i]) + "]";
            }else{
                expected = expected + to_string(expNumOfChangesPerBranch_[nodeIds[n]][i])+ "\\";
            }            

        }
        nodeName = nodeName + expected;
        printTree->setNodeName(nodeIds[n], nodeName);

        
        
    }

    return printTree;

}
//****************************************************************************************/
void ComputeChromosomeTransitionsExp::init(){
    jumpProbs_.reserve(alphabet_->getSize());
    vector <int> nodesIds = tree_->getNodesId();
    
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        vector <double> stateTransitionsProbs;
        stateTransitionsProbs.reserve(alphabet_->getSize());
        waitingTimes_.push_back(-1 * (model_->Qij(i, i)));
        for (size_t j = 0; j < alphabet_->getSize(); j++){
            if (i == j){
                stateTransitionsProbs.push_back(0);
            }else{
                stateTransitionsProbs.push_back(model_->Qij(i,j)/waitingTimes_[i]);
            }            
        }
        jumpProbs_.push_back(stateTransitionsProbs);
    }
    for (size_t n = 0; n < nodesIds.size(); n++){
        if (nodesIds[n] == tree_->getRootId()){
            continue;
        }
        const Node node = *(tree_->getNode(nodesIds[n]));
        branchOrder_.push_back(node);
        for (size_t i = 0; i < alphabet_->getSize(); i++){
            for (size_t j = 0; j < alphabet_->getSize(); j ++){
                std::pair<int,int> ancestralTerminals;
                ancestralTerminals.first = (int)i;
                ancestralTerminals.second = (int)j;
                ancestralTerminalsCounts_[nodesIds[n]][ancestralTerminals] = 0;
                for (size_t k = 0; k < ChromosomeSubstitutionModel::NUMTYPES; k++){
                    branchTransitionsExp_[nodesIds[n]][ancestralTerminals].push_back(0);

                }

            }
        }
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
            expNumOfChangesPerBranch_[nodesIds[n]][i] = 0;
        }
    }
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
        expNumOfChanges_[i] = 0;
    }
    sort(branchOrder_.begin(), branchOrder_.end(), compareBranches);
}

/************************************************************************************/
bool ComputeChromosomeTransitionsExp::compareBranches(Node& node1, Node& node2){
    return (node1.getDistanceToFather() < node2.getDistanceToFather());
}
/************************************************************************************/

void ComputeChromosomeTransitionsExp::runSimulations(int numOfSimulations){
    init();
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        for (size_t j = 0; j < (size_t)numOfSimulations; j++){
            runIteration((int)i);
        }
    }
    computeExpectationAndPosterior();
}
/*************************************************************************************/
void ComputeChromosomeTransitionsExp::runIteration(int beginState, map <int, vector<pair<int,int>>>* unAccountedNodesAndTerminals){
    double totalTimeTillJump = 0;
    double maxBranch = branchOrder_[branchOrder_.size()-1].getDistanceToFather();
    int currentState = beginState;
    vector <std::pair<int,int>> jumpsUntilNow;
    int indexOfSmallestNotUpdatedBranch = 0;
    while (totalTimeTillJump < maxBranch){
        double averageWaitingTime = 1 / waitingTimes_[currentState];
        totalTimeTillJump +=  RandomTools::randExponential(averageWaitingTime);
        for (size_t i = indexOfSmallestNotUpdatedBranch; i < branchOrder_.size(); i++){
            double branchLength = branchOrder_[i].getDistanceToFather();
            if (branchLength > totalTimeTillJump){
                indexOfSmallestNotUpdatedBranch = (int)i;
                break;
            }
            std::pair<int, int> ancestralTerminals;
            ancestralTerminals.first = beginState;
            ancestralTerminals.second = currentState;
            if (unAccountedNodesAndTerminals){
                if (find((*unAccountedNodesAndTerminals)[branchOrder_[i].getId()].begin(), (*unAccountedNodesAndTerminals)[branchOrder_[i].getId()].end(), ancestralTerminals) == (*unAccountedNodesAndTerminals)[branchOrder_[i].getId()].end()){
                    continue;
                }

            }


            ancestralTerminalsCounts_[branchOrder_[i].getId()][ancestralTerminals] += 1;
            for (size_t j = 0; j < jumpsUntilNow.size(); j++){
                updateExpectationsPerBranch(branchOrder_[i].getId(), ancestralTerminals, jumpsUntilNow[j]);

            }
        }

        int nextState = getRandomState(currentState);
        if ((nextState == (int)alphabet_->getSize()-1) && (!isMaxStateValid(currentState))){
            break;
        }
        std::pair<int,int> combOcurredStates;
        combOcurredStates.first = currentState;
        combOcurredStates.second = nextState;
        jumpsUntilNow.push_back(combOcurredStates);
        updateMapOfJumps(currentState, nextState);
        currentState = nextState;
    }

}
/********************************************************************************/
bool ComputeChromosomeTransitionsExp::isMaxStateValid(int prevState) const{
    int valid = false;
    int maxState = model_->getMax();
    int initState = model_->getMin() + prevState;
    // gain
    if (maxState == initState + 1){
        valid = true;
        return valid;
    }
    // dupl
    if ((model_->getConstDupl() != IgnoreParam) || (model_->getChangeRateDupl() != IgnoreParam)){
        if (maxState == 2 * initState){
            valid = true;
            return valid;
        }

    }
    // demi dupl
    if (model_->getDemiDupl() != IgnoreParam){
        if (initState % 2 == 0){
            if ((int)(initState * 1.5) == maxState){
                valid = true;
            }
        }else{
            if ((maxState == (int)ceil(initState * 1.5)) || (maxState == (int)floor(initState * 1.5))){
                valid = true;
            }
        }
        if (valid){
            return valid;
        }
    }
    // base number
    int baseNumber = model_->getBaseNumber();
    if (baseNumber != IgnoreParam){
        if (maxState > initState){
            if (((maxState - initState) % baseNumber == 0) && ((maxState - initState) <= (int)(model_->getMaxChrRange()))){
                valid = true;
                return valid;
            }
        }        
    }
    return valid;  
    
}
/********************************************************************************/
void ComputeChromosomeTransitionsExp::updateMapOfJumps(int startState, int endState){
    bool legalMove = false;
    pair <int, int> jumpStates;
    jumpStates.first = startState;
    jumpStates.second = endState;

    // convert from state index to real chromsome number
    int chrStart = startState + alphabet_->getMin();
    int chrEnd = endState + alphabet_->getMin();
    double sumOfRates = 0; //for normalization of weights

    if (stateJumpTypeProb_.find(jumpStates) == stateJumpTypeProb_.end()){
        //gain
        if (chrStart + 1 == chrEnd){
            if ((model_->getConstGain() != IgnoreParam) || (model_->getChangeRateGain() != IgnoreParam)){
                if (chrStart > 3){
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::GAIN_T] = 1;
                    return;
                }else{
                    double gainRate = model_->getRate(chrStart, model_->getConstGain(), model_->getChangeRateGain());
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::GAIN_T] = gainRate;
                    sumOfRates += gainRate;
                    legalMove = true;
                }

            }
        }
        //loss
        if (chrStart - 1 == chrEnd){
            if ((model_->getConstLoss() != IgnoreParam) || (model_->getChangeRateLoss() != IgnoreParam)){
                stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::LOSS_T] = 1;
                return;
            }
        }
        //baseNumber transitions
        int baseNumber = model_->getBaseNumber();
        if (baseNumber != IgnoreParam){
            if (chrEnd > chrStart){
                if (((chrEnd - chrStart) % baseNumber == 0) && ((chrEnd - chrStart) <= (int)(model_->getMaxChrRange()))){
                    legalMove = true;
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::BASENUM_T] = model_->getBaseNumR();
                    sumOfRates += model_->getBaseNumR();                                 
                }
            }        
        }
        //duplication
        if (chrEnd == 2 * chrStart){
            if ((model_->getConstDupl() != IgnoreParam) || (model_->getChangeRateDupl() != IgnoreParam)){
                legalMove = true;
                double duplRate = model_->getRate(chrStart, model_->getConstDupl(), model_->getChangeRateDupl());
                stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::DUPL_T] = duplRate;
                sumOfRates += duplRate;
            }

        }
        //demi-duplication
        if (model_->getDemiDupl() != IgnoreParam){
            if (chrStart % 2 == 0){
                if (chrEnd == chrStart * 1.5){
                    legalMove = true;
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::DEMIDUPL_T] = model_->getDemiDupl();
                    sumOfRates += model_->getDemiDupl();
                }
            }else{
                if ((chrEnd == (int)ceil(chrStart * 1.5)) || (chrEnd == (int)floor(chrStart * 1.5))){
                    legalMove = true;
                    double demiDupRate;
                    if (chrStart == 1){
                        demiDupRate =  model_->getDemiDupl();
                    }else{
                        demiDupRate = model_->getDemiDupl()/2;
                    }
                    stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::DEMIDUPL_T] = demiDupRate;
                    sumOfRates += demiDupRate;
                }
            }

        }
        //maxChr not assigned to any of the possible transitions
        // if ((chrEnd  == alphabet_->getMax()) && (!legalMove)){
        //     legalMove = true;
        //     stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::MAXCHR_T] = 1;
        //     return;
        
        // }
        if (chrEnd  == alphabet_->getMax()){
            legalMove = true;
            //stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::MAXCHR_T] = 1;
            double toMaxRate = model_->Qij(startState, endState)-sumOfRates;
            stateJumpTypeProb_[jumpStates][ChromosomeSubstitutionModel::MAXCHR_T] = toMaxRate;
            sumOfRates += toMaxRate;
            //return;
        }
        //if nothing fits
        if(!legalMove){
            throw Exception ("ERROR: ComputeChromosomeTransitionsExp::updateMapOfJumps(): Illegal transition!");
            return;
        }
        //DEBUG!!!
        // if (sumOfRates != model_->Qij(startState, endState)){
        //     cout <<"sum of rates is: "<< sumOfRates <<endl;
        //     cout << "entry in matrix is: " << model_->Qij(startState, endState)<< endl;
        //     throw Exception ("ERROR: ComputeChromosomeTransitionsExp::updateMapOfJumps(): sumOfRates does not equal its supposed value!");
        // }
        // normalize according to weights
        map <int, double>::iterator it = stateJumpTypeProb_[jumpStates].begin();
        while (it != stateJumpTypeProb_[jumpStates].end()){
            int typeOfJump = it->first;
            double rate = stateJumpTypeProb_[jumpStates][typeOfJump];
            stateJumpTypeProb_[jumpStates][typeOfJump] = rate / sumOfRates;
            it ++;
        }
            
    }


}
/********************************************************************************/
void ComputeChromosomeTransitionsExp::updateExpectationsPerBranch(int nodeId, pair<int, int> ancestralTerminals, pair<int, int> jumpStates){
    map <int, double>::iterator it = stateJumpTypeProb_[jumpStates].begin();
    while (it != stateJumpTypeProb_[jumpStates].end()){
        int typeOfJump = it->first;
        double prob = stateJumpTypeProb_[jumpStates][typeOfJump];
        branchTransitionsExp_[nodeId][ancestralTerminals][typeOfJump] += prob;
        it ++;
    }
    
}

/********************************************************************************/
int ComputeChromosomeTransitionsExp::getRandomState(int currentState){
    double prob = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.0);
    double cumulativeProb = 0;
    int nextState = currentState;
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        cumulativeProb += jumpProbs_[currentState][(int)i];
        if (prob < cumulativeProb){
            nextState = (int)i;
            return nextState;
        }
    }
    throw Exception("ERROR: ComputeChangesExpectations::getRandomState(): could not sample new state.");
    return 1;
}
/*************************************************************************************/
double ComputeChromosomeTransitionsExp::getCumulativeProbability(int nodeId, vector <pair<int, int>>* terminalsToAccount){
    double cumulativeProb  = 0;
    std::map<std::pair<int, int>, int> :: iterator iterTerminalStates = ancestralTerminalsCounts_[nodeId].begin();
    while(iterTerminalStates != ancestralTerminalsCounts_[nodeId].end()){
        std::pair <int, int> currentPairOfAncestralTerminals = iterTerminalStates->first;
        int countForPairOfTerminals = ancestralTerminalsCounts_[nodeId][currentPairOfAncestralTerminals];
        if (countForPairOfTerminals == 0){
            iterTerminalStates ++;
            if (terminalsToAccount){
                terminalsToAccount->push_back(currentPairOfAncestralTerminals);
                
            }
            continue;

        }
        cumulativeProb += jointProbabilitiesFatherSon_[nodeId][0][currentPairOfAncestralTerminals.second][currentPairOfAncestralTerminals.first];
        iterTerminalStates ++;
    }
    return cumulativeProb;

}
/*************************************************************************************/
double ComputeChromosomeTransitionsExp::isNeededHeuristics(int nodeId, map <int, vector<pair<int,int>>>* unAccountedNodesAndTerminals){
    vector <pair<int, int>> terminalsToAccount;
    double cumulativeProb = getCumulativeProbability(nodeId, &terminalsToAccount);
    if (cumulativeProb < THRESHOLD_HEURISTIC){
        (*unAccountedNodesAndTerminals)[nodeId] = terminalsToAccount;       
        

    }else{
        if (unAccountedNodesAndTerminals->find(nodeId) != unAccountedNodesAndTerminals->end()){
            unAccountedNodesAndTerminals->erase(nodeId);
        }

    }    
    return cumulativeProb;

}
/*************************************************************************************/
void ComputeChromosomeTransitionsExp::updateNumNonAccountedBranches(map <int, vector<pair<int,int>>>* unAccountedNodesAndTerminals, int iteration, const string filePath){
    vector <double> cumulativeProbs;
    vector<Node> branches;
    for (size_t n = 0; n < branchOrder_.size(); n++){
        branches.push_back(branchOrder_[n]);
    }
   
    for (size_t i = 0; i < branches.size(); i++){
        int nodeId = branches[i].getId();
        double thresholdHeuristics = isNeededHeuristics(nodeId, unAccountedNodesAndTerminals);
        if (thresholdHeuristics >= THRESHOLD_HEURISTIC){
            Node node = branches[i];
            branchOrder_.erase(remove(branchOrder_.begin(), branchOrder_.end(), node), branchOrder_.end());
            //if (iteration > 0){
                //branchMultiplier->erase(nodeId);
            //}

        }else{
            if (iteration == 0){
                cumulativeProbs.push_back(thresholdHeuristics);
            }

        }
        

    }
    if (iteration == 0){
        
        if (filePath == "none"){
            return;
        }
        ofstream outFile;
        outFile.open(filePath);
        outFile << "Non Accounted for transitions:"<<endl;
        for (size_t k = 0; k < cumulativeProbs.size(); k++){
            string nodeName;
            int nodeId = branchOrder_[k].getId();
            if (tree_->isLeaf(nodeId)){
                nodeName = tree_->getNodeName(nodeId);

            }else{
                nodeName = "N" + std::to_string(nodeId);
            }
            outFile << nodeName <<": "<<cumulativeProbs[k]<<endl;

        }
        outFile.close();


    }


}
/*************************************************************************************/
void ComputeChromosomeTransitionsExp::runHeuristics(const string FilePath){
    map <int, vector<pair<int, int>>> nonAccountedForBranchesFromFirstRun;
    //map <int, double> branchMultiplier;
    map <int, vector<pair<int,int>>> unAccountedNodesAndTerminals;
    map <int, double> ratesPerState;
    for (int i = 0; i < MAX_ITER_HEURISTICS; i++){
        updateNumNonAccountedBranches(&unAccountedNodesAndTerminals, i, FilePath);
        if (i == 0){
            nonAccountedForBranchesFromFirstRun = unAccountedNodesAndTerminals;
        }
        if (branchOrder_.size() == 0){
            break;
        }
        vector<int> initStates = setVectorOfInitStatesForHeuristics(unAccountedNodesAndTerminals);
      
        for (int k = 0; k < (int)initStates.size(); k++){
            updateBranchLengths(initStates[k], i, &ratesPerState);
            for (int j = 0; j < MAX_SIM_HEURISTICS; j++){
                runIteration(k, &unAccountedNodesAndTerminals);
            }
        }

    }
    getPosteriorAndExpForNonAccountedFor(nonAccountedForBranchesFromFirstRun); 
    computeExpPerTypeHeuristics(nonAccountedForBranchesFromFirstRun);
}
/*************************************************************************************/
void ComputeChromosomeTransitionsExp::updateBranchLengths(int initState, int iteration, map <int, double>* ratesPerState){
    double sumOfRates = 0;
    if (iteration == 0){
        vector<double> rates;
        //dupl
        rates.push_back(model_->getRate(initState, model_->getConstDupl(), model_->getChangeRateDupl()));
        //demi-dupl
        rates.push_back(model_->getDemiDupl());
        //base number
        double baseNumRate = 0;
        int baseNumber = model_->getBaseNumber();
        // we would like to take into consideration the overall rate for a base number transition from the current state
        if (baseNumber != IgnoreParam){
            for (int i = initState + baseNumber; i < model_->getMax()-model_->getMin(); i+=baseNumber){
                if ((i - initState) > (int)(model_->getMaxChrRange())){
                    break;
                }
                baseNumRate += model_->getBaseNumR();

            }
        }
        (baseNumber == IgnoreParam) ? (rates.push_back(model_->getBaseNumR())) : (rates.push_back(baseNumRate));
        //gain
        rates.push_back(model_->getRate(initState, model_->getConstGain(), model_->getChangeRateGain()));
        //loss
        rates.push_back(model_->getRate(initState, model_->getConstLoss(), model_->getChangeRateLoss()));
        for (size_t i = 0; i < rates.size(); i++){
            if (rates[i] == IgnoreParam){
                continue;
            }
            sumOfRates += rates[i];

        }
        (*ratesPerState)[initState] = sumOfRates;

    }


    //multiplying by factor
    for (size_t k = 0; k < branchOrder_.size(); k++){
        if (iteration == 0){
            if ((*ratesPerState)[initState] * branchOrder_[k].getDistanceToFather() < 1){
                branchOrder_[k].setDistanceToFather(1/(*ratesPerState)[initState]);
            }
        }else{
            if (tree_->getDistanceToFather(branchOrder_[k].getId()) * (*ratesPerState)[initState] >= 1){
                branchOrder_[k].setDistanceToFather(tree_->getDistanceToFather(branchOrder_[k].getId()) * BRANCH_MULTIPLIER_FACTOR * iteration);

            }else{
                branchOrder_[k].setDistanceToFather((1/(*ratesPerState)[initState]) * (BRANCH_MULTIPLIER_FACTOR * iteration));
            }
            
        }
         

    }
    sort(branchOrder_.begin(), branchOrder_.end(), compareBranches);

}
/*************************************************************************************/
vector <int> ComputeChromosomeTransitionsExp::setVectorOfInitStatesForHeuristics(map <int, vector<pair<int,int>>>& unAccountedNodesAndTerminals) const{
    vector <int> initStates;
    std::map <int, vector<pair<int, int>>>::iterator it = unAccountedNodesAndTerminals.begin();
    while (it != unAccountedNodesAndTerminals.end()){
        int nodeId = it->first;
        vector <pair<int, int>> terminalsPerNode = unAccountedNodesAndTerminals[nodeId];
        for (size_t i = 0; i < terminalsPerNode.size(); i ++){
            int initState = terminalsPerNode[i].first;
            if (!std::count(initStates.begin(), initStates.end(), initState)){
                initStates.push_back(initState);
            }

        }
        it ++;
    }
    return initStates;
}
/*************************************************************************************/
void ComputeChromosomeTransitionsExp::getPosteriorAndExpForNonAccountedFor(map <int, vector<pair<int, int>>>& nonAccountedForBranchesFromFirstRun){
    map <int, vector<pair<int, int>>>::iterator it = nonAccountedForBranchesFromFirstRun.begin();
    while (it != nonAccountedForBranchesFromFirstRun.end()){
        int nodeId = it->first;
        for (size_t i = 0; i < nonAccountedForBranchesFromFirstRun[nodeId].size(); i++){
            std::pair <int, int> currentPairOfAncestralTerminals = nonAccountedForBranchesFromFirstRun[nodeId][i];
            int countForPairOfTerminals = ancestralTerminalsCounts_[nodeId][currentPairOfAncestralTerminals];
            if (countForPairOfTerminals == 0){
                continue;
            }
            for (size_t j = 0; j < ChromosomeSubstitutionModel::NUMTYPES; j++){
                branchTransitionsExp_[nodeId][currentPairOfAncestralTerminals][j] /= ancestralTerminalsCounts_[nodeId][currentPairOfAncestralTerminals];
            }

        }
        it++;
    }

}

/*************************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationAndPosterior(){
    std::map <int, std::map<std::pair<int, int>, int>>::iterator it = ancestralTerminalsCounts_.begin();
    while (it != ancestralTerminalsCounts_.end()){
        int nodeId = it->first;
        std::map<std::pair<int, int>, int> :: iterator iterTerminalStates = ancestralTerminalsCounts_[nodeId].begin();
        while(iterTerminalStates != ancestralTerminalsCounts_[nodeId].end()){
            std::pair <int, int> currentPairOfAncestralTerminals = iterTerminalStates->first;
            int countForPairOfTerminals = ancestralTerminalsCounts_[nodeId][currentPairOfAncestralTerminals];
            if (countForPairOfTerminals == 0){
                iterTerminalStates ++;
                continue;
            }
            for (size_t i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
                branchTransitionsExp_[nodeId][currentPairOfAncestralTerminals][i] /= ancestralTerminalsCounts_[nodeId][currentPairOfAncestralTerminals];
            }

            iterTerminalStates ++;
        }        
        it ++;
    }
}
/***************************************************************************************/
double ComputeChromosomeTransitionsExp::getExpectation(int nodeId, int startAncestral, int endAncestral, int typeOfChange){
    std::pair <int, int> ancestralTerminals;
    ancestralTerminals.first = startAncestral;
    ancestralTerminals.second = endAncestral;
    return branchTransitionsExp_[nodeId][ancestralTerminals][typeOfChange];

}

