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
void ComputeChromosomeTransitionsExp::printResults() {
    std::cout << "**********************\n" << endl;
    std::map <int, string> jumpTypeToString;
    std::vector<int> nodesIds = tree_->getNodesId();
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
        std::cout << nodeName <<":" << endl;
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
            string jumpType;
            if (i == ChromosomeSubstitutionModel::GAIN_T){
                jumpType = "\tGain Expectation: ";
            }else if (i == ChromosomeSubstitutionModel::LOSS_T){
                jumpType = "\tLoss Expectation: ";
            }else if (i == ChromosomeSubstitutionModel::DUPL_T){
                jumpType = "\tDupl Expectation: ";
            }else if (i == ChromosomeSubstitutionModel::DEMIDUPL_T){
                jumpType = "\tDemi-Dupl Expectation: ";
            }else if (i == ChromosomeSubstitutionModel::BASENUM_T){
                jumpType = "\tBaseNumber Expectation: ";
            }else{
                jumpType = "\tMaxChr Expectation: ";
            }
            std::cout << jumpType << expNumOfChangesPerBranch_[nodesIds[n]][i] <<endl;
            std::cout << "+++++" <<endl;
            jumpTypeToString[i] = jumpType;
        }
        
    }
    std::cout << "********************************************"<<endl;
    std::cout <<"Total Expectations:" << endl;
    // print total expectations
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i++){
        std::cout <<jumpTypeToString[i] << expNumOfChanges_[i] <<endl;
        if ((i == ChromosomeSubstitutionModel::MAXCHR_T) && (expNumOfChanges_[i] > 0)){
            std::cout <<"Note: Max chr transitions exist-> consider to increase the maximum possible chromosome mumber!"<<endl;
        }
    }

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
void ComputeChromosomeTransitionsExp::runIteration(int beginState){
    double totalTimeTillJump = 0;
    double maxBranch = branchOrder_[branchOrder_.size()-1].getDistanceToFather();
    int currentState = beginState;
    vector <std::pair<int,int>> jumpsUntilNow;
    int indexOfSmallestNotUpdatedBranch = 0;
    while (totalTimeTillJump < maxBranch){
        double averageWaitingTime = 1 / waitingTimes_[currentState];
        totalTimeTillJump +=  RandomTools::randExponential(averageWaitingTime);
        for (size_t i = indexOfSmallestNotUpdatedBranch; i < branchOrder_.size(); i++){
            if (branchOrder_[i].getDistanceToFather() > totalTimeTillJump){
                indexOfSmallestNotUpdatedBranch = (int)i;
                break;
            }
            std::pair<int, int> ancestralTerminals;
            ancestralTerminals.first = beginState;
            ancestralTerminals.second = currentState;
            ancestralTerminalsCounts_[branchOrder_[i].getId()][ancestralTerminals] += 1;
            for (size_t j = 0; j < jumpsUntilNow.size(); j++){
                //ChromosomeSubstitutionModel::typeOfTransition typeOfJump;
                updateExpectationsPerBranch(branchOrder_[i].getId(), ancestralTerminals, jumpsUntilNow[j]);
                // if (jumpTypeMethod_ == 0){
                //     typeOfJump = getTypeOfTransition(jumpsUntilNow[j].first, jumpsUntilNow[j].second);
                // }else{
                //     typeOfJump = getTypeOfTransitionWithProb(jumpsUntilNow[j].first, jumpsUntilNow[j].second);
                // }
                
                // branchTransitionsExp_[branchOrder_[i].getId()][ancestralTerminals][typeOfJump] += 1;
            }
        }

        int nextState = getRandomState(currentState);
        std::pair<int,int> combOcurredStates;
        combOcurredStates.first = currentState;
        combOcurredStates.second = nextState;
        jumpsUntilNow.push_back(combOcurredStates);
        updateMapOfJumps(currentState, nextState);
        currentState = nextState;
    }

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
            return;
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
        branchTransitionsExp_[nodeId][ancestralTerminals][typeOfJump] += stateJumpTypeProb_[jumpStates][typeOfJump];
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

