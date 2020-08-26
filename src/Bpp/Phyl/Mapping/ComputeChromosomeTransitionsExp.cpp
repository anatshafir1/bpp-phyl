#include "ComputeChromosomeTransitionsExp.h"
using namespace bpp;

/****************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationOfChangePerBranch(int nodeId, int jumpStateStart, int jumpStateEnd, VVdouble &jointProbFatherNode, int transitionType){
    int jumpType;
    if (transitionType == -1){
        if (jumpTypeMethod_ == 0){
            jumpType = getTypeOfTransition(jumpStateStart, jumpStateEnd);
        }else{
            jumpType = getTypeOfTransitionWithProb(jumpStateStart, jumpStateEnd);
        }
        if (jumpType == ChromosomeSubstitutionModel::ILLEGAL){
            throw Exception("Ilegal transition found!");
            return;
        }

    }else{
        jumpType = transitionType;
    }

    double expectation = 0;
    for (size_t x = 0; x < alphabet_->getSize(); x ++){
        for (size_t y = 0; y < alphabet_->getSize(); y++){
            expectation += jointProbFatherNode[x][y] * simExpectations_->getExpectation(nodeId, (int)y, (int)x, jumpStateStart, jumpStateEnd);

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
        if ((chrStart % 2 == 0) && (chrEnd == chrStart * 1.5)){
            legalMove = true;
            jumpType.push_back(ChromosomeSubstitutionModel::DEMIDUPL_T);
        }else if ((chrStart % 2 != 0) && (chrEnd == (int)ceil(chrStart * 1.5))){
            legalMove = true;
            jumpType.push_back(ChromosomeSubstitutionModel::DEMIDUPL_T);
        }else if ((chrStart % 2 != 0) && (chrEnd == (int)floor(chrStart * 1.5))){
            legalMove = true;
            jumpType.push_back(ChromosomeSubstitutionModel::DEMIDUPL_T);
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
    std::vector<double> weights;
    std::vector<size_t> indices;
    std::vector<size_t> result;
    result.push_back(0);
    double sumOfRates = 0;
    if (jumpType.size() == 1){
        return jumpType[0];
    }else{
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

        }
        for (size_t j = 0; j < weights.size(); j++){
            weights[j] /= sumOfRates;
        }
        
        RandomTools::getSample(indices, weights, result);
        
        
    }
    return jumpType[result[0]];    
}
/**********************************************************************************/
void ComputeChromosomeTransitionsExp::computeExpectationPerType(){
    init();
    vector<int> nodeIds = tree_->getNodesId();
    for (size_t n = 0; n < nodeIds.size(); n++){
        if (tree_->getRootId() == nodeIds[n]){
            continue;
        }
        VVdouble jointProbFatherNode = ancr_->getAllJointFatherNodeProbabilities()[nodeIds[n]][0];
        for (int jumpStateStart = alphabet_->getMin(); jumpStateStart <= alphabet_->getMax(); jumpStateStart++){
            //base number transitions
            //Remember to account for an option where dupl is set to be ignored!
            computeBaseNumExpectation(nodeIds[n], jumpStateStart, jointProbFatherNode);
            //gain transitions
            if (jumpStateStart < alphabet_->getMax()){
                if ((model_->getConstGain() != IgnoreParam) || (model_->getChangeRateGain() != IgnoreParam)){
                    computeExpectationOfChangePerBranch(nodeIds[n], jumpStateStart - alphabet_->getMin(), (jumpStateStart +1) - alphabet_->getMin(), jointProbFatherNode);
                }
            }
            //loss transitions
            if (jumpStateStart > alphabet_->getMin()){
                if ((model_->getConstLoss() != IgnoreParam) || (model_->getChangeRateLoss() != IgnoreParam)){
                    computeExpectationOfChangePerBranch(nodeIds[n], jumpStateStart - alphabet_->getMin(), (jumpStateStart - 1) - alphabet_->getMin(), jointProbFatherNode);
                }

            }
            // duplication transitions
            if ((jumpStateStart * 2 <= alphabet_->getMax()) && (jumpStateStart > 1)){
                if ((model_->getConstDupl() != IgnoreParam) || (model_->getChangeRateDupl() != IgnoreParam)){
                    computeExpectationOfChangePerBranch(nodeIds[n], jumpStateStart - alphabet_->getMin(), (jumpStateStart * 2) - alphabet_->getMin(), jointProbFatherNode);
                }

            }
            // demi-duplication transitions
            computeDemiDuplExpectation(nodeIds[n], jumpStateStart, jointProbFatherNode);
            // transitions to max chromosome number
            computeExpectationOfChangePerBranch(nodeIds[n], jumpStateStart - alphabet_->getMin(), alphabet_->getMax() - alphabet_->getMin(), jointProbFatherNode, ChromosomeSubstitutionModel::MAXCHR_T);

            // for (size_t jumpStateEnd = 0; jumpStateEnd < alphabet_->getSize(); jumpStateEnd++){
            //     if (jumpStateStart == jumpStateEnd){
            //         continue;
            //     }
            //     if (model_->Qij(jumpStateStart, jumpStateEnd) == 0){
            //         continue;
            //     }
            //     computeExpectationOfChangePerBranch(nodeIds[n], (int)jumpStateStart, (int)jumpStateEnd, jointProbFatherNode);

            // }
        }
        
    }

}
/**************************************************************************************/
void ComputeChromosomeTransitionsExp::computeBaseNumExpectation(int nodeId, int jumpStateStart, VVdouble jointProbFatherNode){
    int baseNumber = model_->getBaseNumber();
    if (baseNumber != IgnoreParam){
        int jumpStateEnd = jumpStateStart + baseNumber;
        while(jumpStateEnd <= alphabet_->getMax()){
            if ((double)jumpStateEnd/2 == (double)jumpStateStart){
                jumpStateEnd += baseNumber;
                continue;
            }
            computeExpectationOfChangePerBranch(nodeId, jumpStateStart - alphabet_->getMin(), jumpStateEnd - alphabet_->getMin(), jointProbFatherNode);
            jumpStateEnd += baseNumber;

        }

    }
}
/**************************************************************************************/
void ComputeChromosomeTransitionsExp::computeDemiDuplExpectation(int nodeId, int jumpStateStart, VVdouble jointProbFatherNode){
    if ((model_->getDemiDupl() != IgnoreParam) && (jumpStateStart > 2)){
        if ((jumpStateStart % 2 == 0) && (jumpStateStart * 1.5 <= alphabet_->getMax())){
            computeExpectationOfChangePerBranch(nodeId, jumpStateStart - alphabet_->getMin(), (int)(jumpStateStart * 1.5) - alphabet_->getMin(), jointProbFatherNode);
        }else if (jumpStateStart % 2 != 0){
            if ((int)ceil(jumpStateStart * 1.5) <= alphabet_->getMax()){
                computeExpectationOfChangePerBranch(nodeId, jumpStateStart - alphabet_->getMin(), (int)ceil(jumpStateStart * 1.5) - alphabet_->getMin(), jointProbFatherNode);
            }
            if (((int)floor(jumpStateStart * 1.5) <= alphabet_->getMax()) && (jumpStateStart != 3)){
                computeExpectationOfChangePerBranch(nodeId, jumpStateStart - alphabet_->getMin(), (int)floor(jumpStateStart * 1.5) - alphabet_->getMin(), jointProbFatherNode);

            }

        }
    }
}
/**************************************************************************************/
void ComputeChromosomeTransitionsExp::init(){
    vector<int> nodeIds = tree_->getNodesId();
    for (size_t n = 0; n < nodeIds.size(); n++){
        if (tree_->getRootId() == nodeIds[n]){
            continue;
        }
        for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
            expNumOfChangesPerBranch_[nodeIds[n]][i] = 0;
        }
        
    }
    for (int i = 0; i < ChromosomeSubstitutionModel::NUMTYPES; i ++){
        expNumOfChanges_[i] = 0;
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
    }

}

