#include "ComputeChangesExpectations.h"

using namespace bpp;

/*******************************************************************/
void ComputeChangesExpectations::init(){
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
                for (size_t k = 0; k < alphabet_->getSize(); k++){
                    for (size_t l = 0; l < alphabet_->getSize(); l++){
                        if (jumpProbs_[k][l] == 0){
                            continue; //No need to account for impossible transitions
                        }
                        pairOfpairs terminalsAndJumpsComb;
                        terminalsAndJumpsComb.first = ancestralTerminals;
                        std::pair <int, int> jumpsPair;
                        jumpsPair.first = (int)k;
                        jumpsPair.second = (int)l;
                        terminalsAndJumpsComb.second = jumpsPair;
                        branchTransitionsExp_[nodesIds[n]][terminalsAndJumpsComb] = 0;
                    }
                }

            }
        }
    }
    sort(branchOrder_.begin(), branchOrder_.end(), compareBranches);
}

/************************************************************************************/
bool ComputeChangesExpectations::compareBranches(Node& node1, Node& node2){
    return (node1.getDistanceToFather() < node2.getDistanceToFather());
}
/************************************************************************************/

void ComputeChangesExpectations::runSimulations(int numOfSimulations){
    init();
    for (size_t i = 0; i < alphabet_->getSize(); i++){
        for (size_t j = 0; j < (size_t)numOfSimulations; j++){
            runIteration((int)i);
        }
    }
    computeExpectationAndPosterior();
}
/*************************************************************************************/
void ComputeChangesExpectations::runIteration(int beginState){
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
                pairOfpairs terminalsAndJumpsComb;
                terminalsAndJumpsComb.first = ancestralTerminals;
                terminalsAndJumpsComb.second = jumpsUntilNow[j];
                branchTransitionsExp_[branchOrder_[i].getId()][terminalsAndJumpsComb] += 1;
            }
        }

        int nextState = getRandomState(currentState);
        std::pair<int,int> combOcurredStates;
        combOcurredStates.first = currentState;
        combOcurredStates.second = nextState;
        jumpsUntilNow.push_back(combOcurredStates);
        currentState = nextState;
    }

}
/********************************************************************************/
int ComputeChangesExpectations::getRandomState(int currentState){
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
void ComputeChangesExpectations::computeExpectationAndPosterior(){
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
            for (size_t i = 0; i < alphabet_->getSize(); i++){
                for (size_t j = 0; j < alphabet_->getSize(); j ++){
                    if (jumpProbs_[i][j] == 0){
                        continue; //no need to account for impossible jumps
                    }
                    pairOfpairs terminalsAndJumpsComb;
                    terminalsAndJumpsComb.first = currentPairOfAncestralTerminals;
                    std::pair <int, int> jumpStates;
                    jumpStates.first = (int)i;
                    jumpStates.second = (int)j;
                    terminalsAndJumpsComb.second = jumpStates;
                    branchTransitionsExp_[nodeId][terminalsAndJumpsComb] /= ancestralTerminalsCounts_[nodeId][currentPairOfAncestralTerminals];


                }
            }
            iterTerminalStates ++;
        }        
        it ++;
    }
}
/***************************************************************************************/
double ComputeChangesExpectations::getExpectation(int nodeId, int startAncestral, int endAncestral, int jumpStart, int jumpEnd){
    std::pair <int, int> ancestralTerminals;
    ancestralTerminals.first = startAncestral;
    ancestralTerminals.second = endAncestral;
    std::pair <int, int> jumpStates;
    jumpStates.first = jumpStart;
    jumpStates.second = jumpEnd;
    pairOfpairs terminalsAndJumpStates;
    terminalsAndJumpStates.first = ancestralTerminals;
    terminalsAndJumpStates.second = jumpStates;
    return branchTransitionsExp_[nodeId][terminalsAndJumpStates];

}


    

