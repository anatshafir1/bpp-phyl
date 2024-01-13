// Tester for StochasticMapping implementation

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>

// From bpp-core:
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h> 
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/AlphabetIndex/UserAlphabetIndex1.h>
#include <Bpp/Seq/Alphabet/NumericAlphabet.h>

// From bpp-phyl
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Legacy/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/G2001.h>
#include <Bpp/Phyl/Model/TwoParameterBinarySubstitutionModel.h>
#include <Bpp/Phyl/Legacy/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Legacy/OptimizationTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Simulation/DetailedSiteSimulator.h>
//#include <Bpp/Phyl/Simulation/NonHomogeneousSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>
#include <Bpp/Phyl/Mapping/RewardMappingTools.h>
#include <Bpp/Phyl/Mapping/Reward.h>
#include <Bpp/Phyl/Mapping/DecompositionReward.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>




using namespace bpp;
using namespace std;



/******************************************************************************************************/

std::vector<size_t> decimalToBinaryPowers(int decimalNumber) {
    std::vector<size_t> powers;
    int power = 0;
    while (decimalNumber > 0) {
        if (decimalNumber % 2 == 1) {
            powers.push_back((size_t)(power));
        }

        decimalNumber /= 2;
        ++power;
    }

    return powers;
}
/******************************************************************************************************/
std::map<size_t, vector<size_t>> createEdges(std::map<size_t, double> &vertices, std::map<std::pair<size_t, size_t>, double> &transitions){
  std::map<size_t, vector<size_t>> edges;
  auto it = transitions.begin();
  while (it != transitions.end()){
    auto transition = it->first;
    if (transitions[transition] <= 0){
      it ++;
      continue;

    }
    size_t outgoing = transition.first;
    if (vertices[outgoing] < EPSILON_THRESHOLD){
        it ++;
        continue;
    }
    size_t incoming = transition.second;
    if (vertices[incoming] < EPSILON_THRESHOLD){
        it ++;
        continue;
    }
    if (edges.find(outgoing) == edges.end()){
      edges[outgoing];
      
    }
    edges[outgoing].push_back(incoming);
    it ++;
  }
  return edges;
}
/******************************************************************************/
void findBestPath(std::pair<size_t,size_t> &bestCandidatePathId, std::map<std::pair<size_t, size_t>, double> &paths, size_t desiredPathId, std::map<size_t, vector<size_t>> &edges, std::map<std::pair<size_t, size_t>, double> &transitions, size_t end, bool &foundPath){
  std::pair<size_t,size_t> candidatePathId;
  double weightBest = 0;
  auto itPath = paths.begin();
  vector<size_t> neighbors;
  while (itPath != paths.end()){
    size_t pathId = (itPath->first).first;
    if (pathId != desiredPathId){
      itPath++;
      continue;
    }
    size_t lastMid = (itPath->first).second;
    neighbors = edges[lastMid];
    double weight;
    
    if (std::find(neighbors.begin(), neighbors.end(), end) != neighbors.end()){
      candidatePathId = itPath->first;
      std::pair<size_t, size_t> lastEdge(lastMid, end);
      weight = paths[candidatePathId] + transitions[lastEdge];
      foundPath = true;
      if (weight > weightBest){
        bestCandidatePathId = itPath->first;
        weightBest = weight;

      }
    }
    itPath++;
  }
}
/******************************************************************************/
void reconstructBestPath(std::vector<size_t> &bestPath, size_t lengthOfPath, std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> &pathReconstruction, std::pair<size_t,size_t> bestCandidatePathId, size_t start, size_t end){
  bestPath.resize(lengthOfPath);
  std::pair<size_t, size_t> fatherPath;
  size_t currState;
  for (size_t k = 1; k <= lengthOfPath; k++){
    size_t reverseIndex = lengthOfPath-k;
    if (k == 1){
      bestPath[reverseIndex] = end;
    }else if (k == lengthOfPath){
      bestPath[reverseIndex] = start;
    }else{
      if (k  == 2){
        fatherPath = bestCandidatePathId;
      }else{
        fatherPath = pathReconstruction[fatherPath];
      }
      currState = fatherPath.second;
      bestPath[reverseIndex] = currState;

    }
  }

}

/******************************************************************************/
vector<size_t> findExpectedMappingPathForEachNode(size_t start, size_t end, std::map<std::pair<size_t, size_t>, double> &transitions, vector<double> &dwellingTimes, double totalDurationTime){
  std::map<size_t, double> relativeTimeDuration;
  std::vector<size_t> bestPath;
  size_t desiredPathId = 0;
  
  size_t numOfNotAdded = 0;
  for (size_t i = 0; i < dwellingTimes.size(); i++){
    relativeTimeDuration[i] = dwellingTimes[i]/totalDurationTime;
    if (relativeTimeDuration[i] < EPSILON_THRESHOLD){
      numOfNotAdded ++;
    }else{
      if (i != end){
        desiredPathId += std::pow(2, i);
      }else{
        if (start == end){
          desiredPathId += std::pow(2, i);
        }
      }
    }
  }
  auto edges = createEdges(relativeTimeDuration, transitions);
  size_t pathLength = dwellingTimes.size()-numOfNotAdded-1; // we don't include the final end node.
  
  std::map<std::pair<size_t, size_t>, double> paths;
  std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> pathReconstruction;
  auto neighbors = edges[start];
  if (neighbors.size() == 0){  
    //bestPath.push_back(end);
    return bestPath;
  }
  if ((pathLength == 1) && (start != end)){
    if (std::find(neighbors.begin(), neighbors.end(), end) != neighbors.end()){
      bestPath.push_back(start);
      bestPath.push_back(end);
      return bestPath;
    }
    throw Exception("StochasticMapping::findExpectedMappingPathForEachNode: No such path!");
  }
  for (size_t i = 0; i < neighbors.size(); i++){
    if (neighbors[i] == end){
      continue;
    }
    size_t pathId = std::pow(2, start) + std::pow(2, neighbors[i]);
    size_t fatherPathId = 0;
    size_t endOfFather = start;
    std::pair<size_t, size_t> fatherPathWithEnd(fatherPathId, endOfFather);

    
    std::pair<size_t, size_t> pathWithEnd(pathId, neighbors[i]);
    std::pair<size_t, size_t> edge(start, neighbors[i]);
    paths[pathWithEnd] = transitions[edge];
    pathReconstruction[pathWithEnd] = fatherPathWithEnd;
  }
  if (start == end){
    pathLength -= 1; // path length excludes start and an additional neighbor but not the dest node, since start and end are the same.

  }else{
    pathLength -= 2; // path length excludes start and an additional neighbor in addition to the dest node

  }
  vector<std::pair<size_t, size_t>> pathsIdsOfInitialLength;
  auto it = paths.begin();
  while (it != paths.end()){
    size_t id = (it->first).first;
    if (id != 0){
        pathsIdsOfInitialLength.push_back(it->first);
    }

    it ++;
  }
  auto pathsIdsOfGivenLength = pathsIdsOfInitialLength;
  
  for (size_t i = 0; i < pathLength; i++){
    vector<std::pair<size_t, size_t>> pathsIdsOfCurrentLength;
    for (size_t j = 0; j < pathsIdsOfGivenLength.size(); j++){
      size_t fatherIdPath = pathsIdsOfGivenLength[j].first;
      size_t newStart = pathsIdsOfGivenLength[j].second;
      neighbors = edges[newStart];
      for (size_t k = 0; k< neighbors.size(); k++){
        if (neighbors[k] == end){
          continue;
        }
        auto used = decimalToBinaryPowers(fatherIdPath);
        if (std::find(used.begin(), used.end(), neighbors[k]) != used.end()){
          continue;
        }
        size_t currentPathId = fatherIdPath + std::pow(2, neighbors[k]);

        std::pair<size_t, size_t> currentPathWithEnd(currentPathId, neighbors[k]);
        std::pair<size_t, size_t> edge(newStart, neighbors[k]);
        double weight = paths[pathsIdsOfGivenLength[j]] + transitions[edge];
        if (paths.find(currentPathWithEnd) != paths.end()){
          if (weight <= paths[currentPathWithEnd]){
            continue;
          }
        }else{
          paths[currentPathWithEnd] = weight;
          pathReconstruction[currentPathWithEnd] = pathsIdsOfGivenLength[j];
          pathsIdsOfCurrentLength.push_back(currentPathWithEnd);
        }
        
      }

    }
    pathsIdsOfGivenLength = pathsIdsOfCurrentLength;
  }
  // Find the best path
  
  bool foundPath = false;
  std::pair<size_t,size_t> bestCandidatePathId; 
  findBestPath(bestCandidatePathId, paths, desiredPathId, edges, transitions, end, foundPath);

  // reconstruct the best path
  size_t lengthOfPath = dwellingTimes.size()-numOfNotAdded;
  if (start == end){
    lengthOfPath ++;
  }
  if (!foundPath){
    throw Exception("Path does not exist!");
  }
  reconstructBestPath(bestPath, lengthOfPath, pathReconstruction, bestCandidatePathId, start, end);
  return bestPath;


}
void test3(){
  size_t start = 1;
  size_t end = 1;
  size_t numOfStates = 2;
  map<pair<size_t, size_t>, double> transitions;
  pair<size_t, size_t> transition1(0,1);
  transitions[transition1] = 0.9;
  pair<size_t, size_t> transition2(1,0);
  transitions[transition2] = 0.6;
  vector<double> dwellingTimes;
  dwellingTimes.resize(numOfStates);
  dwellingTimes[0] = 0.3;
  dwellingTimes[1] = 0.2;


  double totalDurationTime = 0.5;
  auto path = findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }
}

void test4(){
  size_t start = 1;
  size_t end = 1;
  size_t numOfStates = 2;
  map<pair<size_t, size_t>, double> transitions;
  pair<size_t, size_t> transition1(0,1);
  transitions[transition1] = 0.9;
  pair<size_t, size_t> transition2(1,0);
  transitions[transition2] = 0.6;
  vector<double> dwellingTimes;
  dwellingTimes.resize(numOfStates);
  dwellingTimes[0] = 0.00001;
  dwellingTimes[1] = 0.2;


  double totalDurationTime = 0.20001;
  auto path = findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }
}

void test1(){
  size_t start = 1;
  size_t end = 3;
  size_t numOfStates = 4;
  map<pair<size_t, size_t>, double> transitions;
  for (size_t i = 0; i < 4; i++){
    for (size_t j = 0; j < 4; j++){
      if (i == j){
        continue;
      }
      pair<size_t, size_t> transition(i,j);
      if ((i == 0) &&(j == 1)){
        transitions[transition] = 1.2;
      }
      else if ((i == 0) && (j == 2)){
        transitions[transition] = 0.6;
      }
      else if ((i == 0) && (j == 3)){
        transitions[transition] = 0.3;
      }
      else if ((i == 1) && (j == 0)){
        transitions[transition] = 0.2;
      }
      else if ((i == 1) && (j == 2)){
        transitions[transition] = 0.4;
      }
      else if ((i == 1) && (j == 3)){
        transitions[transition] = 1.2;
      }
      else if ((i == 2) && (j == 0)){
        transitions[transition] = 0.4;

      }
      else if ((i == 2) && (j == 3)){
        transitions[transition] = 0.9;
      }
      else{
        transitions[transition] = 0;

      }
          
    }
  }

  vector<double> dwellingTimes;
  dwellingTimes.resize(numOfStates);
  dwellingTimes[0] = 0.3;
  dwellingTimes[1] = 0.2;
  dwellingTimes[2] = 0.1;
  dwellingTimes[3] = 0.5;

  double totalDurationTime = 1.1;
  auto path = findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }

}


void test2(){
  size_t start = 1;
  size_t end = 1;
  size_t numOfStates = 4;
  map<pair<size_t, size_t>, double> transitions;
  for (size_t i = 0; i < 4; i++){
    for (size_t j = 0; j < 4; j++){
      if (i == j){
        continue;
      }
      pair<size_t, size_t> transition(i,j);
      if ((i == 0) &&(j == 1)){
        transitions[transition] = 1.2;
      }
      else if ((i == 0) && (j == 2)){
        transitions[transition] = 0.6;
      }
      else if ((i == 0) && (j == 3)){
        transitions[transition] = 0.3;
      }
      else if ((i == 1) && (j == 0)){
        transitions[transition] = 0.2;
      }
      else if ((i == 1) && (j == 2)){
        transitions[transition] = 0.4;
      }
      else if ((i == 1) && (j == 3)){
        transitions[transition] = 1.2;
      }
      else if ((i == 2) && (j == 0)){
        transitions[transition] = 0.4;

      }
      else if ((i == 2) && (j == 3)){
        transitions[transition] = 0.9;
      }
      else if ((i == 3) && (j == 1)){
        transitions[transition] = 0.8;
      }
      else if ((i ==3) && (j == 0)){
        transitions[transition] = 0.1;
      }
      else{
        transitions[transition] = 0;

      }
          
    }
  }

  vector<double> dwellingTimes;
  dwellingTimes.resize(numOfStates);
  dwellingTimes[0] = 0.3;
  dwellingTimes[1] = 0.2;
  dwellingTimes[2] = 0.1;
  dwellingTimes[3] = 0.5;

  double totalDurationTime = 1.1;
  auto path = findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }

}


int main(){
  std::cout << "Test 1" << std::endl;
  test1();
  std::cout << "Test 2" << std::endl;
  test2();
  std::cout << "Test 3" << std::endl;
  test3();
  std::cout << "Test 4" << std::endl;
  test4();

  return 0;



}

