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


int main() {
    return 0;
}

// std::vector<std::pair<double, size_t>> dijkstra(std::map<size_t, vector<size_t>> &graph, std::map<std::pair<size_t, size_t>, double> &weights, size_t numberOfNodes, size_t startNode){
//     std::vector<std::pair<double, size_t>> shortestPathInfo;
//     shortestPathInfo.resize(numberOfNodes);
//     size_t nodeMin = startNode;
//     std::vector<bool> usedVertices;
//     size_t numOfUpdated = 0;
//     for (size_t i = 0; i < numberOfNodes; i++){
//         shortestPathInfo[i].second = startNode;
//         if (i == startNode){
//             shortestPathInfo[i].first = 0;
//             continue;
//         }
//         shortestPathInfo[i].first = std::numeric_limits<double>::infinity();
//         usedVertices[i] = false;       
//     }
//     while (numOfUpdated < numOfNodes){
//         auto neighbors = graph[nodeMin];
//         for (size_t j = 0; j < neighbors.size(); j++){
//             auto edge = std::pair<size_t, size_t>(nodeMin, neighbors[j])
//             double candidatePathDistance = shortestPathInfo[nodeMin].first + weights[edge];
//             if (candidatePathDistance < shortestPathInfo[neighbors[j]].first){
//                 shortestPathInfo[neighbors[j]].first = candidatePathDistance;
//                 shortestPathInfo[neighbors[j]].second = nodeMin;
//             }

//         }
//         numOfUpdated ++;
//     }


// }
