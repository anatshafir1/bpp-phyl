#ifndef BPP_PHYL_MAPPING_MULTISTATEMAPPINGPATH_H
#define BPP_PHYL_MAPPING_MULTISTATEMAPPINGPATH_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>
#include <stack>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <functional>
#include <queue>
#include <cmath>
#include <Bpp/Exceptions.h>


using namespace std;
namespace bpp{
  class MultiStateMappingPath{
    public:
     #define EPSILON_THRESHOLD 0.005
      MultiStateMappingPath() {}
      virtual ~MultiStateMappingPath() {}

    public:
      // get in and out degrees of each vertex
      static std::vector<std::pair<size_t, size_t>> getDegreesForEachVertex(vector<vector<size_t>>& graph);
      // checks if an euler circuit can be obtained from the graph (euler circuit iff each vertex has an equal degree of in and out edges)
      static bool isEulerCurcuit(vector<vector<size_t>>& graph);
      // checks if an euler path from start to end can be obtained from the graph (same conditions as for euler circuit with the exception of the source and destination node
      // where the out degree should ecceed the in degree by one and in degree should ecceed the out degree by one, respectively).
      static bool isEulerPath(vector<vector<size_t>>& graph, size_t start, size_t end);
      // get total number of edges in the graph, where a graph G is represented by a vector, where each index represents a node, and at each index there is a vector of neighbor nodes.
      static size_t getNumberOfEdges(vector<vector<size_t>>& graph);
      // For each node get in_degree-out_degree. Positive delta values represent nodes which have more incoming edges. Negative delta values represent those with more outcoming edges.
      static std::vector<int> getDeltaDegrees(std::vector<std::pair<size_t, size_t>> &inAndOutDeg);
      // Given delta degrees values for each node, fill the inbalanced nodes in the relevant vectors inVertices (with more incoming edges) and outVertices (with more outcoming edges).
      static void fillInAndOutOddVertices(std::vector<int>& deltaDegrees, std::vector<size_t> &inVertices, std::vector<size_t> &outVertices);
      // The main function which finds euler path in the graph from source to destination node
      static std::vector <std::pair<size_t, size_t>> getEulerPath(vector<vector<size_t>>& graph, size_t start, size_t end);
      // A function that finds an euler circuit according to Hierholzer's algorithm. Called by getEulerPath()
      static vector<std::pair<size_t, size_t>> findEulerCircuit(vector<vector<size_t>>& graph, size_t start, size_t end);
      // Acording to Hierholzer's algorithm: finds a subcircuit, i.e., a path starting from node v to v, but not necessarily covers the all the edges in the graph
      static std::vector<pair<size_t, size_t>> findSubCircuit(vector<vector<size_t>>& graph, size_t startNode, vector<size_t> &unusedEdges, size_t &numOfEdges, std::map<pair<size_t,size_t>,size_t> &edgesOccurrences);
      // A function called by constructFinalEulerPathFromCircuits() to recursively reconstruct the euler circuit by concatenating all the subcircuits
      static void constructEulerRec(vector<pair<size_t, size_t>> &finalPath, size_t j, vector<size_t> &usedCircuitsIndices, std::vector<vector<pair<size_t,size_t>>> &circuits);
      // A function that aims to construct an euler circuit from all the found subcircuits
      static vector<pair<size_t, size_t>> constructFinalEulerPathFromCircuits(std::vector<vector<pair<size_t,size_t>>> &circuits, size_t start, size_t end);
      
      // a function that finds a path, such that each edge of the graph is visited at least once
      // The algorithm works as following:
      // 1. We calculate the delta degrees of each node (in-out).
      // 2. Odd vertices (unbalanced nodes with unequal in and out degrees) are used to construct a bipartite graph
      //    where possible edges can be put from delta + to delta -. The weights of the edges are the shortest paths between the two nodes
      //    calculated by dijkstra.
      // 3. We find the best matching (combination of edges) in the bipartite graph, such that the addition of
      //    these edges results in delta degree of 0 for all the odd vertices, and also represents the shortest path
      //    if all the found paths are infinite, then we choose the combination of edges with minimal number of infinite edges.
      //    In such cases, where there is at least one infinite edge, we run dijkstra again, but this time not excluding
      //    The artificial edge end->start, so that now we the infinite edges can be reconstructed to a finite path.
      // 4. Now since all the delta degrees are 0, we can find an euler path. All the added edges are then reconstructed to their
      //    respective path, and we can get either one chinese postman path, or several paths, if there were any disjoint paths in the graph.
      
      static vector<vector<pair<size_t, size_t>>> chinesePostman(std::vector<vector<size_t>> &graph, size_t start, size_t end, bool &foundPath);
      // A function that calculates the weights of the edges between the odd vertices
      static std::map<std::pair<size_t, size_t>, double> findWeightsForOddVerticesEdges(std::vector<vector<size_t>> graph, size_t start, size_t end, std::vector<std::pair<size_t, size_t>> &indexToEdges, bool circuitCloseEdgeAdded, std::unordered_map<size_t, unordered_map<size_t, size_t>> &trackPathsFromEachStart);
      // finds the shortest path between a start node and the rest of the vertices
      static std::map<std::pair<size_t, size_t>, double> dijkstra(std::vector<vector<size_t>> &graph, size_t start, std::unordered_map<size_t, size_t> &track);
      // finds the perfect matching between + odd vertices and - odd vertices.
      static std::vector<pair<size_t, size_t>> findPerfectMatching(std::map<std::pair<size_t, size_t>, size_t> &edgesToIndices, std::vector<std::pair<size_t, size_t>> &indexToEdge, std::vector<double> &weightsPerIndices, std::vector<int> &deltaDegrees, bool &isPahInfinite);
      static void findPerfectMatchingRec(vector<vector<size_t>> &includedEdges, std::vector<size_t> &combinationOfBits, size_t startIndex, std::vector<std::pair<size_t, size_t>> &indexToEdge, std::vector<int> &deltaDegrees, std::map<size_t, int> &leftOptions);
      // calculates the sum of weights of the added vertices
      static double calculateSumOfPath(vector<size_t> &path, std::vector<double> &weightsPerIndices);
      // counts the number of infinite edges.
      static void countNumberOfInfiniteEdges(std::vector <size_t> &infinitePath, size_t &numOfInfiniteEdges, std::vector<double> &weightsPerIndices);
      static void fillRelativeTimeDuration(std::unordered_map<size_t, double> &relativeTimeDuration, vector<double> &dwellingTimes, double totalDurationTime);
      static vector<vector<size_t>> createGraphFromPath(vector<pair<size_t, size_t>> &disjointPath, std::unordered_map<size_t, size_t> &indicesToNodesSub, std::unordered_map<size_t, size_t> &nodesToIndicesSub, size_t originalGraphSize);
      static vector<vector<size_t>> createGraphForChinesePostman(std::unordered_map<size_t, size_t> &indicesToNodes, std::unordered_map<size_t, size_t> &nodesToIndices, std::unordered_map<size_t, double> &relativeTimeDuration, std::map<std::pair<size_t, size_t>, double> &transitions);
      // TSP implementation
      static vector<size_t> TSP(size_t start, size_t end, std::map<std::pair<size_t, size_t>, double> &transitions, double totalDurationTime, bool &validPath, std::unordered_map<size_t, double> &relativeTimeDuration);
      static vector<size_t> findExpectedMappingPathForEachNode(size_t start, size_t end, std::map<std::pair<size_t, size_t>, double> &transitions, vector<double> &dwellingTimes, double totalDurationTime, bool &foundPath);
      static std::unordered_map<size_t, vector<size_t>> createEdges(std::unordered_map<size_t, double> &vertices, std::map<std::pair<size_t, size_t>, double> &transitions);
      static std::vector<size_t> decimalToBinaryPowers(int decimalNumber);
      static void findBestPath(std::pair<size_t,size_t> &bestCandidatePathId, std::map<std::pair<size_t, size_t>, double> &paths, size_t desiredPathId, std::unordered_map<size_t, vector<size_t>> &edges, std::map<std::pair<size_t, size_t>, double> &transitions, size_t end, bool &foundPath);
      static void reconstructBestPath(std::vector<size_t> &bestPath, size_t lengthOfPath, std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> &pathReconstruction, std::pair<size_t,size_t> bestCandidatePathId, size_t start, size_t end);
      static void printGraph(vector<vector<size_t>>& graph);


  };

}
#endif // BPP_PHYL_MAPPING_MULTISTATEMAPPINGPATH_H