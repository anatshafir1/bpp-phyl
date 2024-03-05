#include <Bpp/Phyl/Mapping/MultiStateMappingPath.h>



using namespace std;
using namespace bpp;

// vector<size_t> dfs(vector<vector<size_t>>& graph, size_t startNode) {
//     std::vector<size_t> result;
//     size_t numNodes = graph.size();
//     vector<bool> visited(numNodes, false);
//     stack<size_t> s;

//     // Push the start node onto the stack
//     s.push(startNode);

//     while (!s.empty()) {
//         // Pop a node from the stack
//         size_t currentNode = s.top();
//         s.pop();

//         // Visit the node if not already visited
//         if (!visited[currentNode]) {
//             result.push_back(currentNode);
//             std::cout << " " << currentNode;
//             visited[currentNode] = true;
//         }

//         // Push adjacent nodes onto the stack
//         for (size_t neighbor : graph[currentNode]) {
//             if (!visited[neighbor]) {
//                 s.push(neighbor);
//             }
//         }
//     }
//     std::cout << endl;
//     return result;
// }
// /**********************************************************************************************************/
// vector<vector<size_t>> transposeGraph(vector<vector<size_t>>& graph){
//     std::vector<vector<size_t>> transposed;
//     size_t graphSize = graph.size();
//     for (size_t i = 0; i < graphSize; i++){
//         std::vector<size_t> edges;
//         transposed.push_back(edges);
//     }
//     for (size_t i = 0; i < graphSize; i++){
//         auto& neighbors = graph[i];
//         for (size_t j = 0; j < neighbors.size(); j++){
//             transposed[neighbors[j]].push_back(i);
//         }
//     }
//     return transposed;

// }
// /**********************************************************************************************************/

// /**********************************************************************************************************/
// bool isStronglyConnected(vector<vector<size_t>>& graph){
//     size_t graphSize = graph.size();
//     auto foundVertices = dfs(graph, 0);
//     if (foundVertices.size() < graphSize){
//         return false;
//     }
//     auto reversedGraph = transposeGraph(graph);
//     auto foundVerticeInRev = dfs(reversedGraph, 0);
//     if (foundVerticeInRev.size() < graphSize){
//         return false;
//     }
//     return true;

// }
// /**********************************************************************************************************/
// std::vector<std::pair<size_t, size_t>> getDegreesForEachVertex(vector<vector<size_t>>& graph){
//     size_t graphSize = graph.size();
//     std::vector<std::pair<size_t, size_t>> degrees(graphSize, {0, 0});
//     for (size_t i = 0; i < graphSize; i++){
//         for (size_t j = 0; j < graph[i].size(); j++){
//             degrees[i].second += 1; // out degree
//             degrees[graph[i][j]].first += 1; // in degree
//         }
//     }
//     return degrees;
// }
// /**********************************************************************************************************/
// bool isEulerCurcuit(vector<vector<size_t>>& graph){
//     std::vector<std::pair<size_t, size_t>> degrees = getDegreesForEachVertex(graph);
//     for (size_t i = 0; i < degrees.size(); i++){
//         if (degrees[i].first != degrees[i].second){
//             return false;
//         }
//     }
//     return true;
// }
// /**********************************************************************************************************/
// bool isEulerPath(vector<vector<size_t>>& graph, size_t start, size_t end){
//     std::vector<std::pair<size_t, size_t>> degrees = getDegreesForEachVertex(graph);
//     for (size_t i = 0; i < degrees.size(); i++){
//         if (start == i){
//             if (degrees[i].second != degrees[i].first + 1){
//                 return false;

//             }
//         }else if (end == i){
//             if (degrees[i].first != degrees[i].second + 1){
//                 return false;

//             }
//         }else{
//             if (degrees[i].first != degrees[i].second){
//                 return false;
//             }   
//         }
//     }
//     return true;
// }
// /**********************************************************************************************************/
// size_t getNumberOfEdges(vector<vector<size_t>>& graph){
//     size_t numOfEdges = 0;
//     size_t graphSize = graph.size();
//     for (size_t i = 0; i < graphSize; i++){
//         numOfEdges += graph[i].size();
//     }
//     return numOfEdges;

// }
// /**********************************************************************************************************/
// std::vector<pair<size_t, size_t>> findSubCircuit(vector<vector<size_t>>& graph, size_t startNode, vector<size_t> &unusedEdges, size_t &numOfEdges, std::map<pair<size_t,size_t>,size_t> &edgesOccurrences){
//     size_t currentNode = startNode; 
//     size_t nextNode = graph[startNode][0];
//     vector<pair<size_t, size_t>> circuit;    
//     while (nextNode != startNode){
//         auto neighbors = graph[currentNode];
//         for (size_t i = 0; i < neighbors.size(); i++){
//             std::pair<size_t, size_t> edge(currentNode, neighbors[i]);
//             if (edgesOccurrences[edge] == 0){
//                 continue;
//             }
//             unusedEdges[currentNode]--;
//             numOfEdges--;
//             circuit.push_back(std::pair<size_t,size_t>(currentNode, neighbors[i]));
//             edgesOccurrences[edge] --;
//             currentNode = neighbors[i];
            
//             break; // we choose only one of the neighbors to continue with
            

//         }
//         nextNode = currentNode;
//     }
//     return circuit;

// }
// /**********************************************************************************************************/
// std::vector<int> getDeltaDegrees(std::vector<std::pair<size_t, size_t>> &inAndOutDeg){
//     std::vector<int> deltaDegrees;
//     deltaDegrees.resize(inAndOutDeg.size());
//     size_t inDeg;
//     size_t outDeg;
//     for (size_t i = 0; i < inAndOutDeg.size(); i++){
//         inDeg = inAndOutDeg[i].first;
//         outDeg = inAndOutDeg[i].second;
//         deltaDegrees[i] = static_cast<int>(inDeg) - static_cast<int>(outDeg);

//     }
//     return deltaDegrees;

// }
// /**********************************************************************************************************/
// void fillInAndOutOddVertices(std::vector<int>& deltaDegrees, std::vector<size_t> &inVertices, std::vector<size_t> &outVertices){
//     for (size_t i = 0; i < deltaDegrees.size(); i++){
//         if (deltaDegrees[i] < 0){
//             outVertices.push_back(i);
//         }else if (deltaDegrees[i] > 0){
//             inVertices.push_back(i);
//         }
//     }
//     return;
// }

// /**********************************************************************************************************/
// std::map<std::pair<size_t, size_t>, double> dijkstra(std::vector<vector<size_t>> &graph, size_t start, std::unordered_map<size_t, size_t> &track){
//     vector<pair<size_t, double>> nodesWithdistances;
//     nodesWithdistances.resize(graph.size());
//     for (size_t i = 0; i < nodesWithdistances.size(); i++){
//         nodesWithdistances[i].first = i;
//         if (i == start){
//             nodesWithdistances[i].second = 0;
//         }else{
//             nodesWithdistances[i].second = std::numeric_limits<double>::infinity();
//         }
//     }

//     vector<double> distances(graph.size(), std::numeric_limits<double>::infinity());
//     distances[start] = 0;
//     std::map<std::pair<size_t, size_t>, double> shortestPaths;
//     std::priority_queue<std::pair<size_t, double>, std::vector<std::pair<size_t, double>>, std::greater<std::pair<size_t, double>>> minHeap;
//     minHeap.push({start, 0});
    
//     while (!minHeap.empty()) {
//         size_t node = minHeap.top().first;
//         double nodeDist = minHeap.top().second;
//         minHeap.pop();
//         if (nodeDist > distances[node]) {
//             continue;
//         }
//         for (auto &neighbor : graph[node]) {
//             if (distances[node] + 1 < distances[neighbor]) {
//                 distances[neighbor] = distances[node] + 1;
//                 minHeap.push({neighbor, distances[neighbor]});
//                 track[neighbor] = node;
//             }
//         }

//     }

//     for (size_t i = 0; i < distances.size(); i++){
//         shortestPaths[std::pair<size_t, size_t>(start, i)] = distances[i];
//     }
//     return shortestPaths;

// }
// /**********************************************************************************************************/
// // the parameter graph is by value on purpose, because I'm going to change the object here
// std::map<std::pair<size_t, size_t>, double> findWeightsForOddVerticesEdges(std::vector<vector<size_t>> graph, size_t start, size_t end, std::vector<std::pair<size_t, size_t>> &indexToEdges, bool circuitCloseEdgeAdded, std::unordered_map<size_t, unordered_map<size_t, size_t>> &trackPathsFromEachStart){
//     std::map<std::pair<size_t, size_t>, double> weights;
//     std::vector<size_t> used;
//     if (circuitCloseEdgeAdded){
//         graph[end].pop_back(); // the last artificially added edge was the start node
//     }
//     for (size_t i = 0; i < indexToEdges.size(); i++){
//         if (std::find(used.begin(), used.end(), indexToEdges[i].first) == used.end()){
//             used.push_back(indexToEdges[i].first);
//             std::unordered_map<size_t, size_t> trackPaths;
//             std::map<std::pair<size_t, size_t>, double> shortestPaths= dijkstra(graph, indexToEdges[i].first, trackPaths);
//             trackPathsFromEachStart[indexToEdges[i].first] = trackPaths;
//             auto it = shortestPaths.begin();
//             while (it != shortestPaths.end()){
//                 weights[std::pair<size_t, size_t>(it->first)] = shortestPaths[it->first];
//                 it ++;
//             }
//         }
//     }
//     return weights;

// }
// /**********************************************************************************************************/
// void findPerfectMatchingRec(vector<vector<size_t>> &includedEdges, std::vector<size_t> &combinationOfBits, size_t startIndex, std::vector<std::pair<size_t, size_t>> &indexToEdge, std::vector<int> &deltaDegrees, std::map<size_t, int> &leftOptions){
//     auto edge = indexToEdge[startIndex];
//     auto outGoing = edge.first;
//     auto incoming = edge.second;
//     auto degreeIncoming = deltaDegrees[incoming];
//     auto degreeOutGoing = deltaDegrees[outGoing];
//     auto updatedLeftOptions = leftOptions;
//     updatedLeftOptions[outGoing] --;
//     updatedLeftOptions[incoming] --;
//     if (startIndex == indexToEdge.size()-1){
//         std::vector<size_t> newCombination = combinationOfBits;
//         if ((degreeIncoming == 0) && (degreeOutGoing == 0)){
//             newCombination.push_back(0);
//             includedEdges.push_back(newCombination);
//             return;

//         }else{
//             newCombination.push_back(1);
//             includedEdges.push_back(newCombination);
//             return;
//         }

//     }else if ((degreeOutGoing > 0) && (degreeIncoming < 0)){
//         std::vector<size_t> newCombination = combinationOfBits;
//         newCombination.push_back(1);
//         vector<int> updatedDegrees = deltaDegrees;
//         updatedDegrees[outGoing]--;
//         updatedDegrees[incoming] ++;
//         findPerfectMatchingRec(includedEdges, newCombination, startIndex+1, indexToEdge, updatedDegrees, updatedLeftOptions);
//     }
//     // check if it is possible additionally to fix it to zero
//     if ((leftOptions[outGoing] > deltaDegrees[outGoing]) && (leftOptions[incoming] > std::abs(deltaDegrees[incoming]))){
//         // I can use 0, and there will be still enough options in the future
//         std::vector<size_t> newCombination = combinationOfBits;
//         newCombination.push_back(0);
//         vector<int> updatedDegrees = deltaDegrees;
//         findPerfectMatchingRec(includedEdges, newCombination, startIndex+1, indexToEdge, updatedDegrees, updatedLeftOptions);


//     }






// }
// /**********************************************************************************************************/
// double calculateSumOfPath(vector<size_t> &path, std::vector<double> &weightsPerIndices){
//     double sumOfWeights = 0;
//     for (size_t i = 0; i < path.size(); i++){
//         if (path[i] == 1){
//            sumOfWeights +=  weightsPerIndices[i];

//         }
//     }
//     return sumOfWeights;

// }
// /**********************************************************************************************************/
// void countNumberOfInfiniteEdges(std::vector <size_t> &infinitePath, size_t &numOfInfiniteEdges, std::vector<double> &weightsPerIndices){
//     for (size_t i = 0; i < infinitePath.size(); i++){
//         if (infinitePath[i] == 1){
//            if (weightsPerIndices[i] == std::numeric_limits<double>::infinity()){
//                 numOfInfiniteEdges ++;
//            }
//         }
//     }
//     return;
// }
// /**********************************************************************************************************/
// std::vector<pair<size_t, size_t>> findPerfectMatching(std::map<std::pair<size_t, size_t>, size_t> &edgesToIndices, std::vector<std::pair<size_t, size_t>> &indexToEdge, std::vector<double> &weightsPerIndices, std::vector<int> &deltaDegrees, bool &isPahInfinite){
//     vector<vector<size_t>> includedEdges;
//     std::vector <size_t> bestPath;
//     std::vector<pair<size_t, size_t>> selectedEdges;
//     std::map<size_t, int> leftOptions;
//     auto it = edgesToIndices.begin();
//     while(it != edgesToIndices.end()){
//         auto from = it->first.first;
//         auto to = it->first.second;
//         if (leftOptions.find(from) != leftOptions.end()){
//             leftOptions[from] ++;
//         }else{
//             leftOptions[from] = 1;
//         }
//         if (leftOptions.find(to) != leftOptions.end()){
//             leftOptions[to] ++;
//         }else{
//             leftOptions[to] = 1;
//         }

//         it ++;
//     }
//     vector<size_t> combinationOfBits;
//     findPerfectMatchingRec(includedEdges, combinationOfBits, 0, indexToEdge, deltaDegrees, leftOptions);
//     double minSumShortestPath = std::numeric_limits<double>::infinity();
//     for (auto &path : includedEdges){
//         auto pathSum = calculateSumOfPath(path, weightsPerIndices);
//         if (pathSum < minSumShortestPath){
//             minSumShortestPath = pathSum;
//             bestPath = path;
//         }
//     }
//     // handling a case when no path can be found (i.e., disjoint paths)
//     if (minSumShortestPath ==  std::numeric_limits<double>::infinity()){
//         isPahInfinite = true;
//         //try to search for a path with minimal number of edges with inifinity
//         size_t numOfInfiniteEdges = std::numeric_limits<double>::infinity();
//         size_t minNumberOfInfiniteEdges = std::numeric_limits<double>::infinity();
//         for (auto &infinitePath : includedEdges){
//             numOfInfiniteEdges = 0;
//             // here we do the check
//             countNumberOfInfiniteEdges(infinitePath, numOfInfiniteEdges, weightsPerIndices);
//             //
//             if (numOfInfiniteEdges <  minNumberOfInfiniteEdges){
//                 bestPath = infinitePath;
//             }
//         }
//     }

//     // finding the best path
//     for (size_t i = 0; i < bestPath.size(); i++){
//         if (bestPath[i] == 0){
//             continue;
//         }
//         selectedEdges.push_back(indexToEdge[i]);
//     }
//     return selectedEdges;

// }


// /*********************************************************************************************************/
// void constructEulerRec(vector<pair<size_t, size_t>> &finalPath, size_t j, vector<size_t> &usedCircuitsIndices, std::vector<vector<pair<size_t,size_t>>> &circuits){
//     usedCircuitsIndices.push_back(j);
//     for (size_t i = 0; i < circuits[j].size();i++){
//         finalPath.push_back(circuits[j][i]);
//         auto outgoingVertex = circuits[j][i].second;
//         for (size_t k = 1; k < circuits.size(); k++){
//             if (std::find(usedCircuitsIndices.begin(), usedCircuitsIndices.end(), k) != usedCircuitsIndices.end()){
//                 continue;
//             }
//             if (outgoingVertex == circuits[k][0].first){
//                 constructEulerRec(finalPath, k, usedCircuitsIndices, circuits);
//             }
//         }
//     }

// }
// /*********************************************************************************************************/
// vector<pair<size_t, size_t>> constructFinalEulerPathFromCircuits(std::vector<vector<pair<size_t,size_t>>> &circuits, size_t start, size_t end){
//     auto firstCurcuit = circuits[0];
//     vector<pair<size_t, size_t>> finalPath;
//     vector<size_t> usedCircuitsIndices;
//     usedCircuitsIndices.push_back(0);
//     for (size_t i = 0; i < firstCurcuit.size();i++){
//         finalPath.push_back(firstCurcuit[i]);
//         auto outgoingVertex = firstCurcuit[i].second;
//         for (size_t j = 1; j < circuits.size(); j++){
//             if (outgoingVertex == circuits[j][0].first){
//                 constructEulerRec(finalPath, j, usedCircuitsIndices, circuits);
//             }
//         }
//     }
//     return finalPath;

// }
// /**********************************************************************************************************/
// vector<std::pair<size_t, size_t>> findEulerCircuit(vector<vector<size_t>>& graph, size_t start, size_t end){
//     size_t graphSize = graph.size();
//     vector<size_t> unusedEdges;
//     unusedEdges.resize(graphSize);
//     size_t numOfEdges = getNumberOfEdges(graph);
//     for (size_t i = 0; i < graphSize; i++){
//         unusedEdges[i] = graph[i].size(); // for each vertex we have its respective number of edges
//     }
//     // some edges might present twice
//     std::map<pair<size_t,size_t>,size_t> edgesOccurrences;
//     for (size_t i = 0; i < graph.size(); i++){
//         for (size_t j = 0; j < graph[i].size(); j++){
//             pair<size_t, size_t> graphEdge(i, graph[i][j]);
//             if (edgesOccurrences.find(graphEdge) != edgesOccurrences.end()){
//                 edgesOccurrences[graphEdge] += 1;
//             }else{
//                 edgesOccurrences[graphEdge] = 1;
//             }
//         }
//     }

//     size_t startNode = start;
//     std::vector<vector<pair<size_t, size_t>>> curcuits;
//     vector<size_t> usedVertices;
//     while (numOfEdges > 0){
//         auto foundCircuit = findSubCircuit(graph, startNode, unusedEdges, numOfEdges, edgesOccurrences);
//         curcuits.push_back(foundCircuit);
//         for (size_t j = 0; j < graphSize; j++){
//             if (unusedEdges[j] > 0){ // not all the edges of this node are covered
//                 startNode = j;
//                 break;
//             }
//         }
        
//         //firstVector.insert(firstVector.begin() + 2, secondVector.begin(), secondVector.end());
//     }
//     vector<std::pair<size_t, size_t>> concatenatedPath = constructFinalEulerPathFromCircuits(curcuits, start, end);
//     return concatenatedPath;


// }
// /**********************************************************************************************************/
// std::vector <std::pair<size_t, size_t>> getEulerPath(vector<vector<size_t>>& graph, size_t start, size_t end){
//     std::vector <std::pair<size_t, size_t>> eulerPath;
//     bool hasEulerPath = false;
//     bool hasEulerCircuit = false;
//     if (start == end){
//         hasEulerCircuit = isEulerCurcuit(graph);

//     }else{
//         hasEulerPath = isEulerPath(graph, start, end);

//     }
//     if (hasEulerPath || hasEulerCircuit){
//         if (hasEulerPath){
//           graph[end].push_back(start);  
//         }
//         eulerPath = findEulerCircuit(graph, start, end);
//         if (hasEulerPath){
//             graph[end].pop_back();
//             std::pair<size_t, size_t> artificialEdge(end, start);
//             auto it = std::find(eulerPath.begin(), eulerPath.end(), artificialEdge);
//             size_t index = std::distance(eulerPath.begin(), it);
//             size_t startPath = index+1;
//             vector<pair<size_t, size_t>> eulerPathCorrected;
//             // two scenarios:
//             // v-->uv
//             // v-->uv-->v
//             if (startPath < eulerPath.size()){
//                 for (size_t i = startPath; i < eulerPath.size(); i++){
//                     eulerPathCorrected.push_back(eulerPath[i]);
//                 }

//             }
//             for (size_t i = 0; i < index; i++){
//                 eulerPathCorrected.push_back(eulerPath[i]);
//             }
//             eulerPath = eulerPathCorrected;


//         }

        
//     }
    
//     return eulerPath;  

// } 
// /**********************************************************************************************************/
// // bit for each edge. For each vertex we have its degree, which gets updated.
// vector<vector<pair<size_t, size_t>>> chinesePostman(std::vector<vector<size_t>> &graph, size_t start, size_t end){
//     vector<pair<size_t,size_t>> bestPath;
//     bool circuitCloseEdgeAdded = false;
//     if (std::find(graph[end].begin(), graph[end].end(), start) == graph[end].end()){
//         graph[end].push_back(start);
//         circuitCloseEdgeAdded = true;

//     }
//     std::vector<std::pair<size_t, size_t>> inAndOut = getDegreesForEachVertex(graph);

//     // now find the delta degrees. Note that in each pair, first is the in-degree, while second is the out degree.
//     std::vector<int> deltaDegrees = getDeltaDegrees(inAndOut);
//     std::vector<size_t> inVertices; // those with extra incoming edges (should be first in the bipartite graph)
//     std::vector<size_t> outVertices; // those with extra outcoming edges (should be second in the bipartite graph)
//     fillInAndOutOddVertices(deltaDegrees, inVertices, outVertices);
//     // match between indices and edges in the bipartite graph
//     std::vector<std::pair<size_t, size_t>> indexToEdge;
//     for (size_t i = 0; i < inVertices.size(); i++){
//         for (size_t j = 0; j < outVertices.size(); j++){

//             indexToEdge.push_back(std::pair<size_t, size_t>(inVertices[i], outVertices[j]));
//         }
//     }
//     // for each edge find its edge (u, v), which is the shortest path between u and v.
//     std::unordered_map<size_t, unordered_map<size_t, size_t>> trackPathsFromEachStart;
//     std::map<std::pair<size_t, size_t>, double> weights = findWeightsForOddVerticesEdges(graph, start, end, indexToEdge, circuitCloseEdgeAdded,trackPathsFromEachStart);

//     std::map<std::pair<size_t, size_t>, size_t> edgesToIndices;
//     std::vector<double> weightsPerIndices;
//     for (size_t i = 0; i < indexToEdge.size(); i++){
//         weightsPerIndices.push_back(weights[indexToEdge[i]]);
//         edgesToIndices[indexToEdge[i]] = i;
//     }
//     // now find the edges between the odd vertices
//     bool infinitePath = false;
//     bestPath  = findPerfectMatching(edgesToIndices, indexToEdge, weightsPerIndices, deltaDegrees, infinitePath);
//     std::vector<pair<size_t, size_t>> infiniteEdges;
//     // if there are some infinite edges - find them
//     if (infinitePath){
//         for (size_t i = 0; i < bestPath.size(); i++){
//             if (weights[bestPath[i]] == std::numeric_limits<double>::infinity()){
//                 infiniteEdges.push_back(bestPath[i]);
//             }
//         }
//     }
//     std::map<pair<size_t,size_t>, std::vector<pair<size_t, size_t>>> edgeWithPath;
//     for (auto &edge : bestPath){
//         if ((infinitePath) && (std::find(infiniteEdges.begin(), infiniteEdges.end(), edge) != infiniteEdges.end())){
//             continue;
//         }
//         if (std::find(graph[edge.first].begin(), graph[edge.first].end(), edge.second) == graph[edge.first].end()){
//             // this edge does not present in the graph, so it represents a path that should be restored from dijkstra
//             auto source = edge.first;
//             auto dst = edge.second;
//             size_t currentNode = dst;
//             auto &trackPath = trackPathsFromEachStart[source];
//             edgeWithPath[edge];
//             while (currentNode != source){
//                 std::pair<size_t,size_t> edgeInPath(trackPath[currentNode], currentNode);
//                 edgeWithPath[edge].push_back(edgeInPath); //putting it in a reverse order
//                 currentNode = trackPath[currentNode];
                
//             }
//         }
//     }
//     if (infinitePath){
//         // we apply dijkstra on a modified graph, that includes the <dst, source> edge
//         std::unordered_map<size_t, unordered_map<size_t, size_t>> trackPathsUpdated;
//         std::map<std::pair<size_t, size_t>, double> weightsUpdated = findWeightsForOddVerticesEdges(graph, start, end, indexToEdge, false, trackPathsUpdated);
//         for (auto &infiniteEdge : infiniteEdges){
//             auto source = infiniteEdge.first;
//             auto dst = infiniteEdge.second;
//             size_t currentNode = dst;
//             auto &trackInfinitePath = trackPathsUpdated[source];
//             edgeWithPath[infiniteEdge];
//             while (currentNode != source){
//                 std::pair<size_t,size_t> edgeInPath(trackInfinitePath[currentNode], currentNode);
//                 edgeWithPath[infiniteEdge].push_back(edgeInPath); //putting it in a reverse order
//                 currentNode = trackInfinitePath[currentNode];
                
//             }

//         }
//     } 
//     if (circuitCloseEdgeAdded){
//         graph[end].pop_back(); // we do it because this edge is added in the eiuler tour function
//     }

//     // here we add the additional edges that were constructed for the odd vertices
//     for (size_t i = 0; i < bestPath.size(); i++){
//         graph[bestPath[i].first].push_back(bestPath[i].second);
//     }

//     // now we can apply the euler path algorithm on the modified graph with extra edges.
//     std::vector <std::pair<size_t, size_t>> euler = getEulerPath(graph, start, end);
//     std::vector<pair<size_t, size_t>> chinese;
//     for (size_t i = 0; i < euler.size(); i++){
//         if (edgeWithPath.find(euler[i]) == edgeWithPath.end()){
//             // edge exists in the graph
//             chinese.push_back(euler[i]);
//         }else{ //its a path
//             auto &repeatedEdges = edgeWithPath[euler[i]];
//             for (int j = repeatedEdges.size()-1; j >=0; j--){
//                 chinese.push_back(repeatedEdges[j]);
//             }

//         }

//     }
//     // if there are disjoint paths, i.e., infinitePath is true, we need to find all the subpaths
//     vector<vector<pair<size_t,size_t>>> disjointPaths;
//     if (infinitePath){
//         std::pair<size_t,size_t> artificialEdge(end, start);   
//         vector<pair<size_t,size_t>> subPath;
//         for (size_t i = 0; i < chinese.size();i++){
//             if (chinese[i] != artificialEdge){
//                 subPath.push_back(chinese[i]);

//             }else{
//                 disjointPaths.push_back(subPath);
//                 subPath.clear();
//             }
//             if(i == chinese.size()-1){
//                 disjointPaths.push_back(subPath);
//             }


//         }
//     }else{
//         // there are no disjoint paths, so we return the path results from the euler function without any modifications
//         disjointPaths.push_back(chinese);
//     }
//     return disjointPaths;
    

// }
/**********************************************************************************************************/
void printGraph(vector<vector<size_t>>& graph){
    std::cout << "Printing Graph!" << std::endl;
    for (size_t i = 0; i < graph.size(); i++){
        std::cout << "Current Vertex is " << i << std::endl;
        for (size_t j = 0; j < graph[i].size(); j++){
            std::cout << "\t" << graph[i][j] << std::endl;
        }
    }
}




/**********************************************************************************************************/
void testEulerPath(){
    vector<vector<size_t>> G;
    G.resize(5);
    G[0].push_back(1);
    G[0].push_back(2);
    G[0].push_back(3);

    G[1].push_back(2);

    G[2].push_back(0);
    G[2].push_back(3);

    G[3].push_back(0);
    G[3].push_back(4);

    G[4].push_back(0);
    bool foundPath = true;

    auto path = MultiStateMappingPath::getEulerPath(G, 0, 0);
    std::cout << "Euler path is:" << std::endl;
    for (size_t i = 0; i < path.size(); i++){
        std::cout << path[i].first << "," << path[i].second << "\t";
    }



}
/**********************************************************************************************************/
void testChinese2(){
    vector<vector<size_t>> G;
    G.resize(6);
    G[0].push_back(1);
    G[0].push_back(3);

    G[1].push_back(3);
    G[1].push_back(2);




    G[3].push_back(2);
    G[3].push_back(4);
    G[3].push_back(5);

    G[4].push_back(2);
    G[5].push_back(0);
    bool foundPath = true;
    auto paths = MultiStateMappingPath::chinesePostman(G, 0, 2, foundPath);
    std::cout << "chinese path is:" << std::endl;
    for (auto &path : paths){
        std::cout << "\tsubpath is:" << std::endl;
        for (size_t i = 0; i < path.size(); i++){
            std::cout << path[i].first << "," << path[i].second << "\t";
        }

    }

    std::cout << endl;


}
/**********************************************************************************************************/
void testChinese3(){
    vector<vector<size_t>> G;
    G.resize(4);
    G[1].push_back(2);
    G[1].push_back(3);
    G[2].push_back(0);
    G[3].push_back(0);
    bool foundPath = true;


    auto paths = MultiStateMappingPath::chinesePostman(G, 1, 0, foundPath);
    std::cout << "chinese path is:" << std::endl;
    for (auto &path : paths){
        std::cout << "\tsubpath is:" << std::endl;
        for (size_t i = 0; i < path.size(); i++){
            std::cout << path[i].first << "," << path[i].second << "\t";
        }

    }

    std::cout << endl;


}
/**********************************************************************************************************/
void testTSP3(){
  std::cout << "test TSP 3:" <<std::endl;
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
  bool foundPath = true;
  auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime, foundPath);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }
}

void testTSP4(){
  std::cout << "test TSP 4:" << std::endl;
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
  bool foundPath = true;
  auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime, foundPath);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }
}

void testTSP1(){
  std::cout << "test TSP 1:" << std::endl;
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
  bool foundPath = true;
  auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime, foundPath);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }

}


void testTSP2(){
  std::cout << "test TSP 2:" << std::endl;
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
  bool foundPath = true;
  auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(start, end, transitions, dwellingTimes, totalDurationTime, foundPath);
  for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
  }

}
/**********************************************************************************************************/
void testTSP5(){
    std::cout << "test TSP 5: finding TSP and not CPP" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0,1}, {1,3}, {2,1}, {2,4}, {3,2}, {3,4},{4,1}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(5);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    dwellingTimes[3] = 0.25;
    dwellingTimes[4] = 0.25;
    double totalDurationTime = 1.1;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 4, transitions, dwellingTimes, totalDurationTime, foundPath);
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }

    // G.resize(5);
    // G[0].push_back(1);

    // G[1].push_back(3);

    // G[2].push_back(1);
    // G[2].push_back(4);

    // G[3].push_back(2);
    // G[3].push_back(4);

    // G[4].push_back(1);


}
/**********************************************************************************************************/
void testCPP1(){
    std::cout << "test CPP 1: Euler" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0,1}, {1,0}, {0, 2}, {2,0}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(3);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    double totalDurationTime = 0.6;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 0, transitions, dwellingTimes, totalDurationTime, foundPath);
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }

}
/**********************************************************************************************************/
void testCPP2(){
    std::cout << "test CPP 2: chinese" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0,3}, {3,0}, {0, 2}, {2,3}, {0, 1}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(4);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    dwellingTimes[3] = 0.4;
    double totalDurationTime = 1.0;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 1, transitions, dwellingTimes, totalDurationTime, foundPath);
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }

}
/**********************************************************************************************************/
void testCPP3(){
    std::cout << "test CPP 3: chinese indexing" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0,4}, {4,0}, {0, 2}, {2,4}, {0, 1}, {0, 3}, {3,0}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(5);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    dwellingTimes[3] = 0.0000001;
    dwellingTimes[4] = 0.1;
    double totalDurationTime = 0.7000001;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 1, transitions, dwellingTimes, totalDurationTime, foundPath);
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }

}
/****************************************************************************************/
void testCPP4(){
    std::cout << "test CPP 4: chinese disjoint" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0,1}, {0,3}, {1,3}, {1, 2}, {3, 2}, {3,4}, {3, 5}, {4, 2}, {5,0}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(6);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    dwellingTimes[3] = 0.1;
    dwellingTimes[4] = 0.1;
    dwellingTimes[5] = 0.1;
    double totalDurationTime = 0.9;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 2, transitions, dwellingTimes, totalDurationTime, foundPath);
    if (!foundPath){
      std::cout << "path not found!!" << std::endl;
    }
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }

}
void testCPP5(){
    std::cout << "test CPP 5: chinese" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0, 1}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(3);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    double totalDurationTime = 0.6;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 2, transitions, dwellingTimes, totalDurationTime, foundPath);
    if (!foundPath){
      std::cout << "not found path!" << std::endl;
    }
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }
}

void testCPP6(){
    std::cout << "test CPP 6: chinese" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0, 1}, {1,0}, {0,2}, {2,0}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(3);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    double totalDurationTime = 0.6;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 2, transitions, dwellingTimes, totalDurationTime, foundPath);
    if (!foundPath){
      std::cout << "not found path!" << std::endl;
    }
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }
}
void testCPP7(){
    std::cout << "test CPP 7: chinese" << std::endl;
    std::vector<std::pair<size_t, size_t>> edges = {{0, 1}, {1,3}, {0,2}, {2,3}};
    std::map<std::pair<size_t, size_t>, double> transitions;
    for (auto &edge : edges){
        transitions[edge] = 1.0;
    }
    vector<double> dwellingTimes;
    dwellingTimes.resize(4);
    dwellingTimes[0] = 0.3;
    dwellingTimes[1] = 0.2;
    dwellingTimes[2] = 0.1;
    dwellingTimes[3] = 0.1;
    double totalDurationTime = 0.7;
    bool foundPath = true;
    auto path = MultiStateMappingPath::findExpectedMappingPathForEachNode(0, 3, transitions, dwellingTimes, totalDurationTime, foundPath);
    if (!foundPath){
      std::cout << "not found path!" << std::endl;
    }
    for (size_t i = 0; i < path.size(); i++){
      std::cout << path[i] << std::endl;
    }
}
void testChinese1(){
    vector<vector<size_t>> G;
    G.resize(5);
    G[0].push_back(1);

    G[1].push_back(3);

    G[2].push_back(1);
    G[2].push_back(4);

    G[3].push_back(2);
    G[3].push_back(4);

    G[4].push_back(1);
    bool foundPath = true;
    auto paths = MultiStateMappingPath::chinesePostman(G, 0, 4, foundPath);
    std::cout << "chinese path is:" << std::endl;
    for (auto &path : paths){
        std::cout << "\tsubpath is:" << std::endl;
        for (size_t i = 0; i < path.size(); i++){
            std::cout << path[i].first << "," << path[i].second << "\t";
        }

    }

    std::cout << endl;


}

/**********************************************************************************************************/
// vector<vector<vector<size_t>>> findDisJointPaths(vector<vector<size_t>> &graph, size_t start, size_t end){
//     // if there are node x and y that are not reachable from each other they are not together in the path
//     // if at least one of them is reachable from another, it means that there is a path which contains them both
    
//     vector<vector<size_t>> nodesComeTogether;
//     std::map<std::pair<size_t, size_t>, bool> pathExists;
//     for (size_t i = 0; i < graph.size(); i++){
//         std::unordered_map<size_t, size_t> trackPaths;
//         std::map<std::pair<size_t, size_t>, double> shortestPaths= dijkstra(graph, i, trackPaths);
//         auto it = shortestPaths.begin();
//         while(it != shortestPaths.end()){
//             if (shortestPaths[it->first] < std::numeric_limits<double>::infinity()){
//                 pathExists[it->first] = true;
//             }
//             it++;
//         }
//     }
//     std::vector<bool> strangerNodes(graph.size(), false);
//     for (size_t i = 0; i < graph.size()-1; i++){
//         for (size_t j = i+1; j < graph.size(); j++){
//             if ((i == start) || (i == end)){
//                 continue;

//             }
//             std::pair<size_t, size_t> path(i, j);
//             std::pair<size_t, size_t> oppositePath(j, i);

//             if ((!pathExists[path]) && (!pathExists[oppositePath])){
//                 strangerNodes[i] = true;
//                 strangerNodes[j] = true;

//             }
//         }
//     }
//     vector<bool> used(graph.size(), false);
//     vector<vector<size_t>> nodesInDisjointPaths;
//     for (size_t i = 0; i < strangerNodes.size(); i++){
//         if (strangerNodes[i]){
//             if (!used[i]){
//                 std::vector<size_t> possiblePath;
//                 possiblePath.push_back(i);
//                 for (size_t j = 0; j < graph.size(); j++){
//                     if (j == i){
//                         continue;
//                     }
//                     if (!strangerNodes[j]){
//                         possiblePath.push_back(j);
//                     }else{
//                         std::pair<size_t, size_t> path(i, j);
//                         std::pair<size_t, size_t> oppositePath(j, i);
//                         if ((pathExists[path]) || (pathExists[oppositePath])){
//                             possiblePath.push_back(j)


//                         }

//                     }

//                 }

//                 used[i] = true;

//             }
//         }

//     }
// }
/**********************************************************************************************************/
void testDijkstra(){
    vector<vector<size_t>> G;
    G.resize(9);
    G[0].push_back(1);
    G[0].push_back(3);
    G[1].push_back(2);
    G[1].push_back(4);
    G[2].push_back(3);
    G[2].push_back(5);
    G[3].push_back(4);
    G[4].push_back(8);
    G[5].push_back(6);
    G[6].push_back(7);
    G[6].push_back(0);
    G[7].push_back(0);
    G[8].push_back(7);
    G[8].push_back(4);
    std::unordered_map<size_t, size_t> trackPaths;
    std::map<std::pair<size_t, size_t>, double> shortestPaths= MultiStateMappingPath::dijkstra(G, 0, trackPaths);
    auto it = shortestPaths.begin();
    std::cout << "Dijkstra: " << std::endl;
    while(it != shortestPaths.end()){
        std::cout << (it->first).first << "-->" << (it->first).second << " = " << it->second << std::endl;
        it ++;
    }

}
void testPerfectMatching(){
    std::vector<pair<size_t, size_t>> indexToEdge = {{0,4}, {0,5}, {1,4}, {1,5}, {2,4}, {2,5}};
    std::map<std::pair<size_t, size_t>, size_t> edgesToIndices;
    std::vector<double> weightsPerIndices = {1.0, 3.0, 2.0, 2.0, 2.0, 1.0};
    for (size_t i = 0; i < indexToEdge.size(); i++){
        edgesToIndices[indexToEdge[i]] = i;
    }
    vector<int> deltaDegrees = {1, 1, 1, 0, -2, -1};
    bool infinitePath = false;

    auto bestPath  = MultiStateMappingPath::findPerfectMatching(edgesToIndices, indexToEdge, weightsPerIndices, deltaDegrees, infinitePath);
    std::cout << endl;
    std::cout << "Perfect Matching is:" << std::endl;
    for (size_t i = 0; i < bestPath.size(); i++){   
        std::cout << bestPath[i].first << "," << bestPath[i].second << "\t";
        
    }
    std::cout << endl;
}
int main(){
    //test_dfs();
    testTSP1();
    testTSP2();
    testTSP3();
    testTSP4();
    testTSP5();
    testCPP1();
    testCPP2();
    testCPP3();
    testCPP4();
    testCPP5();
    testCPP6();
    testCPP7();

    return 0;
}