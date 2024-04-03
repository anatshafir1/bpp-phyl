#include "MultiStateMappingPath.h"
using namespace bpp;


std::vector<std::pair<size_t, size_t>> MultiStateMappingPath::getDegreesForEachVertex(vector<vector<size_t>>& graph){
    size_t graphSize = graph.size();
    std::vector<std::pair<size_t, size_t>> degrees(graphSize, {0, 0});
    for (size_t i = 0; i < graphSize; i++){
        for (size_t j = 0; j < graph[i].size(); j++){
            degrees[i].second += 1; // out degree
            degrees[graph[i][j]].first += 1; // in degree
        }
    }
    return degrees;
}
/**********************************************************************************************************/
bool MultiStateMappingPath::isEulerCurcuit(vector<vector<size_t>>& graph){
    std::vector<std::pair<size_t, size_t>> degrees = getDegreesForEachVertex(graph);
    for (size_t i = 0; i < degrees.size(); i++){
        if (degrees[i].first != degrees[i].second){
            return false;
        }
    }
    return true;
}
/**********************************************************************************************************/
bool MultiStateMappingPath::isEulerPath(vector<vector<size_t>>& graph, size_t start, size_t end){
    std::vector<std::pair<size_t, size_t>> degrees = getDegreesForEachVertex(graph);
    for (size_t i = 0; i < degrees.size(); i++){
        if (start == i){
            if (degrees[i].second != degrees[i].first + 1){
                return false;

            }
        }else if (end == i){
            if (degrees[i].first != degrees[i].second + 1){
                return false;

            }
        }else{
            if (degrees[i].first != degrees[i].second){
                return false;
            }   
        }
    }
    return true;
}
/**********************************************************************************************************/
size_t MultiStateMappingPath::getNumberOfEdges(vector<vector<size_t>>& graph){
    size_t numOfEdges = 0;
    size_t graphSize = graph.size();
    for (size_t i = 0; i < graphSize; i++){
        numOfEdges += graph[i].size();
    }
    return numOfEdges;

}
/**********************************************************************************************************/
std::vector<pair<size_t, size_t>> MultiStateMappingPath::findSubCircuit(vector<vector<size_t>>& graph, size_t startNode, vector<size_t> &unusedEdges, size_t &numOfEdges, std::map<pair<size_t,size_t>,size_t> &edgesOccurrences){
    size_t currentNode = startNode; 
    size_t nextNode = graph[startNode][0];
    vector<pair<size_t, size_t>> circuit;    
    while (nextNode != startNode){
        auto neighbors = graph[currentNode];
        for (size_t i = 0; i < neighbors.size(); i++){
            std::pair<size_t, size_t> edge(currentNode, neighbors[i]);
            if (edgesOccurrences[edge] == 0){
                continue;
            }
            unusedEdges[currentNode]--;
            numOfEdges--;
            circuit.push_back(std::pair<size_t,size_t>(currentNode, neighbors[i]));
            edgesOccurrences[edge] --;
            currentNode = neighbors[i];
            
            break; // we choose only one of the neighbors to continue with
            

        }
        nextNode = currentNode;
    }
    return circuit;

}
/**********************************************************************************************************/
std::vector<int> MultiStateMappingPath::getDeltaDegrees(std::vector<std::pair<size_t, size_t>> &inAndOutDeg){
    std::vector<int> deltaDegrees;
    deltaDegrees.resize(inAndOutDeg.size());
    size_t inDeg;
    size_t outDeg;
    for (size_t i = 0; i < inAndOutDeg.size(); i++){
        inDeg = inAndOutDeg[i].first;
        outDeg = inAndOutDeg[i].second;
        deltaDegrees[i] = static_cast<int>(inDeg) - static_cast<int>(outDeg);

    }
    return deltaDegrees;

}
/**********************************************************************************************************/
void MultiStateMappingPath::fillInAndOutOddVertices(std::vector<int>& deltaDegrees, std::vector<size_t> &inVertices, std::vector<size_t> &outVertices){
    for (size_t i = 0; i < deltaDegrees.size(); i++){
        if (deltaDegrees[i] < 0){
            outVertices.push_back(i);
        }else if (deltaDegrees[i] > 0){
            inVertices.push_back(i);
        }
    }
    return;
}

/**********************************************************************************************************/
std::map<std::pair<size_t, size_t>, double> MultiStateMappingPath::dijkstra(std::vector<vector<size_t>> &graph, size_t start, std::unordered_map<size_t, size_t> &track){
    vector<pair<size_t, double>> nodesWithdistances;
    nodesWithdistances.resize(graph.size());
    for (size_t i = 0; i < nodesWithdistances.size(); i++){
        nodesWithdistances[i].first = i;
        if (i == start){
            nodesWithdistances[i].second = 0;
        }else{
            nodesWithdistances[i].second = std::numeric_limits<double>::infinity();
        }
    }

    vector<double> distances(graph.size(), std::numeric_limits<double>::infinity());
    distances[start] = 0;
    std::map<std::pair<size_t, size_t>, double> shortestPaths;
    std::priority_queue<std::pair<size_t, double>, std::vector<std::pair<size_t, double>>, std::greater<std::pair<size_t, double>>> minHeap;
    minHeap.push({start, 0});
    
    while (!minHeap.empty()) {
        size_t node = minHeap.top().first;
        double nodeDist = minHeap.top().second;
        minHeap.pop();
        if (nodeDist > distances[node]) {
            continue;
        }
        for (auto &neighbor : graph[node]) {
            if (distances[node] + 1 < distances[neighbor]) {
                distances[neighbor] = distances[node] + 1;
                minHeap.push({neighbor, distances[neighbor]});
                track[neighbor] = node;
            }
        }

    }

    for (size_t i = 0; i < distances.size(); i++){
        shortestPaths[std::pair<size_t, size_t>(start, i)] = distances[i];
    }
    return shortestPaths;

}
/**********************************************************************************************************/
// the parameter graph is by value on purpose, because I'm going to change the object here
std::map<std::pair<size_t, size_t>, double> MultiStateMappingPath::findWeightsForOddVerticesEdges(std::vector<vector<size_t>> graph, size_t start, size_t end, std::vector<std::pair<size_t, size_t>> &indexToEdges, bool circuitCloseEdgeAdded, std::unordered_map<size_t, unordered_map<size_t, size_t>> &trackPathsFromEachStart){
    std::map<std::pair<size_t, size_t>, double> weights;
    std::vector<size_t> used;
    if (circuitCloseEdgeAdded){
        graph[end].pop_back(); // the last artificially added edge was the start node
    }
    for (size_t i = 0; i < indexToEdges.size(); i++){
        if (std::find(used.begin(), used.end(), indexToEdges[i].first) == used.end()){
            used.push_back(indexToEdges[i].first);
            std::unordered_map<size_t, size_t> trackPaths;
            std::map<std::pair<size_t, size_t>, double> shortestPaths= dijkstra(graph, indexToEdges[i].first, trackPaths);
            trackPathsFromEachStart[indexToEdges[i].first] = trackPaths;
            auto it = shortestPaths.begin();
            while (it != shortestPaths.end()){
                weights[std::pair<size_t, size_t>(it->first)] = shortestPaths[it->first];
                it ++;
            }
        }
    }
    return weights;

}
/**********************************************************************************************************/
void MultiStateMappingPath::findPerfectMatchingRec(vector<vector<size_t>> &includedEdges, std::vector<size_t> &combinationOfBits, size_t startIndex, std::vector<std::pair<size_t, size_t>> &indexToEdge, std::vector<int> &deltaDegrees, std::map<size_t, int> &leftOptions){
    auto edge = indexToEdge[startIndex];
    auto outGoing = edge.first;
    auto incoming = edge.second;
    auto degreeIncoming = deltaDegrees[incoming];
    auto degreeOutGoing = deltaDegrees[outGoing];
    auto updatedLeftOptions = leftOptions;
    updatedLeftOptions[outGoing] --;
    updatedLeftOptions[incoming] --;
    if (startIndex == indexToEdge.size()-1){
        std::vector<size_t> newCombination = combinationOfBits;
        if ((degreeIncoming == 0) && (degreeOutGoing == 0)){
            newCombination.push_back(0);
            includedEdges.push_back(newCombination);
            return;

        }else{
            newCombination.push_back(1);
            includedEdges.push_back(newCombination);
            return;
        }

    }else if ((degreeOutGoing > 0) && (degreeIncoming < 0)){
        std::vector<size_t> newCombination = combinationOfBits;
        newCombination.push_back(1);
        vector<int> updatedDegrees = deltaDegrees;
        updatedDegrees[outGoing]--;
        updatedDegrees[incoming] ++;
        findPerfectMatchingRec(includedEdges, newCombination, startIndex+1, indexToEdge, updatedDegrees, updatedLeftOptions);
    }
    // check if it is possible additionally to fix it to zero
    if ((leftOptions[outGoing] > deltaDegrees[outGoing]) && (leftOptions[incoming] > std::abs(deltaDegrees[incoming]))){
        // I can use 0, and there will be still enough options in the future
        std::vector<size_t> newCombination = combinationOfBits;
        newCombination.push_back(0);
        vector<int> updatedDegrees = deltaDegrees;
        findPerfectMatchingRec(includedEdges, newCombination, startIndex+1, indexToEdge, updatedDegrees, updatedLeftOptions);


    }






}
/**********************************************************************************************************/
double MultiStateMappingPath::calculateSumOfPath(vector<size_t> &path, std::vector<double> &weightsPerIndices){
    double sumOfWeights = 0;
    for (size_t i = 0; i < path.size(); i++){
        if (path[i] == 1){
           sumOfWeights +=  weightsPerIndices[i];

        }
    }
    return sumOfWeights;

}
/**********************************************************************************************************/
void MultiStateMappingPath::countNumberOfInfiniteEdges(std::vector <size_t> &infinitePath, size_t &numOfInfiniteEdges, std::vector<double> &weightsPerIndices){
    for (size_t i = 0; i < infinitePath.size(); i++){
        if (infinitePath[i] == 1){
           if (weightsPerIndices[i] == std::numeric_limits<double>::infinity()){
                numOfInfiniteEdges ++;
           }
        }
    }
    return;
}
/**********************************************************************************************************/
std::vector<pair<size_t, size_t>> MultiStateMappingPath::findPerfectMatching(std::map<std::pair<size_t, size_t>, size_t> &edgesToIndices, std::vector<std::pair<size_t, size_t>> &indexToEdge, std::vector<double> &weightsPerIndices, std::vector<int> &deltaDegrees, bool &isPahInfinite){
    vector<vector<size_t>> includedEdges;
    std::vector <size_t> bestPath;
    std::vector<pair<size_t, size_t>> selectedEdges;
    std::map<size_t, int> leftOptions;
    auto it = edgesToIndices.begin();
    while(it != edgesToIndices.end()){
        auto from = it->first.first;
        auto to = it->first.second;
        if (leftOptions.find(from) != leftOptions.end()){
            leftOptions[from] ++;
        }else{
            leftOptions[from] = 1;
        }
        if (leftOptions.find(to) != leftOptions.end()){
            leftOptions[to] ++;
        }else{
            leftOptions[to] = 1;
        }

        it ++;
    }
    vector<size_t> combinationOfBits;
    findPerfectMatchingRec(includedEdges, combinationOfBits, 0, indexToEdge, deltaDegrees, leftOptions);
    double minSumShortestPath = std::numeric_limits<double>::infinity();
    for (auto &path : includedEdges){
        auto pathSum = calculateSumOfPath(path, weightsPerIndices);
        if (pathSum < minSumShortestPath){
            minSumShortestPath = pathSum;
            bestPath = path;
        }
    }
    // handling a case when no path can be found (i.e., disjoint paths)
    if (minSumShortestPath ==  std::numeric_limits<double>::infinity()){
        isPahInfinite = true;
        //try to search for a path with minimal number of edges with inifinity
        size_t numOfInfiniteEdges = std::numeric_limits<double>::infinity();
        size_t minNumberOfInfiniteEdges = std::numeric_limits<double>::infinity();
        for (auto &infinitePath : includedEdges){
            numOfInfiniteEdges = 0;
            // here we do the check
            countNumberOfInfiniteEdges(infinitePath, numOfInfiniteEdges, weightsPerIndices);
            //
            if (numOfInfiniteEdges <  minNumberOfInfiniteEdges){
                bestPath = infinitePath;
            }
        }
    }

    // finding the best path
    for (size_t i = 0; i < bestPath.size(); i++){
        if (bestPath[i] == 0){
            continue;
        }
        selectedEdges.push_back(indexToEdge[i]);
    }
    return selectedEdges;

}


/*********************************************************************************************************/
void MultiStateMappingPath::constructEulerRec(vector<pair<size_t, size_t>> &finalPath, size_t j, vector<size_t> &usedCircuitsIndices, std::vector<vector<pair<size_t,size_t>>> &circuits){
    usedCircuitsIndices.push_back(j);
    for (size_t i = 0; i < circuits[j].size();i++){
        finalPath.push_back(circuits[j][i]);
        auto outgoingVertex = circuits[j][i].second;
        for (size_t k = 1; k < circuits.size(); k++){
            if (std::find(usedCircuitsIndices.begin(), usedCircuitsIndices.end(), k) != usedCircuitsIndices.end()){
                continue;
            }
            if (outgoingVertex == circuits[k][0].first){
                constructEulerRec(finalPath, k, usedCircuitsIndices, circuits);
            }
        }
    }

}
/*********************************************************************************************************/
vector<pair<size_t, size_t>> MultiStateMappingPath::constructFinalEulerPathFromCircuits(std::vector<vector<pair<size_t,size_t>>> &circuits, size_t start, size_t end){
    auto firstCurcuit = circuits[0];
    vector<pair<size_t, size_t>> finalPath;
    vector<size_t> usedCircuitsIndices;
    usedCircuitsIndices.push_back(0);
    for (size_t i = 0; i < firstCurcuit.size();i++){
        finalPath.push_back(firstCurcuit[i]);
        auto outgoingVertex = firstCurcuit[i].second;
        for (size_t j = 1; j < circuits.size(); j++){
            if (outgoingVertex == circuits[j][0].first){
                if (std::find(usedCircuitsIndices.begin(), usedCircuitsIndices.end(), j) != usedCircuitsIndices.end()){
                    continue;
                }
                constructEulerRec(finalPath, j, usedCircuitsIndices, circuits);
            }
        }
    }
    return finalPath;

}
/**********************************************************************************************************/
vector<std::pair<size_t, size_t>> MultiStateMappingPath::findEulerCircuit(vector<vector<size_t>>& graph, size_t start, size_t end){
    size_t graphSize = graph.size();
    vector<size_t> unusedEdges;
    unusedEdges.resize(graphSize);
    size_t numOfEdges = getNumberOfEdges(graph);
    for (size_t i = 0; i < graphSize; i++){
        unusedEdges[i] = graph[i].size(); // for each vertex we have its respective number of edges
    }
    // some edges might present twice
    std::map<pair<size_t,size_t>,size_t> edgesOccurrences;
    for (size_t i = 0; i < graph.size(); i++){
        for (size_t j = 0; j < graph[i].size(); j++){
            pair<size_t, size_t> graphEdge(i, graph[i][j]);
            if (edgesOccurrences.find(graphEdge) != edgesOccurrences.end()){
                edgesOccurrences[graphEdge] += 1;
            }else{
                edgesOccurrences[graphEdge] = 1;
            }
        }
    }

    size_t startNode = start;
    std::vector<vector<pair<size_t, size_t>>> curcuits;
    vector<size_t> usedVertices;
    while (numOfEdges > 0){
        auto foundCircuit = findSubCircuit(graph, startNode, unusedEdges, numOfEdges, edgesOccurrences);
        curcuits.push_back(foundCircuit);
        for (size_t j = 0; j < graphSize; j++){
            if (unusedEdges[j] > 0){ // not all the edges of this node are covered
                startNode = j;
                break;
            }
        }
        
        //firstVector.insert(firstVector.begin() + 2, secondVector.begin(), secondVector.end());
    }
    vector<std::pair<size_t, size_t>> concatenatedPath = constructFinalEulerPathFromCircuits(curcuits, start, end);
    return concatenatedPath;


}
/**********************************************************************************************************/
std::vector <std::pair<size_t, size_t>> MultiStateMappingPath::getEulerPath(vector<vector<size_t>>& graph, size_t start, size_t end){
    std::vector <std::pair<size_t, size_t>> eulerPath;
    bool hasEulerPath = false;
    bool hasEulerCircuit = false;
    if (start == end){
        hasEulerCircuit = isEulerCurcuit(graph);

    }else{
        hasEulerPath = isEulerPath(graph, start, end);

    }
    if (hasEulerPath || hasEulerCircuit){
        if (hasEulerPath){
          graph[end].push_back(start);  
        }
        eulerPath = findEulerCircuit(graph, start, end);
        if (hasEulerPath){
            graph[end].pop_back();
            std::pair<size_t, size_t> artificialEdge(end, start);
            auto it = std::find(eulerPath.begin(), eulerPath.end(), artificialEdge);
            size_t index = std::distance(eulerPath.begin(), it);
            size_t startPath = index+1;
            vector<pair<size_t, size_t>> eulerPathCorrected;
            // two scenarios:
            // v-->uv
            // v-->uv-->v
            if (startPath < eulerPath.size()){
                for (size_t i = startPath; i < eulerPath.size(); i++){
                    eulerPathCorrected.push_back(eulerPath[i]);
                }

            }
            for (size_t i = 0; i < index; i++){
                eulerPathCorrected.push_back(eulerPath[i]);
            }
            eulerPath = eulerPathCorrected;


        }

        
    }
    
    return eulerPath;  

} 
/**********************************************************************************************************/
// bit for each edge. For each vertex we have its degree, which gets updated.
vector<vector<pair<size_t, size_t>>> MultiStateMappingPath::chinesePostman(std::vector<vector<size_t>> &graph, size_t start, size_t end, bool &foundPath){
    vector<vector<pair<size_t,size_t>>> disjointPaths; // this object will store the resulted path or optional paths
    vector<pair<size_t,size_t>> bestPath;
    bool circuitCloseEdgeAdded = false;
    if (start != end){
        //if (std::find(graph[end].begin(), graph[end].end(), start) == graph[end].end()){
        graph[end].push_back(start);
        circuitCloseEdgeAdded = true;

        //}

    }

    std::vector<std::pair<size_t, size_t>> inAndOut = getDegreesForEachVertex(graph);

    // now find the delta degrees. Note that in each pair, first is the in-degree, while second is the out degree.
    std::vector<int> deltaDegrees = getDeltaDegrees(inAndOut);
    std::vector<size_t> inVertices; // those with extra incoming edges (should be first in the bipartite graph)
    std::vector<size_t> outVertices; // those with extra outcoming edges (should be second in the bipartite graph)
    fillInAndOutOddVertices(deltaDegrees, inVertices, outVertices);
    // match between indices and edges in the bipartite graph
    std::vector<std::pair<size_t, size_t>> indexToEdge;
    for (size_t i = 0; i < inVertices.size(); i++){
        for (size_t j = 0; j < outVertices.size(); j++){

            indexToEdge.push_back(std::pair<size_t, size_t>(inVertices[i], outVertices[j]));
        }
    }
    // for each edge find its edge (u, v), which is the shortest path between u and v.
    std::unordered_map<size_t, unordered_map<size_t, size_t>> trackPathsFromEachStart;
    std::map<std::pair<size_t, size_t>, double> weights = findWeightsForOddVerticesEdges(graph, start, end, indexToEdge, circuitCloseEdgeAdded,trackPathsFromEachStart);

    std::map<std::pair<size_t, size_t>, size_t> edgesToIndices;
    std::vector<double> weightsPerIndices;
    for (size_t i = 0; i < indexToEdge.size(); i++){
        weightsPerIndices.push_back(weights[indexToEdge[i]]);
        edgesToIndices[indexToEdge[i]] = i;
    }
    // now find the edges between the odd vertices
    bool infinitePath = false;
    if ((inVertices.size() > 0) || (outVertices.size() > 0)){
        bestPath  = findPerfectMatching(edgesToIndices, indexToEdge, weightsPerIndices, deltaDegrees, infinitePath);

    }
    
    std::vector<pair<size_t, size_t>> infiniteEdges;
    // if there are some infinite edges - find them
    if (infinitePath){
        foundPath = false;
        printGraph(graph);
        return disjointPaths;
        // for (size_t i = 0; i < bestPath.size(); i++){
        //     if (weights[bestPath[i]] == std::numeric_limits<double>::infinity()){
        //         infiniteEdges.push_back(bestPath[i]);
        //     }
        // }
    }

    std::map<pair<size_t,size_t>, std::vector<pair<size_t, size_t>>> edgeWithPath;
    for (auto &edge : bestPath){
        // if ((infinitePath) && (std::find(infiniteEdges.begin(), infiniteEdges.end(), edge) != infiniteEdges.end())){
        //     continue;
        // }
        if (std::find(graph[edge.first].begin(), graph[edge.first].end(), edge.second) == graph[edge.first].end()){
            // this edge does not present in the graph, so it represents a path that should be restored from dijkstra
            auto source = edge.first;
            auto dst = edge.second;
            size_t currentNode = dst;
            auto &trackPath = trackPathsFromEachStart[source];
            edgeWithPath[edge];
            while (currentNode != source){
                std::pair<size_t,size_t> edgeInPath(trackPath[currentNode], currentNode);
                edgeWithPath[edge].push_back(edgeInPath); //putting it in a reverse order
                currentNode = trackPath[currentNode];
                
            }
        }
    }
    // if (infinitePath){
    //     // we apply dijkstra on a modified graph, that includes the <dst, source> edge
    //     std::unordered_map<size_t, unordered_map<size_t, size_t>> trackPathsUpdated;
    //     bool isReachable = true;
    //     std::map<std::pair<size_t, size_t>, double> weightsUpdated = findWeightsForOddVerticesEdges(graph, start, end, indexToEdge, false, trackPathsUpdated, isReachable, false);
    //     for (auto &infiniteEdge : infiniteEdges){
    //         auto source = infiniteEdge.first;
    //         auto dst = infiniteEdge.second;
    //         size_t currentNode = dst;
    //         auto &trackInfinitePath = trackPathsUpdated[source];
    //         edgeWithPath[infiniteEdge];
    //         while (currentNode != source){
    //             std::pair<size_t,size_t> edgeInPath(trackInfinitePath[currentNode], currentNode);
    //             edgeWithPath[infiniteEdge].push_back(edgeInPath); //putting it in a reverse order
    //             currentNode = trackInfinitePath[currentNode];
                
    //         }

    //     }
    // } 
    if (circuitCloseEdgeAdded){
        graph[end].pop_back(); // we do it because this edge is added in the eiuler tour function
    }

    // here we add the additional edges that were constructed for the odd vertices
    for (size_t i = 0; i < bestPath.size(); i++){
        graph[bestPath[i].first].push_back(bestPath[i].second);
    }

    // now we can apply the euler path algorithm on the modified graph with extra edges.
    std::vector <std::pair<size_t, size_t>> euler = getEulerPath(graph, start, end);
    std::vector<pair<size_t, size_t>> chinese;
    for (size_t i = 0; i < euler.size(); i++){
        if (edgeWithPath.find(euler[i]) == edgeWithPath.end()){
            // edge exists in the graph
            chinese.push_back(euler[i]);
        }else{ //its a path
            auto &repeatedEdges = edgeWithPath[euler[i]];
            for (int j = repeatedEdges.size()-1; j >=0; j--){
                chinese.push_back(repeatedEdges[j]);
            }

        }

    }
    // if there are disjoint paths, i.e., infinitePath is true, we need to find all the subpaths
    //vector<vector<pair<size_t,size_t>>> disjointPaths;
    // if (infinitePath){
    //     std::pair<size_t,size_t> artificialEdge(end, start);   
    //     vector<pair<size_t,size_t>> subPath;
    //     for (size_t i = 0; i < chinese.size();i++){
    //         if (chinese[i] != artificialEdge){
    //             subPath.push_back(chinese[i]);

    //         }else{
    //             disjointPaths.push_back(subPath);
    //             subPath.clear();
    //         }
    //         if(i == chinese.size()-1){
    //             disjointPaths.push_back(subPath);
    //         }


    //     }
    //}else{
        // there are no disjoint paths, so we return the path results from the euler function without any modifications
    disjointPaths.push_back(chinese);
    //}
    return disjointPaths;
    

}
/**********************************************************************************/
void MultiStateMappingPath::fillRelativeTimeDuration(std::unordered_map<size_t, double> &relativeTimeDuration, vector<double> &dwellingTimes, double totalDurationTime){
    for (size_t i = 0; i < dwellingTimes.size(); i++){
        relativeTimeDuration[i] = dwellingTimes[i]/totalDurationTime;

    }
    
}
/**********************************************************************************/

void MultiStateMappingPath::printGraph(vector<vector<size_t>>& graph){
    std::cout << "Printing Graph!" << std::endl;
    for (size_t i = 0; i < graph.size(); i++){
        std::cout << "Current Vertex is " << i << std::endl;
        for (size_t j = 0; j < graph[i].size(); j++){
            std::cout << "\t" << graph[i][j] << std::endl;
        }
    }
}
/**********************************************************************************/
vector<vector<size_t>> MultiStateMappingPath::createGraphForChinesePostman(std::unordered_map<size_t, size_t> &indicesToNodes, std::unordered_map<size_t, size_t> &nodesToIndices, std::unordered_map<size_t, double> &relativeTimeDuration, std::map<std::pair<size_t, size_t>, double> &transitions, bool &isValidGraph, size_t &start, size_t &end){
    size_t index = 0;
    size_t numberOfVertices = 0;
    for (size_t i = 0; i < relativeTimeDuration.size(); i++){
        if (relativeTimeDuration[i] < EPSILON_THRESHOLD){
            continue;
        }
        indicesToNodes[index] = i;
        nodesToIndices[i] = index;
        index ++;
        numberOfVertices ++;
    }
    vector<vector<size_t>> graph;
    graph.resize(numberOfVertices);
    auto edges = createEdges(relativeTimeDuration, transitions);
    auto it = edges.begin();

    // should check if the source and destination nodes exist as terminals in at least one of the edges
    bool startFound = false;
    bool endFound = false;

    while (it != edges.end()){
        auto &node = it->first;
        if (node == start){
            startFound = true;
        }
        auto &neighbors = edges[node];
        for (size_t i = 0; i < neighbors.size(); i++){
            if (neighbors[i] == end){
                endFound = true;
            }
            graph[nodesToIndices[node]].push_back(nodesToIndices[neighbors[i]]);

        }
        it++;

    }
    isValidGraph = (startFound) && (endFound);
    return graph;

}
/**********************************************************************************/
vector<vector<size_t>> MultiStateMappingPath::createGraphFromPath(vector<pair<size_t, size_t>> &disjointPath, std::unordered_map<size_t, size_t> &indicesToNodesSub, std::unordered_map<size_t, size_t> &nodesToIndicesSub, size_t originalGraphSize){
    vector<vector<size_t>> graph;
    vector<bool> verticesToCount(originalGraphSize, false);

    for (size_t i = 0; i < disjointPath.size(); i++){
        if (i == 0){
            verticesToCount[disjointPath[i].first] = true;
        }
        verticesToCount[disjointPath[i].second] = true;

    }
    int subGraphSize = std::count(verticesToCount.begin(), verticesToCount.end(), true);
    graph.resize(static_cast<size_t>(subGraphSize));

    // fill the maps
    size_t index = 0;
    for (size_t i = 0; i < verticesToCount.size(); i++){
        if (verticesToCount[i]){
            indicesToNodesSub[index] = i;
            nodesToIndicesSub[i] = index;
            index++;
        }
    }
    // create the graph
    for (size_t i = 0; i < disjointPath.size();i++){
        auto &currentNode = nodesToIndicesSub[disjointPath[i].first];
        auto &toNode = nodesToIndicesSub[disjointPath[i].second];
        if (std::find(graph[currentNode].begin(), graph[currentNode].end(), toNode) != graph[currentNode].end()){
            continue;
        }
        graph[currentNode].push_back(toNode);
    }
    return graph;
}

/**********************************************************************************/

vector<size_t> MultiStateMappingPath::findExpectedMappingPathForEachNode(size_t start, size_t end, std::map<std::pair<size_t, size_t>, double> &transitions, vector<double> &dwellingTimes, double totalDurationTime, bool &foundPath){
    vector<size_t> finalPath;
    bool validTSP = false;
    std::unordered_map<size_t, double> relativeTimeDuration;
    fillRelativeTimeDuration(relativeTimeDuration, dwellingTimes, totalDurationTime);
    vector<size_t> pathTSP = TSP(start, end, transitions, totalDurationTime, validTSP, relativeTimeDuration);
    if (!validTSP){
        vector<pair<size_t, size_t>> bestChinesePath;
        std::unordered_map<size_t, size_t> indicesToNodes;
        std::unordered_map<size_t, size_t> nodesToIndices;
        bool isValidGraph = true;
        vector<vector<size_t>> graph = createGraphForChinesePostman(indicesToNodes, nodesToIndices, relativeTimeDuration, transitions, isValidGraph, start, end);
        if (!isValidGraph){
            foundPath = false;
            return finalPath;
        }
        vector<vector<pair<size_t, size_t>>> disjointPaths = chinesePostman(graph, nodesToIndices[start], nodesToIndices[end], foundPath);
        if (!(foundPath)){
            return finalPath;
        }
        if (disjointPaths.size() == 0){
            return finalPath;
        }
        if (disjointPaths.size() == 1){
            bestChinesePath = disjointPaths[0];
            for (size_t i = 0; i < bestChinesePath.size();i++){
                if (i == 0){
                    finalPath.push_back(indicesToNodes[bestChinesePath[i].first]);
                }
                finalPath.push_back(indicesToNodes[bestChinesePath[i].second]);            
            }
        }else{
            double totalExpectedNumberOfTransitions = 0;
            auto itTransitions = transitions.begin();
            while (itTransitions != transitions.end()){
                totalExpectedNumberOfTransitions += transitions[itTransitions->first];
                itTransitions++;
            }
            double heaviestWeightPath = 0;
            for (auto &disjointPath : disjointPaths){
                vector<size_t> nodesInOrder;
                std::unordered_map<size_t, size_t> indicesToNodesSub;
                std::unordered_map<size_t, size_t> nodesToIndicesSub;
                vector<vector<size_t>> subGraph = createGraphFromPath(disjointPath, indicesToNodesSub, nodesToIndicesSub, graph.size());
                vector<vector<pair<size_t, size_t>>> disjointSubPaths = chinesePostman(subGraph, nodesToIndicesSub[nodesToIndices[start]], nodesToIndicesSub[nodesToIndices[end]], foundPath);
                auto currChinesePath = disjointSubPaths[0];
                double currentWeight = 0;
                nodesInOrder.push_back(start);
                std::map<std::pair<size_t, size_t>, double> edgeContribution;

                for (auto &edge : currChinesePath){
                    std::pair<size_t,size_t> tranformedEdge(indicesToNodes[indicesToNodesSub[edge.first]], indicesToNodes[indicesToNodesSub[edge.second]]);
                    edgeContribution[tranformedEdge] = transitions[tranformedEdge]/totalExpectedNumberOfTransitions;
                    nodesInOrder.push_back(tranformedEdge.second);
                }
                auto itPresentedEdges = edgeContribution.begin();
                while (itPresentedEdges != edgeContribution.end()){
                    currentWeight += edgeContribution[itPresentedEdges->first];
                    itPresentedEdges++;

                }
                if (currentWeight > heaviestWeightPath){
                    heaviestWeightPath = currentWeight;
                    bestChinesePath = currChinesePath;
                    finalPath = nodesInOrder;
                }

            }
        }
    }else{
        finalPath = pathTSP;
    }
    return finalPath;

}
/**********************************************************************************/
vector<size_t> MultiStateMappingPath::TSP(size_t start, size_t end, std::map<std::pair<size_t, size_t>, double> &transitions, double totalDurationTime, bool &validPath, std::unordered_map<size_t, double> &relativeTimeDuration){
  std::vector<size_t> bestPath;
  size_t desiredPathId = 0;
  
  size_t numOfNotAdded = 0;
  for (size_t i = 0; i < relativeTimeDuration.size(); i++){
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
  if (edges.size() == 0){
    if (start == end){
        validPath = true;
    }
    return bestPath;
  }
  size_t pathLength = relativeTimeDuration.size()-numOfNotAdded-1; // we don't include the final end node.
  
  std::map<std::pair<size_t, size_t>, double> paths;
  std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> pathReconstruction;
  auto neighbors = edges[start];
  if (neighbors.size() == 0){  
    //bestPath.push_back(end);
    if (start == end){
        validPath = true;
    }
    return bestPath;
  }
  if ((pathLength == 1) && (start != end)){
    if (std::find(neighbors.begin(), neighbors.end(), end) != neighbors.end()){
      bestPath.push_back(start);
      bestPath.push_back(end);
      validPath = true;
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
  
  std::pair<size_t,size_t> bestCandidatePathId; 
  findBestPath(bestCandidatePathId, paths, desiredPathId, edges, transitions, end, validPath);

  // reconstruct the best path
  size_t lengthOfPath = relativeTimeDuration.size()-numOfNotAdded;
  if (start == end){
    lengthOfPath ++;
  }
  if (!validPath){
    return bestPath;
  }
  reconstructBestPath(bestPath, lengthOfPath, pathReconstruction, bestCandidatePathId, start, end);
  return bestPath;


}
/**********************************************************************************/
std::unordered_map<size_t, vector<size_t>> MultiStateMappingPath::createEdges(std::unordered_map<size_t, double> &vertices, std::map<std::pair<size_t, size_t>, double> &transitions){
  std::unordered_map<size_t, vector<size_t>> edges;
  auto it = transitions.begin();
  while (it != transitions.end()){
    auto transition = it->first;
    if ((transitions[transition] <= 0) || (transition.first == transition.second)){
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
std::vector<size_t> MultiStateMappingPath::decimalToBinaryPowers(int decimalNumber) {
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
/******************************************************************************/
void MultiStateMappingPath::findBestPath(std::pair<size_t,size_t> &bestCandidatePathId, std::map<std::pair<size_t, size_t>, double> &paths, size_t desiredPathId, std::unordered_map<size_t, vector<size_t>> &edges, std::map<std::pair<size_t, size_t>, double> &transitions, size_t end, bool &foundPath){
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
void MultiStateMappingPath::reconstructBestPath(std::vector<size_t> &bestPath, size_t lengthOfPath, std::map<std::pair<size_t,size_t>, std::pair<size_t, size_t>> &pathReconstruction, std::pair<size_t,size_t> bestCandidatePathId, size_t start, size_t end){
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