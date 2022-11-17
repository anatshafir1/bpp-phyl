#include <Bpp/Phyl/Likelihood/UndirectedGraph.h>
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;

void printClusters(std::map<uint, vector<Vertex*>> &clusters){
    auto it = clusters.begin();
    while (it != clusters.end()){
        std::cout << "Root is " << it->first << std::endl;
        for (size_t i = 0; i < clusters[it->first].size(); i++){
            std::cout << "\t" << clusters[it->first][i]->getId() << std::endl;
        }
        it ++;
    }
}
void printBestEdges(UndirectedGraph* G, std::map<uint, vector<pair<uint,uint>>> &bestEdges){
    auto it = bestEdges.begin();
    while (it != bestEdges.end()){
        std::cout << "Root is " << it->first << std::endl;
        for (size_t i = 0; i < bestEdges[it->first].size(); i++){
            std::cout << "\t" << bestEdges[it->first][i].first << "->" <<  bestEdges[it->first][i].second << " value is: "<<  G->getEdgeValue(bestEdges[it->first][i]) << std::endl;
        }
        it ++;
    }

}

int main() {
    Vertex* v1 = new Vertex(1);
    Vertex* v2 = new Vertex(2);
    Vertex* v3 = new Vertex(3);
    Vertex* v4 = new Vertex(4);
    Vertex* v5 = new Vertex(5);
    Vertex* v6 = new Vertex(6);
    Vertex* v7 = new Vertex(7);
    Vertex* v8 = new Vertex(8);

    UndirectedGraph* G = new UndirectedGraph(v1);
    G->addNewOrphanVertex(v2);
    G->addNewOrphanVertex(v3);
    G->addNewOrphanVertex(v4);
    G->addNewOrphanVertex(v5);
    G->addNewOrphanVertex(v6);
    G->addNewOrphanVertex(v7);
    G->addNewOrphanVertex(v8);
    G->addEdgeBetweenTwoNodes(v1, v2, 1.5);
    G->addEdgeBetweenTwoNodes(v1, v3, 2.5);
    G->addEdgeBetweenTwoNodes(v3, v4, 3.5);
    G->addEdgeBetweenTwoNodes(v4, v2, 4.5);
    G->addEdgeBetweenTwoNodes(v3, v6, 5.5);
    G->addEdgeBetweenTwoNodes(v3, v5, 6.5);
    G->addEdgeBetweenTwoNodes(v5, v6, 7.5);
    G->addEdgeBetweenTwoNodes(v7, v8, 8.5);
    std::map<uint, vector<Vertex*>> clusters;
    std::map<uint, vector<pair<uint,uint>>> bestEdges;
    G->countNumOfVerticesInClusters(clusters, bestEdges);
    printClusters(clusters);
    printBestEdges(G, bestEdges);
    auto fullyConected = G->findFullyConnectedClusters(clusters);
    auto it = fullyConected.begin();
    while (it != fullyConected.end()){
        if (fullyConected[it->first]){
            std::cout << "Root: " << it->first << "is fully conncted " << std::endl;

        }else{
           std::cout << "Root: " << it->first << "is not fully conncted " << std::endl; 
        }
        

        it ++;
    }


    G->printGraph();
    delete G;
    return 0;





}
