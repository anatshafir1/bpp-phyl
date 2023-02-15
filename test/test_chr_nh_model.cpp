#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <iostream>

using namespace bpp;
using namespace std;

void printNode(uint fatherIndex, std::shared_ptr<PhyloNode> node, std::shared_ptr<PhyloTree> tree){
    shared_ptr<PhyloBranch> branch = tree->getEdgeToFather(node);
    double branch_length = branch->getLength();
    if (tree->isLeaf(node)){
        std::cout << "Node is " << node->getName() << " Father is: N" <<  fatherIndex;
        std::cout << " Length is " << branch_length << std::endl;
    }else{
        std::cout << "Node is N" << tree->getNodeIndex(node) << " Father is: N" <<  fatherIndex;
        std::cout << " Length is " << branch_length << std::endl;
        auto sons = tree->getSons(node);
        for (size_t i = 0; i < sons.size(); i++){
            printNode(tree->getNodeIndex(node), sons[i], tree);
        }
    }
}

void printEditedTree(std::shared_ptr<PhyloTree> tree){
    uint rootIndex = tree->getRootIndex();
    auto sons = tree->getSons(tree->getNode(rootIndex));
    for (size_t i = 0; i < sons.size(); i++){
        printNode(rootIndex, sons[i], tree);
    }

}

int main(){
    Newick reader;

    std::shared_ptr<PhyloTree> tree = std::shared_ptr<PhyloTree>(reader.parenthesisToPhyloTree("(((Sida_cordata:0.005588911,Sida_cordifolia:0.005588911):0.01625967,(Sida_cerradoensis:0.006791234,Sida_angustissima:0.006791234):0.01505734):0.018021629,(Sida_aggregata:0.012374647,Sida_salviifolia:0.01237465):0.02749556);"));
    tree->scaleTree(23.7447);
    std::shared_ptr<PhyloNode> son= tree->getNode(9);
    auto edgeIndex =  tree->getIncomingEdges(tree->getNodeIndex(son))[0]; 
    auto fatherIndex = tree->getFatherOfEdge(edgeIndex);
    std::shared_ptr<PhyloNode> father = tree->getNode(fatherIndex);
    shared_ptr<PhyloBranch> branch = tree->getEdgeToFather(son);
    double original_edge_length = branch->getLength();
    tree->createNodeOnEdge(tree->getEdgeIndex(branch), original_edge_length/3);
    printEditedTree(tree);
    std::cout << "******" << std::endl;
    std::vector<uint> nodes = tree->getNodeIndexes(tree->getAllNodes());
    std::vector<uint> edges = tree->getEdgeIndexes(tree->getAllEdges());
    for (size_t i = 0; i < nodes.size(); i++){
        if (tree->isLeaf(tree->getNode(nodes[i]))){
            std::cout << "Node is " << tree->getNode(nodes[i])->getName() << std::endl;
            std::cout << "\tnode id is N" << nodes[i] << std::endl;
        }else{
            std::cout << "Node is N" << nodes[i] << std::endl;
            auto sons = tree->getSons(tree->getNode(nodes[i]));
            for (size_t j = 0; j < sons.size(); j++){
                std::cout << "\tSon is N" << tree->getNodeIndex(sons[j]) << std::endl;
            }
        }
    }
    std::cout << "##########################" << std::endl;
    for (size_t i = 0; i < edges.size(); i++){
        auto edge = tree->getEdge(edges[i]);
        std::cout << "Edge is " << edges[i] << std::endl;
        std::cout << "\tBranch length is " << edge->getLength() << std::endl;
        std::cout << "\tFather of edge is N" << tree->getNodeIndex(tree->getFatherOfEdge(edge)) << std::endl;
        auto sonNode = tree->getSon(edge);
        if (tree->isLeaf(sonNode)){
            std::cout << "\tSon of edge is " << sonNode->getName() << std::endl;
        }else{
            std::cout << "\tSon of edge is N" << tree->getNodeIndex(sonNode) << std::endl;
        }
    }
    std::cout << "Root index is " << tree->getRootIndex() << std::endl;
    return 0;

}
