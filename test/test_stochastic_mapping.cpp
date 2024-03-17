// Tester for StochasticMapping implementation

// From the STL:
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

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

#define STATE "state"
/******************************************************************************/



void checkIfMappingLegal(const StochasticMapping* stocMapping, const Tree* mapping, const TreeTemplate<Node>* baseTree, const TreeLikelihood* tl)
{
/*
    vector<const Node*> origNodes = baseTree->getNodes();
    vector<const Node*> mappingNodes = dynamic_cast<const TreeTemplate<Node>*>(mapping)->getNodes();
    const SiteContainer* leafsStates = tl->getData();
    for (size_t i=0; i<origNodes.size(); ++i)
    {
        if (origNodes[i]->isLeaf())
        {
            string nodeName = origNodes[i]->getName();
            size_t origNodeState = static_cast<size_t>(tl->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(nodeName);
            size_t nodeStateInMapping = StochasticMapping::getNodeState(nodeInMapping);
            if ((nodeStateInMapping < 2) & (origNodeState != nodeStateInMapping)) // if the node's state corresaponds to a concrete character state (not unknown character or a combination of state and rate in the case of a markov modulated model)
            {
                throw Exception("Leafs states not maintained in mapping");
            }
        }
    }
            
    // make sure that branch lengths are maintained in the mapping

    for (size_t i=0; i<origNodes.size(); ++i)
    {
        if (origNodes[i]->hasFather())
        {
            double origBranchLength = origNodes[i]->getDistanceToFather();
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origNodes[i]->getName());
            const Node* origNodeFatherInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origNodes[i]->getFather()->getName());
            double branchLengthInMapping = 0;
            const Node* curNode = nodeInMapping;
            while (curNode != origNodeFatherInMapping)
            {
                if (curNode->hasFather())
                {
                    branchLengthInMapping += curNode->getDistanceToFather(); // failes on node with id 9, but returns exception on node with id 6 (which is the root)
                    curNode = curNode->getFather();
                }
                else

                {
                    branchLengthInMapping = origBranchLength;
                    curNode = origNodeFatherInMapping;
                }
            }
            if (abs(branchLengthInMapping - origBranchLength) > 0.0001) // expected history fails in the root - but the branch length of the branch coming from the root is insignificant as it isn't included in the tree. what happened here?
            {
                throw Exception("branch lengths not maintained in mapping");
            }
        }
    }
    // make sure there are no two nodes in a row such that both don't exist in the base tree and both recieve the same state (indicator of illegal transition)
    for (size_t i=0; i<origNodes.size(); ++i)
    {
        if (origNodes[i]->hasFather())
        {
            const Node* nodeInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origNodes[i]->getName());
            string origFatherName = origNodes[i]->getFather()->getName();
            const Node* origNodeFatherInMapping = (dynamic_cast<const TreeTemplate<Node>*>(mapping))->getNode(origFatherName);
            const Node* curNode = nodeInMapping;
            while (curNode->getFather() != origNodeFatherInMapping)
            {
                string curNodeFatherName = curNode->getFather()->getName();
                size_t curNodeState = StochasticMapping::getNodeState(curNode);
                size_t nextNodeState = StochasticMapping::getNodeState(curNode->getFather());
                if (curNodeState == nextNodeState && ((curNode->getName()).find("mapping") != std::string::npos) && ((curNode->getFather()->getName()).find("mapping") != std::string::npos)) // such transition is permitted in case the father is the root

                {
                    throw Exception("illegal transitions in the mapping");
                }
                curNode = curNode->getFather();
            }
        }
    }
*/
}

void giveNamesToInternalNodes(Tree* tree)
{
    // TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(tree);
    // vector<Node*> nodes = ttree->getNodes();
    // for (size_t i=0; i<nodes.size(); ++i) {
    //     if (!nodes[i]->hasName())
    //         nodes[i]->setName("_baseInternal_" + TextTools::toString(nodes[i]->getId()));
    // }  
}

void setNodeState(Node* node, size_t state)
{
    // BppInteger* stateProperty = new BppInteger(static_cast<int>(state));
    // node->setNodeProperty(STATE, *stateProperty);
    // delete stateProperty;
}



void computePosteriors(VVDouble& posteriorProbabilities, Tree* baseTree, RHomogeneousTreeLikelihood* tl)
{
/*    // some auxiliiary variables

    size_t statesNum = tl->getNumberOfStates();
    const TransitionModel* model = tl->getModelForSite(0,0); // this calls assumes that all the sites and all the branches are assoiacted with the same node

    const SiteContainer* leafsStates = tl->getData();
    TreeTemplate<Node>* ttree = dynamic_cast<TreeTemplate<Node>*>(baseTree); 
    vector<Node*> nodes = ttree->getNodes();
    // compute the posterior probabilities according to Felsenstein prunnig algorithm: for each node nodes[i] and state s compute: P(Data[leafs under node[i]]|node[i] has state s] 
    for (size_t i=0; i<nodes.size(); ++i) // traverse the tree in post-order

    {
        int nodeId = nodes[i]->getId();
        string nodeName = nodes[i]->getName();
        if (nodes[i]->isLeaf()) // if the node is a leaf, set the posterior probability of its state to 1, and the rest ot 0

        {
            size_t leafState = static_cast<int>(tl->getAlphabetStateAsInt(leafsStates->getSequence(nodeName).getValue(0)));
            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                if (nodeState != leafState)
                {
                    posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] = 0;
                }
                else {
                    posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] = 1; 
                }
            }   
        }
        else                   // if the node is internal, follow the Felesenstein computation rule to compute the posterior probability

        {   
            double dataProb = 0;
            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                double fullProb = 1;
                for (size_t j=0; j<(nodes[i]->getNumberOfSons()); ++j) // for each son of the node, sum over the probabilities of all its assignments given its father's state (i.e, nodeState)
                {
                    double sonProb = 0;
                    double bl = nodes[i]->getSon(j)->getDistanceToFather();
                    for(size_t sonState=0; sonState<statesNum; ++sonState)
                    {
                        sonProb += model->Pij_t(nodeState, sonState, bl) * posteriorProbabilities[nodeIdToIndex[nodes[i]->getSon(j)->getId()]][sonState];
                    }
                    fullProb *= sonProb;
                }
                posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] = fullProb;
                dataProb += posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState];
            }
            // now, compute from the so far compued partial likelihoods the posterior probabilities by dividing by the probability of the data (prior(data)=1 in ML world)
            // because the sum of partial likelihoods is in fact the probablity of the data, it is sufficient to standardize the vector

            for (size_t nodeState=0; nodeState<statesNum; ++nodeState)
		    {
                posteriorProbabilities[nodeId][nodeState] = posteriorProbabilities[nodeIdToIndex[nodeId]][nodeState] / dataProb;
            }    
		}
    }
*/}
/*
std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*tree);
    string fixedRootFreqPath = ChromEvolOptions::fixedFrequenciesFilePath_;
    bool weightedRootFreqs;
    std::map<uint, std::map<int, vector<string>>> mapOfParamsNamesPerModelType;
    LikelihoodUtils::setParamsNameInForMultiProcess(mapOfParamsNamesPerModelType, modelParams);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[1].second, modelParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    if (fixedRootFreqPath == "none"){
        weightedRootFreqs = true;
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        weightedRootFreqs = false;
        vector <double> rootFreqs = LikelihoodUtils::setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
        std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);
    }
*/
/******************************************************/
void printNodeRec(std::shared_ptr<PhyloTree> tree, uint nodeId, bool names){
    if ((!(tree->isLeaf(nodeId)))){
        auto sons = tree->getSons(nodeId);
        if (names){
            std::cout << "current node is " << (tree->getNode(nodeId))->getName();
            std::cout << " id is " << nodeId << std::endl;

        }else{
            std::cout << "current node is " << nodeId << std::endl;
        }
        
        std::cout << "\tSons are: " << std::endl;
        for (size_t i = 0; i < sons.size(); i++){
            shared_ptr<PhyloBranch> branch=tree->getEdgeToFather(sons[i]);
            if (names){
                std::cout << "\t\t" << (tree->getNode(sons[i]))->getName() << " branch length is " << branch->getLength();
                std::cout << " id is " << sons[i] << std::endl;

            }else{
                std::cout << "\t\t" <<sons[i] << " branch length is " << branch->getLength() << std::endl;
            }
            
        }
        for (size_t i = 0; i < sons.size(); i++){
            printNodeRec(tree, sons[i], names);
        }
    }
}
/******************************************************/
void printTree(std::shared_ptr<PhyloTree> tree, bool names){
    uint rootNodeId = tree->getRootIndex();
    printNodeRec(tree, rootNodeId, names);

}
/******************************************************/
void testClone(std::shared_ptr<PhyloTree> tree){
    std::shared_ptr<PhyloTree> clonedTree = std::shared_ptr<PhyloTree>(tree->clone());
    auto nodes = clonedTree->getAllNodes();
    uint nodeId;
    for (size_t i = 0; i < nodes.size(); i++){
        if (clonedTree->getNodeIndex(nodes[i]) != clonedTree->getRootIndex()){
            nodeId = clonedTree->getNodeIndex(nodes[i]);
            break;
        }
    }
    auto edge_to_fragment = clonedTree->getEdgeToFather(nodeId);      
    double weightFather = 0.5;
    double branchLength = edge_to_fragment->getLength();
    uint newNodeId = clonedTree->createNodeOnEdge(clonedTree->getEdgeIndex(edge_to_fragment), weightFather* branchLength);
    printTree(clonedTree, true);
    std::cout << "***** ***** *****" << std::endl;
    printTree(tree, true);


}

/******************************************************/
void stochasticMapping(){
    //fix seed for debugging purposes
    double seedUb = 10000000;
    //double mySeed = RandomTools::giveRandomNumberBetweenZeroAndEntry(seedUb);
    RandomTools::setSeed(static_cast<long int>(seedUb));
    Newick reader;
    shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
    std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);
    auto clonedTree = pTree->clone();
    delete clonedTree;
        
    // create a binary model
    const BinaryAlphabet* alphabet = new BinaryAlphabet();
    double mu = 1.;
    double pi0 = 0.5;
    auto twoParamModel = std::make_shared<TwoParameterBinarySubstitutionModel>(alphabet,mu,pi0);
    //std::shared_ptr<ReversibleSubstitutionModel> nestedModel = std::dynamic_pointer_cast<ReversibleSubstitutionModel>(twoParamModel);
    std::shared_ptr<DiscreteDistribution> rdist = std::shared_ptr<DiscreteDistribution>(new GammaDiscreteRateDistribution(1, 1.0));
    vector <double> rootFreqs;
    rootFreqs.push_back(0.25);
    rootFreqs.push_back(0.75);
    std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(twoParamModel->getStateMap(), false)), rootFreqs);
    std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies);        

    // process character data
    std::shared_ptr<VectorSiteContainer> sites = std::make_shared<VectorSiteContainer>(alphabet);
	sites->addSequence(BasicSequence("S1", "1", alphabet));
	sites->addSequence(BasicSequence("S2", "1", alphabet));
	sites->addSequence(BasicSequence("S3", "0", alphabet));
	sites->addSequence(BasicSequence("S4", "0", alphabet));
	sites->addSequence(BasicSequence("S5", "1", alphabet));

    auto nodes = pTree->getAllNodes();
    std::vector<uint> modelNodes;
    for (size_t i = 0; i < nodes.size(); i++){
        auto nodeId = pTree->getNodeIndex(nodes[i]);
        if (nodeId == pTree->getRootIndex()){
            continue;
        }
        modelNodes.push_back(nodeId);
    }
    subProSim->addModel(std::shared_ptr<TwoParameterBinarySubstitutionModel>(twoParamModel->clone()), modelNodes);
		        
    // create tree likelihood function
    Context context;
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *sites->clone(), *subProSim->clone());
    SingleProcessPhyloLikelihood ntl(context, lik);
    auto lik_val = ntl.getValue();
    std::cout << "Likelihood is: " << lik_val << std::endl;


    unsigned int mappingsNum = 3;
    StochasticMapping* stm = new StochasticMapping(lik, mappingsNum);//ChromEvolOptions::NumOfSimulations_);
    stm->generateStochasticMapping();
    std::vector<std::shared_ptr<PhyloTree>> mappings;
    for (size_t i = 0; i < mappingsNum; i++){
        auto mappingTree = stm->createMappingHistoryTree(i);
        mappings.push_back(mappingTree);
        printTree(mappingTree, true);
        std::cout << "****************" << std::endl;



    }
    std::cout << "Expected mapping tree: " << std::endl;
        
    auto expected_mapping = stm->generateExpectedMapping();
    printTree(expected_mapping, true);
    // print the trees
    delete stm;

}


/******************************************************/


int main() 
{
   try

    {
        Newick reader;
        shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
        testClone(pTree);


 
    }catch (exception & e){
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}
