//from bpp-core
#include <Bpp/Version.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/App/BppApplication.h>


//from bpp-seq
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/chrFasta.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

//from bpp-phyl
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/MLAncestralStateReconstruction.h>
#include <Bpp/Phyl/Likelihood/MarginalNonRevAncestralStateReconstruction.h>
#include <Bpp/Phyl/Likelihood/ChromosomeNumberOptimizer.h>
#include <Bpp/Phyl/Mapping/ComputeChangesExpectations.h>
#include <Bpp/Phyl/Mapping/ComputeChromosomeTransitionsExp.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/App/ChromEvolOptions.h>

//standard libraries
#include <string>
#include <vector>
#include <iostream>

using namespace bpp;
using namespace std;


//Functions for initialization of the model
VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* initialAlpha);
VectorSiteContainer* getCharacterData(const std :: string &path, unsigned int* numberOfUniqueStates, unsigned int* chrRange);
void setMaxChrNum(unsigned int maxNumberOfChr);
void setMinChrNum(unsigned int minNumberOfChr);
void rescale_tree(TreeTemplate<Node>* tree, unsigned int chrRange);
TreeTemplate<Node>* getTree(const std :: string &path, unsigned int numOfUniqueCharacterStates);
void getMaxParsimonyUpperBound(VectorSiteContainer* vsc, TreeTemplate<Node>* tree, double* parsimonyScore);

//core functions of ChromEvol
void runChromEvol(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange);
void optimizeLikelihoodMultiStartPoints(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange);
ChromosomeNumberOptimizer* optimizeLikelihoodMultiStartPoints(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange);
void getJointMLAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik);
std::map<int, std::map<size_t, VVdouble>> getMarginalAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik);
void computeExpectations(DRNonHomogeneousTreeLikelihood* lik, std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int numOfSimulations);

// functions to print the tree with ancestral reconstruction
void printTreeWithStates(DRNonHomogeneousTreeLikelihood* lik, TreeTemplate<Node> tree, std::map<int, std::vector<size_t> > ancestors, std::map<int, map<size_t, std::vector<double>>>* probs = 0);
string printTree(const TreeTemplate<Node>& tree);
string nodeToParenthesis(const Node& node);

/******************************************************************************/
VectorSiteContainer* getCharacterData (const string& path, unsigned int* numberOfUniqueCharacterStates, unsigned int* chrRange){
    ChromosomeAlphabet* alphaInitial = new ChromosomeAlphabet(ChromEvolOptions::minAlpha_, ChromEvolOptions::maxAlpha_);
    VectorSequenceContainer* initialSetOfSequences = chrFasta::readSequencesFromFile(path, alphaInitial);
    size_t numOfSequences = initialSetOfSequences->getNumberOfSequences();
    vector <string> sequenceNames = initialSetOfSequences->getSequencesNames();

    unsigned int maxNumberOfChr = 1; //the minimal number of chromosomes cannot be zero
    unsigned int minNumOfChr = ChromEvolOptions::maxAlpha_;

    std::vector <int> UniqueCharacterStates;
    cout<<"vector size is "<< UniqueCharacterStates.size()<<endl;
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = initialSetOfSequences->getSequence(sequenceNames[i]);
        int character = seq.getValue(0);
        if (character == -1){
            continue;
        }
        if (character == static_cast<int>(ChromEvolOptions::maxAlpha_)+1){
            continue;
        }
        // if it is a composite state
        if (character > static_cast<int>(ChromEvolOptions::maxAlpha_) +1){
            const std::vector<int> compositeCharacters = alphaInitial->getSetOfStatesForAComposite(character);
            for (size_t j = 0; j < compositeCharacters.size(); j++){
                if ((unsigned int) compositeCharacters[j] > maxNumberOfChr){
                    maxNumberOfChr = compositeCharacters[j];
                }
                if ((unsigned int) compositeCharacters[j] < minNumOfChr){
                    minNumOfChr = compositeCharacters[j];
                }
                
            }
            continue;
        }

        if (!std::count(UniqueCharacterStates.begin(), UniqueCharacterStates.end(), character)){
            UniqueCharacterStates.push_back(character);
        }
        if ((unsigned int) character > maxNumberOfChr){
            maxNumberOfChr = character;
        }
        if ((unsigned int) character < minNumOfChr){
            minNumOfChr = character;
        }

    }
    *numberOfUniqueCharacterStates = (unsigned int)UniqueCharacterStates.size() + alphaInitial->getNumberOfCompositeStates();
    *chrRange = maxNumberOfChr - minNumOfChr;
    cout <<"Number of unique states is " << *numberOfUniqueCharacterStates <<endl;

    setMaxChrNum(maxNumberOfChr);
    setMinChrNum(minNumOfChr);

    VectorSiteContainer* vsc = resizeAlphabetForSequenceContainer(initialSetOfSequences, alphaInitial);
    delete initialSetOfSequences;
    delete alphaInitial;
    return vsc;
}

/****************************************************************************/
void setMaxChrNum(unsigned int maxNumberOfChr){
    if (ChromEvolOptions::maxChrNum_ < 0){
        ChromEvolOptions::maxChrNum_ = maxNumberOfChr + abs(ChromEvolOptions::maxChrNum_);
    }else{
        if ((int)maxNumberOfChr > ChromEvolOptions::maxChrNum_){
            ChromEvolOptions::maxChrNum_ = maxNumberOfChr;
        }
    }

}
/****************************************************************************/
void setMinChrNum(unsigned int minNumberOfChr){
    if (ChromEvolOptions::minChrNum_ < 0){
        ChromEvolOptions::minChrNum_ = minNumberOfChr - abs(ChromEvolOptions::minChrNum_);
    }else{
        if ((int)minNumberOfChr < ChromEvolOptions::minChrNum_){
            ChromEvolOptions::minChrNum_ = minNumberOfChr;
        }

    }
}

/*****************************************************************************/
VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc, ChromosomeAlphabet* alphaInitial){
    size_t numOfSequences = vsc->getNumberOfSequences();
    vector <string> sequenceNames = vsc->getSequencesNames();
    ChromosomeAlphabet* new_alphabet = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
        // fill with composite values
    if (alphaInitial->getNumberOfCompositeStates() > 0){
        const std::map <int, std::map<int, double>> compositeStates = alphaInitial->getCompositeStatesMap();
        std::map <int, std::map<int, double>>::const_iterator it = compositeStates.begin();
        while (it != compositeStates.end()){
            int compositeState = it->first;
            std::string charComposite = alphaInitial->intToChar(compositeState);
            new_alphabet->setCompositeState(charComposite);
            it++;
        }
    }
    VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(new_alphabet);
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = vsc->getSequence(sequenceNames[i]);
        BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), new_alphabet);
        resized_alphabet_site_container->addSequence(new_seq);

    }
    return resized_alphabet_site_container;
}
/****************************************************************************/
TreeTemplate<Node>* getTree(const string& path, unsigned int NumOfUniqueCharacterStates){
    Newick newick;
    TreeTemplate<Node>* tree = newick.readTree(path);
    rescale_tree(tree, NumOfUniqueCharacterStates);
    return tree;

}
/****************************************************************************/
void rescale_tree(TreeTemplate<Node>* tree, unsigned int chrRange){
    double scale_tree_factor = 1.0;
    //string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    //std :: cout << tree_str << endl;
    bool rooted = tree->isRooted();
    if (!rooted){
        throw UnrootedTreeException("The given input tree is unrooted. Tree must be rooted!", tree);
    }
    if (ChromEvolOptions::branchMul_ == 1.0){
        return;
    }else{
        //tree must be rescaled
        double treeLength = tree->getTotalLength();
        if (ChromEvolOptions::branchMul_ == 999){
            scale_tree_factor = (double)chrRange/treeLength;
        }
        tree->scaleTree(scale_tree_factor);

    }

}

/******************************************************************************/
void getMaxParsimonyUpperBound(VectorSiteContainer* vsc, TreeTemplate<Node>* tree, double* parsimonyBound){

    TreeTemplate<Node>* treeForParsimonyBound = tree->clone();
    DRTreeParsimonyScore maxParsimonyObject = DRTreeParsimonyScore(*treeForParsimonyBound, *vsc);
    *parsimonyBound = (maxParsimonyObject.getScore())/(tree->getTotalLength());
    delete treeForParsimonyBound;
    return;   

}

/******************************************************************************/
void runChromEvol(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange){
    // optimize likelihood
    ChromosomeNumberOptimizer* chrOptimizer = optimizeLikelihoodMultiStartPoints(alpha, vsc, tree, chrRange);
    std::vector <DRNonHomogeneousTreeLikelihood> lik_vec = chrOptimizer->getVectorOfLikelihoods();
    // get joint ML ancestral reconstruction
    getJointMLAncestralReconstruction(&lik_vec[0]);
    //get Marginal ML ancestral reconstruction, and with the help of them- calculate expectations of transitions
    std::map<int, std::map<size_t, VVdouble>>  jointProbabilitiesFatherSon = getMarginalAncestralReconstruction(&lik_vec[0]);
    //compute expectations
    computeExpectations(&lik_vec[0], jointProbabilitiesFatherSon, ChromEvolOptions::NumOfSimulations_);
    //Clear the vector of likelihoods entirely
    delete chrOptimizer;


}
/******************************************************************************/
ChromosomeNumberOptimizer* optimizeLikelihoodMultiStartPoints(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange){
    
    vector<double> modelParams;
    modelParams.reserve(ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS);
    ChromEvolOptions::initVectorOfChrNumParameters(modelParams);
    double parsimonyBound = 0;
    if (ChromEvolOptions::maxParsimonyBound_){
        getMaxParsimonyUpperBound(vsc, tree, &parsimonyBound);
    }
    bool calculateDerivatives = true;
    if (ChromEvolOptions::optimizationMethod_ == "Brent"){
        calculateDerivatives  = false;
    }

    ChromosomeNumberOptimizer* opt = new ChromosomeNumberOptimizer(tree, alpha, vsc, chrRange);
    opt->initModels(modelParams, parsimonyBound, ChromEvolOptions::optimizeBaseNumber_, ChromEvolOptions::rateChangeType_, ChromEvolOptions::seed_, ChromEvolOptions::OptPointsNum_[0], calculateDerivatives, ChromEvolOptions::fixedFrequenciesFilePath_);

    //initialize all the optimization specific parameters
    opt->initOptimizer(ChromEvolOptions::OptPointsNum_, ChromEvolOptions::OptIterNum_, ChromEvolOptions::optimizationMethod_, ChromEvolOptions::baseNumOptimizationMethod_,
        ChromEvolOptions::tolerance_, ChromEvolOptions::standardOptimization_, ChromEvolOptions::BrentBracketing_, 
        ChromEvolOptions::probsForMixedOptimization_, ChromEvolOptions::fixedParams_);
    //optimize models
    opt->optimize();
    // it is safe to delete the chrOptimizer, because the destructor doesn't delete nothing associated with the vector of likelihoods
    return opt;
       
}

/******************************************************************************/
std::map<int, std::map<size_t, VVdouble>> getMarginalAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik){
    std::cout << "Marginal Ancestral Reconstruction"<<endl;
    MarginalNonRevAncestralStateReconstruction* ancr = new MarginalNonRevAncestralStateReconstruction(lik);
    ancr->computePosteriorProbabilitiesOfNodesForEachStatePerSite();
    std::map<int, std::vector<size_t> > ancestors = ancr->getAllAncestralStates();
    std::map<int, map<size_t, std::vector<double>>>* probs = ancr->getPosteriorProbForAllNodesAndStatesPerSite();
    printTreeWithStates(lik, lik->getTree(), ancestors, probs);
    std::map<int, std::map<size_t, VVdouble>> jointProbabilitiesFatherSon = ancr->getAllJointFatherNodeProbabilities();
    delete ancr;
    return jointProbabilitiesFatherSon;

}
/*****************************************************************************/
void computeExpectations(DRNonHomogeneousTreeLikelihood* lik, std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int numOfSimulations){

    ComputeChromosomeTransitionsExp* expCalculator = new ComputeChromosomeTransitionsExp(lik, jointProbabilitiesFatherSon, ChromEvolOptions::jumpTypeMethod_);
    expCalculator->runSimulations(numOfSimulations);
    expCalculator->computeExpectationPerType();
    expCalculator->printResults();

    //delete
    delete expCalculator;
}
/******************************************************************************/
void getJointMLAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik){
    std::cout << "ML Joint Ancestral Reconstruction"<<endl;
    TransitionModel* model = lik->getSubstitutionModelSet()->getModel(0);
    std::vector<double> rootFreqs = lik->getRootFrequencies(0);
    std::map<int, VVVdouble> Pijt;
    TreeTemplate<Node> tree = lik->getTree();
    vector <int> nodesIds = tree.getNodesId();
    for (size_t n = 0; n < nodesIds.size(); n++){
        Pijt[nodesIds[n]] = lik->getTransitionProbabilitiesPerRateClass(nodesIds[n], 0);
    }
    MLAncestralStateReconstruction* ancr = new MLAncestralStateReconstruction(lik, model, rootFreqs, &Pijt);
    ancr->computeJointLikelihood();
    std::map<int, std::vector<size_t> > ancestors = ancr->getAllAncestralStates();
    printTreeWithStates(lik, lik->getTree(), ancestors);

    delete ancr;

}
/******************************************************************************/
void printTreeWithStates(DRNonHomogeneousTreeLikelihood* lik, TreeTemplate<Node> tree, std::map<int, std::vector<size_t> > ancestors, std::map<int, map<size_t, std::vector<double>>>* probs){
    std::vector<int> nodesIds = tree.getNodesId();
    for (size_t n= 0; n < nodesIds.size(); n++){
        for (size_t i = 0; i < ancestors[nodesIds[0]].size(); i++){
            size_t state = ancestors[nodesIds[n]][i] + (dynamic_cast<const ChromosomeAlphabet*>(lik->getAlphabet()))->getMin();
            std::string prevName;
            
            if (tree.isLeaf(nodesIds[n])){
                prevName = tree.getNodeName(nodesIds[n]);
                const std::string newName = (prevName + "-"+ std::to_string(state));
                tree.setNodeName(nodesIds[n], newName);
                if (probs){
                    double postProb = (*probs)[nodesIds[n]][i][ancestors[nodesIds[n]][i]];
                    std::cout <<newName<< " state index: " << state << " : " << postProb <<endl;

                }else{
                    std::cout <<newName<< " state index: " << state <<endl;
                }
            }else{
                prevName = "N" + std::to_string(nodesIds[n]);
                const std::string newName = (prevName + "-"+ std::to_string(state));
                tree.setNodeName(nodesIds[n], newName);
                if (probs){
                    double postProb = (*probs)[nodesIds[n]][i][ancestors[nodesIds[n]][i]];
                    std::cout <<newName<< " state index: " << state << " : " << postProb <<endl;

                }else{
                    std::cout <<newName<< " state index: " << state <<endl;
                }
            }
            
            
        }
        
    }
    string tree_str = printTree(tree);
    std:: cout << tree_str << endl;

}
/******************************************************************************/

string printTree(const TreeTemplate<Node>& tree)
{
  ostringstream s;
  s << "(";
  const Node* node = tree.getRootNode();
  if (node->isLeaf() && node->hasName()) // In case we have a tree like ((A:1.0)); where the root node is an unamed leaf!
  {
    s << node->getName();
    for (size_t i = 0; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i));
    }
  }
  else
  {
    s << nodeToParenthesis(*node->getSon(0));
    for (size_t i = 1; i < node->getNumberOfSons(); ++i)
    {
      s << "," << nodeToParenthesis(*node->getSon(i));
    }
  }
  s << ")";
  if (node->hasDistanceToFather())
    s << ":" << node->getDistanceToFather();
 
  s << tree.getRootNode()->getName();
  s << ";" << endl;
  return s.str();
}

/******************************************************************************/
string nodeToParenthesis(const Node& node)
{
  ostringstream s;
  if (node.isLeaf())
  {
    s << node.getName();
  }
  else
  {
    s << "(";
    s << nodeToParenthesis(*node[0]);
    for (int i = 1; i < static_cast<int>(node.getNumberOfSons()); i++)
    {
      s << "," << nodeToParenthesis(*node[i]);
    }
    s << ")";
  }
  if (! node.isLeaf()){
      s << node.getName();
  }
  
  if (node.hasDistanceToFather())
    s << ":" << node.getDistanceToFather();
  return s.str();
}
/***********************************************************************************************************/


int main(int args, char **argv) {

    if (args == 1){
        std::cout << "No arguments provided"<<endl;
        return 0;
    }
    try{
        BppApplication ChromEvol(args, argv, "ChromEvol");
        ChromEvolOptions::initAllParameters(ChromEvol);
        unsigned int numberOfUniqueStates = 0;
        unsigned int chrRange = 0;
        VectorSiteContainer* vsc = getCharacterData(ChromEvolOptions::characterFilePath_, &numberOfUniqueStates, &chrRange);
        const ChromosomeAlphabet* alpha = dynamic_cast<const ChromosomeAlphabet*>(vsc->getAlphabet());
        TreeTemplate<Node>* tree = getTree(ChromEvolOptions::treeFilePath_, numberOfUniqueStates);
        std::cout << "****** Max allowed chromosome number: "<< ChromEvolOptions::maxChrNum_ <<endl;
        std::cout << "****** Max observed chromosome range: " << chrRange <<endl;
        runChromEvol(alpha, vsc, tree, chrRange);
        delete vsc;
        delete tree;
        delete alpha;

    }
    catch (exception& e)
    {
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}