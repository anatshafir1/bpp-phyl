
#include "ChromEvolOptions.h"
using namespace bpp;
using namespace std;

TreeTemplate<Node>* ChromEvolOptions::tree_;
const VectorSiteContainer* ChromEvolOptions::vsc_;
const ChromosomeAlphabet* ChromEvolOptions::alpha_;
int ChromEvolOptions::maxChrNum_;
int ChromEvolOptions::minChrNum_;
double ChromEvolOptions::branchMul_;
std::vector <unsigned int> ChromEvolOptions::OptPointsNum_;
std::vector <unsigned int> ChromEvolOptions::OptIterNum_;
double ChromEvolOptions::constGain_;
double ChromEvolOptions::constLoss_;
double ChromEvolOptions::constDupl_;
double ChromEvolOptions::constDemiDupl_;
double ChromEvolOptions::tolerance_;
unsigned int ChromEvolOptions::maxIterations_;
bool ChromEvolOptions::maxParsimonyBound_;
unsigned int ChromEvolOptions::maxAlpha_;
unsigned int ChromEvolOptions::minAlpha_;

/*************************************************************************/
void ChromEvolOptions::initAllParameters(BppApplication& ChromEvol){
    initDefaultParameters();
    initParametersFromFile(ChromEvol);

}
/*************************************************************************/
ChromEvolOptions::~ChromEvolOptions(){
    delete tree_;
    delete vsc_;
    delete alpha_;
}
/*************************************************************************/
void ChromEvolOptions::initDefaultParameters(){
    maxAlpha_ = 500;
    minAlpha_ = 1;
}
/*************************************************************************/
void ChromEvolOptions::initParametersFromFile(BppApplication& ChromEvol){
    maxChrNum_ = ApplicationTools::getIntParameter("_maxChrNum", ChromEvol.getParams(), -10, "", true, 0);
    minChrNum_ = ApplicationTools::getIntParameter("_minChrNum", ChromEvol.getParams(), 1, "", true, 0);
    string pathForCharacterData = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    string pathForTree = ApplicationTools::getAFilePath("_treeFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    unsigned int numberOfUniqueCharacterStates = 0;
    vsc_ = getCharacterData(pathForCharacterData, &numberOfUniqueCharacterStates);
    alpha_ = dynamic_cast<const ChromosomeAlphabet*>(vsc_->getAlphabet());
    branchMul_ = ApplicationTools::getDoubleParameter("_branchMul", ChromEvol.getParams(), 999, "", true, 0);
    tree_ = getTree(pathForTree, numberOfUniqueCharacterStates);
    maxIterations_ = (unsigned int)ApplicationTools::getIntParameter("_maxOptimizationItarations", ChromEvol.getParams(), 5, "", true, 0);
    tolerance_ = ApplicationTools::getDoubleParameter("_tolParamOptimization", ChromEvol.getParams(), 0.01, "", true, 0);
    constGain_ = ApplicationTools::getDoubleParameter("_gainConstR", ChromEvol.getParams(), IgnoreParam, "", true, 0);
    constLoss_ = ApplicationTools::getDoubleParameter("_lossConstR", ChromEvol.getParams(), IgnoreParam, "", true, 0);
    constDupl_ = ApplicationTools::getDoubleParameter("_duplConstR", ChromEvol.getParams(), IgnoreParam, "", true, 0);
    constDemiDupl_ = ApplicationTools::getDoubleParameter("_demiPloidyR", ChromEvol.getParams(), IgnoreParam, "", true, 0);
    maxParsimonyBound_ = ApplicationTools::getBooleanParameter("_maxParsimonyBound", ChromEvol.getParams(), false, "", true, 0);
    OptPointsNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizePointsNum", ChromEvol.getParams(), ',', "10,3,1", "", true, 0);
    OptIterNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizeIterNum", ChromEvol.getParams(), ',', "0,2,5", "", true, 0);

}
/*************************************************************************/
VectorSiteContainer* ChromEvolOptions::getCharacterData (const string& path, unsigned int* numberOfUniqueCharacterStates){
    Fasta fasta;
    ChromosomeAlphabet* alphaInitial = new ChromosomeAlphabet(minAlpha_, maxAlpha_);
    VectorSequenceContainer* initialSetOfSequences = fasta.readSequences(path, alphaInitial);
    size_t numOfSequences = initialSetOfSequences->getNumberOfSequences();
    vector <string> sequenceNames = initialSetOfSequences->getSequencesNames();

    unsigned int maxNumberOfChr = 1; //the minimal number of chromosomes cannot be zero
    unsigned int minNumOfChr = maxAlpha_;

    std::vector <int> UniqueCharacterStates;
    cout<<"vector size is "<< UniqueCharacterStates.size()<<endl;
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = initialSetOfSequences->getSequence(sequenceNames[i]);
        int character = seq.getValue(0);
        if (character == -1){
            continue;
        }

        if (!std::count(UniqueCharacterStates.begin(), UniqueCharacterStates.end(), character)){
            UniqueCharacterStates.push_back(character);

        }
        if (character == static_cast<int>(maxAlpha_)+1){
            continue;
        }
        if ((unsigned int) character > maxNumberOfChr){
            maxNumberOfChr = character;
        }
        if ((unsigned int) character < minNumOfChr){
            minNumOfChr = character;
        }

    }
    *numberOfUniqueCharacterStates = (unsigned int)UniqueCharacterStates.size();
    cout <<"Number of unique states is " << *numberOfUniqueCharacterStates <<endl;

    setMaxChrNum(maxNumberOfChr);
    setMinChrNum(minNumOfChr);

    VectorSiteContainer* vsc = resizeAlphabetForSequenceContainer(initialSetOfSequences);
    delete initialSetOfSequences;
    delete alphaInitial;
    return vsc;
}
/****************************************************************************/
void ChromEvolOptions::setMaxChrNum(unsigned int maxNumberOfChr){
    if (maxChrNum_ < 0){
        maxChrNum_ = maxNumberOfChr + abs(maxChrNum_);
    }else{
        if ((int)maxNumberOfChr > maxChrNum_){
            maxChrNum_ = maxNumberOfChr;
        }
    }

}
/****************************************************************************/
void ChromEvolOptions::setMinChrNum(unsigned int minNumberOfChr){
    if (minChrNum_ < 0){
        minChrNum_ = minNumberOfChr - abs(minChrNum_);
    }else{
        if ((int)minNumberOfChr < minChrNum_){
            minChrNum_ = minNumberOfChr;
        }

    }
}
/****************************************************************************/
VectorSiteContainer* ChromEvolOptions::resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc){
    size_t numOfSequences = vsc->getNumberOfSequences();
    vector <string> sequenceNames = vsc->getSequencesNames();
    ChromosomeAlphabet* new_alphabet = new ChromosomeAlphabet(minChrNum_,maxChrNum_);
    VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(new_alphabet);
    for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = vsc->getSequence(sequenceNames[i]);
        BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), new_alphabet);
        resized_alphabet_site_container->addSequence(new_seq);

    }
    return resized_alphabet_site_container;
}
/****************************************************************************/
TreeTemplate<Node>* ChromEvolOptions::getTree(const string& path, unsigned int NumOfUniqueCharacterStates){
    Newick newick;
    TreeTemplate<Node>* tree = newick.readTree(path);
    rescale_tree(tree, NumOfUniqueCharacterStates);
    return tree;

}
/****************************************************************************/
void ChromEvolOptions::rescale_tree(TreeTemplate<Node>* tree, unsigned int chrRange){
    double scale_tree_factor = 1.0;
    string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    std :: cout << tree_str << endl;
    bool rooted = tree->isRooted();
    if (!rooted){
        throw UnrootedTreeException("The given input tree is unrooted. Tree must be rooted!", tree);
    }
    if (branchMul_ == 1.0){
        return;
    }else{
        //tree must be rescaled
        double treeLength = tree->getTotalLength();
        if (branchMul_ == 999){
            scale_tree_factor = (double)chrRange/treeLength;
        }
        tree->scaleTree(scale_tree_factor);

    }
    string tree_str_final = TreeTemplateTools::treeToParenthesis(*tree);
    std :: cout << tree_str_final << endl;
}
