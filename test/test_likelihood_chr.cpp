//from bpp-core
#include <Bpp/Version.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
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
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

//from bpp-phyl
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
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


//Functions for initialization of Data
VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc);
VectorSiteContainer* getCharacterData(const std :: string &path, unsigned int* numberOfUniqueStates);
void setMaxChrNum(unsigned int maxNumberOfChr);
void setMinChrNum(unsigned int minNumberOfChr);
void rescale_tree(TreeTemplate<Node>* tree, unsigned int chrRange);
TreeTemplate<Node>* getTree(const std :: string &path, unsigned int numOfUniqueCharacterStates);

//Likelihood Optimization Functions
void testInitialLikelihood(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree);
void testOptimizeLikelihood(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree);
void OptimizeMultiStartingPoints(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree);
//Axillary functions for likelihood optimization
bool compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2);
void printLikParameters(DRNonHomogeneousTreeLikelihood &lik, unsigned int optimized);
void printRootFrequencies(DRNonHomogeneousTreeLikelihood &lik);
void printLikelihoodVectorValues(std::vector <DRNonHomogeneousTreeLikelihood> lik_vec);
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood &lik);
void clearVectorOfLikelihoods(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, size_t new_size);
ChromosomeSubstitutionModel* initModel(const ChromosomeAlphabet* alpha);
ChromosomeSubstitutionModel* initRandomModel(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int pointNum);
void initLikelihoodVector(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec);

/******************************************************************************/
//Delete items from the vector of likelihood until reaching the required size new_size
void clearVectorOfLikelihoods(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, size_t new_size){
    while(lik_vec.size() > new_size){
        deleteTreeLikAssociatedAttributes(lik_vec[lik_vec.size()-1]);
        lik_vec.pop_back();        
    }
}

/******************************************************************************/
VectorSiteContainer* getCharacterData (const string& path, unsigned int* numberOfUniqueCharacterStates){
    Fasta fasta;
    ChromosomeAlphabet* alphaInitial = new ChromosomeAlphabet(ChromEvolOptions::minAlpha_, ChromEvolOptions::maxAlpha_);
    VectorSequenceContainer* initialSetOfSequences = fasta.readSequences(path, alphaInitial);
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

        if (!std::count(UniqueCharacterStates.begin(), UniqueCharacterStates.end(), character)){
            UniqueCharacterStates.push_back(character);

        }
        if (character == static_cast<int>(ChromEvolOptions::maxAlpha_)+1){
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
/****************************************************************************/
VectorSiteContainer* resizeAlphabetForSequenceContainer(VectorSequenceContainer* vsc){
    size_t numOfSequences = vsc->getNumberOfSequences();
    vector <string> sequenceNames = vsc->getSequencesNames();
    ChromosomeAlphabet* new_alphabet = new ChromosomeAlphabet(ChromEvolOptions::minChrNum_,ChromEvolOptions::maxChrNum_);
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
    string tree_str = TreeTemplateTools::treeToParenthesis(*tree);
    std :: cout << tree_str << endl;
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
    string tree_str_final = TreeTemplateTools::treeToParenthesis(*tree);
    std :: cout << tree_str_final << endl;
}

/******************************************************************************/
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood &lik){
    const SubstitutionModelSet* modelSet = lik.getSubstitutionModelSet();
    const DiscreteDistribution* rateDist = lik.getRateDistribution();
    delete modelSet;
    delete rateDist;
}
/******************************************************************************/

void printLikelihoodVectorValues(std::vector <DRNonHomogeneousTreeLikelihood> lik_vec){
    std :: cout <<"The likelihoods at the end of cycle are :"<<endl;
    for (size_t i = 0; i < lik_vec.size(); i++){
        std :: cout << lik_vec[i].getValue() << endl;
    }
}

/******************************************************************************/
void printRootFrequencies(DRNonHomogeneousTreeLikelihood &lik){
    std :: vector <double> rootFreq = lik.getRootFrequencies(0);
    for (int i = 0; i < (int)(rootFreq.size()); i++){
        std :: cout << "root freq "<< i << " "<< rootFreq[i] <<endl;
    }

}
/******************************************************************************/

void printLikParameters(DRNonHomogeneousTreeLikelihood &lik, unsigned int optimized){
    //double res = lik.getLikelihood();
    if (optimized == 0){
        std:: cout << "Initial likelihood is : "<< lik.getValue() << endl;
    }else{
        std:: cout << "Optimized likelihood is : "<< lik.getValue() << endl;
    }
    //std:: cout << "likelihood is : " << res << endl;
    
    std:: cout << "Parameters are:" << endl;
    ParameterList params = lik.getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = params.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        std::cout << paramsNames[i] << "value is "<< params.getParameterValue(paramsNames[i])<<endl;
    }
    std::cout <<"***"<<endl;

}
/******************************************************************************/
//For likelihood sorting
bool compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2){
    return (lik1.getValue() < lik2.getValue());
}

/******************************************************************************/
ChromosomeSubstitutionModel* initModel(const ChromosomeAlphabet* alpha){
    double gain = ChromEvolOptions::constGain_;
    double loss = ChromEvolOptions::constLoss_;
    double dupl = ChromEvolOptions::constDupl_;
    double demiDupl = ChromEvolOptions::constDemiDupl_;
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(alpha, gain, loss, dupl,  demiDupl, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL);
    return chrModel;

}
/******************************************************************************/
ChromosomeSubstitutionModel* initRandomModel(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int pointNum){
    double gain, loss, dupl, demiDupl;
    double upperBound = upperBoundOfRateParam;
    if (ChromEvolOptions::maxParsimonyBound_){
        TreeTemplate<Node>* treeForParsimonyBound = tree->clone();
        DRTreeParsimonyScore maxParsimonyObject = DRTreeParsimonyScore(*treeForParsimonyBound, *vsc);
        double parsimonyBound = pointNum * ((maxParsimonyObject.getScore())/(tree->getTotalLength()));
        upperBound = std::min(upperBoundOfRateParam, parsimonyBound);
        delete treeForParsimonyBound;

    }
    if (ChromEvolOptions::constGain_ != IgnoreParam){
        gain = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);    
    }else{
        gain = IgnoreParam;
    }
    if (ChromEvolOptions::constLoss_ != IgnoreParam){
        loss = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
    }else{
        loss = IgnoreParam;
    }
    if (ChromEvolOptions::constDupl_ != IgnoreParam){
        dupl = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
    }else{
        dupl = IgnoreParam;
    }
    if ((ChromEvolOptions::constDemiDupl_ != IgnoreParam) & (ChromEvolOptions::constDemiDupl_ != DemiEqualDupl)){
        demiDupl = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);

    }else{
        demiDupl = ChromEvolOptions::constDemiDupl_;
    }
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(alpha, gain, loss, dupl,  demiDupl, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL);
    return chrModel;

}

/******************************************************************************/
void testInitialLikelihood(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree){
    //const ChromosomeAlphabet* chr_alpha = dynamic_cast <const ChromosomeAlphabet*>(vsc->getAlphabet());
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    SubstitutionModelSet* modelSet = new SubstitutionModelSet(alpha);
    ChromosomeSubstitutionModel* chrModel = initModel(alpha);
    vector <int> nodeIds = tree->getNodesId();
    nodeIds.pop_back();
    modelSet->addModel(chrModel, nodeIds);
    DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, modelSet, rdist);
    lik.initialize();
    printLikParameters(lik, 0);
    printRootFrequencies(lik);
    //delete chrModel;
    delete modelSet;
    delete rdist;
    //delete chr_alpha;

}

/******************************************************************************/
void testOptimizeLikelihood(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree){

    //const ChromosomeAlphabet* chr_alpha = dynamic_cast < const ChromosomeAlphabet*>(vsc->getAlphabet());
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    SubstitutionModelSet* modelSet = new SubstitutionModelSet(alpha);
    ChromosomeSubstitutionModel* chrModel = initModel(alpha);
    // setting the likelihood instance (the root is not attached to any model in the modelset)
    vector <int> nodeIds = tree->getNodesId();
    nodeIds.pop_back();
    modelSet->addModel(chrModel, nodeIds);
    DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, modelSet, rdist);
    //initializing and printing the initial likelihood value and the associated parameters
    lik.initialize();
    printLikParameters(lik, 0);
    printRootFrequencies(lik);

    // optimize likelihood
    ParameterList params = lik.getSubstitutionModelParameters();//the model is already set to include only the unignored parameters
    unsigned int optimization_res = OptimizationTools::optimizeNumericalParameters(&lik, params, 0, 1, ChromEvolOptions::tolerance_, ChromEvolOptions::maxIterations_, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 1, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, 1);
    std::cout <<"optimization iterations : "<< optimization_res<<endl;
    printLikParameters(lik, 1);
    //free pointers
    delete modelSet;
    delete rdist;


}
/******************************************************************************/
//Setting the initial starting points
void initLikelihoodVector(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree){
    //Preparing tree nodes for the creation of likelihood instance (excluding the root)
    vector <int> nodeIds = tree->getNodesId();
    nodeIds.pop_back();
    
    //create starting models
    for (size_t n = 0; n < ChromEvolOptions::OptPointsNum_[0]; n++){
        DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
        SubstitutionModelSet* modelSet = new SubstitutionModelSet(alpha);
        ChromosomeSubstitutionModel* chrModel;
        if (n == 0){
            chrModel = initModel(alpha);
        }else{
            chrModel = initRandomModel(alpha, vsc, tree, static_cast<unsigned int>(n));
        }       
        modelSet->addModel(chrModel, nodeIds);
        DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, modelSet, rdist);

        //initializing the likelihood instance
        lik.initialize();
        lik_vec.push_back(lik);//add to vector of likelihoods

        
    }

}

/******************************************************************************/
void OptimizeMultiStartingPoints(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree){
    //initializing the stating points at cycle 0
    std::vector <DRNonHomogeneousTreeLikelihood> lik_vec;
    lik_vec.reserve(ChromEvolOptions::OptPointsNum_[0]);
    initLikelihoodVector(lik_vec, alpha, vsc, tree);

    //Go over each cycle
    for (size_t i = 0; i < ChromEvolOptions::OptIterNum_.size(); i++){
        clearVectorOfLikelihoods(lik_vec, ChromEvolOptions::OptPointsNum_[i]);
        std::cout <<"##################################" << endl;
        std:: cout << "*********  cycle "<< i <<"  **************"<<endl;     
        //Go over each point at cycle i 
        for (size_t j = 0; j < ChromEvolOptions::OptPointsNum_[i]; j++){
            std::cout << "Starting cycle with Point #" << j <<"...."<<endl;;
            printLikParameters(lik_vec[j], 0);
            //If the number of optimization iterations is larger than zero, optimize the number of times as specified
            double prevLogLik;
            if (ChromEvolOptions::OptIterNum_[i] > 0){
                ParameterList params;
                for (size_t it = 0; it < ChromEvolOptions::OptIterNum_[i]; it++){
                    std::cout << "Iteration #"<<it <<endl;
                    prevLogLik = lik_vec[j].getValue();
                    params = lik_vec[j].getSubstitutionModelParameters();
                    OptimizationTools::optimizeNumericalParameters(&lik_vec[j], params, 0, 1,ChromEvolOptions::tolerance_, 2, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, 1);
                    printLikParameters(lik_vec[j], 1);
                    if (abs(lik_vec[j].getValue() - prevLogLik) < ChromEvolOptions::tolerance_){
                        break;
                    }
                }
                std::cout <<"..."<<endl;                                
            }            
        }
        //sort the vector of likelihoods, such that the worst likelihood is at the end
        sort(lik_vec.begin(), lik_vec.end(), compareLikValues);
        printLikelihoodVectorValues(lik_vec);
        
    }
    printRootFrequencies(lik_vec[0]);
    std::cout <<"*****  Final Optimized -logL *********"  <<endl;
    printLikParameters(lik_vec[0], 1);
    //Clear the vector of likelihoods entirely
    clearVectorOfLikelihoods(lik_vec, 0);
    //std::cout << "Length of vector is: "<< lik_vec.size() << endl;
    //delete chr_alpha;

}
/******************************************************************************/


int main(int args, char **argv) {

    if (args == 1){
        std::cout << "No arguments provided"<<endl;
        return 0;
    }
    try{
        BppApplication ChromEvol(args, argv, "ChromEvol");
        ChromEvolOptions::initAllParameters(ChromEvol);
        unsigned int numberOfUniqueStates = 0;
        VectorSiteContainer* vsc = getCharacterData(ChromEvolOptions::characterFilePath_, &numberOfUniqueStates);
        const ChromosomeAlphabet* alpha = dynamic_cast<const ChromosomeAlphabet*>(vsc->getAlphabet());
        TreeTemplate<Node>* tree = getTree(ChromEvolOptions::treeFilePath_, numberOfUniqueStates);
        OptimizeMultiStartingPoints(alpha, vsc, tree);
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