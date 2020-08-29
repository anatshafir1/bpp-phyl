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
#include <Bpp/Numeric/Function/BrentOneDimension.h>
#include <Bpp/Numeric/Function/ConjugateGradientMultiDimensions.h>
#include <Bpp/Numeric/Function/AbstractNumericalDerivative.h>
#include <Bpp/Numeric/Function/TwoPointsNumericalDerivative.h>

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
void setFixedRootFrequencies(const std::string &path, SubstitutionModelSet* modelSet);
ChromosomeSubstitutionModel* initModel(const ChromosomeAlphabet* alpha, unsigned int chrRange);
ChromosomeSubstitutionModel* initRandomModel(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int pointNum, unsigned int chrRange);

//core functions of ChromEvol
void runChromEvol(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange);
void optimizeLikelihoodMultiStartPoints(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange);
void getJointMLAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik);
std::map<int, std::map<size_t, VVdouble>> getMarginalAncestralReconstruction(DRNonHomogeneousTreeLikelihood* lik);
void computeExpectations(DRNonHomogeneousTreeLikelihood* lik, std::map<int, std::map<size_t, VVdouble>>& jointProbabilitiesFatherSon, int numOfSimulations);

// functions used for likelihood optimization
unsigned int optimizeModelParameters(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector<unsigned int> &baseNumCandidates);//, unsigned int inwardBracketing, bool standardOptimization);
unsigned int optimizeModelParametersOneDimension(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector<unsigned int> &baseNumCandidates, bool mixed = false, unsigned int currentIterNum = 0);
unsigned int optimizeMultiDimensions(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, bool mixed = false, unsigned int currentIterNum = 0);
unsigned int useMixedOptimizers(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates);
DRNonHomogeneousTreeLikelihood getLikelihoodFunction(TreeTemplate<Node>* tree, VectorSiteContainer* vsc, SubstitutionModelSet* modelSet, DiscreteDistribution* rdist);
void optimizeBaseNum(ParameterList &params, size_t j, DRNonHomogeneousTreeLikelihood* tl, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, unsigned int* numOfBaseNumEval, unsigned int* numOfLikEvaluations, double lowerBound, double upperBound);
std::vector <string> getNonFixedParams(std::vector <unsigned int> fixedParams, ParameterList &allParams);
void fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, VectorSiteContainer* vsc, double lowerBound, double upperBound);
void getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates, VectorSiteContainer* vsc);
void setNewBounds(ParameterList params, Parameter &param, double* lowerBound);
bool compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2);
void printLikParameters(DRNonHomogeneousTreeLikelihood &lik, unsigned int optimized);
void printRootFrequencies(DRNonHomogeneousTreeLikelihood &lik);
void printLikelihoodVectorValues(std::vector <DRNonHomogeneousTreeLikelihood> lik_vec);
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood &lik);
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood* lik);
void clearVectorOfLikelihoods(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, size_t new_size);
void initLikelihoodVector(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, unsigned int chrRange);

// functions to print the tree with ancestral reconstruction
void printTreeWithStates(DRNonHomogeneousTreeLikelihood* lik, TreeTemplate<Node> tree, std::map<int, std::vector<size_t> > ancestors, std::map<int, map<size_t, std::vector<double>>>* probs = 0);
string printTree(const TreeTemplate<Node>& tree);
string nodeToParenthesis(const Node& node);



/******************************************************************************/
//Delete items from the vector of likelihood until reaching the required size new_size
void clearVectorOfLikelihoods(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, size_t new_size){
    while(lik_vec.size() > new_size){
        deleteTreeLikAssociatedAttributes(lik_vec[lik_vec.size()-1]);
        lik_vec.pop_back();        
    }
}

/******************************************************************************/
DRNonHomogeneousTreeLikelihood getLikelihoodFunction(TreeTemplate<Node>* tree, VectorSiteContainer* vsc, SubstitutionModelSet* modelSet, DiscreteDistribution* rdist){
    bool calculateDerivatives = true;
    if (ChromEvolOptions::optimizationMethod_ == "Brent"){
        calculateDerivatives  = false;
    }
    if (ChromEvolOptions::fixedFrequenciesFilePath_ != "none"){
        setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, modelSet);
        DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, false, calculateDerivatives, modelSet, rdist);
        return lik;
    }else{
        DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, true, calculateDerivatives, modelSet, rdist);
        return lik;

    }
}
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
void getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates, VectorSiteContainer* vsc){
    size_t numOfSequences = vsc->getNumberOfSequences();
    unsigned int minRange = 0;
    vector <string> sequenceNames = vsc->getSequencesNames();
    for (size_t i = 0; i < numOfSequences; i++){
        if (i == numOfSequences-1){
            continue;
        }
        BasicSequence seq1 = vsc->getSequence(sequenceNames[i]);
        int chrNum1 = seq1.getValue(0);
        if (chrNum1 == -1){
            continue;
        }
        for (size_t j = i + 1; j < numOfSequences; j++){
            BasicSequence seq2 = vsc->getSequence(sequenceNames[j]);
            int chrNum2 = seq2.getValue(0);
            if (chrNum2 == -1){
                continue;
            }
            unsigned int chrRange = (unsigned int)(abs(chrNum1 - chrNum2));
            if (chrRange == 0 || chrRange == 1){
                continue;
            }
            else if (chrRange == 2){
                continue;
            }
            if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
                if (minRange == 0){
                    minRange = chrRange;
                }else{
                    if (chrRange < minRange){
                        minRange = chrRange;
                    }
                }
                baseNumCandidates.push_back(chrRange);

            }

        }
    }
    if (minRange > 3){
        for (unsigned int i = 3; i < minRange; i++){
            baseNumCandidates.push_back(i);
        }

    }

}
/****************************************************************************/
std::vector <string> getNonFixedParams(std::vector <unsigned int> fixedParams, ParameterList &allParams){
    std:: vector <string> paramsNames;
    for (size_t i = 0; i < allParams.size(); i++){
        if (!(fixedParams[i])){
            paramsNames.push_back(allParams[i].getName());
        }

    }
    return paramsNames;
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
void setNewBounds(ParameterList params, Parameter &param, double* lowerBound){
    if (param.getName() == "Chromosome.baseNum_1"){
        return;
    }
    if (ChromEvolOptions::rateChangeType_ != ChromosomeSubstitutionModel::LINEAR){
        if (*lowerBound == 0){
            *lowerBound = lowerBoundOfRateParam;
        }
    }else{
        string paramName = param.getName();
        if (paramName == "Chromosome.gainR_1" ){
            if (ChromEvolOptions::constGain_ != IgnoreParam){
                double gain = params.getParameterValue("Chromosome.gain_1");
                *lowerBound = -gain/(ChromEvolOptions::maxChrNum_-1); 
            }else{
                *lowerBound = lowerBoundOfRateParam; 
            }
           
        }
        else if (paramName == "Chromosome.lossR_1"){
            if (ChromEvolOptions::constLoss_ != IgnoreParam){
                double loss = params.getParameterValue("Chromosome.loss_1");
                *lowerBound = -loss/(ChromEvolOptions::maxChrNum_-1);
            }else{
                *lowerBound = lowerBoundOfRateParam;
            }


        }else if (paramName == "Chromosome.duplR_1"){
            if (ChromEvolOptions::constDupl_ != IgnoreParam){
                double dupl = params.getParameterValue("Chromosome.dupl_1");
                *lowerBound = -dupl/(ChromEvolOptions::maxChrNum_-1);

            }else{
                *lowerBound = lowerBoundOfRateParam;
            }

       
        }else{
            if (paramName != "Chromosome.baseNumR_1"){
                if ((paramName == "Chromosome.gain_1") && (ChromEvolOptions::gainR_ != IgnoreParam)){
                    double gainR = params.getParameterValue("Chromosome.gainR_1");
                    *lowerBound = std::max(lowerBoundOfRateParam, -gainR * (ChromEvolOptions::maxChrNum_-1));

                }else if ((paramName == "Chromosome.loss_1")&&(ChromEvolOptions::lossR_ != IgnoreParam)){
                    double lossR = params.getParameterValue("Chromosome.lossR_1");
                    *lowerBound = std::max(lowerBoundOfRateParam, -lossR * (ChromEvolOptions::maxChrNum_-1));

                }else if ((paramName == "Chromosome.dupl_1")&&(ChromEvolOptions::duplR_ != IgnoreParam)){
                    double duplR = params.getParameterValue("Chromosome.duplR_1");
                    *lowerBound = std::max(lowerBoundOfRateParam, -duplR * (ChromEvolOptions::maxChrNum_-1));
                
                }else{
                    if (*lowerBound == 0){
                        *lowerBound = lowerBoundOfRateParam;
                    }
                }
            }
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
    //string tree_str_final = TreeTemplateTools::treeToParenthesis(*tree);
    //std :: cout << tree_str_final << endl;
}

/******************************************************************************/
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood &lik){
    const SubstitutionModelSet* modelSet = lik.getSubstitutionModelSet();
    const DiscreteDistribution* rateDist = lik.getRateDistribution();
    delete modelSet;
    delete rateDist;
}
/******************************************************************************/
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood* lik){
    const SubstitutionModelSet* modelSet = lik->getSubstitutionModelSet();
    const DiscreteDistribution* rateDist = lik->getRateDistribution();
    delete modelSet;
    delete rateDist;
    delete lik;
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
    //unsigned int numOfLikEvaluations = lik.getNumberOfLikelihoodEvaluations();
    //std:: cout << "number of likelihood evaluations is : " << numOfLikEvaluations << endl;
    
    std:: cout << "Parameters are:" << endl;
    ParameterList params = lik.getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = params.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i] == "Chromosome.baseNum_1"){
            std::cout << paramsNames[i] << "value is "<< (int)(params.getParameterValue(paramsNames[i]))<<endl;
        }else{
            std::cout << paramsNames[i] << "value is "<< params.getParameterValue(paramsNames[i])<<endl;
        }
        
    }
    std::cout <<"***"<<endl;

}
/******************************************************************************/
//For likelihood sorting
bool compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2){
    return (lik1.getValue() < lik2.getValue());
}

/******************************************************************************/
ChromosomeSubstitutionModel* initModel(const ChromosomeAlphabet* alpha, unsigned int chrRange){
    double gain = ChromEvolOptions::constGain_;
    double loss = ChromEvolOptions::constLoss_;
    double dupl = ChromEvolOptions::constDupl_;
    double demiDupl = ChromEvolOptions::constDemiDupl_;
    double gainR = ChromEvolOptions::gainR_;
    double lossR = ChromEvolOptions::lossR_;
    int baseNum = ChromEvolOptions::baseNum_;
    double baseNumR = ChromEvolOptions::baseNumR_;
    double duplR = ChromEvolOptions::duplR_;
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(alpha, gain, loss, dupl, demiDupl, gainR, lossR, baseNum, baseNumR, duplR, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, ChromEvolOptions::optimizeBaseNumber_);
    return chrModel;

}
/******************************************************************************/
ChromosomeSubstitutionModel* initRandomModel(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int pointNum, unsigned int chrRange){
    double gain, loss, dupl, demiDupl, gainR, lossR, baseNumR, duplR;
    int baseNum;
    double upperBound = upperBoundOfRateParam;
    double upperBoundLinear = upperBoundLinearRateParam;
    if (ChromEvolOptions::maxParsimonyBound_){
        TreeTemplate<Node>* treeForParsimonyBound = tree->clone();
        DRTreeParsimonyScore maxParsimonyObject = DRTreeParsimonyScore(*treeForParsimonyBound, *vsc);
        double parsimonyBound = pointNum * ((maxParsimonyObject.getScore())/(tree->getTotalLength()));
        upperBound = std::min(upperBoundOfRateParam, parsimonyBound);
        upperBoundLinear = std::min(upperBoundLinearRateParam, parsimonyBound);
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
    //gainR
    if ((ChromEvolOptions::gainR_ != IgnoreParam)){
        if (ChromEvolOptions::constGain_ != IgnoreParam){
            if (ChromEvolOptions::rateChangeType_ == ChromosomeSubstitutionModel::LINEAR){
                double lowerBoundForGainR = -gain/(ChromEvolOptions::maxChrNum_-1);
                gainR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundForGainR, upperBoundLinear);
            }else{
                gainR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfExpParam, upperBound);

            }
        }else{
            gainR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
        }

    }else{
        gainR = ChromEvolOptions::gainR_;
    }
    //lossR
    if ((ChromEvolOptions::lossR_ != IgnoreParam)){
        if (ChromEvolOptions::constLoss_ != IgnoreParam){
            if (ChromEvolOptions::rateChangeType_ == ChromosomeSubstitutionModel::LINEAR){
                double lowerBoundForlossR = -loss/(ChromEvolOptions::maxChrNum_-1);
                lossR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundForlossR, upperBoundLinear);
            }else{
                lossR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfExpParam, upperBound);
            }
        }else{
            lossR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
        }
    }else{
        lossR = ChromEvolOptions::lossR_;
    }
    //////
    //duplR
    if ((ChromEvolOptions::duplR_ != IgnoreParam)){
        if (ChromEvolOptions::constDupl_ != IgnoreParam){
            if (ChromEvolOptions::rateChangeType_ == ChromosomeSubstitutionModel::LINEAR){
                double lowerBoundForDuplR = -dupl/(ChromEvolOptions::maxChrNum_-1);
                duplR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundForDuplR, upperBoundLinear);
            }else{
                duplR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfExpParam, upperBound);
            }
        }else{
            duplR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
        }

    }else{
        duplR = ChromEvolOptions::duplR_;
    }
    //////
    if ((ChromEvolOptions::baseNum_ != IgnoreParam) && (ChromEvolOptions::optimizeBaseNumber_)){
        int lowerBoundBaseNum = lowerBoundBaseNumber;
        int upperBoundBaseNum = std::max((int)chrRange, lowerBoundBaseNumber+1);
        baseNum = lowerBoundBaseNum + RandomTools::giveIntRandomNumberBetweenZeroAndEntry(upperBoundBaseNum-lowerBoundBaseNum);
    }else{
        baseNum = ChromEvolOptions::baseNum_;
    }
    if (ChromEvolOptions::baseNumR_ != IgnoreParam){
        baseNumR = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
    }else{
        baseNumR = ChromEvolOptions::baseNumR_;
    }
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(alpha, gain, loss, dupl, demiDupl, gainR, lossR, baseNum, baseNumR, duplR, chrRange, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, ChromEvolOptions::optimizeBaseNumber_);
    return chrModel;

}
/******************************************************************************/
unsigned int optimizeModelParameters(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates){
    unsigned int numOfEvaluations = 0;
    if (ChromEvolOptions::standardOptimization_){
        double prevLogLik;
        ParameterList params;
        for (size_t i = 0; i < maxNumOfIterations; i++){
            std::cout << "Iteration #"<<i <<endl;
            prevLogLik = tl->getValue();
            params = tl->getSubstitutionModelParameters();
            numOfEvaluations += OptimizationTools::optimizeNumericalParameters(tl, params, 0, 1, tol, 2, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, (unsigned int)(ChromEvolOptions::BrentBracketing_));
            printLikParameters(*tl, 1);
            if (abs(tl->getValue() - prevLogLik) < tol){
                break;
            }
        }
        std::cout <<"..."<<endl;
    }else{
        if (ChromEvolOptions::optimizationMethod_ == "Brent"){
            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, maxNumOfIterations, baseNumCandidates);
        }else if (ChromEvolOptions::optimizationMethod_ == "gradient"){

            numOfEvaluations += optimizeMultiDimensions(tl, tol, maxNumOfIterations);

        }else{
            numOfEvaluations += useMixedOptimizers(tl, tol, maxNumOfIterations, baseNumCandidates);
        }
        
    }
    
    return numOfEvaluations;
    

}
/*******************************************************************************/
void fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, VectorSiteContainer* vsc, double lowerBound, double upperBound){
    if (ChromEvolOptions::baseNumOptimizationMethod_ == "Ranges"){
        getAllPossibleChrRanges(baseNumCandidates, vsc);

    }
    else if ((ChromEvolOptions::baseNumOptimizationMethod_ == "Sequential") || (baseNumCandidates.size() == 0)){

        for (unsigned int chr = (unsigned int)lowerBound; chr <= (unsigned int)upperBound; chr++){
            baseNumCandidates.push_back(chr);
        }

    }

}
/*******************************************************************************/
void optimizeBaseNum(ParameterList &params, size_t j, DRNonHomogeneousTreeLikelihood* tl, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, unsigned int* numOfBaseNumEval, unsigned int* numOfLikEvaluations, double lowerBound, double upperBound){

    *numOfBaseNumEval += (unsigned int)baseNumCandidates.size();
    Function* func = tl;
    size_t best_i = (size_t)(params[j].getValue());
    //double f_value = func->f(params);
    double f_value = *currentLikelihood;
    for (size_t i = 0; i < baseNumCandidates.size(); i++){
        unsigned int baseNum = baseNumCandidates[i];
        params[j].setValue((double)baseNum);
        double f_i = func->f(params);
        if (f_i < f_value){
            best_i = baseNum;
            f_value = f_i;
        }
    }
    params[j].setValue((double)best_i);
    func->f(params);
    *currentLikelihood = f_value;

}
/******************************************************************************/
unsigned int optimizeModelParametersOneDimension(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, bool mixed, unsigned curentIterNum){

    // Initialize optimizer
    BrentOneDimension* optimizer = new BrentOneDimension(tl);
    optimizer->setVerbose(1);
    //optimizer->setProfiler(ApplicationTools::message.get());
    //optimizer->setMessageHandler(ApplicationTools::message.get());
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    std::cout <<"max chromosome number: " << ChromEvolOptions::maxChrNum_ << endl;
    if (ChromEvolOptions::BrentBracketing_ == 1){
        optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);

    }else if (ChromEvolOptions::BrentBracketing_ == 2){
        optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);
    }else{
        optimizer->setBracketing(BrentOneDimension::BRACKET_OUTWARD);
    }
    
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    unsigned int numOfEvaluations = 0;
    unsigned int numOfBaseNumEval = 0;

    for (size_t i = 0; i < maxNumOfIterations; i++){
        if (mixed){
            std::cout << "Iteration #"<<curentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<<i <<endl;
        }
        //ParameterList params = tl->getSubstitutionModelParameters();
        size_t nbParams = (tl->getSubstitutionModelParameters()).size();
        prevLikelihood = currentLikelihood;
        unsigned int numOfLikEvaluations;
        
        for (size_t j = 0; j < nbParams; j ++){
            numOfLikEvaluations = tl->getNumberOfLikelihoodEvaluations();
            ParameterList params = tl->getSubstitutionModelParameters();
            if (ChromEvolOptions::fixedParams_[j]){
                continue;
            }
            const string nameOfParam = params[j].getName();
            std::cout << "Parameter name is: "<< nameOfParam <<endl;

            std::shared_ptr<Constraint> parameterConstraints = params[j].getConstraint();
            std::shared_ptr<IntervalConstraint> parameterIntervals = dynamic_pointer_cast<IntervalConstraint>(parameterConstraints);
            double lowerBound = parameterIntervals->getLowerBound();
            double upperBound = parameterIntervals->getUpperBound();
            std::cout << "old lower bound: "<< lowerBound <<endl;

            //std::cout << "Parameter lowerBound is: " << lowerBound<< endl;

            //need to set bounds- the bounds of the parameters of the model are updated in the model
            //itself but not in the likelihood function, therefore I need to set the bounds here.
            //It also means that when a linear model is used, I can use only Brent.
            setNewBounds(params, params[j], &lowerBound);
            std::cout <<"new lower bound: " << lowerBound << endl;
            std::cout <<"parameter value before optimization is "<< params[j].getValue() << endl;
            //std::cout << "Parameter new lowerBound is: " << lowerBound<< endl;
            //std::cout << "Parameter upperBound is: " << upperBound <<endl;
            //std::cout << "Parameter value before optimization: "<< params[j].getValue() <<endl;
            if ((ChromEvolOptions::baseNumOptimizationMethod_ != "Brent") && (params[j].getName() == "Chromosome.baseNum_1")){
                if (ChromEvolOptions::optimizeBaseNumber_){
                    optimizeBaseNum(params, j, tl, baseNumCandidates, &currentLikelihood, &numOfBaseNumEval, &numOfLikEvaluations, lowerBound, upperBound);
                    std::cout << "parameter value after optimization "<< params[j].getValue() << endl;
                    continue;

                }
            }
            
            if ((i == 1) & (maxNumOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(tol* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(tol);
            }
            if (params[j].getName() != "Chromosome.baseNum_1"){
                optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
            }else{
                optimizer->setInitialInterval(lowerBound, upperBound);
            }
 
            
            optimizer->init(params.createSubList(j));
            currentLikelihood = optimizer->optimize();
            std::cout <<"Parameter value after optimization: "<< tl->getSubstitutionModelParameters()[j].getValue()<<endl;
            std::cout << "***"<<endl;
            numOfLikEvaluations = tl->getNumberOfLikelihoodEvaluations() - numOfLikEvaluations;
            if (params[j].getName() == "Chromosome.baseNum_1"){
                numOfBaseNumEval += numOfLikEvaluations;
            }
            
            
        }
        printLikParameters(*tl, 1);
        
        if (abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();
       
    }
    std::cout << "Number of likelihood evaluations per parameter is "<< numOfBaseNumEval << endl;
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    //std:: cout << "final number of evaluations is : " << numOfEvaluations << endl;
    delete optimizer;
    return numOfEvaluations;
}

/*******************************************************************************/
unsigned int optimizeMultiDimensions(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, bool mixed, unsigned int currentIterNum){
    DerivableSecondOrder* f = tl;
    unique_ptr<AbstractNumericalDerivative> fnum;
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0000001);
    ConjugateGradientMultiDimensions* optimizer = new ConjugateGradientMultiDimensions(fnum.get());
    ParameterList tmp = tl->getNonDerivableParameters();
    fnum->setParametersToDerivate(tmp.getParameterNames());
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->getStopCondition()->setTolerance(tol* 0.1);
    optimizer->setMaximumNumberOfEvaluations(1000);

    unsigned int numOfEvaluations = 0;
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        if(mixed){
            std::cout << "Iteration #"<< currentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<< i <<endl;
        }
        
        ParameterList paramsFull = tl->getSubstitutionModelParameters();
        std::vector <string> paramsNames = getNonFixedParams(ChromEvolOptions::fixedParams_, paramsFull);
        ParameterList params = paramsFull.createSubList(paramsNames);
        std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBoundOfRateParam + 0.0000000001, upperBoundOfRateParam, true, true);
        for (size_t j = 0; j < params.size(); j++){
            params[j].setConstraint(interval);
        }
        prevLikelihood = currentLikelihood;
        optimizer->init(params);
        currentLikelihood = optimizer->optimize();
        printLikParameters(*tl, 1);
        if (abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        
        
    }
    numOfEvaluations += optimizer->getNumberOfEvaluations();
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    //std::cout << "The final number of evaluations is: "<< numOfEvaluations << endl;
    delete optimizer;
    return numOfEvaluations;

}
/*******************************************************************************/
//Setting the initial starting points
void initLikelihoodVector(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange){
    //Preparing tree nodes for the creation of likelihood instance (excluding the root)
    vector <int> nodeIds = tree->getNodesId();
    nodeIds.pop_back();
    if (ChromEvolOptions::seed_ != 0){
        RandomTools::setSeed(static_cast<long>(ChromEvolOptions::seed_));
    }
    
    //create starting models
    for (size_t n = 0; n < ChromEvolOptions::OptPointsNum_[0]; n++){
        DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
        SubstitutionModelSet* modelSet = new SubstitutionModelSet(alpha);
        ChromosomeSubstitutionModel* chrModel;
        if (n == 0){
            chrModel = initModel(alpha, chrRange);
        }else{
            chrModel = initRandomModel(alpha, vsc, tree, static_cast<unsigned int>(n), chrRange);
        }       
        modelSet->addModel(chrModel, nodeIds);
        DRNonHomogeneousTreeLikelihood lik = getLikelihoodFunction(tree, vsc, modelSet, rdist);
        
        lik.initialize();
        if (std::isnan(lik.getValue())){
            std::cout << "value is nan"<<endl;
        }
        
        while (std::isinf(lik.getValue()))
        {
            deleteTreeLikAssociatedAttributes(lik);
            rdist = new GammaDiscreteRateDistribution(1, 1.0);
            modelSet = new SubstitutionModelSet(alpha);
            chrModel = initRandomModel(alpha, vsc, tree, static_cast<unsigned int>(n), chrRange);
            modelSet->addModel(chrModel, nodeIds);
            lik = getLikelihoodFunction(tree, vsc, modelSet, rdist);
            
            lik.initialize();

        }
            
        

        //initializing the likelihood instance
        //lik.initialize();
        lik_vec.push_back(lik);//add to vector of likelihoods
        
    }

}

/******************************************************************************/
void runChromEvol(const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange){
    // initialize the vector of likelihoods
    std::vector <DRNonHomogeneousTreeLikelihood> lik_vec;

    // optimize likelihood
    optimizeLikelihoodMultiStartPoints(lik_vec, alpha, vsc, tree, chrRange);
    // get joint ML ancestral reconstruction
    getJointMLAncestralReconstruction(&lik_vec[0]);
    //get Marginal ML ancestral reconstruction, and with the help of them- calculate expectations of transitions
    std::map<int, std::map<size_t, VVdouble>>  jointProbabilitiesFatherSon = getMarginalAncestralReconstruction(&lik_vec[0]);
    //compute expectations
    computeExpectations(&lik_vec[0], jointProbabilitiesFatherSon, ChromEvolOptions::NumOfSimulations_);
    //Clear the vector of likelihoods entirely
    clearVectorOfLikelihoods(lik_vec, 0);


}
/******************************************************************************/
void optimizeLikelihoodMultiStartPoints(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, const ChromosomeAlphabet* alpha, VectorSiteContainer* vsc, TreeTemplate<Node>* tree, unsigned int chrRange){
    lik_vec.reserve(ChromEvolOptions::OptPointsNum_[0]);
    initLikelihoodVector(lik_vec, alpha, vsc, tree, chrRange);
    unsigned int totalNumOfEvaluations = 0;
    unsigned int numOfEvaluations;
    unsigned int numOfEvaluationsPerCycle;
    std::vector <unsigned int> baseNumCandidates;

    // If base number is one of the parameters
    if ((ChromEvolOptions::baseNumOptimizationMethod_ != "Brent") && (ChromEvolOptions::optimizeBaseNumber_)){
        fillVectorOfBaseNumCandidates(baseNumCandidates, vsc, lowerBoundBaseNumber, chrRange);
        std::cout << "**** Chromosome ranges are: " << endl;
        for (size_t ch = 0; ch < baseNumCandidates.size(); ch++){
            std::cout <<baseNumCandidates[ch] << endl;
        }
        std::cout << "****"<< endl;

    }

    //Go over each cycle
    for (size_t i = 0; i < ChromEvolOptions::OptIterNum_.size(); i++){
        numOfEvaluationsPerCycle = 0;
        clearVectorOfLikelihoods(lik_vec, ChromEvolOptions::OptPointsNum_[i]);
        std::cout <<"##################################" << endl;
        std:: cout << "*********  cycle "<< i <<"  **************"<<endl;     
        //Go over each point at cycle i 
        for (size_t j = 0; j < ChromEvolOptions::OptPointsNum_[i]; j++){
            numOfEvaluations = 0;
            std::cout << "Starting cycle with Point #" << j <<"...."<<endl;;
            printLikParameters(lik_vec[j], 0);
            //If the number of optimization iterations is larger than zero, optimize the number of times as specified
            if (ChromEvolOptions::OptIterNum_[i] > 0){
                numOfEvaluations = optimizeModelParameters(&lik_vec[j], ChromEvolOptions::tolerance_, ChromEvolOptions::OptIterNum_[i], baseNumCandidates);
                          
            }
            std:: cout << "Number of evaluations per point is : " << numOfEvaluations << endl;
            numOfEvaluationsPerCycle += numOfEvaluations;
            std:: cout <<"*****************************" << endl;            
        }
        totalNumOfEvaluations += numOfEvaluationsPerCycle;
        //sort the vector of likelihoods, such that the worst likelihood is at the end
        sort(lik_vec.begin(), lik_vec.end(), compareLikValues);
        printLikelihoodVectorValues(lik_vec);
        
    }
    printRootFrequencies(lik_vec[0]);
    std::cout <<"*****  Final Optimized -logL *********"  <<endl;
    printLikParameters(lik_vec[0], 1);
    std:: cout << "final number of evaluations is : " << totalNumOfEvaluations << endl;

}
/*****************************************************************************/
unsigned int useMixedOptimizers(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates){
    std::vector<size_t> optimization = RandomTools::randMultinomial(maxNumOfIterations, ChromEvolOptions::probsForMixedOptimization_);
    unsigned int numOfEvaluations = 0;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        double prevLikelihood = tl->getValue();
        if (optimization[i] == 0){
            std::cout << "Optimizing with Brent" <<endl;
            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, 1, baseNumCandidates, true, (unsigned int)i);
        }else{
            std::cout << "Optimizing with Gradient Descent" <<endl;
            numOfEvaluations += optimizeMultiDimensions(tl, tol, 1, true, (unsigned int)i);
        }
        double currentLikValue = tl->getValue();
        if (abs(prevLikelihood-currentLikValue) < tol){
            break;
        }


    }
    return numOfEvaluations;

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
    //delete sim;
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
void setFixedRootFrequencies(const std::string &path, SubstitutionModelSet* modelSet){
    ifstream stream;
    stream.open(path.c_str());
    vector <double> freqs;
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    for (size_t i = 0; i < lines.size(); i++){
        string freq_i_str = TextTools::removeSurroundingWhiteSpaces(lines[i]);
        if (freq_i_str == ""){
            continue;
        }
        double freq_i = TextTools::toDouble(freq_i_str);
        freqs.push_back(freq_i);
    }
    size_t nbStates = modelSet->getNumberOfStates();
    if (nbStates != freqs.size()){
        throw "Invalid fixed frequencies file!";
    }
    FixedFrequencySet* rootFreqs = new FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(modelSet->getStateMap(), false)), freqs);
    modelSet->setRootFrequencies(static_cast<FrequencySet*>(rootFreqs));
   

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