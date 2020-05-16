//from bpp-core
#include <Bpp/Version.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/App/BppApplication.h>

//from bpp-seq
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Io/AbstractISequence.h>
#include <Bpp/Seq/Io/ISequence.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

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

//Likelihood Optimization Functions
void testInitialLikelihood();
void testOptimizeLikelihood();
void OptimizeMultiStartingPoints();
//Axillary functions for likelihood optimization
bool compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2);
void printLikParameters(DRNonHomogeneousTreeLikelihood &lik, unsigned int optimized);
void printRootFrequencies(DRNonHomogeneousTreeLikelihood &lik);
void printLikelihoodVectorValues(std::vector <DRNonHomogeneousTreeLikelihood> lik_vec);
void deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood &lik);
void clearVectorOfLikelihoods(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, size_t cycle);
ChromosomeSubstitutionModel* initModel();
ChromosomeSubstitutionModel* initRandomModel();
void initLikelihoodVector(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec);

/******************************************************************************/
void clearVectorOfLikelihoods(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec, size_t cycle){
    size_t new_size = 0;
    if ((cycle < ChromEvolOptions::OptPointsNum_.size()) & (cycle != ChromEvolOptions::OptPointsNum_.size()-1)){
        new_size = ChromEvolOptions::OptPointsNum_[cycle + 1];
    }
    if (cycle != ChromEvolOptions::OptPointsNum_.size()-1){
        while(lik_vec.size() > new_size){
            deleteTreeLikAssociatedAttributes(lik_vec[lik_vec.size()-1]);
            lik_vec.pop_back();
        }
    }
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
    std :: cout <<"The likelihoods at rhe end of cycle are :"<<endl;
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
    double res = lik.getLikelihood();
    if (optimized == 0){
        std:: cout << "Initial likelihood is : "<< lik.getValue() << endl;
    }else{
        std:: cout << "Optimized likelihood is : "<< lik.getValue() << endl;
    }
    std:: cout << "likelihood is : " << res << endl;
    
    std:: cout << "Parameters are:" << endl;
    ParameterList params = lik.getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = params.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        std::cout << paramsNames[i] << "value is "<< params.getParameterValue(paramsNames[i])<<endl;
    }
    std::cout <<"***"<<endl;

}
/******************************************************************************/

bool compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2){
    return (lik1.getValue() < lik2.getValue());
}

/******************************************************************************/
ChromosomeSubstitutionModel* initModel(){
    double gain = ChromEvolOptions::constGain_;
    double loss = ChromEvolOptions::constLoss_;
    double dupl = ChromEvolOptions::constDupl_;
    double demiDupl = ChromEvolOptions::constDemiDupl_;
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(ChromEvolOptions::alpha_, gain, loss, dupl,  demiDupl, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL);
    return chrModel;

}
/******************************************************************************/
ChromosomeSubstitutionModel* initRandomModel(){
    double gain, loss, dupl, demiDupl;
    if (ChromEvolOptions::constGain_ != IgnoreParam){
        gain = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBoundOfRateParam);    
    }else{
        gain = IgnoreParam;
    }
    if (ChromEvolOptions::constLoss_ != IgnoreParam){
        loss = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBoundOfRateParam);
    }else{
        loss = IgnoreParam;
    }
    if (ChromEvolOptions::constDupl_ != IgnoreParam){
        dupl = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBoundOfRateParam);
    }else{
        dupl = IgnoreParam;
    }
    if ((ChromEvolOptions::constDemiDupl_ != IgnoreParam) & (ChromEvolOptions::constDemiDupl_ != DemiEqualDupl)){
        demiDupl = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBoundOfRateParam);

    }else{
        demiDupl = ChromEvolOptions::constDemiDupl_;
    }
    ChromosomeSubstitutionModel* chrModel = new ChromosomeSubstitutionModel(ChromEvolOptions::alpha_, gain, loss, dupl,  demiDupl, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL);
    return chrModel;

}

/******************************************************************************/
void testInitialLikelihood(){
    //const ChromosomeAlphabet* chr_alpha = dynamic_cast <const ChromosomeAlphabet*>(vsc->getAlphabet());
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    SubstitutionModelSet* modelSet = new SubstitutionModelSet(ChromEvolOptions::alpha_);
    ChromosomeSubstitutionModel* chrModel = initModel();
    vector <int> nodeIds = ChromEvolOptions::tree_->getNodesId();
    nodeIds.pop_back();
    modelSet->addModel(chrModel, nodeIds);
    DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*ChromEvolOptions::tree_, *ChromEvolOptions::vsc_, modelSet, rdist);
    lik.initialize();
    printLikParameters(lik, 0);
    printRootFrequencies(lik);
    //delete chrModel;
    delete modelSet;
    delete rdist;
    //delete chr_alpha;

}

/******************************************************************************/
void testOptimizeLikelihood(){

    //const ChromosomeAlphabet* chr_alpha = dynamic_cast < const ChromosomeAlphabet*>(vsc->getAlphabet());
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    SubstitutionModelSet* modelSet = new SubstitutionModelSet(ChromEvolOptions::alpha_);
    ChromosomeSubstitutionModel* chrModel = initModel();
    vector <int> nodeIds = ChromEvolOptions::tree_->getNodesId();
    nodeIds.pop_back();
    modelSet->addModel(chrModel, nodeIds);
    DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*ChromEvolOptions::tree_, *ChromEvolOptions::vsc_, modelSet, rdist);
    lik.initialize();
    printLikParameters(lik, 0);
    printRootFrequencies(lik);

    // optimize likelihood
    ParameterList params = lik.getSubstitutionModelParameters();
    unsigned int optimization_res = OptimizationTools::optimizeNumericalParameters(&lik, params, 0, 1, ChromEvolOptions::tolerance_, ChromEvolOptions::maxIterations_, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 1, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, 1);
    std::cout <<"optimization iterations : "<< optimization_res<<endl;
    printLikParameters(lik, 1);
    delete modelSet;
    delete rdist;
    //delete chr_alpha;

}
/******************************************************************************/
void initLikelihoodVector(std::vector <DRNonHomogeneousTreeLikelihood> &lik_vec){
    vector <int> nodeIds = ChromEvolOptions::tree_->getNodesId();
    nodeIds.pop_back();
    
    //create starting models
    for (size_t n = 0; n < ChromEvolOptions::OptPointsNum_[0]; n++){
        std::cout <<"##################################" << endl;
        std:: cout << "*********  cycle 0  **************"<<endl;
        DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
        SubstitutionModelSet* modelSet = new SubstitutionModelSet(ChromEvolOptions::alpha_);
        ChromosomeSubstitutionModel* chrModel;
        if (n == 0){
            chrModel = initModel();
        }else{
            chrModel = initRandomModel();
        }       
        modelSet->addModel(chrModel, nodeIds);
        DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*ChromEvolOptions::tree_, *ChromEvolOptions::vsc_, modelSet, rdist);

        lik.initialize();
        lik_vec.push_back(lik);
        printLikParameters(lik, 0);
        
    }

}

/******************************************************************************/
void OptimizeMultiStartingPoints(){
    std::vector <DRNonHomogeneousTreeLikelihood> lik_vec;
    lik_vec.reserve(ChromEvolOptions::OptPointsNum_[0]);
    initLikelihoodVector(lik_vec);
 
    for (size_t i = 0; i < ChromEvolOptions::OptIterNum_.size(); i++){ 
        std::cout <<"##################################" << endl;
        std:: cout << "*********  cycle "<< i <<"  **************"<<endl;     
      
        for (size_t j = 0; j < ChromEvolOptions::OptPointsNum_[i]; j++){
            std::cout << "Optimizing Point #" << j <<"...."<<endl;;
            printLikParameters(lik_vec[j], 0);
            if (ChromEvolOptions::OptIterNum_[i] > 0){
                ParameterList params;
                for (size_t it = 0; it < ChromEvolOptions::OptIterNum_[i]; it++){
                    std::cout << "Iteration #"<<it <<endl;
                    params = lik_vec[j].getSubstitutionModelParameters();
                    OptimizationTools::optimizeNumericalParameters(&lik_vec[j], params, 0, 1,ChromEvolOptions::tolerance_, 2, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, 1);
                    printLikParameters(lik_vec[j], 1);
                }
                std::cout <<"..."<<endl;                                
            }            
        }
        sort(lik_vec.begin(), lik_vec.end(), compareLikValues);
        printLikelihoodVectorValues(lik_vec);
        clearVectorOfLikelihoods(lik_vec, i);
        
    }
    printRootFrequencies(lik_vec[0]);
    std::cout <<"*****  Final Optimized -logL *********"  <<endl;
    printLikParameters(lik_vec[0], 1);
    clearVectorOfLikelihoods(lik_vec, ChromEvolOptions::OptPointsNum_.size());
    std::cout << "Length of vector is: "<< lik_vec.size() << endl;
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
        OptimizeMultiStartingPoints();

    }
    catch (exception& e)
    {
        cout << e.what() << endl;
        return 1;
    }


    return 0;
}