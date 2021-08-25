
#include "ChromEvolOptions.h"
using namespace bpp;
using namespace std;

string ChromEvolOptions::treeFilePath_;
string ChromEvolOptions::characterFilePath_;
int ChromEvolOptions::maxChrNum_;
int ChromEvolOptions::minChrNum_;
double ChromEvolOptions::branchMul_;
std::vector <unsigned int> ChromEvolOptions::OptPointsNum_;
std::vector <unsigned int> ChromEvolOptions::OptIterNum_;
vector<double> ChromEvolOptions::gain_;
vector<double> ChromEvolOptions::loss_;
vector<double> ChromEvolOptions::dupl_;
vector<double> ChromEvolOptions::demiDupl_;
int ChromEvolOptions::baseNum_;
vector<double> ChromEvolOptions::baseNumR_;
double ChromEvolOptions::tolerance_;
unsigned int ChromEvolOptions::maxIterations_;
bool ChromEvolOptions::maxParsimonyBound_;
unsigned int ChromEvolOptions::maxAlpha_;
unsigned int ChromEvolOptions::minAlpha_;
int ChromEvolOptions::BrentBracketing_;
bool ChromEvolOptions::standardOptimization_;
string ChromEvolOptions::optimizationMethod_;
int ChromEvolOptions::seed_;
std::vector <double> ChromEvolOptions::probsForMixedOptimization_;
string ChromEvolOptions::rootFreqs_;
string ChromEvolOptions::fixedFrequenciesFilePath_;
std::vector<int> ChromEvolOptions::rateChangeType_;
//bool ChromEvolOptions::optimizeBaseNumber_;
string ChromEvolOptions::baseNumOptimizationMethod_;
std::vector<int> ChromEvolOptions::fixedParams_;
int ChromEvolOptions::NumOfSimulations_;
int ChromEvolOptions::jumpTypeMethod_;
bool ChromEvolOptions::simulateData_;
int ChromEvolOptions::numOfDataToSimulate_;
string ChromEvolOptions::resultsPathDir_;
int ChromEvolOptions::maxBaseNumTransition_;
double ChromEvolOptions::treeLength_;
int ChromEvolOptions::maxNumOfTrials_;
/*************************************************************************/
void ChromEvolOptions::initAllParameters(BppApplication& ChromEvol){
    initDefaultParameters();
    initParametersFromFile(ChromEvol);

}
/*************************************************************************/
void ChromEvolOptions::initDefaultParameters(){
    maxAlpha_ = 500;
    minAlpha_ = 1;
    maxChrNum_ = -10;
    minChrNum_ = 1;
    maxIterations_ = 5;
    tolerance_ = 0.01;
    branchMul_ = 999;
    baseNum_ = IgnoreParam;
    maxParsimonyBound_ = false;
    standardOptimization_ = false;
    BrentBracketing_ = 2;
    optimizationMethod_ = "Brent";
    seed_ = 0;
    rootFreqs_ = "weighted";
    //optimizeBaseNumber_ = false;
    baseNumOptimizationMethod_ = "Brent";
    NumOfSimulations_ = 10000;
    jumpTypeMethod_ = 0;
    simulateData_ = false;
    numOfDataToSimulate_ = 1;
    maxBaseNumTransition_ = 18;
    treeLength_ = 0;
    maxNumOfTrials_ = 100;



}
/*************************************************************************/
std::vector<int> ChromEvolOptions::translateStringParamsToInt(std::vector<string> strParams){
    std::vector <int> params;
    for (size_t i = 0; i < strParams.size(); i++){
        if (strParams[i] == "gain"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::GAIN));   
        }else if (strParams[i] == "loss"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::LOSS));
        }else if (strParams[i] == "dupl"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::DUPL));
        }else if (strParams[i] == "demiPloidyR"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL));
        }else if (strParams[i] == "baseNumR"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::BASENUMR));
        }else if (strParams[i] == "baseNum"){
            params.push_back(static_cast<int>(ChromosomeSubstitutionModel::BASENUM));
        }else{
            throw Exception("ChromEvolOptions::translateStringParamsToInt(): No such parameter!!!");
        }
    }
    return params;
}
/*************************************************************************/
// void ChromEvolOptions::setFixedParams(std::vector<unsigned int> fixedParams){
    
//     //if (optimizeBaseNumber_){
//     for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
//         switch (i)
//         {
//         case ChromosomeSubstitutionModel::BASENUM:
//             if (baseNum_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::BASENUM]);      
//             }
//             break;
//         case ChromosomeSubstitutionModel::BASENUMR:
//             if (baseNumR_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::BASENUMR]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::DUPL:
//             if (constDupl_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::DUPL]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::LOSS:
//             if (constLoss_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::LOSS]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::GAIN:
//             if (constGain_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::GAIN]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::LOSSR:
//             if (lossR_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::LOSSR]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::GAINR:
//             if (gainR_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::GAINR]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::DUPLR:
//             if (duplR_ != IgnoreParam){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::DUPLR]);
//             }
//             break;
//         case ChromosomeSubstitutionModel::DEMIDUPL:
//             if ((constDemiDupl_ != IgnoreParam) && (constDemiDupl_ != DemiEqualDupl)){
//                 fixedParams_.push_back(fixedParams[ChromosomeSubstitutionModel::DEMIDUPL]);
//             }
//             break;
       
//         default:
//             throw Exception("ChromEvolOptions::setFixedParams(): Invalid rate type!");
//             break;
//         }

//     }
// }
/*************************************************************************/
void ChromEvolOptions::initParametersFromFile(BppApplication& ChromEvol){
    maxChrNum_ = ApplicationTools::getIntParameter("_maxChrNum", ChromEvol.getParams(), maxChrNum_, "", true, 0);
    minChrNum_ = ApplicationTools::getIntParameter("_minChrNum", ChromEvol.getParams(), minChrNum_, "", true, 0);
    seed_ = ApplicationTools::getIntParameter("_seed", ChromEvol.getParams(), seed_, "", true, 0);
    simulateData_ = ApplicationTools::getBooleanParameter("_simulateData", ChromEvol.getParams(), simulateData_, "", true, 0);
    if (simulateData_){
        characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), false, true, "", true, "none", 1);
    }else{
        characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    }
    treeFilePath_ = ApplicationTools::getAFilePath("_treeFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    branchMul_ = ApplicationTools::getDoubleParameter("_branchMul", ChromEvol.getParams(), branchMul_, "", true, 0);
    maxIterations_ = (unsigned int)ApplicationTools::getIntParameter("_maxOptimizationItarations", ChromEvol.getParams(), maxIterations_, "", true, 0);
    tolerance_ = ApplicationTools::getDoubleParameter("_tolParamOptimization", ChromEvol.getParams(), tolerance_, "", true, 0);
    gain_ = ApplicationTools::getVectorParameter<double>("_gain", ChromEvol.getParams(), ',', "", "", true, 0);
    loss_ = ApplicationTools::getVectorParameter<double>("_loss", ChromEvol.getParams(), ',', "", "", true, 0);
    dupl_ = ApplicationTools::getVectorParameter<double>("_dupl", ChromEvol.getParams(), ',', "", "", true, 0);
    demiDupl_ = ApplicationTools::getVectorParameter<double>("_demiPloidyR", ChromEvol.getParams(), ',', "", "", true, 0);
    baseNum_ = ApplicationTools::getIntParameter("_baseNum", ChromEvol.getParams(), baseNum_, "", true, 0);
    baseNumR_ = ApplicationTools::getVectorParameter<double>("_baseNumR", ChromEvol.getParams(),',', "",  "", true, 0);
    maxParsimonyBound_ = ApplicationTools::getBooleanParameter("_maxParsimonyBound", ChromEvol.getParams(), maxParsimonyBound_, "", true, 0);
    standardOptimization_ = ApplicationTools::getBooleanParameter("_standardOptimization", ChromEvol.getParams(), standardOptimization_, "", true, 0);
    BrentBracketing_ = ApplicationTools::getIntParameter("_BrentBracketing", ChromEvol.getParams(), BrentBracketing_, "", true, 0);
    optimizationMethod_ = ApplicationTools::getStringParameter("_optimizationMethod", ChromEvol.getParams(), optimizationMethod_, "", true, 0);
    string defaultValForOptPointsNum = "10,3,1";
    string defaultValForOptIterNum = "0,2,5";
    string defaultValForProbsForMixedOpt = "1,0";
    OptPointsNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizePointsNum", ChromEvol.getParams(), ',', defaultValForOptPointsNum, "", true, 0);
    OptIterNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizeIterNum", ChromEvol.getParams(), ',', defaultValForOptIterNum, "", true, 0);
    probsForMixedOptimization_ = ApplicationTools::getVectorParameter<double>("_probsForMixedOptimization", ChromEvol.getParams(), ',', defaultValForProbsForMixedOpt, "", true, 0);
    fixedFrequenciesFilePath_ = ApplicationTools::getAFilePath("_fixedFrequenciesFilePath", ChromEvol.getParams(), false, true, "", true, "none", 0);
    rootFreqs_ = ApplicationTools::getStringParameter("_rootFreqs", ChromEvol.getParams(), rootFreqs_, "", true, 0);
    std::string gainFunc = ApplicationTools::getStringParameter("_gainFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string lossFunc = ApplicationTools::getStringParameter("_lossFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string duplFunc = ApplicationTools::getStringParameter("_duplFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string demiDuplFunc = ApplicationTools::getStringParameter("_demiDuplFunc", ChromEvol.getParams(), "None", "", true, 0);
    std::string baseNumRFunc = ApplicationTools::getStringParameter("_baseNumRFunc", ChromEvol.getParams(), "None", "", true, 0);
    setFunctions(gainFunc, lossFunc, duplFunc, demiDuplFunc, baseNumRFunc);
    //optimizeBaseNumber_ = ApplicationTools::getBooleanParameter("_optimizeBaseNumber", ChromEvol.getParams(), optimizeBaseNumber_, "", true, 0);
    baseNumOptimizationMethod_ = ApplicationTools::getStringParameter("_baseNumOptimizationMethod", ChromEvol.getParams(), baseNumOptimizationMethod_, "", true, 0);
    std::vector<string> fixedParamsStr = ApplicationTools::getVectorParameter<string>("_fixedParams", ChromEvol.getParams(), ',', "", "", true, 0);
    if (fixedParamsStr.size() > 0){
        fixedParams_ = translateStringParamsToInt(fixedParamsStr);
    }else{
        fixedParams_ = std::vector<int>();
    }
    NumOfSimulations_ = ApplicationTools::getIntParameter("_NumOfSimulations", ChromEvol.getParams(), NumOfSimulations_, "", true, 0);
    jumpTypeMethod_ = ApplicationTools::getIntParameter("_jumpTypeMethod", ChromEvol.getParams(), jumpTypeMethod_, "", true, 0);
    numOfDataToSimulate_ = ApplicationTools::getIntParameter("_numOfDataToSimulate", ChromEvol.getParams(), numOfDataToSimulate_, "", true, 0);
    resultsPathDir_ = ApplicationTools::getAFilePath("_resultsPathDir", ChromEvol.getParams(), false, true, "", true, "none", 0);
    maxBaseNumTransition_ = ApplicationTools::getIntParameter("_maxBaseNumTransition", ChromEvol.getParams(), maxBaseNumTransition_, "", true, 0);
    treeLength_ = ApplicationTools::getDoubleParameter("_treeLength", ChromEvol.getParams(), treeLength_, "", true, 0);
    maxNumOfTrials_ = ApplicationTools::getIntParameter("_maxNumOfTrials", ChromEvol.getParams(), maxNumOfTrials_, "", true, 0);

}
/************************************************************************/
void ChromEvolOptions::setFunctions(std::string gainFunc, std::string lossFunc, std::string duplFunc, std::string demiDuplFunc, std::string baseNumRFunc){
    for (size_t i = 0; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
        switch (i)
        {
        case ChromosomeSubstitutionModel::BASENUM:
            break;
        case ChromosomeSubstitutionModel::GAIN:
            rateChangeType_.push_back(getFunctionFromString(gainFunc));
            break;
        case ChromosomeSubstitutionModel::LOSS:
            rateChangeType_.push_back(getFunctionFromString(lossFunc));
            break;
        case ChromosomeSubstitutionModel::DUPL:
            rateChangeType_.push_back(getFunctionFromString(duplFunc));
            break;
        case ChromosomeSubstitutionModel::DEMIDUPL:
            rateChangeType_.push_back(getFunctionFromString(demiDuplFunc));
            break;
        case ChromosomeSubstitutionModel::BASENUMR:
            rateChangeType_.push_back(getFunctionFromString(baseNumRFunc));
            break;
   
        default:
            throw Exception("ChromEvolOptions::setFunctions: parameter not found !!!");
        }
    }

}
/*************************************************************************/
int ChromEvolOptions::getFunctionFromString(string funcStr){
    int func;
    if (funcStr == "CONST"){
        func = static_cast<int>(compositeParameter::CONSTANT);
        
    }else if (funcStr == "LINEAR"){
        func = static_cast<int>(compositeParameter::LINEAR);
    }else if (funcStr == "LINEAR_BD"){
        func = static_cast<int>(compositeParameter::LINEAR_BD);
    }else if (funcStr == "EXP"){
        func = static_cast<int> (compositeParameter::EXP);
    }else if (funcStr == "POLYNOMIAL"){
        func = static_cast<int> (compositeParameter::POLYNOMIAL);
    }else if (funcStr == "LOGNORMAL"){
        func = static_cast<int> (compositeParameter::LOGNORMAL);
    }else if (funcStr == "REVERSE_SIGMOID"){
        func = static_cast<int> (compositeParameter::REVERSE_SIGMOID);
    }else if (funcStr == "IGNORE"){ 
        func = static_cast<int> (compositeParameter::IGNORE);
    }else{
        throw Exception("ChromEvolOptions::getFunctionFromString(): No such function exists!!!");
    }
    return func;
}
/*************************************************************************/
void ChromEvolOptions::getInitialValuesForComplexParams(std::map<int, std::vector<double>> &mapOfParams){
    mapOfParams[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain_;
    mapOfParams[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss_;
    mapOfParams[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl_;
    mapOfParams[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = demiDupl_;
    mapOfParams[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = baseNumR_;

    
}
/*************************************************************************/
// void ChromEvolOptions::initVectorOfChrNumParameters(vector<double>& paramVector){
//     for (size_t i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
//         switch(i){
//             case ChromosomeSubstitutionModel::BASENUM:
//                 paramVector.push_back(baseNum_);
//                 break;
//             case ChromosomeSubstitutionModel::BASENUMR:
//                 paramVector.push_back(baseNumR_);
//                 break;
//             case ChromosomeSubstitutionModel::DUPL:
//                 paramVector.push_back(constDupl_);
//                 break;
//             case ChromosomeSubstitutionModel::LOSS:
//                 paramVector.push_back(constLoss_);
//                 break;
//             case ChromosomeSubstitutionModel::GAIN:
//                 paramVector.push_back(constGain_);
//                 break;
//             case ChromosomeSubstitutionModel::DEMIDUPL:
//                 paramVector.push_back(constDemiDupl_);
//                 break;
//             case ChromosomeSubstitutionModel::LOSSR:
//                 paramVector.push_back(lossR_);
//                 break;
//             case ChromosomeSubstitutionModel::GAINR:
//                 paramVector.push_back(gainR_);
//                 break;
//             case ChromosomeSubstitutionModel::DUPLR:
//                 paramVector.push_back(duplR_);
//                 break;
//             default:
//                 throw Exception("ChromEvolOptions::initVectorOfChrNumParameters(): Invalid rate type!");
//                 break;
//         }

//     }   

// }
