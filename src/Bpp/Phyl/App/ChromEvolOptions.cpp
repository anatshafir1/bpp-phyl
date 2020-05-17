
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
void ChromEvolOptions::initDefaultParameters(){
    maxAlpha_ = 500;
    minAlpha_ = 1;
    maxChrNum_ = -10;
    minChrNum_ = 1;
    maxIterations_ = 5;
    tolerance_ = 0.01;
    branchMul_ = 999;
    constGain_ = IgnoreParam;
    constLoss_ = IgnoreParam;
    constDupl_ = IgnoreParam;
    constDemiDupl_= IgnoreParam;
    maxParsimonyBound_ = false;

}
/*************************************************************************/
void ChromEvolOptions::initParametersFromFile(BppApplication& ChromEvol){
    maxChrNum_ = ApplicationTools::getIntParameter("_maxChrNum", ChromEvol.getParams(), maxChrNum_, "", true, 0);
    minChrNum_ = ApplicationTools::getIntParameter("_minChrNum", ChromEvol.getParams(), minChrNum_, "", true, 0);
    characterFilePath_ = ApplicationTools::getAFilePath("_dataFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    treeFilePath_ = ApplicationTools::getAFilePath("_treeFile", ChromEvol.getParams(), true, true, "", true, "none", 1);
    //unsigned int numberOfUniqueCharacterStates = 0;
    //vsc_ = getCharacterData(pathForCharacterData, &numberOfUniqueCharacterStates);
    //alpha_ = dynamic_cast<const ChromosomeAlphabet*>(vsc_->getAlphabet());
    branchMul_ = ApplicationTools::getDoubleParameter("_branchMul", ChromEvol.getParams(), branchMul_, "", true, 0);
    //tree_ = getTree(pathForTree, numberOfUniqueCharacterStates);
    maxIterations_ = (unsigned int)ApplicationTools::getIntParameter("_maxOptimizationItarations", ChromEvol.getParams(), maxIterations_, "", true, 0);
    tolerance_ = ApplicationTools::getDoubleParameter("_tolParamOptimization", ChromEvol.getParams(), tolerance_, "", true, 0);
    constGain_ = ApplicationTools::getDoubleParameter("_gainConstR", ChromEvol.getParams(), constGain_, "", true, 0);
    constLoss_ = ApplicationTools::getDoubleParameter("_lossConstR", ChromEvol.getParams(), constLoss_, "", true, 0);
    constDupl_ = ApplicationTools::getDoubleParameter("_duplConstR", ChromEvol.getParams(), constDupl_, "", true, 0);
    constDemiDupl_ = ApplicationTools::getDoubleParameter("_demiPloidyR", ChromEvol.getParams(), constDemiDupl_, "", true, 0);
    maxParsimonyBound_ = ApplicationTools::getBooleanParameter("_maxParsimonyBound", ChromEvol.getParams(), maxParsimonyBound_, "", true, 0);
    string defaultValForOptPointsNum = "10,3,1";
    string defaultValForOptIterNum = "0,2,5";
    OptPointsNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizePointsNum", ChromEvol.getParams(), ',', defaultValForOptPointsNum, "", true, 0);
    OptIterNum_ = ApplicationTools::getVectorParameter<unsigned int>("_optimizeIterNum", ChromEvol.getParams(), ',', defaultValForOptIterNum, "", true, 0);

}
/*************************************************************************/

