#include <Bpp/Version.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/FrequencySet/NucleotideFrequencySet.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
#include <Bpp/Phyl/Model/ChromosomeSubstitutionModel.h>

#include <Bpp/Phyl/NewLikelihood/ParametrizablePhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/NewLikelihood/RateAcrossSitesSubstitutionProcess.h>

#include <Bpp/Phyl/NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h>
#include <Bpp/Phyl/App/ChromosomeNumberMng.h>
#include <Bpp/Phyl/App/ChromEvolOptions.h>

#include <iostream>

using namespace bpp;
using namespace std;


unsigned int optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, ChromosomeAlphabet* alpha, std::vector<int> &mapOfRateType, size_t maxNumOfIterations);
void updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType, SingleProcessPhyloLikelihood* tl);
void updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName);
void printLikParameters(SingleProcessPhyloLikelihood* lik);

int main(){
    Newick reader;
    //unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.15):0.25,(D:0.35,(E:0.26,F:0.05):0.12):0.16);", false, "", false, false));
    unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((sp1:1,(sp2:0.5,sp3:0.5):0.5):2,(sp4:1.5,(sp5:0.9,sp6:0.9):0.6):1.5);", false, "", false, false));
    ParametrizablePhyloTree parTree(*pTree);
    ChromosomeAlphabet* alpha = new ChromosomeAlphabet(1,4);
    VectorSiteContainer* vsc = new VectorSiteContainer(alpha);
        
    // setting sequence data
    BasicSequence seq1 = BasicSequence("sp1", "2", alpha);
    BasicSequence seq2 = BasicSequence("sp2", "4", alpha);
    BasicSequence seq3 = BasicSequence("sp3", "3", alpha);
    BasicSequence seq4 = BasicSequence("sp4", "4", alpha);
    BasicSequence seq5 = BasicSequence("sp5", "4", alpha);
    BasicSequence seq6 = BasicSequence("sp6", "3", alpha);
    vsc->addSequence(seq1);
    vsc->addSequence(seq2);
    vsc->addSequence(seq3);
    vsc->addSequence(seq4);
    vsc->addSequence(seq5);
    vsc->addSequence(seq6);

    // setting first model
    std::vector<double> gain1;
    gain1.push_back(2);
    std::vector<double> loss1;
    loss1.push_back(1);
    std::vector<double> dupl1;
    dupl1.push_back(3);
    int baseNumber1 = IgnoreParam;


    std::vector<int> rateFuncType;
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::IGNORE);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::CONSTANT);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::CONSTANT);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::CONSTANT);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::IGNORE);

    std::map<int, std::vector<double>> mapOfParamsModel1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = vector<double>();
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = vector<double>();

    std::shared_ptr<ChromosomeSubstitutionModel> chrModel1 = std::make_shared<ChromosomeSubstitutionModel>(alpha, mapOfParamsModel1, baseNumber1, 0, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateFuncType);
    //std::cout << "*** *** *** Model is:" << std::endl;
    //chrModel1->correctBaseNumForSimulation(20);
    // RowMatrix <double> matrix = chrModel1->getGenerator();
    // for (size_t i = 0; i < matrix.getNumberOfRows(); i++){
    //     for (size_t j = 0; j < matrix.getNumberOfColumns(); j++){
    //         std::cout << matrix(i, j) << "\t";

    //     }
    //     std::cout << endl;
    // }

    // setting the second model
    std::vector<double> gain2;
    gain2.push_back(5.31089);
    std::vector<double> loss2;
    loss2.push_back(1.9216);
    std::vector<double> dupl2;
    dupl2.push_back(0.55745);
    int baseNumber2 = IgnoreParam;


    std::map<int, std::vector<double>> mapOfParamsModel2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = vector<double>();
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = vector<double>();
    std::vector<int> rateFuncType2 = rateFuncType;

    std::shared_ptr<ChromosomeSubstitutionModel> chrModel2 = std::make_shared<ChromosomeSubstitutionModel>(alpha, mapOfParamsModel2, baseNumber2, 0, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateFuncType2);

    vector<shared_ptr<PhyloNode> > nodes = pTree->getAllNodes();
    size_t nbNodes = nodes.size();
    vector<shared_ptr<PhyloNode>> subtreeNodes;
    vector<uint> leavesUnderNode;
    vector<string> subtreeLeaves = {"sp2", "sp3"};
    for (size_t i = 0; i < nbNodes; i++){
        uint nodeId = pTree->getNodeIndex(nodes[i]);
        if (nodeId == pTree->getRootIndex()){
            continue;
        }
        if (pTree->isLeaf(nodeId)){
            continue;
        }
        subtreeNodes = pTree->getSubtreeNodes(pTree->getNode(nodeId));
        if (subtreeNodes.size() == 3){
            leavesUnderNode = pTree->getLeavesUnderNode(nodeId);
            vector<string> leavesUnderNodeNames;
            for (size_t k = 0; k < leavesUnderNode.size(); k++){
                string leafName = pTree->getNode(leavesUnderNode[k])->getName();
                leavesUnderNodeNames.push_back(leafName);
            }
            bool notFound = false;
            for (size_t j = 0; j < subtreeLeaves.size(); j++){
                auto it = std::find(leavesUnderNodeNames.begin(), leavesUnderNodeNames.end(), subtreeLeaves[j]);
                if (it == leavesUnderNodeNames.end()){
                    notFound = true;
                    break;
                }
            }
            if (!(notFound)){
                break;
            }
        }
        
    }
    // split nodes into models
    vector<uint> model2NodeIds = pTree->getNodeIndexes(subtreeNodes);
    vector<uint> model1NodeIds;
    for (size_t i = 0; i < nbNodes; i++){
        if (pTree->getRootIndex() == pTree->getNodeIndex(nodes[i])){
            continue;
        }
        auto it = std::find(model2NodeIds.begin(), model2NodeIds.end(), pTree->getNodeIndex(nodes[i]));
        if (it == model2NodeIds.end()){
            model1NodeIds.push_back(pTree->getNodeIndex(nodes[i]));
        }

    }
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    NonHomogeneousSubstitutionProcess* subProSim = new NonHomogeneousSubstitutionProcess(rdist, &parTree);
    subProSim->addModel(chrModel1, model1NodeIds);
    subProSim->addModel(chrModel2, model2NodeIds);
    //subProSim->aliasParameters("Chromosome.gain0_1","Chromosome.loss0_1");
    Context context;
    SubstitutionProcess* nsubPro=subProSim->clone();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *vsc, *nsubPro, true);
    SingleProcessPhyloLikelihood ntl(context, lik, lik->getParameters());
    std::cout << "likelihood is: " << ntl.getValue() << std::endl;

    // set homogeneous model
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel1Homo = std::make_shared<ChromosomeSubstitutionModel>(alpha, mapOfParamsModel1, baseNumber1, 0, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateFuncType);
    DiscreteDistribution* rdistHomo = new GammaDiscreteRateDistribution(1, 1.0);
    auto treeHomo = pTree->clone();
    ParametrizablePhyloTree parTreeHomo(*treeHomo);
    NonHomogeneousSubstitutionProcess* subProSimHomo = NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(chrModel1Homo, rdistHomo, &parTreeHomo);
    Context contextHomo;
    SubstitutionProcess* nsubProHomo=subProSimHomo->clone();
    auto likHomo = std::make_shared<LikelihoodCalculationSingleProcess>(contextHomo, *vsc, *nsubProHomo, true);
    SingleProcessPhyloLikelihood ntlHomo(contextHomo, likHomo, likHomo->getParameters());
    std::cout << "likelihood is for homogeneous model: " << ntlHomo.getValue() << std::endl;

    // merge two models into one
    subProSim->aliasParameters("Chromosome.gain0_1","Chromosome.gain0_2");
    subProSim->aliasParameters("Chromosome.loss0_1","Chromosome.loss0_2");
    subProSim->aliasParameters("Chromosome.dupl0_1","Chromosome.dupl0_2");
    SubstitutionProcess* nsubProMerged=subProSim->clone();
    auto likMerged = std::make_shared<LikelihoodCalculationSingleProcess>(context, *vsc, *nsubProMerged, true);
    SingleProcessPhyloLikelihood ntlMerged(context, likMerged, likMerged->getParameters());
    std::cout << "likelihood is: " << ntlMerged.getValue() << std::endl;


    //static NonHomogeneousSubstitutionProcess* createHomogeneousSubstitutionProcess(
      //std::shared_ptr<BranchModel> model,
      //DiscreteDistribution* rdist,
      //ParametrizablePhyloTree* tree,
      //std::shared_ptr<FrequencySet> rootFreqs = 0,
      //std::shared_ptr<ModelScenario> scenario = 0

    // ParameterList substitutionModelParams = ntl.getSubstitutionModelParameters();
    // std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    // for (size_t i = 0; i < paramsNames.size(); i++){
    //     std::cout << paramsNames[i] << std::endl;
    // }
    // std::cout << "****  model 1 ****" << std::endl;
    // for (size_t i = 0; i < chrModel1->getGenerator().getNumberOfRows(); i++){
    //     for (size_t j = 0; j < chrModel1->getGenerator().getNumberOfColumns(); j++){
    //         std::cout << chrModel1->Qij(i, j) << " ";
    //     }
    //     std::cout << "\n";
    // }
    // std::cout << "****  model 2 ****" << std::endl;
    // for (size_t i = 0; i < chrModel2->getGenerator().getNumberOfRows(); i++){
    //     for (size_t j = 0; j < chrModel2->getGenerator().getNumberOfColumns(); j++){
    //         std::cout << chrModel2->Qij(i, j) << " ";
    //     }
    //     std::cout << "\n";
    // }


    // //RandomTools::setSeed(static_cast<long>(1));
    // //optimizeModelParametersOneDimension(&ntl, alpha, rateFuncType, 5);

    // RowMatrix<double> pij_0_5_m2 = chrModel2->getPij_t(0.5);
    // RowMatrix<double> pij_0_9_m1 = chrModel1->getPij_t(0.9);
    // RowMatrix<double> pij_1_m1 = chrModel1->getPij_t(1);
    // RowMatrix<double> pij_2_m1 = chrModel1->getPij_t(2);
    // RowMatrix<double> pij_0_6_m1 = chrModel1->getPij_t(0.6);
    // RowMatrix<double> pij_1_5_m1 = chrModel1->getPij_t(1.5);

    // Vdouble L_sp1 = {0,1,0,0};
    // Vdouble L_sp2 = {0,0,0,1};
    // Vdouble L_sp3 = {0,0,1,0};
    // Vdouble L_sp4 = {0,0,0,1};
    // Vdouble L_sp5 = {0,0,0,1};
    // Vdouble L_sp6 = {0,0,1,0};

    // // calculating L_sp56
    // Vdouble L_sp56;
    // for (size_t i = 0; i < L_sp5.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp5.size(); j++){
    //         res += (pij_0_9_m1(i, j) * L_sp5[j]);
    //     }
    //     L_sp56.push_back(res);
    // }
    // for (size_t i = 0; i < L_sp6.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp6.size(); j++){
    //         res += (pij_0_9_m1(i, j) * L_sp6[j]);
    //     }
    //     L_sp56[i] *= res;
    // }

    // // calculating L_sp456
    // Vdouble L_sp456;
    // for (size_t i = 0; i < L_sp4.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp4.size(); j++){
    //         res += (pij_1_5_m1(i, j) * L_sp4[j]);
    //     }
    //     L_sp456.push_back(res);
    // }
    // for (size_t i = 0; i < L_sp56.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp56.size(); j++){
    //         res += (pij_0_6_m1(i, j) * L_sp56[j]);
    //     }
    //     L_sp456[i] *= res;
    // }

    // // calculating L_sp23
    // Vdouble L_sp23;
    // for (size_t i = 0; i < L_sp2.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp2.size(); j++){
    //         res += (pij_0_5_m2(i, j) * L_sp2[j]);
    //     }
    //     L_sp23.push_back(res);
    // }
    // for (size_t i = 0; i < L_sp3.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp3.size(); j++){
    //         res += (pij_0_5_m2(i, j) * L_sp3[j]);
    //     }
    //     L_sp23[i] *= res;
    // }

    // // caculating L_sp123
    // Vdouble L_sp123;
    // for (size_t i = 0; i < L_sp1.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp1.size(); j++){
    //         res += (pij_1_m1(i, j) * L_sp1[j]);
    //     }
    //     L_sp123.push_back(res);
    // }
    // for (size_t i = 0; i < L_sp23.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp23.size(); j++){
    //         res += (pij_0_5_m2(i, j) * L_sp23[j]);
    //     }
    //     L_sp123[i] *= res;
    // }

    // // caculating L_sp123456
    // Vdouble L_sp123456;
    // for (size_t i = 0; i < L_sp123.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp123.size(); j++){
    //         res += (pij_2_m1(i, j) * L_sp123[j]);
    //     }
    //     L_sp123456.push_back(res);
    // }
    // for (size_t i = 0; i < L_sp456.size(); i++){
    //     double res = 0;
    //     for (size_t j = 0; j < L_sp456.size(); j++){
    //         res += (pij_1_5_m1(i, j) * L_sp456[j]);
    //     }
    //     L_sp123456[i] *= res;
    // }
    // // get root frequencies
    // Vdouble rootFreqs;
    // double sumOfRootFreqs = 0;
    // for (size_t i = 0; i < L_sp123456.size(); i++){
    //     sumOfRootFreqs += L_sp123456[i];
    // }
    // for (size_t i = 0; i < L_sp123456.size(); i++){
    //     rootFreqs.push_back(L_sp123456[i]/sumOfRootFreqs);
    // }
    // double likelihood = 0;
    // for (size_t i = 0; i < rootFreqs.size(); i++){
    //     likelihood += (rootFreqs[i]*L_sp123456[i]);
    // }
    // double logLikelihood = std::log(likelihood);
    // std::cout << "Manually calculated log likelihood: " << logLikelihood << std::endl;



}

unsigned int optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, ChromosomeAlphabet* alpha, std::vector<int> &modelRateType, size_t maxNumOfIterations){

    // Initialize optimizer

    DerivableSecondOrder* f = tl;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    std::cout <<"max chromosome number: " << alpha->getMax() << endl;
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();

    // setting bracketing for Brent optimization

    optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);

    // initializing the likelihood values
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    unsigned int numOfEvaluations = 0;
    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;
    std::map<string, std::pair<int, uint>> paramNameAndType;
    updateMapsOfParamTypesAndNames(typeWithParamNames, paramNameAndType, tl);
    ParameterList params;
    // starting iterations of optimization
    for (size_t i = 0; i < maxNumOfIterations; i++){

        std::cout << "Iteration #"<<i <<endl;
        
        //ParameterList params = tl->getParameters();// = tl->getParameters();
        ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
        size_t nbParams = substitutionModelParams.size();
        prevLikelihood = currentLikelihood;
        
        for (size_t j = 0; j < nbParams; j ++){
            params = tl->getParameters();
            const string nameOfParam = substitutionModelParams[j].getName();
            //std::cout << "Previous value of " << nameOfParam << " is: " << f->getParameter(nameOfParam).getValue() << endl;
            std::cout << "Previous value of " << nameOfParam << " is: " << params.getParameter(nameOfParam).getValue() << endl;

            int rateParamType = paramNameAndType[nameOfParam].first;

            
            //int rateCompositeParamType;
            double lowerBound;
            double upperBound;
            // param names corresponding to the parameter type
            std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
            Parameter param = params.getParameter(nameOfParam);
            //Parameter* param = &(f->getParameter(nameOfParam));


            auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
            if (it == paramsNames.end()){
                throw Exception("ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(): index out of range!");
            }
            size_t index = it - paramsNames.begin();
            if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
                //rateCompositeParamType = compositeParameter::getCompositeRateType(rateParamType);
                ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(modelRateType[rateParamType-startCompositeParams]);
                ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
                functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, alpha->getMax());
                functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
                delete functionOp;

                // compositeParameter::updateBounds(params, paramsNames, index, &lowerBound, &upperBound, func, alphabet_->getMax());// 24_08 chck if I can remove this function
                // compositeParameter::updateBounds(f, nameOfParam, lowerBound, upperBound, func);
                std::shared_ptr<IntervalConstraint> intervalFuncUpdated = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(nameOfParam).getConstraint());
                double updated_lowerBound = intervalFuncUpdated->getLowerBound();
                std::cout << "*** ***" << nameOfParam << ": Updated lower bound: " << updated_lowerBound << std::endl;


                std::shared_ptr<IntervalConstraint> intervalFuncUpdatedTL = dynamic_pointer_cast<IntervalConstraint>(tl->getParameter(nameOfParam).getConstraint());
                double updated_lowerBoundTL = intervalFuncUpdatedTL->getLowerBound();
                std::cout << "*** ***" << nameOfParam << ": Updated lower bound TL: " << updated_lowerBoundTL << std::endl;  

            }else{
                // baseNumber parameter
                throw Exception("not testing base number!");
            }        
                 
            cout <<"Parameter name is: "<< nameOfParam << endl; 
            //model->checkParametersBounds();     
            if ((i == 1) & (maxNumOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(0.01* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(0.01);
            }
            if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
                optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
            }else{
                optimizer->setInitialInterval(lowerBound, upperBound);
            }            
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            std::cout <<"Parameter value after optimization: "<< tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue() <<endl;
            //std::cout <<"Parameter value after optimization: "<< tl->getLikelihoodCalculation()->getParameter(param->getName()).getValue() <<endl;
            std::cout << "***"<<endl;                        
        }
        printLikParameters(tl);
        
        if (std::abs(prevLikelihood-currentLikelihood) < 0.01){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();
       
    }
    //std::cout << "Number of likelihood evaluations per parameter is "<< numOfBaseNumEval << endl;
    auto model_1 = tl->getSubstitutionProcess().getModel(1);
    auto model_2 = tl->getSubstitutionProcess().getModel(2);
    auto model1Parameters = (model_1->getParameters()).getParameterNames();
    auto model2Parameters = (model_2->getParameters()).getParameterNames();
    std::cout << "Extracting model 1 parameters:" << std::endl;
    for (size_t i = 0; i < model1Parameters.size(); i++){     
        std::cout << model1Parameters[i] << std::endl;
    }
    std::cout << "Extracting model 2 parameters:" << std::endl;
    for (size_t i = 0; i < model2Parameters.size(); i++){     
        std::cout << model2Parameters[i] << std::endl;
    }
 
    std::cout <<"..."<<endl;
    delete optimizer;
    return numOfEvaluations;
}

void updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>> &paramNameAndType, SingleProcessPhyloLikelihood* tl){
    std::map<std::string, int> typeGeneralName;
    updateWithTypeAndCorrespondingName(typeGeneralName);
    ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
    std::vector<std::string> namesAllParams = substitutionModelParams.getParameterNames();
    std::map<std::string, int>::iterator it = typeGeneralName.begin();
    uint numOfModels = static_cast<uint>(tl->getSubstitutionProcess().getNumberOfModels());
    while (it != typeGeneralName.end()){
        string name = it->first;
        int type = it->second;
        std::map<uint, std::vector<std::string>> modelAndParameterNames = ChromosomeNumberOptimizer::getRelatedParameterNamesForEachModel(substitutionModelParams, name, numOfModels);
        typeWithParamNames[type] = modelAndParameterNames;
        auto paramNamesIt = modelAndParameterNames.begin();
        while(paramNamesIt != modelAndParameterNames.end()){
            std::vector<string> parametersNames = paramNamesIt->second;
            for (size_t i = 0; i < parametersNames.size(); i++){
                paramNameAndType[parametersNames[i]] = std::pair<int, uint>(type, paramNamesIt->first);

            }
            
            paramNamesIt ++;
        }
        it ++;
    }

}
// /*******************************************************************************/
void updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName){
    typeGeneralName["gain"] = static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    typeGeneralName["loss"] = static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    typeGeneralName["dupl"] = static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    typeGeneralName["demi"] = static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    typeGeneralName["baseNumR"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    typeGeneralName["baseNum_"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    
}
void printLikParameters(SingleProcessPhyloLikelihood* lik) {


    std:: cout << "Optimized likelihood is : "<< lik->getValue() << endl;

    std:: cout << "Parameters are:" << endl;
 
    ParameterList substitutionModelParams = lik->getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
            std::cout << paramsNames[i] << "value is "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;

        }else{
            std::cout << paramsNames[i] << "value is "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;

        }
        
    }
    std::cout <<"***"<<endl;

}

