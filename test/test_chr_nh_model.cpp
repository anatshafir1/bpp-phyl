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

int main(){
    Newick reader;
    //unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((A:0.1, B:0.2):0.3,C:0.15):0.25,(D:0.35,(E:0.26,F:0.05):0.12):0.16);", false, "", false, false));
    unique_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("((sp1:1,(sp2:0.5,sp3:0.5):0.5):2,(sp4:1.5,(sp5:0.9,sp6:0.9):0.6):1.5);", false, "", false, false));
    ParametrizablePhyloTree parTree(*pTree);
    ChromosomeAlphabet* alpha = new ChromosomeAlphabet(1, 4);
    VectorSiteContainer* vsc = new VectorSiteContainer(alpha);
    const Alphabet* alphabet = static_cast<const Alphabet*>(alpha);
        
    // setting sequence data
    BasicSequence seq1 = BasicSequence("sp1", "2", alphabet);
    BasicSequence seq2 = BasicSequence("sp2", "4", alphabet);
    BasicSequence seq3 = BasicSequence("sp3", "3", alphabet);
    BasicSequence seq4 = BasicSequence("sp4", "4", alphabet);
    BasicSequence seq5 = BasicSequence("sp5", "4", alphabet);
    BasicSequence seq6 = BasicSequence("sp6", "3", alphabet);
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
    int baseNumber1 = IgnoreParam;

    std::vector<int> rateFuncType;
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::IGNORE);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::IGNORE);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::CONSTANT);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::CONSTANT);
    rateFuncType.push_back(ChromosomeNumberDependencyFunction::IGNORE);

    std::map<int, std::vector<double>> mapOfParamsModel1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss1;
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = vector<double>();
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = vector<double>();
    mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = vector<double>();

    std::shared_ptr<ChromosomeSubstitutionModel> chrModel1 = std::make_shared<ChromosomeSubstitutionModel>(alpha, mapOfParamsModel1, baseNumber1, 0, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateFuncType);

    // setting the second model
    std::vector<double> gain2;
    gain2.push_back(2);
    std::vector<double> loss2;
    loss2.push_back(1.5);
    std::vector<double> dupl2;
    dupl2.push_back(2.5);
    int baseNumber2 = IgnoreParam;

    std::map<int, std::vector<double>> mapOfParamsModel2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl2;
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = vector<double>();
    mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = vector<double>();
    std::vector<int> rateFuncType2 = rateFuncType;
    rateFuncType2[1] = ChromosomeNumberDependencyFunction::CONSTANT;

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
    subProSim->aliasParameters("Chromosome.gain0_1","Chromosome.gain0_2");
    Context context;
    SubstitutionProcess* nsubPro=subProSim->clone();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *vsc, *nsubPro, true);
    SingleProcessPhyloLikelihood ntl(context, lik, lik->getParameters());
    std::cout << "likelihood is: " << ntl.getValue() << std::endl;
    ParameterList substitutionModelParams = ntl.getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    for (size_t i = 0; i < paramsNames.size(); i++){
        std::cout << paramsNames[i] << std::endl;
    }
    std::cout << "****  model 1 ****" << std::endl;
    for (size_t i = 0; i < chrModel1->getGenerator().getNumberOfRows(); i++){
        for (size_t j = 0; j < chrModel1->getGenerator().getNumberOfColumns(); j++){
            std::cout << chrModel1->Qij(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "****  model 2 ****" << std::endl;
    for (size_t i = 0; i < chrModel2->getGenerator().getNumberOfRows(); i++){
        for (size_t j = 0; j < chrModel2->getGenerator().getNumberOfColumns(); j++){
            std::cout << chrModel2->Qij(i, j) << " ";
        }
        std::cout << "\n";
    }
    RowMatrix<double> pij_0_5_m2 = chrModel2->getPij_t(0.5);
    RowMatrix<double> pij_0_9_m1 = chrModel1->getPij_t(0.9);
    RowMatrix<double> pij_1_m1 = chrModel1->getPij_t(1);
    RowMatrix<double> pij_2_m1 = chrModel1->getPij_t(2);
    RowMatrix<double> pij_0_6_m1 = chrModel1->getPij_t(0.6);
    RowMatrix<double> pij_1_5_m1 = chrModel1->getPij_t(1.5);

    Vdouble L_sp1 = {0,1,0,0};
    Vdouble L_sp2 = {0,0,0,1};
    Vdouble L_sp3 = {0,0,1,0};
    Vdouble L_sp4 = {0,0,0,1};
    Vdouble L_sp5 = {0,0,0,1};
    Vdouble L_sp6 = {0,0,1,0};

    // calculating L_sp56
    Vdouble L_sp56;
    for (size_t i = 0; i < L_sp5.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp5.size(); j++){
            res += (pij_0_9_m1(i, j) * L_sp5[j]);
        }
        L_sp56.push_back(res);
    }
    for (size_t i = 0; i < L_sp6.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp6.size(); j++){
            res += (pij_0_9_m1(i, j) * L_sp6[j]);
        }
        L_sp56[i] *= res;
    }

    // calculating L_sp456
    Vdouble L_sp456;
    for (size_t i = 0; i < L_sp4.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp4.size(); j++){
            res += (pij_1_5_m1(i, j) * L_sp4[j]);
        }
        L_sp456.push_back(res);
    }
    for (size_t i = 0; i < L_sp56.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp56.size(); j++){
            res += (pij_0_6_m1(i, j) * L_sp56[j]);
        }
        L_sp456[i] *= res;
    }

    // calculating L_sp23
    Vdouble L_sp23;
    for (size_t i = 0; i < L_sp2.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp2.size(); j++){
            res += (pij_0_5_m2(i, j) * L_sp2[j]);
        }
        L_sp23.push_back(res);
    }
    for (size_t i = 0; i < L_sp3.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp3.size(); j++){
            res += (pij_0_5_m2(i, j) * L_sp3[j]);
        }
        L_sp23[i] *= res;
    }

    // caculating L_sp123
    Vdouble L_sp123;
    for (size_t i = 0; i < L_sp1.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp1.size(); j++){
            res += (pij_1_m1(i, j) * L_sp1[j]);
        }
        L_sp123.push_back(res);
    }
    for (size_t i = 0; i < L_sp23.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp23.size(); j++){
            res += (pij_0_5_m2(i, j) * L_sp23[j]);
        }
        L_sp123[i] *= res;
    }

    // caculating L_sp123456
    Vdouble L_sp123456;
    for (size_t i = 0; i < L_sp123.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp123.size(); j++){
            res += (pij_2_m1(i, j) * L_sp123[j]);
        }
        L_sp123456.push_back(res);
    }
    for (size_t i = 0; i < L_sp456.size(); i++){
        double res = 0;
        for (size_t j = 0; j < L_sp456.size(); j++){
            res += (pij_1_5_m1(i, j) * L_sp456[j]);
        }
        L_sp123456[i] *= res;
    }
    // get root frequencies
    Vdouble rootFreqs;
    double sumOfRootFreqs = 0;
    for (size_t i = 0; i < L_sp123456.size(); i++){
        sumOfRootFreqs += L_sp123456[i];
    }
    for (size_t i = 0; i < L_sp123456.size(); i++){
        rootFreqs.push_back(L_sp123456[i]/sumOfRootFreqs);
    }
    double likelihood = 0;
    for (size_t i = 0; i < rootFreqs.size(); i++){
        likelihood += (rootFreqs[i]*L_sp123456[i]);
    }
    double logLikelihood = std::log(likelihood);
    std::cout << "Manually calculated log likelihood: " << logLikelihood << std::endl;






    // if (args == 1){
    //     std::cout << "No arguments provided"<<endl;
    //     return 0;
    // }
    //try{




        /*
            VectorSiteContainer* resized_alphabet_site_container = new VectorSiteContainer(alphabet_);
        for (size_t i = 0; i < numOfSequences; i++){
        BasicSequence seq = vsc->getSequence(sequenceNames[i]);
        BasicSequence new_seq = BasicSequence(seq.getName(), seq.getChar(0), alphabet_);
        resized_alphabet_site_container->addSequence(new_seq);

    }
        */

        
        // BppApplication ChromEvol(args, argv, "ChromEvol");
        // ChromEvolOptions::initAllParameters(ChromEvol);
        // ChromosomeNumberMng* mng = new ChromosomeNumberMng();
        // if (!ChromEvolOptions::simulateData_){
        //     mng->getCharacterData(ChromEvolOptions::characterFilePath_);
        // }
        // mng->getTree(ChromEvolOptions::treeFilePath_, ChromEvolOptions::treeLength_);
        // // setting parameters for model1
        // std::map<int, std::vector<double>> mapOfParamsModel1;
        // mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = ChromEvolOptions::gain_;
        // mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = ChromEvolOptions::loss_;
        // mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = ChromEvolOptions::dupl_;
        // mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = ChromEvolOptions::demiDupl_;
        // mapOfParamsModel1[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = ChromEvolOptions::baseNumR_;
        // int baseNumberModel1 = ChromEvolOptions::baseNum_;
        // const ChromosomeAlphabet* alpha = mng->getAlphabet();
        // auto seqData = mng->getSeqData();
        // auto chrRange = mng->getChromosomeRange();  
        // unsigned int maxBaseNumTransition = (ChromEvolOptions::simulateData_) ? ChromEvolOptions::maxBaseNumTransition_ : chrRange;
        // std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alpha, mapOfParamsModel1, baseNumberModel1, maxBaseNumTransition, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);   

        // // // setting the second model
        // std::vector<double> gain2;
        // std::vector<double> loss2;
        // std::vector<double> dupl2;
        // std::vector<double> demi2;
        // std::vector<double> baseNumR2;
        // int baseNumber2 = 9;
        // // setting parameters for model2
        // gain2.push_back(1);
        // loss2.push_back(2);
        // dupl2.push_back(2);
        // demi2.push_back(0.01);
        // baseNumR2.push_back(0.5);

        // std::map<int, std::vector<double>> mapOfParamsModel2;
        // mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::GAIN)] = gain2;
        // mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::LOSS)] = loss2;
        // mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::DUPL)] = dupl2;
        // mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)] = demi2;
        // mapOfParamsModel2[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)] = baseNumR2;

        // //setting functions
        // vector<int> functions;
        // for (size_t i = 0; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
        //     functions.push_back(ChromosomeNumberDependencyFunction::CONSTANT);
        // }
        // std::shared_ptr<ChromosomeSubstitutionModel> chrModel2 = std::make_shared<ChromosomeSubstitutionModel>(alpha, mapOfParamsModel2, baseNumber2, maxBaseNumTransition, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, functions);

        // get the node ids to be set to model 2
        // PhyloTree* tree = mng->getPhyloTree()->clone();
    //     vector<shared_ptr<PhyloNode> > nodes = tree->getAllNodes();
    //     size_t nbNodes = nodes.size();
    //     vector<shared_ptr<PhyloNode>> subtreeNodes;
    //     vector<uint> leavesUnderNode;
    //     for (size_t i = 0; i < nbNodes; i++){
    //         uint nodeId = tree->getNodeIndex(nodes[i]);
    //         if (nodeId == tree->getRootIndex()){
    //             continue;
    //         }
    //         if (tree->isLeaf(nodeId)){
    //             continue;
    //         }
    //        leavesUnderNode = tree->getLeavesUnderNode(nodeId);
    //         if ((leavesUnderNode.size() >= 5) && (leavesUnderNode.size() < tree->getAllLeavesNames().size())){
    //             subtreeNodes = tree->getSubtreeNodes(tree->getNode(nodeId));
    //             for (size_t j = 0; j < leavesUnderNode.size(); j++){
    //                 string leafName = tree->getNode(leavesUnderNode[j])->getName();
    //                 std::cout << leafName << std::endl;
    //             }
    //             break;

    //         }

    //     }
    //     vector<uint> model2NodeIds = tree->getNodeIndexes(subtreeNodes);
    //     vector<uint> model1NodeIds;
    //     for (size_t i = 0; i < nbNodes; i++){
    //         if (tree->getRootIndex() == tree->getNodeIndex(nodes[i])){
    //             continue;
    //         }
    //         auto it = std::find(model2NodeIds.begin(), model2NodeIds.end(), tree->getNodeIndex(nodes[i]));
    //         if (it == model2NodeIds.end()){
    //             model1NodeIds.push_back(tree->getNodeIndex(nodes[i]));
    //         }

    //     }
    //     ParametrizablePhyloTree parTree(*tree);
    //     DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    //     NonHomogeneousSubstitutionProcess* subProSim = new NonHomogeneousSubstitutionProcess(rdist, &parTree);
    //     subProSim->addModel(chrModel, model1NodeIds);
    //     subProSim->addModel(chrModel2, model2NodeIds);
    //     Context context;
    //     SubstitutionProcess* nsubPro=subProSim->clone();
    //     auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(context, *seqData, *nsubPro, true);
    //     SingleProcessPhyloLikelihood ntl(context, lik, lik->getParameters());
    //     std::cout << "likelihood is: " << ntl.getValue() << std::endl;
    //     ParameterList substitutionModelParams = ntl.getSubstitutionModelParameters();
    //     std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    //     for (size_t i = 0; i < paramsNames.size(); i++){
    //         std::cout << paramsNames[i] << std::endl;
    //     }




        
        
    //     //delete subProSim;
    //     //delete tree;
    //     //std::cout << "****** Max allowed chromosome number: "<< ChromEvolOptions::maxChrNum_ <<endl;
    //     //delete mng;

    // }
    // catch (exception& e)
    // {
    //     cout << e.what() << endl;
    //     return 1;
    // }

    // return 0;
    return 0;


}