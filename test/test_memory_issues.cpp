#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/JointPhyloLikelihood.h>
#include <Bpp/Phyl/Model/Character/CharacterSubstitutionModel.h>
#include <Bpp/Phyl/Model/Character/SingleRateModel.h>
#include <Bpp/Phyl/Model/Character/RatePerEntryModel.h>
#include <Bpp/Phyl/Model/Character/RatePerExitModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairSymModel.h>
#include <Bpp/Phyl/Model/Character/RatePerPairModel.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <iostream>
#include <fstream>
#include <regex>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Phyl/Model/RateDistribution/GammaDiscreteRateDistribution.h>
using namespace bpp;
using namespace std;


void getTraitNodes(vector<uint> &modelNodes, std::shared_ptr<PhyloTree> tree){
    auto nodes = tree->getAllNodes();
    for (size_t i = 0; i < nodes.size(); i++){
        auto nodeId = tree->getNodeIndex(nodes[i]);
        if (nodeId == tree->getRootIndex()){
            continue;
        }
        modelNodes.push_back(nodeId);
    }
 }


void testMemoryTree(){
  Newick reader;
  shared_ptr<PhyloTree> pTree(reader.parenthesisToPhyloTree("(((S1:0.1,S2:0.1):0.3,S3:0.4):0.2,(S4:0.3,S5:0.3):0.3);", false, "", false, false));
  std::shared_ptr<ParametrizablePhyloTree> parTree = std::make_shared<ParametrizablePhyloTree>(*pTree);
  const IntegerAlphabet* alpha = new IntegerAlphabet(1);
  vector<double> freqVals(2);
  freqVals[0] = 0.25;
  freqVals[1] = 0.75;
  double global_rate = 0.5;
  shared_ptr<IntegerFrequencySet> freqs = make_shared<FullIntegerFrequencySet>(alpha, freqVals);
  shared_ptr<CharacterSubstitutionModel> characterModel = make_shared<SingleRateModel>(alpha, freqs, false);
  characterModel->setFrequencySet(*freqs);
  characterModel->setParameterValue("global_rate", global_rate);
  
  std::shared_ptr<DiscreteDistribution> rdistTrait = make_shared<GammaDiscreteRateDistribution>(1, 1.0);
  std::shared_ptr<NonHomogeneousSubstitutionProcess> subProT = make_shared<NonHomogeneousSubstitutionProcess>(std::shared_ptr<DiscreteDistribution>(rdistTrait->clone()), parTree);
  std::vector<uint> modelNodesTrait;
  getTraitNodes(modelNodesTrait, pTree);
  subProT->addModel(std::shared_ptr<CharacterSubstitutionModel>(characterModel->clone()), modelNodesTrait);
  Context* context = new Context();
  std::shared_ptr<VectorSiteContainer> sites = make_shared<VectorSiteContainer>(alpha);
  sites->addSequence(BasicSequence("S1", "1", alpha));
  sites->addSequence(BasicSequence("S2", "1", alpha));
  sites->addSequence(BasicSequence("S3", "0", alpha));
  sites->addSequence(BasicSequence("S4", "0", alpha));
  sites->addSequence(BasicSequence("S5", "1", alpha));
  auto likT = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *sites->clone(), *(subProT->clone()));
  auto phyloT = new SingleProcessPhyloLikelihood(*context, likT);
  auto lik = phyloT->getValue();
  std::cout << "value is: " << lik << std::endl;
  auto sequenceDataTrait = phyloT->getData();
  auto processTrait = &(phyloT->getSubstitutionProcess());
  auto contextT = &(phyloT->getContext());
  delete sequenceDataTrait;
  delete processTrait;
  delete contextT;
  delete alpha;
  
}
int main(){
    testMemoryTree();
    return 0;

}
