//
// File: FormulaOfPhyloLikelihood.h
// Authors:
//   Laurent GuÃÂ©guen
// Created: jeudi 8 dÃÂ©cembre 2016, ÃÂ  10h 35
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_JOINTOFPHYLOLIKELIHOOD_H
#define BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_JOINTOFPHYLOLIKELIHOOD_H

#include <Bpp/Numeric/Function/Operators/ComputationTree.h>

#include "SetOfAbstractPhyloLikelihood.h"
#include <Bpp/Phyl/Mapping/StochasticMapping.h>
#include "SingleDataPhyloLikelihood.h"
#include "SingleProcessPhyloLikelihood.h"
#include <Bpp/Phyl/Likelihood/NonHomogeneousSubstitutionProcess.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Likelihood/JointMLAncestralReconstruction.h>

namespace bpp
{
/**
 * @brief The JointPhyloLikelihood class is intended for likelihood computation of the following form- given two types of data in the leaves,
 * D1 and D2, if the two evolutionary process are not independent, then if for example, D1= {0,1}, each state dictates a different evolutionary process for D2, therefore the formula is:
 * P(D1,D2| teta) = P(D1|teta)*P(D2|D1,teta). The first term can be calculated via the Felsenstein algorithm, while in order to compute the second term, we have to know the history
 * of the transitions of D1 anlong the phylogeny. This history is unknown, but can be estimated using a large number of stochastic mappings, that can be summarized to generate
 * an expected history of trajectories. Knowing the expected history will allow to assign each model of D2 according to the states of D1, and thus enable to calculate the second term using Felsenstein
 * 
 *
 * WARNING: This formula applies on the log-likelihoods (ie getValues())
 *
 */ 

class JointPhyloLikelihood :
  public SetOfAbstractPhyloLikelihood
{
protected:

  std::shared_ptr<LikelihoodCalculation> likCal_;

  bool traitOptimization_;
  bool expectedHistory_;
  size_t numOfMappings_;
  bool weightedFrequencies_;
  SingleProcessPhyloLikelihood* tempLik_;
  bool firstLikChange_;
  std::shared_ptr<PhyloTree> tempTree_;
  bool ML_;


public:
  JointPhyloLikelihood(Context& context, std::shared_ptr<PhyloLikelihoodContainer> pC, bool expectedHistory, bool weightedFrequencies, size_t numOfMappings, bool ML, bool inCollection = true);

  ~JointPhyloLikelihood() {
    // auto sequenceData = tempLik_->getData();
    // auto process = &(tempLik_->getSubstitutionProcess());
    //auto contextDel = &(tempLik_->getContext());
    // delete process;
    // delete sequenceData;
    // if (getPhyloContainer()->getContext() != contextDel){
    //   delete contextDel;
    // }
    //delete tempLik_;
  }

  JointPhyloLikelihood* clone() const
  {
    return new JointPhyloLikelihood(*this);
  }
  SingleProcessPhyloLikelihood* getPhylo2() const{
    return tempLik_;
  }
  const std::shared_ptr<PhyloTree> getStochasticMappingTree() const{
    return tempTree_;
  }
  void setStochasticMappingTree(std::shared_ptr<PhyloTree> tree){
    tempTree_ = tree;
  }
  std::shared_ptr<PhyloTree> getStochasticMappingTree(){
    return tempTree_;
  }

  JointPhyloLikelihood(const JointPhyloLikelihood& sd);
  // virtual void optimizeTraitModel(double tol, uint numOfIterations);
  // virtual void optimizeSecondModel(double tol, uint numOfIterations);
  // virtual void optimize(double tol, uint numOfIterations, uint numberOfIterationsPerOneOptimization);

public:


  /**
   * @brief Get the logarithm of the likelihood for the whole dataset.
   *
   * @return The logarithm of the likelihood of the dataset.
   */
  std::shared_ptr<LikelihoodCalculation> getLikelihoodCalculation() const
  {
    return likCal_;
  }

  void fireParameterChanged(const ParameterList& params);

  void setParamUpdateMode(bool traitUpdated){
    traitOptimization_ = traitUpdated;
  }

protected:
  /**
   * @brief Build the LikelihoodNode from the computation Tree
   *
   */
  ValueRef<DataLik> makeLikelihoods();
  std::shared_ptr<FrequencySet> copyRootFrequencies(const NonHomogeneousSubstitutionProcess* prevSubstitutionModel, SingleProcessPhyloLikelihood* lik);
  std::map<uint, std::vector<size_t>> getMLAncestralReconstruction(SingleProcessPhyloLikelihood* likProcess);


};
} // end of namespace bpp.
#endif // BPP_PHYL_LIKELIHOOD_PHYLOLIKELIHOODS_JOINTPHYLOLIKELIHOOD_H
