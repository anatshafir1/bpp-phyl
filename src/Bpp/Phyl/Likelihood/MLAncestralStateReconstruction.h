//
// File: MLAncestralStateReconstruction.h
// Created by: Anat Shafir
// Created on: Tuesday June 16 17:52 2020
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _MLANCESTRALSTATERECONSTRUCTION_H_
#define _MLANCESTRALSTATERECONSTRUCTION_H_

#include "../AncestralStateReconstruction.h"
#include "DRTreeLikelihood.h"
#include "DRASRTreeLikelihoodData.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

// From the STL:
#include <vector>
#include <map>
using namespace std;
namespace bpp

{

/**
 * @brief Likelihood ancestral states reconstruction: ML method.
 *
 * Reference:
 * T. Pupko, I. Pe'er, R. Shamir et al (2000), Molecular Biology and Evolution 17(6) 890-896.
 */
class MLAncestralStateReconstruction
{
	private:
    mutable DRASRTreeLikelihoodData *likelihoodData_;
    TreeTemplate<Node> tree_;
		const Alphabet* alphabet_;
    const SiteContainer* data_;
    const TransitionModel* model_;
    std::vector<double> rootFrequiencies_;
    std::map<int, VVVdouble>* pxy_;
		size_t nbSites_;
		size_t nbDistinctSites_;
		size_t nbClasses_;
		size_t nbStates_;
    std::vector<size_t> rootPatternLinks_;
    std::vector<double> r_;
    size_t nbCategories_;
    std::map<int, std::map<size_t, std::vector<size_t>>>* ancestors_;
    //std::vector<double> l_;
		
	public:
		MLAncestralStateReconstruction(const DRTreeLikelihood* drl, const TransitionModel* model, std::vector<double> rootFrequencies, std::map<int, VVVdouble>* Pijt)://map<int, VVVdouble>* Pijt)://, )://,
      likelihoodData_(0),
      tree_(drl->getTree()),
      alphabet_ (drl->getAlphabet()),
      data_(drl->getData()),
      model_(model),
      rootFrequiencies_(rootFrequencies),
      pxy_ (Pijt),
			nbSites_ (drl->getLikelihoodData()->getNumberOfSites()),
			nbDistinctSites_(drl->getLikelihoodData()->getNumberOfDistinctSites()),
			nbClasses_ (drl->getLikelihoodData()->getNumberOfClasses()),
			nbStates_ (drl->getLikelihoodData()->getNumberOfStates()),
			rootPatternLinks_(drl->getLikelihoodData()->getRootArrayPositions()),
      r_ (drl->getRateDistribution()->getProbabilities()),
      nbCategories_(drl->getRateDistribution()->getNumberOfCategories()),
      ancestors_(0)
    {}

    MLAncestralStateReconstruction(const MLAncestralStateReconstruction& masr) :
      likelihoodData_  (masr.likelihoodData_),
      tree_            (masr.tree_),
      alphabet_        (masr.alphabet_),
      data_            (masr.data_),
      model_           (masr.model_),
      rootFrequiencies_(masr.rootFrequiencies_),
      pxy_             (masr.pxy_),
      nbSites_         (masr.nbSites_),
      nbDistinctSites_ (masr.nbDistinctSites_),
      nbClasses_       (masr.nbClasses_),
      nbStates_        (masr.nbStates_),
      rootPatternLinks_(masr.rootPatternLinks_),
      r_               (masr.r_), 
      nbCategories_    (masr.nbCategories_),
      ancestors_       (masr.ancestors_)
    {}

    MLAncestralStateReconstruction& operator=(const MLAncestralStateReconstruction& masr)
    {
      likelihoodData_   = masr.likelihoodData_;
      tree_             = masr.tree_;
      alphabet_         = masr.alphabet_;
      data_             = masr.data_;
      model_            = masr.model_;
      rootFrequiencies_ = masr.rootFrequiencies_;
      pxy_              = masr.pxy_;
      nbSites_          = masr.nbSites_;
      nbDistinctSites_  = masr.nbDistinctSites_;
      nbClasses_        = masr.nbClasses_;
      nbStates_         = masr.nbStates_;
      rootPatternLinks_ = masr.rootPatternLinks_;
      r_                = masr.r_;
      nbCategories_     = masr.nbCategories_;
      ancestors_        = masr.ancestors_;
      return *this;
    }

    MLAncestralStateReconstruction* clone() const { return new MLAncestralStateReconstruction(*this); }

	virtual ~MLAncestralStateReconstruction() {
    delete likelihoodData_;
    delete ancestors_;
  }

	public:
		
        //std::vector<size_t> getAncestralStatesForNode(int nodeId) const;
    void computeJointLikelihood();		
    std::map<int, std::vector<size_t> > getAllAncestralStates() const;		
        //Sequence* getAncestralSequenceForNode(int nodeId) const;
        //SiteContainer* getAncestralSequences() const;

  private:
    void fillRootLikelihoodsArrays(const Node* node, VVVdouble* likNode, size_t nbSites, size_t nbClasses, size_t nbStates, size_t nbNodes);
    void fillLikelihoodsArrays(const Node* node, VVVdouble* likNode, VVVdouble* currNodeProb, size_t nbSites, size_t nbClasses, size_t nbStates, size_t nbNodes);
    void computeSubtreeBestLikelihood(const Node* node);
    void getAllAncestralStatesRec(const Node* node, std::map<int, std::vector<size_t> > &ancestors) const;
		
};

} //end of namespace bpp.

#endif // _MLANCESTRALSTATERECONSTRUCTION_H_

