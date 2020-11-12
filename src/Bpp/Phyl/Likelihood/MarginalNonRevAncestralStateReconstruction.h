//
// File: MarginalNonRevAncestralStateReconstruction.h
// Created by: Anat Shafir
// Created on: Wedensday Jul 01 15:04 2020
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

#ifndef _MARGINALNONREVANCESTRALSTATESRECONSTRUCTION_H_
#define _MARGINALNONREVANCESTRALSTATESRECONSTRUCTION_H_

#include "../AncestralStateReconstruction.h"
#include "DRNonHomogeneousTreeLikelihood.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>

// From the STL:
#include <vector>
#include <map>
using namespace std;

namespace bpp
{

/**
 * @brief Likelihood ancestral states reconstruction: marginal method.
 *
 * Reference:
 * Z Yang, S Kumar and M Nei (1995), _Genetics_ 141(4) 1641-50.
 */
class MarginalNonRevAncestralStateReconstruction
{
	private:
	    DRNonHomogeneousTreeLikelihood* likelihood_;
        TreeTemplate<Node> tree_;
		const Alphabet* alphabet_;
		size_t nbSites_;
		size_t nbDistinctSites_;
		size_t nbClasses_;
		size_t nbStates_;
        std::vector<size_t> rootPatternLinks_; //not used (maybe will be needed for some other purposes)
        std::vector<double> r_;//probability for a specific rate class
        std::vector<double> l_;//likelihood for site
        std::map<int, map<size_t, std::vector<double>>>* postProbNode_; // node->site-> vector of posterior probabilities where the indices are states
        std::map<int, std::map<size_t, VVdouble>>* jointProbabilities_; //node->site:(nodeState,fatherState)-> joint probability of father state and a given node state
		
	public:
		MarginalNonRevAncestralStateReconstruction(DRNonHomogeneousTreeLikelihood* drl) :
            likelihood_      (drl),
            tree_            (drl->getTree()),
		    alphabet_        (drl->getAlphabet()),
	        nbSites_         (drl->getLikelihoodData()->getNumberOfSites()),
		    nbDistinctSites_ (drl->getLikelihoodData()->getNumberOfDistinctSites()),
		    nbClasses_       (drl->getLikelihoodData()->getNumberOfClasses()),
		    nbStates_        (drl->getLikelihoodData()->getNumberOfStates()),
		    rootPatternLinks_(drl->getLikelihoodData()->getRootArrayPositions()),
            r_               (drl->getRateDistribution()->getProbabilities()),
            l_               (drl->getLikelihoodData()->getRootRateSiteLikelihoodArray()),
            postProbNode_    (0),
            jointProbabilities_(0)
        {}

        MarginalNonRevAncestralStateReconstruction(const MarginalNonRevAncestralStateReconstruction& masr) :
            likelihood_      (masr.likelihood_),
            tree_            (masr.tree_),
            alphabet_        (masr.alphabet_),
            nbSites_         (masr.nbSites_),
            nbDistinctSites_ (masr.nbDistinctSites_),
            nbClasses_       (masr.nbClasses_),
            nbStates_        (masr.nbStates_),
            rootPatternLinks_(masr.rootPatternLinks_),
            r_               (masr.r_),
            l_               (masr.l_),
            postProbNode_    (new std::map<int, map<size_t, std::vector<double>>>(*masr.postProbNode_)),
            jointProbabilities_(new std::map<int, std::map<size_t, VVdouble>>(*masr.jointProbabilities_))
        {}

        MarginalNonRevAncestralStateReconstruction& operator=(const MarginalNonRevAncestralStateReconstruction& masr)
        {
            likelihood_       = masr.likelihood_;
            tree_             = masr.tree_;
            alphabet_         = masr.alphabet_;
            nbSites_          = masr.nbSites_;
            nbDistinctSites_  = masr.nbDistinctSites_;
            nbClasses_        = masr.nbClasses_;
            nbStates_         = masr.nbStates_;
            rootPatternLinks_ = masr.rootPatternLinks_;
            r_                = masr.r_;
            l_                = masr.l_;
            postProbNode_     = new std::map<int, map<size_t, std::vector<double>>>(*masr.postProbNode_);
            jointProbabilities_ = new std::map<int, std::map<size_t, VVdouble>>(*masr.jointProbabilities_);
            return *this;
        }

        MarginalNonRevAncestralStateReconstruction* clone() const { return new MarginalNonRevAncestralStateReconstruction(*this); }

	    virtual ~MarginalNonRevAncestralStateReconstruction() {
            delete postProbNode_;
            delete jointProbabilities_;
        }

    private:

    	/**
        * @brief Compute the posterior probabilities consitional on the root state for each state for a given node given its state and the state of the father.
        * @param drl A DR tree likelihood object.
        * @param nodeId The nodeId at which probabilities must be computed.
        * @param nodeState The state of the node
        * @param fatherState The state of the father
        * @param rootState The root state
        * @param site The site at which the calculation should be performed
        * @param nbClasses Number of rate classes
        * @return the joint probability of the father state and the state of the node given the state of the root
        */
    
        double getJointLikelihoodFatherNode(
            DRNonHomogeneousTreeLikelihood* drl, 
            int nodeId, 
            size_t nodeState, 
            size_t fatherState,
            size_t rootState,
            size_t site,
            size_t nbClasses);


        /**
        * @brief Compute the posterior probabilities consitional on the root state for each state for a given node given its state and the state of the father.
        * @param drl A DR tree likelihood object.
        * @param rootState The state of the root on which the posterior probability is conditioned
        * @param nodeIds All node Ids in the  tree
        * @return For each node, for each site- the posterior probability of the node x being assigned to a partcular state given the root state.
        */

        map<int, map<size_t, std::vector<double>>>* getPosteriorProbabilitiesOfNodesForEachRootStatePerSite(
            DRNonHomogeneousTreeLikelihood* drl,
            size_t rootState,
            std::vector<int> nodeIds);



	public:

        /**
        * @brief Compute the posterior probabilities of each node and state
        * @return For each node, for each site- the posterior probability of the node x being assigned to a partcular state
        */


        void computePosteriorProbabilitiesOfNodesForEachStatePerSite();
        
        std::map<int, std::map<size_t, VVdouble>> getAllJointFatherNodeProbabilities(){
            return *jointProbabilities_;
        }


        std::map<int, map<size_t, std::vector<double>>>* getPosteriorProbForAllNodesAndStatesPerSite(){
            if (!postProbNode_){
                computePosteriorProbabilitiesOfNodesForEachStatePerSite();
            }
            return postProbNode_;
        }


        const std::map<int, std::vector<size_t> > getAllAncestralStates() const;
        vector <double> getRootPosteriorProb() const;

        
        

	
 

		
};

} //end of namespace bpp.

#endif // _MARGINALNONREVANCESTRALSTATESRECONSTRUCTION_H_