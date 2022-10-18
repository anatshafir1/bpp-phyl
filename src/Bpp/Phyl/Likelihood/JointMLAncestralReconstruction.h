//
// File: JointMLAncestralReconstruction.h
// Created by: Anat Shafir
// Created on: Monday March 29 11.25 2021
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
#ifndef _JOINTMLANCESTRALRECONSTRUCTION_H_
#define _JOINTMLANCESTRALRECONSTRUCTION_H_

#include "../AncestralStateReconstruction.h"

#include "DataFlow/LikelihoodCalculationSingleProcess.h"
#include "DataFlow/DataFlowCWise.h"

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/AlignedSequenceContainer.h>
#include <Bpp/Seq/Sequence.h>

// From the STL:
#include <vector>

namespace bpp{
    class JointMLAncestralReconstruction:
        public virtual AncestralStateReconstruction
    {
        private:
            std::shared_ptr<LikelihoodCalculationSingleProcess> likelihood_;
            const ParametrizablePhyloTree* tree_;
            const Alphabet* alphabet_;
            size_t nbSites_;
            size_t nbDistinctSites_;
            // size_t nbClasses_;
            size_t nbStates_;
            PatternType rootPatternLinks_;
            std::map<uint, std::vector<size_t>>* ancestors_;
            std::map<uint, uint>* fathers_;

        public:
            JointMLAncestralReconstruction(std::shared_ptr<LikelihoodCalculationSingleProcess> drl) :
                likelihood_      (drl),
                tree_            (&drl->getSubstitutionProcess().getParametrizablePhyloTree()),
                alphabet_        (drl->getStateMap().getAlphabet()),
                nbSites_         (drl->getNumberOfSites()),
                nbDistinctSites_ (drl->getNumberOfDistinctSites()),
                // nbClasses_       (drl->getNumberOfClasses()),
                nbStates_        (drl->getStateMap().getNumberOfModelStates()),
                rootPatternLinks_(drl->getRootArrayPositions()),
                ancestors_(0),
                fathers_(0)
            {}

            JointMLAncestralReconstruction(const JointMLAncestralReconstruction& masr) :
                likelihood_      (masr.likelihood_),
                tree_            (masr.tree_),
                alphabet_        (masr.alphabet_),
                nbSites_         (masr.nbSites_),
                nbDistinctSites_ (masr.nbDistinctSites_),
                // nbClasses_       (masr.nbClasses_),
                nbStates_        (masr.nbStates_),
                rootPatternLinks_(masr.rootPatternLinks_),
                ancestors_(new std::map<uint, std::vector<size_t>>(*masr.ancestors_)),
                fathers_(new std::map<uint, uint>(*masr.fathers_))
            {}

            JointMLAncestralReconstruction& operator=(const JointMLAncestralReconstruction& masr)
            {
                likelihood_       = masr.likelihood_;
                tree_             = masr.tree_;
                alphabet_         = masr.alphabet_;
                nbSites_          = masr.nbSites_;
                nbDistinctSites_  = masr.nbDistinctSites_;
                // nbClasses_        = masr.nbClasses_;
                nbStates_         = masr.nbStates_;
                rootPatternLinks_ = masr.rootPatternLinks_;
                ancestors_ = new std::map<uint, std::vector<size_t>>(*masr.ancestors_);
                fathers_ = new std::map<uint, uint>(*masr.fathers_);
                return *this;
            }


            JointMLAncestralReconstruction* clone() const { return new JointMLAncestralReconstruction(*this); }

            virtual ~JointMLAncestralReconstruction() {
                delete ancestors_;
                delete fathers_;
            }
            public:
                void init();
                std::vector<size_t> getAncestralStatesForNode(uint nodeId) const;
                std::map<uint, std::vector<size_t> > getAllAncestralStates() const;
                Sequence* getAncestralSequenceForNode(uint nodeId) const;
                AlignedSequenceContainer* getAncestralSequences() const;

            private:
                void recursiveJointMLAncestralStates(
                    const std::shared_ptr<PhyloNode> node,
                    AlignedValuesContainer& data) const;

    };

}
#endif // _JOINTMLANCESTRALRECONSTRUCTION_H_