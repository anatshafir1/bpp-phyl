//
// File: StochasticMapping.h
// Created by: Keren Halabi
// Created on: June 2018
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef ___STOCHASTIC_MAPPING_H
#define ___STOCHASTIC_MAPPING_H

#include "../Simulation/SubstitutionProcessSequenceSimulator.h"
#include "../Simulation/MutationProcess.h"

#include "../NewLikelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "../NewLikelihood/DataFlow/DataFlowCWise.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>

using namespace std;

/* Store the countings on a DAG similar to the computing DAG */



typedef vector<vector<vector<double> > > VVVDouble;
typedef vector<vector<double> > VVDouble;
typedef vector<double> VDouble;

/* class for reprenting the framework of Stochastic mapping
 *
 *   A StochasticMapping instance can be used to sample histories of
 *   state transitions along a tree, given a substitution model and the
 *   states at the tip taxa. For more information, see: Nielsen, Rasmus.
 *   "Mapping mutations on phylogenies." Systematic biology 51.5 (2002):
 *   729-739.
 */

namespace bpp
{
  class StochasticMapping
  {
  protected:
    /*
     * @brief The tree likelihood instance is used for computing the
     * the conditional sampling probabilities of the ancestral states as
     * well as the root assignment probabilities.
     *
     */
  
    std::shared_ptr<LikelihoodCalculationSingleProcess> likelihood_;  

    std::shared_ptr<PhyloTree> tree_;

    /*
     * @brief this instance will hold the parameters required for the
     * stochastic mapping procedure, and be used to generate stochastic
     * mappings.
     */
  
//    SimpleSubstitutionProcessSequenceSimulator mappingParameters_;

  
    VVVDouble ConditionalProbabilities_;             // vector that holds the conditionl states assignment probabilities of the nodes in the tree (node*father_states*son_states)
    size_t nodesCounter_;                            // counter of nodes hat allows adding unique names to the generated nodes while breaking branching in a mapping
    size_t numOfMappings_;                           // the number of stochastic mappings to generate
    map<uint, vector<size_t>> ancetralStates_;         // The sampled ancestral states For each node and for each mapping
    map<uint, vector<MutationPath>> mappings_;        // A map of nodes and the mappings of the branch that leads to them
    map<size_t, VVdouble> jumpsProbs_;                  // For each model: Jump probabilities: for each j: Qij/-Qii

  public:
    /* constructors and destructors */

    explicit StochasticMapping(std::shared_ptr<LikelihoodCalculationSingleProcess> drl, size_t numOfMappings = 1000000); // it is a good general practice to use "explicit" keyword on constructors with a single argument: https://stackoverflow.com/questions/121162/what-does-the-explicit-keyword-mean

    ~StochasticMapping();

    StochasticMapping(const StochasticMapping& sm) : 
      likelihood_(sm.likelihood_),
      tree_(sm.tree_),
//      mappingParameters_(likelihood_->getSubstitutionProcess()),
      ConditionalProbabilities_(sm.ConditionalProbabilities_),
      nodesCounter_(0), numOfMappings_(sm.numOfMappings_),
      ancetralStates_(),
      mappings_(),
      jumpsProbs_(sm.jumpsProbs_)
    { 

    }

    /**
     * @brief cloning function used by the copy constructor of
     * ./Likelihood/JointLikelihoodFunction/h
     *
     */
  
    StochasticMapping* clone() const { return new StochasticMapping(*this); }


    /*
     * @brief generates a stochastic mappings based on the sampling
     * parameters
     *
     * @param     Number of histories to sample
     *
     */

    //void generateStochasticMapping(std::vector<std::shared_ptr<PhyloTree>>& mappings);
    void generateStochasticMapping();

    /**
     *@brief Creates a single expected (i.e, average) history based on
     * a given set of mappings steps. Correspond to Nielsen, Rasmus.
     * "Mapping mutations on phylogenies." Systematic biology 51.5
     * (2002): 729-739.
     *
     * the function assumes that there is only one site to simulate history
     *
     * @param mappings          A vector of stochastic mappings to average
     *
     * @param divMethod The method used in the case that the son and
     * father share the same state (either divide the dwelling time of
     * the staed state by 2 for two transitions (method 0) or allocate
     * the entire dwelling time to be adjacent to the son(method 1))
     *
     **/
  
    std::shared_ptr<PhyloTree> generateExpectedMapping(std::vector<std::shared_ptr<PhyloTree>>& mappings, size_t divMethod = 0);

    /**
     *@brief Creates a single expected (i.e, average) history based the rewards provided by te algorithm of Minin and Suchard (2008)
     * the function assumes that there is only one site to simulate history for
     *
     *@param divMethod The method used in the case that the son and
     *father share the same state (either divide the wdelling time of
     *the staed state by 2 for two transitions (method 0) or allocate
     *the entire dwelling time to be adjacent to the son(method 1))
     */
  
    std::shared_ptr<PhyloTree> generateAnalyticExpectedMapping(size_t divMethod = 0);

    /* extracts the state of a node in a mapping
     * @param node              The node to get the state of
     * @return                  Node state is int
     */
    int getNodeState(const PhyloNode* node) const;
    /**
     *@brief Gets the dwelling times for each state in a given mapping. Note: the function assums thay the 
     *        vector dwellingTimes is already resized according to the number of states!
     *
     *@param dwellingTimes A vector where the dewelling times for each state will be stored
     *@param mappingIndex The index of the mapping
     */
    void getDewellingTimesUnderEachStatePerMapping(vector<double> &dwellingTimes, size_t mappingIndex);
    /**
     *@brief Gets the of ocurrences of each transition in a given mapping. 
     *@param mappingIndex The index of the mapping
     *@param transitionOcurrences The output map of occurrences for each transition
     */
    void getNumOfOcuurencesForEachTransitionPerMapping(size_t mappingIndex, std::map<uint, std::map<std::pair<size_t, size_t>, double>> &transitionOcurrences);
    /**
     *@brief Gets the expected dwelling time under each states given all the mappings
     *@return A vector that contains the expected dwelling times under each state
    */
    Vdouble getExpectedDwellingTimesUnderEachState();
    /**
     *@brief Gets the expected number of occurences of each transition in each node (summarizing over all the mapppings)
     *@return A map of the expected number of occurences for each pair of states (for each transition)
    */
    std::map<uint, std::map<pair<size_t, size_t>, double>> getExpectedNumOfOcuurencesForEachTransitionPerNode();
      /**
     *@brief Gets the number of occurences of each transition in each node (summarizing over all the mappings)
     *@return A map of the number of occurences for each pair of states (for each transition)
    */
    std::map<uint, std::map<pair<size_t, size_t>, double>> getNumOfOcuurencesForEachTransitionPerNode();

    /**
     *@brief Gets the expected number of occurences of each transition in the whole tree (summarizing over all the mapppings)
     *@return A map of the expected number of occurences for each pair of states (for each transition)
    */
    std::map<pair<size_t, size_t>, double> getExpectedNumOfOcuurencesForEachTransition(std::map<uint, std::map<pair<size_t, size_t>, double>>* expectedTransitionOccurencesPerNode = 0);

    /**
     *@brief Gets the total number of occurences of each transition in the whole tree (summarizing over all the mapppings)
     *@return A map of the total number of occurences for each pair of states (for each transition)
    */
    std::map<pair<size_t, size_t>, double> getTotalNumOfOcuurencesForEachTransition(std::map<uint, std::map<pair<size_t, size_t>, double>>* expectedTransitionOccurencesPerNode = 0);
    /**
     *@brief Gets the total number of occurences of each transition in the whole tree for each mapping
     *@return A map of mappings, which contains for each encountered pair of states its number of occurrences
    */

    std::map<size_t, std::map<std::pair<size_t, size_t>, double>> getNumOfOccurencesForEachTransitionForEachMapping();
    /**
     *@brief Gets the total the rates of a transition given a certain state. 
     *@return The rates of each possible transition under each state rates. For example,
     * the [i][j] value will represent the rate of (i,j) transition given that we are in state i.
     * These rates are integrated over all the mappings
    */
    VVdouble getExpectedRateOfTransitionGivenState();
    /**
     *@brief Gets the total dwelling time under each state for each mapping
     *@return The total dwelling time under each state for each mapping. The first index represents the 
     * mapping index, while the second index represents the state index.
    */

    VVdouble getDwellingTimeOfStatePerEachMapping();

    
    

  private:
    void initMapOfNumOfOccurences(std::map<uint, std::map<pair<size_t, size_t>, double>> &transitionOcurrences);
    std::map<pair<size_t, size_t>, double> sumTotalOccurences(std::map<uint, std::map<pair<size_t, size_t>, double>>* transitionOcurrencesPerNodePtr);
    /**
     *Gets the of ocurrences of each transition in a given mapping recursively. 
     * @param nodeId      The current node for which the mutation path is inspected
     * @param initialState  The initial state of the branch
     * @param mappingIndex The index of the mapping
     * @param transitionOcurrences A reference to the output map of occurences
     */
    void getNumOfOcuurencesForEachTransitionPerMappingRecursively(uint nodeId, size_t initialState, size_t mappingIndex, std::map<uint, std::map<pair<size_t, size_t>, double>> &transitionOcurrences);
    /* Gets the dwelling times for each state in a given mapping from a given node
     * @param nodeId               Current node in the recursive call.
     * @param initialState               Intital state in the mutation path
     * @param  A vector where the dewelling times for each state will be stored
     * @param mappingIndex               Mapping index
     */
    void getDewellingTimesUnderEachStatePerMappingRecursively(uint nodeId, size_t initialState, vector<double> &dwellingTimes, size_t mappingIndex);

    /* Fills the transition probabilities given that itransition has occured (Qij/-Qii).
     * These probabilities are filled for each model.
    */
    void initJumpProbs();
    /* adds names to the internal nodes, in case of absence.
     * @param tree               The tree whose nodes should be edited if needed.
     */
    void giveNamesToInternalNodes(PhyloTree& tree);
    /* samples random state given the current state
     * @param beginState          The current state
     * @param modelIndex          The model index
     */

    size_t giveRandomState(size_t beginState, size_t modelIndex) const;
    
    /* sets the state of a node in a mapping
     * @param node               The node to get the state of
     * @param state              The state that needs to be assigned to the node
     */
    void setNodeState(PhyloNode* node, size_t state);

    /* set the character states of the leafs as properties of thier nodes instances
     * @param mapping - the tree to sets the properties in
     */
    void setLeafsStates(std::shared_ptr<PhyloTree> mapping);



    /* compute the conditional probabilities of all the nodes assignments.
     * The function assumes that the data is a character data, i.e., one site.
     * @param rootProbabilities  The root frequencies
     */
    void ComputeConditionals();

    /* compute the ancestral frequenceis of character states of all the nodes based on the mappings
     * @param                     A vector of the posterior probabilities probabilities to fill in (node**state combinaion in each entry)
     * @param                     A vector of mappings to base the frequencies on
     */
    void computeStatesFrequencies(VVDouble& ancestralStatesFreuquencies, vector<shared_ptr<PhyloTree>>& mappings);

    /* auxiliary function that samples a state based on a given discrete distribution
     * @param distibution       The distribution to sample states based on
     */
    size_t sampleState(const VDouble& distibution); // k: best by ref
    /* Fill root conditional ancestral probabilities
     * @param conditional fractions in extended float types
    */
    void fillRootConditionals(ExtendedFloatArrayXd &conditionals);
    /* samples ancestral states based on the conditional probabilities at each node in the base (user input) tree and the root assignment probabilities. States will be updated as nodes properties
     * @param mapping               The tree whose nodes names should be updated according to their assigned states.
     */
    //void sampleAncestrals(shared_ptr<PhyloTree> mapping);
    /* samples ancestral states based on the conditional probabilities at each node in the tree
     * @param mappingIndex         The index of the mapping
     */
    void sampleAncestrals(size_t mappingIndex);

    /* samples ancestral states based on the conditional probabilities at each node in the tree recursively
     * @param nodeId               Node Index for which we do the calculation recursively
     * @param mappingIndex         The index of the mapping
     * @param fatherIndex          A pointer to the father index (if father exists)
    */

    void sampleAncestralsRecursively(uint nodeId, size_t mappingIndex, uint *fatherIndex = 0);

    /* set ancestral states in the expected history based on the conditional probabilities at each node in the base (user input) tree and the root assignment probabilities. States will be updated as nodes properties
     * @param expectedMapping           The expected mapping instance whose nodes names should be updated according to their assigned states.
     * @param posteriorProbabilities    Vector of posterior assignment proabilities to inner node to decide on assignments
     */
    void setExpectedAncestrals(shared_ptr<PhyloTree> expectedMapping, VVDouble& posteriorProbabilities);

    /* simulates mutations on phylogeny based the sampled ancestrals, tips data, and the simulation parameters
     * @param mappingIndex               mapping history index
     */
    void sampleMutationsGivenAncestrals(size_t mappingIndex);

    /* adds a branch mapping to the mapping in a tree format by repeatedly braking branches and adding internal nodes with single children
     * @param son                   The node at the bottom of the branch
     * @param branchMapping         The branchMapping of transitions in a MutationProcess format
     */
    void updateBranchMapping(PhyloNode* son, const MutationPath& branchMapping);

    /* sample mutations based on the stochastic mapping parameters, the source and destination state, and the branch length, and updates the simulated history along the branch in the input tree, on the fly
     * @param father                The index of the father of the node of interest
     * @param son                   The index of the node of interest
     * @param mappingIndex               mapping history index
     * @param maxIterNum            Maximal number of imulation trials
     */
    void sampleMutationsGivenAncestralsPerBranch(uint father, uint son, size_t mappingIndex, size_t maxIterNum = 1000000);

    /* converts a vector of dwelling times to a mutation path and then updates the bracnh stemming from the given node */
    /* @param node                      The node at the bottom of the branch
     * @param dwellingTimes             A vector of dwelling times where the value at each entry i corresponds to the dwelling time under the i'th state
     *                                  Note that this function generates a new tree instance, that must be deleted by the calling function.
     * @param posteriorProbabilities    Posterior probaibitlies to divide the time spent in a shared state between father and son into two transitions
     @param divMethod                 The method used in the case that the son and father share the same state (either divide the wdelling time of the staed state by 2 for  two transitions (method 0) or allocate the entire dwelling time to be adjacent to the son(method 1))
    */
    void updateBranchByDwellingTimes(PhyloNode* node, VDouble& dwellingTimes, VVDouble& posteriorProbabilities, size_t divMethod = 0);
  };
}

#endif// ___STOCHASTIC_MAPPING_H
