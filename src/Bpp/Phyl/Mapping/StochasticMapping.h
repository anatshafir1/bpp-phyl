//
// File: StochasticMapping.h
// Authors:
//   Keren Halabi
// Created: 2018-06-08 00:00:00
//

/*
  Copyright or ÃÂ© or Copr. CNRS, (November 16, 2004)
  
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

#ifndef BPP_PHYL_MAPPING_STOCHASTICMAPPING_H
#define BPP_PHYL_MAPPING_STOCHASTICMAPPING_H


#include "../Likelihood/DataFlow/DataFlowCWise.h"
#include "../Likelihood/DataFlow/LikelihoodCalculationSingleProcess.h"
#include "../Simulation/MutationProcess.h"
#include "../Simulation/SubstitutionProcessSequenceSimulator.h"
#include <Bpp/Phyl/Io/Newick.h>
#include "MultiStateMappingPath.h"

// From the STL:
#include <iostream>
#include <iomanip>
#include <map>
#include <regex>


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
    map<uint, vector<size_t>> notRepresentedNodes_;     // a map of nodes that were underrepresented, because the mapping didn't match any possible evolutionary path
    size_t numOfMappingTrials_;
    std::map<uint, std::vector<size_t>>* MLAncr_;         // ML ancestors for the expected mapping
    
public:
  #define EPSILON_THRESHOLD 0.01

  explicit StochasticMapping(std::shared_ptr<LikelihoodCalculationSingleProcess> drl, size_t numOfMappings, size_t numOfMappingTrials = 1000000); // it is a good general practice to use "explicit" keyword on constructors with a single argument: https://stackoverflow.com/questions/121162/what-does-the-explicit-keyword-mean


  ~StochasticMapping();

  StochasticMapping(const StochasticMapping& sm) :
    likelihood_(sm.likelihood_),
    tree_(sm.tree_),
//   mappingParameters_(likelihood_->getSubstitutionProcess()),
    ConditionalProbabilities_(sm.ConditionalProbabilities_),
    nodesCounter_(0), numOfMappings_(sm.numOfMappings_),
    ancetralStates_(),
    mappings_(),
    jumpsProbs_(sm.jumpsProbs_),
    notRepresentedNodes_(sm.notRepresentedNodes_),
    numOfMappingTrials_(sm.numOfMappingTrials_),
    MLAncr_(sm.MLAncr_)
      
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
    /*
     * @brief set ML ancestors. This function should be used prior to finding the expected mapping
     * @param  ML ancestors. 
     */
    void setMLAncestors(std::map<uint, vector<size_t>>* ancestors){
      MLAncr_ = ancestors;

    }
    std::map<uint, vector<size_t>> getMLAncestors(){
      return *MLAncr_;
    }
    std::shared_ptr<PhyloTree> createExpectedMappingHistory(size_t mappingsNum){
      //generateStochasticMapping();
      auto expected_mapping = generateExpectedMapping();
      return expected_mapping;
    }
    



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
  
    std::shared_ptr<PhyloTree> generateExpectedMapping();

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
    int getNodeState(std::shared_ptr<PhyloNode> node) const{
      auto nodeName = node->getName();
      std::smatch match_state;
      std::regex state_rgx("-([\\d]+)");
      regex_search(nodeName, match_state, state_rgx);
      int state = stoi(match_state[1]);
      return state;
    }
    /**
     *@brief Gets the dwelling times for each state in a given mapping. Note: the function assums thay the 
     *        vector dwellingTimes is already resized according to the number of states!
     *
     *@param dwellingTimes A vector where the dewelling times for each state will be stored
     *@param mappingIndex The index of the mapping
     */
    void getDewellingTimesUnderEachStatePerMapping(vector<double> *dwellingTimes, size_t mappingIndex);
    /**
     *@brief Gets the expected dwelling times for each state for each branch.
     *@param dwellingTimes A map of nodeIds and their corresponding expected dewelling times for each state
     */
    void getDewellingTimesUnderEachStatePerNode(std::map<uint, vector<double>> *dwellingTimes);
    /**
     *@brief Gets the of ocurrences of each transition in a given mapping. 
     *@param mappingIndex The index of the mapping
     *@param transitionOcurrences The output map of occurrences for each transition
     */
    void getNumOfOcuurencesForEachTransitionPerMapping(size_t mappingIndex, std::map<uint, std::map<std::pair<size_t, size_t>, double>> &transitionOcurrences);
    /**
     *@brief Gets the expected dwelling time under each states given all the mappings
     *@param expectedDuration     true if we want to calculated the expected tim duration under each state. False
     * if we want to get the total duration of time under each state in all the mappings together
     *@return A vector that contains the expected dwelling times under each state
    */
    Vdouble getDwellingTimesUnderEachState(bool expectedDuration = false);
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
     *@param dwellingTimesPerState  the total time spent under state j.
     *@param numOfOccurencesPerTransition  for each transition, the total number of its occurrences.
     *@return The rates of each possible transition under each state rates. For example,
     * the [i][j] value will represent the rate of (i,j) transition given that we are in state i.
     * These rates are integrated over all the mappings
    */
    VVdouble getExpectedRateOfTransitionGivenState(Vdouble &dwellingTimesPerState, std::map<std::pair<size_t, size_t>, double> &numOfOccurencesPerTransition);
    /**
     *@brief Gets the total dwelling time under each state for each mapping
     *@return The total dwelling time under each state for each mapping. The first index represents the 
     * mapping index, while the second index represents the state index.
    */

    VVdouble getDwellingTimeOfStatePerEachMapping();


    /**
     *@brief Prints the unrepresented leaves
    */
   void printUnrepresentedLeavesWithCorrespondingMappings(ofstream &stream);
   std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> getNumOfOccurrencesFromRootToNode(std::map<uint, std::map<size_t, bool>> &presentMapping);
   void updateFromRootToNodeRecursively(std::map<uint, std::map<size_t, bool>> &presentMapping, std::map<uint, std::map<pair<size_t, size_t>, double>> &occurrencesPerMapping, size_t mappingIndex, uint nodeId, std::map<uint, std::map<size_t, std::map<std::pair<size_t, size_t>, double>>> &occurrencesFromRootToLeaf);
    /*
    *@brief compute posterior probability for each node and state
    */
   void getPosteriorProbabilities(std::map<uint, std::vector<double>> &ancestralStatesFreqs);

    /*
    * get mappings (not const)
    */
   const map<uint, vector<MutationPath>> getMappings() const{
     return mappings_;
   }
   size_t getNumberOfMappings() const{
    return numOfMappings_;
   }
   const std::shared_ptr<PhyloTree> getTree() const{
    return tree_;
   }
   const std::map<uint, std::vector<size_t>> getAncestralStates() const{
    return ancetralStates_;
   }


   map<uint, vector<size_t>> getFailedNodes(){
     return notRepresentedNodes_;
   }
   void removeFailedNodes(uint nodeId, size_t mappingIndex){
     auto &failedMappings = notRepresentedNodes_[nodeId];
     failedMappings.erase(std::remove(failedMappings.begin(), failedMappings.end(), mappingIndex), failedMappings.end());
     if (failedMappings.size() == 0){
       auto it = notRepresentedNodes_.find(nodeId);    
       notRepresentedNodes_.erase(it);
     }
   }
    /*
    * try to fix a mapping for a given node
    */
   bool tryToReplaceMapping(double branchLength, uint nodeId, size_t mappingIndex, vector<MutationPath> &mappings, size_t maxNumOfIterations);
   double getRateToLeaveState(uint nodeId, size_t mapping);
    /**
     *@brief create a tree from a mapping
     *@param mapping index
     *@return A tree with transitions along a branch represented by internal nodes
    */ 
   std::shared_ptr<PhyloTree> createMappingHistoryTree(size_t mappingIndex) const;
    /**
     *@brief create a trees from all stochastic mappings
     *@return A vector of trees with transitions along a branch represented by internal nodes
    */ 
   std::vector<std::shared_ptr<PhyloTree>> createMappingHistoryTrees() const;

  private:
  void assignDewellingTimesUnderEachStatePerMappingPerBranch(uint nodeId, size_t initialState, vector<double> &dwellingTimes, MutationPath &mutationPath);
    void sampleAllAncestals();
    bool sampleEvolutionaryPathForBranch(size_t sonState, size_t fatherState, uint father, uint son, double branchLength, size_t mappingIndex, vector<MutationPath>& mappings, size_t maxIterNum, bool replace = false);
    bool isAccounted(uint nodeId, size_t mappingIndex);
    void clearMapping(size_t mappingIndex);
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
    void getDewellingTimesUnderEachStatePerMappingRecursively(uint nodeId, size_t initialState, vector<double> *dwellingTimes, size_t mappingIndex, std::map<uint, std::vector<double>> *dwellingTimesPerNode=0);

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
    
    // /* sets the state of a node in a mapping
    //  * @param node               The node to get the state of
    //  * @param state              The state that needs to be assigned to the node
    //  */
    // void setNodeState(PhyloNode* node, size_t state);

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
     * @param                     A map where the key is node id, and the value is a vector of state frequencies, i.e.,
     *                            the size of the vector is in the size of number of states.
     */
    void computeStatesFrequencies(std::map<uint, std::vector<double>> &ancestralStatesFreqs);

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
    bool setExpectedAncestrals(shared_ptr<PhyloTree> expectedMapping, std::map<uint, std::vector<double>> &ancestralStatesFrequencies);

    /* simulates mutations on phylogeny based the sampled ancestrals, tips data, and the simulation parameters
     * @param mappingIndex               mapping history index
     * @param failedNodes                 the node ids of the nodes for which the mapping has failed
     * @return                           true if mapping succeeded. False otherwise
     */
    bool sampleMutationsGivenAncestrals(size_t mappingIndex, vector<uint>* failedNode = 0);

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
     * @return: true if the mapping was successful. Otherwise, false.
     */
    bool sampleMutationsGivenAncestralsPerBranch(uint father, uint son, size_t mappingIndex, vector<MutationPath> &mappings, size_t maxIterNum = 10000);

    /* converts a vector of dwelling times to a mutation path and then updates the bracnh stemming from the given node */
    /* @param node                      The node at the bottom of the branch
     * @param dwellingTimes             A vector of dwelling times where the value at each entry i corresponds to the dwelling time under the i'th state
     *                                  Note that this function generates a new tree instance, that must be deleted by the calling function.
     * @param posteriorProbabilities    Posterior probaibitlies to divide the time spent in a shared state between father and son into two transitions
     @param divMethod                 The method used in the case that the son and father share the same state (either divide the wdelling time of the staed state by 2 for  two transitions (method 0) or allocate the entire dwelling time to be adjacent to the son(method 1))
    */
    void updateBranchByDwellingTimes(PhyloNode* node, VDouble& dwellingTimes, VVDouble& posteriorProbabilities, size_t divMethod = 0);
    /* Recursively creates a tree of a mapping history, where each transition is represented by an additional node in the tree.
     * @param nodeId               Node Index for which we do the calculation recursively
     * @param mappingIndex         The index of the mapping
     * @param tree          A pointer to the constructed tree
    */
    void assignTransitionOnHistoryTreeRec(uint nodeId, size_t mappingIndex, std::shared_ptr<PhyloTree> tree) const;
    /*
    Get the expected number of transitions given the expected ancestral states. This function is needed in order to obtain the expected mapping
    history for a multi-state trait
    * @param nodeId         Current node for whcih we want to get the expected number of transitions along the branch
    * @param fatherId       Father id of the current node
    * @param startState     The desired start state of the branch (the expected ancestral state)
    * @param endState       The desired end state of the branch (the expected ancestral state)
    * @param transitionOcurrences   The map which should store the expected number of transitions for each possible transition given the termianl states, per each node
    */
    void getExpectedNumberOfTransitionsPerBranchGivenTerminals(uint nodeId, uint fatherId, size_t startState, size_t endState, std::map<pair<size_t, size_t>, double> &transitionOcurrences, std::map<pair<size_t, size_t>, double> &timeDurations, vector<size_t> &mostFreqPath, vector<MutationPath> &mappings);
    void getExpectedNumberOfTransitionsPerGivenTermianls(std::shared_ptr<PhyloTree> expectedTree, std::map<uint, std::map<pair<size_t, size_t>, double>> &transitionOcurrences, std::map<uint, std::map<pair<size_t, size_t>, double>> &timeDurations, std::unordered_map<uint, vector<size_t>> &mostFreqPaths);

    void findTransitionsAndTimeDurationsForBinary(std::shared_ptr<PhyloTree> expectedMapping, std::map<uint, std::vector<double>> &dwellingTimes, std::map<uint, std::vector<double>> &ancestralStatesFrequencies);
    void findExpectedHistoryTransitionsAndTimeDurationsMultiState(std::shared_ptr<PhyloTree> expectedMapping, std::map<uint, std::vector<double>> &dwellingTimes);
    void getTimeDurationsPerStateGivenAncestrals(std::map<uint, std::map<pair<size_t, size_t>, double>> &timeDurations, std::map<uint, std::vector<double>> &timeDurationsPerState);
    void stringToVector(const std::string& str, vector<size_t> &res);
    void findExpectedPathOnBranch(std::shared_ptr<PhyloTree> expectedMapping, size_t fatherState, size_t sonState, uint fatherId, uint nodeId, std::map<pair<size_t, size_t>, double> &transitionOcurrences, std::map<pair<size_t, size_t>, double> &timeDurations, vector<double> &timeDurationsPerState, vector<MutationPath> &nodeMappings, vector<size_t> &mostFreqPath);



  };
}
#endif // BPP_PHYL_MAPPING_STOCHASTICMAPPING_H
