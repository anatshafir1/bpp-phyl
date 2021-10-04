#include "ChromosomeNumberOptimizer.h"
using namespace bpp;

void ChromosomeNumberOptimizer::initModels(std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelParams, double parsimonyBound, std::vector<int>& rateChange, int seed, unsigned int numOfPoints, const string& fixedRootFreqPath, std::map<uint, vector<int>>& fixedParams, std::map<uint, std::vector<uint>> mapModelNodesIds){
    fixedParams_ = fixedParams;
    sharedParams_ = ChromEvolOptions::sharedParameters_;
    optimizeBaseNumber_ = false;
    // if we should optimize it for at least one model, we will set it to true. Otherwise it is fixed for all the models.
    for (uint i = 1; i <= static_cast<uint>(ChromEvolOptions::numOfModels_); i++){
        if (!(std::count(fixedParams_[i].begin(), fixedParams_[i].end(), ChromosomeSubstitutionModel::BASENUM))){
            optimizeBaseNumber_ = true;
        }
    }
    vectorOfLikelohoods_.reserve(numOfPoints);
    //vectorOfContexts_.reserve(numberOfModels);

    if (seed != 0){
        RandomTools::setSeed(static_cast<long>(seed));
    }
    SingleProcessPhyloLikelihood* lik;
    for (size_t n = 0; n < numOfPoints; n++){
        if (n == 0){
            lik = setHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, ChromEvolOptions::numOfModels_);
        }else{
            lik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, ChromEvolOptions::numOfModels_, parsimonyBound * (double)n, fixedParams_);
    
        }
        
        if (std::isnan(lik->getValue())){
            std::cout << "value is nan"<<endl;
        }
        int countNumOfTrials = 0;
        
        while (((std::isinf(lik->getValue())) || (std::isnan(lik->getValue())))||(lik->getValue() < 0))
        {
            if (countNumOfTrials >= ChromEvolOptions::maxNumOfTrials_){
                break;
            }
            deleteLikObject(lik);
            lik = setRandomHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, ChromEvolOptions::numOfModels_, parsimonyBound * (double)n, fixedParams_);

            countNumOfTrials ++;

        }
        vectorOfLikelohoods_.push_back(lik);//add to vector of likelihoods
        
    }

}
// // /****************************************************************************/
// SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::getLikelihoodFunction(const PhyloTree* tree, const VectorSiteContainer* vsc, std::shared_ptr<ChromosomeSubstitutionModel> &chrModel, DiscreteDistribution* rdist, const string& fixedRootFreqPath){ 
//     bool weightedRootFreqs;
//     std::shared_ptr<SubstitutionModel> model(static_pointer_cast<SubstitutionModel>(chrModel)->clone());
//     NonHomogeneousSubstitutionProcess* subProSim;
//     ParametrizablePhyloTree parTree(*tree_);
//     if (fixedRootFreqPath != "none"){
//         vector <double> rootFreqs = setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
//         std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
//         std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
//         subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone(), shared_ptr<FrequencySet>(rootFrequencies->clone()));
//         weightedRootFreqs = false;
        
//     }else{
//         subProSim= NonHomogeneousSubstitutionProcess::createHomogeneousSubstitutionProcess(model, rdist, parTree.clone());
//         weightedRootFreqs = true;
        
//     }
    
//     SubstitutionProcess* nsubPro=subProSim->clone();
//     std::shared_ptr<Context> context = make_shared<Context>();
//     //vectorOfContexts_.push_back(context);

//     auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc_->clone(), *nsubPro, weightedRootFreqs);
//     //lik->setFactor(factor);
    
//     SingleProcessPhyloLikelihood* ntl = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
//     delete subProSim;
//     return ntl;
// }
// /****************************************************************************/
vector <double> ChromosomeNumberOptimizer::setFixedRootFrequencies(const std::string &path, std::shared_ptr<ChromosomeSubstitutionModel> chrModel){
    ifstream stream;
    stream.open(path.c_str());
    vector <double> freqs;
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    for (size_t i = 0; i < lines.size(); i++){
        string freq_i_str = TextTools::removeSurroundingWhiteSpaces(lines[i]);
        if (freq_i_str == ""){
            continue;
        }
        double freq_i = TextTools::toDouble(freq_i_str);
        if (static_cast<unsigned int>(freqs.size()) >= chrModel->getNumberOfStates()){
            if (freq_i > 0){
                throw Exception("Invalid fixed frequencies file!");
            }

        }else{
            freqs.push_back(freq_i);
        }
        
    }
    size_t nbStates = chrModel->getNumberOfStates();
    if (freqs.size() < nbStates){
        for (size_t s = freqs.size(); s < nbStates; s++){
            freqs.push_back(0);
        }
        
    }
    if (nbStates != freqs.size()){
        throw Exception("Invalid fixed frequencies file!");
    }
    return freqs;
}
/*******************************************************************************/
void ChromosomeNumberOptimizer::optimizeHeterogeneous()
{
    vector <unsigned int> baseNumCandidates;
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        uint maxBaseNumCandidate = getMaxBaseNumAmongModels(baseNumberUpperBound_);
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerBoundBaseNumber, maxBaseNumCandidate);

    }
    int maxNumOfModels;
    (ChromEvolOptions::maxNumOfModels_ == 1) ? (maxNumOfModels = (static_cast<int>((tree_->getAllLeavesNames()).size())-1)) : (maxNumOfModels = ChromEvolOptions::maxNumOfModels_);
    vector<uint> candidateShiftNodesIds;
    getValidCandidatesForShift(candidateShiftNodesIds, ChromEvolOptions::minCladeSize_);
    for (size_t i = 0; i < vectorOfLikelohoods_.size(); i++){
        optimizeSingleHeterogeneousModel(i, maxNumOfModels, candidateShiftNodesIds, baseNumCandidates);
    }
    sort(vectorOfLikelohoods_.begin(), vectorOfLikelohoods_.end(), compareLikValues);
    printRootFrequencies(vectorOfLikelohoods_[0], ChromEvolOptions::resultsPathDir_ + "//" + "inferred_rootFreq.txt");
    cout <<"*****  Final Optimized -logL *********"  <<endl;
    printLikParameters(vectorOfLikelohoods_[0], 1, ChromEvolOptions::resultsPathDir_ + "//" + "likelihood.txt");

}

// /****************************************************************************/

void ChromosomeNumberOptimizer::optimizeHomogeneous()
{

    unsigned int totalNumOfEvaluations = 0;
    unsigned int numOfEvaluations;
    unsigned int numOfEvaluationsPerCycle;
    vector <unsigned int> baseNumCandidates;

    // If base number is one of the parameters
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        uint maxBaseNumCandidate = getMaxBaseNumAmongModels(baseNumberUpperBound_);
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerBoundBaseNumber, maxBaseNumCandidate);

    }

    //Go over each cycle
    for (size_t i = 0; i < numOfIterations_.size(); i++){
        numOfEvaluationsPerCycle = 0;
        clearVectorOfLikelihoods(numOfPoints_[i]);
        cout <<"##################################" << endl;
        cout << "*********  cycle "<< i <<"  **************"<<endl;     
        //Go over each point at cycle i 
        for (size_t j = 0; j < numOfPoints_[i]; j++){
            numOfEvaluations = 0;
            std::cout << "Starting cycle with Point #" << j <<"...."<<endl;;
            printLikParameters(vectorOfLikelohoods_[j], 0);
            //If the number of optimization iterations is larger than zero, optimize the number of times as specified
            if (numOfIterations_[i] > 0){
                numOfEvaluations = optimizeModelParameters(vectorOfLikelohoods_[j], tolerance_, numOfIterations_[i], baseNumCandidates, &(ChromEvolOptions::sharedParameters_), &(ChromEvolOptions::fixedParams_));
                          
            }
            std:: cout << "Number of evaluations per point is : " << numOfEvaluations << endl;
            numOfEvaluationsPerCycle += numOfEvaluations;
            std:: cout <<"*****************************" << endl;            
        }
        totalNumOfEvaluations += numOfEvaluationsPerCycle;
        //sort the vector of likelihoods, such that the worst likelihood is at the end
        sort(vectorOfLikelohoods_.begin(), vectorOfLikelohoods_.end(), compareLikValues);
        printLikelihoodVectorValues(vectorOfLikelohoods_);
        
    }
    const string outPath = (ChromEvolOptions::resultsPathDir_ == "none") ? (ChromEvolOptions::resultsPathDir_) : (ChromEvolOptions::resultsPathDir_ + "//" + "likelihood.txt");
    const string outPathFreq = (ChromEvolOptions::resultsPathDir_ == "none") ? (ChromEvolOptions::resultsPathDir_) : (ChromEvolOptions::resultsPathDir_ + "//" + "inferred_rootFreq.txt");
    
    for (size_t i = 0; i < vectorOfLikelohoods_.size(); i++){
        printRootFrequencies(vectorOfLikelohoods_[i], outPathFreq);
        cout <<"*****  Final Optimized -logL *********"  <<endl;
        printLikParameters(vectorOfLikelohoods_[i], 1, outPath);
        std:: cout << "final number of evaluations is : " << totalNumOfEvaluations << endl;
        if (!(ChromEvolOptions::heterogeneousModel_)){
            break;
        }
    }
   
}
/********************************************************************************/
void ChromosomeNumberOptimizer::deleteLikObject(SingleProcessPhyloLikelihood* lik_to_del){
    auto sequenceData = lik_to_del->getData();
    auto process = &(lik_to_del->getSubstitutionProcess());
    //auto tree = &(lik_to_del->getTree());
    auto context = &(lik_to_del->getContext());
    delete process;
    delete sequenceData;
    //delete tree;
    delete context;
    delete lik_to_del;
    

}

// /********************************************************************************/
void ChromosomeNumberOptimizer::clearVectorOfLikelihoods(size_t new_size){
    while(vectorOfLikelohoods_.size() > new_size){
        //deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
        SingleProcessPhyloLikelihood* lik_to_del = vectorOfLikelohoods_.back(); 
        vectorOfLikelohoods_.pop_back();
        deleteLikObject(lik_to_del);

    }
}

// /***********************************************************************************/
bool ChromosomeNumberOptimizer::compareLikValues(SingleProcessPhyloLikelihood* lik1, SingleProcessPhyloLikelihood* lik2){
    return (lik1->getValue() < lik2->getValue());
}
// /***********************************************************************************/
void ChromosomeNumberOptimizer::printLikParameters(SingleProcessPhyloLikelihood* lik, unsigned int optimized, const string filePath) const{
    ofstream outFile;
    if (filePath != "none"){
        outFile.open(filePath);
    }
    if (optimized == 0){
        std:: cout << "Initial likelihood is : "<< lik->getValue() << endl;
    }else{
        std:: cout << "Optimized likelihood is : "<< lik->getValue() << endl;
        if (filePath != "none"){
            outFile << "Final optimized likelihood is: "<< lik->getValue() << endl;
        }
    }
    
    std:: cout << "Parameters are:" << endl;
    if (filePath != "none"){
        outFile << "Optimized parameters are:"<<endl;
    }
    ParameterList substitutionModelParams = lik->getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = substitutionModelParams.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i].find("Chromosome.baseNum_") != std::string::npos){
            std::cout << paramsNames[i] << "value is "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;
            if (filePath != "none"){
                outFile <<  paramsNames[i] << "value is "<< (int)(lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue()) <<endl;
            }
        }else{
            std::cout << paramsNames[i] << "value is "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;
            if (filePath != "none"){
                outFile << paramsNames[i] << "value is "<< lik->getLikelihoodCalculation()->getParameter(paramsNames[i]).getValue() <<endl;
            }
        }
        
    }
    if (filePath != "none"){
        outFile.close();
    }
    std::cout <<"***"<<endl;

}
// /*************************************************************************************/
void ChromosomeNumberOptimizer::printLikelihoodVectorValues(std::vector <SingleProcessPhyloLikelihood*> lik_vec) const{
    std :: cout <<"The likelihoods at the end of cycle are :"<<endl;
    for (size_t i = 0; i < lik_vec.size(); i++){
        std :: cout << lik_vec[i]->getValue() << endl;
    }
}

// /******************************************************************************/
void ChromosomeNumberOptimizer::printRootFrequencies(SingleProcessPhyloLikelihood* lik, const string filePath) const{
    ofstream outFile;
    if (filePath != "none"){
        outFile.open(filePath);
    }
    
    ValueRef <Eigen::RowVectorXd> rootFreqVector = lik->getLikelihoodCalculationSingleProcess()->getRootFreqs();
    for (size_t s = 0; s < (size_t)rootFreqVector->getTargetValue().size(); s ++){
        cout << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << endl;
        if (filePath != "none"){
            outFile << "F[" << s + alphabet_->getMin() << "] = " << rootFreqVector.get()->getTargetValue()[s] << endl;
        }
    }
    if (filePath != "none"){
        outFile.close();
    }

}
/**************************************************************************************/
uint ChromosomeNumberOptimizer::getMaxBaseNumAmongModels(std::map<uint, uint> baseNumberUpperBound) const{
    auto it = baseNumberUpperBound.begin();
    uint maxBaseNumBound = 0;
    while(it != baseNumberUpperBound.end()){
        if (baseNumberUpperBound[it->first] > maxBaseNumBound){
            maxBaseNumBound = baseNumberUpperBound[it->first];
        }
        it ++;
    }
    return maxBaseNumBound;
}

// /***********************************************************************************/
void ChromosomeNumberOptimizer::fillVectorOfBaseNumCandidates(std::vector <unsigned int> &baseNumCandidates, unsigned int lowerBound, unsigned int upperBound) const{
    if (baseNumOptimizationMethod_ == "Ranges"){
        getAllPossibleChrRanges(baseNumCandidates);

    }
    else if ((baseNumOptimizationMethod_ == "Sequential") || (baseNumCandidates.size() == 0)){

        for (unsigned int chr = (unsigned int)lowerBound; chr <= upperBound; chr++){
            baseNumCandidates.push_back(chr);
        }

    }

}
// /***************************************************************************************/
void ChromosomeNumberOptimizer::getAllPossibleChrRanges(std::vector <unsigned int> &baseNumCandidates) const{
    size_t numOfSequences = vsc_->getNumberOfSequences();
    unsigned int minRange = 0;
    vector <string> sequenceNames = vsc_->getSequencesNames();
    for (size_t i = 0; i < numOfSequences; i++){
        if (i == numOfSequences-1){
            continue;
        }
        BasicSequence seq1 = vsc_->getSequence(sequenceNames[i]);
        int chrNum1 = seq1.getValue(0);
        if (chrNum1 == -1){
            continue;
        }
        for (size_t j = i + 1; j < numOfSequences; j++){
            BasicSequence seq2 = vsc_->getSequence(sequenceNames[j]);
            int chrNum2 = seq2.getValue(0);
            if (chrNum2 == -1){
                continue;
            }
            unsigned int chrRange = (unsigned int)(std::abs(chrNum1 - chrNum2));
            if (chrRange == 0 || chrRange == 1){
                continue;
            }
            else if (chrRange == 2){
                continue;
            }
            if (!std::count(baseNumCandidates.begin(), baseNumCandidates.end(), chrRange)){
                if (minRange == 0){
                    minRange = chrRange;
                }else{
                    if (chrRange < minRange){
                        minRange = chrRange;
                    }
                }
                baseNumCandidates.push_back(chrRange);

            }

        }
    }
    if (minRange > 3){
        for (unsigned int i = 3; i < minRange; i++){
            baseNumCandidates.push_back(i);
        }

    }

}

// /**********************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParameters(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<uint>>* sharedParams, std::map<uint, vector<int>>* fixedParams){
    unsigned int numOfEvaluations = 0;
 
    if (typeOfOptimizer_ == "Brent"){
        numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, maxNumOfIterations, baseNumCandidates, sharedParams, fixedParams);
    }else if (typeOfOptimizer_ == "gradient"){
        checkLegalUseOfGradientOptimization();
        numOfEvaluations += optimizeMultiDimensions(tl, tol, maxNumOfIterations, sharedParams, fixedParams);

    }else{
        checkLegalUseOfGradientOptimization();
        numOfEvaluations += useMixedOptimizers(tl, tol, maxNumOfIterations, baseNumCandidates, sharedParams, fixedParams);
    }
        
    return numOfEvaluations;
    
}
// /***************************************************************************************/
void ChromosomeNumberOptimizer::checkLegalUseOfGradientOptimization(){
    // Only base num is not a composite parameter
    size_t startForComposite = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();
    for (size_t k = startForComposite; k < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; k++){
        if (ChromEvolOptions::rateChangeType_[k-startForComposite] == ChromosomeNumberDependencyFunction::CONSTANT){
            continue;
        }else if (ChromEvolOptions::rateChangeType_[k-startForComposite] == ChromosomeNumberDependencyFunction::IGNORE){
            continue;
        }else if (ChromEvolOptions::rateChangeType_[k-startForComposite] == ChromosomeNumberDependencyFunction::EXP){
            continue;
        }else{
            throw Exception ("ChromosomeNumberOptimizer::checkLegalUseOfGradientOptimization: Cannot use a gradient descent optimization! Use only Brent!!!");

        }
    }
}

// /****************************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeMultiDimensions(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::map<int, std::vector<uint>>* sharedParams, std::map<uint, vector<int>>* fixedParams, bool mixed, unsigned int currentIterNum){
    DerivableSecondOrder* f = tl;
    ParameterList tmp = tl->getSubstitutionModelParameters();
    unique_ptr<AbstractNumericalDerivative> fnum;
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0000001);
    ConjugateGradientMultiDimensions* optimizer = new ConjugateGradientMultiDimensions(fnum.get());
    fnum->setParametersToDerivate(tmp.getParameterNames());
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->getStopCondition()->setTolerance(tol* 0.1);
    optimizer->setMaximumNumberOfEvaluations(1000);
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, tl);
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();

    unsigned int numOfEvaluations = 0;
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        if(mixed){
            std::cout << "Iteration #"<< currentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<< i <<endl;
        }
        
        ParameterList paramsFull = tl->getSubstitutionModelParameters();
        std::vector <string> nonFixedparamsNames = getNonFixedParams(tl, paramsFull, fixedParams);
        ParameterList params = tl->getParameters().createSubList(nonFixedparamsNames);
        int rateParamType;
        double lowerBound;
        double upperBound;
        
        for (size_t j = 0; j < params.size(); j++){
            std::string nameOfParam = params[j].getName();
            rateParamType = paramNameAndType[nameOfParam].first;
            /////////////////////////////////////////////////////
            std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
            auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
            if (it == paramsNames.end()){
                throw Exception("ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(): index out of range!");
            }
            size_t index = it - paramsNames.begin();
            ////////////////////////////////////////////////////////
            if (rateParamType != ChromosomeSubstitutionModel::BASENUM){
                ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(ChromEvolOptions::rateChangeType_[rateParamType-startCompositeParams]);
                ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
                functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, alphabet_->getMax());
                std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(nameOfParam).getConstraint());
                interval->setLowerBound(lowerBound, interval->strictLowerBound());
                functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
                delete functionOp;

            }    

        }
        prevLikelihood = currentLikelihood;
        optimizer->init(params);
        currentLikelihood = optimizer->optimize();
        printLikParameters(tl, 1);
        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        
        
    }
    numOfEvaluations += optimizer->getNumberOfEvaluations();
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    //std::cout << "The final number of evaluations is: "<< numOfEvaluations << endl;
    delete optimizer;
    return numOfEvaluations;

}
// /*******************************************************************************/

unsigned int ChromosomeNumberOptimizer::useMixedOptimizers(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates, std::map<int, std::vector<uint>>* sharedParams, std::map<uint, vector<int>>* fixedParams){
    std::vector<size_t> optimization = RandomTools::randMultinomial(maxNumOfIterations, probsForMixedOptimization_);
    unsigned int numOfEvaluations = 0;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        double prevLikelihood = tl->getValue();
        if (optimization[i] == 0){
            std::cout << "Optimizing with Brent" <<endl;

            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, 1, baseNumCandidates, sharedParams, fixedParams, true, (unsigned int)i);
        }else{
            std::cout << "Optimizing with Gradient Descent" <<endl;
            numOfEvaluations += optimizeMultiDimensions(tl, tol, 1, sharedParams, fixedParams, true, (unsigned int)i);
        }
        double currentLikValue = tl->getValue();
        if (std::abs(prevLikelihood-currentLikValue) < tol){
            break;
        }


    }
    return numOfEvaluations;

}
// /*******************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(SingleProcessPhyloLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, std::map<int, vector<uint>>* sharedParams, std::map<uint, vector<int>>* fixedParams, bool mixed, unsigned curentIterNum){

    // Initialize optimizer
    DerivableSecondOrder* f = tl;
    BrentOneDimension* optimizer = new BrentOneDimension(f);
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    std::cout <<"max chromosome number: " << alphabet_->getMax() << endl;
    size_t startCompositeParams = ChromosomeSubstitutionModel::getNumberOfNonCompositeParams();

    // setting bracketing for Brent optimization
    if (BrentBracketing_ == 1){
        optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);

    }else if (BrentBracketing_ == 2){
        optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);
    }else{
        optimizer->setBracketing(BrentOneDimension::BRACKET_OUTWARD);
    }
    // initializing the likelihood values
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    unsigned int numOfEvaluations = 0;
    // setting maps of parameter type and the corresponding parameters, and vice versa
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, tl);
    ParameterList params;
    // starting iterations of optimization
    for (size_t i = 0; i < maxNumOfIterations; i++){
        if (mixed){
            std::cout << "Iteration #"<<curentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<<i <<endl;
        }
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
            
            if (std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), rateParamType)){
                continue;
            }
            
            //int rateCompositeParamType;
            double lowerBound;
            double upperBound;
            // param names corresponding to the parameter type
            std::vector<string> paramsNames = typeWithParamNames[rateParamType][paramNameAndType[nameOfParam].second];
            Parameter param = params.getParameter(nameOfParam);


            auto it = std::find(paramsNames.begin(), paramsNames.end(), nameOfParam);
            if (it == paramsNames.end()){
                throw Exception("ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(): index out of range!");
            }
            size_t index = it - paramsNames.begin();
            if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
                ChromosomeNumberDependencyFunction::FunctionType funcType = static_cast<ChromosomeNumberDependencyFunction::FunctionType>(ChromEvolOptions::rateChangeType_[rateParamType-startCompositeParams]);
                ChromosomeNumberDependencyFunction* functionOp = compositeParameter::setDependencyFunction(funcType);
                functionOp->updateBounds(params, paramsNames, index, &lowerBound, &upperBound, alphabet_->getMax());
                functionOp->updateBounds(f, nameOfParam, lowerBound, upperBound);
                delete functionOp;

                std::shared_ptr<IntervalConstraint> intervalFuncUpdated = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(nameOfParam).getConstraint());
                double updated_lowerBound = intervalFuncUpdated->getLowerBound();
                std::cout << "*** ***" << nameOfParam << ": Updated lower bound: " << updated_lowerBound << std::endl;


                std::shared_ptr<IntervalConstraint> intervalFuncUpdatedTL = dynamic_pointer_cast<IntervalConstraint>(tl->getParameter(nameOfParam).getConstraint());
                double updated_lowerBoundTL = intervalFuncUpdatedTL->getLowerBound();
                std::cout << "*** ***" << nameOfParam << ": Updated lower bound TL: " << updated_lowerBoundTL << std::endl;  

            }else{
                // baseNumber parameter
                if (baseNumOptimizationMethod_ != "Brent"){
                    if (!std::count((*fixedParams)[paramNameAndType[nameOfParam].second].begin(), (*fixedParams)[paramNameAndType[nameOfParam].second].end(), ChromosomeSubstitutionModel::BASENUM)){
                        optimizeBaseNum(tl, j, baseNumCandidates, &currentLikelihood, lowerBound, upperBound, nameOfParam, params, paramNameAndType[nameOfParam].second);
                        std::cout << "parameter value after optimization "<< tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue() << endl;
                        continue;
                    }
                }
            }        
                 
            cout <<"Parameter name is: "<< nameOfParam << endl; 
            if ((i == 1) & (maxNumOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(tol* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(tol);
            }
            if (rateParamType != static_cast<int>(ChromosomeSubstitutionModel::BASENUM)){
                optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
            }else{
                optimizer->setInitialInterval(lowerBound, upperBound);
            }            
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            std::cout <<"Parameter value after optimization: "<< tl->getLikelihoodCalculation()->getParameter(param.getName()).getValue() <<endl;
            std::cout << "***"<<endl;                        
        }
        printLikParameters(tl, 1);
        
        if (std::abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();
       
    }
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    delete optimizer;
    return numOfEvaluations;
}
// /*******************************************************************************/
void ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(std::map<int, std::map<uint, std::vector<string>>> &typeWithParamNames, std::map<string, std::pair<int, uint>>* paramNameAndType, SingleProcessPhyloLikelihood* tl, std::map<int, std::vector<uint>>* sharedParams){
    std::map<std::string, int> typeGeneralName;
    updateWithTypeAndCorrespondingName(typeGeneralName);
    ParameterList substitutionModelParams = tl->getSubstitutionModelParameters();
    std::vector<std::string> namesAllParams = substitutionModelParams.getParameterNames();
    std::map<std::string, int>::iterator it = typeGeneralName.begin();
    uint numOfModels = static_cast<uint>(tl->getSubstitutionProcess().getNumberOfModels());
    while (it != typeGeneralName.end()){
        string name = it->first;
        int type = it->second;
        std::map<uint, std::vector<std::string>> modelAndParameterNames = ChromosomeNumberOptimizer::getRelatedParameterNamesForEachModel(substitutionModelParams, name, numOfModels, sharedParams);
        typeWithParamNames[type] = modelAndParameterNames;
        auto paramNamesIt = modelAndParameterNames.begin();
        if (paramNameAndType == 0){
            it ++;
            continue;    
        }
        while(paramNamesIt != modelAndParameterNames.end()){
            std::vector<string> parametersNames = paramNamesIt->second;
            for (size_t i = 0; i < parametersNames.size(); i++){
                (*paramNameAndType)[parametersNames[i]] = std::pair<int, uint>(type, paramNamesIt->first);
            }
            paramNamesIt ++;
        }
        it ++;
    }

}
// /*******************************************************************************/
void ChromosomeNumberOptimizer::updateWithTypeAndCorrespondingName(std::map<std::string, int> &typeGeneralName){
    typeGeneralName["gain"] = static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    typeGeneralName["loss"] = static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    typeGeneralName["dupl"] = static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    typeGeneralName["demi"] = static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    typeGeneralName["baseNumR"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    typeGeneralName["baseNum_"] = static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    
}
/*******************************************************************************/
int ChromosomeNumberOptimizer::getEnumOfParamName(std::string pattern){

    if (pattern == "gain"){
        return static_cast<int>(ChromosomeSubstitutionModel::GAIN);
    }else if (pattern == "loss"){
        return static_cast<int>(ChromosomeSubstitutionModel::LOSS);
    }else if (pattern == "dupl"){
        return static_cast<int>(ChromosomeSubstitutionModel::DUPL);
    }else if (pattern ==  "baseNumR"){
        return static_cast<int>(ChromosomeSubstitutionModel::BASENUMR);
    }else if (pattern ==  "baseNum_"){
        return static_cast<int>(ChromosomeSubstitutionModel::BASENUM);
    }else if (pattern ==  "demi"){
        return static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL);
    
    }else{
        throw Exception("ChromosomeNumberOptimizer::getEnumOfParamName(): No such parameter exists!");
    }
    return -1;

}
/*******************************************************************************/
std::map<uint, std::pair<int, std::map<int, vector<double>>>> ChromosomeNumberOptimizer::getModelParameters(SingleProcessPhyloLikelihood* tl){
    uint numOfModels = static_cast<uint>(tl->getSubstitutionProcess().getNumberOfModels());
    std::map<uint, std::pair<int, std::map<int, vector<double>>>> modelsParams;
    std::map<std::string, int> typeGeneralName;
    updateWithTypeAndCorrespondingName(typeGeneralName);
    for (uint i = 1; i <= numOfModels; i++){
        auto params = tl->getSubstitutionProcess().getModel(static_cast<size_t>(i))->getParameters();
        std::vector<string> paramNames = params.getParameterNames();
        for (size_t j = 0; j < params.size(); j ++){
            std::string fullParamName = paramNames[j];
            auto it = typeGeneralName.begin();
            while(it != typeGeneralName.end()){
                string pattern = it->first;
                if (fullParamName.find(pattern) != string::npos){
                    auto type = typeGeneralName[pattern];
                    if (type == ChromosomeSubstitutionModel::BASENUM){
                        modelsParams[i].first = static_cast<int>(params[j].getValue());
                    }else{
                        modelsParams[i].second[type].push_back(params[j].getValue());
                    }

                    break;
                }
                it ++;
            }

        }

    }
    return modelsParams;

}
/*******************************************************************************/
std::map<uint, std::vector<string>> ChromosomeNumberOptimizer::getRelatedParameterNamesForEachModel(ParameterList &params, std::string pattern, uint numOfModels, std::map<int, vector<uint>>* mapSharedParams){
  std::map<uint, std::vector<string>> matchingParamsPerModel;
  for (uint i = 1; i <= numOfModels; i++){
    matchingParamsPerModel[i] = vector<string>();
  }
  std::vector<std::string> paramNames = params.getParameterNames();
  std::regex modelPattern ("_([\\d]+)");
  for (size_t i = 0; i < paramNames.size(); i++){
    std::string fullParamName = paramNames[i];
    if (fullParamName.find(pattern) != string::npos){
      // put to the correct model category
      std::smatch sm;
      std::regex_search(fullParamName, sm, modelPattern);
      std::string modelSuffix = sm[sm.size()-1];
      uint modelId = static_cast<uint>(stoi(modelSuffix));
      if (mapSharedParams){
        int type = getEnumOfParamName(pattern);
        auto it = mapSharedParams->find(type);
        if (it != mapSharedParams->end()){
            for (size_t j = 0; j < (*mapSharedParams)[type].size(); j ++){
                uint model_j = (*mapSharedParams)[type][j];
                matchingParamsPerModel[model_j].push_back(fullParamName);

            }
            continue;
        }
          
      }
      matchingParamsPerModel[modelId].push_back(fullParamName);
    }
  }

  return matchingParamsPerModel;

}

// /*******************************************************************************/
std::vector <string> ChromosomeNumberOptimizer::getNonFixedParams(SingleProcessPhyloLikelihood* tl, ParameterList &allParams, std::map<uint, vector<int>>* fixedParams) const{
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    std::map<string, std::pair<int, uint>> paramNameAndType; // parameter name, its type and number of model
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, &paramNameAndType, tl);
    uint numOfModels = static_cast<uint>(tl->getSubstitutionProcess().getNumberOfModels());
    vector<string> nonFixed;
    for (int i = 0; i < ChromosomeSubstitutionModel::NUM_OF_CHR_PARAMS; i++){
        auto it = typeWithParamNames.find(i);
        if (it == typeWithParamNames.end()){
            continue;
        }
        int type = it->first;
        auto modelAndParameterNames = typeWithParamNames[type];
        for(uint j = 1; j <= numOfModels; j ++){
            if (!(std::count((*fixedParams)[j].begin(), (*fixedParams)[j].end(), type))){
                vector<string> parameterNames = modelAndParameterNames[j];
                for (size_t k = 0; k < parameterNames.size(); k++){
                    nonFixed.push_back(parameterNames[k]);
                }
            }

        }
        
    }
    return nonFixed;
}

// /***************************************************************************************/
void ChromosomeNumberOptimizer::optimizeBaseNum(SingleProcessPhyloLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, double lowerBound, 
                                                double upperBound, const string &paramName, ParameterList& params, uint model){

    Function* func = tl;
    ParameterList substitutionParams = tl->getSubstitutionModelParameters();
    vector <std::string> names;
    for (size_t k = 0; k < substitutionParams.size(); k ++){
        names.push_back(substitutionParams[k].getName());
    }
    ParameterList updatedSubstitutionParams = params.createSubList(names);
    size_t best_i = (size_t)(params.getParameter(paramName).getValue());
    double f_value = *currentLikelihood;
    
    for (size_t i = 0; i < baseNumCandidates.size(); i++){
        unsigned int baseNum = baseNumCandidates[i];
        if (baseNum > baseNumberUpperBound_[model]){
            break;
        }
        params.getParameter(paramName).setValue((double)baseNum);
        double f_i = func->f(params);
        if (f_i < f_value){
            best_i = baseNum;
            f_value = f_i;
        }
    }
    params.getParameter(paramName).setValue((double)best_i);
    updatedSubstitutionParams.getParameter(paramName).setValue((double)best_i);
    func->f(updatedSubstitutionParams);
    //param.setValue((double)best_i);
    *currentLikelihood = f_value;

}
/******************************************
Functions for heterogeneous ChromEvol model
*******************************************/
double ChromosomeNumberOptimizer::calculateAICc(SingleProcessPhyloLikelihood* lik) const{
    // the number of shifts
    size_t numOfModels = lik->getSubstitutionProcess().getNumberOfModels();
    // the number of substitution params takes into account also the backward phase
    auto numOfSubstitutionParams = lik->getSubstitutionModelParameters().size();
    // N (sample size)
    auto sampleSize = tree_->getAllLeavesNames().size();
    // p (number of overall parameters)
    double numOfParams = static_cast<double>(numOfModels) + static_cast<double>(numOfSubstitutionParams);
    // sample size correction term
    double sampleSizeCorrection = (2*numOfParams*(numOfParams + 1))/(static_cast<double>(sampleSize)-numOfParams-1);
    //Calculating AICc
    double AICc = 2*(lik->getValue()) + (2*numOfParams) + sampleSizeCorrection;
    return AICc;

}
/***********************************************/


/***********************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::getNewLikObject(SingleProcessPhyloLikelihood* currentLik, uint nodeToSplit, std::map<int, std::vector<uint>>* sharedParams) const{
    SingleProcessPhyloLikelihood* newLik;
    uint modelNumInCurrentLik = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeToSplit));
    uint numOfModels = static_cast<uint>(currentLik->getSubstitutionProcess().getNumberOfModels());
    std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames;//parameter type, num of model, related parameters
    ChromosomeNumberOptimizer::updateMapsOfParamTypesAndNames(typeWithParamNames, 0, currentLik, sharedParams);
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> modelParams = getMapOfParamsForComplexModel(currentLik, typeWithParamNames, numOfModels);

    // add the new regime to the map (i.e., copy the parameters of the regime to which the node was assigned to previously)
    modelParams[numOfModels+1] = modelParams[modelNumInCurrentLik];
    //getting nodeIds per each model (for the current model) *not the new one
    std::map<uint, std::vector<uint>> mapModelNodesIds;
    getMutableMapOfModelAndNodeIds(mapModelNodesIds, currentLik);

    // getting an updated map of nodes per model (for the new model)
    PhyloTree tree = *tree_;
    vector<std::shared_ptr<PhyloNode>> newSubtree = tree.getSubtreeNodes(tree.getNode(nodeToSplit));
    vector<uint> newSubtreeIds = tree.getNodeIndexes(newSubtree);
    mapModelNodesIds[numOfModels+1] = std::vector<uint>();
    for (size_t i = 0; i < newSubtreeIds.size(); i++){
        uint nodeId = newSubtreeIds[i];
        uint prevModel = static_cast<uint>(currentLik->getSubstitutionProcess().getModelNumberForNode(nodeId));
        if (prevModel == modelNumInCurrentLik){
            mapModelNodesIds[prevModel].erase(std::remove(mapModelNodesIds[prevModel].begin(), mapModelNodesIds[prevModel].end(), nodeId), mapModelNodesIds[prevModel].end());
            mapModelNodesIds[numOfModels+1].push_back(nodeId);

        }
    }
    for (size_t i = 1; i <= numOfModels; i ++){
        auto branchProcess = currentLik->getSubstitutionProcess().getModel(i);
        baseNumberUpperBound_[static_cast<uint>(i)] = dynamic_cast<const ChromosomeSubstitutionModel*>(branchProcess)->getMaxChrRange();
    }
    baseNumberUpperBound_[numOfModels+1] = baseNumberUpperBound_[modelNumInCurrentLik];
    // setting the heterogeneous model
    newLik = setHeterogeneousModel(tree_, vsc_, alphabet_, baseNumberUpperBound_, mapModelNodesIds, modelParams, numOfModels+1);


    return newLik;

}
/****************************************************************/
std::map<uint, pair<int, std::map<int, std::vector<double>>>> ChromosomeNumberOptimizer::getMapOfParamsForComplexModel(SingleProcessPhyloLikelihood* lik, std::map<int, std::map<uint, std::vector<string>>> typeWithParamNames, uint numOfModels) {
    std::map<uint, pair<int, std::map<int, std::vector<double>>>> heterogeneousModelParams;
    for (uint i = 1; i <= numOfModels; i++){
        heterogeneousModelParams[i] = pair<int, std::map<int, std::vector<double>>>();
    }
    auto it = typeWithParamNames.begin();
    while (it != typeWithParamNames.end()){
        int type = it->first;
        auto modelIt = typeWithParamNames[type].begin();
        while(modelIt != typeWithParamNames[type].end()){
            uint model = modelIt->first;
            if (type == ChromosomeSubstitutionModel::BASENUM){
                heterogeneousModelParams[model].first = static_cast<int>(lik->getParameter(typeWithParamNames[type][model][0]).getValue());
            }else{
                vector<string> paramNames = typeWithParamNames[type][model];
                for (size_t i = 0; i < paramNames.size(); i++){
                    double paramValue = lik->getParameter(paramNames[i]).getValue();
                    heterogeneousModelParams[model].second[type].push_back(paramValue);
                }
                
            }
            modelIt ++;
        }
        it ++;
    }
    return heterogeneousModelParams;

}
/***********************************************************************************************/
void ChromosomeNumberOptimizer::getMutableMapOfModelAndNodeIds(std::map<uint, vector<uint>> &mapModelNodesIds, SingleProcessPhyloLikelihood* lik){
    uint numOfModels = static_cast<uint>(lik->getSubstitutionProcess().getNumberOfModels());
    for (uint i = 1; i <= numOfModels; i++){
        auto vectorOfNodes = lik->getSubstitutionProcess().getNodesWithModel(i);
        for (size_t j = 0; j < vectorOfNodes.size(); j++){
            mapModelNodesIds[i].push_back(vectorOfNodes[j]);
        }
    }
}
/**********************************************************************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::setRandomHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels, double parsimonyBound, std::map<uint, vector<int>> &fixedParams){
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    ParametrizablePhyloTree* parTree = new ParametrizablePhyloTree(*tree);
    string fixedRootFreqPath = ChromEvolOptions::fixedFrequenciesFilePath_;
    bool weightedRootFreqs;
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet, modelParams[1].first, modelParams[1].second, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, fixedParams[1], parsimonyBound));
    if (fixedRootFreqPath == "none"){
        weightedRootFreqs = true;
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        weightedRootFreqs = false;
        vector <double> rootFreqs = setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
        std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies.get());
    }

    
    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            //std::make_shared<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet, modelParams[1].first, modelParams[1].second, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, fixedParams_[1], parsimonyBound));
            chrModel = std::shared_ptr<ChromosomeSubstitutionModel>(ChromosomeSubstitutionModel::initRandomModel(alphabet, modelParams[i].first, modelParams[i].second, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_, fixedParams[i], parsimonyBound));
        }  
        subProSim->addModel(chrModel, mapModelNodesIds[i]);
    }


    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    return newLik;

}
/**********************************************************************************************/
SingleProcessPhyloLikelihood* ChromosomeNumberOptimizer::setHeterogeneousModel(const PhyloTree* tree, const VectorSiteContainer* vsc, const ChromosomeAlphabet* alphabet, std::map<uint, uint> baseNumberUpperBound, std::map<uint, vector<uint>> &mapModelNodesIds, std::map<uint, pair<int, std::map<int, std::vector<double>>>> &modelParams, uint numOfModels){
    DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
    ParametrizablePhyloTree* parTree = new ParametrizablePhyloTree(*tree);
    string fixedRootFreqPath = ChromEvolOptions::fixedFrequenciesFilePath_;
    bool weightedRootFreqs;
    std::shared_ptr<NonHomogeneousSubstitutionProcess> subProSim;
    std::shared_ptr<ChromosomeSubstitutionModel> chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[1].second, modelParams[1].first, baseNumberUpperBound[1], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
    if (fixedRootFreqPath == "none"){
        weightedRootFreqs = true;
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree);

    }else{
        weightedRootFreqs = false;
        vector <double> rootFreqs = setFixedRootFrequencies(ChromEvolOptions::fixedFrequenciesFilePath_, chrModel);
        std::shared_ptr<FixedFrequencySet> rootFreqsFixed = std::make_shared<FixedFrequencySet>(std::shared_ptr<const StateMap>(new CanonicalStateMap(chrModel->getStateMap(), false)), rootFreqs);
        std::shared_ptr<FrequencySet> rootFrequencies = static_pointer_cast<FrequencySet>(rootFreqsFixed);
        subProSim = std::make_shared<NonHomogeneousSubstitutionProcess>(rdist, parTree, rootFrequencies.get());
    }
    

    
    // adding models
    for (uint i = 1; i <= numOfModels; i++){
        if (i > 1){
            chrModel = std::make_shared<ChromosomeSubstitutionModel>(alphabet, modelParams[i].second, modelParams[i].first, baseNumberUpperBound[i], ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, ChromEvolOptions::rateChangeType_);
        }   
        subProSim->addModel(chrModel, mapModelNodesIds[i]);
    }


    SubstitutionProcess* nsubPro= subProSim->clone();
    Context* context = new Context();
    auto lik = std::make_shared<LikelihoodCalculationSingleProcess>(*context, *vsc->clone(), *nsubPro, weightedRootFreqs);
    SingleProcessPhyloLikelihood* newLik = new SingleProcessPhyloLikelihood(*context, lik, lik->getParameters());
    /// DEBUG /////////////////////////////////////////////////////
    std::cout << newLik->getValue() <<std::endl;
    ///////////////////////////////////////////////////////////////

    return newLik;

}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::updateSharedParameters(std::map<int, vector<uint>> &sharedParams, uint prevShift, uint numOfShifts) const{
    auto sharedParamIt = sharedParams.begin();
    while(sharedParamIt != sharedParams.end()){
        if (std::count(sharedParams[sharedParamIt->first].begin(), sharedParams[sharedParamIt->first].end(), numOfShifts + 1)){
            sharedParams[sharedParamIt->first].erase(std::remove(sharedParams[sharedParamIt->first].begin(), sharedParams[sharedParamIt->first].end(), numOfShifts + 1), sharedParams[sharedParamIt->first].end());
        }
        sharedParamIt ++;
    }
    sharedParamIt = sharedParams.begin();
    while (sharedParamIt != sharedParams.end()){
        if (std::count(sharedParams[sharedParamIt->first].begin(), sharedParams[sharedParamIt->first].end(), prevShift)){
            sharedParams[sharedParamIt->first].push_back(numOfShifts + 1);
        }
        sharedParamIt ++;
    }

}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::optimizeSingleHeterogeneousModel(size_t index, int maxNumOfModels, std::vector<uint> &candidateShiftNodesIds, vector<uint> &baseNumCandidates){
    auto lik = vectorOfLikelohoods_[index];
    vectorOfLikelohoods_[index] = 0;
    bool deltaAICcImproved = true;
    SingleProcessPhyloLikelihood* minAICcLik = lik;
    uint numOfShifts = ChromEvolOptions::numOfModels_;
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        uint maxBaseNumCandidate = getMaxBaseNumAmongModels(baseNumberUpperBound_);
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerBoundBaseNumber, maxBaseNumCandidate);

    }
    std::map<int, std::vector<uint>> sharedParams = ChromEvolOptions::sharedParameters_;
    std::map<uint, vector<int>> fixedParameters = ChromEvolOptions::fixedParams_;
    while((deltaAICcImproved) && (numOfShifts < (size_t)maxNumOfModels)){
        if (candidateShiftNodesIds.size() == 0){
            break;
        }
        uint minDetaAICcNode;
        bool improvedModelFound = false;
        double initialAICc = calculateAICc(minAICcLik);
        double minAICc = initialAICc;
        for (size_t i = 0; i < candidateShiftNodesIds.size(); i++){
            uint prevShift = static_cast<uint>(lik->getSubstitutionProcess().getModelNumberForNode(candidateShiftNodesIds[i]));
            updateSharedParameters(sharedParams, prevShift, numOfShifts);
            fixedParameters[numOfShifts + 1] = fixedParameters[prevShift];

            SingleProcessPhyloLikelihood* candidateLik = getNewLikObject(lik, candidateShiftNodesIds[i], &ChromEvolOptions::sharedParameters_);
            optimizeModelParameters(candidateLik, ChromEvolOptions::tolerance_, ChromEvolOptions::maxIterations_, baseNumCandidates, &sharedParams, &fixedParameters);
            double AICc_candidate =  calculateAICc(candidateLik);
            SingleProcessPhyloLikelihood* likToDel;
            if ((initialAICc - AICc_candidate > ChromEvolOptions::deltaAICcThreshold_) && (AICc_candidate < minAICc)){
                minDetaAICcNode = candidateShiftNodesIds[i];
                likToDel = minAICcLik;
                minAICcLik = candidateLik;
                sharedParams_ = sharedParams;
                improvedModelFound = true;                           
            }else{
                likToDel = candidateLik;
            }
            deleteLikObject(likToDel);
        }
        if (improvedModelFound){
            deltaAICcImproved = true;
            numOfShifts ++;
        }else{
            deltaAICcImproved = false;
        }   
        
        candidateShiftNodesIds.erase(std::remove(candidateShiftNodesIds.begin(), candidateShiftNodesIds.end(), minDetaAICcNode), candidateShiftNodesIds.end());

    }
    // no iteration
    vectorOfLikelohoods_[index] = minAICcLik;
    
}
/**********************************************************************************************/
void ChromosomeNumberOptimizer::getValidCandidatesForShift(std::vector<uint> &candidateShiftNodesIds, int minCladeSize){
    vector<shared_ptr<PhyloNode>> nodes = tree_->getAllNodes();
    for (size_t i = 0; i < nodes.size(); i++){
        if (tree_->isLeaf(nodes[i])){
            continue;
        }
        if (tree_->getRootIndex() == tree_->getNodeIndex(nodes[i])){
            continue;
        }
        if (std::find(ChromEvolOptions::initialModelNodes_.begin(), ChromEvolOptions::initialModelNodes_.end(), tree_->getNodeIndex(nodes[i])) != ChromEvolOptions::initialModelNodes_.end()){
            continue;
        }
        auto leavesUnderNode = tree_->getLeavesUnderNode(nodes[i]);
        if (leavesUnderNode.size() >= (size_t)minCladeSize){
            candidateShiftNodesIds.push_back(tree_->getNodeIndex(nodes[i]));
        }
    }

}

