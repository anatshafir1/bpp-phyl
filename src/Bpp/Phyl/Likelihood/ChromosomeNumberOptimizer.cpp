#include "ChromosomeNumberOptimizer.h"
using namespace bpp;

void ChromosomeNumberOptimizer::initModels(vector<double> modelParams, double parsimonyBound, ChromosomeSubstitutionModel::rateChangeFunc rateChange, int seed, unsigned int numberOfModels, bool calculateDerivatives, const string& fixedRootFreqPath, vector<unsigned int>& fixedParams){
    //optimizeBaseNumber_ = optimizeBaseNumber;
    fixedParams_ = fixedParams;
    map <int, double> setOfFixedParams;
    ChromosomeSubstitutionModel::getSetOfFixedParameters(modelParams, fixedParams_, setOfFixedParams);
    optimizeBaseNumber_ = setOfFixedParams.count(ChromosomeSubstitutionModel::BASENUM) == 0;
    vectorOfLikelohoods_.reserve(numberOfModels);
    vector <int> nodeIds = tree_->getNodesId();
    nodeIds.pop_back();
    if (seed != 0){
        RandomTools::setSeed(static_cast<long>(seed));
    }
    for (size_t n = 0; n < numberOfModels; n++){
        DiscreteDistribution* rdist = new GammaDiscreteRateDistribution(1, 1.0);
        SubstitutionModelSet* modelSet = new SubstitutionModelSet(alphabet_);
        ChromosomeSubstitutionModel* chrModel;
        if (n == 0){
            chrModel = new ChromosomeSubstitutionModel(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::rootFreqType::ROOT_LL, rateChange);//initModel(alpha, chrRange);
        }else{
            //chrModel = initRandomModel(alpha, chrRange, parsimonyBound * (double)n);
            chrModel = ChromosomeSubstitutionModel::initRandomModel(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::ROOT_LL, rateChange, fixedParams_, parsimonyBound * (double)n);
    
        }       
        modelSet->addModel(chrModel, nodeIds);
        DRNonHomogeneousTreeLikelihood lik = getLikelihoodFunction(tree_, vsc_, modelSet, rdist, calculateDerivatives, fixedRootFreqPath);
        
        lik.initialize();
        if (std::isnan(lik.getValue())){
            std::cout << "value is nan"<<endl;
        }
        
        while (std::isinf(lik.getValue()))
        {
            deleteTreeLikAssociatedAttributes(lik);
            rdist = new GammaDiscreteRateDistribution(1, 1.0);
            modelSet = new SubstitutionModelSet(alphabet_);
            chrModel = ChromosomeSubstitutionModel::initRandomModel(alphabet_, modelParams, baseNumberUpperBound_, ChromosomeSubstitutionModel::ROOT_LL, rateChange, fixedParams_, parsimonyBound * (double)n);
            //chrModel = initRandomModel(alpha, chrRange, parsimonyBound * (double)n);
            modelSet->addModel(chrModel, nodeIds);
            lik = getLikelihoodFunction(tree_, vsc_, modelSet, rdist, calculateDerivatives, fixedRootFreqPath);
            
            lik.initialize();

        }
            
        

        //initializing the likelihood instance
        //lik.initialize();
        vectorOfLikelohoods_.push_back(lik);//add to vector of likelihoods
        
    }

}
/****************************************************************************/
DRNonHomogeneousTreeLikelihood ChromosomeNumberOptimizer::getLikelihoodFunction(const TreeTemplate<Node>* tree, const VectorSiteContainer* vsc, SubstitutionModelSet* modelSet, DiscreteDistribution* rdist, bool calculateDerivatives, const string& fixedRootFreqPath) const{
    // bool calculateDerivatives = true;
    // if (ChromEvolOptions::optimizationMethod_ == "Brent"){
    //     calculateDerivatives  = false;
    // }
    if (fixedRootFreqPath != "none"){
        setFixedRootFrequencies(fixedRootFreqPath, modelSet);
        DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, false, calculateDerivatives, modelSet, rdist);
        return lik;
    }else{
        DRNonHomogeneousTreeLikelihood lik = DRNonHomogeneousTreeLikelihood(*tree, *vsc, true, calculateDerivatives, modelSet, rdist);
        return lik;

    }
}
/****************************************************************************/
void ChromosomeNumberOptimizer::setFixedRootFrequencies(const std::string &path, SubstitutionModelSet* modelSet){
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
        if (static_cast<unsigned int>(freqs.size()) >= modelSet->getNumberOfStates()){
            if (freq_i > 0){
                throw Exception("Invalid fixed frequencies file!");
            }

        }else{
            freqs.push_back(freq_i);
        }
        
    }
    size_t nbStates = modelSet->getNumberOfStates();
    if (freqs.size() < nbStates){
        for (size_t s = freqs.size(); s < nbStates; s++){
            freqs.push_back(0);
        }
        
    }
    if (nbStates != freqs.size()){
        throw Exception("Invalid fixed frequencies file!");
    }
    FixedFrequencySet* rootFreqs = new FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(modelSet->getStateMap(), false)), freqs);
    modelSet->setRootFrequencies(static_cast<FrequencySet*>(rootFreqs));
   

}

/****************************************************************************/

void ChromosomeNumberOptimizer::optimize()
{

    unsigned int totalNumOfEvaluations = 0;
    unsigned int numOfEvaluations;
    unsigned int numOfEvaluationsPerCycle;
    vector <unsigned int> baseNumCandidates;

    // If base number is one of the parameters
    if ((baseNumOptimizationMethod_ != "Brent") && (optimizeBaseNumber_)){
        fillVectorOfBaseNumCandidates(baseNumCandidates, lowerBoundBaseNumber, baseNumberUpperBound_);

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
                numOfEvaluations = optimizeModelParameters(&vectorOfLikelohoods_[j], tolerance_, numOfIterations_[i], baseNumCandidates);
                          
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
   
    printRootFrequencies(vectorOfLikelohoods_[0]);
    cout <<"*****  Final Optimized -logL *********"  <<endl;
    printLikParameters(vectorOfLikelohoods_[0], 1);
    std:: cout << "final number of evaluations is : " << totalNumOfEvaluations << endl;

}

/********************************************************************************/
void ChromosomeNumberOptimizer::clearVectorOfLikelihoods(size_t new_size){
    while(vectorOfLikelohoods_.size() > new_size){
        deleteTreeLikAssociatedAttributes(vectorOfLikelohoods_[vectorOfLikelohoods_.size()-1]);
        vectorOfLikelohoods_.pop_back();        
    }
}
/*********************************************************************************/
void ChromosomeNumberOptimizer::deleteTreeLikAssociatedAttributes(DRNonHomogeneousTreeLikelihood &lik){
    const SubstitutionModelSet* modelSet = lik.getSubstitutionModelSet();
    const DiscreteDistribution* rateDist = lik.getRateDistribution();
    delete modelSet;
    delete rateDist;
}
/***********************************************************************************/
bool ChromosomeNumberOptimizer::compareLikValues(DRNonHomogeneousTreeLikelihood &lik1, DRNonHomogeneousTreeLikelihood &lik2){
    return (lik1.getValue() < lik2.getValue());
}
/***********************************************************************************/
void ChromosomeNumberOptimizer::printLikParameters(DRNonHomogeneousTreeLikelihood &lik, unsigned int optimized) const{
    //double res = lik.getLikelihood();
    if (optimized == 0){
        std:: cout << "Initial likelihood is : "<< lik.getValue() << endl;
    }else{
        std:: cout << "Optimized likelihood is : "<< lik.getValue() << endl;
    }
    
    std:: cout << "Parameters are:" << endl;
    ParameterList params = lik.getSubstitutionModelParameters();
    std::vector<std::string> paramsNames = params.getParameterNames();
    for (int i = 0; i < (int)(paramsNames.size()); i++){
        if (paramsNames[i] == "Chromosome.baseNum_1"){
            std::cout << paramsNames[i] << "value is "<< (int)(params.getParameterValue(paramsNames[i]))<<endl;
        }else{
            std::cout << paramsNames[i] << "value is "<< params.getParameterValue(paramsNames[i])<<endl;
        }
        
    }
    std::cout <<"***"<<endl;

}
/*************************************************************************************/
void ChromosomeNumberOptimizer::printLikelihoodVectorValues(std::vector <DRNonHomogeneousTreeLikelihood> lik_vec) const{
    std :: cout <<"The likelihoods at the end of cycle are :"<<endl;
    for (size_t i = 0; i < lik_vec.size(); i++){
        std :: cout << lik_vec[i].getValue() << endl;
    }
}

/******************************************************************************/
void ChromosomeNumberOptimizer::printRootFrequencies(DRNonHomogeneousTreeLikelihood &lik) const{
    std :: vector <double> rootFreq = lik.getRootFrequencies(0);
    for (int i = 0; i < (int)(rootFreq.size()); i++){
        std :: cout << "F["<< ((int)i + alphabet_->getMin()) << "] = "<< rootFreq[i] <<endl;
    }

}

/***********************************************************************************/
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
/***************************************************************************************/
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
            unsigned int chrRange = (unsigned int)(abs(chrNum1 - chrNum2));
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
/**********************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParameters(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates){
    unsigned int numOfEvaluations = 0;
    if (standardOptimization_){
        double prevLogLik;
        ParameterList params;
        for (size_t i = 0; i < maxNumOfIterations; i++){
            std::cout << "Iteration #"<<i <<endl;
            prevLogLik = tl->getValue();
            params = tl->getSubstitutionModelParameters();
            numOfEvaluations += OptimizationTools::optimizeNumericalParameters(tl, params, 0, 1, tol, 2, ApplicationTools::message.get(), ApplicationTools::message.get(), false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT, (unsigned int)(BrentBracketing_));
            printLikParameters(*tl, 1);
            if (abs(tl->getValue() - prevLogLik) < tol){
                break;
            }
        }
        std::cout <<"..."<<endl;
    }else{
        if (typeOfOptimizer_ == "Brent"){
            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, maxNumOfIterations, baseNumCandidates);
        }else if (typeOfOptimizer_ == "gradient"){

            numOfEvaluations += optimizeMultiDimensions(tl, tol, maxNumOfIterations);

        }else{
            numOfEvaluations += useMixedOptimizers(tl, tol, maxNumOfIterations, baseNumCandidates);
        }
        
    }
    
    return numOfEvaluations;
    
}
/****************************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeMultiDimensions(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, bool mixed, unsigned int currentIterNum){
    DerivableSecondOrder* f = tl;
    unique_ptr<AbstractNumericalDerivative> fnum;
    fnum.reset(new TwoPointsNumericalDerivative(f));
    fnum->setInterval(0.0000001);
    ConjugateGradientMultiDimensions* optimizer = new ConjugateGradientMultiDimensions(fnum.get());
    ParameterList tmp = tl->getNonDerivableParameters();
    fnum->setParametersToDerivate(tmp.getParameterNames());
    optimizer->setVerbose(1);
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->getStopCondition()->setTolerance(tol* 0.1);
    optimizer->setMaximumNumberOfEvaluations(1000);
    optimizer->setBrentOptimizer(static_cast<BrentOneDimension::Bracketing>(BrentBracketing_));

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
        std::vector <string> paramsNames = getNonFixedParams(fixedParams_, paramsFull);
        ParameterList params = paramsFull.createSubList(paramsNames);
        std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBoundOfRateParam + 0.0000000001, upperBoundOfRateParam, true, true);
        for (size_t j = 0; j < params.size(); j++){
            params[j].setConstraint(interval);
        }
        prevLikelihood = currentLikelihood;
        optimizer->init(params);
        currentLikelihood = optimizer->optimize();
        printLikParameters(*tl, 1);
        if (abs(prevLikelihood-currentLikelihood) < tol){
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
/*******************************************************************************/

unsigned int ChromosomeNumberOptimizer::useMixedOptimizers(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, vector <unsigned int> &baseNumCandidates){
    std::vector<size_t> optimization = RandomTools::randMultinomial(maxNumOfIterations, probsForMixedOptimization_);
    unsigned int numOfEvaluations = 0;
    for (size_t i = 0; i < maxNumOfIterations; i++){
        double prevLikelihood = tl->getValue();
        if (optimization[i] == 0){
            std::cout << "Optimizing with Brent" <<endl;
            numOfEvaluations += optimizeModelParametersOneDimension(tl, tol, 1, baseNumCandidates, true, (unsigned int)i);
        }else{
            std::cout << "Optimizing with Gradient Descent" <<endl;
            numOfEvaluations += optimizeMultiDimensions(tl, tol, 1, true, (unsigned int)i);
        }
        double currentLikValue = tl->getValue();
        if (abs(prevLikelihood-currentLikValue) < tol){
            break;
        }


    }
    return numOfEvaluations;

}
/*******************************************************************************/
unsigned int ChromosomeNumberOptimizer::optimizeModelParametersOneDimension(DRNonHomogeneousTreeLikelihood* tl, double tol, unsigned int maxNumOfIterations, std::vector <unsigned int> &baseNumCandidates, bool mixed, unsigned curentIterNum){

    // Initialize optimizer
    BrentOneDimension* optimizer = new BrentOneDimension(tl);
    optimizer->setVerbose(1);
    //optimizer->setProfiler(ApplicationTools::message.get());
    //optimizer->setMessageHandler(ApplicationTools::message.get());
    optimizer->setProfiler(0);
    optimizer->setMessageHandler(0);
    optimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    optimizer->setMaximumNumberOfEvaluations(100);
    std::cout <<"max chromosome number: " << alphabet_->getMax() << endl;
    if (BrentBracketing_ == 1){
        optimizer->setBracketing(BrentOneDimension::BRACKET_INWARD);

    }else if (BrentBracketing_ == 2){
        optimizer->setBracketing(BrentOneDimension::BRACKET_SIMPLE);
    }else{
        optimizer->setBracketing(BrentOneDimension::BRACKET_OUTWARD);
    }
    
    double currentLikelihood = tl->getValue();
    double prevLikelihood;
    unsigned int numOfEvaluations = 0;
    unsigned int numOfBaseNumEval = 0;

    for (size_t i = 0; i < maxNumOfIterations; i++){
        if (mixed){
            std::cout << "Iteration #"<<curentIterNum <<endl;

        }else{
            std::cout << "Iteration #"<<i <<endl;
        }
        ParameterList paramsCopied = tl->getSubstitutionModelParameters();
        const ParameterList params = optimizer->getFunction()->getParameters();
        size_t nbParams = paramsCopied.size();
        prevLikelihood = currentLikelihood;
        unsigned int numOfLikEvaluations;
        
        for (size_t j = 0; j < nbParams; j ++){
            numOfLikEvaluations = tl->getNumberOfLikelihoodEvaluations();
            //ParameterList params = tl->getSubstitutionModelParameters();
            if (fixedParams_[j]){
                continue;
            }
            const string nameOfParam = paramsCopied[j].getName();
            Parameter param = params.getParameter(nameOfParam);
            std::cout << "Parameter name is: "<< nameOfParam <<endl;
            string paramNameInModel = findParameterNameInModel(nameOfParam);
            const ChromosomeSubstitutionModel* model = dynamic_cast <const ChromosomeSubstitutionModel*>(tl->getSubstitutionModelSet()->getModel(0));
            model->setBoundsForEquivalentParameter(param, paramNameInModel);
            model->checkParametersBounds();
            double lowerBound = dynamic_pointer_cast<IntervalConstraint>(param.getConstraint())->getLowerBound();
            double upperBound = dynamic_pointer_cast<IntervalConstraint>(param.getConstraint())->getUpperBound();
            //setNewBounds(params, params[j], &lowerBound);
 
            //model->checkParametersBounds();
            if ((baseNumOptimizationMethod_ != "Brent") && (param.getName() == "Chromosome.baseNum_1")){
                if (optimizeBaseNumber_){
                    optimizeBaseNum(tl, j, baseNumCandidates, &currentLikelihood, &numOfBaseNumEval, &numOfLikEvaluations, lowerBound, upperBound);
                    std::cout << "parameter value after optimization "<< param.getValue() << endl;
                    continue;

                }
            }           
            if ((i == 1) & (maxNumOfIterations > 2)){
                optimizer->getStopCondition()->setTolerance(tol* 2);
            }else{
                optimizer->getStopCondition()->setTolerance(tol);
            }
            if (param.getName() != "Chromosome.baseNum_1"){
                optimizer->setInitialInterval(lowerBound + 1e-10, upperBound);
            }else{
                optimizer->setInitialInterval(lowerBound, upperBound);
            }            
            optimizer->init(params.createSubList(param.getName()));
            currentLikelihood = optimizer->optimize();
            std::cout <<"Parameter value after optimization: "<< tl->getSubstitutionModelParameters()[j].getValue()<<endl;
            std::cout << "***"<<endl;
            numOfLikEvaluations = tl->getNumberOfLikelihoodEvaluations() - numOfLikEvaluations;
            if (param.getName() == "Chromosome.baseNum_1"){
                numOfBaseNumEval += numOfLikEvaluations;
            }
                        
        }
        printLikParameters(*tl, 1);
        
        if (abs(prevLikelihood-currentLikelihood) < tol){
            break;
        }
        numOfEvaluations += optimizer->getNumberOfEvaluations();
       
    }
    std::cout << "Number of likelihood evaluations per parameter is "<< numOfBaseNumEval << endl;
    if (!mixed){
        std::cout <<"..."<<endl;
    }
    delete optimizer;
    return numOfEvaluations;
}
/*******************************************************************************/
std::vector <string> ChromosomeNumberOptimizer::getNonFixedParams(std::vector <unsigned int> fixedParams, ParameterList &allParams) const{
    std:: vector <string> paramsNames;
    for (size_t i = 0; i < allParams.size(); i++){
        if (!(fixedParams[i])){
            paramsNames.push_back(allParams[i].getName());
        }

    }
    return paramsNames;
}
/***********************************************************************************/
string ChromosomeNumberOptimizer::findParameterNameInModel(string fullParameterName) const{
    StringTokenizer st(fullParameterName, "._", false, false);
    string paramName;
    while(st.hasMoreToken()){
        string token = st.nextToken();
        if (token == "Chromosome"){
            continue;
        }else if (token == "1"){
            continue;
        }else{
            paramName = token;
            break;
        }
    }
    return paramName;
}
/***************************************************************************************/
void ChromosomeNumberOptimizer::optimizeBaseNum(DRNonHomogeneousTreeLikelihood* tl, size_t index, std::vector <unsigned int> baseNumCandidates, double* currentLikelihood, unsigned int* numOfBaseNumEval, unsigned int* numOfLikEvaluations, double lowerBound, double upperBound){

    *numOfBaseNumEval += (unsigned int)baseNumCandidates.size();
    Function* func = tl;
    ParameterList params = tl->getSubstitutionModelParameters();
    size_t best_i = (size_t)(params[index].getValue());
    //double f_value = func->f(params);
    double f_value = *currentLikelihood;
    
    for (size_t i = 0; i < baseNumCandidates.size(); i++){
        unsigned int baseNum = baseNumCandidates[i];
        params[index].setValue((double)baseNum);
        //ParameterList params = tl->getSubstitutionModelParameters();
        double f_i = func->f(params);
        if (f_i < f_value){
            best_i = baseNum;
            f_value = f_i;
        }
    }
    params[index].setValue((double)best_i);
    func->f(params);
    //param.setValue((double)best_i);
    *currentLikelihood = f_value;

}