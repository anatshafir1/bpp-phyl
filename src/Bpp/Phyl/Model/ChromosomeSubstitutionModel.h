//
// File: CromosomeSubstitutionModel.h
// Created by: Anat Shafir
// Created on: 2020
//


#ifndef _CHROMOSOMESUBSTITUTIONMODEL_H_
#define _CHROMOSOMESUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>
#include <Bpp/Exceptions.h>
#include <regex>
//#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloatTools.h>

#define lowerBoundOfRateParam 0.0
#define lowerBoundOfExpParam -3.0
#define lowerBoundBaseNumber 3
#define upperBoundOfRateParam 100.0
#define upperBoundLinearRateParam 5.0
#define upperBoundExpParam 4.6
#define IgnoreParam -999
#define DemiEqualDupl -2
#define EPSILON 2.22045e-016
using namespace std;
namespace bpp
{
class ChromosomeSubstitutionModel;


class compositeParameter{
  public:
    enum FunctionType {CONSTANT, LINEAR, LINEAR_BD, EXP, POLYNOMIAL, LOGNORMAL, REVERSE_SIGMOID, FUNC_COUNT};
    enum ParamName {BASENUMR, LOSS, GAIN, DUPL, DEMI_DUPL, PARAMNAME_COUNT};
    typedef void (compositeParameter::*functionOp)(size_t, double*, double*);

  private:
    std::vector<Parameter*> params_;
    FunctionType func_;
    std::string name_;
    functionOp updateParamFunc_;
    size_t size_;
    int maxChrNumber_;
    //vector<double> values_;
    

  public:
    std::vector<double> getParameterValues() const;
    static bool isIgnored(compositeParameter* param){return param == 0;}
    void getBounds(size_t index, double* lowerBound, double* upperBound){return (this->*updateParamFunc_)(index, lowerBound, upperBound);}
    const FunctionType getFuncType() const {return func_;}
    //vector <double>& getValues(size_t index); 
    const std::string getName() const {return name_;}
    const size_t getSize() const {return size_;}
    static size_t getNumOfParameters(FunctionType funcType);
    //const Parameter& getParameter(size_t index){return *(params_[index]);}
    //void setNameStr(ParamName paramName);
    void getParamUpdateFunction(FunctionType funcType);
    double getRate(size_t state) const;
    //not relevant for const, linearBD, and Exp
    static compositeParameter::ParamName getCompositeRateType(int param);
    static void updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, FunctionType funcType, int maxChrNum);
    //static void updateBounds(Function* f, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, FunctionType funcType, int maxChrNum);
    static void updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound, FunctionType funcType);
    static void updateBoundsLinear(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false);
    static void updateBoundsExp(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum);
    static void updateBoundsPolynomial(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false){
      throw Exception("Not impelemented yet!");
    }
    static void updateBoundsLogNormal(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false){
      throw Exception("Not implemented yet!");
    }
    static void updateBoundsReverseSigmoid(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random = false){
      throw Exception("Not implemented yet!");
    }
    static void getBoundsForInitialParams(FunctionType func, size_t index, vector<double> paramValues, double* lowerBound, double* upperBound, int maxChrNumber, bool random = false);
    static void getAbsoluteBounds(FunctionType func, size_t index, double* lowerBound, double* upperBound, int maxChrNumber);

    static std::vector<std::string> getRelatedParameterNames(ParameterList &params, std::string pattern);

  


    compositeParameter(int &maxChromosomeNum, FunctionType func, std::string paramName, vector<Parameter*> &params):
      params_(), func_(func), name_(paramName), updateParamFunc_(0), size_(params.size()), maxChrNumber_(maxChromosomeNum)
    {
      for (size_t i = 0; i < params.size(); i++){
        params_.push_back(params[i]);
      }
      getParamUpdateFunction(func);
    }


    virtual ~compositeParameter(){}



  protected:
    //Parameter* getParameter_(size_t index){return (params_[index]);}
    void setName(ParamName name){name_ = name;}
    void setFunction(FunctionType func){func_ = func;}
    std::vector<Parameter*>& getParams(){return params_;}
    void setParams(std::vector<Parameter*> params){params_ = params;}
    
    

    // functions to update the parameters and their respective intervals
    void getConstBounds(size_t index, double* lowerBound, double* upperBound);
    void getLinearBounds(size_t index, double* lowerBound, double* upperBound);
    void getLinearBDBounds(size_t index, double* lowerBound, double* upperBound);
    void getExpBounds(size_t index, double* lowerBound, double* upperBound);
   

    void getPolynomialBounds(size_t index, double* lowerBound, double* upperBound){
      throw Exception("Not implemented yet!");
    }
    void getLogNormalBounds(size_t index, double* lowerBound, double* upperBound){
      throw Exception("Not implemented yet!");
    }
    void getReverseSigmoidBounds(size_t index, double* lowerBound, double* upperBound){
      throw Exception("Not implemented yet!");
    }
    friend ChromosomeSubstitutionModel;  
};


class ChromosomeSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED};
  enum rateChangeFunc {LINEAR = 0, EXP = 1};
  enum typeOfTransition {GAIN_T = 0, LOSS_T = 1, DUPL_T = 2, DEMIDUPL_T = 3, BASENUM_T = 4, MAXCHR_T = 5, NUMTYPES = 6, ILLEGAL = 7};
  enum paramType {BASENUM = 0, BASENUMR = 1, DUPL = 2, LOSS = 3, GAIN = 4, DEMIDUPL = 5, NUM_OF_CHR_PARAMS = 6};

private:
  compositeParameter* gain_;
  compositeParameter* loss_;
  compositeParameter* dupl_;
  compositeParameter* demiploidy_;
  int baseNum_;
  compositeParameter* baseNumR_;
  unsigned int maxChrRange_;
  rootFreqType freqType_;
  int ChrMinNum_;
  int ChrMaxNum_;
  double firstNormQ_;
  mutable bool pijtCalledFromDeriv_;
  compositeParameter::FunctionType gainFunc_;
  compositeParameter::FunctionType lossFunc_;
  compositeParameter::FunctionType duplFunc_;
  compositeParameter::FunctionType demiFunc_;
  compositeParameter::FunctionType baseNumRFunc_;
 


protected:
  mutable std::vector< RowMatrix<double> > vPowExp_;



public:
  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    vector<double> gain, 
    vector<double> loss, 
    vector<double> dupl, 
    vector<double> demi,
    int baseNum,
    vector<double> baseNumR,
    unsigned int maxChrRange, 
    rootFreqType freqType,
    vector<int> rateChangeType);

  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    std::map<int, vector<double>> mapOfParamValues,
    int baseNum,
    unsigned int maxChrRange, 
    rootFreqType freqType,
    vector<int> rateChangeType);

  //constructor with vector of parameters
  // ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
  //   vector<double> modelParams,
  //   unsigned int maxChrRange,
  //   rootFreqType freqType,
  //   rateChangeFunc rateChangeType);

  virtual ~ChromosomeSubstitutionModel() {
    delete gain_;
    delete loss_;
    delete dupl_;
    delete demiploidy_;
    delete baseNumR_;
  }
  ChromosomeSubstitutionModel(const ChromosomeSubstitutionModel& model):
    AbstractParameterAliasable(model),
    AbstractSubstitutionModel(model),
    gain_(0),
    loss_(0),
    dupl_(0),
    demiploidy_(0),
    baseNum_(model.baseNum_),
    baseNumR_(0),
    maxChrRange_(model.maxChrRange_),
    freqType_(model.freqType_),
    ChrMinNum_(model.ChrMinNum_),
    ChrMaxNum_(model.ChrMaxNum_),
    firstNormQ_(model.firstNormQ_),
    pijtCalledFromDeriv_(model.pijtCalledFromDeriv_),
    gainFunc_(model.gainFunc_),
    lossFunc_(model.lossFunc_),
    duplFunc_(model.duplFunc_),
    demiFunc_(model.demiFunc_),
    baseNumRFunc_(model.baseNumRFunc_),
    vPowExp_(model.vPowExp_)
  {
    std::vector<compositeParameter**> newModelParams = {&gain_, &loss_, &dupl_, &demiploidy_, &baseNumR_};
    std::vector<compositeParameter*> originalModelParams = {model.gain_, model.loss_, model.dupl_, model.demiploidy_, model.baseNumR_};
    for (size_t i = 0; i < compositeParameter::PARAMNAME_COUNT; i++){
      if (originalModelParams[i] == 0){
        continue;
      }
      vector<Parameter*> params = originalModelParams[i]->getParams();
      vector<Parameter*> newParams;
      for (size_t j = 0; j < params.size(); j++){
        std::string name = params[j]->getName();
        string noPrefixName = std::regex_replace(name, std::regex("Chromosome."), "");
        newParams.push_back(&(getParameter_(noPrefixName)));
      }

      *(newModelParams[i]) = new compositeParameter(ChrMaxNum_, originalModelParams[i]->getFuncType(), originalModelParams[i]->getName(), newParams);

    }

  }


  ChromosomeSubstitutionModel* clone() const { return new ChromosomeSubstitutionModel(*this);}
  

  
public:
  static ChromosomeSubstitutionModel* initRandomModel(
    const ChromosomeAlphabet* alpha,
    int &baseNumber,
    std::map<int, vector<double>> initParams,
    unsigned int chrRange,
    rootFreqType rootFrequenciesType,
    std::vector<int> rateChangeType,
    std::vector<int>& fixedParams,
    double parsimonyBound = 0);





  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;
  //The following four functions are just for test.. Finally will bw removed.
  const Matrix<double>& getPij_t_func2(double d) const;
  const Matrix<double>& getPijt_test(double d) const;
  const Matrix<double>& getPij_t_func3(double d) const;
  const Matrix<double>& getPij_t_func4(double d) const;


  std::string getName() const { return "Chromosome"; }
  void setFreq(std::map<int, double>& freqs);
  size_t getNumberOfStates() const { return size_; }
  int getMin() const {return ChrMinNum_;}
  int getMax() const {return ChrMaxNum_;}
  unsigned int getMaxChrRange() const {return maxChrRange_;}
  bool checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const;
  double getInitValue(size_t i, int state) const;
  int getBaseNumber() const {return baseNum_;}
  bool isIgnoredGain() const {return gain_ == 0;}
  bool isIgnoredLoss() const {return loss_ == 0;}
  bool isIgnoredDupl() const {return dupl_ == 0;}
  bool isIgnoredDemiDupl() const {return demiploidy_ == 0;}
  bool isIgnoredBaseNumR() const {return baseNumR_ == 0;}
  const compositeParameter* getDemiDupl() const {return demiploidy_;}
  const compositeParameter* getGain() const {return gain_;}
  const compositeParameter* getLoss() const {return loss_;}
  const compositeParameter* getDupl() const {return dupl_;}
  const compositeParameter* getBaseNumR() const {return baseNumR_;}
  //double getConstDupl () const{return dupl_;}
  //double getChangeRateDupl() const {return duplR_;}
  //double getConstGain() const {return gain_;}
  //double getChangeRateGain() const {return gainR_;}
  //double getConstLoss() const {return loss_;}
  //double getChangeRateLoss() const {return lossR_;}
  //double getBaseNumR() const {return baseNumR_;}
  //double getRate (size_t state, double constRate, double changeRate) const;
  //static void getSetOfFixedParameters(vector<double>& initParams, vector<unsigned int>& fixedParams, map<int, double>& setOfFixedParams);
  //static void getSetOfFixedParameters(vector<double>& initGain, vector<double>& initLoss, vector<double>& initDupl, vector<double>& initDemiDupl, vector<double>& initBaseNumR, int &baseNumber, vector<unsigned int>& fixedParams, map<int, double>& setOfFixedParams);

  
  //These functions should be used from chromsome number optimizer
  void checkParametersBounds() const;


protected:
  void getCompositeParametersValues(std::string paramName, compositeParameter* param);
  void calculatePijtUsingEigenValues(double t) const;
  //static void getRandomParameter(paramType type, double initParamValue, vector<double>& randomParams, double upperBound, double upperBoundLinear, double upperBoundExp, rateChangeFunc rateFunc, int maxChrNum, unsigned int chrRange, map<int, double>& setOfFixedParameters);
  //void updateParameters();
  void updateParameters(vector<double> &gain, vector<double> &loss, vector<double> &dupl, vector<double> &demi, vector<double> &baseNumR);
  //void updateLinearParameters();
  //void updateExpParameters();
  //void updateBaseNumParameters(std::shared_ptr<IntervalConstraint> interval);
  void updateMatrices();
  void updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateQWithGain(size_t currChrNum, size_t minChrNum);
  void updateQWithLoss(size_t currChrNum, size_t minChrNum);
  void updateQWithDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum = 0);
  void updateQWithDemiDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateEigenMatrices();
  void getParametersValues();
  void calculateExp_Qt(size_t pow, double s, size_t m, double v) const;
  void calculateExp_Qt(size_t pow, double* s, double v) const;
  double getFirstNorm() const;
  double get_epsilon() const{ return 0.0001;};
  //void updateConstRateParameter(double paramValueConst, double paramValueChange, string parameterName, std::shared_ptr<IntervalConstraint> interval);
  //void updateLinearChangeParameter(double paramValueConst, double paramValueChange, string parameterName);
  //void setNewBoundsForLinearParameters(double &constRate, double &changeRate, string paramNameConst, string paramNameLinear);
  std::vector<Parameter*> createCompositeParameter(compositeParameter::FunctionType &func, std::string paramName, vector<double> &vectorOfValues);
  friend class compositeParameter;
};
} // end of namespace bpp.

#endif  // _CHROMOSOMESUBSTITUTIONMODEL_H_