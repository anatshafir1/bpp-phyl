//
// File: ChromosomeSubstitutionModel.cpp
// Created by: Anat Shafir
// Created on: 2020
//


#include "ChromosomeSubstitutionModel.h"

// From the STL:
#include <cmath>

using namespace bpp;

#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

using namespace std;
/******************************************************************************/
// void compositeParameter::setNameStr(ParamName paramName){
//   for (size_t i = 0; i < static_cast<int>(PARAMNAME_COUNT); i++){
//     switch(i){
//       case ParamName::BASENUMR:
//         nameStr_ = "baseNumR";
//         break;
//       case ParamName::DEMI_DUPL:
//         nameStr_ = "demi";
//         break;
//       case ParamName::DUPL:
//         nameStr_ = "dupl";
//         break;
//       case ParamName::GAIN:
//         nameStr_ = "gain";
//         break;
//       case ParamName::LOSS:
//         nameStr_ = "loss";
//         break;
//       default:
//         throw Exception("compositeParameter::setNameStr(): Invalid rate type!");
//         break;
//     }
//   }
// } 
/******************************************************************************/
void compositeParameter::getParamUpdateFunction(FunctionType funcType){

  switch(funcType){
    case FunctionType::CONSTANT:
      updateParamFunc_ = &compositeParameter::getConstBounds;
      break;
    case FunctionType::LINEAR:
      updateParamFunc_ = &compositeParameter::getLinearBounds;
      break;
    case FunctionType::LINEAR_BD:
      updateParamFunc_ = &compositeParameter::getLinearBDBounds;
      break;
    case FunctionType::EXP:
      updateParamFunc_ = &compositeParameter::getExpBounds;
      break;
    case FunctionType::POLYNOMIAL:
      updateParamFunc_ = &compositeParameter::getPolynomialBounds;
      break;
    case FunctionType::LOGNORMAL:
      updateParamFunc_ = &compositeParameter::getLogNormalBounds;
      break;
    case FunctionType::REVERSE_SIGMOID:
      updateParamFunc_ = &compositeParameter::getReverseSigmoidBounds;
      break;
    default:
      throw Exception("compositeParameter::getParamUpdateFunction(): Invalid function type!");
      break;

  }
}
/******************************************************************************/
double compositeParameter::getRate(size_t state) const{
    switch(func_){
    case FunctionType::CONSTANT:
      return params_[0]->getValue();
    case FunctionType::LINEAR:
      return params_[0]->getValue() + ((double)(state-1)*params_[1]->getValue());
    case FunctionType::LINEAR_BD:
      return params_[0]->getValue() * (double)state;
    case FunctionType::EXP:
      return params_[0]->getValue() * std::exp((double)(state-1)*params_[1]->getValue());
    case FunctionType::POLYNOMIAL:
      throw Exception("Not implemented yet!");
      break;
    case FunctionType::LOGNORMAL:
      throw Exception("Not implemented yet!");
      break;
    case FunctionType::REVERSE_SIGMOID:
      throw Exception("Not implemented yet!");
      break;
    default:
      throw Exception("compositeParameter::createFunctionsMap(): Invalid function type!");
      break;

  }

} 
/******************************************************************************/
void compositeParameter::getConstBounds(size_t index, double* lowerBound, double* upperBound){
  if (index > 0){
    throw Exception ("compositeParameter::getConstBounds(): Too many parameters!");
  }
  *lowerBound = lowerBoundOfRateParam;
  *upperBound = upperBoundOfRateParam;
  
}
/******************************************************************************/
void compositeParameter::getLinearBounds(size_t index, double* lowerBound, double* upperBound){
  if (index > 1){
    throw Exception("compositeParameter::getLinearBounds(): Too many parameters!!!");
  }
  if (index == 0){
    *lowerBound = std::max(lowerBoundOfRateParam, -params_[index + 1]->getValue() * (maxChrNumber_-1));
    *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    *lowerBound = -params_[index - 1]->getValue()/(maxChrNumber_-1);
    *upperBound = upperBoundLinearRateParam;

  }else{
    throw Exception("compositeParameter::updateLinear(): To many parameters!!");

  }
  
}
/******************************************************************************/
void compositeParameter::getExpBounds(size_t index, double* lowerBound, double* upperBound){
  if (index > 1){
    throw Exception("compositeParameter::getExpBounds(): Too many parameters!!!");
  }
  if (index == 0){
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;

  }else{
    *lowerBound = lowerBoundOfExpParam;
    *upperBound = upperBoundExpParam/(maxChrNumber_-1);

  }


}
/******************************************************************************/
void compositeParameter::getLinearBDBounds(size_t index, double* lowerBound, double* upperBound){
  if (index > 0){
    throw Exception("compositeParameter::getLinearBDBounds(): Too many parameters!!!");
  }
  *lowerBound = lowerBoundOfRateParam;
  *upperBound = upperBoundOfRateParam;
}
/******************************************************************************/
size_t compositeParameter::getNumOfParameters(FunctionType func){
  if (func == FunctionType::CONSTANT || func == FunctionType::LINEAR_BD){
    return 1;
  }else if (func == FunctionType::EXP || func == FunctionType::LINEAR){
    return 2;
  }else if (func == FunctionType::FUNC_COUNT){
    return 0;
  }
  return 3;
}
/******************************************************************************/
compositeParameter::ParamName compositeParameter::getCompositeRateType(int param){
  switch (param)
  {
  case ChromosomeSubstitutionModel::GAIN:
    return compositeParameter::GAIN;
  case ChromosomeSubstitutionModel::LOSS:
    return compositeParameter::LOSS;
  case ChromosomeSubstitutionModel::DUPL:
    return compositeParameter::DUPL;
  case ChromosomeSubstitutionModel::DEMIDUPL:
    return compositeParameter::DEMI_DUPL;
  case ChromosomeSubstitutionModel::BASENUMR:
    return compositeParameter::BASENUMR;

  default:
    throw Exception("compositeParameter::getCompositeRateType(): No such parameter exists!!!");
    return compositeParameter::PARAMNAME_COUNT;
  }
}
/******************************************************************************/
void compositeParameter::updateBounds(ParameterList& params, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, FunctionType funcType, int maxChrNum){
  // if param name is in the following form: gain1 for example
  if ((funcType == FunctionType::CONSTANT || funcType == FunctionType::LINEAR_BD) || (funcType == FunctionType::EXP)){
    std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());
    *lowerBound = interval->getLowerBound();
    *upperBound = interval->getUpperBound();

    return;
  }

  std::vector<double> paramsValues;
  for (size_t i = 0; i < paramsNames.size(); i++){
    paramsValues.push_back(params.getParameter(paramsNames[i]).getValue());

  }
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>(params.getParameter(paramsNames[index]).getConstraint());

  if (funcType == FunctionType::LINEAR){
    updateBoundsLinear(paramsValues, index, lowerBound, upperBound, maxChrNum);
    interval->setLowerBound(*lowerBound, interval->strictLowerBound());
  }else if (funcType == FunctionType::LOGNORMAL){
    updateBoundsLogNormal(paramsValues, index, lowerBound, upperBound, maxChrNum);
  }else if (funcType == FunctionType::POLYNOMIAL){
    updateBoundsPolynomial(paramsValues, index, lowerBound, upperBound, maxChrNum);
  }else if (funcType == FunctionType::REVERSE_SIGMOID){
    updateBoundsReverseSigmoid(paramsValues, index, lowerBound, upperBound, maxChrNum);

  }else{
    throw Exception("compositeParameter::updateBounds(): no match for a given function!");
  }
  *lowerBound = interval->getLowerBound();
  *upperBound = interval->getUpperBound();

}
/******************************************************************************/
// void compositeParameter::updateBounds(Function* f, std::vector<string> paramsNames, size_t index, double* lowerBound, double* upperBound, FunctionType funcType, int maxChrNum){
//   // if param name is in the following form: gain1 for example
//   auto params = f->getParameters();
//   if ((funcType == FunctionType::CONSTANT || funcType == FunctionType::LINEAR_BD) || (funcType == FunctionType::EXP)){
//     std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>((&(f->getParameter(paramsNames[index])))->getConstraint());
//     *lowerBound = interval->getLowerBound();
//     *upperBound = interval->getUpperBound();

//     return;
//   }
  

//   std::vector<double> paramsValues;
//   for (size_t i = 0; i < paramsNames.size(); i++){
//     paramsValues.push_back(params.getParameter(paramsNames[i]).getValue());

//   }
//   std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>((&(f->getParameter(paramsNames[index])))->getConstraint());

//   if (funcType == FunctionType::LINEAR){
//     updateBoundsLinear(paramsValues, index, lowerBound, upperBound, maxChrNum);
//     interval->setLowerBound(*lowerBound, interval->strictLowerBound());
//   }else if (funcType == FunctionType::LOGNORMAL){
//     updateBoundsLogNormal(paramsValues, index, lowerBound, upperBound, maxChrNum);
//   }else if (funcType == FunctionType::POLYNOMIAL){
//     updateBoundsPolynomial(paramsValues, index, lowerBound, upperBound, maxChrNum);
//   }else if (funcType == FunctionType::REVERSE_SIGMOID){
//     updateBoundsReverseSigmoid(paramsValues, index, lowerBound, upperBound, maxChrNum);

//   }else{
//     throw Exception("compositeParameter::updateBounds(): no match for a given function!");
//   }
//   *lowerBound = interval->getLowerBound();
//   *upperBound = interval->getUpperBound();

// }
/******************************************************************************/
void compositeParameter::updateBounds(Function* f, const std::string &paramName, double &lowerBound, double &upperBound, FunctionType funcType){
  // if param name is in the following form: gain1 for example
  
  if ((funcType == FunctionType::CONSTANT || funcType == FunctionType::LINEAR_BD) || (funcType == FunctionType::EXP)){
    return;
  }
  
  std::shared_ptr<IntervalConstraint> interval = dynamic_pointer_cast<IntervalConstraint>((&(f->getParameter(paramName)))->getConstraint());

  if (funcType == FunctionType::LINEAR){
    interval->setLowerBound(lowerBound, interval->strictLowerBound());
  }else if (funcType == FunctionType::LOGNORMAL){
    throw Exception("compositeParameter::updateBounds(): Not implemented yet!");
  }else if (funcType == FunctionType::POLYNOMIAL){
    throw Exception("compositeParameter::updateBounds(): Not implemented yet!");
  }else if (funcType == FunctionType::REVERSE_SIGMOID){
    throw Exception("compositeParameter::updateBounds(): Not implemented yet!");

  }else{
    throw Exception("compositeParameter::updateBounds(): no match for a given function!");
  }

}
/******************************************************************************/
void compositeParameter::getBoundsForInitialParams(FunctionType funcType, size_t index, vector<double> paramsValues, double* lowerBound, double* upperBound, int maxChrNum, bool random){
  switch(funcType){
    case FunctionType::CONSTANT:
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;

      break;
    case FunctionType::LINEAR:
      updateBoundsLinear(paramsValues, index, lowerBound, upperBound, maxChrNum, random);
      break;
    case FunctionType::LINEAR_BD:
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
      break;
    case FunctionType::EXP:
      updateBoundsExp(paramsValues, index, lowerBound, upperBound, maxChrNum);
      break;
    case FunctionType::POLYNOMIAL:
      updateBoundsPolynomial(paramsValues, index, lowerBound, upperBound, maxChrNum);
      break;
    case FunctionType::LOGNORMAL:
      updateBoundsLogNormal(paramsValues, index, lowerBound, upperBound, maxChrNum);
      break;
    case FunctionType::REVERSE_SIGMOID:
      updateBoundsReverseSigmoid(paramsValues, index, lowerBound, upperBound, maxChrNum);
      break;
    default:
      throw Exception("compositeParameter::getBoundsForInitialParams(): Invalid function type!");
      break;

  }

}
/***************************************************************************************/
void compositeParameter::getAbsoluteBounds(FunctionType func, size_t index, double* lowerBound, double* upperBound, int maxChrNumber){
  switch(func){
    case FunctionType::CONSTANT:
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
      break;
    case FunctionType::LINEAR:
      if (index == 0){
        *lowerBound = -upperBoundLinearRateParam*(maxChrNumber - 1);
        *upperBound = upperBoundOfRateParam;
      }else{
        *lowerBound = -upperBoundOfRateParam/(maxChrNumber -1);
        *upperBound = upperBoundLinearRateParam;
      }
      
      break;
    case FunctionType::LINEAR_BD:
      *lowerBound = lowerBoundOfRateParam;
      *upperBound = upperBoundOfRateParam;
      break;
    case FunctionType::EXP:
      if (index == 0){
        *lowerBound = lowerBoundOfRateParam;
        *upperBound = upperBoundOfRateParam;

      }else{
        *lowerBound = lowerBoundOfExpParam;
        *upperBound = upperBoundExpParam/(maxChrNumber-1);

      }
      break;
    case FunctionType::POLYNOMIAL:
      throw Exception("compositeParameter::getAbsoluteBounds(): Not implemented yet!!");
      break;
    case FunctionType::LOGNORMAL:
      throw Exception("compositeParameter::getAbsoluteBounds(): Not implemented yet!!");
      break;
    case FunctionType::REVERSE_SIGMOID:
      throw Exception("compositeParameter::getAbsoluteBounds(): Not implemented yet!!");
      break;
    default:
      throw Exception("compositeParameter::getAbsoluteBounds(): Invalid function type!");
      break;

  }

}
/***************************************************************************************/
void compositeParameter::updateBoundsExp(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum){
  if (index > 1){
    throw Exception("compositeParameter::getExpBounds(): Too many parameters!!!");
  }
  if (index == 0){
    *lowerBound = lowerBoundOfRateParam;
    *upperBound = upperBoundOfRateParam;

  }else if (index == 1){
    *lowerBound = lowerBoundOfExpParam;
    *upperBound = upperBoundExpParam/(maxChrNum-1);

  }


}
/******************************************************************************/
void compositeParameter::updateBoundsLinear(std::vector<double> paramsValues, size_t index, double* lowerBound, double* upperBound, int &maxChrNum, bool random){
  if (index == 0){
    if (random){
      *lowerBound = lowerBoundOfRateParam;
    }else{
      *lowerBound = std::max(lowerBoundOfRateParam, -paramsValues[1]*(maxChrNum-1));
    }
    *upperBound = upperBoundOfRateParam;
    
  }else{
    *lowerBound = -paramsValues[0]/(maxChrNum-1);
    *upperBound = upperBoundLinearRateParam;
        
  }

}

/****************************************************************************/
std::vector<double> compositeParameter::getParameterValues() const{
  std::vector<double> values;
  for (size_t i = 0; i < size_; i++){
    values.push_back(params_[i]->getValue());
  }
  return values;
 
}
/*****************************************************************************/
std::vector<std::string> compositeParameter::getRelatedParameterNames(ParameterList &params, std::string pattern){
  std::vector<std::string> paramNames = params.getParameterNames();
  std::vector<std::string> matchingNames;
  for (size_t i = 0; i < paramNames.size(); i++){
    std::string fullParamName = paramNames[i];
    if (fullParamName.find(pattern) != string::npos){
      matchingNames.push_back(fullParamName);
    }
  }
  return matchingNames;
}

/*****************************************************************************/
// Chromosome model ///////////////////////////////////////////////////////////
/******************************************************************************/
ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(
  const ChromosomeAlphabet* alpha, 
  std::vector<double> gain, 
  std::vector<double> loss, 
  std::vector<double> dupl, 
  std::vector<double> demi,
  int baseNum,
  std::vector<double> baseNumR,
  unsigned int chrRange, 
  rootFreqType freqType,
  std::vector<int> rateChangeType):
    AbstractParameterAliasable("Chromosome."),
    AbstractSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, alpha->getMin(), alpha->getMax(), false)), "Chromosome."),
    gain_(0),
    loss_(0),
    dupl_(0),
    demiploidy_(0),
    baseNum_(baseNum),
    baseNumR_(0),
    maxChrRange_(chrRange),
    freqType_(freqType),
    ChrMinNum_(alpha->getMin()),
    ChrMaxNum_(alpha->getMax()),
    firstNormQ_(0),
    pijtCalledFromDeriv_(false),
    gainFunc_(),
    lossFunc_(),
    duplFunc_(),
    demiFunc_(),
    baseNumRFunc_(),
    vPowExp_()
{
    for (size_t i = 0; i < rateChangeType.size(); i++){
      switch (i)
      {
      case compositeParameter::GAIN:
        gainFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
        break;
      case compositeParameter::LOSS:
        lossFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
        break;
      case compositeParameter::DUPL:
        duplFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
        break;
      case compositeParameter::DEMI_DUPL:
        demiFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
        break;
      case compositeParameter::BASENUMR:
        baseNumRFunc_ =  static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
        break;
      default:
        throw Exception("ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(): No such parameter!");
        break;
      }

    }

    updateParameters(gain, loss, dupl, demi, baseNumR);
    computeFrequencies(false);
    isScalable_ = false;    //in ChromEvol the matrix should be not normalized
    updateMatrices();

}
/******************************************************************************/
ChromosomeSubstitutionModel::ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
  std::map<int, vector<double>> mapOfParamValues,
  int baseNum,
  unsigned int chrRange, 
  rootFreqType freqType,
  vector<int> rateChangeType):
    AbstractParameterAliasable("Chromosome."),
    AbstractSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, alpha->getMin(), alpha->getMax(), false)), "Chromosome."),
    gain_(0),
    loss_(0),
    dupl_(0),
    demiploidy_(0),
    baseNum_(baseNum),
    baseNumR_(0),
    maxChrRange_(chrRange),
    freqType_(freqType),
    ChrMinNum_(alpha->getMin()),
    ChrMaxNum_(alpha->getMax()),
    firstNormQ_(0),
    pijtCalledFromDeriv_(false),
    gainFunc_(),
    lossFunc_(),
    duplFunc_(),
    demiFunc_(),
    baseNumRFunc_(),
    vPowExp_()
{
  for (size_t i = 0; i < rateChangeType.size(); i++){
    switch (i)
    {
    case compositeParameter::GAIN:
      gainFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
      break;
    case compositeParameter::LOSS:
      lossFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
      break;
    case compositeParameter::DUPL:
      duplFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
      break;
    case compositeParameter::DEMI_DUPL:
      demiFunc_ = static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
      break;
    case compositeParameter::BASENUMR:
      baseNumRFunc_ =  static_cast<compositeParameter::FunctionType>(rateChangeType[i]);
      break;
    default:
      throw Exception("ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(): No such parameter!");
      break;
    }

  }

  updateParameters(mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::GAIN)], 
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::LOSS)], 
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DUPL)],
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::DEMIDUPL)],
    mapOfParamValues[static_cast<int>(ChromosomeSubstitutionModel::BASENUMR)]);
  computeFrequencies(false);
  isScalable_ = false;    //in ChromEvol the matrix should be not normalized
  updateMatrices();

}
// ChromosomeSubstitutionModel::ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
//   vector<double> modelParams,
//   unsigned int chrRange,
//   rootFreqType freqType,
//   rateChangeFunc rateChangeType):
//     AbstractParameterAliasable("Chromosome."),
//     AbstractSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, alpha->getMin(), alpha->getMax(), false)), "Chromosome."),
//     gain_(0),
//     loss_(0),
//     dupl_(0),
//     demiploidy_(0),
//     gainR_(0),
//     lossR_(0),
//     baseNum_(0),
//     baseNumR_(0),
//     duplR_(0),
//     maxChrRange_(chrRange),
//     freqType_(freqType),
//     rateChangeFuncType_(rateChangeType),
//     ChrMinNum_(alpha->getMin()),
//     ChrMaxNum_(alpha->getMax()),
//     firstNormQ_(0),
//     pijtCalledFromDeriv_(false),
//     vPowExp_()
//   {
//     //initialize model parameters
//     for (size_t i = 0; i < modelParams.size(); i++){
//       switch(i){
//         case BASENUM:
//           baseNum_ = (int)modelParams[i];
//           break;
//         case BASENUMR:
//           baseNumR_ = modelParams[i];
//           break;
//         case DUPL:
//           dupl_ = modelParams[i];
//           break;
//         case LOSS:
//           loss_ = modelParams[i];
//           break;
//         case GAIN:
//           gain_ = modelParams[i];
//           break;
//         case DEMIDUPL:
//           demiploidy_ = modelParams[i];
//           break;
//         case LOSSR:
//           lossR_ = modelParams[i];
//           break;
//         case GAINR:
//           gainR_ = modelParams[i];
//           break;
//         case DUPLR:
//           duplR_ = modelParams[i];
//           break;
//         default:
//           throw Exception("ChromsomeSubstitutionModel::ChromsomeSubstitutionModel(): Invalid rate type!");
//           break;

//       }
//     }

//     updateParameters();
//     computeFrequencies(false);
//     isScalable_ = false;    //in ChromEvol the matrix should be not normalized
//     updateMatrices();

//   }
/******************************************************************************/
ChromosomeSubstitutionModel* ChromosomeSubstitutionModel::initRandomModel(
  const ChromosomeAlphabet* alpha,
  int &baseNumber,
  map<int, vector<double>> initParams,
  unsigned int chrRange,
  rootFreqType rootFrequenciesType,
  vector<int> rateChangeType,
  vector<int>& fixedParams,
  double parsimonyBound)
{
  std::map<int, vector<double>> mapRandomParams;
  int newBaseNumber = baseNumber;

  //double gain, loss, dupl, demiDupl, gainR, lossR, baseNumR, duplR;
  //int baseNum;

  //vector<double> randomParams;
  //randomParams.reserve(NUM_OF_CHR_PARAMS);
  //map<int, double> setOfFixedParameters;
  //getSetOfFixedParameters(initParams, fixedParams, setOfFixedParameters);
  for (int i = 0; i < ChromosomeSubstitutionModel::paramType::NUM_OF_CHR_PARAMS; i++){
    if (std::find(fixedParams.begin(), fixedParams.end(), i) != fixedParams.end()){
      if (static_cast<ChromosomeSubstitutionModel::paramType>(i) != ChromosomeSubstitutionModel::BASENUM){
        mapRandomParams[i] = initParams[i];

      }
    }else{
      vector<double> paramValues;
      double lowerBound;
      double upperBound;
      if (static_cast<ChromosomeSubstitutionModel::paramType>(i) == ChromosomeSubstitutionModel::BASENUM){
        lowerBound = lowerBoundBaseNumber;
        upperBound = std::max((int)chrRange, lowerBoundBaseNumber+1);
        newBaseNumber = static_cast<int>(lowerBound + RandomTools::giveIntRandomNumberBetweenZeroAndEntry(upperBound-lowerBound));
        continue;

      }else{
        int compositeParamType = compositeParameter::getCompositeRateType(i);
        if (initParams[i].size() != 0){
          compositeParameter::FunctionType func = static_cast<compositeParameter::FunctionType>(rateChangeType[compositeParamType]);
          int numOfParameters = static_cast<int>(compositeParameter::getNumOfParameters(func));
          for (size_t j = 0; j < (size_t)numOfParameters; j++){
            compositeParameter::getBoundsForInitialParams(func, j, paramValues, &lowerBound, &upperBound, alpha->getMax(), true);
            if (parsimonyBound > 0){
              upperBound = std::min(upperBound, parsimonyBound);
            }
            double randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBound, upperBound);
            paramValues.push_back(randomValue);

          }
          
        }
      }
      mapRandomParams[i] = paramValues;

    }
  }
  // double upperBound = upperBoundOfRateParam;
  // double upperBoundLinear = upperBoundLinearRateParam;
  // double upperBoundExp = upperBoundExpParam/(alpha->getMax()-1);
  // if (parsimonyBound > 0){
  //   upperBound = std::min(upperBoundOfRateParam, parsimonyBound);
  //   upperBoundLinear = std::min(upperBoundLinearRateParam, parsimonyBound);
  //   upperBoundExp = std::min(upperBoundExpParam/(alpha->getMax()-1), parsimonyBound);
  // }
  // for (size_t i = 0; i < initParams.size(); i ++){
  //   getRandomParameter(static_cast<paramType>(i), initParams[i], randomParams, upperBound, upperBoundLinear, upperBoundExp, rateChangeType, alpha->getMax(), chrRange, setOfFixedParameters);
  // }

  ChromosomeSubstitutionModel* model = new ChromosomeSubstitutionModel(alpha, mapRandomParams, newBaseNumber, chrRange, rootFrequenciesType, rateChangeType);//, useExtendedFloat);
  return model;

}
/******************************************************************************/
// void ChromosomeSubstitutionModel::getRandomParameter(paramType type, double initParamValue, vector<double>& randomParams, double upperBound, double upperBoundLinear, double upperBoundExp, rateChangeFunc rateFunc, int maxChrNum, unsigned int chrRange, map<int, double>& setOfFixedParameters){
//   // there is an assumption that the rate change parameters are sampled after the const ones, since they are dependent on them
//   double randomValue = initParamValue;
//   if (type == BASENUM){
//     int lowerBoundBaseNum = lowerBoundBaseNumber;
//     int upperBoundBaseNum = std::max((int)chrRange, lowerBoundBaseNumber+1);
//     if ((initParamValue != IgnoreParam) && (setOfFixedParameters.count(BASENUM) == 0)){
//       randomValue = lowerBoundBaseNum + RandomTools::giveIntRandomNumberBetweenZeroAndEntry(upperBoundBaseNum-lowerBoundBaseNum);
//     }

//   }
//   else if (type == BASENUMR){ //|| (type == DUPL)) || ((type == GAIN) || (type == LOSS))){
//     if ((initParamValue != IgnoreParam) && (setOfFixedParameters.count(BASENUMR) == 0)){
//       randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
//     }
//   }
//   else if (type == DEMIDUPL){
//     if ((initParamValue != IgnoreParam) && (initParamValue != DemiEqualDupl)){
//       if (setOfFixedParameters.count(type) == 0){
//         randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
//       }      
//     }
//   }else{
//     paramType typeOfPairedRate;
//     switch (type)
//     {
//       case GAIN:
//         typeOfPairedRate = GAINR;
//         break;
//       case LOSS:
//         typeOfPairedRate = LOSSR;
//         break;
//       case DUPL:
//         typeOfPairedRate = DUPLR;
//         break;
//       case GAINR:
//         typeOfPairedRate = GAIN;
//         break;
//       case LOSSR:
//         typeOfPairedRate = LOSS;
//         break;
//       case DUPLR:
//         typeOfPairedRate = DUPL;
//         break;

//       default:
//         throw Exception("ChromosomeSubstitutionModel::getRandomParameter(): Invalid rate type!");
//         break;
//     }
//     if (((type == GAIN) || (type == LOSS)) || (type == DUPL)){
//       double lowerBound = lowerBoundOfRateParam;
//       if ((initParamValue != IgnoreParam) && (setOfFixedParameters.count(type) == 0)){
//         if ((setOfFixedParameters.count(typeOfPairedRate) > 0) && (setOfFixedParameters[typeOfPairedRate] != IgnoreParam)){
//           if ( rateFunc == LINEAR){
//             lowerBound = std::max(lowerBoundOfRateParam, -setOfFixedParameters[typeOfPairedRate]*(maxChrNum-1));
//           }        
//         }
//         randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBound, upperBound);
//       }
//     }else{
//       if ((initParamValue != IgnoreParam) && (setOfFixedParameters.count(type) == 0)){
//         if (randomParams[typeOfPairedRate] != IgnoreParam){
//           if (rateFunc == LINEAR){
//             double lowerBoundForRate = -randomParams[typeOfPairedRate]/(maxChrNum-1);
//             randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundForRate, upperBoundLinear);
//           }else{
//             randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfExpParam, upperBoundExp);      
//           }
//         }else{
//           randomValue = RandomTools::giveRandomNumberBetweenTwoPoints(lowerBoundOfRateParam, upperBound);
//         }
//       }

//     }


//   }
//   randomParams.push_back(randomValue);

// }
/******************************************************************************/
// void ChromosomeSubstitutionModel::getSetOfFixedParameters(vector<double>& initParams, vector<unsigned int>& fixedParams, map<int, double>& setOfFixedParams){
//   //map<paramType, double> setOfFixedParams;
//   size_t index = 0;
//   for (size_t i = 0; i < NUM_OF_CHR_PARAMS; i++){
//     if (initParams[i] == IgnoreParam){
//       continue;
//     }
//     if (fixedParams[index]){
//       setOfFixedParams[(int)i] = initParams[i];
//     }
//     index++;
        
//   }
//   return;

// }
/******************************************************************************/
std::vector<Parameter*> ChromosomeSubstitutionModel::createCompositeParameter(compositeParameter::FunctionType &func, std::string paramName, vector<double> &vectorOfValues){
  std::vector<Parameter*> params;
  if (func == compositeParameter::FunctionType::FUNC_COUNT){
    return params;
  }
  for (size_t i = 0; i < vectorOfValues.size(); i++){
    auto paramValue = vectorOfValues[i];
    if (paramValue == IgnoreParam){
      throw Exception("ChromosomeSubstitutionModel::createCompositeParameter(): Function is not supposed to be defined as a legal function name if the parameter should be ignored!");
    }
    double lowerBound;
    double upperBound;
    compositeParameter::getAbsoluteBounds(func, i, &lowerBound, &upperBound, ChrMaxNum_);

    //compositeParameter::getBoundsForInitialParams(func, i, vectorOfValues, &lowerBound, &upperBound, ChrMaxNum_);
    std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBound, upperBound, false, true);
    Parameter* param = new Parameter("Chromosome."+ paramName + std::to_string(i), paramValue, interval);
    params.push_back(param);
    //addParameter_(param);
  }
  return params;


}

/******************************************************************************/
void ChromosomeSubstitutionModel::updateParameters(vector<double> &gain, vector<double> &loss, vector<double> &dupl, vector<double> &demi, vector<double> &baseNumR){
  if (baseNum_ != IgnoreParam){
    std::shared_ptr<IntervalConstraint> interval_baseNum = make_shared<IntervalConstraint>(lowerBoundBaseNumber, (int)maxChrRange_, true, true);
    addParameter_(new Parameter("Chromosome.baseNum", baseNum_, interval_baseNum));
    
  }
  std::vector<size_t> numOfParamsVector;
  auto baseNumParams = createCompositeParameter(baseNumRFunc_, "baseNumR", baseNumR);
  numOfParamsVector.push_back(baseNumParams.size());
  auto duplParams = createCompositeParameter(duplFunc_, "dupl", dupl);
  numOfParamsVector.push_back(duplParams.size());
  auto lossParams = createCompositeParameter(lossFunc_, "loss", loss);
  numOfParamsVector.push_back(lossParams.size());
  auto gainParams = createCompositeParameter(gainFunc_, "gain", gain);
  numOfParamsVector.push_back(gainParams.size());
  auto demiParams = createCompositeParameter(demiFunc_, "demi", demi);
  numOfParamsVector.push_back(demiParams.size());
  size_t maxSize = *max_element(numOfParamsVector.begin(), numOfParamsVector.end());
  for (size_t i = 0; i < maxSize; i++){
    if (i < baseNumParams.size()){
      addParameter_(baseNumParams[i]);
    }
    if (i < duplParams.size()){
      addParameter_(duplParams[i]);
    }
    if (i < lossParams.size()){
      addParameter_(lossParams[i]);
    }
    if (i < gainParams.size()){
      addParameter_(gainParams[i]);
    }
    if (i < demiParams.size()){
      addParameter_(demiParams[i]);
    }
  }
  if (gainFunc_ != compositeParameter::FunctionType::FUNC_COUNT){
    gain_ = new compositeParameter(ChrMaxNum_, gainFunc_, "gain", gainParams);

  }
  if (lossFunc_ != compositeParameter::FunctionType::FUNC_COUNT){
    loss_ = new compositeParameter(ChrMaxNum_, lossFunc_, "loss", lossParams);
  }
  if (duplFunc_ != compositeParameter::FunctionType::FUNC_COUNT){
    dupl_ = new compositeParameter(ChrMaxNum_, duplFunc_, "dupl", duplParams);
  }
  if (demiFunc_ != compositeParameter::FunctionType::FUNC_COUNT){
    demiploidy_ = new compositeParameter(ChrMaxNum_, demiFunc_, "demi", demiParams);
  }
  if (baseNumRFunc_ != compositeParameter::FunctionType::FUNC_COUNT){
    baseNumR_ = new compositeParameter(ChrMaxNum_, baseNumRFunc_, "baseNumR", baseNumParams);
  }

    // std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBoundOfRateParam, upperBoundOfRateParam, false, true);
    // if ((baseNum_ != IgnoreParam) && (baseNumR_ != IgnoreParam)){
    //   updateBaseNumParameters(interval); 
    // }
    // updateConstRateParameter(dupl_, duplR_, "Chromosome.dupl", interval);
    // updateConstRateParameter(loss_, lossR_, "Chromosome.loss", interval);
    // updateConstRateParameter(gain_, gainR_, "Chromosome.gain", interval);
 
    // if (rateChangeFuncType_ == rateChangeFunc::LINEAR){
    //   updateLinearParameters();
    // }else if (rateChangeFuncType_ == rateChangeFunc::EXP){
    //   updateExpParameters();
    // }
    // if ((demiploidy_ != IgnoreParam) & (demiploidy_!= DemiEqualDupl)){
    //   addParameter_(new Parameter("Chromosome.demi", demiploidy_, interval));
          
    // }

}
/******************************************************************************/
// void ChromosomeSubstitutionModel::updateConstRateParameter(double paramValueConst, double paramValueChange, string parameterName, std::shared_ptr<IntervalConstraint> interval)
// {
//   if (paramValueConst != IgnoreParam){
//     if (paramValueChange == IgnoreParam){
//       addParameter_(new Parameter(parameterName, paramValueConst, interval));
//     }else{
//       double lowerBound = lowerBoundOfRateParam;
//       if (rateChangeFuncType_ == rateChangeFunc::LINEAR){
//         lowerBound = std::max(lowerBoundOfRateParam, -paramValueChange*(getMax()-1));       
//       }
//       std::shared_ptr<IntervalConstraint> intervalForCompositeRate = make_shared<IntervalConstraint>(lowerBound, upperBoundOfRateParam, false, true);
//       addParameter_(new Parameter(parameterName, paramValueConst, intervalForCompositeRate));
//     }
      
//   }

// }
/******************************************************************************/
// void ChromosomeSubstitutionModel::updateExpParameters(){
//   std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBoundOfExpParam, upperBoundExpParam/(getMax()-1), false, true);
//   if (lossR_ != IgnoreParam){
//     if (loss_ != IgnoreParam){
//       addParameter_(new Parameter("Chromosome.lossR", lossR_, interval));
//     }else{
//       std::shared_ptr<IntervalConstraint> intervalNoConstRateLoss = make_shared<IntervalConstraint>(lowerBoundOfRateParam, upperBoundExpParam/(getMax()-1), false, true);
//       addParameter_(new Parameter("Chromosome.lossR", lossR_, intervalNoConstRateLoss));
//     }

//   }
//   if (gainR_ != IgnoreParam){
//     if (gain_ != IgnoreParam){
//       addParameter_(new Parameter("Chromosome.gainR", gainR_, interval));
//     }else{
//       std::shared_ptr<IntervalConstraint> intervalNoConstRateGain = make_shared<IntervalConstraint>(lowerBoundOfRateParam, upperBoundExpParam/(getMax()-1), false, true);
//       addParameter_(new Parameter("Chromosome.gainR", gainR_, intervalNoConstRateGain));
//     }
    
//   }
//   if (duplR_ != IgnoreParam){
//     if (dupl_ != IgnoreParam){
//       addParameter_(new Parameter("Chromosome.duplR", duplR_, interval));
      
//     }else{
//       std::shared_ptr<IntervalConstraint> intervalNoConstRateDupl = make_shared<IntervalConstraint>(lowerBoundOfRateParam, upperBoundExpParam/(getMax()-1), false, true);
//       addParameter_(new Parameter("Chromosome.duplR", duplR_, intervalNoConstRateDupl));

//     }
    
//   }
// }
/******************************************************************************/
// void ChromosomeSubstitutionModel::updateLinearChangeParameter(double constParam, double linearParam, string paramName){
//   double lowerBound, upperBound;
//   if (linearParam != IgnoreParam){
//     if (constParam != IgnoreParam){
//       lowerBound = -(constParam/(getMax()-1));
//       upperBound = upperBoundLinearRateParam;
//     }else{
//       lowerBound = lowerBoundOfRateParam;
//       upperBound = upperBoundOfRateParam;
//     }
//     std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBound, upperBound, false, true);
//     addParameter_(new Parameter(paramName, linearParam, interval));

//   }


// }
/******************************************************************************/
// void ChromosomeSubstitutionModel::updateLinearParameters(){
//   updateLinearChangeParameter(loss_, lossR_, "Chromosome.lossR");
//   updateLinearChangeParameter(gain_, gainR_, "Chromosome.gainR");
//   updateLinearChangeParameter(dupl_, duplR_, "Chromosome.duplR");

// }
/******************************************************************************/
// void ChromosomeSubstitutionModel::updateBaseNumParameters(std::shared_ptr<IntervalConstraint> interval){
//   std::shared_ptr<IntervalConstraint> interval_baseNum = make_shared<IntervalConstraint>(lowerBoundBaseNumber, (int)maxChrRange_, true, true);
//   addParameter_(new Parameter("Chromosome.baseNum", baseNum_, interval_baseNum));
//   addParameter_(new Parameter("Chromosome.baseNumR", baseNumR_, interval));

// }
/******************************************************************************/
void ChromosomeSubstitutionModel::getCompositeParametersValues(std::string paramName, compositeParameter* param){
  for (size_t i = 0; i < param->getSize(); i++){   
    getParameterValue(paramName + std::to_string(i)); //do I really need it?
    //auto intervals = dynamic_pointer_cast<IntervalConstraint>(param->getParams()[i]->getConstraint());
    //double lowerBound = intervals->getLowerBound();
    //double upperBound = intervals->getUpperBound();
    //param->getBounds(i, &lowerBound, &upperBound);
    //intervals->setLowerBound(lowerBound, false);
    //intervals->setUpperBound(upperBound, true);
  }
}
/******************************************************************************/
void ChromosomeSubstitutionModel::getParametersValues(){
    if (gain_ != 0){
      getCompositeParametersValues("gain", gain_);
      //gain_ = getParameterValue("gain");
    }
    if (loss_ != 0){
      getCompositeParametersValues("loss", loss_);
      //loss_ = getParameterValue("loss");
    }
    if (dupl_ != 0){
      getCompositeParametersValues("dupl", dupl_);
      //dupl_ = getParameterValue("dupl");
    }
    if (demiploidy_ != 0){// && (demiploidy_ != DemiEqualDupl)){
      getCompositeParametersValues("demi", demiploidy_);
      //demiploidy_ = getParameterValue("demi");
    }
    // update intervals
    //setNewBoundsForLinearParameters(gain_, gainR_, "gain", "gainR");
    //setNewBoundsForLinearParameters(loss_, lossR_, "loss", "lossR");
    //setNewBoundsForLinearParameters(dupl_, duplR_, "dupl", "duplR");

    if(baseNum_ != IgnoreParam){
      baseNum_ = (int)getParameterValue("baseNum");
    }
    if (baseNumR_ != 0){
      getCompositeParametersValues("baseNumR", baseNumR_);
    }  
    //checkParametersBounds();
}

/******************************************************************************/
// void ChromosomeSubstitutionModel::setNewBoundsForLinearParameters(double &constRate, double &changeRate, string paramNameConst, string paramNameLinear){

//   if (changeRate != IgnoreParam){
//     changeRate = getParameterValue(paramNameLinear);
//     if ((rateChangeFuncType_ == rateChangeFunc::LINEAR) && (constRate != IgnoreParam)){

//       std::dynamic_pointer_cast<IntervalConstraint>(getParameter(paramNameLinear).getConstraint())->setLowerBound(-constRate/(getMax()-1), true);

//       std::dynamic_pointer_cast<IntervalConstraint>(getParameter(paramNameConst).getConstraint())->setLowerBound(std::max(lowerBoundOfRateParam, -changeRate * (getMax()-1)), true);
        
//     }
//   }

// }

/*******************************************************************************/
void ChromosomeSubstitutionModel::checkParametersBounds() const{
  std::cout << "All bounds" <<endl;
  const ParameterList params = getParameters();
  for (size_t i = 0; i < params.size(); i++){
    std::cout << params[i].getName() << " bound: "<< dynamic_pointer_cast<IntervalConstraint>(params[i].getConstraint())->getLowerBound() <<  ", value: " << params[i].getValue() << std::endl;

  }
}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateMatrices(){
    //update model parameters
    getParametersValues();

    //update generator matrix
    size_t maxChrNum = (size_t)(getMax());
    size_t minChrNum = (size_t)(getMin()); 
    MatrixTools::fill(generator_, 0);
 
    // updating Q matrix
    for (size_t i = minChrNum; i < maxChrNum+1; i++){
        // gain
        if (i + 1 < maxChrNum+1){
            updateQWithGain(i, minChrNum);
        //loss
        }if (i-1 >= minChrNum){
            updateQWithLoss(i, minChrNum);
        //duplication         
        }if (2*i <= maxChrNum){
            updateQWithDupl(i, minChrNum);
        }else if (i != maxChrNum){
            updateQWithDupl(i, minChrNum, maxChrNum);
        }
        //demi-ploidy
        updateQWithDemiDupl(i, minChrNum, maxChrNum);   

        if (i < maxChrNum){
          if (baseNum_ != IgnoreParam){
            updateQWithBaseNumParameters(i, minChrNum, maxChrNum);

          }
        }
        
    }
    setDiagonal();  //sets Qii to -sigma(Qij)
    updateEigenMatrices();
    firstNormQ_ = getFirstNorm();

}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithDemiDupl(size_t i, size_t minChrNum, size_t maxChrNum){
  //double demiploidy;
  
  if (demiploidy_ != 0){
    if (demiploidy_->getRate(i) < 0){
      throw Exception("ChromosomeSubstitutionModel::updateQWithDemiDupl(): Negative demiploidy rate!");
    }
    if (i % 2 == 0 && (double)i * 1.5 <= (double)maxChrNum){

      generator_(i-minChrNum, (size_t)((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i);

                        
    }else if (i % 2 != 0 && (size_t)ceil((double)i*1.5) <= maxChrNum){
      if (i == 1){
        generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i);
      }else{
        generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i)/2;
        generator_(i-minChrNum, (size_t)floor((double)i * 1.5)-minChrNum) += demiploidy_->getRate(i)/2;

      }

    }else{
      if (i != maxChrNum){
        generator_(i-minChrNum, maxChrNum-minChrNum) += demiploidy_->getRate(i);
      }

    }

  }
}
/*******************************************************************************/
// double ChromosomeSubstitutionModel::getRate (size_t state, double constRate, double changeRate) const{
//   if ((constRate == IgnoreParam) && (changeRate == IgnoreParam)){
//     return IgnoreParam;
//   }
//   double totalRate;
//   if (constRate == IgnoreParam){
//     // a birth-death-like model
//     totalRate = changeRate;
//   }else{
//     //const rate is not to be ignored
//     totalRate = constRate;
//   }
//   if (changeRate == IgnoreParam){
//     return totalRate; //only const rate
//   }else{
//     if (rateChangeFuncType_ == rateChangeFunc::LINEAR){
//       totalRate += (changeRate* (double)(state-1));
//     }else if (rateChangeFuncType_ == rateChangeFunc::EXP){
//       totalRate *= (exp(changeRate* (double)(state-1)));
//     }
//   }
//   return totalRate;
// }
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithGain(size_t i, size_t minChrNum){
  if (gain_ == 0){
    return;
  }
  double gainRate = gain_->getRate(i);
  if (gainRate < 0){
    throw Exception ("ChromosomeSubstitutionModel::updateQWithGain(): negative gain rate!");
  }
  generator_(i-minChrNum, i+1-minChrNum) += gainRate;


}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithLoss(size_t i, size_t minChrNum){
  //generator_(i-minChrNum, i-1-minChrNum) = loss_ + (lossR_* i);
  if (loss_ == 0){
    return;
  }
  double lossRate = loss_->getRate(i);
  if (lossRate < 0){
    throw Exception ("ChromosomeSubstitutionModel::updateQWithLoss(): negative loss rate!");
  }
  generator_(i-minChrNum, i-1-minChrNum) += lossRate;


}
/*******************************************************************************/
void ChromosomeSubstitutionModel::updateQWithDupl(size_t i, size_t minChrNum, size_t maxChrNum){
  if (dupl_ == 0){
    return;
  }
  // if the transition is not to maxChr
  double duplRate = dupl_->getRate(i);
  if (duplRate < 0){
    throw Exception("ChromosomeSubstitutionModel::updateQWithDupl(): Negative dupl rate!");
  }

  if (maxChrNum == 0){
    generator_(i-minChrNum, (2 * i)-minChrNum) += duplRate;

  }else{
     generator_(i-minChrNum, maxChrNum-minChrNum) += duplRate;

  }
}


/********************************************************************************/
void ChromosomeSubstitutionModel::updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum){
  if ( baseNumR_->getRate(currChrNum) < 0){
    throw Exception("ChromosomeSubstitutionModel::updateQWithBaseNumParameters():Negative base number rate!");
  }
  for (size_t j = currChrNum + 1; j < maxChrNum + 1; j ++){
    if (j == maxChrNum){
      if ((j-currChrNum) <= maxChrRange_){
        generator_(currChrNum-minChrNum, maxChrNum-minChrNum) += baseNumR_->getRate(currChrNum);
      }
    }else{
      if ((j-currChrNum) % baseNum_ == 0){
        if ((j-currChrNum) <= maxChrRange_){
          generator_(currChrNum - minChrNum, j - minChrNum) += baseNumR_->getRate(currChrNum);
        }
      }
    }
  }
}

/******************************************************************************/
void ChromosomeSubstitutionModel::setFreq(map<int, double>& freqs)
{
  for (size_t i = 0; i < size_; ++i)
  {
    freq_[i] = freqs[static_cast<int>(i)];
  }

}
/******************************************************************************/
void ChromosomeSubstitutionModel::updateEigenMatrices()
{
  // Compute eigen values and vectors:
  if (enableEigenDecomposition())
  {
    // Look for null lines (such as stop lines)
    // ie null diagonal elements

    size_t nbStop=0;
    size_t salph = getNumberOfStates();
    vector<bool> vnull(salph); // vector of the indices of lines with
                               // only zeros

    for (size_t i = 0; i < salph; i++)
    {
      if (abs(generator_(i, i)) < NumConstants::TINY())
      {
        nbStop++;
        vnull[i]=true;
      }
      else
        vnull[i]=false;
    }
        
    if (nbStop != 0)
    {
      size_t salphok=salph - nbStop;
      
      RowMatrix<double> gk(salphok, salphok);
      size_t gi = 0, gj = 0;

      for (size_t i = 0; i < salph; i++)
      {
        if (!vnull[i])
        {
          gj = 0;
          for (size_t j = 0; j < salph; j++)
          {
            if (!vnull[j])
            {
              gk(i - gi, j - gj) = generator_(i, j);
            }
            else
              gj++;
          }
        }
        else
          gi++;
      }

      EigenValue<double> ev(gk);
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();

      for (size_t i = 0; i < nbStop; i++)
      {
        eigenValues_.push_back(0);
        iEigenValues_.push_back(0);
      }

      RowMatrix<double> rev = ev.getV();
      rightEigenVectors_.resize(salph, salph);
      gi = 0;
      for (size_t i = 0; i < salph; i++)
      {
        if (vnull[i])
        {
          gi++;
          for (size_t j = 0; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }

          rightEigenVectors_(i, salphok + gi - 1) = 1;
        }
        else
        {
          for (size_t j = 0; j < salphok; j++)
          {
            rightEigenVectors_(i, j) = rev(i - gi, j);
          }

          for (size_t j = salphok; j < salph; j++)
          {
            rightEigenVectors_(i, j) = 0;
          }
        }
      }
    }
    else
    {
      EigenValue<double> ev(generator_);
      rightEigenVectors_ = ev.getV();
      eigenValues_ = ev.getRealEigenValues();
      iEigenValues_ = ev.getImagEigenValues();
      nbStop = 0;
    }

    /// Now check inversion and diagonalization
    try
    {
      MatrixTools::inv(rightEigenVectors_, leftEigenVectors_);

      // is it diagonalizable ?
      isDiagonalizable_ = true;

      if (!dynamic_cast<ReversibleSubstitutionModel*>(this))
      {
        for (auto& vi : iEigenValues_)
        {
          if (abs(vi) > NumConstants::TINY())
          {
            isDiagonalizable_ = false;
            break;
          }
        }
      }
      
      // looking for the vector of 0 eigenvalues

      vector<size_t> vNullEv;
      for (size_t i = 0; i< salph - nbStop; i++)
        if ((abs(eigenValues_[i]) < NumConstants::SMALL()) && (abs(iEigenValues_[i]) < NumConstants::SMALL()))
          vNullEv.push_back(i);
      

      // pb to find unique null eigenvalue      
      isNonSingular_=(vNullEv.size()==1);

      size_t nulleigen;
      
      double val;
      if (!isNonSingular_)
      {
        //look or check which non-stop right eigen vector elements are
        //equal.
        for (auto cnull : vNullEv)
        {
          size_t i = 0;
          while (vnull[i])
            i++;
          
          val = rightEigenVectors_(i, cnull);
          i++;
          
          while (i < salph)
          {
            if (!vnull[i])
            {
              if (abs(rightEigenVectors_(i, cnull) - val) > NumConstants::SMALL())
                break;
            }
            i++;
          }
          
          if (i >= salph)
          {
            isNonSingular_ = true;
            nulleigen=cnull;
            break;
          }
        }
      }
      else
        nulleigen=vNullEv[0];
      
      if (isNonSingular_)
      {
        eigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches
        iEigenValues_[nulleigen] = 0; // to avoid approximation errors on long long branches


      }
      else
      {
        //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Unable to find eigenvector for eigenvalue 0. Taylor series used instead.");
        isDiagonalizable_ = false;
      }
    }
    // if rightEigenVectors_ is singular
    catch (ZeroDivisionException& e)
    {
      //ApplicationTools::displayMessage("AbstractSubstitutionModel::updateMatrices : Singularity during diagonalization. Taylor series used instead.");
      isNonSingular_ = false;
      isDiagonalizable_ = false;
    }

    if (vPowExp_.size() == 0){
      vPowExp_.resize(30);
    }
    MatrixTools::getId(salph, vPowExp_[0]);
    MatrixTools::Taylor(generator_, 30, vPowExp_);

  }

}
/******************************************************************************/
void ChromosomeSubstitutionModel::calculatePijtUsingEigenValues(double t) const{
  if (isDiagonalizable_){
    MatrixTools::mult<double>(rightEigenVectors_, VectorTools::exp(eigenValues_ * (rate_ * t)), leftEigenVectors_, pijt_);
  }else{
    std::vector<double> vdia(size_);
    std::vector<double> vup(size_ - 1);
    std::vector<double> vlo(size_ - 1);
    double c = 0, s = 0;
    double l = rate_ * t;
    for (size_t i = 0; i < size_; i++){
      vdia[i] = std::exp(eigenValues_[i] * l);
      if (iEigenValues_[i] != 0){
        s = std::sin(iEigenValues_[i] * l);
        c = std::cos(iEigenValues_[i] * l);
        vup[i] = vdia[i] * s;
        vlo[i] = -vup[i];
        vdia[i] *= c;
        vdia[i + 1] = vdia[i]; // trick to avoid computation
        i++;
      }else{
        if (i < size_ - 1){
          vup[i] = 0;
          vlo[i] = 0;
        }
      }
    }
    MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, pijt_);
  }
}


/******************************************************************************/

const Matrix<double>& ChromosomeSubstitutionModel::getPij_t(double t) const
{
  if (t == 0)
  {
    MatrixTools::getId(size_, pijt_);
  }
  else if (isNonSingular_)
  {
    calculatePijtUsingEigenValues(t);
    
  }
  else
  {
    RowMatrix<double> pijt_temp;
    MatrixTools::getId(size_, pijt_temp);
    double s = 1.0;
    double v = rate_ * t;
    double norm = v * firstNormQ_;
    size_t m = 0;
    bool converged = false;
    //while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
    while (norm > 0.5)
    {
      m += 1;
      v /= 2;
      norm /= 2;
    }
    for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
      calculateExp_Qt(iternum, &s, v);

      if (iternum > 2){
        converged = checkIfReachedConvergence(pijt_, pijt_temp);
        if (converged){
          break;
        }
      }
      MatrixTools::copy(pijt_, pijt_temp);
      if (iternum > 250){
        //std :: cout << "ERROR: Pijt did not reach convergence for t = "<< t <<"!"<<endl;
        throw Exception("ChromosomeSubstitutionModel: Taylor series did not reach convergence!");
        break;
      }
      if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
        RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
        MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
        vPowExp_.push_back(new_pow);

      }

    }
    while (m > 0){  // recover the 2^m
      
      MatrixTools::mult(pijt_, pijt_, tmpMat_);
      MatrixTools::copy(tmpMat_, pijt_);

      m--;
    }
  }
  
  //just for test/////////////////////
  // bool correct = true;
  if (!pijtCalledFromDeriv_){
    for (size_t i = 0; i < size_; i++){
      for (size_t j = 0; j < size_; j++){
        if (pijt_(i,j) < 0){
          pijt_(i,j) = NumConstants::VERY_TINY(); // trying to do it exactly as in ChromEvol. Maybe the "nan" problem will be solved
          //pijt_(i,j) = 0;
        }
        else if (pijt_(i, j) > 1){
          pijt_(i,j) = 1.0;
        }

      }
    }

  }

  /////////////////////////////
  return pijt_;
}
double ChromosomeSubstitutionModel::getFirstNorm() const{
  double norm = 0;
  for (size_t i = 0; i < size_; i++){
    for (size_t j = 0; j < size_; j++){
      norm += fabs(generator_(i,j));

    }
  }
  return norm;
}

/******************************************************************************/

bool ChromosomeSubstitutionModel::checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const{
    for (size_t i = 0; i < pijt.getNumberOfRows(); i++){
        for (size_t j = 0; j < pijt.getNumberOfColumns(); j++){
            double diff = fabs(pijt(i,j) - mt_prev(i,j));
            if (diff > get_epsilon()){
                return false;
            }else if ((pijt(i,j) + get_epsilon() < 0) || (pijt(i,j) > 1 + get_epsilon())){
              return false;
            }
        }
    }
    return true;
}

/******************************************************************************/
void ChromosomeSubstitutionModel::calculateExp_Qt(size_t pow, double s, size_t m, double v) const{
  MatrixTools::getId(size_, pijt_);
  for (size_t i = 1; i <= pow; i++){
    s *= v / static_cast<double>(i);// the initial value of v is rt/(2^m)
    MatrixTools::add(pijt_, s, vPowExp_[i]);

  }
  while (m > 0)  // recover the 2^m
  {
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  }

}

/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getdPij_dt  (double d) const{
  pijtCalledFromDeriv_ = true;
  RowMatrix<double> pijt;
  //if (!pijt_calculated_){
  pijt = getPij_t(d);
  MatrixTools::mult(pijt, generator_, dpijt_);
  MatrixTools::scale(dpijt_, rate_);
  pijtCalledFromDeriv_ = false;
  return dpijt_;

}

/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getd2Pij_dt2(double d) const{
  pijtCalledFromDeriv_ = true;
  RowMatrix<double> pijt;
  //if (!pijt_calculated_){
  pijt = getPij_t(d);
  
  MatrixTools::mult(vPowExp_[2], pijt, d2pijt_);
  MatrixTools::scale(d2pijt_, rate_ * rate_);
  pijtCalledFromDeriv_ = false;
  return d2pijt_;
}

/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getPij_t_func2(double d) const{
  RowMatrix<double> pijt_temp;
  MatrixTools::getId(size_, pijt_temp);
  double s = 1.0;
  double v = rate_ * d;
  size_t m = 0;
  bool converged = false;
  while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  {
    m += 1;
    v /= 2;
  }
  for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
    calculateExp_Qt(iternum, s, m, v);

    if (iternum > 1){
      converged = checkIfReachedConvergence(pijt_, pijt_temp);
      if (converged){
        break;
      }
    }
    MatrixTools::copy(pijt_, pijt_temp);
    if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
      RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
      MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
      vPowExp_.push_back(new_pow);

    }

  }
  return pijt_;

}
/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getPij_t_func3(double d) const{
  MatrixTools::getId(size_, pijt_);
  double s = 1.0;
  double v = rate_ * d;
  size_t m = 0;
  while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  {
    m += 1;
    v /= 2;
  }
  for (size_t i = 1; i < vPowExp_.size(); i++)
  {
    s *= v / static_cast<double>(i);
    MatrixTools::add(pijt_, s, vPowExp_[i]);
  }
  while (m > 0)  // recover the 2^m
  {
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);
    m--;
  }

//  MatrixTools::print(pijt_);
  return pijt_;

}

/******************************************************************************/
void ChromosomeSubstitutionModel::calculateExp_Qt(size_t pow, double* s, double v) const{
  if (pow == 2){
    MatrixTools::getId(size_, pijt_);
    for (size_t i = 1; i <= pow; i++){
      *s *= v / static_cast<double>(i);// the initial value of v is rt/(2^m)
      MatrixTools::add(pijt_, *s, vPowExp_[i]);
      
    }

  }else{
    *s *= v / static_cast<double>(pow);
    MatrixTools::add(pijt_, *s, vPowExp_[pow]);
 
  }
  

}
/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getPij_t_func4(double d) const{
  RowMatrix<double> pijt_temp;
  MatrixTools::getId(size_, pijt_temp);
  double s = 1.0;
  double v = rate_ * d;
  size_t m = 0;
  bool converged = false;
  while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  {
    m += 1;
    v /= 2;
  }
  for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
    calculateExp_Qt(iternum, &s, v);

    if (iternum > 2){
      converged = checkIfReachedConvergence(pijt_, pijt_temp);
      if (converged){
        break;
      }
    }
    MatrixTools::copy(pijt_, pijt_temp);
    if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
      RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
      MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
      vPowExp_.push_back(new_pow);

    }

  }
  while (m > 0){  // recover the 2^m
      
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  }
  return pijt_;

}
/*********************************************************************************/
double ChromosomeSubstitutionModel::getInitValue(size_t i, int state) const
{
  if (i >= size_)
    throw IndexOutOfBoundsException("ChromosomeSubstitutionModel::getInitValue", i, 0, size_ - 1);
  if (state < 0 || !alphabet_->isIntInAlphabet(state))
    throw BadIntException(state, "ChromosomeSubstitutionModel::getInitValue. Character " + alphabet_->intToChar(state) + " is not allowed in model.");
  vector<int> states = alphabet_->getAlias(state);
  for (size_t j = 0; j < states.size(); j++)
  {
     if (getAlphabetStateAsInt(i) == states[j]){
       if (dynamic_cast<const ChromosomeAlphabet*>(alphabet_)){
         const ChromosomeAlphabet* alpha = dynamic_cast<const ChromosomeAlphabet*>(alphabet_);
         // it is a composite state
         if (state > alpha->getMax() + 1){
           return alpha->getProbabilityForState(state, states[j]);

         }else{
           return 1.0;
         }

       }else{
         return 1.;
       }

     }
  }
  return 0.;
}

const Matrix<double>& ChromosomeSubstitutionModel::getPijt_test(double t) const {
  RowMatrix<double> pijt_temp;
  MatrixTools::getId(size_, pijt_temp);
  double s = 1.0;
  double v = rate_ * t;
  double norm = v * firstNormQ_;
  size_t m = 0;
  bool converged = false;
  //while (v > 0.5)    // exp(r*t*A)=(exp(r*t/(2^m) A))^(2^m)
  while (norm > 0.5)
  {
    m += 1;
    v /= 2;
    norm /= 2;
  }
  for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
    calculateExp_Qt(iternum, &s, v);

    if (iternum > 2){
      converged = checkIfReachedConvergence(pijt_, pijt_temp);
      if (converged){
        break;
      }
    }
    MatrixTools::copy(pijt_, pijt_temp);
    if (iternum > 250){
      //std :: cout << "ERROR: Pijt did not reach convergence for t = "<< t <<"!"<<endl;
      throw Exception("ChromosomeSubstitutionModel: Taylor series did not reach convergence!");
      break;
    }
    if (iternum == vPowExp_.size()-1 && !converged){  //need to add more powers to the matrix
      RowMatrix<double> new_pow;
      //new_pow.resize(size_, size_);
      MatrixTools :: mult(vPowExp_[vPowExp_.size()-1], generator_, new_pow);
      vPowExp_.push_back(new_pow);

    }

  }
  while (m > 0){  // recover the 2^m
      
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  }
  //just for test/////////////////////
  // bool correct = true;
  if (!pijtCalledFromDeriv_){
    for (size_t i = 0; i < size_; i++){
      for (size_t j = 0; j < size_; j++){
        if (pijt_(i,j) < 0){
          pijt_(i,j) = NumConstants::VERY_TINY(); // trying to do it exactly as in ChromEvol. Maybe the "nan" problem will be solved
          //pijt_(i,j) = 0;
        }
        else if (pijt_(i, j) > 1){
          pijt_(i,j) = 1.0;
        }

      }
    }

  }
  return pijt_;

}

/******************************************************************************/
