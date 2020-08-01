//
// File: CromosomeSubstitutionModel.h
// Created by: Anat Shafir
// Created on: 2020
//


#ifndef _CHROMOSOMESUBSTITUTIONMODEL_H_
#define _CHROMOSOMESUBSTITUTIONMODEL_H_

#include "AbstractSubstitutionModel.h"
#include <Bpp/Seq/Alphabet/ChromosomeAlphabet.h>

#define lowerBoundOfRateParam 0.0
#define lowerBoundOfExpParam -100.0
#define lowerBoundBaseNumber 3
#define upperBoundOfRateParam 100.0
#define upperBoundLinearRateParam 5.0
#define IgnoreParam -999
#define DemiEqualDupl -2
#define EPSILON 2.22045e-016
using namespace std;
namespace bpp
{


class ChromosomeSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED};
  enum rateChangeFunc {LINEAR = 0, EXP = 1};

private:
  double gain_;
  double loss_;
  double dupl_;
  double demiploidy_;
  double gainR_;
  double lossR_;
  int baseNum_;
  double baseNumR_;
  unsigned int maxChrRange_;
  rootFreqType freqType_;
  rateChangeFunc rateChangeFuncType_;
  bool optimizeBaseNum_;
  int ChrMinNum_;
  int ChrMaxNum_;
  double firstNormQ_;
  mutable bool pijtCalledFromDeriv_;

  

protected:
  mutable std::vector< RowMatrix<double> > vPowExp_;



public:
  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, 
    double gain, 
    double loss, 
    double dupl, 
    double demi,
    double gainR,
    double lossR,
    int baseNum,
    double baseNumR,
    unsigned int maxChrRange, 
    rootFreqType freqType,
    rateChangeFunc rateChangeType,
    bool optimizeBaseNumber);

  virtual ~ChromosomeSubstitutionModel() {}

  ChromosomeSubstitutionModel* clone() const { return new ChromosomeSubstitutionModel(*this); }

  
public:
  const Matrix<double>& getPij_t    (double d) const;
  const Matrix<double>& getPij_t_func2(double d) const;
  const Matrix<double>& getPij_t_func3(double d) const;
  const Matrix<double>& getPij_t_func4(double d) const;
  const Matrix<double>& getdPij_dt  (double d) const;
  const Matrix<double>& getd2Pij_dt2(double d) const;

  //void calculateExp_Qt(int pow, double s, size_t m, double v);
  std::string getName() const { return "Chromosome"; }
  void setFreq(std::map<int, double>& freqs);
  size_t getNumberOfStates() const { return size_; }
  int getMin() const {return ChrMinNum_;}
  int getMax() const {return ChrMaxNum_;}
  bool checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const;
  double getInitValue(size_t i, int state) const;
  //size_t getMaxChrNum(const Alphabet* alpha);
  //size_t getMinChrNum(const Alphabet* alpha);    
  //const ChromosomeAlphabet* getChromosomeAlphabet() const { return chromosomeAlpha_; }

protected:
  void updateParameters();
  void updateLinearParameters();
  void updateExpParameters();
  void updateBaseNumParameters(std::shared_ptr<IntervalConstraint> interval);
  void updateMatrices();
  void updateQWithBaseNumParameters(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateQWithGain(size_t currChrNum, size_t minChrNum);
  void updateQWithLoss(size_t currChrNum, size_t minChrNum);
  void updateQWithDemiDupl(size_t currChrNum, size_t minChrNum, size_t maxChrNum);
  void updateEigenMatrices();
  void getParametersValues();
  void calculateExp_Qt(size_t pow, double s, size_t m, double v) const;
  void calculateExp_Qt(size_t pow, double* s, double v) const;
  double getFirstNorm() const;
  double get_epsilon() const{ return 0.0001;};
};
} // end of namespace bpp.

#endif  // _CHROMOSOMESUBSTITUTIONMODEL_H_