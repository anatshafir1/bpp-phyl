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
#define upperBoundOfRateParam 100.0
#define IgnoreParam -999
#define DemiEqualDupl -2
using namespace std;
namespace bpp
{


class ChromosomeSubstitutionModel :
  public AbstractSubstitutionModel
{
public:
  enum rootFreqType {UNIFORM, ROOT_LL, STATIONARY, FIXED};

private:
  double gain_;
  double loss_;
  double dupl_;
  double demiploidy_;
  rootFreqType freqType_;
  int ChrMinNum_;
  int ChrMaxNum_;
  

protected:
  mutable std::vector< RowMatrix<double> > vPowExp_;
  std::vector <int> parameterIncluded_;
  //const ChromosomeAlphabet* chromosomeAlpha_;
  //mutable double lambda_, exp_;
  //mutable RowMatrix<double> p_;
  //bool pijt_calculated_;


public:
  ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, double gain, double loss, double dupl, double demi, rootFreqType freqType);

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
  //size_t getMaxChrNum(const Alphabet* alpha);
  //size_t getMinChrNum(const Alphabet* alpha);    
  //const ChromosomeAlphabet* getChromosomeAlphabet() const { return chromosomeAlpha_; }

protected:
  void updateParameters();
  void updateMatrices();
  void updateEigenMatrices();
  void getParametersValues();
  void calculateExp_Qt(size_t pow, double s, size_t m, double v) const;
  void calculateExp_Qt(size_t pow, double* s, double v) const;
  double get_epsilon() const{ return 0.0001;};
};
} // end of namespace bpp.

#endif  // _CHROMOSOMESUBSTITUTIONMODEL_H_