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
ChromosomeSubstitutionModel :: ChromosomeSubstitutionModel(const ChromosomeAlphabet* alpha, double gain, double loss, double dupl, double demi, rootFreqType freqType):
    AbstractParameterAliasable("Chromosome."),
    AbstractSubstitutionModel(alpha, std::shared_ptr<const StateMap>(new CanonicalStateMap(alpha, alpha->getMin(), alpha->getMax(), false)), "Chromosome."),
    gain_(gain),
    loss_(loss),
    dupl_(dupl),
    demiploidy_(demi),
    freqType_(freqType),
    ChrMinNum_(alpha->getMin()),
    ChrMaxNum_(alpha->getMax()),
    vPowExp_(),
    parameterIncluded_()
{

    updateParameters();
    computeFrequencies(false);
    isScalable_ = false;    //in ChromEvol the matrix should be not normalized
    updateMatrices();

}

/******************************************************************************/
void ChromosomeSubstitutionModel::updateParameters(){
    std::shared_ptr<IntervalConstraint> interval = make_shared<IntervalConstraint>(lowerBoundOfRateParam, upperBoundOfRateParam, true, true);
    if (dupl_ != IgnoreParam){
      addParameter_(new Parameter("Chromosome.dupl", dupl_, interval));
      parameterIncluded_.push_back(1);
    }else{
      dupl_ = 0;
      parameterIncluded_.push_back(0);
    }
    if (loss_ != IgnoreParam){
      addParameter_(new Parameter("Chromosome.loss", loss_, interval));
      parameterIncluded_.push_back(1);
    }else{
      loss_ = 0;
      parameterIncluded_.push_back(0);
    }
    if (gain_ != IgnoreParam){
      addParameter_(new Parameter("Chromosome.gain", gain_, interval));
      parameterIncluded_.push_back(1);
    }else{
      gain_ = 0;
      parameterIncluded_.push_back(0);
    }
    if ((demiploidy_ != IgnoreParam) & (demiploidy_!= DemiEqualDupl)){
      addParameter_(new Parameter("Chromosome.demi", demiploidy_, interval));
      parameterIncluded_.push_back(1);
    }else{
      if (demiploidy_ == IgnoreParam){
        demiploidy_ = 0;
        parameterIncluded_.push_back(0);
      }else{
        demiploidy_ = dupl_;
        parameterIncluded_.push_back(2);
      }
      
    }   

}
/******************************************************************************/
void ChromosomeSubstitutionModel::getParametersValues(){
    if (parameterIncluded_[0]){
      gain_ = getParameterValue("gain");
    }
    if (parameterIncluded_[1]){
      loss_ = getParameterValue("loss");
    }
    if (parameterIncluded_[2]){
      dupl_ = getParameterValue("dupl");
    }
    if (parameterIncluded_[3] == 1){
      demiploidy_ = getParameterValue("demi");
    }else if (parameterIncluded_[3] == 2){
      demiploidy_ = dupl_;
    }   
}

/******************************************************************************/
void ChromosomeSubstitutionModel::updateMatrices(){
    //update model parameters
    getParametersValues();

    //update generator matrix
    size_t maxChrNum = (size_t)(getMax());
    size_t minChrNum = (size_t)(getMin()); 
    MatrixTools::fill(generator_, 0);
 
    // updating Q matrix
    for (size_t i = minChrNum; i < maxChrNum+1; i++){
        //if (i + 1 <= max_chr_num+1){
        if (i + 1 < maxChrNum+1){
            generator_(i-minChrNum, i+1-minChrNum) += gain_;

        }if (i-1 >= 1){
            generator_(i-minChrNum, i-1-minChrNum) = loss_;
          
        }if (2*i <= maxChrNum){

            generator_(i-minChrNum, (2 * i)-minChrNum) += dupl_;
        }else if (i != maxChrNum){
            generator_(i-minChrNum, maxChrNum-minChrNum) += dupl_;

        }
            
        if (i % 2 == 0 && (double)i * 1.5 <= (double)maxChrNum){

            generator_(i-minChrNum, (size_t)((double)i * 1.5)-minChrNum) += demiploidy_;

                        
        }else if (i % 2 != 0 && (size_t)ceil((double)i*1.5) <= maxChrNum){
            if (i == 1){
              generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_;
            }else{
              generator_(i-minChrNum, (size_t)ceil((double)i * 1.5)-minChrNum) += demiploidy_/2;
              generator_(i-minChrNum, (size_t)floor((double)i * 1.5)-minChrNum) += demiploidy_/2;

            }

        }else{
          if (i != maxChrNum){
            generator_(i-minChrNum, maxChrNum-minChrNum) += demiploidy_;
          }

        }
        
    }
    setDiagonal();  //sets Qii to -sigma(Qij)
    updateEigenMatrices();

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

    //if (!isNonSingular_)
    //{
/*       double min = generator_(0, 0);
      for (size_t i = 1; i < salph; i++)
      {
        if (min > generator_(i, i))
          min = generator_(i, i);
      } */

      //setScale(-1 / min);

    if (vPowExp_.size() == 0)
      vPowExp_.resize(30);

      

    MatrixTools::getId(salph, vPowExp_[0]);
    //}

    // normalization
    //normalize();
    
    //if (!isNonSingular_)
    MatrixTools::Taylor(generator_, 30, vPowExp_);
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
    if (isDiagonalizable_)
    {
      MatrixTools::mult<double>(rightEigenVectors_, VectorTools::exp(eigenValues_ * (rate_ * t)), leftEigenVectors_, pijt_);
    }
    else
    {
      std::vector<double> vdia(size_);
      std::vector<double> vup(size_ - 1);
      std::vector<double> vlo(size_ - 1);
      double c = 0, s = 0;
      double l = rate_ * t;
      for (size_t i = 0; i < size_; i++)
      {
        vdia[i] = std::exp(eigenValues_[i] * l);
        if (iEigenValues_[i] != 0)
        {
          s = std::sin(iEigenValues_[i] * l);
          c = std::cos(iEigenValues_[i] * l);
          vup[i] = vdia[i] * s;
          vlo[i] = -vup[i];
          vdia[i] *= c;
          vdia[i + 1] = vdia[i]; // trick to avoid computation
          i++;
        }
        else
        {
          if (i < size_ - 1)
          {
            vup[i] = 0;
            vlo[i] = 0;
          }
        }
      }
      MatrixTools::mult<double>(rightEigenVectors_, vdia, vup, vlo, leftEigenVectors_, pijt_);
    }
  }
  else
  {
    RowMatrix<double> pijt_temp;
    MatrixTools::getId(size_, pijt_temp);
    double s = 1.0;
    double v = rate_ * t;
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
      if (iternum > 250){
        std :: cout << "ERROR: Pijt did not reach convergence for t = "<< t <<"!"<<endl;
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
  return pijt_;
}

  /*     for (size_t iternum = 2; iternum <  vPowExp_.size(); iternum++){
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

    } */
//  MatrixTools::print(pijt_);
  //pijt_calculated_ = true;



/******************************************************************************/

bool ChromosomeSubstitutionModel::checkIfReachedConvergence(const Matrix<double>& pijt, const Matrix<double>& mt_prev) const{
    for (size_t i = 0; i < pijt.getNumberOfRows(); i++){
        for (size_t j = 0; j < pijt.getNumberOfColumns(); j++){
            double diff = fabs(pijt(i,j) - mt_prev(i,j));
            if (diff > get_epsilon()){
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
  RowMatrix<double> pijt;
  //if (!pijt_calculated_){
  pijt = getPij_t(d);
  MatrixTools::mult(pijt, generator_, dpijt_);
  MatrixTools::scale(dpijt_, rate_);
    //pijt_calculated_ = true;
  //}else{
    //mult(pijt_, generator_, dpijt_);
  //}
  
  // dp(t) = p(t)*rate_*Q
  return dpijt_;

}


/******************************************************************************/
const Matrix<double>& ChromosomeSubstitutionModel::getd2Pij_dt2(double d) const{
  RowMatrix<double> pijt;
  //if (!pijt_calculated_){
  pijt = getPij_t(d);
  
  MatrixTools::mult(vPowExp_[2], pijt, d2pijt_);
  MatrixTools::scale(d2pijt_, rate_ * rate_);
    //pijt_calculated_ = true;
  //}else{
    //mult(vPowGen_[2], pijt_, d2pijt_);
  //}
  // ddp(t) = Q^2 * p(t)*rate_^2
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
  

/*   while (m > 0)  // recover the 2^m
  {
    MatrixTools::mult(pijt_, pijt_, tmpMat_);
    MatrixTools::copy(tmpMat_, pijt_);

    m--;
  } */

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

/* size_t ChromosomeSubstitutionModel::getMaxChrNum(const Alphabet* alpha){
  //Alphabet* new_alpha = alpha->clone();
  //AbstractAlphabet* chr_alpha = dynamic_cast <AbstractAlphabet*>(new_alpha);
  
  size_t number_of_alpha_states = alpha->getNumberOfStates();
  AlphabetState state = alpha->getStateAt(number_of_alpha_states-1);
  size_t max_state_num = state.getNum();
  return max_state_num;

} */
/******************************************************************************/
/* 
size_t ChromosomeSubstitutionModel::getMinChrNum(const Alphabet* alpha){
  //Alphabet* new_alpha = alpha->clone();
  //AbstractAlphabet* chr_alpha = dynamic_cast <AbstractAlphabet*>(new_alpha);
  AlphabetState state = alpha->getStateAt(1);
  size_t min_state_num = state.getNum();
  return min_state_num;

} */
