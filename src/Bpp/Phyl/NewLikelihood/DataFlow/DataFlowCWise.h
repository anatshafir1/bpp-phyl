//
// File: DataFlowCwise.h
// Authors: Francois Gindraud (2017), Laurent Guéguen (2019)
// Created: 2018-06-07
// Last modified: 2018-07-11
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef DATAFLOW_CWISE_H
#define DATAFLOW_CWISE_H

#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloat.h>
#include <Eigen/Core>
#include <algorithm>
#include <cassert>
#include <string>
#include <tuple>
#include <list>
#include <type_traits>
#include <iostream>

#include "DataFlowNumeric.h"
#include <Bpp/Numeric/Parameter.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Parameter.h>
#include "Definitions.h"

namespace bpp {

  // Return a reference to the object for component-wise operations
  namespace numeric {
    template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    T & cwise (T & t) {
      return t; // Do nothing for basic types
    }

    inline ExtendedFloat & cwise (ExtendedFloat & t) {
      return t; // Do nothing 
    }

    inline const ExtendedFloat & cwise (const ExtendedFloat & t) {
      return t; // Do nothing 
    }

    template <typename Derived> auto cwise (const Eigen::MatrixBase<Derived> & m) -> decltype (m.array ()) {
      return m.array (); // Use Array API in Eigen
    }
    
    template <typename Derived> auto cwise (Eigen::MatrixBase<Derived> & m) -> decltype (m.array ()) {
      return m.array (); // Use Array API in Eigen
    }

    inline auto cwise (const Eigen::RowVectorXi & m)  -> decltype (m.template cast<double>().array ())
    {
      return m.template cast<double>().array ();
    }

    template<int R, int C>
    inline ExtendedFloatArray<R,C> cwise (const ExtendedFloatMatrix<R,C>& m)
    {
      return ExtendedFloatArray<R,C>(m.float_part().array(), m.exponent_part());
    }


    template<int R, int C>
    inline ExtendedFloatArrayWrapper<R,C> cwise (ExtendedFloatMatrix<R,C>& m)
    {
      return ExtendedFloatArrayWrapper<R,C>(m);
    }

  }
  
  /******************************************************************************
   * Data flow nodes for those numerical functions.
   *
   */

  template <typename Result, typename From> class CWiseFill;
  template <typename Result, typename From> class CWiseMatching;
  template <typename Result, typename From> class CWiseCompound;

  /*************************************************************************
   * @brief build a Value to a Matrix or rowVector filled with
   * values of references.  
   *
   * Node construction should be done with the create static method.
   */
    
  template <typename R, typename T> class CWiseFill : public Value<R> {
  public:
    using Self = CWiseFill;

    /// Build a new CWiseFill node.
    static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
      // Check dependencies
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 1);
      checkNthDependencyIsValue<T> (typeid (Self), deps, 0);

      return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
    }

    CWiseFill (NodeRefVec && deps, const Dimension<R> & dim)
      : Value<R> (std::move (deps)), targetDimension_ (dim)
    {
      this->accessValueMutable().resize(dim.rows,dim.cols);
    }

    std::string debugInfo () const override {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // CWiseFill additional arguments = ().
    bool compareAdditionalArguments (const Node_DF & other) const final {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef derive (Context & c, const Node_DF & node) final {
      if (&node == this) {
        return ConstantOne<R>::create (c, targetDimension_);
      }
      return Self::create (c, {this->dependency(0)->derive (c, node)}, targetDimension_);
    }

    NodeRef recreate (Context & c, NodeRefVec && deps) final {
      return Self::create (c, std::move (deps), targetDimension_);
    }

  private:
    void compute() override { compute<T>();}
      
    template<class U=T>
    typename std::enable_if<std::is_arithmetic<U>::value, void>::type
    compute () {
      using namespace numeric;
      auto & result = this->accessValueMutable ();
      const auto & x0 = accessValueConstCast<T> (*this->dependency (0));
      result=convert(x0,targetDimension_);
    }      

    template<class U=T>
    typename std::enable_if<std::is_same<U,RowLik>::value>::type
    compute ()
    {
      using namespace numeric;
      auto & result = this->accessValueMutable ();
      const auto & x0 = accessValueConstCast<T> (*this->dependency (0));
      result.colwise()=x0.transpose();
    }      

    template<class U=T>
    typename std::enable_if<std::is_same<U,VectorLik>::value>::type
    compute ()
    {
      using namespace numeric;
      auto & result = this->accessValueMutable ();
      const auto & x0 = accessValueConstCast<T> (*this->dependency (0));
      result.colwise() =x0;
    }      

    Dimension<R> targetDimension_;

  };

  
  
  /*************************************************************************
   * @brief build a Value to a Eigen T which columns are accessible
   * through a pattern of positions.
   *
   * Node construction should be done with the create static method.
   */

  typedef Eigen::Matrix<size_t,-1,1> PatternType;

  template <typename R> class CWisePattern : public Value<R> {

    class pattern_functor
    {
      const R& m_arg_;
      const PatternType& pattern_;
    public:
      pattern_functor(const R& arg, const PatternType& pattern) :
        m_arg_(arg), pattern_(pattern) {}

      
      const typename R::Scalar& operator()(Eigen::Index row, Eigen::Index col) const
      {
        return m_arg_(row, Eigen::Index(pattern_[col]));
      }

      const typename R::Scalar& operator()(Eigen::Index col) const
      {
        return m_arg_(Eigen::Index(pattern_[col]));
      }

      // Specific for ExtendedFloatEigen
      template<typename R2=R>
      ExtendedFloat::ExtType exponent_part(typename std::enable_if< std::is_base_of<ExtendedFloatEigenBase<R2>,R2>::value>::type* = 0) const
      {
        return m_arg_.exponent_part();
      }

    };
            
  public:
    using Self = CWisePattern;

    /// Build a new CWisePattern node.
    static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
      // Check dependencies
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyVectorSize (typeid (Self), deps, 2);
      checkNthDependencyIsValue<R> (typeid (Self), deps, 0);
      checkNthDependencyIsValue<PatternType> (typeid(Self), deps, 1);
      // Remove 0s from deps
      if (deps[1]->hasNumericalProperty (NumericalProperty::ConstantOne))
        return convertRef<Value<R>>(deps[0]);
      else
        return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
    }

    CWisePattern (NodeRefVec && deps, const Dimension<R> & dim)
      : Value<R> (deps), targetDimension_ (dim)
    {
    }

    std::string debugInfo () const override {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // CWisePattern additional arguments = ().
    bool compareAdditionalArguments (const Node_DF & other) const final {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef derive (Context & c, const Node_DF & node) final {
      if (&node == this) {
        return ConstantOne<R>::create (c, targetDimension_);
      }
      return Self::create (c, {this->dependency(0)->derive (c, node), this->dependency(1)}, targetDimension_);
    }

    NodeRef recreate (Context & c, NodeRefVec && deps) final {
      return Self::create (c, std::move (deps), targetDimension_);
    }

  private:
    void compute() override {
      const auto& arg=accessValueConstCast<R>(*this->dependency(0));
      const auto& pattern=accessValueConstCast<PatternType>(*this->dependency(1));
      this->accessValueMutable()=R::NullaryExpr(targetDimension_.rows, targetDimension_.cols, pattern_functor(arg, pattern));
    };
      
    Dimension<R> targetDimension_;

  };


  
  /*************************************************************************
   * @brief build a Value to a Matrix R which columns and rows are
   * accessible through a vector of T objects and a function of
   * matching positions from T objects to R object.
   *
   * This class is originally made for partitions of likelihoods.
   *
   * Node construction should be done with the create static method.
   *
   */


  /*
   * Matrix of matching positions : site X (index of T in the vector of Ts, position for corresponding T)
   *
   */
  
  typedef Eigen::Matrix<size_t,-1,2> MatchingType;

  template <typename R, typename T> class CWiseMatching<R, ReductionOf<T>> : public Value<R> {

    class matching_functor
    {
      const std::vector<const T*>& m_arg_;
      const MatchingType& matching_;
    public:
      matching_functor(const std::vector<const T*>& arg, const MatchingType& matching) :
        m_arg_(arg), matching_(matching) {}

      const typename R::Scalar& operator()(Eigen::Index row, Eigen::Index col) const
      {
        return compute<T>(row, col);
      }
      
      template<typename T2 = T>
      const typename R::Scalar& compute(Eigen::Index row, Eigen::Index col,
                                        typename std::enable_if< !std::is_same<T2,typename R::Scalar>::value, T*>::type* = 0) const
      {
        return (*m_arg_[matching_(col,0)])(row, Eigen::Index(matching_(col,1)));
      }


      // Specific case of Eigen::RowVector made from several elements
      template<typename T2 = T>
      const typename R::Scalar& compute(Eigen::Index row, Eigen::Index col,
                                        typename std::enable_if< std::is_same<T2,typename R::Scalar>::value, T*>::type* = 0) const
      {
        return *m_arg_[matching_(col,0)];
      }

      // Specific for ExtendedFloat
      template<typename R2=R>
      ExtendedFloat::ExtType exponent_part(typename std::enable_if< std::is_base_of<ExtendedFloatEigenBase<R2>,R2>::value>::type* = 0) const
      {
        std::vector<ExtendedFloat::ExtType> vexp(m_arg_.size());
        std::transform(m_arg_.begin(),m_arg_.end(),vexp.begin(),[](const T* t){return(t->exponent_part());});

        if (! std::equal(vexp.begin() + 1, vexp.end(), vexp.begin()) )
          throw Exception("DataFlowCwise::CWiseMatching not possible on ExtendedFloatEigen data with different exponents. Ask developpers.");
        
        return vexp[0];
      }
      
    };
            
  public:
    using Self = CWiseMatching;

    /// Build a new CWiseMatching node.
    static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
      // Check dependencies
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyRangeIsValue<T> (typeid (Self), deps, 0, deps.size () -1);
      checkNthDependencyIsValue<MatchingType> (typeid(Self), deps, deps.size()-1);
      
      // Remove 0s from deps

      return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
    }

    CWiseMatching (NodeRefVec && deps, const Dimension<R> & dim)
      : Value<R> (std::move(deps)), targetDimension_ (dim)
    {
    }

    std::string debugInfo () const override {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // CWisePattern additional arguments = ().
    bool compareAdditionalArguments (const Node_DF & other) const final {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef derive (Context & c, const Node_DF & node) final {
      if (&node == this) {
        return ConstantOne<R>::create (c, targetDimension_);
      }
      const auto n = this->nbDependencies ();
      NodeRefVec derivedDeps (n);
      for (std::size_t i = 0; i < n-1; ++i) {
        derivedDeps[i] = this->dependency (i)->derive (c, node);
      }
      derivedDeps[n-1] = this->dependency(n-1);
      
      return Self::create (c, std::move (derivedDeps), targetDimension_);
    }

    NodeRef recreate (Context & c, NodeRefVec && deps) final {
                  return Self::create (c, std::move (deps), targetDimension_);
    }

  private:
    void compute() override {
      const auto n = this->nbDependencies ();
      std::vector<const T*> vR(n-1);
      for (std::size_t i = 0; i < n-1; ++i) {
        vR[i] = &accessValueConstCast<T> (*this->dependency(i));
      }
      
      const auto& matching=accessValueConstCast<MatchingType>(*this->dependency(n-1));
      
      this->accessValueMutable()=R::NullaryExpr(targetDimension_.rows, targetDimension_.cols, matching_functor(vR, matching));
    };
      
    Dimension<R> targetDimension_;

  };
    
  /*************************************************************************
   * @brief build a Value to a Eigen R from a compound of lignes or columns.
   *
   * Node construction should be done with the create static method.
   *
   */

  template <typename R, typename T>  class CWiseCompound<R, ReductionOf<T>> : public Value<R> {

    class compound_functor
    {
      const std::vector<const T*>& m_arg_;
    public:
      compound_functor(const std::vector<const T*>& arg) :
        m_arg_(arg) {}

      const typename R::Scalar& operator()(Eigen::Index row, Eigen::Index col) const
      {
        return compute<T>(row, col);
      }
      
      template<typename T2 = T>
      const typename R::Scalar& compute(Eigen::Index row, Eigen::Index col,
                                        typename std::enable_if< std::is_same<T2,RowLik>::value, T*>::type* = 0) const
      {
        return (*m_arg_[size_t(row)])(col);
      }

      template<typename T2 = T>
      const typename R::Scalar& compute(Eigen::Index row, Eigen::Index col,
                                        typename std::enable_if< std::is_same<T2,VectorLik>::value, T*>::type* = 0) const
      {
        return (*m_arg_[size_t(col)])(row);
      }

      // Specific for ExtendedFloat
      template<typename R2=R>
      ExtendedFloat::ExtType exponent_part(typename std::enable_if< std::is_base_of<ExtendedFloatEigenBase<R2>,R2>::value>::type* = 0) const
      {
        std::vector<ExtendedFloat::ExtType> vexp(m_arg_.size());
        std::transform(m_arg_.begin(),m_arg_.end(),vexp.begin(),[](const T* t){return(t->exponent_part());});

        if (! std::equal(vexp.begin() + 1, vexp.end(), vexp.begin()) )
          throw Exception("DataFlowCwise::CWiseCompound not possible on ExtendedFloatEigen data with different exponents. Ask developpers.");
        
        return vexp[0];
      }

    };
            
  public:
    using Self = CWiseCompound;

    /// Build a new CWiseCompound node.
    static ValueRef<R> create (Context & c, NodeRefVec && deps, const Dimension<R> & dim) {
      // Check dependencies
      checkDependenciesNotNull (typeid (Self), deps);
      checkDependencyRangeIsValue<T> (typeid (Self), deps, 0, deps.size ());
      
      return cachedAs<Value<R>> (c, std::make_shared<Self> (std::move (deps), dim));
    }

    CWiseCompound (NodeRefVec && deps, const Dimension<R> & dim)
      : Value<R> (std::move(deps)), targetDimension_ (dim)
    {
    }

    std::string debugInfo () const override {
      using namespace numeric;
      return debug (this->accessValueConst ()) + " targetDim=" + to_string (targetDimension_);
    }

    // CWiseCompound additional arguments = ().
    bool compareAdditionalArguments (const Node_DF & other) const final {
      return dynamic_cast<const Self *> (&other) != nullptr;
    }

    NodeRef derive (Context & c, const Node_DF & node) final {
      if (&node == this) {
        return ConstantOne<R>::create (c, targetDimension_);
      }
      const auto n = this->nbDependencies ();
      NodeRefVec derivedDeps (n);
      for (std::size_t i = 0; i < n; ++i) {
        derivedDeps[i] = this->dependency (i)->derive (c, node);
      }

      return Self::create (c, std::move (derivedDeps), targetDimension_);
    }

    NodeRef recreate (Context & c, NodeRefVec && deps) final {
      return Self::create (c, std::move (deps), targetDimension_);
    }

  private:
    void compute() override {
      const auto n = this->nbDependencies ();
      std::vector<const T*> vR(n);
      for (std::size_t i = 0; i < n; ++i) {
        vR[i] = &accessValueConstCast<T> (*this->dependency(i));
      }
      
      this->accessValueMutable()=R::NullaryExpr(targetDimension_.rows, targetDimension_.cols, compound_functor(vR));
    };
      
    Dimension<R> targetDimension_;

  };
    
  // Precompiled instantiations
  extern template class CWiseFill<RowLik, double>;
  extern template class CWiseFill<VectorLik, double>;
  extern template class CWiseFill<MatrixLik, VectorLik>;
  extern template class CWiseFill<MatrixLik, RowLik>;

  extern template class CWisePattern<RowLik>;
  extern template class CWisePattern<MatrixLik>;
  
  extern template class CWiseMatching<RowLik, ReductionOf<RowLik>>;
  extern template class CWiseMatching<MatrixLik, ReductionOf<MatrixLik>>;
  extern template class CWiseMatching<MatrixLik, ReductionOf<RowLik>>;
  extern template class CWiseMatching<RowLik, ReductionOf<double>>;

  extern template class CWiseCompound<MatrixLik, ReductionOf<RowLik>>;
  extern template class CWiseCompound<MatrixLik, ReductionOf<VectorLik>>;
  extern template class CWiseCompound<RowLik, ReductionOf<double>>;


} // namespace bpp

#endif // DATAFLOW_CWISE_H