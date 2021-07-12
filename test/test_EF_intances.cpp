#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloat.h>
#include <Eigen/Core>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWise.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Definitions.h>

using namespace bpp;
using namespace std;
using namespace numeric;

int main() {
    MatrixLik mat = MatrixLik::Ones(3,3);
    cout << mat << endl;
    ExtendedFloat ef = ExtendedFloat{1.8e-10};
    //mat(0,0) = ef;
    mat(0,0) = ExtendedFloat::convert(ef);
    mat(0,1) = 1.5;
    mat(0,2) = 1.9e-30;
    mat(1,0) = 1.8e-50;
    mat.normalize();
    MatrixLik mat2 = MatrixLik::Ones(3,1);
    mat2(0,0) = 2.2e-60;
    mat2(1,0) = 0.01;
    mat2(2,0) = 0.0005;
    mat2.normalize();
    auto res = maxProduct(mat, mat2);
    cout << "first matrix: "<< mat << endl;
    cout << "Second matrix: " << mat2 << endl;
    cout << "result: "<< res << endl;

    return 0;
}


      // const auto & x0 = accessValueConstCast<DepT0> (*this->dependency (0));
      // const auto & x1 = accessValueConstCast<DepT1> (*this->dependency (1));
      // size_t nrows = x0.rows();  
      // size_t ncols = x1.cols();
      // result = zero (targetDimension_);
      // if (x0.cols() == 1){
      //   nrows = 1;
      // }
      // for (size_t i = 0; i < nrows; i++){
      //  for (size_t j = 0; j < ncols; j++){
      //    if (nrows == 1){
      //       auto y1 = cwise(x0.col(i).transpose());
      //       auto y2 = cwise(x1.col(j).transpose());
      //       auto prod = y1 * y2;
      //       auto maxRes = ExtendedFloat::convert(prod.maxCoeff());
      //       result (i, j) = maxRes;


      //    }else{
      //       auto y1 = cwise(x0.row(i));
      //       auto y2 = cwise(x1.col(j).transpose());
      //       auto prod = y1 * y2;
      //       auto maxRes = ExtendedFloat::convert(prod.maxCoeff());
      //       result (i, j) = maxRes;

