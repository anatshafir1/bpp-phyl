#include <Bpp/Phyl/NewLikelihood/DataFlow/ExtendedFloat.h>
#include <Eigen/Core>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlowCWise.h>
#include <Bpp/Phyl/NewLikelihood/DataFlow/Definitions.h>

using namespace bpp;
using namespace std;
using namespace numeric;

int main() {
    MatrixLik mat = MatrixLik::Ones(3,2);
    Eigen::MatrixXd Pijt = Eigen::MatrixXd::Ones(3,3);
    Pijt(0,0) = 0.8;
    Pijt(0,1) = 0.19999999999;
    Pijt(0,2) = 9.999945316252479e-12;
    Pijt(1,0) = 0.37;
    Pijt(1,1) = 9.992007221626409e-15;
    Pijt(1,2) = 0.62999999999999;
    Pijt(2,0) = 0.99;
    Pijt(2,1) = 0.005;
    Pijt(2,2) = 0.005;
    cout << mat << endl;
    mat(0,0) = 2.6e-10;
    mat(0,1) = 0.00001;
    mat(1,0) = 0.02;
    mat(1,1) = 0.0009;
    mat(2,0) = 1.8e-7;
    mat(2,1) = 0.8;
    mat.normalize();
    auto fatherStatePijt = numeric::cwise(Pijt.row(1).transpose());
    auto sonLik = numeric::cwise(mat.col(0));
    auto fatherSonJoint = fatherStatePijt * sonLik;
    auto conditionals = fatherSonJoint/fatherSonJoint.sum();
    std::cout << "***** start *******" << std::endl;
    std::cout << "Pi: " << fatherStatePijt << std::endl;
    std::cout << "SonLik: " << sonLik << std::endl;
    std::cout << "The product: " << fatherSonJoint << std::endl;
    std::cout << "The conditionals: " << conditionals << std::endl;
    std::cout << "Conditionals as float: " << std::endl;
    for (size_t i = 0; i< 3; i++){
        auto conditional = ExtendedFloat(conditionals.float_part()(i), conditionals.exponent_part());
        double converted = ExtendedFloat::convert(conditional);
        std::cout << converted << std::endl;

    }

    std::cout << "*****" << std::endl;
    std::cout << mat << std::endl;
    auto firstSiteArray = numeric::cwise(mat.col(0));
    Eigen::RowVectorXd freqs(3);
    freqs(0) = 0.2;
    freqs(1) = 0.5;
    freqs(2) = 0.3;
    //freqs.normalize();
    auto freqsArray = numeric::cwise(freqs.row(0).transpose());
    std::cout << "*****" << std::endl;
    std::cout << freqsArray << std::endl;
    std::cout << "*****" << std::endl;
    std::cout << firstSiteArray << std::endl;
    ExtendedFloatArrayXd res = freqsArray * firstSiteArray;
    std::cout << "checking ...." << std::endl;
    std::cout << res << std::endl;
    for (size_t i = 0; i< 3; i++){
        std::cout << res(i) << std::endl;

    }
    std::cout << "#####" << std::endl;
    auto sumRes = res.sum(); 
    sumRes.normalize();
    ExtendedFloatArrayXd divRes = res/sumRes;
    divRes.normalize();
    auto exponent = divRes.exponent_part();
    for (size_t i = 0; i < 3; i++){
        auto elem = ExtendedFloat(divRes(i), exponent);
        elem.normalize();
        double elemDouble = ExtendedFloat::convert(elem);
        std::cout << elemDouble << std::endl;


    }



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

