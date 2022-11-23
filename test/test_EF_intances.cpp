#include <Bpp/Phyl/Likelihood/DataFlow/ExtendedFloat.h>
#include <Eigen/Core>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowNumeric.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlowCWise.h>
#include <Bpp/Phyl/Likelihood/DataFlow/Definitions.h>
#include<set>
#include<iterator>

// just to test regexes extracted from fasta file
#include <string>
#include <vector>
#include <iostream>
#include <time.h>
#include <sys/stat.h>
#include <regex>

using namespace bpp;
using namespace std;
using namespace numeric;

void printSet(set <uint> set_to_print){
    std::cout << "printing set:" << std::endl;
    for (auto it = set_to_print.begin(); it != set_to_print.end(); ++it)
        std::cout << *it << std::endl;
}

void create_pasta_file(string &file_name, int min, int max, std::map<std::string, std::map<string, double>> &sp_states_map){
    ofstream outFile;
    outFile.open(file_name);
    int min_offset = -1;
    int max_offset = 10;
    for (int i = min+min_offset; i < max+max_offset; i++){
        outFile << i << "\t";
    }
    outFile << max+max_offset << std::endl;
    auto it = sp_states_map.begin();
    while (it != sp_states_map.end()){
        outFile <<">" << it->first << std::endl;
        
        auto &states_probs_map = sp_states_map[it->first];
        if (states_probs_map.find("X") != states_probs_map.end()){
            for (int i = min+min_offset; i < max+max_offset; i++){
                outFile << "1.0\t";
            }
            outFile << "1.0" << std::endl;
        }else{
            auto it_sp_states = states_probs_map.begin();
            vector<int> states;
            while (it_sp_states != states_probs_map.end()){
                states.push_back(std::stoi(it_sp_states->first));
                it_sp_states ++;
            }
            for (int i = min+min_offset; i <= max+max_offset; i++){
                if (std::find(states.begin(), states.end(), i) != states.end()){
                    string state_str = std::to_string(i);
                    if (i < max+max_offset){
                        outFile << states_probs_map[state_str] << "\t";

                    }else{
                        outFile << states_probs_map[state_str] << std::endl;
                    }
                    

                }else{
                    if (i < max+max_offset){
                        outFile <<"0\t";

                    }else{
                        outFile <<"0" << std::endl;
                    }
                    
                }

            }
        }
        it ++;
    }
    outFile.close();
}

std::map<std::string, std::map<string, double>> extract_alphabet_states(string &file_path, int &min, int &max){
    
    ifstream stream;
    min = 500;
    max = 1;
    std::regex rgx_composite("([\\d]+)=[\\d]+");
    std::regex rgx_prob("[\\d]+=([\\d]+\\.*[\\d]*)");
    std::regex rgx_state("([\\d]+)");
    std::regex rgx_species(">([\\S]+)");
    stream.open(file_path.c_str());
    vector <string> lines = FileTools::putStreamIntoVectorOfStrings(stream);
    stream.close();
    std::string species_name;
    std::map<std::string, std::map<string, double>> sp_with_states;
    for (size_t i = 0; i < lines.size(); i ++){
        vector <string> states;
        vector<double> probs;
        if (lines[i] == ""){
            continue;
        }else if (lines[i].rfind(">", 0) == 0){
            std::smatch match_sp_name;
            regex_search(lines[i], match_sp_name, rgx_species);
            species_name = match_sp_name[1];
            continue;
        }
        std::smatch match_composite;
        std::smatch match_prob;
        std::smatch match_single_state;
        string content = lines[i];
        bool composite = false;
        bool single_state = false;
        //uint state;
        while(regex_search(content, match_composite, rgx_composite))
        {
            composite = true;
            int state = stoi(match_composite[1]);
            if (state < min){
                min =  state;
            }
            if (state > max){
                max = state;
            }
            //state = static_cast<uint>(stoi(match[1]));
            
            states.push_back(match_composite[1]);

            regex_search(content, match_prob, rgx_prob);
            probs.push_back(std::stod(match_prob[1]));
            content = match_composite.suffix();
        }
        if (!composite){
            if (lines[i] == "X"){
                states.push_back("X");
                probs.push_back(1.0);
                single_state = true;
            }else if (regex_search(lines[i], match_single_state, rgx_state)){
                if ((size_t)(match_single_state[1].length()) == (size_t)(lines[i].length())){
                    //state = static_cast<uint>(stoi(match_single_state[1]));
                    int state = stoi(match_single_state[1]);
                    if (state < min){
                        min =  state;
                    }
                    if (state > max){
                        max = state;
                    }
                    states.push_back(match_single_state[1]);
                    single_state = true;
                    probs.push_back(1.0);
                }
            }
            if (!single_state){
                throw Exception("Not a state!!!");
            }
        }
        for (size_t j = 0; j < states.size(); j++){
            sp_with_states[species_name][states[j]] = probs[j];

        }

    }
    return sp_with_states;


}

int main() {
    Eigen::MatrixXd Pijt = Eigen::MatrixXd::Ones(3,3);
    Pijt(0,0) = 0.2;
    Pijt(0,1) = 0.3;
    Pijt(0,2) = 0.5;
    Pijt(1,0) = 0.1;
    Pijt(1,1) = 0.8;
    Pijt(1,2) = 0.1;
    Pijt(2,0) = 0.35;
    Pijt(2,1) = 0.6;
    Pijt(2,2) = 0.05;

    MatrixLik Ln;// = MatrixLik::Ones(2,3);

    ExtendedFloat ef1 = ExtendedFloat{0.25,0};
    ef1.normalize();

    ExtendedFloat ef2 = ExtendedFloat{1.2e-30,0};
    ef2.normalize();
    ef2 *= ef2;
    ExtendedFloat ef4 = ExtendedFloat{1.2e-300,0};
    ef4.normalize();
    ef2 *= ef4;
    std::cout << ef2 << std::endl;
    ExtendedFloat ef3 = ExtendedFloat{0.05,0};
    ef3.normalize();
    std::vector<ExtendedFloat> vec = {ef1, ef2, ef3, ExtendedFloat{0,0}};
    ExtendedFloatVectorXd v;
    copyBppToEigen (vec, v);
    std::cout <<"Mid res:" << v <<std::endl;
    std::cout << "printed element " << v(0);
    std::vector<ExtendedFloatVectorXd> vector_v;
    vector_v.push_back(v);
    copyBppToEigen (vector_v, Ln);
    size_t nrows = Pijt.rows();
    size_t ncols = Ln.cols();
      // if (x0.cols() == 1){
      //   nrows = 1;
      // }
    vector<ExtendedFloatVectorXd> result;
    for (size_t i = 0; i < nrows; i++){
        vector <ExtendedFloat> row_col_res;
        for (size_t j = 0; j < ncols; j++){
            auto y1 = Pijt.row(i).array();
            auto y2 = cwise(Ln.col(j));
            auto product = y1*y2;
            auto max = product.maxCoeff();
            row_col_res.push_back(max);
        }
        ExtendedFloatVectorXd res;
        copyBppToEigen(row_col_res, res);
        result.push_back(res);

    }
    MatrixLik matLik = MatrixLik::Zero(nrows, ncols);
    copyBppToEigen(result, matLik);




      //     if (nrows == 1){
      //       auto y1 = x0.col(i).transpose().array();
      //       auto y2 = (x1.col(j).transpose()).array();
      //       auto prod = y1 * y2;
      //       result (i, j) = prod.maxCoeff();
      //     }else{
      //       auto y1 = x0.row(i).array();
      //       auto y2 = (x1.col(j).transpose()).array();
      //       auto prod = y1 * y2;
      //       result (i, j) = prod.maxCoeff();
      //     }

      //   }
      // }


    //Ln.normalize();



    // Ln = mat_double.unaryExpr ([](double d) {
    //     ExtendedFloat ef{d, 0};
    //     return d;
    // });
    // Ln = mat_int.NullaryExpr([&Ln](int i){
    //     Ln()
    // })
    // Ln = Ln.binaryExpr(mat_int),
    //         ([](ExtendedFloat d, int i){
    //             return ExtendedFloat{d.get_float_part(), i};});


    
    // Ln = Pijt.unaryExpr ([](double d) {
    //       ExtendedFloat ef{d};
    //       ef.normalize ();
    //       return d;
    //     });

    // ef3.normalize();
    // Eigen::Matrix<ExtendedFloat, 1, Eigen::Dynamic> temp;
    // temp = Ln.unaryExpr ([](ExtendedFloat d) {
    //     //   ExtendedFloat ef{d};
    //     //   ef.normalize ();
    //       return d;
    //     });




    std::cout << "Pijt:" << std::endl;
    std::cout << Pijt << std::endl;

    std::cout << "Ln: " <<  Ln << std::endl;
    std::cout << "Max result: " << std::endl;
    std::cout << matLik << std::endl;

    // string path = "/home/anat/Docs/Sida/counts_f.fasta";
    // string outFileName = "/home/anat/Docs/Sida/counts_p.pasta";
    // int min;
    // int max;
    // auto map_species_states = extract_alphabet_states(path, min, max);
    // auto it = map_species_states.begin();
    // while (it != map_species_states.end()){
    //     std::cout << "species: " << it->first << std::endl;
    //     auto &map_states_probs = map_species_states[it->first];
    //     auto it_probs = map_states_probs.begin();
    //     while (it_probs != map_states_probs.end()){
    //         std::cout << "\t" << it_probs->first << " " << map_states_probs[it_probs->first] << std::endl;
    //         it_probs ++;
    //     }
    //     it++;
    // }
    // std::cout << "****" << std::endl;
    // std::cout << "min is: " << min << std::endl;
    // std::cout << "max is: " << max << std::endl;
    // create_pasta_file(outFileName, min, max, map_species_states);



    // vector<uint> nodes = {1,2,5,7,8,6,9,7,7,9,11,10};
    // vector <uint> vectorOfNodesM2 = {7,5,5,6,12};
    // set<uint> setNodes(nodes.begin(), nodes.end());
    // printSet(setNodes);
    // set<uint> setNodesM2(vectorOfNodesM2.begin(), vectorOfNodesM2.end());
    // printSet(setNodesM2);
    // set<uint> intersect;
    // set_intersection(setNodes.begin(), setNodes.end(), setNodesM2.begin(), setNodesM2.end(),
    //              std::inserter(intersect, intersect.begin()));
    // printSet(intersect);
    // std::cout << "size: " << intersect.size() << std::endl;
    // std::cout << intersect << std::endl;
    // MatrixLik matrixLik = MatrixLik::Ones(5, 2);
    // int rows = matrixLik.rows();
    // int cols = matrixLik.cols();
    // std::cout << rows << "," << cols << std::endl;
    
    // matrixLik(0, 0) = 1;
    // matrixLik(1, 0) = 0.01;
    // matrixLik(2, 0) = 0;
    // matrixLik(3, 0) = -0.4;
    // matrixLik(4, 0) = 0.1;
    // matrixLik(0, 1) = 0;
    // matrixLik(1, 1) = -0.003;
    // matrixLik(2, 1) = 0.005;
    // matrixLik(3, 1) = -0.0002;
    // matrixLik(4, 1) = 1.2;
    // std::cout << matrixLik << std::endl;
    // auto floatPart = matrixLik.float_part();
    // std::cout << floatPart << std::endl;
    
    // // std::cout << arr.float_part() << std::endl;
    // auto absArr = floatPart.cwiseAbs();
    // std::cout << "abs: " << absArr << std::endl;
    // int minIndex_i;
    // int minIndex_j;

    // int maxIndex_i;
    // int maxIndex_j;
    // auto minValue = absArr.minCoeff(&minIndex_i, &minIndex_j);
    // auto maxValue = absArr.maxCoeff(&maxIndex_i, &maxIndex_j);
    // std::cout << "min index row: " << minIndex_i << " min index col: "<< minIndex_j << " , min value: " << minValue << std::endl;
    // std::cout << "max index row: " << maxIndex_i << " max index col: "<< maxIndex_j << " , max value: " << maxValue << std::endl;
    // //auto arr_no_zeros = (absArr.array() > 0);
    // //std::cout << arr_no_zeros << std::endl;
    // int minIndex_k;
    // int minIndex_l;
    // ((absArr.array() == 0).select(maxValue, absArr)).minCoeff(&minIndex_k, &minIndex_l);
    // std::cout << "After selecting: " << absArr << std::endl;
    // //std::cout <<"minimum is: " <<  result << std::endl;
    // std::cout << "min index row: " << minIndex_k << " max index col: " << minIndex_l << std::endl;

    // auto minNonZeroValue = maxValue;
    // //std::cout << "min non zero value: " << minNonZeroValue << std::endl;
    // int minNonZeroIndex = maxIndex;
    // for (int i = 0; i < (int)absArr.size(); i++){
    //     if ((absArr(i) < minNonZeroValue) && (absArr(i) > 0)){
    //         minNonZeroIndex = i;
    //         minNonZeroValue = absArr(i);

    //     }
    // }
    // std::cout << "index: " << minNonZeroIndex << std::endl;
    // std::cout << "array: " << absArr << std::endl;
    // std::cout << "min non zero value is: " << absArr(minNonZeroIndex) << std::endl;



    // MatrixLik mat = MatrixLik::Ones(3,2);
    // Eigen::MatrixXd Pijt = Eigen::MatrixXd::Ones(3,3);
    // Pijt(0,0) = 0.8;
    // Pijt(0,1) = 0.19999999999;
    // Pijt(0,2) = 9.999945316252479e-12;
    // Pijt(1,0) = 0.37;
    // Pijt(1,1) = 9.992007221626409e-15;
    // Pijt(1,2) = 0.62999999999999;
    // Pijt(2,0) = 0.99;
    // Pijt(2,1) = 0.005;
    // Pijt(2,2) = 0.005;
    // cout << mat << endl;
    // mat(0,0) = 2.6e-10;
    // mat(0,1) = 0.00001;
    // mat(1,0) = 0.02;
    // mat(1,1) = 0.0009;
    // mat(2,0) = 1.8e-7;
    // mat(2,1) = 0.8;
    // mat.normalize();
    // auto fatherStatePijt = numeric::cwise(Pijt.row(1).transpose());
    // auto sonLik = numeric::cwise(mat.col(0));
    // auto fatherSonJoint = fatherStatePijt * sonLik;
    // auto conditionals = fatherSonJoint/fatherSonJoint.sum();
    // std::cout << "***** start *******" << std::endl;
    // std::cout << "Pi: " << fatherStatePijt << std::endl;
    // std::cout << "SonLik: " << sonLik << std::endl;
    // std::cout << "The product: " << fatherSonJoint << std::endl;
    // std::cout << "The conditionals: " << conditionals << std::endl;
    // std::cout << "Conditionals as float: " << std::endl;
    // for (size_t i = 0; i< 3; i++){
    //     auto conditional = ExtendedFloat(conditionals.float_part()(i), conditionals.exponent_part());
    //     double converted = ExtendedFloat::convert(conditional);
    //     std::cout << converted << std::endl;

    // }

    // std::cout << "*****" << std::endl;
    // std::cout << mat << std::endl;
    // auto firstSiteArray = numeric::cwise(mat.col(0));
    // Eigen::RowVectorXd freqs(3);
    // freqs(0) = 0.2;
    // freqs(1) = 0.5;
    // freqs(2) = 0.3;
    // //freqs.normalize();
    // auto freqsArray = numeric::cwise(freqs.row(0).transpose());
    // std::cout << "*****" << std::endl;
    // std::cout << freqsArray << std::endl;
    // std::cout << "*****" << std::endl;
    // std::cout << firstSiteArray << std::endl;
    // ExtendedFloatArrayXd res = freqsArray * firstSiteArray;
    // std::cout << "checking ...." << std::endl;
    // std::cout << res << std::endl;
    // for (size_t i = 0; i< 3; i++){
    //     std::cout << res(i) << std::endl;

    // }
    // std::cout << "#####" << std::endl;
    // auto sumRes = res.sum(); 
    // sumRes.normalize();
    // ExtendedFloatArrayXd divRes = res/sumRes;
    // divRes.normalize();
    // auto exponent = divRes.exponent_part();
    // for (size_t i = 0; i < 3; i++){
    //     auto elem = ExtendedFloat(divRes(i), exponent);
    //     elem.normalize();
    //     double elemDouble = ExtendedFloat::convert(elem);
    //     std::cout << elemDouble << std::endl;


    // }



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

