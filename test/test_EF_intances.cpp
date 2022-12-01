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
    Eigen::MatrixXd Pijt = Eigen::RowVectorXd::Ones(3);
    Pijt(0,0) = 0.2;
    Pijt(0,1) = 0.3;
    Pijt(0,2) = 0.5;
    // Pijt(1,0) = 0.1;
    // Pijt(1,1) = 0.8;
    // Pijt(1,2) = 0.1;
    // Pijt(2,0) = 0.35;
    // Pijt(2,1) = 0.6;
    // Pijt(2,2) = 0.05;

    MatrixLik Ln;// = MatrixLik::Ones(2,3);

    ExtendedFloat ef1 = ExtendedFloat{0.25,0};
    ef1.normalize();

    ExtendedFloat ef2 = ExtendedFloat{0.44,0};
    std::cout << ef2 << std::endl;
    ExtendedFloat ef3 = ExtendedFloat{0.05,0};
    ExtendedFloat ef5 = ExtendedFloat{0.5,1000};
    ef5.normalize();
    ef3.normalize();
    ef2.normalize();
    std::vector<ExtendedFloat> vec = {ef5, ef2, ef3};
    ExtendedFloat ef4 = ExtendedFloat{0.5, -1050};
    std::cout << "converted is: " << ExtendedFloat::convert(ef4) << std::endl;
    std::cout << "converted is: " << ExtendedFloat::convert(ef5) << std::endl;
    ef4.normalize();
    std::vector<ExtendedFloat> vec2 = {ef1*ef1, ef4, ef3*ef3};
    ExtendedFloatVectorXd v;
    ExtendedFloatVectorXd v2;
    copyBppToEigen (vec, v);
    copyBppToEigen (vec2, v2);
    std::cout <<"Mid res:" << v <<std::endl;
    std::cout <<"Mid res:" << v2 <<std::endl;
    std::cout << "printed element " << v(0);
    std::vector<ExtendedFloatVectorXd> vector_v;
    vector_v.push_back(v);
    vector_v.push_back(v2);
    copyBppToEigen (vector_v, Ln);
    std::cout << "results vector ### " << Ln << std::endl;
    std::cout << "Converted is:" << std::endl;
    for (size_t i = 0; i < Ln.rows(); i++){
        for (size_t j = 0; j < Ln.cols(); j++){
            std::cout << ExtendedFloat::convert(Ln(i, j)) << std::endl;
        }
    }
    Ln.normalize();
    size_t nrows = Pijt.rows();
    size_t ncols = Ln.cols();
    vector<ExtendedFloatVectorXd> result;
    for (size_t i = 0; i < nrows; i++){
        vector <ExtendedFloat> row_col_res;
        for (size_t j = 0; j < ncols; j++){
            auto y1 = Pijt.row(i).array().transpose();
            auto y2 = cwise(Ln.col(j));
            std::cout <<"****************" << std::endl;
            std::cout << "y1: " << y1 << std::endl;
            std::cout << "y2: " << y2 << std::endl;
            auto product = y1*y2;
            std::cout << "Product of " << i << "," << j << ": " << product << std::endl;
            auto max = product.maxCoeff();
            row_col_res.push_back(max);
        }
        ExtendedFloatVectorXd res;
        copyBppToEigen(row_col_res, res);
        result.push_back(res);

    }
    MatrixLik matLik;
    MatrixLik MatLik_final;
    copyBppToEigen(result, matLik);
    std::cout << "rows: " << matLik.rows() << std::endl;
    std::cout << "cols: " << matLik.cols() << std::endl;

    std::cout << "Pijt:" << std::endl;
    std::cout << Pijt << std::endl;

    std::cout << "Ln: " <<  Ln << std::endl;
    std::cout << "Max result: " << std::endl;
    MatLik_final = matLik.transpose();
    std::cout << MatLik_final << std::endl;


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

