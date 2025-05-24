#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm> 

using Matrix = std::vector<std::vector<int>>;
using namespace std;
// Function to read a matrix from a file
vector<vector<int>> readMatrix(const string& filename) {
    vector<vector<int>> matrix;
    ifstream file(filename);
    string line;

    if (!file) {
        cerr << "Error opening file: " << filename << endl;
        exit(1);
    }

    while (getline(file, line)) {
        stringstream ss(line);
        vector<int> row;
        int num;
        while (ss >> num) {
            row.push_back(num);
        }
        matrix.push_back(row);
    }

    file.close();
    return matrix;
}

int main() {
    Matrix file1 = readMatrix("Tan_fix_results.txt");
    Matrix file2 = readMatrix("Tan_results.txt");
    Matrix file3 = readMatrix("current_results.txt");

    if (file1.empty() || file2.empty() || file3.empty()) {
        std::cerr << "One or more files could not be read.\n";
        return 1;
    }

    // Check all three matrices are the same size
    size_t rows = file1.size();
    size_t cols = file1[0].size();

    for (size_t i = 0; i < rows; ++i) {
        if (file1[i].size() != cols || file2[i].size() != cols || file3[i].size() != cols) {
            std::cerr << "Mismatch in column sizes at row "<< file3[i].size()<<file2[i].size()<<file1[i].size()<< i << "\n";
            return 1;
        }
    }
    vector<pair<int, string>> res;
    std::cout << "Matching zero positions (i, j):\n";
    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            if (file1[i][j] == 0 && file2[i][j] == file3[i][j]) {
                if(j<i){
                stringstream str;
                str<<"(" << j+2 << ", " << i+2 << ") = " << file2[i][j];
                std::cout << str.str() << "\n";
                //std::cout<<file1[i][j] << file2[i][j] << file3[i][j]<<'\n';
                res.emplace_back(j, str.str());
                }
            }
        }
    }
    cout << endl<<endl;
    sort(res.begin(), res.end());
    for(const auto& elm : res){
        cout<<elm.second<<endl;
    }
    return 0;
}
