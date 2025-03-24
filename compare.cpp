#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
/*
g++ -o compare compare.cpp
compare DYNP MT
*/


using namespace std;

const string RESET = "\033[0m";  // Text Reset

const string RED = "\033[0;31m";     // RED

const string GREEN = "\033[0;32m";   // GREEN

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

// Function to compute element-wise difference of two matrices
vector<vector<int>> subtractMatrices(const vector<vector<int>>& mat1, const vector<vector<int>>& mat2) {
    int rows = max(mat1.size(), mat2.size());
    vector<vector<int>> result(rows);

    for (int i = 0; i < rows; i++) {
        int cols = max(i < mat1.size() ? mat1[i].size() : 0, i < mat2.size() ? mat2[i].size() : 0);
        result[i].resize(cols, 0);

        for (int j = 0; j < cols; j++) {
            int val1 = (i < mat1.size() && j < mat1[i].size()) ? mat1[i][j] : 0;
            int val2 = (i < mat2.size() && j < mat2[i].size()) ? mat2[i][j] : 0;
            result[i][j] = val1 - val2;
        }
    }
    return result;
}

// Function to print matrix to console
void printMatrix(const vector<vector<int>>& matrix) {
    for (const auto& row : matrix) {
        for (int val : row) {
            if(val!=0){
                if(val>0){
                    cout << GREEN << val << RESET << " ";
                }else{
                    cout << RED << val << RESET << " ";
                }
            }else{
                cout << val << " ";
            }
        }
        cout << endl;
    }
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input1.txt> <input2.txt>" << endl;
        return 1;
    }

    string ending = "_results.txt";
    string file1 = argv[1]+ending; // First input file
    string file2 = argv[2]+ending; // Second input file

    vector<vector<int>> matrix1 = readMatrix(file1);
    vector<vector<int>> matrix2 = readMatrix(file2);

    vector<vector<int>> result = subtractMatrices(matrix1, matrix2);

    cout << "Difference Matrix:\n";
    printMatrix(result); // Print to console instead of writing to file

    return 0;
}
