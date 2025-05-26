#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <dirent.h>
#include <vector>
#include <string>
#include <algorithm> 
using namespace std;

/*
g++ -o current_results current_results.cpp
current_results 2 40
*/

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void print_in_latexformat(const vector<vector<int>>& res, int size, const vector<vector<int>>& file1, const vector<vector<int>>& file3){
    ofstream outfile ("tan_latex.txt");
    string start=R"(\begin{center}
\renewcommand{\arraystretch}{1.2}
\setlength{\tabcolsep}{3pt}
\scriptsize
\resizebox{\textwidth}{!}{
\begin{tabular}{r|*{34}{c}} % 34 oszlop n=2..39
$m \backslash n$ & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 & 14 & 15 & 16 & 17 & 18 & 19 & 20 & 21 & 22 & 23 & 24 & 25 & 26 & 27 & 28 & 29 & 30 & 31 & 32 & 33 & 34 & 35 \\
\hline)";
string end= R"(\end{tabular}
}
\vspace{1em}
\captionof{table}{$Z_{2,2}(m,n)$ alsó becslések}
\end{center})";
    outfile<<start<<endl;
    for(int i=0;i<size-16;i++){
        outfile<<i+2;
        for(int j=0;j<size-5;j++){
            int val=res[j][i];
            if (file1[j][i] == 0) {
                if(res[j][i] == file3[j][i]){
                
                stringstream str;
                if(val)
                outfile<<" & \\greenbf{"<<res[j][i]<<'}';
                else outfile<<" &  ";
                //str<<"(" << j+2 << ", " << i+2 << ") = " << res[i][j];
                //std::cout << str.str() << "\n";
                //std::cout<<file1[i][j] << file2[i][j] << file3[i][j]<<'\n';
                }else{
                if(val)
                    outfile<<" & "<<val;
                else
                    outfile<<" &  ";
                }
            }else{
            if(val)
                outfile<<" & \\textbf{"<<val<<'}';
            else
                outfile<<" &  ";
            }
        }
        outfile<<" \\\\"<<endl;
    }
    outfile<<end;
    outfile.close();
}



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
    print_in_latexformat(file2, 39, file1, file3);
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