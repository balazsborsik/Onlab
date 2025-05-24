#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <dirent.h>
using namespace std;

/*
g++ -o current_results current_results.cpp
current_results 2 40
*/

int readBestResult(int m, int n) {
    stringstream filestart;
    filestart << "Z" << m << "_" << n << "_" << 2 << "_" << 2 << "_";
    string prefix = filestart.str();

    regex pattern(prefix + R"((\d+)\.txt)");

    DIR* dir = opendir("output");
    if (!dir) {
        cerr << "Failed to open 'output' directory." << endl;
        return -1;
    }

    struct dirent* entry;
    int maxNum = 0;

    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;

        smatch match;
        if (std::regex_match(filename, match, pattern)) {
            int num = std::stoi(match[1]);
            if (num > maxNum) {
                maxNum = num;
            }
        }
    }

    closedir(dir);
    return maxNum;
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

void print_in_latexformat(const vector<vector<int>>& res, int size){
    stringstream str;
    ofstream outfile ("results_latex.txt");
    string start=R"(\begin{center}
\renewcommand{\arraystretch}{1.2}
\setlength{\tabcolsep}{3pt}
\scriptsize
\resizebox{\textwidth}{!}{
\begin{tabular}{r|*{39}{c}} % 34 oszlop n=2..39
$m \backslash n$ & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 & 12 & 13 & 14 & 15 & 16 & 17 & 18 & 19 & 20 & 21 & 22 & 23 & 24 & 25 & 26 & 27 & 28 & 29 & 30 & 31 & 32 & 33 & 34 & 35 & 36 & 37 & 38 & 39 & 40 \\
\hline)";
string end= R"(\end{tabular}
}
\vspace{1em}
\captionof{table}{$Z_{2,2}(m,n)$ alsó becslések}
\end{center})";
    outfile<<start<<endl;
    for(int i=0;i<size;i++){
        outfile<<i+2;
        for(int j=0;j<size;j++){
            int val=res[j][i];
            if(val)
                outfile<<" & "<<val;
            else
                outfile<<" &  ";
        }
        outfile<<" \\\\"<<endl;
    }
    outfile<<end;
    outfile.close();
}

int main(int argc, char* argv[]){
    if (argc !=1 && (argc != 3 || !is_number(argv[1]) || !is_number(argv[2]))) {
        cerr << "Usage: " << argv[0] << " [<min size> <max size>]" << endl;
        return 1;
    }
    int min = 2;
    int max = 40;
    if(argc==3){
        min = stoi( argv[1] );
        max = stoi( argv[2] );
        if (min>max || min<2 || max>50) {
            cerr << "first parameter can't be greater than second parameter and both parameters must be in the range of: 2-50" << endl;
            return 2;
        }
    }
    int size = max - min+1;
    vector<vector<int>> res(size, vector<int>(size, 0));
    for(int i=0;i<size;i++){
        for(int j=i;j<size;j++){
            res[j][i] = readBestResult(min+i, min+j);
        }
    }
    for(int i=0;i<size;i++){
        for(int j=0;j<i+1;j++){
            cout<<res[i][j]<<" ";
        }
        cout<<endl;
    }
    stringstream str;
    ofstream outfile ("current_results.txt");
    for(int i=0;i<size;i++){
        for(int j=0;j<size;j++){
            outfile<<res[i][j]<<" ";
        }
        outfile<<endl;
    }
    outfile.close();
    print_in_latexformat(res,size);

    
}