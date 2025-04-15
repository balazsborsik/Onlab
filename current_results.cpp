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

int main(int argc, char* argv[]){
    if (argc != 3 || !is_number(argv[1]) || !is_number(argv[2])) {
        cerr << "Usage: " << argv[0] << " <min size> <max size>" << endl;
        return 1;
    }
    int min = stoi( argv[1] );
    int max = stoi( argv[2] );
    int size = max - min+1;
    if (min>max || min<2 || max>50) {
        cerr << "first parameter can't be greater than second parameter and both parameters must be in the range of: 2-50" << endl;
        return 2;
    }
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
        for(int j=0;j<i+1;j++){
            outfile<<res[i][j]<<" ";
        }
        outfile<<endl;
    }
    outfile.close();
}