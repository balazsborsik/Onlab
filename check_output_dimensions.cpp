#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <tuple>
#include <algorithm>
#include <set>
#include <random>
#include <chrono>
#include <regex>
#include <dirent.h>
#include <cstdio>
using namespace std;

vector<string> problematic_filenames;

int linenum(string filename){
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "Could not open the file.\n";
        return 1;
    }

    std::string line;
    int lineCount = 0;

    while (std::getline(file, line)) {
        ++lineCount;
    }
    return lineCount;
}

void evaluateAllGraphs(string directory, int m, int n) {
    stringstream filestart;
    filestart << "Z" << m << "_" << n << "_" << 2 << "_" << 2 << "_";
    string prefix = filestart.str();

    regex pattern(prefix + R"((\d+)\.txt)");

    DIR* dir = opendir(directory.c_str());
    if (!dir) {
        cerr << "Failed to open '"+directory+"' directory." << endl;
        __throw_runtime_error("Failed to open 'output' directory.");
    }

    struct dirent* entry;

    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;

        smatch match;
        if (std::regex_match(filename, match, pattern)) {
            int num = std::stoi(match[1]);
            if(linenum(directory+"\\"+filename)!=m)
                problematic_filenames.emplace_back(directory+"\\"+filename);
        }
    }
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

int main(int argc, char* argv[]) {
    if (argc < 3 || argc > 5 || argc == 4 || ( argc == 5 && ( !is_number(argv[3]) || !is_number(argv[4]) ) ) || ((string)argv[2] != "-l" && (string)argv[2] != "-d")) {
        cerr << "Usage: " << argv[0] << " <directory> <-l (to list) or -d (to delete)> [<min size> <max size>]" <<endl;
        return 1;
    }
    int min = 2;
    int max = 40;
    if(argc==5){
        int min = stoi( argv[3] );
        int max = stoi( argv[4] );
        if (min>max || min<2 || max>50) {
            cerr << "third parameter can't be greater than fourth parameter and both parameters must be in the range of: 2-50" << endl;
            return 2;
        }
    }
    string dir = argv[1];
    for(int i=min;i<=max;i++){
        for(int j=i;j<=max;j++){
            cout<< i<<", "<<j<<endl;
            evaluateAllGraphs(dir ,i, j);
        }
    }
    for(const auto& elm : problematic_filenames){
        cout<< elm << endl;
        if((string)argv[2]=="-d"){
            if (std::remove(elm.c_str()) == 0) {
                std::cout << "File deleted successfully.\n";
            } else {
                std::perror("Error deleting file");
            }
        }
    }

    if(problematic_filenames.empty())
        cout<< "All the files in the given directory are in the defined format";
    return 0;
}