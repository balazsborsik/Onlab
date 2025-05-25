#include <iostream>
#include <fstream>
#include <sstream>
#include <regex>
#include <dirent.h>
using namespace std;


string readBestGraph(int m, int n) {
    stringstream filestart;
    filestart << "Z" << m << "_" << n << "_" << 2 << "_" << 2 << "_";
    string prefix = filestart.str();

    regex pattern(prefix + R"((\d+)\.txt)");

    DIR* dir = opendir("output");
    if (!dir) {
        cerr << "Failed to open 'output' directory." << endl;
        return "";
    }

    struct dirent* entry;
    int maxNum = -1;
    string result;

    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;

        smatch match;
        if (std::regex_match(filename, match, pattern)) {
            int num = std::stoi(match[1]);
            if (num > maxNum) {
                maxNum = num;
                result = "output\\" + filename;
            }
        }
    }

    closedir(dir);
    return result;
}

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
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
    for(int i=0;i<size;i++){
        for(int j=i;j<size;j++){
            string str_src = readBestGraph(min+i, min+j);
            ifstream src(str_src, std::ios::binary);
            ofstream dst("best_"+str_src, std::ios::binary);

            if (!src || !dst) {
                cerr << "Error opening files.\n";
                return 1;
            }

            dst << src.rdbuf();
        }
    }
}