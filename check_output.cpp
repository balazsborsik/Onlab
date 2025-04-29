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
using namespace std;

vector<string> problematic_filenames;

// Function to check if adding edge (u, v) creates a K_{2,2}
bool createsK22(const vector<vector<int>>& adj, int m, int n, int u, int v) {
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            // Now, look for a second common neighbor
            for (int v2 = 0; v2 < n; ++v2) {
                if (v2 == v) continue;
                if (adj[u][v2] && adj[u2][v2]) {
                    return true; // Found a K_{2,2}
                }
            }
        }
    }
    return false;
}


bool has_c4(int m, int n, string filename){
    vector<vector<int>> graph(m, vector<int>(n, 0));
    std::ifstream file(filename);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            int edge=0;
            file>>edge;
            graph[i][j]=edge;
        }
    }
    file.close();
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            if(graph[i][j]){
                graph[i][j]=0;
                if(createsK22(graph, m, n, i, j)){
                    return true;
                }
                graph[i][j]=1;
            }
        }
    }
    return false;
}

void evaluateAllGraphs(int m, int n) {
    stringstream filestart;
    filestart << "Z" << m << "_" << n << "_" << 2 << "_" << 2 << "_";
    string prefix = filestart.str();

    regex pattern(prefix + R"((\d+)\.txt)");

    DIR* dir = opendir("output");
    if (!dir) {
        cerr << "Failed to open 'output' directory." << endl;
        __throw_runtime_error("Failed to open 'output' directory.");
    }

    struct dirent* entry;

    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;

        smatch match;
        if (std::regex_match(filename, match, pattern)) {
            int num = std::stoi(match[1]);
            if(has_c4(m, n, "output\\" + filename))
                problematic_filenames.emplace_back(filename);
        }
    }
}


int main() {
    cout<< "alma";
    for(int i=2;i<=40;i++){
        for(int j=i;j<=40;j++){
            cout<< i<<", "<<j<<endl;
            evaluateAllGraphs(i, j);
        }
    }
    for(const auto& elm : problematic_filenames){
        cout<< elm << endl;
    }
    return 0;
}