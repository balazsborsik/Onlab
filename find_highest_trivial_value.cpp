#include <iostream>
#include <fstream>
#include <sstream>
#include <queue>
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

vector<vector<int>> create_from_file(int m, int n, string filename){
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
    return graph;
}

int CurrentBestValue(int m, int n) {
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
    return maxNum;
}

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

int addVertexToM(vector<vector<int>> &adj, const vector<vector<int>> &inputgraph){
    int m=adj.size();
    int n=adj[0].size();
    for(int i=0;i<m-1;i++){
        for(int j=0;j<n;j++){
            adj[i][j]=inputgraph[i][j];
        }
    }

    vector<int> temp(n, 0);
    for(int i=0;i<n;i++){
        temp[i]=i;
    }
    vector<int> chosenvertices;
    int best_size=0;
    for(int k=0;k<m*n/2+5;k++){
        random_shuffle(begin(temp), end(temp));
        vector<int> vertices;
        vector<int> neighbours(m, 0);
        for(int i=0;i<n;i++){
            bool addvertex=true;
            int vertex=temp[i];
            for(int j=0;j<m;j++){
                if(adj[j][vertex]==1&&neighbours[j]==1){
                    addvertex=false;
                    break;
                }
            }
            if(addvertex){
                for(int j=0;j<m;j++){
                    if(adj[j][vertex]==1)
                        neighbours[j]=1;
                }
                vertices.push_back(vertex);
            }
        }
        if(vertices.size()>best_size){
            best_size=vertices.size();
            chosenvertices=vertices;
        }
    }
    for(const auto& elm : chosenvertices){
        adj[m-1][elm]=1;
    }
    return chosenvertices.size();
}

int addVertexToN(vector<vector<int>> &adj, const vector<vector<int>> &inputgraph){
    int m=adj.size();
    int n=adj[0].size();
    for(int i=0;i<m;i++){
        for(int j=0;j<n-1;j++){
            adj[i][j]=inputgraph[i][j];
        }
    }

    vector<int> temp(m, 0);
    for(int i=0;i<m;i++){
        temp[i]=i;
    }
    vector<int> chosenvertices;
    int best_size=0;
    for(int k=0;k<m*n/2+5;k++){
        random_shuffle(begin(temp), end(temp));
        vector<int> vertices;
        vector<int> neighbours(n, 0);
        for(int i=0;i<m;i++){
            bool addvertex=true;
            int vertex=temp[i];
            for(int j=0;j<n;j++){
                if(adj[vertex][j]==1&&neighbours[j]==1){
                    addvertex=false;
                    break;
                }
            }
            if(addvertex){
                for(int j=0;j<n;j++){
                    if(adj[vertex][j]==1)
                        neighbours[j]=1;
                }
                vertices.push_back(vertex);
            }
        }
        if(vertices.size()>best_size){
            best_size=vertices.size();
            chosenvertices=vertices;
        }
    }
    for(const auto& elm : chosenvertices){
        adj[elm][n-1]=1;
    }
    return chosenvertices.size();
}

int createStartingGraphFromInput(vector<vector<int>> &adj, const vector<vector<int>> &inputgraph){
    int m=adj.size();
    int n=adj[0].size();
    if(m-1 == inputgraph.size()&&n == inputgraph[0].size())
        return addVertexToM(adj, inputgraph);
    else if(m == inputgraph.size()&&n-1 == inputgraph[0].size())
        return addVertexToN(adj, inputgraph);
    else __throw_invalid_argument("az input és az adj gráfok nem megfelelő méretűek");
}

pair<int, int> min_degree(const vector<vector<int>> &inputgraph, int m, int n){ // returns the min_degree like this (min_in_M, min_in_N)
    vector<int> m_degree(m, 0);
    vector<int> n_degree(n, 0);
    for(int i=0; i<m;i++){
        for(int j=0;j<n;j++){
            if(inputgraph[i][j]){
                m_degree[i]++;
                n_degree[j]++;
            }
        }
    }
    int m_min=500;
    int n_min=500;
    for(const auto& elm : m_degree){
        if(elm<m_min) m_min=elm;
    }
    for(const auto& elm : n_degree){
        if(elm<n_min) n_min=elm;
    }
    return make_pair(m_min, n_min);
}

void backwards(){
    int n, m;
    queue<pair<int,int>> n_m_queue;
    n_m_queue.emplace(13,18);
    vector<vector<int>> done(41, vector<int>(41, 0));
    while(!n_m_queue.empty()){
        m=n_m_queue.front().first;
        n=n_m_queue.front().second;
        if(m<=2||n<=2||m>=40||n>=40||done[m][n]){n_m_queue.pop(); continue;}
        done[m][n]=1;
        cout<<"Z("<<m<<","<<n<<")"<<endl;
        vector<vector<int>> graph_m(m+1, vector<int>(n, 0));
        vector<vector<int>> graph_n(m, vector<int>(n+1, 0));
        vector<vector<int>> inputgraph;
        pair<int, int> edges_removed_packed;
        int edges_removed;
        int edges_needed;
        
        edges_needed=CurrentBestValue(m,n) - CurrentBestValue(m-1,n);
        inputgraph=create_from_file(m,n,readBestGraph(m,n));
        edges_removed_packed = min_degree(inputgraph, m, n);
        edges_removed=edges_removed_packed.first;
        if(edges_needed>=edges_removed){
            if(edges_needed>edges_removed)
                cout<<"M_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m-1,n);
        }
        if(m!=n){
        edges_needed=CurrentBestValue(m,n) - CurrentBestValue(m,n-1);
        edges_removed=edges_removed_packed.second;

        if(edges_needed>=edges_removed){
            if(edges_needed>edges_removed)
                cout<<"N_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m,n-1);
        }}
    }
}

int main() {
    backwards(); return 0;
    srand(time(0));
    int n, m;
    queue<pair<int,int>> n_m_queue;
    n_m_queue.emplace(3,3);
    vector<vector<int>> done(41, vector<int>(41, 0));
    while(!n_m_queue.empty()){
        m=n_m_queue.front().first;
        n=n_m_queue.front().second;
        if(m>=40||n>=40||done[m][n]){n_m_queue.pop(); continue;}
        
        done[m][n]=1;
        if(n>10&&m>10)
        cout<<"Z("<<m<<","<<n<<")"<<endl;
        vector<vector<int>> graph_m(m+1, vector<int>(n, 0));
        vector<vector<int>> graph_n(m, vector<int>(n+1, 0));
        vector<vector<int>> inputgraph;
        int edges_added;
        int edges_needed;
        if(m!=n){
        int edges_needed=CurrentBestValue(m+1,n) - CurrentBestValue(m,n);
        inputgraph=create_from_file(m,n,readBestGraph(m,n));
        edges_added = createStartingGraphFromInput(graph_m,inputgraph);

        if(edges_needed<=edges_added){
            if(edges_needed<edges_added)
                cout<<"M_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m+1,n);
        }}

        edges_needed=CurrentBestValue(m,n+1) - CurrentBestValue(m,n);
        inputgraph=create_from_file(m,n,readBestGraph(m,n));
        edges_added = createStartingGraphFromInput(graph_n,inputgraph);

        if(edges_needed<=edges_added){
            if(edges_needed<edges_added)
                cout<<"N_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m,n+1);
        }
    }

    return 0;
}