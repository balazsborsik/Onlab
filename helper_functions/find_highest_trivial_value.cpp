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
    for(int k=0;k<40*40/2+5;k++){
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
    for(int k=0;k<40*40/2+5;k++){
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

pair<pair<int, int>, pair<int, int>> min_degree(const vector<vector<int>> &inputgraph, int m, int n){ // returns the min_degree like this (min_in_M, min_in_N)
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
    int n_vertex=-1;
    int m_vertex=-1;
    for(int i=0;i<m;i++){
        if(m_degree[i]<m_min){ m_min=m_degree[i]; m_vertex=i;}
    }
    for(int i=0;i<n;i++){
        if(n_degree[i]<n_min){ n_min=n_degree[i]; n_vertex=i;}
    }
    return make_pair(make_pair(m_vertex, m_min), make_pair(n_vertex, n_min));
}

void backwards(){
    int n=19, m=13;
    queue<pair<int,int>> n_m_queue;
    n_m_queue.emplace(m,n);
    vector<vector<vector<vector<int>>>> graphs(41, vector<vector<vector<int>>>(41, vector<vector<int>>(0, vector<int>(0, 0))));
    vector<vector<int>>gra=create_from_file(m,n,readBestGraph(m,n));
    graphs[m][n]=gra;
    vector<vector<int>> done(41, vector<int>(41, 0));
    while(!n_m_queue.empty()){
        m=n_m_queue.front().first;
        n=n_m_queue.front().second;
        if(m<=2||n<=2||m>=40||n>=40||done[m][n]){n_m_queue.pop(); continue;}
        done[m][n]=1;
        cout<<"Z("<<m<<","<<n<<")"<<endl;
        vector<vector<int>> inputgraph=graphs[m][n];
        pair<pair<int, int>, pair<int, int>> edges_removed_packed;
        pair<int, int> edges_removed;
        int edges_needed;
        
        edges_needed=CurrentBestValue(m,n) - CurrentBestValue(m-1,n);
        //inputgraph=create_from_file(m,n,readBestGraph(m,n));
        edges_removed_packed = min_degree(inputgraph, m, n);
        edges_removed=edges_removed_packed.first;
        if(edges_needed>=edges_removed.second){
            if(edges_needed>edges_removed.second)
                cout<<"M_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m-1,n);
            vector<vector<int>> res=inputgraph;
            res.erase(res.begin()+edges_removed.first);
            graphs[m-1][n]=res;
        }
        if(m!=n){
        edges_needed=CurrentBestValue(m,n) - CurrentBestValue(m,n-1);
        edges_removed=edges_removed_packed.second;

        if(edges_needed>=edges_removed.second){
            if(edges_needed>edges_removed.second)
                cout<<"N_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m,n-1);
            vector<vector<int>> res=inputgraph;
            for(int i=0;i<m;i++){
                res[i].erase(res[i].begin()+edges_removed.first);
            }
            graphs[m][n-1]=res;
        }}
    }
    for(const auto& elm : done){
        for(const auto& elm2 : elm){
            cout<<elm2<<" ";
        }
        cout<<endl;
    }
}

void forward(int n_n, int m_m){
    int n=n_n, m=m_m;
    queue<pair<int,int>> n_m_queue;
    n_m_queue.emplace(m,n);
    vector<vector<vector<vector<int>>>> graphs(41, vector<vector<vector<int>>>(41, vector<vector<int>>(0, vector<int>(0, 0))));
    graphs[m][n] = create_from_file(m,n,readBestGraph(m,n));
    vector<vector<int>> done(41, vector<int>(41, 0));
    vector<vector<pair<int, int>>> parent(41, vector<pair<int, int>>(41, {0,0}));
    while(!n_m_queue.empty()){
        m=n_m_queue.front().first;
        n=n_m_queue.front().second;
        if(m>=40||n>=40||done[m][n]){n_m_queue.pop(); continue;}
        
        done[m][n]=1;
        if(n>10&&m>10)
        cout<<"Z("<<m<<","<<n<<")"<<endl;
        vector<vector<int>> graph_m(m+1, vector<int>(n, 0));
        vector<vector<int>> graph_n(m, vector<int>(n+1, 0));
        vector<vector<int>> inputgraph=graphs[m][n];
        int edges_added;
        int edges_needed;
        if(m!=n){
        int edges_needed=CurrentBestValue(m+1,n) - CurrentBestValue(m,n);
        edges_added = createStartingGraphFromInput(graph_m,inputgraph);

        if(edges_needed<=edges_added){
            if(edges_needed<edges_added)
                cout<<"M_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m+1,n);
            parent[m+1][n]= {m,n};
            graphs[m+1][n]=graph_m;
        }}

        edges_needed=CurrentBestValue(m,n+1) - CurrentBestValue(m,n);
        //inputgraph=create_from_file(m,n,readBestGraph(m,n));
        edges_added = createStartingGraphFromInput(graph_n,inputgraph);

        if(edges_needed<=edges_added){
            if(edges_needed<edges_added)
                cout<<"N_Nagy a baj!"<<m<<", "<<n<<endl;
            n_m_queue.emplace(m,n+1);
            parent[m][n+1]= {m,n};
            graphs[m][n+1]=graph_n;
        }
    }
    for(const auto& elm : done){
        for(const auto& elm2 : elm){
            cout<<elm2<<" ";
        }
        cout<<endl;
    }
    cout<<endl<<"13,16: "<<endl;
    int _n=16, _m=13;
    while(parent[_m][_n].first){
        cout<<"("<<_m<<","<<_n<<")"<<endl;
        auto par= parent[_m][_n];
        _m=par.first;
        _n=par.second;
    }
}

int main() {
    //backwards(); return 0;
    srand(time(0));
    /*for(int i=0;i<100;i++){
    int n=20 + std::rand() % (20 + 1);
    int m=20 + std::rand() % (20 + 1);
    if(m>n) swap(n, m);*/
    forward(3,3);
   // }

    return 0;
}