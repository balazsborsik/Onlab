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
#include <windows.h>
using namespace std;

std::string readBestGraph(int m, int n) {
    WIN32_FIND_DATAA findData;
    stringstream filestart;
    filestart<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_";
    HANDLE hFind = FindFirstFileA(("output\\"+filestart.str()+"*.txt").c_str(), &findData);
    if (hFind == INVALID_HANDLE_VALUE) {
        std::cerr << "No matching files found.\n";
        return "";
    }

    std::regex pattern(filestart.str() + R"((\d+)\.txt)");
    int maxNum = -1;
    std::string result;

    do {
        std::string filename = findData.cFileName;
        std::smatch match;
        if (std::regex_match(filename, match, pattern)) {
            int num = std::stoi(match[1]);
            if (num > maxNum) {
                maxNum = num;
                result = "output\\" + filename;
            }
        }
    } while (FindNextFileA(hFind, &findData));

    FindClose(hFind);
    return result;
}

struct dynp{
    vector<int> edgenum_m;
    vector<int> edgenum_n;
    int expected_m;
    int expected_n;
    double expected_n_percent;

    dynp(int n, int m,  int upper_bound)
    {   
        expected_m=((double)upper_bound/m+0.5);
        expected_n=((double)upper_bound/n+0.5);
        cout<<n<<" , "<<m<<upper_bound<<expected_m<<" "<<expected_n<<endl;
        expected_n_percent=(double)upper_bound/(n*m);
        edgenum_n.assign(n,0);
        edgenum_m.assign(m,2.0);
    }

    dynp(const vector<vector<int>> &graph, int upper_bound)
    {   
        int m=graph.size();
        int n=graph[0].size();
        expected_m=((double)upper_bound/m+0.5); // zárójel nélkül nincs értelme
        expected_n=((double)upper_bound/n+0.5);
        expected_n_percent=(double)upper_bound/(n*m);
        edgenum_n.assign(n,0);
        edgenum_m.assign(m,2.0); //!! nem így volt hanem valamiért 2 volt
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                if(graph[i][j]){    // ha van él
                    edgenum_m[i]++;
                    edgenum_n[j]++;
                }
            }
        }
    }

    double get_multiplier(int v_m){
        if(edgenum_m[v_m]>=expected_m){
            return 0.15;
        }
        return 1.4-1.1*((double)edgenum_m[v_m]/expected_m);
    }

    double get_p(int v_m, int v_n){
        int n=edgenum_n[v_n];
        return min(4*get_multiplier(v_m)*(expected_n_percent*1.1-(n/expected_n)*expected_n_percent),0.8);
    }

    void delete_edge(int v_m, int v_n){
        edgenum_m[v_m]--;
        edgenum_n[v_n]--;
    }

    void add_edge(int v_m, int v_n){
        edgenum_m[v_m]++;
        edgenum_n[v_n]++;
    }

};

struct c4{
    int u;
    int u2;
    int v;
    int v2;

    c4(int u_, int u2_, int v_, int v2_)
        : u(u_), u2(u2_), v(v_), v2(v2_) {}
};

struct logs{
    int siz = 0;          // Count of numbers
    int maximum = 0;
    double mean = 0.0;  // Running mean
    double M2 = 0.0;    // Sum of squares of differences from the mean
    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    array<int, 5> top5 = {0,0,0,0,0};

    void startTimer(){
        start=std::chrono::steady_clock::now();
    }

    void add(int x) {
        siz += 1;
        double delta = x - mean;
        mean += delta / siz;
        double delta2 = x - mean;
        M2 += delta * delta2;
        if(x>maximum){
            for(int i=x-maximum;i<top5.size();i++){
                top5[i-(x-maximum)]=top5[i];
            }
            for(int i=max((int)top5.size()-(x-maximum),0);i<top5.size();i++){
                top5[i]=0;
            }
            maximum=x;
        }
        if(maximum-x<top5.size()){
            top5[top5.size()-1-(maximum-x)]++;
        }
    }

    double variance() const {
        return (siz > 0) ? M2 / (siz) : 0.0;  // Sample variance
    }

    void print(ostream& out=cout) const {
        out << "Max: " << maximum << "; ";
        out << "Mean: " << mean << "; ";
        out << "Variance: " << variance() << "; ";
        out << "Size: " << siz << "; ";
        out <<"Time: "<< std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-start).count() << " ms" << "; ";
        out << "Top5_Histogram_(relative_to_max): ";
        for (auto count : top5) {
            out << count << " ";
        }
        out <<";";
    }

};

long nCr(int n,int r)
{
    long ans=1;
    r=r>n-r?n-r:r;
    int j=1;
    for(;j<=r;j++,n--)
    {
        if(n%j==0)
        {
            ans*=n/j;
        }else
        if(ans%j==0)
        {
            ans=ans/j*n;
        }else
        {
            ans=(ans*n)/j;
        }
    }
    return ans;
}


int upperBound(int s, int t, int m, int n){
    int p=s-1;
    if(m>n){
        swap(s, t);
        swap(m, n);
    }
    if(m<s)
        return m*s;
    double previous=numeric_limits<double>::max();
    double result=numeric_limits<double>::max()- 1e307;
    while(previous>result){
        previous=result;
        result =
        ( (double)(t - 1) / nCr(p, s - 1) ) * nCr(m, s)
        + n * ((p + 1) * (s - 1)) / (double)s;
        p++;
    }
    return (int)previous;
}

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

// Function to check if adding edge (u, v) creates a K_{2,2}
void storeK22(vector<c4>& circles, vector<vector<int>>& edges_in_circles, const vector<vector<int>>& adj, int m, int n, int u, int v) {
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            // Now, look for a second common neighbor
            for (int v2 = 0; v2 < n; ++v2) {
                if (v2 == v) continue;
                if (adj[u][v2] && adj[u2][v2]) {
                    edges_in_circles[u][v]++;
                    edges_in_circles[u][v2]++;
                    edges_in_circles[u2][v]++;
                    edges_in_circles[u2][v2]++;
                    circles.emplace_back(u, u2, v, v2); // Found a K_{2,2}
                    
                }
            }
        }
    }
}

void reeval_circles(vector<c4> &circles,vector<vector<int>> &edges_in_circles, const vector<vector<int>> &graph){
    circles.erase(std::remove_if(circles.begin(), circles.end(), [&graph, &edges_in_circles](const c4& circle){
        if(graph[circle.u][circle.v]!=1 || graph[circle.u][circle.v2]!=1 || graph[circle.u2][circle.v]!=1 || graph[circle.u2][circle.v2]!=1){
            edges_in_circles[circle.u][circle.v]--;
            edges_in_circles[circle.u][circle.v2]--;
            edges_in_circles[circle.u2][circle.v]--;
            edges_in_circles[circle.u2][circle.v2]--;
            return true;
        }
        return false;
    }), circles.end());
}

void reflip_circle(vector<c4>& new_circles, const vector<vector<int>> &edges_in_circles, const c4 circle, vector<vector<int>> &adj, dynp &dyn_p, int m, int n) {
    pair<int, int> edge = make_pair(circle.u, circle.v);
    int fok = 0;

    for (const auto& elm : new_circles) {
        if (edges_in_circles[elm.u][elm.v] > fok || 
            (edges_in_circles[elm.u][elm.v] == fok && make_pair(elm.u, elm.v) < edge)) {
            fok = edges_in_circles[elm.u][elm.v];
            edge = make_pair(elm.u, elm.v);
        }
        if (edges_in_circles[elm.u][elm.v2] > fok || 
            (edges_in_circles[elm.u][elm.v2] == fok && make_pair(elm.u, elm.v2) < edge)) {
            fok = edges_in_circles[elm.u][elm.v2];
            edge = make_pair(elm.u, elm.v2);
        }
        if (edges_in_circles[elm.u2][elm.v] > fok || 
            (edges_in_circles[elm.u2][elm.v] == fok && make_pair(elm.u2, elm.v) < edge)) {
            fok = edges_in_circles[elm.u2][elm.v];
            edge = make_pair(elm.u2, elm.v);
        }
        if (edges_in_circles[elm.u2][elm.v2] > fok || 
            (edges_in_circles[elm.u2][elm.v2] == fok && make_pair(elm.u2, elm.v2) < edge)) {
            fok = edges_in_circles[elm.u2][elm.v2];
            edge = make_pair(elm.u2, elm.v2);
        }
    }
    
    adj[edge.first][edge.second] = 0;
    dyn_p.delete_edge(edge.first, edge.second);
}

void run_with_p(vector<vector<int>>& adj, dynp &dyn_p, int m, int n){
    vector<c4> circles;
    vector<vector<int>> edges_in_circles(m, vector<int>(n, 0));
    for(int itr=0; itr<200;itr++){
        for (int u = 0; u < m; ++u){
            for (int v = 0; v < n; ++v){
                if(dyn_p.get_p(u,v)>((double)rand()/(double)RAND_MAX)){
                    if(adj[u][v]==0)
                        dyn_p.add_edge(u,v);
                    adj[u][v]=0;
                    storeK22(circles,edges_in_circles,adj,m,n,u,v);
                    
                    adj[u][v]=1;
                }
            }
        }
        while(!circles.empty()){
            reflip_circle(circles, edges_in_circles, circles[0], adj, dyn_p, m, n);
            reeval_circles(circles, edges_in_circles, adj);
        }
    }
}

void print_graph(const vector<vector<int>>& graph, ostream& out=cout){
    for(const auto& vertex : graph){
        for (const auto& edge : vertex)
        {
           out<<edge<<' ';
        }
        out<<endl;
    }
}

vector<vector<int>> create_from_file(int m, int n, string filename){
    vector<vector<int>> graph(m, vector<int>(n, 0));
    std::ifstream file(filename);
    for(int i=0;i<m;i++){
        for(int j=0;j<n;j++){
            int edge;
            file>>edge;
            graph[i][j]=edge;
        }
    }
    file.close();
    return graph;
}

int main() {
    srand(time(0));

    int m, n;
    /*cout << "Enter number of vertices on side U (m): ";
    cin >> m;
    cout << "Enter number of vertices on side V (n): ";
    cin >> n;*/
    vector<vector<int>> results(29, vector<int>(29,0));
    vector<pair<int,int>> n_mqueue;
    /*for(int firstcord=2;firstcord<=30;firstcord++){
        for (int secondcord = firstcord; secondcord <= 30; secondcord++)
        {
            n_mqueue.emplace_back(firstcord,secondcord);
        }
    }*/
    //n_mqueue.emplace_back(14,29);
    //n_mqueue.emplace_back(29,14);
    //n_mqueue.emplace_back(18,28);
    //n_mqueue.emplace_back(19,27);
    n_mqueue.emplace_back(22,30);
    //n_mqueue.emplace_back(21,25);
    auto start = chrono::steady_clock::now();
    for(int iters=0;iters<n_mqueue.size();iters++){
        int m=n_mqueue[iters].first;
        int n=n_mqueue[iters].second;
        cout<<m<<", "<<n<<endl;
        int iterations = 20000;  // Number of trials
        int maxEdges = 0;       // Best lower bound found
        logs stats;
        stats.startTimer();
        vector<vector<int>> graph(m, vector<int>(n, 0));
        vector<vector<int>> inputgraph=create_from_file(21,25,readBestGraph(21,25));

        ///double p=((double)upperBound(2,2,n,m)/(n*m))*0.85;
        for (int iter = 0; iter < iterations; ++iter) {
            vector<vector<int>> adj(m, vector<int>(n, 0));
            for(int i=0;i<21;i++){
                for(int j=0;j<25;j++){
                    adj[i][j]=inputgraph[i][j];
                }
            }
            for(int j=25;j<n;j++){
                adj[21][j]=1;
            }
            vector<pair<int, int>> edges;
            //dynp dyn_p(n, m, upperBound(2,2,n,m));
            dynp dyn_p(adj, upperBound(2,2,n,m));
            run_with_p(adj, dyn_p, m, n);

            // Generate all possible edges
            for (int u = 0; u < m; ++u)
                for (int v = 0; v < n; ++v)
                    edges.emplace_back(u, v);

            // Shuffle edge order randomly
            random_shuffle(edges.begin(), edges.end());

            int edgeCount = 0;

            // Try to add edges one by one, only if they don't create a K_{2,2}
            for (auto& e : edges) {
                int u = e.first;
                int v = e.second;
                if(adj[u][v]==1){
                    ++edgeCount;
                }else if (!createsK22(adj, m, n, u, v)) {
                    adj[u][v] = 1;
                    ++edgeCount;
                }
            }
            stats.add(edgeCount);
            if (edgeCount > maxEdges) {
                maxEdges = edgeCount;
                graph=adj;
            }
        }

        ofstream logfile;
        logfile.open("DYNP_ITER_INPUT_log.txt", std::ios_base::app);
        logfile<<"Z(" << m << ", " << n << "; 2, 2): ";
        stats.print(logfile);
        logfile<<endl;
        ///logfile<<" p: "<<p<<";"<<endl;
        logfile<<" iter: "<<200<<";"<<endl;
        logfile.close();

        stringstream str;
        str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_"<<maxEdges<<".txt";
        ofstream outfile (str.str());
        print_graph( graph, outfile);
        outfile.close();
        results[n-2][m-2]=maxEdges;

        cout<<maxEdges<<endl;
        cout<<readBestGraph(22,30)<<endl;
    }
    
    /*stringstream str;
    str<<"DYNP_ITER_results.txt";
    ofstream resfile (str.str());
    print_graph( results, resfile);
    resfile.close();
    print_graph(results);*/

    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    /*ofstream logfile;
    logfile.open("exec_time.txt", std::ios_base::app);
    logfile<<"DYNP: "<<duration<<"seconds\n";
    logfile.close();*/

    return 0;
}