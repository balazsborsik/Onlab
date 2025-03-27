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
#include <stack>
using namespace std;

struct commmon_neighbours{
    int siz;
    vector<vector<int>> neighbours;
    //TODO!!! számontartani a párokat, ahol cmn>=2;
    commmon_neighbours(int siz_)
        : siz(siz_) { neighbours.assign(siz-1, vector<int>(siz, 0)); }

    int get_cmn_num(int a_, int b_){
        int a = min(a_,b_);
        int b = max(a_,b_);
        if(a>siz-1){
            throw invalid_argument("in get_cmn_num(int a_, int b_), int a = min(a_,b_)>siz-1");
        }
        return neighbours[a][b];
    }

    void add_cmn(int a_, int b_){
        int a = min(a_,b_);
        int b = max(a_,b_);
        if(a>siz-1){
            throw invalid_argument("in add_cmn(int a_, int b_), int a = min(a_,b_)>siz-1");
        }
        neighbours[a][b]++;
    }

    void remove_cmn(int a_, int b_){
        int a = min(a_,b_);
        int b = max(a_,b_);
        if(a>siz-1){
            throw invalid_argument("in add_cmn(int a_, int b_), int a = min(a_,b_)>siz-1");
        }
        neighbours[a][b]--;
    }
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
void storeCMN(commmon_neighbours& cmn, const vector<vector<int>>& adj, int m, int n, int u, int v) {// Its important, that adj[u][i]==0, when calling storeCMN
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            cmn.add_cmn(u,u2);
        }
    }
}

void removeFromCMN(commmon_neighbours& cmn, const vector<vector<int>>& adj, int m, int n, int u, int v) {// Its important, that adj[u][i]==0, when calling removeCMN
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            cmn.remove_cmn(u,u2);
        }
    }
}

void reflip_pair(commmon_neighbours& cmn, const pair<int, int>& pair, vector<vector<int>> &adj, double p, int m, int n){
    int u=pair.first;
    for(int i=0;i<n;i++){
        if(p>((double)rand()/(double)RAND_MAX)){
            if(!adj[u][i]){
                storeCMN(cmn,adj,m,n,u,i);  // Its important, that adj[u][i]==0, when calling storeCMN
                adj[u][i]=1;
            }
        }else{
            if(adj[u][i]){
                adj[u][i]=0;
                removeFromCMN(cmn,adj,m,n,u,i); // Its important, that adj[u][i]==0, when calling removeCMN
            }
        }
    }
    u=pair.second;
    for(int i=0;i<n;i++){
        if(p>((double)rand()/(double)RAND_MAX)){
            if(!adj[u][i]){
                storeCMN(cmn,adj,m,n,u,i);  // Its important, that adj[u][i]==0, when calling storeCMN
                adj[u][i]=1;
            }
        }else{
            if(adj[u][i]){
                adj[u][i]=0;
                removeFromCMN(cmn,adj,m,n,u,i); // Its important, that adj[u][i]==0, when calling removeCMN
            }
        }
    }
}

void run_with_p(vector<vector<int>>& adj, double p, int m, int n){
    commmon_neighbours cmn(m);
    for (int u = 0; u < m; ++u){
        for (int v = 0; v < n; ++v){
            if(p>((double)rand()/(double)RAND_MAX)){
                storeCMN(cmn,adj,m,n,u,v);
                adj[u][v]=1;
            }
        }
    }
    bool hasChanged = true;
    while(hasChanged){
        hasChanged = false;
        for(int i=0;i<cmn.siz-1;i++){
            for(int j=i;j<cmn.siz;j++){
                if(cmn.neighbours[i][j]>=2){
                    reflip_pair(cmn, make_pair(i, j), adj, p, m, n);
                    hasChanged = true;
                }
            }
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

int main() {
    srand(time(0));

    int m, n;
    /*cout << "Enter number of vertices on side U (m): ";
    cin >> m;
    cout << "Enter number of vertices on side V (n): ";
    cin >> n;*/
    
    vector<vector<int>> results(29, vector<int>(29,0));
    vector<pair<int,int>> n_mqueue;
    for(int firstcord=2;firstcord<=30;firstcord++){
        for (int secondcord = firstcord; secondcord <= 30; secondcord++)
        {
            n_mqueue.emplace_back(firstcord,secondcord);
        }
    }
    auto start = chrono::steady_clock::now();
    for(int iters=0;iters<n_mqueue.size();iters++){
        int m=n_mqueue[iters].first;
        int n=n_mqueue[iters].second;
        cout<<m<<", "<<n<<endl;
        int iterations = 50000;  // Number of trials
        int maxEdges = 0;       // Best lower bound found
        logs stats;
        stats.startTimer();
        vector<vector<int>> graph(m, vector<int>(n, 0));
        double p=((double)upperBound(2,2,n,m)/(n*m))*0.8;
        for (int iter = 0; iter < iterations; ++iter) {
            vector<vector<int>> adj(m, vector<int>(n, 0));
            vector<pair<int, int>> edges;

            run_with_p(adj, p, m, n);

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
        logfile.open("REAL_MT_log.txt", std::ios_base::app);
        logfile<<"Z(" << m << ", " << n << "; 2, 2): ";
        stats.print(logfile);
        logfile<<endl;
        logfile<<" p: "<<p<<";"<<endl;
        logfile.close();

        stringstream str;
        str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_"<<maxEdges<<".txt";
        ofstream outfile (str.str());
        print_graph( graph, outfile);
        outfile.close();
        results[n-2][m-2]=maxEdges;
    }

    stringstream str;
    str<<"REAL_MT_results.txt";
    ofstream resfile (str.str());
    print_graph( results, resfile);
    resfile.close();
    print_graph(results);

    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    ofstream logfile;
    logfile.open("exec_time.txt", std::ios_base::app);
    logfile<<"REAL_MT: "<<duration<<"seconds\n";
    logfile.close();

    return 0;
}