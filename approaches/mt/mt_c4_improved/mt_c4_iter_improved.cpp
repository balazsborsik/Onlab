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
using namespace std;

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

double calculate_p_from_LLL(int n, int m){
    double d=nCr(m,2)*nCr(n,2)-(nCr(m-2,2)*nCr(n,2)+2*(m-2)*nCr(n-2,2)+nCr(n-2,2))-1;
    return(pow(1.0/(4*d),1.0/4));
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

void reflip_circle(vector<c4>& new_circles, const vector<vector<int>> &edges_in_circles, const c4 circle, vector<vector<int>> &adj, int m, int n) {
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
}


void run_with_p(vector<vector<int>>& adj, double p, int m, int n){
    vector<c4> circles;
    vector<vector<int>> edges_in_circles(m, vector<int>(n, 0));
    for(int itr=0; itr<20;itr++){
        for (int u = 0; u < m; ++u){
            for (int v = 0; v < n; ++v){
                if(p>((double)rand()/(double)RAND_MAX)){
                    adj[u][v]=0;
                    storeK22(circles,edges_in_circles,adj,m,n,u,v);
                    adj[u][v]=1;
                }
            }
        }
        while(!circles.empty()){
            reflip_circle(circles,edges_in_circles, circles[0], adj, m, n);
            reeval_circles(circles,edges_in_circles, adj);
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
    
    vector<vector<int>> results(39, vector<int>(39,0));
    vector<pair<int,int>> n_mqueue;
    for(int firstcord=2;firstcord<=40;firstcord++){
        for (int secondcord = firstcord; secondcord <= 40; secondcord++)
        {
            n_mqueue.emplace_back(firstcord,secondcord);
        }
    }
    auto start = chrono::steady_clock::now();
    for(int iters=0;iters<n_mqueue.size();iters++){
        int m=n_mqueue[iters].first;
        int n=n_mqueue[iters].second;
        cout<<m<<", "<<n<<endl;
        int iterations = 120;  // Number of trials
        int maxEdges = 0;       // Best lower bound found
        logs stats;
        stats.startTimer();
        vector<vector<int>> graph(m, vector<int>(n, 0));

        vector<pair<int, int>> edges;
        for (int u = 0; u < m; ++u)
                for (int v = 0; v < n; ++v)
                    edges.emplace_back(u, v);
        double p=((double)upperBound(2,2,n,m)/(n*m))*0.85;
        //double p=calculate_p_from_LLL(n,m);
        for (int iter = 0; iter < iterations; ++iter) {
            vector<vector<int>> adj(m, vector<int>(n, 0));

            run_with_p(adj, p, m, n);

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
        logfile.open("MT_ITER_IMPROVED_log.txt", std::ios_base::app);
        logfile<<"Z(" << m << ", " << n << "; 2, 2): ";
        stats.print(logfile);
        logfile<<endl;
        logfile<<" p: "<<p<<";"<<endl;
        logfile<<" iter: "<<20<<";"<<endl;
        logfile.close();

        stringstream str;
        str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_"<<maxEdges<<".txt";
        ofstream outfile (str.str());
        print_graph( graph, outfile);
        outfile.close();
        results[n-2][m-2]=maxEdges;
    }

    stringstream str;
    str<<"MT_ITER_IMPROVED_results.txt";
    ofstream resfile (str.str());
    print_graph( results, resfile);
    resfile.close();
    print_graph(results);

    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    ofstream logfile;
    logfile.open("exec_time.txt", std::ios_base::app);
    logfile<<"MT_ITER_IMPROVED: "<<duration<<"seconds\n";
    logfile.close();

    return 0;
}