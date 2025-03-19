#include <iostream>
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
void storeK22(vector<c4>& circles, const vector<vector<int>>& adj, int m, int n, int u, int v) {
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            // Now, look for a second common neighbor
            for (int v2 = 0; v2 < n; ++v2) {
                if (v2 == v) continue;
                if (adj[u][v2] && adj[u2][v2]) {
                    circles.emplace_back(u, u2, v, v2); // Found a K_{2,2}
                }
            }
        }
    }
}

void reeval_circles(vector<c4> &circles, const vector<vector<int>> &graph){
    circles.erase(std::remove_if(circles.begin(), circles.end(), [&graph](const c4& circle){
        if(graph[circle.u][circle.v]!=1 || graph[circle.u][circle.v2]!=1 || graph[circle.u2][circle.v]!=1 || graph[circle.u2][circle.v2]!=1){
            return true;
        }
        return false;
    }), circles.end());
}

void reflip_circle(vector<c4>& new_circles, const c4& circle, vector<vector<int>> &adj, double p, int m, int n){
    adj[circle.u][circle.v]=0;
    adj[circle.u][circle.v2]=0;
    adj[circle.u2][circle.v]=0;
    adj[circle.u2][circle.v2]=0;
    if(p>((double)rand()/(double)RAND_MAX)){
        storeK22(new_circles,adj,m,n,circle.u,circle.v);
        adj[circle.u][circle.v]=1;
    }
    if(p>((double)rand()/(double)RAND_MAX)){
        storeK22(new_circles,adj,m,n,circle.u,circle.v2);
        adj[circle.u][circle.v2]=1;
    }
    if(p>((double)rand()/(double)RAND_MAX)){
        storeK22(new_circles,adj,m,n,circle.u2,circle.v);
        adj[circle.u2][circle.v]=1;
    }
    if(p>((double)rand()/(double)RAND_MAX)){
        storeK22(new_circles,adj,m,n,circle.u2,circle.v2);
        adj[circle.u2][circle.v2]=1;
    }
}

void run_with_p(vector<vector<int>>& adj, double p, int m, int n){
    vector<c4> circles;
    for (int u = 0; u < m; ++u){
        for (int v = 0; v < n; ++v){
            if(p>((double)rand()/(double)RAND_MAX)){
                storeK22(circles,adj,m,n,u,v);
                adj[u][v]=1;
            }
        }
    }
    while(!circles.empty()){
        vector<c4> new_circles;
        for(const auto& elm : circles){
            reflip_circle(new_circles, elm, adj, p, m, n);
        }
        circles.insert(circles.end(), new_circles.begin(), new_circles.end());
        reeval_circles(circles, adj);
    }
}

int main() {
    srand(time(0));

    int m, n;
    cout << "Enter number of vertices on side U (m): ";
    cin >> m;
    cout << "Enter number of vertices on side V (n): ";
    cin >> n;

    int iterations = 100000;  // Number of trials
    int maxEdges = 0;       // Best lower bound found

    auto start = chrono::steady_clock::now();

    for (int iter = 0; iter < iterations; ++iter) {
        vector<vector<int>> adj(m, vector<int>(n, 0));
        vector<pair<int, int>> edges;

        run_with_p(adj, 0.18, m, n);

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

        if (edgeCount > maxEdges) {
            maxEdges = edgeCount;
            cout << "New lower bound found: " << maxEdges << " edges (iteration " << iter << ")" << endl;
        }
    }

    auto end = chrono::steady_clock::now();

    cout << "\nEstimated lower bound for Z(" << m << ", " << n << "; 2, 2): " << maxEdges << endl;

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    return 0;
}
