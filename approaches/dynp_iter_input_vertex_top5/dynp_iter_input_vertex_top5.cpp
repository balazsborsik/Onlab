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

vector<string> read5BestGraphs(int m, int n) {
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
    int maxNum = -1;
    vector<pair<int, string>> filenames;
    filenames.emplace_back(0, filestart.str()+"0.txt"); //default value so if no graphs are stored it starts form somewhere

    while ((entry = readdir(dir)) != nullptr) {
        string filename = entry->d_name;

        smatch match;
        if (std::regex_match(filename, match, pattern)) {
            int num = std::stoi(match[1]);
            filenames.emplace_back(num, "output\\" + filename);
        }
    }
    sort(filenames.begin(), filenames.end());
    closedir(dir);
    vector<string> results;
    int siz=filenames.size();
    if(siz>=5){
        for(int i=0;i<5;i++){
            results.push_back(filenames[siz-1-i].second);
        }
    }else{
        for(int i=0;i<siz;i++){
            results.push_back(filenames[siz-1-i].second);
        }
    }
    return results;
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
        expected_n_percent=(double)upper_bound/(n*m);
        edgenum_n.assign(n,0);
        edgenum_m.assign(m,0);
    }

    dynp(const vector<vector<int>> &graph, int upper_bound)
    {   
        int m=graph.size();
        int n=graph[0].size();
        expected_m=((double)upper_bound/m+0.5);
        expected_n=((double)upper_bound/n+0.5);
        expected_n_percent=(double)upper_bound/(n*m);
        edgenum_n.assign(n,0);
        edgenum_m.assign(m,0);
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                if(graph[i][j]){
                    edgenum_m[i]++;
                    edgenum_n[j]++;
                }
            }
        }
    }

    double get_multiplier(int v_m){
        if(edgenum_m[v_m]>=expected_m){
            if(edgenum_m[v_m]==expected_m)
                return 0.9;
            return 0.5;
        }
        return 6.7-5.6*((double)edgenum_m[v_m]/expected_m);
    }

    double get_p(int v_m, int v_n){
        int n=edgenum_n[v_n];
        return min(get_multiplier(v_m)*(expected_n_percent*1.1-(n/expected_n)*expected_n_percent),0.8);
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




void run_with_p(vector<vector<int>>& adj, dynp &dyn_p, int insideIterations, int m, int n){
    vector<c4> circles;
    vector<vector<int>> edges_in_circles(m, vector<int>(n, 0));
    for(int itr=0; itr<insideIterations;itr++){
        for (int u = 0; u < m; ++u){
            for (int v = 0; v < n; ++v){
                if(dyn_p.get_p(u,v)>((double)rand()/(double)RAND_MAX)){
                    if(adj[u][v]==0){
                        dyn_p.add_edge(u,v);
                    }
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

void addVertexToM(vector<vector<int>> &adj, const vector<vector<int>> &inputgraph){
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
}

void addVertexToN(vector<vector<int>> &adj, const vector<vector<int>> &inputgraph){
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
}

void createStartingGraphFromInput(vector<vector<int>> &adj, const vector<vector<int>> &inputgraph){
    int m=adj.size();
    int n=adj[0].size();
    if(m-1 == inputgraph.size()&&n == inputgraph[0].size())
        addVertexToM(adj, inputgraph);
    else if(m == inputgraph.size()&&n-1 == inputgraph[0].size())
        addVertexToN(adj, inputgraph);
    else __throw_invalid_argument("az input és az adj gráfok nem megfelelő méretűek");
}

bool top5contains(int value, const vector<pair<int,vector<vector<int>>>> &top5){
    for(const auto &elm: top5){
        if(value==elm.first)
            return true;
    }
    return false;
}

vector<vector<int>> run_in_range(int min, int max, bool appendToM, int runid){
    int m, n;
    vector<vector<int>> results(39, vector<int>(39,0));
    vector<pair<int,int>> n_mqueue;
    for(int firstcord=min;firstcord<=max;firstcord++){
        for (int secondcord = firstcord; secondcord <= max; secondcord++)
        {
            n_mqueue.emplace_back(firstcord,secondcord);
        }
    }

    //n_mqueue.emplace_back(18,28);
    auto start = chrono::steady_clock::now();
    for(int iters=0;iters<n_mqueue.size();iters++){
        int m=n_mqueue[iters].first;
        int n=n_mqueue[iters].second;
        //cout<<m<<", "<<n<<endl;
        int iterations = 24;  // Number of trials
        int insideIterations[] = {1,3,6,12}; // 1, 3, 6, 12
        logs stats;
        stats.startTimer();
        
        vector<pair<int, int>> edges;

        for (int u = 0; u < m; ++u)
                for (int v = 0; v < n; ++v)
                    edges.emplace_back(u, v);
        vector<pair<int,vector<vector<int>>>> top5(5, {0, vector<vector<int>>(m, vector<int>(n, 0))});
        //vector<vector<int>> graph(m, vector<int>(n, 0));

        int input_m;
        int input_n;
        
        if(appendToM || m==n){
            input_m=m-1;
            input_n=n;
        }else{
            input_m=m;
            input_n=n-1;
        }
        vector<string> filenames=read5BestGraphs(input_m, input_n);
        for(const auto& inputfilename : filenames){
            vector<vector<int>> inputgraph=create_from_file(input_m, input_n, inputfilename);
            
            ///double p=((double)upperBound(2,2,n,m)/(n*m))*0.85;
            for (int iter = 0; iter < iterations; ++iter) {
                vector<vector<int>> adj(m, vector<int>(n, 0));
                
                createStartingGraphFromInput(adj, inputgraph);

                //dynp dyn_p(n, m, upperBound(2,2,n,m));
                dynp dyn_p(adj, upperBound(2,2,n,m));
                run_with_p(adj, dyn_p, insideIterations[iters%4], m, n);

                // Generate all possible edges
                

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
                if (edgeCount > top5.back().first && !top5contains(edgeCount, top5)) {
                    int spot;
                    for(int i=0;i<5;i++){
                        if(top5[i].first<edgeCount){
                            spot=i;
                            break;
                        }
                    }
                    top5.pop_back();
                    top5.insert(top5.begin()+spot, {edgeCount, adj});
                }
            }
        }

        stringstream logname;
        if(appendToM)
            logname<<"DYNP_ITER_INPUT_VERTEX_TOP5_M"<<runid<<"_log.txt";
        else
            logname<<"DYNP_ITER_INPUT_VERTEX_TOP5_N"<<runid<<"_log.txt";
        ofstream logfile;
        logfile.open(logname.str(), std::ios_base::app);
        logfile<<"Z(" << m << ", " << n << "; 2, 2): ";
        stats.print(logfile);
        //logfile<<endl;
        ///logfile<<" p: "<<p<<";"<<endl;
        logfile<<" iter: "<<insideIterations<<";"<<endl;
        logfile.close();
        
        //cout<<maxEdges<<endl;
        for(int i=0;i<5;i++){
        stringstream str;
        str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_"<<top5[i].first<<".txt";
        ofstream outfile (str.str());
        print_graph( top5[i].second, outfile);
        outfile.close();
        
        }
        results[n-2][m-2]=top5[0].first;
    }
    
    stringstream str;
     if(appendToM)
            str<<"DYNP_ITER_INPUT_VERTEX_TOP5_M"<<runid<<"_result.txt";
        else
            str<<"DYNP_ITER_INPUT_VERTEX_TOP5_N"<<runid<<"_results.txt";
    ofstream resfile (str.str());
    print_graph( results, resfile);
    resfile.close();
    //print_graph(results);

    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    ofstream logfile;
    logfile.open("exec_time.txt", std::ios_base::app);
    if(appendToM)
        logfile<<"DYNP_ITER_INPUT_VERTEX_TOP5_M"<<runid<<": "<<duration<<"seconds\n";
    else
        logfile<<"DYNP_ITER_INPUT_VERTEX_TOP5_N"<<runid<<": "<<duration<<"seconds\n";
    logfile.close();

    return results;
}

int main() {
    srand(time(0));

    vector<vector<int>> results(39, vector<int>(39,0));
    auto start = chrono::steady_clock::now();
    for(int i=0;i<10;i++){
        vector<vector<int>> temp=run_in_range(2,40,i%2,i/2+1);
        for(int j=0;j<results.size();j++){
            for(int k=0;k<results.size();k++){
                if(temp[j][k]>results[j][k]){
                    results[j][k]=temp[j][k];
                }
            }
        }
    }
    
    stringstream str;
    str<<"DYNP_ITER_INPUT_VERTEX_TOP5_results.txt";
    ofstream resfile (str.str());
    print_graph( results, resfile);
    resfile.close();
    print_graph(results);

    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    ofstream logfile;
    logfile.open("exec_time.txt", std::ios_base::app);
    logfile<<"DYNP_ITER_INPUT_VERTEX_TOP5: "<<duration<<"seconds\n";
    logfile.close();

    return 0;
}