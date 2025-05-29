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

#define VALUE_OF_S 5
#define VALUE_OF_T 5

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

vector<string> read2BestGraphs(int m, int n) {
    stringstream filestart;
    filestart << "Z" << m << "_" << n << "_" << VALUE_OF_S << "_" << VALUE_OF_T << "_";
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
    if(siz>=2){
        for(int i=0;i<2;i++){
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
    vector<int> degree_m;
    vector<int> degree_n;
    int expected_m;
    int expected_n;
    double expected_n_percent;

    dynp(int n, int m,  int upper_bound)
    {   
        expected_m=((double)upper_bound/m+0.5);
        expected_n=((double)upper_bound/n+0.5);
        expected_n_percent=(double)upper_bound/(n*m);
        degree_n.assign(n,0);
        degree_m.assign(m,0);
    }

    dynp(const vector<vector<int>> &graph, int upper_bound)
    {   
        int m=graph.size();
        int n=graph[0].size();
        expected_m=((double)upper_bound/m+0.5);
        expected_n=((double)upper_bound/n+0.5);
        expected_n_percent=(double)upper_bound/(n*m);
        degree_n.assign(n,0);
        degree_m.assign(m,0);
        for(int i=0;i<m;i++){
            for(int j=0;j<n;j++){
                if(graph[i][j]){
                    degree_m[i]++;
                    degree_n[j]++;
                }
            }
        }
    }

    double get_multiplier(int v_m){
        if(degree_m[v_m]>=expected_m){
            if(degree_m[v_m]==expected_m)
                return 0.9;
            return 0.5;
        }
        return 6.7-5.6*((double)degree_m[v_m]/expected_m);
    }

    double get_p(int v_m, int v_n){
        int n=degree_n[v_n];
        return min(get_multiplier(v_m)*(expected_n_percent*1.1-(n/(expected_n+0.0001))*expected_n_percent),0.8);
    }

    void delete_edge(int v_m, int v_n){
        degree_m[v_m]--;
        degree_n[v_n]--;
    }

    void add_edge(int v_m, int v_n){
        degree_m[v_m]++;
        degree_n[v_n]++;
    }

};

struct K55{
    int u[VALUE_OF_S];
    int v[VALUE_OF_T];

    K55(int u_, int u2_, int u3_, int u4_, int u5_, int v_, int v2_, int v3_, int v4_, int v5_){
        u[0]=u_; u[1]=u2_; u[2]=u3_, u[3]=u4_, u[4]=u5_;
        v[0]=v_; v[1]=v2_; v[2]=v3_, v[3]=v4_, v[4]=v5_;
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
bool createsK55(const vector<vector<int>>& adj, int m, int n, int u, int v) {
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            // Now, look for a second common neighbor
            for (int v2 = 0; v2 < n; ++v2) {
                //if (v2 == v) continue;
                if (adj[u][v2] && adj[u2][v2]) {

                    for(int u3 = u2 + 1; u3 < m; ++u3){

                        if(adj[u3][v] && adj[u3][v2]){
                            
                            for(int v3 = v2 + 1; v3 < n; ++v3){

                                if(adj[u][v3] && adj[u2][v3] && adj[u3][v3]){

                                    for(int u4 = u3 + 1; u4 < m; ++u4){

                                        if(adj[u4][v] && adj[u4][v2] && adj[u4][v3]){

                                            for(int v4 = v3 + 1; v4 < n; ++v4){

                                                if(adj[u][v4] && adj[u2][v4] && adj[u3][v4] && adj[u4][v4]){

                                                    for(int u5 = u4 + 1; u5 < m; ++u5){
                                                        
                                                        if(adj[u5][v] && adj[u5][v2] && adj[u5][v3] && adj[u5][v4]){

                                                            for(int v5 = v4 + 1; v5 < n; ++v5){
                                                                
                                                                if(adj[u][v5] && adj[u2][v5] && adj[u3][v5] && adj[u4][v5] && adj[u5][v5]){
                                                                    return true; // Found a K_{5,5}
                                                                }

                                                            }

                                                        }

                                                    }

                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }
  
                }
            }
        }
    }
    return false;
}

// Function to check if adding edge (u, v) creates a K_{2,2}
void storeK55(vector<K55>& circles, vector<vector<int>>& edges_in_circles, const vector<vector<int>>& adj, int m, int n, int u, int v) {
    // Check for another vertex u2 that shares neighbor v
    for (int u2 = 0; u2 < m; ++u2) {
        //if (u2 == u) continue;
        if (adj[u2][v]) {
            // Now, look for a second common neighbor
            for (int v2 = 0; v2 < n; ++v2) {
                //if (v2 == v) continue;
                if (adj[u][v2] && adj[u2][v2]) {

                    for(int u3 = u2 + 1; u3 < m; ++u3){

                        if(adj[u3][v] && adj[u3][v2]){
                            
                            for(int v3 = v2 + 1; v3 < n; ++v3){

                                if(adj[u][v3] && adj[u2][v3] && adj[u3][v3]){

                                    for(int u4 = u3 + 1; u4 < m; ++u4){

                                        if(adj[u4][v] && adj[u4][v2] && adj[u4][v3]){

                                            for(int v4 = v3 + 1; v4 < n; ++v4){

                                                if(adj[u][v4] && adj[u2][v4] && adj[u3][v4] && adj[u4][v4]){

                                                    for(int u5 = u4 + 1; u5 < m; ++u5){
                                                        
                                                        if(adj[u5][v] && adj[u5][v2] && adj[u5][v3] && adj[u5][v4]){

                                                            for(int v5 = v4 + 1; v5 < n; ++v5){
                                                                
                                                                if(adj[u][v5] && adj[u2][v5] && adj[u3][v5] && adj[u4][v5] && adj[u5][v5]){
                                                                    K55 created(u, u2, u3, u4, u5, v, v2, v3, v4, v5);
                                                                    for(int i=0; i < VALUE_OF_S; i++){
                                                                        for(int j=0; j < VALUE_OF_T; j++){
                                                                            edges_in_circles[created.u[i]][created.v[j]]++;
                                                                        }
                                                                    }
                                                                    circles.push_back(created); // Found a K_{5,5}
                                                                }

                                                            }

                                                        }

                                                    }

                                                    
                                                }

                                            }

                                        }

                                    }

                                }

                            }

                        }

                    }
  
                }

            }
        }
    }
}

void reeval_circles(vector<K55> &circles,vector<vector<int>> &edges_in_circles, const vector<vector<int>> &graph){
    circles.erase(std::remove_if(circles.begin(), circles.end(), [&graph, &edges_in_circles](const K55& circle){
        for(int i=0; i < VALUE_OF_S; i++){
            for(int j=0; j < VALUE_OF_T; j++){
                if(!graph[circle.u[i]][circle.v[j]]){
                    for(int i=0; i < VALUE_OF_S; i++)
                        for(int j=0; j < VALUE_OF_T; j++)
                            edges_in_circles[circle.u[i]][circle.v[j]]--;
                    return true;
                }
            }
        }
        return false;
    }), circles.end());
}

void reflip_circle(vector<K55>& new_circles, const vector<vector<int>> &edges_in_circles, const K55 circle, vector<vector<int>> &adj, dynp &dyn_p, int m, int n) {
    pair<int, int> edge = make_pair(circle.u[0], circle.v[0]);
    int fok = 0;
    for (const auto& elm : new_circles) {
        for(int i=0; i < VALUE_OF_S; i++){
            for(int j=0; j < VALUE_OF_T; j++){
                int u = elm.u[i];
                int v = elm.v[j];
                if(edges_in_circles[u][v] > fok || edges_in_circles[u][v] == fok && 
                (u < edge.first || (!(edge.first < u) && v < edge.second))) { 
                    fok = edges_in_circles[u][v];
                    edge.first=u;
                    edge.second=v;
                }
            }
        }
    }
    
    adj[edge.first][edge.second] = 0;
    dyn_p.delete_edge(edge.first, edge.second);
}


void run_with_p(vector<vector<int>>& adj, dynp &dyn_p, int insideIterations, int m, int n){
    vector<K55> circles;
    vector<vector<int>> edges_in_circles(m, vector<int>(n, 0));
    for(int itr=0; itr<insideIterations;itr++){
        for (int u = 0; u < m; ++u){
            for (int v = 0; v < n; ++v){
                if(dyn_p.get_p(u,v)>((double)rand()/(double)RAND_MAX)){
                    if(adj[u][v]==0){
                        dyn_p.add_edge(u,v);
                    }
                    adj[u][v]=0;
                    storeK55(circles,edges_in_circles,adj,m,n,u,v);
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
    int sqrt_n = pow(n,0.5);
    for(int k=0;k<m*sqrt_n/8+5;k++){
        random_shuffle(begin(temp), end(temp));
        vector<int> vertices;
        vector<int> neighbours(m, 0);
        for(int i=0;i<n;i++){
            int vertex=temp[i];
            if(!createsK55(adj,m,n,m-1,vertex)){
                vertices.push_back(vertex);
                adj[m-1][vertex]=1;
            }
        }
        if(vertices.size()>best_size){
            best_size=vertices.size();
            chosenvertices=vertices;
        }
        for(const auto& elm : vertices){
            adj[m-1][elm]=0;
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
    int sqrt_n = pow(n,0.5);
    for(int k=0;k<m*sqrt_n/4+5;k++){
        random_shuffle(begin(temp), end(temp));
        vector<int> vertices;
        for(int i=0;i<m;i++){
                int vertex=temp[i];
                if(!createsK55(adj,m,n,vertex,n-1)){
                    vertices.push_back(vertex);
                    adj[vertex][n-1]=1;
                }
            }
            if(vertices.size()>best_size){
                best_size=vertices.size();
                chosenvertices=vertices;
            }
            for(const auto& elm : vertices){
                adj[elm][n-1]=0;
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

bool top2contains(int value, const vector<pair<int,vector<vector<int>>>> &top2){
    for(const auto &elm: top2){
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
        int iterations = 2;  // Number of trials
        int insideIterations[] = {1,3}; // 1, 3, 6, 12
        logs stats;
        stats.startTimer();
        
        vector<pair<int, int>> edges;

        for (int u = 0; u < m; ++u)
                for (int v = 0; v < n; ++v)
                    edges.emplace_back(u, v);
        vector<pair<int,vector<vector<int>>>> top2(2, {0, vector<vector<int>>(m, vector<int>(n, 0))});
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
        vector<string> filenames=read2BestGraphs(input_m, input_n);
        for(const auto& inputfilename : filenames){
            vector<vector<int>> inputgraph=create_from_file(input_m, input_n, inputfilename);
            
            ///double p=((double)upperBound(2,2,n,m)/(n*m))*0.85;
            for (int iter = 0; iter < iterations; ++iter) {
                vector<vector<int>> adj(m, vector<int>(n, 0));
                
                createStartingGraphFromInput(adj, inputgraph);

                //dynp dyn_p(n, m, upperBound(2,2,n,m));
                dynp dyn_p(adj, upperBound(VALUE_OF_S,VALUE_OF_T,n,m));
                run_with_p(adj, dyn_p, insideIterations[iters%2], m, n);

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
                    }else if (!createsK55(adj, m, n, u, v)) {
                        adj[u][v] = 1;
                        ++edgeCount;
                    }
                }
                stats.add(edgeCount);
                if (edgeCount > top2.back().first && !top2contains(edgeCount, top2)) {
                    int spot;
                    for(int i=0;i<2;i++){
                        if(top2[i].first<edgeCount){
                            spot=i;
                            break;
                        }
                    }
                    top2.pop_back();
                    top2.insert(top2.begin()+spot, {edgeCount, adj});
                }
            }
        }

        stringstream logname;
        if(appendToM)
            logname<<"DYNP_ITER_INPUT_VERTEX_TOP2_M"<<runid<<"_log.txt";
        else
            logname<<"DYNP_ITER_INPUT_VERTEX_TOP2_N"<<runid<<"_log.txt";
        ofstream logfile;
        logfile.open(logname.str(), std::ios_base::app);
        logfile<<"Z(" << m << ", " << n << "; "<<VALUE_OF_S<<", "<<VALUE_OF_T<<"): ";
        stats.print(logfile);
        //logfile<<endl;
        ///logfile<<" p: "<<p<<";"<<endl;
        logfile<<" iter: [1,3];"<<endl;
        logfile.close();
        
        //cout<<maxEdges<<endl;
        for(int i=0;i<2;i++){
        stringstream str;
        str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<VALUE_OF_S<<"_"<<VALUE_OF_T<<"_"<<top2[i].first<<".txt";
        ofstream outfile (str.str());
        print_graph( top2[i].second, outfile);
        outfile.close();
        
        }
        results[n-2][m-2]=top2[0].first;
    }
    
    stringstream str;
     if(appendToM)
            str<<"DYNP_ITER_INPUT_VERTEX_TOP2_M"<<runid<<"_result.txt";
        else
            str<<"DYNP_ITER_INPUT_VERTEX_TOP2_N"<<runid<<"_results.txt";
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
        logfile<<"DYNP_ITER_INPUT_VERTEX_TOP2_M"<<runid<<": "<<duration<<"seconds\n";
    else
        logfile<<"DYNP_ITER_INPUT_VERTEX_TOP2_N"<<runid<<": "<<duration<<"seconds\n";
    logfile.close();

    return results;
}

int main() {
    srand(time(0));

    vector<vector<int>> results(39, vector<int>(39,0));
    auto start = chrono::steady_clock::now();
    for(int i=0;i<10;i++){
        vector<vector<int>> temp=run_in_range(2,30,i%2,i/2+1);
        for(int j=0;j<results.size();j++){
            for(int k=0;k<results.size();k++){
                if(temp[j][k]>results[j][k]){
                    results[j][k]=temp[j][k];
                }
            }
        }
    }
    
    stringstream str;
    str<<"DYNP_ITER_INPUT_VERTEX_TOP2_results.txt";
    ofstream resfile (str.str());
    print_graph( results, resfile);
    resfile.close();
    print_graph(results);

    auto end = chrono::steady_clock::now();

    auto duration = chrono::duration_cast<chrono::seconds>(end - start).count();
    cout << "Execution time: " << duration << " seconds\n";

    ofstream logfile;
    logfile.open("exec_time.txt", std::ios_base::app);
    logfile<<"DYNP_ITER_INPUT_VERTEX_TOP2: "<<duration<<"seconds\n";
    logfile.close();

    return 0;
}