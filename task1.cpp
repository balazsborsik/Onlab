#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <random>
#include <bits/stdc++.h>

//ötlet randommal lefuttatni, elraktározni a seedet, amivel jó lett és azt a seedet/azokat a seedeket használni

int arr[20][20] = {
    {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20},
    {2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21},
    {3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23},
    {4, 5, 7, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
    {5, 6, 8, 10, 12, 14, 15, 17, 18, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30},
    {6, 7, 9, 12, 14, 16, 18, 19, 21, 22, 24, 25, 27, 28, 30, 31, 32, 33, 34, 35},
    {7, 8, 10, 13, 15, 18, 21, 22, 24, 25, 27, 28, 30, 31, 33, 34, 36, 37, 39, 40},
    {8, 9, 11, 14, 17, 19, 22, 24, 26, 28, 30, 32, 33, 35, 36, 38, 39, 41, 42, 44},
    {9, 10, 12, 15, 18, 21, 24, 26, 29, 31, 33, 36, 37, 39, 40, 42, 43, 45, 46, 48},
    {10, 11, 13, 16, 20, 22, 25, 28, 31, 34, 36, 39, 40, 42, 44, 45, 47, 49, 51, 52},
    {11, 12, 14, 17, 21, 24, 27, 30, 33, 36, 39, 42, 44, 45, 47, 49, 51, 53, 54, 56},
    {12, 13, 15, 18, 22, 25, 28, 32, 36, 39, 42, 45, 48, 49, 51, 52, 55, 56, 59, 60},
    {13, 14, 16, 19, 23, 27, 30, 33, 37, 40, 44, 48, 49, 53, 54, 56, 58, 60, 62, 64},
    {14, 15, 17, 20, 24, 28, 31, 35, 39, 42, 45, 49, 53, 55, 57, 60, 62, 63, 66, 67},
    {15, 16, 18, 21, 25, 30, 33, 36, 40, 44, 47, 51, 54, 57, 60, 63, 65, 67, 69, 71},
    {16, 17, 19, 22, 26, 31, 34, 38, 42, 45, 49, 52, 56, 60, 63, 65, 68, 70, 72, 74},
    {17, 18, 20, 23, 27, 32, 36, 39, 43, 47, 51, 55, 58, 62, 65, 68, 70, 73, 76, 78},
    {18, 19, 21, 24, 28, 33, 37, 41, 45, 49, 53, 56, 60, 63, 67, 70, 73, 76, 80, 82},
    {19, 20, 22, 25, 29, 34, 39, 42, 46, 51, 54, 59, 62, 66, 69, 72, 76, 80, 83, 87},
    {20, 21, 23, 26, 30, 35, 40, 44, 48, 52, 56, 60, 64, 67, 71, 74, 78, 82, 87, 88}
};
unsigned const MAX_VERTICES = 101;

using namespace std;

vector<int> vec;

struct c4
{
    int idx1;
    int idx2;
    int idx3;
    int idx4;

    c4(int a, int b, int c, int d) : idx1(a), idx2(b), idx3(c), idx4(d){}
    
    void print(){
        cout<<idx1<<", "<<idx2<<", "<<idx3<<", "<<idx4<<",";
    }

    int get(int idx){
        if(idx<0 || idx>3)
            throw invalid_argument( "received negative value, or greater value than 3" );
        if(idx==0)
            return idx1;
        if(idx==1)
            return idx2;
        if(idx==2)
            return idx3;
        return idx4;
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

/*bool isCircle(const c4& circle){
    if(graph[circle.idx1]!=1 || graph[circle.idx2]!=1 || graph[circle.idx3]!=1 || graph[circle.idx4]!=1){
        return false;
    }
    return true;
}

bool isNotCircle(const c4& circle){
    return !isCircle(circle);
}*/

void reeval_circles(vector<c4> &circles, const vector<int> &graph){
    circles.erase(std::remove_if(circles.begin(), circles.end(), [&graph](const c4& circle){
        if(graph[circle.idx1]!=1 || graph[circle.idx2]!=1 || graph[circle.idx3]!=1 || graph[circle.idx4]!=1){
            return true;
        }
        return false;
    }), circles.end());
}

vector<c4> eval(int idx, const vector<int> &graph, int n, int m){
    vector<c4> res;
    int x = idx/m;
    int y = idx%m;
    int n_neighbours[MAX_VERTICES];
    int m_neighbours[MAX_VERTICES];
    int spot=0;
    for(int i=0;i<m;i++){
        if(graph[x*m+i]>=1){
            n_neighbours[spot++]=i;
        }
    }
    n_neighbours[spot]=-1;
    spot=0;
    for(int i=0;i<n;i++){
        if(graph[i*m+y]>=1){
            m_neighbours[spot++]=i;
        }
    }
    m_neighbours[spot]=-1;

    spot=0;
    while(m_neighbours[spot]!=-1){
        int j=0;
        while(n_neighbours[j]!=-1){
            int x1=m_neighbours[spot];
            int y1=n_neighbours[j];
            if(graph[m*m_neighbours[spot]+n_neighbours[j]]>=1){
                res.push_back(c4(idx, x*m+y1, x1*m+y1, x1*m+y));
                if(graph[x*m+y1]==0 || graph[x1*m+y1]==0 || graph[x1*m+y]==0){
                    throw invalid_argument( "nagy a baj az egyik él 0, de beraktuk c4-nek" );
                }
            }
            j++;
        }
        spot++;
    }
    return res;
}

void print_graph(const vector<int> &graph,int n, int m, ostream& out=cout){
    for(int i = 0; i < n; i++){
        for (int j = 0; j < m; j++)
        {
           out<<graph[i*m+j]<<' ';
        }
        out<<endl;
    }
}

void print_circle(c4 circle, int n, int m){
    cout<<circle.idx1/m<<", "<<circle.idx1%m<<", "<<circle.idx3/m<<", "<<circle.idx3%m<<",";
}

void check__for_free_edges(vector<int> &graph, int n, int m){
    vector<int> temp;
    temp.assign(n*m, {});
    for(int i=0;i<n*m;i++)
        temp[i]=i;
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    auto rng = default_random_engine {seed};
    shuffle(begin(temp), end(temp), rng);

    for(auto i : temp){
        if(graph[i]!=1){
            if(eval(i, graph, n, m).size()==0){
                graph[i]=1;
                //cout<<"siker hozzáadva egy free él["<<i/m<<", "<<i%m<<"]"<<endl;
            }
        }
    }
}

int number_of_edges(const vector<int> &graph){
    int edges=0;
    for(auto element : graph){
        edges+=element;
    }
    return edges;
}

vector<int> run_iterations(int iteration, int n, int m, double p)
{
    vector<int> graph;
    graph.assign(n*m, 0);
    for(int iter=0; iter<iteration;iter++){
        if(iter%(iteration/100)==0)
            cout<<iter/(iteration/100)+1<<"%" <<"completed"<<endl;
        vector<int> gr;
        gr.assign(n*m, 0);
        vector<c4> circs;
        for(int i=0;i<n*m;i++){
            if(p>((double)rand()/(double)RAND_MAX)){
                vector<c4> temp=eval(i, gr, n ,m);
                circs.insert(circs.end(), temp.begin(), temp.end());
                gr[i]=1;
            }
        }
        int iternum=0;
        while(circs.size()!=0){
            iternum++;
            for(int i=0;i<4;i++){
                int edge=circs[0].get(i);
                if(p>((double)rand()/(double)RAND_MAX)){
                    gr[edge]=0;
                    vector<c4> temp=eval(edge, gr, n, m);
                    circs.insert(circs.end(), temp.begin(), temp.end());
                    gr[edge]=1;
                }else{
                    gr[edge]=0;
                }
            }
            reeval_circles(circs, gr);
        }

        double hanyszor = (nCr(n,2)*nCr(m,2))/(nCr(m,2)*nCr(n,2) - (nCr(m,2)*nCr(n-2,2)+2*nCr(m-2,1)*nCr(n-2,2)+nCr(n-2,2))-1);
        //cout<<"n/(d-1) = "<<hanyszor<<endl;
        //cout<<"valojaban = "<<iternum<<endl;
        check__for_free_edges(gr, n, m);
        vec.push_back(number_of_edges(gr));
        if(number_of_edges(gr)>number_of_edges(graph)){
            graph=gr;
        }
    }
    return graph;
}

int main(){
    srand((unsigned)time(NULL));
    int n=-1,m;
    vector<int> n_m_queue;
    cin>>n>>m;
    long long d=nCr(m,2)*nCr(n,2) - (nCr(m,2)*nCr(n-2,2)+2*nCr(m-2,1)*nCr(n-2,2)+nCr(n-2,2));
    double p=sqrt(sqrt((double)1/(4*d)))+0.105;
    //p=0.85*((double)((double)arr[m-1][n-1]/(m*n)));
    auto start = std::chrono::steady_clock::now();

    vector<int> graph = run_iterations(100000, n, m, p);

    auto end = std::chrono::steady_clock::now();
    print_graph(graph, n ,m);
    int edges=0;
    for(auto element : graph){
        edges+=element;
    }
    cout<<endl<<endl<< "edges: "<< edges<< endl;  
    double average=0;
    double variance=0;
    int array[] = {0,0,0,0,0,0,0,0};
    for(const auto& elm : vec){
        if(elm>=80)
            array[elm-80]++;
        average+=elm;
    }
    average/=vec.size();
    for(const auto& elm : vec){
        variance+=(elm - average) * (elm - average);
    }
    variance /= vec.size();
    double D = sqrt(variance);
    auto diff = end - start;
    /*9ofstream myfile;
    myfile.open("logs2.txt", std::ios_base::app);
    myfile<< n<<", "<<m<<endl;
    myfile << "average: " << average << endl;
    myfile << "p: "<<p<<endl;
    myfile <<"time: "<< std::chrono::duration<double, std::milli>(diff).count() << " ms" << std::endl;
    myfile << "variance: " << variance << " szórás: "<< D << endl;
    myfile.close();*/
    stringstream str;
    str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_"<<edges<<".txt";
    ofstream outfile (str.str());
    print_graph(graph, n, m, outfile);
    outfile.close();
    vec.clear();
    return 0;
}


/*
random:
variance: 2.12238 szórás: 1.45684
average: 80.2798
80: 26580 db
81: 26290 db
82: 14735 db
83: 4277 db
84: 564 db
85: 19 db
86: 0 db
87: 0 db

p=0
variance: 2.09465 szórás: 1.44729
average: 80.3511
80: 26334 db
81: 26911 db
82: 15713 db
83: 4575 db
84: 666 db
85: 41 db
86: 0 db
87: 0 db

p=0.09
variance: 1.95605 szórás: 1.39859
average: 80.4829
80: 25993 db
81: 28538 db
82: 17165 db
83: 5210 db
84: 702 db
85: 38 db
86: 0 db
87: 0 db

p=0.105
variance: 1.98126 szórás: 1.40757
average: 80.475
80: 26019 db
81: 28138 db
82: 17080 db
83: 5320 db
84: 712 db
85: 44 db
86: 0 db
87: 0 db
variance: 1.97326 szórás: 1.40473
average: 80.4833
80: 129483 db
81: 141834 db
82: 85641 db
83: 26608 db
84: 3691 db
85: 202 db
86: 4 db
87: 0 db
variance: 1.97242 szórás: 1.40443
average: 80.4816
80: 129390 db
81: 141575 db
82: 85802 db
83: 26464 db
84: 3628 db
85: 221 db
86: 7 db
87: 0 db

p=0.12
variance: 1.98483 szórás: 1.40884
average: 80.4667
80: 26160 db
81: 28184 db
82: 16820 db
83: 5247 db
84: 746 db
85: 30 db
86: 1 db
87: 0 db
variance: 1.98946 szórás: 1.41048
average: 80.4694
80: 129639 db
81: 141115 db
82: 85231 db
83: 26068 db
84: 3599 db
85: 200 db
86: 3 db
87: 0 db
variance: 1.98715 szórás: 1.40966
average: 80.4703
80: 130067 db
81: 140860 db
82: 85158 db
83: 26140 db
84: 3683 db
85: 193 db
86: 4 db
87: 0 db

p=0.2
variance: 2.26752 szórás: 1.50583
average: 80.2393
80: 26013 db
81: 25753 db
82: 14798 db
83: 4266 db
84: 556 db
85: 21 db
86: 1 db
87: 0 db

p=0.3
variance: 3.11217 szórás: 1.76413
average: 79.5868
80: 24056 db
81: 19559 db
82: 9526 db
83: 2403 db
84: 294 db
85: 16 db
86: 1 db
87: 0 db


edges: 85
time: 83694.9 ms
100000 ugye 100 000

variance: 1.96694 szórás: 1.40248
average: 80.4901
80: 25926 db
81: 28346 db
82: 17267 db
83: 5301 db
84: 761 db
85: 46 db
86: 0 db
87: 0 db


time: 27295.4 ms
*/