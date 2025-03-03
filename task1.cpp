#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <random>
#include <bits/stdc++.h>

///!TODO megcsinálni ugyanúgy a random sorrendű sorsolást ne csak a free élekere

unsigned const MAX_VERTICES = 101;

using namespace std;

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
    circles.erase(std::remove_if(circles.begin(), circles.end(), [graph](const c4& circle){
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
    double p=sqrt(sqrt((double)1/(4*d)))+0.12;
    
    
    vector<int> graph=run_iterations(100000, n, m, p);
    /*
        reeval_circles();
        for(auto element : circles){
            print_circle(element);
            cout<<endl;
        }
        print_graph();
        cout<<"Circles reevaluated"<<endl;
        cout<<"size: "<<circles.size()<<endl;
    }
    check__for_free_edges();
    */
    print_graph(graph, n ,m);
    int edges=0;
    for(auto element : graph){
        edges+=element;
    }
    cout<<endl<<endl<< "edges: "<< edges<< endl;
    stringstream str;
    str<<"output/"<<"Z"<<m<<"_"<<n<<"_"<<2<<"_"<<2<<"_"<<edges<<".txt";
    ofstream outfile (str.str());
    print_graph(graph, n, m, outfile);
    outfile.close();
    return 0;
}