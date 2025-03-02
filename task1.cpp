#include <iostream>
#include <map>
#include <vector>
#include <algorithm>
#include <cmath>
#include <time.h>
#include <random>
#include <bits/stdc++.h>

unsigned const MAX_VERTICES = 100;

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


int n, m;
int iternum=0;
vector<int> graph;
vector<c4> circles;

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

bool isCircle(const c4& circle){
    if(graph[circle.idx1]!=1 || graph[circle.idx2]!=1 || graph[circle.idx3]!=1 || graph[circle.idx4]!=1){
        return false;
    }
    return true;
}

bool isNotCircle(const c4& circle){
    return !isCircle(circle);
}

void reeval_circles(vector<c4> circles, const vector<int> &graph){
    circles.erase(std::remove_if(circles.begin(), circles.end(), [graph](const c4& circle){
        if(graph[circle.idx1]!=1 || graph[circle.idx2]!=1 || graph[circle.idx3]!=1 || graph[circle.idx4]!=1){
            return true;
        }
        return false;
    }), circles.end());
}

vector<c4> eval(int idx, const vector<int> &graph){
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
                //circles.push_back(c4(idx, x*m+y1, x1*m+y1, x1*m+y));
                /*El kell tárolni a C4-eket, amiket találtunk, ez most egy C4 befejező éle*/
                if(graph[x*m+y1]==0 || graph[x1*m+y1]==0 || graph[x1*m+y]==0){
                    cout<< "nagy a baj az egyik él 0, de beraktuk c4-nek";
                }
            }
            j++;
        }
        spot++;
    }
    return res;
}

void print_graph(const vector<int> &graph){
    for(int i = 0; i < n; i++){
        for (int j = 0; j < m; j++)
        {
           cout<<graph[i*m+j]<<' ';
        }
        cout<<endl;
    }
}

void print_circle(c4 circle){
    cout<<circle.idx1/m<<", "<<circle.idx1%m<<", "<<circle.idx3/m<<", "<<circle.idx3%m<<",";
}

void check__for_free_edges(vector<int> &graph){
    vector<int> temp;
    temp.assign(n*m, {});
    for(int i=0;i<n*m;i++)
        temp[i]=i;
    auto rng = default_random_engine {};
    shuffle(begin(temp), end(temp), rng);

    for(auto i : temp){
        if(graph[i]!=1){
            if(eval(i, graph).size()==0){
                graph[i]=1;
                cout<<"siker hozzáadva egy free él["<<i/m<<", "<<i%m<<"]"<<endl;
            }
        }
    }
    /*for(int i=0;i<n*m;i++){
        if(graph[i]!=1){
            if(eval(i).size()==0){
                graph[i]=1;
                cout<<"siker hozzáadva egy free él["<<i/m<<", "<<i%m<<"]"<<endl;
            }
        }
    }*/
}

void run_iterations(int iteration, int _n, int _m, double p, ofstream out)
{
    vector<int> gr;
    gr.assign(n*m, 0);
    vector<c4> circs;
    for(int i=0;i<n*m;i++){
        if(p>((double)rand()/(double)RAND_MAX)){
            vector<c4> temp=eval(i, gr);
            circs.insert(circs.end(), temp.begin(), temp.end());
            gr[i]=1;
        }
    }
    while(circs.size()!=0){
        iternum++;
        for(int i=0;i<4;i++){
            int edge=circs[0].get(i);
            if(p>((double)rand()/(double)RAND_MAX)){
                gr[edge]=0;
                eval(edge, gr);
                gr[edge]=1;
            }else{
                gr[edge]=0;
            }
        }
        reeval_circles(circs, gr);
        for(auto element : circs){
            print_circle(element);
            cout<<endl;
        }
        print_graph(gr);
        cout<<"Circles reevaluated"<<endl;
        cout<<"size: "<<circs.size()<<endl;
    }
    check__for_free_edges(graph);
    print_graph(gr);
}

int main(){
    srand((unsigned)time(NULL));
    cin>>n>>m;
    long long d=nCr(m,2)*nCr(n,2) - (nCr(m,2)*nCr(n-2,2)+2*nCr(m-2,1)*nCr(n-2,2)+nCr(n-2,2));
    double p=sqrt(sqrt((double)1/(4*d)));
    graph.assign(n*m, 0);
    /*for(int i=0;i<n*m;i++){
        if(p>((double)rand()/(double)RAND_MAX)){
            vector<c4> temp=eval(i);
            circles.insert(circles.end(), temp.begin(), temp.end());
            graph[i]=1;
        }
    }
    int iternum=0;
    while(circles.size()!=0){
        iternum++;
        for(int i=0;i<4;i++){
            int edge=circles[0].get(i);
            if(p>((double)rand()/(double)RAND_MAX)){
                graph[edge]=0;
                eval(edge);
                graph[edge]=1;
            }else{
                graph[edge]=0;
            }
        }
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
    //print_graph();
    cout<<endl;
    /*for(auto element : circles){
        element.print();
        cout<<endl;
    }*/
   cout<<"circles: "<<circles.size()<<endl;

    int edges=0;
    for(auto element : graph){
        edges+=element;
    }

    cout<<endl<<endl<< "edges: "<< edges<< endl;
    double hanyszor = (nCr(n,2)*nCr(m,2))/(d-1);
    cout<<"n/(d-1) = "<<hanyszor<<endl;
    cout<<"valojaban = "<<iternum<<endl;
    return 0;
}