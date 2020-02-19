//
//  UtilsMethods.cpp
//  
//
//  Created by Rui Shibasaki on 30/11/2018.
//
//

#include "UtilsMethods.hpp"


//==========================================


int
final_feas(const Data* data, const std::vector<int>& idx, const IloNumArray & xsol, const std::vector<double> & y1, std::vector<int> & y_){
    
    const int & nnodes = data->nnodes;
    const int & ndemands = data->ndemands;
    const int & narcs = data->narcs;
    //sort commodities in a decreasing order cijk*(Dk - Ok)
    std::list<HeapCell> heap;
    //std::vector<double> primal(narcs*ndemands,0.0);
    for(int k=0;k<ndemands;++k){
        if(xsol[k]>0){
            heap.push_back(HeapCell(k,xsol[k]));
            //std::cout<<std::setprecision(15)<<"k "<<k<<" "<<xsol[k]<<" ("<<(xsol[k]/double(data->d_k[k].quantity))*100<<"%)"<<std::endl;
        }
    }
    if(heap.empty()) return 0;
    heap.sort(comp());
    
    std::vector<int> grid(nnodes*nnodes,0);
    std::vector<double> wij(narcs,0.0);
    std::vector<double> uij(narcs,0.0);
    std::vector<Pair1> preced(nnodes);
    
    for (int a=0; a<narcs; ++a) {
        uij[a] = data->arcs[a].capa;
        if(y_[a]==1){
            for (int k = 0; k < ndemands; ++k)
            if(idx[k*narcs+a]>=0){
                uij[a] -= xsol[idx[k*narcs+a]];
                //primal[k*narcs+a] = xsol[idx[k*narcs+a]];
            }
        }
    }
    
    int dmand=0;
    double epsP=0.0;
    //------- main loop ----------
    while(!heap.empty()){
        dmand = heap.back().k;
        epsP =  heap.back().rc_;
        //std::cout<<"k "<<dmand<<" "<<epsP<<std::endl;
        for (int a=0; a<narcs; ++a) {
            if(epsP <= uij[a]){
                //wij[a]= 1e-30;
                grid[(data->arcs[a].i-1)*nnodes+data->arcs[a].j-1]=a+1;
                if(y_[a]==0) wij[a]=data->arcs[a].c[dmand]*epsP + data->arcs[a].f*(1-y1[a]);
                else wij[a]=data->arcs[a].c[dmand]*epsP;
            }else grid[(data->arcs[a].i-1)*nnodes+data->arcs[a].j-1]=0;
        }
        dijkstra(data->d_k[dmand].O,data->d_k[dmand].D,grid, preced, wij);//find path
        
        //path retrival
        int i,j, arcij;
        int contnodes=1;
        j=data->d_k[dmand].D;
        
        while (j != data->d_k[dmand].O){
            i = preced[j-1].node;
            arcij = preced[j-1].pos;
            
            if(arcij<0){
                std::cout<<"ERROR 0/0: infeasible heuristic solution !!!!! "<<"commdoity "<<dmand<<" node: "<<i<<" "<<arcij<<"/"<<narcs<<std::endl;
                return -1; //no feasible solution found
            }else if (uij[arcij]<=1e-10 || contnodes>nnodes) { //check if impossible path
                std::cout<<"ERROR 0/1: infeasible heuristic solution !!!!! "<<"commdoity "<<dmand<<" node: "<<i<<" "<<arcij<<"/"<<narcs<<std::endl;
                return -1; //no feasible solution found
            }
            
            
            if(y_[arcij]==0)y_[arcij]=-1;
            //primal[k*narcs+arcij] += epsP;
            uij[arcij]-=epsP;
            // if(dmand==31)std::cout<<data->arcs[arcij].i<<"-"<<data->arcs[arcij].j<<std::endl;
            if(uij[arcij]<1e-8){//arc is saturated
                wij[arcij] = 1e30; //block
                grid[(i-1)*nnodes+j-1]=0;
            }
            
            j = i;
            contnodes++;
        }
        
        heap.pop_back();
    }
    
    
    double solvalue=0.0;
    double ttflow=0.0;
    
    
    grid.clear();
    wij.clear();
    uij.clear();
    preced.clear();
    return 1;
}


//==========================================


void MCND_read_data(std::string fname, Data & data) {
    
    std::ifstream file;
    file.open(fname.c_str());
    if (!file.is_open()) {
        std::cout<<"Failure to open ufl datafile: %s\n "<<fname;
        abort();
    }
    
    int& ndemands = data.ndemands;
    int& narcs = data.narcs;
    int& nnodes = data.nnodes;
    
    
    std::string s;
    std::istringstream ss;
    int style=0;
    getline(file,s);
    if(s == "MULTIGEN.DAT:" ||s == " MULTIGEN.DAT:"){
        style=1;
        getline(file,s);
    }
    ss.str(s);
    
    // read number of locations and number of customers
    
    ss>>nnodes;
    ss>>narcs;
    ss>>ndemands;
    
    ss.clear();
    
    data.d_k.resize(ndemands);
    data.arcs.resize(narcs);
    
    if(style == 1){
        double cost=0;
        for(int i=0; i<narcs; ++i){
            getline(file,s);
            ss.str(s);
            ss>>data.arcs[i].i;
            ss>>data.arcs[i].j;
            ss>>cost;
            ss>>data.arcs[i].capa;
            ss>>data.arcs[i].f;

            ss.clear();
            data.arcs[i].c.resize(ndemands, cost);
            data.arcs[i].b.resize(ndemands, data.arcs[i].capa);
        }
        
        for(int i=0; i<ndemands; ++i){
            getline(file,s);
            ss.str(s);
            ss>> data.d_k[i].O >> data.d_k[i].D >> data.d_k[i].quantity;
            ss.clear();
        }
        //calculate variable bounds
        for(int k=0;k<ndemands;++k){
            for(int e=0; e<narcs; ++e){
                if(data.d_k[k].quantity < data.arcs[e].b[k])
                data.arcs[e].b[k] = data.d_k[k].quantity;
            }
        }
    }else if(style == 0){
        int pk=0;
        for(int i=0; i<narcs; ++i){
            getline(file,s);
            ss.str(s);
            ss>>data.arcs[i].i;
            ss>>data.arcs[i].j;
            ss>>data.arcs[i].f;
            ss>>data.arcs[i].capa;
            ss>>pk;

            ss.clear();
            data.arcs[i].c.resize(ndemands, 0.0);
            data.arcs[i].b.resize(ndemands, 0.0);
            for(int k=0;k<ndemands;++k){
                getline(file,s);
                ss.str(s);
                ss>>pk;
                ss>>data.arcs[i].c[pk-1];
                ss>>data.arcs[i].b[pk-1];
                ss.clear();
            }
        }
        // demands origin b+, destination b-
        
        int node; double b;
        for(int i=0; i<2*ndemands; ++i){
            getline(file,s);
            ss.str(s);
            ss>> pk >> node >> b;
            if(b<0) data.d_k[pk-1].D = node;
            else{data.d_k[pk-1].O = node; data.d_k[pk-1].quantity = b;}
            ss.clear();
        }
    }
    file.close();
    
}


//==========================================


void remove_i(int p, std::list<int> &N){
    
    std::list<int>::iterator it = N.begin();
    std::advance(it,p);
    N.erase(it);
}

Pair1 arg_min(const std::list<int> & nodes, const std::vector<double> & costs){
    Pair1 k;
    k.node = nodes.front();
    k.pos = 0;
    std::list<int>::const_iterator it = nodes.begin();
    for (int i=1; i<nodes.size(); ++i){
        if(costs[k.node-1]>costs[*(++it)-1]){
            k.node = *(it);
            k.pos = i;
        }
    }
    return k;
    
}


double dijkstra(int source, int sink, const std::vector<int> &grid,
                std::vector<Pair1> & preced, const std::vector<double>& cij){
    
    int nnodes=preced.size();
    std::vector<double> potentials(nnodes,0.0);
    std::list<int> nodes(nnodes);
    std::list<int>::iterator it;
    Pair1 i;
    
    it = nodes.begin();
    for (int i=0; i<nnodes; ++i) {
        *(it++) = i+1;
        potentials[i] = __DBL_MAX__;
        preced[i].node = 0;
        preced[i].pos = -1;
    }
    
    potentials[source-1] = 0.0;
    while (nodes.size()>=1) {
        i=arg_min(nodes, potentials);
        remove_i(i.pos, nodes);
        if(i.node == sink) return potentials[sink-1];
        for (int j=nnodes;j--;) {
            if(grid[(i.node-1)*nnodes+j]<=0) continue;
            int a = grid[(i.node-1)*nnodes+j]-1;
            if (potentials[j]-potentials[i.node-1]-cij[a]>1e-30){
                potentials[j]=potentials[i.node-1]+cij[a];
                preced[j].node=i.node;
                preced[j].pos=a;
            }
        }
    }
    
    double ret = potentials[sink-1];
    potentials.clear();
    return ret;
}

//==========================================



void transp(const std::deque<std::deque<int> > & cand,std::vector<int> & y, int pos){
    int opened = cand[pos].size();
    for(int i=opened;i--;){
        y[cand[pos][i]]=1;
    }
}

void transp(const std::deque<int> & cand,std::vector<int> & y){
    int opened = cand.size();
    for(int i=opened;i--;){
        y[cand[i]]=1;
    }
}


int deque_union_front(std::deque<int> & fst,const std::deque<int> & snd, int sz){
    std::vector<int>ink(sz,0);
    int cont=0;
    for(int i=fst.size();i--;){
        ink[fst[i]]=1;
    }
    for(int i=snd.size();i--;){
        if(ink[snd[i]]==0){
            ++cont;
            fst.push_front(snd[i]);
            //std::cout<<"insert fst: "<<snd[i]<<std::endl;
            ink[snd[i]]=1;
        }
    }
    ink.clear();
    return cont;
}

int deque_union_back(std::deque<int> & fst,const std::deque<int> & snd, int sz){
    std::vector<int>ink(sz,0);
    int cont=0;
    for(int i=fst.size();i--;){
        ink[fst[i]]=1;
    }
    for(int i=snd.size();i--;){
        if(ink[snd[i]]==0){
            ++cont;
            fst.push_back(snd[i]);
            ink[snd[i]]=1;
        }
    }
    ink.clear();
    return cont;
}


void deque_to_string(const std::deque<int>& d, std::string & s, int sz) {
    std::vector<bool> y(sz,false);
    for(int i=d.size();i--; ){
        y[d[i]]=true;
    }
    for(int i=0;i<sz;++i){
        if(y[i])s+="1";
        else s+="0";
    }
    y.clear();
}


