#include "MCND_y.hpp"
#include "UtilsMethods.hpp"



#define INF 9999999
ILOSTLBEGIN

void 
FlowY::clear_model(){
	x.end();
    Obj.end();
	flow_row.end();
	capa_row.end();
	cplex->clearModel();
	cplex->end();
	env->end();
    
	delete env;
	delete cplex;
	delete model;
	model=0;
}


FlowY::FlowY() :
	cplex(0),
	model(0), env(0)
{
	
	
}

void 
FlowY::set_data(const std::deque<int>& non0, const Data * d){
	data = d;
	nnodes = data->nnodes;
	ndemands = data->ndemands;
	narcs = data->narcs;
	sznz = non0.size();
	
	yint.assign(narcs,0);
	int arc;
	for(int a=0;a<sznz;++a){
		arc = non0[a];
		yint[arc] = 1;
		//std::cout<<" arc: "<<arc<<std::endl;
	}
	
	idf_row.assign(nnodes*ndemands,-1);
	idc_row.assign(narcs,-1);
	
	idx.assign(narcs*ndemands,-1);
}

void FlowY::create_model(int phase) {
    env = new IloEnv;
    cplex = new IloCplex(*env);
    cplex->setParam(IloCplex::Threads,1);
    cplex->setParam(IloCplex::RootAlg, 1);
    //else cplex->setParam(IloCplex::RootAlg, 3);
    cplex->setParam(IloCplex::ClockType, 1);
    //cplex->setParam(IloCplex::MIPDisplay, 4);
    cplex->setOut(env->getNullStream());
    cplex->setError(env->getNullStream());
    cplex->setParam(IloCplex::TiLim, 36000.0); // Time limit in seconds
    
    
    nx=0;
    IloExpr obj(*env);
    x = IloNumVarArray(*env);
    
    for (int k = 0; k < ndemands; ++k){
        x.add(IloNumVar(*env));
        obj += 1e10*x[nx];
        nx++;
    }
    
    if(phase==1){
        for(int a=0;a<narcs;++a){
            if(yint[a]){
                //std::cout<<"arc: "<<arc<<std::endl;
                for (int k = 0; k < ndemands; ++k){
                    x.add(IloNumVar(*env));
                    idx[k*narcs+a]= nx++;
                }
                
            }
            //std::cout<<a<<" "<<idy[a]<<std::endl;
        }
    }
    
    //addicional FlowYibility vars
    model = new IloModel(*env);
    Obj = IloMinimize(*env, obj);
    model->add(Obj);
    obj.end();
    //constraints
    nfr=0; ncr=0;
    bool flag=false;
    flow_row =IloRangeArray(*env);
    for (int k = 0; k < ndemands; ++k) {
        for (int i = 1; i <= nnodes; i++) {
            
            flag=false;
            IloExpr constraint(*env);
            
            
            
            for(int a=0;a<narcs;++a){
                //std::cout<<"what "<<a<<" "<<yint[a]<<" "<<idx[k*narcs+a]<<std::endl;
                if(yint[a]==0) continue;
                
                if(i == data->arcs[a].i){
                    flag=true;
                    if(idx[k*narcs+a]>=0)constraint -= x[idx[k*narcs+a]];
                }else if(i == data->arcs[a].j){
                    flag=true;
                    if(idx[k*narcs+a]>=0)constraint += x[idx[k*narcs+a]];
                }
                
                
            }
            
            // std::cout<<"what1 "<<std::endl;
            
            
            if( i == data->d_k[k].D){
                flag=true;
                constraint += x[k];
                constraint -=data->d_k[k].quantity;
            }
            if( i == data->d_k[k].O){
                flag=true;
                constraint -= x[k];
                constraint +=data->d_k[k].quantity;
            }
            
            if(flag){
                flow_row.add((constraint == 0));
                model->add(flow_row[nfr]);
                idf_row[k*nnodes+i-1]= nfr++;
            }
            constraint.end();
            
        }
    }
    
    capa_row =IloRangeArray(*env);
    for(int a=0;a<narcs;++a){
        if(yint[a]){
            IloExpr constraint(*env);
            for (int k = 0; k < ndemands; ++k)
                if(idx[k*narcs+a]>=0)constraint -= x[idx[k*narcs+a]];
            constraint+=data->arcs[a].capa;
            
            capa_row.add((constraint >= 0));
            model->add(capa_row[ncr]);
            idc_row[a]= ncr++;
            constraint.end();
        }
    }
    cplex->extract(*model);
}




int FlowY::solve(double &solution, std::deque<int>& non0, const std::vector<double> & y1, std::vector<double> & x1, int phase){
	
	//cplex->setParam(IloCplex::RootAlg, 2);
	//cplex->exportModel("test.lp");
	cplex->solve();
    if(phase==2){
        IloNumArray x_(*env);
        IloNumArray pi_(*env);
        IloNumArray alpha_(*env);
        
        cplex->getDuals(pi_,flow_row);
        cplex->getDuals(alpha_,capa_row);
        
        std::deque<HeapCell> cols_to_add;
        while(price(cols_to_add, pi_, alpha_)){
            add_cols(cols_to_add);
            cplex->solve();
            //std::cout<<"after price "<<cplex->getObjValue()<<std::endl;
            
            cplex->getDuals(pi_,flow_row);
            cplex->getDuals(alpha_,capa_row);
        }
        
        x_.end();
        pi_.end();
        alpha_.end();
    }

    
    //std::cout<<"TEST final "<<cplex->getStatus()<<std::endl;
	if(cplex->getStatus() == IloAlgorithm::Infeasible){	
		 std::cout<<"infeasible "<<std::endl;
		return 2;
	}else if(cplex->getStatus() == IloAlgorithm::Unbounded){
		std::cout<<"Unbounded"<<std::endl;
		return 2;
	}else if(cplex->getStatus() == IloAlgorithm::Optimal){
        if(phase==0) return 1;
		IloNumArray xsol(*env);
		cplex->getValues(x, xsol);
        solution=cplex->getObjValue();
        if(phase==1){
            int retval = final_feas(data,idx,xsol, y1, yint);
            if(retval==1){
                for (int a=0; a<narcs; ++a){
                    if(yint[a]==-1) non0.push_back(a);
                    //if(yint[a]==-1) std::cout<<"put "<<a<<" "<<y1[a]<<std::endl;
                }
                return 0;
            }else if(retval==-1) return 2;
            return 0;
            
        }
        
        //std::cout<<"no put"<<std::endl;
        non0.clear();
        bool used;
        for (int a=0; a<narcs; ++a){
            if(yint[a]==1){
                used = false;
                for (int k=0; k<ndemands; ++k){
                    if(idx[k*narcs+a]>=0){
                        if(xsol[idx[k*narcs+a]]>1e-30){
                            used = true;
                            x1[k*narcs+a] = xsol[idx[k*narcs+a]];
                        }else x1[k*narcs+a] = 0;
                    }else x1[k*narcs+a] = 0;
                }
                if(used){
                    solution+=data->arcs[a].f;
                    non0.push_back(a);
                }//else   std::cout<<"unsed "<<a<<" "<<y1[a]<<std::endl;
            }else for (int k=0; k<ndemands; ++k) x1[k*narcs+a] = 0;
        }
        std::cout<<std::setprecision(10)<<"Solution: "<<solution<<std::endl;
        return 1;
        
	}else{
		//std::cout<<"what "<<std::endl;
		 return 2;
	}
    std::cout<<"N/A "<<std::endl;

    return 3;
}


int
FlowY::price(std::deque<HeapCell>& cols_to_add, const IloNumArray & pi_,const IloNumArray & alpha_ ){
    
    double rc, rcf, totalrc;
    int i,j;
    for(int a=0;a<narcs;++a){
        if(yint[a]==0) continue;
        i = data->arcs[a].i-1;
        j = data->arcs[a].j-1;
        totalrc=0;
        //std::cout<<" arc: "<<arc<<std::endl;
        for (int k = 0; k < ndemands; ++k){
            if(idx[k*narcs+a]<0){
                //std::cout<<"what: "<<idf_row[k*nnodes+i]<<" "<<idf_row[k*nnodes+j]<<" "<<idc_row[a]<<std::endl;
                rc = data->arcs[a].c[k] /*+ data->arcs[a].f/data->arcs[a].capa*/ + alpha_[idc_row[a]];
                rc += pi_[idf_row[k*nnodes+i]]- pi_[idf_row[k*nnodes+j]];
                if(rc<-1e-10){
                    cols_to_add.push_back(HeapCell(k*narcs+a,rc));
                    //std::cout<<"add col k: "<<k<<" id: "<<cols_to_add.back()<<std::endl;
                }
            }
        }
    }
    
    return cols_to_add.size();
}

void
FlowY::add_cols(std::deque<HeapCell>& cols_to_add){
    int k,a,i,j, cont;
    //std::stable_sort(cols_to_add.begin(),cols_to_add.end(),comp());
    //std::cout<<"num cols to add: "<<cols_to_add.size()<<std::endl;
    cont=0;
    for(int c=cols_to_add.size();c--;){
        k = cols_to_add[c].k/narcs;
        a = cols_to_add[c].k%narcs;
        i = data->arcs[a].i-1;
        j = data->arcs[a].j-1;
        //std::cout<<"add col k: "<<k<<" arc: "<<a<<" id: "<<cols_to_add[c].rc_<<std::endl;
        IloNumColumn col =capa_row[idc_row[a]](-1);
        col+=flow_row[idf_row[k*nnodes+i]](-1);
        col+=flow_row[idf_row[k*nnodes+j]](1);
        col+= Obj(data->arcs[a].c[k] /*+ data->arcs[a].f/data->arcs[a].capa*/);
        
        x.add(IloNumVar(col));
        idx[cols_to_add[c].k]= nx++;
        col.end();
        //if(++cont>=10)
        //break;
    }
    //std::cout<<"added: "<<cont<<std::endl;
    cols_to_add.clear();
}




FlowY::~FlowY() {

	try {
		idx.clear();
		idf_row.clear();
		idc_row.clear();
		yint.clear();

	} catch (IloException& e) {
		std::cerr << "ERROR: " << e.getMessage() << std::endl;
	} catch (...) {
		std::cerr << "Error" << std::endl;
	}
}







