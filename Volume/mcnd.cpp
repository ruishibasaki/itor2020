
#include "mcnd.hpp"



MCND::MCND(const Data  *d){
	  data  = d;
	  ndemands = data->ndemands;
	  narcs = data->narcs;
	  nnodes = data->nnodes;
	  h1.resize(narcs);
	  non0.resize(narcs);
	  for(int a=narcs; a--;) non0[a] = a;
	  Iu.resize(ndemands*nnodes,-1);
	  //yl.resize(narcs,0);
	  szopnd=0;
	  szunfxd=0; 
	  sznon0=0;

}




//###### USER HOOKS
// compute reduced costs
int
MCND::compute_rc(const VOL_dvector& u, VOL_dvector& rc)
{
	
	int arc;
   for(int a=szunfxd; a--;){
	 arc = non0.at(a);
   	rc[a] = data->arcs[arc].f;
	for(int k=0; k<ndemands; ++k){
		rc[szunfxd + k*sznon0 + a] = data->arcs[arc].c[k] - u[Iu[k*nnodes + data->arcs[arc].j-1]] + u[Iu[k*nnodes + data->arcs[arc].i-1]];
	}
   }	
   for(int a=szunfxd; a<sznon0;++a){
	arc = non0.at(a);
	for(int k=0; k<ndemands; ++k){
		rc[szunfxd + k*sznon0 + a] = data->arcs[arc].c[k] - u[Iu[k*nnodes + data->arcs[arc].j-1]] + u[Iu[k*nnodes + data->arcs[arc].i-1]];
	}
   }
   return 0;
}



// IN: dual vector u
// OUT: primal solution to the Lagrangian subproblem (x)
//      optimal value of Lagrangian subproblem (lcost)
//      v = difference between the rhs and lhs when substituting
//                  x into the relaxed constraints (v)
//      objective value of x substituted into the original problem (pcost)
//      xrc
// return value: -1 (volume should quit) 0 ow

int 
MCND::solve_subproblem(const VOL_dvector& u, const VOL_dvector& rc,
		      double& lcost, 
		      VOL_dvector& x, VOL_dvector& v, double& pcost)
{
	//std::cout<<"solve knap "<<std::endl;
	int arc;
	double cost_a, ratio;
	pcost= 0;
	lcost =0;
	for(int a=szunfxd; a--;){
		arc = non0.at(a);
		//std::cout<<a<<std::endl;
		cost_a = knapsack(a, rc, x)+ data->arcs[arc].f;
		if( cost_a  <= 0.0){
            ++h1[arc];
			x[a] =1.0;	//calcul of y_ij
			lcost += cost_a;
			pcost += data->arcs[arc].f*x[a];
			for(int k=0; k<ndemands; ++k)
				pcost += data->arcs[arc].c[k] * x[szunfxd + k*sznon0 + a];
		}else{
			x[a] = 0.0;
			for(int k=0; k<ndemands; ++k)
				x[szunfxd + k*sznon0 + a]=0.0;
		}
	}  	
	for(int a=szunfxd; a<sznon0;++a){
		arc = non0.at(a);
        cost_a = knapsack(a, rc, x);
        if(cost_a<=0){
            ++h1[arc];
            lcost += cost_a;
        }
        lcost += data->arcs[arc].f;
		for(int k=0; k<ndemands; ++k){
			pcost += (data->arcs[arc].c[k]) * x[szunfxd + k*sznon0 + a];
			//if(x[szunfxd + k*sznon0 + a]>data->arcs[arc].b[k])
				//std::cout<<"opa "<<arc<<" k: "<<k<<std::endl;
		}
	}
	
	//calcul of v
	v =0;
	for(int k=0; k<ndemands; ++k){
		lcost += data->d_k[k].quantity * ( u[Iu[k*nnodes + data->d_k[k].D-1]] - u[Iu[k*nnodes + data->d_k[k].O-1]]);		
		
		v[Iu[k*nnodes + data->d_k[k].O-1]] -= data->d_k[k].quantity;
		v[Iu[k*nnodes + data->d_k[k].D-1]] += data->d_k[k].quantity;
		for(int a=sznon0; a--;){
			arc = non0.at(a);
			v[Iu[k*nnodes + data->arcs[arc].i-1]] += x[szunfxd+k*sznon0+ a];
			v[Iu[k*nnodes + data->arcs[arc].j-1]] -= x[szunfxd+ k*sznon0+ a];
		}
	}
	
	return 0;
}


// IN:  fractional primal solution (x),
//      best feasible soln value so far (icost)
// OUT: integral primal solution (ix) if better than best so far
//      and primal value (icost)
// returns -1 if Volume should stop, 0/1 if feasible solution wasn't/was
//         found.
// We use randomized rounding. We look at x[i] as the probability
// of opening facility i.
int
MCND::heuristics(const VOL_dvector& p,
		const VOL_dvector& x, double& heur_val)
{
		
		return 0;
}

double 
MCND::knapsack(int a, const VOL_dvector& rc, VOL_dvector& x){
	double kpsack =0;
	double fillUp =0;
	int arc = non0.at(a);
	//get reduced cost for each commodity in arc e
	for(int k=0;k<ndemands;++k){
		if(rc[szunfxd + k*sznon0 + a]<=0.0){
			heap.push_back(HeapCell(k,rc[szunfxd + k*sznon0 + a]));
		}else x[szunfxd + k*sznon0 + a]=0.0;
	}  
	
	heap.sort(comp());
	//std::stable_sort(heap.begin(), heap.end(), comp());
    if(heap.size()==0) return 0.0;
	//solve knap
	//std::cout<<"solve knap "<<szunfxd<<" "<<sznon0<<std::endl;
	while(heap.size()>0){
		
		if(fillUp < data->arcs[arc].capa){
			x[szunfxd + heap.back().k*sznon0 + a] = std::min((data->arcs[arc].capa - fillUp), data->d_k[heap.back().k].quantity);
			fillUp += x[szunfxd + heap.back().k*sznon0 + a];
			kpsack += heap.back().rc_ * x[szunfxd + heap.back().k*sznon0 + a];
			heap.pop_back();
			//std::cout<<x[szunfxd + heap.back().k*sznon0 + a]<<std::endl;
		}else{
			x[szunfxd + heap.back().k*sznon0 + a] = 0.0;
			heap.pop_back();
		}
	}
	return kpsack;
}



int  
MCND::reset_topology(){

	sznon0 = non0.size();
	
	int arc;
	int dsize=0;
	bool flag=false;
	for(int k=0; k<ndemands; ++k){
		for(int i=0; i<nnodes; ++i){
			flag = false;
			Iu[k*nnodes + i]=-1;
			for(int a=sznon0; a--;){
				arc = non0.at(a);
				if((i+1) == data->arcs[arc].i){
					flag = true;
					break;
				}else if((i+1) == data->arcs[arc].j){
					flag = true;
					break;
				}
			}

			if( (i+1) == data->d_k[k].D){
				flag = true;
			}else if( (i+1) == data->d_k[k].O){
				flag = true;
			}
			if(flag){
				Iu[k*nnodes + i] = dsize++;
			 }
		}
	}
	
	return dsize;
}



void 
MCND::reset(VOL_problem & volp, int unfixd, const std::vector<double> & hs, bool hotstart){
	volp.dsize =reset_topology();   
    volp.psize = sznon0 * ndemands + unfixd;
    szunfxd = unfixd;
	szopnd =sznon0-szunfxd;
	//for(int a=sznon0; a--;) std::cout<<non0[a]<<std::endl;
	 //std::cout<<"dsize: "<<volp.dsize<<" sznz: "<<sznon0<<std::endl;
	volp.dsol.clear();
    volp.dsol.allocate(volp.dsize);
    h1.assign(narcs,0);
    if(hotstart){
        for(int k=0; k<ndemands; ++k)
            for(int i=0; i<nnodes; ++i)
                if(Iu[k*nnodes + i]>=0)
                    volp.dsol[Iu[k*nnodes + i]] = hs[k*nnodes + i];
    }else volp.dsol = 0;
}



void
MCND::trans_dualsol(std::vector<double> & dual, const VOL_dvector &dsol){
    dual.assign(dual.size(),0);
    for(int k=0; k<ndemands; ++k)
        for(int i=0; i<nnodes; ++i)
            if(Iu[k*nnodes + i]>=0)
                dual[k*nnodes + i] = dsol[Iu[k*nnodes + i]];
    
}

void
MCND::transp_h(std::vector<double> & y1, double div, int place){
    int arc;
    if(place==0){
        place=sznon0;
        y1.assign(y1.size(),0);
    }
    for(int a=0; a<place;++a){
        arc = non0.at(a);
        y1[arc] = h1[arc]/div;
    }
}


double
MCND::re_rand_fix(const std::vector<int> & ya2, const std::vector<int> & ya1,  const std::vector<double> & y1,  std::deque<std::string> & hist){
    int rr;
	double r;
	double fixc=0;
    bool equal=false;
	non0.clear();
    std::string set = "";
	for(int a=0;a<narcs;++a){
		if(ya2[a] && ya1[a]){
			fixc+=data->arcs[a].f;
			non0.push_front(a);
            set+="1";
		}else if(ya2[a] || ya1[a]){
			 r = (double) ( rand()%101 )/100.0;
			if(r <= y1[a]){
				fixc+=data->arcs[a].f;
				non0.push_front(a);
				set+="1";
			}else set+="0";
		}else set+="0";
		
	}
    //std::cout<<"set: "<<set<<std::endl;

    for(int a=hist.size();a--;){
        //std::cout<<"h: "<<hist[a]<<std::endl;
        if(set==hist[a]){
            equal = true;
            //std::cout<<"equal"<<std::endl;
            break;
        }
    }
    if(equal){
        std::string set2 = "";
        non0.clear();
        fixc=0;
        for(int a=narcs;a--;){
            if(ya2[a] && ya1[a]){
                fixc+=data->arcs[a].f;
                non0.push_front(a);
                set2+="1";
            }else if(ya2[a] || ya1[a]){
                //r = (double) ( rand()%101 )/100.0;
                rr = ( rand()%2 );
                if(rr && !(int(set[a])-48) ){
                    //std::cout<<(int(set[a])-48)<<std::endl;
                    fixc+=data->arcs[a].f;
                    non0.push_front(a);
                    set2+="1";
                }else set2+="1";
            }else set2+="0";
        }
        set = set2;
        set2.clear();
    }
    hist.push_back(set);
    set.clear();
	return fixc;
}


int
MCND::deque_concat(const std::deque<int> & out, const std::deque<int> & from){
    std::vector<int>y(narcs,0);
    for(int n=non0.size();n--;){
        y[non0[n]] = 1;
        // std::cout<<"fix "<<to[n]<<std::endl;
    }
    for(int n=out.size();n--;){
        y[out[n]] = -1;
        // std::cout<<"out "<<out[n]<<std::endl;
    }
    
    int added = 0;
    for(int n=from.size();n--;){
        if(y[from[n]]==0){
            non0.push_front(from[n]);
             //std::cout<<"fix "<<from[n]<<std::endl;
            ++added;
        }
    }
   // std::cout<<"after concat"<<std::endl;
    //for(int n=non0.size();n--;){
      //  std::cout<<"fix "<<non0[n]<<std::endl;
    //}
    
    y.clear();
    return added;
}

void
MCND::insert_in_place(const std::deque<int> & in, int nplace){
    non0.erase(non0.begin(), non0.begin()+nplace);
    non0.insert(non0.begin(),in.begin(), in.end());
}











