/*--------------------------------------------------------------------------*/
/*-------------------------- File TestFi.C ---------------------------------*/
/*--------------------------------------------------------------------------*/
/*--                							  --*/
/*-- TestFi is a simple example of a "concrete" class which implements    --*/
/*-- the interface defined by the abstract base class FiOracle.           --*/
/*--                							  --*/
/*--                            VERSION 0.20                              --*/
/*--                           20 - 02 - 2013                             --*/
/*--                			     				  --*/
/*-- 		     Original Idea and Implementation by:		  --*/
/*--                							  --*/
/*--			       Antonio Frangioni       			  --*/
/*--                							  --*/
/*--   			   Operations Research Group			  --*/
/*--			  Dipartimento di Informatica			  --*/
/*--   			     Universita' di Pisa			  --*/
/*--                                                                      --*/
/*--               Copyright 2001 - 2013 by Antonio Frangioni   

                        KNAPSACK  MULTICOMMODITY                            --*/
/*--                							  --*/  
/*--------------------------------------------------------------------------*/ 
/*--------------------------------------------------------------------------*/
/*--------------------------- IMPLEMENTATION -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- MACROS --------------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------- INCLUDES -------------------------------------*/
/*--------------------------------------------------------------------------*/

#include "TestFi.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

static const Index InINF = Inf<Index>();

/*--------------------------------------------------------------------------*/
/*--------------------------- PUBLIC METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------------------------- CONSTRUCTOR ---------------------------------*/
/*--------------------------------------------------------------------------*/





TestFi::TestFi(const Data  *d)
:
FiOracle()
{
	NumVar = 0;				
    Cntr = 0;
    Lam1 =0;
    LamB = 0;
    LamBd =  0;
   
    
    data  = d;
	ndemands = data->ndemands;
	narcs = data->narcs;
	nnodes = data->nnodes;
    h1.resize(narcs);
    
	non0.resize(narcs);
	for(int a=narcs; a--;) non0[a] = a;
	Iu.resize(ndemands*nnodes,-1);
	szopnd=0;
	szunfxd=0; 
	sznon0=0;
	
	
	 
    
}


/*--------------------------------------------------------------------------*/
/*----------------------- METHODS FOR CHANGING DATA ------------------------*/
/*--------------------------------------------------------------------------*/
void TestFi::SetGiName( cIndex Name )
{
    memx[Name] = x;
    memy[Name] = y;
}

void TestFi::SetLambda( cLMRow Lmbd )
{
    Lam1 = Lmbd;
}

/*--------------------------------------------------------------------------*/

void TestFi::SetLamBase( cIndex_Set LmbdB , cIndex LmbdBD )
{
    LamB = LmbdB;
    LamBd = LmbdBD;
}


/*--------------------------------------------------------------------------*/
/*-------------------- METHODS FOR SOLVING THE PROBLEM ---------------------*/
/*--------------------------------------------------------------------------*/

double 
TestFi::knapsack(int a,  const std::vector<double> & rc, std::vector<double> & x){
	double kpsack =0;
	double fillUp =0;
	int arc = non0.at(a);
	//get reduced cost for each commodity in arc e
	for(int k=0;k<ndemands;++k){
		if(rc[k*sznon0 + a]>=0.0){
			heap.push_back(HeapCell(k,rc[k*sznon0 + a]));
		}else x[k*sznon0 + a]=0.0;
	}  
	
	heap.sort(comp());
	//std::stable_sort(heap.begin(), heap.end(), comp());
		
	//solve knap
	
	while(heap.size()>0){
		
		if(fillUp < data->arcs[arc].capa){
			x[heap.front().k*sznon0 + a] = std::min((data->arcs[arc].capa - fillUp),  data->d_k[heap.front().k].quantity);
			fillUp += x[heap.front().k*sznon0 + a];
			kpsack += heap.front().rc_ * x[heap.front().k*sznon0 + a];
			heap.pop_front();
			//std::cout<<"in"<<std::endl;
		}else{
			x[heap.front().k*sznon0 + a] = 0.0;
			heap.pop_front();
		}	
	}
	return kpsack;
}

/*--------------------------------------------------------------------------*/

void TestFi::calcul_rc(){
    
	int arc;
   for(int a=sznon0; a--;){
	 arc = non0.at(a);
	for(int k=0; k<ndemands; ++k){
		rc[k*sznon0 + a] = -data->arcs[arc].c[k] + Lam1[Iu[k*nnodes + data->arcs[arc].j-1]] - Lam1[Iu[k*nnodes + data->arcs[arc].i-1]];
	}
   }	
    
}

/*--------------------------------------------------------------------------*/


HpNum TestFi::solve (){
   
	int arc;
	double cost_a, ratio;
	double lcost =0;
	calcul_rc();
	
	for(int a=szunfxd; a--;){
		arc = non0.at(a);
		cost_a = knapsack(a, rc, x)- data->arcs[arc].f;
		if( cost_a  >= 0.0){
            ++h1[arc];
			y[a] =1.0;	//calcul of y_ij
			lcost += cost_a;
			
		}else{
			y[a] = 0.0;
			for(int k=0; k<ndemands; ++k)
				x[k*sznon0 + a]=0.0;
		}
	}  	
	
	for(int a=szunfxd; a<sznon0;++a){
		arc = non0.at(a);
        cost_a = knapsack(a, rc, x);
        if(cost_a >= 0.0){
            ++h1[arc];
            lcost += cost_a;
        }
		lcost += - data->arcs[arc].f;

	}
	
	//calcul of v
	for(int k=0; k<ndemands; ++k){
		lcost -= data->d_k[k].quantity * ( Lam1[Iu[k*nnodes + data->d_k[k].D-1]] - Lam1[Iu[k*nnodes + data->d_k[k].O-1]]);
	}
	//std::cout<<"solve knap "<<lcost<<std::endl;
	return lcost;
	
}

/*--------------------------------------------------------------------------*/


HpNum TestFi::Fi( cIndex wFi )
{
    
    HpNum FiValue = 0;
    if(wFi){
        FiValue = solve();
    }
    
    
    return FiValue;
}  // end( TestFi::Fi )



/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR READING RESULTS -----------------------*/
/*--------------------------------------------------------------------------*/

bool TestFi::NewGi( cIndex wFi )
{
	
    if( wFi && Lam1 ){
		//cout<<"lam1 "<<*Lam1<<endl;
        return( true );
    }   
    else
        return( false );
    
}


/*--------------------------------------------------------------------------*/

Index TestFi::GetGi( SgRow SubG , cIndex_Set &SGBse , cIndex Name ,
                    cIndex strt , Index stp )
{
    
    if( Name < MBN )
        throw( NDOException( "GetGi: past information is not recorded" ) );
    
    if( stp > NumVar )
        stp = NumVar;
    
    if( Name == MBN ) {
        SGBse = &InINF;
        return( 0 );
    }
    
    SGBse = 0;
   
    
    
    //calcul of v
    int arc;
	for(int k=0; k<ndemands; ++k){
		for(int i=0; i<nnodes; ++i)		
			if(Iu[k*nnodes + i]>=0)
				SubG[Iu[k*nnodes + i]] = 0.0;
		
		SubG[Iu[k*nnodes + data->d_k[k].O-1]] += data->d_k[k].quantity;
		SubG[Iu[k*nnodes + data->d_k[k].D-1]] -= data->d_k[k].quantity;
		for(int a=0; a<sznon0; ++a){
			arc = non0.at(a);
			SubG[Iu[k*nnodes + data->arcs[arc].i-1]] -= x[k*sznon0+ a];
			SubG[Iu[k*nnodes + data->arcs[arc].j-1]] += x[k*sznon0+ a];
		}
	}
    
    Lam1 = 0;
    
    return( NumVar );
}

/*--------------------------------------------------------------------------*/

HpNum TestFi::GetVal( cIndex Name )
{
    if( Name < MBN )
        throw( NDOException( "GetVal: past information is not recorded" ) );
    
    return( 0 );
}


/*--------------------------------------------------------------------------*/


void TestFi::Aggregate( cHpRow Mlt , cIndex_Set NmSt , cIndex Dm , cIndex NwNm ){
    memx[MBN].assign(x.size(),0.0);
	memy[MBN].assign(y.size(),0.0);            
	if(NmSt){
	///Mlt[ i ] is the multiplier of the dual object with name NmSt[ i ], for i = 0, ..., Dm - 1
		for(int i=0;i<Dm;++i){
			for(int a=szunfxd; a--;){
				memy[MBN][a] += memy[NmSt[i]][a] * Mlt[i];
				for(int k=0;k<ndemands;++k)
					memx[MBN][k*sznon0 + a] += memx[NmSt[i]][k*sznon0 + a] * Mlt[i];
			}
			for(int a=szunfxd; a<sznon0;++a){
				for(int k=0;k<ndemands;++k)
					memx[MBN][k*sznon0 + a] += memx[NmSt[i]][k*sznon0 + a] * Mlt[i];
			}
		}
	}else{
		for(int i=0;i<Dm;++i){
			for(int a=szunfxd; a--;){
				memy[MBN][a] += memy[i][a] * Mlt[i];
				for(int k=0;k<ndemands;++k)
					memx[MBN][k*sznon0 + a] += memx[i][k*sznon0 + a] * Mlt[i];
			}
			for(int a=szunfxd; a<sznon0;++a){
				for(int k=0;k<ndemands;++k)
					memx[MBN][k*sznon0 + a] += memx[i][k*sznon0 + a] * Mlt[i];
			}
		}
	}
	memx[NwNm] = memx[MBN];
	memy[NwNm] = memy[MBN];  
}

/*--------------------------------------------------------------------------*/


void TestFi::formfinalsol(){
	
	y.assign(y.size(),0.0); 
	x.assign(x.size(),0.0); 
	//x.assign(x.size(),0.0); 
    cIndex_Set I;
    Index D; 
    Lam1 = Slvr->ReadBestSol(I,D);
    if(D<NumVar) std::cout<<"Somethings Wrong READBESTSOL"<<std::endl;
    cHpRow Mlt = Slvr->ReadMult(I,D);
    if(I){
	///Mlt[ i ] is the multiplier of the dual object with name I[ i ], for i = 0, ..., D - 1
		for(int i=0;i<D;++i){
			for(int a=szunfxd; a--;){
				y[a] += memy[I[i]][a] * Mlt[i];
				for(int k=0;k<ndemands;++k)
					x[k*sznon0 + a] += memx[I[i]][k*sznon0 + a] * Mlt[i];
				//std::cout<<"a: "<<a+1<<" "<< memy[I[i]][a]<<" * "<< Mlt[i]<<" = "<< y[a]<<std::endl;
			}
			for(int a=szunfxd; a<sznon0;++a){
				for(int k=0;k<ndemands;++k)
					x[k*sznon0 + a] += memx[I[i]][k*sznon0 + a] * Mlt[i];
				//std::cout<<"a: "<<a+1<<" "<< memy[I[i]][a]<<" * "<< Mlt[i]<<" = "<< y[a]<<std::endl;
			}
		}
	}else{
		for(int i=0;i<D;++i){
			for(int a=szunfxd; a--;){
				y[a] += memy[i][a] * Mlt[i];
				for(int k=0;k<ndemands;++k)
					x[k*sznon0 + a] += memx[i][k*sznon0 + a] * Mlt[i];
			}
			for(int a=szunfxd; a<sznon0;++a){
				for(int k=0;k<ndemands;++k)
					x[k*sznon0 + a] += memx[I[i]][k*sznon0 + a] * Mlt[i];
				//std::cout<<"a: "<<a+1<<" "<< memy[I[i]][a]<<" * "<< Mlt[i]<<" = "<< y[a]<<std::endl;
			}
		}
	}
    
}

/*--------------------------------------------------------------------------*/


int
TestFi::reset_topology(){
    
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


/*--------------------------------------------------------------------------*/


void
TestFi::reset( int unfixd, const std::vector<double> & hs, bool hotstart){
    
    NumVar = reset_topology();
    double xsize = sznon0 * ndemands;
    szunfxd = unfixd;
    szopnd =sznon0-szunfxd;
				
    Cntr = 0;
    
    Lam1 =0;
    LamB = 0;
    LamBd =  NumVar;
    
    y.assign(szunfxd,0);
    x.assign(xsize,0);
    rc.assign(xsize,0);
    h1.assign(narcs,0);
    //std::cout<<NumVar<<" : "<<y.size()<<std::endl;
    LMRow lambda0 = new LMNum[NumVar];
    if(hotstart){
        for(int k=0; k<ndemands; ++k)
            for(int i=0; i<nnodes; ++i)
                if(Iu[k*nnodes + i]>=0)
                    lambda0[Iu[k*nnodes + i]] = hs[k*nnodes + i];
    }else for(int i=NumVar; i--;) lambda0[i] = 0;
    
    Slvr->SetFiOracle( this );
    Slvr->SetLambda(lambda0);
    delete[] lambda0;
}


/*--------------------------------------------------------------------------*/
/*------------------------------ DESTRUCTOR --------------------------------*/
/*--------------------------------------------------------------------------*/

TestFi::~TestFi()
{
   
   	x.clear();
   	y.clear();
   	rc.clear();
   	non0.clear();
   	Iu.clear();
   	heap.clear();
   
	for(int n=0;n<=MBN;++n){
		memx[n].clear();
		memy[n].clear();
	}
	memx.clear();
	memy.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- HEURISTIC METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

double
TestFi::re_rand_fix(const std::vector<int> & ya2, const std::vector<int> & ya1,  const std::vector<double> & y1, std::deque<std::string> & hist){
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



/*--------------------------------------------------------------------------*/
/*-------------------------- AUXILIARY METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/





/*--------------------------------------------------------------------------*/


void
TestFi::transp_h(std::vector<double> & y1, double div, int place){
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

/*--------------------------------------------------------------------------*/


int
TestFi::deque_concat(const std::deque<int> & out,const std::deque<int> & from){
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

/*--------------------------------------------------------------------------*/


void
TestFi::insert_in_place(const std::deque<int> & in, int nplace){
    non0.erase(non0.begin(), non0.begin()+nplace);
    non0.insert(non0.begin(),in.begin(), in.end());
}






/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------ End File TestFi.C -------------------------------*/
/*--------------------------------------------------------------------------*/


