
#include "mcnd.hpp"
#include "MCND_y.hpp"
#include <sys/times.h>
#include <algorithm>    // std::sort
#include "UtilsMethods.hpp"
#include "structures.hpp"

void reopt(const Data & data, MCND & volmcnd, VOL_problem & volp, std::vector<double>& h1, std::vector<double> &dual, std::deque<int>& non0, int nnfxdarc, double f2);


double feasible_solve(const std::vector<double>& h1, std::deque<int>& firstnon0, const Data & data, std::vector<double>& x1, int phase=1);


void disturb_best(const Data & data, MCND & volmcnd, VOL_problem & volp,const std::vector<double>& h1, std::vector<double>& dual, const std::deque<int>& bestnon0, const std::deque<int>& basenon0, std::deque<int>& non0, double & best, double ratio);

bool crossover(const Data & data, MCND & volmcnd, VOL_problem & volp,const std::vector<double>& h1, std::vector<double> &dual, const std::deque<std::deque<int> >& cand, std::deque<int>& non0, double & best,int itmax);


//=========================================================================================
//========================        MAIN                =====================================
//=========================================================================================


int main(int argc, char* argv[]) {
    
    //-----------initialize
	srand(1);
	std::string instance(argv[1]);
	Data data;
	MCND_read_data(instance, data);
	MCND  volmcnd(&data);
	VOL_problem volp("mcnd.par", instance);
    
    
    //__________attributes
    std::deque<int> non0;
    std::deque<int> bstfeasnon0;
    std::deque<int> basenon0;

    
	double UB=0;
    double bestfeas ;
    
    std::vector<double> h1(data.narcs,0);
    std::vector<double> h(data.narcs,0);
    std::vector<double> dual(data.nnodes*data.ndemands,0);
    std::vector<double> xbest(data.narcs*data.ndemands,0);

    
    int retval, arc;
    int nnfxdarc=0;
    bool modified=false;
    
    //-----------initialize time
    clock_t t_u;
    struct tms buff;
    times( &buff );
    t_u = buff.tms_utime;
    double t = 0;
    
    //-----------First solve
    volmcnd.reset(volp, data.narcs, dual, 0);
    volp.solve(volmcnd, 1);
    
    //----------log
    std::ofstream file("fileout", std::ios::app);
    file<<std::setprecision(15)<<instance<<" "<<volp.value<<" "<<volp.iter()<<" ";
    std::cout<<std::setprecision(15)<<instance<<std::endl;
    times( &buff );
    t = ( double( buff.tms_utime - t_u ) ) /double( CLK_TCK );
    file<<" t: "<<t<<" ";
    
    //----------------get Solution of First solve
    volmcnd.transp_h(h1, double(volp.iter()));
    volmcnd.trans_dualsol(dual, volp.dsol);
    //volmcnd.trans_sol(y1,volp.psol);
    //volmcnd.trans_sol(rel1,volp.psol);
    //for(int i=data.nnodes*data.ndemands;i--;) dual[i] = volp.dsol[i]; //get solution
    //for(int i=data.narcs;i--;) std::cout<<i<<" "<<volp.psol[i]<<std::endl;
	//std::cout<<std::setprecision(15)<<volp.value<<std::endl;
   // volp.parm.maxsgriters = 100;

   
    //________________________________________________
    //________________________________________________
    //_____________REOPTIMIZE - TRI
    
    nnfxdarc=0;
    for (int a = volmcnd.sznon0; a--; ){
        arc = volmcnd.non0[a];
        if(h1[arc]>=0.3)non0.push_back(arc);
        else if(h1[arc]>=0.001){
            non0.push_front(arc);
            ++nnfxdarc;
        }//else std::cout<<"A: "<<arc<<" h1: "<<h1[arc]<<std::endl;
    }

    reopt(data, volmcnd,volp, h1,dual, non0, nnfxdarc, 0.3);
    basenon0 = volmcnd.non0;
    //std::cout<<"base size: "<<basenon0.size()<<std::endl;
   
    
    /*std::sort(basenon0.begin(), basenon0.end());
     for (int a = basenon0.size(); a--; ){
     arc = basenon0[a];
     std::cout<<"base: "<<arc<<" h1: "<<h1[arc]<<std::endl;
     }
     std::cout<<"base size: "<<basenon0.size()<<std::endl;*/
    
    
    //________________________________________________
    //________________________________________________
    //_____________FIRST FEAS SOL
  
    non0.clear();
    times( &buff );
    t = ( double( buff.tms_utime - t_u ) ) /double( CLK_TCK );
    for (int a = basenon0.size(); a--; ){
        arc = basenon0[a];
        if(h1[arc]>=0.3){
            non0.push_back(arc);
            //std::cout<<"Fst: "<<arc<<" h1: "<<h1[arc]<<std::endl;
        }
    }
    
    UB = feasible_solve(h1, non0, data, xbest);
    deque_union_front(basenon0, non0, data.narcs);
    bstfeasnon0 = non0;
    bestfeas = UB;
    
    times( &buff );
    t = ( double( buff.tms_utime - t_u ) ) /double( CLK_TCK );
    file<<bestfeas<<" t: "<<t<<" ";
    
    /*std::sort(bstfeasnon0.begin(), bstfeasnon0.end());
    for (int a = bstfeasnon0.size(); a--; ){
        arc = bstfeasnon0[a];
        std::cout<<"final: "<<arc<<" h1: "<<h1[arc]<<std::endl;
    }
    std::cout<<"END FIRST FEAS size: "<<bstfeasnon0.size()<<std::endl;*/
    //________________________________________________
    //________________________________________________
    //_____________PHASE 2 ELIMINATION EVALUATION
    std::deque<std::deque<int> >candidates;

    std::deque<int> bestinon0;
    std::deque<int> bestnon0;
    volp.maxlb = bestfeas*100; 
    double best=bestfeas;
    for(int i=0;i<10;++i){
        //std::cout<<"DISTURB"<<std::endl;
        double best = 1e30;
	times(&buff);
	t = (double(buff.tms_utime - t_u) )/double(CLK_TCK);
	if(t>86400)break;
        for(int it=2;it<=10;it+=2){
            //std::cout<<"RATIO "<<it<<"%"<<std::endl;
            double besti = 1e30;
            for(int itt=0;itt<5;itt++){
                disturb_best(data, volmcnd, volp, h1, dual, bstfeasnon0, basenon0, non0, UB, it/100.0);
                if(UB < besti){
                    //std::cout<<"distrub 0.05: "<<UB<<std::endl;
                    besti = UB;
                    bestinon0 = non0;
                }
            }
            //std::cout<<"BEST : "<<besti<<std::endl;
            if(besti<best){ bestnon0 = bestinon0;  best =besti;}
            if(!bestinon0.empty())candidates.push_back(bestinon0);
            bestinon0.clear();
        }
        //std::cout<<"CROSS inside "<<i<<std::endl;
	times(&buff);
	t = (double(buff.tms_utime - t_u) )/double(CLK_TCK);
	if(t>86400)break;        
        crossover(data,volmcnd, volp, h1,dual, candidates, bestnon0, best, 3);
        if(!bestnon0.empty()){
            if(best < bestfeas){
                UB = feasible_solve(h1, bestnon0, data, xbest);
                if(UB>0 && UB<bestfeas){
                    bestfeas = UB;
                    bstfeasnon0 = bestnon0;
                }
            }
        }
        for(int s=candidates.size();s--;)
            candidates[s].clear();
        candidates.clear();
   	times(&buff);
	t = (double(buff.tms_utime - t_u) ) /double(CLK_TCK);
	if(t>86400)break;
    }
    
    times( &buff );
    t = ( double( buff.tms_utime - t_u ) ) /double( CLK_TCK );
    file<<bestfeas<<" t: "<<t<<std::endl;

    //________________________________________________
    //________________________________________________
    //_____________ END ____________________________

       /*std::sort(bstfeasnon0.begin(), bstfeasnon0.end());
    for (int a = bstfeasnon0.size(); a--; ){
        arc = bstfeasnon0[a];
        std::cout<<"final: "<<arc<<" h1: "<<h1[arc]<<std::endl;
    }
    std::cout<<"END FIRST FEAS size: "<<bstfeasnon0.size()<<std::endl;*/


    
    file.close();
    dual.clear();
    h1.clear();
    h.clear();
    non0.clear();
    bestinon0.clear();
    basenon0.clear();
    bstfeasnon0.clear();
    bestnon0.clear();
    xbest.clear();

	return 0;
}

//========================================================================
//
//==================== AUXILIARY FUNCTIONS ===============================
//
//========================================================================


void reopt(const Data & data, MCND & volmcnd, VOL_problem & volp, std::vector<double>& h1, std::vector<double> &dual, std::deque<int>& non0, int nnfxdarc, double f2){
    
    
    volp.parm.maxsgriters = 100;
    
    int arc;
    bool mod = true;
    double f1 = 1/double(volp.iter());
    
    //std::cout<<"1 iter: "<<double(volp.iter())<<" f1: "<<f1<<std::endl;
    volmcnd.non0 = non0;
    volmcnd.reset(volp, nnfxdarc, dual,0);
    volp.solve(volmcnd, 1);
    volmcnd.transp_h(h1, double(volp.iter()),nnfxdarc);
    
    f1 = 0.01;
    //std::cout<<"2 iter: "<<double(volp.iter())<<" f1: "<<f1<<" unfix: "<<nnfxdarc<<" fixed: "<<volmcnd.non0.size()-nnfxdarc<<std::endl;
   // std::cout<<"REOPT:"<<std::endl;
    while(f1 <= 0.05){
        mod = false;
        nnfxdarc=0;
        non0.clear();
        for(int a = volmcnd.szunfxd; a--; ){
            int arc = volmcnd.non0[a];
            if(h1[arc]>=f2){
                //std::cout<<"A1: "<<arc<<" h: "<<h1[arc]<<std::endl;
                non0.push_back(arc);
                mod = true;
            }
            else if(h1[arc]>=f1){
                //std::cout<<"A0: "<<arc<<std::endl;
                non0.push_front(arc);
                ++nnfxdarc;
            }else{
                //std::cout<<"Ot: "<<arc<<std::endl;
                mod = true;
            }
            //std::cout<<"Phi: "<<arc<<" h1: "<<h1[arc]<<std::endl;
        }
        if(mod==false){
            f1 +=0.01;
            //std::cout<<"f1: "<<f1<<std::endl;
        }else{
            volmcnd.insert_in_place(non0, volmcnd.szunfxd);
            volmcnd.reset(volp, nnfxdarc, dual,0);
            volp.solve(volmcnd, 1);
            volmcnd.transp_h(h1, double(volp.iter()), nnfxdarc);
            //std::cout<<"n iter: "<<double(volp.iter())<<" unfix: "<<nnfxdarc<<" fixed: "<<volmcnd.non0.size()-nnfxdarc<<std::endl;
        }
    }
    non0.clear();
}


//========================================================================


bool crossover(const Data & data, MCND & volmcnd, VOL_problem & volp, const std::vector<double>& h1, std::vector<double> &dual, const std::deque<std::deque<int> >& cand, std::deque<int>& non0, double & best, int itmax){
    
    std::vector<int> ya1(data.narcs,0);
    std::vector<int> ya2(data.narcs,0);
    std::deque<std::string> hist;
    volp.parm.maxsgriters = 250;
    int sz = cand.size();
    //makestd::cout<<"size cnad: "<<sz<<std::endl;
    double UB;
    int retval;
    bool to_better = false;
    for(int p=0;p<sz;++p){
        if(cand[p].empty()) continue;
        transp(cand[p],ya2);
        
        for(int pp=p+1;pp<sz;++pp){
            if(cand[pp].empty()) continue;
            transp(cand[pp],ya1);
            int t=0;
            int delta = diff<int>(ya1,ya2);
            //std::cout<<p<<" to "<<pp<<" diff: "<<delta<<std::endl;
            while(t < pow(2,delta) && t < (itmax)){
                
                volmcnd.re_rand_fix(ya1, ya2, h1, hist);
                volmcnd.reset(volp, 0, dual, 1);
                retval = volp.solve(volmcnd, 1);
                UB = volp.value;
                //std::cout<<"cross ub: "<<UB<<std::endl;
                if(UB<best && (retval>=0)){
                    //std::cout<<"cross ub: "<<UB<<std::endl;
                    best = UB;
                    non0=volmcnd.non0;
                    to_better=true;
                }
                ++t;
                //std::cout<<p<<" "<<pp<<" itbest "<<itbest<<" t: "<<t<<std::endl;
                
            }
            ya1.assign(data.narcs,0);
        }
        ya2.assign(data.narcs,0);
    }
    hist.clear();
    ya2.clear();
    ya1.clear();
    return to_better;
}


//========================================================================


void disturb_best(const Data & data, MCND & volmcnd, VOL_problem & volp, const std::vector<double>& h1, std::vector<double>& dual, const std::deque<int>& bestnon0, const std::deque<int>& basenon0, std::deque<int>& non0, double & best, double ratio){
    
    int retval;
    int cont;
    int arc, nnfix;
    int size =bestnon0.size();
    int chngble = size*ratio;
    if(chngble==0) chngble=1;
    
    std::vector<bool> ya1(data.narcs,false);
    std::vector<double> h(data.narcs,0);
    std::deque<int> out;
    std::list<HeapCell> heap;

    for(int i=bestnon0.size();i--;){
        ya1[bestnon0[i]] = true;
    }
    
    //Randomly choose arcs to be closed in the best solution.
   // std::cout<<"best size: "<<bestnon0.size()<<" close "<<chngble<<std::endl;
    cont=0;
    while(cont<chngble){
        arc = rand()%size;
        double r = (rand()%101)/100.0;
        if(ya1[bestnon0[arc]] && r>=h1[bestnon0[arc]]){
            ya1[bestnon0[arc]] = false;
            out.push_front(bestnon0[arc]);
            ++cont;
            //std::cout<<"close: "<<bestnon0[arc]<<" h1: "<<h1[bestnon0[arc]]<<std::endl;
        }
    }
    volmcnd.non0.clear();
    for(int a=data.narcs;a--;)
        if(ya1[a]) volmcnd.non0.push_front(a);
    
    
    //Put the unfixed arcs from the base (do not consider the previously closed ones).
    //std::cout<<"best size: "<<volmcnd.non0.size()<<std::endl;
    volp.parm.maxsgriters = 100;
    nnfix= volmcnd.deque_concat(out, basenon0);
    volmcnd.reset(volp, nnfix, dual, 1);
    retval = volp.solve(volmcnd, 1);
    volmcnd.transp_h(h, double(volp.iter()), nnfix);
    for (int a = nnfix; a--; ){
        arc = volmcnd.non0[a];
        if(h[arc]>=0.3){
            heap.push_back(HeapCell(arc,h[arc]));
            //std::cout<<"disturbd: "<<arc<<" h1: "<<h1[arc]<<std::endl;
        }
    }
    
    //Open the "chngble" highest values
    heap.sort(comp());
    volmcnd.non0.erase(volmcnd.non0.begin(),volmcnd.non0.begin()+nnfix);
    cont=0;
    while(cont<chngble && !heap.empty()){
        arc = heap.front().k;
        volmcnd.non0.push_front(arc);
        ++cont;
       // std::cout<<"add: "<<arc<<" rc "<<heap.front().rc_<<std::endl;
        heap.pop_front();
    }
    
    //Mini local search. try to close some arcs with the Flow problem.
    volp.parm.maxsgriters = 250;
    volmcnd.reset(volp, 0, dual, 1);
    retval = volp.solve(volmcnd, 1);
    best =volp.value;
    non0 = volmcnd.non0;
    volmcnd.non0.clear();
    
    ya1.clear();
    out.clear();
    heap.clear();
    
}



//========================================================================


double feasible_solve(const std::vector<double>& h1, std::deque<int>& firstnon0, const Data & data, std::vector<double>& x1, int phase ){
    
    double UB=-1;
    FlowY sol0;
    int retval;
    if(phase==2){
        sol0.set_data(firstnon0, &data);
        sol0.create_model(2);
        retval = sol0.solve(UB,firstnon0,h1,x1,2);
        sol0.clear_model();
        if(retval == 1) return UB;
        else return -1;
    }else if(phase==1){
        //std::cout<<"check feas"<<std::endl;
        sol0.set_data(firstnon0, &data);
        sol0.create_model(1);
        retval = sol0.solve(UB,firstnon0, h1,x1,1);
        sol0.clear_model();
        
        if(retval==0){
           // std::cout<<"optmisin"<<std::endl;
            sol0.set_data(firstnon0, &data);
            sol0.create_model(2);
            retval = sol0.solve(UB,firstnon0, h1,x1,2);
            sol0.clear_model();
            if(retval == 1) return UB;
            else return -1;
        }else return -1;
    }else  return -2;
    
}

//========================================================================

