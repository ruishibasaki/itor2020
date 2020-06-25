/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "TestFi.h"
#ifndef WIN32
#include <sys/times.h>
#endif


#include "Bundle.h"
#include "QPPnltMP.h"
#include "structures.hpp"
#include "MCND_y.hpp"
#include "UtilsMethods.hpp"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace NDO_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- CONSTANTS ----------------------------------*/
/*--------------------------------------------------------------------------*/

// parameter file name


const char *const ParF = "ParValue.qp";
/*--------------------------------------------------------------------------*/
/*----------------------------- FUNCTIONS ----------------------------------*/
/*--------------------------------------------------------------------------*/

template<class T>
static inline void str2val( const char* const str , T &sthg )
{
    istringstream( str ) >> sthg;
}

void trans_dualsol(const Data & data, const TestFi & Oracle, std::vector<double> & dual, Bundle *s);


void reopt(const Data & data, TestFi & Oracle, Bundle * s, std::vector<double>& h1, std::vector<double> &dual, std::deque<int>& non0, int nnfxdarc, double f2);

void disturb_best(const Data & data, TestFi & Oracle, Bundle * s, const std::vector<double>& h1, std::vector<double>& dual, const std::deque<int>& bestnon0, const std::deque<int>& basenon0, std::deque<int>& non0, double & best, double ratio);

bool crossover(const Data & data, TestFi & Oracle, Bundle * s, const std::vector<double>& h1, std::vector<double> &dual, const std::deque<std::deque<int> >& cand, std::deque<int>& non0, double & best, int itmax);

double feasible_solve(const std::vector<double>& h1, std::deque<int>& firstnon0, const Data & data, std::vector<double>& x1, int phase=1);

/*--------------------------------------------------------------------------*/
/*-------------------------------- main() ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char **argv )
{
    // read the command-line parameters- - - - - - - - - - - - - - - - - - - - -
    if( argc < 2 ) {
        cerr << "Usage: " << argv[ 0 ] << " < instance >."
        << endl;
        
        return( 1 );
    }
    
   //===============================================================================
   //============================ initialize =======================================
   //===============================================================================
   srand(1);
   std::string instance(argv[1]);
   Data data;
   MCND_read_data(instance, data);
    // enter the try-block - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    try {
        // open parameters file - - - - - - - - - - - - - - - - - - - - - - - - - -
        ifstream ParFile( ParF );
        if( ! ParFile.is_open() )
            cerr << "Warning: cannot open parameters file """ << ParF << """" << endl;
        
         // construct the FiOracle object- - - - - - - - - - - - - - - - - - - - - -
        
        TestFi Oracle(&data);
        FiOracle *Fi = &Oracle;
        Bundle *s = new Bundle( &ParFile );
        QPPenaltyMP *MP = new QPPenaltyMP( &ParFile );
 
        // pass the MPSolver to the Bundle    
        s->SetMPSolver( MP );
        Fi->SetNDOSolver( s );
        
        // set the verbosity of log - - - - - - - - - - - - - - - - - - - - - - - -
        // need to read an int, otherwise e.g. '2' in the input file is read as
        // the ASCII code rather than as the number 2
        
        int lvl;
        DfltdSfInpt( &ParFile , lvl , int( 0 ) ); 
        if( lvl ) s->SetNDOLog( &clog , char( lvl ) );
       
        //__________attributes
        std::deque<int> non0;
        std::deque<int> bstfeasnon0;
        std::deque<int> basenon0;
        
        
        double UB=0;
        double bestfeas;
        
        std::vector<double> h1(data.narcs,0);
        std::vector<double> h(data.narcs,0);
        std::vector<double> dual(data.nnodes*data.ndemands,0);
        std::vector<double> xbest(data.narcs*data.ndemands,0);
        
        
        int arc;
        int nnfxdarc=0;
        bool modified=false;
       
        
        //-----------initialize time
        OPTtimers timer;
        timer.Start();
        double tu=0;
        double ts=0;
        
        //-----------First solve
        Oracle.reset(data.narcs, dual,0);
        s->Solve();
        
        //----------log
        std::ofstream file("fileout", std::ios::app);
        file<<std::setprecision(15)<<instance<<" "<<-s->ReadBestFiVal()<<" "<<s->getNumIter()<<" ";
        std::cout<<std::setprecision(15)<<instance<<std::endl;
        
        tu=0;ts=0;
        timer.Read(tu,ts);
        file<<std::setprecision(15)<<tu<<" ";

        //----------------get Solution of First solve
        Oracle.transp_h(h1, double(s->getNumIter()));
        trans_dualsol(data, Oracle, dual, s);

        
        //________________________________________________
        //________________________________________________
        //_____________REOPTIMIZE - TRI
        
        
        nnfxdarc=0;
        for (int a = Oracle.sznon0; a--; ){
            arc = Oracle.non0[a];
            if(h1[arc]>=0.3)non0.push_back(arc);
            else if(h1[arc]>=0.001){
                non0.push_front(arc);
                ++nnfxdarc;
            }
        }
        
        reopt(data, Oracle,s, h1,dual, non0, nnfxdarc, 0.3);
        basenon0 = Oracle.non0;
        
        //________________________________________________
        //________________________________________________
        //_____________FIRST FEAS SOL
        
        non0.clear();
        for (int a = basenon0.size(); a--; ){
            arc = basenon0[a];
            if(h1[arc]>=0.3){
                non0.push_back(arc);
            }
        }
        
        UB = feasible_solve(h1, non0, data, xbest);
        deque_union_front(basenon0, non0, data.narcs);
        bstfeasnon0 = non0;
        bestfeas = UB;
        
        tu=0;ts=0;
        timer.Read(tu,ts);
        file<<bestfeas<<" t: "<<tu<<" ";

        //________________________________________________
        //________________________________________________
        //_____________PHASE 2 ELIMINATION EVALUATION
        std::deque<std::deque<int> >candidates;
        
        std::deque<int> bestinon0;
        std::deque<int> bestnon0;
        s->maxlb = bestfeas*100; 
        double best=bestfeas;
        for(int i=0;i<10;++i){
            double best = 1e30;
            tu=0;ts=0;
            timer.Read(tu,ts);
            if(tu>86400) break;
            for(int it=2;it<=10;it+=2){
                double besti = 1e30;
                for(int itt=0;itt<5;itt++){
                    disturb_best(data, Oracle, s, h1, dual, bstfeasnon0, basenon0, non0, UB, it/100.0);
                    if(UB < besti){
                        besti = UB;
                        bestinon0 = non0;
                    }
                }
                if(besti<best){ bestnon0 = bestinon0;  best =besti;}
                if(!bestinon0.empty())candidates.push_back(bestinon0);
                bestinon0.clear();
            }
            tu=0;ts=0;
            timer.Read(tu,ts);
            if(tu>86400) break;
            crossover(data,Oracle, s, h1,dual, candidates, bestnon0, best, 3);
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
            tu=0;ts=0;
            timer.Read(tu,ts);
            if(tu>86400)break;
        }
        
        tu=0;ts=0;
        timer.Read(tu,ts);
        file<<bestfeas<<" t: "<<tu<<std::endl;
     
        file.close();
        dual.clear();
        h1.clear();
        h.clear();
        non0.clear();
        bestinon0.clear();
        basenon0.clear();
        bstfeasnon0.clear();
		
        //________________________________________________
        //________________________________________________
        //_____________  HEURISTIC END   _________________
	
	
	
        if( lvl ) s->SetNDOLog( NULL , 0 );
        
        delete( s );
        delete( MP );      
        
    }  // end( try-block )
    
 
    catch( exception &e ) {
        cerr << e.what() << endl;
        return( 1 );
    }
    catch(...) {
        cerr << "Error: unknown exception thrown" << endl;
        return( 1 );
    }
    
    return( 0 );
    
}  // end( Main )



//========================================================================
//
//==================== AUXILIARY FUNCTIONS ===============================
//
//========================================================================

void
trans_dualsol(const Data & data, const TestFi & Oracle, std::vector<double> & dual, Bundle *s){
    
    cIndex_Set I ; Index D;
    cLMRow Lm = s->ReadBestSol( I , D);
    dual.assign(dual.size(),0);
    for(int k=0; k<data.ndemands; ++k)
        for(int i=0; i<data.nnodes; ++i)
            if(Oracle.Iu[k*data.nnodes + i]>=0)
                dual[Oracle.Iu[k*data.nnodes + i]] = Lm[Oracle.Iu[k*data.nnodes + i]];
    
    
}

//========================================================================


bool crossover(const Data & data, TestFi & Oracle, Bundle * s, const std::vector<double>& h1, std::vector<double> &dual, const std::deque<std::deque<int> >& cand, std::deque<int>& non0, double & best, int itmax){
    
    std::vector<int> ya1(data.narcs,0);
    std::vector<int> ya2(data.narcs,0);
    std::deque<std::string> hist;
    
    s->SetPar( NDOSolver::kMaxItr , 250);
    int sz = cand.size();
    double UB;
    NDOSolver::NDOStatus retval1;
    
    bool to_better = false;
    for(int p=0;p<sz;++p){
        if(cand[p].empty()) continue;
        transp(cand[p],ya2);
        for(int pp=p+1;pp<sz;++pp){
            if(cand[pp].empty()) continue;
            transp(cand[pp],ya1);
            int t=0;
            int delta;
            
            delta = diff<int>(ya1,ya2);
            while(t < pow(2,delta) && (t < itmax)){
                
                Oracle.re_rand_fix(ya1, ya2, h1, hist);
                Oracle.reset( 0, dual, 1);
                retval1 = s->Solve();
                UB = -s->ReadBestFiVal();
                if(UB<best && (retval1!=2)){
                    best = UB;
                    non0=Oracle.non0;
                    to_better=true;
                }
                ++t;
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


void disturb_best(const Data & data, TestFi & Oracle, Bundle * s, const std::vector<double>& h1, std::vector<double>& dual, const std::deque<int>& bestnon0, const std::deque<int>& basenon0, std::deque<int>& non0, double & best, double ratio){
    
    int retval;
    int cont, it;
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
    cont=0;
    it =0;
    while(cont<chngble && it<100){
        arc = rand()%size;
        ++it;
        double r = rand()%101/100.0;
        if(ya1[bestnon0[arc]] && r>=h1[bestnon0[arc]]){
            ya1[bestnon0[arc]] = false;
            out.push_front(bestnon0[arc]);
            ++cont;
            it=0;
        }
    }
    Oracle.non0.clear();
    for(int a=data.narcs;a--;)
        if(ya1[a]) Oracle.non0.push_front(a);
    
    
    //Put the unfixed arcs from the base (do not consider the previously closed ones).
    s->SetPar( NDOSolver::kMaxItr , 100);
    nnfix= Oracle.deque_concat(out, basenon0);
    Oracle.reset(nnfix, dual, 1);
    s->Solve();
    Oracle.transp_h(h, double(s->getNumIter()), nnfix);
    for (int a = nnfix; a--; ){
        arc = Oracle.non0[a];
        if(h[arc]>=0.3){
            heap.push_back(HeapCell(arc,h[arc]));
        }
    }
    
    //Open the "chngble" highest values
    heap.sort(comp());
    Oracle.non0.erase(Oracle.non0.begin(),Oracle.non0.begin()+nnfix);
    cont=0;
    while(cont<chngble && !heap.empty()){
        arc = heap.front().k;
        Oracle.non0.push_front(arc);
        ++cont;
        heap.pop_front();
    }
    
    //Mini local search. try to close some arcs with the Flow problem.
    s->SetPar( NDOSolver::kMaxItr , 250);
    Oracle.reset(0, dual, 1);
    s->Solve();
    best = -s->ReadBestFiVal();
    non0 = Oracle.non0;
    Oracle.non0.clear();
    
    
    ya1.clear();
    out.clear();
    heap.clear();
    
}

//========================================================================


void reopt(const Data & data, TestFi & Oracle, Bundle * s, std::vector<double>& h1, std::vector<double> &dual, std::deque<int>& non0, int nnfxdarc, double f2){
    
    s->SetPar( NDOSolver::kMaxItr , 100);
    
    int arc;
    bool mod = true;
    double f1 = 1/double(s->getNumIter());
    
    Oracle.non0 = non0;
    Oracle.reset(nnfxdarc, dual,0);
    s->Solve();
    Oracle.transp_h(h1, double(s->getNumIter()),nnfxdarc);
    
    f1 = 0.01;
    while(f1 <= 0.05){
        mod = false;
        nnfxdarc=0;
        non0.clear();
        for(int a = Oracle.szunfxd; a--; ){
            int arc = Oracle.non0[a];
            if(h1[arc]>=f2){
                non0.push_back(arc);
                mod = true;
            }
            else if(h1[arc]>=f1){
                non0.push_front(arc);
                ++nnfxdarc;
            }else{
                mod = true;
            }
        }
        if(mod==false){
            f1 +=0.01;
        }else{
            Oracle.insert_in_place(non0, Oracle.szunfxd);
            Oracle.reset(nnfxdarc, dual,0);
            s->Solve();
            Oracle.transp_h(h1, double(s->getNumIter()), nnfxdarc);
        }
    }
    non0.clear();
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
        sol0.set_data(firstnon0, &data);
        sol0.create_model(1);
        retval = sol0.solve(UB,firstnon0, h1,x1,1);
        sol0.clear_model();
        
        if(retval==0){
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



/*--------------------------------------------------------------------------*/
/*-------------------------- End File Main.C -------------------------------*/
/*--------------------------------------------------------------------------*/
