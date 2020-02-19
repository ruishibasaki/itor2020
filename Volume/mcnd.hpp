#ifndef __PMP_HPP__
#define __PMP_HPP__

#include "VolVolume.hpp"
#include "structures.hpp"
#include <cstdlib>

#include <sys/times.h>




class MCND : public VOL_user_hooks {
public:
	//volume hooks
   int compute_rc(const VOL_dvector& u, VOL_dvector& rc);

   int solve_subproblem(const VOL_dvector& u, const VOL_dvector& rc,
			double& lcost, VOL_dvector& x, VOL_dvector&v,
			double& pcost);

   int heuristics(const VOL_dvector& p, 
		  const VOL_dvector& x, double& heur_val);

   double knapsack(int a,  const VOL_dvector& rc, VOL_dvector& x);

   //heuristic
    //(re)initiate
   int reset_topology();
   void reset(VOL_problem & volp, int unfixd, const std::vector<double> & hs, bool hotstart=true);
   
    //getters translaters
   void trans_dualsol(std::vector<double> & dual, const VOL_dvector &dsol);
   void transp_h(std::vector<double> & y1, double div, int place=0);

    //modifiers
    int deque_concat(const std::deque<int> & out,const std::deque<int> & from);
    void insert_in_place(const std::deque<int> & in, int nplace);

    //randomers
    double re_rand_fix(const std::vector<int> & ya2, const std::vector<int> & ya1,  const std::vector<double> & y1, std::deque<std::string> & hist);
   
   
public: 
  const Data *data;
  int ndemands, narcs, nnodes;
  int szopnd, szunfxd, sznon0;
  std::deque<int> non0;
  std::vector<int> h1;
  
  //double bestlb;
  std::list<HeapCell> heap;  
  std::vector<int> Iu; //index of lag multipliers
  
public:
  MCND(const Data  *d);
  virtual ~MCND(){
	 Iu.clear();
	 non0.clear();
	 heap.clear();
	 //yl.clear();
  }   
};


#endif
