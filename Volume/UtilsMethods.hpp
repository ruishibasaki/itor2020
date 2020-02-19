//
//  UtilsMethods.hpp
//  
//
//  Created by Rui Shibasaki on 30/11/2018.
//
//

#ifndef UtilsMethods_hpp
#define UtilsMethods_hpp

#include <stdio.h>
#include "structures.hpp"
#include <ilcplex/ilocplex.h>

void deque_to_string(const std::deque<int>& d, std::string & s, int sz);
void transp(const std::deque<std::deque<int> > & cand,std::vector<int> & y, int pos);
void transp(const std::deque<int> & cand,std::vector<int> & y);
int deque_union_front(std::deque<int> & fst,const std::deque<int> & snd, int sz);
int deque_union_back(std::deque<int> & fst,const std::deque<int> & snd, int sz);


void remove_i(int p, std::list<int> &N);
Pair1 arg_min(const std::list<int> & nodes, const std::vector<double> & costs);
double dijkstra(int source, int sink, const std::vector<int> &grid, std::vector<Pair1> & preced, const std::vector<double>& cij);


void MCND_read_data(std::string fname, Data & data);
int final_feas(const Data* data, const std::vector<int>& idx, const IloNumArray & xsol, const std::vector<double> & y1, std::vector<int> & y_);




#endif /* UtilsMethods_hpp */
