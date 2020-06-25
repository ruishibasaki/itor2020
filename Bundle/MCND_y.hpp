#ifndef _FlowYBCP_H
#define _FlowYBCP_H

#include <ilcplex/ilocplex.h>

#include "structures.hpp"
#include <vector>
#include <deque>


class FlowY{

public:
	IloEnv* env;
	IloCplex* cplex;
	
	IloModel* model;
    IloObjective Obj;

	IloNumVarArray x;

	std::vector<int> idx;
	int nx;
	
    IloRangeArray flow_row;
	IloRangeArray capa_row;
	
	std::vector<int> idf_row;
	std::vector<int> idc_row;
	int nfr, ncr;
	
	std::vector<int> yint;
	
	const Data * data;
	int nnodes;
	int ndemands;
	int narcs;
	int sznz;
	bool FlowY_init;
	
	
	std::vector<int> grid;
	
	FlowY();
	
	void set_data(const std::deque<int>& non0, const Data * d);
    int solve(double &solution, std::deque<int>& non0, const std::vector<double> & y1, std::vector<double> & x1, int phase);
	void create_model(int phase);
	void clear_model();
    int price(std::deque<HeapCell>& cols_to_add,const IloNumArray & pi_,const IloNumArray & alpha_ );
    void add_cols(std::deque<HeapCell>& cols_to_add);

    ~FlowY();
};

#endif 

