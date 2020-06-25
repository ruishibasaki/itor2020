#ifndef STRUCT_H
#define STRUCT_H

#include <vector>
#include <deque>
#include <list>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>


#define DBLPRES 1e-8
  

struct Pair1{
	int node, pos;
	inline Pair1 &operator=(const Pair1 & source){
		node = source.node;
		pos = source.pos;
		return *this;
	}
	inline Pair1(int n,int p):node(n), pos(p) {};
	inline Pair1():node(0),pos(-1){};
	inline ~Pair1(){};
};


struct Pair2{
	int fst, snd;
	inline Pair2 &operator=(const Pair2 & source){
		fst = source.fst;
		snd = source.snd;
		return *this;
	}
	inline Pair2(int n,int p):fst(n), snd(p) {};
	inline Pair2():fst(-1),snd(-1){};
	inline ~Pair2(){};
};


struct Arc1{
	int i, j; 
	double capa;
	double f;
	std::vector<double> c;
	std::vector<double> b;
	inline Arc1(int ii, int jj, double cap, double ff, int nk, double k):
	i(ii), j(jj), capa(cap), f(ff){
		c.resize(nk,k);
	}
	inline Arc1(int nk, double k){c.resize(nk,k);}
	inline Arc1(int nk){c.resize(nk,0.0);}
	inline Arc1(){}
	inline ~Arc1(){c.clear();b.clear();}
};

struct Demand{
	int O, D;
	double quantity;
};

struct HeapCell{
	int k;
	double rc_;
	
	inline HeapCell(){};
    inline HeapCell(int dk, double rc): k(dk), rc_(rc){}
	inline HeapCell(const HeapCell & copy){
		this->k = copy.k;
		this->rc_ = copy.rc_;
	}
	
	inline HeapCell& operator=(HeapCell other){
		this->k= other.k;
		this->rc_= other.rc_;
		return *this;
	}
	inline ~HeapCell(){};
};

class comp{

public:	
	bool operator()(const HeapCell & x, const HeapCell & y)const{
        return (x.rc_ > y.rc_);
	}
};


class Data{
public:
	std::vector<Arc1> arcs;
	std::vector<Demand> d_k;
	int ndemands, narcs, nnodes;
	inline ~Data(){
		d_k.clear();
		arcs.clear();
	}
	
};

template <typename T>
int diff(const std::vector<T>& y1, const std::vector<T>& y2){
	int sz1 = y1.size(); 
	int sz2 = y2.size();
	if(sz1!=sz2) return -1;
	
	int c=0;
	for(int a=sz1;a--;){
		if(y1[a]!=y2[a]) ++c;
	}
	return c;
}



#endif 

