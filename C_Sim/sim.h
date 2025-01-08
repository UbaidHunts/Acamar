#ifndef SIM_H
#define SIM_H

#include <stdio.h>
#include <iostream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <vector> 
#include <string>
#include <cstring>
#include <unordered_map>
#include <vector>
#include "func.h"
using namespace std;

#define DFX_LATENCY 10
#define SpMVSamples 32
//#define diaDom_unroll_factor 16

void init_config();//read config files
void switch_hist(vector <string> iterator, int size);
void matProp(bool dom,bool sym);
void scale(metrics *perf, int chunks);
void print_perf(metrics perf, string s);
void init_ones(float a[DIM]);
void init_metrics(methods *alg);

void read_values(const string& filepath,  float* values, int *size);
void read_ind (const string& filepath,  int* values, int *size);
void read_off (const string& filepath,  int* values, int *size);

void perf_accum(string s, metrics* rpt, units unit, param* params);

void initialize(matrix* A, float b[DIM], float x[DIM], units unit, metrics *rpt, param *params);
void modSolver(string& solver_, string case_, int act_sol[4], string df[4], bool* break_ );
void solver(matrix A, float b[DIM], float u[DIM], float x[DIM], param* params,
			double error[ITR], string solver, methods alg,  units* unit, metrics* rpt);

void analyzer(string &solver_, methods alg, float error [ITR], param params, int *swh, metrics *perf);
bool diagDom(matrix A, units* unit);
void csr_tocsc(matrix A, matrix *B, units* unit);
bool symmetry(matrix A, units* unit);
void init_solver(matrix A, string& solver_, int act_sol[4], units* unit);

void opt_reconfig_rate(units* unit, int ptr[DIM], int unroll[SpMVSamples]);

template <typename T, int siz>
void vector_print(T x[siz]){
    for (int i=0;i<siz;i++){
    	cout<<"x["<<i<<"]= " <<x[i]<<endl;
    }
}

template <typename T, int siz>
void init_perf(T x[siz]){
    for (int i=0;i<siz;i++){
    	x[i] = 0;
    }
}

void readConfigFile(const string& filename, int ptr[DIM], units& unit) ;
void SpMVconfig(units* unit, int ptr[DIM], int unroll[SpMVSamples]);
void SpMVCSR_rtrace(int ptr[DIM], int unroll [SpMVSamples], string s, float SpMVRU [SpMVSamples]);
void SpMVCSR_ltrace(int MAC_latency, int ptr[DIM], int unroll [SpMVSamples], string s, int SpMVLT[SpMVSamples]);
void SpMV_perf(units* unit, int ptr[DIM], int unroll [SpMVSamples], 
			   float SpMVRU [SpMVSamples],int SpMVLT[SpMVSamples] );
void throughput_acamar(int unroll[SpMVSamples], units unit);

/*BASELINE*/

void perf_accum_bs(string s, metrics* rpt, units baseline, param* params);
void throughput_baseline(int unroll[SpMVSamples], units unit);

#endif // SIM_H
/**************************************************************/