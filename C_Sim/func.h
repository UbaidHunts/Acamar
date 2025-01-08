#include "datastr.h"

#define KB 1024
#define DIM 4*KB
#define ROW_A DIM
#define COL_A DIM
#define ITR 100
#define SAMPLES 10
#define SETUP 2*ITR 
#define th 0.00001

extern int NNZ;
extern float reconfig_th;
extern int MAX_STAGES;
extern int SpMV_unroll_factor_baseline;
//extern int SpMVSamples;

struct matrix{
	float* values=new float [NNZ];
	int* ind=  new int [NNZ];
	int* offsets = new int [DIM];
	int* ptr = new int [DIM];
};

struct param{
	float p[DIM];
	float r[DIM];
	float r_star[DIM];
	float rs_old=0;
	int c=0;
	int d=0;
	int itr=0;
	int e_ptr=0;
	int itr_solver=0;
	int itr_round=0;

	void reset(){
		for (int i=0;i<DIM;i++){
			p[i] = 0;
			r[i] =0;
			r_star[i]=0;
		}
		rs_old=0;
		c=0;
		d=0;
		itr=0;
		e_ptr=0;
		itr_solver=0;
		itr_round=0;
	}


};

void SpMVCSR(
		float *values, 
		int *ind, 
		int ptr[DIM], 
		float vec[DIM],
		float res[DIM]);

template <typename T>
void init_ones(
		T res[DIM]);

void init_zeros(
		float a[ITR]);

void rowlength(
		int offsets[DIM], 
		int ptr[DIM]);

float vector_vector_dot(
		float a[DIM],
		float b[DIM]);

void scalar_vector(
		float a,
		float b[DIM],
		float res[DIM]);

void vector_arith(
		float a[DIM],
		float b[DIM],
		float res[DIM],
		int op);

void vector_copy(
		float a[DIM],
		float res[DIM]);

double norm(
		float a[DIM]);

void vector_print_file(
	 	float x[DIM], 
	 	string file_name);

void split_mult_dia_rand(float* values, int* ind, int offset[DIM],
                        float* T_values, int* T_ind, int T_ptr[DIM], 
						float D_values[DIM], int D_ind[DIM], int D_ptr[DIM], units* unit);

void multiply_diag_rand(float* R_values, int* R_ind, int R_ptr[DIM],
                         float D_values[DIM], int D_ind[DIM], int D_ptr[DIM],
                          float* T_values, int* T_ind, int T_ptr[DIM], int *latency, units* unit);

void conjugate_gradient(
		float* values,
		int* ind,
		int offsets[DIM],
		float b[DIM],
		float u[DIM],
		float x[DIM],
		param* params,
		double error[ITR]);

void bicgstab(
		float* values,
		int* ind,
		int offsets[DIM],
		float b[DIM],
		float u[DIM],
		float x[DIM],
		param* params,
		double error[ITR]);

void jacobi(
		float* values,
		int* ind,
		int offsets[DIM],
		float b[DIM],
		float u[DIM],
		float x[DIM],
		param* params,
		double error[ITR],
		units* unit);

void top(
		float *values,
		int *ind,
		int offsets[DIM],
		float b[DIM],
		float x[DIM]);

/**************************************************************/

