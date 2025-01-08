using namespace std;

struct metrics
{
	double latency=0;
	double SpMV=0;
	int mulOps=0;
	int addOps=0;
	int memX=0;
	int pre_opt_ru=0;
	int pre_opt_lt=0;

	int BRAM=0;
	int DSP=0;
	int FF=0;
	int LUT=0;
	int URAM=0;
	int num_metrics=6;
};

struct methods
{
	metrics CG;
	metrics BiCG;
	metrics JB;
	metrics GS;
};


struct units
{
	int rowlength_latency = 0;
	int vector_arith_latency = 0;
	int vector_copy_latency = 0;
	int vector_vector_dot_latency = 0;
	int norm_latency = 0;
	int scalar_vector_latency = 0;
	int init_ones_latency=0;
	int div_latency=0;
	//int SpMVCSR_unroll=3;
	int DFX_latency = 0;
	int MAC_latency = 0;
	int diaDom_unit_latency = 5;
	//int reconfig_x = 0;
	int diaDom_latency =0;
    int diaDom_unroll_factor = 16;
	int sym_latecny =0;
	int sym_unroll_factor=8;
	//optimizer performance
	float pre_opt_ru =0;
	int pre_opt_lt =0;
	int pre_opt_rr =0;
	//split function
	int split_latency=0;
	int split_unroll_factor=16;

	float post_opt_ru=0;
	int post_opt_lt=0;
	int post_opt_rr=0;

/*Will calculate*/
	int SpMVCSR_latency=0;		
	float SpMVCSR_avg_underutil=0;
	
};