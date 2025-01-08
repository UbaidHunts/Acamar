#include "sim.h"


/*Acamar: A Dynamically Reconfigurable Scientific Computing Accelerator for Robust Convergence and Minimal Resource Utilization*/
string dataset;
int NNZ;
int MAX_STAGES;
float reconfig_th;
int SpMV_unroll_factor_baseline;


//int SpMVSamples;

int main(int argc, char* argv[]){
	if (SpMVSamples>DIM || SpMVSamples == 0){
		cout<<"Illegal SpMV Sampling Rate\nRe-run using correct input"<<endl;
	}
	if (argc!=8){
		cout<<"\nInput Aruments not given properly"<<endl; 
		cout<<"<opt_config_file>.txt <baseline_config_file>.txt <dataset_name> <NNZ> <rOpt_stages> <rOpt_th> <SpMV_unroll_factor_baseline>\n"<<endl;
		cout<<"Dataset name must be present in the dataset directory. Must be given the name only. For Example; wang3 "<<endl;
		return 0;
	}
	
	for (int i=1;i<argc;i++){
		// if (i==1) readConfigFile(argv[i], A.ptr, unit);
		// if (i==2) readConfigFile(argv[i], A.ptr, baseline);
		if (i==3) dataset= argv[3];
		if (i==4) NNZ= atoi(argv[4]);
		if (i==5) MAX_STAGES = atoi (argv[5]);
		if (i==6) reconfig_th = atof(argv[6]);;
		if (i==7) SpMV_unroll_factor_baseline = atoi (argv[7]);
	}

	if (MAX_STAGES>SpMVSamples){
		cout<<"ERROR: Reconfiguration Optimization Stages must be less than " <<SpMVSamples<<endl;
		return 0;
	}
	float freq = 100*10e6;
	
	//cout<<"\n******************* RUNNING "<< dataset <<" *******************"<<endl;
	matrix A;
	matrix A_csc;
	// Static Arrays
	// float b[DIM];
    // float x[DIM];
    // float u[DIM];
    // float error[ITR];
	
	// Dynamic Arrays
	float *b = new float [DIM];
	float *x = new float [DIM];
	float *u = new float [DIM];
	double *error = new double [ITR];
	int* unroll = new int[SpMVSamples];
	int SpMVLT[SpMVSamples];
	float SpMVRU[SpMVSamples];
	init_perf<int,SpMVSamples>(unroll);
	init_perf<int,SpMVSamples>(SpMVLT);
	init_perf<float,SpMVSamples>(SpMVRU);
	bool break_=false;
	int act_sol[]={0,0,0,1};
	methods alg;
	metrics rpt;
	param params;
	units unit,baseline;
	vector<string> iterator;
	bool converged=false;
	int size=0;
	string df[4] = {"BiCG","CG","JB","GS"};

	/****************:-:INITIALIZING VALUES:-:****************/
	init_metrics(&alg);
	init_ones(x);
	init_ones(u);
	
	string f2=dataset;
	string fo = "dataset/offsets";
	string fc = "dataset/colIdx";
	string fv = "dataset/values";
	string f3 = "_ex4k.txt";
	string file_off= fo+f2+f3;
	string file_col= fc+f2+f3;
	string file_val= fv+f2+f3;
		
	read_off(file_off, A.offsets, &size);
    read_ind(file_col, A.ind, &size);
    read_values(file_val, A.values, &size);
	rowlength(A.offsets,A.ptr);

	/****************:-:READING UNIT LATENCY FROM CONFIG FILES:-:****************/
	readConfigFile(argv[1], A.ptr, unit);
	
	SpMV_perf(&unit, A.ptr,unroll,SpMVRU,SpMVLT);
	int sum=0;
	//for (int i=0;i<SpMVSamples;i++)sum=sum+unroll[i];
	//cout<<"Average Unroll Resources: " << sum/double(SpMVSamples)<<endl;
	//printf("\nSpMVCSR PERFORMANCE FOR SAMPLING RATE %d\n", SpMVSamples);
	//printf("SpMVCSR latency = %d\n",unit.SpMVCSR_latency);
	//printf("SpMVCSR Avg. Under RU = %f%% (lower is better)\n\n",unit.SpMVCSR_avg_underutil);
	//vector_print<int,SpMVSamples>(unroll);
	readConfigFile(argv[2], A.ptr, baseline);
	for (int i =0;i<SpMVSamples;i++) unroll[i]=SpMV_unroll_factor_baseline;
	//cout<<"Baseline SpMV Latency = "<<baseline.SpMVCSR_latency <<endl;
	
	//throughput_acamar(unroll,unit);
	// throughput_baseline(unroll,baseline);

	//return 0;


	/****************:-:RUNNING OUR SYSTEM:-:****************/
 
   	string solver_;//="JB";
	init_solver(A,solver_,act_sol, &unit);
	iterator.push_back(solver_);
	int it=0;
	int swh=0;
	// int reset=0;
	initialize(&A, b, x, unit, &rpt, &params);
	while (!converged){
		for (int i=0; i<ITR; i++) error[i]=0;
		solver(A, b, u, x, &params, error, solver_, alg, &unit, &rpt); //u = input, x=output
		//cout<<params.itr<<" Error[" <<params.e_ptr<<"]: "<<error[params.e_ptr]<<endl;
		
		if (params.c){
			// cout<<"\nConverged Achieved in *"<<params.itr<<  
			// 	  "* itr with error = "<<error[params.e_ptr]<<endl;
			converged=true;	
		}
		else {
			if (params.d==1){
				modSolver(solver_, "DIV",act_sol,df,&break_);
				if (break_){
					//cout<<"Solution Diverged in " << params.itr<<"th itr"<<"\nExiting";
					goto BASELINES;
				}
				params.d=0;
				swh=1;
				it++;
				iterator.push_back(solver_);
			}
			if (swh){
				cout<<"SWITHCED to "<<solver_<<endl;
				init_ones(x);
				initialize(&A, b, x, unit, &rpt, &params);
				perf_accum (solver_, &rpt, unit, &params);
				swh=0;
		//		reset=0;
			}
			vector_copy(x,u);
		} 
		//it++;
		//reset=reset+ITR;
	}
	switch_hist(iterator,it);
	print_perf(rpt,"ACAMAR Performance");
	cout<< "Latency to Calculate T matrix in JB "<<unit.split_latency<<endl;
	cout<< "Latency to check symmetry: "<<unit.sym_latecny<<endl;
	cout<< "Latency to check diagonal dominance: "<<unit.diaDom_latency<<endl<<endl;


	/****************:-:RUNNING BASELINE:-:****************/
	
	BASELINES:;
	init_perf<int,SpMVSamples>(unroll);
	init_perf<int,SpMVSamples>(SpMVLT);
	init_perf<float,SpMVSamples>(SpMVRU);

    readConfigFile(argv[2], A.ptr, baseline);
   	string solver__[]={"JB","CG","BiCG"};
	param paramsb;
	cout<<"********************** BASELINES **********************"<<endl;

	SpMVCSR_ltrace(baseline.MAC_latency, A.ptr, unroll,"baseline",SpMVLT);
	for (int i=0;i<SpMVSamples;i++) baseline.SpMVCSR_latency += SpMVLT[i];
	cout<<"Baseline SpMV Latency = "<<baseline.SpMVCSR_latency<<endl;
	throughput_baseline(unroll,baseline);

	SpMVCSR_rtrace(A.ptr, unroll, "baseline", SpMVRU);
	for (int i=0;i<SpMVSamples;i++) baseline.SpMVCSR_avg_underutil += SpMVRU[i];
	baseline.SpMVCSR_avg_underutil /= SpMVSamples;
	cout<<"Baseline SpMV Underutilization = "<<baseline.SpMVCSR_avg_underutil<<"% (lower is better)"<<endl;

	for (int j=0;j<3;j++){
		cout<<"\n___________________________Running Baseline Varient of "<<solver__[j]<<"___________________________"<<endl;
		//string solver__=solver_[j];
		rpt.latency = 0;
		rpt.SpMV =0;
		paramsb.reset();
		init_ones(u);
		init_ones(x);
		converged=false;
		initialize(&A, b, x, baseline, &rpt, &paramsb);
		while (!converged){
			for (int i=0; i<ITR; i++) error[i]=0;
			
			solver(A, b, u, x, &paramsb, error, solver__[j], alg, &baseline, &rpt); //u = input, x=output
			if (paramsb.c){
				//cout<<"after solver"<<paramsb.c<<endl;
				converged=true;	
			}
			else {
				if (paramsb.d==1){
					break;
				}
				vector_copy(x,u);
			} 
		}
    	printf("Latency of (initialize, solver, modSolver): %f ",rpt.latency);
    	cout<<"("<<rpt.latency<<") ="<< rpt.latency * 1/freq <<"ns"<<endl;
		printf("Latency portion for SpMV: %f ",rpt.SpMV);
		cout<<"("<<rpt.SpMV<<") ="<< rpt.SpMV * 1/freq <<"ns"<<endl;

		if(solver__[j]=="JB") rpt.latency+=baseline.split_latency;// cout<< "Latency to Calculate T matrix in JB: "<<baseline.split_latency<<endl;

		
	}

	delete [] b;
	delete [] x;
	delete [] u;
	delete [] error;
	cout<<"_______________________________________________________________________________________________________________________________________\n";
	
}

/**************************************************************/