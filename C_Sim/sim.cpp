#include "sim.h"

void SpMV_perf(units* unit, int ptr[DIM], int unroll [SpMVSamples], 
			   float SpMVRU [SpMVSamples],int SpMVLT[SpMVSamples] ){
			
	SpMVconfig(unit, ptr,unroll);
	SpMVCSR_rtrace(ptr,unroll,"oursys",SpMVRU);
	SpMVCSR_ltrace(unit->MAC_latency, ptr,unroll,"oursys",SpMVLT);
	
	for (int i=0;i<SpMVSamples;i++){
		unit->pre_opt_lt += SpMVLT[i];
		//add dfx latency in SpMVCSR_latency for SpMVCSR as well
		unit->pre_opt_ru += SpMVRU[i];
		if (i==SpMVSamples-1) unit->pre_opt_ru /=SpMVSamples;
	}
	cout<<"Pre-optimization- Stage:LT:RU "<< unit->pre_opt_rr<<"::"<<unit->pre_opt_lt<< "::"<<unit->pre_opt_ru<<"%"<<endl;
	if (MAX_STAGES==0)cout<<"No Optimization in effect. Optimization turned off\n";
	init_perf<int,SpMVSamples>(SpMVLT);
	init_perf<float,SpMVSamples>(SpMVRU);

	opt_reconfig_rate(unit,ptr,unroll);
	SpMVCSR_rtrace(ptr,unroll,"oursys",SpMVRU);
	SpMVCSR_ltrace(unit->MAC_latency, ptr,unroll,"oursys",SpMVLT);
	
	//unit->pre_opt_lt = unit->SpMVCSR_latency;

	for (int i=0;i<SpMVSamples;i++){
		unit->SpMVCSR_latency += SpMVLT[i];
		unit->post_opt_lt += SpMVLT[i];
		unit->SpMVCSR_avg_underutil += SpMVRU[i];
		unit->post_opt_ru += SpMVRU[i];
		if (i==SpMVSamples-1) {
			unit->SpMVCSR_avg_underutil /=SpMVSamples; 
			unit->post_opt_ru /=SpMVSamples;
		}
	}
	cout<<"Post-optimization:- Stage:LT:RU "<< unit->post_opt_rr<<"::"<<unit->post_opt_lt<< "::"<<unit->post_opt_ru<<"%"<<endl;
}
//const int MAX_STAGES = 2;

void opt_reconfig_rate(units* unit, int ptr[DIM], int unroll[SpMVSamples]){
	//optimizing algorithm: 
	int sum_=0;
	int STAGES=MAX_STAGES+1;
    float** stage = new float*[STAGES];
    for (int i = 0; i < STAGES; ++i) {
        stage[i] = new float[SpMVSamples];
    }

	float sum =0;
	int uidx=0;

	for (int i=0; i<SpMVSamples; i++){
		stage[0][i] = unroll[i];
	}

	for (int j=1;j<STAGES;j++){

		for (int i=0;i<j;i++){
			stage[j][i]=stage[j-1][i];
		}

		for (int i =j; i<SpMVSamples; i++){
			float frac=0;
				frac = abs((stage[j-1][i]/stage[j-1][i-1])-1);
				if (frac<=reconfig_th)stage[j][i] = stage[j-1][i-1];
				else stage[j][i] = stage[j-1][i];
				// if (stage[j][i] != stage[j][i-1]) sum_++;
				// if (i==SpMVSamples-1){
				// 	cout<<" Reconfig in stage["<<j<<"]: = " << sum_+1<<endl;
				// 	if (j==STAGES-1)unit->reconfig_x = sum_+1;
				// 	sum_=1;
				// } 
		}
	}
	sum_=1;
	for (int i =0;i<SpMVSamples;i++){
		// cout<< i<<"th sample: "<< stage[0][i]<<" ::: "<<stage[1][i]<<" ::: " 
		// 						<<stage[2][i]<<" ::: "<<stage[3][i]<<endl;
		if (i!=0 && stage[STAGES-1][i] != stage[STAGES-1][i-1]) sum_++;
		unroll[i] = stage[STAGES-1][i];
		//cout<< i<<"th sample: "<< unroll[i]<<endl;
	}
	unit->post_opt_rr=sum_;
	//cout<<"Optimized # of Reconfig Events: = " << unit->post_opt_rr <<endl;
}

void SpMVconfig(units* unit, int ptr[DIM], int unroll[SpMVSamples]){
	
	int sum_=0;
	int numVal = DIM/SpMVSamples;
	float sum =0;
	int uidx=0;
    int ptr_=1;
	for (int i=0; i<DIM; i++){
		sum +=ptr[i];
		if (ptr_==numVal){
            if (sum==0)sum=1;
            unroll[uidx] = floor(sum/numVal);
			//stage[0][uidx] = unroll[uidx];
			if (uidx!=0 && unroll[uidx] != unroll[uidx-1]) sum_++;
			uidx++;
			sum=0;
            ptr_ = 0;
		}
        ptr_++;
	}	
	//for (int i=0;i<SpMVSamples;i++)cout<< i<<"th sample: "<< unroll[i]<<endl;
	unit->pre_opt_rr = sum_+1;
	//cout<<"Un-optimized # of Reconfig Events: = " << unit->pre_opt_rr<<endl;
	
}

void init_solver(matrix A, string& solver_, int act_sol[4],units* unit){
	bool dom = diagDom(A,unit);
	bool sym = symmetry(A,unit);
	matProp(dom,sym);
	if (dom==true) {solver_="JB"; act_sol[2]=1;}
	else if (sym==false) {solver_="BiCG"; act_sol[0]=1;}
	else if (sym==true) {solver_="CG"; act_sol[1] =1;}
}


void initialize(matrix* A, float b[DIM], float x[DIM], units unit, metrics *rpt, param *params){
	//Latency of iniitial operations (static Part)
	float Ax[DIM];
	init_ones(b);
	rowlength(A->offsets,A->ptr);
	SpMVCSR(A->values,A->ind,A->ptr,x,Ax);
	vector_arith(b,Ax,params->r,0);	
	vector_copy(params->r,params->p);
	params->rs_old = vector_vector_dot(params->r,params->r);
	vector_copy(params->r,params->r_star);
	params->itr_solver = 0;
	perf_accum("initialize",rpt,unit,params);
}

// void baseline(matrix A, float b[DIM], float u[DIM], float x[DIM], param* params,
// 	          double error[ITR], string solver, methods alg,  units* unit, metrics* rpt){
	
// 	conjugate_gradient(A.values,
// 					   A.ind,
// 					   A.offsets,
// 					   b,
// 					   u,
// 					   x,
// 					   params,
// 					   error)

// 	perf_accum_bs (solver, rpt, *unit, params);


// }

void solver(matrix A, float b[DIM], float u[DIM], float x[DIM], param* params,
			double error[ITR], string solver, methods alg,  units* unit, metrics* rpt){ //call it using solver (solver,alg,&perf)
		
	if (solver=="CG"){
		
		conjugate_gradient(A.values,
						   A.ind,
						   A.offsets,
						   b,
						   u,
						   x,
						   params,
						   error);
		perf_accum (solver, rpt, *unit, params);
		goto ENDIF;
	}
	if (solver=="BiCG"){
		//perf_accum (alg.BiCG,perf);
		bicgstab		  (A.values,
						   A.ind,
						   A.offsets,
						   b,
						   u,
						   x,
						   params,
						   error);
		perf_accum (solver, rpt, *unit, params);
		//cout<<"Iteration: "<<params->itr_round<<endl;
		goto ENDIF;
	}
	if (solver=="JB"){
		jacobi  		  (A.values,
						   A.ind,
						   A.offsets,
						   b,
						   u,
						   x,
						   params,
						   error,
						   unit);
		
		perf_accum (solver, rpt, *unit, params);
		//cout<<"Iteration: "<<params->itr_round<<endl;
		goto ENDIF;
	}
	ENDIF:;
}

void modSolver(string& solver_, string case_, int act_sol[4], string df[4], bool *break_ ) {
	if (case_ == "DIV")
	{	
		for (int i=0; i<4; i++){
			if (act_sol[i] == 0) {
				solver_ = df[i];
				act_sol[i] = 1;
				goto END_MOD;;
			} 
		}
		cout<<"All the solver have been tested\n"<<endl;
		*break_=true;

    }
    else if (case_=="RATE")
    {  
	    if (solver_ == "JB") {
	        solver_ = "GS";
	    } else if (solver_ == "GS") {
	        solver_ = "CG";
	    } else if (solver_ == "CG") {
	        solver_ = "BiCG"; 
	    } 
	    goto END_MOD;
	}

    END_MOD:;
}

/**************************************************************/
//SpMVconfig() function
/*	int stages = 3;
	int numVal = DIM / SpMVSamples;
    float sum = 0;
    int uidx = 0;
    int ptr_ = 1;
    
    // Compute unroll values
    for (int i = 0; i < DIM; i++) {
        sum += ptr[i];
        if (ptr_ == numVal) {
            if (sum == 0) sum = 1;
            unroll[uidx] = ceil(sum / numVal);
            uidx++;
            sum = 0;
            ptr_ = 0;
        }
        ptr_++;
    }

    // Print unroll values
    std::cout << "Unroll values:" << std::endl;
    for (int i = 0; i < SpMVSamples; i++) {
        std::cout << "unroll[" << i << "]: " << unroll[i] << std::endl;
    }

    // Define stage arrays dynamically based on the number of stages
	//int MAX_STAGES = 3;
	float** stagesArr = new float*[MAX_STAGES];
	for (int i = 0; i < MAX_STAGES; ++i) {
		stagesArr[i] = nullptr;
	}

    //float* stagesArr[MAX_STAGES] = {nullptr};
    for (int i = 0; i < stages; i++) {
        stagesArr[i] = new float[SpMVSamples];
    }
	int sum_;
    // Loop through stages
    for (int s = 0; s < stages; s++) {
        float* stage = stagesArr[s];
        stage[0] = unroll[0];
        sum_ = 0;

        // Compute loop for each stage
        for (int i = 1; i < SpMVSamples; i++) {
            float frac = abs((unroll[i] / unroll[i - 1]) - 1);
            if (frac <= reconfig_th) {
                stage[i] = unroll[i - 1];
            } else {
                stage[i] = unroll[i];
            }
            if (stage[i] != stage[i - 1]) sum_++;
        }

        cout << "Reconfig in stage " << s + 1 << " = " << sum_ << endl;

        // Update unroll values for the next stage
        for (int i = 0; i < SpMVSamples; i++) {
            unroll[i] = stage[i];
        }
    }

    // Clean up dynamic memory
    for (int i = 0; i < stages; i++) {
        delete[] stagesArr[i];
    }

    // Update reconfig_x value in the units struct
    unit->reconfig_x = sum_;
*/
