#include "sim.h"

void switch_hist(vector <string> iterator, int size){
	cout<<"\n___________________________ Switching activity of iterative solvers ___________________________ "<<endl;
	for (int i=0;i<iterator.size();i++){
		cout<< i <<"."<< iterator[i] <<endl;
	}	
}

void matProp(bool dom,bool sym){
    if(dom)cout<<"\n___________________________ Strictly Diagonally Dominant Matrix ___________________________"<<endl;
    if(!dom)cout<<"___________________________ Not Strictly Diagonally Dominant Matrix ___________________________"<<endl;
    if(sym)cout<<"___________________________ Symmetric Matrix ___________________________\n"<<endl;
    if(!sym)cout<<"__________________________ Non-symmetric Matrix ___________________________\n"<<endl;
}

void scale(metrics *perf, int chunks){
	perf->latency=perf->latency*chunks;
}

void print_perf(metrics perf, string s){
    float freq = 100*10e6;
    float clock = 1/freq;
	cout<<"___________________________ "<< s<<" ___________________________ "<<endl;
    printf("Latency of (initialize, solver, modSolver): %f ",perf.latency);
    cout<<"("<<perf.latency<<") ="<< perf.latency * clock <<"ns"<<endl;
	printf("Latency portion for SpMV: %f ",perf.SpMV);
	cout<<"("<<perf.SpMV<<") ="<< perf.SpMV * clock <<"ns"<<endl;
}

void init_ones(float a[DIM]){
	for (int i =0; i<DIM; i++){
		a[i]=1;
	}
}

template <typename T>
void init_zeros(T a[ITR]){
	for (int i =0; i<ITR; i++){
		a[i]=0;
	}
}

void read_values(const string& filepath,  float* values, int *size) {

    ifstream file(filepath);
    if (!file.is_open()) {
        cout << "Error opening file." << endl;
        return;
    }
    
    float value;
    int i=0;
   
    while (file >> value) {
        values[i]=value;
        i++;
    }
    *size=i;
    // Close the file
    file.close();
}
void read_ind (const string& filepath,  int* values, int *size) {

    ifstream file(filepath);
    if (!file.is_open()) {
        cout << "Error opening file." << endl;
        return;
    }

    int value;
    int i=0;

    while (file >> value) {
        values[i]=value;
        i++;
    }
    *size=i;
    // Close the file
    file.close();
}
void read_off (const string& filepath,  int values[DIM], int *size) {

    ifstream file(filepath);
    if (!file.is_open()) {
        cout << "Error opening file." << endl;
        return;
    }

    int value;
    int i=0;

    while (file >> value) {
        values[i]=value;
        i++;
    }
    *size=i;
    // Close the file
    file.close();
}

void init_metrics(methods *alg){
	alg->CG.latency = 1;
	alg->CG.BRAM = 1;
	alg->CG.DSP = 1;
	alg->CG.FF = 1;
	alg->CG.LUT = 1;
	alg->CG.URAM =1;

	alg->BiCG.latency = 1;
	alg->BiCG.BRAM = 1;
	alg->BiCG.DSP = 1;
	alg->BiCG.FF = 1;
	alg->BiCG.LUT = 1;
	alg->BiCG.URAM =1;

	alg->JB.latency = 1;
	alg->JB.BRAM = 1;
	alg->JB.DSP = 1;
	alg->JB.FF = 1;
	alg->JB.LUT = 1;
	alg->JB.URAM =1;

	alg->GS.latency = 1;
	alg->GS.BRAM = 1;
	alg->GS.DSP = 1;
	alg->GS.FF = 1;
	alg->GS.LUT = 1;
	alg->GS.URAM =1;
}

bool diagDom(matrix A, units* unit) {
	int j=0;
	int row_start =0;
    int row_end =A.offsets[0];
	int k;
    for (int i = 0; i < DIM; i++) {
		k=0;
        float diag_element = 0;
        float off_diag_sum = 0;
       	for (j = row_start; j < row_end; j++){
            if (A.ind[j] == i) {
                diag_element = A.values[j];
            } else {
                off_diag_sum += abs(A.values[j]);
            }
			//CALCULATING THE LATENCY 
			if (j-row_start >= unit->diaDom_unroll_factor * k) {
				unit->diaDom_latency = unit->diaDom_latency + unit->diaDom_unit_latency;
				k++;
			}	
        }	
        row_start = row_end;
        if (i<DIM-1)row_end = A.offsets[i+1];
        if (off_diag_sum >= abs(diag_element)){
			//cout<<"non-diagDom...BREAKING!!! in "<<i<<"th row"<<endl;
			return false;
		} 
    }
    return true;
}

void csr_tocsc(matrix A, matrix *B, units* unit){

    int *Ap = new int[DIM+1];
    int *Bp = new int[DIM+1];
    Ap[0]=0; Bp[0]=0;
	int k;
	//cout<< "Latency split 0: "<<unit->sym_latecny<<endl;
    for (int i=0; i<DIM;i++){ //Equals the Dim, iteration latency=2, interval =1. trip count change with unroll factor. Reduce by the unroll factor
        //#pragma HLS unroll factor =16
		Ap[i+1] = A.offsets[i];
        Bp[i+1] =0;
    }
	unit->sym_latecny += (DIM/unit->sym_unroll_factor);
	//cout<< "Latency split 1: "<<unit->sym_latecny<<endl;	

    int n_col = DIM;
    int n_row = n_col;

    for (int i=0; i<NNZ; i++){
        Bp[A.ind[i]]++; //Iteration Interval = 3. Do the Unroll optimization on your own
    }
	unit->sym_latecny += ceil((3*NNZ/unit->sym_unroll_factor));
	//cout<< "Latency split 2: "<<unit->sym_latecny<<endl;

    for(int col = 0, cumsum = 0; col < DIM; col++){ //Equals the Dim, iteration latency=2, interval =1. trip count change with unroll factor. Reduce by the unroll factor
        int temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = NNZ;
	unit->sym_latecny += (DIM/unit->sym_unroll_factor);
	//cout<< "Latency split 3: "<<unit->sym_latecny<<endl;
    for(int row = 0; row < n_row; row++){
		k=0;
        for(int jj = Ap[row]; jj < Ap[row+1]; jj++){//Iteration Interval = 3. Do the Unroll optimization on your own
            int col  = A.ind[jj];
            int dest = Bp[col];
            B->ind[dest] = row;
            B->values[dest] = A.values[jj];
            Bp[col]++;
			//cout<<" jj-Ap[row]:: "<<jj-Ap[row]<<endl;
			if (jj-Ap[row] >= unit->sym_unroll_factor * k) {
				unit->sym_latecny += 3;
				k++;
			}

        }
    }
	//cout<< "Latency split 4: "<<unit->sym_latecny<<endl;
    for(int col = 0, last = 0; col <= n_col; col++){ //Equals the Dim, iteration latency=2, interval =1. trip count change with unroll factor. Reduce by the unroll factor
        int temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }  
	unit->sym_latecny += (DIM/unit->sym_unroll_factor);
	//cout<< "Latency split 5: "<<unit->sym_latecny<<endl;
    for (int i=1; i<DIM+1;i++){
        B->offsets[i-1] = Bp[i];
    }
	unit->sym_latecny += (DIM/unit->sym_unroll_factor);
	//cout<< "Latency split 6: "<<unit->sym_latecny<<endl;
    delete [] Ap;
    delete [] Bp;
}

bool symmetry(matrix A, units* unit){
	
	matrix A_csc;
	csr_tocsc(A,&A_csc, unit);

	for (int i=0;i<DIM;i++){//consider II=1 with full unroll. can break in between
		if (A.offsets[i] != A_csc.offsets[i]) {
			if(i > unit->sym_unroll_factor) unit->sym_latecny += ceil((i/unit->sym_unroll_factor));
			else unit->sym_latecny ++;
			return false;
		}
	}
	//cout<<"latency of Symmetry 0 : "<<unit->sym_latecny<<endl; 
	unit->sym_latecny += ceil(1*DIM/unit->sym_unroll_factor);
	//cout<<"latency of Symmetry 1 : "<<unit->sym_latecny<<endl; 

	for (int i=0;i<NNZ;i++){// iteration latency =5. can use the simulator for unroll 
		if (A.ind[i] != A_csc.ind[i]){
			if(i > unit->sym_unroll_factor) unit->sym_latecny += ceil((5*i/unit->sym_unroll_factor));
			else unit->sym_latecny +=5;			
			return false;
		} 
		if (A.values[i] != A_csc.values[i]){
			if(i > unit->sym_unroll_factor) unit->sym_latecny += ceil((5*i/unit->sym_unroll_factor));
			else unit->sym_latecny +=5;
			return false;
		}	 
		
	}
	unit->sym_latecny += ceil(5*NNZ/unit->sym_unroll_factor);
	//cout<<"latency of Symmetry 2 : "<<unit->sym_latecny<<endl;
	return true;
}

void readConfigFile(const string& filename, int ptr[DIM], units& unit) {
    ifstream configFile(filename);
    if (!configFile) {
        cerr << "Error: Failed to open the configuration file." <<endl;
        return ;
    }

    string line;
    while (getline(configFile, line)) {
        istringstream iss(line);
        string key;
        int value;
        if (iss >> key >> value) {
                 if (key == "rowlength_latency") unit.rowlength_latency = value;
            else if (key == "vector_arith_latency") unit.vector_arith_latency = value;
            else if (key == "vector_copy_latency") unit.vector_copy_latency = value;
            else if (key == "vector_vector_dot_latency") unit.vector_vector_dot_latency = value;
            else if (key == "norm_latency") unit.norm_latency = value;
            else if (key == "scalar_vector_latency") unit.scalar_vector_latency = value;
            else if (key == "init_ones_latency") unit.init_ones_latency = value;
            else if (key == "div_latency") unit.div_latency = value;
            else if (key == "MAC_latency") unit.MAC_latency = value;
            else if (key == "DFX_latency") unit.DFX_latency = value;
            else if (key == "NNZ") continue; // Ignore NNZ
            else {
                cerr << "Warning: Unknown key '" << key << "' in the configuration file." << std::endl;
            }
        }
    }
    return ;
}

void perf_accum(string s, metrics* rpt, units unit, param* params){

	if (s=="initialize"){
	rpt->latency += (unit.init_ones_latency+
					unit.rowlength_latency+ 
					unit.SpMVCSR_latency+
					unit.vector_arith_latency+
					unit.vector_copy_latency+
					unit.vector_vector_dot_latency+
					unit.vector_copy_latency);
	rpt->SpMV += unit.SpMVCSR_latency;				
	}
	else if (s=="init_solver"){
		//Exported Seperately
	}

	else if (s=="CG"){
	rpt->latency =  rpt->latency +
					((unit.SpMVCSR_latency + 
					unit.vector_vector_dot_latency+
					unit.div_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.vector_vector_dot_latency+
					unit.norm_latency+
					unit.div_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency)*(params->itr_round));

	rpt->SpMV += (unit.SpMVCSR_latency * params->itr_round);	

	}
	else if (s=="BiCG"){
	rpt->latency += ((unit.SpMVCSR_latency + 
					unit.vector_vector_dot_latency+
					unit.vector_vector_dot_latency+
					unit.div_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.SpMVCSR_latency + 
					unit.vector_vector_dot_latency+// ignoring because can be parallelized unit.vector_vector_dot_latnecy+
					unit.div_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.vector_vector_dot_latency+
					unit.norm_latency+
					unit.div_latency+
					unit.div_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency+
					unit.scalar_vector_latency+
					unit.vector_arith_latency)* (params->itr_round));
	rpt->SpMV += ((2*unit.SpMVCSR_latency)*(params->itr_round));

	}
	else if (s=="JB"){
	rpt->latency += ((//unit.split_latency+
					unit.SpMVCSR_latency +
					unit.SpMVCSR_latency + 
					unit.vector_arith_latency +
					unit.vector_arith_latency +
					unit.vector_copy_latency +
					unit.norm_latency
					)* (params->itr_round));
	rpt->SpMV += ((2*unit.SpMVCSR_latency)*(params->itr_round));

	}
	else if (s=="modSolver"){
		rpt->latency+=4;
	}
}



void SpMVCSR_ltrace(int MAC_latency, int ptr[DIM], int unroll [SpMVSamples], string s, int SpMVLT[SpMVSamples]) {
	int k;
	int numVal = DIM/SpMVSamples;
	int sum=0;
	int uidx=0;
	int ptr_=1;
	
	if (s!="baseline"){
		for (int i = 0; i < DIM; i++) {
			k = 0;
			for (int j = 0; j < DIM; j++) {
				if (j == ptr[i])break;
				if (j >= unroll[uidx] * k) {
					SpMVLT[uidx] = SpMVLT[uidx] + MAC_latency;
					k++;
				}
			}
			if (ptr_==numVal){
				uidx++;
				ptr_=0;
			}
			ptr_++;
		}
	}
	else {
		uidx=0;
		ptr_=1;
		for (int i=0;i<SpMVSamples;i++)SpMVLT[i]=0;
		for (int i = 0; i < DIM; i++) {
			k = 0;
			for (int j = 0; j < DIM; j++) {
				if (j == ptr[i])break;
				if (j >= SpMV_unroll_factor_baseline* k) {
					SpMVLT[uidx] = SpMVLT[uidx] + MAC_latency;
					k++;
				}
			}
			if (ptr_==numVal){
				uidx++;
				ptr_=0;
			}
			ptr_++;
		}
	}
}

void SpMVCSR_rtrace(int ptr[DIM], int unroll [SpMVSamples], string s, float SpMVRU [SpMVSamples])  {
	
	int numVal=DIM/SpMVSamples;
	float sum=0;
	int ptr_=1;
	int uidx=0;
	if (s!="baseline"){
		for (int i = 0; i < DIM; i++) {
			float util = (ptr[i]%unroll[uidx])*100/unroll[uidx];
			if(util==0)util=100;
			sum += 100-(util);
			if (ptr_==numVal){
				SpMVRU[uidx] = sum / numVal;
				//cout<<"unroll: "<< unroll[uidx]<<endl;
				sum=0;
				ptr_ = 0;
				uidx++;
			}
			ptr_++;
		}
	}
	else {
		sum=0;
		ptr_=1;
		uidx=0;
		for (int i = 0; i < DIM; i++) {
			float util = (ptr[i]%SpMV_unroll_factor_baseline)*100/SpMV_unroll_factor_baseline;
			if(util==0)util=100;
			sum += 100-(util);
			if (ptr_==numVal){
				SpMVRU[uidx] = sum / numVal;
				//cout<<"unroll: "<< SpMV_unroll_factor_baseline<<endl;
				sum=0;
				ptr_ = 0;
				uidx++;
			}
			ptr_++;
		}

	}
	
	//vector_print<float, SpMVSamples>(SpMVRU);
}
//Not called in the host file. Used the cc and NNZ in the excel to find out the peak throughput. 
void throughput_acamar(int unroll[SpMVSamples], units unit){
	double prod=0;
	for (int i=0;i<SpMVSamples;i++){
		prod = prod + (2*unroll[i]/double(unit.MAC_latency));
		//cout<<"Product: "<<prod<<endl;
	}
	cout<<"unit.SpMVCSR_latency: "<<unit.SpMVCSR_latency<<endl;
	float th_flops = prod/double(SpMVSamples);
	float ac_flops = 2*NNZ/double(unit.SpMVCSR_latency);
	float pct	   =ac_flops*100/th_flops;
	cout<<"Average Peak Theoretical FLOPS/cc = "<< th_flops<<endl;
	cout<<"Average Peak Acheived    FLOPS/cc = "<< ac_flops<<endl;
	cout<<"% Throughput Acheived 		  = "<<   pct   <<endl;
}

void throughput_baseline(int* unroll, units baseline){
	double prod=0;
	for (int i=0;i<SpMVSamples;i++){
		prod = prod + (2*SpMV_unroll_factor_baseline/double(6));
		//cout<<"Product: "<<prod<<endl;
	}
	//cout<<"Product: "<<prod<<endl;
	//cout<<"baseline.SpMVCSR_latency: "<<baseline.SpMVCSR_latency<<endl;
	//cout<<baseline.SpMVCSR_latency;
	float th_flops = prod/double(SpMVSamples);
	float ac_flops = 2*NNZ/double(baseline.SpMVCSR_latency);
	float pct	   =ac_flops*100/th_flops;
	cout<<"Average Peak Theoretical FLOPS/cc = "<< th_flops<<endl;
	//cout<<"Average Peak Acheived    FLOPS/cc = "<< ac_flops<<endl;
	//cout<< ac_flops;
	//cout<<"% Throughput Acheived 		  = "<<   pct   <<endl;
}

/**************************************************************/
