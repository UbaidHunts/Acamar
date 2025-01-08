#include "sim.h"

void SpMV(
		float A[DIM][DIM],
		float v[DIM],
		float Av[DIM]){
//***********************************************
//#pragma HLS UNROLL
	for (int i = 0; i < ROW_A; i++) {
			Av[i]=0;
	        for (int j = 0; j < COL_A; j++){
	            Av[i] = Av[i]+ A[i][j] * v[j];
	        }
	    }
}

void SpMVCSR(
		float *values, 
		int *ind, 
		int ptr[DIM], 
		float vec[DIM],
		float res[DIM]){

	int valIdx=0;
	int colIdx=0;
	int k;
	for (int i=0;i<DIM;i++){
		res[i]=0;
		for(int j=0;j<DIM;j++){
			if (j==ptr[i])break;
			int a =valIdx+j;
			int b =colIdx+j;
			float mul=values[a]*vec[ind[b]];
			res[i] = res[i]+mul;
		}
		valIdx+=ptr[i];
		colIdx+=ptr[i];
	}
}

void rowlength(int offsets[DIM], int ptr[DIM])
{	
	ptr[0]=offsets[0];
	for (int i=1;i<DIM;i++){
		ptr[i]= offsets[i]-offsets[i-1];
		}
}

float vector_vector_dot(
		float a[DIM],
		float b[DIM]){
//***********************************************
	float sum=0;
	for (int i=0; i<DIM; i++){
		sum = sum + a[i]*b[i];
	}
	return sum;
}

void scalar_vector(
		float a,
		float b[DIM],
		float res[DIM]){
//***********************************************
	for (int i=0; i<DIM; i++){
		res[i] = a*b[i];
	}
}

void vector_arith(
		float a[DIM],
		float b[DIM],
		float res[DIM],
		int op){

	//op=1 implies add operation
	//op=0 implies sub operation

	if (op==1){
		ADD: for (int i=0; i<DIM;i++){
			res[i]=a[i]+b[i];
		}
		return;
	}
	else if(op==0){
		SUB: for (int i=0; i<DIM;i++){
			 res[i]=a[i]-b[i];
		}
		return;
	}
	printf("op parameter in vector_arith() is not given properly");
	return;

}

void multiply_diag_rand(float* R_values, int* R_ind, int R_ptr[DIM],
                         float D_values[DIM], int D_ind[DIM], int D_ptr[DIM],
                          float* T_values, int* T_ind, int T_ptr[DIM], int *latency, units* unit){
    int temp = 0;
	int mulops=0;
	int k;
    for (int i = 0; i <DIM; i++){
        T_ptr[i] = 0;
		k=0;
        if(R_ptr[i] && D_ptr[i]){
            for (int j = temp; j < temp +R_ptr[i]; j++){
                T_ptr[i] ++;
                T_ind[j] = R_ind[j];
                T_values[j] = D_values[i] * R_values[j];	
				if (j-temp  >= unit->split_unroll_factor * k) {
					*latency = *latency + 3;
					k++;
				}
            }
            temp = temp+ R_ptr[i];
        }
    }
}

void split_mult_dia_rand(float* values, int* ind, int offset[DIM],
                        float* T_values, int* T_ind, int T_ptr[DIM], 
						float D_values[DIM], int D_ind[DIM], int D_ptr[DIM], units* unit){

    int ptr[DIM];
    int diag = 0;
    float R_values[NNZ]; int R_ind[NNZ]; int R_ptr[DIM];
    rowlength(offset,ptr);
    int temp = 0;
    int R_pnt = 0;
    int D_pnt = 0;
	int k;
	int latency=0;
    for (int i = 0; i <DIM; i++){
		k=0;
        D_ptr[i] = 0;
        R_ptr[i] = 0;
        if(ptr[i] > 0){
            for(int j = temp; j < ptr[i]+temp; j++){
                if(i == ind[j]){
                    D_values[D_pnt] = 1/ values[j];
                    D_ind[D_pnt] = ind[j];
                    D_pnt++;
                    D_ptr[i]++;
					if (j-temp  >= unit->split_unroll_factor * k) {
						latency = latency + 11; //div_latency = 10. 1 for add;
						k++;
					}
				} 
                else{
                    R_values[R_pnt] = values[j];
                    R_ind[R_pnt] = ind[j];
                    R_pnt++;
                    R_ptr[i]++; 
					if (j-temp  >= unit->split_unroll_factor * k) {
						latency = latency + 1; //div_latency = 10;
						k++;
					}
                }
            }
            temp = temp + ptr[i];
        }
    }
    multiply_diag_rand(R_values, R_ind, R_ptr, D_values, D_ind, D_ptr,
                         T_values, T_ind, T_ptr, &latency, unit);
	
	unit->split_latency = latency;
}



double norm(
		float a[DIM]){
//***********************************************
	double error=0;
	for (int i=0;i<DIM;i++){
		error = error + a[i]*a[i];
		//cout<< i<<"Error in Function: "<<error<<endl;
	}
	//printf("Error in Func= %f\n",error);
	return sqrt(error);
}

void vector_copy(
		float a[DIM],
		float res[DIM]){
//***********************************************
	for (int i=0;i<DIM;i++){
		res[i]=a[i];
	}
}

void conjugate_gradient(
		float* values,
		int* ind,
		int offsets[DIM],
		float b[DIM],
		float u[DIM],
		float x[DIM],
		param* params,
		double error[ITR]){
//***********************************************
	//float p[DIM];
	//float x[DIM];
	//float Ax[DIM];
	//float r[DIM];
	//float r_star[DIM];
	//float rs_old;
	float Ap[DIM];
	float pAp;
	float alpha_p[DIM];
	float alpha_Ap[DIM];
	float beta_p[DIM];
	int   ptr[DIM];

	float alpha;
	float rs_new;
	float beta;
	rowlength(offsets,ptr);
/*
	//row lengths
	rowlength(offsets,ptr);
	
	//x=u.copy()
	//init_ones(x);
	vector_copy(u,x);
	
	//r=b-A.dot(x)
	//SpMV(A,x,Ax);
	SpMVCSR(values, ind, ptr, x, Ax);
	vector_arith(b,Ax,r,0);
	
	//p=r.copy()
	vector_copy(r,p);
	
	//rs_old=r.dot(r)
	rs_old=vector_vector_dot(r,r);
*/
	params->itr_round=0;
	for(int i=0; i<ITR; i++){
		params->itr_round++;
		params->itr=params->itr+1;
		params->itr_solver=params->itr_solver+1;
		params->e_ptr=i;
		//Ap=A.dot(p)
		//SpMV(A,p,Ap);
						
		SpMVCSR(values, ind, ptr, params->p, Ap);
		//alpha=rs_old/p.dot(Ap)
		pAp=vector_vector_dot(params->p,Ap);
		alpha = params->rs_old/pAp;

		//x=x+alpha*p;
		scalar_vector(alpha,params->p,alpha_p);
		vector_arith(x,alpha_p,x,1);


		//r=r-alpha*Ap
		scalar_vector(alpha,Ap,alpha_Ap);
		vector_arith(params->r,alpha_Ap,params->r,0);

		//rs_new=r.dot(r)
		rs_new=vector_vector_dot(params->r,params->r);

		//error = norm(r)
		error[i] = norm(params->r); 
		if (error[i]<th){ 
			params->c=1;
			cout << "CG Converged in *"<< params->itr<< "th* itr with error: "<<error[i]<<endl;
			//params->e_ptr=i;
			break;
		}
		if (isnan(error[i])||isinf(error[i])){
            cout << "CG Diverged in "<< params->itr<< "th itr with error: "<<error[i]<<endl;
			params->d=1;
            break;			
		}
		if (params->itr_solver > SETUP && error[i] > 10e4){
            cout << "CG Diverged in "<< params->itr<< "th itr with error: "<<error[i]<<endl;
			params->d=1;
            break;
		}
		if (params->itr > 100000){
			cout << "CG Diverged in "<< params->itr<< "th due to bouncing values"<<endl;
			params->d=1;
			break;
		}
		//p = r+beta*p
		beta=rs_new/params->rs_old;
		scalar_vector(beta,params->p,beta_p);
		vector_arith(params->r,beta_p,params->p,1);

		params->rs_old=rs_new;
	}
	//vector_copy(x,u);

}

void bicgstab(
		float* values,
		int* ind,
		int offsets[DIM],
		float b[DIM],
		float u[DIM],
		float x[DIM],
		param* params,
		double error[ITR]){
//***********************************************
	//float p[DIM];
	//float r_star[DIM];
	//float x[DIM];
	//float Ax[DIM];
	//float r[DIM];
	//float rs_old;
	float Ap[DIM];
	float As[DIM];
	float r_starAp;
	float alpha_p[DIM];
	float alpha_Ap[DIM];
	float w_As[DIM];
	float beta_p[DIM];
	float s[DIM];
	float ws[DIM];
	float wAp[DIM];
	float temp[DIM];
	int   ptr[DIM];
	float alpha;
	float rs_new;
	float beta;
	float w;

	rowlength(offsets,ptr);
	/*//x=u.copy()
	//init_ones(x);
	vector_copy(u,x);

	//r=b-A.dot(x)
	//SpMV(A,x,Ax);
	SpMVCSR(values, ind, ptr, x, Ax);
	
	vector_arith(b,Ax,r,0);
	
	//p=r.copy()
	vector_copy(r,p);
	rs_old=vector_vector_dot(r,r);
	vector_copy(r, r_star);
	*/
	params->itr_round=0;
	for(int i=0; i<ITR; i++){
		params->itr_round++;
		params->itr=params->itr+1;
		params->e_ptr=i;
		params->itr_solver = params->itr_solver+1;
		//rs_old=r.dot(r)

		//SpMV(A,p,Ap);
		SpMVCSR(values, ind, ptr, params->p, Ap);
		
		//alpha=rs_old/r_star.dot(Ap)
		r_starAp=vector_vector_dot(params->r_star,Ap);
		alpha = params->rs_old/r_starAp;
		
		//s = r - alpha*Ap;
		scalar_vector(alpha,Ap,alpha_Ap);
		vector_arith(params->r,alpha_Ap,s,0);

		//As = A.dot(s)
		SpMVCSR(values, ind, ptr, s, As);

		//w = As.dot(s)/As.dot(As)
		w = vector_vector_dot(As, s)/vector_vector_dot(As, As);

		//x=x+alha*p
		scalar_vector(alpha,params->p,alpha_p);
		vector_arith(x,alpha_p,x,1);

		//x=x+w*s 
		scalar_vector(w, s, ws);
		vector_arith(x,ws,x,1);

		//r=s - w*As
		scalar_vector(w,As,w_As);
		vector_arith(s,w_As,params->r,0);

		//rs_new=r.dot(r_star)
		rs_new=vector_vector_dot(params->r,params->r_star);

		//error = norm(r)
		error[i] = norm(params->r);

		if (error[i]<th){ 
			cout << "BiCG Converged in *"<< params->itr<< "th* itr with error: "<<error[i]<<endl;
			params->c=1;
			//params->e_ptr=i;	
			break;
		}
		if (isnan(error[i])||isinf(error[i])){
			cout << "BiCG Diverged in "<< params->itr<< "th itr with error: "<<error[i]<<endl;
			params->d=1;
            break;
		}
		if (params->itr_solver > SETUP && error[i] > 10e4){
            cout << "BiCG Diverged in "<< params->itr<< "th itr with error: "<<error[i]<<endl;
			params->d=1;
            break;
		}
		if (params->itr > 100000){
			cout << "BiCG Diverged in "<< params->itr<< "th due to bouncing values"<<endl;
			params->d=1;
			break;
		}
		
		//beta = (rs_new/rs_old)* (alpha/w)
		beta=(rs_new/params->rs_old) * (alpha/w);
		
		//wAp = w*Ap
		scalar_vector(w,Ap,wAp);

		//beta_p = beta*(p-w*Ap)
		vector_arith(params->p,wAp,temp,0);
		scalar_vector(beta,temp,beta_p);


		//p = r+ beta_p	
		vector_arith(params->r,beta_p,params->p,1);
	
		//rs_old=rs_new
		params->rs_old=rs_new;
	}
	//vector_copy(x,u);
	//printf(" Printing after for loop: \n");
//vector_print(x);
//vector_print_file(x, "bicgstab_out2.txt");
}

void jacobi(
		float* values,
		int* ind,
		int offsets[DIM],
		float b[DIM],
		float u[DIM],
		float x[DIM],
		param* params,
		double error[ITR],
		units* unit){

	float x_new[DIM];
	float temp[DIM];
	float res[DIM];
	float norm_inp[DIM];
	float T_values[NNZ]; 
	int T_ind[NNZ]; 
	int T_ptr[DIM];
	float D_values[NNZ]; 
	int D_ind[NNZ]; 
	int D_ptr[DIM];
	params->itr_round=0;

	split_mult_dia_rand(values, ind, offsets,
	 T_values, T_ind, T_ptr, D_values, D_ind, D_ptr, unit);
	vector_copy(u, x);
	for(int i = 0; i < ITR; i++){
		params->itr_round++;
		params->itr=params->itr+1;
		params->itr_solver=params->itr_solver+1;
		params->e_ptr=i;	
		//rs_old=r.dot(r)
		//x_new = D_inv.b - T.x
		SpMVCSR(D_values, D_ind, D_ptr, b, res);
		SpMVCSR(T_values, T_ind, T_ptr, x, temp);
		vector_arith(res, temp, x_new, 0);
		vector_arith(x_new, x, norm_inp, 0);
		vector_copy(x_new, x);
		double er = norm(norm_inp);
		error[i] = er;
		if (error[i]<th){ 
			cout << "JB Converged in *"<< params->itr<< "th* itr with error: "<<error[i]<<endl;
			//cout<<"IN THE JACOBI: "<<params->itr_round<<endl;
			params->c=1;
				break;
		}
		if (isnan(error[i])||isinf(error[i])){
            cout << "JB Diverged in "<< params->itr<< "th itr with error: "<<error[i]<<endl;
			params->d=1;
            break;			
		}
		if (params->itr_solver > SETUP && error[i] > 10e4){
            cout << "JB Diverged in "<< params->itr<< "th itr with error: "<<error[i]<<endl;
			params->d=1;
            break;
		}
		if (params->itr > 100000){
			cout << "JB Diverged in "<< params->itr<< "th due to bouncing values"<<endl;
			params->d=1;
			break;
		}
	}
}

/**************************************************************/