//compute factorial
function print_vec(vec_whatever)
{
	for(i=0; i<vec_whatever.length; i++)
	{
		console.log(i+":\t"+vec_whatever[i])
	}
	return null;
}
function print_mat(mat_whatever)
{
	if(mat_whatever && mat_whatever.length)
	for(i=0; i<mat_whatever.length; i++)
	{
		console.log(i+":\t"+mat_whatever[i]);
	}
	return null;
}

function factrial(num)
{
    var rval=1;
    for (var i = 2; i <= num; i++)	rval = rval * i;
	return rval;
}

function get_mat(num_i, num_j)
{
	var matrix=new Array();
	for(var i=0; i<num_i; i++)
	{
		var vec_buf= new Array(num_j);
		for(var j=0; j<num_j; j++)
		{
			vec_buf[j]=0.0;
		}
		matrix[i]=vec_buf;
	}
	return matrix;
}

function get_mat3d(num_i, num_j, num_k)
{
	var matrix=new Array();
	for(var i=0; i<num_i; i++)
	{
		var matrix_buf=get_mat(num_j, num_k);
		matrix[i]=matrix_buf;
	}
	return matrix;
}

function get_mat_copy(mat_whatever)
	{
		var i_num=mat_whatever.length;
		var j_num=mat_whatever[0].length;
		var copy_mat=get_mat(i_num, j_num);
		for(i=0; i<i_num; i++)
		{
			for(j=0; j<j_num; j++)
			{
				copy_mat[i][j]=mat_whatever[i][j];
			}
		}
		
		return copy_mat;
	}

function get_vec_Xi(vec_spkt, bin_width)
{
	var spkt_dura=vec_spkt[vec_spkt.length-1];
	var bin_num=Math.ceil(spkt_dura/bin_width);
	var vec_Xi= new Array();
	for(i=0; i<bin_num; i++)
	{
			vec_Xi[i]=0;
	}
	for(i=0; i<vec_spkt.length; i++)
	{
		var bin_id=Math.floor(vec_spkt[i]/bin_width);
		if(bin_id<bin_num) vec_Xi[bin_id]++;
	}
	return vec_Xi;
}

//func to get emisson probs.
function get_mat_emission(vec_Xi, vec_lambda)
{
    var mat_emission=new Array();
    for(n=0; n<vec_Xi.length; n++)
    {
        var vec_emission=new Array();
        for(i=0; i<vec_lambda.length; i++)
        {
            var pdf=Math.pow(vec_lambda[i], vec_Xi[n])*Math.exp(-1.0*vec_lambda[i])/factrial(vec_Xi[n])
			vec_emission.push(pdf);
        }
        mat_emission.push(vec_emission);
    }
    return mat_emission;
}	

//func to get alpha and C		
function get_alpha_C(mat_A, vec_pi, mat_emission)
{
	var num_of_states=vec_pi.length;
	var num_of_obs=mat_emission.length;

	var mat_alpha_buf=new Array();
	var vec_C_buf=new Array();

	//n=0
	/////////////////////////////////////////
	var alpha_0=new Array();
	for(var i=0; i<num_of_states; i++)
	{
		var alpha_0_i=mat_emission[0][i]*vec_pi[i];
		alpha_0.push(alpha_0_i);
	}

	var C_0=0.0;
	for(var i=0; i<num_of_states; i++)
	{
		C_0+=alpha_0[i];
	}
	vec_C_buf.push(C_0);

	for(var i=0; i<num_of_states; i++)
	{
		alpha_0[i]=(alpha_0[i]/C_0);
	}
	mat_alpha_buf.push(alpha_0);
	/////////////////////////////////////////

	//n>0
	/////////////////////////////////////////
	for(var n=1; n<num_of_obs; n++)
	{
		var alpha_n= new Array();
		for(var i=0; i<num_of_states; i++)
		{
			var sum_j=0.0;
			for(var j=0; j<num_of_states; j++)
			{
				sum_j+=mat_alpha_buf[n-1][j]*mat_A[j][i];
			}
			var alpha_n_i=mat_emission[n][i]*sum_j;
			alpha_n.push(alpha_n_i);
		}

		var C_n=0.0;
		for(var i=0; i<num_of_states; i++)
		{
			C_n+=alpha_n[i];
		}
		vec_C_buf.push(C_n);

		for(var i=0; i<num_of_states; i++)
		{
			alpha_n[i]=(alpha_n[i]/C_n);
		}
		mat_alpha_buf.push(alpha_n);
	}
	/////////////////////////////////////////

	var res=new Array(2);
	res[0]=vec_C_buf;
	res[1]=mat_alpha_buf;
	
	return res;
}

function get_beta(mat_A, vec_pi, mat_emission, vec_C)
{
    var num_of_states=vec_pi.length;
    var num_of_obs=mat_emission.length;
	
	//initialize
	var mat_beta_buf=get_mat(num_of_obs, num_of_states);
	
	//n=N-1
    /////////////////////////////////////////
	for(var i=0; i<num_of_states; i++) mat_beta_buf[num_of_obs-1][i]=1.0;
 	/////////////////////////////////////////

    //n<N-1
    /////////////////////////////////////////
    for(var m=1; m<=num_of_obs-1; m++)
    {
        var n=(num_of_obs-1-m);
        for(var i=0; i<num_of_states; i++)
        {
            var sum_j=0.0;
            for(var j=0; j<num_of_states; j++)
            {
                sum_j+=mat_beta_buf[n+1][j]*mat_emission[n+1][j]*mat_A[i][j];
            }
            mat_beta_buf[n][i]=(sum_j/vec_C[n+1]);
        }
    }
    /////////////////////////////////////////
	
    return mat_beta_buf;
}
	
function get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C)
{
    var num_of_states=mat_emission[0].length;
    var num_of_obs=mat_emission.length;
    
	//gamma matrix
     ////////////////////////////////////////////////////////////////
     var mat_Gamma_buf=get_mat(num_of_obs, num_of_states);
     for(var n=0; n<num_of_obs; n++)
     {
         for(var i=0; i<num_of_states; i++)
         {
                mat_Gamma_buf[n][i]=mat_alpha[n][i]*mat_beta[n][i];
         }
     }
     ////////////////////////////////////////////////////////////////

     //Xi matrix
     ////////////////////////////////////////////////////////////////
     var mat_Xi_buf=get_mat3d(num_of_obs-1, num_of_states, num_of_states);

     for(var m=0; m<num_of_obs-1; m++)
     {
         for(var i=0; i<num_of_states; i++)
         {
             for(var j=0; j<num_of_states; j++)
             {
                 mat_Xi_buf[m][i][j]=(mat_alpha[m][i]*mat_emission[m+1][j]*mat_A[i][j]*mat_beta[m+1][j])/vec_C[m+1];
             }
         }
     }
     //////////////////////////////////////////////////////////////

	var res=new Array(2);
	res[0]=mat_Gamma_buf;
	res[1]=mat_Xi_buf;
	
	return res;

}
	
function HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi)
{
    var mat_emission=get_mat_emission(vec_Xi,vec_lambda);

    var res_a_C=get_alpha_C(mat_A, vec_pi, mat_emission);
	var vec_C=res_a_C[0];
	var mat_alpha=res_a_C[1];

    var mat_beta=get_beta(mat_A, vec_pi, mat_emission, vec_C);
	
    var res_gamma_xi=get_Gamma_Xi(mat_A, mat_emission, mat_alpha, mat_beta, vec_C);
	var mat_Gamma=res_gamma_xi[0];
	var mat_Xi=res_gamma_xi[1];
	
	var res=new Array(2);
	res[0]=mat_Gamma;
	res[1]=mat_Xi;
	
	return res;
}

function HMM_M_step(vec_Xi,  mat_A,  vec_lambda,  vec_pi,  mat_Gamma,  mat_Xi)
{
    var num_of_states=mat_A.length;
    var num_of_obs=vec_Xi.length;

    //maximize wrt pi vector
    /////////////////////////////////////////////////////
    var pi_denom=0.0;
    for(var k=0; k<num_of_states; k++)
    {
        pi_denom+=mat_Gamma[0][k];
    }
    var vec_pi_new=new Array(num_of_states);
    for(var k=0; k<num_of_states; k++)
    {
        vec_pi_new[k]=mat_Gamma[0][k]/pi_denom;
    }
    /////////////////////////////////////////////////////

    //maximize wrt lambda vector
    /////////////////////////////////////////////////////
    var vec_lambda_new=new Array(num_of_states);
    for(var k=0; k<num_of_states; k++)
    {
        var lambda_denom=0.0;
        for(var n=0; n<num_of_obs; n++)
        {
            lambda_denom+=mat_Gamma[n][k];
        }

        var lambda_nume=0.0;
        for(var n=0; n<num_of_obs; n++)
        {
            lambda_nume+=(mat_Gamma[n][k]*vec_Xi[n]);
        }

        if(lambda_nume==0.0) vec_lambda_new[k]=0.0;
        else vec_lambda_new[k]=lambda_nume/lambda_denom;
    }
    /////////////////////////////////////////////////////

    //maximize wrt A matirx
    /////////////////////////////////////////////////////
    var mat_A_new=get_mat(num_of_states, num_of_states);
    for(var j=0; j<num_of_states; j++)
    {
        var A_denome=0.0;
        for(var n=0; n<num_of_obs-1; n++)
        {
            for(var l=0; l<num_of_states; l++)
            {
                A_denome+=mat_Xi[n][j][l];
            }
        }
        for(var k=0; k<num_of_states; k++)
        {
            var A_nume=0.0;
            for(var n=0; n<num_of_obs-1; n++)
            {
                A_nume+=mat_Xi[n][j][k];
            }
            if(A_nume==0.0) mat_A_new[j][k]=0.0;
            else mat_A_new[j][k]=A_nume/A_denome;
        }
    }
    /////////////////////////////////////////////////////
	res=new Array(3);
	res[0]=vec_pi_new;
    res[1]=vec_lambda_new;
    res[2]=mat_A_new;
	
	return res;
}

	
function HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi)
{
	var mat_emission=get_mat_emission(vec_Xi, vec_lambda);

    var num_of_states=mat_A.length;
    var num_of_obs=vec_Xi.length;

    var mat_hs_seq=get_mat(num_of_states, num_of_obs);
    var vec_logp_seq=new Array(num_of_states);
	
	//n=0
    for(var j=0; j<num_of_states; j++)
    {
        mat_hs_seq[j][0]=j;
        vec_logp_seq[j]=Math.log(vec_pi[j]*mat_emission[0][j])/Math.LN10;
    }
	
	//n>0
	for(var n=1; n<num_of_obs; n++)
	{
		//copy the seq. up to n-1
		/////////////////////////////////////
		var mat_hs_seq_buf=get_mat(num_of_states, num_of_obs);
		for(var k=0; k<num_of_states; k++)
		{
			for(var m=0; m<num_of_obs; m++)
			{
				mat_hs_seq_buf[k][m]=mat_hs_seq[k][m];
			}
		}
		var vec_logp_seq_buf=vec_logp_seq.concat();
		/////////////////////////////////////
		
		// nth node->j
		for(var j=0; j<num_of_states; j++)
		{
			// n-1th node->j
							//compute logp for i->j trans
				var vec_h_logprob_i=new Array(num_of_states);
				for(var i=0; i<num_of_states; i++)
				{
						vec_h_logprob_i[i]=vec_logp_seq[i]+Math.log(mat_emission[n][j]*mat_A[i][j])/Math.LN10;
				}
				//get max logp
				var max_element=Math.max.apply(null, vec_h_logprob_i);
				var max_pos;
				for(var sea=0; sea<vec_h_logprob_i.length; sea++)
				{
					if(vec_h_logprob_i[sea]==max_element) max_pos=sea;
				}
				vec_logp_seq_buf[j]=max_element;
				for(var m=0; m<num_of_obs; m++)
				{
					mat_hs_seq_buf[j][m]=mat_hs_seq[max_pos][m];
				}
				mat_hs_seq_buf[j][n]=j;
			
		}
		//updata the seq. up to n-1
		/////////////////////////////////////
		for(var k=0; k<num_of_states; k++)
		{
			for(var m=0; m<num_of_obs; m++)
			{
				mat_hs_seq[k][m]=mat_hs_seq_buf[k][m];
			}
		}
		vec_logp_seq=vec_logp_seq_buf.concat();
		/////////////////////////////////////
	}
	
	var max_element=Math.max.apply(null, vec_logp_seq);
	var max_pos;
	for(var sea=0; sea<vec_logp_seq.length; sea++)
	{
		if(vec_logp_seq[sea]==max_element) max_pos=sea;
	}
	
	var vec_hs_seq=new Array(num_of_obs);
	for(var n=0; n<num_of_obs; n++)
	{
		vec_hs_seq[n]=mat_hs_seq[max_pos][n];
	}

	return vec_hs_seq;
}	

		
function get_hmm_ratefunc(spike_time, bin_width)
{
    var EMloop_num=5000;		//number of EM itteration
	var mat_A=new Array();
	var mat_A0=[0.999,0.001];	mat_A.push(mat_A0);
	var mat_A1=[0.001,0.999];  mat_A.push(mat_A1);
	var vec_pi=[0.5,0.5];
	var vec_lambda=new Array(vec_pi.length);
	var mean_rate=spike_time.length/(spike_time[spike_time.length-1]-spike_time[0]);
	vec_lambda[0]=(mean_rate*0.75)*bin_width;
	vec_lambda[1]=(mean_rate*1.25)*bin_width;
	
	
	//2D array stores
	//0: begining time of each bins in second
	//1: rate of each bin
	
    //get hmm rate func
    //////////////////////////////////////////////////////////////////////
    var vec_spkt=new Array(spike_time.length)
	for(var i=0; i<spike_time.length; i++) vec_spkt[i]=spike_time[i]-spike_time[0];

    vec_Xi=get_vec_Xi(vec_spkt, bin_width);
    
    var E_res1=HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi);
	var mat_Gamma=E_res1[0];
	var mat_Xi=E_res1[1];
	
    var mat_A_old=get_mat_copy(mat_A);
    var vec_pi_old=vec_pi.concat();
    var vec_lambda_old=vec_lambda.concat();
    var loop=0;
	var flag=0;
	while(loop<=EMloop_num || flag==0)
    {
		var M_res=HMM_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi);
		vec_pi_new=M_res[0];
        vec_lambda_new=M_res[1];
        mat_A_new=M_res[2];
		
		vec_pi=vec_pi_new.concat();
        vec_lambda=vec_lambda_new.concat();
		mat_A=get_mat_copy(mat_A_new);
		
        var sum_check=0.0;
        var num_state=vec_pi.length;
		
        for(var i=0; i<num_state; i++)
        {
            for(var j=0; j<num_state; j++)
            {
                sum_check+=Math.abs(mat_A_old[i][j]-mat_A[i][j]);
            }
            sum_check+=Math.abs(vec_pi_old[i]-vec_pi[i]);
            sum_check+=Math.abs(vec_lambda_old[i]-vec_lambda[i]);
        }
        if(sum_check/(1.0*num_state*(num_state+2))<1.0e-7) flag++;

        mat_A_old=get_mat_copy(mat_A);
        vec_pi_old=vec_pi.concat();
        vec_lambda_old=vec_lambda.concat();

        E_res=HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi);
		mat_Gamma=E_res[0];
		mat_Xi=E_res[1];
		
		loop++;
    }
    
    var vec_hidden=HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi);
    //////////////////////////////////////////////////////////////////////

    var rate_func=get_mat(vec_Xi.length,2);
    
	var c_time=0.0;
    for(var n=0; n<vec_Xi.length; n++)
    {
        var state_id=vec_hidden[n];
        rate_func[n][0]=Math.round(c_time*100)/100.0;
        rate_func[n][1]=Math.round(vec_lambda[state_id]*100)/(bin_width*100.0);
        c_time+=bin_width;
    }
	
    return rate_func;
}
/*
 * 旧関数
function Estimate_rate_hmm(spike_time, bin_width)
{
    var EMloop_num=5000;		//number of EM itteration
	var mat_A=new Array();
	var mat_A0=[0.999,0.001];	mat_A.push(mat_A0);
	var mat_A1=[0.001,0.999];  mat_A.push(mat_A1);
	var vec_pi=[0.5,0.5];
	var vec_lambda=new Array(vec_pi.length);
	var mean_rate=spike_time.length/(spike_time[spike_time.length-1]-spike_time[0]);
	vec_lambda[0]=(mean_rate*0.75)*bin_width;
	vec_lambda[1]=(mean_rate*1.25)*bin_width;
	
	
	//2D array stores
	//0: begining time of each bins in second
	//1: rate of each bin
	
    //get hmm rate func
    //////////////////////////////////////////////////////////////////////
    var vec_spkt=new Array(spike_time.length)
	for(var i=0; i<spike_time.length; i++) vec_spkt[i]=spike_time[i]-spike_time[0];

    vec_Xi=get_vec_Xi(vec_spkt, bin_width);
    
    var E_res1=HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi);
	var mat_Gamma=E_res1[0];
	var mat_Xi=E_res1[1];
	
    var mat_A_old=get_mat_copy(mat_A);
    var vec_pi_old=vec_pi.concat();
    var vec_lambda_old=vec_lambda.concat();
    var loop=0;
	var flag=0;
	while(loop<=EMloop_num || flag==0)
    {
		var M_res=HMM_M_step(vec_Xi, mat_A, vec_lambda, vec_pi, mat_Gamma, mat_Xi);
		vec_pi_new=M_res[0];
        vec_lambda_new=M_res[1];
        mat_A_new=M_res[2];
		
		vec_pi=vec_pi_new.concat();
        vec_lambda=vec_lambda_new.concat();
		mat_A=get_mat_copy(mat_A_new);
		
        var sum_check=0.0;
        var num_state=vec_pi.length;
		
        for(var i=0; i<num_state; i++)
        {
            for(var j=0; j<num_state; j++)
            {
                sum_check+=Math.abs(mat_A_old[i][j]-mat_A[i][j]);
            }
            sum_check+=Math.abs(vec_pi_old[i]-vec_pi[i]);
            sum_check+=Math.abs(vec_lambda_old[i]-vec_lambda[i]);
        }
        if(sum_check/(1.0*num_state*(num_state+2))<1.0e-7) flag++;

        mat_A_old=get_mat_copy(mat_A);
        vec_pi_old=vec_pi.concat();
        vec_lambda_old=vec_lambda.concat();

        E_res=HMM_E_step(vec_Xi, mat_A, vec_lambda, vec_pi);
		mat_Gamma=E_res[0];
		mat_Xi=E_res[1];
		
		loop++;
    }
    
    var vec_hidden=HMM_Viterbi(vec_Xi, mat_A, vec_lambda, vec_pi);
    //////////////////////////////////////////////////////////////////////

    var rate_hmm=get_mat(vec_Xi.length,2);
    
	var c_time=0.0;
    for(var n=0; n<vec_Xi.length; n++)
    {
        var state_id=vec_hidden[n];
        rate_hmm[n][0]=Math.round(c_time*100)/100.0;
        rate_hmm[n][1]=Math.round(vec_lambda[state_id]*100)/(bin_width*100.0);
        c_time+=bin_width;
    }
	
    return rate_hmm;
}
*/
