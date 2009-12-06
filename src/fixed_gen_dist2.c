/*
 * Functions starting with r__ use the SEXP datatype
 */

#include "fixed_gen_dist2.h"
 
void find_low_id(
	int* 		low_id, 
	double* 	mid_arr, 
	int* 		start_arr, 
	int 		start_idx, 
	int 		stop_idx,
	double 		frag_len)
{
	int id = 0, i = 0;
	
	int start_pt = start_arr[start_idx];
	
	for(i=start_idx;i < stop_idx; i++)
	{
		if (mid_arr[i] - start_pt <= frag_len / 2.0)
		{
			id++;
		}
	}
	
	*low_id = id;
}

// wleft is returned as relative to sub_start_arr (start_idx) and a 0-based index
void get_wleft(
	int* wleft,
	int j,
	double sub_mid,
	int* start_arr,
	int start_idx,
	double frag_len)
{
	int i;
	
	for(i=0;i<j;i++)
	{
//		printf("test-(%d)--> %f <= %f \n", i, sub_mid - start_arr[i + start_idx - 1], frag_len / 2.0);
		if (sub_mid - start_arr[i + start_idx - 1] <= frag_len / 2.0)
		{
//			printf("  -> taken! (%d)\n", i);
			*wleft = i;
			return;
		}
	}
	
	*wleft = 0;
}

// wright is returned as relative to sub_start_arr (start_idx) and a 0-based index 
void get_wright(
	int* wright,
	int j,
	double sub_mid,
	int* stop_arr,
	int start_idx,
	int stop_idx,
	double frag_len)
{
	int i, len = stop_idx - start_idx + 1;
	
	*wright = 0;
	
	for(i = j-1;i<len;i++)
	{
		if (stop_arr[i + start_idx - 1] - sub_mid <= frag_len / 2.0)
		{
			if (j - 1 + i > *wright)
			{
				*wright = i;
			}
		}
	}
}
	

void get_wleft_wright_mean(
	double* mn,
	double* logr_arr,
	int start_idx,
	int wleft,
	int wright)
{
	mean(logr_arr, start_idx + wleft - 1, wright - wleft, mn);
}

void c__entry_pt(
	int* 			low_id_ptr,		// return low_id
	dbl_arr_ptr*	ma_low_ret,		// return ma_low
	int_arr_ptr*	w_low_ret,		// return w_low
	
	dbl_arr_ptr		ma_arr,			// current ma array
	int_arr_ptr		w_arr,			// current w array
	
	int*			uniq,
	int				uniq_idx,
	double* 		mid_arr,
	double* 		logr_arr,
	int* 			start_arr,
	int*			stop_arr,
	double 			frag_len)
{
	int i,j;
	
	// Uniq is R construct
	int start_idx = uniq[uniq_idx] - 1;
	int stop_idx = uniq[uniq_idx+1] - 2;

	find_low_id(low_id_ptr, mid_arr, start_arr, start_idx, stop_idx, frag_len);
	
	int low_id = *low_id_ptr;
	
	dbl_arr_ptr ma_low = NULL;
	int_arr_ptr w_low = NULL;
	
	if (low_id > 0)
	{
		alloc_dbl_arr(&ma_low, low_id);
		alloc_int_arr(&w_low, low_id);
	
		///// RAW Array access -- speedup /////
		for(i=0;i<low_id;i++)
		{
			mean(logr_arr,start_idx,i,&(ARRAY_RAW(ma_low)[i]));
	
			ARRAY_RAW(w_low)[i] = i + 1;
		}
		
		ma_low->len = low_id;	// have to assign our own len because we aren't using the macros
		w_low->len = low_id;	// have to assign our own len because we aren't using the macros
		///// End RAW Array access -- speedup /////
	
		
		for(i=0;i<ma_low->len;i++)
		{
			append_dbl(ma_arr, ARRAY_GET(ma_low,i));
		}
		for(i=0;i<w_low->len;i++)
		{
			append_int(w_arr, ARRAY_GET(w_low,i));
		}		
	}
	
	if (low_id < stop_idx - start_idx + 1)
	{
		for(j = low_id + 1; j <= stop_idx - start_idx + 1; j++)
		{
			double mn;
			
			int wleft;
			get_wleft(&wleft, j, mid_arr[start_idx + j - 1], start_arr, start_idx + 1, frag_len);
			
			int wright;
			get_wright(&wright, j, mid_arr[start_idx + j - 1], stop_arr, start_idx + 1, stop_idx + 1, frag_len);			
			
			get_wleft_wright_mean(&mn, logr_arr, start_idx + 1, wleft, wright);
			
			append_dbl(ma_arr,mn);
			
			append_int(w_arr,wright - wleft + 1);
		}
	}	
	
	free_int_arr(w_low);
	
	*ma_low_ret = ma_low;
}

SEXP r__entry_pt_fgd2(
	SEXP s_ma_arr,
	SEXP s_w_arr,
	
	SEXP s_uniq_arr,
	SEXP s_uniq_arr_len,

	SEXP s_mid_arr,
	SEXP s_start_arr,
	SEXP s_stop_arr,
	SEXP s_logr_arr,
	
	SEXP s_frag_len)
{
	int i;
	
	int uniq_arr_len = *INTEGER(s_uniq_arr_len);
	
	dbl_arr_ptr ma_ret;
	int_arr_ptr w_ret;	
	
	alloc_int_arr(&w_ret, 10000);
	alloc_dbl_arr(&ma_ret, 10000);
	
	for(i=0;i<uniq_arr_len-1;i++)
	{
		int low_id;
	
		dbl_arr_ptr ma_low;
		int_arr_ptr w_low;
	
		c__entry_pt(
			&low_id,
			&ma_low,
			&w_low,
			
			ma_ret,
			w_ret,
			
			INTEGER(s_uniq_arr),
			i,
			REAL(s_mid_arr),
			REAL(s_logr_arr),
			INTEGER(s_start_arr), 
			INTEGER(s_stop_arr),
			(double)*INTEGER(s_frag_len));
	}

	SEXP z, s_ma, s_w;
	
	PROTECT(z = NEW_LIST(2));
	PROTECT(s_ma = NEW_NUMERIC(ma_ret->len));
	PROTECT(s_w = NEW_INTEGER(w_ret->len));
	
	for(i=0;i<ma_ret->len;i++)
	{
		REAL(s_ma)[i] = ARRAY_GET(ma_ret,i);
	}
	for(i=0;i<w_ret->len;i++)
	{
		INTEGER(s_w)[i] = ARRAY_GET(w_ret,i);
	}
	
	SET_VECTOR_ELT(z,0,s_ma);
	SET_VECTOR_ELT(z,1,s_w);
	
	UNPROTECT(3);
	
	return z;
}

void mean(double* arr, int start, int len, double* mn)
{
	if (len < 0)
	{
		*mn = 0;
		return;
	}

	double accm = 0;
	int i;
	
	for(i = start;i <= start + len;i++)
	{
		accm += arr[i];
	}
	
	*mn = accm / ((double)len + 1.0);
}
