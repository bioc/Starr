#include "rm_small_peak.h"


//	peak.start=which(bdd.method==1)
//	peak.start=c(min(peak.start),peak.start[(which(peak.start[2:length(peak.start)]-peak.start[1:(length(peak.start)-1)]>1))+1])
// twooff = peak.start[2:length(peak.start)] = peak.start[idx+1]
// lastoff = peak.start[1:(length(peak.start)-1)] = peak.start[idx]
// which(twooff-lastoff>1)+1
void get_peak_start1(
	int_arr_ptr* ret,
	int* bdd_method,
	int bdd_method_len,
	int bdd_method_tgt_val)
{
	int i,mn = -1;

	int_arr_ptr pkst_which;
	
	alloc_int_arr(&pkst_which, bdd_method_len); // larget is when all match
	
	for(i=0;i<bdd_method_len;i++)
	{
		if (bdd_method[i] == bdd_method_tgt_val)
		{
			if (mn == -1)
			{
				mn = i;
			}
			
			append_int(pkst_which, i);
		}
	}
	
	int_arr_ptr pkstr;
	
	alloc_int_arr(&pkstr, pkst_which->len+1);
	
	append_int(pkstr, mn);
	
	int subt = 0;
	
	if (bdd_method_tgt_val == 0)
	{
		subt = -1;
	}
	
	for(i=0;i<pkst_which->len - 1;i++)
	{
		if (ARRAY_GET(pkst_which,i+1) - ARRAY_GET(pkst_which,i) > 1)
		{
			append_int(pkstr, ARRAY_GET(pkst_which,i+1) + subt);
		}
	}
	
	*ret = pkstr;
	
	free_int_arr(pkst_which);
}

SEXP r__rm_small_peak(
	SEXP bdd_method_s,
	SEXP bdd_method_len,
	SEXP minrun_s)
{
	int i,j;
	double minrun;

	int_arr_ptr peak_start, peak_end, peak_size;	
	
	minrun = *REAL(minrun_s);
	int* bdd_method = INTEGER(bdd_method_s);

	get_peak_start1(&peak_start, bdd_method, *INTEGER(bdd_method_len),1);
	get_peak_start1(&peak_end, bdd_method, *INTEGER(bdd_method_len),0);	
	
	// Remove the first entry -- lots of work :(
	if(ARRAY_GET(peak_end,0) == 0)
	{
		int_arr_ptr tmp;
		alloc_int_arr(&tmp, peak_end->len-1);
		
		for(i=1;i<peak_end->len;i++)
		{
			append_int(tmp, ARRAY_GET(peak_end,i));
		}
		
		free_int_arr(peak_end);
		peak_end = tmp;
	}	
	
	if(peak_start->len > peak_end->len)
	{
		peak_start->len = peak_start->len - 1;
	}	
	
	alloc_int_arr(&peak_size, peak_start->len);
	
	for(i=0;i<peak_start->len;i++)
	{
		append_int(peak_size, ARRAY_GET(peak_end,i) - ARRAY_GET(peak_start,i) + 1);
	}
	
	int_arr_ptr id1, id2;
	
	alloc_int_arr(&id1, peak_size->len);
	alloc_int_arr(&id2, peak_size->len);
	
	for(i=0;i<peak_size->len;i++)
	{
		if (ARRAY_GET(peak_size,i) < minrun)
		{
			append_int(id1, ARRAY_GET(peak_start,i));
			append_int(id2, ARRAY_GET(peak_end,i));			
		}
	}
	
	if (id1->len > 0)
	{
		for(i=0;i<id1->len;i++)
		{
			for(j = ARRAY_GET(id1,i); j <= ARRAY_GET(id2,i); j++)
			{
				bdd_method[j] = 0;
			}
		}		
	}
	
	// Reload all data with the new bdd_method
	free_int_arr(peak_start);
	free_int_arr(peak_end);
	free_int_arr(peak_size);
	
	get_peak_start1(&peak_start, bdd_method, *INTEGER(bdd_method_len),1);
	get_peak_start1(&peak_end, bdd_method, *INTEGER(bdd_method_len),0);	
	
	// Remove the first entry -- lots of work :(
	if(ARRAY_GET(peak_end,0) == 0)
	{
		int_arr_ptr tmp;
		alloc_int_arr(&tmp, peak_end->len-1);
		
		for(i=1;i<peak_end->len;i++)
		{
			append_int(tmp, ARRAY_GET(peak_end,i));
		}
		
		free_int_arr(peak_end);
		peak_end = tmp;
	}	
	
	if(peak_start->len > peak_end->len)
	{
		peak_start->len = peak_start->len - 1;
	}	
	
	alloc_int_arr(&peak_size, peak_start->len);
	
	for(i=0;i<peak_start->len;i++)
	{
		append_int(peak_size, ARRAY_GET(peak_end,i) - ARRAY_GET(peak_start,i) + 1);
	}	
	
	SEXP z, s_pkst, s_pken, s_size;
	
	PROTECT(z = NEW_LIST(3));
	PROTECT(s_pkst = NEW_INTEGER(peak_start->len));
	PROTECT(s_pken = NEW_INTEGER(peak_end->len));	
	PROTECT(s_size = NEW_INTEGER(peak_size->len));
	
	for(i=0;i<peak_end->len;i++)
	{
		// +1 for one based arrays
		INTEGER(s_pken)[i] = ARRAY_GET(peak_end,i) + 1;
		INTEGER(s_size)[i] = ARRAY_GET(peak_size,i);
		INTEGER(s_pkst)[i] = ARRAY_GET(peak_start,i) + 1;		
	}
	
	SET_VECTOR_ELT(z,0,s_pkst);
	SET_VECTOR_ELT(z,1,s_pken);
	SET_VECTOR_ELT(z,2,s_size);
	
	UNPROTECT(4);
	
	free_int_arr(peak_start);
	free_int_arr(peak_end);
	free_int_arr(peak_size);
	free_int_arr(id1);
	free_int_arr(id2);
	
	return z;
}































