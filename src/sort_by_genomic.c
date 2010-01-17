/*
 *  sort_by_genomic.c
 *  RPort
 *
 *  Created by hinz on March 24 2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 * @version=1.0.1
 * @bugfix 7/20/2008 AH Fixed counter in the region id method (version updated to 1.0.1
 *
 */

#include "sort_by_genomic.h"


///// TODO
///// Implement a real sorting algo
void sort2(t_rec_ptr** recsx, int len)
{	
	t_rec_ptr* recs = *recsx;
	
	if (len == 0)
	{
		return;
	}
	
	int i,j;
	
	for(i=0;i<len;i++)
	{
		t_rec_ptr value = recs[i];
		j = i - 1;
		
		while(j >= 0 && recs[j]->start > value->start)
		{
			recs[j+1] = recs[j];
			j--;
		}
		recs[j+1] = value;
	}
}

void get_reg_tmp(
	int_arr_ptr* ret_arr, 
	int_arr_ptr* ret_sub_start_sort,
	int_arr_ptr* ret_sub_stop_sort,
	dbl_arr_ptr* ret_sub_logR_sort,
	int_arr_ptr start_arr, 
	int_arr_ptr stop_arr,
	dbl_arr_ptr logR_arr,
	int gap)
{
	int_arr_ptr sub_start_sort;
	int_arr_ptr sub_stop_sort;
	dbl_arr_ptr sub_logR_sort;
	
	int i;
	
	t_rec_ptr* recs = malloc(sizeof(t_rec_ptr)*start_arr->len);
	
	for(i=0;i<start_arr->len;i++)
	{
		t_rec_ptr rec = malloc(sizeof(t_rec));
		rec->start = ARRAY_GET(start_arr,i);
		rec->stop = ARRAY_GET(stop_arr,i);
		rec->logR = ARRAY_GET(logR_arr,i);		
		
		recs[i] = rec;
	}
	
	sort2(&recs,start_arr->len);
	
	alloc_int_arr(&sub_start_sort, start_arr->len);
	alloc_int_arr(&sub_stop_sort, start_arr->len);
	alloc_dbl_arr(&sub_logR_sort, start_arr->len);
	
	///// RAW array access block, length updated manually /////
	for(i=0;i<start_arr->len;i++)
	{
		t_rec_ptr rec = recs[i];
		ARRAY_RAW(sub_start_sort)[i] = rec->start;
		ARRAY_RAW(sub_stop_sort)[i] = rec->stop;
		ARRAY_RAW(sub_logR_sort)[i] = rec->logR;
		
		free(rec);
	}	
	///// END -> RAW array access block, length updated manually /////	
	
	sub_start_sort->len = start_arr->len;
	sub_stop_sort->len = start_arr->len;
	sub_logR_sort->len = start_arr->len;	
	
	free(recs);
	
	int_arr_ptr ret;
	
	alloc_int_arr(&ret, start_arr->len-1);
	
	for(i=0;i<start_arr->len-1;i++)
	{
		if (ARRAY_GET(sub_start_sort,i+1) - ARRAY_GET(sub_stop_sort,i) > gap)
		{
			append_int(ret, i+1);
		}
	}
	
	append_int(ret, stop_arr->len);	
	
	*ret_arr = ret;
	
	*ret_sub_start_sort = sub_start_sort;
	*ret_sub_stop_sort = sub_stop_sort;
	*ret_sub_logR_sort = sub_logR_sort;
}


void get_reg_id(
	int_arr_ptr* ret, 
	int* old_ret_id, 
	int old_ret_id_len, 
	int* k_ptr, 
	int_arr_ptr reg_tmp)
{
	int k = *k_ptr;
	int reg_tmp_len = reg_tmp->len;

	int_arr_ptr repArr;
	
	alloc_int_arr(&repArr, reg_tmp_len);
	
	int i,j;
	
	for(i=k;i<=k+reg_tmp_len-1;i++)
	{
		ARRAY_SET_INT(repArr,i-k,i);		// Protected set
	}
	
	int_arr_ptr numRepsArr;
	
	alloc_int_arr(&numRepsArr, reg_tmp_len+1);
	
	append_int(numRepsArr, ARRAY_GET(reg_tmp,0));
	
	for(i=1;i<reg_tmp_len;i++)
	{
		append_int(numRepsArr, ARRAY_GET(reg_tmp,i) - ARRAY_GET(reg_tmp,i-1));
	}
	
	int_arr_ptr ret_arr;
	
	alloc_int_arr(&ret_arr, numRepsArr->len*REG_ID_ALLOC_TUNE);
	
	for(i=0;i<numRepsArr->len;i++)
	{
		for(j=0; j < ARRAY_GET(numRepsArr,i); j++)
		{
			append_int(ret_arr, ARRAY_GET(repArr,i));
		}
	}
	
	alloc_int_arr(ret, old_ret_id_len + ret_arr->len);
	
	for(i=0;i<old_ret_id_len;i++)
	{
		append_int(*ret, old_ret_id[i]);
	}
	for(i=0;i<ret_arr->len;i++)
	{
		append_int(*ret, ARRAY_GET(ret_arr,i));
	}
	
	*k_ptr = k + reg_tmp->len;
	
	free_int_arr(repArr);
	free_int_arr(numRepsArr);
	free_int_arr(ret_arr);
}

void get_sub_arrays(
	int_arr_ptr* start_ret, 
	int_arr_ptr* stop_ret, 
	int_arr_ptr* chr_ret, 
	dbl_arr_ptr* logR_ret,
	int* start_arr,
	int* stop_arr,
	double* logR_arr,
	int chr_len,
	int* chr_arr,
	int chridx)
{
	int i;
	
	int_arr_ptr start;
	int_arr_ptr stop;
	int_arr_ptr chr; 
	dbl_arr_ptr logR;
	
	alloc_int_arr(&start, chr_len); // max size= chr_len
	alloc_int_arr(&stop, chr_len);
	alloc_int_arr(&chr, chr_len);
	alloc_dbl_arr(&logR, chr_len);
	
	for(i=0;i<chr_len;i++)
	{
		if (chr_arr[i] == chridx)
		{
			append_int(start, start_arr[i]);
			append_int(stop, stop_arr[i]);			
			append_int(chr, chridx);
			append_dbl(logR, logR_arr[i]);
		}
	}
	
	*start_ret = start;
	*stop_ret = stop;
	*chr_ret = chr;
	*logR_ret = logR;
		
}

SEXP r__entry_point(
	SEXP compl_start_arr,
	SEXP compl_stop_arr,
	SEXP compl_start_arr_len,
	SEXP compl_logR_arr,
	SEXP chr_arr,
	SEXP chr_arr_len,
	SEXP chrID,
	SEXP chrID_len,
	SEXP gap,
	SEXP s_reg_id,
	SEXP s_reg_id_len)
{

	int ii=0,i,main_len = *INTEGER(chrID_len);
	
	int_arr_ptr running_start;
	int_arr_ptr running_stop;
	int_arr_ptr running_chr;
	int_arr_ptr running_reg_id;
	
	alloc_int_arr(&running_start, 10000);
	alloc_int_arr(&running_stop, 10000);
	alloc_int_arr(&running_chr, 10000);
	alloc_int_arr(&running_reg_id, 10000);
	
	dbl_arr_ptr running_logR;
	
	alloc_dbl_arr(&running_logR, 10000);
	
	
	int k=0;
	
	for(ii=0;ii<main_len;ii++)
	{
		int_arr_ptr sub_start_sort;
		int_arr_ptr sub_stop_sort;
		dbl_arr_ptr sub_logR_sort;
		int_arr_ptr reg_tmp;
		
		int_arr_ptr start_sub_arr;
		int_arr_ptr stop_sub_arr;
		int_arr_ptr chr_sub_arr;
		dbl_arr_ptr logR_sub_arr;
		
		get_sub_arrays(
			&start_sub_arr,
			&stop_sub_arr,
			&chr_sub_arr,
			&logR_sub_arr,
			INTEGER(compl_start_arr),
			INTEGER(compl_stop_arr),
			REAL(compl_logR_arr),
			*INTEGER(chr_arr_len),
			INTEGER(chr_arr),
			INTEGER(chrID)[ii]);
		
		get_reg_tmp(
			&reg_tmp,
			&sub_start_sort,
			&sub_stop_sort,
			&sub_logR_sort,
			start_sub_arr,
			stop_sub_arr,
			logR_sub_arr,
			*INTEGER(gap));
		
		int_arr_ptr reg_id;
		
		get_reg_id(
			&reg_id,
			INTEGER(s_reg_id), 
			*INTEGER(s_reg_id_len), 
			&k,
			reg_tmp);
		
		for(i=0;i<sub_start_sort->len;i++)
		{
			append_int(running_start, ARRAY_GET(sub_start_sort,i));
			append_int(running_stop, ARRAY_GET(sub_stop_sort,i));
			/// **BUG FIX** 
			//append_int(running_chr, ARRAY_GET(chr_sub_arr,i));
			append_int(running_chr, INTEGER(chrID)[ii]);
			append_dbl(running_logR, ARRAY_GET(sub_logR_sort,i));
		}
		
		for(i=0;i<reg_id->len;i++)
		{
			append_int(running_reg_id, ARRAY_GET(reg_id,i));
		}
		
		
		free_int_arr(reg_tmp);
		free_int_arr(sub_start_sort);
		free_int_arr(sub_stop_sort);
		free_int_arr(reg_id);
	}
	
	SEXP z,s_start,s_stop,s_id,s_logR,s_chr;
	
	PROTECT(s_start = NEW_INTEGER(running_start->len));
	PROTECT(s_stop = NEW_INTEGER(running_start->len));
	PROTECT(s_id = NEW_INTEGER(running_reg_id->len));
	PROTECT(s_logR = NEW_NUMERIC(running_start->len));
	PROTECT(s_chr = NEW_INTEGER(running_start->len));
	
	for(i=0;i<running_reg_id->len;i++)
	{
		INTEGER(s_id)[i] = ARRAY_GET(running_reg_id,i);
	}
	for(i=0;i<running_start->len;i++)
	{
		INTEGER(s_start)[i] = ARRAY_GET(running_start,i);
		INTEGER(s_stop)[i] = ARRAY_GET(running_stop,i);
		INTEGER(s_chr)[i] = ARRAY_GET(running_chr,i);
		REAL(s_logR)[i] = ARRAY_GET(running_logR,i);
	}
	
	PROTECT(z = NEW_LIST(6));
	
	SET_VECTOR_ELT(z,1,s_start);
	SET_VECTOR_ELT(z,2,s_stop);	
	SET_VECTOR_ELT(z,3,s_id);	
	SET_VECTOR_ELT(z,4,s_logR);
	SET_VECTOR_ELT(z,5,s_chr);
	
	UNPROTECT(6);
	
	free_int_arr(running_start);
	free_int_arr(running_stop);
	free_int_arr(running_chr);
	free_dbl_arr(running_logR);
	free_int_arr(running_reg_id);			
	
	return z;
}
