#include "collect.h"

extern void ensure_cap_int(int_arr_ptr p, int sz)
{
	if (p->storage_len < sz)
	{
		grow_int_arr(p, sz*2);
	}
}

extern void ensure_cap_dbl(dbl_arr_ptr p, int sz)
{
	if (p->storage_len < sz)
	{
		grow_dbl_arr(p, sz*2);
	}
}


void alloc_int_arr(int_arr_ptr* arr, int size)
{
	*arr = malloc(sizeof(int_arr));
	ARRAY_RAW((*arr)) = malloc(sizeof(int)*size);
	(*arr)->storage_len = size;
	(*arr)->len = 0;
}

void free_int_arr(int_arr_ptr arr)
{
	if (arr == NULL)
	{
		return;
	}
	
	free(ARRAY_RAW(arr));
	free(arr);
}

void alloc_dbl_arr(dbl_arr_ptr* arr, int size)
{
	*arr = malloc(sizeof(dbl_arr));
	ARRAY_RAW((*arr)) = malloc(sizeof(double)*size);
	(*arr)->storage_len = size;
	(*arr)->len = 0;	
}

void free_dbl_arr(dbl_arr_ptr arr)
{
	if (arr == NULL)
	{
		return;
	}
	
	free(ARRAY_RAW(arr));
	free(arr);
}

void grow_dbl_arr(dbl_arr_ptr arr, int size)
{
	int i;
	double* tmp = malloc(sizeof(double)*(arr->storage_len+size));
	
	for(i=0;i<arr->storage_len;i++)
	{
		tmp[i] = ARRAY_RAW(arr)[i];
	}
	
	free(ARRAY_RAW(arr));
	ARRAY_RAW(arr) = tmp;
	arr->storage_len += size;
}

void grow_int_arr(int_arr_ptr arr, int size)
{
	int i;
	int* tmp = malloc(sizeof(int)*(arr->storage_len+size));
	
	for(i=0;i<arr->storage_len;i++)
	{
		tmp[i] = ARRAY_RAW(arr)[i];
	}
	
	free(ARRAY_RAW(arr));
	ARRAY_RAW(arr) = tmp;
	arr->storage_len += size;
}

void append_dbl(dbl_arr_ptr arr, double ap)
{
	if (arr->storage_len < arr->len + 1)
	{
		grow_dbl_arr(arr, arr->storage_len*2);
	}
	
	ARRAY_RAW(arr)[arr->len] = ap;
	arr->len++;
}

void append_int(int_arr_ptr arr, int ap)
{
	if (arr->storage_len < arr->len + 1)
	{
		grow_int_arr(arr, arr->storage_len*2);
	}
	
	ARRAY_RAW(arr)[arr->len] = ap;
	arr->len++;
}

