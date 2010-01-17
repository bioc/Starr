#ifndef __COLLECT__
#define __COLLECT__

#include <stdlib.h>

/**
 * Struct to hold info about an array
 *
 * @item storage_len current length of the backing array
 * @item len the length of the array (values that have data)
 * @item arr the actual array
 */
typedef struct
{
	int storage_len;
	int len;
	double* arrdata;
} dbl_arr;

/**
 * New type for easy reading and fewer *s
 */
typedef dbl_arr* dbl_arr_ptr;

/**
 * Struct to hold info about an array
 *
 * @item storage_len current length of the backing array
 * @item len the length of the array (values that have data)
 * @item arr the actual array
 */
typedef struct
{
	int storage_len;
	int len;
	int* arrdata;
} int_arr;

/**
 * New type for easy reading and fewer *s
 */
typedef int_arr* int_arr_ptr;

/**
 * Access to the arrays
 * Provides bounds auto-extension and length checking
 */
/** Get the raw array data (not prefered) */
#define ARRAY_RAW(X)           (X)->arrdata

/** Grab an int or double from the data array */
#define ARRAY_GET(X,IDX)       (X)->arrdata[IDX]
#define ARRAY_GET_INT(X,IDX)   (X)->arrdata[IDX]
#define ARRAY_GET_DBL(X,IDX)   (X)->arrdata[IDX]

/** 
 * Add data to the array. Keeps length up to date and
 * Maintains the proper capacity
 */
#define ARRAY_SET_INT(X,IDX,V) 								\
							if(IDX >= (X)->storage_len) 	\
							{ 								\
								ensure_cap_int((X),(IDX));	\
							}								\
							if (IDX >= (X)->len)			\
							{								\
								(X)->len = IDX + 1;			\
							}								\
							(X)->arrdata[IDX] = (V)	

#define ARRAY_SET_DBL(X,IDX,V) 								\
							if(IDX >= (X)->storage_len) 	\
							{ 								\
								ensure_cap_dbl((X),(IDX));	\
							}								\
							if (IDX >= (X)->len)			\
							{								\
								(X)->len = IDX + 1;			\
							}								\
							(X)->arrdata[IDX] = (V)	

/**
 * Make sure that the given array can store an item at a given index
 *
 * @param p the array pointer
 * @param sz the index
 */
extern void ensure_cap_int(int_arr_ptr p, int sz);

/**
 * Make sure that the given array can store an item at a given index
 *
 * @param p the array pointer
 * @param sz the index
 */
extern void ensure_cap_dbl(dbl_arr_ptr p, int sz);

/**
 * Allocate space for a new array struct
 * 
 * @param arr a pointer to an array pointer
 * @param size the size of the backing array to allocate
 */
extern void alloc_int_arr(int_arr_ptr* arr, int size);

/**
 * Free the resources used by an array struct
 *
 * @param arr the array struct to free
 */
extern void free_int_arr(int_arr_ptr arr);

/**
 * Allocate space for a new array struct
 *
 * @param arr a pointer to an array pointer
 * @param size the size of the backing array to allocate
 */
extern void alloc_dbl_arr(dbl_arr_ptr* arr, int size);

/**
 * Free the resources used by an array struct
 *
 * @param arr the array struct to free
 */
extern void free_dbl_arr(dbl_arr_ptr arr);

/**
 * Increase the size of an array
 *
 * @param arr the array to increase
 * @param size the new size of the array
 */
extern void grow_dbl_arr(dbl_arr_ptr arr, int size);

/**
 * Increase the size of an array
 *
 * @param arr the array to increase
 * @param size the new size of the array
 */
extern void grow_int_arr(int_arr_ptr arr, int size);

/**
 * Append a double to an array. The array will grow automatically if needed
 *
 * @param arr the array pointer
 * @param ap the double to add
 */
extern void append_dbl(dbl_arr_ptr arr, double ap);

/**
 * Append an integer to an array. The array will grow automatically if needed
 *
 * @param arr the array pointer
 * @param ap the integer to append
 */
extern void append_int(int_arr_ptr arr, int ap);

#endif

