/*
 *
 * Implementation of inner moving avg algo. in R
 *
 */
 
#include <stdlib.h>
#include <stdio.h>
 
#include <R.h>
#include <R_ext/Memory.h>     /* R_alloc and S_alloc */
#include <Rdefines.h>
#include <Rinternals.h>

#include "collect.h" 

/**
 * General purpose average function
 */
void mean(double* arr, int start, int len, double* mn);


