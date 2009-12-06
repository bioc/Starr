/*
 *  sort_by_genomic.h
 *  RPort
 *
 *  Created by hinz on March 24 2008.
 *  Copyright 2008 __MyCompanyName__. All rights reserved.
 *
 */

#include <stdlib.h>
#include <stdio.h>
 
#include <R.h>
#include <R_ext/Memory.h>     /* R_alloc and S_alloc */
#include <Rdefines.h>
#include <Rinternals.h>

#include "collect.h"

typedef struct { int start; int stop; double logR; } t_rec;
typedef t_rec* t_rec_ptr;

#define REG_ID_ALLOC_TUNE 10
