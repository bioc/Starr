#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)


#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


int DEBUG = 0;
int num_nodes = 0;


/**
  A	T	C	G
  0	1	2	3
**/
struct node{
int hit_id;
struct node *flink;
struct node *links[4];
};
typedef struct node *nodep;
typedef struct node node;



void destroy(nodep p){
  if(p!=NULL){
    destroy(p->links[0]);
    destroy(p->links[1]);
    destroy(p->links[2]);
    destroy(p->links[3]);
    free(p);
  }
  num_nodes=0;
}


int getPos(char c) {
  int pos=-1;
    switch(c)
{
   case 'A':
     pos = 0;
  break;

  case 'T':
     pos = 1;
  break;

  case 'C':
     pos = 2;
  break;

case 'G':
     pos = 3;
  break;
 
case 'a':
     pos = 0;
  break;

  case 't':
     pos = 1;
  break;

  case 'c':
     pos = 2;
  break;

case 'g':
     pos = 3;
  break;

}

return(pos);
}


int get_char(int p) {
char c='-';
    switch(p)
{
   case 0:
     c ='A';
  break;

  case 1:
     c = 'T';
  break;

  case 2:
     c= 'C';
  break;

case 3:
     c = 'G';
  break;
}

return(c);
}


int get_complement(char c) {
  char comp='-';
    switch(c)
{
   case 'A':
     comp = 'T';
  break;

  case 'T':
     comp = 'A';
  break;

  case 'C':
     comp = 'G';
  break;

case 'G':
     comp = 'C';
  break;
}

return(comp);
}


nodep init_tree() {
  nodep root = (nodep)malloc(sizeof(node));
  if(root == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }

  root->flink=NULL;
  root->hit_id=-1;

  int i;
  for(i=0; i<4; i++) {
    root->links[i] = NULL;
  }
  num_nodes++;
  return(root);
}


nodep insert_node(nodep current, char c) {
    nodep newnode;
    newnode=(nodep)malloc(sizeof(node));
     if(newnode == NULL) {
	return NULL;	
    }
    int i;
    for(i=0; i<4; i++) {
      newnode->links[i] = NULL;
    }
    newnode->hit_id = -1;
    newnode->flink = NULL;
    int pos = getPos(c);
  
    if(current->links[pos] == NULL) {
      current->links[pos] = newnode;
      current=newnode;
      num_nodes++;
    }
    else {
      current=current->links[pos];
    }
 
    return current;
}

int insert_word(nodep root, const char w[], int id) {
  int i;
  nodep currnode = root;
 
  i=0;
  while(w[i] != '\0') {
    if(currnode->links[getPos(w[i])] == NULL) {
    	currnode = insert_node(currnode, w[i]);
        if(currnode == NULL) {
		printf("Not enough memory!\nSee prameter nseq for help!\n");
		destroy(root);
		exit(-1);
    	}
    }
    else {
      currnode = currnode->links[getPos(w[i])];
    }
    if(DEBUG) {
      printf("Inserting: %c\n", w[i]);
    }
    i++;
  }
  
  if(currnode->hit_id == -1) {
    currnode->hit_id = id;
  }
  else {
      currnode->hit_id = -2;
  }

  if(DEBUG) {
    printf("pdict_id->%d\n", id);
  }

  return i;
}




void breadth_first_search(nodep root) {
  
  int j = 0;
  int counter = 0;
  int i;
  nodep *queue = (nodep *)malloc(sizeof(nodep)*num_nodes);
   if(queue == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }
  for(i=0; i<num_nodes; i++) {
	queue[i] = NULL;
  } 
  queue[0] = root;
  nodep currnode = queue[j++];
  counter++;

  // initialize flinks in first level
  for(i=0; i<=3; i++) {

      if(root->links[i] != NULL) {
	queue[j++] = currnode->links[i]; 	
	nodep fnode = currnode->links[i];
	fnode->flink = root;
	
	if(DEBUG) {
	  printf("flink: %c -> %c\n",get_char(i), '-');
	}
      }
    }

    if(queue[counter] != NULL) {
	currnode = queue[counter];
    } 

  while(counter++ < num_nodes) {
    
    for(i=0; i<=3; i++) {
      if(currnode->links[i] != NULL) {
	queue[j++] = currnode->links[i];
	nodep a = currnode;
	a = a->flink;
	
	nodep fnode = currnode->links[i];
	int pos = i;
	int stop = 0;
	
	while(stop != 1) {
	  // only follow first flink
	  int first_flink = 1;
	  if(a->links[pos] != NULL) {
	    a = a->links[pos];
	    
	    fnode->flink = a;
	    if(DEBUG) {
	      printf("flink: %c -> %c\n",get_char(i), get_char(pos));
	    }
	    stop=1;
	  }
	  else if(a->flink == NULL) {
	    fnode->flink = a;
	    if(DEBUG) {
	      printf("flink: %c -> %c\n",get_char(i), '-');
	    }
	    stop = 1;
	  }
	  else {
	    if(DEBUG) {
	     /* if(a->flink == NULL) {
		printf("follow flink from - " );
	      }
	      else {
		printf("follow flink from %c ", get_char(i));
	      }*/
	    }
	    a = a->flink;
	    first_flink = 0;
	    if(DEBUG) {
	     // printf("to %c\n", a->info);
	    }
	  }
	  if((a->hit_id >= 0) & (a->flink != NULL) & (first_flink == 1)) {
	int ahits = a->hit_id;
	if(fnode->hit_id == -1) {
	  fnode->hit_id = ahits;
	  /*if(DEBUG) {
	    hitsp fhits = fnode->dict_ids;
	    printf("Creating hit entries:");
	    while(fhits->next != NULL) {
	    printf(" %d", fhits->pdict_id);
	    fhits=fhits->next;
	  }
	  printf(" %d\n", fhits->pdict_id);
	  }*/
	}
         }
	}
	
      }
      
    }

   currnode = queue[counter];
  }
  free(queue);
}



SEXP find_ac2(SEXP dict, SEXP wcount_p, SEXP text, SEXP num_text, SEXP complementary, SEXP reverse, SEXP reverse_complementary, SEXP nseq) {
   
  SEXP text_names;
  PROTECT(text_names = NEW_CHARACTER(INTEGER(num_text)[0]));

  nodep root;
  int comp = INTEGER(complementary)[0];
  int revcomp = INTEGER(reverse_complementary)[0];
  int rev = INTEGER(reverse)[0];
  int index, t;

  int *wlens = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
 if(wlens == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }
  int *nindex = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
   if(nindex == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }
  int *ntext = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
   if(ntext == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }
  int *strand = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
   if(strand == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }
  int *nmatch = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
   if(nmatch == NULL) {
	printf("Not enough memory!\nSee prameter nseq for help!\n");
	destroy(root);
        exit(-1);
  }  


  int i;
  for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    wlens[i] = 0;
    ntext[i] = -1;
    nindex[i] = -1;
    strand[i] = 0;
    nmatch[i] = 0;
  }


  int step;
  for(step=0; step<INTEGER(wcount_p)[0]; step=min(step+INTEGER(nseq)[0], INTEGER(wcount_p)[0])) {
    int currTree=0;
    root = init_tree();

    for(currTree = step; currTree < min(step+INTEGER(nseq)[0], INTEGER(wcount_p)[0]); currTree++)           {
    wlens[currTree] = insert_word(root, CHAR (STRING_ELT (dict, currTree)), currTree);
    
    }
    Rprintf("Number of nodes: %d\n",num_nodes);
    // Initialize failure links
    breadth_first_search(root);
    Rprintf("Searching: ");
    

    for(t=0; t<INTEGER(num_text)[0]; t++) {
	
	i=0;
	gzFile fp;
	kseq_t *seq;

	
	fp = gzopen(CHAR (STRING_ELT (text, t)), "r");
	seq = kseq_init(fp);
	 kseq_read(seq);
		Rprintf("%s ", seq->name.s);
	SET_STRING_ELT(text_names, t, mkChar(seq->name.s));
	
	nodep currnode = root;
	nodep lastNorm = currnode;
	nodep lastComp = currnode;
	index = 0;


      while(seq->seq.s[index] != '\0') {
	int z;
	for(z=0; z<1+comp; z++) {
     	int pos = getPos(seq->seq.s[index]);
	currnode = lastNorm;
	if(z == 1) {
	  pos = getPos(get_complement(seq->seq.s[index]));
	  currnode = lastComp;
	}
	if(DEBUG) {
	  // printf("char=%c\n", get_complement(seq->seq.s[index]));
	}
	int stop = 0;
	while((currnode->links[pos] == NULL) & (stop != 1)) {
	    if(currnode->flink == NULL) {
	      stop = 1;
	    }
	    else {
	      currnode = currnode->flink;
	    }
	}
	if(currnode->links[pos] != NULL) {
	  currnode = currnode->links[pos];
	  if((currnode->hit_id >= 0) | (currnode->hit_id == -2)) {
	    if(DEBUG) {
	      // printf("id=");
	    }
	    nodep currhit = currnode;
	    while((currhit->hit_id >= 0) | (currhit->hit_id == -2)) {
	      if(currhit->hit_id >= 0) {
		int curr_id=currhit->hit_id; 
		ntext[curr_id] = t+1;
		nindex[curr_id] = index-wlens[curr_id]+2;
		strand[curr_id] = 1;
		nmatch[curr_id] = nmatch[curr_id]+1;
	      }
	      currhit=currhit->flink;
	    }

	  }
	}
	if(z == 1) {
	  lastComp = currnode;
	}
	else if(z == 0){
	  lastNorm = currnode;
	}
      }
	index++;
	
      }

      if((rev == 1) | (revcomp == 1)) {
      currnode = root;
      lastNorm = currnode;
      lastComp = currnode;
      int from = 0;
      if(rev == 0) { from=1;}
     index--;
      while(index != -1) {
	int z;
	
	
	for(z=from; z<1+revcomp; z++) {
     	int pos = getPos(seq->seq.s[index]);
	currnode = lastNorm;
	if(z == 1) {
	  pos = getPos(get_complement(seq->seq.s[index]));
	  currnode = lastComp;
	}
	if(DEBUG) {
	  // printf("char=%c\n", seq->seq.s[index]);
	}
	int stop = 0;
	while((currnode->links[pos] == NULL) & (stop != 1)) {
	    if(currnode->flink == NULL) {
	      stop = 1;
	    }
	    else {
	      currnode = currnode->flink;
	    }
	}
	if(currnode->links[pos] != NULL) {
	  currnode = currnode->links[pos];
	  if((currnode->hit_id >= 0) | (currnode->hit_id == -2)) {
	    if(DEBUG) {
	      // printf("id=");
	    }
	    nodep currhit = currnode;
	    while((currhit->hit_id >= 0) | (currhit->hit_id == -2)) {
	      if(currhit->hit_id >= 0) {
		int curr_id=currhit->hit_id; 
		ntext[curr_id] = t+1;
		nindex[curr_id] = index+1;
		strand[curr_id] = 0;
		nmatch[curr_id] = nmatch[curr_id]+1;
	      }
	      currhit=currhit->flink;
	    }

	  }
	}
	if(z == 1) {
	  lastComp = currnode;
	}
	else if(z == 0){
	  lastNorm = currnode;
	}
      }
	index--;
      }
      }
      
      kseq_destroy(seq);
      gzclose(fp);

  }
  destroy(root);
}

  Rprintf("\n");

    SEXP result, word_ind, word_t, wnames, word_strd;
    PROTECT(result = NEW_LIST(4));
	PROTECT(wnames = NEW_CHARACTER(4));
	SET_STRING_ELT(wnames, 0, mkChar("text"));
	SET_STRING_ELT(wnames, 1, mkChar("index"));
	SET_STRING_ELT(wnames, 2, mkChar("strand"));
	SET_STRING_ELT(wnames, 3, mkChar("textnames"));
	SET_NAMES(result, wnames);
	UNPROTECT(1);
    PROTECT(word_ind = NEW_INTEGER(INTEGER(wcount_p)[0]));
    PROTECT(word_t = NEW_INTEGER(INTEGER(wcount_p)[0]));
    PROTECT(word_strd = NEW_INTEGER(INTEGER(wcount_p)[0]));

    for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    
      if(nmatch[i] == 1) {
	INTEGER_POINTER(word_t)[i] = ntext[i];
	INTEGER_POINTER(word_ind)[i] = nindex[i];
	INTEGER_POINTER(word_strd)[i] = strand[i];
	 
	
	}
      else if(nmatch[i] > 1){
	 INTEGER_POINTER(word_t)[i] = -2;
	INTEGER_POINTER(word_ind)[i] = -2;
	INTEGER_POINTER(word_strd)[i] = -2;
    
      }
      else if(nmatch[i] == 0){
	 INTEGER_POINTER(word_t)[i] = -1;
	INTEGER_POINTER(word_ind)[i] = -1;
	INTEGER_POINTER(word_strd)[i] = -1;
    
      }
    }

    SET_ELEMENT(result, 0, word_t);
    SET_ELEMENT(result, 1, word_ind);
    SET_ELEMENT(result, 2, word_strd);
    SET_ELEMENT(result, 3, text_names);
    UNPROTECT(5);
    
    free(nmatch);
    free(nindex);
    free(ntext);
    free(strand);
    free(wlens);
      
    return result;
}
