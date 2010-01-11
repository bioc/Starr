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

/**
  Definition of data structure.
**/

int DEBUG = 0;
int num_nodes = 1;
int c_malloc = 0;
int m_malloc = 0;
int c_free = 0;
int h_malloc = 0;
int l_malloc = 0;
int t_malloc = 0;
int m_free = 0;
int h_free = 0;
int l_free = 0;

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



struct list{
   nodep entry;          /* Storing marks of a node */
    struct list *next;    /* Storing next address */
};
/*****  Redefining struct list as node  *****/
typedef struct list *listp;
typedef struct list list;

struct matches{
   int index;
   int ntext;
    int strand;
    struct matches *next;    /* Storing next address */
};
/*****  Redefining struct list as node  *****/
typedef struct matches *matchp;
typedef struct matches matches;


void destroy(nodep p){
  if(p!=NULL){
     destroy(p->links[0]);
     destroy(p->links[1]);
    destroy(p->links[2]);
      destroy(p->links[3]);
     free(p);
    c_free++;
 }
  num_nodes=1;
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
  t_malloc++;
  nodep root = (nodep)malloc(sizeof(node));
  root->flink=NULL;
  root->hit_id=-1;

  int i;
  for(i=0; i<4; i++) {
    root->links[i] = NULL;
    //printf("Number %d\n", i);
  }
 
  return(root);
}

nodep insert_node(nodep current, char c) {
    nodep newnode;
    newnode=(nodep)malloc(sizeof(node));
    t_malloc++;
     int i;
  for(i=0; i<4; i++) {
    newnode->links[i] = NULL;
    //printf("Number %d\n", i);
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

void insert_word(nodep root, const char w[], int id) {
  int i;
  nodep currnode = root;
 
  i=0;
  while(w[i] != '\0') {
    if(currnode->links[getPos(w[i])] == NULL) {
    	currnode = insert_node(currnode, w[i]);
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
}


void insert_word2(nodep root, char *w, int id) {
  int i;
  nodep currnode = root;
 
  i=0;
  while(w[i] != '\0') {
    if(currnode->links[getPos(w[i])] == NULL) {
    	currnode = insert_node(currnode, w[i]);
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
  printf("w[i]=%d\n", i);
}

listp init_list(nodep node) {
  listp head = (listp)malloc(sizeof(list));
  c_malloc++;
  head->entry = node;
  head->next = NULL;
 
  return(head);
}

// returns new tail of the list
listp add(listp tail, nodep node) {
  listp new = (listp)malloc(sizeof(list));
  l_malloc++;
  new->entry = node;
  new->next = NULL;
  tail->next = new;

  return(new);
}


listp breadth_first_search(nodep root) {
  listp queue = init_list(root);
  listp curr = queue;
  listp tail = queue;

  nodep currnode = curr->entry;
  int i;

  // initialize flinks in first level
  for(i=0; i<=3; i++) {
      if(currnode->links[i] != NULL) {
	tail = add(tail, currnode->links[i]);
	nodep fnode = currnode->links[i];
	fnode->flink = root;
	if(DEBUG) {
	  //printf("flink: %c -> %c\n",fnode->info, root->info);
	}
	//nodep a = currnode->links[i];
	//printf("adding: %c, %d, %d\n", a->info, i, j);
      }
    }
    
curr = curr->next;
    if(curr != NULL) {
      currnode = curr->entry;
    }
  
  while(curr != NULL) {
    
    for(i=0; i<=3; i++) {
      if(currnode->links[i] != NULL) {
	tail = add(tail, currnode->links[i]);
	nodep a = currnode;
	a = a->flink;
	//printf("adding: %c, %d, %d\n", a->info, i, j);
	
	
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
	      //printf("flink: %c -> %c\n",fnode->info, a->info);
	    }
	    stop=1;
	  }
	  else if(a->flink == NULL) {
	    fnode->flink = a;
	    if(DEBUG) {
	      //printf("flink: %c -> %c\n",fnode->info, a->info);
	    }
	    stop = 1;
	  }
	  else {
	    if(DEBUG) {
	      //printf("follow flink from %c ", a->info);
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
	  //printf("a->info: %c", a->info);
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
	else {
	  /*hitsp fhits = fnode->dict_ids;
	  //printf("updating hits:", fhits->pdict_id);
	  while(fhits->next != NULL) {
	    //printf(" %d", fhits->pdict_id);
	    fhits=fhits->next;
	  }
	  //printf(" %d, ", fhits->pdict_id);
	  fhits->next=ahits;
	  if(DEBUG) {
	    printf(" with:", ahits->pdict_id);
	    while(ahits->next != NULL) {
	    printf(" %d", ahits->pdict_id);
	    ahits=ahits->next;
	  }
	  printf(" %d\n", ahits->pdict_id);
	  }*/
	 }
	  }
	}
	
      }
      
    }

    
    
    curr = curr->next;
    if(curr != NULL) {
      currnode = curr->entry;
    }
  }
  return(queue);
}


void breadth_first_search1(nodep root) {
  //listp queue = init_list(root);
  //listp curr = queue;
  //listp tail = queue;

  int j = 0;
  int counter = 0;
  int i;
  nodep *queue = (nodep *)malloc(sizeof(nodep)*num_nodes);
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
	  //printf("a->info: %c", a->info);
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
	else {
	 /* hitsp fhits = fnode->dict_ids;
	  if(DEBUG) {
	    printf("updating hits ");
	  }
	  while(fhits->next != NULL) {
	    //printf(" %d", fhits->pdict_id);
	    fhits=fhits->next;
	  }
	  //printf(" %d, ", fhits->pdict_id);
	  fhits->next=ahits;
	  if(DEBUG) {
	    printf(" with:", ahits->pdict_id);
	    while(ahits->next != NULL) {
	    printf(" %d", ahits->pdict_id);
	    ahits=ahits->next;
	  }
	  printf(" %d\n", ahits->pdict_id);
	  }*/
	}
	  }
	}
	
      }
      
    }

    //printf("2\n");
  
   currnode = queue[counter];
  }
  free(queue);
}

void ac_from_fasta(const char *path, nodep root) {
  int currTree = 0;
  gzFile words;
	  kseq_t *dictionary;
	  int l;
	  words = gzopen(path, "r");
	  dictionary = kseq_init(words);
	  while ((l = kseq_read(dictionary)) >= 0) {
		  insert_word2(root, dictionary->seq.s, currTree++);
	  }
	  kseq_destroy(dictionary);
	  gzclose(words);
}

SEXP find_ac2(SEXP dict, SEXP wcount_p, SEXP text, SEXP num_text, SEXP complementary, SEXP reverse, SEXP reverse_complementary, SEXP nseq) {
  
 
  
  nodep root;
 // listp queue;
  
 // printf("number of nodes: %d\nSearching...", num_nodes);
  
   int comp = INTEGER(complementary)[0];
    int revcomp = INTEGER(reverse_complementary)[0];
   // int fastad = INTEGER(fasta_dict)[0];
int rev = INTEGER(reverse)[0];
    int index, t;

  
	  

  
	

  //matchp matcher[INTEGER(wcount_p)[0]];
  int *nindex = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
  int *ntext = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
  int *strand = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
  int *nmatch = (int*)malloc(INTEGER(wcount_p)[0] * sizeof(int));
  //int nmatch[INTEGER(wcount_p)[0]];
  int i;
  for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    ntext[i] = -1;
    nindex[i] = -1;
    strand[i] = 0;
    nmatch[i] = 0;
  }


int step;
for(step=0; step<INTEGER(wcount_p)[0]; step=min(step+INTEGER(nseq)[0], INTEGER(wcount_p)[0])) {
	int currTree=0;
	root = init_tree();

	//if(fastad == 0) {
	  for(currTree = step; currTree < min(step+INTEGER(nseq)[0], INTEGER(wcount_p)[0]); currTree++)           {
	  insert_word(root, CHAR (STRING_ELT (dict, currTree)), currTree);
	  
	  }
	/*}
	else {
	  ac_from_fasta(CHAR (STRING_ELT (dict, 0)), root);
	  //char path = (char) CHAR (STRING_ELT (dict, 0));
	  
	}*/
	//printf("currTree=%d\n", currTree-1);
	//queue = 
	Rprintf("Number of nodes: %d\n",num_nodes);
      breadth_first_search1(root);
	
	
	
Rprintf("Searching: ");
    for(t=0; t<INTEGER(num_text)[0]; t++) {
	
	i=0;
	gzFile fp;
	kseq_t *seq;

	
	fp = gzopen(CHAR (STRING_ELT (text, t)), "r");
	seq = kseq_init(fp);
	 kseq_read(seq);
		Rprintf("%s ", seq->name.s);
	//printf("return value: %d\n", l);

	
	/*

	for(step=0; step<INTEGER(wcount_p)[0]; step=min(step+INTEGER(nseq)[0], INTEGER(wcount_p)[0])) {
	int currTree;
	//printf("step=%d\n", step);
	root = init_tree();
	//printf("hier2\n");
	for(currTree = step; currTree < min(step+INTEGER(nseq)[0], INTEGER(wcount_p)[0]); currTree++)           {
	insert_word(root, CHAR (STRING_ELT (dict, currTree)), currTree);
	
	} 
	//printf("currTree=%d\n", currTree-1);
	queue = breadth_first_search(root);
	*/

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
		//if(matcher[curr_id] == NULL) {
		      /*matchp new = (matchp)malloc(sizeof(matches));
		      m_malloc++;
		      new->index=index;
		      new->ntext = t;
		      new->strand = 1;
		      new->next = NULL;
		      matcher[curr_id] = new;*/
		      ntext[curr_id] = t;
		      nindex[curr_id] = index+1;
		      strand[curr_id] = 1;
		      nmatch[curr_id] = nmatch[curr_id]+1;
		      //printf("id=%d, index=%d, text=%d\n", curr_id, index, t);
		//}
		/*else if (matcher[curr_id] != NULL){
		    matchp curr_match = matcher[curr_id];
		    while(curr_match->next != NULL) {
			curr_match = curr_match->next;
		    }
		    matchp new = (matchp)malloc(sizeof(matches));
		    m_malloc++;
		    new->index=index;
		    new->ntext = t;
		    new->strand = 1;
		    new->next = NULL;
		    curr_match->next = new;
		    nmatch[curr_id] = nmatch[curr_id]+1;
		    //printf("else, id=%d, index=%d, text=%d\n", curr_id, index, t);
		}*/
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
		//if(matcher[curr_id] == NULL) {
		      /*matchp new = (matchp)malloc(sizeof(matches));
		      m_malloc++;
		      new->index=index;
		      new->ntext = t;
		      new->strand = 1;
		      new->next = NULL;
		      matcher[curr_id] = new;*/
		      ntext[curr_id] = t;
		      nindex[curr_id] = index+1;
		      strand[curr_id] = -1;
		      nmatch[curr_id] = nmatch[curr_id]+1;
		      //printf("id=%d, index=%d, text=%d\n", curr_id, index, t);
		//}
		/*else if (matcher[curr_id] != NULL){
		    matchp curr_match = matcher[curr_id];
		    while(curr_match->next != NULL) {
			curr_match = curr_match->next;
		    }
		    matchp new = (matchp)malloc(sizeof(matches));
		    m_malloc++;
		    new->index=index;
		    new->ntext = t;
		    new->strand = 1;
		    new->next = NULL;
		    curr_match->next = new;
		    nmatch[curr_id] = nmatch[curr_id]+1;
		    //printf("else, id=%d, index=%d, text=%d\n", curr_id, index, t);
		}*/
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

	//printf("i=%d\n", i);
	
     // 
	



 kseq_destroy(seq);
	gzclose(fp);

  }
  destroy(root);

  /*while(queue->next != NULL) {
    listp currel = queue->next;
   // free(queue->entry);
    free(queue);
    l_free++;
    queue = currel;
  }
  free(queue);*/
    }

  Rprintf("\n");

    SEXP result, word_ind, word_t, wnames, word_strd;
    PROTECT(result = NEW_LIST(3));
	PROTECT(wnames = NEW_CHARACTER(3));
	SET_STRING_ELT(wnames, 0, mkChar("text"));
	SET_STRING_ELT(wnames, 1, mkChar("index"));
	SET_STRING_ELT(wnames, 2, mkChar("strand"));
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
    UNPROTECT(4);
    
	
/*for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    if(matcher[i] != NULL){
         matchp curr = matcher[i];
      while(curr->next != NULL) {
	matchp nextmatch = curr->next;
	free(curr);
	m_free++;
	curr = nextmatch;
      }
      free(curr);
      m_free++;
    }
   
  }*/
  free(nmatch);
  free(nindex);
   free(ntext);
   free(strand);
   /* printf("m_malloc: %d\n", m_malloc);
    printf("t_malloc: %d\n", t_malloc);
    printf("l_malloc: %d\n", l_malloc);
    printf("malloc_h: %d\n", h_malloc);
    printf("free: %d\n", c_free);
    printf("m_free: %d\n", m_free);
     printf("h_free: %d\n", h_free);
     printf("l_free: %d\n", l_free);
    
    printf("num_nodes: %d\n", num_nodes);*/
    
    return result;
}







// R CMD SHLIB ahocorasick.c
