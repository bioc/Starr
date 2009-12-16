#include <stdio.h>
#include <R.h>
#include <Rdefines.h>




/**
  Definition of data structure.
**/

int DEBUG = 0;
int num_nodes = 1;

/**
  A	T	C	G
  0	1	2	3
**/


struct hits{
   int pdict_id;          /* Storing marks of a node */
    struct hits *next;    /* Storing next address */
};
/*****  Redefining struct list as node  *****/
typedef struct hits *hitsp;
typedef struct hits hits;


struct node{
char info;
hitsp dict_ids;
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


void preorder(nodep p){
  if(p!=NULL){
     printf("%c ",p->info);
     preorder(p->links[0]);
     preorder(p->links[1]);
    preorder(p->links[2]);
      preorder(p->links[3]);
 }

}


int getPos(char c) {
  int pos;
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


int get_complement(char c) {
  char comp;
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


nodep insert_node(nodep current, char c) {
    nodep newnode;
    newnode=(nodep)malloc(sizeof(node));
     int i;
  for(i=0; i<(sizeof(newnode->links)/4); i++) {
    newnode->links[i] = NULL;
    //printf("Number %d\n", i);
  }
newnode->dict_ids = NULL;
newnode->flink = NULL;
  int pos = getPos(c);
  
      if(current->links[pos] == NULL) {
      current->links[pos] = newnode;
      newnode->info = c;
      current=newnode;
      num_nodes++;
      }
      else {
	current=current->links[pos];
      }
 
    return current;
}


nodep init_tree() {
  nodep root = (nodep)malloc(sizeof(node));
  root->info = '-';
  root->flink=NULL;
  root->dict_ids=NULL;

  int i;
  for(i=0; i<(sizeof(root->links)/4); i++) {
    root->links[i] = NULL;
    //printf("Number %d\n", i);
  }
 
  return(root);
}

void append_hit(hitsp hit, int id) {
  hitsp new_hit = (hitsp)malloc(sizeof(hits));
  new_hit->pdict_id = id;
  new_hit->next=hit;
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
  
  if(currnode->dict_ids == NULL) {
    hitsp new_hit = (hitsp)malloc(sizeof(hits));
    new_hit->pdict_id = id;
    new_hit->next=NULL;
    currnode->dict_ids = new_hit;
  }
  else {
      hitsp currhit = currnode->dict_ids;
     while(currhit->next != NULL) {
	currhit = currhit->next;
      }
      hitsp new_hit = (hitsp)malloc(sizeof(hits));
    new_hit->pdict_id = id;
    new_hit->next=NULL;
    currhit->next = new_hit;
     //append_hit(currnode->dict_ids, id);
  }

  if(DEBUG) {
    printf("pdict_id->%d\n", id);
  }
}


listp init_list(nodep node) {
  listp head = (listp)malloc(sizeof(list));
  head->entry = node;
  head->next = NULL;
 
  return(head);
}

// returns new tail of the list
listp add(listp tail, nodep node) {
  listp new = (listp)malloc(sizeof(list));
  new->entry = node;
  new->next = NULL;
  tail->next = new;

  return(new);
}


listp breadth_first_search(nodep root) {
  listp queue = init_list(root);
  listp curr = queue;
  listp tail = queue;
  int j = 1;
  nodep currnode = curr->entry;
  int i;

  // initialize flinks in first level
  for(i=0; i<=3; i++) {
      if(currnode->links[i] != NULL) {
	tail = add(tail, currnode->links[i]);
	nodep fnode = currnode->links[i];
	fnode->flink = root;
	if(DEBUG) {
	  printf("flink: %c -> %c\n",fnode->info, root->info);
	}
	//nodep a = currnode->links[i];
	//printf("adding: %c, %d, %d\n", a->info, i, j);
      }
    }
    j++;
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
	int pos = getPos(fnode->info);
	int stop = 0;
	
	while(stop != 1) {
	  // only follow first flink
	  int first_flink = 1;
	  if(a->links[pos] != NULL) {
	    a = a->links[pos];
	    
	    fnode->flink = a;
	    if(DEBUG) {
	      printf("flink: %c -> %c\n",fnode->info, a->info);
	    }
	    stop=1;
	  }
	  else if(a->info =='-') {
	    fnode->flink = a;
	    if(DEBUG) {
	      printf("flink: %c -> %c\n",fnode->info, a->info);
	    }
	    stop = 1;
	  }
	  else {
	    if(DEBUG) {
	      printf("follow flink from %c ", a->info);
	    }
	    a = a->flink;
	    first_flink = 0;
	    if(DEBUG) {
	      printf("to %c\n", a->info);
	    }
	  }
	  if(a->dict_ids != NULL & a->info != '-' & first_flink == 1) {
	hitsp ahits = a->dict_ids;
	if(fnode->dict_ids == NULL) {
	  fnode->dict_ids = ahits;
	  //printf("a->info: %c", a->info);
	  if(DEBUG) {
	    hitsp fhits = fnode->dict_ids;
	    printf("Creating hit entries:");
	    while(fhits->next != NULL) {
	    printf(" %d", fhits->pdict_id);
	    fhits=fhits->next;
	  }
	  printf(" %d\n", fhits->pdict_id);
	  }
	}
	else {
	  hitsp fhits = fnode->dict_ids;
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
	  }
	  }
	  }
	}
	
      }
      
    }

    
    j++;
    curr = curr->next;
    if(curr != NULL) {
      currnode = curr->entry;
    }
  }
  return(queue);
}


char getChar(int c) {
  char pos;
    switch(c)
{
   case 0:
     pos = 'A';
  break;

  case 1:
     pos = 'T';
  break;

  case 2:
     pos = 'C';
  break;

case 3:
     pos = 'G';
  break;
}
return(pos);
}



nodep make_tree(char **dict, int *wcount_p, int *numnodes_p) {
  nodep root = init_tree();

 int wcount = wcount_p[0];
  int *numnodes = numnodes_p;
  
  int i,j;
  for(i=0; i<*wcount_p; i++) {
          insert_word(root, dict[i], i);
  }
  breadth_first_search(root);
  *numnodes_p = num_nodes;
  //lengths[2] = 3;
  // printf("number of nodes: %d\n", num_nodes);
  
  return(root);
}


/*void find_ac(char **dict, int *wcount_p, int *numnodes_p, char **text, int *num_text, int matches[]) {
    nodep root = make_tree(dict, wcount_p, numnodes_p);
    
    int index, t;
    for(t=0; t<num_text[0]; t++) {
      nodep currnode = root;
      index = 0;
      while(text[t][index] != '\0') {
	int pos = getPos(text[t][index]);
	if(DEBUG) {
	  printf("char=%c\n", text[t][index]);
	}
	int stop = 0;
	while(currnode->links[pos] == NULL & stop != 1) {
	    if(currnode->info == '-') {
	      stop = 1;
	    }
	    else {
	      currnode = currnode->flink;
	    }
	}
	if(currnode->links[pos] != NULL) {
	  currnode = currnode->links[pos];
	  if(currnode->pdict_id != -1) {
	    if(DEBUG) {
	      printf("id=%d found at index=%d in text=%d\n", currnode->pdict_id, index, t);
	    }
	    if(matches[currnode->pdict_id+(t*wcount_p[0])] > 0) {
	      matches[currnode->pdict_id+(t*wcount_p[0])] = -1;     
	    }
	    else if (matches[currnode->pdict_id+(t*wcount_p[0])] == 0){
	      matches[currnode->pdict_id+(t*wcount_p[0])] = index;
	    }
	  }
	}
	index++;
      }
    }
}*/



SEXP find_ac2(SEXP dict, SEXP wcount_p, SEXP text, SEXP num_text, SEXP complementary, SEXP reverse, SEXP reverse_complementary) {
  matchp matcher[INTEGER(wcount_p)[0]];
  int nmatch[INTEGER(wcount_p)[0]];
  int i;
  for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    matcher[i] = NULL;
    nmatch[i] = 0;
  }

  nodep root = init_tree();

  
  for(i=0; i<INTEGER(wcount_p)[0]; i++) {
          insert_word(root, CHAR (STRING_ELT (dict, i)), i);
  }
  listp queue = breadth_first_search(root);
  
 // printf("number of nodes: %d\nSearching...", num_nodes);
  
   int comp = INTEGER(complementary)[0];
    int revcomp = INTEGER(reverse_complementary)[0];
int rev = INTEGER(reverse)[0];
    int index, t;
    for(t=0; t<INTEGER(num_text)[0]; t++) {
      nodep currnode = root;
      nodep lastNorm = currnode;
      nodep lastComp = currnode;
      index = 0;
      while(CHAR (STRING_ELT (text, t))[index] != '\0') {
	int z;
	for(z=0; z<1+comp; z++) {
     	int pos = getPos(CHAR (STRING_ELT (text, t))[index]);
	currnode = lastNorm;
	if(z == 1) {
	  pos = getPos(get_complement(CHAR (STRING_ELT (text, t))[index]));
	  currnode = lastComp;
	}
	if(DEBUG) {
	  printf("char=%c\n", get_complement(CHAR (STRING_ELT (text, t))[index]));
	}
	int stop = 0;
	while(currnode->links[pos] == NULL & stop != 1) {
	    if(currnode->info == '-') {
	      stop = 1;
	    }
	    else {
	      currnode = currnode->flink;
	    }
	}
	if(currnode->links[pos] != NULL) {
	  currnode = currnode->links[pos];
	  if(currnode->dict_ids != NULL) {
	    if(DEBUG) {
	      printf("id=");
	    }
	    hitsp curr_hits = currnode->dict_ids;
	    while(curr_hits != NULL) {
	      int curr_id=curr_hits->pdict_id; 
	      if(matcher[curr_id] == NULL) {
		    matchp new = (matchp)malloc(sizeof(matches));
		    new->index=index;
		    new->ntext = t;
		    new->strand = 1;
		    new->next = NULL;
		    matcher[curr_id] = new;
		    nmatch[curr_id] = nmatch[curr_id]+1;
		    //printf("id=%d, index=%d, text=%d\n", curr_id, index, t);
	      }
	      else if (matcher[curr_id] != NULL){
		  matchp curr_match = matcher[curr_id];
		  while(curr_match->next != NULL) {
		      curr_match = curr_match->next;
		  }
		  matchp new = (matchp)malloc(sizeof(matches));
		  new->index=index;
		  new->ntext = t;
		  new->strand = 1;
		  new->next = NULL;
		  curr_match->next = new;
		  nmatch[curr_id] = nmatch[curr_id]+1;
		  //printf("else, id=%d, index=%d, text=%d\n", curr_id, index, t);
	      }
	      curr_hits=curr_hits->next;
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

      if(rev == 1 | revcomp == 1) {
      currnode = root;
      lastNorm = currnode;
      lastComp = currnode;
int from = 0;
      if(rev == 0) { from=1;}
     index--;
      while(index != -1) {
	int z;
	
	
	for(z=from; z<1+revcomp; z++) {
     	int pos = getPos(CHAR (STRING_ELT (text, t))[index]);
	currnode = lastNorm;
	if(z == 1) {
	  pos = getPos(get_complement(CHAR (STRING_ELT (text, t))[index]));
	  currnode = lastComp;
	}
	if(DEBUG) {
	  printf("char=%c\n", CHAR (STRING_ELT (text, t))[index]);
	}
	int stop = 0;
	while(currnode->links[pos] == NULL & stop != 1) {
	    if(currnode->info == '-') {
	      stop = 1;
	    }
	    else {
	      currnode = currnode->flink;
	    }
	}
	if(currnode->links[pos] != NULL) {
	  currnode = currnode->links[pos];
	  if(currnode->dict_ids != NULL) {
	    if(DEBUG) {
	      printf("id=");
	    }
	    hitsp curr_hits = currnode->dict_ids;
	    while(curr_hits != NULL) {
	      int curr_id=curr_hits->pdict_id; 
	      if(matcher[curr_id] == NULL) {
		    matchp new = (matchp)malloc(sizeof(matches));
		    new->index=index;
		    new->ntext = t;
		    new->strand = -1;
		    new->next = NULL;
		    matcher[curr_id] = new;
		    nmatch[curr_id] = nmatch[curr_id]+1;
		    //printf("id=%d, index=%d, text=%d\n", curr_id, index, t);
	      }
	      else if (matcher[curr_id] != NULL){
		  matchp curr_match = matcher[curr_id];
		  while(curr_match->next != NULL) {
		      curr_match = curr_match->next;
		  }
		  matchp new = (matchp)malloc(sizeof(matches));
		  new->index=index;
		  new->ntext = t;
		  new->strand = -1;
		  new->next = NULL;
		  curr_match->next = new;
		  nmatch[curr_id] = nmatch[curr_id]+1;
		  //printf("else, id=%d, index=%d, text=%d\n", curr_id, index, t);
	      }
	      curr_hits=curr_hits->next;
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

    }

    

    SEXP result, curr_word, word_ind, word_t, wnames, word_strd;
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
	matchp curr = matcher[i];
	INTEGER_POINTER(word_t)[i] = curr->ntext;
	INTEGER_POINTER(word_ind)[i] = curr->index;
	INTEGER_POINTER(word_strd)[i] = curr->strand;
	 
	
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
    // free all nodes
    while(queue->next != NULL) {
	nodep current_node = queue->entry;
	
	
	if(current_node->dict_ids != NULL) {
		hitsp currHit = current_node->dict_ids;
		while(currHit->next != NULL) {
		  hitsp nextHit = currHit->next;
		  free(currHit);
		  currHit = nextHit;
		}
	}
	
	free(current_node);
	listp before = queue;
	free(before);
	queue = queue->next;
    }
    free(queue->entry);
	
for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    free(matcher[i]);
  }

    return result;
}



SEXP naive_search(SEXP dict, SEXP wcount_p, SEXP text, SEXP num_text) {
  matchp matcher[INTEGER(wcount_p)[0]];
  int nmatch[INTEGER(wcount_p)[0]];
  int i;
  for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    matcher[i] = NULL;
    nmatch[i] = 0;
  }

 
   
    int index, t;
    for(t=0; t<INTEGER(num_text)[0]; t++) {
      int textlen=0;
      const char *t_curr = CHAR (STRING_ELT (text, t));
      while(t_curr[textlen] != '\0') {
	textlen++;
      }
      for(i=0; i<INTEGER(wcount_p)[0]; i++) {
	int wlen=0;
	const char *w_curr=CHAR (STRING_ELT (dict, i));
	while(w_curr[wlen] != '\0') {
	  //printf("w1=%c\n", w_curr[wlen]);
	  wlen++;
	}
	int q,r;
	for(q=0; q<=textlen-wlen; q++) {
	  for(r=0; r<wlen; r++) {
	    
	     if(t_curr[q+r] != w_curr[r]) {
		r=wlen;
	     }
	      if(t_curr[q+r] == w_curr[r] & r == wlen-1) {
		 //printf("t=%c, w=%c\n", t_curr[q+r], w_curr[r]);
		  if(matcher[i] == NULL) {
		    matchp new = (matchp)malloc(sizeof(matches));
		    new->index=q+r;
		    new->ntext = t;
		    new->next = NULL;
		    matcher[i] = new;
		    nmatch[i] = nmatch[i]+1;
		    //printf("id=%d, index=%d, text=%d\n", i, index, t);
	      }
	      else if (matcher[i] != NULL){
		  matchp curr_match = matcher[i];
		  while(curr_match->next != NULL) {
		      curr_match = curr_match->next;
		  }
		  matchp new = (matchp)malloc(sizeof(matches));
		  new->index=q+r;
		  new->ntext = t;
		  new->next = NULL;
		  curr_match->next = new;
		  nmatch[i] = nmatch[i]+1;
		 // printf("else, id=%d, index=%d, text=%d\n", i, index, t);
	      }
	      }
	  }
	}
      }
      
    }

    SEXP result, curr_word, word_ind, word_t, wnames;
    PROTECT(result = NEW_LIST(INTEGER(wcount_p)[0]));
    
    for(i=0; i<INTEGER(wcount_p)[0]; i++) {
    
      if(nmatch[i] > 0) {
	PROTECT(curr_word = NEW_LIST(2));
	PROTECT(wnames = NEW_CHARACTER(2));
	SET_STRING_ELT(wnames, 0, mkChar("text"));
	SET_STRING_ELT(wnames, 1, mkChar("index"));
	SET_NAMES(curr_word, wnames);
		  
	  PROTECT(word_ind = NEW_INTEGER(nmatch[i]));
	  PROTECT(word_t = NEW_INTEGER(nmatch[i]));
	 
	int pos=0;

	matchp curr = matcher[i];
	while(curr != NULL) {
	INTEGER_POINTER(word_t)[pos] = curr->ntext;
	INTEGER_POINTER(word_ind)[pos] = curr->index;
	 
	pos++;
	curr=curr->next;
	}
	  
	  SET_ELEMENT(curr_word, 0, word_t);
	  
	  SET_ELEMENT(curr_word, 1, word_ind);
	  SET_ELEMENT(result, i, curr_word);
	  UNPROTECT(4);
      }
      else {
	 PROTECT(curr_word = NEW_LIST(2));
	
	PROTECT(wnames = NEW_CHARACTER(2));
	SET_STRING_ELT(wnames, 0, mkChar("text"));
	SET_STRING_ELT(wnames, 1, mkChar("index"));
	SET_NAMES(curr_word, wnames);
	//UNPROTECT(1);
	
	  PROTECT(word_t = NEW_INTEGER(1));
	 INTEGER_POINTER(word_t)[0] = -1;
	 SET_ELEMENT(curr_word, 0, word_t);
	 //UNPROTECT(1);
	
	  PROTECT(word_ind = NEW_INTEGER(1));
	 INTEGER_POINTER(word_ind)[0] = -1;
	 SET_ELEMENT(curr_word, 1, word_ind);
	 //UNPROTECT(1);

	  SET_ELEMENT(result, i, curr_word);
	  UNPROTECT(4);
    
      }
    }

    UNPROTECT(1);
    return result;
}




/*SEXP ac2Rlist(SEXP dict, SEXP nentry, SEXP lengths)
{
	const char *dict1 = CHAR (STRING_ELT (dict, 1));
	printf("dict[0]=%s\n", dict1);
	int array1 = INTEGER(lengths)[0];
	printf("array1=%d\n", array1);
	nodep root = init_tree();
	int wcount = INTEGER(nentry)[0];
  
  int i,j;
  for(i=0; i<wcount; i++) {
          insert_word(root, CHAR (STRING_ELT (dict, i)), INTEGER(lengths)[i], i);
  }
  breadth_first_search(root);
  preorder(root);
  	printf("hier\n");
	SEXP r_tree, l_names, element, links, info;
	PROTECT(r_tree = NEW_LIST(num_nodes));
	
	listp queue = init_list(root);
  listp curr = queue;
  listp tail = queue;

  nodep currnode = curr->entry;
  // set the root element 
	PROTECT(element = NEW_LIST(4));
	PROTECT(l_names = NEW_CHARACTER(4));
	SET_STRING_ELT(l_names, 0, mkChar("info"));
	SET_STRING_ELT(l_names, 1, mkChar("id"));
	SET_STRING_ELT(l_names, 2, mkChar("children"));
	SET_STRING_ELT(l_names, 3, mkChar("flink"));
	SET_NAMES(element, l_names);
	UNPROTECT(1);

	PROTECT(info = NEW_INTEGER(1));
       INTEGER_POINTER(info)[0] = -1;
       SET_ELEMENT(element, 0, info);
	UNPROTECT(1);
	SET_ELEMENT(r_tree, 0, element);
	UNPROTECT(1);
	

 int entry_in_Rtree = 1;
  
 // initialize flinks in first level
  nchildren = 0;
  for(i=0; i<=3; i++) {
      if(currnode->links[i] != NULL) {
	nchildren++;
	tail = add(tail, currnode->links[i]);
	nodep fnode = currnode->links[i];
	fnode->flink = root;
	if(DEBUG) {
	  printf("flink: %c -> %c\n",fnode->info, root->info);
	}
	PROTECT(element = NEW_LIST(4));
	PROTECT(l_names = NEW_CHARACTER(4));
	SET_STRING_ELT(l_names, 0, mkChar("info"));
	SET_STRING_ELT(l_names, 1, mkChar("id"));
	SET_STRING_ELT(l_names, 2, mkChar("children"));
	SET_STRING_ELT(l_names, 3, mkChar("flink"));
	SET_NAMES(element, l_names);
	UNPROTECT(1);

	PROTECT(info = NEW_INTEGER(1));
       INTEGER_POINTER(info)[0] = i;
       SET_ELEMENT(element, 0, info);
	UNPROTECT(1);
	SET_ELEMENT(r_tree, entry_in_Rtree++, element);
	UNPROTECT(1);
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
	int pos = getPos(fnode->info);
	int stop = 0;
	
	while(stop != 1) {
	  if(a->links[pos] != NULL) {
	    a = a->links[pos];
	    
	    fnode->flink = a;
	    if(DEBUG) {
	      printf("flink: %c -> %c\n",fnode->info, a->info);
	    }
	    stop=1;
	  }
	  else if(a->info =='-') {
	    fnode->flink = a;
	    if(DEBUG) {
	      printf("flink: %c -> %c\n",fnode->info, a->info);
	    }
	    stop = 1;
	  }
	  else {
	    if(DEBUG) {
	      printf("follow flink from %c ", a->info);
	    }
	    a = a->flink;
	    if(DEBUG) {
	      printf("to %c\n", a->info);
	    }
	  }
	}
	PROTECT(element = NEW_LIST(4));
	PROTECT(l_names = NEW_CHARACTER(4));
	SET_STRING_ELT(l_names, 0, mkChar("info"));
	SET_STRING_ELT(l_names, 1, mkChar("id"));
	SET_STRING_ELT(l_names, 2, mkChar("children"));
	SET_STRING_ELT(l_names, 3, mkChar("flink"));
	SET_NAMES(element, l_names);
	UNPROTECT(1);

	PROTECT(info = NEW_INTEGER(1));
       INTEGER_POINTER(info)[0] = getPos(fnode->info);
       SET_ELEMENT(element, 0, info);
	UNPROTECT(1);
	SET_ELEMENT(r_tree, entry_in_Rtree++, element);
	UNPROTECT(1);
      }
    }
    j++;
    curr = curr->next;
    if(curr != NULL) {
      currnode = curr->entry;
    }
  }
   

	UNPROTECT(1);
	return r_tree;
}*/

// R CMD SHLIB ahocorasick.c
