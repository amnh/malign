/*
Copyright 1992 Ward Wheeler all rights reserved
change back gaps at least as an option
but need to remember the leading and trailing
could change cost or something (to 1/2) perhaps
*/
#include "align3.h"

/*
#define SETME 		diag=(values->gap_cost * (vbs1[4]|vbs2[4]));\
									down=vbs1[5];\
			 						right=vbs2[5];\

#define DOME 			diag+=(values->change_cost * ((vbs1[2]|vbs2[2])+(vbs1[1]|vbs2[1])+(vbs1[0]|vbs2[0])+(vbs1[3]|vbs2[3])-1));\
																		
*/

alignment *nw_real_new_opt (a1, a2,values)
     alignment *a1, *a2;
     parameters *values;
{
  cell **m;
  int i, j, k, l, n, down, right, diag, minimum, a1_offs, a2_offs;
  char holder,had_to_go_back;
  int net_gap,hold_i,hold_j;
  alignment *a;
  register int vbs1, vbs2, vbs1h,vbs2h;
  register cell *pmi,*pmi1;
  char rand_marker,direc_holder;
  int gaps_test=0;
  char *s1, *s2;
  
  GOMAC
    
    say ("allocating")
      if (values->show_mem) fprintf(stderr,"Allocating %d bytes\n",(1+a2->length)*(1+a1->length)*sizeof(cell));
  m = allocate_matrix (1+a1->length, 1+a2->length);
  if (values->iter) values->max_gap=(1+max(a1->length-a2->length,a2->length-a1->length));		
  if (values->pref_direc==RANDOM) rand_marker=1;
  else rand_marker=0;
  
  /*redo if to gappy*/
  had_to_go_back=0;
 redo_the_alignment:;
  
  if (values->iter) values->max_gap=min(values->max_gap,a2->length+1);
  
  if (values->max_gap!=MAX_SEQUENCE_SIZE) {		
    for (i=values->max_gap+1;((i<=a1->length) && (i-values->max_gap-1)<=a2->length);i++) {
      m[i][i-values->max_gap-1].score=BARRIER_COST;
      m[i][i-values->max_gap-1].direction=3;
    }
    for (i=0;((i<=a1->length) && ((i+values->max_gap+1)<=a2->length));i++) {
      m[i][i+values->max_gap+1].score=BARRIER_COST;	
      m[i][i+values->max_gap+1].direction=3;	
    }
  }	
  
  say ("filling in")
    m[0][0].direction = DIAGONAL;
  m[0][0].score = 0;
  m[0][0].length = 0;
  /*the vbs-1 and vbs2 -1 because of teh -tsring alocation*/
  /*leading gaps*/
  for (i = 1; i < min(a1->length+1,values->max_gap+1); i++) {
    m[i][0].direction = DOWN;
    m[i][0].score = ((values->gap_cost-values->leading_gap_cost)*(!(values->b_string1[0][0]&GAP_BIT))) + m[i-1][0].score; /*fix for leading gaps*/;
    m[i][0].length = 1 + m[i-1][0].length;
  }
  pmi=m[0];
  for (j = 1; j < min(a2->length+1,values->max_gap+1); j++) {
    pmi[j].direction = RIGHT;
    pmi[j].score=((values->gap_cost-values->leading_gap_cost)*(!(values->b_string2[0][0]&GAP_BIT))) + pmi[j-1].score; 
    pmi[j].length = 1 + pmi[j-1].length;
  }
  
  say ("main loop")
    
    if ((values->extra_adjustment==-32000) && (values->coding==-32000)) {
      for (i = 1; i< min(a1->length,values->max_gap+1); i++) {
	pmi = m[i];
	pmi1= m[i-1];
	vbs1=values->b_string1[0][i-1];
	for (j = 1; j< min(a2->length,(i+values->max_gap+1)); j++) {
	  vbs2=values->b_string2[0][j-1];
          SETME_NEW_OPT
	    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
	    else if (values->ttr) {DOMENOW_NEW_OPT}
	    else {DOMELATER_NEW_OPT}
          ADDEM
	    CHECKUM
	}
      }
      for (i = values->max_gap+1; i < a1->length; i++) {
	pmi = m[i];
	pmi1= m[i-1];
	vbs1=values->b_string1[0][i-1];
	for (j = i-values->max_gap; j < min(a2->length,(i+values->max_gap+1)); j++) {
	  vbs2=values->b_string2[0][j-1];
          SETME_NEW_OPT
	    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
	    else if (values->ttr) {DOMENOW_NEW_OPT}
	    else {DOMELATER_NEW_OPT}
          ADDEM
	    CHECKUM
	}
      }
    }
    else if (values->extra_adjustment!=-32000) {
      for (i = 1; i< min(a1->length,values->max_gap+1); i++) {
	pmi = m[i];
	pmi1= m[i-1];
	vbs1=values->b_string1[0][i-1];
	for (j = 1; j< min(a2->length,(i+values->max_gap+1)); j++) {
	  vbs2=values->b_string2[0][j-1];
          SETME_NEW_OPT
	    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
	    else if (values->ttr) {DOMENOW_NEW_OPT}
	    else {DOMELATER_NEW_OPT}
          ADDEM
	    if (pmi1[j].direction == DOWN) down -= values->extra_adjustment;
	  if (pmi[j-1].direction == RIGHT) right -= values->extra_adjustment;
          CHECKUM
	}
      }
      for (i = values->max_gap+1; i < a1->length; i++) {
	pmi = m[i];
	pmi1= m[i-1];
	vbs1=values->b_string1[0][i-1];
	for (j = i-values->max_gap; j < min(a2->length,(i+values->max_gap+1)); j++) {
	  vbs2=values->b_string2[0][j-1];
          SETME_NEW_OPT
	    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
	    else if (values->ttr) {DOMENOW_NEW_OPT}
	    else {DOMELATER_NEW_OPT}
	  ADDEM
	    if (pmi1[j].direction == DOWN) down -= values->extra_adjustment;
	  if (pmi[j-1].direction == RIGHT) right -= values->extra_adjustment;
          CHECKUM
	}
      }
    }
    else if (values->coding!=-32000) {
      for (i = 1; i< min(a1->length,values->max_gap+1); i++) {
	pmi = m[i];
	pmi1= m[i-1];
	vbs1=values->b_string1[0][i-1];
	for (j = 1; j< min(a2->length,(i+values->max_gap+1)); j++) {
	  vbs2=values->b_string2[0][j-1];
          SETME_NEW_OPT
	    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
	    else if (values->ttr) {DOMENOW_NEW_OPT}
	    else {DOMELATER_NEW_OPT}
	  ADDEM
	    /*fixes for gap lengths*/
	    holder=1;
	  for (k=i-1;k>=0;k--) {
	    if (m[k][j].direction==DOWN) ++holder;
	    else break;
	  }
	  if ((holder%3)==0) down -= values->coding;
	  holder=1;
	  for (k=i-1;k>=0;k--) {
	    if (pmi[k].direction==RIGHT) ++holder;
	    else break;
	  }
	  if ((holder%3)==0) right -= values->coding;
          CHECKUM
	}
      }
      for (i = values->max_gap+1; i < a1->length; i++) {
	pmi = m[i];
	pmi1= m[i-1];
	vbs1=values->b_string1[0][i-1];
	for (j = i-values->max_gap; j < min(a2->length,(i+values->max_gap+1)); j++) {
	  vbs2=values->b_string2[0][j-1];
          SETME_NEW_OPT
	    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
	    else if (values->ttr) {DOMENOW_NEW_OPT}
	    else {DOMELATER_NEW_OPT}
          ADDEM
	    /*fixes for gap lengths*/
	    holder=1;
	  for (k=i-1;k>=0;k--) {
	    if (m[k][j].direction==DOWN) ++holder;
	    else break;
	  }
	  if ((holder%3)==0) down -= values->coding;
	  holder=1;
	  for (k=i-1;k>=0;k--) {
	    if (pmi[k].direction==RIGHT) ++holder;
	    else break;
	  }
	  if ((holder%3)==0) right -= values->coding;
          CHECKUM
	}
      }
    }
  /* trailing gap fix - moved to optimize i = a1->length;*/
  pmi = m[a1->length];
  pmi1= m[a1->length-1];
  vbs1=values->b_string1[0][a1->length-1];
  for (j = max(1,a1->length-values->max_gap); j < a2->length; j++) {
    vbs2=values->b_string2[0][j-1];
    SETME_NEW_OPT
      if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
      else if (values->ttr) {DOMENOW_NEW_OPT}
      else {DOMELATER_NEW_OPT}
    if (right==values->gap_cost) right -=values->trailing_gap_cost;
    ADDEM
      CHECKUM
  }
  /*j = a2->length;*/
  vbs2=values->b_string2[0][a2->length-1];
  for (i = max(1,a2->length-values->max_gap); i < a1->length; i++) {
    vbs1=values->b_string1[0][i-1];
    SETME_NEW_OPT
      if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
      else if (values->ttr) {DOMENOW_NEW_OPT}
      else {DOMELATER_NEW_OPT}
    if (down==values->gap_cost) down -= values->trailing_gap_cost;
    diag	+= m[i-1][a2->length-1].score;
    down	+= m[i-1][a2->length].score;
    right += m[i][a2->length-1].score;
    minimum = min (min (diag, down), right);
    m[i][a2->length].score = minimum;
    if (diag == minimum) {
      m[i][a2->length].direction = DIAGONAL;
      m[i][a2->length].length = 1 + m[i-1][a2->length-1].length;
    }
    else if (down == minimum) {
      m[i][a2->length].direction = DOWN;
      m[i][a2->length].length = 1 + m[i-1][a2->length].length;
    }
    else {
      m[i][a2->length].direction = RIGHT;
      m[i][a2->length].length = 1 + m[i][a2->length-1].length;
    }
  }
  
  /*	j = a2->length;
	i = a1->length; 
	*/
  vbs1=values->b_string1[0][a1->length-1];
  vbs2=values->b_string2[0][a2->length-1];
  SETME_NEW_OPT
    if ((!values->delta) && (!values->ttr)) {DOME_NEW_OPT}
    else if (values->ttr) {DOMENOW_NEW_OPT}
    else {DOMELATER_NEW_OPT}
  if (right==values->gap_cost) right -= values->trailing_gap_cost;
  if (down==values->gap_cost)  down -= values->trailing_gap_cost;	
  diag	+= m[a1->length-1][a2->length-1].score;
  down	+= m[a1->length-1][a2->length].score;
  right += m[a1->length][a2->length-1].score;
  
  minimum = min (min (diag, down), right);
  m[a1->length][a2->length].score = minimum;
  if (diag == minimum) {
    m[a1->length][a2->length].direction = DIAGONAL;
    m[a1->length][a2->length].length = 1 + m[a1->length-1][a2->length-1].length;
  }
  else if (down == minimum) {
    m[a1->length][a2->length].direction = DOWN;
    m[a1->length][a2->length].length = 1 + m[a1->length-1][a2->length].length;
  }
  else {
    m[a1->length][a2->length].direction = RIGHT;
    m[a1->length][a2->length].length = 1 + m[a1->length][a2->length-1].length;
  }
  say ("trace back")
  i = a1->length;
  j = a2->length;
  k = m[i][j].length;
  n = a1->n_seqs + a2->n_seqs;
  s1=(char *)malloc((1+k)*sizeof(char));
  assert((int) s1);
  s2=(char *)malloc((1+k)*sizeof(char));
  assert((int) s2);
  s1[k]='\0';
  s2[k]='\0';

   net_gap=hold_i=hold_j=0;
  while (k > 0) {
    k--;
    switch (m[i][j].direction) {
    case DIAGONAL:
      s1[k] = a1->s[0][i-1];
      s2[k] = a2->s[0][j-1];
      i--;
      j--;
      /*if (!(a1->s[1][k]&a->s[0][k])) ++gaps_test;*/
      break;
    case DOWN:
      s1[k] = a1->s[0][i-1];
      s2[k] = '-';
      i--;
      ++hold_i;
      break;
    case RIGHT:
      s1[k] = '-';
      s2[k] = a2->s[0][j-1];
      j--;
      ++hold_j;
      break;
    default:
      fprintf (stderr, "Bad direction in matrix: %d\n", m[i][j].direction);
      exit (-1);
    }
  }
  
  if (values->iter) {	
    net_gap=max(hold_i,hold_j);
    if (net_gap >= values->max_gap) {
      values->max_gap=(2*(net_gap+1));
      free(s1);
      free(s2);
      had_to_go_back=1;	
      goto redo_the_alignment;
    }
    else values->max_gap=net_gap+1;
  }
  
 i = a1->length;
  j = a2->length;
  k = m[i][j].length;
  n = a1->n_seqs + a2->n_seqs;
  a = allocate_alignment_and_strings (n, k+1);
  a->score = m[i][j].score;
  a->length = m[i][j].length;
  if (strcmp (a1->name, a2->name) > 0) {
    a1_offs = a2->n_seqs;
    a2_offs = 0;
    a->name = pair_names (a2, a1);
    for (l=0;l<a2->n_seqs;l++) {
      a->taxon_name[l]=(char *)malloc((1+(strlen(a2->taxon_name[l])))*sizeof(char));
      assert((int)a->taxon_name[l]);
      strcpy(a->taxon_name[l],a2->taxon_name[l]);
    }
    for (l=a2->n_seqs;l<(a2->n_seqs+a1->n_seqs);l++) {
      a->taxon_name[l]=(char *)malloc((1+(strlen(a1->taxon_name[l-(a2->n_seqs)])))*sizeof(char));
      assert((int)a->taxon_name[l]);
      strcpy(a->taxon_name[l],a1->taxon_name[l-(a2->n_seqs)]);
    }
  }
  else {
    a1_offs = 0;
    a2_offs = a1->n_seqs;
    a->name = pair_names (a1, a2);
    for (l=0;l<a1->n_seqs;l++) {
      a->taxon_name[l]=(char *)malloc((1+(strlen(a1->taxon_name[l])))*sizeof(char));
      assert((int)a->taxon_name[l]);
      strcpy(a->taxon_name[l],a1->taxon_name[l]);
    }
    for (l=a1->n_seqs;l<(a1->n_seqs+a2->n_seqs);l++) {
      a->taxon_name[l]=(char *)malloc((1+(strlen(a2->taxon_name[l-(a1->n_seqs)])))*sizeof(char));
      assert((int)a->taxon_name[l]);
      strcpy(a->taxon_name[l],a2->taxon_name[l-(a1->n_seqs)]);
    }
  }
  strcpy(a->s[a1_offs],s1);
  strcpy(a->s[a2_offs],s2);
  
say ("freeing")
  free(s1);
  free(s2);
  free_matrix (m, 1+a1->length);
  free(m);
  if (values->phylo_score>0) {
    if ((!values->in_some) || (a->n_seqs==values->actual_num_sequences)) {
      if (!compare_aligns(a,values->previous)) a->score=values->previous->score;
      else {
	a->score = cladogram(a,values);
	if (values->previous) values->previous=dump_align(values->previous);
	values->previous=make_align(a);
      }
      /*fprintf(stderr,"l=%d",a->score);*/
    }
  }
  say ("returning")
    /*printf("gaps %d changes %d ",hold_i+hold_j,gaps_test);*/
    /*	printf("%s %s\n",a1->name, a2->name);
	for (i=0;i<a1->length;i++) printf("(%c)%2d ",a1->s[0][i],values->b_string1[0][i]);
	printf("\n");
	for (i=0;i<a2->length;i++) printf("(%c)%2d ",a2->s[0][i],values->b_string2[0][i]);
	printf("\n");*/
    return a;
  
}
