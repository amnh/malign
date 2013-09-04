/*
Copyright 1992 Ward Wheeler all rights reserved
change back gaps at least as an option
but need to remember the leading and trailing
could change cost or something (to 1/2) perhaps
*/
#include "align3.h"

alignment *nw_real (a1, a2,values)
alignment *a1, *a2;
parameters *values;
{
	cell **m;
	int i, j, k, l, n, down, right, diag, minimum, a1_offs, a2_offs;
	char holder,had_to_go_back;
	int net_gap,hold_i,hold_j;
	alignment *a;
	register cell *pmi,*pmi1;
	register int *vbs1, *vbs2;
	char rand_marker,direc_holder;

	GOMAC
	
		say ("allocating")
		if (values->show_mem) fprintf(stderr,"Allocating %d bytes\n",(1+a2->length)*(1+a1->length)*sizeof(cell));
		/*m = allocate_matrix (1+a1->length, 1+a2->length);*/
		m=(cell **)malloc((1+a1->length)*sizeof(cell *));
		assert((int) m);
		for (i=0;i<1+a1->length;i++) {
			m[i]=(cell *)malloc((1+a2->length)*sizeof(cell));
			assert((int) m[i]);
			}
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

		for (i = 1; i < min(a1->length+1,values->max_gap+1); i++) {
			m[i][0].direction = DOWN;
			vbs1=values->b_string1[i-1];
			if ((!values->delta) && (!values->ttr)) {
      	if (values->gap_must_cost) down=(values->change_cost * (vbs1[0]+vbs1[1]+vbs1[2]+vbs1[3]-1) + values->gap_cost * (vbs1[4]+1));
        else down=(values->change_cost * (vbs1[0]+vbs1[1]+vbs1[2]+vbs1[3]-1) + values->gap_cost);
        }
      else if (values->ttr) {
				if (!(vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3]-1)) {
        	if (values->gap_must_cost) down=(values->gap_cost * (vbs1[4]+1));
          else down=values->gap_cost;
          }
				else {
					down=(values->transition+values->transversion)*((vbs1[0]&vbs1[1])+(vbs1[0]&vbs1[3])+(vbs1[1]&vbs1[2])+(vbs1[2]&vbs1[3]));
					down+=(values->transition*((vbs1[0]&vbs1[2])+(vbs1[1]&vbs1[3])));
					if ((vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3])==3) {
						down << 1;
						down/=3;
						}
					else	if ((vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3])==4) down >> 1;
					if (values->gap_must_cost) down += (values->gap_cost * (vbs1[4]+1));
          else down += values->gap_cost;
					}
				}
			else {
				if (!(vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3]-1)) {
        	if (values->gap_must_cost) down= (values->gap_cost * (vbs1[4]+1));
          else down= values->gap_cost;
          }
				else {
					down+=(values->delta[0][1]*(vbs1[0]&vbs1[1]));
					down+=(values->delta[0][2]*(vbs1[0]&vbs1[2]));
					down+=(values->delta[0][3]*(vbs1[0]&vbs1[3]));
					down+=(values->delta[1][2]*(vbs1[1]&vbs1[2]));
					down+=(values->delta[1][3]*(vbs1[1]&vbs1[3]));
					down+=(values->delta[2][3]*(vbs1[2]&vbs1[3]));
					if ((vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3])==3) {
						down << 1;
						down/=3;
						}
					else	if ((vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3])==4) down >> 1;
					if (values->gap_must_cost) down += (values->gap_cost * (vbs1[4]+1));
          else down += values->gap_cost;
					}
				}
			m[i][0].score = down + m[i-1][0].score - values->leading_gap_cost; /*fix for leading gaps*/;
			m[i][0].length = 1 + m[i-1][0].length;
	 		}
		pmi=m[0];
	 	for (j = 1; j < min(a2->length+1,values->max_gap+1); j++) {
			pmi[j].direction = RIGHT;
			vbs2=values->b_string2[j-1];
			if ((!values->delta) && (!values->ttr)) {
      	if (values->gap_must_cost) right=(values->change_cost * (vbs2[0]+vbs2[1]+vbs2[2]+vbs2[3]-1) + values->gap_cost * (1+vbs2[4]));
       	else right=(values->change_cost * (vbs2[0]+vbs2[1]+vbs2[2]+vbs2[3]-1) + values->gap_cost);
        }
			else if (values->ttr) {
				if (!(vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3]-1)) {
        	if (values->gap_must_cost) right=(values->gap_cost * (1+vbs2[4]));
          else right=values->gap_cost;
          }
				else {
					right=((values->transition+values->transversion)*((vbs2[0]&vbs2[1])+(vbs2[0]&vbs2[3])+(vbs2[1]&vbs2[2])+(vbs2[2]&vbs2[3])));
					right+=(values->transition*((vbs2[0]&vbs2[2])+(vbs2[1]&vbs2[3])));
					if ((vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3])==3) {
						right << 1;
						right/=3;
						}
					else	if ((vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3])==4) right >> 1;
					if (values->gap_must_cost) right += (values->gap_cost * (1+vbs2[4]));
          else right += values->gap_cost;
					}
				}
			else {
				if (!(vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3]-1)) {
        	if (values->gap_must_cost) right= (values->gap_cost * (1+vbs2[4]));
          else right= values->gap_cost;
          }
				else {
					right+=(values->delta[0][1]*(vbs2[0]&vbs2[1]));
					right+=(values->delta[0][2]*(vbs2[0]&vbs2[2]));
					right+=(values->delta[0][3]*(vbs2[0]&vbs2[3]));
					right+=(values->delta[1][2]*(vbs2[1]&vbs2[2]));
					right+=(values->delta[1][3]*(vbs2[1]&vbs2[3]));
					right+=(values->delta[2][3]*(vbs2[2]&vbs2[3]));
					if ((vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3])==3) {
						right << 1;
						right/=3;
						}
					else	if ((vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3])==4) right >> 1;
					if (values->gap_must_cost) right += (values->gap_cost * (1+vbs2[4]));
          else right += values->gap_cost;
					}
				}
			/*
      else column_score (&diag, &down, &right, -1, j-1, values);
      */
			pmi[j].score=right + pmi[j-1].score - values->leading_gap_cost; 
			pmi[j].length = 1 + pmi[j-1].length;
			}
		
		say ("main loop")
		
		if ((values->extra_adjustment==-32000) && (values->coding==-32000)) {
			for (i = 1; i< min(a1->length,values->max_gap+1); i++) {
				pmi = m[i];
				pmi1= m[i-1];
				vbs1=values->b_string1[i-1];
				for (j = 1; j< min(a2->length,(i+values->max_gap+1)); j++) {
					vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM
          CHECKUM
					}
				}
			for (i = values->max_gap+1; i < a1->length; i++) {
				pmi = m[i];
				pmi1= m[i-1];
				vbs1=values->b_string1[i-1];
				for (j = i-values->max_gap; j < min(a2->length,(i+values->max_gap+1)); j++) {
					vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM
          CHECKUM
					}
				}
			}
		else if (values->extra_adjustment!=-32000) {
			for (i = 1; i< min(a1->length,values->max_gap+1); i++) {
				pmi = m[i];
				pmi1= m[i-1];
				vbs1=values->b_string1[i-1];
				for (j = 1; j< min(a2->length,(i+values->max_gap+1)); j++) {
					vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM
					if (pmi1[j].direction == DOWN) down -= values->extra_adjustment;
					if (pmi[j-1].direction == RIGHT) right -= values->extra_adjustment;
          CHECKUM
					}
				}
			for (i = values->max_gap+1; i < a1->length; i++) {
				pmi = m[i];
				pmi1= m[i-1];
				vbs1=values->b_string1[i-1];
				for (j = i-values->max_gap; j < min(a2->length,(i+values->max_gap+1)); j++) {
					vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
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
				vbs1=values->b_string1[i-1];
				for (j = 1; j< min(a2->length,(i+values->max_gap+1)); j++) {
					vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
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
				vbs1=values->b_string1[i-1];
				for (j = i-values->max_gap; j < min(a2->length,(i+values->max_gap+1)); j++) {
					vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
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
		vbs1=values->b_string1[a1->length-1];
		for (j = max(1,a1->length-values->max_gap); j < a2->length; j++) {
			vbs2=values->b_string2[j-1];
      SETME
			if ((!values->delta) && (!values->ttr)) {DOME}
			else if (values->ttr) {DOMENOW}
			else {DOMELATER}
			ADDEM
      right -= values->trailing_gap_cost;
			CHECKUM
     	}
		/*j = a2->length;*/
		vbs2=values->b_string2[a2->length-1];
		for (i = max(1,a2->length-values->max_gap); i < a1->length; i++) {
			vbs1=values->b_string1[i-1];
      SETME
			if ((!values->delta) && (!values->ttr)) {DOME}
			else if (values->ttr) {DOMENOW}
			else {DOMELATER}
			diag	+= m[i-1][a2->length-1].score;
			down	+= m[i-1][a2->length].score;
			right += m[i][a2->length-1].score;
				/*trailing gap fix - move to optimize */
				/*if (i == a1->length) right -= values->trailing_gap_cost;*/
				down -= values->trailing_gap_cost;
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
	vbs1=values->b_string1[a1->length-1];
	vbs2=values->b_string2[a2->length-1];
  SETME
	if ((!values->delta) && (!values->ttr)) {DOME}
	else if (values->ttr) {DOMENOW}
	else {DOMELATER}
	diag	+= m[a1->length-1][a2->length-1].score;
	down	+= m[a1->length-1][a2->length].score;
	right += m[a1->length][a2->length-1].score;
				
				right -= values->trailing_gap_cost;
				down -= values->trailing_gap_cost;
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
		a = allocate_alignment_and_strings (n, k+1);
		a->score = m[i][j].score;
		a->length = k;
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
		for (l = 0; l < n; l++) a->s[l][k] = '\0';
		net_gap=hold_i=hold_j=0;
		while (k > 0) {
			k--;
			switch (m[i][j].direction) {
			case DIAGONAL:
				for (l = 0; l < a1->n_seqs; l++) a->s[l + a1_offs][k] = a1->s[l][i-1];
				for (l = 0; l < a2->n_seqs; l++) a->s[l + a2_offs][k] = a2->s[l][j-1];
				i--;
				j--;
				break;
			case DOWN:
				for (l = 0; l < a1->n_seqs; l++) a->s[l + a1_offs][k] = a1->s[l][i-1];
				for (l = 0; l < a2->n_seqs; l++) a->s[l + a2_offs][k] = '-';
				i--;
				++hold_i;
				break;
			case RIGHT:
				for (l = 0; l < a1->n_seqs; l++) a->s[l + a1_offs][k] = '-';
				for (l = 0; l < a2->n_seqs; l++) a->s[l + a2_offs][k] = a2->s[l][j-1];
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
			a=dump_align(a);
			had_to_go_back=1;	
			goto redo_the_alignment;
			}
	else values->max_gap=net_gap+1;
		}
	
		say ("freeing")
		/*free_matrix (m, 1+a1->length);*/
		for (i=0;i<1+a1->length;i++) free(m[i]);
		free(m);
		if (values->phylo_score>0) {
			if ((!values->in_some) || (a->n_seqs==values->actual_num_sequences)) {
				if (!compare_aligns(a,values->previous)) a->score=values->previous->score;
				else {
					a->score = cladogram(a,values);
					if (values->previous) values->previous=dump_align(values->previous);
					values->previous=make_align(a);
					}
				}
			}

	say ("returning")
		return a;

}

