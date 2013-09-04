/*
Copyright 1992 Ward Wheeler all rights reserved
*/

#include "align3.h"

alignment *nw_low_mem (a1, a2,values)
alignment *a1, *a2;
parameters *values;
{
	cell **m;
	int i, j, k, l, n, down, right, diag, minimum, a1_offs, a2_offs;
	char holder,had_to_go_back;
	int net_gap,hold_gap,hold_i,hold_j;
	alignment *a;
	int size_holder;
	register cell *pmi, *pmi1;
	register int *vbs1, *vbs2;
	register int reg_int_holder;

	GOMAC
	say ("allocating")	
	values->max_gap=(1+max(a1->length-a2->length,a2->length-a1->length));
	
		/*redo if to gappy*/
		 had_to_go_back=0;
		 redo_the_low_alignment:;
		
		values->max_gap=min(values->max_gap,a2->length+1);
		size_holder=((2*values->max_gap)+3);

	if (values->show_mem) fprintf(stderr,"Allocating %d bytes\n",size_holder*(1+a1->length)*sizeof(cell));

	m = allocate_matrix (1+a1->length, size_holder);
	
	for (i=values->max_gap+1;((i<=a1->length) && (i-values->max_gap-1)<=a2->length);i++) {
		m[i][0].score=BARRIER_COST;
		m[i][0].direction=3;
		}
	for (i=0;((i<=a1->length) && ((i+values->max_gap+1)<=a2->length));i++) {
		m[i][values->max_gap+values->max_gap+2].score=BARRIER_COST;	
		m[i][values->max_gap+values->max_gap+2].direction=3;	
		}
			

	say ("filling in")
		 m[0][values->max_gap+1].direction = DIAGONAL;
		 m[0][values->max_gap+1].score = 0;
		 m[0][values->max_gap+1].length = 0;
		 for (i = 1; ((i <= a1->length) && (i<=values->max_gap)); i++) {
      vbs1=values->b_string1[i-1];
			m[i][1-i+values->max_gap].direction = DOWN;
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
          else down+=values->gap_cost;
					}
				}
			else {
				if (!(vbs1[2]+vbs1[1]+vbs1[0]+vbs1[3]-1)) {
        	if (values->gap_must_cost) down= (values->gap_cost * (vbs1[4]+1));
          else down=values->gap_cost;
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
          else down +=values->gap_cost;
					}
				}
			/*
      else column_score (&diag, &down, &right,i-1, -1, values);
      */
      m[i][1-i+values->max_gap].score = down + m[i-1][2-i+values->max_gap].score - values->leading_gap_cost;
			m[i][1-i+values->max_gap].length = 1 + m[i-1][2-i+values->max_gap].length;
	 		}
			pmi=m[0];
		 for (j = 1; ((j<=a2->length) && (j<=values->max_gap)); j++) {
      vbs2=values->b_string2[j-1];
			pmi[j+values->max_gap+1].direction = RIGHT;
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
          else right=values->gap_cost;
					}
				}
			else {
				if (!(vbs2[2]+vbs2[1]+vbs2[0]+vbs2[3]-1)) {
        	if (values->gap_must_cost) right= (values->gap_cost * (1+vbs2[4]));
          else right=values->gap_cost;
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
      pmi[j+values->max_gap+1].score = right + m[0][j+values->max_gap].score - values->leading_gap_cost;
			pmi[j+values->max_gap+1].length = 1 + pmi[j+values->max_gap].length;
			}
		
		say ("main loop")
		
		if ((values->extra_adjustment==-32000) && (values->coding==-32000)) {
		say("First half of main middle")
			for (i = 1; (i < a1->length) && (i<=values->max_gap); i++) {
				pmi=m[i];
	      vbs1=values->b_string1[i-1];
				pmi1=m[i-1];
				for (j = 1; (j < a2->length) && (j<=(i+values->max_gap)); j++) {
		      vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM_LOWMEM
          CHECKUM_LOWMEM
					}
				}
	say("Second half of main middle")
			for (i = values->max_gap+1; i < a1->length; i++) {
				pmi=m[i];
				pmi1=m[i-1];
	      vbs1=values->b_string1[i-1];
				for (j = i-values->max_gap; (j < a2->length) && (j<=(i+values->max_gap)); j++) {
	      	vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM_LOWMEM
          CHECKUM_LOWMEM
					}
				}
			}
		else if (values->extra_adjustment!=-32000) {
			for (i = 1; (i < a1->length) && (i<=values->max_gap); i++) {
			pmi=m[i];
			pmi1=m[i-1];
      vbs1=values->b_string1[i-1];
				for (j = 1; (j < a2->length) && (j<=(i+values->max_gap)); j++) {
     		 vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM_LOWMEM
					if (pmi1[reg_int_holder+1].direction == DOWN) down -= values->extra_adjustment;
					if (pmi[reg_int_holder-1].direction == RIGHT) right -= values->extra_adjustment;
          CHECKUM_LOWMEM
					}
				}
			for (i = values->max_gap+1; i < a1->length; i++) {
				pmi=m[i];
				pmi1=m[i-1];
      	vbs1=values->b_string1[i-1];
				for (j = i-values->max_gap; (j < a2->length) && (j<=(i+values->max_gap)); j++) {
		      vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM_LOWMEM
					if (pmi1[reg_int_holder+1].direction == DOWN) down -= values->extra_adjustment;
					if (pmi[reg_int_holder-1].direction == RIGHT) right -= values->extra_adjustment;
          CHECKUM_LOWMEM
					}
				}
			}
		else if (values->coding!=-32000) {
			for (i = 1; (i < a1->length) && (i<=values->max_gap); i++) {
				pmi=m[i];
				pmi1=m[i-1];
	      vbs1=values->b_string1[i-1];
				for (j = 1; (j < a2->length) && (j<=(i+values->max_gap)); j++) {
	      	vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
					ADDEM_LOWMEM
          /*fixes for gap lengths*/
					holder=1;
					for (k=i-1;k>=0;k--) {
						if (m[k][j-k+values->max_gap+1].direction==DOWN) ++holder;
						else break;
						}
					if ((holder%3)==0) down -= values->coding;
					holder=1;
					for (k=i-1;k>=0;k--) {
						if (pmi[k-i+values->max_gap+1].direction==RIGHT) ++holder;
						else break;
						}
					if ((holder%3)==0) right -= values->coding;
					CHECKUM_LOWMEM
          }
				}
			for (i = values->max_gap+1; i < a1->length; i++) {
				pmi=m[i];
				pmi1=m[i-1];
      	vbs1=values->b_string1[i-1];
				for (j = i-values->max_gap; (j < a2->length) && (j<=(i+values->max_gap)); j++) {
      		vbs2=values->b_string2[j-1];
          SETME
					if ((!values->delta) && (!values->ttr)) {DOME}
					else if (values->ttr) {DOMENOW}
					else {DOMELATER}
          ADDEM_LOWMEM
					/*fixes for gap lengths*/
					holder=1;
					for (k=i-1;k>=0;k--) {
						if (m[k][j-k+values->max_gap+1].direction==DOWN) ++holder;
						else break;
						}
					if ((holder%3)==0) down -= values->coding;
					holder=1;
					for (k=i-1;k>=0;k--) {
						if (pmi[k-i+values->max_gap+1].direction==RIGHT) ++holder;
						else break;
						}
					if ((holder%3)==0) right -= values->coding;
					CHECKUM_LOWMEM		
					}
				}
			}

	say("Trailing stuff")
	say("	i")	
	/* trailing gap fix - moved to optimize i = a1->length*/
	pmi=m[a1->length];
	pmi1=m[a1->length-1];
	vbs1=values->b_string1[a1->length-1];
	for (j =max(1,a1->length-values->max_gap); j < a2->length; j++) {
		vbs2=values->b_string2[j-1];
    SETME
		if ((!values->delta) && (!values->ttr)) {DOME}
		else if (values->ttr) {DOMENOW}
		else {DOMELATER}
    ADDEM_LOWMEM
		right -= values->trailing_gap_cost;
    CHECKUM_LOWMEM
		}
	say("	j")
	/*j = a2->length;*/
  vbs2=values->b_string2[a2->length-1];
	for (i = max(1,a2->length-values->max_gap); i < a1->length; i++) {
  	vbs1=values->b_string1[i-1];
    SETME
		if ((!values->delta) && (!values->ttr)) {DOME}
		else if (values->ttr) {DOMENOW}
		else {DOMELATER}
		diag += m[i-1][a2->length-i+1+values->max_gap].score;
		down += m[i-1][a2->length-i+2+values->max_gap].score;
		right += m[i][a2->length-i+values->max_gap].score;
		down -= values->trailing_gap_cost;
		minimum = min (min (diag, down), right);
		m[i][a2->length-i+values->max_gap+1].score = minimum;
		if (diag == minimum) {
			m[i][a2->length-i+values->max_gap+1].direction = DIAGONAL;
			m[i][a2->length-i+values->max_gap+1].length = 1 + m[i-1][a2->length-i+1+values->max_gap].length;
			}
		else if (down == minimum) {
			m[i][a2->length-i+values->max_gap+1].direction = DOWN;
			m[i][a2->length-i+values->max_gap+1].length = 1 + m[i-1][a2->length-i+2+values->max_gap].length;
			}
		else {
			m[i][a2->length-i+values->max_gap+1].direction = RIGHT;
			m[i][a2->length-i+values->max_gap+1].length = 1 + m[i][a2->length-i+values->max_gap].length;
			}	
		}
		say("	Both")	
/*
		j = a2->length;
		i = a1->length; 
*/
  vbs1=values->b_string1[a1->length-1];
  vbs2=values->b_string2[a2->length-1];
  SETME
	if ((!values->delta) && (!values->ttr)) {DOME}
	else if (values->ttr) {DOMENOW}
	else {DOMELATER}
	diag += m[a1->length-1][a2->length-a1->length+1+values->max_gap].score;
	down += m[a1->length-1][a2->length-a1->length+2+values->max_gap].score;
	right += m[a1->length][a2->length-a1->length+values->max_gap].score;
	right -= values->trailing_gap_cost;
	down -= values->trailing_gap_cost;
				
	minimum = min (min (diag, down), right);
	m[a1->length][a2->length-a1->length+values->max_gap+1].score = minimum;
	if (diag == minimum) {
		m[a1->length][a2->length-a1->length+values->max_gap+1].direction = DIAGONAL;
		m[a1->length][a2->length-a1->length+values->max_gap+1].length = 1 + m[a1->length-1][a2->length-a1->length+1+values->max_gap].length;
		}
	else if (down == minimum) {
		m[a1->length][a2->length-a1->length+values->max_gap+1].direction = DOWN;
		m[a1->length][a2->length-a1->length+values->max_gap+1].length = 1 + m[a1->length-1][a2->length-a1->length+2+values->max_gap].length;
		}
	else {
		m[a1->length][a2->length-a1->length+values->max_gap+1].direction = RIGHT;
		m[a1->length][a2->length-a1->length+values->max_gap+1].length = 1 + m[a1->length][a2->length-a1->length+values->max_gap].length;
		}
	
		say ("trace back")
		i = a1->length;
		j = a2->length;
		k = m[i][j-i+values->max_gap+1].length;
		n = a1->n_seqs + a2->n_seqs;
		a = allocate_alignment_and_strings (n, k+1);
		a->score = m[i][j-i+values->max_gap+1].score;
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
			switch (m[i][j-i+values->max_gap+1].direction) {
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
				fprintf (stderr, "Bad direction in matrix: %d\n", m[i][j-i+values->max_gap+1].direction);
			exit (-1);
			}
		}
	if (values->iter) {	
		net_gap=max(hold_i,hold_j);
		if (net_gap >= values->max_gap) {
			values->max_gap=(2*net_gap)+2;
			a=dump_align(a);
			had_to_go_back=1;	
			say ("freeing")
			free_matrix (m, 1+a1->length);
			free(m);
			if (values->max_gap < a2->length ) goto redo_the_low_alignment;
			else goto end_thang_for_bag; /*do regular alignment with "iter" off*/
			}
	else values->max_gap=net_gap+1;
		}
	
		say ("freeing")
		free_matrix (m, 1+a1->length);
		free(m);

		if (values->phylo_score>0) {
		if ((!values->in_some) || (a->n_seqs==values->actual_num_sequences)) {
			if (!compare_aligns(a,values->previous)) a->score=values->previous->score;
			else 	{
				a->score = cladogram(a,values);
				if (values->previous) values->previous=dump_align(values->previous);
				values->previous=make_align(a);
				}
			}
		}
 	say ("returning")
	return a;
	/*if lowmem inappropriate*/
	end_thang_for_bag:;
	values->iter=0;
	a=nw_real(a1,a2,values);
	values->iter=1;
	return a;
}
