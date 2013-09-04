/*
Copyright 1992 Ward Wheeler all rights reserved
change back gaps at least as an option
but need to remember the leading and trailing
could change cost or something (to 1/2) perhaps
*/
#include "align3.h"

int first;
extern int charbit[256];
extern char bitchar[32];

alignment *nw (input_a1,input_a2,v_in)
alignment *input_a1, *input_a2;
parameters *v_in;
{
	int implicit_gap_size,smaller,gap_not_too_big;
	alignment *b,**a, **a1, **a2;
	int check_diffs;
	char iter_holder,max_gap_holder,fix_condition,low_mem_holder;
	char *b_string1,*b_string2;
	int i, number_of_data_sets, *type, j;
	int current_start, total_length, ii, holder;
	int leading_holder, trailing_holder; /*,temp1,temp2;*/
	int a1_mis, a2_mis;
	parameters *values,*valuesI;


	values=v_in;
	if (((!values->in_some) || (((input_a1->n_seqs+input_a2->n_seqs)==values->actual_num_sequences))) && (values->number_of_input_alignments>1) &&
	    (values->phylo_score>0) && (!values->new_optimization)) {
		holder=0;
		values->in_some=1;
	}
	else holder=1;

	if ((values->number_of_input_alignments>1) && (values->chop)) {
		leading_holder=values->leading_gap_cost;
		trailing_holder=values->trailing_gap_cost;
		values->trailing_gap_cost=0;
		values->leading_gap_cost=0;
	}

	/*allocate stuff for the number*/
	if (values->number_of_input_alignments==0) {
		input_a1->type_weight=input_a2->type_weight=0;
		number_of_data_sets=1;
	}
	/*else if ((!values->new_optimization) && (values->chop==0)) {
		number_of_data_sets=1;
		fprintf(stderr,"Here ");
		}*/
	else number_of_data_sets=values->number_of_input_alignments;
	type=(int *)malloc(number_of_data_sets*sizeof(int));
	assert((int) type);
	a=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a);
	a1=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a1);
	a2=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a2);
	/*temp1=temp2=0;*/
	for (i=0;i<number_of_data_sets;i++) {
		a[i]=NULL;
		a1[i]=NULL;
		a2[i]=NULL;
		if (!values->new_optimization || (values->new_optimization && values->in_optimize) || (number_of_data_sets==1)) {
		  a1[i]=make_align(input_a1);
		  a2[i]=make_align(input_a2);
		}
		if (!values->new_optimization || (values->new_optimization && values->in_optimize)) {
		  if (values->number_of_input_alignments>0) {
		    redo_for_type(a1[i],i,values);
		    redo_for_type(a2[i],i,values);
		  }
		}

	}
	if ((values->new_optimization && (!values->in_optimize)) && (values->number_of_input_alignments>0)) {
		redo_for_type_2(a1,input_a1,values,number_of_data_sets);
	       	redo_for_type_2(a2,input_a2,values,number_of_data_sets);
		}

/*set flag to align or no type[]==1 align 0==no*/
	if (values->number_of_input_alignments==0) type[0]=1;
	/*else if ((!values->new_optimization) && (values->chop==0)) type[0]=1;*/
	else {
		for (i=0;i<number_of_data_sets;i++) {
			if ((values->data_sets[i]+1)==1) type[i]=1;
			else if ((values->data_sets[i]+1)==2) type[i]=0;
			else if ((values->data_sets[i]+1)==3) type[i]=0;
			else if ((values->data_sets[i]+1)==4) type[i]=1;
		}
	}
	first=0;
	for (ii=0;ii<number_of_data_sets;ii++) {
		if (!v_in->other_parm) valuesI=v_in;
		else if (!v_in->other_parm[ii]) valuesI=v_in;
		else valuesI=v_in->other_parm[ii];
		valuesI->in_optimize=values->in_optimize;
		valuesI->in_some=values->in_some;
		valuesI->new_optimization=values->new_optimization;
		/*fix for leading and trialing in chop*/
		if ((values->number_of_input_alignments>1) && (values->chop)) {
			if (ii==0) valuesI->leading_gap_cost=leading_holder;
			else valuesI->leading_gap_cost=0;
			if (ii==(number_of_data_sets-1)) valuesI->trailing_gap_cost=trailing_holder;
			else valuesI->trailing_gap_cost=0;
		}
		if (type[ii]) {
			/*check to see if missing*/
			a1_mis=only_missing(a1[ii]);
			a2_mis=only_missing(a2[ii]);
			/*fprintf(stderr,"A1 (%s) %d A2 (%s) %d\n",input_a1->name,a1_mis,input_a2->name,a2_mis);*/
			}
		if (!type[ii]) a[ii]=optimization_alignment(a1[ii],a2[ii],valuesI,values->data_sets[ii]);
		else if ((a1_mis || a2_mis)) {
			if (a1_mis && (!a2_mis)) a[ii]=missing_alignment(a1[ii],a2[ii],valuesI,0);
			else if ((!a1_mis) && a2_mis) a[ii]=missing_alignment(a1[ii],a2[ii],valuesI,1);
			else a[ii]=missing_alignment(a1[ii],a2[ii],valuesI,2);
			a[ii]->score=0;
			}
		else {
			gap_not_too_big=1;
			fix_condition=0;
			if (((5*max (a1[ii]->length-a2[ii]->length,a2[ii]->length-a1[ii]->length))/2) >= (min(a1[ii]->length,a2[ii]->length))) {
				iter_holder=values->iter;
				max_gap_holder=valuesI->max_gap;
				low_mem_holder=valuesI->low_mem;
				valuesI->iter=0;
				valuesI->low_mem=0;
				valuesI->max_gap=MAX_SEQUENCE_SIZE;
				fix_condition=1;
			}
			if (valuesI->low_mem) {
				implicit_gap_size=max(a1[ii]->length - a2[ii]->length,a2[ii]->length - a1[ii]->length);
				implicit_gap_size*=2;
				implicit_gap_size=max(valuesI->max_gap,implicit_gap_size);
				if (a1[ii]->length<a2[ii]->length) smaller=a1[ii]->length;
				else smaller=a2[ii]->length;
				if ((smaller*((2*implicit_gap_size)+3)) > ((a2[ii]->length+1)*(a1[ii]->length+1))) gap_not_too_big=0;
				if ((a1[ii]->length<100) || (a2[ii]->length<100)) gap_not_too_big=0;
				if (gap_not_too_big) {
					if (a1[ii]->length > a2[ii]->length) {
						if (!values->new_optimization) {
					    	valuesI->b_string1=do_bit_thing(a2[ii],valuesI);
    						valuesI->b_string2=do_bit_thing(a1[ii],valuesI);
							a[ii]=nw_low_mem(a2[ii], a1[ii],valuesI);
							for (i=0;i<a2[ii]->length;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a1[ii]->length;i++) free(valuesI->b_string2[i]);
						}
						else {
						    valuesI->b_string1=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string1);
                            valuesI->b_string1[0]=(int *)malloc(a2[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string1[0]);
                            for (i=0;i<a2[ii]->length;i++)  valuesI->b_string1[0][i]=charbit[a2[ii]->s[0][i]];
						    valuesI->b_string2=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string2);
                            valuesI->b_string2[0]=(int *)malloc(a1[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string2[0]);
                            for (i=0;i<a1[ii]->length;i++)  valuesI->b_string2[0][i]=charbit[a1[ii]->s[0][i]];
							a[ii]=nw_low_mem_new_opt(a2[ii], a1[ii],valuesI);
							for (i=0;i<a2[ii]->n_seqs;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a1[ii]->n_seqs;i++) free(valuesI->b_string2[i]);
						}

					}
					else {
						if (!values->new_optimization) {
	    					valuesI->b_string1=do_bit_thing(a1[ii],valuesI);
    						valuesI->b_string2=do_bit_thing(a2[ii],valuesI);
							a[ii]=nw_low_mem(a1[ii], a2[ii],valuesI);
							for (i=0;i<a1[ii]->length;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a2[ii]->length;i++) free(valuesI->b_string2[i]);
						}
						else {
						    valuesI->b_string1=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string1);
                            valuesI->b_string1[0]=(int *)malloc(a1[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string1[0]);
							valuesI->b_string2=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string2);
                            valuesI->b_string2[0]=(int *)malloc(a2[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string2[0]);
                            for (i=0;i<a1[ii]->length;i++)  valuesI->b_string1[0][i]=charbit[a1[ii]->s[0][i]];
                            for (i=0;i<a2[ii]->length;i++)  valuesI->b_string2[0][i]=charbit[a2[ii]->s[0][i]];
							a[ii]=nw_low_mem_new_opt(a1[ii], a2[ii],valuesI);
							for (i=0;i<a1[ii]->n_seqs;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a2[ii]->n_seqs;i++) free(valuesI->b_string2[i]);
						}
					}
				}
				else {
					if (a1[ii]->length > a2[ii]->length) {
						if (!values->new_optimization) {
			    			valuesI->b_string1=do_bit_thing(a2[ii],valuesI);
		    				valuesI->b_string2=do_bit_thing(a1[ii],valuesI);
							a[ii]=nw_real(a2[ii], a1[ii],valuesI);
							for (i=0;i<a2[ii]->length;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a1[ii]->length;i++) free(valuesI->b_string2[i]);
						}
						else {
							if (!ii) first=1;
						    valuesI->b_string1=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string1);
                            valuesI->b_string1[0]=(int *)malloc(a2[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string1[0]);
                            for (i=0;i<a2[ii]->length;i++)  valuesI->b_string1[0][i]=charbit[a2[ii]->s[0][i]];
						    valuesI->b_string2=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string2);
                            valuesI->b_string2[0]=(int *)malloc(a1[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string2[0]);
                            for (i=0;i<a1[ii]->length;i++)  valuesI->b_string2[0][i]=charbit[a1[ii]->s[0][i]];
							a[ii]=nw_real_new_opt(a2[ii], a1[ii],valuesI);
							for (i=0;i<a2[ii]->n_seqs;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a1[ii]->n_seqs;i++) free(valuesI->b_string2[i]);
						}
					}
					else {
						if (!values->new_optimization) {
					    	valuesI->b_string1=do_bit_thing(a1[ii],valuesI);
				    		valuesI->b_string2=do_bit_thing(a2[ii],valuesI);
							a[ii]=nw_real(a1[ii], a2[ii],valuesI);
							for (i=0;i<a1[ii]->length;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a2[ii]->length;i++) free(valuesI->b_string2[i]);
						}
						else {
						    valuesI->b_string1=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string1);
                            valuesI->b_string1[0]=(int *)malloc(a1[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string1[0]);
							valuesI->b_string2=(int **)malloc(sizeof(int *));
                            assert((int)valuesI->b_string2);
                            valuesI->b_string2[0]=(int *)malloc(a2[ii]->length*sizeof(int));
                            assert((int)valuesI->b_string2[0]);
                            for (i=0;i<a1[ii]->length;i++)  valuesI->b_string1[0][i]=charbit[a1[ii]->s[0][i]];
                            for (i=0;i<a2[ii]->length;i++)  valuesI->b_string2[0][i]=charbit[a2[ii]->s[0][i]];
							a[ii]=nw_real_new_opt(a1[ii], a2[ii],valuesI);
							for (i=0;i<a1[ii]->n_seqs;i++) free(valuesI->b_string1[i]);
							for (i=0;i<a2[ii]->n_seqs;i++) free(valuesI->b_string2[i]);
						}
					}
				}
			}
			else {
				if (a1[ii]->length > a2[ii]->length) {
					if (!values->new_optimization) {
	    				valuesI->b_string1=do_bit_thing(a2[ii],valuesI);
    					valuesI->b_string2=do_bit_thing(a1[ii],valuesI);
						a[ii]=nw_real(a2[ii], a1[ii],valuesI);
						for (i=0;i<a2[ii]->length;i++) free(valuesI->b_string1[i]);
						for (i=0;i<a1[ii]->length;i++) free(valuesI->b_string2[i]);
					}
					else {
						if (!ii) first=1;
						valuesI->b_string1=(int **)malloc(sizeof(int *));
                        assert((int)valuesI->b_string1);
                        valuesI->b_string1[0]=(int *)malloc(a2[ii]->length*sizeof(int));
                        assert((int)valuesI->b_string1[0]);
                        for (i=0;i<a2[ii]->length;i++)  valuesI->b_string1[0][i]=charbit[a2[ii]->s[0][i]];
						valuesI->b_string2=(int **)malloc(sizeof(int *));
                        assert((int)valuesI->b_string2);
                        valuesI->b_string2[0]=(int *)malloc(a1[ii]->length*sizeof(int));
                        assert((int)valuesI->b_string2[0]);
                        for (i=0;i<a1[ii]->length;i++)  valuesI->b_string2[0][i]=charbit[a1[ii]->s[0][i]];
						a[ii]=nw_real_new_opt(a2[ii], a1[ii],valuesI);
						for (i=0;i<a2[ii]->n_seqs;i++) free(valuesI->b_string1[i]);
						for (i=0;i<a1[ii]->n_seqs;i++) free(valuesI->b_string2[i]);
					}
				}
				else {
					if (!values->new_optimization) {
		    			valuesI->b_string1=do_bit_thing(a1[ii],valuesI);
		    			valuesI->b_string2=do_bit_thing(a2[ii],valuesI);
						a[ii]=nw_real(a1[ii], a2[ii],valuesI);
						for (i=0;i<a1[ii]->length;i++) free(valuesI->b_string1[i]);
						for (i=0;i<a2[ii]->length;i++) free(valuesI->b_string2[i]);
					}
					else {
						valuesI->b_string1=(int **)malloc(sizeof(int *));
                        assert((int)valuesI->b_string1);
                        valuesI->b_string1[0]=(int *)malloc(a1[ii]->length*sizeof(int));
                        assert((int)valuesI->b_string1[0]);
						valuesI->b_string2=(int **)malloc(sizeof(int *));
                        assert((int)valuesI->b_string2);
                        valuesI->b_string2[0]=(int *)malloc(a2[ii]->length*sizeof(int));
                        assert((int)valuesI->b_string2[0]);
                        for (i=0;i<a1[ii]->length;i++)  valuesI->b_string1[0][i]=charbit[a1[ii]->s[0][i]];
                        for (i=0;i<a2[ii]->length;i++)  valuesI->b_string2[0][i]=charbit[a2[ii]->s[0][i]];
						a[ii]=nw_real_new_opt(a1[ii], a2[ii],valuesI);
						for (i=0;i<a1[ii]->n_seqs;i++) free(valuesI->b_string1[i]);
						for (i=0;i<a2[ii]->n_seqs;i++) free(valuesI->b_string2[i]);
					}
				}
			}
			if (fix_condition) {
				valuesI->iter=iter_holder;
				valuesI->low_mem=low_mem_holder;
				valuesI->max_gap=max_gap_holder;
			}

			free(valuesI->b_string1);
			free(valuesI->b_string2);
			/*set to zero because NW only deals with one data set at a time*/
			a[ii]->type_weight=0;
			if (values->new_optimization) if (!values->in_optimize) a[ii] = make_ambig(a[ii],valuesI);
			if (a1_mis || a2_mis) a[ii]->score=0;
		}
		a[ii]->type_weight=0;
		if (type[ii]) {
		    for (i=0;i<a[ii]->n_seqs;i++) if (only_missing_string(a[ii]->s[i])) for (j=0;j<a[ii]->length;j++) a[ii]->s[i][j]='X';
		    if (values->expand_X) make_contig_X(a[ii],valuesI);
        	}
		if ((values->phylo_score == 0) && (!values->new_optimization) && type[ii]) a[ii]->score=no_align_column_score(a[ii],valuesI);
		}
	/*fprintf(stderr,"\n");*/
	/*put them back together*/
	if (values->number_of_input_alignments==0) b=make_align(a[0]);
	else if (!values->new_optimization) b=glue_back_after_chop(a,values);
	else if (values->new_optimization && values->in_optimize)  b=glue_back_after_chop(a,values);
	else {
		/*stick back together and make b*/
		/*iteratively reallocate and copy a[0]*/
		b=(alignment *)malloc(sizeof(alignment));
		assert((int) b);
		if (strcmp (input_a1->name, input_a2->name) > 0) b->name = pair_names (input_a2, input_a1);
		else b->name = pair_names (input_a1, input_a2);
		b->taxon_name=(char **)malloc(sizeof(char *));
		assert((int) b->taxon_name);
		b->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
		assert ((int) b->taxon_name[0]);
		b->taxon_name[0][0]='H';
		b->taxon_name[0][1]='y';
		b->taxon_name[0][2]='p';
		b->taxon_name[0][3]='A';
		b->taxon_name[0][4]='n';
		b->taxon_name[0][5]='c';
		b->taxon_name[0][6]='\0';
		total_length=0;
		for (i=0;i<number_of_data_sets;i++) total_length+=a[i]->length;
		b->length=total_length;
		b->score=0;
		for (i=0;i<number_of_data_sets;i++) {
			b->score+=(values->data_set_weights[i]*a[i]->score);
		}
		b->n_seqs=1;
		b->type_weight=1;
		b->s=(char **)malloc(2*sizeof(char));
		assert((int) b->s);
		for (i=0;i<2;i++) {
			b->s[i]=(char *)malloc((1+b->length)*sizeof(char));
			assert((int) b->s[i]);
			b->s[i][b->length]='\0';
		}
		current_start=0;
		for (i=0;i<number_of_data_sets;i++) {
			for (j=0;j<a[i]->length;j++) {
				b->s[0][j+current_start]=a[i]->s[0][j];
				b->s[1][j+current_start]=values->data_sets[i]+1;
				if (j==(a[i]->length-1))        b->s[1][j+current_start] += (64);
			}
			current_start+=a[i]->length;
		}
	}
	if (values->new_optimization && b->type_weight && values->in_optimize) {
		b->score=0;
		for (i=0;i<number_of_data_sets;i++) b->score+=(values->data_set_weights[i]*a[i]->score);
		}
	for (i=0;i<number_of_data_sets;i++) {
		if (a[i]) a[i]=dump_align(a[i]);
		if (a1[i]) a1[i]=dump_align(a1[i]);
		if (a2[i]) a2[i]=dump_align(a2[i]);
	}
	free(a);
	free(a1);
	free(a2);
	free(type);
	/*fprintf(stderr,"%s %s\n",b->name,b->s[0]);
    if (!values->new_optimization) {
        printf("Just after NW\n");
        for (ii=0;ii< b->n_seqs;ii++) {
            for (j=0;j<b->length;j++) printf("%4d",b->s[ii][j]);
            printf("\n");
    }}*/
	if ((values->number_of_input_alignments>1) && (values->chop)) {
		values->trailing_gap_cost=trailing_holder;
		values->leading_gap_cost=leading_holder;
	}
	if (!holder) {
		values->in_some=0;
		b->score=cladogram(b,values);
	}
    if (b->type_weight!=1) b->type_weight=0;
	return b;
}

void make_contig_X(a,values)
alignment *a;
parameters *values;
{
/*function makes '-' which hit 'X' before other bases into 'X' */
int i,hit_X,j,k;

for (k=0;k<a->n_seqs;k++) {
    for (i=0;i<a->length;i++) if (a->s[k][i]=='-') {
    	hit_X=0;
    	for (j=(i-1);j>=0;j--) {
    		if (a->s[k][j]=='X') { hit_X=1;break; }
    		else if (a->s[k][j]!='-') break;
    		}
    	if (hit_X) a->s[k][i]='X';
    	else {
    		for (j=(i+1);j<a->length;j++) {
    			if (a->s[k][j]=='X') { hit_X=1;break; }
    			else if (a->s[k][j]!='-') break;
    			}
    		if (hit_X) a->s[k][i]='X';
    		}
    	if (hit_X && values->new_optimization) { /*but oposite must != - or more*/
    	    if (!(is_there_a_gap(a->s[!k][i]))) {
    	        if (is_internal(a->s[k],i,a->length)) a->score-=values->gap_cost;
    	        else a->score -=(values->gap_cost-values->leading_gap_cost);
    	        }
    	    }
        }
    }
}

int is_there_a_gap(a)
char a;
{
int in_a;
switch (a) {
	case 'A': return 0;
	case 'C': return 0;
	case 'G': return 0;
	case 'T': return 0;
	case 'N': return 0;
	case 'U': return 0;
	case 'R': return 0;
	case 'Y': return 0;
	case 'M': return 0;
	case 'W': return 0;
	case 'S': return 0;
	case 'K': return 0;
	case 'B': return 0;
	case 'D': return 0;
	case 'H': return 0;
	case 'V': return 0;
	default:  return 1;
	}
}

alignment *optimization_alignment(a1,a2,values,type)
alignment *a1, *a2;
parameters *values;
int type; /* if 1 => bases; if 2=> characters*/
{
	alignment *a;
	int i,j,leave_two;



	if (values->in_optimize || (!values->new_optimization)) leave_two=1;
	else leave_two=0;

	a=NULL;
	if (a1->length!=a2->length) {
		fprintf(stderr,"Cannot optimize alignments of different length\n");
		print_inter_dig(a1,values);
		print_inter_dig(a2,values);
		fprintf(stderr,"%s length %d %s length %d\n",a1->name,a1->length,a2->name,a2->length);
		exit(-1);
	}
	/*allocate*/
	a=(alignment *)malloc(sizeof(alignment));
	assert((int) a);
	a->length=a1->length;
	a->type_weight=0;
	if (strcmp (a1->name, a2->name) > 0) a->name = pair_names (a2, a1);
	else a->name = pair_names (a1, a2);

	if (values->new_optimization) {
		a->n_seqs=1+leave_two;
		a->s=(char **)malloc((1+leave_two)*sizeof(char *));
		assert((int) a->s);
		a->s[0]=(char *)malloc((a->length+1)*sizeof(char));
		assert ((int) a->s[0]);
		a->s[0][a->length]='\0';
		if (leave_two) {
			a->s[1]=(char *)malloc((a->length+1)*sizeof(char));
			assert ((int) a->s[1]);
			a->s[1][a->length]='\0';
			}
		a->taxon_name=(char **)malloc((1+leave_two)*sizeof(char *));
		assert((int) a->taxon_name);
		if (!leave_two) {
			a->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
			assert ((int) a->taxon_name[0]);
			a->taxon_name[0][0]='H';
			a->taxon_name[0][1]='y';
			a->taxon_name[0][2]='p';
			a->taxon_name[0][3]='A';
			a->taxon_name[0][4]='n';
			a->taxon_name[0][5]='c';
			a->taxon_name[0][6]='\0';
			}
		else {
			a->taxon_name[0]=(char *)malloc((1+strlen(a1->taxon_name[0]))*sizeof(char));
			assert ((int) a->taxon_name[0]);
			a->taxon_name[1]=(char *)malloc((1+strlen(a2->taxon_name[0]))*sizeof(char));
			assert ((int) a->taxon_name[1]);
			a->taxon_name[0]=(char *)strcpy(a->taxon_name[0],a1->taxon_name[0]);
			a->taxon_name[1]=(char *)strcpy(a->taxon_name[1],a2->taxon_name[0]);
			}
		a->score=0;

		/*optimize*/
		if (leave_two) {
			if (!values->jackboot) {
			    if (type==1) {
			        for (i=0;i<a->length;i++) {
				        a->s[0][i]=a1->s[0][i];
					    a->s[1][i]=a2->s[0][i];
					    a->score+=values->lookup[a1->s[0][i]][a2->s[0][i]].cost;
					    }
					}
				else {
    				for (i=0;i<a->length;i++) {
    					a->s[0][i]=a1->s[0][i];
    					a->s[1][i]=a2->s[0][i];
                        if (!(a1->s[0][i]&a2->s[0][i])) ++a->score;
                        }
					}
				}
			else {
                if (type==1) {
			        for (i=0;i<a->length;i++) {
				        a->s[0][i]=a1->s[0][i];
					    a->s[1][i]=a2->s[0][i];
					    if (values->jack_array[(i*1000)/a->length]) a->score+=values->lookup[a1->s[0][i]][a2->s[1][i]].cost;
					    }
					}
				else {
    				for (i=0;i<a->length;i++) {
    					a->s[0][i]=a1->s[0][i];
    					a->s[1][i]=a2->s[0][i];
                        if (!(a1->s[0][i]&a2->s[0][i])) if (values->jack_array[(i*1000)/a->length]) ++a->score;
                        }
					}

				}
			}
		else {
    		if (!values->jackboot) {
			    if (type==1) {
			        for (i=0;i<a->length;i++) {
				        a->s[0][i]=values->lookup[a1->s[0][i]][a2->s[0][i]].base;
					    a->score+=values->lookup[a1->s[0][i]][a2->s[0][i]].cost;
					    }
					}
				else {
    				for (i=0;i<a->length;i++) {
					a->s[0][i]=(a1->s[0][i]&a2->s[0][i]);
					if (!a->s[0][i]) {
						a->s[0][i]=(a1->s[0][i]|a2->s[0][i]);
						++a->score;
						}
                        }
					}
				}
			else {
                if (type==1) {
			        for (i=0;i<a->length;i++) {
				        a->s[0][i]=values->lookup[a1->s[0][i]][a2->s[0][i]].base;
					    a->score+=values->lookup[a1->s[0][i]][a2->s[0][i]].cost;
					    }
					}
				else {
    				for (i=0;i<a->length;i++) {
       					a->s[0][i]=(a1->s[0][i]&a2->s[0][i]);
	    				if (!a->s[0][i]) {
		    				a->s[0][i]=(a1->s[0][i]|a2->s[0][i]);
    			    		if (values->jack_array[(i*1000)/a->length]) ++a->score;
				    		}
                        }
					}

				}

			}
		}
	else /*if (values->chop>0)*/ {
	    a->n_seqs=a1->n_seqs+a2->n_seqs;
		a->s=(char **)malloc(a->n_seqs*sizeof(char *));
		assert((int) a->s);
		a->taxon_name=(char **)malloc(a->n_seqs*sizeof(char *));
		assert((int) a->taxon_name);
		for (j=0;j<a->n_seqs;j++) {
			a->s[j]=(char *)malloc((a->length+1)*sizeof(char));
			assert ((int) a->s[j]);
			a->s[j][a->length]='\0';
		}

		for (j=0;j<a1->n_seqs;j++) {
				a->taxon_name[j]=(char *)malloc((1+strlen(a1->taxon_name[j]))*sizeof(char));
				assert ((int) a->taxon_name[j]);
				a->taxon_name[j]=strcpy(a->taxon_name[j],a1->taxon_name[j]);
				a->s[j]=strcpy(a->s[j],a1->s[j]);
		}
		for (j=0;j<a2->n_seqs;j++) {
				a->taxon_name[j+a1->n_seqs]=(char *)malloc((1+strlen(a2->taxon_name[j]))*sizeof(char));
				assert ((int) a->taxon_name[j+a1->n_seqs]);
				a->taxon_name[j+a1->n_seqs]=strcpy(a->taxon_name[j+a1->n_seqs],a2->taxon_name[j]);
				a->s[j+a1->n_seqs]=strcpy(a->s[j+a1->n_seqs],a2->s[j]);
				}

		/*this will be wrong if score==0 becuaes of the singleton variants*/
		a->score=0;
	    /*printf("In NW\n");
	    for (j=0;j<a->n_seqs;j++) {
	        for (i=0;i<a->length;i++) printf("%4d",a->s[j][i]);
	        printf("\n");
	    }*/
	}
	if (a->score && (type==2)) a->score*=values->change_cost;
	return a;
}


void    redo_for_type_2(a,input,values,number_of_data_sets)
alignment **a,*input;
int number_of_data_sets;
parameters *values;
{
  int i, *is_the_one, internal_segment, j,k,l;
  int length_of_segment;
  char *segment_holder;
  int segment,previous_length;

  if (!(values->new_optimization && (!values->in_optimize)) ){
    fprintf(stderr,"Wrong redo_for_type\n");
    exit(-1);
  }

  is_the_one=(int *)malloc(input->length*sizeof(int));
  assert((int) is_the_one);
  segment_holder=(char *)malloc(input->length*sizeof(char));
  assert((int) segment_holder);
  previous_length=0;
for (segment=0;segment<number_of_data_sets;segment++) {
  /*Allocate the alignments and copy name if segment==0*/
  a[segment]=(alignment *)malloc(sizeof(alignment));
  assert((int) a[segment]);
  a[segment]->taxon_name=(char **)malloc(sizeof(char *));
  assert((int) a[segment]->taxon_name);
  a[segment]->taxon_name[0]=(char *)malloc(2*sizeof(char));
  assert((int) a[segment]->taxon_name[0]);
  a[segment]->taxon_name[0][0]='h';
  a[segment]->taxon_name[0][1]='\0';
  a[segment]->n_seqs=1;
  if (segment==0) {
    a[segment]->name=(char *)malloc((1+(strlen(input->name)))*sizeof(char));
    assert((int) a[segment]->name);
    a[segment]->name=(char *)strcpy(a[segment]->name,input->name);
  }
  else {
    a[segment]->name=(char *)malloc(2*sizeof(char));
    assert((int) a[segment]->name);
    a[segment]->name[0]='h';
    a[segment]->name[1]='\0';
  }
  if (segment==0) is_the_one[0]=1;
  else is_the_one[0]=0;
  if (!(input->s[1][0]>63)) internal_segment=0;
  else internal_segment=1;
  for (i=1;i<input->length;i++) {
    if (internal_segment==segment) is_the_one[i]=1;
    else is_the_one[i]=0;
    if (input->s[1][i]>63) ++internal_segment;
    if ((is_the_one[i]==0) && (is_the_one[i-1]==1)) break;
  }
  length_of_segment=0;
  for (i=previous_length;i<input->length;i++) {
    if (is_the_one[i]) ++length_of_segment;
    else break;
  }
  a[segment]->s=(char **)malloc(sizeof(char *));
  assert((int) a[segment]->s);
  a[segment]->length=length_of_segment;
  a[segment]->s[0]=(char *)malloc((1+a[segment]->length)*sizeof(char));
  assert((int) a[segment]->s[0]);
  j=0;
  for (i=previous_length;i<previous_length+length_of_segment;i++) a[segment]->s[0][j++]=input->s[0][i];
  a[segment]->s[0][a[segment]->length]='\0';
  a[segment]->type_weight=0;
  previous_length+=length_of_segment;
  /*fprintf(stderr,"%d ",length_of_segment);*/
}
	free(is_the_one);
	free(segment_holder);
}

void    redo_for_type(a,segment,values)
alignment *a;
int segment;
parameters *values;
{
	int i, *is_the_one, internal_segment, j,k;
	int length_of_segment;
	char *segment_holder, **segment_holder2;

	is_the_one=(int *)malloc(a->length*sizeof(int));
	assert((int) is_the_one);

	if (segment==0) is_the_one[0]=1;
	else is_the_one[0]=0;
	if (!(a->s[a->n_seqs][0]>63)) internal_segment=0;
	else internal_segment=1;
	for (i=1;i<a->length;i++) {
		if (internal_segment==segment) is_the_one[i]=1;
		else is_the_one[i]=0;
		if (a->s[a->n_seqs][i]>63) ++internal_segment;
	}
	length_of_segment=0;
	for (i=0;i<a->length;i++) if (is_the_one[i]) ++length_of_segment;
	/*redo here down for chopped regular alignments*/
	if ((values->new_optimization) && (!values->in_optimize)) {
		segment_holder=(char *)malloc(length_of_segment*sizeof(char));
		assert((int) segment_holder);
		j=0;
		for (i=0;i<a->length;i++) if (is_the_one[i]) segment_holder[j++]=a->s[0][i];
		free(a->s[1]);
		free(a->s[0]);
		free(a->s);
		a->s=(char **)malloc(sizeof(char *));
		assert((int) a->s);
		a->length=length_of_segment;
		a->s[0]=(char *)malloc((1+a->length)*sizeof(char));
		assert((int) a->s[0]);
		for (i=0;i<a->length;i++) a->s[0][i]=segment_holder[i];
		a->s[0][a->length]='\0';
		/*a->s[0]=(char *)realloc(a->s[0],(1+strlen(a->s[0]))*sizeof(char));
		assert((int) a->s[0]);*/
		a->type_weight=0;
		free(is_the_one);
		free(segment_holder);
		}
	else /*if (values->chop>0)*/ {
		segment_holder2=(char **)malloc(a->n_seqs*sizeof(char *));
		assert((int) segment_holder2);
		for (i=0;i<a->n_seqs;i++) {
			segment_holder2[i]=(char *)malloc(length_of_segment*sizeof(char));
			assert((int) segment_holder2[i]);
			j=0;
			for (k=0;k<a->length;k++) if (is_the_one[k]) segment_holder2[i][j++]=a->s[i][k];
			free(a->s[i]);
		}
		free(a->s[a->n_seqs]);
		free(a->s);
		a->s=(char **)malloc(a->n_seqs*sizeof(char *));
		assert((int) a->s);
		a->length=length_of_segment;
		for (k=0;k<a->n_seqs;k++) {
			a->s[k]=(char *)malloc((1+a->length)*sizeof(char));
			assert((int) a->s[k]);
			j=0;
			for (i=0;i<a->length;i++) a->s[k][i]=segment_holder2[k][i];
			a->s[k][a->length]='\0';
		}
		a->type_weight=0;
		free(is_the_one);
		for (i=0;i<a->n_seqs;i++) free(segment_holder2[i]);
		free(segment_holder2);
	}
}

/*remember type wirght stuff*/
alignment *glue_back_after_chop(a,values)
alignment **a;
parameters *values;
{
	alignment *b;
	int i,j,k,l,current_start;


	b=make_align(a[0]);
	b->score*=values->data_set_weights[0];
	b->type_weight=1;
	for (i=0;i<b->n_seqs;i++) free(b->s[i]);
	free(b->s);
	b->s=(char **)malloc((1+b->n_seqs)*sizeof(char *));
	assert((int) b->s);
	for (i=1;i<values->number_of_input_alignments;i++) {
		b->length+=a[i]->length;
		b->score+=(a[i]->score*values->data_set_weights[i]);
	}
	for (i=0;i<=b->n_seqs;i++) {
		b->s[i]=(char *)malloc((1+b->length)*sizeof(char));
		assert((int) b->s[i]);
		b->s[i][b->length]='\0';
	}
/*need to reorder names why?*/

	current_start=0;
	for (i=0;i<values->number_of_input_alignments;i++) {
		for (l=0;l<b->n_seqs;l++)  {
			for (j=0;j<a[i]->length;j++) b->s[l][j+current_start]=a[i]->s[l][j];
			}

		for (j=0;j<a[i]->length;j++) {
			b->s[b->n_seqs][j+current_start]=values->data_sets[i]+1;
			if (j==(a[i]->length-1))        b->s[b->n_seqs][j+current_start] += (64);
			}
		current_start+=a[i]->length;
	}
	return b;
}

/*
int do_two_score(a,values)
alignment * a;
parameters *values;
{
int i, score;

score=0;
for (i=0;i<a->length;i++) if (a->s[0][i]!=a->s[1][i]) {
	if ((a->s[0][i]=='-') || (a->s[1][i]=='-')) score+=values->gap_cost;
	else score += values->change_cost;
	if (!((is_internal(a->s[0],i,a->length)) || (is_internal(a->s[1],i,a->length)))) score -= values->leading_gap_cost;
	}


return  score;
}

int do_three_score(a,values)
alignment * a;
parameters *values;
{
int i,j,score;

for (i=0;i<a->length;i++) {
	nstates=0;
	ngaps=0;
	for (j=0;j<3;j++) if (a->s[j][i]=='-') ++ngaps;
	if      (a->s[j][i]!='-') ++ngaps;

	}

return score;
}
*/

int only_missing(a)
alignment *a;
{
int i,j;

for (i=0;i<a->n_seqs;i++) for (j=0;j<a->length;j++) if ((a->s[i][j]!='X') && (a->s[i][j]!='N') && (a->s[i][j]!='-') && (a->s[i][j]!='?')) return 0;
return 1;
}

int only_missing_string(s)
char *s;
{
    int i;
    for (i=0;i<strlen(s);i++) if ((s[i]!='X') && (s[i]!='N') && (s[i]!='-') && (s[i]!='?')) return 0;
    return 1;
}

alignment *missing_alignment(a1,a2,values,type)
alignment *a1, *a2;
parameters *values;
int type;
{
/* type == 0 a1 missing
   type == 1 a2 missing
   type == 2 both missing*/
/* for now only new_optimization*/
/*remember==> in_optimize leave two*/

alignment *a;
int i,j;
/*fprintf(stderr,"m(%d) ",type);*/
if (!values->new_optimization) {
	a=(alignment *)malloc(sizeof(alignment));
	assert((int) a);
	if (type==0) a->length=a2->length;
	else a->length=a1->length;
	a->type_weight=0;
	if (strcmp (a1->name, a2->name) > 0) a->name = pair_names (a2, a1);
	else a->name = pair_names (a1, a2);
	a->n_seqs=a1->n_seqs+a2->n_seqs;
	a->s=(char **)malloc(a->n_seqs*sizeof(char *));
	assert((int) a->s);
	for (i=0;i<a->n_seqs;i++) {
		a->s[i]=(char *)malloc((a->length+1)*sizeof(char));
		assert ((int) a->s[i]);
		a->s[i][a->length]='\0';
		}
	a->taxon_name=(char **)malloc(a->n_seqs*sizeof(char *));
	assert((int) a->taxon_name);
	for (i=0;i<a1->n_seqs;i++) {
		a->taxon_name[i]=(char *)malloc((1+strlen(a1->taxon_name[i]))*sizeof(char));
		assert ((int) a->taxon_name[i]);
		a->taxon_name[i]=(char *)strcpy(a->taxon_name[i],a1->taxon_name[i]);
		}
	for (i=0;i<a2->n_seqs;i++) {
		a->taxon_name[i+a1->n_seqs]=(char *)malloc((1+strlen(a2->taxon_name[i]))*sizeof(char));
		assert ((int) a->taxon_name[i+a1->n_seqs]);
		a->taxon_name[i+a1->n_seqs]=(char *)strcpy(a->taxon_name[i+a1->n_seqs],a2->taxon_name[i]);
		}
	/*add the actual bases*/
	if (type==0) {
		for (j=0;j<a1->n_seqs;j++) for (i=0;i<a->length;i++) a->s[j][i]='X';
		for (j=0;j<a2->n_seqs;j++) for (i=0;i<a->length;i++) a->s[j+a1->n_seqs][i]=a2->s[j][i];
		}
	else if (type==1) {
		for (j=0;j<a1->n_seqs;j++) for (i=0;i<a->length;i++) a->s[j][i]=a1->s[j][i];
		for (j=0;j<a2->n_seqs;j++) for (i=0;i<a->length;i++) a->s[j+a1->n_seqs][i]='X';
		}
	else for (j=0;j<a->n_seqs;j++) for (i=0;i<a->length;i++) a->s[j][i]='X';
	}
else {
	a=(alignment *)malloc(sizeof(alignment));
	assert((int) a);
	if (type==0) a->length=a2->length;
	else a->length=a1->length;
	a->type_weight=0;
	if (strcmp (a1->name, a2->name) > 0) a->name = pair_names (a2, a1);
	else a->name = pair_names (a1, a2);
	a->n_seqs=1+values->in_optimize;
	a->s=(char **)malloc((1+values->in_optimize)*sizeof(char *));
	assert((int) a->s);
	a->s[0]=(char *)malloc((a->length+1)*sizeof(char));
	assert ((int) a->s[0]);
	a->s[0][a->length]='\0';
	if (values->in_optimize) {
		a->s[1]=(char *)malloc((a->length+1)*sizeof(char));
		assert ((int) a->s[1]);
		a->s[1][a->length]='\0';
		}
	a->taxon_name=(char **)malloc((1+values->in_optimize)*sizeof(char *));
	assert((int) a->taxon_name);
	if (!values->in_optimize) {
		a->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
		assert ((int) a->taxon_name[0]);
		a->taxon_name[0][0]='H';
		a->taxon_name[0][1]='y';
		a->taxon_name[0][2]='p';
		a->taxon_name[0][3]='A';
		a->taxon_name[0][4]='n';
		a->taxon_name[0][5]='c';
		a->taxon_name[0][6]='\0';
		}
	else {
		a->taxon_name[0]=(char *)malloc((1+strlen(a1->taxon_name[0]))*sizeof(char));
		assert ((int) a->taxon_name[0]);
		a->taxon_name[1]=(char *)malloc((1+strlen(a2->taxon_name[0]))*sizeof(char));
		assert ((int) a->taxon_name[1]);
		a->taxon_name[0]=(char *)strcpy(a->taxon_name[0],a1->taxon_name[0]);
		a->taxon_name[1]=(char *)strcpy(a->taxon_name[1],a2->taxon_name[0]);
		}
	/*set the actual bases*/
	if (type==0) {
		if (!values->in_optimize) for (i=0;i<a->length;i++) a->s[0][i]=a2->s[0][i];
		else {
			for (i=0;i<a->length;i++) {
				a->s[0][i]=a2->s[0][i];
				a->s[1][i]='X';
				}
			}
		}
	else if (type==1) {
		if (!values->in_optimize) for (i=0;i<a->length;i++) a->s[0][i]=a1->s[0][i];
		else {
			for (i=0;i<a->length;i++) {
				a->s[0][i]='X';
				a->s[1][i]=a1->s[0][i];
				}
			}
		}
	else {
		if (!values->in_optimize) for (i=0;i<a->length;i++) a->s[0][i]='X';
		else for (i=0;i<a->length;i++) a->s[0][i]=a->s[1][i]='X';
		}
	}
a->score=0;
/*fprintf(stderr,"%s %s / ",a->name,a->s[0]);
if (values->in_optimize) fprintf(stderr,"%s\n",a->s[1]);*/

/*if (values->in_optimize) fprintf(stderr,"Uping");*/
return a;
}

int no_align_column_score(a,values)
parameters *values;
alignment *a;
{
    int i,j,A,C,G,T,GAP;
    a->score=0;
    for (i=0;i<a->length;i++) {
        A=C=G=T=GAP=0;
        for (j=0;j<a->n_seqs;j++) {
                switch (a->s[j][i]) {
                    case 'A' : A=1; break;
                    case 'C' : C=1; break;
                    case 'G' : G=1; break;
                    case 'T' : T=1; break;
                    case '-' : GAP=1; break;
                    case '0' : A=1; break;
                    case '1' : C=1; break;
                    case '2' : G=1; break;
                    case '3' : T=1; break;
                }
        }
        a->score+=(((A+C+G+T-1)*values->change_cost)+(GAP*values->gap_cost));
    }
return a->score;
}


