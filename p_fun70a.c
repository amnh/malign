/*
Copyright 1992 Ward Wheeler all rights reserved

function to get tree lengths
need to redo the Sankoff stuff
change not to redo like for alignments
     1) nodes pointer or something keeps tree description above
  2) need to make leading and trailing gaps count
     at least in an option
    currently screwing up BandB
    two new classes of change or a weight vector if easier to code
*/

#include "align3.h"

int aut;
change_counter_structure *changes;
change_counter_structure *new_changes;
char *pre_tree_made;
int *num_seqs;
int *new_num_seqs;
int *cc_weights;
/*int return_value;*/


#define THANG holder=0;\
	      ls1=strlen(changes[i].name);\
	      if (ls1 != strlen(new_changes[j].name)) holder=1;\
	      else for (ii=0;ii<ls1;ii++) if (changes[i].name[ii] != new_changes[j].name[ii]) {\
	      holder=1;\
	      break;}\


int cladogram(a,values)
alignment *a;
parameters *values;
{
	int i,m,j,l,ntaxa,k;
	int **nodes;
	look_up_thang **sequence;
	int best_tree,cur_b_tree;
	int *tree_rep,n_gaps=0;
	char *p,*igaps;
	int holder,h1,h2,h3,h4,h5,h6;
	int alpha,tree_rand_order,temp_tree;
	char name_buff[10];
	alignment **a1;
	int **pieces,total_length,*weights,same,yet_newer_length;
	parameters *valuesI;
	int *symaut,state[7],nstates,length_holder,single,new_length;
	int m_holder,all_same,k2;


	GOMAC
	say("Recoding data")
	if (a->n_seqs==1) return 0;
        if (a->n_seqs < 4) if (values->number_of_input_alignments<2) return a->score;
	if (values->phylo_time>-1) alpha=(int) time(NULL);
	sequence=NULL;
	pieces=NULL;
	weights=NULL;
	valuesI=NULL;
	symaut=NULL;
	values->tree_made=NULL;

	if (values->phylo_gap==0) {
		holder=values->gap_cost;
		values->gap_cost=0;
		}

	/*reading data*/
	ntaxa=a->n_seqs;
	cc_weights=NULL;


aut=0;
if (values->number_of_input_alignments<2) {
	symaut=(int *)malloc(a->length*sizeof(int));
	assert((int) symaut);
	m=0;
	for (i=0;i<a->length;i++) {
	    all_same=1;
	    for (j=1;j<a->n_seqs;j++) if (a->s[j][i]!=a->s[0][i]) {all_same=0;break;}
	    if (!all_same) {++m;symaut[i]=0;}
	    else symaut[i]=1;
	}
        if (m==0) {
                best_tree=0;
                goto end_of_cladogram;
                }
	sequence  = (look_up_thang **) malloc(a->n_seqs * sizeof (look_up_thang *));
	assert((int)sequence);
	for (i=0;i<a->n_seqs;i++) {
	    sequence[i]=(look_up_thang *)malloc(m*sizeof(look_up_thang));
	    assert((int) sequence[i]);
	    m_holder=0;
	    for (j=0;j<a->length;j++) if (!symaut[j]) {
	            switch (a->s[i][j]) {
						case 'A': sequence[i][m_holder].base=A_BIT; break;
						case 'C': sequence[i][m_holder].base=C_BIT; break;
						case 'G': sequence[i][m_holder].base=G_BIT; break;
						case 'T': sequence[i][m_holder].base=T_BIT; break;
						case '-': {
							if (values->phylo_gap) {
								if (is_internal(a->s[i],j,a->length)) sequence[i][m_holder].base=GAP_BIT;
								else sequence[i][m_holder].base=X_BIT;
								}
							else sequence[i][m_holder].base=X_BIT;
							break;
							}
						case 'N': sequence[i][m_holder].base=N_BIT; break;
						case 'X': sequence[i][m_holder].base=X_BIT; break;
						case 'U': sequence[i][m_holder].base=U_BIT; break;
						case 'R': sequence[i][m_holder].base=R_BIT; break;
						case 'Y': sequence[i][m_holder].base=Y_BIT; break;
						case 'M': sequence[i][m_holder].base=M_BIT; break;
						case 'W': sequence[i][m_holder].base=W_BIT; break;
						case 'S': sequence[i][m_holder].base=S_BIT; break;
						case 'K': sequence[i][m_holder].base=K_BIT; break;
						case 'B': sequence[i][m_holder].base=B_BIT; break;
						case 'D': sequence[i][m_holder].base=D_BIT; break;
						case 'H': sequence[i][m_holder].base=H_BIT; break;
						case 'V': sequence[i][m_holder].base=V_BIT; break;
						default:
							fprintf (stderr, "Bad character %c in sequence[%d][%d] of input %d\n",a->s[i][j], i,j);
							print_alignment (a);
							exit(-1);
						}
				/*sequence[i][m_holder].cost=sequence[i][m_holder].union_bit=0;*/
        	    m_holder++;
                }
  	}


	best_tree=HUGE_COST;
	if ((m)==0) {
		best_tree=0;
		goto end_of_cladogram;
	}
	 weights=(int *)malloc((m)*sizeof(int));
	assert((int) weights);
	for (i=0;i<m;i++) weights[i]=1;
    if (values->jackboot) for (j=0;j<m;j++) weights[j]*=values->jack_array[(j*1000)/m];
	/*crunch*/
	if (1) {
		for (j=0;j<m;j++) if (weights[j]>0) {
			for (k=j+1;k<m;k++) if (weights[k]>0) {
				same=1;
				for (i=0;i<a->n_seqs;i++) if (sequence[i][j].base!=sequence[i][k].base) { same=0; break;}
				if (same) {
					weights[j]+=weights[k];
					weights[k]=0;
					}
				}
			}

		new_length=0;
		for (j=0;j<m;j++) if (weights[j]>0) ++new_length;
		}
	else new_length=m;
	/*remove autapomorphies*/
	for (j=0;j<m;j++) if (weights[j]>0) {
	    for (i=0;i<7;i++) state[i]=0;
	    for (i=0;i<a->n_seqs;i++) {
	        if (sequence[i][j].base==A_BIT) state[0]+=1;
	        else if (sequence[i][j].base==C_BIT) state[1]+=1;
	        else if (sequence[i][j].base==G_BIT) state[2]+=1;
	        else if (sequence[i][j].base==T_BIT) state[3]+=1;
	        else if (sequence[i][j].base==GAP_BIT) state[4]+=1;
	    }
	    nstates=0;
	    for (i=0;i<5;i++) if (state[i]>0) ++nstates;
	    if (nstates==2) {
	        same=0;
	        for (i=0;i<5;i++) if (state[i]==1) same=1;
	        if (same) {/*is autoapomorphy*/
	            symaut[j]=1;
	            h1=h2= (-1);
	            for (i=0;i<5;i++) if (state[i]>0) {
	                if (h1<0) h1=i;
	                else h2=i;
	            }
                if (h1==0) h1=A_BIT;
                else if (h1==1) h1=C_BIT;
                else if (h1==2) h1=G_BIT;
                else if (h1==3) h1=T_BIT;
                else if (h1==4) h1=GAP_BIT;
                if (h2==0) h2=A_BIT;
                else if (h2==1) h2=C_BIT;
                else if (h2==2) h2=G_BIT;
                else if (h2==3) h2=T_BIT;
                else if (h2==4) h2=GAP_BIT;
                aut+= weights[j]*values->lookup[h1][h2].cost;
	            weights[j]=0;
	        }
	    }
	}
	new_length=0;
    for (j=0;j<m;j++) if (weights[j]>0) ++new_length;
        if (new_length==0) {
                best_tree=0;
                free(weights);
                goto end_of_cladogram;
                }
	cc_weights=(int *)malloc((new_length)*sizeof(int));
	assert((int) cc_weights);
	m_holder=0;
	for (j=0;j<m;j++) if (weights[j]>0) cc_weights[m_holder++]=weights[j];
	for (i=0;i<a->n_seqs;i++) {
	    m_holder=0;
	    for (j=0;j<m;j++) if (weights[j]>0) sequence[i][m_holder++]=sequence[i][j];
	    sequence[i]=(look_up_thang *)realloc(sequence[i],new_length*sizeof(look_up_thang));
	    assert((int) sequence[i]);
	}
	free(weights);
    m=new_length;
    values->start=(int *)malloc(sizeof(int));
    assert((int) values->start);
    values->stop=(int *)malloc(sizeof(int));
    assert((int) values->stop);
    values->start[0]=0;
    values->stop[0]=m;
    /*fprintf(stderr,"length=%d vs. %d\n",m,a->length);
    for (j=0;j<m;j++) fprintf(stderr,"(%d %d)",j,cc_weights[j]);*/
	}
else {
	values->start=(int *)malloc(values->number_of_input_alignments*sizeof(int));
    assert((int) values->start);
    values->stop=(int *)malloc(values->number_of_input_alignments*sizeof(int));
    assert((int) values->stop);
    /*split*/
	a1=(alignment **)malloc(values->number_of_input_alignments*sizeof(alignment *));
	assert((int) a1);
	pieces=(int **)malloc(a->n_seqs*sizeof (int *));
	assert((int)pieces);
	for (i=0;i<a->n_seqs;i++) {
		pieces[i]=(int *)malloc(sizeof(int));
		assert((int) pieces[i]);
		}
	weights=(int *)malloc(a->length*sizeof(int));
	assert((int) weights);
	total_length=0;
	for (i=0;i<values->number_of_input_alignments;i++) {
	    values->start[i]=total_length;
	    if (!values->other_parm) valuesI=values;
		else if (!values->other_parm[i]) valuesI=values;
		else valuesI=values->other_parm[i];
		a1[i]=NULL;
		a1[i]=make_align(a);
		redo_for_type(a1[i],i,values);
	   /*assign raw states*/
		if (values->data_sets[i]==2) { /*morph*/
			for (j=0;j<a->n_seqs;j++) {
				pieces[j]=(int *)realloc(pieces[j],(total_length+a1[i]->length)*sizeof(int));
				assert((int) pieces[j]);
				}
			weights=(int *)realloc(weights,(total_length+a1[i]->length)*sizeof(int));
			assert((int) weights);
			for (j=total_length;j<total_length+a1[i]->length;j++) weights[j]=values->data_set_weights[i];
			for (j=0;j<a1[i]->n_seqs;j++) for (k=0;k<a1[i]->length;k++) {
				if (a1[i]->s[j][k]=='0') pieces[j][k+total_length]=1;
				else if (a1[i]->s[j][k]=='1') pieces[j][k+total_length]=2;
				else if (a1[i]->s[j][k]=='2') pieces[j][k+total_length]=4;
				else if (a1[i]->s[j][k]=='3') pieces[j][k+total_length]=8;
				else if (a1[i]->s[j][k]=='4') pieces[j][k+total_length]=16;
				else if (a1[i]->s[j][k]=='5') pieces[j][k+total_length]=32;
				else if (a1[i]->s[j][k]=='6') pieces[j][k+total_length]=64;
				else if (a1[i]->s[j][k]=='?') pieces[j][k+total_length]=127;
				else if (a1[i]->s[j][k]=='-') pieces[j][k+total_length]=127;
				else pieces[j][k+total_length]=127;
				}
			/*crunch it down*/
			if (symaut) free(symaut);
			symaut=(int *)malloc(a1[i]->length*sizeof(int));
	        assert((int) symaut);
	        for (j=0;j<a1[i]->length;j++) {
	                symaut[j]=1;
	                for (k=1;k<a1[i]->n_seqs;k++) if (a1[i]->s[0][j]!=a1[i]->s[k][j]) {symaut[j]=0;break;}
	        }
	        /*jackboot stuff*/
	        if (values->jackboot) for (j=0;j<a1[i]->length;j++) {
	                weights[j]*=values->jack_array[(j*1000)/a1[i]->length];
	                if (!weights[j]) symaut[j]=1;
	        }

	        /*for (j=0;j<a1[i]->length;j++) if (!symaut[j]) {
	            for (l=j+1;l<a1[i]->length;l++) if (!symaut[l]) {
	                same=1;
	                for (k=0;k<a1[i]->n_seqs;k++) if (a1[i]->s[k][j]!=a1[i]->s[k][l]) {same=0;break;}
	                if (same) {
	                        weights[total_length+j]+=weights[total_length+l];
	                        weights[total_length+l]=0;
	                        symaut[l]=1;
	                }
	            }
	        }
	        new_length=0;
	        for (j=0;j<a1[i]->length;j++)  if (!symaut[j]) ++new_length;
	        for (j=0;j<a1[i]->n_seqs;j++) {
	                m_holder=0;
	                for (k=0;k<a1[i]->length;k++) if (!symaut[k]) pieces[j][total_length+(m_holder++)]=pieces[j][total_length+k];
	                pieces[j]=(int *)realloc(pieces[j],(total_length+new_length)*sizeof(int));
	                assert((int) pieces[j]);
	        }
             m_holder=0;
	         for (k=0;k<a1[i]->length;k++) if (!symaut[k]) weights[total_length+(m_holder++)]=weights[total_length+k];
	         weights=(int *)realloc(weights,(total_length+new_length)*sizeof(int));
	         assert((int) weights);*/
	         free(symaut);
	         symaut=NULL;
			total_length+=new_length;
	    }
		else { /*Molecular Data*/
		    /*crunch it down*/
			if (symaut) free(symaut);
			symaut=(int *)malloc(a1[i]->length*sizeof(int));
	        assert((int) symaut);
	        for (j=0;j<a1[i]->length;j++) {
	                symaut[j]=1;
	                for (k=1;k<a1[i]->n_seqs;k++) if (a1[i]->s[0][j]!=a1[i]->s[k][j]) {symaut[j]=0;break;}
	        }
	        /*jackboot stuff*/
	        if (values->jackboot) for (j=0;j<a1[i]->length;j++) if (!symaut[j]) symaut[j]=(!(values->jack_array[(j*1000)/a1[i]->length]));
	        weights=(int *)realloc(weights,(total_length+(a1[i]->length))*sizeof(int));
			assert((int) weights);
			for (j=total_length;j<total_length+a1[i]->length;j++) weights[j]=values->data_set_weights[i];
			for (j=0;j<a1[i]->length;j++) if (!symaut[j]) {
	            for (l=j+1;l<a1[i]->length;l++) if (!symaut[l]) {
	                same=1;
	                for (k=0;k<a1[i]->n_seqs;k++) if (a1[i]->s[k][j]!=a1[i]->s[k][l]) {same=0;break;}
	                if (same) {
	                        weights[total_length+j]+=weights[total_length+l];
	                        weights[total_length+l]=0;
	                        symaut[l]=1;
	                }
	            }
	        }
	        new_length=0;
	        for (j=0;j<a1[i]->length;j++)  if (!symaut[j]) ++new_length;

            for (j=0;j<a->n_seqs;j++) {
					pieces[j]=(int *)realloc(pieces[j],(total_length+new_length)*sizeof(int));
					assert((int) pieces[j]);
					}
			for (j=0;j<a1[i]->n_seqs;j++) {
			    k2=0;
			    for (k=0;k<a1[i]->length;k++) if (!symaut[k]) {
				switch (a1[i]->s[j][k]) {
					case 'A': pieces[j][k2+total_length]=A_BIT; break;
					case 'C': pieces[j][k2+total_length]=C_BIT; break;
					case 'G': pieces[j][k2+total_length]=G_BIT; break;
					case 'T': pieces[j][k2+total_length]=T_BIT; break;
					case '-': {
						if (values->phylo_gap) {
							if (is_internal(a1[i]->s[j],k,a1[i]->length)) pieces[j][k2+total_length]=GAP_BIT;
							else pieces[j][k2+total_length]=X_BIT;
							}
						else pieces[j][k2+total_length]=X_BIT;
						break;
						}
					case 'N': pieces[j][k2+total_length]=N_BIT; break;
					case 'X': pieces[j][k2+total_length]=X_BIT; break;
					case 'U': pieces[j][k2+total_length]=U_BIT; break;
					case 'R': pieces[j][k2+total_length]=R_BIT; break;
					case 'Y': pieces[j][k2+total_length]=Y_BIT; break;
					case 'M': pieces[j][k2+total_length]=M_BIT; break;
					case 'W': pieces[j][k2+total_length]=W_BIT; break;
					case 'S': pieces[j][k2+total_length]=S_BIT; break;
					case 'K': pieces[j][k2+total_length]=K_BIT; break;
					case 'B': pieces[j][k2+total_length]=B_BIT; break;
					case 'D': pieces[j][k2+total_length]=D_BIT; break;
					case 'H': pieces[j][k2+total_length]=H_BIT; break;
					case 'V': pieces[j][k2+total_length]=V_BIT; break;
					default:
						fprintf (stderr, "Bad character %c in sequence[%d][%d] of input %d\n",a1[i]->s[j][k], j, k,i);
						print_alignment (a1[i]);
						exit(-1);
					}
				k2++;
				}
			}
        	/*remove autapomorphies*/
        	for (j=0;j<new_length;j++) if (!symaut[j]) {
        	    for (k=0;k<7;k++) state[k]=0;
        	    for (k=0;k<a1[i]->n_seqs;k++) {
        	        if (pieces[k][j+total_length]==A_BIT) state[0]+=1;
        	        else if (pieces[k][j+total_length]==C_BIT) state[1]+=1;
        	        else if (pieces[k][j+total_length]==G_BIT) state[2]+=1;
        	        else if (pieces[k][j+total_length]==T_BIT) state[3]+=1;
        	        else if (pieces[k][j+total_length]==GAP_BIT) state[4]+=1;
        	    }
        	    nstates=0;
        	    for (k=0;k<5;k++) if (state[k]>0) ++nstates;
        	    if (nstates==2) {
        	        same=0;
        	        for (k=0;k<5;k++) if (state[k]==1) same=1;
        	        if (same) {/*is autoapomorphy*/
        	            symaut[j]=1;
        	            h1=h2= (-1);
        	            for (k=0;k<5;k++) if (state[k]>0) {
        	                if (h1<0) h1=k;
        	                else h2=k;
        	            }
                        if (h1==0) h1=A_BIT;
                        else if (h1==1) h1=C_BIT;
                        else if (h1==2) h1=G_BIT;
                        else if (h1==3) h1=T_BIT;
                        else if (h1==4) h1=GAP_BIT;
                        if (h2==0) h2=A_BIT;
                        else if (h2==1) h2=C_BIT;
                        else if (h2==2) h2=G_BIT;
                        else if (h2==3) h2=T_BIT;
                        else if (h2==4) h2=GAP_BIT;
                        aut+= weights[total_length+j]*values->lookup[h1][h2].cost;
        	            weights[total_length+j]=0;
        	        }
        	    }
        	}
        	new_length=0;
            for (j=0;j<a1[i]->length;j++) if (!symaut[j]) ++new_length;

            /*fprintf(stderr,"length=%d vs. %d\n",new_length,a1[i]->length);*/
	        m_holder=0;
	        for (k=0;k<a1[i]->length;k++) if (!symaut[k]) weights[total_length+(m_holder++)]=weights[total_length+k];
	        weights=(int *)realloc(weights,(total_length+new_length)*sizeof(int));
	        free(symaut);
	        symaut=NULL;
			total_length+=new_length;
    		}
	    values->stop[i]=total_length;
	    }

	/*allocate sequences and copy over*/
	sequence  = (look_up_thang **) malloc(a->n_seqs * sizeof (look_up_thang *));
	assert((int)sequence);
	for (i=0;i<a->n_seqs;i++) {
		sequence[i]=(look_up_thang *)malloc(total_length*sizeof(look_up_thang));
		assert((int) sequence[i]);
		for (j=0;j<total_length;j++) sequence[i][j].base=pieces[i][j];
		}
	/*create new weights and copy over*/
	/*cc_weights=(int *)malloc(total_length*sizeof(int));
	assert((int) cc_weights);
	for (j=0;j<total_length;j++) cc_weights[j]=weights[j];*/
	cc_weights=weights;
	weights=NULL;


	/*deallocate local stuff*/
	if (weights) free(weights);
	for (i=0;i<a->n_seqs;i++) free(pieces[i]);
	free(pieces);
	if (symaut) free(symaut);
	for (i=0;i<values->number_of_input_alignments;i++) if (a1[i]) a1[i]=dump_align(a1[i]);
	free(a1);
	m=total_length;
	}
/*for (i=0;i<values->number_of_input_alignments;i++) fprintf(stderr,"(%d->%d)",values->start[i],values->stop[i]);*/
best_tree=HUGE_COST;
if (a->n_seqs<4) {
	if (a->n_seqs==2) best_tree=0;
        else if (cc_weights) {
		best_tree=0;
		for (i=0;i<m;i++) best_tree+=(2*cc_weights[i]);
		}
        else best_tree=0;
	goto end_of_cladogram;
	}
if (values->phylo_score > 6) {
		if (!values->ttr && (!values->delta)) values->n_transv=0;
		for (tree_rand_order=0;tree_rand_order<max(1,values->tree_rand_order_max);tree_rand_order++) {
			if (values->tree_rand_order_max>0) if (tree_rand_order>0) reorder_tree(sequence,ntaxa,m,n_gaps);
			if (m==0) temp_tree=0;
			else temp_tree=daves_cladogram(a->n_seqs,values,sequence,m,cc_weights);
			if (best_tree > temp_tree) best_tree=temp_tree;
			}
		}
else {fprintf(stderr,"Other tree construction procedures under revision\n");}
end_of_cladogram:;
if (values->phylo_time>-1)    values->phylo_time+=(((int)time(NULL))-alpha);
	if (cc_weights) free(cc_weights);
	for (i = 0; i < ((a->n_seqs)); i++) if (sequence[i]) free(sequence[i]);
	if (sequence) free(sequence);
	if (values->tree_made) free(values->tree_made);
    free(values->start);
    free(values->stop);

	if (num_seqs) free(num_seqs);
	if (pre_tree_made) free(pre_tree_made);
	if (new_num_seqs) free(new_num_seqs);
	if (values->phylo_gap==0) values->gap_cost=holder;
	return (best_tree+aut);
}/*end cladogram*/

void loop(tree_rep,element,ntax,nodes,nbases,sequence,best_tree,l,cur_b_tree,n_gaps,values)
int element,ntax,nbases;
int *tree_rep;
int **nodes;
int *best_tree,n_gaps;
int **sequence,*l;
int *cur_b_tree;
parameters *values;
{
	int i,cur_tree;

	for (i=((2*(element+4))-5);i>=1;--i)
	{
		tree_rep[element]=i;
		if (element<(ntax-4)) {
			GOMAC
			    make_nodes(tree_rep,nodes,ntax,nbases,sequence,element+4,n_gaps,values,&cur_tree);
			if (cur_tree>(*cur_b_tree)) goto end_of_loop;/*because of autoapomaorphy removal, otherwise >= */
			loop(tree_rep,element+1,ntax,nodes,nbases,sequence,best_tree,l,cur_b_tree,n_gaps,values);
		}
		else {
			tree_rep[element]=i;
			make_nodes(tree_rep,nodes,ntax,nbases,sequence,ntax,n_gaps,values,&cur_tree);
			/*      cur_tree=return_value;*/
			if (cur_tree < (*cur_b_tree)) {
				*cur_b_tree = *best_tree=cur_tree;
				/*
     values->length_best_clad=cur_tree;
	       values->number_best_clad=1;
	       for (k=0;k<ntax-3;k++) values->best_rep[0][k]=tree_rep[k];
	       */
			}
			/*
	  if (cur_tree == (*cur_b_tree)) {
		    if (values->number_best_clad < ((values->keep_trees)-1)) {
		    for (k=0;k<ntax-3;k++) values->best_rep[values->number_best_clad][k]=tree_rep[k];
		    ++values->number_best_clad;
		    }
	       else say("Overflow in exact solution")
	       }
     */
		}
end_of_loop:
		;
	}
}


int make_nodes_weights(tree_rep,node,ntaxa,nbases,sequence,nent,n_gaps,values)
int **node;
int *tree_rep;
int ntaxa,nbases,nent;
int **sequence,n_gaps;
parameters *values;
{
	int i,dummy1,dummy2,tree_length;

	node[0][0]=0;
	node[0][1]=ntaxa+1;
	node[1][0]=1;
	node[1][1]=2;

	/*node names start after the taxa*/
	for (i=3;i<nent;i++)
	{
		node[i-1][0]=i;
		dummy1=(tree_rep[i-3])/2;
		dummy2=(tree_rep[i-3])%2;
		node[i-1][1]=node[dummy1][dummy2];
		node[dummy1][dummy2]=ntaxa+i-1;
	}
	/* order do here*/
	tree_length=diagnose_tree_and_get_weights(sequence,node,nbases,ntaxa,nent,n_gaps,values);
	return tree_length;
}

int get_bound(sequence,ntax,nbases,n_gaps,nodes,tree_rep,randtrees,values)
int **sequence,randtrees;
int ntax,nbases,n_gaps;
int *tree_rep;
int **nodes;
parameters *values;
{
	int i,j,k,l,rand_tree,best;
	int swap,best_swap;
	int *in_rep,counter;
	int counter2,temp1,temp2;
	b_tree **b_swap,**b_swap2;

	in_rep=(int *)malloc((ntax-3)*sizeof(int));
	assert((int)in_rep);
	b_swap = (b_tree **) malloc(values->keep_trees * sizeof (b_tree *));
	assert((int)b_swap);
	b_swap2 = (b_tree **) malloc(values->keep_trees * sizeof (b_tree *));
	assert((int)b_swap2);
	for (i = 0; i < values->keep_trees; i++) {
		b_swap[i] = (b_tree *)malloc (sizeof(b_tree));
		assert((int)b_swap[i]);
		b_swap[i]->rep=(int *)malloc((ntax-3)*sizeof(int));
		assert((int)b_swap[i]->rep);
		b_swap2[i] = (b_tree *)malloc (sizeof(b_tree));
		assert((int)b_swap2[i]);
		b_swap2[i]->rep=(int *)malloc((ntax-3)*sizeof(int));
		assert((int)b_swap2[i]->rep);
	}

	best=HUGE_COST;
	counter=0;
	/*get random trees and lengths*/
	rand_tree=HUGE_COST;
	for (j=0;j<randtrees;j++){
		for (i=0;i<=ntax-4;i++){
			temp1=(int) rand()/(2*(i+4));
			temp2=(int) max_rand/(2*(i+4));
			tree_rep[i]=(int) ((((2*(i+4))-5)*temp1)/temp2);
			if (tree_rep[i]!=((2*(i+4))-5)) ++tree_rep[i];
		}

		make_nodes(tree_rep,nodes,ntax,nbases,sequence,ntax,n_gaps,values,&rand_tree);
		/*    rand_tree=return_value;*/
		if ((rand_tree==best) && (!already_clad(values,in_rep,ntax))) {
			if (counter==(values->keep_trees-1)) if (values->rep_error) fprintf(stderr,"Overflow in cladogram construction -- but that's OK.\n");
			else {
				++counter;
				for (i=0;i<ntax-3;i++) b_swap[counter]->rep[i]=values->best_rep[values->number_best_clad][i]=tree_rep[i];
				++values->number_best_clad;
			}
		}
		else if (rand_tree<best) {
			values->length_best_clad=best=rand_tree;
			counter=0;
			for (i=0;i<ntax-3;i++) b_swap[counter]->rep[i]=values->best_rep[0][i]=tree_rep[i];
			values->number_best_clad=1;
		}
	}
	if (values->tree_swap==1) {
		say("Swapping")
		    /*swap on rand*/
		GOMAC
		    best_swap=HUGE_COST;
		if (randtrees > 0) {
SWAP:
			counter2=0;
			for (k=0;k<=counter;k++) {
				for (i=0;i<ntax-3;i++) in_rep[i]=b_swap[k]->rep[i];
				for (i=0;i<ntax-3;i++){
					for (j=1;j<=((2*(i+4))-5);j++){
						in_rep[i]=j;
						make_nodes(in_rep,nodes,ntax,nbases,sequence,ntax,n_gaps,values,&swap);
						/*          swap=return_value;*/
						if ((swap==best_swap) && (!already_clad(values,in_rep,ntax)) && (counter2<(values->keep_trees-1))) {
							if (counter2==(values->keep_trees-1)) if (values->rep_error) fprintf(stderr,"Overflow in cladogram construction -- but that's OK.");
							else {
								++counter2;
								for (l=0;l<ntax-3;l++) b_swap2[counter2]->rep[l]=values->best_rep[values->number_best_clad][l]=in_rep[l];
								++values->number_best_clad;
							}
						}
						else if (swap<best_swap) {
							values->length_best_clad=best_swap=swap;
							counter2=0;
							for (l=0;l<ntax-3;l++) b_swap2[counter2]->rep[l]=values->best_rep[0][l]=in_rep[l];
							values->number_best_clad=1;
						}
					}
					in_rep[i]=b_swap[k]->rep[i];
				}
			}
			/*fprintf(stderr,"swap %d ",best_swap);*/
			if (best_swap<best) {
				best=best_swap;
				for (k=0;k<=counter2;k++) {
					for (i=0;i<ntax-3;i++) b_swap[k]->rep[i]=b_swap2[k]->rep[i];
				}
				counter=counter2;
				goto SWAP;
			}
		}
	}/*end swap*/
	/*fprintf(stderr,"After swap %d\n",best_swap);*/
	/*Build tree sensibly*/

	for (i=0;i<values->keep_trees;i++) {
		free(b_swap[i]->rep);
		free(b_swap2[i]->rep);
		free(b_swap[i]);
		free(b_swap2[i]);
	}
	free(b_swap);
	free(b_swap2);
	free(in_rep);
	return best;
}

void shuffle(sequence,ntaxa,m,ngaps)
int **sequence;
int ntaxa,m,ngaps;
{
	int i;
	int *temp;

	temp=(int *)malloc((m+ngaps)*sizeof(int));
	for (i=0;i<(m+ngaps);i++) *(temp+i)=(*(sequence[1]+i));
	for (i=0;i<(m+ngaps);i++) *(sequence[1]+i) = (*(sequence[ntaxa-1]+i));
	for (i=0;i<(m+ngaps);i++) *(sequence[ntaxa-1]+i) = (*(temp+i));

	if (ntaxa>5) {
		for (i=0;i<(m+ngaps);i++) *(temp+i)=(*(sequence[2]+i));
		for (i=0;i<(m+ngaps);i++) *(sequence[2]+i) = (*(sequence[ntaxa/2]+i));
		for (i=0;i<(m+ngaps);i++) *(sequence[ntaxa/2]+i) = (*(temp+i));
	}

	free(temp);
}


int **deep_swap_tree(cur_b_tree,sequence,n_bases,num_to_swap,old_reps,num_best,ntaxa,depth,values,n_gaps)
int num_to_swap,**old_reps,*num_best,ntaxa,depth,n_gaps;
int **sequence,n_bases,*cur_b_tree;
parameters *values;
{
	int i,j,new_ones=0,overflow=0;
	char *to_swap,*tuple_holder,**already_swapped;
	int **best_nodes,**nodes;
	realloc_holder *passer;

	GOMAC
	    say("Swapping")

	passer=(realloc_holder *)malloc(sizeof(realloc_holder));
	assert((int)passer);
	to_swap=(char *)malloc((ntaxa-3)*sizeof(char));
	assert((int)to_swap);
	tuple_holder=(char *)malloc((ntaxa-3)*sizeof(char));
	assert((int)tuple_holder);
	already_swapped=(char **)malloc(values->keep_trees*sizeof(char *));
	assert((int)already_swapped);
	for (i=0;i<values->keep_trees;i++) {
		already_swapped[i]=(char *)malloc((ntaxa-3)*sizeof(char));
		assert((int)already_swapped[i]);
	}

	nodes=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)nodes);
	best_nodes=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)best_nodes);
	for (i=0;i<(ntaxa-1);i++) {
		nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int)nodes[i]);
		best_nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int)best_nodes[i]);
	}
	nodes[0][0]=0;
	nodes[0][1]=ntaxa+1;
	nodes[1][0]=1;
	nodes[1][1]=2;

	for (i=0;i<num_to_swap;i++) {
		for (j=0;j<ntaxa-4;j++) *(already_swapped[i]+j)=0;
		*(already_swapped[i]+ntaxa-4)=1;
	}
	/*need to add previous swap info i.e. from so if tree_rep!=from etc*/
	for (i=0;i<num_to_swap+new_ones;i++) {
		for (j=0;j<ntaxa-3;j++) *(tuple_holder+j)=0;
		passer=get_tuples_and_loop_tree(&overflow,i,already_swapped,best_nodes,nodes,cur_b_tree,sequence,n_bases,num_best,old_reps[i],&new_ones,num_to_swap,old_reps,0,depth,tuple_holder,ntaxa,values,n_gaps,passer);
		old_reps=passer->rep;
		already_swapped=passer->swap;
	}

	free(to_swap);
	free(tuple_holder);
	for (i=0;i<values->keep_trees;i++) if (already_swapped[i]) free(already_swapped[i]);
	free(already_swapped);
	for (i=0;i<(ntaxa-1);i++) {
		free(nodes[i]);
		free(best_nodes[i]);
	}
	free(nodes);
	free(best_nodes);
	free(passer);

	return old_reps;
}

realloc_holder *get_tuples_and_loop_tree(overflow,alignment_swapping,already_swapped,best_nodes,nodes,cur_b_tree,sequence,n_bases,num_best,current_rep,new_ones,num_to_swap,old_reps,current_depth,depth,tuple_holder,ntaxa,values,n_gaps,passer)
int *num_best,*current_rep,*new_ones,num_to_swap,**old_reps,current_depth,depth,ntaxa;
char *tuple_holder;
parameters *values;
int *cur_b_tree,alignment_swapping,**sequence,n_bases,n_gaps;
int **best_nodes,**nodes;
char **already_swapped;
int *overflow;
realloc_holder *passer;
{
	int i,j,num_tuples;
	int *tree_rep;
	char doit;

	/*this needs to be checked as far as clearing and maintaining previus values*/
	for (i=current_depth;i<ntaxa-3;i++) {
		*(tuple_holder+i)=1;
		if (current_depth < (depth-1)) {
			passer=get_tuples_and_loop_tree(overflow,alignment_swapping,already_swapped,best_nodes,nodes,cur_b_tree,sequence,n_bases,num_best,current_rep,new_ones,num_to_swap,old_reps,current_depth+1,depth,tuple_holder,ntaxa,values,n_gaps,passer);
			old_reps=passer->rep;
			already_swapped=passer->swap;
		}
		else {
			doit=num_tuples=0;
			for (j=0;j<ntaxa-3;j++) {
				if (*(tuple_holder+j)!=(*(already_swapped[alignment_swapping]+j))) doit=1;
				num_tuples+=(*(tuple_holder+j));
			}
			if (num_tuples!=depth) doit=0;
			if (doit==1) {
				tree_rep=(int *)malloc((ntaxa-3)*sizeof(int));
				assert((int)tree_rep);
				for (j=0;j<ntaxa-3;j++) tree_rep[j]=*(current_rep+j);
				passer=special_all_loop_tree(overflow,alignment_swapping,already_swapped,old_reps,num_to_swap,new_ones,tuple_holder,sequence,n_bases,tree_rep,0,ntaxa,nodes,num_best,cur_b_tree,best_nodes,values,n_gaps,passer);
				old_reps=passer->rep;
				already_swapped=passer->swap;
				free(tree_rep);
			}
		}
		*(tuple_holder+i)=0;
	}
	passer->rep=old_reps;
	passer->swap=already_swapped;
	return passer;
}

realloc_holder *special_all_loop_tree(overflow,alignment_swapping,already_swapped,old_reps,num_to_swap,new_ones,tuple_holder,sequence,n_bases,tree_rep,element,ntax,nodes,l,cur_b_tree,best_nodes,values,n_gaps,passer)
int element,ntax,*tree_rep,n_gaps;
int **best_nodes,**nodes;
int **old_reps,alignment_swapping;
int *l,*cur_b_tree,num_to_swap,*new_ones;
parameters *values;
char *tuple_holder,**already_swapped;
int **sequence,n_bases,*overflow;
realloc_holder *passer;
{
	int i,j,cur_tree,doit;

	for (i=((2*(element+4))-5);i>=1;--i) {
		if (*(tuple_holder+element)==1) tree_rep[element]=i;
		else i=1;
		if (element<(ntax-4)) {
			passer=special_all_loop_tree(overflow,alignment_swapping,already_swapped,old_reps,num_to_swap,new_ones,tuple_holder,sequence,n_bases,tree_rep,element+1,ntax,nodes,l,cur_b_tree,best_nodes,values,n_gaps,passer);
			old_reps=passer->rep;
			already_swapped=passer->swap;
		}
		else {
			doit=0;
			for (j=0;j<ntax-3;j++) if (tree_rep[j]!=old_reps[alignment_swapping][j]) doit=1;
			if (doit==0) goto all_end_of_all_loop_special_tree;
			make_nodes(tree_rep,nodes,ntax,n_bases,sequence,ntax,n_gaps,values,&cur_tree);
			/*      cur_tree=return_value;*/
			if (cur_tree < (*cur_b_tree)) {
				values->length_best_clad=(*cur_b_tree) = cur_tree;
				*l=1;
				/*update the ones to swap*/
				if ((num_to_swap+(*new_ones))<values->keep_trees) *new_ones+=1;
				else if (values->max_out_trees) {
					if (((old_reps=(int **)realloc(old_reps,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((already_swapped=(char **)realloc(already_swapped,(values->keep_trees+1)*sizeof(char *)))!=NULL) && ((values->best_rep=(int **)realloc(values->best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL)) {
						++(*new_ones);
						++values->keep_trees;
						old_reps[num_to_swap+(*new_ones)-1]=(int *)malloc((ntax-3)*sizeof(int));
						assert((int)old_reps[num_to_swap+(*new_ones)-1]);
						values->best_rep[num_to_swap+(*new_ones)-1]=(int *)malloc((values->all_done-3)*sizeof(int));
						assert((int)values->best_rep[num_to_swap+(*new_ones)-1]);
						already_swapped[num_to_swap+(*new_ones)-1]=(char *)malloc((ntax-3)*sizeof(char));
						assert((int)already_swapped[num_to_swap+(*new_ones)-1]);
					}
					else if (*overflow==0) {
						if (values->rep_error) fprintf(stderr,"1-Overflow in cladogram swapping -- but that's OK.\n");
						*overflow=1;
					}
				}
				else if (*overflow==0) {
					if (values->rep_error) fprintf(stderr,"2-Overflow in cladogram swapping -- but that's OK.\n");
					*overflow=1;
				}
				for (j=0;j<ntax-3;j++) {
					*(old_reps[num_to_swap+(*new_ones)-1]+j)=values->best_rep[0][j]=tree_rep[j];
					*(already_swapped[num_to_swap+(*new_ones)-1]+j)=*(tuple_holder+j);
				}
				values->number_best_clad=1;
			}
			else if ((cur_tree == (*cur_b_tree)) && (!already_clad(values,tree_rep,ntax))) {
				if ((num_to_swap+(*new_ones))<values->keep_trees) {
					++values->number_best_clad;
					*new_ones+=1;
				}
				else if (values->max_out_trees) {
					if (((old_reps=(int **)realloc(old_reps,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((already_swapped=(char **)realloc(already_swapped,(values->keep_trees+1)*sizeof(char *)))!=NULL) && ((values->best_rep=(int **)realloc(values->best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL)) {
						++(*new_ones);
						++values->number_best_clad;
						++values->keep_trees;
						old_reps[num_to_swap+(*new_ones)-1]=(int *)malloc((ntax-3)*sizeof(int));
						assert((int)old_reps[num_to_swap+(*new_ones)-1]);
						values->best_rep[num_to_swap+(*new_ones)-1]=(int *)malloc((values->all_done-3)*sizeof(int));
						assert((int)values->best_rep[num_to_swap+(*new_ones)-1]);
						already_swapped[num_to_swap+(*new_ones)-1]=(char *)malloc((ntax-3)*sizeof(char));
						assert((int)already_swapped[num_to_swap+(*new_ones)-1]);
					}
					else if (*overflow==0) {
						if (values->rep_error) fprintf(stderr,"4-Overflow in cladogram swapping -- but that's OK.\n");
						*overflow=1;
					}
				}
				else if (*overflow==0) {
					if (values->rep_error) fprintf(stderr,"5-Overflow in cladogram swapping -- but that's OK.\n");
					*overflow=1;
				}
				for (j=0;j<ntax-3;j++) {
					already_swapped[num_to_swap+(*new_ones)-1][j]=tuple_holder[j];
					values->best_rep[values->number_best_clad-1][j]=old_reps[num_to_swap+(*new_ones)-1][j]=tree_rep[j];
				}
				if ((*l)==values->keep_trees) {
					if (*overflow==0) {
						if (values->rep_error) fprintf(stderr,"3-Overflow in cladogram swapping -- but that's OK.\n");
						*overflow=1;
					}
					goto all_end_of_all_loop_special_tree;
				}
				++(*l);
			}
		}
all_end_of_all_loop_special_tree:
		continue;
	}
	passer->rep=old_reps;
	passer->swap=already_swapped;
	return passer;
}




/*this is now Sankoff only*/
int diagnose_tree(taxa,nodes,nbases,ntaxa,nent,n_gaps,values,gap_counts)
int **taxa;
int **nodes,gap_counts;
int ntaxa,nbases,nent,n_gaps;
parameters *values;
{
	int tree_length,transition_counts,transversion_counts,change_counts,m;
	int n,d1,d2,i,j;
	int temp,k,l;
	int indel_temp,indel_best;
	int numanc,numd1,numd2;
	int  anc_pos[4],d1_pos[4],d2_pos[4];
	int cheapest_cost,*ambi_holder;
	char u_counter;

	/*
add here for matrix of costs
*/
	if ((values->delta) && (values->matquick)) {
		/*optimize root from outgroup (first taxon--ok if outgroup/symmetrical matrix)*/
		for (j=0;j<(nbases);j++) if (!unique(taxa[ntaxa][j])) {
			if ((taxa[ntaxa][j])&(taxa[0][j])) taxa[ntaxa][j]=(taxa[0][j])&(taxa[ntaxa][j]);
			else taxa[ntaxa][j]=(taxa[ntaxa][j])|(taxa[0][j]);
		}
		for (n=ntaxa+1;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;
		values->tree_made[ntaxa]=values->tree_made[nodes[0][1]]=1;

up_pass:
		for (n=0;n<=(nent-2);n++){
			if (values->tree_made[ntaxa+n]==1) {
				/*if descendant not unique intersection/union */
				for (i=0;i<2;i++) {
					d1=nodes[n][i];
					/*if (d1>ntaxa) {*/
					values->tree_made[d1]=1;
					for (j=0;j<(nbases);j++) if (!unique(taxa[d1][j])) {
						if ((taxa[ntaxa+n][j])&(taxa[d1][j])) taxa[d1][j]=(taxa[ntaxa+n][j])&(taxa[d1][j]);
						else if (d1>ntaxa) taxa[d1][j]=(taxa[ntaxa+n][j])|(taxa[d1][j]);
					}
					/*}*/
				}
			}
		}
		for (n=ntaxa;n<((ntaxa+nent)-1);++n) if (values->tree_made[n]==0) goto up_pass;

		/*here add different down pass to optimize ala Sankov and Cedergren*/
		temp=(gap_counts*values->gap_cost);
		for (n=ntaxa;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;

down_pass_delta:
		for (n=(nent-2);n>=0;--n){
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->tree_made[d1]==1) && (values->tree_made[d2]==1) && (values->tree_made[ntaxa+n]==0)) {
				values->tree_made[ntaxa+n]=1;
				for (m=0;m<nbases;m++) {
					if (taxa[d1][m]!=taxa[d2][m]) {
						if ((unique(taxa[ntaxa+n][m]))&&(unique(taxa[d1][m]))&&(unique(taxa[d2][m]))) temp+=values->delta[here_log2(taxa[ntaxa+n][m])][here_log2(taxa[d1][m])]+values->delta[here_log2(taxa[ntaxa+n][m])][here_log2(taxa[d2][m])];
						else {
							indel_temp=indel_best=HUGE_COST;
							/*loop through optimizations*/
							numanc=how_many(taxa[ntaxa+n][m],anc_pos);
							numd1=how_many(taxa[d1][m],d1_pos);
							numd2=how_many(taxa[d2][m],d2_pos);
							for (i=0;i<numanc;i++) {
								for (j=0;j<numd1;j++) {
									for (k=0;k<numd2;k++) {
										/*keep best and assign to ancestor*/
										indel_temp=(values->delta[here_log2(anc_pos[i])][here_log2(d1_pos[j])]+values->delta[here_log2(anc_pos[i])][here_log2(d2_pos[k])]);
										if (indel_temp<indel_best) {
											indel_best=indel_temp;
											taxa[ntaxa+n][m]=anc_pos[i];
											taxa[d1][m]=d1_pos[j];
											taxa[d2][m]=d2_pos[k];
										}
									}
								}
							}
							/*temp gets incremented with delta*/
							temp+=indel_best;
							/*
			      taxa[ntaxa+n][m]=anc_holder;
			      taxa[d1][m]=d1_holder;
			      taxa[d2][m]=d2_holder;
			      */
						}/*else not unique*/
					}
				}
			}
		}
		for (n=ntaxa;n<((ntaxa+nent)-1);++n) if (values->tree_made[n]==0) goto down_pass_delta;
		/*for (i=0;i<nent-1;i++) for (k=0;k<2;k++) for (j=nbases;j<(nbases+n_gaps);j++) if (taxa[ntaxa+i][j]!=taxa[nodes[i][k]][j]) temp+=values->gap_cost;*/
		tree_length=temp;
	}/*end matrix quick*/

	if ((values->delta) && (!values->matquick)) {
		/*optimize root from outgroup (first taxon--ok if outgroup/symmetrical matrix)*/
		for (j=0;j<(nbases);j++) if (!unique(taxa[ntaxa][j])) {
			if ((taxa[ntaxa][j])&(taxa[0][j])) taxa[ntaxa][j]=(taxa[0][j])&(taxa[ntaxa][j]);
			else taxa[ntaxa][j]=(taxa[ntaxa][j])|(taxa[0][j]);
		}
		for (n=ntaxa+1;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;
		values->tree_made[ntaxa]=values->tree_made[nodes[0][1]]=1;
up_pass_exhaustive:
		for (n=0;n<=(nent-2);n++){
			if (values->tree_made[ntaxa+n]==1) {
				/*if descendant not unique intersection/union */
				for (i=0;i<2;i++) {
					d1=nodes[n][i];
					/*if (d1>ntaxa) {*/
					values->tree_made[d1]=1;
					for (j=0;j<(nbases);j++) if (!unique(taxa[d1][j])) {
						if ((taxa[ntaxa+n][j])&(taxa[d1][j])) taxa[d1][j]=(taxa[ntaxa+n][j])&(taxa[d1][j]);
						else if (d1>ntaxa) taxa[d1][j]=(taxa[ntaxa+n][j])|(taxa[d1][j]);
					}
					/*}*/
				}
			}
		}
		for (n=ntaxa;n<((ntaxa+nent)-1);++n) if (values->tree_made[n]==0) goto up_pass_exhaustive;

		temp=(gap_counts*values->gap_cost);
		/*temp=0;*/
		ambi_holder=(int *)malloc((ntaxa+nent-1)*sizeof(int));
		assert((int)ambi_holder);

		/*optimize unique*/
		for (i=0;i<nent-1;i++) {
			for (k=0;k<2;k++) {
				for (j=0;j<nbases;j++) {
					u_counter=1;
					for (l=0;l<ntaxa+nent-1;l++) u_counter*=unique(taxa[l][j]);
					if (u_counter) {
						if (taxa[ntaxa+i][j]!=taxa[nodes[i][k]][j]) {
							temp+=values->delta[here_log2(taxa[ntaxa+i][j])][here_log2(taxa[nodes[i][k]][j])];
						}
					}
				}
			}
		}
		/*optimize non-unique*/
		for (j=0;j<nbases;j++) {
			u_counter=1;
			for (l=0;l<ntaxa+nent-1;l++) u_counter*=unique(taxa[l][j]);
			if (!u_counter) {
				cheapest_cost=HUGE_COST;
				cheap_loop(0,&cheapest_cost,ambi_holder,ntaxa,nent,nodes,taxa,j,values);
				temp+=cheapest_cost;
			}
		}
		/*add in gap characters*/
		/*
     for (i=0;i<nent-1;i++) for (k=0;k<2;k++) for (j=nbases;j<(nbases+n_gaps);j++) if (taxa[ntaxa+i][j]!=taxa[nodes[i][k]][j]) temp+=values->gap_cost;
     */
		/*
     printf("%d down steps %d up steps\n",tree_length,temp);
     */
		free(ambi_holder);
		tree_length=temp;
	}/*end matrix exhaustive*/

	return(tree_length);
}

int diagnose_tree_and_get_weights(taxa,nodes,nbases,ntaxa,nent,n_gaps,values)
int **taxa;
int **nodes;
int ntaxa,nbases,nent,n_gaps;
parameters *values;
{
	int tree_length,transition_counts,transversion_counts,gap_counts,change_counts,m;
	int n,d1,d2,i,j;
	int temp,k,l;
	int indel_temp,indel_best;
	int numanc,numd1,numd2;
	int  anc_pos[4],d1_pos[4],d2_pos[4];
	int cheapest_cost,*ambi_holder;
	int temp1,temp2,a,c,g,t;
	int *temp3, *temp4,*temp5;
	int new_transitions,old_transitions;
	int *taxon_holder,n_what_holder,inf,nstates;
	char u_counter,A,C,G,T,*sym;

	tree_length=transversion_counts=transition_counts=gap_counts=change_counts=0;
	/*need to insure that we start with nodes which connect to terminal taxa and then nodes
     which have been defined - otherwise we'll be using the optimizations of
     previous topology's nodes*/



	/*
     recode transition characters for appropriate Farris/Goloboff weights
     0-old_transiitons-1 are A/G
     old_transitions-new_transitions1 are C/T
     */
	/*should transversion cost be upped by transition?*/
	if (values->ttr) {
		old_transitions=nbases-values->n_transv;
		new_transitions=2*old_transitions;
		for (i=0;i<ntaxa;i++) {
			taxa[i]=(int *)realloc(taxa[i],(new_transitions+values->n_transv+n_gaps)*sizeof(int));
			assert((int)taxa[i]);
			taxon_holder=(int *)malloc((nbases+n_gaps)*sizeof(int));
			assert((int)taxon_holder);
			for (j=0;j<nbases+n_gaps;j++) taxon_holder[j]=taxa[i][j];
			for (j=0;j<old_transitions;j++) taxa[i][j]=taxon_holder[j];
			for (j=old_transitions;j<n_gaps+nbases+old_transitions;j++) taxa[i][j]=taxon_holder[j-old_transitions];
			free(taxon_holder);
			for (j=0;j<old_transitions;j++) if ((taxa[i][j]!=2) && (taxa[i][j]!=8)) taxa[i][j]=missing;
			for (j=old_transitions;j<new_transitions;j++) if ((taxa[i][j]!=4) && (taxa[i][j]!=16)) taxa[i][j]=missing;
		}
		nbases+=old_transitions;
		/*recode to get rid of uninformative*/
		sym=(char *)malloc((nbases-values->n_transv)*sizeof(char));
		assert((int)sym);
		for (j=0;j<nbases-values->n_transv;j++) {
			nstates=A=C=G=T=0;
			for (i=0;i<ntaxa;i++){
				if (taxa[i][j]==2) ++A;
				else if (taxa[i][j]==4) ++C;
				else if (taxa[i][j]==8) ++G;
				else if (taxa[i][j]==16) ++T;
			}
			if (A>0) ++nstates;
			if (C>0) ++nstates;
			if (G>0) ++nstates;
			if (T>0) ++nstates;
			if (nstates<2) sym[j]=1;
			else if ((A==1)||(C==1)||(G==1)||(T==1)) sym[j]=1;
			else sym[j]=0;
		}
		inf=0;
		for (j=0;j<nbases-values->n_transv;j++) if (!sym[j]) ++inf;
		for (i=0;i<ntaxa;i++) {
			taxon_holder=(int *)malloc((nbases+n_gaps)*sizeof(int));
			assert((int)taxon_holder);
			for (j=0;j<nbases+n_gaps;j++) taxon_holder[j]=taxa[i][j];
			taxa[i]=(int *)realloc(taxa[i],(inf+values->n_transv+n_gaps)*sizeof(int));
			assert((int)taxa[i]);
			u_counter=0;
			for (j=0;j<nbases-values->n_transv;j++) {
				if (!sym[j]) {
					taxa[i][u_counter]=taxon_holder[j];
					++u_counter;
				}
			}
			for (j=nbases-values->n_transv;j<nbases+n_gaps;j++) {
				taxa[i][u_counter]=taxon_holder[j];
				++u_counter;
			}
			free(taxon_holder);
		}
		free(sym);
		nbases=inf+values->n_transv;
	}

	/*
fprintf(stderr,"recoded nb %d nt %d ot %d nv %d ng %d\n",nbases,new_transitions,old_transitions,values->n_transv,n_gaps);
for (i=0;i<ntaxa;i++) {
     for (j=0;j<nbases+n_gaps;j++) fprintf(stderr,"%3d:%2d ",j,taxa[i][j]);
     fprintf(stderr,"\n");
     }
*/

	temp3=NULL;
	temp4=NULL;
	temp5=NULL;
	if (n_gaps==0) {
		fprintf(stderr,"There are no informative gaps upon which to base weights.\n");
		values->saw_goloboff=0;
		values->saw_farris=0;
	}
	if (values->saw_goloboff > 0) {
		temp5=(int *)malloc(n_gaps*sizeof(int));
		assert((int)temp5);
		for (i=0;i<n_gaps;i++) temp5[i]=0;
		if (!values->ttr) {
			temp3=(int *)malloc(nbases*sizeof(int));
			assert((int)temp3);
			for (i=0;i<nbases;i++) temp3[i]=0;
		}
		else {
			temp3=(int *)malloc((nbases-values->n_transv)*sizeof(int));
			assert((int)temp3);
			for (i=0;i<(nbases-values->n_transv);i++) temp3[i]=0;
			temp4=(int *)malloc(values->n_transv*sizeof(int));
			assert((int)temp4);
			for (i=0;i<values->n_transv;i++ ) temp4[i]=0;
		}
	}

	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;

make_the_node:
	for (n=(nent-2);n>=0;--n){
		d1=nodes[n][0];
		d2=nodes[n][1];
		if ((values->tree_made[d1]==1) && (values->tree_made[d2]==1) && (values->tree_made[ntaxa+n]==0)) {
			values->tree_made[ntaxa+n]=1;
			if (!values->ttr) {
				for (m=0;m<nbases;m++) {
					if (taxa[d1][m]&taxa[d2][m]) taxa[ntaxa+n][m]=(taxa[d1][m]&taxa[d2][m]);
					else {
						taxa[ntaxa+n][m]=(taxa[d1][m]|taxa[d2][m]);
						/*tree_length+=values->change_cost;*/
						++change_counts;
						if (temp3) ++temp3[m];
					}
				}
			}
			else {
				for (m=0;m<nbases-values->n_transv;m++) {
					if (taxa[d1][m]&taxa[d2][m]) taxa[ntaxa+n][m]=(taxa[d1][m]&taxa[d2][m]);
					else {
						taxa[ntaxa+n][m]=(taxa[d1][m]|taxa[d2][m]);
						/*tree_length+=values->transition;
				   ++tree_length;*/
						++transition_counts;
						if (temp3) ++temp3[m];
					}
				}
				for (m=nbases-values->n_transv;m<nbases;m++) {
					if (taxa[d1][m]&taxa[d2][m]) taxa[ntaxa+n][m]=(taxa[d1][m]&taxa[d2][m]);
					else {
						taxa[ntaxa+n][m]=(taxa[d1][m]|taxa[d2][m]);
						/*tree_length+=values->transversion;
				   ++tree_length;*/
						++transversion_counts;
						if (temp4) ++temp4[m-(nbases-values->n_transv)];
					}
				}
			}

			for (m=nbases;m<(nbases+n_gaps);m++) {
				if (taxa[d1][m]&taxa[d2][m]) taxa[ntaxa+n][m]=(taxa[d1][m]&taxa[d2][m]);
				else {
					taxa[ntaxa+n][m]=(taxa[d1][m]|taxa[d2][m]);
					/*tree_length+=values->gap_cost;*/
					++gap_counts;
					if (temp5) ++temp5[m-nbases];
				}
			}

		}
	}
	for (n=ntaxa;n<((ntaxa+nent)-1);++n) if (values->tree_made[n]==0) goto make_the_node;
	tree_length=(change_counts*values->change_cost)+(gap_counts*values->gap_cost)+(transition_counts*values->transition)+(transversion_counts*values->transversion);

	if (!values->ttr) {
		temp1=temp2=0;
		for (i=0;i<nbases;i++) {
			n_what_holder=a=c=g=t=A=C=G=T=0;
			for (j=0;j<ntaxa;j++) {
				switch(taxa[j][i]) {
				case 2 :
					++a;
					break;
				case 4 :
					++c;
					break;
				case 8 :
					++g;
					break;
				case 16 :
					++t;
					break;
				default :
					break;
				}
			}
			temp1+=(a+c+g+t);
			if (a>0) A=1;
			if (c>0) C=1;
			if (g>0) G=1;
			if (t>0) T=1;
			temp2+=(A+C+G+T);
			if (values->saw_goloboff > 0) values->temp_changes_g+=(values->saw_goloboff*values->weight_range/(values->saw_goloboff+temp3[i]-(A+C+G+T-1)));
		}
		temp1-=nbases;
		temp2-=nbases;

		/* this is ci*ri*range*/
		values->temp_changes=((temp1-change_counts)*temp2*values->weight_range)/(change_counts*(temp1-temp2));
		if (values->temp_changes < 0) values->temp_changes=0;
		if (values->saw_goloboff > 0) values->temp_changes_g/=nbases;
	}
	else {
		/*transitions*/
		temp1=temp2=0;
		for (i=0;i<nbases-values->n_transv;i++) {
			n_what_holder=a=c=g=t=A=C=G=T=0;
			for (j=0;j<ntaxa;j++) {
				switch(taxa[j][i]) {
				case 2 :
					++a;
					break;
				case 4 :
					++c;
					break;
				case 8 :
					++g;
					break;
				case 16 :
					++t;
					break;
				default :
					break;
				}

			}
			temp1+=(a+c+g+t);
			if (a>0) A=1;
			if (c>0) C=1;
			if (g>0) G=1;
			if (t>0) T=1;
			temp2+=(A+C+G+T);
			if (values->saw_goloboff > 0) values->temp_ti_g+=(values->saw_goloboff*values->weight_range/(values->saw_goloboff+temp3[i]-(A+C+G+T-1)));
		}
		temp1-=nbases-values->n_transv;
		temp2-=nbases-values->n_transv;

		values->temp_ti=((temp1-transition_counts)*temp2*values->weight_range)/(transition_counts*(temp1-temp2));
		if (values->temp_ti < 0) values->temp_ti=0;
		if (values->saw_goloboff > 0) values->temp_ti_g/=(nbases-values->n_transv);

		/*transversions*/
		/*something is wrong with max min counts here*/
		temp1=temp2=0;
		for (i=nbases-values->n_transv;i<nbases;i++) {
			n_what_holder=a=c=g=t=A=C=G=T=0;
			for (j=0;j<ntaxa;j++) {
				switch(taxa[j][i]) {
				case 2 :
					++a;
					break;
				case 4 :
					++c;
					break;
				case 8 :
					++g;
					break;
				case 16 :
					++t;
					break;
				default :
					break;
				}
			}
			temp1+=(a+c+g+t);
			if (a>0) A=1;
			if (c>0) C=1;
			if (g>0) G=1;
			if (t>0) T=1;
			temp2+=(A+C+G+T);
			if (values->saw_goloboff > 0) values->temp_tv_g+=(values->saw_goloboff*values->weight_range/(values->saw_goloboff+temp4[i-(nbases-values->n_transv)]-(A+C+G+T-1)));
		}
		temp1-=values->n_transv;
		temp2-=values->n_transv;

		values->temp_tv=((temp1-transversion_counts)*temp2*values->weight_range)/(transversion_counts*(temp1-temp2));
		if (values->temp_tv < 0) values->temp_tv=0;
		if (values->saw_goloboff > 0) values->temp_tv_g/=values->n_transv;
	}

	/*getting weights from gaps*/
	temp1=temp2=0;
	for (i=nbases;i<nbases+n_gaps;i++) {
		n_what_holder=a=c=A=C=0;
		for (j=0;j<ntaxa;j++) {
			switch(taxa[j][i]) {
			case 2 :
				++a;
				break;
			case 4 :
				++c;
				break;
			default :
				break;
			}
		}
		temp1+=(a+c);
		if (a>0) A=1;
		if (c>0) C=1;
		temp2+=(A+C);
		if (values->saw_goloboff > 0) values->temp_gaps_g+=(values->saw_goloboff*values->weight_range/(values->saw_goloboff+temp5[i-nbases]-(A+C-1)));
	}
	temp1-=n_gaps;
	temp2-=n_gaps;

	values->temp_gaps=((temp1-gap_counts)*temp2*values->weight_range)/(gap_counts*(temp1-temp2));
	if (values->temp_gaps < 0) values->temp_gaps=0;
	if (values->saw_goloboff > 0) values->temp_gaps_g/=n_gaps;

	/*
add here for matrix of costs
*/
	if (temp3) free(temp3);
	if (temp4) free(temp4);
	if (temp5) free(temp5);
	return(tree_length);

}

int here_log2(x)
int x;
{
	if (x==2) return 0;
	else if (x==4) return 1;
	else if (x==8) return 2;
	else if (x==16) return 3;

	return 4;
}



int how_many(x,there)
int x,there[4];
{
	int number_of_resolutions;

	there[0]=there[1]=there[2]=there[3]=0;

	switch (x) {
	case 2    :
		number_of_resolutions=1;
		there[0]= 2;
		break;

	case 4    :
		number_of_resolutions=1;
		there[0]= 4;
		break;

	case 8    :
		number_of_resolutions=1;
		there[0]= 8;
		break;

	case 16 :
		number_of_resolutions=1;
		there[0]=16;
		break;

	case 6    :
		number_of_resolutions=2;
		there[0]= 2;
		there[1]= 4;
		break;

	case 10 :
		number_of_resolutions=2;
		there[0]= 2;
		there[1]= 8;
		break;

	case 18 :
		number_of_resolutions=2;
		there[0]= 2;
		there[1]=16;
		break;

	case 12 :
		number_of_resolutions=2;
		there[0]= 4;
		there[1]= 8;
		break;

	case 20 :
		number_of_resolutions=2;
		there[0]= 4;
		there[1]=16;
		break;

	case 24 :
		number_of_resolutions=2;
		there[0]= 8;
		there[1]=16;
		break;

	case 28 :
		number_of_resolutions=3;
		there[0]= 4;
		there[1]= 8;
		there[2]=16;
		break;

	case 26 :
		number_of_resolutions=3;
		there[0]= 2;
		there[1]= 8;
		there[2]=16;
		break;

	case 22 :
		number_of_resolutions=3;
		there[0]= 2;
		there[1]= 4;
		there[2]=16;
		break;

	case 14 :
		number_of_resolutions=3;
		there[0]= 2;
		there[1]= 4;
		there[2]= 8;
		break;

	case 30 :
		number_of_resolutions=4;
		there[0]= 2;
		there[1]= 4;
		there[2]= 8;
		there[3]=16;
		break;

	case missing :
		number_of_resolutions=4;
		there[0]= 2;
		there[1]= 4;
		there[2]= 8;
		there[3]=16;
		break;

	default :
		fprintf(stderr,"Unrecognized ambiguity--error.\n");
		exit(1);
	}

	return number_of_resolutions;
}


char unique(x)
int x;
{
	if (x==2) return 1;
	else if (x==4) return 1;
	else if (x==8) return 1;
	else if (x==16) return 1;
	/*
else if (x==missing) return 1;
*/

	return 0;
}

char already_clad(values,current_rep,ntax)
parameters *values;
int *current_rep,ntax;
{
	int i,j,temp,holder;

	holder=0;
	for (i=0;i<values->number_best_clad;i++) {
		temp=0;
		for (j=0;j<ntax-3;j++) if (values->best_rep[i][j]!=current_rep[j]) temp=1;
		holder += temp;
	}
	if (holder==values->number_best_clad) return 0;
	else return 1;
}

void cheap_loop(to_optimize,cost,ambi_holder,ntaxa,nent,nodes,taxa,j,values)
int *cost,*ambi_holder,ntaxa,nent,**nodes,**taxa,j;
parameters *values;
{
	int k,cur_cost,i,l,number;
	int possible[4];


	/*not looping through internal stuff*/
	number=how_many(taxa[to_optimize][j],possible);
	for (k=0;k<number;k++) {
		ambi_holder[to_optimize]=possible[k];
		if (to_optimize<(ntaxa+nent-2)) cheap_loop(to_optimize+1,cost,ambi_holder,ntaxa,nent,nodes,taxa,j,values);
		else {
			cur_cost=0;
			for(i=0;i<(nent-1);i++) {
				for (l=0;l<2;l++) {
					if (ambi_holder[ntaxa+i]!=ambi_holder[nodes[i][l]])
						cur_cost+=values->delta[here_log2(ambi_holder[ntaxa+i])][here_log2(ambi_holder[nodes[i][l]])];
				}
			}
			if (cur_cost < (*cost)) *cost=cur_cost;
		}/*else*/
	}/*k*/

}

int get_quick(tree_rep,nodes,ntax,nbases,sequence,n_gaps,values)
int *tree_rep, **nodes, ntax,nbases,**sequence,n_gaps;
parameters *values;
{
	int i,j,value,*temp_rep, cur_val;

	temp_rep=(int *)malloc((ntax-3)*sizeof(int));
	assert((int)temp_rep);

	/*build tree*/
	value=HUGE_COST;
	for (i=0;i<ntax-3;i++) {
		value=HUGE_COST;
		for (j=1;j<=((2*i)+3);j++) {
			temp_rep[i]=j;
			make_nodes(temp_rep,nodes,ntax,nbases,sequence,i+4,n_gaps,values,&cur_val);
			if (cur_val<value) {
				tree_rep[i]=j;
				value=cur_val;
			}
		}
		temp_rep[i]=tree_rep[i];
		if ((values->clade_swap_while_add) && (i>0) && ((i+4)<ntax))     {
			for (j=0;j<(ntax-3);j++) tree_rep[j]=temp_rep[j];
			clade_addition_swap_1(tree_rep,&value,i,nodes,ntax,nbases,sequence,n_gaps,values);
		}
		if ((values->tree_swap) && ((i+4)==ntax))    {
			for (j=0;j<(ntax-3);j++) tree_rep[j]=temp_rep[j];
			clade_addition_swap_1(tree_rep,&value,i,nodes,ntax,nbases,sequence,n_gaps,values);
		}
	}

	if (values->get_weights) {
		for (i=0;i<ntax-3;i++) values->best_rep[0][i]=temp_rep[i];
		values->length_best_clad=value+aut;
		values->number_best_clad=1;
	}

	free(temp_rep);

	return value;
}
void clade_addition_swap_1(tree_rep,cur_best,current_adding,nodes,ntax,nbases,sequence,n_gaps,values)
int *cur_best,current_adding;
int *tree_rep, **nodes, ntax,nbases,**sequence,n_gaps;
parameters *values;
{
	int i,j,k,cur_val;
	int *temp_rep,*best_rep;
	char found_better;

	temp_rep=(int *)malloc((current_adding+1)*sizeof(int));
	assert((int)temp_rep);
	best_rep=(int *)malloc((current_adding+1)*sizeof(int));
	assert((int)best_rep);

top:
	;
	found_better=0;

	for (i=0;i<=current_adding;i++) temp_rep[i]=tree_rep[i];
	for (i=0;i<current_adding;i++) {
		for (j=1;j<=((2*i)+3);j++) {
			temp_rep[i]=j;
			make_nodes(temp_rep,nodes,ntax,nbases,sequence,current_adding+4,n_gaps,values,&cur_val);
			/*    cur_val=return_value;*/
			if (cur_val<(*cur_best)) {
				for (k=0;k<=current_adding;k++) best_rep[k]=temp_rep[k];
				(*cur_best)=cur_val;
				found_better=1;
			}
		}
		temp_rep[i]=tree_rep[i];
	}

	for (i=0;i<=current_adding;i++) tree_rep[i]=best_rep[i];

	if (found_better) goto top;

	free(best_rep);
	free(temp_rep);
}


int get_quick2(nodes,ntax,nbases,sequence,n_gaps,values)
int ntax,nbases,n_gaps,**sequence;
int **nodes;
parameters *values;
{
	int i,j,k,z,y,counter,value,*temp_rep;
	int cur_val,**new_rep,**old_rep,l,temp_tree_holder,old_keep_trees_holder;
	char **from,overflow;

	temp_rep=(int *)malloc((ntax-3)*sizeof(int));
	assert((int)temp_rep);
	new_rep=(int **)malloc(values->keep_trees*sizeof(int *));
	assert((int)new_rep);
	old_rep=(int **)malloc(values->keep_trees*sizeof(int *));
	assert((int)old_rep);
	for (i=0;i<values->keep_trees;i++) {
		new_rep[i]=(int *)malloc((ntax-3)*sizeof(int));
		assert((int)new_rep[i]);
		old_rep[i]=(int *)malloc((ntax-3)*sizeof(int));
		assert((int)old_rep[i]);
	}

	y=overflow=counter=0;
	for (i=0;i<ntax-3;i++) {
		value=HUGE_COST;
		for (z=0;z<=counter;z++) {
			for (j=0;j<ntax-3;j++) temp_rep[j]=old_rep[z][j];
			for (j=1;j<=((2*i)+3);j++) {
				temp_rep[i]=j;
				make_nodes(temp_rep,nodes,ntax,nbases,sequence,i+4,n_gaps,values,&cur_val);
				/*
	       cur_val=return_value;
*/
				if (cur_val<value) {
					y=0;
					for (k=0;k<ntax-3;k++) *(new_rep[y]+k)=values->best_rep[0][k]=temp_rep[k];
					value=values->length_best_clad=cur_val;
					values->number_best_clad=1;
				}
				else if ((cur_val==value) && (already_clad(values,temp_rep,ntax))) {
					if (y==(values->keep_trees-1)) {
						if ((values->max_out_trees) && ((old_rep=(int **)realloc(old_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((from=(char **)realloc(from,(values->keep_trees+1)*sizeof(char *)))!=NULL) && ((values->best_rep=(int **)realloc(values->best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) &&((new_rep=(int
						**)realloc(new_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL)) {
							++y;
							++values->keep_trees;
							old_rep[y]=(int *)malloc((ntax-3)*sizeof(int));
							assert((int)old_rep[y]);
							new_rep[y]=(int *)malloc((ntax-3)*sizeof(int));
							assert((int)new_rep[y]);
							values->best_rep[y]=(int *)malloc((values->all_done-3)*sizeof(int));
							assert((int)values->best_rep[y]);
							for (k=0;k<(ntax-3);k++) new_rep[y][k]=values->best_rep[values->number_best_clad][k]=temp_rep[k];
							++values->number_best_clad;
						}
						else if (overflow==0) {
							if (values->rep_error) fprintf(stderr,"Overflow in cladogram construction -- but that's OK.\n");
							overflow=1;
							--y;
						}
					}
					else {
						++y;
						for (k=0;k<ntax-3;k++) *(new_rep[y]+k)=values->best_rep[values->number_best_clad][k]=temp_rep[k];
						++values->number_best_clad;
					}
				}
			}/*j*/
			for (k=0;k<=y;k++) {
				for (j=0;j<ntax-3;j++) old_rep[k][j]=new_rep[k][j];
			}
		}/*z*/
		counter=y;
		y=0;
		if ((values->clade_swap_while_add) && (i>0) && ((i+4)<ntax))     {
			/*put in holder here for increments*/
			temp_tree_holder=values->keep_trees;
			old_rep=clade_addition_swap_2(&counter,old_rep,&value,i,nodes,ntax,nbases,sequence,n_gaps,values,&overflow);
			/*can reallocate others here if keep trees has changed*/
			if (values->keep_trees>temp_tree_holder) {
				new_rep=(int **)realloc(new_rep,values->keep_trees*sizeof(int *));
				for (j=temp_tree_holder-1;j<values->keep_trees;j++) {
					new_rep[j]=(int *)malloc((ntax-3)*sizeof(int));
					assert((int)new_rep[j]);
				}
			}
		}
	}/*i*/
	/*insert better swap for trees here*/
	l=counter+1;
	if (values->tree_multi_swap>1) {
		old_keep_trees_holder=values->keep_trees;
		old_rep=deep_swap_tree(&value,sequence,nbases,counter+1,old_rep,&l,ntax,values->tree_multi_swap,values,n_gaps);
		if (values->keep_trees > old_keep_trees_holder) {
			new_rep=(int **)realloc(new_rep,values->keep_trees*sizeof(int *));
			for (j=old_keep_trees_holder-1;j<values->keep_trees;j++) {
				new_rep[j]=(int *)malloc((ntax-3)*sizeof(int));
				assert((int)new_rep[j]);
			}
		}
	}
	else if ((values->tree_swap!=0)|| (values->tree_multi_swap==1)) {
		say("Swapping")
		    GOMAC
		    from=(char **)malloc(values->keep_trees*sizeof(char *));
		assert((int)from);
		for (i=0;i<values->keep_trees;i++) {
			from[i]=(char *)malloc((ntax-3)*sizeof(char));
			assert((int)from[i]);
			for (j=0;j<(ntax-3);j++) *(from[i]+j)=0;
		}
		for (z=0;z<=counter;z++) *(from[z]+(ntax-4))=1;
		for (z=0;z<=counter;z++) {
			for (i=0;i<ntax-3;i++) temp_rep[i]=old_rep[z][i];
			for (i=0;i<ntax-3;i++) {
				if ((*(from[z]+i))==0) {
					for (j=1;j<=((2*i)+3);j++) {
						if (j!=old_rep[z][i]) {
							temp_rep[i]=j;
							make_nodes(temp_rep,nodes,ntax,nbases,sequence,ntax,n_gaps,values,&cur_val);
							/*         cur_val=return_value;*/
							if (cur_val<value) {
								values->length_best_clad=value=cur_val;
								values->number_best_clad=1;
								if (counter==(values->keep_trees-1)) {
									if ((values->max_out_trees) && ((old_rep=(int **)realloc(old_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((from=(char **)realloc(from,(values->keep_trees+1)*sizeof(char *)))!=NULL) &&
									    ((values->best_rep=(int **)realloc(values->best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) &&((new_rep=(int **)realloc(new_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL)) {
										++counter;
										++values->keep_trees;
										old_rep[counter]=(int *)malloc((ntax-3)*sizeof(int));
										assert((int)old_rep[counter]);
										new_rep[counter]=(int *)malloc((ntax-3)*sizeof(int));
										assert((int)new_rep[counter]);
										values->best_rep[counter]=(int *)malloc((values->all_done-3)*sizeof(int));
										assert((int)values->best_rep[counter]);
										from[counter]=(char *)malloc((ntax-3)*sizeof(char));
										assert((int)from[counter]);
										for (k=0;k<(ntax-3);k++) {
											old_rep[counter][k]=values->best_rep[0][k]=temp_rep[k];
											from[counter][k]=0;
										}
										from[counter][i]=1;
									}
									else if (overflow==0) {
										if (values->rep_error) fprintf(stderr,"Overflow in cladogram construction -- but that's OK.\n");
										overflow=1;
										--counter;
										for (k=0;k<ntax-3;k++) *(from[counter+1]+k)=0;
									}
								}
								else {
									++counter;
									from[counter][i]=1;
									for (k=0;k<(ntax-3);k++) old_rep[counter][k]=values->best_rep[0][k]=temp_rep[k];
								}
							}
							else if ((cur_val==value) && (!already_clad(values,temp_rep,ntax))) {
								if (counter==(values->keep_trees-1)) {
									if ((values->max_out_trees==1) && ((old_rep=(int **)realloc(old_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((from=(char **)realloc(from,(values->keep_trees+1)*sizeof(char *)))!=NULL) &&
									    ((values->best_rep=(int **)realloc(values->best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) &&((new_rep=(int **)realloc(new_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL)) {
										++counter;
										++values->keep_trees;
										old_rep[counter]=(int *)malloc((ntax-3)*sizeof(int));
										assert((int)old_rep[counter]);
										new_rep[counter]=(int *)malloc((ntax-3)*sizeof(int));
										assert((int)new_rep[counter]);
										values->best_rep[counter]=(int *)malloc((values->all_done-3)*sizeof(int));
										assert((int)values->best_rep[counter]);
										from[counter]=(char *)malloc((ntax-3)*sizeof(char));
										assert((int)from[counter]);
										for (k=0;k<(ntax-3);k++) {
											old_rep[counter][k]=values->best_rep[values->number_best_clad][k]=temp_rep[k];
											from[counter][k]=0;
										}
										from[counter][i]=1;
										++values->number_best_clad;
									}
									else if (overflow==0) {
										if (values->rep_error) fprintf(stderr,"Overflow in cladogram construction -- but that's OK.\n");
										overflow=1;
									}
								}
								else {
									++counter;
									from[counter][i]=1;
									for (k=0;k<(ntax-3);k++) old_rep[counter][k]=values->best_rep[values->number_best_clad][k]=temp_rep[k];
									++values->number_best_clad;
								}
							} /*else if*/
						}/*j!=*/
					}/*j*/
					temp_rep[i]=old_rep[z][i];
				}/*check i for swap*/
			}/*i*/
		}/*z*/
		for (i=0;i<values->keep_trees;i++) if (from[i]) free(from[i]);
		free(from);
	}/*if*/

	for (i=0;i<values->keep_trees;i++) {
		if (old_rep[i]) free(old_rep[i]);
		if (new_rep[i]) free(new_rep[i]);
	}

	free(new_rep);
	free(old_rep);
	free(temp_rep);
	return value;
}

int **clade_addition_swap_2(counter,old_rep,cur_best,current_adding,nodes,ntax,nbases,sequence,n_gaps,values,overflow)
int *cur_best,current_adding,*counter;
int **old_rep, **nodes, ntax,nbases,**sequence,n_gaps;
char *overflow;
parameters *values;
{
	int i,j,k,m,cur_val,holder;
	int *temp_rep,**best_rep,sub_taxa;
	char found_better,**from;

	sub_taxa=current_adding+4;
	temp_rep=(int *)malloc((sub_taxa-3)*sizeof(int));
	assert((int)temp_rep);

	best_rep=(int **)malloc(values->keep_trees*sizeof(int *));
	assert((int)best_rep);
	for (i=0;i<values->keep_trees;i++) {
		best_rep[i]=(int *)malloc((sub_taxa-3)*sizeof(int));
		assert((int)best_rep[i]);
	}

	from=(char **)malloc(values->keep_trees *sizeof(char *));
	assert((int)from);
	for (i=0;i<values->keep_trees;i++) {
		from[i]=(char *)malloc((sub_taxa-3)*sizeof(char));
		assert((int)from[i]);
		for (j=0;j<(sub_taxa-3);j++) from[i][j]=0;
	}

top:
	;
	found_better=0;

	holder=(*counter);
	for (j=0;j<=(*counter);j++) {
		for (i=0;i<(sub_taxa-3);i++) best_rep[j][i]=old_rep[j][i];
	}


	for (m=0;m<=holder;m++) {
		for (i=0;i<(sub_taxa-3);i++) temp_rep[i]=old_rep[m][i];
		for (i=0;i<(sub_taxa-4);i++) {
			if (!from[m][i]) {
				for (j=1;j<=((2*i)+3);j++) {
					temp_rep[i]=j;
					make_nodes(temp_rep,nodes,ntax,nbases,sequence,sub_taxa,n_gaps,values,&cur_val);
					/*        cur_val=return_value;*/
					if (cur_val<(*cur_best)) {
						found_better=1;
						/*
			 fprintf(stderr,"Found better\n");
			 */
						for (k=0;k<(sub_taxa-3);k++) best_rep[0][k]=values->best_rep[0][k]=temp_rep[k];
						values->length_best_clad=(*cur_best)=cur_val;
						(*counter)=0;
						from[0][i]=1;
						values->number_best_clad=1;
					}
					if ((cur_val==(*cur_best)) && (!already_clad(values,temp_rep,sub_taxa))) {
						if ((*counter)<(values->keep_trees-1)) {
							found_better=1;
							++(*counter);
							for (k=0;k<(sub_taxa-3);k++) best_rep[(*counter)][k]=values->best_rep[values->number_best_clad][k]=temp_rep[k];
							from[*counter][i]=1;/*problem here*/
							++values->number_best_clad;
						}
						else if (values->max_out_trees) {
							if (((old_rep=(int**)realloc(old_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((best_rep=(int**)realloc(best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) && ((values->best_rep=(int**)realloc(values->best_rep,(values->keep_trees+1)*sizeof(int *)))!=NULL) &&((from=(char **)realloc(from,(values->keep_trees+1)*sizeof(char*)))!=NULL)) {
								found_better=1;
								++(*counter);
								++values->keep_trees;
								old_rep[values->keep_trees-1]=(int *)malloc((ntax-3)*sizeof(int));
								assert((int)old_rep[values->keep_trees-1]);
								best_rep[values->keep_trees-1]=(int *)malloc((sub_taxa-3)*sizeof(int));
								assert((int)best_rep[values->keep_trees-1]);
								values->best_rep[values->keep_trees-1]=(int *)malloc((values->all_done-3)*sizeof(int));
								assert((int)values->best_rep[values->keep_trees-1]);
								from[values->keep_trees-1]=(char *)malloc((sub_taxa-3)*sizeof(int));
								assert((int)from[values->keep_trees-1]);
								for (k=0;k<(sub_taxa-3);k++) {
									best_rep[(*counter)][k]=values->best_rep[values->number_best_clad][k]=temp_rep[k];
									from[*counter][k]=0;
								}
								from[*counter][i]=1;/*problem here*/
								++values->number_best_clad;
							}
							else {
								if ((*overflow)==0) {
									if (values->rep_error) fprintf(stderr,"Overflow in cladogram addition swap -- but that's OK.\n");
									(*overflow)=1;
								}
							}
						}
						else {
							if ((*overflow)==0) {
								if (values->rep_error) fprintf(stderr,"Overflow in cladogram addition swap -- but that's OK.\n");
								(*overflow)=1;
							}
						}
					}
				}
				temp_rep[i]=old_rep[m][i];
			}
		}
	}

	for (j=0;j<=(*counter);j++) {
		for (i=0;i<(sub_taxa-3);i++) old_rep[j][i]=best_rep[j][i];
	}

	if (found_better) goto top;

	for (i=0;i<values->keep_trees;i++) {
		if (best_rep[i]) free(best_rep[i]);
		if (from[i]) free(from[i]);
	}
	free(from);
	free(best_rep);
	free(temp_rep);

	return(old_rep);
}

void pre_diagnose_tree(tree_rep,nodes,ntaxa,nbases,nent,n_gaps,values,sequence)
int **sequence;
int ntaxa,nbases,n_gaps,nent;
int *tree_rep;
int **nodes;
parameters *values;
{
	int i,d1,d2,all_nodes_done,j;
	int holder;
	register int ii;
	int ls1;

	if (!changes[ntaxa].name) {
		for (i=ntaxa;i<((nent+ntaxa)-1);i++) changes[i].changes=changes[i].gaps=changes[i].transitions=changes[i].transversions=0;
		for (i=ntaxa;i<(ntaxa+nent-1);i++) num_seqs[i]=pre_tree_made[i]=0;
		all_nodes_done=0;
		for (;;){
			for (i=(nent-2);i>=0;--i){
				d1=nodes[i][0];
				d2=nodes[i][1];
				if ((pre_tree_made[d1]) && (pre_tree_made[d2]) && (!pre_tree_made[ntaxa+i])) {
					pre_tree_made[ntaxa+i]=1;
					++all_nodes_done;
					/*make names here*/
					if (strcmp (changes[d1].name,changes[d2].name) > 0) changes[ntaxa+i].name=(char *)other_pair_names(changes[d2].name,changes[d1].name);
					else changes[ntaxa+i].name=(char *)other_pair_names(changes[d1].name,changes[d2].name);
					num_seqs[ntaxa+i]=num_seqs[d1]+num_seqs[d2];
				}
			}
			if (all_nodes_done>nent-2) break;
		}
	}

	else {
		/*try names for identity first
     is overhead too high?
  maybe only do for some threshold # of taxa
  would we want a cache?
  volatile or no?*/

		/*create new names*/
		for (i=ntaxa;i<(ntaxa+nent-1);i++) new_num_seqs[i]=pre_tree_made[i]=0;
		all_nodes_done=0;
		for (;;){
			for (i=(nent-2);i>=0;--i){
				d1=nodes[i][0];
				d2=nodes[i][1];
				if ((pre_tree_made[d1]) && (pre_tree_made[d2]) && (!pre_tree_made[ntaxa+i])) {
					pre_tree_made[ntaxa+i]=1;
					++all_nodes_done;
					/*make names here*/
					if (new_changes[ntaxa+i].name) free(new_changes[ntaxa+i].name);
					if (strcmp (new_changes[d1].name,new_changes[d2].name) > 0) new_changes[ntaxa+i].name=(char *)other_pair_names(new_changes[d2].name,new_changes[d1].name);
					else new_changes[ntaxa+i].name=(char *)other_pair_names(new_changes[d1].name,new_changes[d2].name);
					new_num_seqs[ntaxa+i]=new_num_seqs[d1]+new_num_seqs[d2];
				}
			}
			if (all_nodes_done>nent-2) break;
		}

		/*compare names*/
		for (i=ntaxa;i<(ntaxa+nent-1);i++) pre_tree_made[i]=0;
		for (i=ntaxa;i<(ntaxa+nent-1);i++) if (changes[i].name) {
			for (j=ntaxa;j<(ntaxa+nent-1);j++) {
				if (num_seqs[i]==new_num_seqs[j]) if (!pre_tree_made[j]) {
					/*
      holder=strcmp_clade(changes[i].name,new_changes[j].name);
      */
					THANG
					    if (!holder) {
						values->tree_made[j]=1;
						pre_tree_made[j]=1;
						new_changes[j].gaps=changes[i].gaps;
						if (!values->ttr) new_changes[j].changes=changes[i].changes;
						else {
							new_changes[j].transitions=changes[i].transitions;
							new_changes[j].transversions=changes[i].transversions;
						}
					}
				}
			}
		}
		/* mark descendents ?*/
		/*copy over new names and lengths*/
		for (i=ntaxa;i<(ntaxa+nent-1);i++) {
			if (changes[i].name) free(changes[i].name);
			changes[i].name=(char *)malloc((1+strlen(new_changes[i].name))*(sizeof(char)));
			assert((int)changes[i].name);
			changes[i].name=(char *)strcpy(changes[i].name,new_changes[i].name);
			num_seqs[i]=new_num_seqs[i];
			changes[i].gaps=new_changes[i].gaps;
			if (!values->ttr) changes[i].changes=new_changes[i].changes;
			else {
				changes[i].transitions=new_changes[i].transitions;
				changes[i].transversions=new_changes[i].transversions;
			}
		}
	}
}


void make_nodes(tree_rep,nodes,ntaxa,nbases,taxa,nent,n_gaps,values,return_value)
int *return_value;
int **nodes;
int *tree_rep;
int ntaxa,nbases,nent;
int **taxa,n_gaps;
parameters *values;
{
	int ii,dummy1,dummy2;
	int nn,d1,d2,all_nodes_done;
	int change_counts,mm;
	register int *td1, *td2, *td3;

	nodes[0][0]=0;
	nodes[0][1]=ntaxa+1;
	nodes[1][0]=1;
	nodes[1][1]=2;
	for (ii=3;ii<nent;ii++){
		nodes[ii-1][0]=ii;
		dummy1=(tree_rep[ii-3])/2;
		dummy2=(tree_rep[ii-3])%2;
		nodes[ii-1][1]=nodes[dummy1][dummy2];
		nodes[dummy1][dummy2]=ntaxa+ii-1;
	}
	for (nn=ntaxa;nn<(ntaxa+nent-1);++nn) values->tree_made[nn]=0;
	if (!values->delta) {
		if ((n_gaps+nbases)>200) pre_diagnose_tree(tree_rep,nodes,ntaxa,nbases,nent,n_gaps,values,taxa);
	}
	/*else all_nodes_done=0;*/
make_the_node:
	/*for (;;){*/
	for (nn=(nent-2);nn>=0;--nn){
		d1=nodes[nn][0];
		d2=nodes[nn][1];
		if ((values->tree_made[d1]) && (values->tree_made[d2]) && (!values->tree_made[ntaxa+nn])) {
			td1=taxa[d1];
			td2=taxa[d2];
			td3=taxa[ntaxa+nn];
			values->tree_made[ntaxa+nn]=1;
			++all_nodes_done;
			if (!values->ttr) {
				changes[ntaxa+nn].changes=changes[d1].changes+changes[d2].changes;
				for (mm=0;mm<nbases;mm++) {
					td3[mm]=(td1[mm]&td2[mm]);
					if (!td3[mm]) {
						td3[mm]=(td1[mm]|td2[mm]);
						++changes[ntaxa+nn].changes;
					}
				}
			}
			else {
				changes[ntaxa+nn].transitions=changes[d1].transitions+changes[d2].transitions;
				for (mm=0;mm<nbases-values->n_transv;mm++) {
					td3[mm]=(td1[mm]&td2[mm]);
					if (!td3[mm]) {
						td3[mm]=(td1[mm]|td2[mm]);
						++changes[ntaxa+nn].transitions;
					}
				}
				changes[ntaxa+nn].transversions=changes[d1].transversions+changes[d2].transversions;
				for (mm=nbases-values->n_transv;mm<nbases;mm++) {
					td3[mm]=(td1[mm]&td2[mm]);
					if (!td3[mm]) {
						td3[mm]=(td1[mm]|td2[mm]);
						++changes[ntaxa+nn].transversions;
					}
				}
			}
			changes[ntaxa+nn].gaps=changes[d1].gaps+changes[d2].gaps;
			for (mm=nbases;mm<(nbases+n_gaps);mm++) {
				td3[mm]=(td1[mm]&td2[mm]);
				if (!td3[mm]) {
					td3[mm]=(td1[mm]|td2[mm]);
					++changes[ntaxa+nn].gaps;
				}
			}
		}
	}
	for (nn=ntaxa;nn<((ntaxa+nent)-1);++nn) if (values->tree_made[nn]==0) goto make_the_node;
	/*
		 if (all_nodes_done>nent-2) break;
		    }
*/

	if (!values->delta) {
		if (values->ttr) (*return_value)=((changes[ntaxa].gaps*values->gap_cost)+(changes[ntaxa].transitions*values->transition)+(changes[ntaxa].transversions*values->transversion));
		else (*return_value)=(changes[ntaxa].changes*values->change_cost)+(changes[ntaxa].gaps*values->gap_cost);
	}
	else *return_value=(diagnose_tree(taxa,nodes,nbases,ntaxa,nent,n_gaps,values,changes[ntaxa].gaps));

}

/*
int make_nodes(tree_rep,nodes,ntaxa,nbases,taxa,nent,n_gaps,values)
int **nodes;
int *tree_rep;
int ntaxa,nbases,nent;
int **taxa,n_gaps;
parameters *values;
{
     int i,dummy1,dummy2;
  int n,d1,d2;
  int transition_counts,transversion_counts,gap_counts,change_counts,m;
  register int *td1, *td2, *td3;


     nodes[0][0]=0;
     nodes[0][1]=ntaxa+1;
     nodes[1][0]=1;
     nodes[1][1]=2;

for (i=3;i<nent;i++)
	  {
	  nodes[i-1][0]=i;
	  dummy1=(tree_rep[i-3])/2;
	  dummy2=(tree_rep[i-3])%2;
	  nodes[i-1][1]=nodes[dummy1][dummy2];
	  nodes[dummy1][dummy2]=ntaxa+i-1;
	  }

transversion_counts=transition_counts=gap_counts=change_counts=0;

     for (n=ntaxa;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;
	  make_the_node:
	  for (n=(nent-2);n>=0;--n){
	       d1=nodes[n][0];
	       d2=nodes[n][1];
	       if ((values->tree_made[d1]==1) && (values->tree_made[d2]==1) && (values->tree_made[ntaxa+n]==0)) {
	   td1=taxa[d1];
	 td2=taxa[d2];
	td3=taxa[ntaxa+n];
		    values->tree_made[ntaxa+n]=1;
		    if (!values->ttr) {
			 for (m=0;m<nbases;m++) {
			      if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
				   else {
					td3[m]=(td1[m]|td2[m]);
					++change_counts;
					}
				   }
			  }
	       else {
		    for (m=0;m<nbases-values->n_transv;m++) {
			 if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
			      else {
				   td3[m]=(td1[m]|td2[m]);
				   ++transition_counts;
				   }
			      }
		    for (m=nbases-values->n_transv;m<nbases;m++) {
			 if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
			      else {
				   td3[m]=(td1[m]|td2[m]);
				   ++transversion_counts;
				   }
			      }
		}

		for (m=nbases;m<(nbases+n_gaps);m++) {
		    if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
			 else {
			      td3[m]=(td1[m]|td2[m]);
			      ++gap_counts;
			      }
			 }

		    }
	  }
for (n=ntaxa;n<((ntaxa+nent)-1);++n) if (values->tree_made[n]==0) goto make_the_node;

if (!values->delta) return ((change_counts*values->change_cost)+(gap_counts*values->gap_cost)+(transition_counts*values->transition)+(transversion_counts*values->transversion));
if (values->delta)  return (diagnose_tree(taxa,nodes,nbases,ntaxa,nent,n_gaps,values,gap_counts));
}
*/

void reorder_tree(sequence,ntaxa,m,n_gaps)
look_up_thang **sequence;
int ntaxa,m,n_gaps;
{
	int i,j,temp,did_something;
	int **seq_holder;
	int *new_order,*ran_nums;

	ran_nums=(int *)malloc(ntaxa*sizeof(int));
	assert((int)ran_nums);
	new_order=(int *)malloc(ntaxa*sizeof(int));
	assert((int)new_order);
	for (i=0;i<ntaxa;i++) {
		new_order[i]=i;
		ran_nums[i]=rand();
	}

	seq_holder=(int **)malloc(ntaxa*sizeof(int *));
	assert((int)seq_holder);
	/*copy over*/
	for (i=0;i<ntaxa;i++) {
		seq_holder[i]=(int *)malloc((m+n_gaps)*sizeof(int));
		assert((int)seq_holder[i]);
		for (j=0;j<(m+n_gaps);j++) seq_holder[i][j]=sequence[i][j].base;
	}

	/*get new order*/
do_reorder_thing_tree:
	;
	did_something=0;
	for (i=0;i<ntaxa-1;i++) {
		if (ran_nums[i] < ran_nums[i+1]) {
			temp=ran_nums[i];
			ran_nums[i]=ran_nums[i+1];
			ran_nums[i+1]=temp;
			temp=new_order[i];
			new_order[i]=new_order[i+1];
			new_order[i+1]=temp;
			did_something=1;
		}
	}
	if (did_something) goto do_reorder_thing_tree;

	/*copy back*/
	for (i=0;i<ntaxa;i++) for (j=0;j<(m+n_gaps);j++) sequence[i][j].base=seq_holder[new_order[i]][j];
	/*free*/
	for (i=0;i<ntaxa;i++) free(seq_holder[i]);
	free(seq_holder);
	free(new_order);
	free(ran_nums);
}


