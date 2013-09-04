/*
Copyright 1992 Ward Wheeler all rights reserved

gets heuristic alignment by building
*/

#include "align3.h"

alignment **make_heur_align(a,ntax,best_aligns,values,previous_solutions,num_best_aligns,parallel_modify)
alignment **a,**best_aligns;
int ntax,previous_solutions;
int *num_best_aligns,parallel_modify;
parameters *values;
{
	int i,j,temp,l,ntaxa;
	int **nodes,extra;
	alignment **old_best_aligns;

	ntaxa=ntax+1; /*to fool into doing all rootings keeping '0' the root just never do the last alignment*/

	nodes=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)nodes);
	for (i=0;i<(ntaxa-1);i++) {
		nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int)nodes[i]);
	}
	nodes[0][0]=0;
	nodes[0][1]=ntaxa+1;
	nodes[1][0]=1;
	nodes[1][1]=2;

	for (i=0;i<previous_solutions;i++) if (!best_aligns[i]) {
		previous_solutions=i+1;
		break;
	}

	if (best_aligns[0]) {
		old_best_aligns=(alignment **)malloc(previous_solutions*sizeof(alignment *));
		assert((int)old_best_aligns);
		for (i=0;i<previous_solutions;i++) {
			old_best_aligns[i]=make_align(best_aligns[i]);
			assert((int)old_best_aligns[i]);
		}
	}
	else old_best_aligns=NULL;
	
	/*added to get rid of best_aligns previous solutions don't seem to need them*/
	for (i=0;i<previous_solutions;i++) if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
	values->number=l=0;
	if (values->do_dave_align) best_aligns=make_dave_align(a,ntaxa,nodes,&l,best_aligns,values,parallel_modify);
	else {
		if (!values->get_heur3) best_aligns=get_align_bound(a,ntaxa,nodes,&l,best_aligns,values);
		else best_aligns=make_heur_align_new(a,ntaxa,nodes,&l,best_aligns,values,parallel_modify);
		}

	for (i=0;i<(ntaxa-1);i++) free(nodes[i]);
	free(nodes);
	extra=0;

	if (old_best_aligns) {
		fprintf(stderr,"Old best score=%d new=%d\n",old_best_aligns[0]->score,best_aligns[0]->score);
		if (old_best_aligns[0]->score < best_aligns[0]->score) {
			/*save old ones*/
			for (i=0;i<l;i++) {
				if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
				best_aligns[i]=NULL;
			}
			for (i=0;i<previous_solutions;i++) best_aligns[i]=make_align(old_best_aligns[i]);
			if (values->VERBOSE) fprintf(stderr,"	Previously found alignments are better.\n");
			l=previous_solutions;
		}
		else if (old_best_aligns[0]->score==best_aligns[0]->score) {
			/*get all unique*/
			/*if temp!=0 not there*/
			for (i=0;i<previous_solutions;i++) {
				temp=1;
				for (j=0;j<l;j++) temp*=compare_aligns(old_best_aligns[i],best_aligns[j]);
				if (temp) {
					if (l==values->keep_aligns) {
						if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,
						    (values->keep_aligns+1)*sizeof(alignment *)))!=NULL)) {
							++values->keep_aligns;
							best_aligns[values->keep_aligns-1]=NULL;
						}
					}
					else if (l<values->keep_aligns) {
						best_aligns[(l+extra)]=make_align(old_best_aligns[i]);
						++extra;
					}
					else if (values->VERBOSE) fprintf(stderr,"0-Overflow in heuristic alignment -- but that's OK.\n");
				}
			}
			if (values->VERBOSE) fprintf(stderr,"	Found %d new alignments.\n",l+extra-previous_solutions);
		}
		else if (values->VERBOSE) fprintf(stderr,"	Found %d better alignments.\n",l);
		for (i=0;i<previous_solutions;i++) old_best_aligns[i]=dump_align(old_best_aligns[i]);
		free(old_best_aligns);
	}
	if (!values->jackboot) if (values->VERBOSE) fprintf(stderr,"Heuristic build yielded %d alignments of cost %d\n",l+extra,best_aligns[0]->score);
	(*num_best_aligns)=l+extra;
	/*print_alignment(best_aligns[0]);*/
	return best_aligns;
}/*end*/

alignment **get_align_bound(a,ntaxa,nodes,l,best_aligns,values)
alignment **a,**best_aligns;
int ntaxa,*l;
int **nodes;
parameters *values;
{
	int i,j,k,m,z,y,counter,value,*temp_rep;
	int cur_val,temp,**tree_rep,**old_rep;
	char **from,overflow;
	int **n_nodes,keep_aligns_holder;
	realloc_align_holder *align_passer;
	int *intermediate_scores,intermediate_best_score;
	alignment **intermediate_align;
	int num_best_aligns, sub_taxa, next_to_do, num_left, info;
	int bufid, bytes, type, source, max_to_do, index, grain, jj;
	int i_thang;

	align_passer=(realloc_align_holder *)malloc(sizeof(realloc_align_holder));
	assert((int)align_passer);
	n_nodes=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)n_nodes);
	for (i=0;i<(ntaxa-1);i++) {
		n_nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int)n_nodes[i]);
	}
	tree_rep=(int **)malloc(values->keep_aligns *sizeof(int *));
	assert((int)tree_rep);
	old_rep=(int **)malloc(values->keep_aligns*sizeof(int *));
	assert((int)old_rep);
	for (i=0;i<values->keep_aligns;i++) {
		tree_rep[i]=(int *)malloc((ntaxa-3)*sizeof(int));
		assert((int)tree_rep[i]);
		old_rep[i]=(int *)malloc((ntaxa-3)*sizeof(int));
		assert((int)old_rep[i]);
	}
	temp_rep=(int *)malloc((ntaxa-3)*sizeof(int));
	assert((int)temp_rep);

	/*these are paralell formulations to be modified later but better work in sequential*/
	intermediate_scores=(int *)malloc(((2*(ntaxa-3))+3+1)*sizeof(int));
	assert((int) intermediate_scores);
	intermediate_align=(alignment **)malloc(((2*(ntaxa-3))+3+1)*sizeof(alignment *));
	assert((int) intermediate_align);
	for (i=0;i<((2*(ntaxa-3))+3+1);i++) intermediate_align[i]=NULL;

	/*
	change to keep track of were swaps are from to avoid
	reconstructing alignments
	also do something to keep equal and best within a swap not all
*/

	overflow=y=counter=0;
	if (values->VERBOSE) fprintf(stderr,"\nAdding sequence ");
	for (i=0;i<ntaxa-3;i++) {
		value=intermediate_best_score=HUGE_COST;
		if (values->VERBOSE) {
			fprintf(stderr,"%d ",i+3);
			fflush(stderr);
		}
		if (i==(ntaxa-4)) if (values->VERBOSE) fprintf(stderr,"\n");
		/*iterate through all equally costly intermediates*/
		(*l)=0;
		for (z=0;z<=counter;z++) {
			for (j=0;j<i;j++) temp_rep[j]=old_rep[z][j];
			for (j=1;j<=((2*i)+3);j++) {
				intermediate_scores[j]=HUGE_COST;
				if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
				intermediate_align[j]=NULL;
			}
			/*Can add paralellization loops in heren dby splitting the loop functions into three*/
if (PARALLEL) do_parallel_shtick(a,values,&i,&ntaxa,&j,intermediate_scores,&intermediate_best_score,intermediate_align,temp_rep);
			/*sequential*/
else {
	for (j=1;j<=((2*i)+3);j++) {
		if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
		temp_rep[i]=j;
		/*fprintf(stderr,"%s->%d\n",a[i]->name,a[i]->n_seqs);*/
		intermediate_scores[j]=all_make_nodes(a,temp_rep,nodes,ntaxa,i+4,values);
		if (intermediate_scores[j]<=intermediate_best_score) {
			if (intermediate_scores[j]<intermediate_best_score) {
				intermediate_best_score=intermediate_scores[j];
				all_other_make_nodes2(temp_rep,nodes,ntaxa,i+4,n_nodes);
				if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
				intermediate_align[j]=make_align(a[n_nodes[0][1]]);
				}
			else if ((!values->aquick) && (intermediate_scores[j]!=HUGE_COST)) {
				all_other_make_nodes2(temp_rep,nodes,ntaxa,i+4,n_nodes);
				if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
				intermediate_align[j]=make_align(a[n_nodes[0][1]]);
			}
		}
		else if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
	}
}
value=intermediate_best_score;
num_best_aligns=0;
for (j=1;j<=((2*i)+3);j++) if (intermediate_scores[j]==intermediate_best_score) ++num_best_aligns;

if (!(*l)) assert ((int) num_best_aligns);
if (z==0) {
	for (m=0;m<(*l);m++) best_aligns[m]=dump_align(best_aligns[m]);
	(*l)=0;
	}
else if (z>0) {
	if (intermediate_best_score<best_aligns[0]->score) {
		for (m=0;m<(*l);m++) best_aligns[m]=dump_align(best_aligns[m]);
		(*l)=0;
		}
	}
for (j=1;j<=((2*i)+3);j++) {
	if ((intermediate_scores[j]==intermediate_best_score) && (intermediate_align[j])) {
		if ((*l)==0) {
			(*l)=1;
			overflow=y=0;
			for (k=0;k<ntaxa-3;k++) tree_rep[y][k]=temp_rep[k];
			tree_rep[y][i]=j;
			if ((i==(ntaxa-4)) && (values->VERBOSE)) fprintf(stderr,"Found one (better) at cost %d\n",intermediate_best_score);
			if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
			best_aligns[0]=make_align(intermediate_align[j]);
		}
		else if (!values->aquick) {
			temp=1;
			for (m=0;m<(*l);m++) if (best_aligns[m]) temp*=compare_aligns(intermediate_align[j],best_aligns[m]);
			if (temp) {
				if (((*l)>values->keep_aligns-2)||((y+1)>values->keep_aligns-1)) {
					if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment
					    *)))!=NULL) && ((tree_rep=(int **)realloc(tree_rep,(values->keep_aligns+1)*sizeof(int
					*)))!=NULL) && ((old_rep=(int **)realloc(old_rep,(values->keep_aligns+1)*sizeof(int
					*)))!=NULL)) {
						++values->keep_aligns;
						best_aligns[values->keep_aligns-1]=NULL;
						tree_rep[values->keep_aligns-1]=(int *)malloc((ntaxa-3)*sizeof(int));
						assert((int)tree_rep[values->keep_aligns-1]);
						old_rep[values->keep_aligns-1]=(int *)malloc((ntaxa-3)*sizeof(int));
						assert((int)old_rep[values->keep_aligns-1]);
					}
				}
				if ((y+1)<values->keep_aligns) {
					++y;
					for (k=0;k<ntaxa-3;k++) tree_rep[y][k]=temp_rep[k];
					tree_rep[y][i]=j;
				}
				else if (!overflow) {
					if (values->rep_error) fprintf(stderr,"Overflow in heuristic align -- but that's OK.\n");
					overflow=1;
				}
				if ((*l)>(values->keep_aligns-1)) {
					if (!overflow) {
						if (values->rep_error) fprintf(stderr,"Overflow in heuristic align -- but that's OK.\n");
						overflow=1;
					}
				}
				else {
					if (i==(ntaxa-4)) if (values->VERBOSE) fprintf(stderr,"Found another for %d alignments\n",
					    (*l)+1);
					if (best_aligns[(*l)]) best_aligns[(*l)]=dump_align(best_aligns[(*l)]);
					best_aligns[(*l)]=make_align(intermediate_align[j]);
					++(*l);
				}
			}
		}
	}
}/*j*/
		}/*z*/
		/*copy tree to old*/
		/*fprintf(stderr,"\ncopying\n");*/
		for (k=0;k<=y;k++) for (j=0;j<ntaxa-3;j++) old_rep[k][j]=tree_rep[k][j];
		counter=y;
		y=overflow=0;

		/*Swap while building here*/
		if ((values->swap_while_add) && (i<(ntaxa-4)) && (i>0)) {
			keep_aligns_holder=values->keep_aligns;
			align_passer=do_addition_swap(&counter,old_rep,a,best_aligns,nodes,n_nodes,l,values,i,ntaxa,align_passer);
			if (values->aquick) counter=0;
			best_aligns=align_passer->best;
			old_rep=align_passer->rep;
			if (values->keep_aligns>keep_aligns_holder) {
				tree_rep=(int **)realloc(tree_rep,values->keep_aligns*sizeof(int *));
				for (j=keep_aligns_holder-1;j<values->keep_aligns;j++) {
					tree_rep[j]=(int *)malloc((ntaxa-3)*sizeof(int));
					assert((int)tree_rep[j]);
				}
			}
		}
		if ((values->align_swap) && (i==(ntaxa-4))) {
			keep_aligns_holder=values->keep_aligns;
			align_passer=do_addition_swap(&counter,old_rep,a,best_aligns,nodes,n_nodes,l,values,i,ntaxa,align_passer);
			if (values->aquick) counter=0;
			best_aligns=align_passer->best;
			old_rep=align_passer->rep;
			if (values->keep_aligns>keep_aligns_holder) {
				tree_rep=(int **)realloc(tree_rep,values->keep_aligns*sizeof(int *));
				for (j=keep_aligns_holder-1;j<values->keep_aligns;j++) {
					tree_rep[j]=(int *)malloc((ntaxa-3)*sizeof(int));
					assert((int)tree_rep[j]);
				}
			}
		}
}/*i*/

/*swapping part*/
	if (values->align_multi_swap>1) {
		keep_aligns_holder=values->keep_aligns;
		align_passer=deep_swap_align(nodes,n_nodes,a,counter+1,old_rep,best_aligns,l,ntaxa,values->align_multi_swap,values,align_passer);
		best_aligns=align_passer->best;
		old_rep=align_passer->rep;
		if (values->keep_aligns>keep_aligns_holder) {
			tree_rep=(int **)realloc(tree_rep,values->keep_aligns*sizeof(int *));
			for (j=keep_aligns_holder-1;j<values->keep_aligns;j++) {
				tree_rep[j]=(int *)malloc((ntaxa-3)*sizeof(int));
				assert((int)tree_rep[j]);
			}
		}
	}

if (values->align_node_swap) {
		keep_aligns_holder=values->keep_aligns;
		align_passer=do_node_swap(&counter,old_rep,a,best_aligns,nodes,n_nodes,l,values,ntaxa-4,ntaxa,align_passer);
		if (values->aquick) counter=0;
		best_aligns=align_passer->best;
		old_rep=align_passer->rep;
		if (values->keep_aligns>keep_aligns_holder) {
			tree_rep=(int **)realloc(tree_rep,values->keep_aligns*sizeof(int *));
			for (j=keep_aligns_holder-1;j<values->keep_aligns;j++) {
				tree_rep[j]=(int *)malloc((ntaxa-3)*sizeof(int));
				assert((int)tree_rep[j]);
			}
		}
	}


	for (i=0;i<values->keep_aligns;i++) {
		if (tree_rep[i]) free(tree_rep[i]);
		if (old_rep[i]) free(old_rep[i]);
	}
	free(tree_rep);
	free(temp_rep);
	free(old_rep);
	for (i=0;i<(ntaxa-1);i++) free(n_nodes[i]);
	free(n_nodes);
	free(align_passer);

	free(intermediate_scores);
	for (i=0;i<((2*(ntaxa-3))+3+1);i++) if (intermediate_align[i]) intermediate_align[i]=dump_align(intermediate_align[i]);
	free(intermediate_align);
	return best_aligns;
}

realloc_align_holder *deep_swap_align(nodes,best_nodes,a,num_to_swap,old_reps,best_aligns,num_best,ntaxa,depth,values,align_passer)
int num_to_swap,**old_reps,*num_best,ntaxa,depth;
int **best_nodes,**nodes;
alignment **best_aligns,**a;
parameters *values;
realloc_align_holder *align_passer;
{
	int i,j,new_ones=0,cur_b_tree;
	char *to_swap,*tuple_holder,**already_swapped;
	int overflow=0;

	to_swap=(char *)malloc((ntaxa-3)*sizeof(char));
	assert((int)to_swap);
	tuple_holder=(char *)malloc((ntaxa-3)*sizeof(char));
	assert((int)tuple_holder);
	already_swapped=(char **)malloc(values->keep_aligns*sizeof(char *));
	assert((int)already_swapped);
	for (i=0;i<values->keep_aligns;i++) {
		already_swapped[i]=(char *)malloc((ntaxa-3)*sizeof(char));
		assert((int)already_swapped[i]);
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
		if (i>values->keep_aligns) break;
		if (values->VERBOSE) fprintf(stderr,"Swapping alignment %d of %d\n",i+1,num_to_swap+new_ones);
		for (j=0;j<ntaxa-3;j++) *(tuple_holder+j)=0;
		cur_b_tree=best_aligns[0]->score;
		align_passer=get_tuples_and_loop(&overflow,i,already_swapped,best_nodes,nodes,&cur_b_tree,a,best_aligns,num_best,
		    old_reps[i],&new_ones,num_to_swap,old_reps,0,depth,tuple_holder,ntaxa,values,align_passer);
		old_reps=align_passer->rep;
		already_swapped=align_passer->swap;
		best_aligns=align_passer->best;
	}

	free(to_swap);
	free(tuple_holder);
	for (i=0;i<values->keep_aligns;i++) free(already_swapped[i]);
	free(already_swapped);

	align_passer->best=best_aligns;
	align_passer->rep=old_reps;
	return align_passer;
}

realloc_align_holder *get_tuples_and_loop(overflow,alignment_swapping,already_swapped,best_nodes,nodes,cur_b_tree,a,best_aligns,
num_best,current_rep,new_ones,num_to_swap,old_reps,current_depth,depth,tuple_holder,ntaxa,values,align_passer)
alignment **best_aligns,**a;
int *num_best,*current_rep,*new_ones,num_to_swap,**old_reps,current_depth,depth,ntaxa;
char *tuple_holder;
parameters *values;
int *cur_b_tree,alignment_swapping,*overflow;
int **best_nodes,**nodes;
char **already_swapped;
realloc_align_holder *align_passer;
{
	int i,j,num_tuples;
	int *tree_rep;
	char doit;

	/*this needs to be checked as far as clearing and maintaining previus values*/
	for (i=current_depth;i<ntaxa-3;i++) {
		*(tuple_holder+i)=1;
		if (current_depth < (depth-1)) {
			align_passer=get_tuples_and_loop(overflow,alignment_swapping,already_swapped,best_nodes,nodes,cur_b_tree,
			    a,best_aligns,num_best,current_rep,new_ones,num_to_swap,old_reps,current_depth+1,depth,tuple_holder,
			    ntaxa,values,align_passer);
			old_reps=align_passer->rep;
			already_swapped=align_passer->swap;
			best_aligns=align_passer->best;
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
				if (values->VERBOSE) fprintf(stderr," swapping sequences ");
				for (j=0;j<ntaxa-3;j++) if ((*(tuple_holder+j)==1) && (values->VERBOSE)) fprintf(stderr,"%d ",
				    j+3);
				if (values->VERBOSE) fprintf(stderr,"simultaneously\n");
				align_passer=special_all_loop(overflow,alignment_swapping,already_swapped,old_reps,num_to_swap,
				    new_ones,tuple_holder,a,tree_rep,0,ntaxa,nodes,num_best,cur_b_tree,best_nodes,best_aligns,
				    values,align_passer);
				old_reps=align_passer->rep;
				already_swapped=align_passer->swap;
				best_aligns=align_passer->best;
				free(tree_rep);
			}
		}
		*(tuple_holder+i)=0;
	}
	align_passer->rep=old_reps;
	align_passer->swap=already_swapped;
	align_passer->best=best_aligns;
	return align_passer;
}

realloc_align_holder *special_all_loop(overflow,alignment_swapping,already_swapped,old_reps,num_to_swap,new_ones,tuple_holder,
a,tree_rep,element,ntax,nodes,l,cur_b_tree,best_nodes,best_aligns,values,align_passer)
int element,ntax,*tree_rep;
int **best_nodes,**nodes;
int **old_reps,alignment_swapping;
int *l,*cur_b_tree,num_to_swap,*new_ones,*overflow;
alignment **a,**best_aligns;
parameters *values;
char *tuple_holder,**already_swapped;
realloc_align_holder *align_passer;
{
	int i,j,cur_tree,temp,doit;

	for (i=((2*(element+4))-5);i>=1;--i)
	{
		if (*(tuple_holder+element)==1) tree_rep[element]=i;
		else i=1;
		if (element<(ntax-4)) {
			align_passer=special_all_loop(overflow,alignment_swapping,already_swapped,old_reps,num_to_swap,new_ones,
			    tuple_holder,a,tree_rep,element+1,ntax,nodes,l,cur_b_tree,best_nodes,best_aligns,values,align_passer);
			old_reps=align_passer->rep;
			already_swapped=align_passer->swap;
			best_aligns=align_passer->best;
		}
		else {
			doit=0;
			for (j=0;j<ntax-3;j++) if (tree_rep[j]!=(*(old_reps[alignment_swapping]+j))) doit=1;
			if (doit==0) goto all_end_of_all_loop_special;
			cur_tree=all_make_nodes(a,tree_rep,nodes,ntax,ntax,values);
			if (cur_tree==HUGE_COST) goto all_end_of_all_loop_special;
			if (cur_tree < (*cur_b_tree)) {
				if (values->VERBOSE) fprintf(stderr,"Found better...");
				(*cur_b_tree) = cur_tree;
				if (values->VERBOSE) fprintf(stderr,"cost %d\n",*cur_b_tree);
				if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
				all_other_make_nodes(tree_rep,nodes,ntax,best_nodes);
				best_aligns[0]=make_align(a[best_nodes[0][1]]);
				/*free the rest*/
				for (j=1;j<(*l);j++) best_aligns[j]=dump_align(best_aligns[j]);
				*l=1;
				/*update the ones to swap*/
				if ((num_to_swap+(*new_ones))==values->keep_aligns) {
					if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment
					    *)))!=NULL) && ((old_reps=(int **)realloc(old_reps,(values->keep_aligns+1)*sizeof(int
					*)))!=NULL) && ((already_swapped=(char **)realloc(already_swapped,(values->keep_aligns+1)*sizeof(char
					*)))!=NULL)) {
						++values->keep_aligns;
						best_aligns[values->keep_aligns-1]=NULL;
						old_reps[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
						assert((int)old_reps[values->keep_aligns-1]);
						already_swapped[values->keep_aligns-1]=(char *)malloc((ntax-3)*sizeof(char));
						assert((int)already_swapped[values->keep_aligns-1]);
						for (j=0;j<ntax-3;j++) already_swapped[values->keep_aligns-1][j]=0;
					}
				}
				if ((num_to_swap+(*new_ones))<values->keep_aligns) *new_ones+=1;
				else if (*overflow==0) {
					if (values->rep_error) fprintf(stderr,"Overflow in alignment swapping--but that's OK.\n");
					*overflow=1;
				}
				for (j=0;j<ntax-3;j++) {
					*(old_reps[num_to_swap+(*new_ones)-1]+j)=tree_rep[j];
					*(already_swapped[num_to_swap+(*new_ones)-1]+j)=*(tuple_holder+j);
				}
			}
			else if (cur_tree == (*cur_b_tree)) {
				all_other_make_nodes(tree_rep,nodes,ntax,best_nodes);
				temp=1;
				for (j=0;j<(*l);j++) temp*=compare_aligns(a[best_nodes[0][1]],best_aligns[j]);
				if ((temp!=0) || ((*l)==0)) {
					if ((*l)==values->keep_aligns) {
						if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,
						    (values->keep_aligns+1)*sizeof(alignment *)))!=NULL) && ((old_reps=(int **)realloc(old_reps,
						    (values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((already_swapped=(char
						**)realloc(already_swapped,(values->keep_aligns+1)*sizeof(char *)))!=NULL)) {
							++values->keep_aligns;
							best_aligns[values->keep_aligns-1]=NULL;
							old_reps[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
							assert((int)old_reps[values->keep_aligns-1]);
							already_swapped[values->keep_aligns-1]=(char *)malloc((ntax-3)*sizeof(char));
							assert((int)already_swapped[values->keep_aligns-1]);
							for (j=0;j<ntax-3;j++) already_swapped[values->keep_aligns-1][j]=0;
						}
					}
					if ((*l)==values->keep_aligns) {
						if (*overflow==0) {
							if (values->rep_error) fprintf(stderr,"Overflow in alignment swapping--but that's OK.\n");
							*overflow=1;
						}
						goto all_end_of_all_loop_special;
					}
					if ((best_aligns[0]) && ((*l)==0)) best_aligns[0]=dump_align(best_aligns[0]);
					if (values->VERBOSE) fprintf(stderr,"Found another for %d alignments\n",(*l)+1);
					best_aligns[(*l)]=make_align(a[best_nodes[0][1]]);
					++(*l);
					if ((num_to_swap+(*new_ones))==values->keep_aligns) {
						if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,
						    (values->keep_aligns+1)*sizeof(alignment *)))!=NULL) && ((old_reps=(int **)realloc(old_reps,
						    (values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((already_swapped=(char
						**)realloc(already_swapped,(values->keep_aligns+1)*sizeof(char *)))!=NULL)) {
							++values->keep_aligns;
							best_aligns[values->keep_aligns-1]=NULL;
							old_reps[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
							assert((int)old_reps[values->keep_aligns-1]);
							already_swapped[values->keep_aligns-1]=(char *)malloc((ntax-3)*sizeof(char));
							assert((int)already_swapped[values->keep_aligns-1]);
							for (j=0;j<ntax-3;j++) already_swapped[values->keep_aligns-1][j]=0;
						}
					}
					if ((num_to_swap+(*new_ones))<values->keep_aligns) {
						*new_ones+=1;
						for (j=0;j<ntax-3;j++) {
							*(old_reps[num_to_swap+(*new_ones)-1]+j)=tree_rep[j];
							*(already_swapped[num_to_swap+(*new_ones)-1]+j)=*(tuple_holder+j);
						}
					}
					else if (*overflow==0) {
						if (values->rep_error) fprintf(stderr,"Overflow in alignment swapping--but that's OK.\n");
						*overflow=1;
					}
				}
			}
		}
all_end_of_all_loop_special:
		continue;
	}
	align_passer->rep=old_reps;
	align_passer->swap=already_swapped;
	align_passer->best=best_aligns;
	return align_passer;
}



void all_other_make_nodes2(tree_rep,onode,ntaxa,nent,node)
int **onode,**node;
int *tree_rep;
int ntaxa,nent;
{
	int i,dummy1,dummy2;

	/*node names start after the taxa*/
	for (i=0;i<(nent-1);++i)
	{
		node[i][0]=onode[i][0];
		node[i][1]=onode[i][1];
	}
	for (i=3;i<nent;++i)
	{
		node[i-1][0]=i;
		dummy1=(tree_rep[i-3])/2;
		dummy2=(tree_rep[i-3])%2;
		node[i-1][1]=node[dummy1][dummy2];
		node[dummy1][dummy2]=ntaxa+i-1;
	}
}

