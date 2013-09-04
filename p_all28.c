/*
Copyright 1992 Ward Wheeler all rights reserved

	program for getting and doing ALL alignment orders
*/

#include "align3.h"

alignment **do_all(a,ntax,bound,best_aligns,values,previous_solutions,num_best_sols)
alignment **a,**best_aligns;
int ntax,bound,previous_solutions,*num_best_sols;
parameters *values;
{
	int i,j,extra,temp,l,cur_b_tree;
	int **nodes;
	int **best_nodes;
	alignment **old_best_aligns;
	int *tree_rep,ntaxa;
	realloc_align_holder *align_passer;

	align_passer=(realloc_align_holder *)malloc(sizeof(realloc_align_holder));
	assert((int)align_passer);
	extra=0;
	ntaxa=ntax+1; /*to fool into doing all rootings keeping '0' the root just never do the last alignment*/
	cur_b_tree=bound;

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
	tree_rep=(int *)malloc((ntaxa-3)*sizeof(int));
	assert((int)tree_rep);

	best_nodes=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)best_nodes);
	for (i=0;i<(ntaxa-1);i++) {
		best_nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int)best_nodes[i]);
		}

	if (values->VERBOSE) fprintf(stderr,"Initial bound = %d\n",cur_b_tree);
	l=values->number=0;

	if (best_aligns[0]) {
		old_best_aligns=(alignment **)malloc(previous_solutions*sizeof(alignment *));
		assert((int)old_best_aligns);
		for (i=0;i<previous_solutions;i++) old_best_aligns[i]=make_align(best_aligns[i]);
		}
	else {
		old_best_aligns=NULL;
		previous_solutions=0;
		}
	values->in_bandb_loop=1;
	align_passer=all_loop(a,tree_rep,0,ntaxa,nodes,&l,&cur_b_tree,best_nodes,best_aligns,values,align_passer);
	best_aligns=align_passer->best;
	values->in_bandb_loop=0;
	if (values->VERBOSE) fprintf(stderr,"	100%% complete\n");

	if (old_best_aligns) {
			if (old_best_aligns[0]->score < best_aligns[0]->score) {
		/*save old ones*/
		for (i=0;i<l;i++) {
			if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
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
				if (temp!=0) {
					if (l==values->keep_aligns) {
						if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment *)))!=NULL)) {
						++values->keep_aligns;
							best_aligns[values->keep_aligns-1]=NULL;
						}
					}
					if (l<values->keep_aligns) {
						best_aligns[(l+extra)]=make_align(old_best_aligns[i]);
						++extra;
						}
					else if (values->VERBOSE) fprintf(stderr,"Overflow in heuristic alignment -- but that's OK.\n");
					}
				}
			if (values->VERBOSE) fprintf(stderr,"	Found %d new alignments.\n",l+extra-previous_solutions);
			}
		else if (values->VERBOSE) fprintf(stderr,"	Found %d better alignments.\n",l);
		if (old_best_aligns) {
			for (i=0;i<previous_solutions;i++) if (old_best_aligns[i]) old_best_aligns[i]=dump_align(old_best_aligns[i]);
			free(old_best_aligns);
			}
		}
	if (values->VERBOSE) fprintf(stderr,"Found %d best alignments",l+extra);
	if (cur_b_tree!=HUGE_COST) if (values->VERBOSE) fprintf(stderr," of cost %d\n",cur_b_tree);
	else if (values->VERBOSE) fprintf(stderr,"\n");

free(tree_rep);
for (i=0;i<(ntaxa-1);i++) {
	free(nodes[i]);
	}
free(nodes);
for (i=0;i<(ntaxa-1);i++) {
	free(best_nodes[i]);
	}
free(best_nodes);
free(align_passer);

*num_best_sols=l+extra;
return best_aligns;
}/*end do_all*/

realloc_align_holder *all_loop(a,tree_rep,element,ntax,nodes,l,cur_b_tree,best_nodes,best_aligns,values,align_passer)
int element,ntax,*tree_rep;
int **best_nodes,**nodes;
int *l,*cur_b_tree;
alignment **a;
alignment **best_aligns;
parameters *values;
realloc_align_holder *align_passer;
{
	int i,j,k,cur_tree;
	int done_trees;
	int temp;
	int *intermediate_scores,intermediate_best_score;
	alignment **intermediate_align;
	int num_best_aligns, sub_taxa, next_to_do, num_left, info;
	int bufid, bytes, type, source, max_to_do, index, grain, jj;
	int old_best_score;

	/*to parallelize keeping all processors busy and minimizing comunication
	     **remember quick must exist stuff**
		1) do all tree_reps for this round
		2) assign costs to all 
		3) for each 
			a) if >= bound skip
			b) if < bound and !all taxa recurse
			c) if all taxa 
				1) keep best
				2) reassign bound 
	*/
	intermediate_scores=(int *)malloc(((2*(element+4))-5+1)*sizeof(int));
	assert((int) intermediate_scores);
	
	if (element == (ntax-4)) {
		intermediate_align=(alignment **)malloc(((2*(element+4))-5+1)*sizeof(alignment *));
		assert((int) intermediate_align);
		for (i=0;i<((2*(element+4))-5+1);i++)	intermediate_align[i]=NULL;
		}
	
/*parallel tag send 5 recieve 7 for no aignment needed*/
if (PARALLEL) {
	grain=values->grain_size;
	if (element < (ntax-4)) {
		/*scores only when alignments not needed send 5 recieve 7*/
			num_left=max_to_do=(2*element)+3;
			for (i=1;((i<values->num_hosts) && (i<max_to_do));i++) {
				tree_rep[element]=i;
				/*fprintf(stderr," %d (grain=%d ntaxa=%d i=%d temp_rep[%d]=%d) ",i,grain,(*ntaxa),element,element,temp_rep[element]);*/
				pvm_initsend( PvmDataDefault );
				pvm_pkint(&grain,1,1);
				pvm_pkint(&(*cur_b_tree),1,1);
				pvm_pkint(&ntax,1,1);
				pvm_pkint(&element,1,1);
				pvm_pkint(tree_rep,element+1,1);
				pvm_pkint(&i,1,1);
				pvm_send(values->tids[i],PVM_DO_ALIGNMENT_SCORE_ONLY_BANDB);
				/*fprintf(stderr," sent\n");*/
			}
			while (num_left) {
				/*fprintf(stderr,"num_left=%d ",num_left);*/
				bufid=pvm_recv(-1,7);
				if (bufid) {
					info=pvm_bufinfo(bufid, &bytes, &type, &source);
					if (i<=max_to_do) {
						tree_rep[element]=i;
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&grain,1,1);
						pvm_pkint(&(*cur_b_tree),1,1);
						pvm_pkint(&ntax,1,1);
						pvm_pkint(&element,1,1);
						pvm_pkint(tree_rep,element+1,1);
						pvm_pkint(&i,1,1);
						pvm_send(source,PVM_DO_ALIGNMENT_SCORE_ONLY_BANDB);
						++i;
					}
					--num_left;
					pvm_upkint(&grain,1,1);
					pvm_upkint(&index,1,1);
					pvm_upkint(&(intermediate_scores[index]),1,1);
					/*if (intermediate_scores[index]<(*cur_b_tree)) (*cur_b_tree)=intermediate_scores[index];*/
				}
			}
		}
	/*unpack only if less than or equal to*/
	else {
		/*fprintf(stderr,"In end loop ");*/
		old_best_score=(*cur_b_tree);
		num_left=max_to_do=(2*element)+3;
			for (i=1;((i<values->num_hosts) && (i<max_to_do));i++) {
				tree_rep[element]=i;
				/*fprintf(stderr," %d (grain=%d ntaxa=%d i=%d temp_rep[%d]=%d) ",i,grain,(*ntaxa),element,element,temp_rep[element]);*/
				pvm_initsend( PvmDataDefault );
				pvm_pkint(&grain,1,1);
				pvm_pkint(&(*cur_b_tree),1,1);
				pvm_pkint(&ntax,1,1);
				pvm_pkint(&element,1,1);
				pvm_pkint(tree_rep,element+1,1);
				pvm_pkint(&i,1,1);
				pvm_send(values->tids[i],PVM_DO_ALIGNMENT_NOW_BANDB);
				/*fprintf(stderr," sent\n");*/
			}
			while (num_left) {
				/*fprintf(stderr,"num_left=%d ",num_left);*/
				bufid=pvm_recv(-1,4);
				if (bufid) {
					info=pvm_bufinfo(bufid, &bytes, &type, &source);
					if (i<=max_to_do) {
						tree_rep[element]=i;
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&grain,1,1);
						pvm_pkint(&(*cur_b_tree),1,1);
						pvm_pkint(&ntax,1,1);
						pvm_pkint(&element,1,1);
						pvm_pkint(tree_rep,element+1,1);
						pvm_pkint(&i,1,1);
						pvm_send(source,PVM_DO_ALIGNMENT_NOW_BANDB);
						++i;
					}
					--num_left;
					pvm_upkint(&grain,1,1);
					pvm_upkint(&index,1,1);
					pvm_upkint(&(intermediate_scores[index]),1,1);
					if (intermediate_scores[index]<=(*cur_b_tree)) {
						if (intermediate_scores[index] < (*cur_b_tree)) (*cur_b_tree)=intermediate_scores[index];
						if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
						intermediate_align[index]=unpack_align(intermediate_scores[index]);
						}/*end unpack alignment*/
					else pvm_freebuf(bufid);
				}	/*bufid*/
			}/*num_left*/		
		}/*all taxa*/
	}/*parallel*/
/*sequential*/
else {
	for (i=((2*(element+4))-5);i>=1;--i) {
		tree_rep[element]=i;
		intermediate_scores[i]=all_make_nodes(a,tree_rep,nodes,ntax,element+4,values);
		if ((element==(ntax-4)) && (intermediate_scores[i]<=(*cur_b_tree))) {
			if (intermediate_scores[i]<(*cur_b_tree)) {
				values->current_bound=(*cur_b_tree)=intermediate_scores[i];
				for (j=((2*(element+4))-5);j>=1;--j) if ((intermediate_align[j]) && (intermediate_scores[j]>(*cur_b_tree))) intermediate_align[j]=dump_align(intermediate_align[j]);
				}
			if (intermediate_align[i]) intermediate_align[i]=dump_align(intermediate_align[i]);
			all_other_make_nodes(tree_rep,nodes,ntax,best_nodes);
			intermediate_align[i]=make_align(a[best_nodes[0][1]]);
			}
		}
	}
/*for (i=((2*(element+4))-5);i>=1;--i) fprintf(stderr,"E%d+S%d ",element, intermediate_scores[i]);*/
for (i=((2*(element+4))-5);i>=1;--i) {
		tree_rep[element]=i;
		cur_tree=intermediate_scores[i];
		if (values->VERBOSE) if (i<((2*(element+4))-5)) {
			if (element==0) fprintf(stderr,"	%3d%% complete.\n",(100*(3-i))/3);
			if (element==1) if (ntax>5) fprintf(stderr,"	%3d%% complete.\n",(100*((5*(3-tree_rep[0]))+(5-i)))/15);
			if (element==2) if (ntax>6) fprintf(stderr,"	%3d%% complete.\n",(100*((35*(3-tree_rep[0]))+(7*(5-tree_rep[1]))+(7-i)))/105);
			}
		if (element<(ntax-4)) {
			/*if no groups file but no gap cost counting*/
			if ((!values->groups)&&(!values->phylo_gap)) {
				align_passer=all_loop(a,tree_rep,element+1,ntax,nodes,l,cur_b_tree,best_nodes,best_aligns,values,align_passer);
				best_aligns=align_passer->best;
				}
			/*groups file and gaps not counting*/
			else if (!values->phylo_gap) {
				if (cur_tree==HUGE_COST) goto all_end_of_all_loop;
				align_passer=all_loop(a,tree_rep,element+1,ntax,nodes,l,cur_b_tree,best_nodes,best_aligns,values,align_passer);
				best_aligns=align_passer->best;
				}
			/*gaps counting*/
			else 	{
				if (cur_tree==HUGE_COST) goto all_end_of_all_loop;
				if (cur_tree > (*cur_b_tree)) goto all_end_of_all_loop;
				else {
					align_passer=all_loop(a,tree_rep,element+1,ntax,nodes,l,cur_b_tree,best_nodes,best_aligns,values,align_passer);
					best_aligns=align_passer->best;
					}
				}
			}
		else if (intermediate_align[i]) {
			if (cur_tree==HUGE_COST) goto all_end_of_all_loop;
			if (cur_tree < old_best_score) {
				old_best_score=cur_tree;
				if (values->VERBOSE) fprintf(stderr,"Found better...");
				if (values->VERBOSE) fprintf(stderr,"cost %d\n",*cur_b_tree);
				if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
				best_aligns[0]=make_align(intermediate_align[i]);
				/*free the rest*/
				for (j=1;j<(*l);j++) best_aligns[j]=dump_align(best_aligns[j]);
				*l=1;
				}
			else if (cur_tree == (*cur_b_tree)) {
				temp=1;
				for (j=0;j<(*l);j++) temp*=compare_aligns(intermediate_align[i],best_aligns[j]);
				if ((temp!=0) || ((*l)==0)) {
					if ((*l)==values->keep_aligns) {
						if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment *)))!=NULL)) {
							++values->keep_aligns;
							best_aligns[values->keep_aligns-1]=NULL;
							}
						}
					if ((*l)==values->keep_aligns) {
						if (values->rep_error) fprintf(stderr,"Overflow in exact alignment construction!\n");
						goto all_end_of_all_loop;
						}
					if ((best_aligns[0]) && ((*l)==0)) best_aligns[0]=dump_align(best_aligns[0]);
					if (values->VERBOSE) fprintf(stderr,"Found another for %d alignments\n",(*l)+1);
					best_aligns[(*l)]=make_align(intermediate_align[i]);
					++(*l);
					}
				}
			}
		all_end_of_all_loop:;
		}


	free(intermediate_scores);
	if (element == (ntax-4)) {
		for (i=0;i<((2*(element+4))-5+1);i++) if (intermediate_align[i]) intermediate_align[i]=dump_align(intermediate_align[i]);
		free(intermediate_align);
		}


align_passer->best=best_aligns;
return align_passer;
}


/*how many nodes??*/
int all_make_nodes(a,tree_rep,onode,ntaxa,nent,values)
int **onode;
int *tree_rep;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int i,dummy1,dummy2,ngroups2,return_value;
	int **node;
	int **cur_groups;
	char test;

	node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)node);
	for (i=0;i<(ntaxa-1);i++) {
		node[i]=(int *)malloc(2*sizeof(int));
		assert((int)node[i]);
		}

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
	/* insert topology test here */
	if (!values->groups) {
		return_value=all_diagnose_tree(a,node,ntaxa,nent,values);
		for (i=0;i<(ntaxa-1);i++) free(node[i]);
		free(node);
		return return_value;
		}
	else {
		cur_groups=get_cur_groups(node,ntaxa,nent,&ngroups2);
		test=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntaxa-1,nent-1);
		/*dump groups*/
		for (i=0;i<ngroups2;i++) free(cur_groups[i]);
		free(cur_groups);
		if (!test) {
			/*
			if (!a[node[0][1]]) a[node[0][1]]=make_align(a[0]);
			a[node[0][1]]->score=HUGE_COST;
			*/
			for (i=0;i<(ntaxa-1);i++) free(node[i]);
			free(node);
			return HUGE_COST;
			}
		else {
			return_value=all_diagnose_tree(a,node,ntaxa,nent,values);
			for (i=0;i<(ntaxa-1);i++) free(node[i]);
			free(node);
			return return_value;
			}
		}
}

void all_other_make_nodes(tree_rep,onode,ntaxa,new_nodes)
int **onode;
int *tree_rep;
int **new_nodes;
int ntaxa;
{
	int i,dummy1,dummy2;

	/*node names start after the taxa*/
	for (i=0;i<(ntaxa-1);++i)
	{
		new_nodes[i][0]=onode[i][0];
		new_nodes[i][1]=onode[i][1];
	}
	for (i=3;i<ntaxa;++i)
	{
		new_nodes[i-1][0]=i;
		dummy1=(tree_rep[i-3])/2;
		dummy2=(tree_rep[i-3])%2;
		new_nodes[i-1][1]=new_nodes[dummy1][dummy2];
		new_nodes[dummy1][dummy2]=ntaxa+i-1;
	}
}


int compare_aligns(a,b)
alignment *a, *b;
/*parameters *values;*/
{
	int i,j,same=0;

	if (!b) return 1;
	if (!a) return 1;
	/*new optimization check*/
	if ((a->n_seqs==1) && (b->n_seqs==1)) {
		if (!strcmp(a->name,b->name)) return 0;
		else return 1;
		}
		
	if (a->length!=b->length) return 1;
	if (a->n_seqs!=b->n_seqs) return 1;
	for (i=0;i<a->n_seqs;i++) {
		for (j=0;j<b->n_seqs;j++) {
		if (!strcmp(a->taxon_name[i],b->taxon_name[j])) {
				if (!strcmp(a->s[i],b->s[j])) ++same;
	}
			}
		}
	if (same!=a->n_seqs) return 1;
	else return 0;
}

int **get_cur_groups(node,ntaxa,nent,ngroups2)
int **node;
int ntaxa,nent,*ngroups2;
{
/*there is a group for each node of ntaxa
	remember that the ntaxa are shifted up one*/

int **groups;
int i,j,k;

/*Remeber ntaxa is ntax+1 for rooting*/
/*first determine number of groups*/
*ngroups2=nent-3;
groups=(int **)malloc((*ngroups2)*sizeof(int *));
assert((int)groups);
for (i=0;i<(*ngroups2);i++) {
		groups[i]=(int *)malloc(ntaxa*sizeof(int));
		assert((int)groups[i]);
	for (j=0;j<ntaxa;j++) *(groups[i]+j)=0;
	}

k=0;
for (i=1;i<nent-1;i++){
	/*for (j=0;j<2;j++) fprintf(stderr,"node[%d][%d]=%d\n",i,j,node[i][j]);*/
	if (i!=node[0][1]-ntaxa) {
		recurse_groups(node,ntaxa,nent,i,k,groups);
		++k;
		}
	}

for (i=0;i<*ngroups2;i++) {
	/*for (j=nent-2;j<ntaxa-1;j++) groups[i][j]=1; fix for nent < ntax */
	for (j=0;j<ntaxa-1;j++) *(groups[i]+ntaxa-1)+=(*(groups[i]+j));
	 }
/*fprintf(stderr,"\n");
for (i=0;i<*ngroups2;i++) {
	for (j=0;j<ntaxa;j++) fprintf(stderr,"%d ",(*(groups[i]+j)));
	fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");*/
return groups;
}

void recurse_groups(node,ntaxa,nent,i,k,groups)
int **node;
int ntaxa,nent,k,i;
int **groups;
{
int ii;
/*fprintf(stderr,"R ");*/
for (ii=0;ii<2;ii++){
	if (node[i][ii]<ntaxa) *(groups[k]+node[i][ii]-1)=1;
	else if (node[i][ii]>ntaxa) recurse_groups(node,ntaxa,nent,node[i][ii]-ntaxa,k,groups);
	}

}
