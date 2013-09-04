/*
Copyright 1994 Ward Wheeler all rights reserved
*/

#include "align3.h"

int do_swap_new(nodes,ntax,nbases,n_gaps,initial_length,sequence,values,mult)
int **nodes, ntax,nbases,n_gaps,initial_length,mult;
int **sequence;
parameters *values;
{
int i,j,k,l,m,num_clipped,check,tot_best,best,d1,d2;
int l_clip, l_rest,cur_best,found_better,node_removed;
int *anc, **cur_nodes, *blocked, *is_it_or_desc;
int **hold_nodes,**temp_nodes;
int clipped_base,***reroot_array,cur_num_trees;
int old_length, tree_counter;
int ***best_nodes_new,***best_nodes_old,num_found;
int new_start,all_others;
int where_it_was_anc,where_it_was_desc;
int *td1,*td2,*td3;
int **seq_buffer,*thing_to_add,*reopt_array;
int *pablo_m,*td4,*td5,*td6,*td7,*td8;
int **check_seq_buffer,**temp_nodes_h,whole_length;

/*allocations*/
temp_nodes_h=(int **)malloc((ntax-1)*sizeof(int *));
assert((int)temp_nodes_h);
for (j=0;j<ntax-1;j++) {
	temp_nodes_h[j]=(int *)malloc(2*sizeof(int));
	assert((int)temp_nodes_h);
	}
if (values->shortcut_tree) {
	pablo_m=(int *)malloc((nbases+n_gaps)*sizeof(int));
	assert((int)pablo_m);
	}
best_nodes_new=(int ***)malloc(values->keep_trees*sizeof(int **));
assert((int) best_nodes_new);
best_nodes_old=(int ***)malloc(values->keep_trees*sizeof(int **));
assert((int) best_nodes_old);
for (i=0;i<values->keep_trees;i++) {
	best_nodes_new[i]=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int) best_nodes_new[i]);
	best_nodes_old[i]=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int) best_nodes_old[i]);
	for (j=0;j<ntax-1;j++) {
		best_nodes_new[i][j]=(int *)malloc(2*sizeof(int));
		assert((int) best_nodes_new[i][j]);
		best_nodes_old[i][j]=(int *)malloc(2*sizeof(int));
		assert((int) best_nodes_old[i][j]);
		}
	}       


hold_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) hold_nodes);
cur_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) cur_nodes);
temp_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) temp_nodes);
for (k=0;k<ntax-1;k++) {
	cur_nodes[k]=(int *)malloc(2*sizeof(int));
	assert((int) cur_nodes[k]);
	hold_nodes[k]=(int *)malloc(2*sizeof(int));
	assert((int) hold_nodes[k]);
	temp_nodes[k]=(int *)malloc(2*sizeof(int));
	assert((int) temp_nodes[k]);
	}
blocked=(int *)malloc((2*(ntax-1))*sizeof(int));
assert((int) blocked);
anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);
is_it_or_desc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) is_it_or_desc);
if (values->tbr) {
	reroot_array=(int ***)malloc((2*ntax-3)*sizeof(int **));
	assert((int) reroot_array);
	for (i=0;i<2*ntax-3;i++) {
		reroot_array[i]=(int **)malloc((ntax-1)*sizeof(int *));
		assert((int) reroot_array[i]);
		for (j=0;j<ntax-1;j++) {
			reroot_array[i][j]=(int *)malloc(2*sizeof(int));
			assert((int) reroot_array[i][j]);
			}
		}
	}
reopt_array=(int *)malloc((nbases+n_gaps)*sizeof(int));
assert((int) reopt_array);
thing_to_add=(int *)malloc((nbases+n_gaps)*sizeof(int));
assert((int) thing_to_add);

/*allocate sequence buffers*/
seq_buffer=(int **)malloc(((2*ntax)-1)*sizeof(int *));
assert((int) seq_buffer);
for (k=0;k<(2*ntax)-1;k++) {
	seq_buffer[k]=(int *)malloc((nbases+n_gaps)*sizeof(int));
	assert((int) seq_buffer[k]);
	}
check_seq_buffer=(int **)malloc(((2*ntax)-1)*sizeof(int *));
assert((int) check_seq_buffer);
for (k=0;k<(2*ntax)-1;k++) {
	check_seq_buffer[k]=(int *)malloc((nbases+n_gaps)*sizeof(int));
	assert((int) check_seq_buffer[k]);
	}
for (k=0;k<ntax;k++) for (i=0;i<nbases+n_gaps;i++) seq_buffer[k][i]=check_seq_buffer[k][i]=sequence[k][i];

for (k=0;k<ntax-1;k++) {
	best_nodes_new[0][k][0]=best_nodes_old[0][k][0]=nodes[k][0];
	best_nodes_new[0][k][1]=best_nodes_old[0][k][1]=nodes[k][1];
	}
new_start=num_found=0;
tot_best=old_length=initial_length;
cur_num_trees=1;
found_better=1;
tree_counter=1;
while (found_better) {
	found_better=0;
	old_length=tot_best;
	for (m=new_start;m<cur_num_trees;m++) {
		/*copy nodes*/
		for (k=0;k<ntax-1;k++) {
			hold_nodes[k][0]=best_nodes_old[m][k][0];
			hold_nodes[k][1]=best_nodes_old[m][k][1];
			}
		for (k=0;k<ntax-1;k++) anc[hold_nodes[k][0]]=anc[hold_nodes[k][1]]=k+ntax;     
		for (k=0;k<nbases+n_gaps;k++) reopt_array[k]=1;
		check=do_optimize_nodes(values,seq_buffer,hold_nodes,ntax,ntax,nbases,n_gaps,anc,reopt_array);/*need to save down pass for nodes also*/
		for (i=2;i<2*(ntax-1);i++) {
			for (k=ntax;k<2*ntax-1;k++) for (l=0;l<nbases+n_gaps;l++) if (reopt_array[l]) sequence[k][l]=seq_buffer[k][l];
			/*copy nodes*/
			for (k=0;k<ntax-1;k++) {
				cur_nodes[k][0]=hold_nodes[k][0];
				cur_nodes[k][1]=hold_nodes[k][1];
				anc[hold_nodes[k][0]]=anc[hold_nodes[k][1]]=k+ntax;     
				}
			anc[ntax]=0;
			/*choose*/
			d1=i/2;
			d2=i%2;
			/*get size*/
			num_clipped=get_size(ntax,cur_nodes,d1,d2,blocked,i,is_it_or_desc,anc);
			if (num_clipped<=(ntax-3)) {
				/*get reopt stuff*/
				if (values->shortcut_tree) {
					if (num_clipped>1) { 
						td1=pablo_m;
						td2=sequence[cur_nodes[cur_nodes[d1][d2]-ntax][0]];
						td3=sequence[cur_nodes[cur_nodes[d1][d2]-ntax][1]];
						for (k=0;k<nbases+n_gaps;k++) {
							td1[k]=td2[k]&td3[k];
							if (!td1[k]) td1[k]=td2[k]|td3[k];
							}
						
						for (k=0;k<nbases+n_gaps;k++) {
							if (pablo_m[k]!=sequence[node_removed][k]) reopt_array[k]=1;
							else reopt_array[k]=0;
							}
						}
					else for (k=0;k<nbases+n_gaps;k++) reopt_array[k]=0;
					
					td1=sequence[cur_nodes[d1][d2]];
					if (cur_nodes[d1][!d2]<ntax) {
						td3=sequence[cur_nodes[d1][!d2]];
						for (k=0;k<nbases+n_gaps;k++) if (!reopt_array[k]) {
							if (td1[k]!=td3[k]) reopt_array[k]=1;
							}
						}
					else {
						td3=sequence[cur_nodes[cur_nodes[d1][!d2]-ntax][0]];
						td4=sequence[cur_nodes[cur_nodes[d1][!d2]-ntax][1]];
						for (k=0;k<nbases+n_gaps;k++) if (!reopt_array[k]) {
							if (td1[k]!=td3[k]) reopt_array[k]=1;
							else if (td1[k]!=td4[k]) reopt_array[k]=1;
							}
						}
					
					td2=sequence[d1+ntax];
					td5=sequence[anc[d1+ntax]];
					for (k=0;k<nbases+n_gaps;k++) if (!reopt_array[k]) if (td2[k]!=td5[k]) reopt_array[k]=1;
					}/*shortcut end*/
				if (values->shortcut_tree2) {
					td1=sequence[anc[d1+ntax]];
					td2=sequence[cur_nodes[d1][!d2]];
					for (k=0;k<nbases+n_gaps;k++) {
						if (td1[k]!=td2[k]) reopt_array[k]=1;
						else reopt_array[k]=0;
						}
					}
				/*redo nodes*/
				node_removed=anc[cur_nodes[d1][d2]];
				where_it_was_anc=anc[node_removed];
				if (cur_nodes[anc[node_removed]-ntax][0]==node_removed) {
					cur_nodes[anc[node_removed]-ntax][0]=cur_nodes[d1][!d2];
					blocked[2*(anc[node_removed]-ntax)]=1;
					where_it_was_desc=cur_nodes[anc[node_removed]-ntax][0];
					}
				else {
					cur_nodes[anc[node_removed]-ntax][1]=cur_nodes[d1][!d2];
					blocked[2*(anc[node_removed]-ntax)+1]=1;
					where_it_was_desc=cur_nodes[anc[node_removed]-ntax][1];
					}
				anc[cur_nodes[d1][!d2]]=anc[node_removed];
				cur_nodes[d1][!d2]=cur_nodes[d1][d2];
				anc[node_removed]=node_removed;
				l_rest=optimize_non_blocked2(ntax,cur_nodes,sequence,nbases,n_gaps,values,anc,is_it_or_desc,reopt_array);
				/*fprintf(stderr,"got lrest ");*/
				/*when one or two taxa*/
				if ((num_clipped < 3) || (!values->tbr)) { /*if no alternative roots or sbr*/
					/* get "best" set sequences[cur_nodes[d1][d2]] as XUY */
					if (num_clipped>1) {
						td2=seq_buffer[cur_nodes[cur_nodes[d1][d2]-ntax][0]];
						td3=seq_buffer[cur_nodes[cur_nodes[d1][d2]-ntax][1]];
						for (k=0;k<nbases+n_gaps;k++) {
							thing_to_add[k]=(td2[k]&td3[k]);
							if (!thing_to_add[k]) thing_to_add[k]=(td2[k]|td3[k]);
							}
						}
					else {
						td2=seq_buffer[cur_nodes[d1][d2]];
						for (k=0;k<nbases+n_gaps;k++) thing_to_add[k]=td2[k];
						}
					/*if (where_it_was_desc<ntax) best=get_interval(thing_to_add,sequence[where_it_was_anc],sequence[where_it_was_desc],values,nbases,n_gaps,HUGE_COST);
					else */
					best=get_interval2(thing_to_add,sequence[where_it_was_anc],sequence[where_it_was_desc],values,nbases,n_gaps,HUGE_COST);
					all_others=old_length-best;
					/*fprintf(stderr,"\nB %d ",best);*/
					if ((best>0) || ((best==0) && (mult))) {
						/*loop through places to put back*/
						blocked[1]=0;
						for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
							/* get_cost*/
							/*if (cur_nodes[j/2][j%2]<ntax) cur_best=get_interval2(thing_to_add,sequence[cur_nodes[j/2][j%2]],sequence[(j/2)+ntax],values,nbases,n_gaps,best+mult);
							else */
							cur_best=get_interval2(thing_to_add,sequence[cur_nodes[j/2][j%2]],sequence[(j/2)+ntax],values,nbases,n_gaps,best+mult);
							whole_length=cur_best+all_others;
							/*fprintf(stderr,"C %d ",cur_best);*/
							/*the following is a filter for missing data heavy data sets--they screw up shortcuts*/
							if ((whole_length < tot_best) || ((whole_length == tot_best) && mult && (tree_counter< values->keep_trees))) {
								for (k=0;k<ntax-1;k++) {
									temp_nodes_h[k][0]=cur_nodes[k][0];
									temp_nodes_h[k][1]=cur_nodes[k][1];
									}
								check=temp_nodes_h[j/2][j%2];
								temp_nodes_h[j/2][j%2]=node_removed;
								temp_nodes_h[d1][!d2]=check;
								check=check_length(temp_nodes_h,ntax,nbases,check_seq_buffer,ntax,n_gaps,values,tot_best+1);
								/*fprintf(stderr,"(%d vs. %d)",check,whole_length);*/
								if (check>whole_length) whole_length=HUGE_COST;
								}
							if (whole_length < tot_best) {
								/*fprintf(stderr,"FO%d ",whole_length);*/
								best=cur_best;
								tot_best=whole_length;
								found_better=1;
								for (k=0;k<ntax-1;k++) {
									best_nodes_new[0][k][0]=temp_nodes_h[k][0];
									best_nodes_new[0][k][1]=temp_nodes_h[k][1];
									}
								/*check=best_nodes_new[0][j/2][j%2];
								best_nodes_new[0][j/2][j%2]=node_removed;
								best_nodes_new[0][d1][!d2]=check;*/
								tree_counter=1;
								} 
							else if ((whole_length == tot_best) && mult && (tree_counter< values->keep_trees)) {
								for (k=0;k<ntax-1;k++) {
									best_nodes_new[tree_counter][k][0]=temp_nodes_h[k][0];
									best_nodes_new[tree_counter][k][1]=temp_nodes_h[k][1];
									}
								/*check=best_nodes_new[tree_counter][j/2][j%2];
								best_nodes_new[tree_counter][j/2][j%2]=node_removed;
								best_nodes_new[tree_counter][d1][!d2]=check;
								fprintf(stderr,"(%d)",check);*/
								if (is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax)) {
									/*fprintf(stderr,"FA%d ",tree_counter);*/
									found_better=1;
									++tree_counter;
									}
								}
							} /*SBR loop--regrafted*/
						} /*BEST split size ok*/
					}/*too few clipped or SBR*/
				else { /*TBR with at least three taxa*/
					check=get_root_array(reroot_array,num_clipped,node_removed,cur_nodes,d1,d2,ntax,is_it_or_desc,anc);
					for (l=0;l<2*num_clipped-3;l++) {/*TBR loop*/
						td1=thing_to_add;
						td2=seq_buffer[reroot_array[l][reroot_array[l][d1][d2]-ntax][0]];
						td3=seq_buffer[reroot_array[l][reroot_array[l][d1][d2]-ntax][1]];
						for (k=0;k<nbases+n_gaps;k++) {
							td1[k]=(td2[k]&td3[k]);
							if (!td1[k]) td1[k]=(td2[k]|td3[k]);
							}
						/*if (where_it_was_desc<ntax) best=get_interval(thing_to_add,sequence[where_it_was_anc],sequence[where_it_was_desc],values,nbases,n_gaps,HUGE_COST);
						else */
						best=get_interval2(thing_to_add,sequence[where_it_was_anc],sequence[where_it_was_desc],values,nbases,n_gaps,HUGE_COST);
						all_others=old_length-best;
						if ((best>0) && ((best==0) && (mult))) {
							/*loop through places to put back*/
							blocked[1]=0;
							for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
								/* get_cost*/
								/*if (cur_nodes[j/2][j%2]<ntax) cur_best=get_interval(thing_to_add,sequence[reroot_array[l][j/2][j%2]],sequence[(j/2)+ntax],values,nbases,n_gaps,best+mult);
								else*/
								cur_best=get_interval2(thing_to_add,sequence[reroot_array[l][j/2][j%2]],sequence[(j/2)+ntax],values,nbases,n_gaps,best+mult);	
								whole_length=cur_best+all_others;
								/*put in check here later copy nodes if OK*/
								if ((whole_length < tot_best) || ((whole_length == tot_best) && mult && (tree_counter< values->keep_trees))) {
									for (k=0;k<ntax-1;k++) {
										temp_nodes_h[k][0]=reroot_array[l][k][0];
										temp_nodes_h[k][1]=reroot_array[l][k][1];
										}
									check=temp_nodes_h[j/2][j%2];
									temp_nodes_h[j/2][j%2]=node_removed;
									temp_nodes_h[d1][!d2]=check;
									check=check_length(temp_nodes_h,ntax,nbases,check_seq_buffer,ntax,n_gaps,values,tot_best+1);
									/*fprintf(stderr,"(%d vs %d)",check,whole_length);*/
									if (check>whole_length) whole_length=HUGE_COST;
									}
								if (whole_length < tot_best) {
									best=cur_best;
									tot_best=whole_length;
									found_better=1;
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[0][k][0]=temp_nodes_h[k][0];
										best_nodes_new[0][k][1]=temp_nodes_h[k][1];
										}
									/*check=best_nodes_new[0][j/2][j%2];
									best_nodes_new[0][j/2][j%2]=node_removed;
									best_nodes_new[0][d1][!d2]=check;*/
									tree_counter=1;
									} 
								else if ((whole_length == tot_best) && mult && (tree_counter< values->keep_trees)) {
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[tree_counter][k][0]=temp_nodes_h[k][0];
										best_nodes_new[tree_counter][k][1]=temp_nodes_h[k][1];
										}
									/*check=best_nodes_new[tree_counter][j/2][j%2];
									best_nodes_new[tree_counter][j/2][j%2]=node_removed;
									best_nodes_new[tree_counter][d1][!d2]=check;*/
									if (is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax)) {
										found_better=1;
										++tree_counter;
										/*fprintf(stderr,"FA%d ",tree_counter);*/
										}
									} 
								} /*SBR loop--regrafted*/
							} /*BEST split size ok*/
						} /*TBR loop*/
					} /*reroot TBR/SBR choice*/
				} /*size of plucked OK*/
			} /*i clades to pick*/
		/*copy after best is found*/
		}/*num trees*/
	/*copy over new set of best_nodes*/
	if (tree_counter>0) {
		if (tot_best < old_length) { /*new trees are better*/
			num_found=tree_counter;
			new_start=0;
			old_length=tot_best;
			for (k=0;k<tree_counter;k++) {
				for (l=0;l<ntax-1;l++) {
					best_nodes_old[k][l][0]=best_nodes_new[k][l][0];
					best_nodes_old[k][l][1]=best_nodes_new[k][l][1];
					}
				}
			} /*from zero*/
		else {          /*New trees are same length as old from cur_num_trees*/
			num_found=tree_counter;
			new_start=cur_num_trees;
			for (k=new_start;k<tree_counter;k++) {
				for (l=0;l<ntax-1;l++) {
					best_nodes_old[k][l][0]=best_nodes_new[k][l][0];
					best_nodes_old[k][l][1]=best_nodes_new[k][l][1];
					}
				}
			}
		}
	/*reset num trees*/
	cur_num_trees=tree_counter;
	}/*found better reswap if found a better one*/
for (k=0;k<ntax-1;k++) {
	nodes[k][0]=best_nodes_new[0][k][0];
	nodes[k][1]=best_nodes_new[0][k][1];
	} 

for (i=0;i<ntax-1;i++) {
	free(cur_nodes[i]);
	free(hold_nodes[i]);
	free(temp_nodes[i]);
	free(temp_nodes_h[i]);
	}
free(cur_nodes);
free(hold_nodes);
free(temp_nodes);
free(temp_nodes_h);
free(anc);
free(is_it_or_desc);
free(blocked);

if (values->tbr) {
	for (i=0;i<2*ntax-3;i++) {
		for (j=0;j<ntax-1;j++) free(reroot_array[i][j]);
		free(reroot_array[i]);
		}
	free(reroot_array);
	}

for (i=0;i<values->keep_trees;i++) {
	for (j=0;j<ntax-1;j++) {
		free(best_nodes_new[i][j]);
		free(best_nodes_old[i][j]);
		}
	free(best_nodes_new[i]);
	free(best_nodes_old[i]);
	}       
free(best_nodes_new);
free(best_nodes_old);
for (k=0;k<(2*ntax)-1;k++) {
	free(seq_buffer[k]);
	free(check_seq_buffer[k]);
	}
free(seq_buffer);
free(check_seq_buffer);
free(reopt_array);	
free(thing_to_add);
if (values->shortcut_tree) free(pablo_m);

/*fprintf(stderr,"Found %d trees\n",num_found);*/
return tot_best;
}

int get_size(ntax,nodes,d1,d2,blocked,place,desc,anc) 
int ntax,d1,d2,**nodes,*blocked,place,*desc,*anc;
{
int num_in,i,found_one;

num_in=0;
desc[0]=desc[1]=0;
blocked[0]=blocked[1]=1;
for (i=2;i<2*(ntax-1);i++) blocked[i]=desc[i]=0;
desc[2*ntax-2]=0;

 /*block the node chosen and its sister so not added later to non-existant branch*/
blocked[place]=1;
if (d2) blocked[place-1]=1;
else blocked[place+1]=1;

found_one=1;
/*for (i=2;i<(2*ntax)-2;i++) fprintf(stderr,"B(%d)%d ",i,blocked[i]);*/
/*count and block descendents of choice*/
/*
while (found_one) {
	found_one=0;
	for (i=2;i<2*(ntax-1);i++) {
		if ((blocked[i]) && (!blocked[nodes[i/2][0]]) && (!blocked[nodes[i/2][1]])) {
			blocked[nodes[i/2][0]]=blocked[nodes[i/2][1]]=1;
			found_one=1;
			}
		}                       
	}*/
found_one=1;
desc[nodes[d1][d2]]=1;
desc[anc[nodes[d1][d2]]]=1;
while (found_one) {
	found_one=0;
	for (i=0;i<ntax-1;i++) if ((desc[i+ntax]) && (!desc[nodes[i][0]]) && (!desc[nodes[i][1]])) {
		desc[nodes[i][0]]=desc[nodes[i][1]]=1;
		found_one=1;
		}
	}
for (i=1;i<ntax-1;i++) if (desc[i+ntax]) {
	if (desc[nodes[i][0]]) blocked[2*i]=1;
	if (desc[nodes[i][1]]) blocked[(2*i)+1]=1;
	}
for (i=1;i<ntax;i++) if (desc[i]) ++num_in;
return num_in;
}

int optimize_blocked2(ntax,nodes,sequence,basal_node,nbases,n_gaps,values,desc,reopt)
int ntax, **nodes, **sequence,nbases,n_gaps,basal_node,*desc,*reopt;
parameters *values;
{
int *done,found_one,i,j;
int *td1,*td2,*td3,length,ch,ti,tv,gap;

length=0;

done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[1]=1;
for (i=ntax;i<(2*ntax)-1;i++) {
	if (desc[i]) done[i]=0;
	else done[i]=1;
	}
found_one=1;
while (found_one) {
	found_one=0;
	for (i=0;i<ntax-1;i++) {
		if ((!done[ntax+i]) && (done[nodes[i][0]]) && (done[nodes[i][1]])) {
			done[ntax+i]=1;
			found_one=1;
			td3=sequence[ntax+i];
			td2=sequence[nodes[i][0]];
			td1=sequence[nodes[i][1]];
			if (nodes[i][0]==nodes[i][1]) for (j=0;j<nbases+n_gaps;j++) td3[j]=td2[j];
			else {
				 for (j=0;j<(nbases+n_gaps);j++) {
					td3[j]=(td1[j]&td2[j]);
					if (!td3[j]) { td3[j]=(td1[j]|td2[j]);}
					}                       
				}
			}
		}
	}
free(done);
return length;
}

int optimize_blocked(ntax,nodes,sequence,basal_node,nbases,n_gaps,values,desc)
int ntax, **nodes, **sequence,nbases,n_gaps,basal_node,*desc;
parameters *values;
{
int *done,found_one,i,j;
int *td1,*td2,*td3,length,ch,ti,tv,gap;

length=0;

ch=values->change_cost;
ti=values->transition;
tv=values->transversion;
gap=values->gap_cost;

done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[1]=1;
for (i=ntax;i<(2*ntax)-1;i++) {
	if (desc[i]) done[i]=0;
	else done[i]=1;
	}
found_one=1;
while (found_one) {
	found_one=0;
	for (i=0;i<ntax-1;i++) {
		if ((!done[ntax+i]) && (done[nodes[i][0]]) && (done[nodes[i][1]])) {
			done[ntax+i]=1;
			found_one=1;
			td3=sequence[ntax+i];
			td2=sequence[nodes[i][0]];
			td1=sequence[nodes[i][1]];
			if (nodes[i][0]==nodes[i][1]) for (j=0;j<nbases+n_gaps;j++) td3[j]=td2[j];
			else {
				if (!values->ttr) {
					for (j=0;j<nbases;j++) {
						td3[j]=(td1[j]&td2[j]);
						if (!td3[j]) { td3[j]=(td1[j]|td2[j]); ++ch;}
						}
					}
				else {
					for (j=0;j<nbases-values->n_transv;j++) {
						td3[j]=(td1[j]&td2[j]);
						if (!td3[j]) { td3[j]=(td1[j]|td2[j]); ++ti;}
						}
					for (j=nbases-values->n_transv;j<nbases;j++) {
						td3[j]=(td1[j]&td2[j]);
						if (!td3[j]) { td3[j]=(td1[j]|td2[j]); ++tv;}
						}            
					}
				 for (j=nbases;j<(nbases+n_gaps);j++) {
					td3[j]=(td1[j]&td2[j]);
					if (!td3[j]) { td3[j]=(td1[j]|td2[j]); ++gap;}
					}                       
				}
			}
		}
	}
free(done);
return length;
}

int optimize_non_blocked(ntax,nodes,taxa,nbases,n_gaps,values,anc,desc)
int ntax, **nodes,**taxa,nbases,n_gaps,*anc,*desc;
parameters *values;
{
int *td1,*td2,*td3,*td4,*td5;
int n,m,d1,d2;
int ***up_nodes,i_am,i, length=0;
int ch,ti,tv,gap,found_one,*reopt;


ch=values->change_cost;
ti=values->transition;
tv=values->transversion;
gap=values->gap_cost;

/*remove this when do reopt shortcut*/
reopt=(int *)malloc((nbases+n_gaps)*sizeof(int));
assert((int) reopt);
for (i=0;i<(nbases+n_gaps);i++) reopt[i]=1;

up_nodes=(int ***)malloc(((2*ntax)-1)*sizeof(int**));
assert((int) up_nodes);
for (i=0;i<(2*ntax)-1;i++) {
	up_nodes[i]=(int **)malloc((3)*sizeof(int *));
	assert((int) up_nodes[i]);
	up_nodes[i][0]=(int *)malloc(((nbases+n_gaps))*sizeof(int)); /*holds down pass*/
	assert((int) up_nodes[i][0]);
	up_nodes[i][1]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*holds up pass to d1*/
	assert((int) up_nodes[i][1]);
	up_nodes[i][2]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*hold uppass to d2*/
	assert((int) up_nodes[i][2]);
	}

/*down pass*/
for (n=0;n<ntax;n++) {
	for (i=0;i<(nbases+n_gaps);i++) up_nodes[n][0][i]=taxa[n][i];
	values->tree_made[n]=1;
	}
/*add tree_made for values->tree_made things*/
for (i=ntax;i<(2*ntax)-1;i++) {
	if (desc[i]) values->tree_made[i]=1;
	else values->tree_made[i]=0;
	}
found_one=1;
while (found_one) {
	found_one=0;
	for (n=(ntax-2);n>=0;--n){
		d1=nodes[n][0];
		d2=nodes[n][1];
		if ((values->tree_made[d1]) && (values->tree_made[d2]) && (!values->tree_made[ntax+n])) {
			found_one=1;
			td1=up_nodes[d1][0];
			td2=up_nodes[d2][0];
			td3=up_nodes[ntax+n][0];
			values->tree_made[ntax+n]=1;
			if (!values->ttr) {
				for (m=0;m<nbases;m++) /* if (reopt[m]) */{
					td3[m]=(td1[m]&td2[m]);
					if (!td3[m]) { td3[m]=(td1[m]|td2[m]); ++ch;}
					}
				}
			else {
				for (m=0;m<nbases-values->n_transv;m++) /* if (reopt[m]) */{
					td3[m]=(td1[m]&td2[m]);
					if (!td3[m]) { td3[m]=(td1[m]|td2[m]); ++ti;}
					}
				for (m=nbases-values->n_transv;m<nbases;m++) /* if (reopt[m]) */{
					td3[m]=(td1[m]&td2[m]);
					if (!td3[m]) { td3[m]=(td1[m]|td2[m]); ++tv;}
					}            
				}
			 for (m=nbases;m<(nbases+n_gaps);m++) /* if (reopt[m]) */{
				td3[m]=(td1[m]&td2[m]);
				if (!td3[m]) { td3[m]=(td1[m]|td2[m]); ++gap;}
				}
			/*fprintf(stderr,"Length from %d + %d -> %d is %d\n",d1,d2,ntax+n,length);*/
			}
		}
	}

/*up pass*/     
/*remember to remove extraneous up pass (ie when leads to a terminal*/
td1=up_nodes[ntax][0];
td2=up_nodes[ntax][1];
td3=up_nodes[ntax][2];
for (m=0;m<(nbases+n_gaps);m++) /* if (reopt[m]) */td1[m]=td2[m]=td3[m]=taxa[ntax][m]=taxa[0][m]; /*set basal node = to outgroup*/

for (i=ntax+1;i<(2*ntax)-1;i++) {
	if (desc[i]) values->tree_made[i]=1;
	else values->tree_made[i]=0;
	}
values->tree_made[ntax]=1;
found_one=1;
while (found_one) {
	found_one=0;
	for (n=1;n<ntax-1;n++){
		if (!values->tree_made[ntax+n]) if (values->tree_made[anc[ntax+n]]) {
			found_one=1;
			values->tree_made[ntax+n]=1;
			if (nodes[anc[ntax+n]-ntax][0]==(n+ntax)) i_am=1;
			else i_am=2;
			td1=up_nodes[n+ntax][1];
			td2=up_nodes[n+ntax][2];
			td3=up_nodes[anc[ntax+n]][i_am];
			td4=up_nodes[nodes[n][1]][0];
			td5=up_nodes[nodes[n][0]][0];
			i=get_best_states_opt(taxa[ntax+n],td3,up_nodes[ntax+n][0],td5,td4,nbases+n_gaps,reopt); /*add reopt as an option after fix*/
			/*get up pass states for this node*/
			if (nodes[n][0]>(ntax-1)) {
				for (m=0;m<(nbases+n_gaps);m++)  /* if (reopt[m]) */{
					td1[m]=(td4[m]&td3[m]);
					if (!td1[m]) td1[m]=(td4[m]|td3[m]);
					}
				}
			if (nodes[n][1]>(ntax-1)) {
				for (m=0;m<(nbases+n_gaps);m++)  /* if (reopt[m]) */{
					td2[m]=(td5[m]&td3[m]);
					if (!td2[m]) td2[m]=(td5[m]|td3[m]);
					}
				}
			}
		}
	}
for (i=0;i<(2*ntax)-1;i++) {
	free(up_nodes[i][0]);
	free(up_nodes[i][1]);
	free(up_nodes[i][2]);
	free(up_nodes[i]);
	}
free(up_nodes);

/*remove when shotcut*/
free(reopt);

return length;
}


int check_length(nodes,ntaxa,nbases,taxa,nent,n_gaps,values,bound)
int **nodes,bound;
int ntaxa,nbases,nent;
int **taxa,n_gaps;
parameters *values;
{
	int ii,dummy1,dummy2;
	int nn,d1,d2,all_nodes_done;
	int change_counts,mm,ch,ti,tv,gap,found_one;
	register int *td1, *td2, *td3;
	int length;
	
length=0;

	for (nn=0;nn<ntaxa;nn++) values->tree_made[nn]=1;
	for (nn=ntaxa;nn<(ntaxa+nent-1);++nn) values->tree_made[nn]=0;
	ch=values->change_cost;
	ti=values->transition;
	tv=values->transversion;
	gap=values->gap_cost;  
	
found_one=1;
while (found_one) {
	found_one=0;
	for (nn=(nent-2);nn>=0;--nn){
		d1=nodes[nn][0];
		d2=nodes[nn][1];
		if ((values->tree_made[d1]) && (values->tree_made[d2]) && (!values->tree_made[ntaxa+nn])) {
			found_one=1;
			td1=taxa[d1];
			td2=taxa[d2];
			td3=taxa[ntaxa+nn];
			values->tree_made[ntaxa+nn]=1;
			if (!values->ttr) {
				for (mm=0;mm<nbases;mm++) {
					td3[mm]=(td1[mm]&td2[mm]);
					if (!td3[mm]) {
						td3[mm]=(td1[mm]|td2[mm]);
						length+=ch;
						if (length >= bound) return HUGE_COST;
					}
				}
			}
			else {
				for (mm=0;mm<nbases-values->n_transv;mm++) {
					td3[mm]=(td1[mm]&td2[mm]);
					if (!td3[mm]) {
						td3[mm]=(td1[mm]|td2[mm]);
						length+=ti;
						if (length >= bound) return HUGE_COST;
					}
				}
				for (mm=nbases-values->n_transv;mm<nbases;mm++) {
					td3[mm]=(td1[mm]&td2[mm]);
					if (!td3[mm]) {
						td3[mm]=(td1[mm]|td2[mm]);
						length+=tv;
						if (length >= bound) return HUGE_COST;
					}
				}
			}
			for (mm=nbases;mm<(nbases+n_gaps);mm++) {
				td3[mm]=(td1[mm]&td2[mm]);
				if (!td3[mm]) {
					td3[mm]=(td1[mm]|td2[mm]);
					length+=gap;
					if (length >= bound) return HUGE_COST;
					}
				}
			}
		}
	}
return length;
}

int do_simul_optimize_nodes(values,taxa,nodes,ntaxa,nent,nbases,n_gaps,anc,lclip,lrest,desc)
int **nodes, **taxa, ntaxa, nent,*anc,*desc;
int *lclip,*lrest;
parameters *values;
{
int *td1,*td2,*td3,*td4,*td5;
int n,m,d1,d2;
int ch,gap,ti,tv;
int ***up_nodes,i_am,i;
int found_one;


up_nodes=(int ***)malloc(((2*ntaxa)-1)*sizeof(int**));
assert((int) up_nodes);
for (i=0;i<(2*ntaxa)-1;i++) {
	up_nodes[i]=(int **)malloc((3)*sizeof(int *));
	assert((int) up_nodes[i]);
	up_nodes[i][0]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*holds down pass*/
	assert((int) up_nodes[i][0]);
	up_nodes[i][1]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*holds up pass to d1*/
	assert((int) up_nodes[i][1]);
	up_nodes[i][2]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*hold uppass to d2*/
	assert((int) up_nodes[i][2]);
	}
	
ch=values->change_cost;
gap=values->gap_cost;
ti=values->transition;
tv=values->transversion;

*lclip=0;
*lrest=0;
/*down pass*/
for (n=0;n<ntaxa;n++) {
	for (i=0;i<nbases+n_gaps;i++) up_nodes[n][0][i]=taxa[n][i];
	values->tree_made[n]=1;
	}
for (n=ntaxa;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;
found_one=1;
while (found_one) {
	found_one=0;
	for (n=(nent-2);n>=0;--n){
		d1=nodes[n][0];
		d2=nodes[n][1];
		if ((values->tree_made[d1]) && (values->tree_made[d2]) && (!values->tree_made[ntaxa+n])) {
			found_one=1;
			td1=up_nodes[d1][0];
			td2=up_nodes[d2][0];
			td3=up_nodes[ntaxa+n][0];
			/*fprintf(stderr,"values->tree_made %d ",ntaxa+n);*/
			values->tree_made[ntaxa+n]=1;
			/*only need this to check*/
			if (desc[n+ntaxa]) { /*in clipped*/
				if (!values->ttr) {
					for (m=0;m<nbases;m++) {
						if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
						else {
							td3[m]=(td1[m]|td2[m]);
						(*lclip)+=ch;
						}
						}
					}
				else {
				  for (m=0;m<nbases-values->n_transv;m++) {
				    if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
				    else {
				      td3[m]=(td1[m]|td2[m]);
				      (*lclip)+=ti;
				 }
					}
				  for (m=nbases-values->n_transv;m<nbases;m++) {
				      if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
					else {
					 td3[m]=(td1[m]|td2[m]);
					 (*lclip)+=tv;
					 }
				       }
				     }            
				 for (m=nbases;m<(nbases+n_gaps);m++) {
				  if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
				  else {
					td3[m]=(td1[m]|td2[m]);
					(*lclip)+=gap;
					}
				 }      
			}/*clipped*/
			else { /*in rest*/
				if (!values->ttr) {
					for (m=0;m<nbases;m++) {
						if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
						else {
							td3[m]=(td1[m]|td2[m]);
						(*lrest)+=ch;
						}
						}
					}
				else {
				  for (m=0;m<nbases-values->n_transv;m++) {
				    if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
				    else {
				      td3[m]=(td1[m]|td2[m]);
				      (*lrest)+=ti;
				 }
					}
				  for (m=nbases-values->n_transv;m<nbases;m++) {
				      if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
					else {
					 td3[m]=(td1[m]|td2[m]);
					 (*lrest)+=tv;
					 }
				       }
				     }            
				 for (m=nbases;m<(nbases+n_gaps);m++) {
				  if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
				  else {
					td3[m]=(td1[m]|td2[m]);
					(*lrest)+=gap;
					}
				 }      
			}/*rest*/
		}
	}
}

/*up pass*/     
/*remember to remove extraneous up pass (ie when leads to a terminal*/
td1=up_nodes[ntaxa][0];
td2=up_nodes[ntaxa][1];
td3=up_nodes[ntaxa][2];
for (m=0;m<nbases+n_gaps;m++) td1[m]=td2[m]=td3[m]=taxa[ntaxa][m]=taxa[0][m]; /*set basal node = to outgroup*/

for (i=ntaxa+1;i<(2*ntaxa)-1;i++) values->tree_made[i]=0;
values->tree_made[ntaxa]=1;
found_one=1;
while (found_one) {
	found_one=0;
for (n=1;n<nent-1;n++){
	if ((values->tree_made[anc[ntaxa+n]]) && (!values->tree_made[ntaxa+n])) {
		found_one=1;
		values->tree_made[ntaxa+n]=1;
		if (nodes[anc[ntaxa+n]-ntaxa][0]==(n+ntaxa)) i_am=1;
		else i_am=2;
		td1=up_nodes[n+ntaxa][1];
		td2=up_nodes[n+ntaxa][2];
		td3=up_nodes[anc[ntaxa+n]][i_am];
		td4=up_nodes[nodes[n][1]][0];
		td5=up_nodes[nodes[n][0]][0];
		i=get_best_states_opt(taxa[ntaxa+n],td3,up_nodes[ntaxa+n][0],td5,td4,nbases+n_gaps);
		/*get up pass states for this node*/
		if (nodes[n][0]>(ntaxa-1)) {
			for (m=0;m<nbases+n_gaps;m++) {
				td1[m]=(td4[m]&td3[m]);
				if (!td1[m]) td1[m]=(td4[m]|td3[m]);
				}
			}
		if (nodes[n][1]>(ntaxa-1)) {
			for (m=0;m<nbases+n_gaps;m++) {
				td2[m]=(td5[m]&td3[m]);
				if (!td2[m]) td2[m]=(td5[m]|td3[m]);
				}
			}
			
		}
	}
}

for (i=0;i<(2*ntaxa)-1;i++) {
	free(up_nodes[i][0]);
	free(up_nodes[i][1]);
	free(up_nodes[i][2]);
	free(up_nodes[i]);
	}
free(up_nodes);
return ((*lclip)+(*lrest));
}

/*remember only deal with desc[] and eep abse same name*/
int get_root_array(root_array,num_clipped,node_removed,cur_nodes,d1,d2,ntax,desc,anc)
int ***root_array,num_clipped,node_removed,**cur_nodes,d1,d2,ntax,*desc;
int *anc;
{
int **t_reroot_nodes,**reroot_nodes,b1,b2;
int i,j,k,l,m,n,clipped_base;
int **block_nodes,rerootings,did_one,*rooted;


/*for (i=0;i<ntax-1;i++) {
	for (j=0;j<2;j++) fprintf(stderr,"cnb[%d][%d]=%d ",i,j,cur_nodes[i][j]);
	fprintf(stderr,"\n");
	}
*/

rooted=(int *)malloc((ntax-1)*sizeof(int));
assert((int)rooted);
clipped_base=cur_nodes[d1][d2];
block_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) block_nodes);
reroot_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) reroot_nodes);
t_reroot_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) t_reroot_nodes);
for (i=0;i<ntax-1;i++) {
	block_nodes[i]=(int *)malloc(3*sizeof(int));
	assert((int) block_nodes[i]);
	reroot_nodes[i]=(int *)malloc(3*sizeof(int));
	assert((int) reroot_nodes[i]);
	t_reroot_nodes[i]=(int *)malloc(3*sizeof(int));
	assert((int) t_reroot_nodes[i]);
	reroot_nodes[i][0]=cur_nodes[i][0]; 
	reroot_nodes[i][1]=cur_nodes[i][1]; 
	if (i>0) reroot_nodes[i][2]=anc[i+ntax];	
	}

/*deroot tree*/
if (reroot_nodes[clipped_base-ntax][0]>=ntax) {
	b1=reroot_nodes[clipped_base-ntax][0];
	b2=reroot_nodes[clipped_base-ntax][1];
	}
else  {
	b1=reroot_nodes[clipped_base-ntax][1];
	b2=reroot_nodes[clipped_base-ntax][0];
	}
reroot_nodes[b1-ntax][2]=b2;
if (b2>=ntax) reroot_nodes[b2-ntax][2]=b1; 

/*get blocked and redundancies*/
for (i=0;i<ntax-1;i++) {
	for (j=0;j<3;j++) block_nodes[i][j]=0;
	if (!desc[i+ntax]) {for (j=0;j<3;j++) block_nodes[i][j]=1;}
	else if ((i==(node_removed-ntax)) || (i==(clipped_base-ntax))) {for (j=0;j<3;j++) block_nodes[i][j]=1;}
	}

for (i=0;i<ntax-1;i++)  for (j=0;j<3;j++) if (!block_nodes[i][j]) {
	for (k=i+1;k<ntax-1;k++) for (l=0;l<3;l++) if (!block_nodes[k][l]) {
		if ((reroot_nodes[i][j]==(k+ntax)) && (reroot_nodes[k][l]==(i+ntax))) block_nodes[k][l]=1;
		}
	}
for (i=0;i<ntax-1;i++)  for (j=0;j<3;j++) if (!block_nodes[i][j]) {
	for (k=i+1;k<ntax-1;k++) for (l=0;l<3;l++) if (!block_nodes[k][l]) {
		if ((reroot_nodes[i][j]==(k+ntax)) && (reroot_nodes[k][l]==(i+ntax))) block_nodes[k][l]=1;
		}
	}

/*pick new root and reroot*/
rerootings=0;
for (i=0;i<ntax-1;i++) for (j=0;j<3;j++) if (!block_nodes[i][j]) {
	for (k=0;k<ntax-1;k++) for (l=0;l<3;l++) t_reroot_nodes[k][l]=reroot_nodes[k][l];
	/*fprintf(stderr,"Before CB %d NR %d\n",clipped_base,node_removed);
	for (k=0;k<ntax-1;k++) fprintf(stderr,"trn[%d] %d %d %d\n",k,t_reroot_nodes[k][0],t_reroot_nodes[k][1],t_reroot_nodes[k][2]);*/
	/*set new invasion*/
	t_reroot_nodes[clipped_base-ntax][0]=i+ntax;
	t_reroot_nodes[clipped_base-ntax][1]=reroot_nodes[i][j];
	t_reroot_nodes[clipped_base-ntax][2]=node_removed;
	t_reroot_nodes[i][j]=clipped_base;
	if (reroot_nodes[i][j]>ntax) for (l=0;l<3;l++) if (reroot_nodes[reroot_nodes[i][j]-ntax][l]==(i+ntax)) t_reroot_nodes[reroot_nodes[i][j]-ntax][l]=clipped_base;
	/*fprintf(stderr,"After\n");
	for (k=0;k<ntax-1;k++) fprintf(stderr,"trn[%d] %d %d %d\n",k,t_reroot_nodes[k][0],t_reroot_nodes[k][1],t_reroot_nodes[k][2]);*/
	for (l=0;l<ntax-1;l++) {
		if (desc[l+ntax]) rooted[l]=0;
		else rooted[l]=1;
		}
	rooted[clipped_base-ntax]=0;
	rooted[node_removed-ntax]=1;
	did_one=1;
	while (did_one) {
		did_one=0;
		for (k=0;k<ntax-1;k++) if (desc[k+ntax]) if (!rooted[k]) {
			for (l=0;l<ntax-1;l++) if (desc[l+ntax]) if (rooted[l]) {
				for (m=0;m<2;m++) if (t_reroot_nodes[l][m]==(k+ntax)) {
					did_one=1;
					rooted[k]=1;
					for (n=0;n<2;n++) if (t_reroot_nodes[k][n]==(l+ntax)){
						t_reroot_nodes[k][n]=t_reroot_nodes[k][2];
						t_reroot_nodes[k][2]=l+ntax;
						}
					}
				}
			
			}
		}
	
	for (k=0;k<ntax-1;k++) {
		root_array[rerootings][k][0]=t_reroot_nodes[k][0];
		root_array[rerootings][k][1]=t_reroot_nodes[k][1];              
		}
	++rerootings;
	}
if (rerootings!=(2*num_clipped-3)) {fprintf(stderr,"Ouch %d vs. %d rerootings\n",rerootings,2*num_clipped-3);}

free(rooted);
for (i=0;i<ntax-1;i++) {
	free(reroot_nodes[i]);
	free(t_reroot_nodes[i]);
	free(block_nodes[i]);
	}
free(reroot_nodes);
free(t_reroot_nodes);
free(block_nodes);

/*if (num_clipped > 3) exit(-1);*/
/*for (k=0;k<rerootings;k++ ) {
	for (i=0;i<ntax-1;i++) {
		for (j=0;j<2;j++) fprintf(stderr,"cna[%d][%d]=%d ",i,j,root_array[k][i][j]);
		fprintf(stderr,"\n");
		}
	fprintf(stderr,"\n");
	}*/

return 1;
}

int is_unique_tree(cur_nodes,new_nodes,num_new,old_nodes,num_old,old_length,new_length,ntax) 
int **cur_nodes;
int ***new_nodes, ***old_nodes;
int num_new, num_old, new_length, old_length, ntax;
{
int i,j,test;
int **cur_groups, **cur_groups2;

cur_groups=(int **)malloc((ntax-2)*sizeof(int *));
assert((int)cur_groups);
for (i=0;i<(ntax-2);i++) {
	cur_groups[i]=(int *)malloc(ntax*sizeof(int));/*only ntaxa because 0 not in any of the groups (always the outtaxon)*/
	assert((int) cur_groups[i]);
	for (j=0;j<ntax;j++) cur_groups[i][j]=0;
	}
get_cur_groups_tree(cur_groups,cur_nodes,ntax,ntax,ntax-2);
cur_groups2=(int **)malloc((ntax-2)*sizeof(int *));
assert((int)cur_groups2);
for (i=0;i<(ntax-2);i++) {
	cur_groups2[i]=(int *)malloc(ntax*sizeof(int));/*only ntaxa because 0 not in any of the groups (always the outtaxon)*/
	assert((int) cur_groups2[i]);
	for (j=0;j<ntax;j++) cur_groups2[i][j]=0;
	}
	

for (i=0;i<num_new;i++) {
	get_cur_groups_tree(cur_groups2,new_nodes[i],ntax,ntax,ntax-2);
	test=compare_groups_tree(cur_groups,cur_groups2,ntax-2,ntax-2,ntax-1,ntax-1);
	if (test) {
		for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
		free(cur_groups2);
		for (j=0;j<ntax-2;j++) free(cur_groups[j]);
		free(cur_groups);
		return 0;
		}
	}
	

/*if (new_length==old_length) for (i=0;i<num_old;i++) {
	get_cur_groups_tree(cur_groups2,old_nodes[i],ntax,ntax,ntax-2);
	test=compare_groups_tree(cur_groups,cur_groups2,ntax-2,ntax-2,ntax-1,ntax-1);
	if (test) {
		for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
		free(cur_groups2);
		for (j=0;j<ntax-2;j++) free(cur_groups[j]);
		free(cur_groups);
		return 0;
		}
	}*/

for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
free(cur_groups2);
for (i=0;i<ntax-2;i++) free(cur_groups[i]);
free(cur_groups);
return 1;
}
									
void get_cur_groups_tree(groups,node,ntaxa,nent,ngroups2)
int **node,**groups;
int ntaxa,nent,ngroups2;
{
int i,j,k;
int found_one,*done;

/*Remeber ntaxa is NOT!!! ntax+1 for rooting*/ 
/*first determine number of groups--one for each node*/
done=(int *)malloc(((2*ntaxa)-1)*sizeof(int));
assert((int) done);
/*for (i=1;i<nent-1;i++) recurse_groups(node,ntaxa,nent,i,i-1,groups);*/
for (i=1;i<ntaxa-1;i++) {
	found_one=1;
	for (j=0;j<(2*ntaxa)-1;j++) done[j]=0;
	done[i+ntaxa]=1;
	while (found_one) {
		found_one=0;
		for (j=1;j<ntaxa-1;j++) if ((done[j+ntaxa]) && (!done[node[j][0]]) && (!done[node[j][1]])) {
		found_one=1;
			done[node[j][0]]=done[node[j][1]]=1;
			if (node[j][0]<ntaxa) groups[i-1][node[j][0]]=1;
			if (node[j][1]<ntaxa) groups[i-1][node[j][1]]=1;
			}
		}
	}

free(done);
for (i=0;i<ngroups2;i++) for (j=0;j<ntaxa-1;j++) groups[i][ntaxa-1]+=groups[i][j];
}

int compare_groups_tree(g1,g2,ngroups1,ngroups2,ntax,nent)
int **g1,**g2;
int ngroups1,ngroups2,ntax,nent;
{
int i,j,k,holder,temp,disjunct,are_both_one, are_not_same;
int g1_tot,g2_tot;

for (k=0;k<ngroups1;k++) {
	for (i=0;i<ngroups2;i++) {
		/*check make sure each group is consisten with others*/
		/*are they completely disjunct*/
		are_both_one=are_not_same=g1_tot=g2_tot=0;
		for (j=0;j<ntax-1;j++) {
			if ((g1[k][j]+g2[i][j])==2) are_both_one+=1;
			else if ((g1[k][j]+g2[i][j])==1) are_not_same=1;
			g1_tot+=g1[k][j];
			g2_tot+=g2[i][j];
			}
		if (are_not_same && are_both_one) {
			if ((are_both_one != g1_tot) && (are_both_one != g2_tot)) return 0;
			}
		/*old way if (g1[k][ntax]==g2[i][ntax]) {
			holder=0;
			for (j=0;j<ngroups2;j++ ) {
				if (g1[k][j]!=g2[i][j]) {
					holder=1;
					break;
					}
				}
			if (!holder) temp=1;
			}*/
		}
	
	}
return 1;
}

void recurse_groups_tree(node,ntaxa,nent,i,k,groups)
int **node;
int ntaxa,nent,k,i;
int **groups;
{
int ii;
fprintf(stderr,"R ");
for (ii=0;ii<2;ii++){
	if (node[i][ii]<ntaxa) groups[k][node[i][ii]]=1;
	else if (node[i][ii]>=ntaxa) recurse_groups(node,ntaxa,nent,node[i][ii]-ntaxa,k,groups);
	}

}

int get_interval2(td1,td2,td3,values,nbases,n_gaps,value)/*add cut out for gap*/
int *td1,*td2,*td3;
int nbases,n_gaps,value;
parameters *values;
{
int m,ch,gap,ti,tv,ret_val,q;

ch=values->change_cost;
gap=values->gap_cost;
ti=values->transition;
tv=values->transversion;

ret_val=0;

/*Down pass if remove > can keep ties*/
/*fprintf(stderr,"{");*/
if (!values->ttr) {
	for (m=0;m<nbases;m++) {
		if (!(td1[m]&(td2[m]|td3[m]))) ret_val+=ch;
		if (ret_val>=value) return HUGE_COST;
		}
	}
else {
	for (m=0;m<nbases-values->n_transv;m++) {
		if (!(td1[m]&(td2[m]|td3[m]))) ret_val+=ti;
		if (ret_val>=value) return HUGE_COST;
		}
	for (m=nbases-values->n_transv;m<nbases;m++) {
		if (!(td1[m]&(td2[m]|td3[m]))) ret_val+=tv;
		if (ret_val>=value) return HUGE_COST;
		}            
	}
 for (m=nbases;m<(nbases+n_gaps);m++) {
	if (!(td1[m]&(td2[m]|td3[m]))) ret_val+=gap;
	if (ret_val>=value) return HUGE_COST;
	}
/*fprintf(stderr,"}");*/
return ret_val;
}


int optimize_non_blocked2(ntax,nodes,taxa,nbases,n_gaps,values,anc,desc,reopt)
int ntax, **nodes,**taxa,nbases,n_gaps,*anc,*desc,*reopt;
parameters *values;
{
int *td1,*td2,*td3,*td4,*td5;
int n,m,d1,d2;
int ***up_nodes,i_am,i;
int found_one;

up_nodes=(int ***)malloc(((2*ntax)-1)*sizeof(int**));
assert((int) up_nodes);
for (i=0;i<(2*ntax)-1;i++) {
	up_nodes[i]=(int **)malloc((3)*sizeof(int *));
	assert((int) up_nodes[i]);
	up_nodes[i][0]=(int *)malloc(((nbases+n_gaps))*sizeof(int)); /*holds down pass*/
	assert((int) up_nodes[i][0]);
	up_nodes[i][1]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*holds up pass to d1*/
	assert((int) up_nodes[i][1]);
	up_nodes[i][2]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*hold uppass to d2*/
	assert((int) up_nodes[i][2]);
	}

/*down pass*/
for (n=0;n<ntax;n++) {
	for (i=0;i<(nbases+n_gaps);i++) up_nodes[n][0][i]=taxa[n][i];
	values->tree_made[n]=1;
	}
/*add tree_made for done things*/
for (i=ntax;i<(2*ntax)-1;i++) {
	if (desc[i]) values->tree_made[i]=1;
	else values->tree_made[i]=0;
	}
found_one=1;
while (found_one) {
	found_one=0;
	for (n=(ntax-2);n>=0;--n){
		d1=nodes[n][0];
		d2=nodes[n][1];
		if ((values->tree_made[d1]) && (values->tree_made[d2]) && (!values->tree_made[ntax+n])) {
			found_one=1;
			td1=up_nodes[d1][0];
			td2=up_nodes[d2][0];
			td3=up_nodes[ntax+n][0];
			values->tree_made[ntax+n]=1;
			if (!values->ttr) {
				for (m=0;m<nbases+n_gaps;m++) /* if (reopt[m]) */{
					td3[m]=(td1[m]&td2[m]);
					if (!td3[m]) td3[m]=(td1[m]|td2[m]);
					}
				}
			}
		}
	}

/*up pass*/     
/*remember to remove extraneous up pass (ie when leads to a terminal*/
td1=up_nodes[ntax][0];
td2=up_nodes[ntax][1];
td3=up_nodes[ntax][2];
for (m=0;m<(nbases+n_gaps);m++) /* if (reopt[m]) */td1[m]=td2[m]=td3[m]=taxa[ntax][m]=taxa[0][m]; /*set basal node = to outgroup*/

for (i=ntax+1;i<(2*ntax)-1;i++) {
	if (desc[i]) values->tree_made[i]=1;
	else values->tree_made[i]=0;
	}
values->tree_made[ntax]=1;
found_one=1;
while (found_one) {
	found_one=0;
	for (n=1;n<ntax-1;n++) {
		if (!values->tree_made[ntax+n]) if (values->tree_made[anc[ntax+n]]) {
			found_one=1;
			values->tree_made[ntax+n]=1;
			if (nodes[anc[ntax+n]-ntax][0]==(n+ntax)) i_am=1;
			else i_am=2;
			td1=up_nodes[n+ntax][1];
			td2=up_nodes[n+ntax][2];
			td3=up_nodes[anc[ntax+n]][i_am];
			td4=up_nodes[nodes[n][1]][0];
			td5=up_nodes[nodes[n][0]][0];
			i=get_best_states_opt(taxa[ntax+n],td3,up_nodes[ntax+n][0],td5,td4,nbases+n_gaps,reopt); /*add reopt as an option after fix*/
			/*get up pass states for this node*/
			if (nodes[n][0]>(ntax-1)) {
				for (m=0;m<(nbases+n_gaps);m++)  /* if (reopt[m]) */{
					td1[m]=(td4[m]&td3[m]);
					if (!td1[m]) td1[m]=(td4[m]|td3[m]);
					}
				}
			if (nodes[n][1]>(ntax-1)) {
				for (m=0;m<(nbases+n_gaps);m++)  /* if (reopt[m]) */{
					td2[m]=(td5[m]&td3[m]);
					if (!td2[m]) td2[m]=(td5[m]|td3[m]);
					}
				}
			}
		}
	}
for (i=0;i<(2*ntax)-1;i++) {
	free(up_nodes[i][0]);
	free(up_nodes[i][1]);
	free(up_nodes[i][2]);
	free(up_nodes[i]);
	}
free(up_nodes);

return 1;
}

