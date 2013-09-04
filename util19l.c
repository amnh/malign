/* fprintfCopyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Goloboff*/

#include "align3.h"

alignment **do_swap_thang_align(a,ntax,nodes,l_in,best_aligns,values,a_up,anc,mult,score_holder,parallel_modify)
alignment **a, **best_aligns,***a_up;
int ntax, *l_in, **nodes,parallel_modify;
parameters *values;
int *anc,mult,score_holder;
{
int i,j,k,l,m,kk,num_clipped,check,tot_best,best,d1,d2;
int l_clip, l_rest,cur_best,found_better,node_removed;
int **cur_nodes, *blocked, *is_it_or_desc;
int **hold_nodes,**temp_nodes;
int clipped_base,***reroot_array,cur_num_trees;
int old_length, tree_counter;
int ***best_nodes_new,***best_nodes_old,num_found;
int new_start,all_others;
int where_it_was_anc,where_it_was_desc;
int **temp_nodes_h,whole_length;
int initial_length,d3,d4;
alignment *thing_to_add,*temp_align,*temp_align2;
int where_it_was_branch,paren_count;
int is_unique,sent,received;
int ii,jj,this_tree_counter,this_best;
int bytes,info,bufid,source,type;
int will_be_OK,**cur_groups,ngroups2;
int num_best, ***corres,all_rest;
alignment **best_aligns_holder,**a_down;
alignment **a_down_buf, **a_buf, ***a_up_buf;
int ***corres_buf,x,y,up_node,down_node;
alignment *up_temp, **check_buf;
int int_holder,in_trees,g_test;
char *collapse_name;

up_temp=NULL;
check_buf=NULL;
collapse_name=NULL;
values->support=NULL;
/*allocations*/
if (values->new_optimization) { 
	a_down=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) a_down);
	for (i=0;i<2*ntax-1;i++) a_down[i]=NULL;
	corres=(int ***)malloc((2*(ntax))*sizeof(int **));
	assert((int) corres);
	for (i=0;i<2*ntax;i++) corres[i]=NULL;
	for (i=0;i<ntax-1;i++) a_down[i]=make_align(a[i]);
	corres_buf=(int ***)malloc((2*(ntax))*sizeof(int **));
	assert((int) corres_buf);
	for (i=0;i<2*ntax;i++) corres_buf[i]=NULL;
	check_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) check_buf);
	for (i=0;i<2*ntax-1;i++) check_buf[i]=NULL;
	for (i=0;i<ntax-1;i++) check_buf[i]=make_align(a[i]);
	if (mult && (values->collapse != 2)) {
		collapse_name=(char *)malloc(MAX_SEQUENCE_SIZE*sizeof(char ));
		assert((int) collapse_name);
		values->support=(int **)malloc(values->keep_aligns*sizeof(int *));
		assert((int) values->support);
		for (i=0;i<values->keep_aligns;i++) {
			values->support[i]=(int *)malloc(ntax*sizeof(int));
			assert((int) values->support[i]);
			for (j=0;j<ntax;j++) values->support[0][j]=0;
			}
		}
	}
a_down_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_down_buf);
for (i=0;i<2*ntax-1;i++) a_down_buf[i]=NULL;
for (i=0;i<ntax-1;i++) a_down_buf[i]=make_align(a[i]);
if (values->tbr_align_shortcut) {
	a_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) a_buf);
	for (i=0;i<2*ntax-1;i++) a_buf[i]=NULL;
	for (i=0;i<ntax-1;i++) a_buf[i]=make_align(a[i]);
	a_up_buf=(alignment ***)malloc(((2*ntax)-1)*sizeof(alignment **));
	assert((int) a_up_buf);
	for (i=0;i<2*ntax-1;i++) {
		a_up_buf[i]=(alignment **)malloc(2*sizeof(alignment *));
		assert((int) a_up_buf[i]);
		a_up_buf[i][0]=NULL;
		a_up_buf[i][1]=NULL;
		}
	}
temp_nodes_h=(int **)malloc((ntax-1)*sizeof(int *));	
assert((int)temp_nodes_h);
for (j=0;j<ntax-1;j++) {
	temp_nodes_h[j]=(int *)malloc(2*sizeof(int));
	assert((int)temp_nodes_h);
	}
best_nodes_new=(int ***)malloc(values->keep_aligns*sizeof(int **));
assert((int) best_nodes_new);
best_nodes_old=(int ***)malloc(values->keep_aligns*sizeof(int **));
assert((int) best_nodes_old);
for (i=0;i<values->keep_aligns;i++) {
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
is_it_or_desc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) is_it_or_desc);
if ((values->atbr) || (values->arrt)) {
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
will_be_OK=1;
if (!mult) if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Initial build yielded an alignment at length %d\n",best_aligns[0]->score);
paren_count=0;
thing_to_add=NULL;
temp_align=NULL;
temp_align2=NULL;
for (k=0;k<ntax-1;k++) {
	best_nodes_new[0][k][0]=best_nodes_old[0][k][0]=nodes[k][0];
	best_nodes_new[0][k][1]=best_nodes_old[0][k][1]=nodes[k][1];
	}
new_start=num_found=0;
tot_best=old_length=initial_length=best_aligns[0]->score;
cur_num_trees=1;
found_better=1;
tree_counter=(*l_in);/*shoulf always be 1*/
while (found_better) {
	if (!(PARALLEL && values->jackboot)) fprintf(stderr,"(swapping");
	found_better=0;
	old_length=tot_best;
	for (m=new_start;m<cur_num_trees;m++) {
		/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"S on %d (%d->%d) ",m,new_start,cur_num_trees);*/
		/*copy nodes*/
		for (k=0;k<ntax-1;k++) {
			hold_nodes[k][0]=best_nodes_old[m][k][0];
			hold_nodes[k][1]=best_nodes_old[m][k][1];
			anc[hold_nodes[k][0]]=anc[hold_nodes[k][1]]=k+ntax;     
			}
		anc[ntax]=ntax;
		/*initialize buffer with down passes should be parallel if can*/
		if (values->new_optimization) {
			if (!(PARALLEL*parallel_modify)) {
				check=all_diagnose_tree_here(a_down_buf,hold_nodes,ntax,ntax,values,HUGE_COST);
			 	if (values->tbr_align_shortcut) get_up_pass_new_opt(a_buf,a_down_buf,a_up_buf,hold_nodes,ntax,values,ntax,anc);
				}
			else {
				check=all_diagnose_tree_here_parallel(a_down_buf,hold_nodes,ntax,ntax,values,0,corres_buf);
				if (values->tbr_align_shortcut) get_up_pass_new_opt_parallel(a_buf,a_down_buf,a_up_buf,hold_nodes,ntax,values,ntax,anc);
				}
			if (mult && (values->collapse!=2)) {
				g_test=0;
				for (k=0;k<ntax;k++) g_test+=values->support[0][k];
				/*collapse and get groups*/
				if (!g_test) check=get_collapsed_thang_with_groups_old(best_aligns[0],hold_nodes,ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
				}
			}
		else {
			if (!(PARALLEL*parallel_modify)) check=all_diagnose_tree_local2(a_down_buf,hold_nodes,ntax,ntax,values);/*get down passes*/
			else  check=all_diagnose_tree_local_parallel2(a_down_buf,hold_nodes,ntax,ntax,values);
			}
		sent=received=0;
		for (i=2;i<2*(ntax-1);i++) {
			/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"i%d ",i);*/
			/*copy nodes*/
			for (k=0;k<ntax-1;k++) {
				cur_nodes[k][0]=hold_nodes[k][0];
				cur_nodes[k][1]=hold_nodes[k][1];
				anc[hold_nodes[k][0]]=anc[hold_nodes[k][1]]=k+ntax;     
				}
			anc[ntax]=ntax;
			/*choose*/
			d1=i/2;
			d2=i%2;
			/*get size*/
			num_clipped=get_size(ntax,cur_nodes,d1,d2,blocked,i,is_it_or_desc,anc);
			if (!values->atbr) values->number+=((2*(ntax-num_clipped))-3);
			else {
				if (num_clipped>1) values->number+=(((2*num_clipped)-3)*((2*(ntax-num_clipped))-3));
				else values->number+=((2*(ntax-num_clipped))-3);
				}
			if (num_clipped<=(ntax-3)) {
				/*redo nodes*/
				node_removed=anc[cur_nodes[d1][d2]];
				where_it_was_anc=anc[node_removed];
				if (cur_nodes[anc[node_removed]-ntax][0]==node_removed) {
					cur_nodes[anc[node_removed]-ntax][0]=cur_nodes[d1][!d2];
					blocked[2*(anc[node_removed]-ntax)]=1;
					where_it_was_desc=cur_nodes[anc[node_removed]-ntax][0];
					where_it_was_branch=0;
					}
				else {
					cur_nodes[anc[node_removed]-ntax][1]=cur_nodes[d1][!d2];
					blocked[2*(anc[node_removed]-ntax)+1]=1;
					where_it_was_desc=cur_nodes[anc[node_removed]-ntax][1];
					where_it_was_branch=1;
					}
				anc[cur_nodes[d1][!d2]]=anc[node_removed];
				cur_nodes[d1][!d2]=cur_nodes[d1][d2];
				anc[node_removed]=node_removed;
				if (PARALLEL*parallel_modify) {
					if (sent < (values->num_hosts-1)) {
						pvm_initsend( PvmDataDefault );
			          		/*pack stuff*/
						pvm_pkint(&num_clipped,1,1);
						for (ii=0;ii<ntax-1;ii++) pvm_pkint(cur_nodes[ii],2,1);
						pvm_pkint(anc,((2*ntax)-1),1);
						pvm_pkint(&d1,1,1);
						pvm_pkint(&d2,1,1);
						pvm_pkint(&where_it_was_desc,1,1);
						pvm_pkint(&where_it_was_anc,1,1);
						pvm_pkint(&where_it_was_branch,1,1);
						pvm_pkint(&score_holder,1,1);
						pvm_pkint(&old_length,1,1);
						pvm_pkint(blocked,(2*(ntax-1)),1);
						pvm_pkint(&node_removed,1,1);
						pvm_pkint(&mult,1,1);
						pvm_pkint(is_it_or_desc,((2*ntax)-1),1);
						for (ii=ntax+1;ii<2*ntax-1;ii++) pack_align(a_down_buf[ii]);
						if (values->tbr_align_shortcut) for (ii=ntax+1;ii<2*ntax-1;ii++) pack_align(a_buf[ii]);
						pvm_send(values->tids[++sent],PVM_NEW_SWAP);
						/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"S->%d ",sent);*/
			   		     	}
			   		  else {
			   			bufid=pvm_recv(-1,PVM_NEW_SWAP_DONE);
		 	         		if (bufid) {
		        	       			info=pvm_bufinfo(bufid, &bytes, &type, &source);
		        	       			/*unpack stuff*/
		        	       			pvm_upkint(&this_best,1,1); /*pack cost*/
		        	       			if ((this_best<tot_best) || ((this_best==tot_best) && mult && (tree_counter< values->keep_aligns))) {
								pvm_upkint(&this_tree_counter,1,1);/*pack #solns*/
								if (!mult) {
									if (this_best<tot_best) this_tree_counter=1;
									else this_tree_counter=0;
									}
								if (this_best<tot_best) {
									if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",this_best);
									found_better=1;
									tot_best=this_best;
									/*unpack the first and then treat the rest like additional solutions*/
									for (ii=0;ii<tree_counter;ii++) if (best_aligns[ii]) best_aligns[ii]=dump_align(best_aligns[ii]);
									for (jj=0;jj<ntax-1;jj++) pvm_upkint(best_nodes_new[0][jj],2,1);
									best_aligns[0]=unpack_align_and_score(values);
										if (values->new_optimization && (values->collapse!=2)){
											if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values);
											else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values,values->support[0]);
											best_aligns[0]->score=tot_best;
											}
									tree_counter=1;
									/*rest 1->this_tree_counter*/
									for (ii=1;ii<this_tree_counter;ii++) { /*make sure OK when no solutions found*/
										for (jj=0;jj<ntax-1;jj++) pvm_upkint(temp_nodes_h[jj],2,1);
										if (temp_align) temp_align=dump_align(temp_align);
										temp_align=unpack_align_and_score(values);
										is_unique=1;
										if (tree_counter < values->keep_aligns) {
											if (!values->new_optimization) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
											else {
												if (values->collapse==2) is_unique=is_unique_tree(temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
												else is_unique=is_unique_tree_align(best_aligns,temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
												}
											if (is_unique) {
												if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
												found_better=1;
												best_aligns[tree_counter]=make_align(temp_align);
												for (k=0;k<ntax-1;k++) {
													best_nodes_new[tree_counter][k][0]=temp_nodes_h[k][0];
													best_nodes_new[tree_counter][k][1]=temp_nodes_h[k][1];
													}
												if (values->new_optimization && (values->collapse!=2)) {
													free(best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
													assert((int) best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
													}
												++tree_counter;
												}
											}
										}									
									}
								else {
									for (ii=0;ii<this_tree_counter;ii++) { /*make sure OK when no solutions found*/
										for (jj=0;jj<ntax-1;jj++) pvm_upkint(temp_nodes_h[jj],2,1);
										if (temp_align) temp_align=dump_align(temp_align);
										temp_align=unpack_align_and_score(values);
										is_unique=1;
										if (tree_counter < values->keep_aligns) {
											if (!values->new_optimization) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
											else {
												if (values->collapse==2) is_unique=is_unique_tree(temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
												else is_unique=is_unique_tree_align(best_aligns,temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
												}
											if (is_unique) {
												if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
												found_better=1;
												best_aligns[tree_counter]=make_align(temp_align);
												for (k=0;k<ntax-1;k++) {
													best_nodes_new[tree_counter][k][0]=temp_nodes_h[k][0];
													best_nodes_new[tree_counter][k][1]=temp_nodes_h[k][1];
													}
												if (values->new_optimization && (values->collapse!=2)) {
													free(best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
													assert((int) best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
													}
												++tree_counter;
												}
											}
										}
									}/*else equal*/
								}
							++received;
		        	       			}
		        	       		pvm_initsend( PvmDataDefault );
			          		/*pack stuff*/
			   			pvm_pkint(&num_clipped,1,1);
						for (ii=0;ii<ntax-1;ii++) pvm_pkint(cur_nodes[ii],2,1);
						pvm_pkint(anc,((2*ntax)-1),1);
						pvm_pkint(&d1,1,1);
						pvm_pkint(&d2,1,1);
						pvm_pkint(&where_it_was_desc,1,1);
						pvm_pkint(&where_it_was_anc,1,1);
						pvm_pkint(&where_it_was_branch,1,1);
						pvm_pkint(&score_holder,1,1);
						pvm_pkint(&old_length,1,1);
						pvm_pkint(blocked,(2*(ntax-1)),1);
						pvm_pkint(&node_removed,1,1);
						pvm_pkint(&mult,1,1);
						pvm_pkint(is_it_or_desc,((2*ntax)-1),1);
						for (ii=ntax+1;ii<2*ntax-1;ii++) pack_align(a_down_buf[ii]);
						if (values->tbr_align_shortcut) for (ii=ntax+1;ii<2*ntax-1;ii++) pack_align(a_buf[ii]);
						pvm_send(source,PVM_NEW_SWAP);
			   			++sent;
			   			} 
					}
				else {
				/*when one or two taxa*/
				if ((num_clipped < 3) || (!values->atbr)) { /*if no alternative roots or sbr*/
					if (!values->new_optimization) {
						check=all_diagnose_tree_local2(a,cur_nodes,ntax,ntax,values);
						/*check=all_diagnose_tree_local3(a,cur_nodes,ntax,ntax,values,is_it_or_desc);
						for (k=ntax+1;k<2*ntax-1;k++) if (is_it_or_desc[k]) {
							if (a[k]) a[k]=dump_align(a[k]);
							a[k]=make_align(a_down_buf[k]);
							}*/
						get_up_pass2(a,cur_nodes,ntax,values,a_up,ntax,anc);
						}
					else {
						check=all_diagnose_tree_local_new_opt2(a_down,cur_nodes,ntax,ntax,values,node_removed,corres,is_it_or_desc);
						/*get up_pass alignments*/
						get_up_pass_new_opt2(a,a_down,a_up,cur_nodes,ntax,values,ntax,anc,corres,node_removed,is_it_or_desc);
						}
					/* get "best" set sequences[cur_nodes[d1][d2]] as XUY */
					d3=cur_nodes[d1][d2];
					if (d3<ntax) d3-=1;
					d4=where_it_was_desc;
					if (d4<ntax) d4-=1;
					if (thing_to_add) thing_to_add=dump_align(thing_to_add);
					thing_to_add=make_align(a_down_buf[d3]);
					if (!values->new_optimization) {
						if (temp_align) temp_align=dump_align(temp_align);
						if (where_it_was_anc>ntax) temp_align=nw(a_up[where_it_was_anc][where_it_was_branch],a[d4],values);
						else temp_align=make_align(a[d4]);
						if (score_holder>0) {
							values->phylo_score=score_holder;
							temp_align->score=cladogram(temp_align,values);
							values->phylo_score=0;
							}
						/*redo score with phyloshit*/
						best=old_length-temp_align->score;
						}
					else {
						all_rest=thing_to_add->score + a_down[cur_nodes[0][1]]->score;
						best=old_length - all_rest;
						/*values->in_optimize=1;
						temp_align2=nw(a[where_it_was_anc],a[d4],values);
						values->in_optimize=0;
						temp_align2=newer_make_ambig(temp_align2,values);
						temp_align=nw(temp_align2,thing_to_add,values);
						if (!(PARALLEL*parallel_modify)) fprintf(stderr,"{OB %d NB ",best);
						best=temp_align2->score;
						if (!(PARALLEL*parallel_modify)) fprintf(stderr,"%d}",best);
						all_rest=old_length-best;
						if (temp_align) temp_align=dump_align(temp_align);
						if (temp_align2) temp_align2=dump_align(temp_align2);*/
						/*values->in_optimize=1;
						temp_align=nw(a[where_it_was_anc],a[cur_nodes[where_it_was_anc - ntax][where_it_was_branch]],values);
						values->in_optimize=0;
						temp_align=newer_make_ambig(temp_align,values);
						temp_align2=nw(thing_to_add,temp_align,values);
						if (!(PARALLEL*parallel_modify)) fprintf(stderr,"OB %d NB ",best);
						best=temp_align2->score;
						if (!(PARALLEL*parallel_modify)) fprintf(stderr,"%d ",best);
						all_rest=old_length-best;
						if (!(PARALLEL*parallel_modify)) fprintf(stderr,"{%d %d %d %d}",best,thing_to_add->score,a_down[cur_nodes[0][1]]->score,old_length);
						if (temp_align) temp_align=dump_align(temp_align);
						if (temp_align2) temp_align2=dump_align(temp_align2);*/
						}
					if (((best==0) && (mult)) || (best>0)) {/*to make sure that that sequence matters ie adds length*/
						/*loop through places to put back*/
						blocked[1]=0;
						for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
							if (values->groups) {
								for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=cur_nodes[ii][jj];
								check=temp_nodes[j/2][j%2];
								temp_nodes[j/2][j%2]=node_removed;
								temp_nodes[d1][!d2]=check;
								cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
								will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
								for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
								free(cur_groups);
								}
							if (will_be_OK) {
								/* get_cost*/
								if (temp_align) temp_align=dump_align(temp_align);
								if (temp_align2) temp_align2=dump_align(temp_align2);
								d3=cur_nodes[j/2][j%2];
								if (d3<ntax) d3-=1;
								if (!values->new_optimization) {
									temp_align=nw(thing_to_add,a[d3],values);
									if (j>1) {
										/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"{");*/
										d4=(j/2)+ntax;
										temp_align2=nw(temp_align,a_up[d4][j%2],values);
										if (score_holder>0) {
											if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
											else {
												values->phylo_score=score_holder;
												temp_align2->score=cladogram(temp_align2,values);
												values->phylo_score=0;
												if (values->previous) values->previous=dump_align(values->previous);
												values->previous=make_align(temp_align2);
												}
											}
										cur_best=temp_align2->score;
										/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"}");*/
										}
									else {
										/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"(");*/
										if (score_holder>0) {
											if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
											else {
												values->phylo_score=score_holder;
												temp_align->score=cladogram(temp_align,values);
												values->phylo_score=0;
												if (values->previous) values->previous=dump_align(values->previous);
												values->previous=make_align(temp_align);
												}
											}
										cur_best=temp_align->score;
										/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,")");*/
										}
									}
								else {
									d4=(j/2)+ntax;
									if (j>1) {
										values->in_optimize=1;
										temp_align2=nw(a[d3],a[d4],values);
										values->in_optimize=0;
										temp_align2=newer_make_ambig(temp_align2,values);
										temp_align=nw(temp_align2,thing_to_add,values);
										}
									else temp_align=nw(thing_to_add,a_down[cur_nodes[0][1]],values);
									cur_best=temp_align->score+all_rest;
									}
								/*filter for type wieight problem for now*/								
								if (values->new_optimization) {
									/*if (1) {*/
									if ((cur_best < tot_best) || ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns))) {
										for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=cur_nodes[ii][jj];
										check=temp_nodes[j/2][j%2];
										temp_nodes[j/2][j%2]=node_removed;
										temp_nodes[d1][!d2]=check;
										check=cur_best;
										cur_best=all_diagnose_tree_local_new_opt(check_buf,temp_nodes,ntax,ntax,values,0,corres);
										/*if (cur_best != check) if (!(PARALLEL*parallel_modify)) fprintf(stderr,".(%d %d)",check,cur_best);
										else if (!(PARALLEL*parallel_modify)) fprintf(stderr,"+");*/
										}
									}
								if (cur_best < tot_best) {
									if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",cur_best);
									tot_best=cur_best;
									found_better=1;
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[0][k][0]=cur_nodes[k][0];
										best_nodes_new[0][k][1]=cur_nodes[k][1];
										}
									check=best_nodes_new[0][j/2][j%2];
									best_nodes_new[0][j/2][j%2]=node_removed;
									best_nodes_new[0][d1][!d2]=check;
									/*add best_aligns stuf free first*/
									for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
									if (temp_align2) best_aligns[0]=make_align(temp_align2);
									else best_aligns[0]=make_align(temp_align);
									if (values->new_optimization && (values->collapse != 2)) {
										/*redo name for collapse*/
										if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values);
										else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
										best_aligns[0]->score=tot_best;
										}
									tree_counter=1;
									} 
								else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
									/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"F");*/
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[tree_counter][k][0]=cur_nodes[k][0];
										best_nodes_new[tree_counter][k][1]=cur_nodes[k][1];
										}
									check=best_nodes_new[tree_counter][j/2][j%2];
									best_nodes_new[tree_counter][j/2][j%2]=node_removed;
									best_nodes_new[tree_counter][d1][!d2]=check;
									/*chec to see if novel*/
									is_unique=1;
									if (!values->new_optimization) {
										if (temp_align2) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align2,best_aligns[k]);
										else for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
										}
									else {
										if (values->collapse==2) is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax);
										else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*1,values,collapse_name);
										}
									if (is_unique) { 
										if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
										found_better=1;
										if (temp_align2) best_aligns[tree_counter]=make_align(temp_align2);
										else best_aligns[tree_counter]=make_align(temp_align);
										if (values->new_optimization && (values->collapse!=2)) {
											free(best_aligns[tree_counter]->name);
											best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
											assert((int) best_aligns[tree_counter]->name);
											best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
											}
										++tree_counter;
										}
									/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"O");*/
									}
								} /*SBR loop--regrafted*/
							}/*will it be OK*/
						} /*BEST split size ok*/
					}/*too few clipped or SBR*/
				else { /*TBR with at least three taxa*/
					check=get_root_array(reroot_array,num_clipped,node_removed,cur_nodes,d1,d2,ntax,is_it_or_desc,anc);
					if (!values->new_optimization) {
						/*fix so no do if in group*/
						check=all_diagnose_tree_local2(a,cur_nodes,ntax,ntax,values);
						/*check=all_diagnose_tree_local3(a,cur_nodes,ntax,ntax,values,is_it_or_desc);
						for (k=ntax+1;k<2*ntax-1;k++) if (is_it_or_desc[k]) {
							if (a[k]) a[k]=dump_align(a[k]);
							a[k]=make_align(a_down_buf[k]);
							}*/
						get_up_pass2(a,cur_nodes,ntax,values,a_up,ntax,anc);
						}
					else {
						/*only down pass if ! in picked*/
						check=all_diagnose_tree_local_new_opt2(a_down,cur_nodes,ntax,ntax,values,node_removed,corres,is_it_or_desc);
						/*get up_pass alignments*/
						get_up_pass_new_opt2(a,a_down,a_up,cur_nodes,ntax,values,ntax,anc,corres,node_removed,is_it_or_desc);
						}
					for (l=0;l<2*num_clipped-3;l++) {
						/*only down pass if in picked*/
						d3=reroot_array[l][d1][d2];
						if (d3<ntax) d3-=1;
						d4=where_it_was_desc;
						if (d4<ntax) d4-=1;
						if (thing_to_add) thing_to_add=dump_align(thing_to_add);
						/*if same rooting as original use the buffer*/
						if ((cur_nodes[cur_nodes[d1][d2]-ntax][0]==reroot_array[l][reroot_array[l][d1][d2]-ntax][0]) &&(cur_nodes[cur_nodes[d1][d2]-ntax][1]==reroot_array[l][reroot_array[l][d1][d2]-ntax][1])) thing_to_add=make_align(a_down_buf[d3]);
						else if (values->tbr_align_shortcut) { /*X U Y*/
							if (values->new_optimization) {
								x=reroot_array[l][reroot_array[l][d1][d2]-ntax][0];
								y=reroot_array[l][reroot_array[l][d1][d2]-ntax][1];
								if (x<ntax) x-=1;
								if (y<ntax) y-=1;
								thing_to_add=nw(a_buf[x],a_buf[y],values);
								}
							else {/*currently disabled*/
								/*get X and Y*/
								x=reroot_array[l][reroot_array[l][d1][d2]-ntax][0];
								y=reroot_array[l][reroot_array[l][d1][d2]-ntax][1];
								if (x<ntax) {
									down_node=x;
									up_node=y;
																			}
								else if (y<ntax) {
									down_node=y;
									up_node=x;
									}
								else {
									if (cur_nodes[y-ntax][0]==x) { down_node=x; up_node=y;}
									else if (cur_nodes[y-ntax][1]==x) { down_node=x; up_node=y; }
									else if (cur_nodes[x-ntax][0]==y) { down_node=y; up_node=x; }
									else { down_node=y; up_node=x; }

									}
								/*get the up*/
								if (up_temp) up_temp=dump_align(up_temp);
								up_temp=get_single_up(a_down_buf[down_node],a_down_buf[cur_nodes[d1][d2]],values);
								/*make it*/
								if (down_node<ntax) x-=1;
								thing_to_add=nw(a_down_buf[down_node],up_temp,values);
								if (up_temp) up_temp=dump_align(up_temp);

								}

							
							}
						else {
							if (!values->new_optimization) check=all_diagnose_tree_local3(a,reroot_array[l],ntax,ntax,values,is_it_or_desc); /*,is_it_or_desc); get down passes*/
							else check=all_diagnose_tree_local_new_opt3(a_down,reroot_array[l],ntax,ntax,values,node_removed,corres,is_it_or_desc);
							if (!values->new_optimization) thing_to_add=make_align(a[d3]);
							else thing_to_add=make_align(a_down[d3]);
							}
						if (!values->new_optimization) {
							if (temp_align) temp_align=dump_align(temp_align);
							if (where_it_was_anc>ntax) temp_align=nw(a_up[where_it_was_anc][where_it_was_branch],a[d4],values);
							else temp_align=make_align(a[d4]);
							if (score_holder>0) {
								values->phylo_score=score_holder;
								temp_align->score=cladogram(temp_align,values);
								values->phylo_score=0;
								}
							/*redo score with phyloshit*/
							best=old_length-temp_align->score;
							}
						else {
							all_rest=thing_to_add->score + a_down[reroot_array[l][0][1]]->score;
							best=old_length - all_rest;
							}
						if ((best>0) || ((best==0) && (mult))) {/*to make sure that that sequence matters ie adds length*/
							/*loop through places to put back*/
							blocked[1]=0;
							for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
								if (values->groups) {
									for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=reroot_array[l][ii][jj];
									check=temp_nodes[j/2][j%2];
									temp_nodes[j/2][j%2]=node_removed;
									temp_nodes[d1][!d2]=check;
									cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
									will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
									for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
									free(cur_groups);
									}
								if (will_be_OK) {
									/* get_cost*/
									if (temp_align) temp_align=dump_align(temp_align);
									if (temp_align2) temp_align2=dump_align(temp_align2);
									d3=reroot_array[l][j/2][j%2];
									if (d3<ntax) d3-=1;
									if (!values->new_optimization) {
										temp_align=nw(thing_to_add,a[d3],values);
										if (j>1) {
											/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"{");*/
											d4=(j/2)+ntax;
											temp_align2=nw(temp_align,a_up[d4][j%2],values);
											if (score_holder>0) {
												if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
												else {
													values->phylo_score=score_holder;
													temp_align2->score=cladogram(temp_align2,values);
													values->phylo_score=0;
													if (values->previous) values->previous=dump_align(values->previous);
													values->previous=make_align(temp_align2);
													}
												}
											cur_best=temp_align2->score;
											/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"}");*/
											}
										else {
											/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"(");*/
											if (score_holder>0) {
												if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
												else {
													values->phylo_score=score_holder;
													temp_align->score=cladogram(temp_align,values);
													values->phylo_score=0;
													if (values->previous) values->previous=dump_align(values->previous);
													values->previous=make_align(temp_align);
													}
												}
											cur_best=temp_align->score;
											/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,")");*/
											}
										}
									else {
										d4=(j/2)+ntax;
										if (j>1) {
											values->in_optimize=1;
											temp_align2=nw(a[d3],a[d4],values);
											values->in_optimize=0;
											temp_align2=newer_make_ambig(temp_align2,values);
											temp_align=nw(temp_align2,thing_to_add,values);
											}
										else temp_align=nw(thing_to_add,a_down[reroot_array[l][0][1]],values);
										cur_best=temp_align->score+all_rest;
										}
									if (values->new_optimization) {
										if ((cur_best < tot_best) || ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns))) {
											for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=reroot_array[l][ii][jj];
											check=temp_nodes[j/2][j%2];
											temp_nodes[j/2][j%2]=node_removed;
											temp_nodes[d1][!d2]=check;
											check=cur_best;
											cur_best=all_diagnose_tree_local_new_opt(check_buf,temp_nodes,ntax,ntax,values,0,corres);
											/*if (cur_best != check) if (!(PARALLEL*parallel_modify)) fprintf(stderr,".");
											else if (!(PARALLEL*parallel_modify)) fprintf(stderr,"+");*/
											}
										}
									if (cur_best < tot_best) {
										if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",cur_best);
										tot_best=cur_best;
										found_better=1;
										for (k=0;k<ntax-1;k++) {
											best_nodes_new[0][k][0]=reroot_array[l][k][0];
											best_nodes_new[0][k][1]=reroot_array[l][k][1];
											}
										check=best_nodes_new[0][j/2][j%2];
										best_nodes_new[0][j/2][j%2]=node_removed;
										best_nodes_new[0][d1][!d2]=check;
										/*add best_aligns stuf free first*/
										for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
										if (temp_align2) best_aligns[0]=make_align(temp_align2);
										else best_aligns[0]=make_align(temp_align);
										if (values->new_optimization && (values->collapse != 2)) {
											/*redo name for collapse*/
											if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values);
											else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
											best_aligns[0]->score=tot_best;											
											}
										tree_counter=1;
										} 
									else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
										for (k=0;k<ntax-1;k++) {
											best_nodes_new[tree_counter][k][0]=reroot_array[l][k][0];
											best_nodes_new[tree_counter][k][1]=reroot_array[l][k][1];
											}
										check=best_nodes_new[tree_counter][j/2][j%2];
										best_nodes_new[tree_counter][j/2][j%2]=node_removed;
										best_nodes_new[tree_counter][d1][!d2]=check;
										/*chec to see if novel*/
										is_unique=1;
										if (!values->new_optimization) {
											if (temp_align2) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align2,best_aligns[k]);
											else for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
											}
										else {
											if (values->collapse==2) is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax);
											else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*1,values,collapse_name);
											}
										if (is_unique) {
											if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
											found_better=1;
											if (temp_align2) best_aligns[tree_counter]=make_align(temp_align2);
											else best_aligns[tree_counter]=make_align(temp_align);
											if (values->new_optimization && (values->collapse!=2)) {
												free(best_aligns[tree_counter]->name);
												best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
												assert((int) best_aligns[tree_counter]->name);
												best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
												}
											++tree_counter;
											}
										}
									} /*SBR loop--regrafted*/
								}
							} /*BEST split size ok*/	
						} 
					} /*reroot TBR/SBR choice*/
					}/*PARALLEL*parallel_modify*/
				} /*size of plucked OK*/
			} /*i clades to pick*/
		if (PARALLEL*parallel_modify) { /*receive the rest after all I's are sent out and some received */
			/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"In cleanup ");*/
			while (received<sent) {
				bufid=pvm_recv(-1,PVM_NEW_SWAP_DONE);
			   		if (bufid) {
			    			info=pvm_bufinfo(bufid, &bytes, &type, &source);
			    			/*unpack stuff*/
			    			pvm_upkint(&this_best,1,1); /*pack cost*/
			        		if ((this_best<tot_best) || ((this_best==tot_best) && mult && (tree_counter< values->keep_aligns))) {
						pvm_upkint(&this_tree_counter,1,1);/*pack #solns*/
						if (!mult) {
							if (this_best<tot_best) this_tree_counter=1;
							else this_tree_counter=0;
							}
						if (this_best<tot_best) {
									if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",this_best);
									found_better=1;
									tot_best=this_best;
									/*unpack the first and then treat the rest like additional solutions*/
									for (ii=0;ii<tree_counter;ii++) if (best_aligns[ii]) best_aligns[ii]=dump_align(best_aligns[ii]);
									for (jj=0;jj<ntax-1;jj++) pvm_upkint(best_nodes_new[0][jj],2,1);
										best_aligns[0]=unpack_align_and_score(values);
										if (values->new_optimization && (values->collapse!=2)){
											if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values);
											else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values,values->support[0]);
											best_aligns[0]->score=tot_best;
											}
									tree_counter=1;
									/*rest 1->this_tree_counter*/
									for (ii=1;ii<this_tree_counter;ii++) { /*make sure OK when no solutions found*/
										for (jj=0;jj<ntax-1;jj++) pvm_upkint(temp_nodes_h[jj],2,1);
										if (temp_align) temp_align=dump_align(temp_align);
										temp_align=unpack_align_and_score(values);
										is_unique=1;
										if (tree_counter < values->keep_aligns) {
											if (!values->new_optimization) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
											else {
												if (values->collapse==2) is_unique=is_unique_tree(temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
												else is_unique=is_unique_tree_align(best_aligns,temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
												}
											if (is_unique) {
												if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
												found_better=1;
												best_aligns[tree_counter]=make_align(temp_align);
												for (k=0;k<ntax-1;k++) {
													best_nodes_new[tree_counter][k][0]=temp_nodes_h[k][0];
													best_nodes_new[tree_counter][k][1]=temp_nodes_h[k][1];
													}
												if (values->new_optimization && (values->collapse!=2)) {
													free(best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
													assert((int) best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
													}
												++tree_counter;
												}
											}
										}									
									}
					/*old way of gathering better solutions
					if (this_best<tot_best) {
							if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",this_best);
							found_better=1;
							for (ii=0;ii<tree_counter;ii++) if (best_aligns[ii]) best_aligns[ii]=dump_align(best_aligns[ii]);
							tot_best=this_best;
							tree_counter=this_tree_counter;
							for (ii=0;ii<tree_counter;ii++){
								for (jj=0;jj<ntax-1;jj++) pvm_upkint(best_nodes_new[ii][jj],2,1);
								best_aligns[ii]=unpack_align_and_score(values);
								if (values->new_optimization && (values->collapse!=2)) {
									if (!mult) check=get_collapsed_thang(best_aligns[ii],best_nodes_new[ii],ntax,check_buf,PARALLEL*parallel_modify*1,values);
									else check=get_collapsed_thang_with_groups_old(best_aligns[ii],best_nodes_new[ii],ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
									best_aligns[ii]->score=tot_best;
									}
								}
							}*/
						else {
							for (ii=0;ii<this_tree_counter;ii++) { /*make sure OK when no solutions found*/
								for (jj=0;jj<ntax-1;jj++) pvm_upkint(temp_nodes_h[jj],2,1);
								if (temp_align) temp_align=dump_align(temp_align);
								temp_align=unpack_align_and_score(values);
								is_unique=1;
								if (!values->new_optimization) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
								else {
									if (values->collapse==2) is_unique=is_unique_tree(temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
									else is_unique=is_unique_tree_align(best_aligns,temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*1,values,collapse_name);
									}
								if (is_unique) {
									if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
									found_better=1;
									best_aligns[tree_counter]=make_align(temp_align);
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[tree_counter][k][0]=temp_nodes_h[k][0];
										best_nodes_new[tree_counter][k][1]=temp_nodes_h[k][1];
										}
									if (values->new_optimization && (values->collapse!=2)) {
											free(best_aligns[tree_counter]->name);
											best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
											assert((int) best_aligns[tree_counter]->name);
											best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
											}
									++tree_counter;
									}/*unique*/
								}/*ii*/
							}
						}/*if worth dealing*/
					}/*bufid*/
				++received;
			    }/*while*/
			} /*PARALLEL*parallel_modify*/
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
	/*reset num_trees*/
	cur_num_trees=tree_counter;
	/*reroot swapping*/
	if (values->arrt) {
		if ((found_better) || (!mult && (new_start==0))) { /*found new or original*/
			in_trees=tree_counter;
			for (kk=new_start;kk<in_trees;kk++) {
				/*Get rerootings for each NEW or original solutions*/
				for (l=1;l<2*ntax-1;l++) is_it_or_desc[l]=1;
				is_it_or_desc[ntax]=0;
				for (x=0;x<ntax-1;x++) {
					anc[best_nodes_old[kk][x][0]]=x+ntax;
					anc[best_nodes_old[kk][x][1]]=x+ntax;
					}
				check=get_root_array(reroot_array,ntax-1,best_nodes_old[kk][0][1],best_nodes_old[kk],0,1,ntax,is_it_or_desc,anc); /*check these for basal reroot*/
				/*get complete down pass length*/
				values->number+=((2*(ntax-1))-3);
				for (l=0;l<(2*(ntax-1))-3;l++) {
					if (!values->new_optimization) {
						if (!(PARALLEL*parallel_modify)) cur_best=all_diagnose_tree_local2(a_down_buf,reroot_array[l],ntax,ntax,values);
						else cur_best=all_diagnose_tree_local_parallel2(a_down_buf,reroot_array[l],ntax,ntax,values);
						}
					else {
						if (!(PARALLEL*parallel_modify)) cur_best=all_diagnose_tree_local_new_opt(a_down_buf,reroot_array[l],ntax,ntax,values,0,corres);
						else cur_best=all_diagnose_tree_local_parallel_new_opt(a_down_buf,reroot_array[l],ntax,ntax,values,0,corres);
						}
					if (!values->new_optimization && (score_holder>0)) {
						values->phylo_score=score_holder;
						cur_best=a_down_buf[reroot_array[l][0][1]]->score=cladogram(a_down_buf[reroot_array[l][0][1]],values);
						values->phylo_score=0;
						}
					/*for (x=0;x<ntax-1;x++) if (!(PARALLEL*parallel_modify)) fprintf(stderr,"%d %d,",reroot_array[l][x][0],reroot_array[l][x][1]);
					if (!(PARALLEL*parallel_modify)) fprintf(stderr,"%d\n",cur_best);*/
					/*compare etc*/
					if (cur_best < tot_best) {
						if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d (R)\n",cur_best);
						tot_best=cur_best;
						found_better=1;
						for (k=0;k<ntax-1;k++) {
							best_nodes_new[0][k][0]=reroot_array[l][k][0];
							best_nodes_new[0][k][1]=reroot_array[l][k][1];
							}
						/*add best_aligns stuf free first*/
						for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
						best_aligns[0]=make_align(a_down_buf[reroot_array[l][0][1]]);
						if (values->new_optimization && (values->collapse != 2)) {
							/*redo name for collapse*/
							if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values);
							else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
							best_aligns[0]->score=tot_best;											
							}				
						tree_counter=1;
						} 
					else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
						for (k=0;k<ntax-1;k++) {
							best_nodes_new[tree_counter][k][0]=reroot_array[l][k][0];
							best_nodes_new[tree_counter][k][1]=reroot_array[l][k][1];
							}
						/*chec to see if novel*/
						is_unique=1;
						if (!values->new_optimization) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(a_down_buf[reroot_array[l][0][1]],best_aligns[k]);
						else {
							if (values->collapse==2) is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax);
							else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*1,values,collapse_name);
							}
						if (is_unique) {
							if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d (R)\n",tree_counter+1);
							found_better=1;
							best_aligns[tree_counter]=make_align(a_down_buf[reroot_array[l][0][1]]);
							if (values->new_optimization && (values->collapse!=2)) {
											free(best_aligns[tree_counter]->name);
											best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
											assert((int) best_aligns[tree_counter]->name);
											best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name); 
											}
							++tree_counter;
							
							}
						}
					}
				}
			/*copy back*/
				/*copy over new set of best_nodes*/
			/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"about to copy back ");*/
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
			}
		}
	/*reset num trees*/
	cur_num_trees=tree_counter;
	++paren_count;
	}/*found better reswap if found a better one*/
for (i=0;i<paren_count;i++) 	if (!(PARALLEL && values->jackboot)) fprintf(stderr,")");
if (!(PARALLEL && values->jackboot)) fprintf(stderr,"\n");
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
free(is_it_or_desc);
free(blocked);

if ((values->arrt) || (values->atbr)) {
	for (i=0;i<2*ntax-3;i++) {
		for (j=0;j<ntax-1;j++) free(reroot_array[i][j]);
		free(reroot_array[i]);
		}
	free(reroot_array);
	}

/*filter*/
if (values->new_optimization) {
	if (values->aquick || (!values->aquick && mult)) { /*only filter at end*/
		if (values->collapse==2) {
			for (i=0;i<tree_counter;i++) {
				if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
				if (!(PARALLEL*parallel_modify)) check=all_diagnose_tree_local_new_opt(a_down,best_nodes_new[i],ntax,ntax,values,0,corres); 
				else check=all_diagnose_tree_local_parallel_new_opt(a_down,best_nodes_new[i],ntax,ntax,values,0,corres);
				best_aligns[i]=make_align(a_down[best_nodes_new[i][0][1]]);
				if (best_aligns[i]->score != tot_best) if (!(PARALLEL*parallel_modify)) fprintf(stderr,"Score problems!\n");
				}
			}
		}
	}
else {
	if (values->aquick || (!values->aquick && mult)) { /*only filter at end*/
		tot_best=HUGE_COST;
		num_best=0;
		for (i=0;i<tree_counter;i++) {
			if (!(PARALLEL*parallel_modify)) check=all_diagnose_tree_local2(a,best_nodes_new[i],ntax,ntax,values);
			else  check=all_diagnose_tree_local_parallel2(a,best_nodes_new[i],ntax,ntax,values);
			if (compare_aligns(best_aligns[i],a[best_nodes_new[i][0][1]])) {
				if (score_holder>0) {
					values->phylo_score=score_holder;
					a[best_nodes_new[i][0][1]]->score=cladogram(a[best_nodes_new[i][0][1]],values);
					values->phylo_score=0;
					}
				if (a[best_nodes_new[i][0][1]]->score < best_aligns[i]->score) {
					if (!(PARALLEL && values->jackboot)) fprintf(stderr,"	Filtering Alignment %d => New Cost %d\n",i,a[best_nodes_new[i][0][1]]->score);
					if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
					best_aligns[i]=make_align(a[best_nodes_new[i][0][1]]);
					print_paup(a[best_nodes_new[i][0][1]],values);exit(-1);
					}
				}
			/*at any rate fix name so matches order for OPTALI too*/
			free(best_aligns[i]->name);
			best_aligns[i]->name=(char *)malloc((1+strlen(a[nodes[0][1]]->name))*sizeof(char));
			assert((int) best_aligns[i]->name);
			strcpy(best_aligns[i]->name,a[nodes[0][1]]->name);
			if (best_aligns[i]->score < tot_best) tot_best=best_aligns[i]->score;
			}
		if (mult && (tree_counter>1)) {
			best_aligns_holder=(alignment **)malloc(tree_counter*sizeof(alignment *));
			assert((int) best_aligns_holder);
			for (i=0;i<tree_counter;i++) {
				best_aligns_holder[i]=make_align(best_aligns[i]);
				best_aligns[i]=dump_align(best_aligns[i]);
				}
			num_best=0;
			for (i=0;i<tree_counter;i++) if (best_aligns_holder[i]->score==tot_best) {
				is_unique=1;
				for (j=0;j<i;j++) if (best_aligns_holder[j]->score==tot_best) is_unique*=compare_aligns(best_aligns_holder[i],best_aligns_holder[j]);
				if (is_unique) best_aligns[num_best++]=make_align(best_aligns_holder[i]);
				}
			for (i=0;i<tree_counter;i++) best_aligns_holder[i]=dump_align(best_aligns_holder[i]);
			free(best_aligns_holder);
			tree_counter=num_best;
			}
		values->number+=((tree_counter*(ntax-1))/2);
		}
	}
for (i=0;i<values->keep_aligns;i++) {
	for (j=0;j<ntax-1;j++) {
		free(best_nodes_new[i][j]);
		free(best_nodes_old[i][j]);
		}
	free(best_nodes_new[i]);
	free(best_nodes_old[i]);
	}       
free(best_nodes_new);
free(best_nodes_old);
for (i=0;i<2*ntax-1;i++) if (a_down_buf[i]) a_down_buf[i]=dump_align(a_down_buf[i]);
free(a_down_buf);
if (values->new_optimization) {
	for (i=0;i<2*ntax-1;i++) if (a_down[i]) a_down[i]=dump_align(a_down[i]);
	free(a_down);
	for (i=0;i<2*ntax;i++) {
		if (corres[i]) {
			if (corres[i][0]) free(corres[i][0]);
			if (corres[i][1]) free(corres[i][1]);
			free(corres[i]);
			}
		}
	free(corres);
	for (i=0;i<2*ntax;i++) {
		if (corres_buf[i]) {
			if (corres_buf[i][0]) free(corres_buf[i][0]);
			if (corres_buf[i][1]) free(corres_buf[i][1]);
			free(corres_buf[i]);
			}
		}
	free(corres_buf);
	for (i=0;i<2*ntax-1;i++) if (check_buf[i]) check_buf[i]=dump_align(check_buf[i]);
	free(check_buf);
	if (collapse_name) free(collapse_name);
	}
if (values->support) {
	for (i=0;i<values->keep_aligns;i++) free(values->support[i]);
	free(values->support);
	values->support=NULL;
	}
if (values->tbr_align_shortcut) {
	for (i=0;i<2*ntax-1;i++) if (a_buf[i]) a_buf[i]=dump_align(a_buf[i]);
	free(a_buf);
	if (a_up_buf) {
		for (i=0;i<2*ntax-1;i++) {
			if (a_up_buf[i][0]) a_up_buf[i][0]=dump_align(a_up_buf[i][0]);
			if (a_up_buf[i][1]) a_up_buf[i][1]=dump_align(a_up_buf[i][1]);
			free(a_up_buf[i]);
			}
		free(a_up_buf);
		}
	}
if (thing_to_add) thing_to_add=dump_align(thing_to_add);
if (temp_align) temp_align=dump_align(temp_align);
if (temp_align2) temp_align2=dump_align(temp_align2);
*l_in=tree_counter;
/*for (i=0;i<tree_counter;i++) best_aligns[i]->score=tot_best;*/
/*if (!(PARALLEL*parallel_modify)) fprintf(stderr,"Found %d trees\n",num_found);*/
return best_aligns;
}



