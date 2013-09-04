/*Copyright 1994 Ward Wheeler all rights reserved*/
/*ned to add something which will pic individual taxa and move them as well*/
#include "align3.h"

realloc_align_holder *do_node_swap(counter,old_rep,a,best_aligns,nodes,n_nodes,l,values,current,ntaxa,align_passer)
int *counter, **old_rep, **nodes, **n_nodes,*l;
int ntaxa,current;
parameters *values;
alignment **a,**best_aligns;
realloc_align_holder *align_passer;
{
	int i,j,k,m,sub_taxa,found_how_many,ii;
	int cur_val,best_value,**tree_rep;
	char overflow,temp,found_better_or_more;
	int *intermediate_scores,intermediate_best_score, ***intermediate_nodes;
	int ***best_nodes;
	alignment **intermediate_align;
	int num_intermediates_allocated;
	int num_best_aligns, next_to_do, num_actually_done;
	int bufid, bytes, type, source, max_to_do, index, grain=1, num_left, info, n_sub_taxa;
	int ***from_nodes, **new_node, current_from[2], *is_descendent, **yet_newer_node, kk;
	int kkk, found_new, ***node_counter, old_number_best, n_left, n_plucked;
	int **remain_node, **prune_node, **unrooted_prune_node, **block_node, **new_prune_node;
	int root, *is_descendent2, root_node, root_desc, num_spots,jj,lll,ll;
	int *is_descendent3, root_node_root,root_desc_root, num_spots_root;
	int **remain_node_root, **prune_node_root, **unrooted_prune_node_root, **block_node_root;
	int **new_prune_node_root;
	int **new_node_root,**yet_newer_node_root, checker, test_score,jjj;
	int **cur_groups, test, ngroups2;
	int number_sent, number_recieved,i1,j1,i2;
	int doit,single_anc;

	sub_taxa=current+4;
	best_value=intermediate_best_score=best_aligns[0]->score;

	new_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)new_node);
	for (i=0;i<(ntaxa-1);i++) {
		new_node[i]=(int *)malloc(2*sizeof(int));
		assert((int) new_node[i]);
	}
	new_node_root=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)new_node_root);
	for (i=0;i<(ntaxa-1);i++) {
		new_node_root[i]=(int *)malloc(2*sizeof(int));
		assert((int) new_node_root[i]);
	}
	yet_newer_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)yet_newer_node);
	for (i=0;i<(ntaxa-1);i++) {
		yet_newer_node[i]=(int *)malloc(2*sizeof(int));
		assert((int) yet_newer_node[i]);
	}
	remain_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)remain_node);
	prune_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)prune_node);
	unrooted_prune_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)unrooted_prune_node);
	block_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)block_node);
	new_prune_node=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)new_prune_node);
	for (i=0;i<(ntaxa-1);i++) {
		remain_node[i]=(int *)malloc(2*sizeof(int));
		assert((int) remain_node[i]);
		prune_node[i]=(int *)malloc(2*sizeof(int));
		assert((int) prune_node[i]);
		unrooted_prune_node[i]=(int *)malloc(3*sizeof(int));
		assert((int) unrooted_prune_node[i]);
		block_node[i]=(int *)malloc(3*sizeof(int));
		assert((int) block_node[i]);
		new_prune_node[i]=(int *)malloc(2*sizeof(int));
		assert((int) new_prune_node[i]);
	}

	remain_node_root=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)remain_node_root);
	prune_node_root=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)prune_node_root);
	unrooted_prune_node_root=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)unrooted_prune_node_root);
	block_node_root=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)block_node_root);
	new_prune_node_root=(int **)malloc((ntaxa-1)*sizeof(int *));
	assert((int)new_prune_node_root);
	for (i=0;i<(ntaxa-1);i++) {
		remain_node_root[i]=(int *)malloc(2*sizeof(int));
		assert((int) remain_node_root[i]);
		prune_node_root[i]=(int *)malloc(2*sizeof(int));
		assert((int) prune_node_root[i]);
		unrooted_prune_node_root[i]=(int *)malloc(3*sizeof(int));
		assert((int) unrooted_prune_node_root[i]);
		block_node_root[i]=(int *)malloc(3*sizeof(int));
		assert((int) block_node_root[i]);
		new_prune_node_root[i]=(int *)malloc(2*sizeof(int));
		assert((int) new_prune_node_root[i]);
	}

	is_descendent=(int *)malloc(((2*ntaxa)-1)*sizeof(int));
	assert((int) is_descendent);
	is_descendent2=(int *)malloc(((2*ntaxa)-1)*sizeof(int));
	assert((int) is_descendent2);
	is_descendent3=(int *)malloc(((2*ntaxa)-1)*sizeof(int));
	assert((int) is_descendent3);
	for (i=1;i<((2*ntaxa)-1);i++) is_descendent3[i]=1;
	is_descendent3[0]=0;

	num_intermediates_allocated=1;
	intermediate_scores=(int *)malloc(num_intermediates_allocated*sizeof(int));
	assert((int) intermediate_scores);
	intermediate_align=(alignment **)malloc(num_intermediates_allocated*sizeof(alignment *));
	assert((int) intermediate_align);
	intermediate_nodes=(int ***)malloc(num_intermediates_allocated*sizeof(int **));
	assert((int) intermediate_nodes);
	for (i=0;i<num_intermediates_allocated;i++) {
		intermediate_align[i]=NULL;
		intermediate_nodes[i]=(int **)malloc((ntaxa-1)*sizeof(int *));
		assert((int) intermediate_nodes[i]);
		for (j=0;j<(ntaxa-1);j++) {
			intermediate_nodes[i][j]=(int *)malloc(2*sizeof(int));
			assert((int) intermediate_nodes[i][j]);
		}
	}

	best_nodes=(int ***)malloc(values->keep_aligns*sizeof(int **));
	assert((int) best_nodes);
	for (i=0;i<values->keep_aligns;i++) {
		best_nodes[i]=(int **)malloc((ntaxa-1)*sizeof(int *));
		assert((int) best_nodes[i]);
		for (j=0;j<(ntaxa-1);j++) {
			best_nodes[i][j]=(int *)malloc(2*sizeof(int));
			assert((int) best_nodes[i][j]);
		}
	}

	node_counter=(int ***)malloc(values->keep_aligns*sizeof(int **));
	assert((int) node_counter);
	for (i=0;i<values->keep_aligns;i++) {
		node_counter[i]=(int **)malloc((ntaxa-1)*sizeof(int *));
		assert((int) node_counter[i]);
		for (j=0;j<(ntaxa-1);j++) {
			node_counter[i][j]=(int *)malloc(2*sizeof(int));
			assert((int) node_counter[i][j]);
		}
	}

	if ((values->ngroups>0)	&& (values->groups_as_start)) {
		for (i=0;i<values->ngroups;i++) free(values->groups[i]);
		free(values->groups);
		values->groups=NULL;
		values->ngroups=0;
		fprintf(stderr,"	Freeing groups\n");
		if (PARALLEL) pvm_mcast(&(values->tids[1]),values->num_hosts-1, PVM_DUMP_INITIAL_GROUPS);
	}

	for (i=0;i<=(*counter);i++) all_other_make_nodes3(old_rep[i],ntaxa,sub_taxa,best_nodes[i]);
	old_number_best=0;
	found_how_many=0;
node_top:
	;
	if (values->VERBOSE) {
		fprintf(stderr,"(node swapping");
		fflush(stderr);
	}
	overflow=found_better_or_more=found_new=0;
	/*nodes_to_swap=0;*/
	for (i=0;i<=(*counter);i++) copy_nodes(ntaxa,sub_taxa,best_nodes[i],node_counter[i]);
	for (i=old_number_best;i<=(*counter);i++) {/*swap on input tree_reps*/
		copy_nodes(ntaxa,sub_taxa,node_counter[i],n_nodes);
		for (j=0;j<num_intermediates_allocated;j++) if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
		ii=-1;
		if (PARALLEL) number_sent=number_recieved=0;
		for (j=1;j<sub_taxa;j++) {
			n_left=sub_taxa-1;
			n_plucked=1;
			/*for each node > 0 with remaing taxa >1 pick node and create new nodes description*/
			/*these next two yield TBR root after adding the swanstucker (maybe get all any way and more efficient)
			could do rerootings of remainder - better to do after addition to get full roots of pruned bit
			could do rerootings of pruned bit*/
			for (k=0;k<sub_taxa-1;k++) {
				if (n_nodes[k][0]==j) single_anc=k;
				else if (n_nodes[k][1]==j) single_anc=k;
			}
			pick_node_and_remake_with_from(ntaxa,sub_taxa,n_nodes,j,new_node,current_from);
			for (k=0;k<sub_taxa-1;k++)  for (kk=0;kk<2;kk++) {
				/*generate new node insertions*/
				if (!((k==0) && (kk==0))) if (!((k==current_from[0]) && (kk==current_from[1]))) if (k!=single_anc) {
					/*fprintf(stderr,"New node:\n"); print_nodes(ntaxa,sub_taxa,new_node);
					fprintf(stderr,"CF0=%d CF1=%d\n",current_from[0],current_from[1]);*/
					get_new_topo( ntaxa,sub_taxa,j,new_node,yet_newer_node,k, kk);
					/*fprintf(stderr,"- ");
					fprintf(stderr,"Yet newer:\n"); print_nodes(ntaxa,sub_taxa,yet_newer_node);*/
					checker=1;
					if (!checker) {
						fprintf(stderr,"AAAH! problems!\n");
						printf("AAAH! problems!\n");
						printf("Unrooted prune node root\n");
						print_nodes3_stdout(ntaxa,sub_taxa,unrooted_prune_node_root);
						printf("New prune node root\n");
						print_nodes3_stdout(ntaxa,sub_taxa,new_prune_node_root);
						printf("New node root\n");
						print_nodes_stdout(ntaxa,sub_taxa,new_node_root);
						printf("Yet Newer node\n");
						print_nodes_stdout(ntaxa,sub_taxa,yet_newer_node);
					}
					else {
						/*get alignment*/
						if (PARALLEL && !SPECIAL) {
							/*insert topology test here if (!values->groups) */
							doit=1;
							if (doit) {
								if (number_sent < (values->num_hosts-1)) {
									/* code folded from here */
									/*send out job*/
									pvm_initsend(PvmDataDefault);
									pvm_pkint(&grain,1,1);
									pvm_pkint(&ntaxa,1,1);
									for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_pkint(&(yet_newer_node[i1][j1]),1,1);
									pvm_send(values->tids[number_sent+1],PVM_REQUEST_DIAGNOSE_ALIGNMENT);
									++number_sent;

									/* unfolding */
								}
								else {
									/* code folded from here */
									/*blocking recieve*/
									bufid=pvm_recv(-1,PVM_SEND_DIAGNOSE_ALIGNMENT);
									if (bufid) {
										info=pvm_bufinfo(bufid, &bytes, &type, &source);
										/*send out next job*/
										pvm_initsend(PvmDataDefault);
										pvm_pkint(&grain,1,1);
										pvm_pkint(&ntaxa,1,1);
										for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_pkint(&(yet_newer_node[i1][j1]),1,1);
										pvm_send(source,PVM_REQUEST_DIAGNOSE_ALIGNMENT);
										++number_sent;
										/*unpack*/
										++number_recieved;
										pvm_upkint(&grain,1,1);
										pvm_upkint(&test_score,1,1);
										/*evaluate reallocate ii etc*/
										if ((test_score<=intermediate_best_score) && (test_score!=HUGE_COST)) {
											if (test_score<intermediate_best_score) {
												if (ii>=num_intermediates_allocated-1) {
													++num_intermediates_allocated;
													intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
													assert((int) intermediate_scores);
													intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
													assert((int) intermediate_align);
													intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
													assert((int) intermediate_nodes);
													intermediate_align[num_intermediates_allocated-1]=NULL;
													intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
													assert((int) intermediate_nodes[num_intermediates_allocated-1]);
													for (jjj=0;jjj<(ntaxa-1);jjj++) {
														intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
														assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
													}
												}
												intermediate_scores[++ii]=test_score;
												intermediate_best_score=intermediate_scores[ii];
												if (values->VERBOSE) fprintf(stderr,"Found better at cost %d.\n",intermediate_best_score);
												if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
												intermediate_align[ii]=unpack_align(test_score,values);
												/*unpack nodes*/
												for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(intermediate_nodes[ii][i1][j1]),1,1);
											}
											else if (!values->aquick) {
												if (ii>=num_intermediates_allocated-1) {
													++num_intermediates_allocated;
													intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
													assert((int) intermediate_scores);
													intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
													assert((int) intermediate_align);
													intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
													assert((int) intermediate_nodes);
													intermediate_align[num_intermediates_allocated-1]=NULL;
													intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
													assert((int) intermediate_nodes[num_intermediates_allocated-1]);
													for (jjj=0;jjj<(ntaxa-1);jjj++) {
														intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
														assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
													}
												}
												intermediate_scores[++ii]=test_score;
												if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
												intermediate_align[ii]=unpack_align(test_score,values);
												/*unpack_nodes*/
												for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(intermediate_nodes[ii][i1][j1]),1,1);
											}
										}
									}
									/* unfolding */
								}
							}
						}
						else { /*not parallel*/
							/* insert topology test here */
							if (!values->groups) test_score=all_diagnose_tree(a,yet_newer_node,
							    ntaxa,sub_taxa,values);
							else {
								cur_groups=get_cur_groups(yet_newer_node,ntaxa,sub_taxa,&ngroups2);
								test=compare_groups(values->groups,cur_groups,values->ngroups,
								    ngroups2,ntaxa-1,sub_taxa-1);
								/*dump groups*/
								for (i=0;i<ngroups2;i++) free(cur_groups[i]);
								free(cur_groups);
								if (!test) test_score=HUGE_COST;
								else test_score=all_diagnose_tree(a,yet_newer_node,ntaxa,
								    sub_taxa,values);
							}
							/*fprintf(stderr,"\n");*/
							if (test_score<=intermediate_best_score) {
								if (test_score<intermediate_best_score) {
									/* code folded from here */
									if (ii>=num_intermediates_allocated-1) {
										++num_intermediates_allocated;
										intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
										assert((int) intermediate_scores);
										intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
										assert((int) intermediate_align);
										intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
										assert((int) intermediate_nodes);
										intermediate_align[num_intermediates_allocated-1]=NULL;
										intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
										assert((int) intermediate_nodes[num_intermediates_allocated-1]);
										for (jjj=0;jjj<(ntaxa-1);jjj++) {
											intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
											assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
										}
									}
									intermediate_scores[++ii]=test_score;
									intermediate_best_score=intermediate_scores[ii];
									if (values->VERBOSE) fprintf(stderr,"Found better at cost %d.\n",intermediate_best_score);
									if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
									intermediate_align[ii]=make_align(a[yet_newer_node[0][1]]);
									copy_nodes(ntaxa,sub_taxa,yet_newer_node,intermediate_nodes[ii]);
									/* unfolding */
								}
								else if (!values->aquick) {
									/* code folded from here */
									if (ii>=num_intermediates_allocated-1) {
										++num_intermediates_allocated;
										intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
										assert((int) intermediate_scores);
										intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
										assert((int) intermediate_align);
										intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
										assert((int) intermediate_nodes);
										intermediate_align[num_intermediates_allocated-1]=NULL;
										intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
										assert((int) intermediate_nodes[num_intermediates_allocated-1]);
										for (jjj=0;jjj<(ntaxa-1);jjj++) {
											intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
											assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
										}
									}
									intermediate_scores[++ii]=test_score;
									if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
									intermediate_align[ii]=make_align(a[yet_newer_node[0][1]]);
									copy_nodes(ntaxa,sub_taxa,yet_newer_node,intermediate_nodes[ii]);
									/* unfolding */
								}
							}
						}
					}
					/* unfolding */
					/* unfolding */
					/* unfolding */
					/* unfolding */

				}
			}/*TBR loop*/
		}/*picked j*/
		/*nodes with more than a single descendent taxon*/
		for (j=1+ntaxa;j<sub_taxa-1+ntaxa;j++) {
			n_left=remainder(n_nodes,ntaxa,sub_taxa,j-ntaxa,is_descendent,is_descendent2);
			n_plucked=sub_taxa-n_left;
			if (n_left>1) {
				/*for each node > 0 with remaing taxa >1 pick node and create new nodes description*/
				/*these next two yield TBR root after adding the swanstucker (maybe get all any way and more efficient)
			could do rerootings of remainder - better to do after addition to get full roots of pruned bit
			could do rerootings of pruned bit*/
				pick_node_and_remake_with_from(ntaxa,sub_taxa,n_nodes,j,new_node,current_from);
				if ((values->align_root_swap) && (n_plucked>2)) {
					split_nodes(ntaxa,sub_taxa,new_node,remain_node,prune_node,is_descendent2);
					unroot_pruned(ntaxa,sub_taxa,prune_node,unrooted_prune_node,is_descendent2,&root_node,
					    &root_desc);
					block_redundant_additions(ntaxa,sub_taxa,unrooted_prune_node,block_node,is_descendent2,
					    root_node,root_desc);
					num_spots=0;
					for (jj=0;jj<sub_taxa-1;jj++) for (k=0;k<3;k++) if (!block_node[jj][k]) ++num_spots;
					if (num_spots==0) num_spots=1;
				}
				else num_spots=1;
				for (jj=0;jj<num_spots;jj++) {
					if ((values->align_root_swap) && (n_plucked>2) && (num_spots>1)) {
						pick_new_root(ntaxa,sub_taxa,unrooted_prune_node,block_node,jj+1,new_prune_node,
						    is_descendent2,root_node,root_desc);
						/*copy_nodes(ntaxa,sub_taxa,prune_node,new_prune_node);*/
						remake_nodes(ntaxa,sub_taxa,remain_node,new_prune_node,new_node,is_descendent2);
						/*fprintf(stderr,"Remade:\n");
						print_nodes(ntaxa,sub_taxa,new_node);*/
					}
					for (k=0;k<sub_taxa-1;k++)  {
						for (kk=0;kk<2;kk++) if ((!is_descendent[new_node[k][kk]])) {
							/*generate new node insertions*/
							if (!((k==0) && (kk==0))) if (!((k==current_from[0]) && (current_from[1]!=kk))) {
								/*fprintf(stderr,"New node:\n");
								print_nodes(ntaxa,sub_taxa,new_node);*/
								get_new_topo( ntaxa,sub_taxa,j,new_node,yet_newer_node,k,
								    kk);
								/*fprintf(stderr,"Yet newer:\n");
								print_nodes(ntaxa,sub_taxa,yet_newer_node);*/
								if ((values->align_complete_root_swap) || (values->align_partial_root_swap)) {
									/* code folded from here */
									split_nodes(ntaxa,sub_taxa,yet_newer_node,remain_node_root,
									    prune_node_root,is_descendent3);
									unroot_pruned(ntaxa,sub_taxa,prune_node_root,unrooted_prune_node_root,is_descendent3,&root_node_root,&root_desc_root);
									if (!values->align_partial_root_swap) block_redundant_additions(ntaxa,sub_taxa,unrooted_prune_node_root,block_node_root,is_descendent3,root_node_root,root_desc_root);
									else block_redundant_additions_and_non_pruned(ntaxa,sub_taxa,unrooted_prune_node_root,block_node_root,is_descendent3,root_node_root,root_desc_root,is_descendent);
									num_spots_root=0;
									for (ll=0;ll<sub_taxa-1;ll++) for (lll=0;lll<3;lll++) if (block_node_root[ll][lll]==0) ++num_spots_root;
									if (num_spots_root==0) num_spots_root=1;
									/* unfolding */
								}
								else num_spots_root=1;
								if (values->align_partial_root_swap) ++num_spots_root;/*this so do original root + the rootings in the pruned group*/
								/*if (j<ntaxa) num_spots_root=0;*/
								for (ll=0;ll<num_spots_root;ll++) {
									/* code folded from here */
									if (!((values->align_partial_root_swap) && (ll==(num_spots_root-1)))) {
										/*put loop here to reroot each addition=rooted TBR*/
										if (values->align_complete_root_swap && (num_spots_root>1)) {
											pick_new_root(ntaxa,sub_taxa,unrooted_prune_node_root,block_node_root,ll+1,new_prune_node_root,is_descendent3,root_node_root,root_desc_root);
											remake_nodes(ntaxa,sub_taxa,remain_node_root,new_prune_node_root,new_node_root,is_descendent3);
											/*print_nodes(ntaxa,sub_taxa,new_node_root);*/
											copy_nodes(ntaxa,sub_taxa,new_node_root,yet_newer_node);
											yet_newer_node[0][0]=0;
										}
									}


									checker=1;

									if (!checker) {
										fprintf(stderr,"AAAH! problems!\n");
										printf("AAAH! problems!\n");
										printf("Unrooted prune node root\n");
										print_nodes3_stdout(ntaxa,sub_taxa,unrooted_prune_node_root);
										printf("New prune node root\n");
										print_nodes3_stdout(ntaxa,sub_taxa,new_prune_node_root);
										printf("New node root\n");
										print_nodes_stdout(ntaxa,sub_taxa,new_node_root);
										printf("Yet Newer node\n");
										print_nodes_stdout(ntaxa,sub_taxa,yet_newer_node);
									}
									else {
										/*get alignment*/
										if (PARALLEL) {
											/* code folded from here */
											doit=1;
											/* code folded from here */
											if (doit) {
												if (number_sent < (values->num_hosts-1)) {
													/*send out job*/
													pvm_initsend(PvmDataDefault);
													pvm_pkint(&grain,1,1);
													pvm_pkint(&ntaxa,1,1);
													for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_pkint(&(yet_newer_node[i1][j1]),1,1);
													pvm_send(values->tids[number_sent+1],PVM_REQUEST_DIAGNOSE_ALIGNMENT);
													++number_sent;

												}
												else {
													/*blocking recieve*/
													bufid=pvm_recv(-1,PVM_SEND_DIAGNOSE_ALIGNMENT);
													if (bufid) {
														info=pvm_bufinfo(bufid, &bytes, &type, &source);
														/*send out next job*/
														pvm_initsend(PvmDataDefault);
														pvm_pkint(&grain,1,1);
														pvm_pkint(&ntaxa,1,1);
														for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_pkint(&(yet_newer_node[i1][j1]),1,1);
														pvm_send(source,PVM_REQUEST_DIAGNOSE_ALIGNMENT);
														++number_sent;
														/*unpack*/
														++number_recieved;
														pvm_upkint(&grain,1,1);
														pvm_upkint(&test_score,1,1);
														/*evaluate reallocate ii etc*/
														if ((test_score<=intermediate_best_score) && (test_score!=HUGE_COST)) {
															if (test_score<intermediate_best_score) {
																if (ii>=num_intermediates_allocated-1) {
																	++num_intermediates_allocated;
																	intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
																	assert((int) intermediate_scores);
																	intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
																	assert((int) intermediate_align);
																	intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
																	assert((int) intermediate_nodes);
																	intermediate_align[num_intermediates_allocated-1]=NULL;
																	intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
																	assert((int) intermediate_nodes[num_intermediates_allocated-1]);
																	for (jjj=0;jjj<(ntaxa-1);jjj++) {
																		intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
																		assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
																	}
																}
																intermediate_scores[++ii]=test_score;
																intermediate_best_score=intermediate_scores[ii];
																if (values->VERBOSE) fprintf(stderr,"Found better at cost %d.\n",intermediate_best_score);
																if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
																intermediate_align[ii]=unpack_align(test_score,values);
																/*unpack nodes*/
																for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(intermediate_nodes[ii][i1][j1]),1,1);
															}
															else if (!values->aquick) {
																if (ii>=num_intermediates_allocated-1) {
																	++num_intermediates_allocated;
																	intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
																	assert((int) intermediate_scores);
																	intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
																	assert((int) intermediate_align);
																	intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
																	assert((int) intermediate_nodes);
																	intermediate_align[num_intermediates_allocated-1]=NULL;
																	intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
																	assert((int) intermediate_nodes[num_intermediates_allocated-1]);
																	for (jjj=0;jjj<(ntaxa-1);jjj++) {
																		intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
																		assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
																	}
																}
																intermediate_scores[++ii]=test_score;
																if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
																intermediate_align[ii]=unpack_align(test_score,values);
																/*unpack_nodes*/
																for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(intermediate_nodes[ii][i1][j1]),1,1);
															}
														}
													}
												}
											}
										}
										else {
											/* insert topology test here */
											if (!values->groups) test_score=all_diagnose_tree(a,yet_newer_node,ntaxa,sub_taxa,values);
											else {
												cur_groups=get_cur_groups(yet_newer_node,ntaxa,sub_taxa,&ngroups2);
												test=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntaxa-1,sub_taxa-1);
												/*dump groups*/
												for (i=0;i<ngroups2;i++) free(cur_groups[i]);
												free(cur_groups);
												if (!test) test_score=HUGE_COST;
												else test_score=all_diagnose_tree(a,yet_newer_node,ntaxa,sub_taxa,values);
											}
											/*fprintf(stderr,"\n");*/
											if (test_score<=intermediate_best_score) {
												if (test_score<intermediate_best_score) {
													if (ii>=num_intermediates_allocated-1) {
														++num_intermediates_allocated;
														intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
														assert((int) intermediate_scores);
														intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
														assert((int) intermediate_align);
														intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
														assert((int) intermediate_nodes);
														intermediate_align[num_intermediates_allocated-1]=NULL;
														intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
														assert((int) intermediate_nodes[num_intermediates_allocated-1]);
														for (jjj=0;jjj<(ntaxa-1);jjj++) {
															intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
															assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
														}
													}
													intermediate_scores[++ii]=test_score;
													intermediate_best_score=intermediate_scores[ii];
													if (values->VERBOSE) fprintf(stderr,"Found better at cost %d.\n",intermediate_best_score);
													if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
													intermediate_align[ii]=make_align(a[yet_newer_node[0][1]]);
													copy_nodes(ntaxa,sub_taxa,yet_newer_node,intermediate_nodes[ii]);
												}
												else if (!values->aquick) {
													if (ii>=num_intermediates_allocated-1) {
														++num_intermediates_allocated;
														intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
														assert((int) intermediate_scores);
														intermediate_align=(alignment **)realloc(intermediate_align,num_intermediates_allocated*sizeof(alignment *));
														assert((int) intermediate_align);
														intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int **));
														assert((int) intermediate_nodes);
														intermediate_align[num_intermediates_allocated-1]=NULL;
														intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
														assert((int) intermediate_nodes[num_intermediates_allocated-1]);
														for (jjj=0;jjj<(ntaxa-1);jjj++) {
															intermediate_nodes[num_intermediates_allocated-1][jjj]=(int *)malloc(2*sizeof(int));
															assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
														}
													}
													intermediate_scores[++ii]=test_score;
													if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
													intermediate_align[ii]=make_align(a[yet_newer_node[0][1]]);
													copy_nodes(ntaxa,sub_taxa,yet_newer_node,intermediate_nodes[ii]);
												}
											}
										}
									}
									/* unfolding */
									/* unfolding */
									/* unfolding */
									/* unfolding */
									/* unfolding */
								}
							}
						} /*if (!descendent)*/
					}/*TBR loop*/
				}
			}
		}/*picked j*/
		if (PARALLEL) {
			/*get the reaminder of jobs not recieved (this should always be values->num_hosts-2)*/
			while (number_recieved< number_sent) {
				/*blocking recieve*/
				bufid=pvm_recv(-1,PVM_SEND_DIAGNOSE_ALIGNMENT);
				if (bufid) {
					info=pvm_bufinfo(bufid, &bytes, &type, &source);
					/*unpack*/
					++number_recieved;
					pvm_upkint(&grain,1,1);
					pvm_upkint(&test_score,1,1);
					/*evaluate reallocate ii etc*/
					if ((test_score<=intermediate_best_score) && (test_score!=HUGE_COST)) {
						if (test_score<intermediate_best_score) {
							if (ii>=num_intermediates_allocated-1) {
								++num_intermediates_allocated;
								intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
								assert((int) intermediate_scores);
								intermediate_align=(alignment **)realloc(intermediate_align,
								    num_intermediates_allocated*sizeof(alignment *));
								assert((int) intermediate_align);
								intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int
								**));
								assert((int) intermediate_nodes);
								intermediate_align[num_intermediates_allocated-1]=NULL;
								intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int
								*));
								assert((int) intermediate_nodes[num_intermediates_allocated-1]);
								for (jjj=0;jjj<(ntaxa-1);jjj++) {
									/* code folded from here */
									/* code folded from here */
									intermediate_nodes[num_intermediates_allocated-1][jjj]=(int
									*)malloc(2*sizeof(int));
									assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
									/* unfolding */
									/* unfolding */
								}
							}
							intermediate_scores[++ii]=test_score;
							intermediate_best_score=intermediate_scores[ii];
							if (values->VERBOSE) fprintf(stderr,"Found better at cost %d.\n",
							    intermediate_best_score);
							if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
							intermediate_align[ii]=unpack_align(test_score,values);
							/*unpack nodes*/
							for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(intermediate_nodes[ii][i1][j1]),
							    1,1);
						}
						else if (!values->aquick) {
							if (ii>=num_intermediates_allocated-1) {
								++num_intermediates_allocated;
								intermediate_scores=(int *)realloc(intermediate_scores,num_intermediates_allocated*sizeof(int));
								assert((int) intermediate_scores);
								intermediate_align=(alignment **)realloc(intermediate_align,
								    num_intermediates_allocated*sizeof(alignment *));
								assert((int) intermediate_align);
								intermediate_nodes=(int ***)realloc(intermediate_nodes,num_intermediates_allocated*sizeof(int
								**));
								assert((int) intermediate_nodes);
								intermediate_align[num_intermediates_allocated-1]=NULL;
								intermediate_nodes[num_intermediates_allocated-1]=(int **)malloc((ntaxa-1)*sizeof(int
								*));
								assert((int) intermediate_nodes[num_intermediates_allocated-1]);
								for (jjj=0;jjj<(ntaxa-1);jjj++) {
									/* code folded from here */
									/* code folded from here */
									intermediate_nodes[num_intermediates_allocated-1][jjj]=(int
									*)malloc(2*sizeof(int));
									assert((int) intermediate_nodes[num_intermediates_allocated-1][jjj]);
									/* unfolding */
									/* unfolding */
								}
							}
							intermediate_scores[++ii]=test_score;
							if (intermediate_align[ii]) intermediate_align[ii]=dump_align(intermediate_align[ii]);
							intermediate_align[ii]=unpack_align(test_score,values);
							/*unpack_nodes*/
							for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(intermediate_nodes[ii][i1][j1]),
							    1,1);
						}
					}
				}/*recieve*/

			}
		}


		/*loop through all alignments and scores to find good ones*/
		num_best_aligns=0;
		num_actually_done=ii+1;
		/*fprintf(stderr,"Num actually done %d\n",num_actually_done);*/
		for (j=0;j<num_actually_done;j++) if (intermediate_scores[j]==intermediate_best_score) ++num_best_aligns;
		if ((best_value>intermediate_best_score) || (num_best_aligns>0)) { /*see if swap did anything*/
			if (intermediate_best_score < best_value) {
				best_value=intermediate_best_score;
				for (k=0;k<(*l);k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
				(*l)=0;
				found_new=1;
			}
			/*change all this rep shit to nodes*/
			for (j=0;j<num_actually_done;j++) if ((intermediate_scores[j]==intermediate_best_score) && (intermediate_align[j])) {
				if (!(*l)) {
					found_better_or_more=1;
					for (k=0;k<(*l);k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
					(*l)=1;
					best_aligns[0]=make_align(intermediate_align[j]);
					/*copy best and from*/
					copy_nodes(ntaxa,sub_taxa,intermediate_nodes[j],best_nodes[0]);
				}/*better*/
else if (!values->aquick) {
	temp=1;
	/*someday change compare aligns to supoprted nodes only difference*/
	for (k=0;k<(*l);k++) if (best_aligns[k]) temp*=compare_aligns(intermediate_align[j],best_aligns[k]);
	if (temp) {
		if ((*l)>(values->keep_aligns-2)) {
			if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment
			    *)))!=NULL) && ((best_nodes=(int ***)realloc(best_nodes,(values->keep_aligns+1)*sizeof(int **)))!=NULL) &&
			    ((node_counter=(int ***)realloc(node_counter,(values->keep_aligns+1)*sizeof(int **)))!=NULL)) {
				++values->keep_aligns;
				best_aligns[values->keep_aligns-1]=NULL;
				best_nodes[values->keep_aligns-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
				assert((int)best_nodes[values->keep_aligns-1]);
				for (k=0;k<ntaxa-1;k++) {
					best_nodes[values->keep_aligns-1][k]=(int *)malloc(2*sizeof(int));
					assert((int) best_nodes[values->keep_aligns-1][k]);
				}
				node_counter[values->keep_aligns-1]=(int **)malloc((ntaxa-1)*sizeof(int *));
				assert((int)node_counter[values->keep_aligns-1]);
				for (k=0;k<ntaxa-1;k++) {
					node_counter[values->keep_aligns-1][k]=(int *)malloc(2*sizeof(int));
					assert((int) node_counter[values->keep_aligns-1][k]);
				}
			}
		}
		if ((*l)>(values->keep_aligns-2)) {
			if (!overflow) {
				if (values->rep_error) fprintf(stderr,"Overflow in alignment node swap -- but that's OK.\n");
				overflow=1;
			}
		}
		else {
			copy_nodes(ntaxa,sub_taxa,intermediate_nodes[j],best_nodes[(*l)]);
			best_aligns[(*l)]=make_align(intermediate_align[j]);
			++(*l);
			found_better_or_more=1;
			if (values->VERBOSE) fprintf(stderr,"Found another for %d alignments\n",(*l));
		}
	}
}/*additional*/
			}
		}
	}
	if (found_new) old_number_best=0;
	else old_number_best=(*counter);
	(*counter)=(*l)-1;
	found_how_many+=found_better_or_more;
	if (found_better_or_more) goto node_top;

	for (i=0;i<(ntaxa-1);i++) {
		free(new_prune_node[i]);
		free(new_node[i]);
		free(yet_newer_node[i]);
		free(remain_node[i]);
		free(prune_node[i]);
		free(unrooted_prune_node[i]);
		free(block_node[i]);
	}
	free(yet_newer_node);
	free(new_node);
	free(remain_node);
	free(prune_node);
	free(unrooted_prune_node);
	free(block_node);
	free(new_prune_node);
	for (i=0;i<(ntaxa-1);i++) {
		free(new_node_root[i]);
		free(remain_node_root[i]);
		free(prune_node_root[i]);
		free(unrooted_prune_node_root[i]);
		free(block_node_root[i]);
		free(new_prune_node_root[i]);
	}
	free(new_node_root);
	free(remain_node_root);
	free(prune_node_root);
	free(unrooted_prune_node_root);
	free(block_node_root);
	free(new_prune_node_root);


	free(is_descendent);
	free(is_descendent2);
	free(is_descendent3);

	free(intermediate_scores);
	for (i=0;i<num_intermediates_allocated;i++) {
		if (intermediate_align[i]) intermediate_align[i]=dump_align(intermediate_align[i]);
		for (j=0;j<(ntaxa-1);j++) free(intermediate_nodes[i][j]);
		free(intermediate_nodes[i]);
	}
	free(intermediate_align);
	free(intermediate_nodes);
	for (i=0;i<values->keep_aligns;i++) for (j=0;j<ntaxa-1;j++) {
		free(best_nodes[i][j]);
		free(node_counter[i][j]);
	}
	free(best_nodes);
	free(node_counter);
	if (values->VERBOSE) {
		for (i=0;i<found_how_many;i++) fprintf(stderr,")");
		fprintf(stderr,") ");
		fflush(stderr);
	}
	align_passer->best=best_aligns;
	align_passer->rep=old_rep;

	return align_passer;
}

void all_other_make_nodes3(tree_rep,ntaxa,nent,node)
int **node;
int *tree_rep;
int ntaxa,nent;
{
	int i,dummy1,dummy2;

	/*node names start after the taxa*/
	/*initialize nodes*/
	node[0][0]=0;
	node[0][1]=ntaxa+1;
	node[1][0]=1;
	node[1][1]=2;

	for (i=3;i<nent;++i)
	{
		node[i-1][0]=i;
		dummy1=(tree_rep[i-3])/2;
		dummy2=(tree_rep[i-3])%2;
		node[i-1][1]=node[dummy1][dummy2];
		node[dummy1][dummy2]=ntaxa+i-1;
	}
}

/*function to determine the number of taxa left after a pruning*/
int remainder(nodes,ntaxa,sub_taxa,node_taken,is_descendent,is_descendent2)
int node_taken, ntaxa, **nodes, *is_descendent, *is_descendent2;
{
	int pruned,i;
	int got_one, anc;

	/*fprintf(stderr,"Getting remainder ");*/
	/*how many descendents of node_taken?*/
	for (i=0;i<((2*ntaxa)-1);i++) is_descendent[i]=0;
	for (i=0;i<sub_taxa-1;i++) {
		if (nodes[i][0]==(node_taken+ntaxa)) anc=i;
		else if (nodes[i][1]==(node_taken+ntaxa)) anc=i;
	}
	is_descendent[node_taken+ntaxa]=1;
	is_descendent[nodes[node_taken][0]]=1;
	is_descendent[nodes[node_taken][1]]=1;
	/*fprintf(stderr,"got fisrt ");*/
	got_one=1;
	while (got_one) {
		got_one=0;
		for (i=0;i<sub_taxa-1;i++){
			/*fprintf(stderr,"%d ",i);*/
			if (is_descendent[i+ntaxa]) {
				if (!is_descendent[nodes[i][0]]) {
					got_one=1;
					is_descendent[nodes[i][0]]=1;
					is_descendent[nodes[i][1]]=1;
				}
			}
		}
	}
	/*fprintf(stderr,"done getting ");*/

	/*for tbr stuff thats why not snc*/
	is_descendent[anc+ntaxa]=1;
	for (i=0;i<((2*ntaxa)-1);i++) is_descendent2[i]=is_descendent[i];
	is_descendent[nodes[anc][0]]=1;
	is_descendent[nodes[anc][1]]=1;

	pruned=0;
	for (i=0;i<sub_taxa;i++) pruned += is_descendent[i];
	/*fprintf(stderr,"p-%d ",pruned);*/
	return (sub_taxa-pruned);
}

void nodes_to_tree_rep(tree_rep,ntaxa,nent,node)
int **node;
int *tree_rep;
int ntaxa,nent;
{
	int i,cur_i,cur_j,j,k;
	int anc_i,anc_j, rep_count;
	int temp_0,temp_1;

	/*node names start after the taxa*/

	/*for (i=3;i<nent;i++) {
		node[i-1][0]=i;
		dummy1=(tree_rep[i-3])/2;
		dummy2=(tree_rep[i-3])%2;
		node[i-1][1]=node[dummy1][dummy2];
		node[dummy1][dummy2]=ntaxa+i-1;
		}
	*/

	for (i=1;i<nent-1;i++) {
		temp_0=node[i][0];
		temp_1=node[i][1];
		if (((temp_0<ntaxa) && (temp_1<ntaxa)) || ((temp_0>=ntaxa) && (temp_1>=ntaxa))) {
			node[i][1]=min(temp_0,temp_1);
			node[i][0]=max(temp_0,temp_1);
		}
	}
	/*by taxa*/
	for (i=nent-1;i>=3;i--) {
		for (j=0;j<nent-1;j++) for (k=0;k<2;k++) if (node[j][k]==i) {
			cur_i=j;
			cur_j=k;
			goto next;
		}
next:
		;
		for (j=0;j<nent-1;j++) for (k=0;k<2;k++) if (node[j][k]==(ntaxa+cur_i)) {
			anc_i=j;
			anc_j=k;
			goto still_next;
		}
still_next:
		;
		node[anc_i][anc_j]=node[cur_i][!cur_j];
		node[cur_i][cur_j]=HUGE_COST;
		node[cur_i][!cur_j]=HUGE_COST;
		rep_count=0;
		for (j=0;j<nent-1;j++) for (k=0;k<2;k++) if (!((j==0) && (k==0))) {
			if (node[j][k]<HUGE_COST) {
				++rep_count;
				if ((j==anc_i) && (k==anc_j)) goto still_more_next;
			}
		}
still_more_next:
		;
		tree_rep[i-3]=rep_count;
		fprintf(stderr,"Iter %d\n",i);
		print_nodes(ntaxa,nent,node);
	}
	for (i=0;i<nent-3;i++) fprintf(stderr,"tr[%d]=%d ",i,tree_rep[i]);
	/*print_nodes(ntaxa,nent,node);*/
	fprintf(stderr,"\n");
}

void	pick_node_and_remake_with_from(ntaxa,sub_taxa,nodes,node_taken,new_node,from)
int ntaxa, sub_taxa, **nodes, node_taken, **new_node, *from;
{
	int i, one_taken, anc_to_taken, one_not_taken;

	/*copy to ne_node*/
	for (i=0;i<(sub_taxa-1);i++) {
		new_node[i][0]=nodes[i][0];
		new_node[i][1]=nodes[i][1];
	}

	for (i=0;i<sub_taxa-1;i++) {
		if ((new_node[i][0]==(node_taken)) || (new_node[i][1]==(node_taken))) {
			if (new_node[i][0]==(node_taken)) {
				one_taken=0;
				one_not_taken=1;
			}
			else if (new_node[i][1]==(node_taken)) {
				one_taken=1;
				one_not_taken=0;
			}
			anc_to_taken=i;
		}
	}

	/*move the sister of the node taken to its grandparent*/
	/*establish where the pruned came from*/
	for (i=0;i<sub_taxa-1;i++) {
		if (new_node[i][0]==(anc_to_taken+ntaxa)) {
			new_node[i][0]=new_node[anc_to_taken][one_not_taken];
			from[0]=i;
			from[1]=0;
		}
		else if (new_node[i][1]==(anc_to_taken+ntaxa)) {
			new_node[i][1]=new_node[anc_to_taken][one_not_taken];
			from[0]=i;
			from[1]=1;
		}
	}

}

void get_new_topo(ntaxa,sub_taxa,j,new_node,yet_newer_node,node_to_be_altered,branch)
int ntaxa, sub_taxa,j,**new_node,**yet_newer_node,node_to_be_altered,branch;
{
	int i,k, pruned_ancestor;

	/*fprintf(stderr,"In topo ");*/
	/*copy over*/
	for (i=0;i<sub_taxa-1;i++) for (k=0;k<2;k++) yet_newer_node[i][k]=new_node[i][k];
	/*get ancestor node to the one pruned*/
	for (i=0;i<sub_taxa-1;i++) for (k=0;k<2;k++) if (new_node[i][k]==(j)) pruned_ancestor=i;
	/*modify the new nodes*/
	/*fprintf(stderr,"Copied and got anc(%d) of node[%d] ",pruned_ancestor,j);*/
	yet_newer_node[pruned_ancestor][0]=j;
	yet_newer_node[pruned_ancestor][1]=yet_newer_node[node_to_be_altered][branch];
	yet_newer_node[node_to_be_altered][branch]=pruned_ancestor+ntaxa;
	/*
fprintf(stderr,"Old ");
print_nodes(ntaxa,sub_taxa,new_node);
fprintf(stderr,"New ");
print_nodes(ntaxa,sub_taxa,yet_newer_node);
*/
}

void 	copy_nodes(ntaxa,sub_taxa,old_nodes,new_nodes)
int ntaxa,sub_taxa,**old_nodes,**new_nodes;
{
	int i;

	for (i=0;i<sub_taxa-1;i++) {
		new_nodes[i][0]=old_nodes[i][0];
		new_nodes[i][1]=old_nodes[i][1];
	}

}

void print_nodes(ntaxa,sub_taxa,nodes)
int ntaxa,sub_taxa,**nodes;
{
	int i,j;

	j=0;
	for (i=0;i<sub_taxa-1;i++) {
		if (nodes[i][0]<HUGE_COST) {
			fprintf(stderr,"N[%d][0]=%d N[%d][1]=%d ",i,nodes[i][0],i,nodes[i][1]);
			++j;
			if (!(j%3)) fprintf(stderr,"\n");
		}
	}
	if (j%3) fprintf(stderr,"\n");

}

void print_nodes3(ntaxa,sub_taxa,nodes)
int ntaxa,sub_taxa,**nodes;
{
	int i,j;

	j=0;
	for (i=0;i<sub_taxa-1;i++) {
		if (nodes[i][0]<HUGE_COST) {
			fprintf(stderr,"N[%d][0]=%d N[%d][1]=%d N[%d][2]=%d ",i,nodes[i][0],i,nodes[i][1],i,nodes[i][2]);
			++j;
			if (!(j%2)) fprintf(stderr,"\n");
		}
	}
	if (j%2) fprintf(stderr,"\n");

}

void print_nodes_stdout(ntaxa,sub_taxa,nodes)
int ntaxa,sub_taxa,**nodes;
{
	int i,j;

	j=0;
	for (i=0;i<sub_taxa-1;i++) {
		if (nodes[i][0]<HUGE_COST) {
			printf("N[%d][0]=%d N[%d][1]=%d ",i,nodes[i][0],i,nodes[i][1]);
			++j;
			if (!(j%3)) printf("\n");
		}
	}
	if (j%3) printf("\n");

}

void print_nodes3_stdout(ntaxa,sub_taxa,nodes)
int ntaxa,sub_taxa,**nodes;
{
	int i,j;

	j=0;
	for (i=0;i<sub_taxa-1;i++) {
		if (nodes[i][0]<HUGE_COST) {
			printf("N[%d][0]=%d N[%d][1]=%d N[%d][2]=%d ",i,nodes[i][0],i,nodes[i][1],i,nodes[i][2]);
			++j;
			if (!(j%2)) printf("\n");
		}
	}
	if (j%2) printf("\n");

}

void split_nodes(ntaxa,sub_taxa,new_node,remain_node,prune_node,is_descendent)
int ntaxa, sub_taxa, *is_descendent;
int **new_node, **remain_node, **prune_node;
{
	int i,k;

	for (i=0;i<sub_taxa-1;i++) for (k=0;k<2;k++) remain_node[i][k]=prune_node[i][k]=-1;

	for (i=0;i<sub_taxa-1;i++) {
		if (is_descendent[i+ntaxa]) {
			prune_node[i][0]=new_node[i][0];
			prune_node[i][1]=new_node[i][1];
		}
		else {
			remain_node[i][0]=new_node[i][0];
			remain_node[i][1]=new_node[i][1];
		}
	}
	/*
fprintf(stderr,"prune:\n");
print_nodes(ntaxa,sub_taxa,prune_node);
fprintf(stderr,"remain:\n");
print_nodes(ntaxa,sub_taxa,remain_node);
*/
}

void unroot_pruned(ntaxa,sub_taxa,prune_node,unrooted_prune_node,is_descendent,root_node,root_desc)
int ntaxa, sub_taxa,*is_descendent, *root_node, *root_desc;
int **unrooted_prune_node, **prune_node;
{
	int i,j,k, num_terminals, thang[2];

	/*basal connection is left in and blocked later so no root attempt is made*/
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<2;j++) unrooted_prune_node[i][j]=prune_node[i][j];
		unrooted_prune_node[i][2]=HUGE_COST;
	}
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<sub_taxa-1;j++)  if (is_descendent[j+ntaxa]) {
			for (k=0;k<2;k++) if (prune_node[j][k]==(i+ntaxa)) unrooted_prune_node[i][2]=(j+ntaxa);
		}
	}
	/*fprintf(stderr,"Bunrooted:\n");
for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
	for (j=0;j<3;j++) fprintf(stderr,"node[%d][%d]=%d ",i,j,unrooted_prune_node[i][j]);
	fprintf(stderr,"\n");
	}*/

	/*move basal node thing*/
	/*forget about root node*/
	/*get_name of descendent*/
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (k=0;k<3;k++) if (unrooted_prune_node[i][k]==HUGE_COST) (*root_node)=i;
	}
	for (k=0;k<3;k++) if (unrooted_prune_node[(*root_node)][k]!=HUGE_COST) if (is_descendent[unrooted_prune_node[(*root_node)][k]]) (*root_desc)=unrooted_prune_node[(*root_node)][k]-ntaxa;
	/*fprintf(stderr,"Root descendent is %d ",(*root_desc));*/
	/*transfer root desc taxa to others*/
	num_terminals=0;
	i=0;
	for (k=0;k<3;k++) if (unrooted_prune_node[(*root_desc)][k]!=(*root_node)+ntaxa) thang[i++]=unrooted_prune_node[(*root_desc)][k];
	for (k=0;k<3;k++) if (unrooted_prune_node[(*root_desc)][k]<ntaxa) ++num_terminals;
	if (num_terminals==0){
		for (k=0;k<3;k++) if (unrooted_prune_node[thang[0]-ntaxa][k]==((*root_desc)+ntaxa)) unrooted_prune_node[thang[0]-ntaxa][k]=thang[1];
		for (k=0;k<3;k++) if (unrooted_prune_node[thang[1]-ntaxa][k]==((*root_desc)+ntaxa)) unrooted_prune_node[thang[1]-ntaxa][k]=thang[0];
	}
	else {
		if (thang[0]>=ntaxa) for (k=0;k<3;k++) if (unrooted_prune_node[thang[0]-ntaxa][k]==((*root_desc)+ntaxa)) unrooted_prune_node[thang[0]-ntaxa][k]=thang[1];
		if (thang[1]>=ntaxa) for (k=0;k<3;k++) if (unrooted_prune_node[thang[1]-ntaxa][k]==((*root_desc)+ntaxa)) unrooted_prune_node[thang[1]-ntaxa][k]=thang[0];
	}
	/*fprintf(stderr,"Th1 %d Th2 %d ",thang[0],thang[1]);*/

	/*for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
	for (k=0;k<3;k++) if (unrooted_prune_node[i][k]==((*root_node)+ntaxa)) unrooted_prune_node[i][k]=HUGE_COST;
	}*/


	/*fprintf(stderr,"unrooted:\n");
for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
	for (j=0;j<3;j++) fprintf(stderr,"ur_node[%d][%d]=%d ",i,j,unrooted_prune_node[i][j]);
	fprintf(stderr,"\n");
	}
						fprintf(stderr,"root node %d\n",*root_node);*/}

/*this function is to only allow rerooting the entire cladogram within the pruned bit*/
void block_redundant_additions_and_non_pruned (ntaxa,sub_taxa,unrooted_prune_node,block_node,is_descendent,root_node,root_desc,
is_descendent2)
int ntaxa, sub_taxa, *is_descendent, root_node, root_desc;
int **block_node, **unrooted_prune_node, *is_descendent2;
{
	int i,j,k,l, is_a_match;

	for (i=0;i<sub_taxa-1;i++) for (j=0;j<3;j++) block_node[i][j]=-1;
	for (j=0;j<3;j++) block_node[root_node][j]=1;
	for (j=0;j<3;j++) block_node[root_desc][j]=1;
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) if (i!=root_node) if (i!=root_desc) {
		for (j=0;j<3;j++) {
			if (unrooted_prune_node[i][j]==HUGE_COST) block_node[i][j]=1;
			else block_node[i][j]=0;
		}
	}
	/*check for overcounted internal nodes
	block the second if the first is not blocked */


	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<3;j++) {
			if (block_node[i][j]<1) {
				if (unrooted_prune_node[i][j]>=ntaxa) {
					for (k=0;k<3;k++) {
						is_a_match=0;
						if (unrooted_prune_node[(unrooted_prune_node[i][j])-ntaxa][k]==(i+ntaxa)) is_a_match=1;
						if (is_a_match) block_node[(unrooted_prune_node[i][j])-ntaxa][k]=1;
					}
				}
			}
		}
	}

	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<3;j++)	if (block_node[i][j]<1) {
			for (k=0;k<sub_taxa-1;k++) if (is_descendent[k+ntaxa]) {
				for (l=0;l<3;l++)	if (block_node[k][l]<1) {
					if ((i!=k) && (j!=l)) {
						if (unrooted_prune_node[i][j]==unrooted_prune_node[k][l]) block_node[k][l]=1;
					}
				}
			}
		}
	}
	for (i=0;i<sub_taxa-1;i++) if (!is_descendent2[i+ntaxa]) for (j=0;j<3;j++) block_node[i][j]=1;
	for (i=0;i<sub_taxa-1;i++) for (j=0;j<3;j++) if (block_node[i][j]==-1) block_node[i][j]=1;
}

void block_redundant_additions(ntaxa,sub_taxa,unrooted_prune_node,block_node,is_descendent,root_node,root_desc)
int ntaxa, sub_taxa, *is_descendent, root_node, root_desc;
int **block_node, **unrooted_prune_node;
{
	int i,j,k,l, is_a_match;

	for (i=0;i<sub_taxa-1;i++) for (j=0;j<3;j++) block_node[i][j]=-1;
	for (j=0;j<3;j++) block_node[root_node][j]=1;
	for (j=0;j<3;j++) block_node[root_desc][j]=1;
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) if (i!=root_node) if (i!=root_desc) {
		for (j=0;j<3;j++) {
			if (unrooted_prune_node[i][j]==HUGE_COST) block_node[i][j]=1;
			else block_node[i][j]=0;
		}
	}
	/*check for overcounted internal nodes
	block the second if the first is not blocked */


	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<3;j++) {
			if (block_node[i][j]<1) {
				if (unrooted_prune_node[i][j]>=ntaxa) {
					for (k=0;k<3;k++) {
						is_a_match=0;
						if (unrooted_prune_node[(unrooted_prune_node[i][j])-ntaxa][k]==(i+ntaxa)) is_a_match=1;
						if (is_a_match) block_node[(unrooted_prune_node[i][j])-ntaxa][k]=1;
					}
				}
			}
		}
	}

	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<3;j++)	if (block_node[i][j]<1) {
			for (k=0;k<sub_taxa-1;k++) if (is_descendent[k+ntaxa]) {
				for (l=0;l<3;l++)	if (block_node[k][l]<1) {
					if ((i!=k) && (j!=l)) {
						if (unrooted_prune_node[i][j]==unrooted_prune_node[k][l]) block_node[k][l]=1;
					}
				}
			}
		}
	}
	for (i=0;i<sub_taxa-1;i++) for (j=0;j<3;j++) if (block_node[i][j]==-1) block_node[i][j]=1;
}

void pick_new_root(ntaxa,sub_taxa,unrooted_prune_node_holder,block_node,root,new_prune_node,is_descendent,root_node,root_desc)
int ntaxa, sub_taxa, root, *is_descendent, root_node,root_desc;
int **unrooted_prune_node_holder, **block_node, **new_prune_node;
{
	int i,j,l,row, column, replace_holder, skip, got_one;
	int num_terminals,k,root_counter, **unrooted_prune_node;

	unrooted_prune_node=(int **)malloc((sub_taxa-1)*sizeof(int *));
	assert((int) unrooted_prune_node);
	for (i=0;i<sub_taxa-1;i++) {
		unrooted_prune_node[i]=(int *)malloc(3*sizeof(int));
		assert((int) unrooted_prune_node[i]);
		for (j=0;j<2;j++) {
			unrooted_prune_node[i][j]=unrooted_prune_node_holder[i][j];
			new_prune_node[i][j]=0;
		}
		unrooted_prune_node[i][2]=unrooted_prune_node_holder[i][2];
	}
	/*get_place*/
	root_counter=0;
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		for (j=0;j<3;j++)	{
			if (!block_node[i][j]) {
				++root_counter;
				if (root_counter==root) {
					row=i;
					column=j;
				}
			}
		}
	}
	/*fprintf(stderr,"New root %d %d ",row, column);*/
	/*set new root*/
	replace_holder=unrooted_prune_node[row][column];
	unrooted_prune_node[root_desc][0]=replace_holder;
	unrooted_prune_node[root_desc][1]=row+ntaxa;
	unrooted_prune_node[root_desc][2]=-1;
	unrooted_prune_node[root_node][0]=root_desc+ntaxa;
	unrooted_prune_node[root_node][1]=root_desc+ntaxa;
	unrooted_prune_node[root_node][2]=-1;
	unrooted_prune_node[row][column]=-1;

	/*set roots to -1*/
	for (i=0;i<sub_taxa-1;i++) if (i!=root_node) if (is_descendent[i+ntaxa]) {
		skip=0;
		for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]==-1) skip=1;
		if (!skip) {
			num_terminals=0;
			for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]<ntaxa) ++num_terminals;
			if (num_terminals==2) {
				for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]>=ntaxa) unrooted_prune_node[i][j]=-1;
			}
		}
	}
	for (i=0;i<sub_taxa-1;i++) if (i!=root_node) if (is_descendent[i+ntaxa]) {
		skip=0;
		for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]==-1) skip=1;
		if (skip) {
			for (k=0;k<sub_taxa-1;k++) if (k!=root_node) if (is_descendent[k+ntaxa]) if (k!=i) {
				for (j=0;j<3;j++) for (l=0;l<3;l++) if (unrooted_prune_node[i][j]==unrooted_prune_node[k][l]) unrooted_prune_node[k][l]=-1;
			}
		}
	}
	/*
fprintf(stderr,"Irooted:\n");
for (i=0;i<sub_taxa-1;i++)  if (is_descendent[i+ntaxa]) {
	for (j=0;j<3;j++) fprintf(stderr,"up_node[%d][%d]=%d ",i,j,unrooted_prune_node[i][j]);
	fprintf(stderr,"\n");
	}*/

	/*root descendents of rooted ancestors*/
	/*zero terminal error here?*/
do_it_again:
	;
	got_one=0;
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		skip=1;
		for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]==-1) skip=0;
		if (!skip) {
			for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]>=ntaxa) {
				for (k=0;k<3;k++) {
					if ((unrooted_prune_node[unrooted_prune_node[i][j]-ntaxa][k]) == i+ntaxa) {
						unrooted_prune_node[unrooted_prune_node[i][j]-ntaxa][k]=-1;
						got_one=1;
					}
				}
			}
		}
	}
	if (got_one) goto do_it_again;

	/*copy over rooted version*/
	for (i=0;i<sub_taxa-1;i++) if (is_descendent[i+ntaxa]) {
		k=0;
		for (j=0;j<3;j++)	if (unrooted_prune_node[i][j]>0) new_prune_node[i][k++]=unrooted_prune_node[i][j];
		if (k>2) {
			fprintf(stderr,"UH OH! ");
			for (j=0;j<3;j++)	fprintf(stderr,"UPN[%d][%d]=%d ",i,j,unrooted_prune_node[i][j]);
			fprintf(stderr,"\n");
		}
	}
	/*fprintf(stderr,"rooted:\n");
for (i=0;i<sub_taxa-1;i++)  if (is_descendent[i+ntaxa]) {
	for (j=0;j<3;j++) fprintf(stderr,"up_node[%d][%d]=%d ",i,j,unrooted_prune_node[i][j]);
	fprintf(stderr,"\n");
	}
fprintf(stderr,"noded:\n");
print_nodes(ntaxa,sub_taxa,new_prune_node);*/
	for (i=0;i<sub_taxa-1;i++) free(unrooted_prune_node[i]);
	free(unrooted_prune_node);
}

void remake_nodes(ntaxa,sub_taxa,remain_node,new_prune_node,new_node,is_descendent)
int ntaxa, sub_taxa, *is_descendent;
int **new_node, **remain_node, **new_prune_node;
{
	int i,k;

	for (i=0;i<sub_taxa-1;i++) {
		if (is_descendent[i+ntaxa]) {
			new_node[i][0]=new_prune_node[i][0];
			new_node[i][1]=new_prune_node[i][1];
		}
		else {
			new_node[i][0]=remain_node[i][0];
			new_node[i][1]=remain_node[i][1];
		}
	}
	/*fprintf(stderr,"new new:\n");
print_nodes(ntaxa,sub_taxa,new_node);*/

}

int dumb_shit(ntaxa,sub_taxa,yet_newer_node)
int ntaxa, sub_taxa, **yet_newer_node;
{
	int i,j, k,l;

	for (i=0;i<sub_taxa-1;i++) for (j=0;j<2;j++) if (yet_newer_node[i][j]>(2*ntaxa)) {
		/*fprintf(stderr,"AAAH! problems!\n");
	print_nodes(ntaxa,sub_taxa,yet_newer_node);
	exit(-1);*/
		return 0;
	}
	for (i=0;i<sub_taxa-1;i++) for (j=0;j<2;j++) {
		for (k=i+1;k<sub_taxa-1;k++) for (l=0;l<2;l++) if (yet_newer_node[i][j]==yet_newer_node[k][l]) {
			/*fprintf(stderr,"AAAH! problems!\n");
		print_nodes(ntaxa,sub_taxa,yet_newer_node);
		exit(-1);*/
			return 0;
		}
	}

	return 1;
}
