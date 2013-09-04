/*Copyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Goloboff*/

#include "align3.h"

alignment **make_heur_align_new(a,ntax,nodes,l,best_aligns,values,parallel_modify)
alignment **a, **best_aligns;
int ntax, *l, **nodes, parallel_modify;
parameters *values;
{
	int i,j,k,ll,value,cur_val,n,m;
	int initial_length,d1,d2;
	int *anc,check;
	int score_holder;
	alignment ***a_up,*temp_align,*temp_align2;
	alignment *cur_align,*cur_align2;
	int d3,best_j;
	int bufid,type,bytes,source,index,info,num_to_do,num_to_receive;
	int **temp_nodes,will_be_OK,**cur_groups,ngroups2;
	int ii,jj,***corres;
	alignment **a_down;
	int d1_thang, d2_thang;

	if (values->groups) {
		temp_nodes=(int **)malloc((ntax-1)*sizeof(int *));
		assert((int)temp_nodes);
		for (i=0;i<ntax-1;i++) {
			temp_nodes[i]=(int *)malloc(2*sizeof(int));
			assert((int) temp_nodes[i]);
		}
	}
	anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
	assert((int) anc);
	a_up=(alignment ***)malloc(((2*ntax)-1)*sizeof(alignment **));
	assert((int) a_up);
	for (i=0;i<2*ntax-1;i++) {
		a_up[i]=(alignment **)malloc(2*sizeof(alignment *));
		assert((int) a_up[i]);
		a_up[i][0]=NULL;
		a_up[i][1]=NULL;
	}
	if (values->new_optimization) { 
		a_down=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
		assert((int) a_down);
		for (i=0;i<2*ntax-1;i++) a_down[i]=NULL;
		corres=(int ***)malloc((2*(ntax))*sizeof(int **));
		assert((int) corres);
		for (i=0;i<2*ntax;i++) corres[i]=NULL;
		}
	cur_align=NULL;
	cur_align2=NULL;
	temp_align=NULL;
	temp_align2=NULL;
	if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Adding sequence ");
	/*set so no automatic cladogram*/  
	score_holder=values->phylo_score;
	values->phylo_score=0;
	values->in_some=1;

	/*get initial tree REMEMBER all < ntax are shifted one*/
	nodes[0][0]=0;/*this should never be referenced*/
	nodes[0][1]=ntax+1; 	/*real base ofmultiple alignment or whatever*/
     nodes[1][0]=1;
	nodes[1][1]=2;

	anc[0]=ntax;
	anc[1]=ntax+1;
	anc[2]=ntax+1;
	anc[ntax+1]=ntax;
	anc[ntax]=ntax;
	will_be_OK=1;

	/* and optimized nodes*/
	if (a[ntax+1]) a[ntax+1]=dump_align(a[ntax+1]);
	if (!values->new_optimization) {
		a[ntax+1]=nw(a[0],a[1],values); /*these are hand redone bit really 1-1 AND 2-1*/
		initial_length=a[ntax+1]->score;
		}
	else {
		values->in_optimize=1;
		a[ntax+1]=nw(a[0],a[1],values);  
		initial_length=a[ntax+1]->score;
		/*print_paup(a[ntax+1],values);*/
		if (temp_align) temp_align=dump_align(temp_align);
		temp_align=make_align(a[ntax+1]);
		temp_align=new_make_ambig(temp_align,values);/*anc with gaps*/
		/*get correspondance arrays*/
		if (strcmp(a[0]->name,a[1]->name)>0) {
			d2_thang=0;
			d1_thang=1;
			}
		else {
			d2_thang=1;
			d1_thang=0;
			}
		if (corres[1]) { /*these numbers are because of the taxon # jack up*/
			free(corres[1][0]);
			free(corres[1][1]);
			free(corres[1]);
			}
		if (corres[2]) {
			free(corres[2][0]);
			free(corres[2][1]);
			free(corres[2]);
			}
		corres[1]=get_correspondances(temp_align->s[0],a[ntax+1]->s[d1_thang],temp_align->length);
		corres[2]=get_correspondances(temp_align->s[0],a[ntax+1]->s[d2_thang],temp_align->length);
		if (temp_align->type_weight) {
			modify_corres(corres[1],temp_align);
			modify_corres(corres[2],temp_align);
			}
		a[ntax+1]=make_ambig(a[ntax+1],values);/*anc wo gaps*/
		values->in_optimize=0;
		free(temp_align);
		}
	/*get up passes (a's are the down passes)
	the uppass notation is for the up-pass which would be needed if the second array were hit
	 that is up pass [x][0] is the uppass needed if the addition comes between node[x] and its left [0] descendent
	 it is in essence the entire alignment minus node[x][0]	*/
	/*there are no up-passes for the root node (since never used or the terminals (again never used) */
	if (!values->new_optimization) {
		a_up[ntax+1][0]=make_align(a[1]);
		a_up[ntax+1][1]=make_align(a[0]);
	}
	else {
		for (i=0;i<ntax-1;i++) a_down[i]=make_align(a[i]);
		a_down[ntax+1]=make_align(a[ntax+1]);
		}

	for (i=0;i<ntax-3;i++) {/*add taxa in turn but numbers are jacked up by one*/
		if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d ",i+3);
		if (i==(ntax-4)) if (!(PARALLEL*values->jackboot)) fprintf(stderr,"\n");
		value=HUGE_COST;
		if (PARALLEL*parallel_modify) {
			num_to_receive=0;
			for (j=1;((num_to_receive<(values->num_hosts-1)) && (j<=((2*i)+3)));j++) {
				if (values->groups) {
					for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=nodes[ii][jj];
					temp_nodes[i+2][0]=i+3;
					temp_nodes[i+2][1]=nodes[j/2][j%2];
					temp_nodes[j/2][j%2]=ntax+i+2;
					cur_groups=get_cur_groups(temp_nodes,ntax,i+4,&ngroups2);
					will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,i+4-1);
					for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
					free(cur_groups);
				}
				if (will_be_OK) {
					/*send out first num_host jobs*/
					d1=j/2;
					d2=j%2;
					d3=nodes[d1][d2];
					if (d3 < ntax) d3-=1;
					pvm_initsend( PvmDataDefault );
					/*pack info*/
					pvm_pkint(&score_holder,1,1);
					pvm_pkint(&j,1,1);
					pvm_pkint(&d1,1,1);
					pvm_pkint(&d2,1,1);
					pvm_pkint(&d3,1,1);
					pvm_pkint(&i,1,1);
					pack_align(a[d3]);
					if (d1>0) {
						if (!values->new_optimization) pack_align(a_up[d1+ntax][d2]);
						else pack_align(a[d1+ntax]);
						}
					pvm_send(values->tids[++num_to_receive],PVM_NEW_BUILD_ADD_FUNCTION);
					
				}
			}
			while (num_to_receive>0) { /*rest as get 'em send*/
				bufid=pvm_recv(-1,PVM_NEW_BUILD_ALIGNMENTS_ONLY_SEND);
				if (bufid) {
					info=pvm_bufinfo(bufid, &bytes, &type, &source);
					/*if more to do send out next*/
					if (j<=((2*i)+3)) {
						if (values->groups) {
							for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=nodes[ii][jj];
							temp_nodes[i+2][0]=i+3;
							temp_nodes[i+2][1]=nodes[j/2][j%2];
							temp_nodes[j/2][j%2]=ntax+i+2;
							cur_groups=get_cur_groups(temp_nodes,ntax,i+4,&ngroups2);
							will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,i+4-1);
							for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
							free(cur_groups);
						}
						if (will_be_OK) {
							d1=j/2;
							d2=j%2;
							d3=nodes[d1][d2];
							if (d3 < ntax) d3-=1;
							pvm_initsend( PvmDataDefault );
							/*pack info*/
							pvm_pkint(&score_holder,1,1);
							pvm_pkint(&j,1,1);
							pvm_pkint(&d1,1,1);
							pvm_pkint(&d2,1,1);
							pvm_pkint(&d3,1,1);
							pvm_pkint(&i,1,1);
							pack_align(a[d3]);
							if (d1>0)  {
								if (!values->new_optimization) pack_align(a_up[d1+ntax][d2]);
								else pack_align(a[d1+ntax]);
								}
							pvm_send(source,PVM_NEW_BUILD_ADD_FUNCTION);
							++num_to_receive;
						}
						++j;
					}
					/*unpack this result*/
					pvm_upkint(&index,1,1);
					pvm_upkint(&cur_val,1,1);
					
					if (cur_val <value) {
						/*if new value set and keep alignment (2 if last i) */
						if (!values->new_optimization) {
								if (cur_align) cur_align=dump_align(cur_align);
								cur_align=unpack_align_and_score(values);
													}
						best_j=index;
						value=cur_val;
					}
					else if (!values->new_optimization) pvm_freebuf(bufid);
					--num_to_receive;
					
				}
			}
		}
		else {
			for (j=1;j<=((2*i)+3);j++) {/*places to go*/
				if (values->groups) {
					for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=nodes[ii][jj];
					temp_nodes[i+2][0]=i+3;
					temp_nodes[i+2][1]=nodes[j/2][j%2];
					temp_nodes[j/2][j%2]=ntax+i+2;
					cur_groups=get_cur_groups(temp_nodes,ntax,i+4,&ngroups2);
					will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,i+4-1);
					for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
					free(cur_groups);
				}
				if (will_be_OK) {
					d1=j/2;
					d2=j%2;
					d3=nodes[d1][d2];
					if (d3 < ntax) d3-=1;
					if (!values->new_optimization) {
						temp_align=nw(a[i+2],a[d3],values);/*would be i+3 but the adjustment for jack-up*/
						if (d1>0) {
							temp_align2=nw(temp_align,a_up[d1+ntax][d2],values);
							if (score_holder>0) {
								if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
								else {
									values->phylo_score=score_holder;
									cur_val=temp_align2->score=cladogram(temp_align2,values);
									values->phylo_score=0;
									if (values->previous) values->previous=dump_align(values->previous);
									values->previous=make_align(temp_align2);
								}
							}
							else cur_val=temp_align2->score;
						}
						else {
							if (score_holder>0) {
								if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
								else {
									values->phylo_score=score_holder;
									cur_val=temp_align->score=cladogram(temp_align,values);
									values->phylo_score=0;
									if (values->previous) values->previous=dump_align(values->previous);
									values->previous=make_align(temp_align);
								}
							}
							else cur_val=temp_align->score;
						}
					}
					else {
						/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"(");*/
						if (d1>0) { /*like interval in phylo stuff*/
							values->in_optimize=1;
							if (temp_align2) temp_align2=dump_align(temp_align2);
							temp_align2=nw(a[d1+ntax],a[d3],values);/*would be i+3 but the adjustment for jack-up*/
							values->in_optimize=0;
							temp_align2=newer_make_ambig(temp_align2,values);
							temp_align=nw(a[i+2],temp_align2,values);
							}
						else temp_align=nw(a[i+2],a[d3],values); /*if below existant root*/
						cur_val=temp_align->score;
						}
					/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d",cur_val);*/
					if (cur_val<value) {
						value=cur_val;
						best_j=j;
						if (!values->new_optimization) {
							if (cur_align) cur_align=dump_align(cur_align);
							cur_align=make_align(temp_align);
							if (cur_align2) cur_align2=dump_align(cur_align2);
							if (temp_align2) cur_align2=make_align(temp_align2);
							}
						}
					temp_align=dump_align(temp_align);
					if (temp_align2) temp_align2=dump_align(temp_align2);
					/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,")");*/
					if (value==0) break;
				}
			}
		}
		
		/*make new nodes*/
		nodes[i+2][0]=i+3; /*i+2 above because of the 1 offset when < ntax*/
		d1=best_j/2;
		d2=best_j%2;
		nodes[i+2][1]=nodes[d1][d2];
		nodes[d1][d2]=ntax+i+2;
		/*get new anc and its anc*/
		anc[nodes[i+2][0]]=i+2+ntax;
		anc[nodes[i+2][1]]=i+2+ntax;
		anc[i+2+ntax]=d1+ntax;
		initial_length+=value;
		if (!values->new_optimization) {
			/*copy the best alignment*/
			if (a[i+2+ntax]) a[i+2+ntax]=dump_align(a[i+2+ntax]);
			a[i+2+ntax]=make_align(cur_align);
			}
		if (i<ntax-4) { /*still building*/
			if (!values->new_optimization) {
				if (!(PARALLEL*parallel_modify)) check=all_diagnose_tree_local(a,nodes,ntax,i+4,values,i+2+ntax);
				else check=all_diagnose_tree_local_parallel(a,nodes,ntax,i+4,values,i+2+ntax);
				/*get up_pass alignments*/
			get_up_pass2(a,nodes,ntax,values,a_up,i+4,anc);
				}
			else {
				if (!(PARALLEL*parallel_modify)) {
					check=all_diagnose_tree_local_new_opt(a_down,nodes,ntax,i+4,values,nodes[i+2][1],corres);
					get_up_pass_new_opt(a,a_down,a_up,nodes,ntax,values,i+4,anc);
					}
				else {
					check=all_diagnose_tree_local_parallel_new_opt(a_down,nodes,ntax,i+4,values,nodes[i+2][1],corres);
					/*get up_pass alignments*/
					get_up_pass_new_opt_parallel(a,a_down,a_up,nodes,ntax,values,i+4,anc);
					}
				}
			}
		else { /*added last*/ 
			if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
			if (!values->new_optimization) {
				if (PARALLEL*parallel_modify) best_aligns[0]=make_align(cur_align);
				else {
					if (cur_align2) best_aligns[0]=make_align(cur_align2);
					else best_aligns[0]=make_align(cur_align);
					}
				}
			else { 	/*this is unnessesary really only need the names the rest is bullshit*/
				if (!(PARALLEL*parallel_modify)) check=all_diagnose_tree_local_new_opt(a_down,nodes,ntax,ntax,values,nodes[i+2][1],corres); 
				else check=all_diagnose_tree_local_parallel_new_opt(a_down,nodes,ntax,ntax,values,nodes[i+2][1],corres);
				best_aligns[0]=make_align(a_down[nodes[0][1]]);
				initial_length=best_aligns[0]->score;
				}
			*l=1;
			}
		}
	values->number+=(2*ntax-5);
	/*if (values->new_optimization) best_aligns[0]->score=initial_length;			*/
	/*do_swap before freeeing*/
	if (values->new_optimization && (values->collapse!=2)) {
		/*redo name for collapse*/
		if (!values->asbr || (values->asbr && values->aquick)) check=get_collapsed_thang(best_aligns[0],nodes,ntax,a_down,PARALLEL*parallel_modify*1,values);
		}
	if (!(PARALLEL && values->jackboot)) if (values->asbr) fprintf(stderr,"Heuristic build yielded a single alignment at cost %d\n",best_aligns[0]->score);
	if (values->asbr) {
		best_aligns=do_swap_thang_align(a,ntax,nodes,l,best_aligns,values,a_up,anc,0,score_holder,parallel_modify);
		if (!values->aquick) best_aligns=do_swap_thang_align(a,ntax,nodes,l,best_aligns,values,a_up,anc,1,score_holder,parallel_modify);
	}
	else if (!values->new_optimization) { /*filter for alignments*/
		if (!(PARALLEL*parallel_modify)) check=all_diagnose_tree_local2(a,nodes,ntax,ntax,values);
		else check=all_diagnose_tree_local_parallel2(a,nodes,ntax,ntax,values);
		if (compare_aligns(best_aligns[0],a[nodes[0][1]])) {
			if (score_holder>0) {
				values->phylo_score=score_holder;
				a[nodes[0][1]]->score=cladogram(a[nodes[0][1]],values);
				values->phylo_score=0;
			}
			if (a[nodes[0][1]]->score < best_aligns[0]->score) {
				if (!(PARALLEL && values->jackboot)) fprintf(stderr,"	Filtering Alignment 0 => New Cost %d\n",a[nodes[0][1]]->score);
				if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
				best_aligns[0]=make_align(a[nodes[0][1]]);
			}
		}
		/*at any rate fix name so matches order for OPTALI too*/
		free(best_aligns[0]->name);
		best_aligns[0]->name=(char *)malloc((1+strlen(a[nodes[0][1]]->name))*sizeof(char));
		assert((int) best_aligns[0]->name);
		strcpy(best_aligns[0]->name,a[nodes[0][1]]->name);
		values->number+=((ntax-1)/2);
	}
	
	if (anc) free(anc);
	if (a_up) {
		for (i=0;i<2*ntax-1;i++) {
			if (a_up[i][0]) a_up[i][0]=dump_align(a_up[i][0]);
			if (a_up[i][1]) a_up[i][1]=dump_align(a_up[i][1]);
			free(a_up[i]);
		}
		free(a_up);
	}
	if (temp_align) temp_align=dump_align(temp_align);
	if (temp_align2) temp_align2=dump_align(temp_align2);
	if (cur_align) cur_align=dump_align(cur_align);
	if (cur_align2) cur_align2=dump_align(cur_align2);
	if (values->groups) {
		for (i=0;i<ntax-1;i++) free(temp_nodes[i]);
		free(temp_nodes);
	}
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
		}
	
	values->phylo_score=score_holder;
	values->in_some=0; /*return initial length ?thang for new opt filter for names???*/
	return best_aligns;
}

void get_up_pass(a,nodes,ntax,values,a_up,nent,anc)
alignment **a,***a_up;
int **nodes,ntax,nent,*anc;
parameters *values;
{
	int i,j,found_one,d1,d2;
	int branch;

	for (i=1;i<nent-1;i++) {
		if (a_up[i+ntax][0]) a_up[i+ntax][0]=dump_align(a_up[i+ntax][0]);
		if (a_up[i+ntax][1]) a_up[i+ntax][1]=dump_align(a_up[i+ntax][1]);
	}
	/*get the first two*/
	d1=nodes[nodes[0][1]-ntax][1];
	d2=nodes[nodes[0][1]-ntax][0];
	/*if terminal*/
	if (d1<ntax) d1-=1;
	if (d2<ntax) d2-=1;
	a_up[nodes[0][1]][0]=make_align(a[d1]);
	a_up[nodes[0][1]][1]=make_align(a[d2]);

	found_one=1;
	while (found_one) {
		found_one=0;
		for (i=1;i<nent-1;i++) if ((!a_up[i+ntax][0]) && (a_up[anc[i+ntax]][0])){
			found_one=1;
			d1=nodes[i][0];
			d2=nodes[i][1];
			if (nodes[anc[i+ntax]-ntax][0]==(i+ntax)) branch=0;
			else branch=1;
			if (d1<ntax) d1-=1;
			if (d2<ntax) d2-=1;
			a_up[i+ntax][0]=nw(a_up[anc[i+ntax]][branch],a[d2],values);
			a_up[i+ntax][1]=nw(a_up[anc[i+ntax]][branch],a[d1],values);
		}
	}
}


void get_up_pass2(a,nodes,ntax,values,a_up,nent,anc)
alignment **a,***a_up;
int **nodes,ntax,nent,*anc;
parameters *values;
{
	int i,j,k,l;
	int n_seqs,length;
	char **names, **bases;
	int found,seq_counter,*all_gaps;
	int new_length,l_counter,d1;
	alignment *an,*auij;

		/*deallocate 'em all*/
	for (i=1;i<nent-1;i++) {
		if (a_up[i+ntax][0]) a_up[i+ntax][0]=dump_align(a_up[i+ntax][0]);
		if (a_up[i+ntax][1]) a_up[i+ntax][1]=dump_align(a_up[i+ntax][1]);
	}

	/*BASIC IDEA => a_up[i][j]=a[nodes[0][1]]-a[nodes[i-ntax][j]]*/
	an=a[nodes[0][1]];
	/*alloacte holders*/
	n_seqs=an->n_seqs;
	length=an->length;
	names=(char **)malloc(n_seqs*sizeof(char *));
	assert((int)names);
	bases=(char **)malloc((n_seqs+an->type_weight)*sizeof(char *));
	assert((int)bases);
	for (i=0;i<(n_seqs+an->type_weight);i++) {
		if (i<n_seqs) {
			names[i]=(char *)malloc(100*sizeof(char));
			assert((int)names[i]);
			}
		bases[i]=(char *)malloc((1+length)*sizeof(char));
		assert((int)bases[i]);
	}
	all_gaps=(int *)malloc(length*sizeof(int));
	assert((int)all_gaps);

	/*doit*/
	for (i=ntax+1;i<ntax+nent-1;i++) for (j=0;j<2;j++){
		d1=nodes[i-ntax][j];
		if (d1 < ntax) d1-=1;
		seq_counter=0;
		for (k=0;k<n_seqs;k++) {
			found=0;
			for (l=0;l<a[d1]->n_seqs;l++) if (!strcmp(an->taxon_name[k],a[d1]->taxon_name[l])) found=1;
			if (!found) { /*want these*/
				strcpy(names[seq_counter],an->taxon_name[k]);
				strcpy(bases[seq_counter],an->s[k]);
				++seq_counter;
			}
		}/*k*/
		if (an->type_weight) {
			strcpy(bases[seq_counter],an->s[an->n_seqs]);/*type_weight stuff*/
			}
		/*check seqs copied are correct in number*/
		if (seq_counter==(n_seqs-a[d1]->n_seqs)) { /*else in sbr/tbr non-needed up-passes*/
			/*scan for all gaps in remainder*/
			for (k=0;k<length;k++) {
				all_gaps[k]=1;
				for (l=0;l<seq_counter+an->type_weight;l++) {
					if (bases[l][k]!='-') all_gaps[k]=0;
					if ((l==seq_counter) && an->type_weight && (an->s[l][k]!=0)) all_gaps[k]=0;
					}
			}
			new_length=0;
			for (k=0;k<length;k++) if (!all_gaps[k]) ++new_length;

			/*allocate and copy to a_up[i][j]*/
			a_up[i][j]=NULL;
			a_up[i][j]=(alignment *)malloc(sizeof(alignment));
			assert((int)a_up[i][j]);
			auij=a_up[i][j];
			auij->n_seqs=seq_counter;
			auij->score=0;
			auij->length=new_length;
			auij->type_weight=an->type_weight;
			auij->name=(char *)malloc(8*sizeof(char));
			assert((int) auij->name);
			auij->name[0]='U';
			auij->name[1]='p';
			auij->name[2]='-';
			auij->name[3]='P';
			auij->name[4]='a';
			auij->name[5]='s';
			auij->name[6]='s';
			auij->name[7]='\0';
			auij->s=(char **)malloc((auij->n_seqs+auij->type_weight)*sizeof(char *));
			assert((int)auij->s);
			for (k=0;k<(auij->n_seqs+auij->type_weight);k++) {
				auij->s[k]=(char *)malloc((auij->length+1)*sizeof(char));
				assert((int)auij->s[k]);
				l_counter=0;
				for (l=0;l<length;l++) {
					if (!all_gaps[l]) auij->s[k][l_counter++]=bases[k][l];
					else if (auij->type_weight && (k==auij->n_seqs)) {
							if (bases[k][l]>60) auij->s[k][l_counter-1]=bases[k][l];
							}
					}
				auij->s[k][new_length]='\0';
			}

			auij->taxon_name=(char **)malloc((auij->n_seqs)*sizeof(char *));
			assert((int)auij->taxon_name);
			for (k=0;k<auij->n_seqs;k++) {
				auij->taxon_name[k]=(char *)malloc((strlen(names[k])+1)*sizeof(char));
				assert((int)auij->taxon_name[k]);
				auij->taxon_name[k]=(char *)strcpy(auij->taxon_name[k],names[k]);
			}
			/*print_alignment(auij);
		printf("\n");*/
		}
	}

	/*freeing*/
	for (i=0;i<n_seqs+an->type_weight;i++) {
		if (i<n_seqs) free(names[i]);
		free(bases[i]);
	}
	free(names);
	free(bases);
	free(all_gaps);
}

void get_up_pass3(a,nodes,ntax,values,a_up,nent,anc,desc)
alignment **a,***a_up;
int **nodes,ntax,nent,*anc,*desc;
parameters *values;
{
	int i,j,k,l;
	int n_seqs,length;
	char **names, **bases;
	int found,seq_counter,*all_gaps;
	int new_length,l_counter,d1;

	/*deallocate 'em all*/
	for (i=1;i<nent-1;i++) {
		if (a_up[i+ntax][0]) a_up[i+ntax][0]=dump_align(a_up[i+ntax][0]);
		if (a_up[i+ntax][1]) a_up[i+ntax][1]=dump_align(a_up[i+ntax][1]);
	}

	/*BASIC IDEA => a_up[i][j]=a[nodes[0][1]]-a[nodes[i-ntax][j]]*/

	/*alloacte holders*/
	n_seqs=a[nodes[0][1]]->n_seqs;
	length=a[nodes[0][1]]->length;
	names=(char **)malloc(n_seqs*sizeof(char *));
	assert((int)names);
	bases=(char **)malloc(n_seqs*sizeof(char *));
	assert((int)bases);
	for (i=0;i<n_seqs;i++) {
		names[i]=(char *)malloc(100*sizeof(char));
		assert((int)names[i]);
		bases[i]=(char *)malloc((1+length)*sizeof(char));
		assert((int)bases[i]);
	}
	all_gaps=(int *)malloc(length*sizeof(int));
	assert((int)all_gaps);

	/*doit*/
	for (i=ntax+1;i<ntax+nent-1;i++) for (j=0;j<2;j++){
		d1=nodes[i-ntax][j];
		if (d1 < ntax) d1-=1;
		seq_counter=0;
		for (k=0;k<n_seqs;k++) {
			found=0;
			for (l=0;l<a[d1]->n_seqs;l++) if (!strcmp(a[nodes[0][1]]->taxon_name[k],a[d1]->taxon_name[l])) found=1;
			if (!found) { /*want these*/
				strcpy(names[seq_counter],a[nodes[0][1]]->taxon_name[k]);
				strcpy(bases[seq_counter],a[nodes[0][1]]->s[k]);
				++seq_counter;
			}
		}/*k*/
		/*check seqs copied are correct in number*/
		if (seq_counter==(n_seqs-a[d1]->n_seqs)) { /*else in sbr/tbr non-needed up-passes*/
			/*scan for all gaps in remainder*/
			for (k=0;k<length;k++) {
				all_gaps[k]=1;
				for (l=0;l<seq_counter;l++) if (bases[l][k]!='-') all_gaps[k]=0;
			}
			new_length=0;
			for (k=0;k<length;k++) if (!all_gaps[k]) ++new_length;

			/*allocate and copy to a_up[i][j]*/
			a_up[i][j]=(alignment *)malloc(sizeof(alignment));
			assert((int)a_up[i][j]);
			a_up[i][j]->n_seqs=seq_counter;
			a_up[i][j]->score=0;
			a_up[i][j]->length=new_length;
			a_up[i][j]->type_weight=0;
			a_up[i][j]->name=(char *)malloc(8*sizeof(char));
			assert((int)a_up[i][j]->name);
			a_up[i][j]->name[0]='U';
			a_up[i][j]->name[1]='p';
			a_up[i][j]->name[2]='-';
			a_up[i][j]->name[3]='P';
			a_up[i][j]->name[4]='a';
			a_up[i][j]->name[5]='s';
			a_up[i][j]->name[6]='s';
			a_up[i][j]->name[7]='\0';
			a_up[i][j]->s=(char **)malloc((a_up[i][j]->n_seqs+a_up[i][j]->type_weight)*sizeof(char *));
			assert((int)a_up[i][j]->s);
			for (k=0;k<(a_up[i][j]->n_seqs+a_up[i][j]->type_weight);k++) {
				a_up[i][j]->s[k]=(char *)malloc((a_up[i][j]->length+1)*sizeof(char));
				assert((int)a_up[i][j]->s[k]);
				l_counter=0;
				for (l=0;l<length;l++) if (!all_gaps[l]) a_up[i][j]->s[k][l_counter++]=bases[k][l];
				a_up[i][j]->s[k][new_length]='\0';
			}

			a_up[i][j]->taxon_name=(char **)malloc((a_up[i][j]->n_seqs)*sizeof(char *));
			assert((int)a_up[i][j]->taxon_name);
			for (k=0;k<a_up[i][j]->n_seqs;k++) {
				a_up[i][j]->taxon_name[k]=(char *)malloc((strlen(names[k])+1)*sizeof(char));
				assert((int)a_up[i][j]->taxon_name[k]);
				a_up[i][j]->taxon_name[k]=(char *)strcpy(a_up[i][j]->taxon_name[k],names[k]);
			}
			/*print_alignment(a_up[i][j]);
		printf("\n");*/
		}
	}

	/*freeing*/
	for (i=0;i<n_seqs;i++) {
		free(names[i]);
		free(bases[i]);
	}
	free(names);
	free(bases);
	free(all_gaps);
	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"Done up_pass\n");exit(-1);*/
}

int all_diagnose_tree_local(a,nodes,ntaxa,nent,values,node_made)
int **nodes,node_made;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,found_one;


	num_seqs=(int *)malloc((ntaxa+ntaxa-1)*sizeof(int));
	assert((int)num_seqs);

	values->in_some=1;
	values->all_made=(char *)malloc(((2*ntaxa)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (values->number_of_input_alignments<1) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization)) nent_minus_one_new=nent-1;
	else for (n=0;n<=(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs);

	values->all_made[node_made]=1;

	found_one=1;
	while (found_one) {
		found_one=0;
		for (n=0;n<nent_minus_one_new;n++) if (values->all_made[n+ntaxa] && ((!values->all_made[nodes[n][0]]) || (!values->all_made[nodes[n][1]]))) {
			found_one=1;
			values->all_made[nodes[n][0]]=1;
			values->all_made[nodes[n][1]]=1;
		}
	}


	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u-",_memavl());*/
	found_one=1;
	while (found_one) {
		found_one=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				found_one=1;
				values->all_made[ntaxa+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				/*avoid unnesesary trees*/
				if ((a[d1]->n_seqs+a[d2]->n_seqs)<(nent_minus_one_new)) {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
					if (values->in_bandb_loop && (a[ntaxa+n]->score>values->current_bound)) goto second_end;
				}
				else {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
				}
			}
		}
	}

	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	values->in_some=0;
	return(a[nodes[0][1]]->score);

	/*if early exceed bound*/
second_end:
	free(num_seqs);
	free(values->all_made);
	values->in_some=0;
	values->all_made=NULL;
	return HUGE_COST;
}

int all_diagnose_tree_local_parallel2(a,nodes,ntaxa,nent,values)
int **nodes;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,found_one;
	int sent,received,score_holder,node_to_be_sent;
	int node_received,bufid,info,source,tag,bytes,type;
	int ii;

	num_seqs=(int *)malloc((ntaxa+ntaxa-1)*sizeof(int));
	assert((int)num_seqs);
	values->in_some=1;
	values->all_made=(char *)malloc(((2*ntaxa)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (values->number_of_input_alignments<1) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization)) nent_minus_one_new=nent-1;
	else for (n=0;n<=(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	/*mark stuff that's already done*/
	pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs);
	score_holder=0;
	found_one=1;
	while (found_one) {
		found_one=0;
		sent=received=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				found_one=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				node_to_be_sent=ntaxa+n;
				if (sent < (values->num_hosts-1)) {
					pvm_initsend( PvmDataDefault );
					pvm_pkint(&score_holder,1,1);
					pvm_pkint(&node_to_be_sent,1,1);
					pack_align(a[d1]);
					pack_align(a[d2]);
					pvm_send(values->tids[++sent],PVM_PARALLEL_NW);
				}
				else {
					bufid=pvm_recv(-1,PVM_PARALLEL_NW_DONE);
					if (bufid) {
						info=pvm_bufinfo(bufid, &bytes, &type, &source);
						pvm_upkint(&node_received,1,1);
						if (a[node_received]) a[node_received]=dump_align(a[node_received]);
						a[node_received]=unpack_align_and_score(values);
						values->all_made[node_received]=1;
						++received;
					}
					pvm_initsend( PvmDataDefault );
					pvm_pkint(&score_holder,1,1);
					pvm_pkint(&node_to_be_sent,1,1);
					pack_align(a[d1]);
					pack_align(a[d2]);
					pvm_send(source,PVM_PARALLEL_NW);
					++sent;
				}
			}
		}
		/*receive the remainder*/
		/*while (received<sent) {*/
		for (ii=received;ii<sent;ii++) {
			bufid=pvm_recv(-1,PVM_PARALLEL_NW_DONE);
			if (bufid) {
				info=pvm_bufinfo(bufid, &bytes, &type, &source);
				pvm_upkint(&node_received,1,1);
				if (a[node_received]) a[node_received]=dump_align(a[node_received]);
				a[node_received]=unpack_align_and_score(values);
				values->all_made[node_received]=1;
				/*++received;*/
			}
		}

	}
	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	values->in_some=0;
	return(a[nodes[0][1]]->score);

	/*if early exceed bound*/
second_end:
	free(num_seqs);
	free(values->all_made);
	values->in_some=0;
	values->all_made=NULL;
	return HUGE_COST;
}

int all_diagnose_tree_local_parallel(a,nodes,ntaxa,nent,values,node_made)
int **nodes,node_made;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,found_one;
	int sent,received,score_holder,node_to_be_sent;
	int node_received,bufid,info,source,tag,bytes,type;
	int ii;

	num_seqs=(int *)malloc((ntaxa+ntaxa-1)*sizeof(int));
	assert((int)num_seqs);
	values->in_some=1;
	values->all_made=(char *)malloc(((2*ntaxa)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (values->number_of_input_alignments<1) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization)) nent_minus_one_new=nent-1;
	else for (n=0;n<=(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	/*mark stuff that's already done*/
	pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs);
	values->all_made[node_made]=1;
	found_one=1;
	while (found_one) {
		found_one=0;
		for (n=0;n<nent_minus_one_new;n++) if (values->all_made[n+ntaxa] && ((!values->all_made[nodes[n][0]]) || (!values->all_made[nodes[n][1]]))) {
			found_one=1;
			values->all_made[nodes[n][0]]=1;
			values->all_made[nodes[n][1]]=1;
		}
	}

	score_holder=0;
	found_one=1;
	while (found_one) {
		found_one=0;
		sent=received=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				found_one=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				node_to_be_sent=ntaxa+n;
				if (sent < (values->num_hosts-1)) {
					pvm_initsend( PvmDataDefault );
					pvm_pkint(&score_holder,1,1);
					pvm_pkint(&node_to_be_sent,1,1);
					pack_align(a[d1]);
					pack_align(a[d2]);
					pvm_send(values->tids[++sent],PVM_PARALLEL_NW);
				}
				else {
					bufid=pvm_recv(-1,PVM_PARALLEL_NW_DONE);
					if (bufid) {
						info=pvm_bufinfo(bufid, &bytes, &type, &source);
						pvm_upkint(&node_received,1,1);
						if (a[node_received]) a[node_received]=dump_align(a[node_received]);
						a[node_received]=unpack_align_and_score(values);
						values->all_made[node_received]=1;
						++received;
					}
					pvm_initsend( PvmDataDefault );
					pvm_pkint(&score_holder,1,1);
					pvm_pkint(&node_to_be_sent,1,1);
					pack_align(a[d1]);
					pack_align(a[d2]);
					pvm_send(source,PVM_PARALLEL_NW);
					++sent;
				}
			}
		}
		/*receive the remainder*/
		/*while (received<sent) {*/
		for (ii=received;ii<sent;ii++) {
			bufid=pvm_recv(-1,PVM_PARALLEL_NW_DONE);
			if (bufid) {
				info=pvm_bufinfo(bufid, &bytes, &type, &source);
				pvm_upkint(&node_received,1,1);
				if (a[node_received]) a[node_received]=dump_align(a[node_received]);
				a[node_received]=unpack_align_and_score(values);
				values->all_made[node_received]=1;
				/*++received;*/
			}
		}

	}

	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	values->in_some=0;
	return(a[nodes[0][1]]->score);

	/*if early exceed bound*/
second_end:
	free(num_seqs);
	free(values->all_made);
	values->in_some=0;
	values->all_made=NULL;
	return HUGE_COST;
}

int all_diagnose_tree_local2(a,nodes,ntaxa,nent,values)
int **nodes;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,found_one;


	num_seqs=(int *)malloc((ntaxa+ntaxa-1)*sizeof(int));
	assert((int)num_seqs);

	values->in_some=1;
	values->all_made=(char *)malloc(((2*ntaxa)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (values->number_of_input_alignments<1) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization)) nent_minus_one_new=nent-1;
	else for (n=0;n<=(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs);


	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u-",_memavl());*/
	found_one=1;
	while (found_one) {
		found_one=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				found_one=1;
				values->all_made[ntaxa+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				/*avoid unnesesary trees*/
				if ((a[d1]->n_seqs+a[d2]->n_seqs)<(nent_minus_one_new)) {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
					if (values->in_bandb_loop && (a[ntaxa+n]->score>values->current_bound)) goto second_end;
				}
				else {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
				}
			}
		}
	}

	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	values->in_some=0;
	return(a[nodes[0][1]]->score);

	/*if early exceed bound*/
second_end:
	free(num_seqs);
	free(values->all_made);
	values->in_some=0;
	values->all_made=NULL;
	return HUGE_COST;
}

int all_diagnose_tree_local3(a,nodes,ntaxa,nent,values,desc)
int **nodes;
int ntaxa,nent,*desc;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,found_one;

	num_seqs=(int *)malloc((ntaxa+ntaxa-1)*sizeof(int));
	assert((int)num_seqs);

	values->in_some=1;
	values->all_made=(char *)malloc(((2*ntaxa)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (values->number_of_input_alignments<1) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization)) nent_minus_one_new=nent-1;
	else for (n=0;n<=(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs);

	for (n=ntaxa+1;n<ntaxa+nent-1;n++) if (!desc[n]) values->all_made[n]=1;

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"L ");
for (n=0;n<2*ntaxa-1;n++) if (!(PARALLEL*values->jackboot)) fprintf(stderr,"s[%d] %d %d, ",n,values->all_made[n],desc[n]);
if (!(PARALLEL*values->jackboot)) fprintf(stderr,"\n");		*/

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u-",_memavl());*/
	found_one=1;
	while (found_one) {
		found_one=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				found_one=1;
				values->all_made[ntaxa+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				/*avoid unnesesary trees*/
				if ((a[d1]->n_seqs+a[d2]->n_seqs)<(nent_minus_one_new)) {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
					if (values->in_bandb_loop && (a[ntaxa+n]->score>values->current_bound)) goto second_end;
				}
				else {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
				}
			}
		}
	}

	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL*values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	values->in_some=0;
	return(a[nodes[0][1]]->score);

	/*if early exceed bound*/
second_end:
	free(num_seqs);
	free(values->all_made);
	values->in_some=0;
	values->all_made=NULL;
	return HUGE_COST;
}





