/*Copyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Gladstein*/
/*In TBR should reroot at top ! in seperate SPR loop*/

#include "align3.h"

extern int ***best_nodes;

alignment **make_dave_align(a,ntax,nodes,l,best_aligns,values,parallel_modify)
alignment **a, **best_aligns;
int ntax, *l, **nodes,parallel_modify;
parameters *values;
{
	int i,j,k,ll,value,cur_val,n,m;
	int initial_length,d1,d2;
	int *anc,check;
	int d3,best_j;
	int bufid,type,bytes,source,index,info,num_to_do,num_to_receive;
	int **temp_nodes,will_be_OK,**cur_groups,ngroups2;
	int ii,jj,same,kk,done;
	int d1_thang, d2_thang;
	alignment **a_best;
	int *order,set,old_cur_val,temp,temp2;
	int *dumb_sup;
	alignment *temp_align,**check_buf;
	int best_k,local_grain,grain_really_sent;
	int found_better,smeg;

	temp_align=NULL;
	check=HUGE_COST;
	if (values->best_order & !values->reorder) {
		for (i=1;i<ntax-1;i++) {
			temp_align=nw(a[0],a[i],values);
			if (temp_align->score < check) {
				check=temp_align->score;
				best_j=i;
			}
			temp_align=dump_align(temp_align);
		}
		if (best_j != 1) {
			temp_align=make_align(a[1]);
			a[1]=dump_align(a[1]);
			a[1]=make_align(a[best_j]);
			a[best_j]=dump_align(a[best_j]);
			a[best_j]=make_align(temp_align);
			temp_align=dump_align(temp_align);
			}
		}

best_nodes=(int ***)malloc(values->keep_aligns*sizeof(int **));
assert((int) best_nodes);
for (j=0;j<values->keep_aligns;j++) {
    best_nodes[j]=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int)best_nodes[j]);
	for (i=0;i<ntax-1;i++) {
		best_nodes[j][i]=(int *)malloc(2*sizeof(int));
		assert((int) best_nodes[j][i]);
	}
	}
	temp_nodes=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int)temp_nodes);
	for (i=0;i<ntax-1;i++) {
		temp_nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int) temp_nodes[i]);
	}
	anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
	assert((int) anc);
	if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Adding sequence ");
	/*set so no automatic cladogram*/
	values->phylo_score=0;
	values->in_some=1;

	/*get initial tree REMEMBER all < ntax are shifted one*/
	nodes[0][0]=0;/*this should never be referenced*/
	nodes[0][1]=ntax+1;     /*real base ofmultiple alignment or whatever*/
	nodes[1][0]=1;
	nodes[1][1]=2;
	if (a[ntax+1]) a[ntax+1]=dump_align(a[ntax+1]);
	a[ntax+1]=nw(a[0],a[1],values); /*first node*/
	anc[0]=ntax;
	anc[1]=ntax+1;
	anc[2]=ntax+1;
	anc[ntax+1]=ntax;
	anc[ntax]=ntax;
	will_be_OK=1;

	a_best=(alignment **)malloc((ntax-1)*sizeof(alignment *));
	assert((int) a_best);
	for (i=0;i<ntax-1;i++) a_best[i]=NULL;

	/* and optimized nodes*/
	for (i=0;i<ntax-3;i++) {/*add taxa in turn but numbers are jacked up by one*/
		best_k=i+2;
		if (!(PARALLEL && values->jackboot)) if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d ",i+3);
		if (i==(ntax-4)) if (!(PARALLEL && values->jackboot)) fprintf(stderr,"\n");
		value=HUGE_COST;
		if (i==0) old_cur_val=HUGE_COST;
		if (PARALLEL*parallel_modify) {
			/*send all i and nodes */
			pvm_initsend( PvmDataDefault );
			pvm_pkint(&i,1,1);
			if (values->best_order) for (ii=0;ii<ntax-1;ii++) pack_align(a[ii]);
			for (j=1;j<i+2;j++) pack_align(a[ntax+j]); /*i+2 because only previous i is done*/
			for (ii=0;ii<i+3;ii++) for (jj=0;jj<2;jj++) pvm_pkint(&(nodes[ii][jj]),1,1);
			pvm_mcast(&(values->tids[1]),min( values->num_hosts-1,((2*i)+3)), PVM_SEND_ANODES);
			num_to_receive=0;
			if (values->grain_size !=1) {
				local_grain=(((2*i)+3)/(values->num_hosts-1));
				if (!local_grain) local_grain=1;
				local_grain=min(local_grain,values->grain_size);
			}
			else local_grain=1;
			grain_really_sent=0;
			j=kk=1;
			/*j is the host kk is th placement ll is the internal counter NB if grain 1 k=j
			fprintf(stderr,"Local grain=%d\n",local_grain);*/
			for (j=1;((num_to_receive<(values->num_hosts-1)) && ((kk+local_grain-1)<=((2*i)+3)));j++) {
				pvm_initsend( PvmDataDefault );
				/*pack info basics*/
				pvm_pkint(&i,1,1);
				/*if (values->best_order) if (i>0) for (ii=i-1;ii<ntax-1;ii++) pack_align(a[ii]);*/
				pvm_pkint(&j,1,1);
				pvm_pkint(&value,1,1);
				pvm_pkint(&local_grain,1,1);
				pvm_pkint(&kk,1,1);
				/*for (smeg=0;smeg<local_grain;smeg++) {
					for (ii=0;ii<i+3;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=nodes[ii][jj];
					temp_nodes[i+2][0]=i+3;
					temp_nodes[i+2][1]=nodes[kk/2][kk%2];
					temp_nodes[kk/2][kk%2]=ntax+i+2;
					for (ii=0;ii<i+3;ii++) for (jj=0;jj<2;jj++) pvm_pkint(&(temp_nodes[ii][jj]),1,1);
					kk++;
					}*/
				pvm_send(values->tids[++num_to_receive],PVM_DAVE_BUILD);
				kk+=local_grain;
				/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"S%d ",num_to_receive);*/
			}
			/*fprintf(stderr,"End of send topos %d parts %d\n",kk-1,num_to_receive);*/
			/*local_grain=1; to get straglers*/
			j=kk;
			/*exit(-1);*/
			/*receive/send*/
			while (num_to_receive>0) { /*rest as get 'em send for now--only one each no matter what the grain*/
				bufid=pvm_recv(-1,PVM_DAVE_BUILD_DONE);
				if (bufid) {
					info=pvm_bufinfo(bufid, &bytes, &type, &source);
					pvm_upkint(&index,1,1);
					pvm_upkint(&cur_val,1,1);
					if (cur_val<value) {value=cur_val;found_better=1;}
					else found_better=0;
					/*if more to do send out next*/
					/*if ((kk+local_grain) > ((2*i)+3)) local_grain=((2*i)+3)-kk;*/
					if ((kk + local_grain -1) > ((2*i)+3)) local_grain=( ((2*i)+3)-kk+1);
					if (kk <=((2*i)+3)) {
						pvm_initsend( PvmDataDefault );
						/*pack info*/
						pvm_pkint(&i,1,1);
						/*if (values->best_order) if (i>0) for (ii=i-1;ii<ntax-1;ii++) pack_align(a[ii]);*/
						pvm_pkint(&j,1,1);
						pvm_pkint(&value,1,1);
						pvm_pkint(&local_grain,1,1);
						pvm_pkint(&kk,1,1);
						/*for (smeg=0;smeg<local_grain;smeg++) {
							for (ii=0;ii<i+3;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=nodes[ii][jj];
							temp_nodes[i+2][0]=i+3;
							temp_nodes[i+2][1]=nodes[kk/2][kk%2];
							temp_nodes[kk/2][kk%2]=ntax+i+2;
							for (ii=0;ii<i+3;ii++) for (jj=0;jj<2;jj++) pvm_pkint(&(temp_nodes[ii][jj]),1,1);
							++kk;
							}*/
						pvm_send(source,PVM_DAVE_BUILD);
						++num_to_receive;
						++j;
						kk+=local_grain;
					}/*j*/
					/*unpack results still only recieving one fomr each--only the best sent*/
					if (found_better) {
						/*value=cur_val;*/
						best_j=index;
						pvm_upkint(&best_k,1,1);
						for (k=0;k<i+2;k++) {
							if (a_best[k]) a_best[k]=dump_align(a_best[k]);
							a_best[k]=unpack_align_and_score(values);
							}
						if (value==old_cur_val) j=HUGE_COST;
					}
					else pvm_freebuf(bufid);
					--num_to_receive;
				}/*bufid*/
			}/*while*/
			/*send all kill signal to dump nodes*/
			pvm_initsend( PvmDataDefault );
			pvm_mcast(&(values->tids[1]),min( values->num_hosts-1,((2*i)+3)), PVM_DAVE_DUMP_NODES);
		}/*PARALLEL*parallel_modify*/
else {
	order=(int *)malloc(ntax*sizeof(int));
	assert((int)order);
	for (j=1;j<=((2*i)+3);j++) {/*places to go*/
		for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=nodes[ii][jj];
		temp_nodes[i+2][0]=i+3;
		temp_nodes[i+2][1]=nodes[j/2][j%2];
		temp_nodes[j/2][j%2]=ntax+i+2;
		anc[i+3]=i+2+ntax;
		anc[temp_nodes[i+2][1]]=i+2+ntax;
		anc[i+2+ntax]=(j/2)+ntax;
		if (values->groups) {
			cur_groups=get_cur_groups(temp_nodes,ntax,i+4,&ngroups2);
			will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,i+4-1);
			for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
			free(cur_groups);
		}
		if (will_be_OK) {
			/*set up nodes to do and order*/
			for (k=0;k<ntax;k++) order[k]=(-1);
			order[0]=i+2;
			num_to_do=1;
			done=0;
			while (!done) {
				done=1;
				for (k=1;k<i+3;k++) {
					same=0;
					for (kk=0;kk<num_to_do;kk++) if (order[kk]==k) same=1;
					if (!same) {
						for (kk=0;kk<num_to_do;kk++) {
							if (temp_nodes[k][0]==(order[kk]+ntax)) {
								order[num_to_do++]=k;
								done=0;
							}
							else if (temp_nodes[k][1]==(order[kk]+ntax)) {
								order[num_to_do++]=k;
								done=0;
							}
						}
					}
				}
			}
			/*do_alignment*/
			/*fprintf(stderr,"%d ",values->new_optimization);*/
			/*could put loop in here for best adddiotn sequence*/
			if (values->best_order) {
				for (k=i+2;k<ntax-1;k++) {
					if (k!=(i+2)) {
						temp_align=make_align(a[i+2]);
						a[i+2]=dump_align(a[i+2]);
						a[i+2]=make_align(a[k]);
					}
					cur_val=all_diagnose_tree_cut_off(a,temp_nodes,ntax,i+4,values,order,num_to_do,value,a_best);
					if (cur_val<value) {
						value=cur_val;
						best_j=j;
						best_k=k;
						if (value==old_cur_val) goto end_of_loop_zero;
					}
					if (k!=(i+2)) {
						a[i+2]=dump_align(a[i+2]);
						a[i+2]=make_align(temp_align);
						temp_align=dump_align(temp_align);
					}
					/*fprintf(stderr,".");*/
				}
			}
			else {
				cur_val=all_diagnose_tree_cut_off(a,temp_nodes,ntax,i+4,values,order,num_to_do,value,a_best);
				if (cur_val<value) {
					value=cur_val;
					best_j=j;
					best_k=i+2;
					if (value==old_cur_val) goto end_of_loop_zero;
				}
			}
		}
	}/*j*/
end_of_loop_zero:
	;
	free(order);
}/*!PARRALEL*/
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
/*value=all_diagnose_tree_here(a,nodes,ntax,i+4,values,HUGE_COST);*/
/*copy bestnodes over for next round*/
if (values->best_order) {
	if (best_k != (i+2))  {
		temp_align=make_align(a[i+2]);
		a[i+2]=dump_align(a[i+2]);
		a[i+2]=make_align(a[best_k]);
		a[best_k]=dump_align(a[best_k]);
		a[best_k]=make_align(temp_align);
		temp_align=dump_align(temp_align);
	}
}
for (k=0;k<i+2;k++) {
	if (a[k+ntax+1]) a[k+ntax+1]=dump_align(a[k+ntax+1]);
	a[k+ntax+1]=a_best[k];
	a_best[k]=NULL;
}
if (i==(ntax-4)) {
	/*added last*/
	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"adding last\n");*/
	if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
	best_aligns[0]=make_align(a[nodes[0][1]]);
	best_aligns[0]->score=value;
	*l=1;
}
else old_cur_val=value;
	}/*i loop*/
	/*do_swap before freeeing*/
	/*fprintf(stderr,"%d length",all_diagnose_tree_here3(a,nodes,ntax,ntax,values,HUGE_COST));*/
	check_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) check_buf);
	for (i=0;i<2*ntax-1;i++) check_buf[i]=NULL;
	for (i=0;i<ntax-1;i++) check_buf[i]=make_align(a[i]);


	if (values->new_optimization && (values->collapse==2)) {
		free(best_aligns[0]->name);
		dumb_sup=(int *)malloc(ntax*sizeof(int));
		assert((int) dumb_sup);
		for (i=0;i<ntax;i++) dumb_sup[i]=1;
		best_aligns[0]->name=(char *)get_names_with_collapse(a,nodes,ntax,ntax,values,dumb_sup);
		free(dumb_sup);
	}
	else if (values->new_optimization && (values->collapse!=2)) {
		/*redo name for collapse*/
		if (!values->asbr || (values->asbr && values->aquick)) check=get_collapsed_thang(best_aligns[0],nodes,ntax,check_buf,PARALLEL*parallel_modify*1,values);
	}
	if (values->poopsy && (!values->jackboot)) fprintf(stderr,"%s\n",best_aligns[0]->name);
	for (i=0;i<ntax-1;i++) for (j=0;j<2;j++) best_nodes[0][i][j]=nodes[i][j];
	if (values->asbr) {
		if (values->VERBOSE && (!(PARALLEL && values->jackboot))) fprintf(stderr,"Heuristic build yielded 1 alignment at length %d\n",best_aligns[0]->score);
		if (!values->tbr_first) {
			temp=0;
			if (values->atbr) {
				temp=1;
				values->atbr=0;
				temp2=values->aquick;
				values->aquick=1;
			}
			if (!values->aquick) best_aligns=do_swap_thang_align_dave(a,ntax,nodes,l,best_aligns,values,anc,1,parallel_modify);
			else best_aligns=do_swap_thang_align_dave(a,ntax,nodes,l,best_aligns,values,anc,0,parallel_modify);
			if (temp) {
				values->atbr=1;
				values->aquick=temp2;
				if (!values->aquick) best_aligns=do_swap_thang_align_dave(a,ntax,nodes,l,best_aligns,values,anc,1,parallel_modify);
				else best_aligns=do_swap_thang_align_dave(a,ntax,nodes,l,best_aligns,values,anc,0,parallel_modify);
			}
		}
		else {
			if (!values->aquick) best_aligns=do_swap_thang_align_dave(a,ntax,nodes,l,best_aligns,values,anc,1,parallel_modify);
			else best_aligns=do_swap_thang_align_dave(a,ntax,nodes,l,best_aligns,values,anc,0,parallel_modify);
		}
	}
	/*for (i=0;i<ntax-1;i++) {
    for (j=i+1;j<ntax-1;j++) {
        temp_align=nw(a[i],a[j],values);
        printf("(%s %s + %d)",a[i]->name,a[j]->name,temp_align->score);
        temp_align=dump_align(temp_align);
    }
    printf("\n");
}
fprintf(stderr,"%d length",all_diagnose_tree_here3(a,nodes,ntax,ntax,values,HUGE_COST));*/

	for (i=0;i<2*ntax-1;i++) if (check_buf[i]) check_buf[i]=dump_align(check_buf[i]);
	free(check_buf);

	if (anc) free(anc);
	for (i=0;i<ntax-1;i++) free(temp_nodes[i]);
	free(temp_nodes);
	for (i=0;i<ntax-1;i++) if (a_best[i]) a_best[i]=dump_align(a_best[i]);
	free(a_best);
	return best_aligns;
}

int all_diagnose_tree_cut_off(a,nodes,ntax,nent,values,order,num_to_do,bound,a_best)
int **nodes;
int ntax,nent;
alignment **a,**a_best;
parameters *values;
int *order, num_to_do,bound;
{
	int n,d1,d2,nent_minus_one_new;
	int new_opt_score,i,j;
	alignment *temp_align;
	alignment **backup;
	int *changed,same,num_seqs;
	int nodes_actual=0;
	char *temp_name;

	if (nent==ntax) ++values->number;

	backup=(alignment **)malloc(num_to_do*sizeof(alignment *));
	assert((int) backup);
	for (i=0;i<num_to_do;i++) backup[i]=NULL;
	changed=(int *)malloc(num_to_do*sizeof(int));
	assert((int)changed);
	for (i=0;i<num_to_do;i++) {
		if (a[ntax+order[i]]) backup[i]=make_align(a[ntax+order[i]]);
		changed[i]=0;
	}
	nent_minus_one_new=nent-1;
	new_opt_score=0;
	for (i=1;i<nent-1;i++) {
		same=0;
		for (j=0;j<num_to_do;j++) if (i==order[j]) same=1;
		if (!same) new_opt_score+=a[ntax+i]->score;
	}
	temp_align=NULL;

	for (i=0;i<num_to_do;i++){
		n=order[i];
		d1=nodes[n][0];
		d2=nodes[n][1];
		/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Making %d from %d and %d\n",order[i],d1,d2);*/
		/*fix for rooting more taxa*/
		if (d1<ntax) d1-=1;
		if (d2<ntax) d2-=1;
		if (temp_align) temp_align=dump_align(temp_align);
		if (values->align_cache) {
			/*make names*/
			if (strcmp (a[d1]->name,a[d2]->name) > 0) temp_name=other_pair_names(a[d2]->name,a[d1]->name);
			else temp_name=other_pair_names(a[d1]->name,a[d2]->name);
			num_seqs=1;
			for (j=0;j<strlen(temp_name);j++) if (temp_name[j]=='(') ++num_seqs;
			temp_align=retrieve_align_cache(temp_name,values,values->align_cache[num_seqs-2]);
			free(temp_name);
		}
		if (!temp_align) {
			temp_align = nw(a[d1],a[d2],values);
			if (values->align_cache) store_align_cache(temp_align,values,values->align_cache[num_seqs-2]);
		}
		/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"D");*/
		++nodes_actual;
		if (!a[ntax+n]) same=0;
		else {
			same=1;
			/*if misssing same*/
			if ((!only_missing(a[ntax+n])) && (!only_missing(temp_align))) for (j=0;j<temp_align->length;j++) if (a[ntax+n]->s[0][j]!=temp_align->s[0][j]) {
				same=0;
				break;
			}
		}
		changed[i]=(!same);
		if (!same) {
			if (a[ntax+n]) a[ntax+n]=dump_align(a[ntax+n]);
			a[ntax+n]=temp_align;
			temp_align=NULL;
			new_opt_score+=a[ntax+n]->score;
			if (new_opt_score>=bound) break;
		}
		else {
			for (j=i;j<num_to_do;j++) new_opt_score+=a[ntax+order[j]]->score;
			break;
		}
	}
	if (new_opt_score < bound) for (i=0;i<nent-2;i++) {
		/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"C%d",i);*/
		if (a_best[i]) a_best[i]=dump_align(a_best[i]);
		a_best[i]=make_align(a[i+ntax+1]);
		/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"P%d",i);*/
	}
	/*else if ((new_opt_score==bound) && PARALLEL*parallel_modify) {
		if (!(PARALLEL && values->jackboot)) fprintf(stderr,"C%d",i);
		if (a_best[i]) a_best[i]=dump_align(a_best[i]);
		a_best[i]=make_align(a[i+ntax+1]);
		if (!(PARALLEL && values->jackboot)) fprintf(stderr,"P%d",i);
		}*/
	if (temp_align) temp_align=dump_align(temp_align);
	/*
if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Nodes to do ");
	for (i=0;i<num_to_do;i++) if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d ",order[i]);
	if (!(PARALLEL && values->jackboot)) fprintf(stderr,"\n");
	for (i=1;i<nent-1;i++) {
		if (temp_align) temp_align=dump_align(temp_align);
		d1=nodes[i][0];
		d2=nodes[i][1];
		if (d1<ntax) d1-=1;
		if (d2<ntax) d2-=1;
		temp_align=nw(a[d1],a[d2],values);
		if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Node[%d] %d vs. %d\n",i,a[ntax+i]->score,temp_align->score);
		temp_align=dump_align(temp_align);
		}*/
	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d/%d ",nodes_actual,num_to_do);*/
	for (i=0;i<num_to_do;i++) {
		if (changed[i] && backup[i]) {
			a[order[i]+ntax]=dump_align(a[order[i]+ntax]);
			a[order[i]+ntax]=backup[i];
			backup[i]=NULL;
		}
		else if (backup[i]) backup[i]=dump_align(backup[i]);
	}
	free(backup);
	free(changed);
	return(new_opt_score);

}

int all_diagnose_tree_here3(a,nodes,ntax,nent,values,bound)
int **nodes;
int ntax,nent,bound;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,done,dont_do;
	char *check_names;
	int n_seqs,j;


	/*need to insure that we start with nodes which connect to terminal taxa and then nodes
	which have been defined - otherwise we'll be using the optimizations of
	previous topology's nodes*/
	/*add in a count for the first nent and nent-1 alignments and number of taxa for entering sequences with !=1 size*/

	check_names=NULL;
	if (nent==ntax) ++values->number;
	values->all_made=(char *)malloc(((2*ntax)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntax;n++) values->all_made[n]=1;
	for (n=ntax;n<(ntax+nent-1);++n) values->all_made[n]=0;

	new_opt_score=0;
	done=1;
	while (done) {
		done=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntax+n])) {
				values->all_made[ntax+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				/*insert name checkto make sure its not there already
				if (strcmp (a[d1]->name,a[d2]->name) > 0) check_names=other_pair_names(a[d2]->name,a[d1]->name);
				else check_names=other_pair_names(a[d1]->name,a[d2]->name);
				if ((a[ntax+n]) && (!strcmp(check_names,a[ntax+n]->name))) dont_do=1;
												else */dont_do=0;
				if (!dont_do) {
					if (a[ntax+n]) a[ntax+n]=dump_align(a[ntax+n]);
					if (values->align_cache) {
						n_seqs=1;
						for (j=0;j<strlen(check_names);j++) if (check_names[j]=='(') ++n_seqs;
						a[ntax+n]=retrieve_align_cache(check_names,values,values->align_cache[n_seqs-2]);
						free(temp_name);
					}
					if (!a[ntax+n]) {
						a[ntax+n] = nw(a[d1],a[d2],values);
						if (values->align_cache) store_align_cache(a[ntax+n],values,values->align_cache[n_seqs-2]);
					}
				}
				done=1;
				new_opt_score+=a[ntax+n]->score;
				free(check_names);
				if (new_opt_score>bound) goto out;
			}
		}
	}
out:
	;
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	return(new_opt_score);
}

int all_diagnose_tree_here(a,nodes,ntax,nent,values,bound)
int **nodes;
int ntax,nent,bound;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,done,dont_do;
	char *check_names;
	int n_seqs,j;


	/*need to insure that we start with nodes which connect to terminal taxa and then nodes
	which have been defined - otherwise we'll be using the optimizations of
	previous topology's nodes*/
	/*add in a count for the first nent and nent-1 alignments and number of taxa for entering sequences with !=1 size*/

	check_names=NULL;
	if (nent==ntax) ++values->number;
	values->all_made=(char *)malloc(((2*ntax)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntax;n++) values->all_made[n]=1;
	for (n=ntax;n<(ntax+nent-1);++n) values->all_made[n]=0;

	new_opt_score=0;
	done=1;
	while (done) {
		done=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntax+n])) {
				values->all_made[ntax+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				/*insert name checkto make sure its not there already*/
				if (strcmp (a[d1]->name,a[d2]->name) > 0) check_names=other_pair_names(a[d2]->name,a[d1]->name);
				else check_names=other_pair_names(a[d1]->name,a[d2]->name);
				if ((a[ntax+n]) && (!strcmp(check_names,a[ntax+n]->name))) dont_do=1;
				else dont_do=0;
				if (!dont_do) {
					if (a[ntax+n]) a[ntax+n]=dump_align(a[ntax+n]);
					if (values->align_cache) {
						n_seqs=1;
						for (j=0;j<strlen(check_names);j++) if (check_names[j]=='(') ++n_seqs;
						a[ntax+n]=retrieve_align_cache(check_names,values,values->align_cache[n_seqs-2]);
						free(temp_name);
					}
					if (!a[ntax+n]) {
						a[ntax+n] = nw(a[d1],a[d2],values);
						if (values->align_cache) store_align_cache(a[ntax+n],values,values->align_cache[n_seqs-2]);
					}
				}
				done=1;
				new_opt_score+=a[ntax+n]->score;
				free(check_names);
				if (new_opt_score>bound) goto out;
			}
		}
	}
out:
	;
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	return(new_opt_score);
}

alignment **do_swap_thang_align_dave(a,ntax,nodes,l_in,best_aligns,values,anc,mult,parallel_modify)
alignment **a, **best_aligns;
int ntax, *l_in, **nodes, parallel_modify;
parameters *values;
int *anc,mult;
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
	alignment *temp_align;
	int where_it_was_branch,paren_count;
	int is_unique,sent,received;
	int ii,jj,this_tree_counter,this_best;
	int bytes,info,bufid,source,type;
	int will_be_OK,**cur_groups,ngroups2;
	int num_best, all_rest;
	alignment **best_aligns_holder,**check_buf;
	int x,y,up_node,down_node;
	int int_holder,in_trees,g_test;
	char *collapse_name;
	int done,same,num_to_do;
	int *order,*dumb_sup,temp_length,initial;
	int can_move,a_group_same,second_group_same,found,match,*temp_group;
	int anc_counter,rem_counter;

	initial=1;

	if ((values->ngroups>0)	&& (values->groups_as_start)) {
		for (i=0;i<values->ngroups;i++) free(values->groups[i]);
		free(values->groups);
		values->groups=NULL;
		values->ngroups=0;
		fprintf(stderr,"	Freeing groups\n");
		if (PARALLEL) pvm_mcast(&(values->tids[1]),values->num_hosts-1, PVM_DUMP_INITIAL_GROUPS);
	}


	check_buf=NULL;
	collapse_name=NULL;
	values->support=NULL;
	/*allocations*/
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
	paren_count=0;
	temp_align=NULL;
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
			/*copy nodes*/
			for (k=0;k<ntax-1;k++) {
				hold_nodes[k][0]=best_nodes_old[m][k][0];
				hold_nodes[k][1]=best_nodes_old[m][k][1];
				anc[hold_nodes[k][0]]=anc[hold_nodes[k][1]]=k+ntax;
			}
			anc[ntax]=ntax;
			/*initialize buffer with down passes should be PARALLEL*parallel_modify if can*/
			if (values->new_optimization) {
				/*deresolves incomming tree on mult pass of swap*/
				if (mult && (values->collapse!=2)) {
					g_test=0;
					for (k=0;k<ntax;k++) g_test+=values->support[0][k];
					/*collapse and get groups*/
					if (!g_test) check=get_collapsed_thang_with_groups_old(best_aligns[0],hold_nodes,ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
				}
			}
			sent=received=0;
			for (i=2;i<2*(ntax-1);i++) {
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
				blocked[1]=0;
				node_removed=anc[cur_nodes[d1][d2]];
				/*redo nodes*/
				where_it_was_anc=anc[node_removed];
				if (cur_nodes[anc[node_removed]-ntax][0]==node_removed) {
					cur_nodes[anc[node_removed]-ntax][0]=cur_nodes[d1][!d2];
					/*blocked[2*(anc[node_removed]-ntax)]=1;*/
					where_it_was_desc=cur_nodes[anc[node_removed]-ntax][0];
					where_it_was_branch=0;
				}
				else {
					cur_nodes[anc[node_removed]-ntax][1]=cur_nodes[d1][!d2];
					/*blocked[2*(anc[node_removed]-ntax)+1]=1;*/
					where_it_was_desc=cur_nodes[anc[node_removed]-ntax][1];
					where_it_was_branch=1;
				}
				anc[cur_nodes[d1][!d2]]=anc[node_removed];
				cur_nodes[d1][!d2]=cur_nodes[d1][d2];
				anc[node_removed]=node_removed;
				if (num_clipped<=(ntax-3)) can_move=1;
				else can_move=0;
				/*if (0) { */
				if (values->groups && can_move) {
					/*fprintf(stderr,"Checking...");*/
					/*if taxon can't move--don't swap it*/
					/*if all descendents of ancestor are a group AND complement of group a goup (or single) ==> no swap);
			    /*get groups of base tree*/
					can_move=1;
					cur_groups=(int **)malloc(2*sizeof(int *));
					assert((int) cur_groups);
					for (ii=0;ii<2;ii++) {
						cur_groups[ii]=(int *)malloc(ntax*sizeof(int));
						assert((int) cur_groups[ii]);
						for (j=0;j<ntax;j++) cur_groups[ii][j]=0;
					}
					/* get group of anc*/
					a_group_same=0;
					recurse_groups2(hold_nodes,ntax,ntax,d1,cur_groups[0]);
					match=-1;
					for (k=0;k< values->ngroups;k++) {
						found=1;
						for (j=0;j<ntax-1;j++) if (cur_groups[0][j]!=values->groups[k][j]) {
							found=0;
							break;
						}
						if (found==1) {
							a_group_same=1;
							match=k;
							break;
						}
					}
					/*fprintf(stderr,"Gone through first %d(%d)...",match,a_group_same);*/
					if (a_group_same) {
						if (hold_nodes[d1][!d2]<ntax) can_move=0; /* if sister singleton=>can't move*/
else {
	recurse_groups2(hold_nodes,ntax,ntax,hold_nodes[d1][!d2]-ntax,cur_groups[1]);
	second_group_same=0;
	for (k=0;k< values->ngroups;k++) {
		found=1;
		for (j=0;j<ntax-1;j++) if (cur_groups[1][j]!=values->groups[k][j]) {
			found=0;
			break;
		}
		if (found==1) {
			second_group_same=1;
			break;
		}
	}
	if (second_group_same) can_move=0;
}
					}
					for (ii=0;ii<2;ii++) free(cur_groups[ii]);
					free(cur_groups);
					if (!can_move) fprintf(stderr,"X");
					/*if (values->atbr && (num_clipped > 2)) { if more than 2 and tbr--leave reroot possibility internally REALLY ALREADY DONE IN PREVIOUS SWAP
        	            if (!can_move) {
        	                can_move=1;
                            for (j=1;j<(2*(ntax-1));j++) blocked[j]=1;
                            if (d2) blocked[i-1]=0;
                            else blocked[i+1]=0;
        	                }
    				    }*/
				}
				if (can_move) {
					order=(int *)malloc(ntax*sizeof(int));
					assert((int)order);
					if (PARALLEL*parallel_modify) {
						if (sent < (values->num_hosts-1)) {
							pvm_initsend( PvmDataDefault );
							/*pack stuff*/
							pvm_pkint(&num_clipped,1,1);
							for (ii=0;ii<ntax-1;ii++) pvm_pkint(cur_nodes[ii],2,1);
							pvm_pkint(anc,((2*ntax)-1),1);
							pvm_pkint(&d1,1,1);
							pvm_pkint(&d2,1,1);
							pvm_pkint(&old_length,1,1);
							pvm_pkint(blocked,(2*(ntax-1)),1);
							pvm_pkint(&node_removed,1,1);
							pvm_pkint(&mult,1,1);
							pvm_pkint(is_it_or_desc,((2*ntax)-1),1);
							for (ii=ntax+1;ii<2*ntax-1;ii++) pack_align(a[ii]);
							pvm_send(values->tids[++sent],PVM_NEW_SWAP_DAVE);
						}
						else {
							bufid=pvm_recv(-1,PVM_NEW_SWAP_DONE_DAVE);
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
										if (values->poopsy && (!values->jackboot)) fprintf(stderr,"%s\n",best_aligns[0]->name);
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
												if (values->collapse==2) is_unique=is_unique_tree(temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
												else is_unique=is_unique_tree_align(best_aligns,temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
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
							pvm_pkint(&old_length,1,1);
							pvm_pkint(blocked,(2*(ntax-1)),1);
							pvm_pkint(&node_removed,1,1);
							pvm_pkint(&mult,1,1);
							pvm_pkint(is_it_or_desc,((2*ntax)-1),1);
							for (ii=ntax+1;ii<2*ntax-1;ii++) pack_align(a[ii]);
							pvm_send(source,PVM_NEW_SWAP_DAVE);
							++sent;
						}
					}
					else {
						/*when one or two taxa*/
						if ((num_clipped < 3) || (!values->atbr)) { /*if no alternative roots or sbr*/
							/*get nodes of clipped bit and rest*/
							all_diagnose_tree_here_pieces(a,cur_nodes,ntax,ntax,values,node_removed);
							/*then do the adding thing loop through places to put back*/
							best=0;
							for (j=1;j<ntax-1;j++) if ((j!=(node_removed-ntax)) && (a[ntax+j])) best+=a[ntax+j]->score;
							best=old_length-best;
							if ((best>0) || ((best==0) && (mult))) {

								for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
									for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=cur_nodes[ii][jj];
									check=temp_nodes[j/2][j%2];
									temp_nodes[j/2][j%2]=node_removed;
									temp_nodes[d1][!d2]=check;
									if (values->groups) {
										cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
										will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
										for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
										free(cur_groups);
									}
									if (will_be_OK) {
										/*set up nodes to do and order*/
										for (k=0;k<ntax;k++) order[k]=(-1);
										order[0]=node_removed-ntax;
										num_to_do=1;
										/*if (j>1) {
								order[1]=(j/2);
								num_to_do=2;
								}*/
										done=0;
										while (!done) {
											done=1;
											for (k=1;k<ntax-1;k++) {
												same=0;
												for (kk=0;kk<num_to_do;kk++) if (order[kk]==k) same=1;
												if (!same) {
													for (kk=0;kk<num_to_do;kk++) {
														if (temp_nodes[k][0]==(order[kk]+ntax)) {
															order[num_to_do++]=k;
															done=0;
														}
														else if (temp_nodes[k][1]==(order[kk]+ntax)) {
															order[num_to_do++]=k;
															done=0;
														}
													}
												}
											}
										}
										/*do_alignment*/
										cur_best=all_diagnose_tree_cut_off_spr(a,temp_nodes,ntax,ntax,values,order,num_to_do,tot_best);
										/*if (cur_best <= tot_best) {
							    temp_length=all_diagnose_tree_here3(check_buf,temp_nodes,ntax,ntax,values,HUGE_COST);
    							if (temp_length!=cur_best) {
    							        fprintf(stderr,"(%d,%d)",temp_length,cur_best);
    							        cur_best=temp_length;
    							        }
    						}*/
										if (cur_best < tot_best) {
											if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",cur_best);
											tot_best=cur_best;
											found_better=1;
											for (k=0;k<ntax-1;k++) {
												best_nodes_new[0][k][0]=temp_nodes[k][0];
												best_nodes_new[0][k][1]=temp_nodes[k][1];
												/*fprintf(stderr,"([%d] %d %d)",k,temp_nodes[k][0],temp_nodes[k][1]);*/
											}
											/*add best_aligns stuf free first*/
											for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
											best_aligns[0]=make_align(a[best_nodes_new[0][0][1]]);
											best_aligns[0]->score=tot_best;
											if (values->new_optimization && (values->collapse != 2)) {
												/*redo name for collapse*/
												if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values);
												else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values,values->support[0]);
												best_aligns[0]->score=tot_best;
											}
											tree_counter=1;
											if (values->poopsy && (!values->jackboot)) fprintf(stderr,"%s\n",best_aligns[0]->name);
										}
										else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
											/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"F");*/
											for (k=0;k<ntax-1;k++) {
												best_nodes_new[tree_counter][k][0]=temp_nodes[k][0];
												best_nodes_new[tree_counter][k][1]=temp_nodes[k][1];
											}
											/*chec to see if novel*/
											is_unique=1;
											if (values->collapse==2) is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax);
											else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
											if (is_unique) {
												if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
												found_better=1;
												best_aligns[tree_counter]=make_align(a[best_nodes_new[tree_counter][0][1]]);
												best_aligns[tree_counter]->score=tot_best;
												if (values->new_optimization && (values->collapse!=2)) {
													free(best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
													assert((int) best_aligns[tree_counter]->name);
													best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name);
												}
												++tree_counter;
											}
											/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"O");*/
										}
									} /*OK*/
								} /*regraft 'j'*/
							}/*best*/
							/*else {values->number+=((2*(ntax-num_clipped))-3);}*/
						}/*too fewipped or SBR*/
else { /*TBR*/
	check=get_root_array(reroot_array,num_clipped,node_removed,cur_nodes,d1,d2,ntax,is_it_or_desc,anc);
	for (l=0;l<2*num_clipped-3;l++) {
		/*get nodes of clipped bit and rest*/
		all_diagnose_tree_here_pieces(a,reroot_array[l],ntax,ntax,values,node_removed);
		/*then do the adding thing loop through places to put back*/
		best=0;
		for (j=1;j<ntax-1;j++) if ((j!=(node_removed-ntax)) && (a[ntax+j])) best+=a[ntax+j]->score;
		best=old_length-best;
		if ((best>0) || ((best==0) && (mult))) {
			for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
				for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=reroot_array[l][ii][jj];
				check=temp_nodes[j/2][j%2];
				temp_nodes[j/2][j%2]=node_removed;
				temp_nodes[d1][!d2]=check;
				if (values->groups) {
					cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
					will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
					for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
					free(cur_groups);
				}
				if (will_be_OK) {
					/*set up nodes to do and order*/
					for (k=0;k<ntax;k++) order[k]=(-1);
					order[0]=node_removed-ntax;
					num_to_do=1;
					/*if (j>1) {
										order[1]=(j/2);
										num_to_do=2;
										}*/
					done=0;
					while (!done) {
						done=1;
						for (k=1;k<ntax-1;k++) {
							same=0;
							for (kk=0;kk<num_to_do;kk++) if (order[kk]==k) same=1;
							if (!same) {
								for (kk=0;kk<num_to_do;kk++) {
									if (temp_nodes[k][0]==(order[kk]+ntax)) {
										order[num_to_do++]=k;
										done=0;
									}
									else if (temp_nodes[k][1]==(order[kk]+ntax)) {
										order[num_to_do++]=k;
										done=0;
									}
								}
							}
						}
					}
					/*do_alignment*/
					cur_best=all_diagnose_tree_cut_off_spr(a,temp_nodes,ntax,ntax,values,order,num_to_do,tot_best);
					/*if (cur_best <= tot_best) {
							            temp_length=all_diagnose_tree_here3(check_buf,temp_nodes,ntax,ntax,values,HUGE_COST);
    							        if (temp_length!=cur_best) {
    							            fprintf(stderr,"(%d,%d)",temp_length,cur_best);
    							            cur_best=temp_length;
    							            }
    						            }*/
					if (cur_best < tot_best) {
						if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found better at %d\n",cur_best);
						tot_best=cur_best;
						found_better=1;
						for (k=0;k<ntax-1;k++) {
							best_nodes_new[0][k][0]=temp_nodes[k][0];
							best_nodes_new[0][k][1]=temp_nodes[k][1];
						}
						/*add best_aligns stuf free first*/
						for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
						best_aligns[0]=make_align(a[best_nodes_new[0][0][1]]);
						best_aligns[0]->score=tot_best;
						if (values->new_optimization && (values->collapse != 2)) {
							/*redo name for collapse*/
							if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values);
							else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*0,values,values->support[0]);
							best_aligns[0]->score=tot_best;
						}
						tree_counter=1;
						if (values->poopsy && (!values->jackboot)) fprintf(stderr,"%s\n",best_aligns[0]->name);
					}
					else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
						/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"F");*/
						for (k=0;k<ntax-1;k++) {
							best_nodes_new[tree_counter][k][0]=temp_nodes[k][0];
							best_nodes_new[tree_counter][k][1]=temp_nodes[k][1];
						}
						/*chec to see if novel*/
						is_unique=1;
						if (values->collapse==2) is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax);
						else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
						if (is_unique) {
							if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d\n",tree_counter+1);
							found_better=1;
							best_aligns[tree_counter]=make_align(a[best_nodes_new[tree_counter][0][1]]);
							best_aligns[tree_counter]->score=tot_best;
							if (values->new_optimization && (values->collapse!=2)) {
								free(best_aligns[tree_counter]->name);
								best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
								assert((int) best_aligns[tree_counter]->name);
								best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name);
							}
							++tree_counter;
						}
						/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"O");*/
					}
				} /*OK*/
			} /*regraft 'j'*/
		}/*best*/
		/*else {values->number+=((2*(ntax-num_clipped))-3);}*/
	}/*rerootings*/
}/*TBR*/
					}/*PARALLEL*parallel_modify*/
					if (order) free(order);
				} /*size of plucked OK*/
			} /*i clades to pick*/
			if (PARALLEL*parallel_modify) { /*receive the rest after all I's are sent out and some received */
				/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"In cleanup ");*/
				while (received<sent) {
					bufid=pvm_recv(-1,PVM_NEW_SWAP_DONE_DAVE);
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
								if (values->poopsy && (!values->jackboot)) fprintf(stderr,"%s\n",best_aligns[0]->name);
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
									if (!mult) check=get_collapsed_thang(best_aligns[ii],best_nodes_new[ii],ntax,check_buf,PARALLEL*parallel_modify*0,values);
									else check=get_collapsed_thang_with_groups_old(best_aligns[ii],best_nodes_new[ii],ntax,check_buf,PARALLEL*parallel_modify*0,values,values->support[0]);
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
			else is_unique=is_unique_tree_align(best_aligns,temp_nodes_h,best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*0,values,collapse_name);
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
		if (!mult) new_start=0;
		/*reroot swapping*/
		if (values->arrt) {
			if (found_better || initial) { /*found new or original*/
				initial=0;
				in_trees=tree_counter;
				if (!(PARALLEL && values->jackboot)) fprintf(stderr,"*");
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
					for (l=0;l<(2*(ntax-1))-3;l++) {
						if (!PARALLEL*parallel_modify) cur_best=all_diagnose_tree_here3(a,reroot_array[l],ntax,ntax,values,tot_best);
						else cur_best=all_diagnose_tree_here_parallel3(a,reroot_array[l],ntax,ntax,values,tot_best);
						/*for (x=0;x<ntax-1;x++) if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d %d,",reroot_array[l][x][0],reroot_array[l][x][1]);*/
						/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d ",cur_best);*/
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
							best_aligns[0]=make_align(a[reroot_array[l][0][1]]);
							best_aligns[0]->score=tot_best;
							if (values->new_optimization && (values->collapse != 2)) {
								/*redo name for collapse*/
								if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values);
								else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*parallel_modify*1,values,values->support[0]);
								best_aligns[0]->score=tot_best;
							}
							tree_counter=1;
							if (values->poopsy && (!values->jackboot)) fprintf(stderr,"%s\n",best_aligns[0]->name);
						}
						else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
							for (k=0;k<ntax-1;k++) {
								best_nodes_new[tree_counter][k][0]=reroot_array[l][k][0];
								best_nodes_new[tree_counter][k][1]=reroot_array[l][k][1];
							}
							/*chec to see if novel*/
							is_unique=1;
							if (values->collapse==2) is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax);
							else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length,tot_best,ntax,check_buf,PARALLEL*parallel_modify*1,values,collapse_name);
							if (is_unique) {
								if (!(PARALLEL && values->jackboot)) fprintf(stderr," Found another for %d (R)\n",tree_counter+1);
								found_better=1;
								best_aligns[tree_counter]=make_align(a[reroot_array[l][0][1]]);
								best_aligns[tree_counter]->score=tot_best;
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
				/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"about to copy back ");*/

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
	/*do names if collapse 2*/
	if (values->new_optimization && (values->collapse==2)) {
		for (k=0;k<tree_counter;k++) {
			free(best_aligns[k]->name);
			dumb_sup=(int *)malloc(ntax*sizeof(int));
			assert((int) dumb_sup);
			for (l=0;l<ntax;l++) dumb_sup[l]=1;
			best_aligns[k]->name=(char *)get_names_with_collapse(a,best_nodes_new[k],ntax,ntax,values,dumb_sup);
			free(dumb_sup);
		}
	}
	for (i=0;i<paren_count;i++)     if (!(PARALLEL && values->jackboot)) fprintf(stderr,")");
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
	for (i=0;i<cur_num_trees;i++) {
		for (j=0;j<ntax-1;j++) {
		    for (k=0;k<2;k++) {
		      fprintf(stderr,"BNN[%d][%d][%d]=%d\n",i,j,k,best_nodes_old[i][j][k]);
		      best_nodes[i][j][k]=best_nodes_old[i][j][k];
		    }
		  }}
	

	/*filter*/
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
	if (values->support) {
		for (i=0;i<values->keep_aligns;i++) free(values->support[i]);
		free(values->support);
		values->support=NULL;
	}
	if (temp_align) temp_align=dump_align(temp_align);
	*l_in=tree_counter;
	for (i=0;i<2*ntax-1;i++) if (check_buf[i]) check_buf[i]=dump_align(check_buf[i]);
	free(check_buf);
	return best_aligns;
}

void all_diagnose_tree_here_pieces(a,nodes,ntax,nent,values,base)
int **nodes;
int ntax,nent,base;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,did_one,dont_do;
	char *check_names;

	values->all_made=(char *)malloc(((2*ntax)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntax;n++) values->all_made[n]=1;
	for (n=ntax;n<(ntax+nent-1);++n) values->all_made[n]=0;
	/*values->all_made[base]=1;*/
	/*if (a[base]) a[base]->score=0;*/

	/*make sure the in-some check works for multiple input alignments*/
	if (!values->number_of_input_alignments) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization) && (values->chop>0)) nent_minus_one_new=nent-1;
	else for (n=0;n<(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	new_opt_score=0;
	did_one=1;
	while (did_one) {
		did_one=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntax+n])) {
				values->all_made[ntax+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				/*if (strcmp (a[d1]->name,a[d2]->name) > 0) check_names=other_pair_names(a[d2]->name,a[d1]->name);
				else check_names=other_pair_names(a[d1]->name,a[d2]->name);
				if (!strcmp(check_names,a[ntax+n]->name)) dont_do=1;
				else dont_do=0;
				free(check_names);*/
				dont_do=0;
				if (!dont_do) {
					if (a[ntax+n]) a[ntax+n]=dump_align(a[ntax+n]);
					a[ntax+n] = nw(a[d1],a[d2],values);
				}
				did_one=1;
			}
		}
	}
	free(values->all_made);
	values->all_made=NULL;
}

int all_diagnose_tree_cut_off_spr(a,nodes,ntax,nent,values,order,num_to_do,bound)
int **nodes;
int ntax,nent;
alignment **a;
parameters *values;
int *order, num_to_do,bound;
{
	int n,d1,d2,nent_minus_one_new;
	int new_opt_score,i,j;
	alignment *temp_align;
	alignment **backup;
	int *changed,same;
	int nodes_actual=0;
	int num_seqs;
	char *temp_name;

	if (nent==ntax) ++values->number;
	/*could move these allocations out/up*/
	backup=(alignment **)malloc(num_to_do*sizeof(alignment *));
	assert((int) backup);
	for (i=0;i<num_to_do;i++) backup[i]=NULL;
	changed=(int *)malloc(num_to_do*sizeof(int));
	assert((int)changed);
	for (i=0;i<num_to_do;i++) {
		if (a[ntax+order[i]]) backup[i]=make_align(a[ntax+order[i]]);
		changed[i]=0;
	}
	nent_minus_one_new=nent-1;
	new_opt_score=0;
	for (i=1;i<nent-1;i++) {
		same=0;
		for (j=0;j<num_to_do;j++) if (i==order[j]) same=1;
		if (!same) new_opt_score+=a[ntax+i]->score;
	}
	temp_align=NULL;
	for (i=0;i<num_to_do;i++){
		n=order[i];
		d1=nodes[n][0];
		d2=nodes[n][1];
		/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Making %d from %d and %d\n",order[i],d1,d2);*/
		/*fix for rooting more taxa*/
		if (d1<ntax) d1-=1;
		if (d2<ntax) d2-=1;
		if (temp_align) dump_align(temp_align);
		if (values->align_cache) {
			/*make names*/
			if (strcmp (a[d1]->name,a[d2]->name) > 0) temp_name=other_pair_names(a[d2]->name,a[d1]->name);
			else temp_name=other_pair_names(a[d1]->name,a[d2]->name);
			num_seqs=1;
			for (j=0;j<strlen(temp_name);j++) if (temp_name[j]=='(') ++num_seqs;
			temp_align=retrieve_align_cache(temp_name,values,values->align_cache[num_seqs-2]);
			free(temp_name);
		}
		if (!temp_align) {
			temp_align = nw(a[d1],a[d2],values);
			if (values->align_cache) store_align_cache(temp_align,values,values->align_cache[num_seqs-2]);
		}
		++nodes_actual;
		/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"a[ntax+%d]->score=%d ",n,temp_align->score);*/
		if (!a[ntax+n]) same=0;
		else {
			same=1;
			if ((!only_missing(a[ntax+n])) && (!only_missing(temp_align)))  for (j=0;j<temp_align->length;j++) if (a[ntax+n]->s[0][j]!=temp_align->s[0][j]) {
				same=0;
				break;
			}
		}
		changed[i]=(!same);
		if (!same) {
			if (a[ntax+n]) a[ntax+n]=dump_align(a[ntax+n]);
			a[ntax+n]=temp_align;
			temp_align=NULL;
			new_opt_score+=a[ntax+n]->score;
			if (new_opt_score>bound) break;
		}
		else {
			for (j=i;j<num_to_do;j++) new_opt_score+=a[ntax+order[j]]->score;
			break;
		}
	}
	if (temp_align) temp_align=dump_align(temp_align);
	/*
if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Nodes to do ");
	for (i=0;i<num_to_do;i++) if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%d ",order[i]);
	if (!(PARALLEL && values->jackboot)) fprintf(stderr,"\n");
	for (i=1;i<nent-1;i++) {
		if (temp_align) temp_align=dump_align(temp_align);
		d1=nodes[i][0];
		d2=nodes[i][1];
		if (d1<ntax) d1-=1;
		if (d2<ntax) d2-=1;
		temp_align=nw(a[d1],a[d2],values);
		if (!(PARALLEL && values->jackboot)) fprintf(stderr,"Node[%d] %d vs. %d\n",i,a[ntax+i]->score,temp_align->score);
		temp_align=dump_align(temp_align);
		}*/
	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"(%d/%d)",nodes_actual,num_to_do);*/
	for (i=0;i<num_to_do;i++) {
		if (changed[i] && backup[i]) {
			a[order[i]+ntax]=dump_align(a[order[i]+ntax]);
			a[order[i]+ntax]=backup[i];
			backup[i]=NULL;
		}
		else if (backup[i]) backup[i]=dump_align(backup[i]);
	}
	free(backup);
	free(changed);

	return(new_opt_score);

}

int all_diagnose_tree_here_parallel(a_down,nodes,ntax,nent,values,best_so_far)
alignment **a_down;
int **nodes, ntax,nent,best_so_far;
parameters *values;
{
	int n,i,found_one;
	alignment *temp_align;
	int d1,d2,d1_thang,d2_thang;
	int nent_minus_one_new,od1,od2;
	int *done;
	int sent,received,score_holder,node_to_be_sent;
	int node_received,bufid,info,source,tag,bytes,type;
	int ii,td1,td2,new_opt_score=0;
	int do_it;
	char *check_names=NULL;

	temp_align=NULL;
	done=(int *)malloc(((2*ntax)-1)*sizeof(int));
	assert((int) done);
	for (i=0;i<ntax;i++) done[i]=1;
	for (i=ntax;i<ntax+nent-1;i++) done[i]=0;

	if (nent==ntax) values->number++;

	/*down pass redone no pre-diagnose because the nodes are optimized!!*/
	values->in_optimize=0;
	found_one=1;
	while (found_one) {
		found_one=0;
		sent=received=0;
		for (n=(nent-2);n>0;--n){
			od1=d1=nodes[n][0];
			od2=d2=nodes[n][1];
			if ((done[d1]) && (done[d2]) && (!done[ntax+n])) {
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				node_to_be_sent=ntax+n;
				do_it=1;
				if (new_opt_score>best_so_far) do_it=0;
				else {
					if (a_down[ntax+n]) {
						if (strcmp(a_down[d1]->name,a_down[d2]->name)>0) check_names=other_pair_names(a_down[d2]->name,a_down[d1]->name);
						else check_names=other_pair_names(a_down[d1]->name,a_down[d2]->name);
						if (!strcmp(a_down[node_to_be_sent]->name,check_names)) {
							do_it=0;
							found_one=1;
							new_opt_score+=a_down[ntax+n]->score;
						}
						free(check_names);
					}
				}
				if (do_it) {
					found_one=1;
					if (sent < (values->num_hosts-1)) {
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&values->in_optimize,1,1);
						pvm_pkint(&node_to_be_sent,1,1);
						pack_align(a_down[d1]);
						pack_align(a_down[d2]);
						pvm_send(values->tids[++sent],PVM_PARALLEL_NW_NEW_OPT);
					}
					else {
						bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_DONE);
						if (bufid) {
							info=pvm_bufinfo(bufid, &bytes, &type, &source);
							pvm_upkint(&node_received,1,1);
							if (a_down[node_received]) a_down[node_received]=dump_align(a_down[node_received]);
							a_down[node_received]=unpack_align_and_score(values);
							done[node_received]=1;
							++received;
							new_opt_score+=a_down[node_received]->score;
						}
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&values->in_optimize,1,1);
						pvm_pkint(&node_to_be_sent,1,1);
						pack_align(a_down[d1]);
						pack_align(a_down[d2]);
						pvm_send(source,PVM_PARALLEL_NW_NEW_OPT);
						++sent;
					}
				}
			} /*if !done*/
		}/*n*/
		/*receive the remainder*/
		for (ii=received;ii<sent;ii++) {
			bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_DONE);
			if (bufid) {
				info=pvm_bufinfo(bufid, &bytes, &type, &source);
				pvm_upkint(&node_received,1,1);
				if (a_down[node_received]) a_down[node_received]=dump_align(a_down[node_received]);
				a_down[node_received]=unpack_align_and_score(values);
				done[node_received]=1;
				new_opt_score+=a_down[node_received]->score;
			}
		}
	}/*found one*/
	values->in_optimize=0;
	free(done);
	a_down[nodes[0][1]]->score=new_opt_score;
	return a_down[nodes[0][1]]->score;
}

int all_diagnose_tree_here_parallel3(a_down,nodes,ntax,nent,values,best_so_far)
alignment **a_down;
int **nodes, ntax,nent,best_so_far;
parameters *values;
{
	int n,i,found_one;
	alignment *temp_align;
	int d1,d2,d1_thang,d2_thang;
	int nent_minus_one_new,od1,od2;
	int *done;
	int sent,received,score_holder,node_to_be_sent;
	int node_received,bufid,info,source,tag,bytes,type;
	int ii,td1,td2,new_opt_score=0;
	int do_it;
	char *check_names=NULL;

	temp_align=NULL;
	done=(int *)malloc(((2*ntax)-1)*sizeof(int));
	assert((int) done);
	for (i=0;i<ntax;i++) done[i]=1;
	for (i=ntax;i<ntax+nent-1;i++) done[i]=0;

	if (nent==ntax) values->number++;

	/*down pass redone no pre-diagnose because the nodes are optimized!!*/
	values->in_optimize=0;
	found_one=1;
	while (found_one) {
		found_one=0;
		sent=received=0;
		for (n=(nent-2);n>0;--n){
			od1=d1=nodes[n][0];
			od2=d2=nodes[n][1];
			if ((done[d1]) && (done[d2]) && (!done[ntax+n])) {
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				node_to_be_sent=ntax+n;
				do_it=1;
				if (new_opt_score>best_so_far) do_it=0;

				if (do_it) {
					found_one=1;
					if (sent < (values->num_hosts-1)) {
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&values->in_optimize,1,1);
						pvm_pkint(&node_to_be_sent,1,1);
						pack_align(a_down[d1]);
						pack_align(a_down[d2]);
						pvm_send(values->tids[++sent],PVM_PARALLEL_NW_NEW_OPT);
					}
					else {
						bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_DONE);
						if (bufid) {
							info=pvm_bufinfo(bufid, &bytes, &type, &source);
							pvm_upkint(&node_received,1,1);
							if (a_down[node_received]) a_down[node_received]=dump_align(a_down[node_received]);
							a_down[node_received]=unpack_align_and_score(values);
							done[node_received]=1;
							++received;
							new_opt_score+=a_down[node_received]->score;
						}
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&values->in_optimize,1,1);
						pvm_pkint(&node_to_be_sent,1,1);
						pack_align(a_down[d1]);
						pack_align(a_down[d2]);
						pvm_send(source,PVM_PARALLEL_NW_NEW_OPT);
						++sent;
					}
				}
			} /*if !done*/
		}/*n*/
		/*receive the remainder*/
		for (ii=received;ii<sent;ii++) {
			bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_DONE);
			if (bufid) {
				info=pvm_bufinfo(bufid, &bytes, &type, &source);
				pvm_upkint(&node_received,1,1);
				if (a_down[node_received]) a_down[node_received]=dump_align(a_down[node_received]);
				a_down[node_received]=unpack_align_and_score(values);
				done[node_received]=1;
				new_opt_score+=a_down[node_received]->score;
			}
		}
	}/*found one*/
	values->in_optimize=0;
	free(done);
	a_down[nodes[0][1]]->score=new_opt_score;
	return a_down[nodes[0][1]]->score;
}

int all_diagnose_tree_here_parallel2(a_down,nodes,ntax,nent,values,best_so_far,a_down2)
alignment **a_down,**a_down2;
int **nodes, ntax,nent,best_so_far;
parameters *values;
{
	int n,i,found_one;
	alignment *temp_align;
	int d1,d2,d1_thang,d2_thang;
	int nent_minus_one_new,od1,od2;
	int *done;
	int sent,received,score_holder,node_to_be_sent;
	int node_received,bufid,info,source,tag,bytes,type;
	int ii,td1,td2,new_opt_score=0;
	int do_it;
	char *check_names=NULL;

	temp_align=NULL;
	done=(int *)malloc(((2*ntax)-1)*sizeof(int));
	assert((int) done);
	for (i=0;i<ntax;i++) done[i]=1;
	for (i=ntax;i<ntax+nent-1;i++) done[i]=0;

	if (nent==ntax) values->number++;

	/*down pass redone no pre-diagnose because the nodes are optimized!!*/
	values->in_optimize=1;
	found_one=1;
	while (found_one) {
		found_one=0;
		sent=received=0;
		for (n=(nent-2);n>0;--n){
			od1=d1=nodes[n][0];
			od2=d2=nodes[n][1];
			if ((done[d1]) && (done[d2]) && (!done[ntax+n])) {
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				node_to_be_sent=ntax+n;
				do_it=1;
				if (do_it) {
					found_one=1;
					if (sent < (values->num_hosts-1)) {
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&values->in_optimize,1,1);
						pvm_pkint(&node_to_be_sent,1,1);
						pack_align(a_down[d1]);
						pack_align(a_down[d2]);
						pvm_send(values->tids[++sent],PVM_PARALLEL_NW_NEW_OPT);
					}
					else {
						bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_DONE);
						if (bufid) {
							info=pvm_bufinfo(bufid, &bytes, &type, &source);
							pvm_upkint(&node_received,1,1);
							if (a_down[node_received]) a_down[node_received]=dump_align(a_down[node_received]);
							if (a_down2[node_received]) a_down2[node_received]=dump_align(a_down2[node_received]);
							a_down2[node_received]=unpack_align_and_score(values);
							a_down[node_received]=make_align(a_down2[node_received]);
							a_down[node_received]=make_ambig(a_down[node_received],values);
							a_down2[node_received] = make_three(a_down2[node_received] ,values);
							done[node_received]=1;
							++received;
							new_opt_score+=a_down[node_received]->score;
						}
						pvm_initsend( PvmDataDefault );
						pvm_pkint(&values->in_optimize,1,1);
						pvm_pkint(&node_to_be_sent,1,1);
						pack_align(a_down[d1]);
						pack_align(a_down[d2]);
						pvm_send(source,PVM_PARALLEL_NW_NEW_OPT);
						++sent;
					}
				}
			} /*if !done*/
		}/*n*/
		/*receive the remainder*/
		for (ii=received;ii<sent;ii++) {
			bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_DONE);
			if (bufid) {
				info=pvm_bufinfo(bufid, &bytes, &type, &source);
				pvm_upkint(&node_received,1,1);
				if (a_down[node_received]) a_down[node_received]=dump_align(a_down[node_received]);
				if (a_down2[node_received]) a_down2[node_received]=dump_align(a_down2[node_received]);
				a_down2[node_received]=unpack_align_and_score(values);
				a_down[node_received]=make_align(a_down2[node_received]);
				a_down[node_received]=make_ambig(a_down[node_received],values);
				a_down2[node_received] = make_three(a_down2[node_received] ,values);
				done[node_received]=1;
				new_opt_score+=a_down[node_received]->score;
			}
		}
	}/*found one*/
	values->in_optimize=0;
	free(done);
	a_down[nodes[0][1]]->score=new_opt_score;
	return a_down[nodes[0][1]]->score;
}

int all_diagnose_tree_here2(a,nodes,ntax,nent,values,bound,a2)
int **nodes;
int ntax,nent,bound;
alignment **a,**a2;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs,done,dont_do;
	char *check_names;
	int n_seqs,j;


	/*need to insure that we start with nodes which connect to terminal taxa and then nodes
	which have been defined - otherwise we'll be using the optimizations of
	previous topology's nodes*/
	/*add in a count for the first nent and nent-1 alignments and number of taxa for entering sequences with !=1 size*/

	check_names=NULL;
	if (nent==ntax) ++values->number;
	if (values->all_made) free (values->all_made);
	values->all_made=(char *)malloc(((2*ntax)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntax;n++) values->all_made[n]=1;
	for (n=ntax;n<(ntax+nent-1);++n) values->all_made[n]=0; 

	/*make sure the in-some check works for multiple input alignments*/
	if (values->number_of_input_alignments==0) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization) && (values->chop>0)) nent_minus_one_new=nent-1;
	else for (n=0;n<(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%u-",_memavl());*/

	new_opt_score=0;
	done=1;
	while (done) {
		done=0;
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntax+n])) {
				values->all_made[ntax+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntax) d1-=1;
				if (d2<ntax) d2-=1;
				/*insert name checkto make sure its not there already*/
				if (a[ntax+n]) a[ntax+n]=dump_align(a[ntax+n]);
				if (a2[ntax+n]) a2[ntax+n]=dump_align(a2[ntax+n]);
				values->in_optimize=1;
				a2[ntax+n] = nw(a[d1],a[d2],values); /*need to add the node result to avoid NW later in up pass*/
				values->in_optimize=0;
				a[ntax+n]=make_align(a2[ntax+n]);
				a[ntax+n]=make_ambig(a[ntax+n],values);
				a2[ntax+n] = make_three(a2[ntax+n],values);
				done=1;
				new_opt_score+=a[ntax+n]->score;
				free(check_names);
				if (new_opt_score>bound) goto out;
			}
		}
	}
out:
	;
	free(values->all_made);
	values->all_made=NULL;

	/*if (!(PARALLEL && values->jackboot)) fprintf(stderr,"%u)",_memavl());*/
	return(new_opt_score);
}

void recurse_groups2(node,ntaxa,nent,i,groups)
int **node;
int ntaxa,nent,i;
int *groups;
{
	int ii;
	/*fprintf(stderr,"R ");*/
	for (ii=0;ii<2;ii++){
		if (node[i][ii]<ntaxa) groups[node[i][ii]-1]=1;
		else if (node[i][ii]>ntaxa) recurse_groups2(node,ntaxa,nent,node[i][ii]-ntaxa,groups);
	}

}


