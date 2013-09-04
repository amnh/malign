/*Copyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Goloboff*/

#include "align3.h"

extern int first;
extern int charbit[256];
extern char bitchar[32];

int all_diagnose_tree_local_new_opt(a_down,nodes,ntax,nent,values,stuff_done,corres)
alignment **a_down;
int **nodes, ntax,nent,stuff_done,***corres;
parameters *values;
{
int n,i,found_one;
alignment *temp_align;
int d1,d2,d1_thang,d2_thang;
int nent_minus_one_new,od1,od2;
int *done;

temp_align=NULL;
done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;

/*down pass redone*/
values->in_optimize=1;
found_one=1;
while (found_one) {
	found_one=0;
	for (n=(nent-2);n>0;--n){
		od1=d1=nodes[n][0];
		od2=d2=nodes[n][1];
		if ((done[d1]) && (done[d2]) && (!done[ntax+n])) {
			if (d1 < ntax) d1-=1;
			if (d2 < ntax) d2 -=1;
			found_one=1;
			done[ntax+n]=1;
			if (a_down[ntax+n]) a_down[ntax+n]=dump_align(a_down[ntax+n]);
			/*get alignment
			fprintf(stderr,"{%d %d=>%d}",d1,d2,ntax+n);
			for (i=0;i<a_down[d1]->length;i++) fprintf(stderr,"%d",a_down[d1]->s[1][i]);
			fprintf(stderr,"\n");
			for (i=0;i<a_down[d2]->length;i++) fprintf(stderr,"%d",a_down[d2]->s[1][i]);
			fprintf(stderr,"\n");*/
			a_down[ntax+n] = nw(a_down[d1],a_down[d2],values);
			if (temp_align) temp_align=dump_align(temp_align);
			temp_align=make_align(a_down[ntax+n]);
			/*get ancestor with gaps*/
			temp_align=new_make_ambig(temp_align,values);
			/*get correspondance arrays*/
			if (od1!=od2) {
				if (strcmp(a_down[d1]->name,a_down[d2]->name)>0) {
					d2_thang=0;
					d1_thang=1;
					}
				else {
					d2_thang=1;
					d1_thang=0;
					}
				if (corres[od1]) {
					free(corres[od1][0]);
					free(corres[od1][1]);
					free(corres[od1]);
					}
				if (corres[od2]) {
					free(corres[od2][0]);
					free(corres[od2][1]);
					free(corres[od2]);
					}
				corres[od1]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[d1_thang],temp_align->length);
				corres[od2]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[d2_thang],temp_align->length);
				if (temp_align->type_weight) {
					modify_corres(corres[od1],temp_align);
					modify_corres(corres[od2],temp_align);
					}
				}
			else {
				if (corres[od1]) {
					free(corres[od1][0]);
					free(corres[od1][1]);
					free(corres[od1]);
					}
				corres[od1]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[0],temp_align->length);
				if (temp_align->type_weight) {
					modify_corres(corres[od1],temp_align);
					}
			}
			/*get ancestor without gaps*/
		    	a_down[ntax+n]=make_ambig(a_down[ntax+n],values);
		  	a_down[ntax+n]->score += (a_down[d1]->score+a_down[d2]->score);
		  	if (temp_align) temp_align=dump_align(temp_align);
	    		}
		}
	}
values->in_optimize=0;
if (temp_align) temp_align=dump_align(temp_align);
free(done);
return a_down[nodes[0][1]]->score;
}


int all_diagnose_tree_local_parallel_new_opt(a_down,nodes,ntax,nent,values,stuff_done,corres)
alignment **a_down;
int **nodes, ntax,nent,stuff_done,***corres;
parameters *values;
{
int n,i,found_one;
alignment *temp_align;
int d1,d2,d1_thang,d2_thang;
int nent_minus_one_new,od1,od2;
int *done;
int sent,received,score_holder,node_to_be_sent;
int node_received,bufid,info,source,tag,bytes,type;
int ii,td1,td2;

temp_align=NULL;
done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;

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
			found_one=1;
			if (d1 < ntax) d1-=1;
			if (d2 < ntax) d2 -=1;
			node_to_be_sent=ntax+n;
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
					if (temp_align) temp_align=dump_align(temp_align);
					temp_align=make_align(a_down[node_received]);
					/*get ancestor with gaps*/
					temp_align=new_make_ambig(temp_align,values);
					/*get correspondance arrays*/
					od1=td1=nodes[node_received-ntax][0];
					od2=td2=nodes[node_received-ntax][1];
					if (td1 < ntax) td1-=1;
					if (td2 < ntax) td2-=1;
					if (strcmp(a_down[td1]->name,a_down[td2]->name)>0) {
						d2_thang=0;
						d1_thang=1;
						}
					else {
						d2_thang=1;
						d1_thang=0;
						}
					if (corres[od1]) {
						free(corres[od1][0]);
						free(corres[od1][1]);
						free(corres[od1]);
						}
					if (corres[od2]) {
						free(corres[od2][0]);
						free(corres[od2][1]);
						free(corres[od2]);
						}
					corres[od1]=get_correspondances(temp_align->s[0],a_down[node_received]->s[d1_thang],temp_align->length);
					corres[od2]=get_correspondances(temp_align->s[0],a_down[node_received]->s[d2_thang],temp_align->length);
					if (temp_align->type_weight) {
						modify_corres(corres[od1],temp_align);
						modify_corres(corres[od2],temp_align);
						}
					/*get ancestor without gaps*/
				    	a_down[node_received]=make_ambig(a_down[node_received],values);
				  	a_down[node_received]->score += (a_down[td1]->score+a_down[td2]->score);
				  	if (temp_align) temp_align=dump_align(temp_align);
					}
				pvm_initsend( PvmDataDefault );
				pvm_pkint(&values->in_optimize,1,1);
				pvm_pkint(&node_to_be_sent,1,1);
				pack_align(a_down[d1]);
				pack_align(a_down[d2]);
				pvm_send(source,PVM_PARALLEL_NW_NEW_OPT);
				++sent;
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
					if (temp_align) temp_align=dump_align(temp_align);
					temp_align=make_align(a_down[node_received]);
					/*get ancestor with gaps*/
					temp_align=new_make_ambig(temp_align,values);
					/*get correspondance arrays*/
					od1=td1=nodes[node_received-ntax][0];
					od2=td2=nodes[node_received-ntax][1];
					if (td1 < ntax) td1-=1;
					if (td2 < ntax) td2-=1;
					if (strcmp(a_down[td1]->name,a_down[td2]->name)>0) {
						d2_thang=0;
						d1_thang=1;
						}
					else {
						d2_thang=1;
						d1_thang=0;
						}
					if (corres[od1]) {
						free(corres[od1][0]);
						free(corres[od1][1]);
						free(corres[od1]);
						}
					if (corres[od2]) {
						free(corres[od2][0]);
						free(corres[od2][1]);
						free(corres[od2]);
						}
					corres[od1]=get_correspondances(temp_align->s[0],a_down[node_received]->s[d1_thang],temp_align->length);
					corres[od2]=get_correspondances(temp_align->s[0],a_down[node_received]->s[d2_thang],temp_align->length);
					if (temp_align->type_weight) {
						modify_corres(corres[od1],temp_align);
						modify_corres(corres[od2],temp_align);
						}
					/*get ancestor without gaps*/
				    	a_down[node_received]=make_ambig(a_down[node_received],values);
				  	a_down[node_received]->score += (a_down[td1]->score+a_down[td2]->score);
				  	if (temp_align) temp_align=dump_align(temp_align);
				}
			}
	}/*found one*/
values->in_optimize=0;
if (temp_align) temp_align=dump_align(temp_align);
free(done);
return a_down[nodes[0][1]]->score;
}

void get_up_pass_new_opt2(a_final,a_down,a_up,nodes,ntax,values,nent,anc,corres,node_removed,desc)
alignment **a_down, ***a_up, **a_final;
int **nodes, ntax,nent, *anc,***corres,node_removed,*desc;
parameters *values;
{
int i,d1,d2,j,k,found_one,n;
alignment *temp1,*temp2,*temp3,*temp4;
int node_is_from;
char *tempstr,*t0,*t1,*s0,*s1;
int nent_minus_one_new,od1,od2;
int *done;

temp1=NULL;
temp2=NULL;
temp3=NULL;
temp4=NULL;


done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;
done[node_removed]=1;
for (i=0;i<2*ntax-1;i++) if (desc[i]) done[i]=1;

/*deallocate 'em all*/
for (i=1;i<nent-1;i++) {
	if (a_up[i+ntax][0]) a_up[i+ntax][0]=dump_align(a_up[i+ntax][0]);
	if (a_up[i+ntax][1]) a_up[i+ntax][1]=dump_align(a_up[i+ntax][1]);
	}

d1=nodes[nodes[0][1]-ntax][1];
d2=nodes[nodes[0][1]-ntax][0];
if (d1 < ntax) d1-=1;
if (d2 < ntax) d2-=1;
a_up[nodes[0][1]][0]=make_align(a_down[d1]);
a_up[nodes[0][1]][1]=make_align(a_down[d2]);

if (a_final[nodes[0][1]]) dump_align(a_final[nodes[0][1]]);
a_final[nodes[0][1]]=make_align(a_down[nodes[0][1]]);

done[nodes[0][1]]=1;

found_one=1;
while (found_one) {
	found_one=0;
	for (n=1;n<nent-1;n++) if ((!done[n+ntax]) && (done[anc[n+ntax]])) {
		found_one=1;
		done[n+ntax]=1;
		if (nodes[anc[ntax+n]-ntax][0]==(ntax+n)) node_is_from=0;
		else node_is_from=1;
		od1=d1=nodes[n][0];
		od2=d2=nodes[n][1];
		if (d1 < ntax) d1-=1;
		if (d2 < ntax) d2-=1;
		/*get 4-thang of node,anc-up,d1-down, and d2-down for up-pass*/
		if (temp1) temp1=dump_align(temp1);
		if (temp2) temp2=dump_align(temp2);
		if (temp3) temp3=dump_align(temp3);
		if (temp4) temp4=dump_align(temp4);
		values->in_optimize=1;
		/*need to make sure the correct taxon order for make neighbor later*/
		/*can this be changed to Fake NW?? with make_neigbor and a convert*/
		temp1=nw(a_down[ntax+n],a_up[anc[ntax+n]][node_is_from],values); /* no corres here because it would have to go through the ancestral down pass node*/
		/*switch if 0 and 1 not correct should be descendent then ancestor in each one*/
		/*if (strcmp(temp1->taxon_name[0],a_down[ntax+n]->taxon_name[0])) {*/
		if ((strcmp(a_down[ntax+n]->name,a_up[anc[ntax+n]][node_is_from]->name)>0)) {
			/*fprintf(stderr,"Switching..");*/
			tempstr=(char *)malloc((1+temp1->length)*sizeof(char));
			assert((int) tempstr);
			tempstr=(char *)strcpy(tempstr,temp1->s[0]);
			temp1->s[0]=(char *)strcpy(temp1->s[0],temp1->s[1]);
			temp1->s[1]=(char *)strcpy(temp1->s[1],tempstr);
			free(tempstr);
			/*fprintf(stderr,"Done ");*/
			}
		temp2=fake_nw(a_down[d1],a_down[ntax+n],values,corres[od1]);
		temp3=fake_nw(a_down[d2],a_down[ntax+n],values,corres[od2]);
		if (temp1->type_weight) temp4=big_make_neighbor(temp1,temp2,temp3,values);
		else temp4=make_neighbor(temp1,temp2,temp3,values);
		values->in_optimize=0;
		/*makem and reallocate to get rid of non-relavent strings then make-ambig*/
		a_up[ntax+n][0]=make_align(temp4);
		convert_to_anc_and_desc(a_up[ntax+n][0],1,2);
		a_up[ntax+n][0]=make_ambig(a_up[ntax+n][0],values);
		a_up[ntax+n][1]=make_align(temp4);
		convert_to_anc_and_desc(a_up[ntax+n][1],1,3);
		a_up[ntax+n][1]=make_ambig(a_up[ntax+n][1],values);
		/*make final for a also*/
		if (a_final[ntax+n]) a_final[ntax+n]=dump_align(a_final[ntax+n]);
		/* final = 1) anc_up & desc_down 2) if not best (anc_up, desc1_down, desc2_down) */
		a_final[ntax+n]=make_align(temp4);
		convert_to_final(a_final[ntax+n],values);
		/*printf("Optimizing node %d\n",n);
		for (i=0;i<4;i++) {
			for (j=0;j<temp4->length;j++) if (temp4->type_weight) {
				if ((temp4->s[temp4->n_seqs][j] != 1) && (temp4->s[temp4->n_seqs][j] != 65)) printf("%2d ",temp4->s[i][j]);
				else  printf("%c",temp4->s[i][j]);
				}
			printf("\n");
			}
		printf("Final node %d\n",n);
		for (j=0;j<a_final[ntax+n]->length;j++) if (a_final[ntax+n]->type_weight) {
			if ((a_final[ntax+n]->s[a_final[ntax+n]->n_seqs][j] != 1) && (a_final[ntax+n]->s[a_final[ntax+n]->n_seqs][j] != 65)) printf("%2d ",a_final[ntax+n]->s[0][j]);
			else  printf("%c",a_final[ntax+n]->s[0][j]);
			}
		printf("\n");
		*/
		}
	}
free(done);
if (temp1) temp1=dump_align(temp1);
if (temp2) temp2=dump_align(temp2);
if (temp3) temp3=dump_align(temp3);
if (temp4) temp4=dump_align(temp4);
}


void convert_to_anc_and_desc(thang,anc,desc)
alignment *thang;
int anc,desc;
{
char *s0, *s1, *s2;
int i, old_seqs;

s0=(char *)malloc((1+thang->length)*sizeof(char));
assert((int) s0);
s0=(char *)strcpy(s0,thang->s[anc]);

s1=(char *)malloc((1+thang->length)*sizeof(char));
assert((int) s1);
s1=(char *)strcpy(s1,thang->s[desc]);
/*fprintf(stderr,"new\nanc\n%s desc \n%s ",s0,s1);*/

if (thang->type_weight) {
	s2=(char *)malloc((1+thang->length)*sizeof(char));
	assert((int) s2);
	for (i=0;i<thang->length;i++) s2[i]=thang->s[4][i];
	s2[thang->length]='\0';
	}
else s2=NULL;
/*reallocate type weight too*/
thang->n_seqs=2;
for (i=0;i<(4+thang->type_weight);i++) free(thang->s[i]);
free(thang->s);
for (i=0;i<4;i++) free(thang->taxon_name[i]);
thang->taxon_name=(char **)realloc(thang->taxon_name,2*sizeof(char *));
assert((int) thang->taxon_name);
thang->s=(char **)malloc((thang->n_seqs+thang->type_weight)*sizeof(char *));
assert((int) thang->s);
for (i=0;i<(thang->n_seqs+thang->type_weight);i++) {
	thang->s[i]=(char *)malloc((1+thang->length)*sizeof(char));
	assert((int) thang->s[i]);
	}
thang->s[0]=(char *)strcpy(thang->s[0],s0);
thang->s[1]=(char *)strcpy(thang->s[1],s1);
if (thang->type_weight) thang->s[2]=(char *)strcpy(thang->s[2],s2);

thang->taxon_name[0]=(char *)malloc(2*sizeof(char));
assert((int) thang->taxon_name[0]);
thang->taxon_name[0][0]='F';
thang->taxon_name[0][1]='\0';
thang->taxon_name[1]=(char *)malloc(2*sizeof(char));
assert((int) thang->taxon_name[1]);
thang->taxon_name[1][0]='F';
thang->taxon_name[1][1]='\0';
/*free*/
free(s0);
free(s1);
if (s2) free(s2);
/*fprintf(stderr,"After Anc \n%s desc \n%s ",thang->s[0],thang->s[1]);*/
}

int all_diagnose_tree_local_new_opt2(a_down,nodes,ntax,nent,values,stuff_done,corres,desc)
alignment **a_down;
int **nodes, ntax,nent,stuff_done,***corres, *desc;
parameters *values;
{
int n,i,found_one;
alignment *temp_align;
int d1,d2,d1_thang,d2_thang;
int nent_minus_one_new,od1,od2;
int *done;


temp_align=NULL;
done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;
for (i=0;i<ntax+nent-1;i++) if (desc[i]) done[i]=1;
/*down pass redone*/
values->in_optimize=1;
found_one=1;
while (found_one) {
	found_one=0;
	for (n=(nent-2);n>0;--n){
		od1=d1=nodes[n][0];
		od2=d2=nodes[n][1];
		if ((done[d1]) && (done[d2]) && (!done[ntax+n])) {
			if (d1 < ntax) d1-=1;
			if (d2 < ntax) d2 -=1;
			found_one=1;
			done[ntax+n]=1;
			if (a_down[ntax+n]) a_down[ntax+n]=dump_align(a_down[ntax+n]);
			/*get alignment*/
			a_down[ntax+n] = nw(a_down[d1],a_down[d2],values);
			if (temp_align) temp_align=dump_align(temp_align);
			temp_align=make_align(a_down[ntax+n]);
			/*get ancestor with gaps*/
			temp_align=new_make_ambig(temp_align,values);
			/*get correspondance arrays*/
			if (od1!=od2) {
				if (strcmp(a_down[d1]->name,a_down[d2]->name)>0) {
					d2_thang=0;
					d1_thang=1;
					}
				else {
					d2_thang=1;
					d1_thang=0;
					}
				if (corres[od1]) {
					free(corres[od1][0]);
					free(corres[od1][1]);
					free(corres[od1]);
					}
				if (corres[od2]) {
					free(corres[od2][0]);
					free(corres[od2][1]);
					free(corres[od2]);
					}
				corres[od1]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[d1_thang],temp_align->length);
				corres[od2]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[d2_thang],temp_align->length);
				if (temp_align->type_weight) {
					modify_corres(corres[od1],temp_align);
					modify_corres(corres[od2],temp_align);
					}
				}
			else {
				if (corres[od1]) {
					free(corres[od1][0]);
					free(corres[od1][1]);
					free(corres[od1]);
					}
				corres[od1]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[0],temp_align->length);
				if (temp_align->type_weight) modify_corres(corres[od1],temp_align);
				}
			/*get ancestor without gaps*/
		    	a_down[ntax+n]=make_ambig(a_down[ntax+n],values);
		  	a_down[ntax+n]->score += (a_down[d1]->score+a_down[d2]->score);
		  	if (temp_align) temp_align=dump_align(temp_align);
	    		}
		}
	}
values->in_optimize=0;
if (temp_align) temp_align=dump_align(temp_align);
free(done);
return a_down[nodes[0][1]]->score;
}

int all_diagnose_tree_local_new_opt3(a_down,nodes,ntax,nent,values,stuff_done,corres,desc)
alignment **a_down;
int **nodes, ntax,nent,stuff_done,***corres, *desc;
parameters *values;
{
int n,i,found_one;
alignment *temp_align;
int d1,d2,d1_thang,d2_thang;
int nent_minus_one_new,od1,od2;
int *done;


temp_align=NULL;
done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;
for (i=0;i<ntax+nent-1;i++) if ((!done[i]) && (!desc[i])) done[i]=1;

/*down pass redone*/
values->in_optimize=1;
found_one=1;
while (found_one) {
	found_one=0;
	for (n=(nent-2);n>0;--n){
		od1=d1=nodes[n][0];
		od2=d2=nodes[n][1];
		if ((done[d1]) && (done[d2]) && (!done[ntax+n])) {
			if (d1 < ntax) d1-=1;
			if (d2 < ntax) d2 -=1;
			found_one=1;
			done[ntax+n]=1;
			if (a_down[ntax+n]) a_down[ntax+n]=dump_align(a_down[ntax+n]);
			/*get alignment*/
			a_down[ntax+n] = nw(a_down[d1],a_down[d2],values);
			if (temp_align) temp_align=dump_align(temp_align);
			temp_align=make_align(a_down[ntax+n]);
			/*get ancestor with gaps*/
			temp_align=new_make_ambig(temp_align,values);
			/*get correspondance arrays*/
			if (od1!=od2) {
				if (strcmp(a_down[d1]->name,a_down[d2]->name)>0) {
					d2_thang=0;
					d1_thang=1;
					}
				else {
					d2_thang=1;
					d1_thang=0;
					}
				if (corres[od1]) {
					free(corres[od1][0]);
					free(corres[od1][1]);
					free(corres[od1]);
					}
				if (corres[od2]) {
					free(corres[od2][0]);
					free(corres[od2][1]);
					free(corres[od2]);
					}
				corres[od1]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[d1_thang],temp_align->length);
				corres[od2]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[d2_thang],temp_align->length);
				if (temp_align->type_weight) {
					modify_corres(corres[od1],temp_align);
					modify_corres(corres[od2],temp_align);
					}
				}
			else {
				if (corres[od1]) {
					free(corres[od1][0]);
					free(corres[od1][1]);
					free(corres[od1]);
					}
				corres[od1]=get_correspondances(temp_align->s[0],a_down[ntax+n]->s[0],temp_align->length);
				if (temp_align->type_weight) modify_corres(corres[od1],temp_align);
				}
			/*get ancestor without gaps*/
		    	a_down[ntax+n]=make_ambig(a_down[ntax+n],values);
		  	a_down[ntax+n]->score += (a_down[d1]->score+a_down[d2]->score);
		  	if (temp_align) temp_align=dump_align(temp_align);
	    		}
		}
	}
values->in_optimize=0;
if (temp_align) temp_align=dump_align(temp_align);
free(done);
return a_down[nodes[0][1]]->score;
}



/*void convert_to_final(a)
alignment *a;
{
char *t,*w;
int i,j,new_length,old_length;;

t=(char *)malloc((1+a->length)*sizeof(char));
assert((int) t);
if (a->type_weight) {
	w=(char *)malloc((1+a->length)*sizeof(char));
	  assert((int) w);
	for (i=0;i<a->length;i++) w[i]=a->s[a->n_seqs][i];
	w[a->length]='\0';
	}

for (i=0;i<a->n_seqs;i++) free(a->taxon_name[i]);
free(a->taxon_name);
a->taxon_name=(char **)malloc(sizeof(char *));
assert((int) a->taxon_name);
a->taxon_name[0]=(char *)malloc(7*sizeof(char));
assert((int)a->taxon_name[0]);
a->taxon_name[0][0]='h';
a->taxon_name[0][1]='y';
a->taxon_name[0][2]='p';
a->taxon_name[0][3]='a';
a->taxon_name[0][4]='n';
a->taxon_name[0][5]='c';
a->taxon_name[0][5]='\0';
new_length=0;
for (i=0;i<a->length;i++) {
	if ((a->type_weight) && (w[i]!=1) && (w[i]!=65)) {
		t[i]=(a->s[0][i] & a->s[1][i]);
		if (!t[i]) t[i]=(char)get_best_states_morph(a->s[1][i],a->s[2][i],a->s[3][i]);
		++new_length;
		}
	else {
		t[i]=overlap_best(a->s[0][i],a->s[1][i],a->s[2][i],a->s[3][i]);
		if (t[i]!='-') ++new_length;
		}
	}
t[a->length]='\0';
old_length=a->length;
a->length=new_length;
for (i=0;i<a->n_seqs+a->type_weight;i++) free(a->s[i]);
free(a->s);
a->s=(char **)malloc((1+a->type_weight)*sizeof(char *));
assert((int) a->s);
for (i=0;i<1+a->type_weight;i++) {
	a->s[i]=(char *)malloc((1+a->length)*sizeof(char));
	assert((int) a->s[i]);
	}
j=0;
for (i=0;i<old_length;i++) {
	if ((a->type_weight) && (w[i]!=1) && (w[i]!=65)) {
		a->s[0][j]=t[i];
		a->s[1][j]=w[i];
		++j;
		}
	else {
		if (t[i]!='-') {
			a->s[0][j]=t[i];
			if (a->type_weight) a->s[1][j]=w[i];
			j++;
			}
		else if (a->type_weight) {
			if (w[i]>=64) a->s[1][j-1]=w[i];
			}
		}
	}
a->s[0][a->length]='\0';
if (a->type_weight) a->s[1][a->length]='\0';
a->n_seqs=1;
free(t);
if (a->type_weight) free(w);
}
*/


int get_best_states_morph(a,b,c,d)
char a,b,c,d;
{
/*c=(a&b);
if (!c) c= (a|b);*/

if (((c & d) == d) && (c != d)) return d;
else if (!(a & b)) return ( c | d);
else  return (c | ((a | b) & d));
}

void convert_to_final(a,values)
alignment *a;
parameters *values;
{
char *t,*w;
int i,j,k,new_length,old_length;
int in_c, in_b, in_a, in_d, out,best;
int c16,c8,c4,c2,c1;
int allgaps;

if (a->n_seqs!=4) {
	fprintf(stderr,"Should be four sequences in 'a' but are %d\n",a->n_seqs);
	exit(-1);
	}
t=(char *)malloc((1+a->length)*sizeof(char));
assert((int) t);
if (a->type_weight) {
	w=(char *)malloc((1+a->length)*sizeof(char));
	  assert((int) w);
	for (i=0;i<a->length;i++) w[i]=a->s[4][i];
	w[a->length]='\0';
	}

for (i=0;i<a->n_seqs;i++) free(a->taxon_name[i]);
free(a->taxon_name);
a->taxon_name=(char **)malloc(sizeof(char *));
assert((int) a->taxon_name);
a->taxon_name[0]=(char *)malloc(7*sizeof(char));
assert((int)a->taxon_name[0]);
a->taxon_name[0][0]='h';
a->taxon_name[0][1]='y';
a->taxon_name[0][2]='p';
a->taxon_name[0][3]='a';
a->taxon_name[0][4]='n';
a->taxon_name[0][5]='c';
a->taxon_name[0][5]='\0';
new_length=0;
for (i=0;i<a->length;i++) {
	if ((a->type_weight) && (w[i]!=1) && (w[i]!=65)) {
		/*t[i]=get_best_states_morph(a->s[0][i],a->s[1][i],a->s[2][i],a->s[3][i]);*/
		if (((a->s[2][i] & a->s[3][i]) == a->s[3][i]) && (a->s[2][i] != a->s[3][i])) t[i]= a->s[3][i];
        else if (!(a->s[0][i] & a->s[1][i])) t[i]=( a->s[2][i] | a->s[3][i]);
        else  t[i]= (a->s[2][i] | ((a->s[0][i] | a->s[1][i]) & a->s[3][i]));
		++new_length;
		}
	else {
		/*t[i]=new_overlap_best(a->s[0][i],a->s[1][i],a->s[2][i],a->s[3][i],values);*/
        in_a=charbit[(int)a->s[0][i]];
        in_b=charbit[(int)a->s[1][i]];
        in_c=charbit[(int)a->s[2][i]];
        in_d=charbit[(int)a->s[3][i]];
        /*if (!values->delta && !values->ttr && (values->gap_cost==values->change_cost)) {
            if (((in_c & in_d) == in_d) && (in_c != in_d)) out= in_d;
            else if (!(in_a & in_b)) out= ( in_c | in_d);
            else  out= (in_c | ((in_a | in_b) & in_d));
            }
        else {*/
            out= (in_a & in_b & in_d);
            if (!out) {
                c1=c2=c4=c8=c16=0;
                c1=(values->lookup[A_BIT][in_a].cost+values->lookup[A_BIT][in_b].cost+values->lookup[A_BIT][in_d].cost);
                c2=(values->lookup[C_BIT][in_a].cost+values->lookup[C_BIT][in_b].cost+values->lookup[C_BIT][in_d].cost);
                c4=(values->lookup[G_BIT][in_a].cost+values->lookup[G_BIT][in_b].cost+values->lookup[G_BIT][in_d].cost);
                c8=(values->lookup[T_BIT][in_a].cost+values->lookup[T_BIT][in_b].cost+values->lookup[T_BIT][in_d].cost);
                c16=(values->lookup[GAP_BIT][in_a].cost+values->lookup[GAP_BIT][in_b].cost+values->lookup[GAP_BIT][in_d].cost);
                best=min(c1,min(c2,min(c4,min(c8,c16))));
                if (c1==best) out+=A_BIT;
                if (c2==best) out+=C_BIT;
                if (c4==best) out+=G_BIT;
                if (c8==best) out+=T_BIT;
                if (c16==best) out+=GAP_BIT;
                }
            /*}*/
        t[i]=bitchar[out];
		if (t[i]!='-') ++new_length;
		}
	}
t[a->length]='\0';
old_length=a->length;
/*prefilter for zero length pieces*/
if (a->type_weight) {
    if (t[0]=='-') allgaps=1;
    else allgaps=0;
    for (i=1;i<old_length;i++) {
        if (w[i]==67) allgaps=1;
        else if (w[i]>63){
            if (allgaps) {
                if (t[i]=='-') {
                    t[i]='X';
                    ++new_length;
                }
            }
            else allgaps=1;
        }
        else if ((w[i]==3) || (t[i]!='-')) allgaps=0;
    }
}
a->length=new_length;
for (i=0;i<a->n_seqs+a->type_weight;i++) free(a->s[i]);
free(a->s);
a->s=(char **)malloc((1+a->type_weight)*sizeof(char *));
assert((int) a->s);
for (i=0;i<1+a->type_weight;i++) {
	a->s[i]=(char *)malloc((1+a->length)*sizeof(char));
	assert((int) a->s[i]);
	}
j=0;

for (i=0;i<old_length;i++) {
	if ((a->type_weight) && (w[i]!=1) && (w[i]!=65)) {
		a->s[0][j]=t[i];
		a->s[1][j]=w[i];
		++j;
		}
	else {
		if (t[i]!='-') {
			a->s[0][j]=t[i];
			if (a->type_weight) a->s[1][j]=w[i];
			j++;
			}
		else {
		    if (a->type_weight && (w[i]>=64)) a->s[1][j-1]=w[i];
		    }
		}
	}
a->s[0][a->length]='\0';
if (a->type_weight) a->s[1][a->length]='\0';
/*if segments have zero length give them an 'X'
if (a->type_weight) {
    for (i=1;i<a->length;i++) {
        if ((a->s[1][i-1]>63) && (a->s[1][i]>63)) {
            ++a->length;
            a->s[0]=(char *)realloc(a->s[0],(1+a->length)*sizeof(char));
            assert((int) a->s[0]);
            a->s[0][a->length]='\0';
            a->s[1]=(char *)realloc(a->s[1],(1+a->length)*sizeof(char));
            assert((int) a->s[1]);
            a->s[1][a->length]='\0';
            for (j=a->length-1;j>i;j--) {
                 a->s[0][j]=a->s[0][j-1];
                 a->s[1][j]=a->s[1][j-1];
                }
            if (a->s[1][i]==65) a->s[0][i]='X';
            else a->s[0][i]=127;
            a->s[1][i]-=64;
            }
        }
    }*/

a->n_seqs=1;
free(t);
if (a->type_weight) free(w);

}


void get_up_pass_new_opt(a_final,a_down,nodes,ntax,values,nent,anc,a_down2)
alignment **a_down, **a_final,**a_down2;
int **nodes, ntax,nent, *anc;
parameters *values;
{
int i,d1,d2,j,found_one,n;
int node_is_from;
int *done;

done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;


if (a_final[nodes[0][1]]) dump_align(a_final[nodes[0][1]]);
a_final[nodes[0][1]]=make_align(a_down[nodes[0][1]]);

done[nodes[0][1]]=1;

found_one=1;
while (found_one) {
	found_one=0;
	for (n=1;n<nent-1;n++) if ((!done[n+ntax]) && (done[anc[n+ntax]])) {
		found_one=1;
		done[n+ntax]=1;
		if (a_final[ntax+n]) a_final[ntax+n]=dump_align(a_final[ntax+n]);
		values->in_optimize=1;
		values->new_optimization=0;
		/*getting order predictable in 4 align*/
		/*Make order s0=D1down s1=D2down s2=NodeDown s3=AncFinal*/
		free(a_down2[ntax+n]->name);
		a_down2[ntax+n]->name=(char *)malloc(2*sizeof(char));
		assert((int)a_down2[ntax+n]->name);
		a_down2[ntax+n]->name[0]='a'; 
	        a_down2[ntax+n]->name[1]='\0';
		a_final[anc[ntax+n]]->name=(char *)malloc(2*sizeof(char));
		assert((int) a_final[anc[ntax+n]]->name);
		a_final[anc[ntax+n]]->name[0]='b'; a_final[anc[ntax+n]]->name[1]='\0';;
		a_final[ntax+n]=nw(a_down2[ntax+n],a_final[anc[ntax+n]],values);
		values->new_optimization=1;
		values->in_optimize=0;
		convert_to_final(a_final[ntax+n],values);
		}
	}
free(done);
}

