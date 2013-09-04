/*Copyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Goloboff*/

#include "align3.h"

int is_unique_tree_align(best_aligns,cur_nodes,new_nodes,num_new,old_nodes,num_old,old_length,new_length,ntax,a,do_in_parallel,values,collapse_name)
int **cur_nodes,do_in_parallel;
int ***new_nodes, ***old_nodes;
int num_new, num_old, new_length, old_length, ntax;
alignment **a, **best_aligns;
parameters *values;
char *collapse_name;
{
int i,j,k,test;
int **cur_groups, **cur_groups2;
int *anc,check;
alignment *temp_align,**a_buf,**a_down2;
char *temp_name;

temp_align=NULL;

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


a_down2=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_down2);
for (i=0;i<2*ntax-1;i++) a_down2[i]=NULL;
a_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_buf);
for (i=0;i<2*ntax-1;i++) a_buf[i]=NULL;
for (i=0;i<ntax-1;i++) a_buf[i]=make_align(a[i]);
anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;
if (!do_in_parallel) {
	check=all_diagnose_tree_here2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
 	get_up_pass_new_opt(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
else {
	check=all_diagnose_tree_here_parallel2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
	get_up_pass_new_opt_parallel(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}

/*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=0;i<ntax-1;i++) values->support[num_new][i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
					if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {values->support[num_new][i-1]=1;break;}
					} /* k */
					
				} /* collapse ==1 */
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) values->support[num_new][i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /*else*/
			} /* else */
		}/* i */

/* get new nAME*/
 temp_name=(char *)get_names_with_collapse(a,cur_nodes,ntax,ntax,values,values->support[num_new]);
 collapse_name=(char *)strcpy(collapse_name,temp_name);
 free(temp_name);

for (i=0;i<num_new;i++) {
	get_cur_groups_tree(cur_groups2,new_nodes[i],ntax,ntax,ntax-2);
	test=compare_groups_tree_align(cur_groups,cur_groups2,ntax-2,ntax-2,ntax-1,ntax-1,values->support[num_new],values->support[i]);
	if (test) {
		for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
		free(cur_groups2);
		for (j=0;j<ntax-2;j++) free(cur_groups[j]);
		free(cur_groups);
		free(anc);
		for (j=0;j<2*ntax-1;j++) if (a_buf[j]) a_buf[j]=dump_align(a_buf[j]);
		free(a_buf);
		return 0;
		}
	}


for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
free(cur_groups2);
for (j=0;j<ntax-2;j++) free(cur_groups[j]);
free(cur_groups);
free(anc);

for (j=0;j<2*ntax-1;j++) {
	if (a_buf[j]) a_buf[j]=dump_align(a_buf[j]);
	if (a_down2[j]) a_down2[j]=dump_align(a_down2[j]);
	}
free(a_down2);
free(a_buf);

return 1;
}

int compare_groups_tree_align(g1,g2,ngroups1,ngroups2,ntax,nent,g1_sup,g2_sup)
int **g1,**g2,*g1_sup,*g2_sup;
int ngroups1,ngroups2,ntax,nent;
{
int i,j,k,holder,temp,disjunct,are_both_one, are_not_same;
int g1_tot, g2_tot;

for (k=0;k<ngroups1;k++) if (g1_sup[k]) {
	for (i=0;i<ngroups2;i++) if (g2_sup[i]) {
		/*for (j=0;j<ntax;j++) fprintf(stderr,"%d:%d %d ",j,g1[k][j],g2[i][j]);*/
		/*check make sure each group is consisten with others*/
		/*are they completely disjunct*/
		are_both_one=are_not_same=g1_tot=g2_tot=0;
		for (j=0;j<ntax-1;j++) {
			/*fprintf(stderr,"{%d+%d,%d=%d%d}",k,i,j,g1[k][j],g2[i][j]);*/
			if ((g1[k][j]+g2[i][j])==2) are_both_one+=1;
			else if ((g1[k][j]+g2[i][j])==1) are_not_same=1;
			g1_tot+=g1[k][j];
			g2_tot+=g2[i][j];
			}
		if (are_not_same && are_both_one) {
			if ((are_both_one != g1_tot) && (are_both_one != g2_tot)) {
				/*fprintf(stderr,"\nans %d abo %d g1 %d g2 %d ",are_not_same,are_both_one,g1_tot,g2_tot);
				for (j=0;j<ntax-1;j++) fprintf(stderr,"{%d+%d,%d=%d%d}",k,i,j,g1[k][j],g2[i][j]);*/
				return 0;
				}
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


char *get_names_with_collapse(a,nodes,ntaxa,nent,values,support)
int **nodes,*support;
int ntaxa,nent;
alignment **a;
parameters *values;
{
char **names,*pre_made,*found;
int i,n,d1,d2,j,those_to_check;
char *new_name;
int found_one;

/*allocate names*/
names=(char **)malloc((ntaxa+ntaxa-1)*sizeof(char *));
assert((int)names);
pre_made=(char *)malloc((ntaxa+ntaxa-1)*sizeof(char));
assert((int)pre_made);

/*optimize names*/
for (i=0;i<ntaxa-1;i++) {
	pre_made[i]=1;
  	names[i]=(char *)malloc((1+strlen(a[i]->name))*sizeof(char));
	assert((int)names[i]);
	names[i]=strcpy(names[i],a[i]->name);
	}
pre_made[ntaxa-1]=1;
names[ntaxa-1]=NULL;
for (i=ntaxa;i<(ntaxa+ntaxa-1);i++) {
  pre_made[i]=0;
  names[i]=NULL;
  }
found_one=1;
while(found_one) {
	found_one=0;
	for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings*/
		d1=nodes[n][0];
		d2=nodes[n][1];
		if ((pre_made[d1]) && (pre_made[d2]) && (!pre_made[ntaxa+n])) {
			pre_made[ntaxa+n]=1;
			/*fix for rooting more taxa*/
			if (d1<ntaxa) d1-=1;
			if (d2<ntaxa) d2-=1;
			found_one=1;
			/*make names*/
			if (support[n-1] || (n==(nodes[0][1]-ntaxa))) {
				if (strcmp (names[d1],names[d2]) > 0) names[ntaxa+n]=other_pair_names(names[d2],names[d1]);
				else names[ntaxa+n]=other_pair_names(names[d1],names[d2]);
				}
			else {
				if (strcmp (names[d1],names[d2]) > 0) names[ntaxa+n]=third_pair_names(names[d2],names[d1]);
				else names[ntaxa+n]=third_pair_names(names[d1],names[d2]);
				}

      }
    }
}

new_name=(char *)malloc((1+strlen(names[nodes[0][1]]))*sizeof(char));
assert((int)new_name);
new_name=(char *)strcpy(new_name,names[nodes[0][1]]);

free(pre_made);
for (i=0;i<(ntaxa+ntaxa-1);i++) if (names[i]) free(names[i]);
free(names);
return(new_name);
}

char *third_pair_names (a1, a2)
char *a1, *a2;
{
	char *p;

	p = (char *) allocate (strlen (a1) + strlen (a2) +3);
	sprintf (p, "%s %s", a1, a2);
	return p;
}


int get_collapsed_thang(best_align,cur_nodes,ntax,a,do_in_parallel,values)
int **cur_nodes,do_in_parallel;
int ntax;
alignment **a, *best_align;
parameters *values;
{
int i,j,k,test, *g1_sup,*anc,check;
alignment *temp_align,**a_buf,**a_down2;
char *temp_name;
int d1,d2;

temp_align=NULL;

/*allocate array for supported groups
0 if group ! supported 1 if is */
g1_sup=(int *)malloc(ntax*sizeof(int));
assert((int) g1_sup);

a_down2=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_down2);
for (i=0;i<2*ntax-1;i++) a_down2[i]=NULL;
a_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_buf);
for (i=0;i<2*ntax-1;i++) a_buf[i]=NULL;
for (i=0;i<ntax-1;i++) a_buf[i]=make_align(a[i]);

anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

/*get support for the new tree*/
anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;
if (!do_in_parallel) {
	check=all_diagnose_tree_here2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
 	get_up_pass_new_opt(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
else {
	check=all_diagnose_tree_here_parallel2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
	get_up_pass_new_opt_parallel(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
/*if NW->score for internode group !sup
the node number -1 is the group number so no need to check for basal support really*/
g1_sup[0]=1;
for (i=1;i<ntax-1;i++) g1_sup[i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
				if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {g1_sup[i-1]=1;break;}
					} /* k */
					
				} /* collapse ==1 */
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) g1_sup[i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /*else*/
			} /* else */
		}/* i */

/*get collapsed name of new thang*/
free(best_align->name);
best_align->name=(char *)get_names_with_collapse(a,cur_nodes,ntax,ntax,values,g1_sup);
/*check on nodes etc
for (i=1;i<ntax-1;i++)  {
    d1=cur_nodes[i][0];
    d2=cur_nodes[i][1];
    if (d1 < ntax) d1-=1;
    if (d2 < ntax) d2-=1;
    printf("Node %d\n",i);
    for (j=0;j<a[0]->length;j++) printf("%4d",a_down2[ntax+i]->s[0][j]);
    printf("\n");
    for (j=0;j<a[0]->length;j++) printf("%4d",a_down2[ntax+i]->s[1][j]);
    printf("\n");
    for (j=0;j<a[0]->length;j++) printf("%4d",a[ntax+i]->s[0][j]);
    printf("\n");
    if (anc[ntax+i]!=ntax) {
        for (j=0;j<a[0]->length;j++) printf("%4d",a_buf[anc[ntax+i]]->s[0][j]);
        printf("\n");
    }
    else printf("BASAL\n");
    for (j=0;j<a[0]->length;j++) printf("%4d",a_buf[i]->s[0][j]);
    printf("\n");

}*/


free(g1_sup);
free(anc);
for (j=0;j<2*ntax-1;j++) if (a_down2[j]) a_down2[j]=dump_align(a_down2[j]);
free(a_down2);

for (j=0;j<2*ntax-1;j++) if (a_buf[j]) a_buf[j]=dump_align(a_buf[j]);
free(a_buf);

return 1;
}

int get_support_and_name(best_align,anc,cur_nodes,a_buf,ntax,g1_sup,values)
int *anc,**cur_nodes;
int ntax,*g1_sup;
alignment *best_align,**a_buf;
parameters *values;
{
int i,k;
alignment *temp_align;

temp_align=NULL;

for (i=0;i<ntax-1;i++) g1_sup[i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
					if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {g1_sup[i-1]=1;break;}
					} /* k */
					
				} /* collapse ==1 */
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) g1_sup[i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /*else*/
			} /* else */
		}/* i */

/*get collapsed name of new thang*/
free(best_align->name);
best_align->name=(char *)get_names_with_collapse(a_buf,cur_nodes,ntax,ntax,values,g1_sup);
return 1;
}

int is_unique_tree_align_new(best_aligns,cur_nodes,new_nodes,num_new,old_nodes,num_old,old_length,new_length,ntax,a,do_in_parallel,values,collapse_name,corres,a_buf,a_up_buf,a_down2)
int **cur_nodes,do_in_parallel;
int ***new_nodes, ***old_nodes;
int num_new, num_old, new_length, old_length, ntax;
alignment **a, **best_aligns,**a_buf,***a_up_buf,**a_down2;
parameters *values;
int ***corres;
char *collapse_name;
{
int i,j,k,test;
int **cur_groups, **cur_groups2;
int *anc,check;
alignment *temp_align;
char *temp_name;
int g1_test;

temp_align=NULL;

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


anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

/*get support for the new tree*/
anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;
if (!do_in_parallel) get_up_pass_new_opt(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
else get_up_pass_new_opt_parallel(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
/*if NW->score for internode group !sup
the node number -1 is the group number so no need to check for basal support really*/

for (i=0;i<ntax-1;i++) values->support[num_new][i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
					if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {values->support[num_new][i-1]=1;break;}
					} /* k */
					
				} /* collapse ==1 */
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) values->support[num_new][i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /*else*/
			} /* else */
		}/* i */
/*get collapsed name of new thang*/
 temp_name=(char *)get_names_with_collapse(a,cur_nodes,ntax,ntax,values,values->support[num_new]);
 collapse_name=(char *)strcpy(collapse_name,temp_name);
 free(temp_name);
 /*fprintf(stderr,"CN %s ",collapse_name);*/

for (i=0;i<num_new;i++) {
	get_cur_groups_tree(cur_groups2,new_nodes[i],ntax,ntax,ntax-2);
	test=compare_groups_tree_align(cur_groups,cur_groups2,ntax-2,ntax-2,ntax-1,ntax-1,values->support[num_new],values->support[i]);
	if (test) {
		for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
		free(cur_groups2);
		for (j=0;j<ntax-2;j++) free(cur_groups[j]);
		free(cur_groups);
		free(anc);
		return 0;
		}
	}

/*is this part necessary??
if (new_length==old_length) for (i=0;i<num_old;i++) {
	get_cur_groups_tree(cur_groups2,old_nodes[i],ntax,ntax,ntax-2);
	test=compare_groups_tree_align(cur_groups,cur_groups2,ntax-2,ntax-2,ntax-1,ntax-1,g1_sup,g2_sup);
	if (test) {
		for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
		free(cur_groups2);
		for (j=0;j<ntax-2;j++) free(cur_groups[j]);
		free(cur_groups);
		free(anc);
		return 0;
		}
	}*/

for (j=0;j<ntax-2;j++) free(cur_groups2[j]);
free(cur_groups2);
for (j=0;j<ntax-2;j++) free(cur_groups[j]);
free(cur_groups);
free(anc);
return 1;
}

int get_collapsed_thang_with_groups(best_align,cur_nodes,ntax,a,do_in_parallel,values,g1_sup,a_buf,a_up_buf,corres)
int **cur_nodes,***corres,do_in_parallel;
int ntax;
int *g1_sup;
alignment **a, *best_align,**a_buf,***a_up_buf;
parameters *values;
{
int i,j,k,test;
int *anc,check;
alignment *temp_align,**a_down2;
char *temp_name;

temp_align=NULL;

a_down2=(alignment **)malloc(((2*ntax)-1)*sizeof(alignment *));
assert((int) a_down2);
for (i=0;i<(2*ntax)-1;i++) a_down2[i]=NULL;
/*allocate array for supported groups
0 if group ! supported 1 if is */



anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

/*get support for the new tree*/
anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;
if (!do_in_parallel) {
	check=all_diagnose_tree_here2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
 	get_up_pass_new_opt(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
else {
	check=all_diagnose_tree_here_parallel2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
	get_up_pass_new_opt_parallel(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
/*if NW->score for internode group !sup
the node number -1 is the group number so no need to check for basal support really*/
for (i=0;i<ntax-1;i++) g1_sup[i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
					if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {g1_sup[i-1]=1;break;}
					} /* k */
					
				} /* collapse ==1 */
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) g1_sup[i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /*else*/
			} /* else */
		}/* i */
/*get collapsed name of new thang*/
free(best_align->name);
best_align->name=(char *)get_names_with_collapse(a,cur_nodes,ntax,ntax,values,g1_sup);
for (i=0;i<(2*ntax)-1;i++) if (a_down2[i]) a_down2[i]=dump_align(a_down2[i]);
free(a_down2);
free(anc);
return 1;
}

int get_collapsed_thang_after_up(best_align,cur_nodes,ntax,a,do_in_parallel,values,a_buf)
int **cur_nodes,do_in_parallel;
int ntax;
alignment **a,**a_buf, *best_align;
parameters *values;
{
int i,j,k,test;
int *g1_sup;
int *anc,***corres,check;
alignment *temp_align;
char *temp_name;

temp_align=NULL;


/*allocate array for supported groups
0 if group ! supported 1 if is */
g1_sup=(int *)malloc(ntax*sizeof(int));
assert((int) g1_sup);

anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

/*get support for the new tree*/
anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;
/*if NW->score for internode group !sup
the node number -1 is the group number so no need to check for basal support really*/
for (i=0;i<ntax-1;i++) g1_sup[i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
					if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {g1_sup[i-1]=1;break;}
					} /* k */
					
				} /* collapse ==1 */
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) g1_sup[i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /*else*/
			} /* else */
		}/* i */
/*get collapsed name of new thang*/
free(best_align->name);
best_align->name=(char *)get_names_with_collapse(a,cur_nodes,ntax,ntax,values,g1_sup);

free(g1_sup);
free(anc);
return 1;
}

int get_collapsed_thang_with_groups_old(best_align,cur_nodes,ntax,a,do_in_parallel,values,g1_sup)
int **cur_nodes,do_in_parallel;
int ntax;
int *g1_sup;
alignment **a, *best_align;
parameters *values;
{
int i,j,k,test;
int *anc,check;
alignment *temp_align,**a_buf,**a_down2;
char *temp_name;

temp_align=NULL;


/*allocate array for supported groups
0 if group ! supported 1 if is */

a_down2=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_down2);
for (i=0;i<2*ntax-1;i++) a_down2[i]=NULL;

a_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_buf);
for (i=0;i<2*ntax-1;i++) a_buf[i]=NULL;
for (i=0;i<ntax-1;i++) a_buf[i]=make_align(a[i]);


anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

/*get support for the new tree*/
anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;
if (!do_in_parallel) {
	check=all_diagnose_tree_here2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
 	get_up_pass_new_opt(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
else {
	check=all_diagnose_tree_here_parallel2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
	get_up_pass_new_opt_parallel(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
	}
/*if NW->score for internode group !sup
the node number -1 is the group number so no need to check for basal support really*/

for (i=0;i<ntax-1;i++) g1_sup[i]=0; /*no need to check the basal most group (X+all) and sets the next up also (ALL) */
for (i=1;i<ntax-1;i++) if (anc[i+ntax]!=ntax) {
	if ((!only_missing(a_buf[i+ntax])) && (!only_missing (a_buf[anc[i+ntax]]))) {
			if (values->collapse==1) {
				if (a_buf[i+ntax]->length != a_buf[anc[i+ntax]]->length) g1_sup[i-1]=1;
				else {
					for (k=0;k<a_buf[i+ntax]->length;k++) if ((a_buf[i+ntax]->s[0][k]!='X') && (a_buf[anc[i+ntax]]->s[0][k]!='X')) {
						if (a_buf[i+ntax]->s[0][k]!=a_buf[anc[i+ntax]]->s[0][k]) {g1_sup[i-1]=1;break;}
						}
					}
				} /* collapse == 1*/
			else {
				if (temp_align) temp_align=dump_align(temp_align);
				temp_align=nw(a_buf[i+ntax],a_buf[anc[i+ntax]],values);
				if (temp_align->score!=0) g1_sup[i-1]=1;
				if (temp_align) temp_align=dump_align(temp_align);
				} /* else */
			} /* else */
		}/* i* /

/*get collapsed name of new thang*/
free(best_align->name);
best_align->name=(char *)get_names_with_collapse(a,cur_nodes,ntax,ntax,values,g1_sup);

for (j=0;j<2*ntax;j++) {
	for (j=0;j<2*ntax-1;j++) if (a_buf[j]) a_buf[j]=dump_align(a_buf[j]);
	free(a_buf);
	}
for (i=0;i<(2*ntax)-1;i++) if (a_down2[i]) a_down2[i]=dump_align(a_down2[i]);
free(a_down2);
free(anc);
return 1;
}

