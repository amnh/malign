/*Copyright 1994 Ward Wheeler all rights reserved*/
/*Make sure that the reconstructed alignments == the initial ones + can count gaps and changes from down pass as a check*/

#include "align3.h"
/*this assumes 2 sequences
    need to add stuff for type weight in make ambig to allow for differnt union procedures
    this could be direct from input characters to lookup?? would it be faster avoid conversion entirely for bases but not morph or whatever*/

alignment *fake_nw(desc,anc,values,corres) /*now does type_weight shit too*/
alignment *desc, *anc;
int **corres;
parameters *values;
{
alignment *new_desc;
int i,j,position,length,piece;
char *s0;
char *s1;
char *s2;

i=length=0;
while(corres[0][i]>-1) {i++,length++;}

s0=(char *)malloc((length+1)*sizeof(char *));
assert((int) s0);
s1=(char *)malloc((length+1)*sizeof(char *));
assert((int) s1);
s0[length]='\0';
s1[length]='\0';
if (anc->type_weight) {
	s2=(char *)malloc((length+1)*sizeof(char *));
	assert((int) s2);
	s2[length]='\0';
	piece=0;
	for (i=0;i<length;i++) {
		/*fprintf(stderr,"%d->",corres[0][i]);*/
		if (corres[0][i]<64) s2[i]=0;
		else s2[i]=64;
		s2[i]+=(values->data_sets[piece]+1);
		if (corres[0][i]>=64) ++piece;
		/*fprintf(stderr,"%d ",s2[i]);*/
		}
	}

position=0;
for (i=0;i<length;i++) {
	if ((!corres[0][i]) || (corres[0][i]==64)) s0[i]=desc->s[0][position++];
	else s0[i]='-';
	}
position=0;
for (i=0;i<length;i++) {
	if ((!corres[1][i]) || (corres[1][i]==64)) s1[i]=anc->s[0][position++];
	else s1[i]='-';
	}
s0[length]=s1[length]='\0';

/*ruthless fix*/
if (strlen(s0)!=strlen(s1)) fprintf(stderr,"Problem in fake_nw--fixing %d %d %d (\n%s\n%s)\n",length,strlen(s0),strlen(s1),s0,s1);
if (strlen(s0)!=strlen(s1)) {
	i=0;
	while (corres[0][i]>-1) fprintf(stderr,"%d",corres[0][i++]);
	fprintf(stderr,"\n");
	i=0;
	while (corres[1][i]>-1) fprintf(stderr,"%d",corres[1][i++]);
	fprintf(stderr,"\n");
	}

/*terminal gaps in ancestor*/
	position=length;
	new_desc=(alignment *)malloc(sizeof(alignment));
	assert ((int) new_desc);
	new_desc->n_seqs=2;
	new_desc->type_weight=desc->type_weight;
	new_desc->length=position;
	new_desc->score=desc->score; /*score_test; b->score;*/
	new_desc->name=(char *)malloc((strlen(desc->name)+1)*sizeof(char));
	assert((int) new_desc->name);
	new_desc->name=(char *)strcpy(new_desc->name,desc->name);
	new_desc->s=(char **)malloc((2+new_desc->type_weight)*sizeof(char *));
	assert((int) new_desc);
	for (i=0;i<(new_desc->n_seqs+new_desc->type_weight);i++) {
		new_desc->s[i]=(char *)malloc((position+1)*sizeof(char));
		assert ((int) new_desc->s[i]);
		new_desc->s[i][position]='\0';
		}
	new_desc->s[0]=strcpy(new_desc->s[0],s0);
	new_desc->s[1]=strcpy(new_desc->s[1],s1);
	if (new_desc->type_weight) new_desc->s[2]=strcpy(new_desc->s[2],s2);
	new_desc->taxon_name=(char **)malloc(2*sizeof(char *));
	assert((int) new_desc->taxon_name);
	new_desc->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) new_desc->taxon_name[0]);
	new_desc->taxon_name[0][0]='H';
	new_desc->taxon_name[0][1]='y';
	new_desc->taxon_name[0][2]='p';
	new_desc->taxon_name[0][3]='A';
	new_desc->taxon_name[0][4]='n';
	new_desc->taxon_name[0][5]='c';
	new_desc->taxon_name[0][6]='\0';
	new_desc->taxon_name[1]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) new_desc->taxon_name[1]);
	new_desc->taxon_name[1][0]='H';
	new_desc->taxon_name[1][1]='y';
	new_desc->taxon_name[1][2]='p';
	new_desc->taxon_name[1][3]='A';
	new_desc->taxon_name[1][4]='n';
	new_desc->taxon_name[1][5]='c';
	new_desc->taxon_name[1][6]='\0';
free(s0);
free(s1);
if (anc->type_weight) free(s2);

/*printf("NW\n");
printf("%s\n%s\n",new_desc->s[0],new_desc->s[1]);*/
return new_desc;
}

/* Need to take two descendent alignments which are basically equall and reconcile them with ancestor alignment
in each case 0 is desc 1 is ancestor
a1 is current node and ancestor
a2 is desc1 and cur
a3 is desc2 and cur
*/
alignment *make_neighbor(a1,a2,a3,values)
parameters *values;
alignment *a1, *a2, *a3;
{
alignment *four;
int **gap_holder,end0,end1,new_length,new_length2;
int i,j,k;
char *tw0, *tw1;

/*this holds the extra gaps the pieces need*/
gap_holder=(int **)malloc(2*sizeof(int *));
assert((int) gap_holder);
for (i=0;i<2;i++) {
	gap_holder[i]=(int *)malloc((a1->length+a2->length)*sizeof(int));
	assert((int) gap_holder[i]);
	}
for (i=0;i<(a1->length+a2->length);i++) gap_holder[0][i]=gap_holder[1][i]=0;
j=k=0;
for (i=0;i<(a1->length+a2->length);i++) if ((j<a1->length) && (k<a2->length)) {
	if (a1->s[0][j]!=a2->s[1][k]) {
		if (a1->s[0][j]=='-') {
			++gap_holder[1][k];
			j++;
			}
		else {
			++gap_holder[0][j];
			k++;
			}
		}
	else {
		k++;
		j++;
		}
	}
/*terminal gaps*/
end0=end1=0;
end0=a1->length-j;
end1=a2->length-k;
new_length=a1->length;
new_length2=a2->length;
for (i=0;i<a1->length;i++) new_length+=gap_holder[0][i];
new_length+=end1;
for (i=0;i<a2->length;i++) new_length2+=gap_holder[1][i];
new_length2+=end0;
if (new_length !=new_length2) {
	fprintf(stderr,"Error in make_four (%d vs. %d)\n",new_length,new_length2);
	exit(-1);
	}
/*make the four alignment*/
four=(alignment *)malloc(sizeof(alignment));
assert((int) four);
four->type_weight=a1->type_weight;
four->n_seqs=4;
four->length=new_length;
four->score=0;
four->name=(char *)malloc(2*sizeof(char));
assert((int) four->name);
four->name[0]='F';
four->name[1]='\0';
four->taxon_name=(char **)malloc(four->n_seqs*sizeof(char *));
assert((int) four->taxon_name);
for (i=0;i<four->n_seqs;i++) {
	four->taxon_name[i]=(char *)malloc(2*sizeof(char));
	assert((int) four->taxon_name[i]);
	four->taxon_name[i][0]='F';
	four->taxon_name[i][1]='\0';
	}
four->s=(char **)malloc((four->type_weight+four->n_seqs)*sizeof(char *));
assert((int) four->s);
for (i=0;i<(four->type_weight+four->n_seqs);i++) {
	four->s[i]=(char *)malloc((1+new_length)*sizeof(char));
	assert((int) four->s[i]);
	four->s[i][new_length]='\0';
	}
/*to make s's
	0 = current node a1->s[0],a2->s[1],a3->s[1];
	1 = its ancestor a1->s[1]
	2 = desc0		  a2->s[0]
	3 = desc1        a3->s[0]
*/
/*current node*/
j=0;
for (i=0;i<a1->length;i++) {
	if (gap_holder[0][i]==0)	four->s[0][j++]=a1->s[0][i];
	else {
		for (k=j;k<(j+gap_holder[0][i]);k++) four->s[0][k]='-';
		four->s[0][k]=a1->s[0][i];
		j=(k+1);
		}
	}
while (j<new_length) four->s[0][j++]='-';

/*ancestral node*/
j=0;
for (i=0;i<a1->length;i++) {
	if (gap_holder[0][i]==0)	four->s[1][j++]=a1->s[1][i];
	else {
		for (k=j;k<(j+gap_holder[0][i]);k++) four->s[1][k]='-';
		four->s[1][k]=a1->s[1][i];
		j=(k+1);
		}
	}
while (j<new_length) four->s[1][j++]='-';

/*desc0 node*/
j=0;
for (i=0;i<a2->length;i++) {
	if (gap_holder[1][i]==0)	four->s[2][j++]=a2->s[0][i];
	else {
		for (k=j;k<(j+gap_holder[1][i]);k++) four->s[2][k]='-';
		four->s[2][k]=a2->s[0][i];
		j=(k+1);
		}
	}
while (j<new_length) four->s[2][j++]='-';

/*desc1 node*/
j=0;
for (i=0;i<a3->length;i++) {
	if (gap_holder[1][i]==0)	four->s[3][j++]=a3->s[0][i];
	else {
		for (k=j;k<(j+gap_holder[1][i]);k++) four->s[3][k]='-';
		four->s[3][k]=a3->s[0][i];
		j=(k+1);
		}
	}
while (j<new_length) four->s[3][j++]='-';

/*type weight stuff*/
if (four->type_weight) {
	tw0=(char *)malloc((1+four->length)*sizeof(char));
	assert((int) tw0);
	tw1=(char *)malloc((1+four->length)*sizeof(char));
	assert((int) tw1);
	i=0;
	for (j=0;j<four->length;j++) {
		if (a1->s[0][i]==four->s[0][j]) {
			tw0[j]=a1->s[2][i];
			++i;
			}
		else tw0[j]=(100);
		}
	i=0;
	for (j=0;j<four->length;j++) {
		if (a2->s[0][i]==four->s[0][j]) {
			tw1[j]=a2->s[2][i];
			++i;
			}
		else tw1[j]=(100);
		}
	for (i=0;i<four->length;i++) {
		if (tw0[i]==tw1[i]) four->s[4][i]=tw0[i];
		else if ((tw0[i]== (100)) && (tw1[i]!= (100))) four->s[4][i]=tw1[i];
		else if ((tw1[i]== (100)) && (tw0[i]!= (100))) four->s[4][i]=tw0[i];
		else if (tw0[i]==(tw1[i]+64)) four->s[4][i]=tw0[i];
		else if (tw1[i]==(tw0[i]+64)) four->s[4][i]=tw1[i];
		}
	for (i=0;i<four->length-1;i++) if ((four->s[4][i] >= 64) && (four->s[4][i+1] >= 64)) four->s[4][i]-=64;
	free(tw0);
	free(tw1);
	}

free(gap_holder[0]);
free(gap_holder[1]);
free(gap_holder);
return four;
}

alignment *make_neighbor_non_align(a1,a2,a3,values)
parameters *values;
alignment *a1, *a2, *a3;
{
alignment *four;
int i,j,k;

if (a1->length!=a2->length) {fprintf(stderr,"Cannot neighbor sequences of unequal length\n"); exit(-1);}
if (a1->length!=a3->length) {fprintf(stderr,"Cannot neighbor sequences of unequal length\n"); exit(-1);}
/*this holds the extra gaps the pieces need*/

/*make the four alignment*/
four=(alignment *)malloc(sizeof(alignment));
assert((int) four);
four->type_weight=a1->type_weight;
four->n_seqs=4;
four->length=a1->length;
four->score=0;
four->name=(char *)malloc(2*sizeof(char));
assert((int) four->name);
four->name[0]='F';
four->name[1]='\0';
four->taxon_name=(char **)malloc(four->n_seqs*sizeof(char *));
assert((int) four->taxon_name);
for (i=0;i<four->n_seqs;i++) {
	four->taxon_name[i]=(char *)malloc(2*sizeof(char));
	assert((int) four->taxon_name[i]);
	four->taxon_name[i][0]='F';
	four->taxon_name[i][1]='\0';
	}
four->s=(char **)malloc((four->type_weight+four->n_seqs)*sizeof(char *));
assert((int) four->s);
for (i=0;i<(four->type_weight+four->n_seqs);i++) {
	four->s[i]=(char *)malloc((1+a1->length)*sizeof(char));
	assert((int) four->s[i]);
	four->s[i][four->length]='\0';
	}
/*to make s's
	0 = current node a1->s[0],a2->s[1],a3->s[1];
	1 = its ancestor a1->s[1]
	2 = desc0		  a2->s[0]
	3 = desc1        a3->s[0]
*/
for (i=0;i<a1->length;i++) {
	four->s[0][i]=a1->s[0][i];
	four->s[1][i]=a1->s[1][i];
	four->s[2][i]=a2->s[0][i];
	four->s[3][i]=a3->s[0][i];
	}
/*type weight stuff*/
if (four->type_weight) for (i=0;i<four->length;i++) four->s[4][i]=a1->s[1][i];
return four;
}


int **get_correspondances(anc,desc,nbases)
char *anc, *desc;
int nbases;
{
int i,j,k;
int **corres;
int gaps=0;

corres=(int **)malloc(2*sizeof(int *));
assert((int)corres);
for (i=0;i<2;i++) {
	corres[i]=(int *)malloc((nbases+1)*sizeof(int));
	assert((int) corres[i]);
	}

j=0;
/*what of leading gaps*/
/*k=0;
for (i=0;i<nbases;i++)  {
	if ((desc[i]=='-') && (anc[i]=='-')) ++k;
	else break;
	}
corres[0][0]=corres[1][0]=k;*/
for (i=0;i<nbases;i++) { /*if (!((desc[i]=='-') && (anc[i]=='-'))) {*/
	if (desc[i]!='-') corres[0][i]=0;
	else corres[0][i]=1;
	if (anc[i]!='-') corres[1][i]=0;
	else corres[1][i]=1;
	}
corres[0][nbases]=corres[1][nbases]=-1;
/*
gaps=0;
for (i=0;i<nbases;i++) if (corres[0][i]!=corres[1][i]) {
	if ((corres[0][i]==1) || (corres[1][i])==1) ++gaps;
	}
printf("\n%s\n%s",anc,desc);
printf("corres gaps %d ",gaps);
for (i=0;i<nbases;i++) printf("%c",desc[i]);
printf("\n");
for (i=0;i<nbases;i++) printf("%d",corres[0][i]);
printf("\n");
for (i=0;i<nbases;i++) printf("%c",anc[i]);
printf("\n");
for (i=0;i<nbases;i++) printf("%d",corres[1][i]);
printf("\n");*/
return corres;
}

void modify_corres(corres,tw)
alignment *tw;
int **corres;
{
int i,j;

for (i=0;i<tw->length;i++) {
	if ((tw->s[1][i]!=1) && (tw->s[1][i]!=65)) corres[0][i]=corres[1][i]=0;
	if (tw->s[1][i]>=64) {
		corres[0][i]+=64;
		corres[1][i]+=64;
		}
	}
}

alignment *big_make_neighbor (input_a1,input_a2,input_a3,values)
alignment *input_a1, *input_a2, *input_a3;
parameters *values;
{
	int implicit_gap_size,smaller,gap_not_too_big;
	alignment *b,**a, **a1, **a2, **a3;
	int check_diffs;
	char iter_holder,max_gap_holder,fix_condition,low_mem_holder;
	char *b_string1,*b_string2;
	int i, number_of_data_sets,j;
	int current_start, total_length, ii, holder;
	int leading_holder, trailing_holder,thang;

	if ((values->number_of_input_alignments>1) && (values->chop)) {
		leading_holder=values->leading_gap_cost;
		trailing_holder=values->trailing_gap_cost;
		values->trailing_gap_cost=0;
		values->leading_gap_cost=0;
	}

	/*allocate stuff for the number*/
	if (values->number_of_input_alignments==0) {
		input_a1->type_weight=input_a2->type_weight=input_a3->type_weight=0;
		number_of_data_sets=1;
	}
	else number_of_data_sets=values->number_of_input_alignments;
	a=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a);
	a1=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a1);
	a2=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a2);
	a3=(alignment **)malloc(number_of_data_sets*sizeof(alignment *));
	assert((int) a3);
	for (i=0;i<number_of_data_sets;i++) {
		a[i]=NULL;
		a1[i]=NULL;
		a2[i]=NULL;
		a3[i]=NULL;
		a1[i]=make_align(input_a1);
		a2[i]=make_align(input_a2);
		a3[i]=make_align(input_a3);
		/*split alignments*/
		if (values->number_of_input_alignments>0) if ((values->new_optimization) || (values->chop>0)) {
			redo_for_type(a1[i],i,values);
			redo_for_type(a2[i],i,values);
			redo_for_type(a3[i],i,values);
			}
	}
	/*fprintf(stderr,"{");
	for (ii=0;ii<number_of_data_sets;ii++) {
		print_inter_dig(a1[ii],values);
		print_inter_dig(a2[ii],values);
		print_inter_dig(a3[ii],values);
		fprintf(stderr,"%d=>%d %d %d ",ii,a1[ii]->length,a2[ii]->length,a3[ii]->length);
		}
	fprintf(stderr,"}\n");*/

	/*exit(-1);*/
	for (ii=0;ii<number_of_data_sets;ii++) {
		if (values->data_sets[ii]==0) a[ii]=make_neighbor(a1[ii],a2[ii],a3[ii],values);
		else a[ii]=make_neighbor_non_align(a1[ii],a2[ii],a3[ii],values);
		}

	/*fprintf(stderr,"\n");*/
	/*put them back together*/
	if (values->number_of_input_alignments==0) b=make_align(a[0]);
	else b=glue_back_after_chop(a,values);

	thang=0;
	for (i=0;i<number_of_data_sets;i++) thang+=a[i]->length;
	fprintf(stderr,"BMN2{%d(%d)->%d %d %d %d}",b->length,thang,strlen(b->s[0]),strlen(b->s[1]),strlen(b->s[2]),strlen(b->s[3]));

	if (values->new_optimization && b->type_weight) {
		b->score=0;
		for (i=0;i<number_of_data_sets;i++) b->score+=(values->data_set_weights[i]*a[i]->score);
		}
	fprintf(stderr," F ");
	for (i=0;i<number_of_data_sets;i++) {
		fprintf(stderr,"%d ",i);
		if (a[i]) a[i]=dump_align(a[i]);
		if (a1[i]) a1[i]=dump_align(a1[i]);
		if (a2[i]) a2[i]=dump_align(a2[i]);
		if (a3[i]) a3[i]=dump_align(a3[i]);
		fprintf(stderr,".");
	}
	free(a);
	free(a1);
	free(a2);
	free(a3);
	fprintf(stderr,"D");
	if ((values->number_of_input_alignments>1) && (values->chop)) {
		values->trailing_gap_cost=trailing_holder;
		values->leading_gap_cost=leading_holder;
		}
	if (b->type_weight!=1) b->type_weight=0;
	fprintf(stderr,"BMN2{%d->%d %d %d %d}",b->length,strlen(b->s[0]),strlen(b->s[1]),strlen(b->s[2]),strlen(b->s[3]));
	return b;
}
