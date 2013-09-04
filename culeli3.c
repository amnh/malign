/*
Copyright 1992 Ward Wheeler all rights reserved

routines for cull and elision

check and make sure that the taxa are reordered so correct cull and elision can be done
*/

#include "align3.h"

alignment **do_elision(a,values,best_aligns,num_best_aligns)
alignment **a,**best_aligns;
int *num_best_aligns;
parameters *values;
{
int i,j,k, total_length, old_length;
int ***weight_info, *weights;
int change_ratio, gap_ratio;
alignment *temp_align;

weight_info=(int ***)malloc(values->number_of_input_alignments*sizeof(int **));
assert((int) weight_info);
for (i=0;i<values->number_of_input_alignments;i++) {
	weight_info[i]=(int **)malloc(2*sizeof(int *));
	assert((int) weight_info[i]);
	for (j=0;j<2;j++) {
		weight_info[i][j]=(int *)malloc(2*sizeof(int));
		assert((int) weight_info[i][j]);
		}
	}

/*add info from files
if ((gap_ratio/change_ratio)<1) gap_ratio=change_ratio=1;
else {
	gap_ratio=gap_ratio/change_ratio;
	change_ratio=1;
	}
*/
change_ratio=1;
gap_ratio=2;
for (i=0;i<values->number_of_input_alignments;i++) {
	weight_info[i][0][0]=a[i]->length;
	weight_info[i][0][1]=change_ratio;
	weight_info[i][1][1]=gap_ratio;	
	/*fprintf(stderr,"%d->%d %d %d\n",i,	weight_info[i][0][0],	weight_info[i][0][1],	weight_info[i][1][1]);*/
	}

if (values->VERBOSE) fprintf(stderr,"Eliding...");
total_length=0;
for (i=0;i<values->number_of_input_alignments;i++) 	{
	recode_with_gaps(a[i]);
	total_length+=a[i]->length;
	}
weights=(int *)malloc(total_length*sizeof(int));
assert((int) weights);	
total_length=0;
for (i=0;i<values->number_of_input_alignments;i++) {
	old_length=total_length;
	total_length+=a[i]->length;
	for (j=old_length;j<(old_length+weight_info[i][0][0]);j++) weights[j]=weight_info[i][0][1];
	for (j=(old_length+weight_info[i][0][0]);j<total_length;j++) weights[j]=weight_info[i][1][1];
	}
for (i=1;i<values->number_of_input_alignments;i++) {
	for (j=0;j<a[i]->n_seqs;j++) {
			a[0]->s[j]=(char *)realloc(a[0]->s[j],(1+a[0]->length+a[i]->length)*sizeof(char));
			assert((int) a[0]->s[j]);
			for (k=a[0]->length;k<(a[0]->length+a[i]->length);k++) {
				a[0]->s[j][k]=a[i]->s[j][k-a[0]->length];
				}
			a[0]->s[j][(a[0]->length+a[i]->length)]='\0';
			}
		a[0]->length=a[0]->length+a[i]->length;
		}
if (values->print_intermediates) {
	printf("Elision 'Grand alignment':\n");
	printem(a,1,values);
	}
/* compact*/
if (values->VERBOSE) fprintf(stderr,"compacting...");
for (i=0;i<a[0]->length;i++) {
	for (j=i+1;j<a[0]->length;j++) {
		if (weights[j]!=0) {
			if (same_position(i,j,a[0])) {
				weights[i]+=weights[j];
				weights[j]=0;
				}
			}		
		}
	}
total_length=0;
for (i=0;i<a[0]->length;i++) if (weights[i]>0) ++total_length;
values->ce_weights=(int *)malloc(total_length*sizeof(int));
assert((int)values->ce_weights);
temp_align=make_align(a[0]);
temp_align->length=total_length;
for (i=0;i<a[0]->n_seqs;i++) {
	temp_align->s[i]=(char *)realloc(temp_align->s[i],(1+temp_align->length)*sizeof(char));
	assert((int) temp_align);
	k=0;
	for (j=0;j<a[0]->length;j++) if (weights[j]>0) temp_align->s[i][k++]=a[0]->s[i][j];
	temp_align->s[i][temp_align->length]='\0';
	}
	
/*copy back*/
if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
best_aligns[0]=make_align(temp_align);
*num_best_aligns=1;
k=0;
for (i=0;i<a[0]->length;i++) if (weights[i]>0) values->ce_weights[k++]=weights[i];


/*free*/
temp_align=dump_align(temp_align);
free(weights);
for (i=0;i<values->number_of_input_alignments;i++){
	 for (j=0;j<2;j++) free(weight_info[i][j]);
	 free(weight_info[i]);
	 }
if (values->VERBOSE) fprintf(stderr,"done 'elision'\n");
return best_aligns;
}

alignment **do_cull(a,values,best_aligns,num_best_aligns)
alignment **a,**best_aligns;
int *num_best_aligns;
parameters *values;
{
int i,j,k;
char *agree;
alignment *temp_align;

if (*num_best_aligns>1) {
	fprintf(stderr,"Performing 'cull' on only the first of %d alignment-alignments.\n");
	fprintf(stderr,"	In order to 'cull' all of these alignment-alignments, print intermediate results and reinput them to 'cull' in a separate run.\n");
	for (i=1;i<(*num_best_aligns);i++) if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
	}

if (values->VERBOSE) fprintf(stderr,"Culling...");
agree=(char *)malloc(best_aligns[0]->length*sizeof(char));
assert((int) agree);
get_robust_positions(agree,best_aligns[0],a[0]->n_seqs,values->number_of_input_alignments);

/*print out accepted positions file*/
j=0;
for (i=0;i<best_aligns[0]->length;i++) if (agree[i]) ++j;

if (values->VERBOSE) fprintf(stderr,"Gatesy-DeSalle residuum...");
temp_align=make_align(a[0]);
temp_align->length=j;
for (i=0;i<temp_align->n_seqs;i++) {
	free(temp_align->s[i]);
	temp_align->s[i]=(char *)malloc((temp_align->length+1)*sizeof(char));
	assert((int) temp_align->s[i]);
 	}
if (values->VERBOSE) fprintf(stderr,"%d taxa %d characters remain\n",temp_align->n_seqs,temp_align->length);
	for (i=0;i<temp_align->n_seqs;i++) {
		k=0;
		for (j=0;j<best_aligns[0]->length;j++) if (agree[j]) {
			temp_align->s[i][k++]=best_aligns[0]->s[i][j];
			}
		temp_align->s[i][temp_align->length]='\0';
		}
	if (values->VERBOSE) fprintf(stderr,"finished 'cull'\n");
	
	/*This stuff to include gap costs
	gap_holder=values->gap_cost;
	change_holder=values->change_cost;
	values->gap_cost=average_gap_ratio;
	values->change_cost=average_change_ratio;
	*/

best_aligns[0]=dump_align(best_aligns[0]);
best_aligns[0]=make_align(temp_align);
temp_align=dump_align(temp_align);
free(agree);
*num_best_aligns=1;
return best_aligns;
}

void new_order_as_input(a,order_of_names,values)
alignment *a;
parameters *values;
char **order_of_names;
{
int i,j;
char **bs,**btaxon_name;
char found_match;

bs=(char **)malloc(a->n_seqs*sizeof(char *));
assert((int) bs);
btaxon_name=(char **)malloc(a->n_seqs*sizeof(char *));
assert((int) btaxon_name);

for (i=0;i<a->n_seqs;i++) {
	found_match=0;
	for (j=0;j<a->n_seqs;j++) {
		if (!stricmp(a->taxon_name[i],order_of_names[j])) {
			found_match=1;
			bs[j]=(char *)malloc((1+(strlen(a->s[i])))*sizeof(char));
			assert((int) bs[j]);
			strcpy(bs[j],a->s[i]);
			btaxon_name[j]=(char *)malloc((1+(strlen(a->taxon_name[i])))*sizeof(char));
			assert((int) btaxon_name[j]);
			strcpy(btaxon_name[j],a->taxon_name[i]);
			}
		}
	if (!found_match) {
		fprintf(stderr,"Found no match for %s. Bye.\n",a->taxon_name[i]);
		exit(-1);
		}
	}

for (i=0;i<a->n_seqs;i++) {
	free(a->s[i]);
	a->s[i]=(char *)malloc((1+a->length)*sizeof(char));
	assert((int) a->s[i]);
	strcpy(a->s[i],bs[i]);
	free(bs[i]);

	free(a->taxon_name[i]);
	a->taxon_name[i]=(char *)malloc((1+(strlen(btaxon_name[i])))*sizeof(char));
	assert((int) a->taxon_name[i]);
	strcpy(a->taxon_name[i],btaxon_name[i]);
	free(btaxon_name[i]);
	}
free(bs);
free(btaxon_name);
}

void get_robust_positions(agree,a,n_seqs,num_aligns)
char *agree;
alignment *a;
int n_seqs,num_aligns;
{
int i,j,k,counter;

for (i=0;i<a->length;i++) {
 	counter=0;
	for (j=1;j<num_aligns;j++) {
		for (k=0;k<n_seqs;k++) {
		if (a->s[k][i]==a->s[k+(j*n_seqs)][i]) ++counter;
	else break;
	}
		}
	if (counter==((num_aligns-1)*n_seqs)) agree[i]=1;
	else agree[i]=0;
	}	
}

void recode_with_gaps(a)
alignment *a;
{
int number_of_gap_positions;
int i,j,k,l,old_length;
char holder,*position_is_gap;
char leading,trailing;

number_of_gap_positions=0;
position_is_gap=(char *)malloc(a->length*sizeof(char));
assert((int) position_is_gap);

/*get number and position of gaps*/
/*add stuff for leading and trailing*/
for (i=0;i<a->length;i++) {
	holder=0;
	for (j=0;j<a->n_seqs;j++) if (a->s[j][i]=='-') holder=1;
	number_of_gap_positions+=holder;
	position_is_gap[i]=holder;
	}

/*reallocate and assign 0/1 gap characters*/
old_length=a->length;
a->length+=number_of_gap_positions;
for (i=0;i<a->n_seqs;i++) {
	a->s[i]=(char *)realloc(a->s[i],(1+a->length)*sizeof(char));
	assert((int) a->s[i]);
	a->s[i][a->length]='\0';
	k=0;
	leading=1;
	for (j=0;j<old_length;j++) {
		if ((leading) && (a->s[i][j]!='-')) leading=0;
	trailing=0;
	for (l=j;l<old_length;l++) {
		if ((a->s[i][l]!='-') && (!leading)) trailing=1;
		}
		if (position_is_gap[j]) {
		if (a->s[i][j]=='-') {
		if ((!leading) && (trailing)) a->s[i][old_length+k]='1';
	else a->s[i][old_length+k]='?';
	}
	else a->s[i][old_length+k]='0';
		k++;
		}
		}
	}

free(position_is_gap);
}

int same_position(first,second,a)
int first,second;
alignment *a;
{
int i;

for (i=0;i<a->n_seqs;i++) {
	if (a->s[i][first]!=a->s[i][second]) return 0;
	}
return 1;
}
