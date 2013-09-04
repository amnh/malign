/*Copyright 1993 Ward Wheeler all rights reserved*/

#include "align3.h"

alignment **read_input_align(values,a,ret_ntax)
alignment **a;
parameters *values;
int *ret_ntax;
{
char buffer[LINESIZE],**sequences;
char **name;
int i,j,ntaxa,nlines,k,best_length;
char temp[LINESIZE];
int max_length=0;

scanf("%d %d",&ntaxa,&nlines);
if ((ntaxa<1) || (nlines<1)) {
	if (values->rep_error) fprintf(stderr,"Problems in numbers of taxa and/or lines\n");
	exit(-1);
	}
name=(char **)malloc(ntaxa*sizeof(char *));
assert((int)name);
sequences=(char **)malloc(ntaxa*sizeof(char *));
assert((int)sequences);
for (i=0;i<ntaxa;i++) {
	scanf("%s %s ",temp,buffer);
	name[i]=(char *)malloc((strlen(temp)+1)*sizeof(char));
	assert((int)name[i]);
	strcpy(name[i],temp);
	sequences[i]=(char *)malloc((strlen(buffer)+1)*sizeof(char));
	assert((int)sequences[i]);
	strcpy(sequences[i],buffer);
	}


for (j=1;j<nlines;j++) {
	for (i=0;i<ntaxa;i++) {
		scanf("%s %s ",temp,buffer);
		if (strcmp(temp,name[i])) fprintf(stderr,"Name mismatch -- could be error in file.\n");
		sequences[i]=(char *)realloc(sequences[i],strlen(sequences[i])+strlen(buffer)+1);
		assert((int)sequences[i]);
	sequences[i]=(char *)strcat(sequences[i],buffer);
		}
	}
for (i=0;i<ntaxa;i++) {if (strlen(sequences[i]) > max_length) max_length=strlen(sequences[i]);}
max_length+=1;

a=(alignment **)malloc(((2*(ntaxa+1))-1)*sizeof(alignment *));
assert((int)a);
for (i=0;i<ntaxa;i++) {
	a[i]=(alignment *)malloc(sizeof(alignment));
	assert((int)a[i]);
	a[i]->s=(char **)malloc(1*sizeof(char *));
	assert((int)a[i]->s);
	a[i]->taxon_name=(char **)malloc(1*sizeof(char *));
	assert((int)a[i]->taxon_name);
	a[i]->s[0]=(char *)malloc(max_length*sizeof(char));
	assert((int)a[i]->s[0]);
	strcpy(a[i]->s[0],sequences[i]);
	for (j=strlen(sequences[i]);j<max_length-1;j++) {a[i]->s[0][j]='-';}
	a[i]->s[0][max_length-1]='\0';
	a[i]->taxon_name[0]=(char *)malloc((strlen(name[i])+1)*sizeof(char));
	assert((int)a[i]->taxon_name[0]);
	strcpy(a[i]->taxon_name[0],name[i]);
	a[i]->name=(char *)malloc((strlen(name[i])+1)*sizeof(char));
	assert((int)a[i]->name);
	strcpy(a[i]->name,name[i]);
	/*printf("%s %s\n\n",name[i],a[i]->s[0]);*/
	a[i]->n_seqs=1;
	a[i]->length=strlen(a[i]->s[0]);
	a[i]->score=0;
	}
a[ntaxa]=(alignment *)malloc(sizeof(alignment));
assert((int)a[ntaxa]);
a[ntaxa]->s=(char **)malloc(ntaxa*sizeof(char *));
assert((int)a[ntaxa]->s);
a[ntaxa]->taxon_name=(char **)malloc(ntaxa*sizeof(char *));
assert((int)a[ntaxa]->taxon_name);
for (i=0;i< ntaxa; i++) {
	a[ntaxa]->s[i]=(char *)malloc((strlen(sequences[i])+1)*sizeof(char));
	assert((int)a[ntaxa]->s[i]);
	strcpy(a[ntaxa]->s[i],sequences[i]);
	a[ntaxa]->taxon_name[i]=(char *)malloc((strlen(name[i])+1)*sizeof(char));
	assert((int)a[ntaxa]->taxon_name[i]);
	strcpy(a[ntaxa]->taxon_name[i],name[i]);
	/*printf("%s %s\n",name[i],a[ntaxa]->s[i]);*/
	a[ntaxa]->name=(char *)malloc(30*sizeof(char));
	assert((int)a[ntaxa]->name);
	sprintf(a[ntaxa]->name,"(%s)","Alignment order unknown");
	}
a[ntaxa]->score=0;
a[ntaxa]->n_seqs=ntaxa;
a[ntaxa]->length=strlen(a[ntaxa]->s[0]);
if (values->VERBOSE) fprintf(stderr,"%d taxa and %d bases entered\n",ntaxa,a[ntaxa]->length);

values->best_rep=(int **)malloc(values->keep_trees*sizeof(int *));
assert((int)values->best_rep);
for (i=0;i<values->keep_trees;i++) {
	values->best_rep[i]=(int *)malloc((ntaxa-3)*sizeof(int));
	assert((int)values->best_rep[i]);
	}
values->number_best_clad=0;
values->length_best_clad=HUGE_COST;
values->all_done=ntaxa;

/*if not cost only
	move alignment to temp
	copy taxa to new alignments
	strip gaps
*/
for (i=0;i<ntaxa;i++) {
	free(name[i]);
	free(sequences[i]);
	}
free(name);
free(sequences);

values->all_done=ntaxa;
for (i=1;i<ntaxa;i++) for (j=0;j<a[i]->length;j++) if (a[i]->s[0][j]=='.') a[i]->s[0][j]=a[0]->s[0][j];
for (i=ntaxa+1;i<(2*(values->all_done+1))-1;i++) a[i]=NULL;
*ret_ntax=ntaxa;
if ((!values->cost_only) && (values->chop==0)) {
	if (values->VERBOSE) fprintf(stderr,"Stripping out gaps before realignment.\n");
	for (i=0;i<ntaxa;i++) a[i]=strip_gaps(a[i],values);
	}
return a;
}


alignment *strip_gaps(a,values)
parameters *values;
alignment *a;
{
int i,num_holder;
char *base_holder;

num_holder=0;
base_holder=(char *)malloc(a->length*sizeof(char));
assert((int)base_holder);
for (i=0;i<a->length;i++) {
	if(a->s[0][i]!='-') {
	base_holder[num_holder]=a->s[0][i];
		++num_holder;
		}
	}
a->length=num_holder;
a->s[0]=(char *)realloc(a->s[0],(num_holder+1)*sizeof(char));
assert((int)a->s[0]);
for (i=0;i<num_holder;i++) a->s[0][i]=base_holder[i];
a->s[0][num_holder]='\0';
free(base_holder);
/*printf("%s %s %d %s\n\n",a->taxon_name[0],a->name,a->length,a->s[0]);*/
return a;
}

/*
void alter(min_val,diag,down,right,pmi,pmi1)
int min_val,*diag,*down,*right;
cell pmi,pmi1;
{
int num_equal=0;

if (min_val==(*diag)) ++num_equal;
if (min_val==(*down)) ++num_equal;
if (min_val==(*right)) ++num_equal;

if (num_equal >1) {
	if ((*diag)==min_val) {
	 	if (pmi1.direction==DOWN) ++(*diag);
	 	else if (pmi.direction==RIGHT) ++(*diag);
	 	}
	}
}
*/

void alter_low_mem(min_val,diag,down,right,pmi,pmi1)
int min_val,*diag,*down,*right;
cell pmi,pmi1;
{
int num_equal=0;

if (min_val==(*diag)) ++num_equal;
if (min_val==(*down)) ++num_equal;
if (min_val==(*right)) ++num_equal;

if (num_equal >1) {
	if ((*diag)==min_val) {
	 	if (pmi1.direction==DOWN) ++(*diag);
	 	else if (pmi.direction==RIGHT) ++(*diag);
	 	}
	}
}

/*must add in reordering of groups file!!!*/
alignment **randomize_taxon_order(a,values,groups)
alignment **a;
parameters *values;
int **groups;
{
int i,j,*ran_nums,*new_order;
alignment **a_holder;
int temp;
char did_something,temp0,temp1;
int *groups_holder;

if (values->VERBOSE) {
	fprintf(stderr,"	Randomizing order...");
	if (values->jackboot) fprintf(stderr,"Save 0...");
	}
ran_nums=(int *)malloc(values->all_done*sizeof(int));
assert((int)ran_nums);
new_order=(int *)malloc(values->all_done*sizeof(int));
assert((int)new_order);
for (i=0;i<values->all_done;i++) {
	new_order[i]=i;
	ran_nums[i]=rand();
	}
a_holder=(alignment **)malloc((values->all_done)*sizeof(alignment *));
assert((int)a_holder);
for (i=0;i<values->all_done;i++) {
	a_holder[i]=make_align(a[i]);
	a[i]=dump_align(a[i]);
	}

/*get new order*/
if (values->jackboot) ran_nums[0]=-1;
do_reorder_thing:;
did_something=0;
for (i=0;i<values->all_done-1;i++) {
	if (ran_nums[i] > ran_nums[i+1]) {
		temp=ran_nums[i];
		ran_nums[i]=ran_nums[i+1];
		ran_nums[i+1]=temp;
		temp=new_order[i];
		new_order[i]=new_order[i+1];
		new_order[i+1]=temp;
		did_something=1;
		}
	}
if (did_something) goto do_reorder_thing;

for (i=0;i<values->all_done;i++) {
	a[i]=make_align(a_holder[new_order[i]]);
	a_holder[new_order[i]]=dump_align(a_holder[new_order[i]]);
	}

/*reorder groups if there*/
if (groups) {
	groups_holder=(int *)malloc(values->all_done*sizeof(int));
	assert((int)groups_holder);
	for (i=0;i<values->ngroups;i++){
		for (j=0;j<values->all_done;j++) groups_holder[j]=groups[i][new_order[j]];
		for (j=0;j<values->all_done;j++) groups[i][j]=groups_holder[j];
		}
	free(groups_holder);
	}

free(new_order);
free(ran_nums);
free(a_holder);
if (values->VERBOSE) fprintf(stderr,"done_rr.\n");
return a;
}

/*
alignment **reorder(a,ntax,type,advancement_scores,groups,values)
alignment **a;
int ntax,type,**advancement_scores;
char **groups;
parameters *values;
{
int best_score,*new_order,*done;
int i,j,first[2],tempij,tot_dist;
char **name,**input_sequence,*groups_holder;
alignment *temp_align;
char temp0,temp1;


if (values->VERBOSE) fprintf(stderr,"Ordering sequences...");
fflush(stderr);

new_order=(int *)malloc(ntax*sizeof(int));
assert((int)new_order);
done=(int *)malloc(ntax*sizeof(int));
assert((int)done);

tot_dist=0;
for (i=0;i<ntax;i++) for (j=i+1;j<ntax;j++) tot_dist+=advancement_scores[i][j];
if (tot_dist==0) {
	for (i=0;i<ntax;i++) {
		fprintf(stderr,"%d ",i+1);
		fflush(stderr);
		for (j=(i+1);j<ntax;j++) {
			temp_align=nw(a[i],a[j],values);
			advancement_scores[i][j]=advancement_scores[j][i]=temp_align->score;
			if (temp_align->score==0) fprintf(stderr,"Sequences %d and %d are identical -- this may cause problems later.\n",i,j);
			temp_align=dump_align(temp_align);
			}
		}
	}
name=(char **)malloc(ntax*sizeof(char *));
assert((int)name);
input_sequence=(char **)malloc(ntax*sizeof(char *));
assert((int)input_sequence);

new_order[0]=0;
new_order[1]=1;

if (type==0) {
	best_score=HUGE_COST;
	for (i=0;i<ntax;i++) {
		for (j=(i+1);j<ntax;j++) {
			if (*(advancement_scores[i]+j)<=best_score) {
				best_score=(*(advancement_scores[i]+j));
				first[0]=i;
				first[1]=j;
				}
			}
		}
	if (strcmp(a[first[0]]->s[0],a[first[1]]->s[0])<0) {
		tempij=first[0];
		first[0]=first[1];
		first[1]=tempij;
		}
	for (i=0;i<2;i++) {
		name[i]=(char *)malloc((strlen(a[first[i]]->name)+1)*sizeof(char *));
		assert((int)name[i]);
		strcpy(name[i],a[first[i]]->name);
		input_sequence[i]=(char *)malloc((strlen(a[first[i]]->s[0])+1)*sizeof(char *));
		assert((int)input_sequence[i]);
		strcpy(input_sequence[i],a[first[i]]->s[0]);
		}
	for (i=0;i<2;i++) {
		free(a[first[i]]->name);
		a[first[i]]->name=(char *)malloc((strlen(a[i]->name)+1)*sizeof(char));
		assert((int)a[first[i]]->name);
		strcpy(a[first[i]]->name,a[i]->name);
		a[first[i]]->s[0]=(char *)malloc((strlen(a[i]->s[0])+1)*sizeof(char));
		assert((int)a[first[i]]->s[0]);
		strcpy(a[first[i]]->s[0],a[i]->s[0]);
		a[first[i]]->length=strlen(a[i]->s[0]);
		}
	for (i=0;i<2;i++) {
		free(a[i]->name);
		a[i]->name=(char *)malloc((strlen(name[i])+1)*sizeof(char));
		assert((int)a[i]->name);
		strcpy(a[i]->name,name[i]);
		free(name[i]);
		free(a[i]->s[0]);
		a[i]->s[0]=(char *)malloc((strlen(input_sequence[i])+1)*sizeof(char));
		assert((int)a[i]->s[0]);
		strcpy(a[i]->s[0],input_sequence[i]);
		a[i]->length=strlen(input_sequence[i]);
		free(input_sequence[i]);
		}
	best_score=HUGE_COST;
	for (i=2;i<ntax;i++) *(done+i)=0;
	for (j=2;j<ntax;j++) {
		for (i=2;i<ntax;i++){
			if ((*(advancement_scores[first[0]]+i)<=best_score)&&(*(done+i)==0)) {
				*(new_order+j)=i;
				best_score=(*(advancement_scores[first[0]]+i));
				}
			}
		best_score=HUGE_COST;
		*(done+(*(new_order+j)))=1;
		}
	}
else if (type==1) {
	best_score=0;
	for (i=0;i<ntax;i++) {
		for (j=(i+1);j<ntax;j++) {
			if (*(advancement_scores[i]+j)>=best_score) {
				best_score=(*(advancement_scores[i]+j));
				first[0]=i;
				first[1]=j;
				}
			}
		}
	if (strcmp(a[first[0]]->s[0],a[first[1]]->s[0])<0) {
		tempij=first[0];
		first[0]=first[1];
		first[1]=tempij;
		}
	for (i=0;i<2;i++) {
		name[i]=(char *)malloc((strlen(a[first[i]]->name)+1)*sizeof(char *));
		assert((int)name[i]);
		strcpy(name[i],a[first[i]]->name);
		input_sequence[i]=(char *)malloc((strlen(a[first[i]]->s[0])+1)*sizeof(char *));
		assert((int)input_sequence[i]);
		strcpy(input_sequence[i],a[first[i]]->s[0]);
		}
	for (i=0;i<2;i++) {
		free(a[first[i]]->name);
		a[first[i]]->name=(char *)malloc((strlen(a[i]->name)+1)*sizeof(char));
		assert((int)a[first[i]]->name);
		strcpy(a[first[i]]->name,a[i]->name);
		a[first[i]]->s[0]=(char *)malloc((strlen(a[i]->s[0])+1)*sizeof(char));
		assert((int)a[first[i]]->s[0]);
		strcpy(a[first[i]]->s[0],a[i]->s[0]);
		a[first[i]]->length=strlen(a[i]->s[0]);
		}
	for (i=0;i<2;i++) {
		free(a[i]->name);
		a[i]->name=(char *)malloc((strlen(name[i])+1)*sizeof(char));
		assert((int)a[i]->name);
		strcpy(a[i]->name,name[i]);
		free(name[i]);
		free(a[i]->s[0]);
		a[i]->s[0]=(char *)malloc((strlen(input_sequence[i])+1)*sizeof(char));
		assert((int)a[i]->s[0]);
		strcpy(a[i]->s[0],input_sequence[i]);
		a[i]->length=strlen(input_sequence[i]);
		free(input_sequence[i]);
		}
	best_score=0;
	for (i=2;i<ntax;i++) *(done+i)=0;
	for (j=2;j<ntax;j++) {
		for (i=2;i<ntax;i++){
			if ((*(advancement_scores[first[0]]+i)>=best_score)&&(*(done+i)==0)) {
				*(new_order+j)=i;
				best_score=(*(advancement_scores[first[0]]+i));
				}
			}
		best_score=0;
		*(done+(*(new_order+j)))=1;
		}
	}
if (groups) {
	groups_holder=(char *)malloc(ntax*sizeof(char));
	assert((int)groups_holder);
	for (i=0;i<values->ngroups;i++){
		temp0=(*(groups[i]+0));
		temp1=(*(groups[i]+1));
		*(groups[i]+0)=(*(groups[i]+first[0]));
		*(groups[i]+1)=(*(groups[i]+first[1]));
		*(groups[i]+first[0])=temp0;
		*(groups[i]+first[1])=temp1;
		for (j=2;j<ntax;j++) *(groups_holder+j)=(*(groups[i]+(*(new_order+j))));
		for (j=2;j<ntax;j++) *(groups[i]+j)=(*(groups_holder+j));
		}
	free(groups_holder);
	}

for (i=2;i<ntax;i++){
	name[i]=(char *)malloc((strlen(a[(*(new_order+i))]->name)+1)*sizeof(char *));
	assert((int)name[i]);
	strcpy(name[i],a[(*(new_order+i))]->name);
	input_sequence[i]=(char *)malloc((strlen(a[(*(new_order+i))]->s[0])+1)*sizeof(char *));
	assert((int)input_sequence[i]);
	strcpy(input_sequence[i],a[(*(new_order+i))]->s[0]);
	}
for (i=2;i<ntax;i++){
	free(a[i]->name);
	a[i]->name=(char *)malloc((strlen(name[i])+1)*sizeof(char));
	assert((int)a[i]->name);
	strcpy(a[i]->name,name[i]);
	free(a[i]->s[0]);
	a[i]->s[0]=(char *)malloc((strlen(input_sequence[i])+1)*sizeof(char));
	assert((int)a[i]->s[0]);
	strcpy(a[i]->s[0],input_sequence[i]);
	a[i]->length=strlen(input_sequence[i]);
	}
for (i=2;i<ntax;i++){
	free(name[i]);
	free(input_sequence[i]);
	}
free(name);
free(input_sequence);
free(done);
free(new_order);


for (i=0;i<ntax;i++) {
	free(a[i]->taxon_name[0]);
	a[i]->taxon_name[0]=(char *)malloc((1+(strlen(a[i]->name)))*sizeof(char));
	assert((int)a[i]->taxon_name[0]);
	strcpy(a[i]->taxon_name[0],a[i]->name);
	}

if (values->VERBOSE) fprintf(stderr,"\n");
return a;
}
*/

alignment **new_reorder(a,ntax,type,advancement_scores,groups,values)
alignment **a;
int ntax,type,**advancement_scores;
int **groups;
parameters *values;
{
int i,j,*ran_nums,*new_order;
alignment **a_holder;
int temp,best_score;
char did_something,temp0,temp1;
int *groups_holder;
char *done;
alignment *temp_align;
int tot_dist,tempij,first[2];
int lowest,lowest_num,this_one;

if (values->VERBOSE) fprintf(stderr,"Reordering ");
ran_nums=(int *)malloc(values->all_done*sizeof(int));
assert((int)ran_nums);
new_order=(int *)malloc(values->all_done*sizeof(int));
assert((int)new_order);
done=(char *)malloc(values->all_done*sizeof(char));
assert((int)new_order);

a_holder=(alignment **)malloc((values->all_done)*sizeof(alignment *));
assert((int)a_holder);
for (i=0;i<values->all_done;i++) {
	a_holder[i]=make_align(a[i]);
	a[i]=dump_align(a[i]);
	}

/*if done before don't again*/
tot_dist=0;
for (i=0;i<ntax;i++) for (j=i+1;j<ntax;j++) tot_dist+=advancement_scores[i][j];
if (!tot_dist) {
	for (i=0;i<ntax;i++) {
		if (values->VERBOSE) fprintf(stderr,"%d ",i+1);
		if (values->VERBOSE) fflush(stderr);
		for (j=(i+1);j<ntax;j++) {
			temp_align=nw(a_holder[i],a_holder[j],values);
			advancement_scores[i][j]=advancement_scores[j][i]=temp_align->score;
			if (0) {
				printf("%d + %d => %d\n",i,j,temp_align->score);
				}
			if (!temp_align->score) if (values->rep_error) fprintf(stderr,"Sequences %d and %d are identical -- this may cause problems later.\n",i,j);
			temp_align=dump_align(temp_align);
			}
		}
	}

if (!type) {
	/*find best pair*/
	best_score=HUGE_COST;
	for (i=0;i<ntax;i++) {
		for (j=(i+1);j<ntax;j++) {
			if (advancement_scores[i][j]<=best_score) {
				best_score=(advancement_scores[i][j]);
				first[0]=i;
				first[1]=j;
				}
			}
		}
  }
else {
	/*take furthest to start them the column of the appropriate one*/
	best_score=0;
	for (i=0;i<ntax;i++) {
		for (j=(i+1);j<ntax;j++) {
			if (advancement_scores[i][j]>=best_score) {
				best_score=(advancement_scores[i][j]);
				first[0]=i;
				first[1]=j;
				}
			}
		}
	}
/*switch these i and j */
if (strcmp(a_holder[first[0]]->s[0],a_holder[first[1]]->s[0])<0) {
	tempij=first[0];
	first[0]=first[1];
	first[1]=tempij;
	}
new_order[0]=first[0];
new_order[1]=first[1];

/*get new order*/
for (i=0;i<ntax;i++) done[i]=0;
done[first[0]]=1;
done[first[1]]=1;
if (!type) {
	for (j=2;j<ntax;j++) {
		lowest=HUGE_COST;
		for (i=0;i<ntax;i++) if (!done[i]) {
      	if (advancement_scores[first[0]][i]<=lowest) {
	    	lowest=advancement_scores[first[0]][i];
	      lowest_num=i;
	    	}
      }
    this_one=0;
    for (i=0;i<ntax;i++) this_one+=done[i];
    new_order[this_one]=lowest_num;
    done[lowest_num]=1;
    }
	}
else {
	for (j=2;j<ntax;j++) {
		lowest=0;
		for (i=0;i<ntax;i++) if (!done[i]) {
	  	if (advancement_scores[first[0]][i]>=lowest) {
	    	lowest=advancement_scores[first[0]][i];
	      lowest_num=i;
	    	}
      }
    this_one=0;
    for (i=0;i<ntax;i++) this_one+=done[i];
    new_order[this_one]=lowest_num;
    done[lowest_num]=1;
  	}
	}

/*copy back*/
for (i=0;i<values->all_done;i++) {
	a[i]=make_align(a_holder[new_order[i]]);
	a_holder[new_order[i]]=dump_align(a_holder[new_order[i]]);
	}

/*reorder groups if there*/
if (groups) {
	groups_holder=(int *)malloc(values->all_done*sizeof(int));
	assert((int)groups_holder);
	for (i=0;i<values->ngroups;i++){
		for (j=0;j<values->all_done;j++) groups_holder[j]=groups[i][new_order[j]];
		for (j=0;j<values->all_done;j++) groups[i][j]=groups_holder[j];
		}
	free(groups_holder);
	}

free(new_order);
free(ran_nums);
free(a_holder);
free(done);
if (values->VERBOSE) fprintf(stderr,"done.\n");
return a;
}

