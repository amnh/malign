/*Copyright 1993 Ward Wheeler all rights reserved
need to add data-sets etc*/

#include "align3.h"
extern int charbit[256];
alignment **new_get_prealigned_sequences_chop(values)
parameters *values;
{
	char buffer[LINESIZE],**sequences;
	int **poly_seq;
	char **name;
	int i,j,ntaxa,nlines,k,best_length,number_of_sequences,*num_bs,done;
	char temp[LINESIZE];
	FILE *fp_in, *fopen();
	alignment **a,***b,**c;
	int old_length,total_length,current_type,original_0;
	int *start, *stop, old_num_alignments;
	int something_to_do,num_poly;
fprintf(stderr,"reading from here\n");
	/*establish*/
	a=NULL;
	c=NULL;
	b=NULL;
    poly_seq=NULL;

	old_num_alignments=values->number_of_input_alignments;
	if (values->chop>0) {
		if (!values->number_of_input_alignments) values->number_of_input_alignments=1;
		if (!values->data_sets) {
			values->data_sets=(int *)realloc(values->data_sets,(values->number_of_input_alignments)*sizeof(int));
			assert((int) values->data_sets);
			values->data_sets[0]=0;
			}
		if (!values->data_set_weights) {
			values->data_set_weights=(int *)realloc(values->data_set_weights,(values->number_of_input_alignments)*sizeof(int));
			assert((int) values->data_set_weights);
			values->data_set_weights[0]=1;
			}
		}
	values->all_done = values->number_of_input_alignments;

	/*allocate*/
	a=(alignment **)malloc(values->all_done*sizeof(alignment *));
	assert((int) a);
	for (j=0;j<values->all_done;j++) a[j]=NULL;
	b=(alignment ***)malloc(values->number_of_input_alignments*sizeof(alignment **));
	assert((int) b);
	for (i=0;i<values->all_done;i++) b[i]=NULL;
	number_of_sequences=0;
	num_bs=(int *)malloc(values->all_done*sizeof(int));
	assert((int) num_bs);

	/*get sequences from files here*/
	/*2=morph*/
	for (k=0;k<values->number_of_input_alignments;k++) {
		if (values->data_sets[k]<3) fp_in=(FILE *)fopen(values->input_file_name[k],"r");
		else fp_in=freopen(values->input_file_name[k],"r",stdin);
		if (!fp_in) {
			fprintf(stderr,"File %s not found\n",values->input_file_name[k]);
			exit(-1);
			}
		else if (values->VERBOSE) fprintf(stderr,"Reading sequences from file %s\n",values->input_file_name[k]);
		if (values->data_sets[k]<3) {
			fscanf(fp_in,"%d %d",&ntaxa,&nlines);
			if (!number_of_sequences) number_of_sequences=ntaxa;
			if ((ntaxa<1) || (nlines<1)) {
				if (values->rep_error) fprintf(stderr,"Problems in numbers of taxa and/or lines in alignment file %s\n",
			    values->input_file_name[k]);
					exit(-1);
					}
			name=(char **)malloc(ntaxa*sizeof(char *));
			assert((int)name);
			sequences=(char **)malloc(ntaxa*sizeof(char *));
			assert((int)sequences);
			for (i=0;i<ntaxa;i++) {
				fscanf(fp_in,"%s %s ",temp,buffer);
				name[i]=(char *)malloc((strlen(temp)+1)*sizeof(char));
				assert((int)name[i]);
				strcpy(name[i],temp);
				sequences[i]=(char *)malloc((strlen(buffer)+1)*sizeof(char));
				assert((int)sequences[i]);
				strcpy(sequences[i],buffer);
			}
			fprintf(stderr,"	");
			for (j=1;j<nlines;j++) {
				for (i=0;i<ntaxa;i++) {
					fscanf(fp_in,"%s %s ",temp,buffer);
					if (strcmp(temp,name[i])) fprintf(stderr,"Name mismatch -- could be error in file %s.\n",
					    values->input_file_name[k]);
					sequences[i]=(char *)realloc(sequences[i],strlen(sequences[i])+strlen(buffer)+1);
					assert((int)sequences[i]);
					sequences[i]=(char *)strcat(sequences[i],buffer);
				}
			}
			/*do ordering compensation*/
			if (values->data_sets[k]==2) {
			    make_all_unordered(ntaxa,sequences,fp_in);
			}
			/*copy to a[k]*/
			a[k]=convert_to_alignment(ntaxa,name,sequences,values->input_file_name[k]);
			if (values->VERBOSE) fprintf(stderr,"%d taxa and %d characters entered from file %s\n",ntaxa,a[k]->length,values->input_file_name[k]);
			/*free temp stuff*/
			for (i=0;i<ntaxa;i++) {
				free(name[i]);
				free(sequences[i]);
			}
			free(name);
			free(sequences);
		}
		else {
			/*genbank format style*/
			b[k]=(alignment **)new_read_old_file(b[k],values,&(num_bs[k]));
			/*values->all_done=num_bs[k];*/
			/*do_modifications(b[k],values);*/
			}
		fclose(fp_in);
		}
	if (values->VERBOSE) fprintf(stderr,"\n");

	/*check num_sequences*/
	for (i=0;i<values->number_of_input_alignments;i++) {
		if (a[i]) {
			if (a[i]->n_seqs!=number_of_sequences) {
				fprintf(stderr,"Input files contain different numbers of taxa/sequences.\n");
				exit(-1);
			}
		}
		else if (b[i]) {
			if (num_bs[i]!=number_of_sequences) {
				fprintf(stderr,"Input files contain different numbers of taxa/sequences.\n");
				exit(-1);
			}
		}
	}

	/*convert b's to a's*/
	for (i=0;i<values->number_of_input_alignments;i++)  if (b[i]) a[i]=convert_b_to_a(b[i],num_bs[i]);
	values->all_done=number_of_sequences;

	/*here add the chop stuff and return a new a vector + increment number of input and set types and weights etc*/
	if (values->chop>0) {
		something_to_do=0;
		for (i=0;i<values->number_of_input_alignments;i++) if (values->data_sets[i]==0) something_to_do=1;
		if (something_to_do) a=do_chop_thing(a,values);
		}

	/*check to make sure they all have the same number*/
	if (values->new_optimization) for (i=0;i<values->number_of_input_alignments;i++)  {
		/*convert morph to numbers*/
		if (values->data_sets[i]==2) {
			for (j=0;j<a[i]->n_seqs;j++) {
				for (k=0;k<a[i]->length;k++) {
					if (a[i]->s[j][k]=='0') a[i]->s[j][k]=1;
					else if (a[i]->s[j][k]=='1') a[i]->s[j][k]=2;
					else if (a[i]->s[j][k]=='2') a[i]->s[j][k]=4;
					else if (a[i]->s[j][k]=='3') a[i]->s[j][k]=8;
					else if (a[i]->s[j][k]=='4') a[i]->s[j][k]=16;
					else if (a[i]->s[j][k]=='5') a[i]->s[j][k]=32;
					else if (a[i]->s[j][k]=='6') a[i]->s[j][k]=64;
					else if (a[i]->s[j][k]=='7') a[i]->s[j][k]=128;
					else if (a[i]->s[j][k]=='?') a[i]->s[j][k]=127;
					else if (a[i]->s[j][k]=='-') a[i]->s[j][k]=127;
					else if (a[i]->s[j][k]&128)  a[i]->s[j][k]=a[i]->s[j][k]-128;
					/*else do nothing -- its a polymorphism preset negative values for polymorphsm*/
				}
			}
		}
		/*convert noalign to bits as well '-' converted to GAP when non-align*/
		else if (values->data_sets[i]==1) for (j=0;j<a[i]->n_seqs;j++) for (k=0;k<a[i]->length;k++) a[i]->s[j][k]=charbit[a[i]->s[j][k]];
		}

	/*reorder*/
	for (i=1;i<values->number_of_input_alignments;i++) new_order_as_input(a[i],a[0]->taxon_name,values);
	original_0=a[0]->length;

	/*elide inputs into a single set */
	total_length=0;
	for (i=0;i<values->number_of_input_alignments;i++) {
		old_length=total_length;
		total_length+=a[i]->length;
	}
	for (i=1;i<values->number_of_input_alignments;i++) {
		for (j=0;j<a[i]->n_seqs;j++) {
			a[0]->s[j]=(char *)realloc(a[0]->s[j],(1+a[0]->length+a[i]->length)*sizeof(char));
			assert((int) a[0]->s[j]);
			for (k=a[0]->length;k<(a[0]->length+a[i]->length);k++) a[0]->s[j][k]=a[i]->s[j][k-a[0]->length];
			a[0]->s[j][(a[0]->length+a[i]->length)]='\0';
		}
		a[0]->length=a[0]->length+a[i]->length;
	}

	/*split out to c's*/
	values->all_done=a[0]->n_seqs;
	c=(alignment **)malloc(((2*(values->all_done+1))-1)*sizeof(alignment *));
	assert((int) c);
	for (i=0;i<((2*(values->all_done+1))-1);i++) c[i]=NULL;
	for (i=0;i<values->all_done;i++) {
		c[i]=(alignment *)malloc(sizeof(alignment));
		assert((int) c[i]);
		c[i]->length=a[0]->length;
		c[i]->n_seqs=1;
		c[i]->score=0;
		c[i]->name=(char *)malloc((1+strlen(a[0]->taxon_name[i]))*sizeof(char));
		assert((int) c[i]->name);
		strcpy(c[i]->name,a[0]->taxon_name[i]);
		c[i]->type_weight=1;
		c[i]->s=(char **)malloc((1+c[i]->type_weight)*sizeof(char *));
		assert((int) c[i]->s);
		c[i]->s[0]=(char *)malloc((1+a[0]->length)*sizeof(char));
		assert((int) c[i]->s[0]);
		c[i]->s[0][c[i]->length]='\0';
		/*strcpy(c[i]->s[0],a[0]->s[i]);*/
		for (j=0;j<c[i]->length;j++) c[i]->s[0][j]=a[0]->s[i][j];
		c[i]->s[1]=(char *)malloc((1+a[0]->length)*sizeof(char));
		/*c[i]->s[1]=(char *)malloc((MAX_SEQUENCE_SIZE)*sizeof(char));*/
		assert((int) c[i]->s[1]);
		c[i]->taxon_name=(char **)malloc(1*sizeof(char *));
		assert((int) c[i]->taxon_name);
		c[i]->taxon_name[0]=(char *)malloc((1+strlen(a[0]->taxon_name[i]))*sizeof(char));
		assert((int) c[i]->taxon_name[0]);
		strcpy(c[i]->taxon_name[0],a[0]->taxon_name[i]);
		for (j=0;j<c[i]->length;j++) {
			if (isalpha(c[i]->s[0][j])) if (islower(c[i]->s[0][j])) c[i]->s[0][j]=toupper(c[i]->s[0][j]);
			if (c[i]->s[0][j]=='.') c[i]->s[0][j]=c[0]->s[0][j];
			}
		}
	/*establish data types and set*/
	start=(int *)malloc(values->number_of_input_alignments*sizeof(int));
	assert((int) start);
	stop=(int *)malloc(values->number_of_input_alignments*sizeof(int));
	assert((int) stop);
	start[0]=0;
	stop[0]=original_0;
	for (i=1;i<values->number_of_input_alignments;i++) {
		start[i]=stop[i-1];
		stop[i]=start[i]+a[i]->length;
		}

for (i=0;i<values->number_of_input_alignments;i++) if (a[i]) a[i]=dump_align(a[i]);
free(a);

for (i=0;i<values->all_done;i++) {
		for (j=0;j<c[i]->length;j++) {
			for (k=0;k<values->number_of_input_alignments;k++) if ((j>=start[k]) && (j<stop[k])) {
				current_type=k;
				c[i]->s[1][j]=values->data_sets[current_type]+1;
				}
			}
		c[i]->s[1][c[i]->length]='\0';
		c[i]->s[1]=(char *)realloc(c[i]->s[1],(1+strlen(c[i]->s[1]))*sizeof(char));
		assert((int) c[i]->s[1]);
		}
/*set end points*/
for (i=0;i<values->all_done;i++) {
	for (k=0;k<values->number_of_input_alignments;k++) {
		c[i]->s[c[i]->n_seqs][stop[k]-1] += (64);
		}
	}

	/*strip our gapif to be realigned and mark begginings and ends*/
	for (i=0;i<values->all_done;i++) {
		strip_out_gaps_and_recount(c[i]);
		}
	/*dump b's*/

if (b) {
	for (i=0;i<old_num_alignments;i++) {
		if (b[i]) {
			for (j=0;j<num_bs[i];j++) if (b[i][j]) b[i][j]=dump_align(b[i][j]);
			free(b[i]);
			}
		}
	free(b);
	free(num_bs);
	}
free(start);
free(stop);
if (values->number_of_input_alignments==1) 	values->number_of_input_alignments=0;
return c;
}

alignment **new_read_old_file(a,values,num_bs)
alignment **a;
parameters *values;
int *num_bs;
{
	int i,j;

	a=(alignment **)malloc(sizeof(alignment *));
	assert((int)a);

	if (values->VERBOSE) {
		fprintf(stderr,"	");
		fflush(stderr);
	}
	for (i = 0; ((a[i] = (alignment *)read_alignment (values))!=0); i++) {
		if (values->VERBOSE){
			fprintf(stderr,"%d ",i+1);
			fflush(stderr);
			a=(alignment **)realloc(a,(i+1)*sizeof(alignment *));
			assert((int) a);
		}
	}

	a=(alignment **)realloc(a,(i)*sizeof(alignment *));
	assert((int)a);
	for (j=0;j<i;j++) a[j]->type_weight=0;
	if (values->VERBOSE) fprintf(stderr,"\n");
	*num_bs=i;
	return a;

}

alignment *convert_b_to_a(b,num)
alignment **b;
int num;
{
	alignment *a;
	int i,j,longest;


	longest=0;
	for (i=0;i<num;i++) if (b[i]->length>longest) longest=b[i]->length;

	a=(alignment *)malloc(sizeof(alignment));
	assert((int) a);
	a->length=longest;
	a->n_seqs=num;
	a->s=(char **)malloc(num*sizeof(char *));
	assert((int) a->s);
	for (i=0;i<num;i++) {
		a->s[i]=(char *)malloc((a->length+1)*sizeof(char));
		assert ((int) a->s[i]);
		a->s[i][a->length]='\0';
		for (j=0;j<b[i]->length;j++) {
			a->s[i][j]=b[i]->s[0][j];
			if (islower(a->s[i][j])) a->s[i][j]=toupper(a->s[i][j]);
			}
		for (j=b[i]->length;j<a->length;j++) a->s[i][j]='-';
	}
	a->taxon_name=(char **)malloc(num*sizeof(char *));
	assert((int) a->taxon_name);
	for (i=0;i<num;i++) {
		a->taxon_name[i]=(char *)malloc((strlen(b[i]->name)+1)*sizeof(char));
		assert ((int) a->taxon_name[i]);
		a->taxon_name[i]=(char *)strcpy(a->taxon_name[i],b[i]->name);
	}
	a->name=(char *)malloc((strlen(b[0]->name)+1)*sizeof(char));
	assert ((int) a->name);
	a->name=(char *)strcpy(a->name,b[0]->name);
	a->type_weight=0;
	return a;
}

void strip_out_gaps_and_recount(a)
alignment *a;
{
int i,j,k,ii;
int *type, *length;
char *holder;
int old_length, new_length;
int all_gaps;

holder=(char *)malloc(a->length*sizeof(char));
assert((int) holder);
type=(int *)malloc(a->length*sizeof(int));
assert((int) type);

/*if segemnt has no length make sure it has at least one gap if all gaps mae boundry +64 then reduce later*/
all_gaps=1;
for (i=0;i<a->length;i++) {
	if (a->s[0][i]!='-') all_gaps=0;
	if (a->s[1][i]>63) {
		if (all_gaps) a->s[1][i]+=64;
		else all_gaps=1;
		}
	}

old_length=a->length;
k=0;
for (i=0;i<a->length;i++) {
	if (!((a->s[0][i]=='-') && ((a->s[1][i]==1) || (a->s[1][i]==4) || (a->s[1][i]==(1+64)) || (a->s[1][i]==(4+64))))) {
		holder[k]=a->s[0][i];
		type[k]=a->s[1][i];
		k++;
		}
	else if (a->s[1][i]>63) type[k-1] += (64);
	}
new_length=k;
for (i=0;i<new_length;i++) if (type[i]>127) type[i]-=64;

a->length=new_length;
free(a->s[0]);
free(a->s[1]);

a->s[1]=(char *)malloc((1+new_length)*sizeof(char));
assert((int) a->s[1]);
a->s[0]=(char *)malloc((1+new_length)*sizeof(char));
assert((int) a->s[0]);

for (i=0;i<a->length;i++) {
	a->s[0][i]=holder[i];
	a->s[1][i]=type[i];
	}

a->s[0][a->length]='\0';
a->s[1][a->length]='\0';

free(holder);
free(type);
}


alignment *convert_to_alignment(ntaxa,name,sequences,file_name)
int ntaxa;
char **name, **sequences;
char *file_name;
{
int i,j;
alignment *a;
int max_length;

max_length=0;
for (i=0;i<ntaxa;i++) {if (strlen(sequences[i]) > max_length) max_length=strlen(sequences[i]);}
max_length+=1;
a=allocate_alignment_and_strings(ntaxa,max_length);
a->score=0;
a->length=max_length-1;
for (i=0;i< ntaxa; i++) {
	strcpy(a->s[i],sequences[i]);
	for (j=strlen(sequences[i]);j<max_length-1;j++) {a->s[i][j]='-';}
	a->s[i][max_length-1]='\0';
	a->taxon_name[i]=(char *)malloc((strlen(name[i])+1)*sizeof(char));
	assert((int)a->taxon_name[i]);
	strcpy(a->taxon_name[i],name[i]);
	}
a->name=(char *)malloc(LINESIZE*sizeof(char));
assert((int)a->name);
sprintf(a->name,"(%s %s)%c","Alignment from",file_name,'\0');
a->name=(char *)realloc(a->name,(1+strlen(a->name))*sizeof(char));
assert((int)a->name);
a->type_weight=0;
return a;
}

alignment **do_chop_thing(a,values)
alignment **a;
parameters *values;
{
int i,j, k,num_to_chop,num_pieces,total_pieces;
alignment ***b,**c;
int *align_info, *num_in_b, *set_holder, *weight_holder;
int place_counter,multi_counter, old_number;

/*stuff for non-infile entry*/
if (!values->number_of_input_alignments) {
	values->number_of_input_alignments=1;
  values->data_sets=(int *)malloc((values->number_of_input_alignments)*sizeof(int));
  assert((int) values->data_sets);
  values->data_set_weights=(int *)malloc((values->number_of_input_alignments)*sizeof(int));
  assert((int) values->data_set_weights);
  values->data_sets[0]=0;
  values->data_set_weights[0]=1;
	}

align_info=(int *)malloc(FIFTEEN_BITS_MAX*sizeof(int));
assert((int) align_info);
b=NULL;
num_to_chop=0;
for (i=0;i<values->number_of_input_alignments;i++) if (values->data_sets[i]==0) ++num_to_chop;
b=(alignment ***)malloc(num_to_chop*sizeof(alignment **));
assert((int) b);
num_in_b=(int *)malloc(num_to_chop*sizeof(int *));
assert((int) num_in_b);
old_number=values->number_of_input_alignments;

total_pieces=j=0;
for (i=0;i<values->number_of_input_alignments;i++) {
	if (values->data_sets[i]==0) {
		num_pieces=0;
		a[i]->type_weight=0;
		b[j]=get_chopped(a[i],values,&num_pieces,align_info+total_pieces);
		++num_pieces;
		num_in_b[j++]=num_pieces;
		total_pieces+=num_pieces;
		}
	else {
		align_info[total_pieces]=values->data_sets[i];/*changed here*/
		total_pieces+=1;
		}
	}
/*allocate new thangs*/
c=(alignment **)malloc(total_pieces*sizeof(alignment *));
assert((int) c);
set_holder=(int *)malloc(total_pieces*sizeof(int));
assert((int) set_holder);
weight_holder=(int *)malloc(total_pieces*sizeof(int));
assert((int) weight_holder);

/*copy over to c after chopping*/
place_counter=0;
multi_counter=0;
for (i=0;i<values->number_of_input_alignments;i++) {
	if (values->data_sets[i]==0) {
		for (j=0;j<num_in_b[multi_counter];j++) {
			c[place_counter]=make_align(b[multi_counter][j]);
			/*set_holder[place_counter]=values->data_sets[i];*/
			weight_holder[place_counter]=values->data_set_weights[i];
			++place_counter;
			}
		++multi_counter;
		}
	else {
		c[place_counter]=make_align(a[i]);
		/*set_holder[place_counter]=values->data_sets[i];*/
		weight_holder[place_counter]=values->data_set_weights[i];
		++place_counter;
		}
}

/*values->all_done
values->number_of_input
values->data_sets, data_set weights*/

values->data_sets=(int *)realloc(values->data_sets,(total_pieces)*sizeof(int));
assert((int) values->data_sets);
values->data_set_weights=(int *)realloc(values->data_set_weights,(total_pieces)*sizeof(int));
assert((int) values->data_set_weights);
values->number_of_input_alignments=total_pieces;

for (i=0;i<total_pieces;i++) {
	values->data_sets[i]=align_info[i];
	values->data_set_weights[i]=weight_holder[i];
	}
for (i=0;i<num_to_chop;i++) {
	for (j=0;j<num_in_b[i];j++) {
		b[i][j]=dump_align(b[i][j]);
		}
	free(b[i]);
	}
free(b);
free(num_in_b);
free(align_info);
free(set_holder);
free(weight_holder);

/*realloc not dump*/
for (i=0;i<old_number;i++) {
	if (a[i]) a[i]=dump_align(a[i]);
	}
free(a);
a=(alignment **)malloc(total_pieces*sizeof(alignment *));
assert((int) a);

for (i=0;i<total_pieces;i++) {
	a[i]=make_align(c[i]);
	c[i]=dump_align(c[i]);
	}
free(c);

return a;
}

alignment **get_chopped(a,values,num_barriers,align_info)
alignment *a;
int *num_barriers, *align_info;
parameters *values;
{
int i,j,k;
alignment **b,**temp_align;
int nbases,*same;
int **piece_info, size_holder;
int piece_holder,threshold, new_cost,speedup;
int num_forced,chop_counter,base_counter,num_aster;

threshold=values->chop;
/*find contiguous areas*/
num_forced=0;
nbases=a->length;
for (j=0;j<nbases;j++) if (a->s[0][j]=='*') {
	++num_forced;
	for (k=1;k<a->n_seqs;k++) if (a->s[k][j]!='*') {
		fprintf(stderr,"'*' characters do not line up in taxon %s bye.\n",a->taxon_name[k]);
		exit(-1);
		}
	}

same=(int *)malloc(nbases*sizeof(int));
assert((int)same);
for (j=0;j<nbases;j++) {
	if (a->s[0][j]!='-') same[j]=1;
	else if (a->s[0][j]=='*') same[j]=1;
	else same[j]=0;
	for (i=1;i<a->n_seqs;i++) if (a->s[i][j]=='.') a->s[i][j]=a->s[0][j];
	for (i=1;i<a->n_seqs;i++) if (a->s[i][j]!=a->s[0][j]) same[j]=0;
	}

/*get size and number of constant areas*/
(*num_barriers)=0;
for (j=1;j<nbases;j++) if (same[j]!=same[j-1]) ++(*num_barriers);

/*fprintf(stderr,"\nOriginal Number of Segments = %d\n",(*num_barriers)+1);*/
/*piece_info[0]=type [1]=beggining [2]=length*/
piece_info=(int **)malloc(((*num_barriers)+1)*sizeof(int *));
assert((int) piece_info);
for (i=0;i<=(*num_barriers);i++) {
	piece_info[i]=(int *)malloc(3*sizeof(int));
	assert((int)piece_info[i]);
	}
for (i=0;i<=(*num_barriers);i++) piece_info[i][2]=0;
piece_info[0][0]=same[0];
piece_info[0][1]=0;
piece_info[0][2]=1;
piece_holder=0;
for (j=1;j<nbases;j++) {
	if (same[j]!=same[j-1]) ++piece_holder;
	++piece_info[piece_holder][2];
	}
for (i=1;i<=(*num_barriers);i++) piece_info[i][1]=piece_info[i-1][1]+piece_info[i-1][2];
for (i=0;i<=(*num_barriers);i++) piece_info[i][0]=same[piece_info[i][1]];

/*set the too small barriers to 0*/
for (i=0;i<=(*num_barriers);i++) {
  if (piece_info[i][0]==1){
    if (piece_info[i][2]<threshold) for (j=piece_info[i][1];j<(piece_info[i][1]+piece_info[i][2]);j++) same[j]=0;
  }
}
for (j=1;j<nbases;j++) if (a->s[0][j]=='*') same[j]=1;
for (i=0;i<=(*num_barriers);i++) free(piece_info[i]);
free(piece_info);

(*num_barriers)=0;
for (j=1;j<nbases;j++) if (same[j]!=same[j-1]) ++(*num_barriers);
/*fprintf(stderr,"\nNewer Number Segments = %d\n",(*num_barriers)+1);*/

/*return if no chop possible*/
if ((*num_barriers)==0) {
	free(same);
	align_info[0]=0;
	b=(alignment**)malloc((1+(*num_barriers))*sizeof(alignment *));
	assert((int) b);
	for (i=0;i<=(*num_barriers);i++) {
		b[i]=make_align(a);
		b[i]->type_weight=0;
		}
	if (values->number_of_input_alignments==1) {
		if (values->chop<MAX_SEQUENCE_SIZE) fprintf(stderr,"No chopping possible in this segment.  On we go.\n");
		}
	return b;
	}


piece_info=(int **)malloc(((*num_barriers)+1)*sizeof(int *));
assert((int) piece_info);
for (i=0;i<=(*num_barriers);i++) {
	piece_info[i]=(int *)malloc(3*sizeof(int));
	assert((int)piece_info[i]);
	}

for (i=0;i<=(*num_barriers);i++) piece_info[i][2]=1;
piece_info[0][1]=0;
piece_info[0][0]=same[0];
piece_holder=0;
for (j=1;j<nbases;j++) {
	if (same[j]!=same[j-1]) ++piece_holder;
	else ++piece_info[piece_holder][2];
	}

for (i=1;i<=(*num_barriers);i++) piece_info[i][1]=piece_info[i-1][1]+piece_info[i-1][2];
for (i=0;i<=(*num_barriers);i++) piece_info[i][0]=same[piece_info[i][1]];
/*for (i=0;i<=(*num_barriers);i++) fprintf(stderr,"Piece %d type %d begin %d size %d\n",i,piece_info[i][0],piece_info[i][1],piece_info[i][2]);*/

/*for (i=0;i<nbases;i++) {
	for (j=0;j<=(*num_barriers);j++) if (i==piece_info[j][1]) printf("\n");
	printf("%d",same[i]);
	}
printf("\n");
for (i=0;i<=(*num_barriers);i++) printf("Frag %d Type %d Beg %d length %d end %d \n",i,piece_info[i][0],piece_info[i][1],piece_info[i][2],piece_info[i][1]+piece_info[i][2]-1);
*/
new_cost=0;
for  (i=0;i<=(*num_barriers);i++) if (!piece_info[i][0]) new_cost+= (piece_info[i][2]*piece_info[i][2]);
speedup=(100*(nbases*nbases))/new_cost;
fprintf(stderr,"	Potential Speedup of %d%% in this segment\n",speedup);

/*output new alignment files*/
/*output malign batch file*/
/*remember to modify malign to ad a gap if lenght ==0
	remember to add elision option to only produce the printermediate files*/
chop_counter=0;
/*determine the number already included*/
for (i=0;i<=(*num_barriers);i++) if (((piece_info[i][0]==1) && (piece_info[i][2]<threshold))) ++chop_counter;
num_forced=chop_counter;
b=(alignment**)malloc((1+(*num_barriers))*sizeof(alignment *));
assert((int) b);
for (i=0;i<=(*num_barriers);i++) {
	b[i]=make_align(a);
	b[i]->type_weight=0;
	b[i]->length=piece_info[i][2];
	for (j=0;j<b[i]->n_seqs;j++) {
		free(b[i]->s[j]);
		b[i]->s[j]=(char *)malloc((1+b[i]->length)*sizeof(char));
		assert((int) b[i]->s[j]);
		for (k=0;k<b[i]->length;k++) b[i]->s[j][k]=a->s[j][k+piece_info[i][1]];
		b[i]->s[j][b[i]->length]='\0';
		}
	align_info[i]=piece_info[i][0];
	free(piece_info[i]);
	}
/*if length 1 and * remove alignment*/
base_counter=num_aster=0;
for (i=0;i<=(*num_barriers);i++) {
	if (b[i]->length!=1) {
		for (j=0;j<b[i]->length;j++) if (b[i]->s[0][j]=='*') {
			++num_aster;
			for (k=0;k<b[i]->n_seqs;k++) for (base_counter=j;base_counter<b[i]->length;base_counter++) {
				b[i]->s[k][base_counter]=b[i]->s[k][base_counter+1];
				}
			}
		if (num_aster>0) {
			b[i]->length-=num_aster;
			for (k=0;k<b[i]->n_seqs;k++) {
				b[i]->s[k]=(char *)realloc(b[i]->s[k],(1+b[i]->length)*sizeof(char));
				assert((int) b[i]->s[k]);
				b[i]->s[k][b[i]->length]='\0';
				}
			}
		}
	num_aster=0;
	base_counter=0;
	}
temp_align=(alignment **)malloc(((*num_barriers)+1-chop_counter)*sizeof(alignment *));
assert((int)temp_align);
chop_counter=0;
for (i=0;i<=(*num_barriers);i++) {
	if (!((b[i]->length==1) && (b[i]->s[0][0]=='*'))) {
		temp_align[chop_counter]=make_align(b[i]);
		align_info[chop_counter]=align_info[i];
		++chop_counter;
		}
	}

for (i=0;i<=(*num_barriers);i++) 	b[i]=dump_align(b[i]);
b=(alignment **)realloc(b,chop_counter*sizeof(alignment *));
assert((int) b);
for (i=0;i<chop_counter;i++) {
	b[i]=make_align(temp_align[i]);
	temp_align[i]=dump_align(temp_align[i]);
	}
free(temp_align);

/*if other pieces --which are not to be aligned-- have * then remove the character and reallocate etc*/

free(same);
free(piece_info);
(*num_barriers)=(chop_counter-1);
/*realloc align_info to reflect new num pieces*/
return b;
}

void make_all_unordered(ntaxa,sequences,fp_in)
int ntaxa;
char **sequences;
FILE *fp_in;
{
int i, j, l, k,num_states,construct,*ordered, nchar, we_are_at,place;
char *buf, *beg, *end;
int num_ordered,new_num,highest_state,start,stop;
char **comp,*t_seq;
char *temp;
/*fix for polymorphism AND reallocate*/
for (i=0;i<ntaxa;i++) {
    nchar=0;
    temp=(char *)malloc((strlen(sequences[i]))*sizeof(char));
    assert((int) temp);
    for (j=0;j<strlen(sequences[i]);j++) {
        /*if open bracket*/
        if (sequences[i][j]=='[') {
            fprintf(stderr,"Poly");
            num_states=construct=0;
            /*read states until ']' */
                for (k=j+1;k<strlen(sequences[i]);k++) {
                    if (sequences[i][k]==']') break;
                    else {
                        construct+=make_bit(sequences[i][k]);
                        num_states+=1;
                    }
                }
            /*put in bitrep*/
            temp[nchar++]=(128 + construct);
            j=k;
        }
        else temp[nchar++]=sequences[i][j];
    }
    /*reallocate and rewritw*/
    sequences[i]=realloc(sequences[i],(1+nchar)*sizeof(char));
    assert((int) sequences[i]);
    for (j=0;j<nchar;j++) sequences[i][j]=temp[j];
    sequences[i][nchar]='\0';
    free(temp);
}

nchar = strlen(sequences[0]);
buf=(char *)malloc(100*sizeof(char *));
assert((int) buf);
beg=(char *)malloc(100*sizeof(char *));
assert((int) beg);
end=(char *)malloc(100*sizeof(char *));
assert((int) end);
ordered=(int *)malloc(nchar*sizeof(int));
assert((int) ordered);
for (i=0;i<nchar;i++) ordered[i]=0;
/*read coding statements*/
	/*read till cc or eof but for now just the unordered multistate*/
while (!feof(fp_in)) {
	fscanf(fp_in,"%s ",buf);
	place=got_dash_or_dot(buf);
	if (place<0) ordered[atoi(buf)]=1; /*fix for begginig and end*/
	else {
		if (place==0) start=0;
		else {
			for (i=0;i<place;i++) beg[i]=buf[i];
			beg[place]='\0';
			start=atoi(beg);
			}
		if (place==(strlen(buf)-1)) stop=nchar-1;
		else {
			for (i=place+1;i<strlen(buf);i++) end[i-(place+1)]=buf[i];
			end[strlen(buf)-(place+1)]='\0';
			stop=atoi(end);
			}
		for (i=start;i<=stop;i++) ordered[i]=1;
		}
	}

/*recode those characters*/
num_ordered=new_num=0; /*new_num id the new number of extra chars 0/1/2=1 0-3=2 etc.*/
for (i=0;i<nchar;i++) if (ordered[i]) {
	++num_ordered;
	highest_state=0;
	for (j=0;j<ntaxa;j++) if ((sequences[j][i]!='?') && (sequences[j][i]!='-')) if ((sequences[j][i]-('0'))>highest_state) highest_state=(sequences[j][i]-('0'));
	if (highest_state < 2) {
		ordered[i]=0;
		--num_ordered;
		}
	else new_num+=highest_state;
	if (highest_state > 6 ) {
		fprintf(stderr,"Character data must in the range {0-6}\n	Recode or delete the offending data. (OK if ordered)\n");
		/*exit(-1);*/
		}
	}

/*make sure ther are some to recode*/
if (new_num>0) {
	comp=(char **)malloc(ntaxa*sizeof(char *));
	assert((int) comp);
	for (i=0;i<ntaxa;i++) {
		comp[i]=(char *)malloc(new_num*sizeof(char));
		assert((int) comp[i]);
		}

	we_are_at=0;
	for (i=0;i<nchar;i++) if (ordered[i]) {
		highest_state=0;
		for (j=0;j<ntaxa;j++) if ((sequences[j][i]!='?') && (sequences[j][i]!='-')) if ((sequences[j][i]-('0'))>highest_state) highest_state=(sequences[j][i]-('0'));
		for (l=0;l<ntaxa;l++) {
			for (j=0;j<highest_state;j++) {
				if ((sequences[l][i]=='?') || (sequences[l][i]=='-')) comp[l][we_are_at+j]='?';
				else if ((sequences[l][i]-('0'))>j) comp[l][we_are_at+j]=1+('0');
				else comp[l][we_are_at+j]=0+('0');
				}
			}
		we_are_at+=highest_state;
		}
	/*add to end of sequences by reallocation and block out the originals*/
	t_seq=(char *)malloc((1+nchar)*sizeof(char));
	assert((int) t_seq);
	for (i=0;i<ntaxa;i++) {
		/*copy to temp*/
		t_seq=strcpy(t_seq,sequences[i]);
		free(sequences[i]);
		sequences[i]=(char *)malloc((1+nchar-num_ordered+new_num)*sizeof(char));
		assert((int) sequences[i]);
		l=0;
		for (j=0;j<nchar;j++) if (!ordered[j]) sequences[i][l++]=t_seq[j];
		for (j=0;j<new_num;j++) sequences[i][l++]=comp[i][j];
		sequences[i][nchar-num_ordered+new_num]='\0';

		}
	free(t_seq);
	for (i=0;i<ntaxa;i++) free(comp[i]);
	free(comp);
	fprintf(stderr,"	%d ordered multi-state characters were recoded as %d binary characters.\n",num_ordered,new_num);
}
/*free*/
free(ordered);
free(buf);
free(beg);
free(end);
}

int got_dash_or_dot(s)
char *s;
{
int i;
for (i=0;i<strlen(s);i++) if ((s[i]=='.') || (s[i]=='-')) return i;
return -1;
}

int make_bit(a)
char a;
{
    if (a=='0') return 1;
    if (a=='1') return 2;
    if (a=='2') return 4;
    if (a=='3') return 8;
    if (a=='4') return 16;
    if (a=='5') return 32;
    if (a=='6') return 64;
    if (a=='7') return 128;
    return 0;
}