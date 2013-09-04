/*Copyright 1993 Ward Wheeler all rights reserved*/

#include "align3.h"

alignment **read_old_file(a,values)
alignment **a;
parameters *values;
{
int i,j;

a=(alignment **)malloc(1*sizeof(alignment *));
assert((int)a);
a[0]=NULL;

if (values->VERBOSE) {
	fprintf(stderr,"Reading ");
	fflush(stderr);
	}
for (i = 0; ((a[i] = (alignment *)read_alignment (values))!=0); i++) {
	if (values->VERBOSE){ fprintf(stderr,"%d ",i+1); fflush(stderr);}
	a=(alignment **)realloc(a,(i+2)*sizeof(alignment *));
	assert((int)a);
	a[i]->type_weight=0;
	a[i+1]=NULL;
	}
values->all_done = i;
a=(alignment **)realloc(a,((2*(values->all_done+1))-1)*sizeof(alignment *));
assert((int)a);
for (j=values->all_done;j<((2*(values->all_done+1))-1);j++) a[j]=NULL;
if (values->VERBOSE) fprintf(stderr,"\n");

return a;
}

void do_pre_stuff(values,a)
parameters *values;
alignment **a;
{
 int i,j,k;

	fix_all_after_values(a,values);

	if ((values->how_many>0) || (values->randtrees>0)) {
		values->rand_seed=(long int) time(NULL);
		srand((int) values->rand_seed);
	 	}

	if ((values->VERBOSE) && (!values->elision) && (values->align_multi_swap==(values->all_done-2))) fprintf(stderr,"This will be an exhaustive alignment search.\n");
	if ((values->VERBOSE) && (!values->elision) && (values->tree_multi_swap==(values->all_done-3))) fprintf(stderr,"This will use an exhaustive cladogram search.\n");
	if (values->all_done<2) {
		fprintf(stderr,"What are you trying to do with only a single sequence?\n");
		exit(-1);
		}
	if ((!values->elision) && (values->all_done==2)) {
		values->aquick=values->get_heur2=values->best=values->rand_align=0;
		fprintf(stderr,"When the input file contains two sequences or alignments only the option 'pair' can be used.\n");
		}
	if (values->all_done<4) values->phylo_score=0;
	/*Wget real number of sequences*/
	for (i=0;i<values->all_done;i++) values->actual_num_sequences+=a[i]->n_seqs;
	if (values->VERBOSE) if (values->number_of_input_alignments>0) fprintf(stderr,"Entered %d sequences from %d alignments\n",values->actual_num_sequences,values->number_of_input_alignments);
	/*for (i=0;i<values->all_done;i++) for (j=0;j<a[i]->n_seqs;j++) for (k=0;k<a[i]->length;k++) if (a[i]->s[j][k]=='X') a[i]->s[j][k]='N'; */

}

void end_free(values)
parameters *values;
{
int i,j, classes;

if (values->groups) {
	for (i=0;i<values->ngroups;i++) if (values->groups[i]) free(values->groups[i]);
	free(values->groups);
	}
if (values->align_file) free(values->align_file);

if (values->delta) {
	for (i=0;i<4;i++) if (values->delta[i]) free(values->delta[i]);
	free(values->delta);
	}
if (values->new_codes) {
	for (i=0;i<values->n_codes;i++ ) if (values->new_codes[i]) free(values->new_codes[i]);
	free(values->new_codes);
	}
if (values->input_names) {
	for (i=0;i<values->all_done;i++) if (values->input_names[i]) free(values->input_names[i]);
	free(values->input_names);
	}
if (values->best_rep) {
	for (i=0;i<values->keep_trees;i++) if (values->best_rep[i]) free(values->best_rep[i]);
	free(values->best_rep);
	}
if (values->ce_weights) free(values->ce_weights);

classes=values->all_done-1;
if ((values->cache_info) && (values->align_cache) && (values->VERBOSE)) {
	fprintf(stderr,"Alignment cache :\n");
	/*if (values->check_cache_for_scores) values->all_done+=1;*/
	for (i=0;i<classes;i++) fprintf(stderr,"	Size[%d]=> %d hits and %d misses -- %d%% yield\n",i+2,values->align_cache[i]->hit,values->align_cache[i]->miss,100*values->align_cache[i]->hit/(values->align_cache[i]->hit+values->align_cache[i]->miss));
	for (i=1;i<classes;i++) {
  	if (!values->new_optimization) values->align_cache[0]->hit+=((i+1)*values->align_cache[i]->hit);
  	else values->align_cache[0]->hit+=(values->align_cache[i]->hit);
    values->align_cache[0]->miss+=values->align_cache[i]->miss;
    }
  fprintf(stderr,"	Total (weighted) => %d hits and %d misses -- %d%% yield\n",values->align_cache[0]->hit,values->align_cache[0]->miss,100*values->align_cache[0]->hit/(values->align_cache[0]->hit+values->align_cache[0]->miss));
 	/*if (values->check_cache_for_scores) values->all_done-=1; */
  }

if (values->align_cache) {
	if (values->VERBOSE) fprintf(stderr,"Freeing cache...");
	for (j=0;j<classes;j++) {
	 	for (i=0;i<values->align_cache[j]->max_num;i++) {
    			if (values->align_cache[j]) {
	   			if (values->align_cache[j]->ba[i]) values->align_cache[j]->ba[i]=dump_align(values->align_cache[j]->ba[i]);
        		}
	   	}
	   free(values->align_cache[j]->age);
	   free(values->align_cache[j]->hack);
	   free(values->align_cache[j]->ba);
	   free(values->align_cache[j]);
	   }
  free(values->align_cache);
  if (values->VERBOSE) fprintf(stderr,"Done.\n");
  /*if (values->check_cache_for_scores) values->all_done-=1;*/
 	}
 if (values->other_parm) {
 	for (i=0;i<values->number_of_input_alignments;i++) if (values->other_parm[i] && (values->other_parm[i]!=values)) end_free(values->other_parm[i]);
 	free(values->other_parm);
	}
if (values->jack_array) free(values->jack_array);

}


void do_estimated_weights(best_aligns,values,count)
alignment **best_aligns;
parameters *values;
int count;
{
int i;
int **combinatorial_weight,w_ti,w_tv,w_gap;

  if (values->saw_wheeler) {
  	combinatorial_weight=get_combinatorial_weights(best_aligns,values,count);
    w_gap=0;
    for (i=0;i<4;i++) w_gap+=combinatorial_weight[i][4];
    w_gap/=4;
    if (values->ttr) {
    	w_ti=(combinatorial_weight[0][2]+combinatorial_weight[1][3])/2;
      w_tv=(combinatorial_weight[0][1]+combinatorial_weight[0][3]+combinatorial_weight[1][2]+combinatorial_weight[2][3])/4;
    	}
    if ((!values->ttr) && (!values->delta)) {
    	w_ti=(combinatorial_weight[0][2]+combinatorial_weight[1][3]+combinatorial_weight[0][1]+combinatorial_weight[0][3]+combinatorial_weight[1][2]+combinatorial_weight[2][3])/6;
    	}
    }
	for (i=0;i< count;i++) {
    if (values->get_weights_2) {
    	values->get_weights=1;
  		cladogram(best_aligns[i],values);
      if (!values->ttr) {
      	values->overall_changes+=values->temp_changes;
      	values->overall_changes_g+=values->temp_changes_g;
        }
      else {
      	values->overall_ti+=values->temp_ti;
      	values->overall_tv+=values->temp_tv;
      	values->overall_ti_g+=values->temp_ti_g;
      	values->overall_tv_g+=values->temp_tv_g;
      	}
      values->overall_gaps+=values->temp_gaps;
      values->overall_gaps_g+=values->temp_gaps_g;
      values->get_weights=0;
    	}
     if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
    }

   /*print weights*/
   if (values->get_weights_2) {
   	if (values->delta) {
    	fprintf(stderr,"\n");
    	if (values->saw_farris) fprintf(stderr,"Farris estimations cannot be performed on matrix characters.\n");
    	if (values->saw_farris) fprintf(stderr,"Goloboff estimations cannot be performed on matrix characters.\n");
    	}
    else {
	   	printf("Relative weights implied by the alignment(s):\n");
	   	fprintf(stderr,"Relative weights implied by the alignment(s):\n");
	    if (values->saw_farris) {
	    	fprintf(stderr,"	Farris [range*ci*ri] (Hennig86, 1988):\n");
	    	printf("	Farris [range*ci*ri] (Hennig86, 1988):\n");
		   	if (!values->ttr) {
		    	printf("		changes		%d\n",values->overall_changes/count);
		    	fprintf(stderr,"		changes		%d\n",values->overall_changes/count);
		      }
		     else {
		    	printf("		transitions	%d\n",values->overall_ti/count);
		    	fprintf(stderr,"		transitions	%d\n",values->overall_ti/count);
		    	printf("		transversions	%d\n",((values->overall_tv)/count));
		    	fprintf(stderr,"		transversions	%d\n",((values->overall_tv)/count));
		     	}
		  	 printf("		gaps		%d\n",values->overall_gaps/count);
		  	 fprintf(stderr,"		gaps		%d\n",values->overall_gaps/count);
	       fprintf(stderr,"\n");
	       printf("\n");
		     }
	    if (values->saw_goloboff) {
	    	fprintf(stderr,"	Goloboff [range*k/(k+extra steps) for k=%d] (1993):\n",values->saw_goloboff);
	    	printf("	Goloboff [range*k/(k+extra steps)] (1993):\n");
		   	if (!values->ttr) {
		    	printf("		changes		%d\n",values->overall_changes_g/count);
		    	fprintf(stderr,"		changes		%d\n",values->overall_changes_g/count);
		      }
		    else {
		    	printf("		transitions	%d\n",values->overall_ti_g/count);
		    	fprintf(stderr,"		transitions	%d\n",values->overall_ti_g/count);
		    	printf("		transversions	%d\n",((values->overall_tv_g)/count));
		    	fprintf(stderr,"		transversions	%d\n",((values->overall_tv_g)/count));
		     	}
		  	printf("		gaps		%d\n",values->overall_gaps_g/count);
		  	fprintf(stderr,"		gaps		%d\n",values->overall_gaps_g/count);
	      fprintf(stderr,"\n");
	      printf("\n");
		    }
     	}
	  if (values->saw_wheeler) {
    	if (values->delta) {
		   	printf("Relative weights implied by the alignment(s):\n");
		   	fprintf(stderr,"Relative weights implied by the alignment(s):\n");
      	}
	  	fprintf(stderr,"	Wheeler combinatorial [rescaled and symmetrical (i->j=j->i)] (1990):\n");
	  	printf("	Wheeler combinatorial [rescaled and symmetrical (i->j=j->i)] (1990):\n");
      /*print weights*/
      if (values->delta) {
      	fprintf(stderr,"		    %3s %3s %3s %3s\n","A","C","G","T");
        fprintf(stderr,"		%3s ","A");
        for (i=0;i<4;i++) fprintf(stderr,"%3d ",combinatorial_weight[0][i]);
        fprintf(stderr,"\n");
        fprintf(stderr,"		%3s ","C");
        for (i=0;i<4;i++) fprintf(stderr,"%3d ",combinatorial_weight[1][i]);
        fprintf(stderr,"\n");
        fprintf(stderr,"		%3s ","G");
        for (i=0;i<4;i++) fprintf(stderr,"%3d ",combinatorial_weight[2][i]);
        fprintf(stderr,"\n");
        fprintf(stderr,"		%3s ","T");
        for (i=0;i<4;i++) fprintf(stderr,"%3d ",combinatorial_weight[3][i]);
        fprintf(stderr,"\n");

      	printf("		    %3s %3s %3s %3s\n","A","C","G","T");
        printf("		%3s ","A");
        for (i=0;i<4;i++) printf("%3d ",combinatorial_weight[0][i]);
        printf("\n");
        printf("		%3s ","C");
        for (i=0;i<4;i++) printf("%3d ",combinatorial_weight[1][i]);
        printf("\n");
        printf("		%3s ","G");
        for (i=0;i<4;i++) printf("%3d ",combinatorial_weight[2][i]);
        printf("\n");
        printf("		%3s ","T");
        for (i=0;i<4;i++) printf("%3d ",combinatorial_weight[3][i]);
        printf("\n");
      	}
      else if (values->ttr) {
		    printf("		transitions	%d\n",w_ti);
		    fprintf(stderr,"		transitions	%d\n",w_ti);
		    printf("		transversions	%d\n",w_tv);
		    fprintf(stderr,"		transversions	%d\n",w_tv);
      	}
      else {
		   	printf("		changes		%d\n",w_ti);
		   	fprintf(stderr,"		changes		%d\n",w_ti);
      	}
		  printf("		gaps		%d\n",w_gap);
		  fprintf(stderr,"		gaps		%d\n",w_gap);
	    fprintf(stderr,"\n");
	    printf("\n");
      /*free weights*/
      for (i=0;i<5;i++) free(combinatorial_weight[i]);
      free(combinatorial_weight);
	    fprintf(stderr,"\n");
	    printf("\n");
	    }
    }
	if (best_aligns) free(best_aligns);
	if (values->previous) values->previous=dump_align(values->previous);
}

void do_modifications(a,values)
alignment **a;
parameters *values;
{
int i,j,is_protein;

/*make sure all are upper case*/
for (j=0;j<values->all_done;j++) for (i=0;i<a[j]->length;i++) if (islower(a[j]->s[0][i])) a[j]->s[0][i]=toupper(a[j]->s[0][i]);

/*protein to nucleotide conversion here*/
for (j=0;j<values->all_done;j++){
	is_protein=0;
	if (a[j]->name[strlen(a[j]->name)-1]=='*') {
		is_protein=1;
		a[j]->name[strlen(a[j]->name)-1]='\0';
		}
	/*else {
		for (i=0;i<a[j]->length;i++) if ((a[j]->s[0][i]=='Q')||(a[j]->s[0][i]=='E')||(a[j]->s[0][i]=='I')||(a[j]->s[0][i]=='L')||(a[j]->s[0][i]=='F')||(a[j]->s[0][i]=='P')) is_protein=1;
		}*/
	if (is_protein) {
		if (values->VERBOSE) fprintf(stderr,"\nSequence %s appears to contain amino acids \nand is being converted to nucleotide ambiguities.\n",a[j]->name);
		protein_to_ambiguity(a[j],values);
		}
  }

/*set up taxon names*/
for (j=0;j<values->all_done;j++){
  if (a[j]->name[strlen(a[j]->name)-1]=='!') 				a[j]=reverse_complement_and_rename(a[j],values);
  else if (a[j]->name[strlen(a[j]->name)-1]=='=') 	a[j]=reverse_and_rename(a[j],values);
  else if (a[j]->name[strlen(a[j]->name)-1]=='+') 	a[j]=complement_and_rename(a[j],values);
	}

/*set up taxon names*/
for (j=0;j<values->all_done;j++){
	if (a[j]->taxon_name[0]) free(a[j]->taxon_name[0]);
	a[j]->taxon_name[0]=(char *)malloc((strlen(a[j]->name)+1)*sizeof(char));
	assert((int)a[j]->taxon_name[0]);
	strcpy(a[j]->taxon_name[0],a[j]->name);
	}
}

alignment **get_prealigned_sequences(values)
parameters *values;
{
char buffer[LINESIZE],**sequences;
char **name;
int i,j,ntaxa,nlines,k,best_length;
char temp[LINESIZE];
FILE *fp_in, *fopen();
alignment **a;

a=NULL;
if (!a) {
	a=(alignment **)malloc(1*sizeof(alignment *));
	assert((int)a);
	a[0]=NULL;
	}

values->all_done = values->number_of_input_alignments;
a=(alignment **)realloc(a,((2*(values->all_done+1))-1)*sizeof(alignment *));
assert((int)a);
for (j=0;j<((2*(values->all_done+1))-1);j++) a[j]=NULL;
/*get shit here*/
for (k=0;k<values->number_of_input_alignments;k++) {
	fp_in=(FILE *)fopen(values->input_file_name[k],"r");
	if (!fp_in) {
		fprintf(stderr,"File %s not found\n",values->input_file_name[k]);
		exit(-1);
		}
	else if (values->VERBOSE) fprintf(stderr,"Reading sequences from file %s\n",values->input_file_name[k]);
	fscanf(fp_in,"%d %d",&ntaxa,&nlines);
	if ((ntaxa<1) || (nlines<1)) {
		if (values->rep_error) fprintf(stderr,"Problems in numbers of taxa and/or lines in alignment file %s\n",values->input_file_name[k]);
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
			if (strcmp(temp,name[i])) fprintf(stderr,"Name mismatch -- could be error in file %s.\n",values->input_file_name[k]);
			sequences[i]=(char *)realloc(sequences[i],strlen(sequences[i])+strlen(buffer)+1);
			assert((int)sequences[i]);
			sequences[i]=(char *)strcat(sequences[i],buffer);
			}
		}
	/*copy to a[k]*/
	a[k]=(alignment *)malloc(sizeof(alignment));
	assert((int)a[k]);
	a[k]->s=(char **)malloc(ntaxa*sizeof(char *));
	assert((int)a[k]->s);
	a[k]->taxon_name=(char **)malloc(ntaxa*sizeof(char *));
	assert((int)a[k]->taxon_name);
	for (i=0;i< ntaxa; i++) {
		a[k]->s[i]=(char *)malloc((strlen(sequences[i])+1)*sizeof(char));
		assert((int)a[k]->s[i]);
		strcpy(a[k]->s[i],sequences[i]);
		a[k]->taxon_name[i]=(char *)malloc((strlen(name[i])+1)*sizeof(char));
		assert((int)a[k]->taxon_name[i]);
		strcpy(a[k]->taxon_name[i],name[i]);
		/*printf("%s %s\n",name[i],a[k]->s[i]);*/
		a[k]->name=(char *)malloc(LINESIZE*sizeof(char));
		assert((int)a[k]->name);
		sprintf(a[k]->name,"(%s %s)","Alignment from",values->input_file_name[k]);
		a[k]->name=(char *)realloc(a[k]->name,(1+strlen(a[k]->name))*sizeof(char));
		assert((int)a[k]->name);
		}
	a[k]->score=0;
	a[k]->n_seqs=ntaxa;
	a[k]->length=strlen(a[k]->s[0]);
	a[k]->type_weight=0;
	if (values->VERBOSE) fprintf(stderr,"%d taxa and %d characters entered from file %s\n",ntaxa,a[k]->length,values->input_file_name[k]);
	/*free temp stuff*/
	for (i=0;i<ntaxa;i++) {
		free(name[i]);
		free(sequences[i]);
		}
	free(name);
	free(sequences);
	fclose(fp_in);
	}
if (values->VERBOSE) fprintf(stderr,"\n");

if (values->new_optimization) {
	/*elide and reset values->all_done etc*/
	}

return a;
}

