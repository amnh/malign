/*
heuristic (better) for getting alignments
*/

#include "align3.h"

alignment **do_some(a,ntax,how_many_rand_align,bound,best_aligns,values,count)
alignment **a,**best_aligns;
int ntax,how_many_rand_align,bound,*count;
parameters *values;
{
int i,l,best_tree;
int **nodes;
int ntaxa;

if (values->groups) bound=HUGE_COST;
if (values->VERBOSE) fprintf(stderr,"\nGenerating random alignments...\n");
ntaxa=ntax+1; /*to fool into doing all rootings keeping '0' the root just never do the last alignment*/

nodes=(int **)malloc((ntaxa-1)*sizeof(int *));
assert((int)nodes);
for (i=0;i<(ntaxa-1);i++) {
	nodes[i]=(int *)malloc(2*sizeof(int));
	assert((int)nodes[i]);
	}

nodes[0][0]=0;
nodes[0][1]=ntaxa+1;
nodes[1][0]=1;
nodes[1][1]=2;

l=0;
/*change this*/
best_aligns=get_random(a,ntaxa,nodes,how_many_rand_align,best_aligns,&l,bound,values);
best_tree=best_aligns[0]->score;

if (values->VERBOSE) fprintf(stderr,"Found a total of %d alignments",l);
if (best_tree!=HUGE_COST) if (values->VERBOSE) fprintf(stderr," of cost %d\n",best_aligns[0]->score);
else if (values->VERBOSE) fprintf(stderr,"\n");
for (i=0;i<(ntaxa-1);i++) {
	free(nodes[i]);
	}
free(nodes);
if ((*count)!=0) (*count)=l;
else (*count)=1;
return best_aligns;
}/*end do_all*/


alignment **get_random(a,ntax,nodes,how_many_rand_align,best_aligns,l,bound,values)
alignment **a,**best_aligns;
int ntax,*l,bound;
int how_many_rand_align;
int **nodes;
parameters *values;
{
int i,j,k,z,m,cur_tree;
int temp,temp1,temp2,cur_b_tree,counter;
int **best_nodes;
char **from,overflow=0;
int *temp_rep,**tree_rep,**old_rep,keep_aligns_holder;
alignment *old_best_align;
realloc_align_holder *align_passer;

align_passer=(realloc_align_holder *)malloc(sizeof(realloc_align_holder));
assert((int)align_passer);
best_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int)best_nodes);
for (i=0;i<(ntax-1);i++) {
	best_nodes[i]=(int *)malloc(2*sizeof(int));
	assert((int)best_nodes[i]);
	}

tree_rep=(int **)malloc(values->keep_aligns*sizeof(int *));
	assert((int)tree_rep);
old_rep=(int **)malloc(values->keep_aligns*sizeof(int *));
	assert((int)old_rep);
for (i=0;i<values->keep_aligns;i++) {
	tree_rep[i]=(int *)malloc((ntax-3)*sizeof(int));
	assert((int)tree_rep[i]);
	old_rep[i]=(int *)malloc((ntax-3)*sizeof(int));
	assert((int)old_rep[i]);
	}
temp_rep=(int *)malloc((ntax-3)*sizeof(int));
assert((int)temp_rep);
cur_b_tree=HUGE_COST;
if (best_aligns[0]) old_best_align=make_align(best_aligns[0]);
else old_best_align=NULL;


(*l)=0;
/*get random trees and lengths could change to only increment when actually makes one*/
for (k=0;k<how_many_rand_align;){
	for (i=0;i<=ntax-4;i++){
		temp1=rand()/(2*(i+4));
		temp2=max_rand/(2*(i+4));
		temp_rep[i]=(int) ((((2*(i+4))-5)*temp1)/temp2);
		if (temp_rep[i]!=((2*(i+4))-5)) ++temp_rep[i];
		}
	cur_tree=all_make_nodes(a,temp_rep,nodes,ntax,ntax,values);
	if (cur_tree==HUGE_COST) goto end_of_it_all;
	++k;
	if (cur_tree < cur_b_tree) {
		if (values->VERBOSE) fprintf(stderr,"Found better...");
		cur_b_tree = cur_tree;
		if (values->VERBOSE) fprintf(stderr,"cost %d\n",cur_b_tree);
		all_other_make_nodes(temp_rep,nodes,ntax,best_nodes);
		if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
		best_aligns[0]=make_align(a[best_nodes[0][1]]);
		/*free the rest*/
		for (j=1;j<(*l);j++) best_aligns[j]=dump_align(best_aligns[j]);
		for (j=0;j<=ntax-4;j++) *(tree_rep[0]+j)=temp_rep[j];
		*l=1;
		}
	else if (cur_tree == cur_b_tree) {
		all_other_make_nodes(temp_rep,nodes,ntax,best_nodes);
		temp=1;
		for (j=0;j<(*l);j++) temp*=compare_aligns(a[best_nodes[0][1]],best_aligns[j]);
		if ((temp!=0) || ((*l)==0)) {
			if (values->VERBOSE) fprintf(stderr,"Found another\n");
			if ((*l)==values->keep_aligns) {
      	if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment *)))!=NULL) && ((tree_rep=(int **)realloc(tree_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((old_rep=(int **)realloc(old_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL)) {
        	++values->keep_aligns;
					best_aligns[(*l)]=make_align(a[best_nodes[0][1]]);
          tree_rep[*l]=(int *)malloc((ntax-3)*sizeof(int));
          assert((int)tree_rep[*l]);
          old_rep[*l]=(int *)malloc((ntax-3)*sizeof(int));
          assert((int)old_rep[*l]);
					for (j=0;j<(ntax-3);j++) tree_rep[*l][j]=temp_rep[j];
					++(*l);
        	}
				else if ((values->rep_error) && (overflow==0)){
					fprintf(stderr,"Overflow in rand aligns -- but that's OK.\n");
					overflow=1;
          }
				}
			else {
				best_aligns[(*l)]=make_align(a[best_nodes[0][1]]);
				for (j=0;j<(ntax-3);j++) tree_rep[*l][j]=temp_rep[j];
				++(*l);
				}
			}
		}/*if else stuff*/
	end_of_it_all:;
	}/*k*/

counter=(*l)-1;
for (i=0;i<=counter;i++) {
	for (j=0;j<(ntax-3);j++) old_rep[i][j]=tree_rep[i][j];
	}

if (values->align_multi_swap>1) {
	keep_aligns_holder=values->keep_aligns;
	align_passer=deep_swap_align(nodes,best_nodes,a,counter+1,old_rep,best_aligns,l,ntax,values->align_multi_swap,values,align_passer);
  best_aligns=align_passer->best;
  old_rep=align_passer->rep;
  if (values->keep_aligns>keep_aligns_holder) {
  	tree_rep=(int **)realloc(tree_rep,values->keep_aligns*sizeof(int *));
    for (j=keep_aligns_holder-1;j<values->keep_aligns;j++) {
    	tree_rep[j]=(int *)malloc((ntax-3)*sizeof(int));
      assert((int)tree_rep[j]);
    	}
  	}
  }
else if ((values->align_swap!=0)|| (values->align_multi_swap==1)) {
	from=(char **)malloc(values->keep_aligns*sizeof(char *));
	assert((int)from);
	for (i=0;i<values->keep_aligns;i++) {
		from[i]=(char *)malloc((ntax-3)*sizeof(char));
		assert((int)from[i]);
		for (j=0;j<(ntax-3);j++) from[i][j]=0;
		}
	if (values->VERBOSE) fprintf(stderr,"Random alignment yeilded %d alignment(s) of cost %d to be swapped\n",counter+1,best_aligns[0]->score);
	/*change to loop through all tree_reps*/
	for (z=0;z<=counter;z++) {
		if (values->VERBOSE) fprintf(stderr,"Swapping on alignment %d of %d\n",z+1,counter+1);
		for (i=0;i<ntax-3;i++) temp_rep[i]=(*(old_rep[z]+i));
		for (i=0;i<ntax-3;i++) {
		if (values->VERBOSE) fprintf(stderr,"	swapping sequence %d\n",i+3);
		if ((*(from[z]+i)==0)) {
			for (j=1;j<=((2*i)+3);j++) {
				if (j!=(*(old_rep[z]+i))) {
					temp_rep[i]=j;
					cur_tree=all_make_nodes(a,temp_rep,nodes,ntax,ntax,values);
					if (cur_tree<cur_b_tree) {
						if (values->VERBOSE) fprintf(stderr,"Found better at cost %d\n",cur_tree);
						cur_b_tree=cur_tree;
						if (counter==(values->keep_aligns-1)) {
             	if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment *)))!=NULL) && ((tree_rep=(int **)realloc(tree_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((old_rep=(int **)realloc(old_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((from=(char **)realloc(from,(values->keep_aligns+1)*sizeof(char *)))!=NULL)) {
        				++values->keep_aligns;
			          	best_aligns[values->keep_aligns-1]=NULL;
								tree_rep[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
			          assert((int)tree_rep[values->keep_aligns-1]);
			          old_rep[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
			          assert((int)old_rep[values->keep_aligns-1]);
			          from[values->keep_aligns-1]=(char *)malloc((ntax-3)*sizeof(char));
			          assert((int)from[values->keep_aligns-1]);
                for (k=0;k<(ntax-3);k++) from[values->keep_aligns-1][k]=0;
			        	}
							else if ((values->rep_error) && (overflow==0)){
								fprintf(stderr,"Overflow in rand aligns -- but that's OK.\n");
								overflow=1;
								--counter;
			          }
							for (k=0;k<ntax-3;k++) *(from[counter]+k)=0;
							}
						++counter;
						from[counter][i]=1;
						for (k=0;k<ntax-3;k++) *(old_rep[counter]+k)=temp_rep[k];
						/*keep alignment*/
						if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
						all_other_make_nodes(temp_rep,nodes,ntax,best_nodes);
						best_aligns[0]=make_align(a[best_nodes[0][1]]);
						/*free the rest*/
						for (m=1;m<(*l);m++) if (best_aligns[m]) best_aligns[m]=dump_align(best_aligns[m]);
						(*l)=1;
						}
					else if ((cur_tree==cur_b_tree) && (cur_tree!=HUGE_COST)) {
						if (values->aquick==0) {
							/*keep alignment*/
							all_other_make_nodes(temp_rep,nodes,ntax,best_nodes);
							temp=1;
							for (m=0;m<(*l);m++) temp*=compare_aligns(a[best_nodes[0][1]],best_aligns[m]);
							if ((temp!=0) || ((*l)==0)) {
              	if (((*l)==values->keep_aligns)||(counter==(values->keep_aligns-1))) {
             			if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment *)))!=NULL) && ((tree_rep=(int **)realloc(tree_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((old_rep=(int **)realloc(old_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((from=(char **)realloc(from,(values->keep_aligns+1)*sizeof(char *)))!=NULL)) {
        						++values->keep_aligns;
			          			best_aligns[values->keep_aligns-1]=NULL;
										tree_rep[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
			          		assert((int)tree_rep[values->keep_aligns-1]);
			          		old_rep[values->keep_aligns-1]=(int *)malloc((ntax-3)*sizeof(int));
			          		assert((int)old_rep[values->keep_aligns-1]);
			          		from[values->keep_aligns-1]=(char *)malloc((ntax-3)*sizeof(char));
			          		assert((int)from[values->keep_aligns-1]);
                    for (k=0;k<ntax-3;k++) from[values->keep_aligns-1][k]=0;
                		}
                  }
								if ((*l)==values->keep_aligns) {
									if (overflow==0) {
										if (values->rep_error) fprintf(stderr,"Overflow in rand aligns -- but that's OK.\n");
										overflow=1;
										}
									}
								else {
									if ((best_aligns[0]) && ((*l)==0)) best_aligns[0]=dump_align(best_aligns[0]);
									if (values->VERBOSE) fprintf(stderr,"Found another for %d alignments\n",(*l)+1);
									best_aligns[(*l)]=make_align(a[best_nodes[0][1]]);
									++(*l);
									if (counter==(values->keep_aligns-1)) {
										if (overflow==0) {
											if (values->rep_error) fprintf(stderr,"Overflow in rand aligns -- but that's OK.\n");
											overflow=1;
											}
										}
									else {
										++counter;
										from[counter][i]=1;
										for (k=0;k<ntax-3;k++) old_rep[counter][k]=temp_rep[k];
										}
									}/*else*/
								}/*if temp*/
							}/*quick*/
						}/*else if*/
					}/*j!=*/
				}/*j*/
			 temp_rep[i]=old_rep[z][i];
			 }/*check swap from*/
			}/*i*/
		}/*z*/
	if (values->VERBOSE) fprintf(stderr,"\n");
	for (i=0;i<values->keep_aligns;i++)	if (from[i]) free(from[i]);
	free(from);
	}/*swap*/

for (i=0;i<values->keep_aligns;i++) {
	if (tree_rep[i]) free(tree_rep[i]);
	if (old_rep[i]) free(old_rep[i]);
	}
free(tree_rep);
free(old_rep);
free(temp_rep);
for (i=0;i<(ntax-1);i++) free(best_nodes[i]);
free(best_nodes);

if (best_aligns[0]->score>bound) {
	for (j=0;j<(*l);j++) best_aligns[j]=dump_align(best_aligns[j]);
	 if (old_best_align) best_aligns[0]=make_align(old_best_align);
	 cur_b_tree=bound;
	 *l=0;
	 if (values->VERBOSE) fprintf(stderr,"	Previously found alignment is better.\n");
	 old_best_align=dump_align(old_best_align);
	 return best_aligns;
	 }
else if (best_aligns[0]->score==bound) {
		/*get all unique*/
		temp=1;
		for (i=0;i<(*l);i++) temp*=compare_aligns(old_best_align,best_aligns[i]);
		if (temp!=0) {
    	if ((*l)==values->keep_aligns) {
				if (values->max_out_aligns) {
        	if ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment *)))!=NULL) {
        		++values->keep_aligns;
        			best_aligns[values->keep_aligns-1]=NULL;
						}
          }
      	}
			if ((*l) < values->keep_aligns) {
				best_aligns[(*l)]=make_align(old_best_align);
				++(*l);
				}
			else if (values->VERBOSE) fprintf(stderr,"Overflow in random alignment generation -- but that's OK.\n");
			}
		if (values->VERBOSE) fprintf(stderr,"	Found %d new alignments.\n",(*l)-1);
		}
else if (values->VERBOSE) fprintf(stderr,"	Found %d better alignments.\n",*l);
if (old_best_align) old_best_align=dump_align(old_best_align);
free(align_passer);
return best_aligns;
}

alignment *make_align(a)
alignment *a;
{
int i;
alignment *b;

if (!a) {fprintf(stderr,"Alignment is not.\n");exit(-1);}
if (!a->n_seqs) {fprintf(stderr,"Alignment with no sequnces.\n");exit(-1);}
if (!a->length) {fprintf(stderr,"Alignment with no bases.\n");exit(-1);}
b=(alignment *)malloc(sizeof(alignment));
assert((int)b);
b->n_seqs=a->n_seqs;
b->length=a->length;
b->score=a->score;
if (a->type_weight!=1)  a->type_weight=0;
b->type_weight=a->type_weight;

b->name=(char *)malloc((strlen(a->name)+1)*sizeof(char));
	assert((int)b->name);
b->name=(char *)strcpy(b->name,a->name);

b->s=(char **)malloc((b->n_seqs+b->type_weight)*sizeof(char *));
assert((int)b->s);
for (i=0;i<(b->n_seqs+b->type_weight);i++) {
	b->s[i]=(char *)malloc((a->length+1)*sizeof(char));
	assert((int)b->s[i]);
	b->s[i]=(char *)strcpy(b->s[i],a->s[i]);
	}

b->taxon_name=(char **)malloc((b->n_seqs)*sizeof(char *));
assert((int)b->taxon_name);
for (i=0;i<b->n_seqs;i++) {
	b->taxon_name[i]=(char *)malloc((strlen(a->taxon_name[i])+1)*sizeof(char));
	assert((int)b->taxon_name[i]);
	b->taxon_name[i]=(char *)strcpy(b->taxon_name[i],a->taxon_name[i]);
	}

return b;
}

