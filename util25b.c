/*Copyright 1992 Ward Wheeler all rights reserved*/

#include "align3.h"

void print_parameters(v)
parameters *v;
{
int i,j;

printf("Parameter values used in MALIGN %s:\n",malign_version);
printf("Internal gap cost %d Leading gap cost %d Trailing gap cost %d ",v->gap_cost,v->gap_cost-v->leading_gap_cost,v->gap_cost-v->trailing_gap_cost);
if (v->extra_adjustment>0) printf("Extragaps %d ",v->gap_cost-v->extra_adjustment);
else if (v->coding >0) printf("Coding %d ",v->gap_cost-v->coding);
printf("\n");
if (!v->ttr && !v->delta) printf("Change cost %d ",v->change_cost);
else if (v->ttr) printf("Transition cost %d Transversion cost %d ",v->transition,v->transition+v->transversion);
else {
	printf("Cost matrix specified:\n");
	for (i=0;i<4;i++) {
		printf("	");
		for (j=0;j<4;j++) printf("%3d ",v->delta[i][j]);
		printf("\n");
		}
	}
printf("\n");
if (v->groups) printf("Groups were specified ");
if (v->groups_as_start) printf("but abandoned after initial build ");
printf("Phylogenetic scoring %d ,",v->phylo_score);
if (v->randtrees>0) printf("-- %d random trees generated ",v->randtrees);
if (v->clade_swap_while_add) printf("with cladogram addition swapping ");
if (v->tree_swap) printf("with cladogram swapping ");
if (v->tree_multi_swap) printf("with tree multi-swapping at level %d ",v->tree_multi_swap);
printf("\n");
if (v->tree_rand_order_max>0) printf("Taxa added in random order %d times\n",v->tree_rand_order_max);
if (v->phylo_gap) printf("Gaps treated as missing data in phylogenetic reconstruction\n");
printf("Maximum number of trees held %d Maximum number of cladograms held %d ",v->keep_trees,v->keep_aligns);
if (v->max_out_trees) printf("with tree buffer increased as needed ");
if (v->max_out_aligns) printf("with alignment buffer increased as needed ");
printf("\n");
if (v->get_heur) printf("Pair alignment peformed ");
if (v->get_heur2) if (v->aquick) printf("Quick alignment peformed ");
if (v->get_heur2) if (!v->aquick) printf("Build alignment peformed ");
if (v->best) printf("Exact (Branch and Bound) alignment peformed ");
if (v->how_many>0) printf("%d Random alignments peformed ",v->how_many);
printf("\n");
if (v->swap_while_add) printf("Alignments swapped while added ");
if (v->align_swap) printf("Simple sequence swapping after addition of all taxa ");
if (v->align_node_swap) printf("Sequence node swapping after addition of all taxa ");
if (v->align_root_swap) printf("with rerooting of pruned taxa ");
if (v->align_partial_root_swap) printf("with partial rerooting of swapped arrangements ");
if (v->align_complete_root_swap) printf("with complete rerooting of swapped arrangements ");
if (v->align_multi_swap) printf("with alignment multi-swapping at level %d ",v->align_multi_swap);
if (v->rand_order>0) printf("Sequences added in random order %d times\n",v->rand_order);
printf("\n");
printf("Output formats ");
if (v->acgt) printf("ascii ");
if (v->farris) printf("Hennig86/CLADOS ");
if (v->inter_dig) printf("interleaved ");
if (v->dot) printf("interleaved with dots ");
if (v->paup) printf("PAUP ");
if (v->paup_dot) printf("PAUP with dots ");
if (v->hen_gap) printf("Hennig86/CLADOS with gap characters ");
printf("with linelength %d ", v->line_length);
if (v->output_order) printf("and output order = input order ");
printf("\n");
if (v->print_intermediates) printf("intermediate solutions printed ");
if (v->VERBOSE) printf("program in verbose mode");
printf("\n");
if (v->new_codes) {
	printf("Amino Acid codes (re)assigned:\n");
	for (i=0;i<v->n_codes;i++) {
		printf("	");
		printf("%c->%c%c%c\n",v->new_codes[i][0],v->new_codes[i][1],v->new_codes[i][2],v->new_codes[i][3]);
		}
	}
if (v->iter) printf("Iter on ");
if (v->low_mem) printf("Low memory use on ");
if (v->show_mem) printf("Showing memory requirement before alignment ");
printf("\n");
if (v->cache_size>0) printf("Alignment cache set at %dk\n",v->cache_size);
if (v->pref_direc==0) printf("Matches favored\n");
else if (v->pref_direc==1) printf("Gaps in the shorter sequences favored\n");
else if (v->pref_direc==2) printf("Gaps in the longer sequences favored\n");
else if (v->pref_direc==3) printf("Random direction chosen\n");
if (v->reorder) printf("Sequences added in order of decreasing similarity\n");
if (v->rand_order>0) printf("Sequences added in random order %d times\n",v->rand_order);
if (v->cull) printf("Culled multiple alignments ");
if (v->elision) printf("Elided multiple alignments ");
if (PARALLEL) {
	printf("Process factor %d and Grain size %d\n",v->process_factor,v->grain_size);
	}
}

void print_alignment (a)
alignment *a;
{
	int i;
	char *p;

	for (i = 0; i < a->n_seqs; i++) {
		printf("%s\n",a->taxon_name[i]);
		for (p = a->s[i]; *p; p++) putchar (*p);
		putchar ('\n');
	}
}


void print_hennig (a,values)
alignment *a;
parameters *values;
{
	int i, j;
	char *p;

if (!values->new_optimization) {
	printf("xread\n");
	if ((!values->delta)&&(!values->ttr)) printf("'%s order %d score %d gap cost %d change cost \n",a->name,a->score,values->gap_cost,values->change_cost);
	else if (values->delta) {
		printf("'%s order %d score %d gap cost \n",a->name,a->score,values->gap_cost);
	printf("cost matrix:\n	A	C	G	T \n");
	printf("A	%3d %3d %3d %3d\n",values->delta[0][0],values->delta[0][1],values->delta[0][2],values->delta[0][3]);
	printf("C	%3d %3d %3d %3d\n",values->delta[1][0],values->delta[1][1],values->delta[1][2],values->delta[1][3]);
	printf("G	%3d %3d %3d %3d\n",values->delta[2][0],values->delta[2][1],values->delta[2][2],values->delta[2][3]);
	printf("T	%3d %3d %3d %3d\n",values->delta[3][0],values->delta[3][1],values->delta[3][2],values->delta[3][3]);
	printf("\n");
	}
	else if (values->ttr) {
		printf("'%s order %d score %d gap cost \n",a->name,a->score,values->gap_cost);
	printf("	transition cost %d transversion cost %d\n",values->transition,values->transversion + values->transition);
	printf("\n");
	}
	if (values->dump_parameters) print_parameters(values);
	printf("'\n");

	if (!values->new_optimization) {
		printf("%d %d\n",a->length,a->n_seqs);
		for (i = 0; i < a->n_seqs; i++) {
			printf("%s\n",a->taxon_name[i]);
			for (j = 0, p = a->s[i]; *p; j++, p++) {
				if (j >= values->line_length) {
					putchar ('\n');
					j = 0;
				}
				if ((*p)=='A') putchar ('0');
				else if ((*p)=='C') putchar ('1');
				else if ((*p)=='G') putchar ('2');
				else if ((*p)=='T') putchar ('3');
				else if ((*p)=='U') putchar ('3');
				else if ((*p)=='-') putchar ('?');
				else putchar ('?');
			}
			putchar ('\n');
		}
		printf(";\n");
		printf("cc-.;\n");
		}
	else {
		printf("1 %d\n",a->n_seqs);
		for (i = 0; i < a->n_seqs; i++) {
			printf("%s 0\n",a->taxon_name[i]);
			}
		printf(";\n");
		printf("cc-.;\n");
		printf("tread\n");
		printf("%s;\n",a->name);
		}
	if (values->ce_weights) {
		for (i=0;i<a->length;i++) {
			printf("cc /%d %d;",values->ce_weights[i],i);
			if (i<(a->length)) printf(" ");
			if (((i+1) % 25)==0) if (i<(a->length-1)) printf("\n");
			}
		}
	printf("\nproc /	;\n");
	}
}

void print_hennig_gap_code (a,values)
alignment *a;
parameters *values;
{
	int i, j,k,n_gaps;
	char *p,m1,m0,*igaps;
	int boundries[2],**sequence,gaps_counter;
	int **pieces, *weights,total_length,multiplier;
	alignment **a1;
	parameters *valuesI;

if (!values->new_optimization) {
	if (values->number_of_input_alignments < 2) {
		igaps=(char *)malloc(a->length*sizeof(char));
		assert((int)igaps);
		n_gaps=0;
		for (j = 0; j<a->length;j++) {
			*(igaps+j)=0;
			for (i = 0; i < a->n_seqs; i++) {
				p=a->s[i];
				if ((*(p+j))=='-') *(igaps+j)=1;
			}
			n_gaps+=*(igaps+j);
		}

		sequence = (int **) malloc(((2*a->n_seqs)-1) * sizeof (int *));
		assert((int)sequence);
		for (i = 0; i < ((2*a->n_seqs)-1); i++) {
			sequence[i] = (int *)malloc ((a->length+n_gaps)*sizeof(int));
			assert((int)sequence[i]);
		}

		for (i = 0; i < a->n_seqs; i++) {
			boundries[0]=0;
			boundries[1]=a->length;
			m0=m1=0;
			for (j = 0, p = a->s[i]; j<a->length; j++) {
				if ((*(p+j))!='-') m1=1;
				if (m1!=m0) {
					boundries[0]=j;
					break;
				}
				else m0=m1;
			}
			m0=m1=0;
			for (j =((a->length)-1), p = a->s[i]; j>=0; j--) {
				if ((*(p+j))!='-') m1=1;
				if (m1!=m0) {
					boundries[1]=j;
					break;
				}
				else m0=m1;
			}
			/*fprintf(stderr,"b[0]=%d b[1]=%d\n",boundries[0],boundries[1]);*/
			gaps_counter=-1;
			for (j = 0, p = a->s[i]; j<a->length; j++) {
				if (*(igaps+j)==1) ++gaps_counter;
				if ((*(p+j))=='A') { *(sequence[i]+j)=2; *(sequence[i]+(gaps_counter+a->length))=2; }
				else if ((*(p+j))=='C') { *(sequence[i]+j)=4; *(sequence[i]+(gaps_counter+a->length))=2; }
				else if ((*(p+j))=='G') { *(sequence[i]+j)=8; *(sequence[i]+(gaps_counter+a->length))=2; }
				else if ((*(p+j))=='T') { *(sequence[i]+j)=16; *(sequence[i]+(gaps_counter+a->length))=2; }
				else if ((*(p+j))=='U') { *(sequence[i]+j)=16; *(sequence[i]+(gaps_counter+a->length))=2; }
				else if ((*(p+j))=='?') { *(sequence[i]+j)=30; *(sequence[i]+(gaps_counter+a->length))=30; }
				else if ((*(p+j))=='N') { *(sequence[i]+j)=30; *(sequence[i]+(gaps_counter+a->length))=2; }
				else if ((*(p+j))=='-') {
					*(sequence[i]+j)=30;
					if ((j>=boundries[0]) && (j<=boundries[1])) *(sequence[i]+(gaps_counter+a->length))=4;
					else *(sequence[i]+(gaps_counter+a->length))=30;
					}
				else  {
					*(sequence[i]+j)=30;
					*(sequence[i]+(gaps_counter+a->length))=30;
					}
		    }
		}

		printf("xread\n");
		if ((!values->delta)&&(!values->ttr)) printf("'%s order %d score %d gap cost %d change cost \n",a->name,a->score,values->gap_cost,values->change_cost);
		else if (values->delta) {
			printf("'%s order %d score %d gap cost \n",a->name,a->score,values->gap_cost);
		printf("cost matrix:\n	A	C	G	T \n");
		printf("A	%3d %3d %3d %3d\n",values->delta[0][0],values->delta[0][1],values->delta[0][2],values->delta[0][3]);
		printf("C	%3d %3d %3d %3d\n",values->delta[1][0],values->delta[1][1],values->delta[1][2],values->delta[1][3]);
		printf("G	%3d %3d %3d %3d\n",values->delta[2][0],values->delta[2][1],values->delta[2][2],values->delta[2][3]);
		printf("T	%3d %3d %3d %3d\n",values->delta[3][0],values->delta[3][1],values->delta[3][2],values->delta[3][3]);
		printf("\n");
		}
		else if (values->ttr) {
			printf("'%s order %d score %d gap cost \n",a->name,a->score,values->gap_cost);
		printf("	transition cost %d transversion cost %d\n",values->transition,values->transversion + values->transition);
		printf("\n");
		}
		if (values->dump_parameters) print_parameters(values);
		printf("'\n");

		printf("%d %d\n",a->length+n_gaps,a->n_seqs);
		for (i = 0; i < a->n_seqs; i++) {
			printf("%s\n",a->taxon_name[i]);
			for (j = 0;j<(a->length+n_gaps); j++) {
				if (*(sequence[i]+j)==2) printf("0");
				else if (*(sequence[i]+j)==4) printf("1");
				else if (*(sequence[i]+j)==8) printf("2");
				else if (*(sequence[i]+j)==16) printf("3");
				else if (*(sequence[i]+j)==30) printf("?");
				if (((j+1) % values->line_length)==0) putchar ('\n');
			}
			putchar ('\n');
		}
		printf(" ;\n");
		if ((!values->delta)&&(!values->ttr)) printf("cc-.;cc/%d 0.%d; cc/%d %d.;\n",values->change_cost,(a->length-1),values->gap_cost,a->length);
		else printf("cc-.;\ncc/%d %d.;\n",values->gap_cost,a->length);
		printf("proc /	;\n");
		free(igaps);
		for (i = 0; i < ((2*a->n_seqs)-1); i++) {
			free(sequence[i]);
		}
		free(sequence);
		}
	else {/*print total evidence shit*/
		a1=(alignment **)malloc(values->number_of_input_alignments*sizeof(alignment *));
		assert((int) a1);
		pieces=(int **)malloc(a->n_seqs*sizeof (int *));
		assert((int)pieces);
		for (i=0;i<a->n_seqs;i++) {
			pieces[i]=(int *)malloc(sizeof(int));
			assert((int) pieces[i]);
			}
		weights=(int *)malloc(a->length*sizeof(int));
		assert((int) weights);
		total_length=0;
		for (i=0;i<values->number_of_input_alignments;i++) {
			if (!values->other_parm) valuesI=values;
			else if (!values->other_parm[i]) valuesI=values;
			else valuesI=values->other_parm[i];
			multiplier=2;
			if (valuesI->ttr) ++multiplier;
			else if (values->data_sets[i]==2) --multiplier; /*morph*/
			a1[i]=NULL;
			a1[i]=make_align(a);
			redo_for_type(a1[i],i,values);
			/*assign raw states*/
			if (values->data_sets[i]==2) { /*morph*/
				for (j=0;j<a->n_seqs;j++) {
					pieces[j]=(int *)realloc(pieces[j],(total_length+a1[i]->length)*sizeof(int));
					assert((int) pieces[j]);
					}
				weights=(int *)realloc(weights,(total_length+a1[i]->length)*sizeof(int));
				assert((int) weights);
				for (j=total_length;j<total_length+a1[i]->length;j++) weights[j]=valuesI->change_cost*values->data_set_weights[i];
				for (j=0;j<a1[i]->n_seqs;j++) for (k=0;k<a1[i]->length;k++) {
					if (a1[i]->s[j][k]=='0') pieces[j][k+total_length]='0';
					else if (a1[i]->s[j][k]=='1') pieces[j][k+total_length]='1';
					else if (a1[i]->s[j][k]=='2') pieces[j][k+total_length]='2';
					else if (a1[i]->s[j][k]=='3') pieces[j][k+total_length]='3';
					else if (a1[i]->s[j][k]=='4') pieces[j][k+total_length]='4';
					else if (a1[i]->s[j][k]=='5') pieces[j][k+total_length]='5';
					else if (a1[i]->s[j][k]=='6') pieces[j][k+total_length]='6';
					else if (a1[i]->s[j][k]=='?') pieces[j][k+total_length]='?';
					else if (a1[i]->s[j][k]=='-') pieces[j][k+total_length]='?';
					else pieces[j][k+total_length]='?';
					}
				total_length+=a1[i]->length;
				}
			else {
				if (!valuesI->ttr) {      /*regular*/
					for (j=0;j<a->n_seqs;j++) {
						pieces[j]=(int *)realloc(pieces[j],(total_length+(2*a1[i]->length))*sizeof(int));
						assert((int) pieces[j]);
						}
					weights=(int *)realloc(weights,(total_length+(2*a1[i]->length))*sizeof(int));
					assert((int) weights);
					for (j=total_length;j<total_length+a1[i]->length;j++) weights[j]=valuesI->change_cost*values->data_set_weights[i];
					for (j=total_length+a1[i]->length;j<total_length+(2*a1[i]->length);j++) weights[j]=valuesI->gap_cost*values->data_set_weights[i];
					for (j=0;j<a1[i]->n_seqs;j++) for (k=0;k<a1[i]->length;k++) {
						switch (a1[i]->s[j][k]) {
							case 'A': pieces[j][k+total_length]='0';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'C': pieces[j][k+total_length]='1';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'G': pieces[j][k+total_length]='2';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'T': pieces[j][k+total_length]='3';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case '-': {
								if (values->phylo_gap) {
									if (!(is_internal(a1[i]->s[j],k,a1[i]->length))) {pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; }
									else {pieces[j][k+total_length]='4';	pieces[j][k+total_length+a1[i]->length]='?'; }
									}
								else {pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; }
								break;
								}
							case 'N': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'X': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; break;
							case 'U': pieces[j][k+total_length]='3';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'R': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'Y': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'M': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'W': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'S': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'K': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'B': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'D': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'H': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							case 'V': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; break;
							default:
								fprintf (stderr, "Bad character %c in sequence[%d][%d] of input %d\n",a1[i]->s[j][k], j, k,i);
								print_alignment (a1[i]);
								exit(-1);
							}
						}
					total_length+=(2*a1[i]->length);
					}
				else {				 /*transition transversion*/
					for (j=0;j<a->n_seqs;j++) {
						pieces[j]=(int *)realloc(pieces[j],(total_length+(3*a1[i]->length))*sizeof(int));
						assert((int) pieces[j]);
						}
					weights=(int *)realloc(weights,(total_length+(3*a1[i]->length))*sizeof(int));
					assert((int) weights);
					for (j=total_length;j<total_length+a1[i]->length;j++) weights[j]=valuesI->transition*values->data_set_weights[i];
					for (j=total_length+a1[i]->length;j<total_length+(2*a1[i]->length);j++) weights[j]=valuesI->transversion*values->data_set_weights[i];
					for (j=total_length+(2*a1[i]->length);j<total_length+(3*a1[i]->length);j++) weights[j]=valuesI->gap_cost*values->data_set_weights[i];
					for (j=0;j<a1[i]->n_seqs;j++) for (k=0;k<a1[i]->length;k++) {
						switch (a1[i]->s[j][k]) {
							case 'A': pieces[j][k+total_length]='0';	pieces[j][k+total_length+a1[i]->length]='0'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'C': pieces[j][k+total_length]='1';	pieces[j][k+total_length+a1[i]->length]='1'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'G': pieces[j][k+total_length]='2';	pieces[j][k+total_length+a1[i]->length]='0'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'T': pieces[j][k+total_length]='3';	pieces[j][k+total_length+a1[i]->length]='1'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case '-': {
								if (values->phylo_gap) {
									if (is_internal(a1[i]->s[j],k,a1[i]->length)) {pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='1';}
									else {pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='?';}
									}
								else {pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='?';}
								break;
							}
							case 'N': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'X': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='?'; break;
							case 'U': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='1'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'R': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='0'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'Y': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='1'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'M': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'W': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'S': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'K': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'B': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'D': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'H': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							case 'V': pieces[j][k+total_length]='?';	pieces[j][k+total_length+a1[i]->length]='?'; pieces[j][k+total_length+(2*a1[i]->length)]='0'; break;
							default:
								fprintf (stderr, "Bad character %c in sequence[%d][%d] of input %d\n",a1[i]->s[j][k], j, k,i);
								print_alignment (a1[i]);
								exit(-1);
							}
						}
					total_length+=(3*a1[i]->length);
					}
				}
			}
		/*print file*/
		printf("xread\n");
		printf("'Data from %d input files\nAlignment order %s'\n",values->number_of_input_alignments,a->name);
		printf("%d %d\n",total_length,a->n_seqs);
		for (i=0;i<a->n_seqs;i++) {
			printf("%-20s\n",a->taxon_name[i]);
			for (j=0;j<total_length;j++) {
				printf("%c",pieces[i][j]);
				if (!((j+1)%values->line_length)) printf("\n");
				}
			printf("\n");
			}
		printf("\n;\ncc-.;\n");
		printf("cc/%d 0",weights[0]);
		for (j=1;j<total_length-1;j++) {
			if (!(weights[j]==weights[j-1])) {
				printf(".%d; ",j-1);
				printf("cc/%d %d",weights[j],j);
				}
			/*if (!((j+1)%(values->line_length/5))) printf("\n"); */
			}

		if (weights[total_length-1]!=weights[total_length-2]) printf(".%d; cc/%d %d",total_length-2,weights[total_length-1],total_length-1);
		else printf(".%d;\n",total_length-1);
		printf("proc/;\n");
		/*deallocate local stuff*/
		free(weights);
		for (i=0;i<a->n_seqs;i++) free(pieces[i]);
		free(pieces);
		for (i=0;i<values->number_of_input_alignments;i++) if (a1[i]) a1[i]=dump_align(a1[i]);
		free(a1);
		}
	}
}

void print_inter_dig (a,values)
alignment *a;
parameters *values;
{
	int i, j, k, l, m, n;
	char *p;
if (!values->new_optimization) {
	printf("%d ",a->n_seqs);
	printf("%d\n",(a->length/values->line_length)+1);
	for (l=0; l<=(a->length/values->line_length);l++){
		for (i = 0; i < a->n_seqs; i++) {
		printf("%-15s",a->taxon_name[i]);
			for (j = l*values->line_length, p = (a->s[i]+j); ((j<((l+1)*values->line_length)) && (*p)); j++, p++) {
				putchar (*p);
			}
			putchar('\n');
		}
		putchar('\n');
	}
}
printf("%s\n",a->name);
if (values->dump_parameters) print_parameters(values);
}

void print_paup (a,values)
alignment *a;
parameters *values;
{
	int i, j, l, n;
	char *p;

if (!values->new_optimization) {
	printf("#NEXUS\n\n");
	if ((!values->delta)&&(!values->ttr))	printf("[!order %s change cost = %d gap cost = %d (U's converted to T's)",a->name,values->change_cost,values->gap_cost);
	else {
		printf("[!order %s gap cost = %d (U's converted to T's)\n",a->name,values->change_cost,values->gap_cost);
	if (values->delta) {
		printf("cost matrix:\n	A	C	G	T \n");
		printf("A	%3d %3d %3d %3d\n",values->delta[0][0],values->delta[0][1],values->delta[0][2],values->delta[0][3]);
		printf("C	%3d %3d %3d %3d\n",values->delta[1][0],values->delta[1][1],values->delta[1][2],values->delta[1][3]);
		printf("G	%3d %3d %3d %3d\n",values->delta[2][0],values->delta[2][1],values->delta[2][2],values->delta[2][3]);
		printf("T	%3d %3d %3d %3d\n",values->delta[3][0],values->delta[3][1],values->delta[3][2],values->delta[3][3]);
		printf("\n");
		}
		else if (values->ttr) {
		printf("	transition cost %d transversion cost %d\n",values->transition,values->transversion + values->transition);
		 	}
		}
	if (values->dump_parameters) print_parameters(values);
	printf("]\n\n");

	printf("begin data;\n");
	printf("	dimensions ntax=%d nchar=%d;\n",a->n_seqs,a->length);
	printf("	format datatype=dna matchchar=. interleave missing=-;\n");
	printf("	matrix\n");

	for (l=0; l<=(a->length/values->line_length);l++){
		for (i = 0; i < a->n_seqs; i++) {
		printf("%-15s",a->taxon_name[i]);
			for (j = l*values->line_length, p = (a->s[i]+j); ((j<((l+1)*values->line_length)) && (*p)); j++, p++) {
				if (*p=='U') putchar('T');
				else putchar (*p);
			}
			putchar('\n');
		}
		putchar('\n');
	}
	printf(";\n");
	printf("end;\n");
	if (values->ce_weights) {
		printf("begin assumptions;\n");
		printf("	wtset Ellision = ");
		for (i=0;i<a->length;i++) {
			printf("%d: %d",values->ce_weights[i],i+1);
			if (i<(a->length)) printf(", ");
			if (((i+1) % 25)==0) printf("\n		");
			}
	printf(";\n");
	printf("end;\n\n");
	}
}
}

void print_paup_dot (a,values)
alignment *a;
parameters *values;
{
	int i, j,	l, n;
	char *p,*q;

if (!values->new_optimization){
	printf("#NEXUS\n");
	if ((!values->delta)&&(!values->ttr))	printf("[!order %s change cost = %d gap cost = %d (U's converted to T's)",a->name,values->change_cost,values->gap_cost);
	else {
		printf("[!order %s gap cost = %d (U's converted to T's)\n",a->name,values->change_cost,values->gap_cost);
	if (values->delta) {
		printf("cost matrix:\n	A	C	G	T \n");
		printf("A	%3d %3d %3d %3d\n",values->delta[0][0],values->delta[0][1],values->delta[0][2],values->delta[0][3]);
		printf("C	%3d %3d %3d %3d\n",values->delta[1][0],values->delta[1][1],values->delta[1][2],values->delta[1][3]);
		printf("G	%3d %3d %3d %3d\n",values->delta[2][0],values->delta[2][1],values->delta[2][2],values->delta[2][3]);
		printf("T	%3d %3d %3d %3d\n",values->delta[3][0],values->delta[3][1],values->delta[3][2],values->delta[3][3]);
		printf("\n");
		}
		else if (values->ttr) {
		printf("	transition cost %d transversion cost %d\n",values->transition,values->transversion + values->transition);
		 	}
	}
		if (values->dump_parameters) print_parameters(values);
		printf("]\n\n");
	printf("begin data;\n");
	printf("	dimensions ntax=%d nchar=%d;\n",a->n_seqs,a->length);
	printf("	format datatype=dna matchchar=. interleave missing=-;\n");
	printf("	matrix\n");


	for (l=0; l<=(a->length/values->line_length);l++){
		for (i = 0; i < a->n_seqs; i++) {
		printf("%-15s",a->taxon_name[i]);
			for (j = l*values->line_length, q = (a->s[0]+j), p = (a->s[i]+j); ((j<((l+1)*values->line_length)) && (*p)); j++, p++,q++) {
				if (*p=='U') *p='T';
				if (*q=='U') *q='T';
				if ((i==0) || ((*p)!=(*q))) putchar (*p);
				else putchar('.');
			}
			putchar('\n');
		}
		putchar('\n');
	}

	printf(";\n");
	printf("end;\n");

	if (values->ce_weights) {
		printf("begin assumptions;\n");
		printf("	wtset Ellision = ");
		for (i=0;i<a->length;i++) {
			printf("%d: %d",values->ce_weights[i],i+1);
			if (i<(a->length)) printf(", ");
			if (((i+1) % 25)==0) printf("\n		");
			}
	printf(";\n");
	printf("end;\n\n");
	}
}
}

void print_dot (a,values)
alignment *a;
parameters *values;
{
	int i, j, l, n;
	char *p,*q;

if (!values->new_optimization){
	for (l=0; l<=(a->length/values->line_length);l++){
		for (i = 0; i < a->n_seqs; i++) {
		printf("%-15s",a->taxon_name[i]);
			for (j = l*values->line_length, q = (a->s[0]+j), p = (a->s[i]+j); ((j<((l+1)*values->line_length)) && (*p)); j++, p++,q++) {
				if ((i==0) || ((*p)!=(*q))) putchar (*p);
				else putchar('.');
			}
			putchar('\n');
		}
		putchar('\n');
	}
}
printf("%s\n",a->name);
	if (values->dump_parameters) print_parameters(values);

}


void printem(best_aligns,count,values,a)
alignment **best_aligns, **a;
int count;
parameters *values;
{
int i,j,counter,pos;
char *buf;

if (values->farris || values->hen_gap) if (values->new_optimization) {

	printf("xread\n");
	printf("'Under new optimization there were %d equally costly rooted topologies \n	this is a dummy file to examine and print out trees",count);
	if ((!values->delta)&&(!values->ttr)) printf(" %d score %d gap cost %d change cost \n",best_aligns[0]->score,values->gap_cost,values->change_cost);
	else if (values->delta) {
		printf(" %d score %d gap cost \n",best_aligns[0]->score,values->gap_cost);
		printf("cost matrix:\n	A	C	G	T \n");
		printf("A	%3d %3d %3d %3d\n",values->delta[0][0],values->delta[0][1],values->delta[0][2],values->delta[0][3]);
		printf("C	%3d %3d %3d %3d\n",values->delta[1][0],values->delta[1][1],values->delta[1][2],values->delta[1][3]);
		printf("G	%3d %3d %3d %3d\n",values->delta[2][0],values->delta[2][1],values->delta[2][2],values->delta[2][3]);
		printf("T	%3d %3d %3d %3d\n",values->delta[3][0],values->delta[3][1],values->delta[3][2],values->delta[3][3]);
		printf("\n");
		}
	else if (values->ttr) {
		printf(" %d score %d gap cost \n",best_aligns[0]->score,values->gap_cost);
		printf("	transition cost %d transversion cost %d\n",values->transition,values->transversion + values->transition);
		printf("\n");
		}
	if (values->dump_parameters) print_parameters(values);
	printf("'\n");

	/*fprintf(stderr,"Max name length=%d\n",values->max_name_length);*/

	printf("1 %d\n",values->all_done);
	if (!values->nona) {
		for (i = 0; i < 2; i++) printf("%s 0\n",values->input_names[i]);
		for (i = 2; i < values->all_done; i++) printf("%s 1\n",values->input_names[i]);
		}
	else {
		for (i = 0; i < 2; i++) {
			for (j=0;j< min(strlen(values->input_names[i]),values->max_name_length);j++) printf("%c",values->input_names[i][j]);
			printf(" 0\n");
			}
		for (i = 2; i < values->all_done; i++) {
			for (j=0;j<min(strlen(values->input_names[i]),values->max_name_length);j++) printf("%c",values->input_names[i][j]);
			printf(" 1\n");
			}
		}
	printf(";\n");
	printf("cc-.;\n");
	printf("tread\n");
	for (i=0;i<count;i++) {
		buf=(char *)malloc((1+strlen(best_aligns[i]->name))*sizeof(char));
		assert((int) buf);
		counter=0;
		pos=0;
		for (j=0;j<strlen(best_aligns[i]->name);j++) {
			if (best_aligns[i]->name[j]=='(') buf[pos++]=best_aligns[i]->name[j];
			else if ((best_aligns[i]->name[j]!=')') && (best_aligns[i]->name[j]!=' ')) {
				if (counter < values->max_name_length) buf[pos++]=best_aligns[i]->name[j];
				counter++;
				}
			else if ((best_aligns[i]->name[j]==')') || (best_aligns[i]->name[j]==' ')) {
				buf[pos++]=best_aligns[i]->name[j];
				counter=0;
				}
			}
		buf[pos]='\0';
		printf("%s",buf);
		if (i<(count-1)) printf("*\n");
		else printf(";\n");
		free(buf);
		}
	printf("\nproc /	;\n");
}


for (i=0;i<count;i++) {
	if ((!values->cull) && (!values->elision) && (!values->new_optimization)) if (values->output_order==1) order_as_input(best_aligns[i],values);
	if (values->acgt) {
			print_alignment (best_aligns[i]);
			printf("\n");
		 }
		if (values->farris) {
			if (!values->new_optimization) print_hennig (best_aligns[i], values);
			printf("\n");
		}
		if (values->inter_dig) {
			print_inter_dig (best_aligns[i], values);
			printf("\n");
		}
	if (values->dot) {
			print_dot (best_aligns[i], values);
			printf("\n");
		}
		if (values->hen_gap) {
			print_hennig_gap_code (best_aligns[i], values);
			printf("\n");
			}
	if (values->paup) {
			print_paup (best_aligns[i], values);
			printf("\n");
			}
	if (values->paup_dot) {
			print_paup_dot (best_aligns[i], values);
			printf("\n");
		}
	}
}

