/*Copyright 1992 Ward Wheeler all rights reserved*/

#include "align3.h"

void get_guess_max_gap(a,ntax,values)
alignment **a;
int ntax;
parameters *values;
{
	int i;

	/*get initial gaps max*/
	if (values->max_gap == MAX_SEQUENCE_SIZE) {
		values->max_gap=0;
		for (i=0;i<ntax;i++) if (a[i]->length > values->max_gap) values->max_gap=a[i]->length;
		values->max_gap=values->max_gap/10;
	}
}

/*
int max(a,b)
int a,b;
{
if (a>b) return a;
return b;
}
*/

realloc_align_holder *do_addition_swap(counter,old_rep,a,best_aligns,nodes,n_nodes,l,values,current,ntaxa,align_passer)
int *counter, **old_rep, **nodes, **n_nodes,*l;
int ntaxa,current;
parameters *values;
alignment **a,**best_aligns;
realloc_align_holder *align_passer;
{
	int i,j,k,m,sub_taxa,found_how_many;
	int cur_val,best_value,**tree_rep,*temp_rep;
	char overflow,temp,found_better_or_more,**from;
	int *intermediate_scores,intermediate_best_score;
	alignment **intermediate_align;
	int num_best_aligns, next_to_do, old_best,found_new;
	int bufid, bytes, type, source, max_to_do, index, grain, num_left, info, n_sub_taxa;


	best_value=intermediate_best_score=best_aligns[0]->score;
	/*allocating*/
	sub_taxa=current+4;

	tree_rep=(int **)malloc(values->keep_aligns *sizeof(int *));
	assert((int)tree_rep);
	for (i=0;i<values->keep_aligns;i++) {
		tree_rep[i]=(int *)malloc((sub_taxa-3)*sizeof(int));
		assert((int)tree_rep[i]);
	}
	from=(char **)malloc(values->keep_aligns *sizeof(char *));
	assert((int)from);
	for (i=0;i<values->keep_aligns;i++) {
		from[i]=(char *)malloc((sub_taxa-3)*sizeof(char));
		assert((int)from[i]);
		for (j=0;j<(sub_taxa-3);j++) from[i][j]=0;
	}
	temp_rep=(int *)malloc((sub_taxa-3)*sizeof(int));
	assert((int)temp_rep);
	/*these are paralell formulations to be modified later but better work in sequential*/
	intermediate_scores=(int *)malloc(((2*(ntaxa-3))+3+1)*sizeof(int));
	assert((int) intermediate_scores);
	intermediate_align=(alignment **)malloc(((2*(ntaxa-3))+3+1)*sizeof(alignment *));
	assert((int) intermediate_align);
	for (i=0;i<((2*(ntaxa-3))+3+1);i++) intermediate_align[i]=NULL;

	found_how_many=old_best=0;
top:
	;
	if (values->VERBOSE) {
		fprintf(stderr,"(swapping");
		fflush(stderr);
	}
	/*Initialization*/
	overflow=found_better_or_more=found_new=0;
	if ((*counter)>(values->keep_aligns-1)) (*counter)=values->keep_aligns-1;
	for (i=0;i<=(*counter);i++) for (j=0;j<(sub_taxa-3);j++) tree_rep[i][j]=old_rep[i][j];/*add if!old rep*/
	for (m=old_best;m<=(*counter);m++) {
		for (i=0;i<(sub_taxa-3);i++) temp_rep[i]=old_rep[m][i];
		for (i=0;i<(sub_taxa-4);i++) if (!from[m][i]) {
			/*redoing for paralell structure here*/
			for (j=1;j<=((2*i)+3);j++) {
				intermediate_scores[j]=HUGE_COST;
				if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
				intermediate_align[j]=NULL;
			}
			if (PARALLEL) {
				grain=values->grain_size;
				/*fprintf(stderr,"Sending jobs to");*/
				num_left=max_to_do=(2*i)+3;
				for (j=1;((j<values->num_hosts) && (j<max_to_do));j++) {
					temp_rep[i]=j;
					pvm_initsend( PvmDataDefault );
					pvm_pkint(&grain,1,1);
					pvm_pkint(&ntaxa,1,1);
					n_sub_taxa=sub_taxa-4;
					pvm_pkint(&n_sub_taxa,1,1);
					pvm_pkint(temp_rep,n_sub_taxa+1,1);
					pvm_pkint(&j,1,1);
					pvm_send(values->tids[j],3);
					/*fprintf(stderr," sent\n");*/
				}
				while (num_left) {
					/*fprintf(stderr,"num_left=%d ",num_left);*/
					bufid=pvm_recv(-1,4);
					if (bufid) {
						info=pvm_bufinfo(bufid, &bytes, &type, &source);
						if (j<=max_to_do) {
							temp_rep[i]=j;
							pvm_initsend( PvmDataDefault );
							pvm_pkint(&grain,1,1);
							pvm_pkint(&ntaxa,1,1);
							n_sub_taxa=sub_taxa-4;
							pvm_pkint(&n_sub_taxa,1,1);
							pvm_pkint(temp_rep,n_sub_taxa+1,1);
							pvm_pkint(&j,1,1);
							pvm_send(source,3);
							++j;
						}
						--num_left;
						pvm_upkint(&grain,1,1);
						pvm_upkint(&index,1,1);
						pvm_upkint(&(intermediate_scores[index]),1,1);
						if (intermediate_scores[index]<=intermediate_best_score) {
							if (intermediate_scores[index]<intermediate_best_score) {
								/* code folded from here */
								intermediate_best_score=intermediate_scores[index];
								if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
								intermediate_align[index]=unpack_align(intermediate_scores[index],values);
							}
							else if (!values->aquick) {
								/* code folded from here */
								if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
								intermediate_align[index]=unpack_align(intermediate_scores[index],values);
							}
							else pvm_freebuf(bufid);
						}/*j<*/
					}/*if message recieved*/
				}/*num_left*/
			}
			/*sequential*/
else {
	for (j=1;j<=((2*i)+3);j++) {
		/*if (intermediate_align[j]) dump_align(intermediate_align[j]);*/
		temp_rep[i]=j;
		intermediate_scores[j]=all_make_nodes(a,temp_rep,nodes,ntaxa,sub_taxa,values);
		if (intermediate_scores[j]<=intermediate_best_score) {
			if (intermediate_scores[j]<intermediate_best_score) {
				intermediate_best_score=intermediate_scores[j];
				all_other_make_nodes2(temp_rep,nodes,ntaxa,sub_taxa,n_nodes);
				if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
				intermediate_align[j]=make_align(a[n_nodes[0][1]]);
				/*for (k=1;k<j;k++) if (intermediate_align[k]) dump_align(intermediate_align[k]);*/
			}
			else if ((!values->aquick) && (intermediate_scores[j]!=HUGE_COST)) {
				all_other_make_nodes2(temp_rep,nodes,ntaxa,sub_taxa,n_nodes);
				if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
				intermediate_align[j]=make_align(a[n_nodes[0][1]]);
			}
		}
	else if (intermediate_align[j]) intermediate_align[j]=dump_align(intermediate_align[j]);
	}
}
num_best_aligns=0;
for (j=1;j<=((2*i)+3);j++) if (intermediate_scores[j]==intermediate_best_score) ++num_best_aligns;
if ((best_value>intermediate_best_score) || (num_best_aligns>0)) { /*see if swap did anything*/
	if (intermediate_best_score < best_value) {
		best_value=intermediate_best_score;
		for (k=0;k<(*l);k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
		(*l)=0;
		found_new=0;
		if (sub_taxa==ntaxa) if (values->VERBOSE) fprintf(stderr,"Found better alignment at %d\n",intermediate_best_score);				
	}
	for (j=1;j<=((2*i)+3);j++) if ((intermediate_scores[j]==intermediate_best_score) && (intermediate_align[j])) {
		if (!(*l)) {
			found_better_or_more=1;
			(*l)=1;
			if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
			best_aligns[0]=make_align(intermediate_align[j]);
			for (k=0;k<(sub_taxa-3);k++) tree_rep[0][k]=temp_rep[k];
			tree_rep[0][i]=j;
			from[0][i]=1;
		}/*better*/
else if (!values->aquick) {
	temp=1;
	for (k=0;k<(*l);k++) if (best_aligns[k]) temp*=compare_aligns(intermediate_align[j],best_aligns[k]);
	if (temp) {
		if ((*l)>(values->keep_aligns-2)) {
			if ((values->max_out_aligns) && ((best_aligns=(alignment **)realloc(best_aligns,(values->keep_aligns+1)*sizeof(alignment
			    *)))!=NULL) && ((tree_rep=(int **)realloc(tree_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) &&
			    ((old_rep=(int **)realloc(old_rep,(values->keep_aligns+1)*sizeof(int *)))!=NULL) && ((from=(char
			**)realloc(from,(values->keep_aligns+1)*sizeof(char *)))!=NULL)) {
				++values->keep_aligns;
				best_aligns[values->keep_aligns-1]=NULL;
				tree_rep[values->keep_aligns-1]=(int *)malloc((ntaxa-3)*sizeof(int));
				assert((int)tree_rep[values->keep_aligns-1]);
				old_rep[values->keep_aligns-1]=(int *)malloc((ntaxa-3)*sizeof(int));
				assert((int)old_rep[values->keep_aligns-1]);
				from[values->keep_aligns-1]=(char *)malloc((ntaxa-3)*sizeof(char));
				assert((int)from[values->keep_aligns-1]);
				for (k=0;k<ntaxa-3;k++) from[values->keep_aligns-1][k]=0;
			}
		}
		if ((*l)>(values->keep_aligns-2)) {
			if (!overflow) {
				if (values->rep_error) fprintf(stderr,"Overflow in alignment addition swap -- but that's OK.\n");
				overflow=1;
				}
		}
		else {
			for (k=0;k<sub_taxa-3;k++) tree_rep[(*l)][k]=temp_rep[k];
			tree_rep[(*l)][i]=j;
			best_aligns[(*l)]=make_align(intermediate_align[j]);
			from[(*l)][i]=1;
			++(*l);
			found_better_or_more=1;
			if (sub_taxa==ntaxa) if (values->VERBOSE) fprintf(stderr,"Found another alignment for %d\n",*l);

		}
	}
}/*additional*/
	}	/*j*/

}/*swap did something*/
temp_rep[i]=old_rep[m][i];
		}/*i*/
	}/*m counter stuff*/
	for (j=0;j<(*l);j++) for (i=0;i<(sub_taxa-3);i++) old_rep[j][i]=tree_rep[j][i];
	if (found_new) old_best=0;
	else old_best=(*counter);
	(*counter)=(*l)-1;

	found_how_many+=found_better_or_more;
	if (found_better_or_more) goto top;
	free(temp_rep);
	for (i=0;i<values->keep_aligns;i++) {
		free(tree_rep[i]);
		free(from[i]);
	}
	free(tree_rep);
	free(from);
	free(intermediate_scores);
	for (i=0;i<((2*(ntaxa-3))+3+1);i++) if (intermediate_align[i]) intermediate_align[i]=dump_align(intermediate_align[i]);
	free(intermediate_align);

	if (values->VERBOSE) {
		for (i=0;i<found_how_many;i++) fprintf(stderr,")");
		fprintf(stderr,") ");
		fflush(stderr);
	}

	align_passer->best=best_aligns;
	align_passer->rep=old_rep;
	return align_passer;
}

void order_as_input(a,values)
alignment *a;
parameters *values;
{
	int i,j;
	char **bs,**btaxon_name;

	bs=(char **)malloc(a->n_seqs*sizeof(char *));
	assert((int) bs);
	btaxon_name=(char **)malloc(a->n_seqs*sizeof(char *));
	assert((int) btaxon_name);

	for (i=0;i<a->n_seqs;i++) {
		for (j=0;j<a->n_seqs;j++) {
			if (!strcmp(values->input_names[i],a->taxon_name[j])) {
				bs[i]=(char *)malloc((1+(strlen(a->s[j])))*sizeof(char));
				assert((int)bs[i]);
				strcpy(bs[i],a->s[j]);
				btaxon_name[i]=(char *)malloc((1+(strlen(a->taxon_name[j])))*sizeof(char));
				assert((int)btaxon_name[i]);
				strcpy(btaxon_name[i],a->taxon_name[j]);
			}
		}
	}

	for (i=0;i<a->n_seqs;i++) {
		free(a->s[i]);
		a->s[i]=(char *)malloc((1+a->length)*sizeof(char));
		assert((int)a->s[i]);
		strcpy(a->s[i],bs[i]);
		free(bs[i]);

		free(a->taxon_name[i]);
		a->taxon_name[i]=(char *)malloc((1+(strlen(btaxon_name[i])))*sizeof(char));
		assert((int)a->taxon_name[i]);
		strcpy(a->taxon_name[i],btaxon_name[i]);
		free(btaxon_name[i]);
	}
	free(bs);
	free(btaxon_name);
}

int **get_combinatorial_weights(best_aligns,values,count)
alignment **best_aligns;
int count;
parameters *values;
{
	int **w_matrix,**t_matrix,i,j,k,biggest;
	/*allocating*/
	w_matrix=(int **)malloc(5*sizeof(int *));
	assert((int)w_matrix);
	t_matrix=(int **)malloc(5*sizeof(int *));
	assert((int)t_matrix);
	for (i=0;i<5;i++) {
		w_matrix[i]=(int *)malloc(5*sizeof(int));
		assert((int)w_matrix[i]);
		t_matrix[i]=(int *)malloc(5*sizeof(int));
		assert((int)t_matrix[i]);
	}

	for (j=0;j<5;j++) for (k=0;k<5;k++) w_matrix[j][k]=0;
	for (i=0;i<count;i++) {
		get_local_matrix(t_matrix,best_aligns[i],values);
		for (j=0;j<5;j++) for (k=0;k<5;k++) if (j!=k) w_matrix[j][k]+=t_matrix[j][k];
	}

	/*normalize*/
	biggest=0;
	for (i=0;i<5;i++) for (j=0;j<5;j++) if (w_matrix[i][j]>biggest) biggest=w_matrix[i][j];
	for (i=0;i<5;i++) for (j=0;j<5;j++) w_matrix[i][j]=(w_matrix[i][j]*values->weight_range/biggest);
	/*freeing*/
	for (i=0;i<5;i++) free(t_matrix[i]);
	free(t_matrix);

	/*returning*/
	/*
for (i=0;i<5;i++) {
	for (j=0;j<5;j++) fprintf(stderr,"%3d ",w_matrix[i][j]);
  fprintf(stderr,"\n");
  }
*/
	return w_matrix;
}

void get_local_matrix(m,align,values)
int **m;
alignment *align;
parameters *values;
{
	int i,j,a,c,g,t,gap,h[5],amount;

	for (i=0;i<5;i++) for (j=0;j<5;j++) m[i][j]=0;
	for (i=0;i<align->length;i++) {
		a=c=g=t=gap=0;
		for (j=0;j<align->n_seqs;j++) {
			switch(align->s[j][i]) {
			case 'A':
				a=1;
				break;
			case 'C':
				c=1;
				break;
			case 'G':
				g=1;
				break;
			case 'T':
				t=1;
				break;
			case '-':
				{
					if (is_internal(align->s[j],i,align->length)) gap=1;
					break;
				}
			default:
				break;
			}
		}
		if ((a+c+g+t+gap) > 1) {
			/* 2 x 30 so all even after division */
			amount=(60)/(a+c+g+t+gap);
			if ((a==1) && (c==1))   m[0][1]+=amount;
			if ((a==1) && (g==1))   m[0][2]+=amount;
			if ((a==1) && (t==1))   m[0][3]+=amount;
			if ((a==1) && (gap==1)) m[0][4]+=amount;
			if ((c==1) && (g==1))   m[1][2]+=amount;
			if ((c==1) && (t==1))   m[1][3]+=amount;
			if ((c==1) && (gap==1)) m[1][4]+=amount;
			if ((g==1) && (t==1))   m[2][3]+=amount;
			if ((g==1) && (gap==1)) m[2][4]+=amount;
			if ((t==1) && (gap==1)) m[3][4]+=amount;
		}
	}
	for (i=1;i<4;i++) m[0][4]+=m[i][4];
	m[0][4]/=4;
	for (i=1;i<4;i++) m[i][4]=m[0][4];
	for (i=0;i<5;i++) for (j=0;j<5;j++) if (m[i][j]==0) if (i!=j) m[i][j]=1;/*avoids infinities*/
	for (i=0;i<5;i++) for (j=0;j<i;j++) m[i][j]=m[j][i];
	for (i=0;i<5;i++) {
		h[i]=0;
		for (j=0;j<5;j++) h[i]+=m[j][i];
	}
	for (i=0;i<5;i++) for (j=0;j<5;j++) if (i!=j) m[i][j]=h[j]/m[i][j];
	for (i=0;i<5;i++) for (j=0;j<i;j++) if (i!=j) m[i][j]=m[j][i]=((m[i][j]+m[j][i])/2);

	/*
for (i=0;i<5;i++) {
	for (j=0;j<5;j++) fprintf(stderr,"%d ",m[i][j]);
  fprintf(stderr,"\n");
  }
*/
}

char is_internal(seq,position,length)
char *seq;
int position,length;
{
	int i;
	char before,after;

	if (position==0) return 0;
	if (position==length) return 0;

	before=after=0;
	for (i=position-1;i>=0;i--) {
		if (seq[i]!='-') {
			before=1;
			break;
		}
	}
	for (i=position+1;i<length;i++) {
		if (seq[i]!='-') {
			after=1;
			break;
		}
	}

	if ((before) && (after)) return 1;
	return 0;
}

char is_leading(seq,position,length)
char *seq;
int position,length;
{
	int i;
	char before,after;

	if (position==0) return 0;
	if (position==length) return 0;

	before=after=0;
	for (i=0;i<position;i++) {
		if (seq[i]!='-') {
			before=1;
			break;
		}
	}
	for (i=position+1;i<length;i++) {
		if (seq[i]!='-') {
			after=1;
			break;
		}
	}

	if ((before) && (!after)) return 1;
	return 0;
}

char is_trailing(seq,position,length)
char *seq;
int position,length;
{
	int i;
	char before,after;

	if (position==0) return 0;
	if (position==length) return 0;

	before=after=0;
	for (i=0;i<position;i++) {
		if (seq[i]!='-') {
			before=1;
			break;
		}
	}
	for (i=position+1;i<length;i++) {
		if (seq[i]!='-') {
			after=1;
			break;
		}
	}

	if ((!before) && (after)) return 1;
	return 0;
}
