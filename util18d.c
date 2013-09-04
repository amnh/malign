/*Copyright 1993 Ward Wheeler all rights reserved*/

/*
	To try--
	1)perhaps add size specific caching 
	2) make cache work for cull
	3) cache at ntaxa for cladecache
 */

#include "align3.h"
int local_hack;

/*bump this up when clade cache in to see in makes a difference--remember to free the last as well at end*/
align_buffer **allocate_align_cache(size,values,a)
int size;
alignment **a;
parameters *values;
{
	int i,j,number,temp,cache_size,temp_ho;
	align_buffer **p;
	int cache_classes;

	/*cache n-taxa when clade cacheing is on*/
	/*if (!values->check_cache_for_scores) cache_classes=values->all_done;
else cache_classes=values->all_done+1;*/
	cache_classes=values->all_done;

	size*=1000;
	/*estimate appropriate cache numbers from memory info*/
	number=temp=0;
	for (i=0;i<values->all_done;i++ ) if (a[i]->length>temp) temp=a[i]->length;
	temp*=11;
	temp/=10;
	if (!values->new_optimization) for (i=2;i<cache_classes;i++) number+=i;
	else for (i=2;i<cache_classes;i++) number+=1;
	cache_size=(size/(number*temp*sizeof(char)));
	/*soring only one sequnce not two*/

	/*
fprintf(stderr,"%d cache classes, %d size %d bumber %d temp %d sizeof(char) %d size per\n",cache_classes,size,number, temp, sizeof(char),number*temp*sizeof(char));
*/
	if (cache_size<1) {
		fprintf(stderr,"Specified cache memory is too low--cache cannot be created--increase cache memory or remove the cache from program options.\n");
		exit(-1);
	}
	if (values->VERBOSE) fprintf(stderr,"Cache allocated to hold %d alignments.\n",(cache_classes-2)*cache_size);

	p=(align_buffer **)malloc((values->all_done-1)*sizeof(align_buffer *));
	assert((int)p);
	for (i=0;i<values->all_done-1;i++) {
		p[i]=(align_buffer *)malloc(sizeof(align_buffer));
		assert((int)p[i]);
		p[i]->max_num=cache_size;
		p[i]->num_in=0;
		p[i]->ba=(alignment **)malloc(p[i]->max_num*sizeof(alignment *));
		assert((int)p[i]->ba);
		for (j=0;j<p[i]->max_num;j++) p[i]->ba[j]=NULL;
		p[i]->age=(int *)malloc(p[i]->max_num*sizeof(int));
		assert((int)p[i]->age);
		p[i]->hack=(int *)malloc(p[i]->max_num*sizeof(int));
		assert((int)p[i]->hack);
	}
	for (i=0;i<cache_classes-1;i++) p[i]->hit=p[i]->miss=0;
	return p;
}

alignment *retrieve_align_cache(name,values,holder)
char *name;
align_buffer *holder;
parameters *values;
{
	int i,j,do_it;
	char *p;

	local_hack=0;
	for (i=1;i<strlen(name)-1;i++) local_hack+=(i*name[i]);
	for (i=0;i<holder->num_in;i++) if (local_hack==holder->hack[i]) {
		do_it=1;
		p=holder->ba[i]->name;
		for (j=1;j<strlen(name)-1;j++) if (name[j]!=(*(p+j))) {do_it=0;break;}
		if (do_it) { /*(!strcmp_here(name,holder->ba[i]->name)) {*/
			if (values->cache_info) ++holder->hit;
			for (j=0;j<holder->num_in;j++) ++holder->age[j];
			holder->age[i]=0;
			/*if (values->new_optimization) holder->ba[i]=make_ambig(holder->ba[i],values);*/
			return(make_align(holder->ba[i]));
		}
	}
	if (values->cache_info) ++holder->miss;
	return(NULL);
}


int retrieve_score_cache(a,holder,values)
align_buffer *holder;
alignment *a;
parameters *values;
{
	int i;
	for (i=0;i<holder->num_in;i++) if (!compare_aligns(a,holder->ba[i])) return holder->ba[i]->score;
	return -1;
}

void store_align_cache(a,values,holder)
alignment *a;
align_buffer *holder;
parameters *values;
{
	int i,oldest_i,oldest_age;

	/*if space*/
	if (holder->num_in < holder->max_num ) {
		holder->ba[holder->num_in]=make_align(a);
		for (i=0;i<holder->num_in;i++) ++holder->age[i];
		holder->age[holder->num_in]=0;
		holder->hack[holder->num_in]=0;
		/*for (i=1;i<strlen(a->name)-1;i++) holder->hack[holder->num_in]+=(i*a->name[i]);*/
		 holder->hack[holder->num_in]=local_hack;
		++holder->num_in;
	}
	/*if no space*/
else {
	/*find oldest*/
	oldest_age=0;
	for (i=0;i<holder->num_in;i++) {
		if (holder->age[i]>oldest_age) {
			oldest_age=holder->age[i];
			oldest_i=i;
		}
	}
	/*dump_oldest*/
	holder->ba[oldest_i]=dump_align(holder->ba[oldest_i]);
	/*assign new*/
	holder->ba[oldest_i]=make_align(a);
	/*new age*/
	for (i=0;i<holder->num_in;i++) ++holder->age[i];
	holder->age[oldest_i]=0;
	holder->hack[oldest_i]=0;
	/*for (i=1;i<strlen(a->name)-1;i++) holder->hack[oldest_i]+=(i*a->name[i]);*/
	holder->hack[oldest_i]=local_hack;
	}	
}

void pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs)
int **nodes, *num_seqs;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	char **names,*pre_made,*found;
	int i,n,d1,d2,j,those_to_check;
	alignment **temp;

	/*allocate names*/
	names=(char **)malloc((ntaxa+ntaxa-1)*sizeof(char *));
	assert((int)names);
	pre_made=(char *)malloc((ntaxa+ntaxa-1)*sizeof(char));
	assert((int)pre_made);

	/*optimize names*/
	for (i=0;i<(ntaxa-1);i++) {
		num_seqs[i]=pre_made[i]=1;
		names[i]=(char *)malloc((1+strlen(a[i]->name))*sizeof(char));
		assert((int)names[i]);
		names[i]=strcpy(names[i],a[i]->name);
	}
	names[ntaxa-1]=NULL;
	num_seqs[ntaxa-1]=pre_made[ntaxa-1]=1;
	for (i=ntaxa;i<(ntaxa+ntaxa-1);i++) {
		num_seqs[i]=pre_made[i]=0;
		names[i]=NULL;
	}

tip_top_mt_olympus:
	for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings*/
		d1=nodes[n][0];
		d2=nodes[n][1];
		if ((pre_made[d1]) && (pre_made[d2]) && (!pre_made[ntaxa+n])) {
			pre_made[ntaxa+n]=1;
			/*fix for rooting more taxa*/
			if (d1<ntaxa) d1-=1;
			if (d2<ntaxa) d2-=1;
			/*make names*/
			if (strcmp (names[d1],names[d2]) > 0) names[ntaxa+n]=other_pair_names(names[d2],names[d1]);
			else names[ntaxa+n]=other_pair_names(names[d1],names[d2]);
			num_seqs[ntaxa+n]=num_seqs[d1]+num_seqs[d2];
		}
	}
	for (n=ntaxa+1;n<((ntaxa+nent)-1);++n) if (!pre_made[n]) goto tip_top_mt_olympus;
	/*if (values->new_optimization) for (n=(nent-2);n>0;--n) num_seqs[ntaxa+n]=2;*/
	/*
fprintf(stderr,"%s\n",names[nodes[0][1]]);
*/

	/*check last time*/
	temp=(alignment **)malloc((nent-2)*sizeof(alignment *));
	assert((int)temp);
	found=(char *)malloc((nent-2)*sizeof(char));
	assert((int)found);
	for (i=0;i<nent-2;i++) {
		temp[i]=NULL;
		if (a[ntaxa+i+1]) temp[i]=make_align(a[ntaxa+i+1]);
		found[i]=0;
	}
	for (n=1;n<(nent-1);n++) {
		if (!values->all_made[ntaxa+n]){
			/*check preexiting*/
			for (j=1;j<(nent-1);j++) if (!found[j-1]) {
				if (temp[j-1]) if (!strcmp(temp[j-1]->name,names[ntaxa+n])) {
					found[j-1]=1;
					if (j!=n) {
						a[ntaxa+n]=dump_align(a[ntaxa+n]);
						a[ntaxa+n]=make_align(temp[j-1]);
					}
					values->all_made[ntaxa+n]=1;
					/*convert descendents to "1"*/
					mark_descendants(nodes,ntaxa,nent,values);
					j=nent-1;
				}
			}
		}
	}
	for (i=0;i<nent-2;i++) if (temp[i]) temp[i]=dump_align(temp[i]);
	free(temp);
	free(found);

	/*check cache*/
	/*if (!values->check_cache_for_scores) those_to_check=ntaxa-1;
else those_to_check=ntaxa;
*/
	those_to_check=ntaxa-1;

	/*if (values->new_optimization) for (n=1;n<(nent-1);n++) num_seqs[ntaxa+n]=2;*/


	if (values->align_cache) {
		for (n=1;n<(nent-1);n++) { /*if (num_seqs[ntaxa+n]<those_to_check) {*/
			if (!values->all_made[ntaxa+n]){
				if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
				a[ntaxa+n]=retrieve_align_cache(names[ntaxa+n],values,values->align_cache[num_seqs[ntaxa+n]-2]);
				if (a[ntaxa+n]) {
					values->all_made[ntaxa+n]=1;
					mark_descendants(nodes,ntaxa,nent,values);
				}
			}
		}
	}

	/*deallocate*/
	free(pre_made);
	for (i=0;i<(ntaxa+ntaxa-1);i++) if (names[i]) free(names[i]);
	free(names);
}

void mark_descendants(nodes,ntaxa,nent,values)
int ntaxa,nent;
int **nodes;
parameters *values;
{
	int found, i;

	found=1;

	while (found) {
		found=0;
		for (i=1;i<(nent-1);i++) {
			if (values->all_made[i+ntaxa]) {
				if (!values->all_made[nodes[i][0]]) if (nodes[i][0]>ntaxa) {
					values->all_made[nodes[i][0]]=1;
					found=1;
				}
				if (!values->all_made[nodes[i][1]]) if (nodes[i][1]>ntaxa) {
					values->all_made[nodes[i][1]]=1;
					found=1;
				}
			}
		}
	}

}

int all_diagnose_tree(a,nodes,ntaxa,nent,values)
int **nodes;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0;
	char *temp_name;
	int new_opt_score;
	int *num_seqs;


	/*need to insure that we start with nodes which connect to terminal taxa and then nodes
	which have been defined - otherwise we'll be using the optimizations of 
	previous topology's nodes*/
	/*add in a count for the first nent and nent-1 alignments and number of taxa for entering sequences with !=1 size*/

	num_seqs=(int *)malloc((ntaxa+ntaxa-1)*sizeof(int));
	assert((int)num_seqs);


	if (nent==ntaxa) {
		++values->number;
		/*if ((values->VERBOSE) && (!PARALLEL)) fprintf(stderr,"Perform multiple alignment --> %ld\n",values->number);*/
	}
	values->all_made=(char *)malloc(((2*ntaxa)-1)*sizeof(char));
	assert((int)values->all_made);

	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent-1);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (!values->number_of_input_alignments) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization)) nent_minus_one_new=nent-1;
	else for (n=0;n<=(nent-1);n++) nent_minus_one_new+=a[n]->n_seqs;

	pre_diagnose_align(a,nodes,ntaxa,nent,values,num_seqs);

	/*fprintf(stderr,"%u-",_memavl());*/

	if (!values->new_optimization) {
all_make_the_node:
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				values->all_made[ntaxa+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				/*avoid unnesesary trees*/
				if ((a[d1]->n_seqs+a[d2]->n_seqs)<(nent_minus_one_new)) {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					values->in_some=1;
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*a[ntaxa+n]->type_weight=0;*/
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
					values->in_some=0;
					if (values->in_bandb_loop && (a[ntaxa+n]->score>values->current_bound)) goto second_end;
				}
				else {
					if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
					a[ntaxa+n] = nw(a[d1],a[d2],values);
					/*++temp_done;*/
					if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[a[ntaxa+n]->n_seqs-2]);
				}
			}
		}
		/*fprintf(stderr,"%u ",_memavl());*/
		for (n=ntaxa+1;n<((ntaxa+nent)-1);++n) if (!values->all_made[n]) goto all_make_the_node;
	}
	else {
		new_opt_score=0;
new_all_make_the_node:
		for (n=(nent-2);n>0;--n){ /* don't do to zero because jack up sequence number to get all rootings (as opposed to cladograms)*/
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				values->all_made[ntaxa+n]=1;
				/*fix for rooting more taxa*/
				if (d1<ntaxa) d1-=1;
				if (d2<ntaxa) d2-=1;
				if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
				a[ntaxa+n] = nw(a[d1],a[d2],values);
				/*a[ntaxa+n] = make_ambig(a[ntaxa+n],values);*/
				a[ntaxa+n]->score += (a[d1]->score+a[d2]->score);
				if (values->align_cache) store_align_cache(a[ntaxa+n],values,values->align_cache[num_seqs[ntaxa+n]-2]);
				if (values->in_bandb_loop && (a[ntaxa+n]->score>values->current_bound)) goto second_end;
			}
		}
		for (n=ntaxa+1;n<((ntaxa+nent)-1);++n) if (!values->all_made[n]) goto new_all_make_the_node;
		/*for (n=ntaxa+1;n<((ntaxa+nent)-1);++n) new_opt_score+=a[n]->score;
		a[nodes[0][1]]->score=new_opt_score;*/
	}

	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;

	/*fprintf(stderr,"%u)",_memavl());*/

	return(a[nodes[0][1]]->score);

	/*if early exceed bound*/
second_end:
	free(num_seqs);
	free(values->all_made);
	values->all_made=NULL;
	return HUGE_COST;
}

char *other_pair_names (a1, a2)
char *a1, *a2;
{
	char *p;

	p = (char *) allocate (strlen (a1) + strlen (a2) +5);
	sprintf (p, "(%s %s)", a1, a2);
	return p;
}

int strcmp_here(a,b)
char *a, *b;
{
/*
while (a && b) {
	if ((*a)!=(*b)) return 1;
	a++;b++;
	}
	*/
	int i;
	for (i=0;i<strlen(a);i++) if (a[i]!=b[i]) return 1;
return 0;
}
