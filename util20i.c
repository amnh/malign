/*Copyright 1993 Ward Wheeler all rights reserved*/
/*Make sure that the reconstructed alignments == the initial ones + can count gaps and changes from down pass as a check*/

#include "align3.h"
extern int charbit[256];
extern char bitchar[32];

int num_states(k)
int k;
{
/*int rv;

rv=0;
if (A_BIT & k) ++rv;
if (C_BIT & k)  ++rv;
if (G_BIT & k) ++rv;
if (T_BIT & k) ++rv;
if (GAP_BIT & k) ++rv;
return rv;*/

switch (k) {
	case A_BIT: return 1;
	case C_BIT: return 1;
	case G_BIT: return 1;
	case T_BIT: return 1;
	case R_BIT: return 2;
	case Y_BIT: return 2;
	case M_BIT: return 2;
	case K_BIT: return 2;
	case W_BIT: return 2;
	case S_BIT: return 2;
	case B_BIT: return 3;
	case D_BIT: return 3;
	case H_BIT: return 3;
	case V_BIT: return 3;
	}
return 4;
}

int **get_groups3(buffer,ntax,a,values)
char *buffer;
int ntax;
alignment **a;
parameters *values;
{
int i,j,*buffer2,paren_check_left,paren_check_right;
int rep_length,groups_marker,groups_number;
int **groups,ngroups;

printf(" %s\n",buffer);

paren_check_left=paren_check_right=0;
for (i=0;i<strlen(buffer);i++) {
	if (buffer[i]=='(') ++paren_check_left;
	else if (buffer[i]==')') ++paren_check_right;
	}
if (paren_check_left!=paren_check_right) {
	fprintf(stderr,"Parentheses do not match in groups specification.\n");
	exit(-1);
	}

buffer2=(int *)malloc((strlen(buffer)+1)*sizeof(int));
assert((int)buffer2);
for (i=0;i<(strlen(buffer)+1);i++) buffer2[i]=-32000;
j=0;
for (i=0;i<strlen(buffer);i++)
	{
	if (buffer[i]=='(') buffer2[j++]=-1;
	else if (buffer[i]==')') buffer2[j++]=-2;
	else if (buffer[i]==',') buffer2[j++]=-3;
	else if (isspace(buffer[i])) buffer2[j++]=-3; /*can get rid of extra commas and terminal spaces laterif need to */
	else if (isalnum(buffer[i])) {
		get_taxon_number(&i,&j,buffer,buffer2,a,ntax);
		--i;
		}
	else {
		fprintf(stderr,"Bad character |%c| in position %d of groups description.\n",buffer[i],i);
		exit(-1);
		}
	}
buffer2[i]=-4;
i=0;
while (buffer2[i]>-4) {
	if (buffer2[i]==-2) rep_length=i+1;
	i++;
	}
/*for (i=0;i<rep_length;i++) fprintf(stderr,"%d,",buffer2[i]);
fprintf(stderr,"\n");*/

/*the real stuff*/
ngroups=paren_check_left;
if (ngroups!=ntax-1) {
	fprintf(stderr,"Something is wrong in the name %s\n",buffer);
	exit(-1);
	}
groups=(int **)malloc(ngroups*sizeof(int *));
assert((int)groups);
for (i=0;i<paren_check_left;i++) {
	/*groups[i]=(int *)malloc((ntax+1)*sizeof(int));*/
	groups[i]=(int *)malloc(((2*ntax))*sizeof(int));
	assert((int)groups[i]);
	for (j=0;j<ntax+1;j++) groups[i][j]=0;
	}

groups_number=0;
for (i=0;i<rep_length;i++) {
	if (buffer2[i]==-1) {
		groups_marker=0;
		for (j=i;j<rep_length;j++) {
			if (buffer2[j]==-1) ++groups_marker;
			if ((buffer2[j]>=0)&&(buffer2[j]<ntax)) groups[groups_number][buffer2[j]]=1;
			if (buffer2[j]==-2) {
				--groups_marker;
				if (groups_marker==0) break;
				}
			}
		++groups_number;
		}
	}

for (i=0;i<ngroups;i++) {
	for (j=0;j<ntax;j++) groups[i][ntax]+=groups[i][j];
	if (groups[i][ntax]==0) {
		fprintf(stderr,"Error in groups specification -- parenthesis pair which contains no sequences (%d).\n",j);
		exit(-1);
		}
/*	for (j=0;j<=ntax;j++) fprintf(stderr," %d ",groups[i][j]);
	fprintf(stderr,"\n");*/
	}
/*fprintf(stderr,"\n");*/

free(buffer2);
return groups;
}

void do_and_print_optimization(a,values,count,best_aligns)
alignment **a;
parameters *values;
int count;
alignment **best_aligns;
{
int i,j,k,ntax, new_score;
int **cur_groups, found_a_two,desc1,desc2;
int **nodes,current_node;

ntax=values->all_done;
nodes=(int **)malloc(ntax*sizeof(int *));
assert((int) nodes);
for (i=0;i<ntax;i++) {
	nodes[i]=(int *)malloc(2*sizeof(int));
	assert((int) nodes[i]);
	nodes[i][0]=nodes[i][1]=-1;
	}
fprintf(stderr,"	Topology");
for (i=0;i<count;i++) {
	if (i>0) fprintf(stderr,"	        ");
	fprintf(stderr," %d\n",i);
	printf("   \nTopology %d:",i);

	/*get initial groups*/
	cur_groups=get_groups3(best_aligns[i]->name,ntax,a,values);

	/*reallocate to include nodes*/
	for (j=0;j<ntax-1;j++) {
		/*cur_groups[j]=(int *)realloc(cur_groups[j],(2*ntax)*sizeof(int));
		assert((int) cur_groups[j]);*/
		cur_groups[j][(2*ntax)-1]=cur_groups[j][ntax];
		for (k=ntax;k<((2*ntax)-1);k++) cur_groups[j][k]=0;
		/*for (k=0;k<(2*ntax);k++) fprintf(stderr," %d ",cur_groups[j][k]);
		fprintf(stderr,"\n");*/
	}

	/*look for twos/ remove/ return*/
	found_a_two=1;
	current_node=0;
	while (found_a_two) {
		found_a_two=0;
		for (j=0;j<ntax-1;j++) {
			if (cur_groups[j][2*ntax-1]==2) {
				/*is a node*/
				found_a_two=1;
				++current_node;
				desc1=-1;
				desc2=-1;
				/*get descendents*/
				for (k=0;k<((2*ntax)-1);k++) {
					if (cur_groups[j][k]==1) {
						if (desc1<0) desc1=k;
						else desc2=k;
						}
					}
				/*mark nodes*/
				nodes[current_node][0]=desc1;
				nodes[current_node][1]=desc2;
				cur_groups[j][desc1]=0;
				cur_groups[j][desc2]=0;
				cur_groups[j][2*ntax-1]=0;
				/*update others*/
				for (k=0;k<ntax-1;k++) if (k!=j) {
					if ((cur_groups[k][desc1]==1) && (cur_groups[k][desc2]==1)) {
						cur_groups[k][desc1]=0;
						cur_groups[k][desc2]=0;
						cur_groups[k][current_node+ntax]=1;
						cur_groups[k][2*ntax-1]-=1;
						}
					}
				}
			}
		}

	/*add root and finish nodes*/
	nodes[0][1]=ntax+ntax-1;
	/*for (j=0;j<ntax;j++) fprintf(stderr,"Nodes %d %d\n",nodes[j][0],nodes[j][1]);
	fprintf(stderr,"\n");*/

	/*dump groups*/
	for (j=0;j<ntax-1;j++) free(cur_groups[j]);
	free(cur_groups);

	/*up_pass*/
	new_score=all_diagnose_tree_w_uppass(a,nodes,ntax,ntax,values);
	if (new_score!=best_aligns[i]->score) {
		fprintf(stderr,"Error optimizing topology %d\n",i);
		}
	}/*next topology*/
fprintf(stderr,"\n");
/*free*/
for (i=0;i<ntax;i++) free(nodes[i]);
free(nodes);
}

int all_diagnose_tree_w_uppass(a,nodes,ntaxa,nent,values)
int **nodes;
int ntaxa,nent;
alignment **a;
parameters *values;
{
	int n,d1,d2,nent_minus_one_new=0,i,j;
	char *temp_name;
	int new_opt_score, **branch_lengths;
	int *num_seqs, *anc;
	alignment *temp_align,*temp_align2;
	int ancestor, descendent;
	apo_thang ***apomorphy;
	apolist_holder *return_holder;
	int check_score,n_gaps,n_changes,n_lt;
	int ***corres, d1_thang, d2_thang;
	alignment **a_holder,*temp_align3,*temp_align4,*temp_align5;
	int c_h,g_h,desc1,desc2;

	/*need to insure that we start with nodes which connect to terminal taxa and then nodes
	which have been defined - otherwise we'll be using the optimizations of
	previous topology's nodes*/
	/*add in a count for the first nent and nent-1 alignments and number of taxa for entering sequences with !=1 size*/

	a_holder=(alignment **)malloc((2*ntaxa)*sizeof(alignment *));
	assert((int) a_holder);
	corres=(int ***)malloc((2*ntaxa)*sizeof(int **));
	assert((int) corres);
	return_holder=(apolist_holder *)malloc(sizeof(apolist_holder));
	assert((int) return_holder);
	if (values->all_made) free (values->all_made);
	values->all_made=(char *)malloc(((2*ntaxa))*sizeof(char));
	assert((int)values->all_made);
	anc=(int *)malloc(((2*ntaxa))*sizeof(int));
	assert((int)anc);
	branch_lengths=(int **)malloc(((2*ntaxa))*sizeof(int *));
	assert((int)branch_lengths);
	for (n=0;n<((2*ntaxa));n++) {
		branch_lengths[n]=(int *)malloc(4*3*sizeof(int)); /*this is 3*3 lt/i/c/cost X min/max/reg because of four types of transformation--could be increased later*/
		assert((int) branch_lengths[n]);
		}
	if (values->apolist) {
		apomorphy=(apo_thang ***)malloc((2*ntaxa)*sizeof(apo_thang **));
		assert((int)apomorphy);
		for (n=0;n<(2*ntaxa);n++) {
			apomorphy[n]=(apo_thang **)malloc(3*sizeof(apo_thang *));/*leading/rtainling internal changes*/
			assert((int)apomorphy[n]);
			for (i=0;i<3;i++) apomorphy[n][i]=NULL;
			}
		}
	else apomorphy=NULL;


	for (n=0;n<ntaxa;n++) values->all_made[n]=1;
	for (n=ntaxa;n<(ntaxa+nent);++n) values->all_made[n]=0;

	/*make sure the in-some check works for multiple input alignments*/
	if (!values->number_of_input_alignments) nent_minus_one_new=nent-1;
	else if ((!values->new_optimization) && (values->chop>0)) nent_minus_one_new=nent-1;
	else for (n=0;n<nent;n++) nent_minus_one_new+=a[n]->n_seqs;
	/*down pass redone*/
		values->in_optimize=1;
		new_opt_score=0;
		new_all_make_the_node2:
		for (n=(nent-1);n>0;--n){
			d1=nodes[n][0];
			d2=nodes[n][1];
			if ((values->all_made[d1]) && (values->all_made[d2]) && (!values->all_made[ntaxa+n])) {
				values->all_made[ntaxa+n]=1;
				if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
				/*get alignment*/
				a[ntaxa+n] = nw(a[d1],a[d2],values);
				temp_align=make_align(a[ntaxa+n]);
				/*get ancestor with gaps*/
				temp_align=new_make_ambig(temp_align,values);
				/*print changes and gaps to check the uppass*/
				g_h=c_h=0;
				for (i=0;i<a[ntaxa+n]->length;i++) if (!overlap(a[ntaxa+n]->s[0][i],a[ntaxa+n]->s[1][i])) {
					if ((a[ntaxa+n]->s[0][i]=='-') && ((a[ntaxa+n]->s[1][i]!='-'))) ++g_h;
					else if ((a[ntaxa+n]->s[1][i]=='-') && ((a[ntaxa+n]->s[0][i]!='-'))) ++g_h;
					else ++c_h;
					}
				/*fprintf(stderr,"Changes from %d: %d gaps %d changes ",ntaxa+n,g_h,c_h);*/
				/*get correspondance arrays*/
				if (strcmp(a[d1]->name,a[d2]->name)>0) {
					d2_thang=0;
					d1_thang=1;
					}
				else {
					d2_thang=1;
					d1_thang=0;
					}
				corres[d1]=get_correspondances(temp_align->s[0],a[ntaxa+n]->s[d1_thang],temp_align->length);
				corres[d2]=get_correspondances(temp_align->s[0],a[ntaxa+n]->s[d2_thang],temp_align->length);
				/*get ancestor without gaps*/
			    	a[ntaxa+n]=make_ambig(a[ntaxa+n],values);
			  	a[ntaxa+n]->score += (a[d1]->score+a[d2]->score);
			  	a_holder[ntaxa+n]=make_align(a[ntaxa+n]);
			  	temp_align=dump_align(temp_align);
	     		}
			}
		for (n=ntaxa+1;n<((ntaxa+nent));++n) if (!values->all_made[n]) goto new_all_make_the_node2;

/*assigned branch lengths first*/
/*up_pass part*/
	/*get_ancestors*/
	for (n=1;n<nent;n++) anc[nodes[n][0]]=anc[nodes[n][1]]=n;
	anc[nodes[0][1]]=0;/*fix this for outtaxa?*/
	printf("Reconstructed Ancestral Nodes\n");

	for (n=ntaxa+1;n<(ntaxa+nent);++n) values->all_made[n]=0;
	/*get ancestors of each node*/
	/*reset basal node first*/
	values->all_made[ntaxa]=1;
	values->all_made[nodes[0][1]]=1;

	up_pass_mm:
	for (n=1;n<nent;n++){
	if ((!values->all_made[n+ntaxa]) && (values->all_made[anc[n+ntaxa]+ntaxa])) {
			values->all_made[ntaxa+n]=1;
		/*add fake nw here to regenerate the initale alignment correspondances*/
			temp_align=fake_nw(a[ntaxa+n],a[anc[n+ntaxa]+ntaxa],values,corres[ntaxa+n]);
			ancestor=1;
			descendent=0;
			if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
			/*change make_anc_up to reflect new way with conditions for counting
				and leave in the damn gaps or the correspondaces will not hold for the next path*/
			return_holder=make_anc_up(temp_align,values,branch_lengths[ntaxa+n]+4,ancestor,descendent,0,apomorphy[ntaxa+n],ntaxa+n,return_holder);
			a[ntaxa+n]=return_holder->z;
			return_holder=make_anc_up(temp_align,values,branch_lengths[ntaxa+n]+8,ancestor,descendent,3,apomorphy[ntaxa+n],ntaxa+n,return_holder);
			temp_align2=return_holder->z;
			if (temp_align) temp_align=dump_align(temp_align);
			if (temp_align2) temp_align2=dump_align(temp_align2);
			}
		}
	for (n=ntaxa+1;n<=((ntaxa+nent)-1);++n) if (!values->all_made[n]) goto up_pass_mm;
/*then do nodes to terminals*/
for (n=0;n<ntaxa;n++) {
	ancestor=1;
	descendent=0;
	temp_align=fake_nw(a[n],a[anc[n]+ntaxa],values,corres[n]);
	return_holder=make_anc_up(temp_align,values,branch_lengths[n]+4,ancestor,descendent,0,apomorphy[n],n,return_holder); /*min lengths*/
	temp_align2=return_holder->z;
	return_holder=make_anc_up(temp_align,values,branch_lengths[n]+8,ancestor,descendent,3,apomorphy[n],n,return_holder);
	temp_align2=return_holder->z;
	if (temp_align) temp_align=dump_align(temp_align);
	if (temp_align2) temp_align2=dump_align(temp_align2);
	}
printf("With possible states as ambiguities\n");
for (n=ntaxa+1;n<=((ntaxa+nent)-1);++n) {
	printf("	Node %3d: ",n);
	for (i=0;i<a[n]->length;i++) if (a[n]->s[0][i]!='-') printf("%c",a[n]->s[0][i]);
	printf("\n");
	}
printf("\n");

/*reset nodes*/
for (n=(nent-1);n>0;--n) {
	if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
	a[ntaxa+n]=make_align(a_holder[ntaxa+n]);
	a_holder[ntaxa+n]=dump_align(a_holder[ntaxa+n]);
	}
free(a_holder);
for (n=ntaxa+1;n<(ntaxa+nent);++n) values->all_made[n]=0;

/*set each each node*/
	/*reset basal node first*/
	values->all_made[ntaxa]=1;
	values->all_made[nodes[0][1]]=1;
	/*set outgroup*/
	/*A second way*/
		/*for (n=0;n<a[nodes[0][1]]->length;n++) a[nodes[0][1]]->s[0][n]=overlap(a[nodes[0][1]]->s[0][n],a[nodes[nodes[0][1]-ntaxa][0]]->s[0][n]);*/
		for (n=0;n<a[nodes[0][1]]->length;n++) a[nodes[0][1]]->s[0][n]=choose_one2(a[nodes[0][1]]->s[0][n],values->rand_apo);

	up_pass_mm_1:
	for (n=1;n<nent;n++){
	if ((!values->all_made[n+ntaxa]) && (values->all_made[anc[n+ntaxa]+ntaxa])) {
			values->all_made[ntaxa+n]=1;
		/*add fake nw here to regenerate the initale alignment correspondances*/
			temp_align=fake_nw(a[ntaxa+n],a[anc[n+ntaxa]+ntaxa],values,corres[ntaxa+n]);
			ancestor=1;
			descendent=0;
			temp_align3=fake_nw(a[nodes[n][0]],a[ntaxa+n],values,corres[nodes[n][0]]);
			temp_align4=fake_nw(a[nodes[n][1]],a[ntaxa+n],values,corres[nodes[n][1]]);
			temp_align5=make_neighbor(temp_align,temp_align3,temp_align4,values);
			if (a[ntaxa+n]) a[ntaxa+n]=dump_align(a[ntaxa+n]);
			/*change make_anc_up to reflect new way with conditions for counting
				and leave in the damn gaps or the correspondaces will not hold for the next path*/
			desc1=1;
			desc2=2;
			return_holder=make_anc_up2(temp_align5,values,branch_lengths[ntaxa+n],ancestor,descendent,desc1,desc2,1,apomorphy[ntaxa+n],ntaxa+n,return_holder);
			a[ntaxa+n]=return_holder->z;
			if (temp_align) temp_align=dump_align(temp_align);
			if (temp_align3) temp_align3=dump_align(temp_align3);
			if (temp_align4) temp_align4=dump_align(temp_align4);
			if (temp_align5) temp_align5=dump_align(temp_align5);
			}
		}
	for (n=ntaxa+1;n<=((ntaxa+nent)-1);++n) if (!values->all_made[n]) goto up_pass_mm_1;

printf("With states assigned\n");
for (n=ntaxa+1;n<=((ntaxa+nent)-1);++n) {
	printf("	Node %3d: ",n);
	for (i=0;i<a[n]->length;i++) if (a[n]->s[0][i]!='-') printf("%c",a[n]->s[0][i]);
	printf("\n");
	}
printf("\n");

/*Now count changes*/
/*then do nodes to terminals*/
for (n=0;n<ntaxa;n++) {
	ancestor=1;
	descendent=0;
	temp_align=fake_nw(a[n],a[anc[n]+ntaxa],values,corres[n]);
	return_holder=make_anc_up(temp_align,values,branch_lengths[n],ancestor,descendent,1,apomorphy[n],n,return_holder); /*min lengths*/
	temp_align2=return_holder->z;
	temp_align=dump_align(temp_align);
	temp_align2=dump_align(temp_align2);
	}

for (n=(nent-1);n>0;--n){
	d1=nodes[n][0];
	d2=nodes[n][1];
	free(corres[d1][0]);
	free(corres[d2][0]);
	free(corres[d1][1]);
	free(corres[d2][1]);
	}
free(corres);
check_score=n_gaps=n_changes=n_lt=0;
for (n=1;n<(ntaxa);n++) if ((anc[n+ntaxa]+ntaxa)!=ntaxa) {
	check_score+=(branch_lengths[n+ntaxa][0]*(values->gap_cost - values->leading_gap_cost));
	check_score+=(branch_lengths[n+ntaxa][1]*(values->gap_cost));
	check_score+=(branch_lengths[n+ntaxa][2]*(values->change_cost));
	n_gaps+=branch_lengths[n+ntaxa][1];
	n_changes+=branch_lengths[n+ntaxa][2];
	n_lt+=branch_lengths[n+ntaxa][0];
	}
for (n=0;n<(ntaxa);n++) {
	check_score+=(branch_lengths[n][0]*(values->gap_cost - values->leading_gap_cost));
	check_score+=(branch_lengths[n][1]*(values->gap_cost));
	check_score+=(branch_lengths[n][2]*(values->change_cost));
	n_gaps+=branch_lengths[n][1];
	n_changes+=branch_lengths[n][2];
	n_lt+=branch_lengths[n][0];
	}
printf("	Branch Lengths (This = One Optimization of (potentially) many)\n");
printf("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
printf("				%-20s						%-20s\n"," ","Insertions/Deletions");
printf("				%-20s			----------------------------------------------------------\n"," ");
printf("	  	%-10s		%-20s	%-20s						%-10s							%-10s							%-10s\n","Branch"," ","Leading/Trailing","Internal","Changes","Cost");
printf("	-----------------%-20s	-------------------------------------		-------------------------------------		-------------------------------------		-------------------------------------\n"," ");
printf("	%-10s	%-20s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s	%-10s\n","From","To","This","Minimum","Maximum","This","Minimum","Maximum","This","Minimum","Maximum","This","Minimum","Maximum");
printf("-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
for (n=1;n<(ntaxa);n++) if ((anc[n+ntaxa]+ntaxa)!=ntaxa) {
	if (anc[n+ntaxa]+ntaxa!=((2*ntaxa)-1))
printf("	%-10d	%-20d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d\n",anc[n+ntaxa]+ntaxa,n+ntaxa,branch_lengths[n+ntaxa][0],branch_lengths[n+ntaxa][4],branch_lengths[n+ntaxa][8],branch_lengths[n+ntaxa][1],branch_lengths[n+ntaxa][5],branch_lengths[n+ntaxa][9],branch_lengths[n+ntaxa][2],branch_lengths[n+ntaxa][6],branch_lengths[n+ntaxa][10],branch_lengths[n+ntaxa][3],branch_lengths[n+ntaxa][7],branch_lengths[n+ntaxa][11]);
	else
printf("	%-10s	%-20d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d\n","root",n+ntaxa,branch_lengths[n+ntaxa][0],branch_lengths[n+ntaxa][4],branch_lengths[n+ntaxa][8],branch_lengths[n+ntaxa][1],branch_lengths[n+ntaxa][5],branch_lengths[n+ntaxa][9],branch_lengths[n+ntaxa][2],branch_lengths[n+ntaxa][6],branch_lengths[n+ntaxa][10],branch_lengths[n+ntaxa][3],branch_lengths[n+ntaxa][7],branch_lengths[n+ntaxa][11]);
	}
for (n=0;n<ntaxa;n++) {
	if (anc[n]+ntaxa!=((2*ntaxa)-1))
printf("	%-10d	%-20s	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d\n",anc[n]+ntaxa,a[n]->name,branch_lengths[n][0],branch_lengths[n][4],branch_lengths[n][8],branch_lengths[n][1],branch_lengths[n][5],branch_lengths[n][9],branch_lengths[n][2],branch_lengths[n][6],branch_lengths[n][10],branch_lengths[n][3],branch_lengths[n][7],branch_lengths[n][11]);
	else
printf("	%-10s	%-20s	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d	%-10d\n","root",a[n]->name,branch_lengths[n][0],branch_lengths[n][4],branch_lengths[n][8],branch_lengths[n][1],branch_lengths[n][5],branch_lengths[n][9],branch_lengths[n][2],branch_lengths[n][6],branch_lengths[n][10],branch_lengths[n][3],branch_lengths[n][7],branch_lengths[n][11]);
	}
printf("\nTotal length %d with %d leading/trailing gaps, %d internal gaps, and %d changes\n",check_score,n_lt,n_gaps,n_changes);
if (values->apolist) {
	printf("\nApomorphy list (This is only one optimization of (potentially) many)\n");
	printf("----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	printf("		[ancestral position] ancestral character state -> [descendant position] descendant character state\n");
	for (n=1;n<(ntaxa);n++) if ((anc[n+ntaxa]+ntaxa)!=ntaxa) {
		printf("	Branch leading to %d\n",ntaxa+n);
		printf("		Leading/trailing indels (%d):\n",branch_lengths[ntaxa+n][0]);
		for (i=0;i<branch_lengths[ntaxa+n][0];i++) {
			printf("			[%3d] %c -> [%3d] %c ",apomorphy[ntaxa+n][0][i].anc_pos,apomorphy[ntaxa+n][0][i].anc_char,apomorphy[ntaxa+n][0][i].desc_pos,apomorphy[ntaxa+n][0][i].desc_char);
			if (apomorphy[ntaxa+n][0][i].desc_char=='-') printf("(deletion)\n");
			else printf("(insertion)\n");
			}
		printf("		Internal indels (%d):\n",branch_lengths[ntaxa+n][1]);
		for (i=0;i<branch_lengths[ntaxa+n][1];i++) {
			printf("			[%3d] %c -> [%3d] %c ",apomorphy[ntaxa+n][1][i].anc_pos,apomorphy[ntaxa+n][1][i].anc_char,apomorphy[ntaxa+n][1][i].desc_pos,apomorphy[ntaxa+n][1][i].desc_char);
			if (apomorphy[ntaxa+n][1][i].desc_char=='-') printf("(deletion)\n");
			else printf("(insertion)\n");
			}
		printf("		Changes (%d):\n",branch_lengths[ntaxa+n][2]);
		for (i=0;i<branch_lengths[ntaxa+n][2];i++) {
			printf("			[%3d] %c -> [%3d] %c ",apomorphy[ntaxa+n][2][i].anc_pos,apomorphy[ntaxa+n][2][i].anc_char,apomorphy[ntaxa+n][2][i].desc_pos,apomorphy[ntaxa+n][2][i].desc_char);
			if (((apomorphy[ntaxa+n][2][i].anc_char=='A') || (apomorphy[ntaxa+n][2][i].anc_char=='G')) && ((apomorphy[ntaxa+n][2][i].desc_char=='C') || (apomorphy[ntaxa+n][2][i].desc_char=='T'))) printf("(transversion)");
			else if (((apomorphy[ntaxa+n][2][i].desc_char=='A') || (apomorphy[ntaxa+n][2][i].desc_char=='G')) && ((apomorphy[ntaxa+n][2][i].anc_char=='C') || (apomorphy[ntaxa+n][2][i].anc_char=='T'))) printf("(transversion)");
			else if (((apomorphy[ntaxa+n][2][i].desc_char=='A') || (apomorphy[ntaxa+n][2][i].desc_char=='G')) && ((apomorphy[ntaxa+n][2][i].anc_char=='A') || (apomorphy[ntaxa+n][2][i].anc_char=='G'))) printf("(transition)");
			else if (((apomorphy[ntaxa+n][2][i].desc_char=='C') || (apomorphy[ntaxa+n][2][i].desc_char=='T')) && ((apomorphy[ntaxa+n][2][i].anc_char=='C') || (apomorphy[ntaxa+n][2][i].anc_char=='T'))) printf("(transition)");
			printf("\n");
			}
		}
	for (n=0;n<ntaxa;n++) {
		printf("	Branch leading to %s\n",a[n]->name);
		printf("		Leading/trailing indels (%d):\n",branch_lengths[n][0]);
		for (i=0;i<branch_lengths[n][0];i++) {
			printf("			[%3d] %c -> [%3d] %c ",apomorphy[n][0][i].anc_pos,apomorphy[n][0][i].anc_char,apomorphy[n][0][i].desc_pos,apomorphy[n][0][i].desc_char);
			if (apomorphy[n][0][i].desc_char=='-') printf("(deletion)\n");
			else printf("(insertion)\n");
			}
		printf("		Internal indels (%d):\n",branch_lengths[n][1]);
		for (i=0;i<branch_lengths[n][1];i++) {
			printf("			[%3d] %c -> [%3d] %c ",apomorphy[n][1][i].anc_pos,apomorphy[n][1][i].anc_char,apomorphy[n][1][i].desc_pos,apomorphy[n][1][i].desc_char);
			if (apomorphy[n][1][i].desc_char=='-') printf("(deletion)\n");
			else printf("(insertion)\n");
			}
		printf("		Changes (%d):\n",branch_lengths[n][2]);
		for (i=0;i<branch_lengths[n][2];i++) {
			printf("			[%3d] %c -> [%3d] %c ",apomorphy[n][2][i].anc_pos,apomorphy[n][2][i].anc_char,apomorphy[n][2][i].desc_pos,apomorphy[n][2][i].desc_char);
			if (((apomorphy[n][2][i].anc_char=='A') || (apomorphy[n][2][i].anc_char=='G')) && ((apomorphy[n][2][i].desc_char=='C') || (apomorphy[n][2][i].desc_char=='T'))) printf("(transversion)");
			else if (((apomorphy[n][2][i].desc_char=='A') || (apomorphy[n][2][i].desc_char=='G')) && ((apomorphy[n][2][i].anc_char=='C') || (apomorphy[n][2][i].anc_char=='T'))) printf("(transversion)");
			else if (((apomorphy[n][2][i].desc_char=='A') || (apomorphy[n][2][i].desc_char=='G')) && ((apomorphy[n][2][i].anc_char=='A') || (apomorphy[n][2][i].anc_char=='G'))) printf("(transition)");
			else if (((apomorphy[n][2][i].desc_char=='C') || (apomorphy[n][2][i].desc_char=='T')) && ((apomorphy[n][2][i].anc_char=='C') || (apomorphy[n][2][i].anc_char=='T'))) printf("(transition)");
			printf("\n");
			}
		}
	}

free(anc);
free(values->all_made);
values->all_made=NULL;
for (n=0;n<((2*ntaxa));n++) free(branch_lengths[n]);
free(branch_lengths);
if (values->apolist) if (apomorphy) {
	for (n=0;n<(2*ntaxa);n++) if (apomorphy[n]) {
		for (j=0;j<3;j++) if (apomorphy[n][j]) free(apomorphy[n][j]);
		free(apomorphy[n]);
		}
	free(apomorphy);
	}
free(return_holder);
return check_score;
return(a[nodes[0][1]]->score);
}



/*this will always produce a base even if a gap is possible because of the assumption of gap cost > change cost*/
int choose_one(anc_state,rand,values)
parameters *values;
int anc_state,rand;
{
int choice,num_states,counter;

if (!rand) {
	if (values->gap_cost < values->change_cost) if (anc_state&GAP_BIT) return GAP_BIT;
	if (anc_state&A_BIT) return A_BIT;
	if (anc_state&C_BIT) return C_BIT;
	if (anc_state&G_BIT) return G_BIT;
	if (anc_state&T_BIT) return T_BIT;
	if (anc_state&GAP_BIT) return GAP_BIT;
	fprintf(stderr,"This cannot happen in optimization\n");
	exit(-1);
	}
else {
	num_states=0;
	if (anc_state&A_BIT) ++num_states;
	if (anc_state&C_BIT) ++num_states;
	if (anc_state&G_BIT) ++num_states;
	if (anc_state&T_BIT) ++num_states;
	if (anc_state&GAP_BIT) ++num_states;
	if (num_states==0) {
		fprintf(stderr,"This cannot happen in optimization\n");
		exit(-1);
		}
	else if (num_states==1) return anc_state;
	choice=my_randomize(num_states-1); /*to include gaps change to unm_states*/
	counter=0;
	if (anc_state&A_BIT) {
		if (counter==choice) return A_BIT;
		else ++counter;
		}
	if (anc_state&C_BIT) {
		if (counter==choice) return C_BIT;
		else ++counter;
		}
	if (anc_state&G_BIT) {
		if (counter==choice) return G_BIT;
		else ++counter;
		}
	if (anc_state&T_BIT) {
		if (counter==choice) return T_BIT;
		else ++counter;
		}
	if (anc_state&GAP_BIT) {
		if (counter==choice) return GAP_BIT;
		else ++counter;
		}
	}
return -1;
}


/*return random 0 -> n-1*/
int my_randomize(n)
int n;
{
int rnum,divisor,i;

if (n==0) {
	fprintf(stderr,"Dividing by zero in my_randomize\n");
	fprintf(stderr,"This could not have happened\n");
	exit(-1);
	}
rnum=rand();
divisor=max_rand/n;

if (rnum <= divisor) return 0;
for (i=1;i<n-1;i++) {
	if ((rnum > (i*divisor)) && (rnum <= ((i+1)*divisor))) return i;
	}
return (n-1);
}


int do_and_print_optimization2(a,values,best_aligns,cur_nodes)
alignment **a;
parameters *values;
alignment *best_aligns;
int **cur_nodes;
{
int i,j,k,ntax,*anc,check,d1,d2;
alignment **a_buf,**a_down2,*temp_align;
int min_changes[4];/* subs ti tv indels*/
int max_changes[4],c1,c2;

ntax=values->all_done+1;



temp_align=NULL;
/*down and up pass */
a_down2=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_down2);
for (i=0;i<2*ntax-1;i++) a_down2[i]=NULL;
a_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_buf);
for (i=0;i<2*ntax-1;i++) a_buf[i]=NULL;
for (i=0;i<ntax-1;i++) a_buf[i]=make_align(a[i]);

anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);
/*
fprintf(stderr,"Before cur\n");
for (i=0;i<ntax-1;i++) for (j=0;j<2;j++) fprintf(stderr,"CN[%d][%d]=%d\n",i,j,cur_nodes[i][j]);*/
/*get support for the new tree*/
anc[ntax]=ntax;
for (j=0;j<ntax-1;j++) anc[cur_nodes[j][0]]=anc[cur_nodes[j][1]]=j+ntax;

check=all_diagnose_tree_here2(a,cur_nodes,ntax,ntax,values,HUGE_COST,a_down2);
get_up_pass_new_opt(a_buf,a,cur_nodes,ntax,values,ntax,anc,a_down2);
printf("HYANC sequences:\n");
for (i=1;i<ntax-1;i++) printf("Node %d %s %s\n",i+ntax,a[i+ntax]->name,a_buf[i+ntax]->s[0]);
for (i=1;i<ntax-1;i++) {
   values->new_optimization=0;
  d1=cur_nodes[i][0];
  d2=cur_nodes[i][1];
  if (d1 < ntax) d1 -=1;
  if (d2 < ntax) d2 -=1;
 /* if ((d1 > ntax) && (d2 > ntax)) fprintf(stderr,"%d -> %d %d ",i+ntax,d1,d2);
  if ((d1 < ntax) && (d2 < ntax)) fprintf(stderr,"%d -> %s %s ",i+ntax,a[d1]->name,a[d2]->name);
  else if (d1 < ntax) fprintf(stderr,"%d -> %s %d ",i+ntax,a[d1]->name,d2);
  else if (d2 < ntax) fprintf(stderr,"%d -> %d %s ",i+ntax,d1,a[d2]->name);*/

  temp_align=nw(a_buf[ntax+i],a_buf[d1],values);
  for (j=0;j<4;j++) min_changes[j]=max_changes[j]=0;
 /*for (j=0;j<2;j++) printf("%s\n",temp_align->s[j]);*/
  printf("To ");
  if (1) /*(d1 < ntax )*/  printf("%s ",a[d1]->name);
  else printf("%d ", d1);
  /*get min changes*/
  for (k=0;k<temp_align->length;k++) if (! (charbit[temp_align->s[0][k]] & charbit[temp_align->s[1][k]])) {
    /*printf("(%c %c)",temp_align->s[0][k],temp_align->s[1][k]);*/
    if ((temp_align->s[0][k]=='-') && (!(charbit[temp_align->s[1][k]] & GAP_BIT))) min_changes[3]+=1;
    else if ((!(charbit[temp_align->s[0][k]] & GAP_BIT))  && (temp_align->s[1][k]=='-')) min_changes[3]+=1;
    else {
      /*printf("(%c %c)",temp_align->s[0][k],temp_align->s[1][k]);*/
      min_changes[0]+=1; /*if not indel then must have a base change if must be change*/
      if ((!(charbit[temp_align->s[0][k]] & R_BIT)) && (!(charbit[temp_align->s[1][k]] & Y_BIT))) min_changes[2]+=1;
      else if ((!(charbit[temp_align->s[0][k]] & Y_BIT)) && (!(charbit[temp_align->s[1][k]] & R_BIT)))  min_changes[2]+=1;
      else  min_changes[1]+=1; /*transition*/
    }
  }
  printf("Min Changes %d %d %d %d ", min_changes[0],min_changes[1],min_changes[2],min_changes[3]);
  /*Get max changes*/
  for (k=0;k<temp_align->length;k++) if ((temp_align->s[0][k]!='X') && (temp_align->s[1][k]!='X') && (temp_align->s[0][k] != temp_align->s[1][k])) {
    if  ((charbit[temp_align->s[0][k]] & GAP_BIT)  && (charbit[temp_align->s[1][k]] & GAP_BIT)) ;
    else if ((charbit[temp_align->s[0][k]] & GAP_BIT)  && (!(charbit[temp_align->s[1][k]] & GAP_BIT))) max_changes[3]+=1;
    else if ((!(charbit[temp_align->s[0][k]] & GAP_BIT))  && (charbit[temp_align->s[1][k]] & GAP_BIT)) max_changes[3]+=1;
    else { 
      /*printf("(%c %c)",temp_align->s[0][k],temp_align->s[1][k]);*/
      max_changes[0]+=1; /*if not indel then must have a base change if must be change*/
      if (((charbit[temp_align->s[0][k]] & R_BIT)) && ((charbit[temp_align->s[1][k]] & Y_BIT))) max_changes[2]+=1;
      else if (((charbit[temp_align->s[0][k]] & Y_BIT)) && ((charbit[temp_align->s[1][k]] & R_BIT)))  max_changes[2]+=1;
      if (((charbit[temp_align->s[0][k]] & A_BIT)) && ((charbit[temp_align->s[1][k]] & G_BIT))) max_changes[1]+=1;
      else if (((charbit[temp_align->s[0][k]] & G_BIT)) && ((charbit[temp_align->s[1][k]] & A_BIT))) max_changes[1]+=1;
      else if (((charbit[temp_align->s[0][k]] & C_BIT)) && ((charbit[temp_align->s[1][k]] & T_BIT))) max_changes[1]+=1;
      else if (((charbit[temp_align->s[0][k]] & T_BIT)) && ((charbit[temp_align->s[1][k]] & C_BIT))) max_changes[1]+=1;
     /*transition*/
    }
  }
  printf("Max Changes %d %d %d %d\n", max_changes[0],max_changes[1],max_changes[2],max_changes[3]);
 
  temp_align=dump_align(temp_align);
  temp_align=nw(a_buf[ntax+i],a_buf[d2],values);
  /*for (j=0;j<2;j++) printf("%s\n",temp_align->s[j]);*/
   for (j=0;j<4;j++) min_changes[j]=max_changes[j]=0;
 /*for (j=0;j<2;j++) printf("%s\n",temp_align->s[j]);*/
  printf("To ");
  if (1) /*(d2 < ntax )*/ printf("%s ",a[d2]->name);
  else printf("%d ", d2);
  /*get min changes*/
  for (k=0;k<temp_align->length;k++) if (!( charbit[temp_align->s[0][k]] & charbit[temp_align->s[1][k]])) {
    if ((temp_align->s[0][k]=='-') && (!(charbit[temp_align->s[1][k]] & GAP_BIT))) min_changes[3]+=1;
    else if ((!(charbit[temp_align->s[0][k]] & GAP_BIT))  && (temp_align->s[1][k]=='-')) min_changes[3]+=1;
    else {
      min_changes[0]+=1; /*if not indel then must have a base change if must be change*/
      if ((!(charbit[temp_align->s[0][k]] & R_BIT)) && (!(charbit[temp_align->s[1][k]] & Y_BIT))) min_changes[2]+=1;
      else if ((!(charbit[temp_align->s[0][k]] & Y_BIT)) && (!(charbit[temp_align->s[1][k]] & R_BIT)))  min_changes[2]+=1;
      else  min_changes[1]+=1; /*transition*/
    }
  }
  printf("Min Changes %d %d %d %d ", min_changes[0],min_changes[1],min_changes[2],min_changes[3]);
  /*Get max changes*/
  for (k=0;k<temp_align->length;k++) if ((temp_align->s[0][k]!='X') && (temp_align->s[1][k]!='X')  && (temp_align->s[0][k] != temp_align->s[1][k])) {
    if  ((charbit[temp_align->s[0][k]] & GAP_BIT)  && (charbit[temp_align->s[1][k]] & GAP_BIT)) ;
    else if ((charbit[temp_align->s[0][k]] & GAP_BIT)  && (!(charbit[temp_align->s[1][k]] & GAP_BIT))) max_changes[3]+=1;
    else if ((!(charbit[temp_align->s[0][k]] & GAP_BIT))  && (charbit[temp_align->s[1][k]] & GAP_BIT)) max_changes[3]+=1;
    else { 
      /*printf("(%c %c)",temp_align->s[0][k],temp_align->s[1][k]);*/
      max_changes[0]+=1; /*if not indel then must have a base change if must be change*/
      if (((charbit[temp_align->s[0][k]] & R_BIT)) && ((charbit[temp_align->s[1][k]] & Y_BIT))) max_changes[2]+=1;
      else if (((charbit[temp_align->s[0][k]] & Y_BIT)) && ((charbit[temp_align->s[1][k]] & R_BIT)))  max_changes[2]+=1;
      if (((charbit[temp_align->s[0][k]] & A_BIT)) && ((charbit[temp_align->s[1][k]] & G_BIT))) max_changes[1]+=1;
      else if (((charbit[temp_align->s[0][k]] & G_BIT)) && ((charbit[temp_align->s[1][k]] & A_BIT))) max_changes[1]+=1;
      else if (((charbit[temp_align->s[0][k]] & C_BIT)) && ((charbit[temp_align->s[1][k]] & T_BIT))) max_changes[1]+=1;
      else if (((charbit[temp_align->s[0][k]] & T_BIT)) && ((charbit[temp_align->s[1][k]] & C_BIT))) max_changes[1]+=1;
     /*transition*/
    }
  }
  printf("Max Changes %d %d %d %d\n", max_changes[0],max_changes[1],max_changes[2],max_changes[3]);
 
temp_align=dump_align(temp_align);
  values->new_optimization=1;
}
free(anc);
for (j=0;j<2*ntax-1;j++) if (a_down2[j]) a_down2[j]=dump_align(a_down2[j]);
free(a_down2);

for (j=0;j<2*ntax-1;j++) if (a_buf[j]) a_buf[j]=dump_align(a_buf[j]);
free(a_buf);

return 1;
}
