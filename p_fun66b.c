/*
Copyright 1994 Ward Wheeler all rights reserved
*/

#include "align3.h"

int get_intial_length_and_opt(values,taxa, ntaxa,nbases,n_gaps)
int **taxa;
int ntaxa,nbases,n_gaps;
parameters *values;
{
int m;
int *td1,*td2,*td3,*td4,*td5;
int change_counts,transition_counts,transversion_counts,gap_counts;


/*Down pass*/
transversion_counts=transition_counts=gap_counts=change_counts=0;
td1=taxa[1];
td2=taxa[2];
td3=taxa[ntaxa+1];
if (!values->ttr) {
	for (m=0;m<nbases;m++) {
		/*fprintf(stderr,"(%d vs. %d)",td1[m],td2[m]);*/
		if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
		else {
			td3[m]=(td1[m]|td2[m]);
			++change_counts;
			}
		}
	}
else {
  for (m=0;m<nbases-values->n_transv;m++) {
    if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
    else {
      td3[m]=(td1[m]|td2[m]);
      ++transition_counts;
      }
    }
   for (m=nbases-values->n_transv;m<nbases;m++) {
      if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
	else {
	 td3[m]=(td1[m]|td2[m]);
	 ++transversion_counts;
	 }
       }
     }            
 for (m=nbases;m<(nbases+n_gaps);m++) {
  if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
  else {
   td3[m]=(td1[m]|td2[m]);
   ++gap_counts;
   }
 }
/*  fprintf(stderr,"%d change counts\n",change_counts);  */
td1=taxa[0];
td2=taxa[ntaxa+1];
td3=taxa[ntaxa];
if (!values->ttr) {
	for (m=0;m<nbases;m++) {
	if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
    else {
	td3[m]=(td1[m]|td2[m]);
      ++change_counts;
      }
    }
  }
else {
  for (m=0;m<nbases-values->n_transv;m++) {
    if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
    else {
      td3[m]=(td1[m]|td2[m]);
      ++transition_counts;
      }
    }
   for (m=nbases-values->n_transv;m<nbases;m++) {
      if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
	else {
	 td3[m]=(td1[m]|td2[m]);
	 ++transversion_counts;
	 }
       }
     }            
   for (m=nbases;m<(nbases+n_gaps);m++) {
     if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
       else {
	td3[m]=(td1[m]|td2[m]);
	++gap_counts;
	}
      }
    
/*up pass*/
td4=taxa[1];
td5=taxa[2];
for (m=0;m<nbases+n_gaps;m++) td3[m]=td1[m]; /*set base to outstates*/

m=get_best_states_init(td2,td3,taxa[1],taxa[2],nbases+n_gaps);
return ((change_counts*values->change_cost)+(gap_counts*values->gap_cost)+(transition_counts*values->transition)+(transversion_counts*values->transversion));
}

int get_interval(td1,td2,td3,values,nbases,n_gaps,value)/*add cut out for gap*/
int *td1,*td2,*td3;
int nbases,n_gaps,value;
parameters *values;
{
int m,ch,gap,ti,tv,ret_val,q;

ch=values->change_cost;
gap=values->gap_cost;
ti=values->transition;
tv=values->transversion;

ret_val=0;

/*Down pass if remove > can keep ties*/
/*fprintf(stderr,"{");*/
if (!values->ttr) {
	for (m=0;m<nbases;m++) {
		q=td2[m]|td3[m];
		if (q>little_missing_minus_one) q=td2[m]&td3[m];
		if (!(td1[m]&q)) { ret_val+=ch; /*fprintf(stderr,"(%d %d %d)",td1[m],td2[m],td3[m]);*/}
		/*else fprintf(stderr,"(%d %d %d)",td1[m],td2[m],td3[m]);*/
	if (ret_val>=value) return HUGE_COST;
	}
	}
else {
	for (m=0;m<nbases-values->n_transv;m++) {
		q=td2[m]|td3[m];
		if (q>little_missing_minus_one) q=td2[m]&td3[m];
		if (!(td1[m]&q)) ret_val+=ti;
		if (ret_val>=value) return HUGE_COST;
		}
	for (m=nbases-values->n_transv;m<nbases;m++) {
		q=td2[m]|td3[m];
		if (q>little_missing_minus_one) q=td2[m]&td3[m];
		if (!(td1[m]&q)) ret_val+=tv;
	if (ret_val>=value) return HUGE_COST;
		}            
	}
 for (m=nbases;m<(nbases+n_gaps);m++) {
	q=td2[m]|td3[m];
	if (q>little_missing_minus_one) q=td2[m]&td3[m];
	if (!(td1[m]&q)) ret_val+=gap;
	if (ret_val>=value) return HUGE_COST;
	}
/*fprintf(stderr,"}");*/
return ret_val;
}

void get_reopt(td1,td2,td3,values,nbases,n_gaps,value,reopt_array)/*add cut out for gap*/
int *td1,*td2,*td3,*reopt_array;
int nbases,n_gaps,value;
parameters *values;
{
int m,ch,gap,ti,tv,ret_val,q;

/*Down pass*/
for (m=0;m<nbases+n_gaps;m++) {
	q=td2[m]|td3[m];
	if (q>little_missing_minus_one) q=td2[m]&td3[m];
	if (!(td1[m]&q)) reopt_array[m]=1;
	else reopt_array[m]=0;
	}
}


int get_best_states_init(a,b,c,d,npos)
int *a,*b,*c,*d,npos;
{
int c16,c8,c4,c2,c1;
int best, ret_val,i;

for (i=0;i<npos;i++) {
	if (b[i]&c[i]&d[i]) a[i]=(b[i]&c[i]&d[i]);
	else {
		c16=c8=c4=c2=c1=0;
		a[i]=0;
		if (b[i]&1) ++c1;
		if (c[i]&1) ++c1;
		if (d[i]&1) ++c1;
		if (b[i]&2) ++c2;
		if (c[i]&2) ++c2;
		if (d[i]&2) ++c2;
		if (b[i]&4) ++c4;
		if (c[i]&4) ++c4;
		if (d[i]&4) ++c4;
		if (b[i]&8) ++c8;
		if (c[i]&8) ++c8;
		if (d[i]&8) ++c8;
		if (b[i]&16) ++c16;
		if (c[i]&16) ++c16;
		if (d[i]&16) ++c16;
		best=max(c1,max(c2,max(c4,max(c8,c16))));
		if (c1==best) a[i]+=1;
		if (c2==best) a[i]+=2;
		if (c4==best) a[i]+=4;
		if (c8==best) a[i]+=8;
		if (c16==best) a[i]+=16;
		}
	}

return 1;
}

int get_best_states_init2(a,b,c,d,npos,reopt)
int *a,*b,*c,*d,npos,*reopt;
{
int c16,c8,c4,c2,c1;
int best, ret_val,i;

for (i=0;i<npos;i++) if (!reopt[i]) {
	if (b[i]&c[i]&d[i]) a[i]=(b[i]&c[i]&d[i]);
	else {
		c16=c8=c4=c2=c1=0;
		a[i]=0;
		if (b[i]&1) ++c1;
		if (c[i]&1) ++c1;
		if (d[i]&1) ++c1;
		if (b[i]&2) ++c2;
		if (c[i]&2) ++c2;
		if (d[i]&2) ++c2;
		if (b[i]&4) ++c4;
		if (c[i]&4) ++c4;
		if (d[i]&4) ++c4;
		if (b[i]&8) ++c8;
		if (c[i]&8) ++c8;
		if (d[i]&8) ++c8;
		if (b[i]&16) ++c16;
		if (c[i]&16) ++c16;
		if (d[i]&16) ++c16;
		best=max(c1,max(c2,max(c4,max(c8,c16))));
		if (c1==best) a[i]+=1;
		if (c2==best) a[i]+=2;
		if (c4==best) a[i]+=4;
		if (c8==best) a[i]+=8;
		if (c16==best) a[i]+=16;
		}
	}

return 1;
}

int get_best_states_opt(a,b,c,d,e,npos,reopt)
int *a,*b,*c,*d,*e,npos,*reopt;
{
int c16,c8,c4,c2,c1;
int best, ret_val,i;

for (i=0;i<npos;i++) if (reopt[i]) {
	if (b[i]&c[i]) a[i]=(b[i]&c[i]);
	else {
		c16=c8=c4=c2=c1=0;
		a[i]=0;
		if (b[i]&1) ++c1;
		if (d[i]&1) ++c1;
		if (e[i]&1) ++c1;
		if (b[i]&2) ++c2;
		if (d[i]&2) ++c2;
		if (e[i]&2) ++c2;
		if (b[i]&4) ++c4;
		if (d[i]&4) ++c4;
		if (e[i]&4) ++c4;
		if (b[i]&8) ++c8;
		if (d[i]&8) ++c8;
		if (e[i]&8) ++c8;
		if (b[i]&16) ++c16;
		if (d[i]&16) ++c16;
		if (e[i]&16) ++c16;
		best=max(c1,max(c2,max(c4,max(c8,c16))));
		if (c1==best) a[i]+=1;
		if (c2==best) a[i]+=2;
		if (c4==best) a[i]+=4;
		if (c8==best) a[i]+=8;
		if (c16==best) a[i]+=16;
		}
	}

return 1;
}

int get_best_states(a,b,c)
int a,b,c;
{
int c16,c8,c4,c2,c1;
int best, ret_val;

c16=c8=c4=c2=c1=0;
ret_val=0;

if (a&1) ++c1;
if (b&1) ++c1;
if (c&1) ++c1;
if (a&2) ++c2;
if (b&2) ++c2;
if (c&2) ++c2;
if (a&4) ++c4;
if (b&4) ++c4;
if (c&4) ++c4;
if (a&8) ++c8;
if (b&8) ++c8;
if (c&8) ++c8;
if (a&16) ++c16;
if (b&16) ++c16;
if (c&16) ++c16;
best=max(c1,max(c2,max(c4,max(c8,c16))));
if (c1==best) ret_val+=1;
if (c2==best) ret_val+=2;
if (c4==best) ret_val+=4;
if (c8==best) ret_val+=8;
if (c16==best) ret_val+=16;
/*fprintf(stderr,"(%d %d %d=>%d)",a,b,c,ret_val);*/
return ret_val;
}

int get_new_quick(ntax,nbases,sequence,n_gaps,values)
int ntax,nbases,**sequence,n_gaps;
parameters *values;
{
	int i,j,value,*temp_rep, cur_val,n,m;
	int initial_length,d1,d2;
	int *anc,check,*reopt_array;
	int *td1,*td2,*td3,*td4,q;
	int *tree_rep, **nodes;


	nodes=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int)nodes);
	for (i=0;i<ntax-1;i++) {
		nodes[i]=(int *)malloc(2*sizeof(int));
		assert((int) nodes[i]);
		}
	tree_rep=(int *)malloc((ntax-3)*sizeof(int));
	assert((int) tree_rep);
	anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
	assert((int) anc);
	temp_rep=(int *)malloc((ntax-3)*sizeof(int));
	assert((int) temp_rep);
	reopt_array=(int *)malloc((nbases+n_gaps)*sizeof(int));
	assert((int) reopt_array);
	
	/*get initial tree*/
	nodes[0][0]=0;
	nodes[0][1]=ntax+1;
	nodes[1][0]=1;
	nodes[1][1]=2;
	anc[0]=ntax;
	anc[1]=ntax+1;
	anc[2]=ntax+1;
	anc[ntax+1]=ntax;
	/* and optimized nodes*/
	initial_length=get_intial_length_and_opt(values,sequence, ntax,nbases,n_gaps);
	for (i=0;i<ntax-3;i++) {/*add taxa in turn*/
		value=HUGE_COST;
		for (j=1;j<=((2*i)+3);j++) {/*places to go*/
			d1=j/2;
			d2=j%2;
			temp_rep[i]=j;
			/*if (nodes[d1][d2]<ntax) cur_val=get_interval(sequence[i+3],sequence[nodes[d1][d2]],sequence[d1+ntax],values,nbases,n_gaps,value);
			else */
			cur_val=get_interval2(sequence[i+3],sequence[nodes[d1][d2]],sequence[d1+ntax],values,nbases,n_gaps,value);
			if (cur_val<value) {
				tree_rep[i]=j;
				value=cur_val;
				}
			if (value==0) break;
			}
		temp_rep[i]=tree_rep[i];
		/*make new nodes*/
		nodes[i+2][0]=i+3;
		d1=tree_rep[i]/2;
		d2=tree_rep[i]%2;
		nodes[i+2][1]=nodes[d1][d2];
		nodes[d1][d2]=ntax+i+2;
		/*get new anc and its anc*/
		anc[nodes[i+2][0]]=i+2+ntax;
		anc[nodes[i+2][1]]=i+2+ntax;
		anc[i+2+ntax]=d1+ntax;
		td1=sequence[i+3];
		td2=sequence[anc[i+2+ntax]];
		td3=sequence[i+2+ntax]; 
		td4=sequence[nodes[i+2][1]];
		for (j=0;j<nbases+n_gaps;j++) {
			q=td4[j]|td2[j];
			/*q=td4[j]&td2[j];*/
			if (nodes[i+2][1]<ntax) if (q>little_missing_minus_one) q=td4[j]&td2[j];
			td3[j]=td1[j]&q;
			if (td3[j]) reopt_array[j]=0;
			else reopt_array[j]=1;
			}
		if ((d1+ntax)>ntax) check=get_best_states_init2(sequence[d1+ntax],sequence[nodes[d1][0]],sequence[nodes[d1][1]],sequence[anc[d1+ntax]],nbases+n_gaps,reopt_array);
		/*for (n=0;n<i+3;n++) anc[nodes[n][0]]=anc[nodes[n][1]]=ntax+n;*/
		/*reoptimize internal*/
		check=do_optimize_nodes(values,sequence,nodes,ntax,i+4,nbases,n_gaps,anc,reopt_array);
		initial_length+=value;  
	}
free(anc);
free(temp_rep);
free(reopt_array);
make_nodes(tree_rep,nodes,ntax,nbases,sequence,ntax,n_gaps,values,&initial_length);
if (values->sbr) {
	initial_length=do_swap_new(nodes,ntax,nbases,n_gaps,initial_length,sequence,values,0);
	if (values->phylo_score==6) initial_length=do_swap_new(nodes,ntax,nbases,n_gaps,initial_length,sequence,values,1);
	}
for (i=0;i<ntax-1;i++) free(nodes[i]);
free(nodes);
free(tree_rep);
return initial_length;
}

int do_optimize_nodes(values,taxa,nodes,ntaxa,nent,nbases,n_gaps,anc,reopt)
int **nodes, **taxa, ntaxa, nent,*anc,*reopt;
parameters *values;
{
int *td1,*td2,*td3,*td4,*td5;
int n,m,d1,d2;
int length=0;
int ch,gap,ti,tv;
int ***up_nodes,i_am,i;
int found_one;

up_nodes=(int ***)malloc(((2*ntaxa)-1)*sizeof(int**));
assert((int) up_nodes);
for (i=0;i<(2*ntaxa)-1;i++) {
	up_nodes[i]=(int **)malloc((3)*sizeof(int *));
	assert((int) up_nodes[i]);
	up_nodes[i][0]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*holds down pass*/
	assert((int) up_nodes[i][0]);
	up_nodes[i][1]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*holds up pass to d1*/
	assert((int) up_nodes[i][1]);
	up_nodes[i][2]=(int *)malloc((nbases+n_gaps)*sizeof(int)); /*hold uppass to d2*/
	assert((int) up_nodes[i][2]);
	}
	
ch=values->change_cost;
gap=values->gap_cost;
ti=values->transition;
tv=values->transversion;


/*down pass*/
for (n=0;n<ntaxa;n++) {
	for (i=0;i<nbases+n_gaps;i++) if (reopt[i]) up_nodes[n][0][i]=taxa[n][i];
	values->tree_made[n]=1;
	}
for (n=ntaxa;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;
found_one=1;
while (found_one) {
	found_one=0;
for (n=(nent-2);n>=0;--n){
	d1=nodes[n][0];
	d2=nodes[n][1];
	if ((values->tree_made[d1]) && (values->tree_made[d2]) && (!values->tree_made[ntaxa+n])) {
		found_one=1;
		td1=up_nodes[d1][0];
		td2=up_nodes[d2][0];
		td3=up_nodes[ntaxa+n][0];
		/*fprintf(stderr,"done %d ",ntaxa+n);*/
		values->tree_made[ntaxa+n]=1;
		if (1) {
			for (m=0;m<nbases+n_gaps;m++) if (reopt[m]) {
				td3[m]=(td1[m]&td2[m]);
				if (!td3[m]) td3[m]=(td1[m]|td2[m]);
				}
			}
		else {
			/*only need this to check*/
			if (!values->ttr) {
				for (m=0;m<nbases;m++) {
					if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
					else {
						td3[m]=(td1[m]|td2[m]);
						length+=ch;
						}
					}
				}
			else {
			  for (m=0;m<nbases-values->n_transv;m++) {
			    if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
			    else {
			      td3[m]=(td1[m]|td2[m]);
			      length+=ti;
			     }
			    }
			   for (m=nbases-values->n_transv;m<nbases;m++) {
			      if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
				else {
				 td3[m]=(td1[m]|td2[m]);
				 length+=tv;
				 }
			       }
			     }            
			 for (m=nbases;m<(nbases+n_gaps);m++) {
			  if (td1[m]&td2[m]) td3[m]=(td1[m]&td2[m]);
			  else {
			   td3[m]=(td1[m]|td2[m]);
			   length+=gap;
			   }
			 }      
			}                                               
		}
	}
}

/*up pass*/     
/*remember to remove extraneous up pass (ie when leads to a terminal*/
td1=up_nodes[ntaxa][0];
td2=up_nodes[ntaxa][1];
td3=up_nodes[ntaxa][2];
for (m=0;m<nbases+n_gaps;m++) if (reopt[m]) td1[m]=td2[m]=td3[m]=taxa[ntaxa][m]=taxa[0][m]; /*set basal node = to outgroup*/

for (n=ntaxa+1;n<(ntaxa+nent-1);++n) values->tree_made[n]=0;
values->tree_made[ntaxa]=1;
found_one=1;
while (found_one) {
	found_one=0;
for (n=1;n<nent-1;n++){
	if ((values->tree_made[anc[ntaxa+n]]) && (!values->tree_made[ntaxa+n])) {
		found_one=1;
		values->tree_made[ntaxa+n]=1;
		if (nodes[anc[ntaxa+n]-ntaxa][0]==(n+ntaxa)) i_am=1;
		else i_am=2;
		td1=up_nodes[n+ntaxa][1];
		td2=up_nodes[n+ntaxa][2];
		td3=up_nodes[anc[ntaxa+n]][i_am];
		td4=up_nodes[nodes[n][1]][0];
		td5=up_nodes[nodes[n][0]][0];
		i=get_best_states_opt(taxa[ntaxa+n],td3,up_nodes[ntaxa+n][0],td5,td4,nbases+n_gaps,reopt);
		/*get up pass states for this node*/
		if (nodes[n][0]>(ntaxa-1)) {
			for (m=0;m<nbases+n_gaps;m++) if (reopt[m]) {
				td1[m]=(td4[m]&td3[m]);
				if (!td1[m]) td1[m]=(td4[m]|td3[m]);
				}
			}
		if (nodes[n][1]>(ntaxa-1)) {
			for (m=0;m<nbases+n_gaps;m++) if (reopt[m]) {
				td2[m]=(td5[m]&td3[m]);
				if (!td2[m]) td2[m]=(td5[m]|td3[m]);
				}
			}
			
		}
	}
}
for (i=0;i<(2*ntaxa)-1;i++) {
	free(up_nodes[i][0]);
	free(up_nodes[i][1]);
	free(up_nodes[i][2]);
	free(up_nodes[i]);
	}
free(up_nodes);
return length;
}

void modify_nodes(t,n,ntaxa,nent)
int *t;
int **n;
int ntaxa,nent;
{
int i,a,b;

n[0][0]=0;
n[0][1]=ntaxa+1;
n[1][0]=1;
n[1][1]=2;
for (i=3;i<nent;i++) {
	a=t[i-3]/2;
	b=t[i-3]%2;
	n[i-1][0]=i;
	n[i-1][1]=n[a][b];
	n[a][b]=ntaxa+i-1;
	}
/*for (i=0;i<nent-2;i++) fprintf(stderr,"n %d 0 %d n %d 1 %d \n",i,n[i][0],i,n[i][1]);*/
	
}
