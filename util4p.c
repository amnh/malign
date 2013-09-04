/*Copyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Goloboff*/

#include "align3.h"

void get_up_pass_new_opt_parallel(a_final,a_down,nodes,ntax,values,nent,anc,a_down2)
alignment **a_down, **a_final,**a_down2;
int **nodes, ntax,nent, *anc;
parameters *values;
{
int i,d1,d2,j,found_one,n,ii;
int node_is_from;
char *tempstr,*t0,*t1,*s0,*s1;
int nent_minus_one_new;
int *done,node_to_be_sent,node_received;
int bufid,type,bytes,source,index,info;
int sent,received;


done=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) done);
for (i=0;i<ntax;i++) done[i]=1;
for (i=ntax;i<ntax+nent-1;i++) done[i]=0;


if (a_final[nodes[0][1]]) dump_align(a_final[nodes[0][1]]);
a_final[nodes[0][1]]=make_align(a_down[nodes[0][1]]);

done[nodes[0][1]]=1;

found_one=1;
while (found_one) {
	found_one=0;
	sent=received=0;
	for (n=1;n<nent-1;n++) if ((!done[n+ntax]) && (done[anc[n+ntax]])) {
		found_one=1;
		node_to_be_sent=ntax+n;
		if (sent < (values->num_hosts-1)) {
			pvm_initsend( PvmDataDefault );
			pvm_pkint(&node_to_be_sent,1,1);
			pack_align(a_final[anc[ntax+n]]); /*a_up[anc[ntax+n]][node_is_from]*/
			pack_align(a_down2[ntax+n]);
			pvm_send(values->tids[++sent],PVM_PARALLEL_NW_NEW_OPT_UP_DO);
			}
		else {
			bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_UP_DONE);
				if (bufid) {
					info=pvm_bufinfo(bufid, &bytes, &type, &source);
					pvm_upkint(&node_received,1,1);
					done[node_received]=1;
					if (a_final[node_received]) a_final[node_received]=dump_align(a_final[node_received]);
					a_final[node_received]=unpack_align_and_score(values);
					++received;
					}
				/*send out next*/
				pvm_initsend( PvmDataDefault );
				pvm_pkint(&node_to_be_sent,1,1);
				pack_align(a_final[anc[ntax+n]]); /*a_up[anc[ntax+n]][node_is_from]*/
			    pack_align(a_down2[ntax+n]);
				pvm_send(source,PVM_PARALLEL_NW_NEW_OPT_UP_DO);
				++sent;
				}

		}/*n loop all sent some received*/
	/*receive rest*/
	for (ii=received;ii<sent;ii++) {
		bufid=pvm_recv(-1,PVM_PARALLEL_NW_NEW_OPT_UP_DONE);
		if (bufid) {
			info=pvm_bufinfo(bufid, &bytes, &type, &source);
			pvm_upkint(&node_received,1,1);
			if (a_final[node_received]) a_final[node_received]=dump_align(a_final[node_received]);
			a_final[node_received]=unpack_align_and_score(values);
			done[node_received]=1;
			}
		}
	}
free(done);
}

int up_node_remote(values)
parameters *values;
{
int node_to_be_done;
alignment *a1,*a2;
alignment *o1;
int i;

a1=NULL;
a2=NULL;
o1=NULL;


pvm_upkint(&node_to_be_done,1,1);
a1=unpack_align_and_score(values); /*final of anc*/
a2=unpack_align_and_score(values); /*d1,d2,and node*/

free(a2->name);
a2->name=(char *)malloc(2*sizeof(char));
assert((int)a2->name);
a2->name[0]='b'; a2->name[1]='\0';
free(a1->name);
a1->name=(char *)malloc(2*sizeof(char));
assert((int)a1->name);
a1->name[0]='a'; a1->name[1]='\0';

values->new_optimization=0;
o1=nw(a2,a1,values);
values->new_optimization=1;
convert_to_final(o1,values);

pvm_initsend( PvmDataDefault );
pvm_pkint(&node_to_be_done,1,1);
pack_align(o1);
pvm_send(values->tids[0],PVM_PARALLEL_NW_NEW_OPT_UP_DONE);
if (a1) a1=dump_align(a1);
if (a2) a2=dump_align(a2);
if (o1) o1=dump_align(o1);
return (0);
}
