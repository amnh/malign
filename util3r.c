/*Copyrite 1995 Ward Wheeler, all rights reserved*/
/*functions to do alignments the new faster way ala Goloboff*/

#include "align3.h"

extern int charbit[256];
extern char bitchar[32];

char char_overlap_best(a,b,c,d)
char a,b,c,d;
{
int in_c, in_b, in_a, in_d, out;
char out_char;
int c16,c8,c4,c2,c1,best;
int c32,c64,c128;

in_a=charbit[(int)a];
in_b=charbit[(int)b];
in_c=charbit[(int)c];
in_d=charbit[(int)d];



out= (in_b & in_c & in_d);
if (!out) {
	c16=c8=c4=c2=c1=0;
	out=0;
	if (in_b&1) ++c1;
	if (in_c&1) ++c1;
	if (in_d&1) ++c1;
	if (in_b&2) ++c2;
	if (in_c&2) ++c2;
	if (in_d&2) ++c2;
	if (in_b&4) ++c4;
	if (in_c&4) ++c4;
	if (in_d&4) ++c4;
	if (in_b&8) ++c8;
	if (in_c&8) ++c8;
	if (in_d&8) ++c8;
	if (in_b&16) ++c16;
	if (in_c&16) ++c16;
	if (in_d&16) ++c16;
	best=max(c1,max(c2,max(c4,max(c8,c16))));
	if (c1==best) out+=1;
	if (c2==best) out+=2;
	if (c4==best) out+=4;
	if (c8==best) out+=8;
	if (c16==best) out+=16;
	}

return ((char) bitchar[out]);
}

char overlap(a,b)
char a,b;
{
int in_a, in_b, out;
char out_char;


in_a=charbit[(int)a];
in_b=charbit[(int)b];
out=(in_a & in_b);
if (!out) return 0;
else  return ((char) bitchar[out]);}

char overlap_best(a,b,c,d)
char a,b,c,d;
{
int in_c, in_b, in_a, in_d, out;
char out_char;
int c16,c8,c4,c2,c1,best;


in_a=charbit[(int)a];
in_b=charbit[(int)b];
in_c=charbit[(int)c];
in_d=charbit[(int)d];


out= (in_b & in_c & in_d);
if (!out) {
	c16=c8=c4=c2=c1=0;
	out=0;
	if (in_b&1) ++c1;
	if (in_c&1) ++c1;
	if (in_d&1) ++c1;
	if (in_b&2) ++c2;
	if (in_c&2) ++c2;
	if (in_d&2) ++c2;
	if (in_b&4) ++c4;
	if (in_c&4) ++c4;
	if (in_d&4) ++c4;
	if (in_b&8) ++c8;
	if (in_c&8) ++c8;
	if (in_d&8) ++c8;
	if (in_b&16) ++c16;
	if (in_c&16) ++c16;
	if (in_d&16) ++c16;
	best=max(c1,max(c2,max(c4,max(c8,c16))));
	if (c1==best) out+=1;
	if (c2==best) out+=2;
	if (c4==best) out+=4;
	if (c8==best) out+=8;
	if (c16==best) out+=16;
	}
return ((char) bitchar[out]);
}

int new_overlap_best(a,b,c,d,values)
char a,b,c,d;
parameters *values;
{
int in_c, in_b, in_a, in_d, out,best;
int out_char;
int c16,c8,c4,c2,c1;

in_a=charbit[(int)a];
in_b=charbit[(int)b];
in_c=charbit[(int)c];
in_d=charbit[(int)d];

if (!values->delta && !values->ttr && (values->gap_cost==values->change_cost)) {
    if (((in_c & in_d) == in_d) && (in_c != in_d)) out= in_d;
    else if (!(in_a & in_b)) out= ( in_c | in_d);
    else  out= (in_c | ((in_a | in_b) & in_d));
    }
else {
    out= (in_a & in_b & in_d);
    if (!out) {
        c1=c2=c4=c8=c16=0;
        c1=(values->lookup[A_BIT][in_a].cost+values->lookup[A_BIT][in_b].cost+values->lookup[A_BIT][in_d].cost);
        c2=(values->lookup[C_BIT][in_a].cost+values->lookup[C_BIT][in_b].cost+values->lookup[C_BIT][in_d].cost);
        c4=(values->lookup[G_BIT][in_a].cost+values->lookup[G_BIT][in_b].cost+values->lookup[G_BIT][in_d].cost);
        c8=(values->lookup[T_BIT][in_a].cost+values->lookup[T_BIT][in_b].cost+values->lookup[T_BIT][in_d].cost);
       c16=(values->lookup[GAP_BIT][in_a].cost+values->lookup[GAP_BIT][in_b].cost+values->lookup[GAP_BIT][in_d].cost);
       best=min(c1,min(c2,min(c4,min(c8,c16))));
       if (c1==best) out+=A_BIT;
       if (c2==best) out+=C_BIT;
       if (c4==best) out+=G_BIT;
       if (c8==best) out+=T_BIT;
       if (c16==best) out+=GAP_BIT;
        }
    }


return ((char) bitchar[out]);
}

/*this is just like the old one but leaves in the optimized gaps so correspondances are ok
but removes the preexisting gaps again for correspondances*/
apolist_holder *make_anc_up(a,values,bl,ancestor,descendent,task,apomorphy,node_number,return_holder)
alignment *a;
parameters *values;
int *bl,ancestor,descendent,task;
apo_thang **apomorphy;
int node_number;
apolist_holder *return_holder;
{
char **nbsa;
int i,j, new_length, old_length, score_test,k0,k1,k2;
alignment *b;
int gaps=0, changes=0, ends=0;

b=NULL;
if  (a->n_seqs>1) {
	if (a->n_seqs!=2) {
		fprintf(stderr,"This is wrong, %d sequences when there should be 2\n",a->n_seqs);
		exit(-1);
		}
	nbsa=(char **)malloc(3*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<3;i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
	/*if not bitted already get_bases*/
	for (i=0;i<2;i++) for (j=0;j<a->length;j++) nbsa[i][j]=charbit[(int)a->s[i][j]];

	/*assume two now and get union/intersection*/
	ends=changes=gaps=new_length=score_test=0;
	bl[0]=bl[1]=bl[2]=bl[3]=0;
	if (task==1) { /*up pass with choose*/
		for (j=0;j<a->length;j++)  {
			if (!is_ambig(nbsa[descendent][j])) nbsa[2][j]=nbsa[descendent][j];
			else nbsa[2][j]=(nbsa[ancestor][j]&nbsa[descendent][j]);
			if (!nbsa[2][j]) nbsa[2][j]=choose_one(nbsa[descendent][j],values->rand_apo,values);
			if (nbsa[ancestor][j]!=nbsa[2][j]) {
				if ((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[2][j]==GAP_BIT))) {
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) ++gaps;
					else ++ends;
					}
				else ++changes;
				}
			if (nbsa[descendent][j]!=GAP_BIT) ++new_length;
			}
		if (values->apolist) {
			apomorphy[0]=(apo_thang *)malloc(sizeof(apo_thang));
			assert((int) apomorphy[0]);
			apomorphy[1]=(apo_thang *)malloc(sizeof(apo_thang));
			assert((int) apomorphy[1]);
			apomorphy[2]=(apo_thang *)malloc(sizeof(apo_thang));
			assert((int) apomorphy[2]);
			k0=k1=k2=0;
			for (j=0;j<a->length;j++)  {
				if (nbsa[ancestor][j]!=nbsa[2][j]) {
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) {
						if ((nbsa[ancestor][j]&GAP_BIT) || (nbsa[2][j]&GAP_BIT)) {
							apomorphy[1]=(apo_thang *)realloc(apomorphy[1],(k1+1)*sizeof(apo_thang));
							assert((int) apomorphy[1]);
							apomorphy[1][k1].anc_pos=real_position(a->s[ancestor],j);
							apomorphy[1][k1].desc_pos=real_position(a->s[descendent],j);
							apomorphy[1][k1].anc_char=a->s[ancestor][j];
							apomorphy[1][k1].desc_char=return_character(nbsa[2][j]);
							++k1;
							}
						else {
							apomorphy[2]=(apo_thang *)realloc(apomorphy[2],(k2+1)*sizeof(apo_thang));
							assert((int) apomorphy[2]);
							apomorphy[2][k2].anc_pos=real_position(a->s[ancestor],j);
							apomorphy[2][k2].desc_pos=real_position(a->s[descendent],j);							apomorphy[2][k2].anc_char=a->s[ancestor][j];
							apomorphy[2][k2].anc_char=a->s[ancestor][j];
							apomorphy[2][k2].desc_char=return_character(nbsa[2][j]);
							++k2;
							}
						}
					else {
						if ((nbsa[ancestor][j]&GAP_BIT) || (nbsa[2][j]&GAP_BIT)) {
							apomorphy[0]=(apo_thang *)realloc(apomorphy[0],(k0+1)*sizeof(apo_thang));
							assert((int) apomorphy[0]);
							apomorphy[0][k0].anc_pos=real_position(a->s[ancestor],j);
							apomorphy[0][k0].desc_pos=real_position(a->s[descendent],j); 							apomorphy[0][k0].anc_char=a->s[ancestor][j];
							apomorphy[0][k0].anc_char=a->s[ancestor][j];
							apomorphy[0][k0].desc_char=return_character(nbsa[2][j]);
							++k0;
							}
						else {
							apomorphy[2]=(apo_thang *)realloc(apomorphy[2],(k2+1)*sizeof(apo_thang));
							assert((int) apomorphy[2]);
							apomorphy[2][k2].anc_pos=real_position(a->s[ancestor],j);
							apomorphy[2][k2].desc_pos=real_position(a->s[descendent],j);
							apomorphy[2][k2].anc_char=a->s[ancestor][j];
							apomorphy[2][k2].desc_char=return_character(nbsa[2][j]);
							++k2;
							}
						}
					}
				}
			if (k0!=ends) fprintf(stderr,"Error in apomorphy list 0 (%d %d)\n",k0,ends);
			/*for (i=0;i<k0;i++) fprintf(stderr,"[%d][0][%d]=%d %c %d %c ",node_number,i,apomorphy[0][i].anc_pos,apomorphy[0][i].anc_char,apomorphy[0][i].desc_pos,apomorphy[0][i].desc_char);*/
			if (k1!=gaps) fprintf(stderr,"Error in apomorphy list 1 (%d %d)\n",k1,gaps);
			/*for (i=0;i<k1;i++) fprintf(stderr,"[%d][1][%d]=%d  %c %d %c ",node_number,i,apomorphy[1][i].anc_pos,apomorphy[1][i].anc_char,apomorphy[1][i].desc_pos,apomorphy[1][i].desc_char);*/
			if (k2!=changes) fprintf(stderr,"Error in apomorphy list 2 (%d %d)\n",k2,changes);
			/*for (i=0;i<k2;i++) fprintf(stderr,"[%d][2][%d]=%d  %c %d %c ",node_number,i,apomorphy[2][i].anc_pos,apomorphy[2][i].anc_char,apomorphy[2][i].desc_pos,apomorphy[2][i].desc_char);
			fprintf(stderr,"end\n");*/
			}/*end of apolist*/
		}
	else if (task==0) { /*min branch lengths*/
		for (j=0;j<a->length;j++)  {
			if (!is_ambig(nbsa[descendent][j])) nbsa[2][j]=nbsa[descendent][j];
			else {
				nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				if (!nbsa[2][j]) nbsa[2][j]=nbsa[descendent][j];
				}
			if (!(nbsa[ancestor][j]&nbsa[2][j])) {/*min gaps*/
				if ((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[2][j]==GAP_BIT))) {
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) ++gaps;
					else ++ends;
					}
				}
			if (!(nbsa[ancestor][j]&nbsa[2][j])) {/*min chnages*/
				if ((nbsa[ancestor][j]&GAP_BIT) || ((nbsa[2][j]&GAP_BIT))) {
					}
				else ++changes;
				}
			if (!(nbsa[ancestor][j]&nbsa[2][j])) {/*min gaps*/
				if (!((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[2][j]==GAP_BIT))))  bl[3]+=(values->change_cost);
				else{
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) bl[3]+=(values->gap_cost);
					else bl[3]+=(values->gap_cost-values->leading_gap_cost);
					}
				}
			if (nbsa[descendent][j]!=GAP_BIT) ++new_length;
			}
		}
	else if (task==3) { /* max branch length*/
		for (j=0;j<a->length;j++)  {
			if (!is_ambig(nbsa[descendent][j])) nbsa[2][j]=nbsa[descendent][j];
			else {
				nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				if (!nbsa[2][j]) nbsa[2][j]=nbsa[descendent][j];
				}
			if (nbsa[ancestor][j]!=nbsa[2][j]) {
				if ((nbsa[ancestor][j]&GAP_BIT) || ((nbsa[2][j]&GAP_BIT))) {/*max gaps*/
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) ++gaps;
					else ++ends;
					}
				}
			if (nbsa[ancestor][j]!=nbsa[2][j]) {/*max changes*/
				if ((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[2][j]==GAP_BIT))) {
					}
				else ++changes;
				}
			if (nbsa[ancestor][j]!=nbsa[2][j]) {
				if ((nbsa[ancestor][j]&GAP_BIT) || ((nbsa[2][j]&GAP_BIT))) {
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) bl[3]+=(values->gap_cost);
					else bl[3]+=(values->gap_cost-values->leading_gap_cost);
					}
				else bl[3]+=(values->change_cost);;
				}
			if (nbsa[descendent][j]!=GAP_BIT) ++new_length;
			}
		}

	a->score=(gaps*values->gap_cost);
	a->score+=(changes*values->change_cost);
	a->score-=(ends*values->leading_gap_cost);
	/*set branch_lengths*/
	bl[0]=ends;
	bl[1]=gaps;
	bl[2]=changes;
	if (task==1) bl[3]=((bl[0]*(values->gap_cost-values->leading_gap_cost))+(bl[1]*values->gap_cost)+(bl[2]*values->change_cost));
	/*copy back*/
	old_length=a->length;
	if (values->gap_in_hypanc) new_length=a->length;
	b=(alignment *)malloc(sizeof(alignment));
	assert ((int) b);
	b->n_seqs=1;
	b->type_weight=a->type_weight;
	b->length=new_length;
	b->score=a->score; /*score_test; b->score;*/
	b->name=(char *)malloc((strlen(a->name)+1)*sizeof(char));
	assert((int) b->name);
	b->name=(char *)strcpy(b->name,a->name);
	b->s=(char **)malloc((1+b->type_weight)*sizeof(char *));
	assert((int) b->s);
	for (i=0;i<(b->n_seqs+b->type_weight);i++) {
		b->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) b->s[i]);
		b->s[i][new_length]='\0';
		}
	if (b->type_weight) b->s[b->n_seqs]=strcpy(b->s[b->n_seqs],a->s[a->n_seqs]);
	b->taxon_name=(char **)malloc(sizeof(char *));
	assert((int) b->taxon_name);
	b->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) b->taxon_name[0]);
	b->taxon_name[0][0]='H';
	b->taxon_name[0][1]='y';
	b->taxon_name[0][2]='p';
	b->taxon_name[0][3]='A';
	b->taxon_name[0][4]='n';
	b->taxon_name[0][5]='c';
	b->taxon_name[0][6]='\0';
	j=0;
	if (!values->gap_in_hypanc) {
		for (i=0;i<old_length;i++) if (nbsa[descendent][i]!=GAP_BIT) { b->s[0][j]=(char)bitchar[nbsa[2][i]];
			j++;
			}
		}
	else {
		for (i=0;i<old_length;i++) { b->s[0][j]=(char)bitchar[nbsa[2][i]];
			j++;
			}
	}
	/*free*/
	for (i=0;i<3;i++) free(nbsa[i]);
	free(nbsa);
	}
return_holder->z=b;
if ((task==1) && (values->apolist)) {
	return_holder->at0=apomorphy[0];
	return_holder->at1=apomorphy[1];
	return_holder->at2=apomorphy[2];
	}
else {
	return_holder->at0=NULL;
	return_holder->at1=NULL;
	return_holder->at2=NULL;
	}
return return_holder;
}

/*like choose one this will only yield gaps if a choice assuming gaps > chnages*/
char choose_one2(inchar,rand)
char inchar,rand;
{
int choice;

/*if unambig then no choice*/
switch (inchar) {
	case 'A': return 'A';	break;
	case 'C': return 'C';	break;
	case 'G': return 'G';	break;
	case 'T': return 'T';	break;
	case '-': return '-';	break;
	case 'U': return 'T';	break;
	}
if (rand){
	switch (inchar) {
		case 'N':
			choice=my_randomize(4);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else if (choice==2)  return 'G';
			else return 'T';
			break;
		case 'X':
			choice=my_randomize(4);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else if (choice==2)  return 'G';
			else return 'T';
			break;
		case 'R':
			choice=my_randomize(2);
			if (choice==0)  return 'A';
			else return 'G';
			break;
		case 'Y':
			choice=my_randomize(2);
			if (choice==0)  return 'C';
			else return 'T';
			break;
		case 'M':
			choice=my_randomize(2);
			if (choice==0)  return 'A';
			else return 'C';
			break;
		case 'W':
			choice=my_randomize(2);
			if (choice==0)  return 'A';
			else return 'T';
			break;
		case 'S':
			choice=my_randomize(2);
			if (choice==0)  return 'G';
			else return 'C';
			break;
		case 'K':
			choice=my_randomize(2);
			if (choice==0)  return 'G';
			else return 'T';
			break;
		case 'B':
			choice=my_randomize(3);
			if (choice==0)  return 'C';
			else if (choice==1)  return 'G';
			else return 'T';
			break;
		case 'D':
			choice=my_randomize(3);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'G';
			else return 'T';
			break;
		case 'H':
			choice=my_randomize(3);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else return 'T';
			break;
		case 'V':
			choice=my_randomize(3);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else return 'G';
			break;

		case 'E':
			/*return 'A'; comment out to preserve gap possibility*/
			choice=my_randomize(2);
			if (choice==0)  return 'A';
			else return '-';
			break;
		case 'F':
			/*return 'C'; comment out to preserve gap possibility*/
			choice=my_randomize(2);
			if (choice==0)  return 'C';
			else return '-';
			break;
		case 'I':
			/*return 'G'; comment out to preserve gap possibility*/
			choice=my_randomize(2);
			if (choice==0)  return 'G';
			else return '-';
			break;
		case 'J':
			/*return 'T'; comment out to preserve gap possibility*/
			choice=my_randomize(2);
			if (choice==0)  return 'T';
			else return '-';
			break;
		case 'L':
			/*return 'T'; comment out to preserve gap possibility*/
			choice=my_randomize(2);
			if (choice==0)  return 'T';
			else return '-';
			break;
		case 'O':
			choice=my_randomize(5);	/*ERMOVE THESE AND FOLLOWING -1 TO ALLOW GAPS TO BE CHOSEN*/
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else if (choice==2)  return 'G';
			else if (choice==3)  return 'T';
			else return '-';
			break;
		case 'P':
			choice=my_randomize(5);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else if (choice==2)  return 'G';
			else if (choice==3)  return 'T';
			else return '-';
			break;
		case 'Q':
			choice=my_randomize(3);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'G';
			else return '-';
			break;
		case 'Z':
			choice=my_randomize(3);
			if (choice==0)  return 'C';
			else if (choice==1)  return 'T';
			else return '-';
			break;
		case '1':
			choice=my_randomize(3);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else return '-';
			break;
		case '2':
			choice=my_randomize(3);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'T';
			else return '-';
			break;
		case '3':
			choice=my_randomize(3);
			if (choice==0)  return 'C';
			else if (choice==1)  return 'G';
			else return '-';
			break;
		case '4':
			choice=my_randomize(3);
			if (choice==0)  return 'G';
			else if (choice==1)  return 'T';
			else return '-';
			break;
		case '5':
			choice=my_randomize(4);
			if (choice==0)  return 'C';
			else if (choice==1)  return 'G';
			else if (choice==2)  return 'T';
			else return '-';
			break;
		case '6':
			choice=my_randomize(4);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'G';
			else if (choice==2)  return 'T';
			else return '-';
			break;
		case '7':
			choice=my_randomize(4);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else if (choice==2)  return 'T';
			else return '-';
			break;
		case '8':
			choice=my_randomize(4);
			if (choice==0)  return 'A';
			else if (choice==1)  return 'C';
			else if (choice==2)  return 'G';
			else return '-';
			break;

		default:
			fprintf (stderr, "Bad character %c\n",inchar);
			exit(-1);
		}
	}
else {
	switch (inchar) {
		case 'N': return 'A';	;	break;
		case 'X': return 'A';	;	break;
		case 'U': return 'T';	;	break;
		case 'R': return 'A';	;	break;
		case 'Y': return 'C';	;	break;
		case 'M': return 'A';	;	break;
		case 'W': return 'A';	;	break;
		case 'S': return 'C';	;	break;
		case 'K': return 'G';	;	break;
		case 'B': return 'C';	;	break;
		case 'D': return 'A';	;	break;
		case 'H': return 'A';	;	break;
		case 'V': return 'A';	;	break;
	/*
		case 'E': return '-';	;	break;
		case 'F': return '-';	;	break;
		case 'I': return '-';	;	break;
		case 'J': return '-';	;	break;
		case 'L': return '-';	;	break;
		case 'O': return '-';	;	break;
		case 'P': return '-';	;	break;
		case 'Q': return '-';	;	break;
		case 'Z': return '-';	;	break;
		case '1': return '-';	;	break;
		case '2': return '-';	;	break;
		case '3': return '-';	;	break;
		case '4': return '-';	;	break;
		case '5': return '-';	;	break;
		case '6': return '-';	;	break;
		case '7': return '-';	;	break;
		case '8': return '-';	;	break;
	*/
		case 'E': return 'A';	;	break;
		case 'F': return 'C';	;	break;
		case 'I': return 'G';	;	break;
		case 'J': return 'T';	;	break;
		case 'L': return 'A';	;	break;
		case 'O': return 'A';	;	break;
		case 'P': return 'T';	;	break;
		case 'Q': return 'A';	;	break;
		case 'Z': return 'C';	;	break;
		case '1': return 'A';	;	break;
		case '2': return 'A';	;	break;
		case '3': return 'C';	;	break;
		case '4': return 'G';	;	break;
		case '5': return 'C';	;	break;
		case '6': return 'A';	;	break;
		case '7': return 'A';	;	break;
		case '8': return 'A';	;	break;

		default:
			fprintf (stderr, "Bad character %c\n",inchar);
			exit(-1);
		}
	}
return '\0';
}

int is_ambig(num)
int num;
{
if (num==A_BIT) return 0;
if (num==C_BIT) return 0;
if (num==G_BIT) return 0;
if (num==T_BIT) return 0;
if (num==GAP_BIT) return 0;
return 1;
}

int real_position(istring,place)
int place;
char *istring;
{
int i,j;

j=0;
for (i=0;i<=place;i++) if (istring[i]!='-') ++j;
return j;
}

char return_character(input)
int input;
{
return ((char) bitchar[input]);
}

apolist_holder *make_anc_up2(a,values,bl,ancestor,cur_node,descendent1,descendent2,task,apomorphy,node_number,return_holder)
alignment *a;
parameters *values;
int *bl,ancestor,cur_node,task,descendent1,descendent2;
apo_thang **apomorphy;
int node_number;
apolist_holder *return_holder;
{
char **nbsa;
int i,j, new_length, old_length, score_test,k0,k1,k2;
alignment *b;
int gaps=0, changes=0, ends=0;

b=NULL;
if  (a->n_seqs>1) {
	if (a->n_seqs!=4) {
		fprintf(stderr,"This is wrong, %d sequences when there should be 4\n",a->n_seqs);
		exit(-1);
		}
	nbsa=(char **)malloc((a->n_seqs+1)*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<(a->n_seqs+1);i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
	/*if not bitted already get_bases*/
	for (i=0;i<a->n_seqs;i++) for (j=0;j<a->length;j++) nbsa[i][j]=charbit[(int) a->s[i][j]];


	/*first set ancestor as best of descendents and ancestor*/
	ends=changes=gaps=new_length=score_test=0;
	bl[0]=bl[1]=bl[2]=bl[3]=0;
	for (j=0;j<a->length;j++)  {
		if (!is_ambig(nbsa[cur_node][j])) nbsa[4][j]=nbsa[cur_node][j];
		else {
			nbsa[4][j]=(nbsa[1][j]&nbsa[2][j]&nbsa[3][j]);
			if (nbsa[4][j]) { if (is_ambig(nbsa[4][j])) nbsa[4][j]=choose_one(nbsa[4][j],values->rand_apo,values); }
			else nbsa[4][j]=get_best_state(nbsa[0][j],nbsa[1][j],nbsa[2][j],nbsa[3][j],values->rand_apo,values);
			}
		if (nbsa[cur_node][j]!=GAP_BIT) ++new_length;
		}
	for (j=0;j<a->length;j++)  {
		if (nbsa[ancestor][j]!=nbsa[4][j]) {
			if ((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[4][j]==GAP_BIT))) {
				if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) ++gaps;
				else ++ends;
				}
			else ++changes;
			}
		}
	if (values->apolist) {
		apomorphy[0]=(apo_thang *)malloc(sizeof(apo_thang));
		assert((int) apomorphy[0]);
		apomorphy[1]=(apo_thang *)malloc(sizeof(apo_thang));
		assert((int) apomorphy[1]);
		apomorphy[2]=(apo_thang *)malloc(sizeof(apo_thang));
		assert((int) apomorphy[2]);
		k0=k1=k2=0;
		for (j=0;j<a->length;j++)  {
			if (nbsa[ancestor][j]!=nbsa[4][j]) {
				if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) {
					if ((nbsa[ancestor][j]==GAP_BIT) || (nbsa[4][j]==GAP_BIT)) {
						apomorphy[1]=(apo_thang *)realloc(apomorphy[1],(k1+1)*sizeof(apo_thang));
						assert((int) apomorphy[1]);
						apomorphy[1][k1].anc_pos=real_position(a->s[ancestor],j);
						apomorphy[1][k1].desc_pos=real_position(a->s[cur_node],j);
						apomorphy[1][k1].anc_char=a->s[ancestor][j];
						apomorphy[1][k1].desc_char=return_character(nbsa[4][j]);
						++k1;
						}
					else {
						apomorphy[2]=(apo_thang *)realloc(apomorphy[2],(k2+1)*sizeof(apo_thang));
						assert((int) apomorphy[2]);
						apomorphy[2][k2].anc_pos=real_position(a->s[ancestor],j);
						apomorphy[2][k2].desc_pos=real_position(a->s[cur_node],j);							apomorphy[2][k2].anc_char=a->s[ancestor][j];
						apomorphy[2][k2].anc_char=a->s[ancestor][j];
						apomorphy[2][k2].desc_char=return_character(nbsa[4][j]);
						++k2;
						}
					}
				else {
					if ((nbsa[ancestor][j]==GAP_BIT) || (nbsa[2][j]==GAP_BIT)) {
						apomorphy[0]=(apo_thang *)realloc(apomorphy[0],(k0+1)*sizeof(apo_thang));
						assert((int) apomorphy[0]);
						apomorphy[0][k0].anc_pos=real_position(a->s[ancestor],j);
						apomorphy[0][k0].desc_pos=real_position(a->s[cur_node],j); 							apomorphy[0][k0].anc_char=a->s[ancestor][j];
						apomorphy[0][k0].anc_char=a->s[ancestor][j];
						apomorphy[0][k0].desc_char=return_character(nbsa[4][j]);
						++k0;
						}
					else {
						apomorphy[2]=(apo_thang *)realloc(apomorphy[2],(k2+1)*sizeof(apo_thang));
						assert((int) apomorphy[2]);
						apomorphy[2][k2].anc_pos=real_position(a->s[ancestor],j);
						apomorphy[2][k2].desc_pos=real_position(a->s[cur_node],j);
						apomorphy[2][k2].anc_char=a->s[ancestor][j];
						apomorphy[2][k2].desc_char=return_character(nbsa[4][j]);
						++k2;
						}
					}
				}
			}
		if (k0!=ends) fprintf(stderr,"Error in apomorphy list 0 (%d %d)\n",k0,ends);
		/*for (i=0;i<k0;i++) fprintf(stderr,"[%d][0][%d]=%d %c %d %c ",node_number,i,apomorphy[0][i].anc_pos,apomorphy[0][i].anc_char,apomorphy[0][i].desc_pos,apomorphy[0][i].desc_char);*/
		if (k1!=gaps) fprintf(stderr,"Error in apomorphy list 1 (%d %d)\n",k1,gaps);
		/*for (i=0;i<k1;i++) fprintf(stderr,"[%d][1][%d]=%d  %c %d %c ",node_number,i,apomorphy[1][i].anc_pos,apomorphy[1][i].anc_char,apomorphy[1][i].desc_pos,apomorphy[1][i].desc_char);*/
		if (k2!=changes) fprintf(stderr,"Error in apomorphy list 2 (%d %d)\n",k2,changes);
		/*for (i=0;i<k2;i++) fprintf(stderr,"[%d][2][%d]=%d  %c %d %c ",node_number,i,apomorphy[2][i].anc_pos,apomorphy[2][i].anc_char,apomorphy[2][i].desc_pos,apomorphy[2][i].desc_char);
		fprintf(stderr,"end\n");*/
		}/*end of apolist*/

	a->score=(gaps*values->gap_cost);
	a->score+=(changes*values->change_cost);
	a->score-=(ends*values->leading_gap_cost);

	/*set branch_lengths*/
	bl[0]=ends;
	bl[1]=gaps;
	bl[2]=changes;
	bl[3]=((bl[0]*(values->gap_cost-values->leading_gap_cost))+(bl[1]*values->gap_cost)+(bl[2]*values->change_cost));
		/*fprintf(stderr,"%d %d %d, ",ends,gaps,changes);*/
	/*copy back*/
	old_length=a->length;
	if (values->gap_in_hypanc) new_length=a->length;
	b=(alignment *)malloc(sizeof(alignment));
	assert ((int) b);
	b->n_seqs=1;
	b->length=new_length;
	b->type_weight=a->type_weight;
	b->score=a->score; /*score_test; b->score;*/
	b->name=(char *)malloc((strlen(a->name)+1)*sizeof(char));
	assert((int) b->name);
	b->name=(char *)strcpy(b->name,a->name);
	b->s=(char **)malloc((1+b->type_weight)*sizeof(char *));
	assert((int) b->s);
	for (i=0;i<(b->n_seqs+b->type_weight);i++) {
		b->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) b->s[i]);
		b->s[i][new_length]='\0';
		}
	if (b->type_weight) b->s[b->n_seqs]=strcpy(b->s[b->n_seqs],a->s[a->n_seqs]);
	b->taxon_name=(char **)malloc(sizeof(char *));
	assert((int) b->taxon_name);
	b->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) b->taxon_name[0]);
	b->taxon_name[0][0]='H';
	b->taxon_name[0][1]='y';
	b->taxon_name[0][2]='p';
	b->taxon_name[0][3]='A';
	b->taxon_name[0][4]='n';
	b->taxon_name[0][5]='c';
	b->taxon_name[0][6]='\0';
	j=0;
	if (!values->gap_in_hypanc) {
		for (i=0;i<old_length;i++) if (nbsa[cur_node][i]!=GAP_BIT) { b->s[0][j]=bitchar[nbsa[4][j]];
			j++;
			}
		}
	else {
		for (i=0;i<old_length;i++) { b->s[0][j]=bitchar[nbsa[4][j]];
			j++;
			}
	}
	/*free*/
	for (i=0;i<a->n_seqs+1;i++) free(nbsa[i]);
	free(nbsa);
	}
return_holder->z=b;
if ((task==1) && (values->apolist)) {
	return_holder->at0=apomorphy[0];
	return_holder->at1=apomorphy[1];
	return_holder->at2=apomorphy[2];
	}
else {
	return_holder->at0=NULL;
	return_holder->at1=NULL;
	return_holder->at2=NULL;
	}
return return_holder;
}

/*branch lengths removed*/
apolist_holder *make_anc_up3(a,values,bl,ancestor,descendent,desc0, desc1,task,apomorphy,node_number,return_holder)
alignment *a;
parameters *values;
int *bl,ancestor,descendent,task;
apo_thang **apomorphy;
int node_number;
apolist_holder *return_holder;
{
char **nbsa;
int i,j, new_length, old_length, score_test,k0,k1,k2;
alignment *b;
int gaps=0, changes=0, ends=0;

b=NULL;
if  (a->n_seqs>1) {
	if (a->n_seqs!=4) {
		fprintf(stderr,"This is wrong, %d sequences when there should be 4\n",a->n_seqs);
		exit(-1);
		}
	nbsa=(char **)malloc((a->n_seqs+1)*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<(a->n_seqs+1);i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
	/*if not bitted already get_bases*/
	for (i=0;i<a->n_seqs;i++) for (j=0;j<a->length;j++)  nbsa[i][j]=bitchar[a->s[i][j]];


	/*assume two now and get union/intersection*/
	ends=changes=gaps=new_length=score_test=0;
	bl[0]=bl[1]=bl[2]=bl[3]=0;
	if (task==0) { /*min branch lengths*//*redo for all MPR brute force*/
		for (j=0;j<a->length;j++)  {
			nbsa[4][j]=(nbsa[1][j]&nbsa[2][j]&nbsa[3][j]);
			if (!nbsa[4][j]) nbsa[4][j]=get_all_states(nbsa[0][j],nbsa[1][j],nbsa[2][j],nbsa[3][j],values->rand_apo,values);
			if (!(nbsa[ancestor][j]&nbsa[4][j])) {/*min gaps*/
				if ((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[4][j]==GAP_BIT))) {
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) ++gaps;
					else ++ends;
					}
				}
			if (!(nbsa[ancestor][j]&nbsa[4][j])) {/*min chnages*/
				if ((nbsa[ancestor][j]&GAP_BIT) || ((nbsa[4][j]&GAP_BIT))) {
					}
				else ++changes;
				}
			/*if (!(nbsa[ancestor][j]&nbsa[4][j])) {
				if (!((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[4][j]==GAP_BIT))))  bl[3]+=(values->change_cost);
				else{
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) bl[3]+=(values->gap_cost);
					else bl[3]+=(values->gap_cost-values->leading_gap_cost);
					}
				}*/
			if (nbsa[descendent][j]!=GAP_BIT) ++new_length;
			}
		}
	else if (task==3) { /* max branch length*/
		for (j=0;j<a->length;j++)  {
			nbsa[4][j]=(nbsa[1][j]&nbsa[2][j]&nbsa[3][j]);
			if (!nbsa[4][j]) nbsa[4][j]=get_all_states(nbsa[0][j],nbsa[1][j],nbsa[2][j],nbsa[3][j],values->rand_apo,values);
			if (nbsa[ancestor][j]!=nbsa[4][j]) {
				if ((nbsa[ancestor][j]&GAP_BIT) || ((nbsa[4][j]&GAP_BIT))) {/*max gaps*/
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) ++gaps;
					else ++ends;
					}
				}
			if (nbsa[ancestor][j]!=nbsa[4][j]) {/*max changes*/
				if ((nbsa[ancestor][j]==GAP_BIT) || ((nbsa[4][j]==GAP_BIT))) {
					}
				else ++changes;
				}
			if (nbsa[ancestor][j]!=nbsa[4][j]) {
				if ((nbsa[ancestor][j]&GAP_BIT) || ((nbsa[4][j]&GAP_BIT))) {
					if ((is_internal(a->s[0],j,a->length)) && (is_internal(a->s[1],j,a->length))) bl[3]+=(values->gap_cost);
					else bl[3]+=(values->gap_cost-values->leading_gap_cost);
					}
				else bl[3]+=(values->change_cost);;
				}
			if (nbsa[descendent][j]!=GAP_BIT) ++new_length;
			}
		}

	/*a->score=(gaps*values->gap_cost);
	a->score+=(changes*values->change_cost);
	a->score-=(ends*values->leading_gap_cost);*/
	/*set branch_lengths*/
/*	bl[0]=ends;
	bl[1]=gaps;
	bl[2]=changes;*/
	if (task==1) bl[3]=((bl[0]*(values->gap_cost-values->leading_gap_cost))+(bl[1]*values->gap_cost)+(bl[2]*values->change_cost));
	/*copy back*/
	old_length=a->length;
	if (values->gap_in_hypanc) new_length=a->length;
	b=(alignment *)malloc(sizeof(alignment));
	assert ((int) b);
	b->n_seqs=1;
	b->type_weight=a->type_weight;
	b->length=new_length;
	b->score=a->score; /*score_test; b->score;*/
	b->name=(char *)malloc((strlen(a->name)+1)*sizeof(char));
	assert((int) b->name);
	b->name=(char *)strcpy(b->name,a->name);
	b->s=(char **)malloc((1+b->type_weight)*sizeof(char *));
	assert((int) b->s);
	for (i=0;i<(b->n_seqs+b->type_weight);i++) {
		b->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) b->s[i]);
		b->s[i][new_length]='\0';
		}
	if (b->type_weight) b->s[b->n_seqs]=strcpy(b->s[b->n_seqs],a->s[a->n_seqs]);
	b->taxon_name=(char **)malloc(sizeof(char *));
	assert((int) b->taxon_name);
	b->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) b->taxon_name[0]);
	b->taxon_name[0][0]='H';
	b->taxon_name[0][1]='y';
	b->taxon_name[0][2]='p';
	b->taxon_name[0][3]='A';
	b->taxon_name[0][4]='n';
	b->taxon_name[0][5]='c';
	b->taxon_name[0][6]='\0';
	j=0;
	if (!values->gap_in_hypanc) {
		for (i=0;i<old_length;i++) if (nbsa[descendent][i]!=GAP_BIT) { b->s[0][j]=bitchar[nbsa[4][i]];
			j++;
			}
		}
	else {
		for (i=0;i<old_length;i++) { b->s[0][j]=bitchar[nbsa[4][i]];
			j++;
			}
	}
	/*free*/
	for (i=0;i<(a->n_seqs+1);i++) free(nbsa[i]);
	free(nbsa);
	}
return_holder->z=b;
if ((task==1) && (values->apolist)) {
	return_holder->at0=apomorphy[0];
	return_holder->at1=apomorphy[1];
	return_holder->at2=apomorphy[2];
	}
else {
	return_holder->at0=NULL;
	return_holder->at1=NULL;
	return_holder->at2=NULL;
	}
return return_holder;
}

alignment *make_ambig(a,values)
alignment *a;
parameters *values;
{
char **nbsa;
int i,j, new_length, old_length, score_test;
alignment *b;
int gaps, changes, ends;
int was_union,allgaps,pieces;

b=NULL;
if  (a) {
	b=make_align(a);
	if (a->n_seqs!=2) {
		fprintf(stderr,"This is wrong, %d sequences when there should be 2\n",a->n_seqs);
		exit(-1);
		}
	nbsa=(char **)malloc(3*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<3;i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
	/*if not bitted already get_bases*/
	for (i=0;i<2;i++) {
		for (j=0;j<a->length;j++) {
			if ((a->type_weight) && (a->s[2][j]!=1) && (a->s[2][j]!=65)) nbsa[i][j]=a->s[i][j];
			else { nbsa[i][j]=charbit[a->s[i][j]];
			}
			}
		}
	/*if (a->type_weight) {
		pieces=0;
		for (i=0;i<a->length;i++) if (a->s[2][i]>63) ++pieces;
		fprintf(stderr,"(%d in",pieces);
		}*/
	a=dump_align(a);
	a=NULL;

	/*assume two now and get union/intersection*/
	new_length=score_test=0;
	if ((!values->no_leading_or_trailing_cost) && (values->phylo_gap)) {
		for (j=0;j<b->length;j++)  {
			if ((b->type_weight) && (b->s[b->n_seqs][j]!=1) && (b->s[b->n_seqs][j]!=65)) {
			    	++new_length;
            			nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
    				if (!nbsa[2][j]) {
	    				was_union=1;
		    			nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
			    		}
			    	else was_union=0;
			   	 }
			else {
			    nbsa[2][j]=values->lookup[nbsa[0][j]][nbsa[1][j]].base;
			    was_union=values->lookup[nbsa[0][j]][nbsa[1][j]].union_bit;
			    if (nbsa[2][j]!=GAP_BIT) ++new_length;
    			}
			if (values->jack_array && was_union) {
				if (!values->jack_array[(j*1000)/b->length]) {
					if ((b->type_weight) && (b->s[b->n_seqs][j]==3) && (b->s[b->n_seqs][j]==67)) b->score-=values->change_cost; /*morph*/
					else b->score-=values->lookup[nbsa[0][j]][nbsa[1][j]].cost;
					}
				}
			}/*J*/
		}
	else {
		gaps=ends=changes=0;
		for (j=0;j<b->length;j++)  {
			if ((b->type_weight) && (b->s[b->n_seqs][j]!=1) && (b->s[b->n_seqs][j]!=65)) {
			    ++new_length;
            	nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
    			if (!nbsa[2][j]) {
	    			was_union=1;
		    		nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
			    	}
			    else was_union=0;
			    }
			else {
			    nbsa[2][j]=values->lookup[nbsa[0][j]][nbsa[1][j]].base;
			    was_union=values->lookup[nbsa[0][j]][nbsa[1][j]].union_bit;
				if (nbsa[2][j]!=GAP_BIT) ++new_length;
    			}
		if (was_union) {
			if (values->jack_array) {
				if (!values->jack_array[(j*1000)/b->length]) {
					if ((b->type_weight) && (b->s[b->n_seqs][j]==3) && (b->s[b->n_seqs][j]==67)) b->score-=values->change_cost; /*morph*/
					else  b->score-=values->lookup[nbsa[0][j]][nbsa[1][j]].cost;
					}
				else {
					if ((is_internal(b->s[0],j,b->length)) && (is_internal(b->s[1],j,b->length))) {
						if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++gaps;
						else ++changes;
						}
					else if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++ends;
					}
				}
			else {
				if ((is_internal(b->s[0],j,b->length)) && (is_internal(b->s[1],j,b->length))) {
					if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++gaps;
					else ++changes;
					}
				else if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++ends;
				}
			}
		}
		if (values->no_leading_or_trailing_cost || (!values->phylo_gap)) b->score-=(ends*(values->gap_cost-values->leading_gap_cost));
		if (!values->phylo_gap) b->score-=(gaps*values->gap_cost);
		}
	
	/*filter for zero length pieces to be*/
	if (b->type_weight) {
    		if (nbsa[2][0]==GAP_BIT) allgaps=1;
    		else allgaps=0;
    		for (i=1;i<b->length;i++) {
	        	if (b->s[2][i]==67) allgaps=1;
 			else if (b->s[2][i]>63){
           			if (allgaps) {
                			if (nbsa[2][i]==GAP_BIT) {
                   				nbsa[2][i]=X_BIT;
                    				++new_length;
                				}
            				}
            			else allgaps=1;
        			}
       			 else if ((b->s[2][i]==3) || (nbsa[2][i]!=GAP_BIT)) allgaps=0;
    		}
}/*copy back*/
	old_length=b->length;
	if (values->gap_in_hypanc) new_length=b->length;
	a=(alignment *)malloc(sizeof(alignment));
	assert ((int) a);
	a->n_seqs=1;
	if (b->type_weight!=1) b->type_weight=0;
	a->type_weight=b->type_weight;
	a->length=new_length;
	a->score=b->score; /*score_test; b->score;*/
	a->name=(char *)malloc((strlen(b->name)+1)*sizeof(char));
	assert((int) a->name);
	a->name=(char *)strcpy(a->name,b->name);
	a->s=(char **)malloc((1+a->type_weight)*sizeof(char *));
	assert((int) a->s);
	for (i=0;i<(a->n_seqs+a->type_weight);i++) {
		a->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) a->s[i]);
		a->s[i][new_length]='\0';
		}
	a->taxon_name=(char **)malloc(sizeof(char *));
	assert((int) a->taxon_name);
	a->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) a->taxon_name[0]);
	a->taxon_name[0][0]='H';
	a->taxon_name[0][1]='y';
	a->taxon_name[0][2]='p';
	a->taxon_name[0][3]='A';
	a->taxon_name[0][4]='n';
	a->taxon_name[0][5]='c';
	a->taxon_name[0][6]='\0';
	j=0;
	if (!values->gap_in_hypanc) {
		for (i=0;i<old_length;i++) {
			if ((a->type_weight) && (b->s[2][i]!=1) && (b->s[2][i]!=65)) {
				a->s[0][j]=nbsa[2][i];
				a->s[1][j]=b->s[2][i];
				++j;
				}
			else if (nbsa[2][i]==GAP_BIT) {
				if (a->type_weight) if (b->s[2][i]>=64) a->s[1][j-1]=b->s[2][i];
				}
			else  {
				if (a->type_weight) a->s[1][j]=b->s[2][i];
				a->s[0][j]=bitchar[nbsa[2][i]];
				j++;
			 	}
			}
		}
	else {
		for (i=0;i<old_length;i++) {
			if ((a->type_weight) && (b->s[2][i]!=1) && (b->s[2][i]!=65)) {
				a->s[0][j]=nbsa[2][i];
				a->s[1][j]=b->s[2][i];
				}
			else {
				if (a->type_weight) a->s[1][j]=b->s[2][i];
				a->s[0][j]=bitchar[nbsa[2][i]];
			}
			j++;
			}
	}
	/*free*/
	for (i=0;i<3;i++) free(nbsa[i]);
	free(nbsa);
	}
else {
	fprintf(stderr,"Alignment is not in make_ambig\n");
	exit(-1);
	}
/*
print_alignment(a);

fprintf(stderr,"score %d vs. %d (%d gaps, %d ends %d changes)\n",a->score, score_test, gaps, ends, changes);
*/
b=dump_align(b);
/*if (a->type_weight) {
	pieces=0;
	for (i=0;i<a->length;i++) if (a->s[1][i]>63) ++pieces;
	fprintf(stderr,"%d out)",pieces);
	}*/
return a;
}

/*this assumes 2 sequences but leaves gaps*/
alignment *new_make_ambig(a,values)/*used in apolist uppass and dead areas of build code*/
alignment *a;
parameters *values;
{
char **nbsa;
int i,j, new_length, old_length, score_test;
alignment *b;
int gaps=0, changes=0, ends=0,was_union;

if (!a) {fprintf(stderr,"Alignment is not in new_make_ambig\n"); exit(-1);}
if (a->n_seqs!=2) {fprintf(stderr,"Alignment with wrong num seqs in new_make_ambig(%d)\n",a->n_seqs); exit(-1);}
b=NULL;
	b=make_align(a);
	nbsa=(char **)malloc(3*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<3;i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
for (i=0;i<2;i++) {
	for (j=0;j<a->length;j++) {
		if (a->type_weight && (a->s[2][j]!=1) && (a->s[2][j]!=65)) nbsa[i][j]=a->s[i][j];
		else {
		    nbsa[i][j]=charbit[a->s[i][j]];
			}
		}
	}

	new_length=score_test=0;
if ((!values->no_leading_or_trailing_cost) && (values->phylo_gap)) {
	for (j=0;j<b->length;j++)  {
		was_union=0;
		nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
		if (!nbsa[2][j]) {
			nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
			was_union=1;
			}
		if (b->type_weight && ((b->s[2][j]==1)||(b->s[2][j]==65))) {
		   if (was_union) nbsa[2][j]=values->lookup[nbsa[0][j]][nbsa[1][j]].base;
    	}
		++new_length;
		}
		}
	else {
		gaps=ends=changes=0;
		for (j=0;j<b->length;j++)  {
			was_union=0;
			nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
			if (!nbsa[2][j]) {
				nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
				was_union=1;
				}
			if ((b->type_weight) && (b->s[b->n_seqs][j]!=1) && (b->s[b->n_seqs][j]!=65)) ++new_length;
			else {
			    if (was_union) nbsa[2][j]=values->lookup[nbsa[0][j]][nbsa[1][j]].base;
			++new_length;
			}
			if (was_union) {
				if ((is_internal(b->s[0],j,b->length)) && (is_internal(b->s[1],j,b->length))) {
					if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++gaps;
					else ++changes;
					}
				else if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++ends;
				}
			}/*j*/
		if (values->no_leading_or_trailing_cost || (!values->phylo_gap)) b->score-=(ends*(values->gap_cost-values->leading_gap_cost));
		if (!values->phylo_gap) b->score-=(gaps*values->gap_cost);
	}
	new_length=b->length;
	old_length=b->length;
	for (i=0;i<a->n_seqs+a->type_weight;i++) free(a->s[i]);
	free(a->s);
	free(a->name);
	for (i=0;i<a->n_seqs;i++) free(a->taxon_name[i]);
	free(a->taxon_name);

	a->n_seqs=1;
	a->type_weight=b->type_weight;
	a->length=new_length;
	a->score=b->score;
	a->name=(char *)malloc((strlen(b->name)+1)*sizeof(char));
	assert((int) a->name);
	a->name=(char *)strcpy(a->name,b->name);
	a->s=(char **)malloc((1+a->type_weight)*sizeof(char *));
	assert((int) a->s);
	for (i=0;i<(a->n_seqs+a->type_weight);i++) {
		a->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) a->s[i]);
		a->s[i][new_length]='\0';
		}
	if (a->type_weight) a->s[a->n_seqs]=strcpy(a->s[a->n_seqs],b->s[b->n_seqs]);
	a->taxon_name=(char **)malloc(sizeof(char *));
	assert((int) a->taxon_name);
	a->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) a->taxon_name[0]);
	a->taxon_name[0][0]='H';
	a->taxon_name[0][1]='y';
	a->taxon_name[0][2]='p';
	a->taxon_name[0][3]='A';
	a->taxon_name[0][4]='n';
	a->taxon_name[0][5]='c';
	a->taxon_name[0][6]='\0';
	j=0;
		for (i=0;i<old_length;i++) {
			if (a->type_weight && (a->s[1][i]!=1) && (a->s[1][i]!=65)) a->s[0][j]=nbsa[2][j];
			else {
			    a->s[0][j]=bitchar[nbsa[2][i]];
				}
			j++;
			}
	for (i=0;i<3;i++) free(nbsa[i]);
	free(nbsa);
	b=dump_align(b);
return a;
}


/*this assumes 2 sequences but his only does unions*/
alignment *newer_make_ambig(a,values)
alignment *a;
parameters *values;
{
char **nbsa;
int i,j, new_length, old_length, score_test;
alignment *b;
int gaps=0, changes=0, ends=0;
int was_union;

if (!a) {fprintf(stderr,"Alignment is not in newer_make_ambig\n"); exit(-1);}
if (a->n_seqs!=2) {fprintf(stderr,"Alignment with wrong num seqs in newer_make_ambig(%d)\n",a->n_seqs); exit(-1);}

b=NULL;
	b=make_align(a);
	nbsa=(char **)malloc(3*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<3;i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
for (i=0;i<2;i++) {
	for (j=0;j<a->length;j++) {
		if (a->type_weight && (a->s[2][j]!=1) && (a->s[2][j]!=65)) nbsa[i][j]=a->s[i][j];
		else {
		    nbsa[i][j]=charbit[a->s[i][j]];
			}
		}
	}

	new_length=score_test=0;
if ((!values->no_leading_or_trailing_cost) && (values->phylo_gap)) {
	for (j=0;j<b->length;j++)  {
		was_union=0;
		/*this needs to be changed for missing in a terminal taxon*/
		if (a->type_weight) {
			if ((a->s[a->n_seqs][j]==1) || (a->s[a->n_seqs][j]==65) || (a->s[a->n_seqs][j]==2) || (a->s[a->n_seqs][j]==66)) {
				nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
				if (nbsa[2][j]==X_BIT) if ((nbsa[0][j]==X_BIT) || (nbsa[1][j]==X_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				else if (nbsa[2][j]==N_BIT) if ((nbsa[0][j]==N_BIT) || (nbsa[1][j]==N_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				}
			else {
				if ((nbsa[0][j]==127) || (nbsa[1][j]==127)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				else if ((nbsa[0][j]==127) || (nbsa[1][j]==127)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				else nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
				}
			}
		else {
			nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
			if (nbsa[2][j]==X_BIT) if ((nbsa[0][j]==X_BIT) || (nbsa[1][j]==X_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				else if (nbsa[2][j]==N_BIT) {
					if ((nbsa[0][j]==N_BIT) || (nbsa[1][j]==N_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
					}
			}
		++new_length;
		}
	}
else {
	for (j=0;j<b->length;j++)  {
		/*this needs to be changed for missing in a terminal taxon*/
		if (a->type_weight) {
			if ((a->s[a->n_seqs][j]==1) || (a->s[a->n_seqs][j]==65) || (a->s[a->n_seqs][j]==2) || (a->s[a->n_seqs][j]==66)) {
				nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
				if (nbsa[2][j]==X_BIT) {
					if ((nbsa[0][j]==X_BIT) || (nbsa[1][j]==X_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
					}
				else if (nbsa[2][j]==N_BIT) {
					if ((nbsa[0][j]==N_BIT) || (nbsa[1][j]==N_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
					}
				}
			else {
				if ((nbsa[0][j]==127) || (nbsa[1][j]==127)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				else if ((nbsa[0][j]==127) || (nbsa[1][j]==127)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				else nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
				}
			}
		else {
			nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
			if (nbsa[2][j]==X_BIT) {
				if ((nbsa[0][j]==X_BIT) || (nbsa[1][j]==X_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				}
			else if (nbsa[2][j]==N_BIT) {
				if ((nbsa[0][j]==N_BIT) || (nbsa[1][j]==N_BIT)) nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
				}
			}
		++new_length;
		if (was_union) {
				if ((is_internal(b->s[0],j,b->length)) && (is_internal(b->s[1],j,b->length))) {
					if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++gaps;
					else ++changes;
					}
				else if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++ends;
				}
		}/*j*/
		if (values->no_leading_or_trailing_cost || (!values->phylo_gap)) b->score-=(ends*(values->gap_cost-values->leading_gap_cost));
		if (!values->phylo_gap) b->score-=(gaps*values->gap_cost);
	}
	old_length=b->length;
	a=dump_align(a);
	a=(alignment *)malloc(sizeof(alignment));
	assert((int) a);
	a->n_seqs=1;
	a->type_weight=b->type_weight;
	a->length=new_length;
	a->score=b->score;
	a->name=(char *)malloc((strlen(b->name)+1)*sizeof(char));
	assert((int) a->name);
	a->name=(char *)strcpy(a->name,b->name);
	a->s=(char **)malloc((1+a->type_weight)*sizeof(char *));
	assert((int) a->s);
	for (i=0;i<(a->n_seqs+a->type_weight);i++) {
		a->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) a->s[i]);
		a->s[i][new_length]='\0';
		}
	if (a->type_weight) a->s[a->n_seqs]=strcpy(a->s[a->n_seqs],b->s[b->n_seqs]);
	a->taxon_name=(char **)malloc(sizeof(char *));
	assert((int) a->taxon_name);
	a->taxon_name[0]=(char *)malloc((6+1)*sizeof(char));
	assert ((int) a->taxon_name[0]);
	a->taxon_name[0][0]='H';
	a->taxon_name[0][1]='y';
	a->taxon_name[0][2]='p';
	a->taxon_name[0][3]='A';
	a->taxon_name[0][4]='n';
	a->taxon_name[0][5]='c';
	a->taxon_name[0][6]='\0';
	j=0;
		for (i=0;i<old_length;i++) {
			if (a->type_weight && (a->s[1][i]!=1) && (a->s[1][i]!=65)) a->s[0][j]=nbsa[2][j];
			else {
			    a->s[0][j]=bitchar[nbsa[2][j]];
				}
			j++;
			}
	for (i=0;i<3;i++) free(nbsa[i]);
	free(nbsa);
	b=dump_align(b);
return a;
}

/*this assumes 2 sequences but adds the result--leaves gaps*/
alignment *make_three(a,values)
alignment *a;
parameters *values;
{
char **nbsa;
int i,j, new_length, old_length, score_test;
alignment *b;
int gaps=0, changes=0, ends=0,was_union;

if (!a) {fprintf(stderr,"Alignment is not in new_make_ambig\n"); exit(-1);}
if (a->n_seqs!=2) {fprintf(stderr,"Alignment with wrong num seqs in new_make_ambig(%d)\n",a->n_seqs); exit(-1);}
b=NULL;
	b=make_align(a);
	nbsa=(char **)malloc(3*sizeof(char *));
	assert((int)nbsa);
	for (i=0;i<3;i++) {
		nbsa[i]=(char *)malloc(a->length*sizeof(char));
		assert((int)nbsa[i]);
		}
for (i=0;i<2;i++) {
	for (j=0;j<a->length;j++) {
		if (a->type_weight && (a->s[2][j]!=1) && (a->s[2][j]!=65)) nbsa[i][j]=a->s[i][j];
		else {
		    nbsa[i][j]=charbit[a->s[i][j]];
			}
		}
	}

	new_length=score_test=0;
if ((!values->no_leading_or_trailing_cost) && (values->phylo_gap)) {
	for (j=0;j<b->length;j++)  {
		was_union=0;
		nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
		if (!nbsa[2][j]) {
			nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
			was_union=1;
			}
		if (b->type_weight && ((b->s[2][j]==1)||(b->s[2][j]==65))) {
		   if (was_union) nbsa[2][j]=values->lookup[nbsa[0][j]][nbsa[1][j]].base;
    	}
		++new_length;
		}
		}
	else {
		gaps=ends=changes=0;
		for (j=0;j<b->length;j++)  {
			was_union=0;
			nbsa[2][j]=(nbsa[0][j]&nbsa[1][j]);
			if (!nbsa[2][j]) {
				nbsa[2][j]=(nbsa[0][j]|nbsa[1][j]);
				was_union=1;
				}
			if ((b->type_weight) && (b->s[b->n_seqs][j]!=1) && (b->s[b->n_seqs][j]!=65)) ++new_length;
			else {
			    if (was_union) nbsa[2][j]=values->lookup[nbsa[0][j]][nbsa[1][j]].base;
			++new_length;
			}
			if (was_union) {
				if ((is_internal(b->s[0],j,b->length)) && (is_internal(b->s[1],j,b->length))) {
					if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++gaps;
					else ++changes;
					}
				else if ((nbsa[0][j]==GAP_BIT) || (nbsa[1][j]==GAP_BIT)) ++ends;
				}
			}/*j*/
		if (values->no_leading_or_trailing_cost || (!values->phylo_gap)) b->score-=(ends*(values->gap_cost-values->leading_gap_cost));
		if (!values->phylo_gap) b->score-=(gaps*values->gap_cost);
	}
	new_length=b->length;
	old_length=b->length;
	for (i=0;i<a->n_seqs+a->type_weight;i++) free(a->s[i]);
	free(a->s);
	free(a->name);
	for (i=0;i<a->n_seqs;i++) free(a->taxon_name[i]);
	free(a->taxon_name);

	a->n_seqs=3;
	a->type_weight=b->type_weight;
	a->length=new_length;
	a->score=b->score;
	a->name=(char *)malloc((strlen(b->name)+1)*sizeof(char));
	assert((int) a->name);
	a->name=(char *)strcpy(a->name,b->name);
	a->s=(char **)malloc((3+a->type_weight)*sizeof(char *));
	assert((int) a->s);
	for (i=0;i<(a->n_seqs+a->type_weight);i++) {
		a->s[i]=(char *)malloc((new_length+1)*sizeof(char));
		assert ((int) a->s[i]);
		a->s[i][new_length]='\0';
		}
	if (a->type_weight) a->s[a->n_seqs]=strcpy(a->s[a->n_seqs],b->s[b->n_seqs]);
	a->taxon_name=(char **)malloc(3*sizeof(char *));
	assert((int) a->taxon_name);
	for (i=0;i<(a->n_seqs);i++) {
    	a->taxon_name[i]=(char *)malloc((6+1)*sizeof(char));
    	assert ((int) a->taxon_name[i]);
    	a->taxon_name[i][0]='H';
    	a->taxon_name[i][1]='y';
    	a->taxon_name[i][2]='p';
    	a->taxon_name[i][3]='A';
    	a->taxon_name[i][4]='n';
    	a->taxon_name[i][5]='c';
    	a->taxon_name[i][6]='\0';
    }
	j=0;
	a->s[0]=(char *)strcpy(a->s[0],b->s[0]);
	a->s[1]=(char *)strcpy(a->s[1],b->s[1]);
	if (a->type_weight)	a->s[3]=(char *)strcpy(a->s[3],b->s[2]);

		for (i=0;i<old_length;i++) {
			if (a->type_weight && (a->s[3][i]!=1) && (a->s[3][i]!=65)) a->s[2][j]=nbsa[2][j];
			else {
			    a->s[2][j]=bitchar[nbsa[2][i]];
				}
			j++;
			}
	for (i=0;i<3;i++) free(nbsa[i]);
	free(nbsa);
	b=dump_align(b);
return a;
}

char *make_ambig_string(a,values,new_length,new_score)
char **a;
parameters *values;
int *new_length, *new_score;
{
char *out, *t_out;
int i,j,old_length;
int gaps, changes, ends;
int was_union,allgaps,pieces;

old_length=(*new_length);
t_out=(char *)malloc((1+old_length)*sizeof(char));
assert((int) t_out);
t_out[old_length]='\0';

(*new_length)=0;
for (i=0;i<old_length;i++) {
  t_out[i]=bitchar[values->lookup[charbit[a[0][i]]][charbit[a[1][i]]].base];
  /*fprintf(stderr,"%c ",t_out[i]);*/
  if (values->jack_array) {
    if (!values->jack_array[(i*1000)/old_length]) {
      (*new_score)-=values->lookup[charbit[a[0][i]]][charbit[a[1][i]]].cost;
    }
  }
  if (t_out[i]!= '-') ++(*new_length);
}
if ((*new_length)==0) allgaps=1;
else allgaps=0;
out=(char *)malloc((1+allgaps+(*new_length))*sizeof(char));
assert((int) out);
out[(*new_length)+allgaps]='\0';
if ((!values->no_leading_or_trailing_cost) && (values->phylo_gap)) { /*Normal run*/
  j=0;
  for (i=0;i<old_length;i++) if (t_out[i]!='-') out[j++]=t_out[i];
  if (allgaps) {
    ++(*new_length);
    out[0]='X';
  }
}
else if (!allgaps) { /*leading trailing etc removal*/
  gaps=ends=0;
  for (i=0;i<old_length;i++) {
   if (values->lookup[charbit[a[0][i]]][charbit[a[1][i]]].union_bit) {
     if (is_internal(a[0],i,old_length) && is_internal(a[1],i,old_length)) ++gaps;
     else ++ends;
     if (values->jack_array) {
       if (values->jack_array[(i*1000)/old_length]) {
	  if (is_internal(a[0],i,old_length) && is_internal(a[1],i,old_length)) --gaps;
	  else --ends;
       }
     }
   }
   *new_score-=(ends*(values->gap_cost-values->leading_gap_cost));
   if (values->phylo_gap) (*new_score)-=(gaps*values->gap_cost);
  }
}

free(t_out);
	
return out;
}

