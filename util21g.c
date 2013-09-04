/*
Copyright 1992 Ward Wheeler all rights reserved
*/
#include "align3.h"

extern int charbit[];
int **do_bit_thing(a,values)
alignment *a;
parameters *values;
{
int **nbsa;
int i,j;

if (!values->new_optimization) {
	/*
	This thing contains
		[0]=A
		[1]=C
		[2]=G
		[3]=T
		[4]=GAP
		[5]=down/right depending
	*/
	nbsa=(int **)malloc(a->length*sizeof(int *));
	assert((int)nbsa);
	for (i=0;i<a->length;i++) {
		nbsa[i]=(int *)malloc(6*sizeof(int));
		assert((int)nbsa[i]);
		}
	/*get_bases*/
	for (i=0;i<a->length;i++) {
		for (j=0;j<5;j++) nbsa[i][j]=0;
		for (j=0;j<a->n_seqs;j++) {
			switch (a->s[j][i]) {
			case 'A': nbsa[i][0]=1;	break;
			case 'C': nbsa[i][1]=1;	break;
			case 'G': nbsa[i][2]=1;	break;
			case 'T': nbsa[i][3]=1;	break;
			case '-': nbsa[i][4]=1;	break;
			case 'N': nbsa[i][0]=nbsa[i][1]=nbsa[i][2]=nbsa[i][3]=1;	break;
			case 'X': nbsa[i][0]=nbsa[i][1]=nbsa[i][2]=nbsa[i][3]=1;	break;
			case 'U': nbsa[i][3]=1;	break;
			case 'R': nbsa[i][0]=nbsa[i][2]=1;	break;
			case 'Y': nbsa[i][1]=nbsa[i][3]=1;	break;
			case 'M': nbsa[i][0]=nbsa[i][1]=1;	break;
			case 'W': nbsa[i][0]=nbsa[i][3]=1;	break;
			case 'S': nbsa[i][1]=nbsa[i][2]=1;	break;
			case 'K': nbsa[i][2]=nbsa[i][3]=1;	break;
			case 'B': nbsa[i][1]=nbsa[i][2]=nbsa[i][3]=1;	break;
			case 'D': nbsa[i][0]=nbsa[i][2]=nbsa[i][3]=1;	break;
			case 'H': nbsa[i][1]=nbsa[i][1]=nbsa[i][3]=1;	break;
			case 'V': nbsa[i][1]=nbsa[i][2]=nbsa[i][0]=1;	break;

			case 'E': nbsa[i][0]=nbsa[i][4]=1;	break;
			case 'F': nbsa[i][1]=nbsa[i][4]=1;	break;
			case 'I': nbsa[i][2]=nbsa[i][4]=1;	break;
			case 'J': nbsa[i][3]=nbsa[i][4]=1;	break;
			case 'L': nbsa[i][3]=nbsa[i][4]=1;	break;
			case 'O': nbsa[i][0]=nbsa[i][1]=nbsa[i][2]=nbsa[i][3]=nbsa[i][4]=1;	break;
			case 'P': nbsa[i][0]=nbsa[i][1]=nbsa[i][2]=nbsa[i][3]=nbsa[i][4]=1;	break;
			case 'Q': nbsa[i][0]=nbsa[i][2]=nbsa[i][4]=1;	break;
			case 'Z': nbsa[i][1]=nbsa[i][3]=nbsa[i][4]=1;	break;
			case '1': nbsa[i][0]=nbsa[i][1]=nbsa[i][4]=1;	break;
			case '2': nbsa[i][0]=nbsa[i][3]=nbsa[i][4]=1;	break;
			case '3': nbsa[i][1]=nbsa[i][2]=nbsa[i][4]=1;	break;
			case '4': nbsa[i][2]=nbsa[i][3]=nbsa[i][4]=1;	break;
		case '5': nbsa[i][1]=nbsa[i][2]=nbsa[i][3]=nbsa[i][4]=1;	break;
		case '6': nbsa[i][0]=nbsa[i][2]=nbsa[i][3]=nbsa[i][4]=1;	break;
		case '7': nbsa[i][1]=nbsa[i][1]=nbsa[i][3]=nbsa[i][4]=1;	break;
		case '8': nbsa[i][1]=nbsa[i][2]=nbsa[i][0]=nbsa[i][4]=1;	break;

		default:
			fprintf (stderr, "Bad character %c in alignment %s base [%d][%d]\n",a->s[j][i], a->name,j, i);
			print_inter_dig (a,values);
			exit(-1);
			}
		}
	}
/*create the down/right thing*/
/*nbsa[5] is the gap cost thing*/
if (values->gap_must_cost) for (i=0;i<a->length;i++) nbsa[i][5]= values->gap_cost * (1+nbsa[i][4]);
else for (i=0;i<a->length;i++) nbsa[i][5]= values->gap_cost;

if ((!values->delta) && (!values->ttr)) for (i=0;i<a->length;i++) nbsa[i][5]+=(values->change_cost * (nbsa[i][2]+nbsa[i][1]+nbsa[i][0]+nbsa[i][3]-1));
else if (values->ttr) for (i=0;i<a->length;i++) {
	nbsa[i][5] +=(values->transition+values->transversion)*((nbsa[i][0]&nbsa[i][1])+(nbsa[i][0]&nbsa[i][3])+(nbsa[i][1]&nbsa[i][2])+(nbsa[i][2]&nbsa[i][3]));
	nbsa[i][5] +=(values->transition*((nbsa[i][0]&nbsa[i][2])+(nbsa[i][1]&nbsa[i][3])));
	if ((nbsa[i][2]+nbsa[i][1]+nbsa[i][0]+nbsa[i][3])==3) {
		nbsa[i][5] <<= 1;
		nbsa[i][5]/=3;
		}
	else if ((nbsa[i][2]+nbsa[i][1]+nbsa[i][0]+nbsa[i][3])==4) nbsa[i][5] >>= 1;
	}
else {
	nbsa[i][5]+=(values->delta[0][1]*(nbsa[i][0]&nbsa[i][1]));
	nbsa[i][5]+=(values->delta[0][2]*(nbsa[i][0]&nbsa[i][2]));
	nbsa[i][5]+=(values->delta[0][3]*(nbsa[i][0]&nbsa[i][3]));
	nbsa[i][5]+=(values->delta[1][2]*(nbsa[i][1]&nbsa[i][2]));
	nbsa[i][5]+=(values->delta[1][3]*(nbsa[i][1]&nbsa[i][3]));
	nbsa[i][5]+=(values->delta[2][3]*(nbsa[i][2]&nbsa[i][3]));
	if ((nbsa[i][2]+nbsa[i][1]+nbsa[i][0]+nbsa[i][3])==3) {
		nbsa[i][5] <<= 1;
		nbsa[i][5]/=3;
		}
	else	if ((nbsa[i][2]+nbsa[i][1]+nbsa[i][0]+nbsa[i][3])==4) nbsa[i][5] >>= 1;
	}
	}
else {
	/*for new optimization really for a single sequence but makes no assumption*/
	nbsa=(int **)malloc(a->n_seqs*sizeof(int *));
	assert((int)nbsa);
	for (i=0;i<a->n_seqs;i++) {
		nbsa[i]=(int *)malloc(a->length*sizeof(int));
		assert((int)nbsa[i]);
		}
	/*get bit reps*/
	for (i=0;i<a->n_seqs;i++) for (j=0;j<a->length;j++)  nbsa[i][j]=charbit[a->s[i][j]];
	}
return nbsa;
}

int **do_bit_thing_string(a,length,values)
char *a;
int length;
parameters *values;
{
int **nbsa;
int i,j;

if (!values->new_optimization) {fprintf(stderr,"Wrong bit thing!"); exit(1);}

/*for new optimization really for a single sequence but makes no assumption*/
nbsa=(int **)malloc(sizeof(int *));
assert((int)nbsa);
nbsa[0]=(int *)malloc(length*sizeof(int));
assert((int)nbsa[0]);

/*get bit reps*/
for (j=0;j<length;j++)  nbsa[0][j]=charbit[a[j]];

return nbsa;
}

