/*Copyright 1992 Ward Wheeler all rights reserved*/

#include "align3.h"

void *allocate (u)
unsigned u;
{
	void *p;

	p = (void *)malloc (u);
	assert ((int) p);
	return p;
}

alignment *allocate_alignment (n)
int n;
{
	alignment *a;

	a = (alignment *) malloc (1*sizeof (alignment));
	assert((int)a);
	a->s = (char **) malloc (n * sizeof (char *));
	assert((int) a->s);
	a->n_seqs = n;
	a->taxon_name = (char**) malloc(n*sizeof(char *));
	assert((int)a->taxon_name);
	a->type_weight=0;
	return a;
}

alignment *allocate_alignment_and_strings (n, l)
int n, l;
{
	alignment *a;
	int i;

	a = allocate_alignment (n);
	for (i = 0; i < n; i++) {
		a->s[i] = (char *) malloc (l);
	}
	return a;
}

alignment *string_to_alignment (s, name)
char *s, *name;
{
	alignment *a;

	a = allocate_alignment (1);
	a->length = strlen (s);
	a->s[0] = s;
	a->name = name;
	a->score = 0;
	return a;
}

void skipspace (s)
char **s;
{
	while (**s && isspace (**s)) (*s)++;
}

void copyword (s, bufp)
char **s, **bufp;
{
	while (**s && !isspace (**s)) *(*bufp)++ = *(*s)++;
}

char linebuf[LINESIZE];

alignment *read_alignment (values)
parameters *values;
{
	char *s, *bufp, *oldbufp;
	alignment *a;
	int w,i;
	char *buf,is_protein=0;

	buf=(char *)malloc(MAX_SEQUENCE_SIZE*sizeof(char));
	assert((int)buf);
	do {
		s = gets (linebuf);
		if (!s) return (alignment *) NULL;
		/* if (VERBOSE) puts (s);*/
	}	while (strlen (s) == 0);
	skipspace (&s);
	bufp = buf;
	copyword (&s, &bufp);
	*bufp = '\0';
	a = allocate_alignment (1);
	a->name = (char *)allocate (strlen (buf) + 1);
	strcpy (a->name, buf);
	bufp = buf;
	s = gets (linebuf);
	/*if (VERBOSE) puts (s);*/
	w = 0;
	while (s && strlen (s) > 0) {
		oldbufp = bufp;
		while (*s) {
			skipspace (&s);
			copyword (&s, &bufp);
			if (w++ == 0) bufp = oldbufp;
		}
		s = gets (linebuf);
		w = 0;
		/*	if (VERBOSE) puts (s);*/
	}
	*bufp = '\0';
	a->s[0] = (char *)allocate (strlen (buf) + 1);
	strcpy (a->s[0], buf);
	a->length = strlen (buf);
	free(buf);
	a->score = 0;
	a->taxon_name[0]=NULL;  
	return a;
}

alignment *reverse_complement_and_rename(a,values)
alignment *a;
parameters *values;
{
int i;
char *temp;

/*redo name*/
a->name[strlen(a->name)-1]='\0';


if (values->VERBOSE) fprintf(stderr,"Reversing and complementing %s\n",a->name);

temp=(char *)malloc(a->length*sizeof(char));
assert((int)temp);
for (i=0;i<a->length;i++) temp[i]=a->s[0][i];

for (i=0;i<a->length;i++) {
	switch (temp[i]) {
		case 'A' :
    	a->s[0][a->length-1-i]='T';
    	break;  	
		case 'C' :
    	a->s[0][a->length-1-i]='G';
    	break;  	
		case 'G' :
    	a->s[0][a->length-1-i]='C';
    	break;  	
		case 'T' :
    	a->s[0][a->length-1-i]='A';
    	break;  	
		case 'U' :
    	a->s[0][a->length-1-i]='A';
    	break;  	
		case 'R' :
    	a->s[0][a->length-1-i]='Y';
    	break;  	
		case 'Y' :
    	a->s[0][a->length-1-i]='R';
    	break;  	
		case 'S' :
    	a->s[0][a->length-1-i]='W';
    	break;  	
		case 'W' :
    	a->s[0][a->length-1-i]='S';
    	break;  	
		case 'M' :
    	a->s[0][a->length-1-i]='K';
    	break;  	
		case 'K' :
    	a->s[0][a->length-1-i]='M';
    	break;  	
		case 'X' :
    	a->s[0][a->length-1-i]='X';
    	break;  	
		case 'N' :
    	a->s[0][a->length-1-i]='N';
    	break;  	
		case 'B' :
    	a->s[0][a->length-1-i]='V';
    	break;  	
		case 'D' :
    	a->s[0][a->length-1-i]='H';
    	break;  	
		case 'H' :
    	a->s[0][a->length-1-i]='D';
    	break;  	
		case 'V' :
    	a->s[0][a->length-1-i]='B';
    	break;  	
		case '-' :
    	a->s[0][a->length-1-i]='-';
    	break;  	
    default:
    	fprintf(stderr,"Unrecognized nucleotide %c in %s!\n",temp[i],a->name);
      exit(-1);
  	}	
	}

free(temp);
return a;
}

alignment *reverse_and_rename(a,values)
alignment *a;
parameters *values;
{
int i;
char *temp;

/*redo name*/


a->name[strlen(a->name)-1]='\0';
if (values->VERBOSE) fprintf(stderr,"Reversing %s\n",a->name);

temp=(char *)malloc(a->length*sizeof(char));
assert((int)temp);
for (i=0;i<a->length;i++) temp[i]=a->s[0][i];

for (i=0;i<a->length;i++) a->s[0][a->length-1-i]=temp[i];

free(temp);
return a;

}
alignment *complement_and_rename(a,values)
alignment *a;
parameters *values;
{
int i;
char *temp;

/*redo name*/
a->name[strlen(a->name)-1]='\0';

if (values->VERBOSE) fprintf(stderr,"Complementing %s\n",a->name);

temp=(char *)malloc(a->length*sizeof(char));
assert((int)temp);
for (i=0;i<a->length;i++) temp[i]=a->s[0][i];

for (i=0;i<a->length;i++) {
	switch (temp[i]) {
		case 'A' :
    	a->s[0][i]='T';
    	break;  	
		case 'C' :
    	a->s[0][i]='G';
    	break;  	
		case 'G' :
    	a->s[0][i]='C';
    	break;  	
		case 'T' :
    	a->s[0][i]='A';
    	break;  	
		case 'U' :
    	a->s[0][i]='A';
    	break;  	
		case 'R' :
    	a->s[0][i]='Y';
    	break;  	
		case 'Y' :
    	a->s[0][i]='R';
    	break;  	
		case 'S' :
    	a->s[0][i]='W';
    	break;  	
		case 'W' :
    	a->s[0][i]='S';
    	break;  	
		case 'M' :
    	a->s[0][i]='K';
    	break;  	
		case 'K' :
    	a->s[0][i]='M';
    	break;  	
		case 'X' :
    	a->s[0][i]='X';
    	break;  	
		case 'N' :
    	a->s[0][i]='N';
    	break;  	
		case 'B' :
    	a->s[0][i]='V';
    	break;  	
		case 'D' :
    	a->s[0][i]='H';
    	break;  	
		case 'H' :
    	a->s[0][i]='D';
    	break;  	
		case 'V' :
    	a->s[0][i]='B';
    	break;  	
		case '-' :
    	a->s[0][i]='-';
    	break;  	
    default:
    	fprintf(stderr,"Unrecognized nucleotide %c in %s!\n",temp[i],a->name);
      exit(-1);
  	}	
	}

free(temp);
return a;
}

void protein_to_ambiguity(a,values)
alignment *a;
parameters *values;
{
char *buf;
int i,holder,j,k,temp;

buf=(char *)malloc(((3*a->length)+1)*sizeof(char));
assert((int)buf);
buf[(3*a->length)]='\0';
for (i=0;i<a->length;i++) {
	holder=i*3;
	temp=0;
	switch (a->s[0][i]) {
		case 'A' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='A') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='G';
					buf[holder+1]='C';
					buf[holder+2]='N';
					}
				break;
		case 'R' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='R') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='M';
					buf[holder+1]='G';
					buf[holder+2]='N';
					}
				break;
		case 'N' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='N') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='A';
					buf[holder+1]='A';
					buf[holder+2]='Y';
					}
				break;
		case 'D' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='D') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='G';
					buf[holder+1]='A';
					buf[holder+2]='Y';
					}
				break;
		case 'C' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='C') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='T';
					buf[holder+1]='G';
					buf[holder+2]='Y';
					}
				break;
		case 'Q' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='Q') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='C';
					buf[holder+1]='A';
					buf[holder+2]='R';
					}
				break;
		case 'E' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='E') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='G';
					buf[holder+1]='A';
					buf[holder+2]='R';
					}
				break;
		case 'G' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='G') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='G';
					buf[holder+1]='G';
					buf[holder+2]='N';
					}
				break;
		case 'H' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='H') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='C';
					buf[holder+1]='A';
					buf[holder+2]='Y';
					}
				break;
		case 'I' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='I') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='A';
					buf[holder+1]='T';
					buf[holder+2]='H';
					}
				break;
		case 'L' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='L') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='Y';
					buf[holder+1]='T';
					buf[holder+2]='N';
					}
				break;
		case 'K' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='K') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='A';
					buf[holder+1]='A';
					buf[holder+2]='R';
					}
				break;
		case 'M' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='M') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='A';
					buf[holder+1]='T';
					buf[holder+2]='G';
					}
				break;
		case 'F' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='F') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='T';
					buf[holder+1]='T';
					buf[holder+2]='Y';
					}
				break;
		case 'P' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='P') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='C';
					buf[holder+1]='C';
					buf[holder+2]='N';
					}
				break;
		case 'S' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='S') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='W';
					buf[holder+1]='S';
					buf[holder+2]='N';
					}
				break;
		case 'T' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='T') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='A';
					buf[holder+1]='C';
					buf[holder+2]='N';
					}
				break;
		case 'W' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='W') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='T';
					buf[holder+1]='G';
					buf[holder+2]='G';
					}
				break;
		case 'Y' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='Y') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='T';
					buf[holder+1]='A';
					buf[holder+2]='Y';
					}
				break;
		case 'V' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='V') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='G';
					buf[holder+1]='T';
					buf[holder+2]='N';
					}
				break;
		case 'X' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='X') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '-' : 
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='-') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='-';
					buf[holder+1]='-';
					buf[holder+2]='-';
          }
				break;
		case '0' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='0') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '1' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='1') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '2' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='2') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '3' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='3') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '4' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='4') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '5' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='5') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '6' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='6') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '7' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='7') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '8' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='8') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		case '9' :
				for (j=0;j<values->n_codes;j++) {
					if (values->new_codes[j][0]=='9') {
						for (k=1;k<4;k++) buf[holder+k-1]=values->new_codes[j][k];
						temp=1;
						}
					}
				if (!temp) {
					buf[holder]	='N';
					buf[holder+1]='N';
					buf[holder+2]='N';
					}
				break;
		default:
			fprintf(stderr,"Unrecognized amino acid code in sequence %s at position %d: %c\n",a->name,i+1,a->s[0][i]);
			exit(-1);
		}
	}
free(a->s[0]);
a->s[0]=(char *)malloc((strlen(buf)+1)*sizeof(char));
assert((int)a->s[0]);
strcpy(a->s[0],buf);
a->length=strlen(buf);
free(buf);
}

char *pair_names (a1, a2)
alignment *a1, *a2;
{
	char *p;

	p = (char *) allocate (strlen (a1->name) + strlen (a2->name) +5);
	sprintf (p, "(%s %s)", a1->name, a2->name);
	return p;
}

cell **allocate_matrix (rows, columns)
int rows, columns;
{
	cell **r;
	int i;

	r = (cell **)allocate (rows * sizeof (cell *));
	for (i = 0; i < rows; i++) {
		r[i] = (cell *)allocate (columns * sizeof (cell));
	}
	return r;
}

void free_matrix (r, rows)
cell **r;
int rows;
{
	int i;

	for (i = 0; i < rows; i++) {
		free (r[i]);
		}
}

/*
int min (i, j)
int i, j;
{
	return j < i? j: i;
}
*/

char compare_groups(g1,g2,ngroups1,ngroups2,ntax,nent)
int **g1,**g2;
int ngroups1,ngroups2,ntax,nent;
{
int i,j,k,holder,temp;

/*
fprintf(stderr,"ntax=%d nent=%d\n",ntax,nent);
for (i=0;i<ngroups1;i++) {
	for (j=0;j<=ntax;j++) fprintf(stderr,"%d ",(*(g1[i]+j)));
	fprintf(stderr,"\n");}

for (i=0;i<ngroups2;i++) {
	for (j=0;j<=ntax;j++) fprintf(stderr,"%d ",(*(g2[i]+j)));
	fprintf(stderr,"\n");}
fprintf(stderr,"\n");
*/

/*when taxon not in g2 but in g1 increment*/
for (k=0;k<ngroups1;k++) {
	for (i=0;i<ngroups2;i++) {
		temp=holder=0;
		for (j=0;j<nent;j++) holder+=(*(g1[k]+j)*(*(g2[i]+j)));
		/*when nent<ntax just match after nent*/
		for (j=nent;j<ntax;j++) temp+=(*(g1[k]+j));
		if (holder!=0) {
			if (((*(g1[k]+ntax))-temp)<=((*(g2[i]+ntax)))) {
				if (holder!=((*(g1[k]+ntax))-temp)) return 0;
					}
			else if (holder!=((*(g2[i]+ntax)))) return 0;
			}
		}
	}

return 1;
}


int **get_groups2(file_in,ntax,a,values)
FILE *file_in;
int ntax;
alignment **a;
parameters *values;
{
int i,j,*buffer2,paren_check_left,paren_check_right;
int rep_length,groups_marker,groups_number;
int **groups;
char *buffer;

i=(values->ngroups)=0;

buffer=(char *)malloc(MAX_SEQUENCE_SIZE*sizeof(char));
assert((int)buffer);

i=0;
while (!feof(file_in)) buffer[i++]= (char) fgetc(file_in);
buffer[i-1]='\0';
/*
fprintf(stderr,"%s\n",buffer);
*/
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
free(buffer);
i=0;
while (buffer2[i]>-4) {
	if (buffer2[i]==-2) rep_length=i+1;
	i++;
	}

/*
for (i=0;i<rep_length;i++) fprintf(stderr,"%d,",buffer2[i]);
fprintf(stderr,"\n");
*/

/*the real stuff*/
values->ngroups=paren_check_left;
groups=(int **)malloc(values->ngroups*sizeof(int *));
assert((int)groups);
for (i=0;i<paren_check_left;i++) {
	groups[i]=(int *)malloc((ntax+1)*sizeof(int));
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

for (i=0;i<values->ngroups;i++) {
	for (j=0;j<ntax;j++) groups[i][ntax]+=groups[i][j];
	if (groups[i][ntax]==0) {
		fprintf(stderr,"Error in groups specification -- parenthesis pair which contains no sequences.\n");
		exit(-1);
		}
	/*
	for (j=0;j<=ntax;j++) fprintf(stderr," %d ",groups[i][j]);
	fprintf(stderr,"\n");*/
	}
/*fprintf(stderr,"\n");*/

fclose(file_in);
if (values->VERBOSE) fprintf(stderr,"%d groups required\n",values->ngroups);
free(buffer2);
return groups;
}

void get_taxon_number(i_start,j_start,in_buffer,out_buffer,a,ntax)
alignment **a;
int *i_start, *j_start;
char *in_buffer;
int *out_buffer,ntax;
{
int i,name_found;
char *name_buffer;

name_buffer=(char *)malloc(50*sizeof(char));
assert((int)name_buffer);

i=0;
while (isalnum(in_buffer[(*i_start)+i])) {
	name_buffer[i]=in_buffer[(*i_start)+i];
	++i;
	}
name_buffer[i]='\0';
*i_start+=i;
if (i>50) {
	fprintf(stderr,"Sequence name too long.\n");
	exit(-1);
	}
name_found=0;
for (i=0;i<ntax;i++) {
	if (!stricmp(name_buffer,a[i]->name)) {
		out_buffer[*j_start]=i;
		*j_start+=1;
		if (name_found==0) name_found=1;
		else {
			fprintf(stderr,"Sequence name (%s) repeated in groups representation.\n",a[i]->name);
			exit(-1);
			}
		}
	}
free(name_buffer);
}


alignment *dump_align(d)
alignment *d;
{
int i;

if (d->type_weight>1) d->type_weight=0;
else if (d->type_weight<0) d->type_weight=0;
for (i=0;i<(d->n_seqs+d->type_weight);i++) 	free(d->s[i]);
for (i=0;i<d->n_seqs;i++)	free(d->taxon_name[i]);
free(d->taxon_name);
free(d->s);
free(d->name);
free(d);
d=NULL;
return NULL;
}

