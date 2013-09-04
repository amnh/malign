/*Copyright 1992 Ward Wheeler all rights reserved*/

#include "align3.h"

alignment **alignments;
int *present;
int *present2;
int *from_i;
int *from_j;
int alignments_extent;

int find_alignment (i, j, values)
int i, j;
parameters *values;
{
	char *name;
	int k;

	/*if (values->VERBOSE) fprintf (stderr, "Seek alignment %d %d\n", i, j);*/
	if (strcmp (alignments[i]->name, alignments[j]->name) > 0) {
		name = pair_names (alignments[j], alignments[i]);
		}
	else {
		name = pair_names (alignments[i], alignments[j]);
		}
	for (k = 0; k < alignments_extent; k++) if (present2[k]) {
		if (strcmp (alignments[k]->name, name) == 0) {
			free (name);
			return k;
		}
	}
	free (name);
	if (values->VERBOSE) fprintf (stderr, "Perform alignment %d %d --> %d\n",
	i, j, alignments_extent);
	alignments[alignments_extent] = nw (alignments[i], alignments[j], values);
	present[alignments_extent] = 0;
	return alignments_extent++;
}

alignment *minimum_spanning_alignment (a, n, values)
alignment **a;
parameters *values;
int n;
{
	int i, j, k, min_cost, min_i, min_j, min_k;
	alignment *b;

	/*fprintf(stderr,"mem=%lu\n",(unsigned long) coreleft());*/
	assert (n > 1);
	alignments = (alignment **)allocate (n * n * sizeof (alignment *));
	present = (int *)allocate (n * n * sizeof (int));
	present2 = (int *)allocate (n * n * sizeof (int));
	from_i = (int *)allocate (n * n * sizeof (int));
	from_j = (int *)allocate (n * n * sizeof (int));
	assert((int)present2);
	for (i = 0; i < n; i++) {
		alignments[i] = make_align(a[i]);
		present[i] = 1;
		}

	for (i=0;i<(n*n);i++) {
		present2[i]=1;
		from_i[i]=-1;
		from_j[i]=-1;
		}
	alignments_extent = n;
	while (n > 1) {
		min_cost = HUGE_COST;
		for (i = 0; i < alignments_extent; i++) if (present[i]) {
			for (j = i + 1; j < alignments_extent; j++) if (present[j]) {
				k = find_alignment (i, j, values);
				if (alignments[k]->score < min_cost) {
				min_cost = alignments[k]->score;
				min_i = i;
				min_j = j;
				min_k = k;
				}
			}
		}
		if (values->VERBOSE) fprintf (stderr, "Select alignment %d, eliminate %d and %d\n",
				min_k, min_i, min_j);

		from_i[min_k]=min_i;
		from_j[min_k]=min_j;

		alignments[min_i]=dump_align(alignments[min_i]);
		alignments[min_j]=dump_align(alignments[min_j]);
		
		present2[min_i]=present2[min_j]=0;

		for (j=0;j<alignments_extent;j++) {
			if ((from_i[j]==min_i) && (present2[min_i])) {
				alignments[j]=dump_align(alignments[j]);
				present2[j]=0;
				}
			else if ((from_i[j]==min_j) && (present2[min_j])) {
				alignments[j]=dump_align(alignments[j]);
				present2[j]=0;
				}
			if ((from_j[j]==min_i) && (present2[min_i])) {
				alignments[j]=dump_align(alignments[j]);
				present2[j]=0;
				}
			else if ((from_j[j]==min_j) && (present2[min_j])) {
				alignments[j]=dump_align(alignments[j]);
				present2[j]=0;
				}
			}
		present[min_i] = present[min_j] = 0;
		present[min_k] = 1;
		n--;
	}
	free (present);

	/*fix this so it makes sense */
	b=make_align(alignments[min_k]);
	for (i=0;i<alignments_extent;i++) if (present2[i]) alignments[i]=dump_align(alignments[i]);
	free (alignments);
	/*changed to make free*/
	free (present2);
	free (from_i);
	free (from_j);

	return b;
}


