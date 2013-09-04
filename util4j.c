/*Copyright 1994 Ward Wheeler all rights reserved*/
/*Make sure that the reconstructed alignments == the initial ones + can count gaps and changes from down pass as a check*/

#include "align3.h"

/* Need to take two descendent alignments which are basically equall and reconcile them with ancestor alignment
in each case 0 is desc 1 is ancestor
a1 is current node and ancestor
a2 is desc1 and cur
a3 is desc2 and cur
*/


int get_best_state(cur,anc,desc0,desc1,rand,values)
int cur, anc, desc0, desc1, rand;
parameters *values;
{
int best_state, holder;
int cost_A,cost_C,cost_G,cost_T,cost_GAP;
int best_cost,num_best,choice;
int could_be_A,could_be_C,could_be_G,could_be_T,could_be_GAP;

could_be_A=could_be_C=could_be_G=could_be_T=could_be_GAP=0;
if (cur & A_BIT) could_be_A=1;
if (cur & C_BIT) could_be_C=1;
if (cur & G_BIT) could_be_G=1;
if (cur & T_BIT) could_be_T=1;
if (cur & GAP_BIT) could_be_GAP=1;

/*exaustive with this
could_be_A=could_be_C=could_be_G=could_be_T=could_be_GAP=1;*/

cost_A=cost_C=cost_G=cost_T=cost_GAP=HUGE_COST;
if (could_be_A) {
	cost_A=0;
	if (!(A_BIT&anc)) {
		if ((anc&C_BIT) || (anc&G_BIT) || (anc&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	if (!(A_BIT&desc0)) {
		if ((desc0&C_BIT) || (desc0&G_BIT) || (desc0&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	if (!(A_BIT&desc1)) {
		if ((desc1&C_BIT) || (desc1&G_BIT) || (desc1&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	}
if (could_be_C) {
	cost_C=0;
	if (!(C_BIT&anc)) {
		if ((anc&A_BIT) || (anc&G_BIT) || (anc&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	if (!(C_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&G_BIT) || (desc0&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	if (!(C_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&G_BIT) || (desc1&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	}
if (could_be_G) {
	cost_G=0;
	if (!(G_BIT&anc)) {
		if ((anc&A_BIT) || (anc&C_BIT) || (anc&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	if (!(G_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&C_BIT) || (desc0&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	if (!(G_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&C_BIT) || (desc1&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	}
if (could_be_T) {
	cost_T=0;
	if (!(T_BIT&anc)) {
		if ((anc&A_BIT) || (anc&C_BIT) || (anc&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	if (!(T_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&C_BIT) || (desc0&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	if (!(T_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&C_BIT) || (desc1&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	}
if (could_be_GAP) {
	cost_GAP=0;
	if (!(GAP_BIT&anc)) cost_GAP+=values->gap_cost;
	if (!(GAP_BIT&desc0)) cost_GAP+=values->gap_cost;
	if (!(GAP_BIT&desc1)) cost_GAP+=values->gap_cost;
	}
if (!rand) {
	best_cost=HUGE_COST;
	if (cost_A<best_cost) {
		best_cost=cost_A;
		best_state=A_BIT;
		}
	if (cost_C<best_cost) {
		best_cost=cost_C;
		best_state=C_BIT;
		}
	if (cost_G<best_cost) {
		best_cost=cost_G;
		best_state=G_BIT;
		}
	if (cost_T<best_cost) {
		best_cost=cost_T;
		best_state=T_BIT;
		}
	if (cost_GAP<best_cost) {
		best_cost=cost_GAP;
		best_state=GAP_BIT;
		}
	}
else {
	best_cost=HUGE_COST;
	if (cost_A<best_cost) {
		best_cost=cost_A;
		best_state=A_BIT;
		}
	if (cost_C<best_cost) {
		best_cost=cost_C;
		best_state=C_BIT;
		}
	if (cost_G<best_cost) {
		best_cost=cost_G;
		best_state=G_BIT;
		}
	if (cost_T<best_cost) {
		best_cost=cost_T;
		best_state=T_BIT;
		}
	if (cost_GAP<best_cost) {
		best_cost=cost_GAP;
		best_state=GAP_BIT;
		}
	num_best=0;
	if (cost_A==best_cost) ++num_best;
	if (cost_C==best_cost) ++num_best;
	if (cost_G==best_cost) ++num_best;
	if (cost_T==best_cost) ++num_best;
	if (cost_GAP==best_cost) ++num_best;
	choice=my_randomize(num_best);
	if ((cost_A==best_cost) && (choice==0)) best_state=A_BIT;
	if ((cost_C==best_cost) && (choice==1)) best_state=C_BIT;
	if ((cost_G==best_cost) && (choice==2)) best_state=G_BIT;
	if ((cost_T==best_cost) && (choice==3)) best_state=T_BIT;
	if ((cost_GAP==best_cost) && (choice==4)) best_state=GAP_BIT;
	}

if ((best_state != A_BIT) && (best_state != C_BIT) && (best_state != G_BIT) && (best_state != T_BIT) && (best_state != GAP_BIT)) {
	fprintf(stderr,"Error in optimization\n");
	exit(-1);
	}
/*fprintf(stderr,"(%d %d %d %d)->%d ",anc,desc0,desc1,cur,best_state);*/
return best_state;
}


int get_all_states(cur,anc,desc0,desc1,rand,values)
int cur, anc, desc0, desc1, rand;
parameters *values;
{
int best_state, holder;
int cost_A,cost_C,cost_G,cost_T,cost_GAP;
int best_cost,num_best,choice;
int could_be_A,could_be_C,could_be_G,could_be_T,could_be_GAP;

/*exaustive with this*/
could_be_A=could_be_C=could_be_G=could_be_T=could_be_GAP=1;

cost_A=cost_C=cost_G=cost_T=cost_GAP=HUGE_COST;
if (could_be_A) {
	cost_A=0;
	if (!(A_BIT&anc)) {
		if ((anc&C_BIT) || (anc&G_BIT) || (anc&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	if (!(A_BIT&desc0)) {
		if ((desc0&C_BIT) || (desc0&G_BIT) || (desc0&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	if (!(A_BIT&desc1)) {
		if ((desc1&C_BIT) || (desc1&G_BIT) || (desc1&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	}
if (could_be_C) {
	cost_C=0;
	if (!(C_BIT&anc)) {
		if ((anc&A_BIT) || (anc&G_BIT) || (anc&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	if (!(C_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&G_BIT) || (desc0&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	if (!(C_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&G_BIT) || (desc1&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	}
if (could_be_G) {
	cost_G=0;
	if (!(G_BIT&anc)) {
		if ((anc&A_BIT) || (anc&C_BIT) || (anc&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	if (!(G_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&C_BIT) || (desc0&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	if (!(G_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&C_BIT) || (desc1&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	}
if (could_be_T) {
	cost_T=0;
	if (!(T_BIT&anc)) {
		if ((anc&A_BIT) || (anc&C_BIT) || (anc&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	if (!(T_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&C_BIT) || (desc0&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	if (!(T_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&C_BIT) || (desc1&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	}
if (could_be_GAP) {
	cost_GAP=0;
	if (!(GAP_BIT&anc)) cost_GAP+=values->gap_cost;
	if (!(GAP_BIT&desc0)) cost_GAP+=values->gap_cost;
	if (!(GAP_BIT&desc1)) cost_GAP+=values->gap_cost;
	}
if (!rand) {
	best_cost=HUGE_COST;
	if (cost_A<best_cost) {
		best_cost=cost_A;
		best_state=A_BIT;
		}
	if (cost_C<best_cost) {
		best_cost=cost_C;
		best_state=C_BIT;
		}
	if (cost_G<best_cost) {
		best_cost=cost_G;
		best_state=G_BIT;
		}
	if (cost_T<best_cost) {
		best_cost=cost_T;
		best_state=T_BIT;
		}
	if (cost_GAP<best_cost) {
		best_cost=cost_GAP;
		best_state=GAP_BIT;
		}
	}
else {
	best_cost=HUGE_COST;
	if (cost_A<best_cost) {
		best_cost=cost_A;
		best_state=A_BIT;
		}
	if (cost_C<best_cost) {
		best_cost=cost_C;
		best_state=C_BIT;
		}
	if (cost_G<best_cost) {
		best_cost=cost_G;
		best_state=G_BIT;
		}
	if (cost_T<best_cost) {
		best_cost=cost_T;
		best_state=T_BIT;
		}
	if (cost_GAP<best_cost) {
		best_cost=cost_GAP;
		best_state=GAP_BIT;
		}
	best_state=0;
	if (cost_A==best_cost) best_state+=A_BIT;
	if (cost_C==best_cost) best_state+=C_BIT;
	if (cost_G==best_cost) best_state+=G_BIT;
	if (cost_T==best_cost) best_state+=T_BIT;
	if (cost_GAP==best_cost) best_state+=GAP_BIT;
	}

return best_state;
}

/*
char get_best_state_root(in_cur,in_desc0,in_desc1,rand,values)
char in_cur,in_desc0,in_desc1;
int rand;
parameters *values;
{
int cur, desc0, desc1, rand;
int best_state, holder;
int cost_A,cost_C,cost_G,cost_T,cost_GAP;
int best_cost,num_best,choice;
int could_be_A,could_be_C,could_be_G,could_be_T,could_be_GAP;

could_be_A=could_be_C=could_be_G=could_be_T=could_be_GAP=0;
if (cur & A_BIT) could_be_A=1;
if (cur & C_BIT) could_be_C=1;
if (cur & G_BIT) could_be_G=1;
if (cur & T_BIT) could_be_T=1;
if (cur & GAP_BIT) could_be_GAP=1;

could_be_A=could_be_C=could_be_G=could_be_T=could_be_GAP=1;

cost_A=cost_C=cost_G=cost_T=cost_GAP=HUGE_COST;
if (could_be_A) {
	cost_A=0;
	if (!(A_BIT&anc)) {
		if ((anc&C_BIT) || (anc&G_BIT) || (anc&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	if (!(A_BIT&desc0)) {
		if ((desc0&C_BIT) || (desc0&G_BIT) || (desc0&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	if (!(A_BIT&desc1)) {
		if ((desc1&C_BIT) || (desc1&G_BIT) || (desc1&T_BIT)) cost_A+=values->change_cost;
		else cost_A+=values->gap_cost;
		}
	}
if (could_be_C) {
	cost_C=0;
	if (!(C_BIT&anc)) {
		if ((anc&A_BIT) || (anc&G_BIT) || (anc&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	if (!(C_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&G_BIT) || (desc0&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	if (!(C_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&G_BIT) || (desc1&T_BIT)) cost_C+=values->change_cost;
		else cost_C+=values->gap_cost;
		}
	}
if (could_be_G) {
	cost_G=0;
	if (!(G_BIT&anc)) {
		if ((anc&A_BIT) || (anc&C_BIT) || (anc&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	if (!(G_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&C_BIT) || (desc0&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	if (!(G_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&C_BIT) || (desc1&T_BIT)) cost_G+=values->change_cost;
		else cost_G+=values->gap_cost;
		}
	}
if (could_be_T) {
	cost_T=0;
	if (!(T_BIT&anc)) {
		if ((anc&A_BIT) || (anc&C_BIT) || (anc&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	if (!(T_BIT&desc0)) {
		if ((desc0&A_BIT) || (desc0&C_BIT) || (desc0&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	if (!(T_BIT&desc1)) {
		if ((desc1&A_BIT) || (desc1&C_BIT) || (desc1&G_BIT)) cost_T+=values->change_cost;
		else cost_T+=values->gap_cost;
		}
	}
if (could_be_GAP) {
	cost_GAP=0;
	if (!(GAP_BIT&anc)) cost_GAP+=values->gap_cost;
	if (!(GAP_BIT&desc0)) cost_GAP+=values->gap_cost;
	if (!(GAP_BIT&desc1)) cost_GAP+=values->gap_cost;
	}
if (!rand) {
	best_cost=HUGE_COST;
	if (cost_A<best_cost) {
		best_cost=cost_A;
		best_state=A_BIT;
		}
	if (cost_C<best_cost) {
		best_cost=cost_C;
		best_state=C_BIT;
		}
	if (cost_G<best_cost) {
		best_cost=cost_G;
		best_state=G_BIT;
		}
	if (cost_T<best_cost) {
		best_cost=cost_T;
		best_state=T_BIT;
		}
	if (cost_GAP<best_cost) {
		best_cost=cost_GAP;
		best_state=GAP_BIT;
		}
	}
else {
	best_cost=HUGE_COST;
	if (cost_A<best_cost) {
		best_cost=cost_A;
		best_state=A_BIT;
		}
	if (cost_C<best_cost) {
		best_cost=cost_C;
		best_state=C_BIT;
		}
	if (cost_G<best_cost) {
		best_cost=cost_G;
		best_state=G_BIT;
		}
	if (cost_T<best_cost) {
		best_cost=cost_T;
		best_state=T_BIT;
		}
	if (cost_GAP<best_cost) {
		best_cost=cost_GAP;
		best_state=GAP_BIT;
		}
	num_best=0;
	if (cost_A==best_cost) ++num_best;
	if (cost_C==best_cost) ++num_best;
	if (cost_G==best_cost) ++num_best;
	if (cost_T==best_cost) ++num_best;
	if (cost_GAP==best_cost) ++num_best;
	choice=my_randomize(num_best);
	if ((cost_A==best_cost) && (choice==0)) best_state=A_BIT;
	if ((cost_C==best_cost) && (choice==1)) best_state=C_BIT;
	if ((cost_G==best_cost) && (choice==2)) best_state=G_BIT;
	if ((cost_T==best_cost) && (choice==3)) best_state=T_BIT;
	if ((cost_GAP==best_cost) && (choice==4)) best_state=GAP_BIT;
	}

if ((best_state != A_BIT) && (best_state != C_BIT) && (best_state != G_BIT) && (best_state != T_BIT) && (best_state != GAP_BIT)) {
	fprintf(stderr,"Error in optimization\n");
	exit(-1);
	}
fprintf(stderr,"(%d %d %d %d)->%d ",anc,desc0,desc1,cur,best_state);
return best_state;
}
*/