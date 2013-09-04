/*Copyright 1992-1996 Ward Wheeler all rights reserved*/

#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <string.h>

#define malign_version "2.7"
#define PARALLEL_VERSION "1.5"

#define MAC   0
#define IBM32 1
#define IBM16 2
#define SPARC 3
#define HPUX  4
#define LINUX 4
#define GCC_USED 4

#define EXE_FILE_NAME "malign_pvmgcc"
#define BIT_SIZE 48 /*total bits of cost + length + up + over*/
#define CURRENT_PLATFORM GCC_USED
#define PARALLEL 0

#define SPECIAL 0

#if PARALLEL==1
#include "pvm3.h"
#endif

#define THIRTY_TWO_BITS_MAX 4294967295
#define THIRTY_ONE_BITS_MAX 2147483647
#define THIRTY_ONE_BITS_MAX_MINUS_FIFTEEN_BITS 2147451000
#define FIFTEEN_BITS_MAX_MINUS_ONE_THOUSAND 31767
#define FIFTEEN_BITS_MAX 32767
#define THIRTEEN_BITS_MAX 8191
#define SEVENTEEN_BITS_MAX 131071

#define MAX_NUM_SEQUENCES FIFTEEN_BITS_MAX

/*for use in lowmem alignments*/
#if BIT_SIZE==16
#define BARRIER_COST 30000
#endif
#if BIT_SIZE==64
#define HUGE_COST THIRTY_ONE_BITS_MAX     /*15 bits maximum*/
#define MAX_SEQUENCE_SIZE THIRTY_ONE_BITS_MAX      /*15 bits maximum*/
#define COST_BITS 31
#define LENGTH_BITS 31
#define BARRIER_COST THIRTY_ONE_BITS_MAX_MINUS_FIFTEEN_BITS
#endif
#if BIT_SIZE==48
#define HUGE_COST THIRTY_ONE_BITS_MAX    /*17 bits maximum*/
#define MAX_SEQUENCE_SIZE FIFTEEN_BITS_MAX /*13 bits maximum*/
#define COST_BITS 31
#define LENGTH_BITS 15
#define BARRIER_COST THIRTY_ONE_BITS_MAX_MINUS_FIFTEEN_BITS
#endif
#if BIT_SIZE==32
#define HUGE_COST FIFTEEN_BITS_MAX      /*15 bits maximum*/
#define MAX_SEQUENCE_SIZE FIFTEEN_BITS_MAX      /*15 bits maximum*/
#define COST_BITS 15
#define LENGTH_BITS 15
#define BARRIER_COST FIFTEEN_BITS_MAX_MINUS_ONE_THOUSAND
#endif

typedef struct
{
    unsigned union_bit :1;
	unsigned base:5;
	signed cost:26;
} look_up_thang;

#include "bit-vect.h"
#include "yapp.h"
#include "utils.h"
#include "pvm.h"

typedef struct
{
	unsigned direction:2;
	unsigned score:COST_BITS;
	unsigned length:LENGTH_BITS;
} cell;

/*PVM Message types*/
#define PVM_SEND_TIDS 0
#define PVM_ALL_MUST_DIE 1
#define PVM_UNPACK_INITIAL_SHIT 2
#define PVM_DO_ALIGNMENT_NOW 3
#define PVM_PACK_ALIGN_AND_SCORE 4
#define PVM_DO_ALIGNMENT_SCORE_ONLY 5
#define PVM_INITIAL_SHIT_RECIEVED 6
#define PVM_SEND_SCORE_ONLY 7
#define PVM_DUMP_INITIAL_GROUPS 8
#define PVM_REQUEST_TIME 10
#define PVM_SEND_TIME 11
#define PVM_REQUEST_NUM_TREES 12
#define PVM_SEND_NUM_TREES 13
#define PVM_REQUEST_DIAGNOSE_ALIGNMENT 14
#define PVM_SEND_DIAGNOSE_ALIGNMENT 15
#define PVM_DO_ALIGNMENT_NOW_BANDB 16
#define PVM_DO_ALIGNMENT_SCORE_ONLY_BANDB 17
#define PVM_NEW_BUILD_START 18
#define PVM_NEW_BUILD_INIT 19
#define PVM_NEW_BUILD_ADD 20
#define PVM_NEW_BUILD_OPTIMIZE 21
#define PVM_NEW_BUILD_DONE 22
#define PVM_NEW_BUILD_ALIGNMENTS_ONLY_DO 23
#define PVM_NEW_BUILD_ALIGNMENTS_ONLY_SEND 24
#define PVM_REDO_INITIAL_ALIGNS 25
#define PVM_NEW_BUILD_ADD_FUNCTION 26
#define PVM_PARALLEL_NW 27
#define PVM_PARALLEL_NW_DONE 28
#define PVM_NEW_SWAP 29
#define PVM_NEW_SWAP_DONE 30
#define PVM_PARALLEL_NW_NEW_OPT 31
#define PVM_PARALLEL_NW_NEW_OPT_DONE 32
#define PVM_PARALLEL_NW_NEW_OPT_UP_DONE 33
#define PVM_PARALLEL_NW_NEW_OPT_UP_DO 34
#define PVM_PARALLEL_NW_NEW_SWAP_START 35
#define PVM_PARALLEL_NW_NEW_SWAP_DO 36
#define PVM_PARALLEL_NW_NEW_SWAP_STOP 37
#define PVM_SEND_ANODES 38
#define PVM_DAVE_BUILD 39
#define PVM_DAVE_DUMP_NODES 40
#define PVM_DAVE_BUILD_DONE 41
#define PVM_NEW_SWAP_DONE_DAVE 42
#define PVM_NEW_SWAP_DAVE 43
#define PVM_JACK_DO 44
#define PVM_JACK_DONE 45
#define PVM_DAVE_BUILD_GRAIN 46
#define PVM_DAVE_BUILD_GRAIN_DONE 47

#if PARALLEL==0
#include "pvm3.h" /*untill get fake one*/
/*crapm out pvm functions*/
#define pvm_pkint(x,y,z)
#define pvm_pkbyte(x,y,z)
#define pvm_pkstr(x)
#define pvm_upkint(x,y,z)
#define pvm_upkbyte(x,y,z)
#define pvm_upkstr(x)
#define pvm_mcast(x,y,z) 1
#define pvm_nrecv(x,y) 1
#define pvm_recv(x,y) 1
#define pvm_send(x,y);
#define pvm_nsend(x,y)
#define pvm_config(x,y,z) 1
#define pvm_initsend(x)
#define pvm_exit()
#define pvm_mytid() 1
#define pvm_advise(x) 11
#define pvm_parent() 1
#define pvm_bufinfo(x,y,z,e) 1
#define pvm_freebuf(x)
#define pvm_spawn(x,y,z,a,b,c) 1
#endif

#if CURRENT_PLATFORM==HPUX
#define min(a,b) ((a<b) ? a : b)
#define max(a,b) ((a>b) ? a : b)
#define strnicmp strnicmp_mine
#define stricmp stricmp_mine
#endif

/*include for mac*/
#if CURRENT_PLATFORM==MAC
#define GOMAC (void)do_mac_shit();
#else
#define GOMAC /*nothing*/
#endif

/*include these for SUN*/
#if CURRENT_PLATFORM==SPARC
#include <malloc.h>
#define max_rand THIRTY_TWO_BITS_MAX
#else
#include <stdlib.h>
#define max_rand RAND_MAX       /* 2**31 -1 on 32 bit machine RAND_MAX if avail*/
#endif


#include <assert.h>

#if CURRENT_PLATFORM==MAC
int min();
int max();
#define strnicmp strnicmp_mine
#define stricmp stricmp_mine
#endif
#define LINESIZE 1000

#define DIAGONAL 0
#define RIGHT 1
#define DOWN 2
#define RANDOM 3

#define new_all_other_make_nodes all_other_make_nodes2

#define say(x) /*fprintf(stderr,"%s\n",x);*/

#define huge_length MAX_SEQUENCE_SIZE
#define missing 2047 /* early 1's       but doesn't matter as much any more*/
#define little_missing_minus_one 29 /* early 1's       but doesn't matter as much any more*/
#define little_missing 31

#define A_BIT 1
#define C_BIT 2
#define G_BIT 4
#define T_BIT 8
#define U_BIT 8
#define N_BIT 15
#define X_BIT /* 15 */ 31
#define GAP_BIT 16
#define R_BIT 5
#define Y_BIT 10
#define M_BIT 3
#define W_BIT 9
#define S_BIT 6
#define K_BIT 12
#define B_BIT 14
#define D_BIT 13
#define H_BIT 11
#define V_BIT 7

#define A_GAP_BIT 17
#define C_GAP_BIT 18
#define G_GAP_BIT 20
#define T_GAP_BIT 24
#define U_GAP_BIT 24
#define N_GAP_BIT 31
#define X_GAP_BIT 31
#define R_GAP_BIT 21
#define Y_GAP_BIT 26
#define M_GAP_BIT 19
#define W_GAP_BIT 25
#define S_GAP_BIT 22
#define K_GAP_BIT 28
#define B_GAP_BIT 30
#define D_GAP_BIT 29
#define H_GAP_BIT 27
#define V_GAP_BIT 23

#define E_BIT A_GAP_BIT
#define F_BIT C_GAP_BIT
#define I_BIT G_GAP_BIT
#define J_BIT T_GAP_BIT
#define L_BIT U_GAP_BIT
#define O_BIT N_GAP_BIT
#define P_BIT X_GAP_BIT
#define Q_BIT R_GAP_BIT
#define Z_BIT Y_GAP_BIT
#define C1_BIT M_GAP_BIT
#define C2_BIT W_GAP_BIT
#define C3_BIT S_GAP_BIT
#define C4_BIT K_GAP_BIT
#define C5_BIT B_GAP_BIT
#define C6_BIT D_GAP_BIT
#define C7_BIT H_GAP_BIT
#define C8_BIT V_GAP_BIT

#define TRANSITION_BIT 64
#define TRANSVERSION_BIT 128

#define SETME_NEW_OPT       if (!(vbs1&GAP_BIT)) down=values->gap_cost;\
					        else down=0;\
						    if (!(vbs2&GAP_BIT)) right=values->gap_cost;\
						    else right=0;

#define DOME_NEW_OPT        diag=values->lookup[vbs1][vbs2].cost;

#define DOMENOW_NEW_OPT     DOME_NEW_OPT

#define DOMELATER_NEW_OPT   DOME_NEW_OPT


/*add previous cell score*/
#define ADDEM               diag  += pmi1[j-1].score;\
							down  += pmi1[j].score;\
							right += pmi[j-1].score;

/*determine best path*/
#define CHECKUM             minimum = min (min (diag, down), right);\
							pmi[j].score = minimum;\
		                    if (values->min_num_gaps>0) {alter(minimum,diag,down,right,pmi[j-1],pmi1[j],values->min_num_gaps)};\
		                    if (values->pref_direc==RANDOM) {\
			                    direc_holder= (rand()/(max_rand/3));\
		                        if (direc_holder==1) values->pref_direc=DOWN;\
		                        else if (direc_holder>1) values->pref_direc=RIGHT;\
			                    }\
		                    if (values->pref_direc==DIAGONAL) {\
								if (diag == minimum) {\
									pmi[j].direction = DIAGONAL;\
									pmi[j].length = 1 + pmi1[j-1].length;\
									}\
								else if (down == minimum) {\
									pmi[j].direction = DOWN;\
									pmi[j].length = 1 + pmi1[j].length;\
									}\
								else {\
									pmi[j].direction = RIGHT;\
									pmi[j].length = 1 + pmi[j-1].length;\
									}\
		                        }\
		                    else if (values->pref_direc==RIGHT) {\
								if (right == minimum) {\
									pmi[j].direction = RIGHT;\
									pmi[j].length = 1 + pmi[j-1].length;\
									}\
								else if (diag == minimum) {\
									pmi[j].direction = DIAGONAL;\
									pmi[j].length = 1 + pmi1[j-1].length;\
									}\
								else  {\
									pmi[j].direction = DOWN;\
									pmi[j].length = 1 + pmi1[j].length;\
									}\
                    		     }\
                    		  else {\
								if (down == minimum) {\
									pmi[j].direction = DOWN;\
									pmi[j].length = 1 + pmi1[j].length;\
									}\
								else if (right == minimum) {\
									pmi[j].direction = RIGHT;\
									pmi[j].length = 1 + pmi[j-1].length;\
									}\
								else  {\
									pmi[j].direction = DIAGONAL;\
									pmi[j].length = 1 + pmi1[j-1].length;\
									}\
                    		     }\
				            if (rand_marker) values->pref_direc=RANDOM;

#define SETME               diag=(values->gap_cost * (vbs1[4]|vbs2[4]));\
							down=vbs1[5];\
						    right=vbs2[5];\

#define DOME                diag+=(values->change_cost * ((vbs1[2]|vbs2[2])+(vbs1[1]|vbs2[1])+(vbs1[0]|vbs2[0])+(vbs1[3]|vbs2[3])-1));\

#define DOMENOW             diag +=(values->transition+values->transversion)*(((vbs1[0]|vbs2[0])&(vbs1[1]|vbs2[1]))+((vbs1[0]|vbs2[0])&(vbs1[3]|vbs2[3]))+((vbs1[1]|vbs2[1])&(vbs1[2]|vbs2[2]))+((vbs1[2]|vbs2[2])&(vbs1[3]|vbs2[3])));\
							diag +=values->transition*(((vbs1[0]|vbs2[0])&(vbs1[2]|vbs2[2]))+((vbs1[1]|vbs2[1])&(vbs1[3]|vbs2[3])));\
							if (((vbs1[2]|vbs2[2])+(vbs1[1]|vbs2[1])+(vbs1[0]|vbs2[0])+(vbs1[3]|vbs2[3]))==3) {\
								diag <<= 1;\
								diag/=3;\
								}\
							else if (((vbs1[2]|vbs2[2])+(vbs1[1]|vbs2[1])+(vbs1[0]|vbs2[0])+(vbs1[3]|vbs2[3]))==4) diag >>= 1;\

#define DOMELATER           diag+=(values->delta[0][1]*((vbs1[0]|vbs2[0])&(vbs1[1]|vbs2[1])));\
							diag+=(values->delta[0][2]*((vbs1[0]|vbs2[0])&(vbs1[2]|vbs2[2])));\
							diag+=(values->delta[0][3]*((vbs1[0]|vbs2[0])&(vbs1[3]|vbs2[3])));\
							diag+=(values->delta[1][2]*((vbs1[1]|vbs2[1])&(vbs1[2]|vbs2[2])));\
							diag+=(values->delta[1][3]*((vbs1[1]|vbs2[1])&(vbs1[3]|vbs2[3])));\
							diag+=(values->delta[2][3]*((vbs1[2]|vbs2[2])&(vbs1[3]|vbs2[3])));\
							if (((vbs1[2]|vbs2[2])+(vbs1[1]|vbs2[1])+(vbs1[0]|vbs2[0])+(vbs1[3]|vbs2[3]))==3) {\
								diag <<= 1;\
								diag/=3;\
								}\
							else if (((vbs1[2]|vbs2[2])+(vbs1[1]|vbs2[1])+(vbs1[0]|vbs2[0])+(vbs1[3]|vbs2[3]))==4) diag >>= 1;\

#define                     alter(minimum_value,diag,down,right,pmi,pmi1,option)\
					        int num_equal=0;\
					        if (minimum_value==(diag)) ++num_equal;\
					        if (minimum_value==(down)) ++num_equal;\
					        if (minimum_value==(right)) ++num_equal;\
					        if (option==1) {if ((num_equal >1) && ((diag)==minimum_value) && (pmi1.direction!=DIAGONAL)) ++(diag);}\
	                        else if ((num_equal >1) && ((diag)==minimum_value) && (pmi1.direction==DIAGONAL)) {++(down);++(right);}


#define ADDEM_LOWMEM        reg_int_holder=j-i+values->max_gap+1;\
							diag  += pmi1[reg_int_holder].score;\
							down  += pmi1[reg_int_holder+1].score;\
							right += pmi[reg_int_holder-1].score;


#define CHECKUM_LOWMEM      minimum = min (min (diag, down), right);\
							pmi[reg_int_holder].score = minimum;\
						    if (values->min_num_gaps==1) {alter(minimum,diag,down,right,pmi[reg_int_holder-1],pmi1[reg_int_holder+1],values->min_num_gaps)};\
						    if (values->pref_direc==DIAGONAL) {\
							    if (diag == minimum) {\
									pmi[reg_int_holder].direction = DIAGONAL;\
									pmi[reg_int_holder].length = 1 + pmi1[reg_int_holder].length;\
									}\
								else if (down == minimum) {\
									pmi[reg_int_holder].direction = DOWN;\
									pmi[reg_int_holder].length = 1 + pmi1[reg_int_holder+1].length;\
									}\
								else {\
									pmi[reg_int_holder].direction = RIGHT;\
									pmi[reg_int_holder].length = 1 + pmi[reg_int_holder-1].length;\
									}\
                			     }\
				    		  else if (values->pref_direc==DOWN) {\
								if (down == minimum) {\
					    			pmi[reg_int_holder].direction = DOWN;\
									pmi[reg_int_holder].length = 1 + pmi1[reg_int_holder+1].length;\
									}\
								else if (diag == minimum) {\
						    		pmi[reg_int_holder].direction = DIAGONAL;\
									pmi[reg_int_holder].length = 1 + pmi1[reg_int_holder].length;\
									}\
								else {\
									pmi[reg_int_holder].direction = RIGHT;\
									pmi[reg_int_holder].length = 1 + pmi[reg_int_holder-1].length;\
									}\
                			     }\
				    		  else {\
								if (right == minimum) {\
									pmi[reg_int_holder].direction = RIGHT;\
									pmi[reg_int_holder].length = 1 + pmi[reg_int_holder-1].length;\
									}\
								else if (down == minimum) {\
									pmi[reg_int_holder].direction = DOWN;\
									pmi[reg_int_holder].length = 1 + pmi1[reg_int_holder+1].length;\
									}\
								else {\
									pmi[reg_int_holder].direction = DIAGONAL;\
									pmi[reg_int_holder].length = 1 + pmi1[reg_int_holder].length;\
					    			}\
                			     }\


typedef struct
{
	int changes;
	int gaps;
	int transitions;
	int transversions;
	char *name;
} change_counter_structure;

typedef struct
{
	int n_seqs;
	int length;
	int score;
	char **s;
	char *name;
	char **taxon_name;
	int type_weight;
} alignment;

typedef struct
{
	alignment **ba;
	int *age;
	int num_in;
	int max_num;
	int hit;
	int miss;
	int *hack;
} align_buffer;

typedef struct parameters
{
	int get_best;
	int get_heur;
	int gap_cost;
	int change_cost;
	int leading_gap_cost;
	int trailing_gap_cost;
	int phylo_score;
	int all_done;
	int rand_seed;
	int how_many;
	int length_at_end;
	int randtrees;
	int hen_gap;
	int acgt;
	int farris;
	int inter_dig;
	int dot;
	int keep_trees;
	int aquick;
	int tree_swap;
	int align_swap;
	int print_intermediates;
	int keep_aligns;
	int line_length;
	int phylo_gap;
	int VERBOSE;
	int rep_error;
	int paup;
	int paup_dot;
	int **groups;
	int ngroups;
	int in_some;
	int number;
	int best;
	int get_heur2;
	int rand_align;
	char *all_made;
	char *tree_made;
	int time;
	int extra_adjustment;
	int coding;
	int align_multi_swap;
	int tree_multi_swap;
	int n_codes;
	int **delta;
	int length_best_clad;
	int number_best_clad;
	int **best_rep;
	char **new_codes;
	int ttr;
	int transversion;
	int transition;
	int n_transv;
	int givem;
	int matquick;
	int max_gap;
	int iter;
	int low_mem;
	int swap_while_add;
	int clade_swap_while_add;
	char **input_names;
	int output_order;
	int show_mem;
	int max_out_aligns;
	int max_out_trees;
	int cost_only;
	int temp_ti;
	int temp_tv;
	int temp_changes;
	int temp_gaps;
	int overall_ti;
	int overall_tv;
	int overall_changes;
	int overall_gaps;
	int temp_ti_g;
	int temp_tv_g;
	int temp_changes_g;
	int temp_gaps_g;
	int overall_ti_g;
	int overall_tv_g;
	int overall_changes_g;
	int overall_gaps_g;
	int get_weights;
	int get_weights_2;
	int weight_range;
	int cache_size;
	align_buffer **align_cache;
	int cache_info;
	int saw_farris;
	int saw_goloboff;
	int saw_wheeler;
	int input_is_align;
	int **b_string1;
	int **b_string2;
  int phylo_time;
  int pref_direc;
  int min_num_gaps;
  int gap_must_cost;
  int rand_order;
  int reorder;
  int tree_rand_order_max;
  int number_of_input_alignments;
  char **input_file_name;
  char *align_file;
  int actual_num_sequences;
  int cull;
  int elision;
  int *ce_weights;
  alignment *previous;
 int grain_size;
 int *tids;
 int num_hosts;
 int process_factor;
 int new_optimization;
 int no_leading_or_trailing_cost;
 int gap_in_hypanc;
 int *data_sets;
 int *data_set_weights;
 int align_node_swap;
 int align_root_swap;
 int align_complete_root_swap;
 int align_partial_root_swap;
 int dump_parameters;
 int groups_as_start;
 int in_bandb_loop;
 int current_bound;
 int chop;
 int optimize_nodes;
 int in_optimize;
 int apolist;
 int rand_apo;
 int freq_cost;
 int del_acc;
 int num_nws;
 int sbr;
 int tbr;
 int shortcut_tree;
 int shortcut_tree2;
 int get_heur3;
 int asbr;
 int atbr;
 int tbr_align_shortcut;
 int collapse;
 int arrt;
 int **support;
 int max_name_length;
 int do_dave_align;
 struct parameters **other_parm;
 int jackboot;
 char **jack_name;
 int *jack_array;
 int nona;
 int clados;
 int cutpoint;
 int make_groups_file;
 char output_groups_name[100];
 look_up_thang lookup[32][32];
 int *start;
 int *stop;
 int expand_X;
 int best_order;
int poopsy;
int tbr_first;
int worst;
 } parameters;

typedef struct
{
	int **rep;
	char **swap;
} realloc_holder;

typedef struct
{
	int anc_pos;
	int desc_pos;
	char anc_char;
	char desc_char;
} apo_thang;


typedef struct
{
	alignment *z;
	apo_thang *at0;
	apo_thang *at1;
	apo_thang *at2;
} apolist_holder;

typedef struct
{
	int **rep;
	char **swap;
	alignment **best;
} realloc_align_holder;

typedef struct
{
	int *rep;
} b_tree;


void *allocate ();
alignment *allocate_alignment ();
alignment *allocate_alignment_and_strings ();
alignment *string_to_alignment ();
void print_alignment ();
void print_hennig ();
void print_inter_dig ();
void print_hennig_gap_code ();
cell **allocate_matrix ();
void free_matrix ();
void sort_alignment ();
int same_alignment ();
char *pair_names ();
alignment *read_alignment ();
int cladogram();
alignment *nw ();
void get_data();
void loop();
int diagnose_tree();
int get_bound();
realloc_align_holder *all_loop();
int all_make_nodes();
int all_diagnose_tree();
void all_other_make_nodes();
int compare_aligns();
void get_interactive();
void get_file();
alignment **do_some();
alignment **get_random();
alignment **do_all();
int **get_groups();
char compare_groups();
alignment *make_align();
int **get_cur_groups();
void recurse_groups_tree();
void get_cur_groups_tree();
int compare_groups_tree();
void recurse_groups();
alignment **make_heur_align();
int get_quick();
int get_new_quick();
void all_other_make_nodes2();
void printem();
void print_dot();
alignment **get_align_bound();
void shuffle();
void print_paup();
void print_paup_dot();
char check_done();
alignment **reorder();
void get_command_line();
realloc_align_holder *deep_swap_align();
realloc_align_holder *get_tuples_and_loop();
realloc_align_holder *special_all_loop();
int **deep_swap_tree();
realloc_holder *special_all_loop_tree();
realloc_holder *get_tuples_and_loop_tree();
int **get_groups2();
void get_taxon_number();
void get_initial_values();
void fix_all_after_values();
void do_mac_shit();
alignment *dump_align();
/*void dump_align();*/
void protein_to_ambiguity();
char unique();
int how_many();
char already_clad();
void cheap_loop();
void get_guess_max_gap();
void print_values();
alignment *nw_low_mem();
alignment *nw_real();
realloc_align_holder *do_addition_swap();
void clade_addition_swap_1();
int **clade_addition_swap_2();
void order_as_input();
int strnicmp_mine();
alignment *reverse_complement_and_rename();
alignment *complement_and_rename();
alignment *reverse_and_rename();
int diagnose_tree_and_get_weights();
int make_nodes_weights();
alignment *retrieve_align_cache();
void store_align_cache();
void pre_diagnose_align();
void pre_diagnose_tree();
char *other_pair_names();
void mark_descendants();
align_buffer **allocate_align_cache();
int **get_combinatorial_weights();
void get_local_matrix();
char is_internal();
int find_alignment();
alignment *minimum_spanning_alignment();
alignment **read_old_file();
void do_pre_stuff();
void end_free();
void do_estimated_weights();
int **do_bit_thing();
alignment **read_input_align();
void column_score();
void make_nodes();
alignment *strip_gaps();
void do_modifications();
alignment **randomize_taxon_order();
void reorder_tree();
alignment **new_reorder();
alignment **get_prealigned_sequences();
void get_align_file();
alignment **do_cull();
alignment **do_elision();
void new_order_as_input();
void get_robust_positions();
void recode_with_gaps();
int same_position();
int retrieve_score_cache();
void do_initialize_thang();
alignment **unpack_init_vals();
void do_alignment_now();
void do_parallel_shtick();
alignment *make_ambig();
alignment *nw_real_new_opt();
alignment *nw_low_mem_new_opt();
void do_alignment_now();
void do_alignment_score_only();
void do_alignment_now_bandb();
void do_alignment_score_only_bandb();
void pack_align();
alignment *unpack_align();
alignment **new_read_old_file();
alignment *convert_b_to_a();
alignment *unpack_align_and_score();
void strip_out_gaps_and_recount();
alignment **old_get_prealigned_sequences();
alignment *convert_to_alignment();
alignment **new_read_input_align();
alignment **new_get_prealigned_sequences();
realloc_align_holder *do_node_swap();
void all_other_make_nodes3();
void nodes_to_tree_rep();
void pick_node_and_remake_with_from();
void copy_nodes();
void print_nodes();
void split_nodes();
void unroot_pruned();
void block_redundant_additions();
void pick_new_root();
void remake_nodes();
void block_redundant_additions_and_non_pruned ();
void print_nodes_stdout();
void print_nodes3_stdout();
void print_parameters();
alignment *optimization_alignment();
void redo_for_type();
void get_new_topo();
alignment **do_chop_thing();
alignment **get_chopped();
alignment **new_get_prealigned_sequences_chop();
alignment *glue_back_after_chop();
void make_all_unordered();
int got_dash_or_dot();
void do_and_print_optimization();
int **get_groups3();
int all_diagnose_tree_w_uppass();
/*alignment *make_anc_up();*/
apolist_holder *make_anc_up();
apolist_holder *make_anc_up2();
char choose_one2();
int choose_one();
int is_ambig();
alignment *make_max();
int real_position();
char return_character();
int my_randomize();
int **get_correspondances();
alignment *fake_nw();
alignment *new_make_ambig();
char overlap();
alignment *make_neighbor();
int get_best_state();
char get_best_state_root();
void get_reopt();
alignment **make_heur_align_new();
void get_up_pass();
alignment **do_swap_thang_align();
void get_up_pass2();
void do_new_build_parallel();
void do_new_build_parallel2();
void do_new_parallel_swap();
void redo_it_all();
void get_up_pass_new_opt();
void get_up_pass_new_opt_parallel();
void get_final_a();
void convert_to_anc_and_desc();
void convert_to_final();
char overlap_best();
void get_up_pass3();
alignment *get_single_up();
alignment *old_make_ambig();
alignment *fake_nw_real();
char char_overlap_best();
void modify_corres();
alignment *big_make_neighbor();
alignment *make_neighbor_non_align();
alignment *newer_make_ambig();
int is_unique_tree_align();
int compare_groups_tree_align();
char *get_names_with_collapse();
char *third_pair_names();
int daves_cladogram();
void ProcessCommandLineArguments();
TreeBufferT *Phylogeny();
int MissingDataCode();
void NewReadData();
void SaveNode();
void RestoreNodes();
NodeT *FindRoot();
int Repair();
void ClearNode();
int AddSib();
NodeT *RemoveAsSib();
void LocalRearrangement();
void ClearMarks();
void MarkTree();
void SPR();
NodeT *Reroot();
void TBR();
double UniformRand (void);
int BitRand (void);
void *AllocateNew (unsigned int);
void *Reallocate (void *p, unsigned int);
int streq (char *s1, char *s2);
int round (double x);
char *CopyString (char *s);
void fail (char *s);
int only_missing();
alignment *missing_alignment();
alignment **make_dave_align();
alignment **do_swap_thang_align_dave();
void all_diagnose_tree_here_pieces();
int do_dave_build();
void do_new_parallel_swap_dave();
void fix_all_after_values2();
void pack_other_vals();
void unpack_other_vals();
int get_jack_groups();
int **get_groups4();
void get_lookup();
void get_lookup2();
void make_contig_X();
int no_align_column_score();
alignment *make_three();
void redo_for_type_2();
void recurse_groups2();
char **nw_string();
int **do_bit_thing_string();
char **optimization_alignment_string();
alignment *nw_tot_evid();
char **missing_alignment_string();
void make_contig_X_string();
void redo_for_type_3();
int only_missing_string2();

/*#define DOME_NEW_OPT        if (!(vbs1&vbs2)) diag=values->change_cost;
#define DOME_NEW_OPT        if (!(vbs1&vbs2)) {\
							    if ((vbs1==GAP_BIT) || (vbs2==GAP_BIT)) diag=values->gap_cost;\
    							else diag=values->change_cost;\
	    						}\
#define DOMENOW_NEW_OPT    if (!(vbs1&vbs2)) diag=values->lookup[vbs1][vbs2].cost;
#define DOMENOW_NEW_OPT    if (!(vbs1&vbs2)) {\
							if ((vbs1==GAP_BIT) || (vbs2==GAP_BIT)) diag=values->gap_cost;\
							else {\
								diag=values->transition;\
								if ((vbs1&A_BIT) && (!(vbs2&G_BIT))) diag+=values->transversion;\
								else if ((vbs1&G_BIT) && (!(vbs2&A_BIT))) diag+=values->transversion;\
								else if ((vbs1&C_BIT) && (!(vbs2&T_BIT))) diag+=values->transversion;\
								else if ((vbs1&T_BIT) && (!(vbs2&C_BIT))) diag+=values->transversion;\
								}\
							}\
takes transitions fist then transversions when ambiguities are involved
#define DOMELATER_NEW_OPT if (!(vbs1&vbs2)) {\
						if ((vbs1==GAP_BIT) || (vbs2==GAP_BIT)) diag=values->gap_cost;\
						else {\
							vbs1h=min(vbs1,vbs2);\
							vbs2h=max(vbs1,vbs2);\
							if ((vbs1h&A_BIT) && (vbs2h&G_BIT)) diag=values->delta[0][2];\
							else    if ((vbs1h&C_BIT) && (vbs2h&T_BIT)) diag=values->delta[1][3];\
							else    if ((vbs1h&A_BIT) && (vbs2h&C_BIT)) diag=values->delta[0][1];\
							else    if ((vbs1h&A_BIT) && (vbs2h&T_BIT)) diag=values->delta[0][3];\
							else    if ((vbs1h&C_BIT) && (vbs2h&G_BIT)) diag=values->delta[1][2];\
							else    if ((vbs1h&G_BIT) && (vbs2h&T_BIT)) diag=values->delta[2][3];\
							}\
						}\
*/



