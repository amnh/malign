/*do_new_build_parallel2 Copyright 1992 Ward Wheeler all rights reserved*/
/*PVM parallel calls of both parent and child so only a single exe*/

#include "align3.h"

void do_initialize_thang(values,a)
alignment **a;
parameters *values;
{
     int i,j,info,the_ok,temp;
     int *reciept,groups_holder;

     reciept=(int *)malloc(values->num_hosts*sizeof(int));
     assert((int) reciept);

     for (i=0;i<values->num_hosts;i++) reciept[i]=0;
     /*initialize*/
     pvm_initsend( PvmDataDefault );
     /*pack values*/
     /*should these be -> or "." because they should be pointers*/
	pack_other_vals(values);
	/*pack other parm files*/
	if (values->number_of_input_alignments > 0) for (i=0;i<values->number_of_input_alignments;i++) {
		if (!values->other_parm[i]) {
			temp=0;
			pvm_pkint(&temp,1,1);
			}
		else {
			temp=1;
			pvm_pkint(&temp,1,1);
			pack_other_vals(values->other_parm[i]);
			}
		}

     /*packing alignments*/
     for (i=0;i<values->all_done;i++) pack_align(a[i]);
     /*send*/
     info=pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_UNPACK_INITIAL_SHIT);

     if (info<0) {
	  fprintf(stderr,"Problems multicasting initial values.\n");
	  /*sending kill message*/
	  pvm_initsend( PvmDataDefault );
	  pvm_pkint(values->tids, values->num_hosts, 1);
	  pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_ALL_MUST_DIE);
	  pvm_exit();
	  exit(-1);
     }

     /*check through recieve*/
     for (i=1;i<values->num_hosts;i++) {
	  reciept[i]=pvm_recv(values->tids[i],PVM_INITIAL_SHIT_RECIEVED);
	  pvm_upkint(&the_ok,1,1);
	  if (!reciept[i]) fprintf(stderr,"Problems with host tid %d.\n",i);
     }
     for (i=1;i<values->num_hosts;i++) if (!reciept[i]) {
	  fprintf(stderr,"Bye.\n");
	  pvm_initsend( PvmDataDefault );
	  pvm_pkint(values->tids, values->num_hosts, 1);
	  pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_ALL_MUST_DIE);
	  pvm_exit();
	  exit(-1);
     }
     free(reciept);
}

alignment **unpack_init_vals(a,values)
alignment **a;
parameters *values;
{
     int i,j,groups_holder,temp;

     get_initial_values(values);
     unpack_other_vals(values);

     values->other_parm=NULL;
	if (values->number_of_input_alignments>0) {
		values->other_parm=(parameters **)malloc(values->number_of_input_alignments*sizeof(parameters *));
		assert((int) values->other_parm);
		for (i=0;i<values->number_of_input_alignments;i++) {
			pvm_upkint(&temp,1,1);
			if (temp) {
				values->other_parm[i]=(parameters *)malloc(sizeof(parameters));
				assert((int) values->other_parm[i]);
				get_initial_values(values->other_parm[i]);
				unpack_other_vals(values->other_parm[i]);
				}
			else values->other_parm[i]=NULL;
			}
		}

     /*unpacking alignments*/
     a=(alignment **)malloc(((2*(values->all_done+1))-1)*sizeof(alignment *));
     assert((int) a);
     for (i=0;i<((2*(values->all_done+1))-1);i++) a[i]=NULL;
     for (i=0;i<values->all_done;i++) a[i]=unpack_align_and_score(values);
     values->previous=NULL;

     if (values->cache_size>0) {
	  values->align_cache=allocate_align_cache(values->cache_size,values,a);
	  fprintf(stderr,"Daughter allocated cache.\n");
     }
     return (a);
}

void unpack_other_vals(values)
parameters *values;
{
int i,j,groups_holder,temp;
     get_initial_values(values);
     pvm_upkint(&(values->get_best), 1, 1);
     pvm_upkint(&(values->get_heur), 1, 1);
     pvm_upkint(&(values->gap_cost), 1, 1);
     pvm_upkint(&(values->change_cost), 1, 1);
     pvm_upkint(&(values->leading_gap_cost), 1, 1);
     pvm_upkint(&(values->trailing_gap_cost), 1, 1);
     pvm_upkint(&(values->phylo_score), 1, 1);
     pvm_upkint(&(values->all_done), 1, 1);
     pvm_upkint(&(values->rand_seed), 1, 1);
     pvm_upkint(&(values->how_many), 1, 1);
     pvm_upkint(&(values->length_at_end), 1, 1);
     pvm_upkint(&(values->randtrees), 1, 1);
     pvm_upkint(&(values->hen_gap), 1, 1);
     pvm_upkint(&(values->acgt), 1, 1);
     pvm_upkint(&(values->farris), 1, 1);
     pvm_upkint(&(values->inter_dig), 1, 1);
     pvm_upkint(&(values->dot), 1, 1);
     pvm_upkint(&(values->keep_trees), 1, 1);
     pvm_upkint(&(values->aquick), 1, 1);
     pvm_upkint(&(values->tree_swap), 1, 1);
     pvm_upkint(&(values->align_swap), 1, 1);
     pvm_upkint(&(values->print_intermediates), 1, 1);
     pvm_upkint(&(values->keep_aligns), 1, 1);
     pvm_upkint(&(values->line_length), 1, 1);
     pvm_upkint(&(values->phylo_gap), 1, 1);
     pvm_upkint(&(values->VERBOSE), 1, 1);
     pvm_upkint(&(values->rep_error), 1, 1);
     pvm_upkint(&(values->paup), 1, 1);
     pvm_upkint(&(values->paup_dot), 1, 1);
     pvm_upkint(&(values->ngroups), 1, 1);
     if (values->ngroups>0) {
	  values->groups=(int **)malloc(values->ngroups*sizeof(int *));
	  assert((int) values->groups);
	  for (i=0;i<values->ngroups;i++) {
	       values->groups[i]=(int *)malloc((values->all_done+1)*sizeof(int));
	       assert((int)values->groups[i]);
	       for (j=0;j<(values->all_done+1);j++) {
		    pvm_upkint(&(values->groups[i][j]),1,1);
	       }
	  }
     }
     pvm_upkint(&(values->in_some), 1, 1);
     pvm_upkint(&(values->number), 1, 1); values->number=0;
     pvm_upkint(&(values->best), 1, 1);
     pvm_upkint(&(values->get_heur2), 1, 1);
     pvm_upkint(&(values->rand_align), 1, 1);
     /*values->all_made where are these things allocated*/
     /*values->tree_made where allocated*/
     pvm_upkint(&(values->time), 1, 1);
     pvm_upkint(&(values->extra_adjustment), 1, 1);
     pvm_upkint(&(values->coding), 1, 1);
     pvm_upkint(&(values->align_multi_swap), 1, 1);
     pvm_upkint(&(values->tree_multi_swap), 1, 1);
     pvm_upkint(&(values->n_codes), 1, 1);
     /*if (values->) pvm_upkint(values->delta,16,1); delta currently diabled anyway*/
     pvm_upkint(&(values->length_best_clad), 1, 1);
     pvm_upkint(&(values->number_best_clad), 1, 1);
     /*values->best_rep where used and allocated*/
     /*values->new_codes not needed*/
     pvm_upkint(&(values->ttr), 1, 1);
     pvm_upkint(&(values->transversion), 1, 1);
     pvm_upkint(&(values->transition), 1, 1);
     pvm_upkint(&(values->n_transv), 1, 1);
     pvm_upkint(&(values->givem), 1, 1);
     pvm_upkint(&(values->matquick), 1, 1);
     pvm_upkint(&(values->max_gap), 1, 1);
     pvm_upkint(&(values->iter), 1, 1);
     pvm_upkint(&(values->low_mem), 1, 1);
     pvm_upkint(&(values->swap_while_add), 1, 1);
     pvm_upkint(&(values->clade_swap_while_add), 1, 1);
     /*values->input names not needed*/
     pvm_upkint(&(values->output_order), 1, 1);
     pvm_upkint(&(values->show_mem), 1, 1);
     pvm_upkint(&(values->max_out_aligns), 1, 1);
     pvm_upkint(&(values->max_out_trees), 1, 1);
     pvm_upkint(&(values->cost_only), 1, 1);
     pvm_upkint(&(values->temp_ti), 1, 1);
     pvm_upkint(&(values->temp_tv), 1, 1);
     pvm_upkint(&(values->temp_changes), 1, 1);
     pvm_upkint(&(values->temp_gaps), 1, 1);
     pvm_upkint(&(values->overall_ti), 1, 1);
     pvm_upkint(&(values->overall_tv), 1, 1);
     pvm_upkint(&(values->overall_changes), 1, 1);
     pvm_upkint(&(values->overall_gaps), 1, 1);
     pvm_upkint(&(values->temp_ti_g), 1, 1);
     pvm_upkint(&(values->temp_tv_g), 1, 1);
     pvm_upkint(&(values->temp_changes_g), 1, 1);
     pvm_upkint(&(values->temp_gaps_g), 1, 1);
     pvm_upkint(&(values->overall_ti_g), 1, 1);
     pvm_upkint(&(values->overall_tv_g), 1, 1);
     pvm_upkint(&(values->overall_changes_g), 1, 1);
     pvm_upkint(&(values->overall_gaps_g), 1, 1);
     pvm_upkint(&(values->get_weights), 1, 1);
     pvm_upkint(&(values->get_weights_2), 1, 1);
     pvm_upkint(&(values->weight_range), 1, 1);
     pvm_upkint(&(values->cache_size), 1, 1);
     /*values->align_cache allocated by daughters*/
     pvm_upkint(&(values->cache_info), 1, 1);
     pvm_upkint(&(values->saw_farris), 1, 1);
     pvm_upkint(&(values->saw_goloboff), 1, 1);
     pvm_upkint(&(values->saw_wheeler), 1, 1);
     pvm_upkint(&(values->input_is_align), 1, 1);
     /*values->b_string1 allocated by daughters*/
     /*values->b_string2 allocated by daughters*/
     pvm_upkint(&(values->phylo_time), 1, 1);
     pvm_upkint(&(values->pref_direc), 1, 1);
     pvm_upkint(&(values->min_num_gaps), 1, 1);
     pvm_upkint(&(values->gap_must_cost), 1, 1);
     pvm_upkint(&(values->rand_order), 1, 1);
     pvm_upkint(&(values->reorder), 1, 1);
     pvm_upkint(&(values->tree_rand_order_max), 1, 1);
     pvm_upkint(&(values->number_of_input_alignments), 1, 1);
     /*values->input_file_names not needed*/
     /*values->align_file not needed*/
     pvm_upkint(&(values->actual_num_sequences), 1, 1);
     pvm_upkint(&(values->cull), 1, 1);
     pvm_upkint(&(values->elision), 1, 1);
     pvm_upkint(&(values->get_best), 1, 1);
     /*ce_weights not needed*/
     /*previous where allocated ?*/
     pvm_upkint(&(values->grain_size), 1, 1);
     pvm_upkint(&(values->process_factor),1,1);
     pvm_upkint(&(values->new_optimization),1,1);
     pvm_upkint(&(values->no_leading_or_trailing_cost),1,1);
     pvm_upkint(&(values->gap_in_hypanc),1,1);
     if (values->number_of_input_alignments) {
	  values->data_sets=(int *)malloc(values->number_of_input_alignments*sizeof(int));
	  assert((int) values->data_sets);
       pvm_upkint((values->data_sets),values->number_of_input_alignments,1);
	  if (1) {
	       values->data_set_weights=(int *)malloc(values->number_of_input_alignments*sizeof(int));
	       assert((int) values->data_set_weights);
	       pvm_upkint((values->data_set_weights),values->number_of_input_alignments,1);
	       }
	  }
     pvm_upkint(&(values->align_node_swap),1,1);
     pvm_upkint(&(values->align_root_swap),1,1);
     pvm_upkint(&(values->align_complete_root_swap),1,1);
     pvm_upkint(&(values->align_partial_root_swap),1,1);
     pvm_upkint(&(values->dump_parameters),1,1);
     pvm_upkint(&(values->groups_as_start),1,1);
     pvm_upkint(&(values->in_bandb_loop),1,1);
     pvm_upkint(&(values->current_bound),1,1);
     pvm_upkint(&(values->chop),1,1);
     pvm_upkint(&(values->optimize_nodes),1,1);
     pvm_upkint(&(values->in_optimize),1,1);
     pvm_upkint(&(values->apolist),1,1);
     pvm_upkint(&(values->rand_apo),1,1);
     pvm_upkint(&(values->sbr),1,1);
     pvm_upkint(&(values->tbr),1,1);
     pvm_upkint(&(values->shortcut_tree),1,1);
     pvm_upkint(&(values->shortcut_tree2),1,1);
     pvm_upkint(&(values->get_heur3),1,1);
     pvm_upkint(&(values->asbr),1,1);
     pvm_upkint(&(values->atbr),1,1);
     pvm_upkint(&(values->tbr_align_shortcut),1,1);
	pvm_upkint(&(values->collapse),1,1);
     pvm_upkint(&(values->arrt),1,1);
     pvm_upkint(&(values->do_dave_align),1,1);
     pvm_upkint(&(values->expand_X),1,1);
     pvm_upkint(&(values->jackboot),1,1);
     if (values->jackboot) {
	     if (values->jack_array) free(values->jack_array);
	     values->jack_array=(int *)malloc(1000*sizeof(int));
	     assert((int) values->jack_array);
     	pvm_upkint(values->jack_array,1000,1);
     	}
     values->other_parm=NULL;
     for (i=0;i<32;i++) for (j=0;j<32;j++) {
     	pvm_upkint(&temp,1,1);
     	values->lookup[i][j].base=temp;
     	pvm_upkint(&temp,1,1);
     	values->lookup[i][j].union_bit=temp;
     	pvm_upkint(&temp,1,1);
     	values->lookup[i][j].cost=temp;
     	}
     pvm_upkint(&(values->best_order),1,1);
}

void pack_other_vals(values)
parameters *values;
{
    int i,j,info,the_ok,temp;
	int groups_holder;

     pvm_pkint(&(values->get_best), 1, 1);
     pvm_pkint(&(values->get_heur), 1, 1);
     pvm_pkint(&(values->gap_cost), 1, 1);
     pvm_pkint(&(values->change_cost), 1, 1);
     pvm_pkint(&(values->leading_gap_cost), 1, 1);
     pvm_pkint(&(values->trailing_gap_cost), 1, 1);
     pvm_pkint(&(values->phylo_score), 1, 1);
     pvm_pkint(&(values->all_done), 1, 1);
     pvm_pkint(&(values->rand_seed), 1, 1);
     pvm_pkint(&(values->how_many), 1, 1);
     pvm_pkint(&(values->length_at_end), 1, 1);
     pvm_pkint(&(values->randtrees), 1, 1);
     pvm_pkint(&(values->hen_gap), 1, 1);
     pvm_pkint(&(values->acgt), 1, 1);
     pvm_pkint(&(values->farris), 1, 1);
     pvm_pkint(&(values->inter_dig), 1, 1);
     pvm_pkint(&(values->dot), 1, 1);
     pvm_pkint(&(values->keep_trees), 1, 1);
     pvm_pkint(&(values->aquick), 1, 1);
     pvm_pkint(&(values->tree_swap), 1, 1);
     pvm_pkint(&(values->align_swap), 1, 1);
     pvm_pkint(&(values->print_intermediates), 1, 1);
     pvm_pkint(&(values->keep_aligns), 1, 1);
     pvm_pkint(&(values->line_length), 1, 1);
     pvm_pkint(&(values->phylo_gap), 1, 1);
     pvm_pkint(&(values->VERBOSE), 1, 1);
     pvm_pkint(&(values->rep_error), 1, 1);
     pvm_pkint(&(values->paup), 1, 1);
     pvm_pkint(&(values->paup_dot), 1, 1);
     pvm_pkint(&(values->ngroups), 1, 1);
     if (values->ngroups>0) {
	  for (i=0;i<values->ngroups;i++) {
	       for (j=0;j<(values->all_done+1);j++) {
		    pvm_pkint(&(values->groups[i][j]),1,1);
	       }
	  }
     }
     pvm_pkint(&(values->in_some), 1, 1);
     pvm_pkint(&(values->number), 1, 1);
     pvm_pkint(&(values->best), 1, 1);
     pvm_pkint(&(values->get_heur2), 1, 1);
     pvm_pkint(&(values->rand_align), 1, 1);
     /*values->all_made*/
     /*values->tree_made*/
     pvm_pkint(&(values->time), 1, 1);
     pvm_pkint(&(values->extra_adjustment), 1, 1);
     pvm_pkint(&(values->coding), 1, 1);
     pvm_pkint(&(values->align_multi_swap), 1, 1);
     pvm_pkint(&(values->tree_multi_swap), 1, 1);
     pvm_pkint(&(values->n_codes), 1, 1);
     /*if (values->delta) pvm_pkint(values->delta,16,1); delta*/
     pvm_pkint(&(values->length_best_clad), 1, 1);
     pvm_pkint(&(values->number_best_clad), 1, 1);
     /*values->best_rep*/
     /*values->new_codes*/
     pvm_pkint(&(values->ttr), 1, 1);
     pvm_pkint(&(values->transversion), 1, 1);
     pvm_pkint(&(values->transition), 1, 1);
     pvm_pkint(&(values->n_transv), 1, 1);
     pvm_pkint(&(values->givem), 1, 1);
     pvm_pkint(&(values->matquick), 1, 1);
     pvm_pkint(&(values->max_gap), 1, 1);
     pvm_pkint(&(values->iter), 1, 1);
     pvm_pkint(&(values->low_mem), 1, 1);
     pvm_pkint(&(values->swap_while_add), 1, 1);
     pvm_pkint(&(values->clade_swap_while_add), 1, 1);
     /*values->input names*/
     pvm_pkint(&(values->output_order), 1, 1);
     pvm_pkint(&(values->show_mem), 1, 1);
     pvm_pkint(&(values->max_out_aligns), 1, 1);
     pvm_pkint(&(values->max_out_trees), 1, 1);
     pvm_pkint(&(values->cost_only), 1, 1);
     pvm_pkint(&(values->temp_ti), 1, 1);
     pvm_pkint(&(values->temp_tv), 1, 1);
     pvm_pkint(&(values->temp_changes), 1, 1);
     pvm_pkint(&(values->temp_gaps), 1, 1);
     pvm_pkint(&(values->overall_ti), 1, 1);
     pvm_pkint(&(values->overall_tv), 1, 1);
     pvm_pkint(&(values->overall_changes), 1, 1);
     pvm_pkint(&(values->overall_gaps), 1, 1);
     pvm_pkint(&(values->temp_ti_g), 1, 1);
     pvm_pkint(&(values->temp_tv_g), 1, 1);
     pvm_pkint(&(values->temp_changes_g), 1, 1);
     pvm_pkint(&(values->temp_gaps_g), 1, 1);
     pvm_pkint(&(values->overall_ti_g), 1, 1);
     pvm_pkint(&(values->overall_tv_g), 1, 1);
     pvm_pkint(&(values->overall_changes_g), 1, 1);
     pvm_pkint(&(values->overall_gaps_g), 1, 1);
     pvm_pkint(&(values->get_weights), 1, 1);
     pvm_pkint(&(values->get_weights_2), 1, 1);
     pvm_pkint(&(values->weight_range), 1, 1);
     pvm_pkint(&(values->cache_size), 1, 1);
     /*values->align_cache*/
     pvm_pkint(&(values->cache_info), 1, 1);
     pvm_pkint(&(values->saw_farris), 1, 1);
     pvm_pkint(&(values->saw_goloboff), 1, 1);
     pvm_pkint(&(values->saw_wheeler), 1, 1);
     pvm_pkint(&(values->input_is_align), 1, 1);
     /*values->b_string1*/
     /*values->b_string2*/
     pvm_pkint(&(values->phylo_time), 1, 1);
     pvm_pkint(&(values->pref_direc), 1, 1);
     pvm_pkint(&(values->min_num_gaps), 1, 1);
     pvm_pkint(&(values->gap_must_cost), 1, 1);
     pvm_pkint(&(values->rand_order), 1, 1);
     pvm_pkint(&(values->reorder), 1, 1);
     pvm_pkint(&(values->tree_rand_order_max), 1, 1);
     pvm_pkint(&(values->number_of_input_alignments), 1, 1);
     /*values->input_file_names*/
     /*values->align_file*/
     pvm_pkint(&(values->actual_num_sequences), 1, 1);
     pvm_pkint(&(values->cull), 1, 1);
     pvm_pkint(&(values->elision), 1, 1);
     pvm_pkint(&(values->get_best), 1, 1);
     /*ce_weights*/
     /*previous*/
     pvm_pkint(&(values->grain_size), 1, 1);
     pvm_pkint(&(values->process_factor),1,1);
     pvm_pkint(&(values->new_optimization),1,1);
     pvm_pkint(&(values->no_leading_or_trailing_cost),1,1);
     pvm_pkint(&(values->gap_in_hypanc),1,1);
     if (values->number_of_input_alignments) {
	  pvm_pkint((values->data_sets),values->number_of_input_alignments,1);
	  if (1) pvm_pkint((values->data_set_weights),values->number_of_input_alignments,1);
	  }
     pvm_pkint(&(values->align_node_swap),1,1);
     pvm_pkint(&(values->align_root_swap),1,1);
     pvm_pkint(&(values->align_complete_root_swap),1,1);
     pvm_pkint(&(values->align_partial_root_swap),1,1);
     pvm_pkint(&(values->dump_parameters),1,1);
     pvm_pkint(&(values->groups_as_start),1,1);
     pvm_pkint(&(values->in_bandb_loop),1,1);
     pvm_pkint(&(values->current_bound),1,1);
     pvm_pkint(&(values->chop),1,1);
     pvm_pkint(&(values->optimize_nodes),1,1);
     pvm_pkint(&(values->in_optimize),1,1);
     pvm_pkint(&(values->apolist),1,1);
     pvm_pkint(&(values->rand_apo),1,1);
     pvm_pkint(&(values->sbr),1,1);
     pvm_pkint(&(values->tbr),1,1);
     pvm_pkint(&(values->shortcut_tree),1,1);
     pvm_pkint(&(values->shortcut_tree2),1,1);
     pvm_pkint(&(values->get_heur3),1,1);
     pvm_pkint(&(values->asbr),1,1);
     pvm_pkint(&(values->atbr),1,1);
     pvm_pkint(&(values->tbr_align_shortcut),1,1);
	pvm_pkint(&(values->collapse),1,1);
     pvm_pkint(&(values->arrt),1,1);
     pvm_pkint(&(values->do_dave_align),1,1);
     pvm_pkint(&(values->expand_X),1,1);
     pvm_pkint(&(values->jackboot),1,1);
     if (values->jackboot) pvm_pkint(values->jack_array,1000,1);
     for (i=0;i<32;i++) for (j=0;j<32;j++) {
     	temp=values->lookup[i][j].base;
     	pvm_pkint(&temp,1,1);
     	temp=values->lookup[i][j].union_bit;
     	pvm_pkint(&temp,1,1);
     	temp=values->lookup[i][j].cost;
     	pvm_pkint(&temp,1,1);

     	}
     pvm_pkint(&(values->best_order),1,1);
}

void do_alignment_now(d_nodes,nodes,d_tree_rep,a,values)
int **nodes, **d_nodes;
int *d_tree_rep;
alignment **a;
parameters *values;
{
     int grain,d_ntaxa,d_i,index;
     int i,j;
     int current_score;
     /*fprintf(stderr,"receieved and starting ");*/
     pvm_upkint(&grain,1,1);
     pvm_upkint(&d_ntaxa,1,1);
     pvm_upkint(&d_i,1,1);
     pvm_upkint(d_tree_rep,d_i+1,1);
     pvm_upkint(&index,1,1);
     nodes[0][0]=0;
     nodes[0][1]=d_ntaxa+1;
     nodes[1][0]=1;
     nodes[1][1]=2;

     d_nodes[0][0]=0;
     d_nodes[0][1]=d_ntaxa+1;
     d_nodes[1][0]=1;
     d_nodes[1][1]=2;

     /*PVM do align*/
     current_score=all_make_nodes(a,d_tree_rep,nodes,d_ntaxa,d_i+4,values);
     /*if only one then can simple make a[d_nodes[0][1]]=a[d_nodes[0][1]];*/
     if (current_score==HUGE_COST) {
	  pvm_initsend(PvmDataDefault);
	  pvm_pkint(&grain,1,1);
	  pvm_pkint(&index,1,1);
	  pvm_pkint(&(current_score),1,1);
     }
     else {
	  all_other_make_nodes2(d_tree_rep,nodes,d_ntaxa,d_i+4,d_nodes);
	  /*if (current_alignment) dump_align(current_alignment);
			 current_alignment=make_align(a[d_nodes[0][1]]);*/
	  /*PVM send back*/
	  pvm_initsend(PvmDataDefault);
	  pvm_pkint(&grain,1,1);
	  pvm_pkint(&index,1,1);
	  pack_align(a[d_nodes[0][1]]);
     }
     pvm_send(values->tids[0],PVM_PACK_ALIGN_AND_SCORE);
     /*pvm dump align but no need if only assigned*/
}

void do_parallel_shtick(a,values,i,ntaxa,j,intermediate_scores,intermediate_best_score,intermediate_align,temp_rep)
int *i, *ntaxa,*j,*intermediate_scores,*intermediate_best_score, *temp_rep;
alignment **a,**intermediate_align;
parameters *values;
{
     int max_to_do,num_left,k,grain,info, index;
     int bufid, bytes,type,source;

     grain=values->grain_size;

     /*fprintf(stderr,"Sending jobs to");*/
     num_left=max_to_do=(2*(*i))+3;
     for ((*j)=1;(((*j)<values->num_hosts) && ((*j)<max_to_do));(*j)++) {
	  temp_rep[(*i)]=(*j);
	  /*fprintf(stderr," %d (grain=%d ntaxa=%d i=%d temp_rep[%d]=%d) ",(*j),grain,(*ntaxa),(*i),(*i),temp_rep[(*i)]);*/
	  pvm_initsend( PvmDataDefault );
	  pvm_pkint(&grain,1,1);
	  pvm_pkint(&(*ntaxa),1,1);
	  pvm_pkint(&(*i),1,1);
	  pvm_pkint(temp_rep,(*i)+1,1);
	  pvm_pkint(&(*j),1,1);
	  pvm_send(values->tids[(*j)],PVM_DO_ALIGNMENT_NOW);
	  /*fprintf(stderr," sent\n");*/
     }
     while (num_left) {
	  /*fprintf(stderr,"num_left=%d ",num_left);*/
	  bufid=pvm_recv(-1,PVM_PACK_ALIGN_AND_SCORE);
	  if (bufid) {
	       info=pvm_bufinfo(bufid, &bytes, &type, &source);
	       if ((*j)<=max_to_do) {
		    temp_rep[(*i)]=(*j);
		    pvm_initsend( PvmDataDefault );
		    pvm_pkint(&grain,1,1);
		    pvm_pkint(&(*ntaxa),1,1);
		    pvm_pkint(&(*i),1,1);
		    pvm_pkint(temp_rep,(*i)+1,1);
		    pvm_pkint(&(*j),1,1);
		    pvm_send(source,PVM_DO_ALIGNMENT_NOW);
		    ++(*j);
	       }
	       --num_left;
	       pvm_upkint(&grain,1,1);
	       pvm_upkint(&index,1,1);
	       pvm_upkint(&(intermediate_scores[index]),1,1);
	       /*
			      fprintf(stderr,"AlignScore[%d]=%d ",index,intermediate_scores[index]);
			      */
	       if (intermediate_scores[index]<=(*intermediate_best_score)) {
		    if (intermediate_scores[index]<(*intermediate_best_score)) {
			 (*intermediate_best_score)=intermediate_scores[index];
			 if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
			 intermediate_align[index]=unpack_align(intermediate_scores[index],values);
		    }
		    else if ((!values->aquick) && (intermediate_scores[index]!=HUGE_COST)) {
			 /* code folded from here */
			 if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
			 intermediate_align[index]=unpack_align(intermediate_scores[index],values);
		    }


		    else pvm_freebuf(bufid);
	       }/*(*j)<*/
	  }/*if message recieved*/
     }/*num_left*/
}

void do_parallel_shtick_swap(a,values,i,ntaxa,j,intermediate_scores,intermediate_best_score,intermediate_align,temp_rep,old_rep,
z)
int *i, *ntaxa,*j,*intermediate_scores,*intermediate_best_score, *temp_rep, **old_rep, *z;
alignment **a,**intermediate_align;
parameters *values;
{
     int max_to_do,num_left,k,grain,info, index;
     int bufid, bytes,type,source;
     int jj,sub_taxa;

     grain=values->grain_size;
     /*fprintf(stderr,"Sending jobs to");*/
     num_left=(2*(*i))+3-1;/*the -1 is for the old_rep thing so swapping on one less than possible*/
     max_to_do=(2*(*i))+3; /*here no since need to go through all j's*/
     jj=1;
     for ((*j)=1;(((*j)<values->num_hosts) && ((*j)<max_to_do));(*j)++) if ((*j)!=old_rep[(*z)][(*i)]) {
	  temp_rep[(*i)]=(*j);
	  /*fprintf(stderr," %d (grain=%d ntaxa=%d i=%d temp_rep[%d]=%d) ",(*j),grain,(*ntaxa),(*i),(*i),temp_rep[(*i)]);*/
	  pvm_initsend( PvmDataDefault );
	  pvm_pkint(&grain,1,1);
	  pvm_pkint(&(*ntaxa),1,1);
	  sub_taxa=(*ntaxa)-4;
	  pvm_pkint(&sub_taxa,1,1);
	  pvm_pkint(temp_rep,sub_taxa+1,1);
	  pvm_pkint(&(temp_rep[(*i)]),1,1);
	  pvm_send(values->tids[jj++],PVM_DO_ALIGNMENT_NOW);/*otherwise will loose hosts*/
	  /*fprintf(stderr," sent\n");*/
     }
     while (num_left) {
	  /*fprintf(stderr,"num_left=%d ",num_left);*/
	  bufid=pvm_recv(-1,PVM_PACK_ALIGN_AND_SCORE);
	  if (bufid) {
	       info=pvm_bufinfo(bufid, &bytes, &type, &source);
	       if ((*j)<=max_to_do) {
		    if ((*j)!=old_rep[(*z)][(*i)]) {
			 temp_rep[(*i)]=(*j);
			 pvm_initsend( PvmDataDefault );
			 pvm_pkint(&grain,1,1);
			 pvm_pkint(&(*ntaxa),1,1);
			 sub_taxa=(*ntaxa)-4;
			 pvm_pkint(&sub_taxa,1,1);
			 pvm_pkint(temp_rep,sub_taxa+1,1);
			 pvm_pkint(&(temp_rep[(*i)]),1,1);
			 pvm_send(source,PVM_DO_ALIGNMENT_NOW);
		    }
		    ++(*j);
	       }
	       --num_left;
	       pvm_upkint(&grain,1,1);
	       pvm_upkint(&index,1,1);
	       pvm_upkint(&(intermediate_scores[index]),1,1);
	       if (intermediate_scores[index]<=(*intermediate_best_score)) {
		    /* code folded from here */
		    if (intermediate_scores[index]<(*intermediate_best_score)) {
			 (*intermediate_best_score)=intermediate_scores[index];
			 if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
			 intermediate_align[index]=unpack_align(intermediate_scores[index],values);
		    }
		    else if ((!values->aquick) && (intermediate_scores[index]!=HUGE_COST)) {
			 if (intermediate_align[index]) intermediate_align[index]=dump_align(intermediate_align[index]);
			 intermediate_align[index]=unpack_align(intermediate_scores[index],values);
		    }
		    /* unfolding */


else pvm_freebuf(bufid);
	       }/*(*j)<*/
	  } /*if receiecved*/
     }/*num_left*/

}

void do_alignment_score_only(d_nodes,nodes,d_tree_rep,a,values)
int **nodes, **d_nodes;
int *d_tree_rep;
alignment **a;
parameters *values;
{
     int grain,d_ntaxa,d_i,index;
     int i,j;
     int current_score;
     /*fprintf(stderr,"receieved and starting ");*/
     pvm_upkint(&grain,1,1);
     pvm_upkint(&d_ntaxa,1,1);
     pvm_upkint(&d_i,1,1);
     pvm_upkint(d_tree_rep,d_i+1,1);
     pvm_upkint(&index,1,1);
     nodes[0][0]=0;
     nodes[0][1]=d_ntaxa+1;
     nodes[1][0]=1;
     nodes[1][1]=2;

     d_nodes[0][0]=0;
     d_nodes[0][1]=d_ntaxa+1;
     d_nodes[1][0]=1;
     d_nodes[1][1]=2;

     /*PVM do align*/
     current_score=all_make_nodes(a,d_tree_rep,nodes,d_ntaxa,d_i+4,values);
     /*if only one then can simple make a[d_nodes[0][1]]=a[d_nodes[0][1]];*/
     pvm_initsend(PvmDataDefault);
     pvm_pkint(&grain,1,1);
     pvm_pkint(&index,1,1);
     pvm_pkint(&current_score,1,1);
     pvm_send(values->tids[0],PVM_SEND_SCORE_ONLY);

}

void pack_align(a)
alignment *a;
{
int i,j;

	  pvm_pkint(&(a->score),1,1);
	  pvm_pkint(&(a->n_seqs),1,1);
	  pvm_pkint(&(a->length),1,1);
	  pvm_pkstr(a->name);
	  if (a->type_weight!=1) a->type_weight=0;
	  pvm_pkint(&(a->type_weight),1,1);
	  for (j=0;j<a->n_seqs;j++) pvm_pkstr(a->taxon_name[j]);
	  for (j=0;j<(a->n_seqs+a->type_weight);j++) pvm_pkstr(a->s[j]);

}

alignment *unpack_align(current_score,values)
int current_score;
parameters *values;
{
alignment *a;
int i,k;

a=(alignment *)malloc(sizeof(alignment));
assert ((int) a);
a->score=current_score;
pvm_upkint(&(a->n_seqs),1,1);
pvm_upkint(&(a->length),1,1);
a->name=(char *)malloc(FIFTEEN_BITS_MAX*sizeof(char));
assert((int)a->name);
pvm_upkstr(a->name);
a->name=(char *)realloc(a->name,(1+strlen(a->name))*sizeof(char));
assert((int)a->name);
pvm_upkint(&(a->type_weight),1,1);
if (a->type_weight!=1) a->type_weight=0;
a->taxon_name=(char **)malloc((a->n_seqs)*sizeof(char *));
assert((int)a->taxon_name);
for (k=0;k<a->n_seqs;k++) {
     a->taxon_name[k]=(char *)malloc(FIFTEEN_BITS_MAX*sizeof(char *));
     assert((int) a->taxon_name[k]);
     pvm_upkstr(a->taxon_name[k]);
     a->taxon_name[k]=(char *)realloc(a->taxon_name[k],(1+strlen(a->taxon_name[k]))*sizeof(char));
     assert((int) a->taxon_name[k]);
}
a->s=(char **)malloc((a->n_seqs+a->type_weight)*sizeof(char *));
assert((int)a->s);
for (k=0;k<(a->n_seqs+a->type_weight);k++) {
     a->s[k]=(char *)malloc(MAX_SEQUENCE_SIZE*sizeof(char *));
     assert((int) a->s[k]);
     pvm_upkstr(a->s[k]);
     a->s[k]=(char *)realloc(a->s[k],(1+strlen(a->s[k]))*sizeof(char));
     assert((int) a->s[k]);
}


return a;
}

alignment *unpack_align_and_score(values)
parameters *values;
{
alignment *a;
int i,k;

a=(alignment *)malloc(sizeof(alignment));
assert ((int) a);
pvm_upkint(&(a->score),1,1);
pvm_upkint(&(a->n_seqs),1,1);
if (!a->n_seqs) {fprintf(stderr,"UnPacking 0 sequence alignment\n");exit(-1);}
pvm_upkint(&(a->length),1,1);
if (!a->length) {fprintf(stderr,"UnPacking 0 length alignment\n");exit(-1);}
a->name=(char *)malloc(FIFTEEN_BITS_MAX*sizeof(char));
assert((int)a->name);
pvm_upkstr(a->name);
a->name=(char *)realloc(a->name,(1+strlen(a->name))*sizeof(char));
assert((int)a->name);
pvm_upkint(&(a->type_weight),1,1);
if (a->type_weight!=1) a->type_weight=0;
a->taxon_name=(char **)malloc((a->n_seqs)*sizeof(char *));
assert((int)a->taxon_name);
for (k=0;k<a->n_seqs;k++) {
     a->taxon_name[k]=(char *)malloc(FIFTEEN_BITS_MAX*sizeof(char *));
     assert((int) a->taxon_name[k]);
     pvm_upkstr(a->taxon_name[k]);
     a->taxon_name[k]=(char *)realloc(a->taxon_name[k],(1+strlen(a->taxon_name[k]))*sizeof(char));
     assert((int) a->taxon_name[k]);
}
a->s=(char **)malloc((a->n_seqs+a->type_weight)*sizeof(char *));
assert((int)a->s);
for (k=0;k<(a->n_seqs+a->type_weight);k++) {
     a->s[k]=(char *)malloc(MAX_SEQUENCE_SIZE*sizeof(char *));
     assert((int) a->s[k]);
     pvm_upkstr(a->s[k]);
     a->s[k]=(char *)realloc(a->s[k],(1+strlen(a->s[k]))*sizeof(char));
     assert((int) a->s[k]);
	}
return a;
}

void do_diagnose_alignment(nodes,a,values)
int **nodes;
alignment **a;
parameters *values;
{
int i1,j1;
int grain, ntaxa, score;
int **cur_groups, test=1, ngroups2;

/*unpack*/
pvm_upkint(&grain,1,1);
pvm_upkint(&ntaxa,1,1);
for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_upkint(&(nodes[i1][j1]),1,1);

if (!values->groups) score=all_diagnose_tree(a,nodes,ntaxa,ntaxa,values);
else {
	cur_groups=get_cur_groups(nodes,ntaxa,ntaxa,&ngroups2);
	test=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntaxa-1,ntaxa-1);
	for (i1=0;i1<ngroups2;i1++) free(cur_groups[i1]);
	free(cur_groups);
	if (test) score=all_diagnose_tree(a,nodes,ntaxa,ntaxa,values);
	else score=HUGE_COST;
	}

/*packresults*/
pvm_initsend(PvmDataDefault);
pvm_pkint(&grain,1,1);
if (test) {
	pack_align(a[nodes[0][1]]);
	/*pack nodes for results later*/
	for (i1=0;i1<ntaxa-1;i1++) for (j1=0;j1<2;j1++) pvm_pkint(&(nodes[i1][j1]),1,1);
	}
else pvm_pkint(&score,1,1);

/*send*/
pvm_send(values->tids[0],PVM_SEND_DIAGNOSE_ALIGNMENT);

}

void do_alignment_score_only_bandb(d_nodes,nodes,d_tree_rep,a,values)
int **nodes, **d_nodes;
int *d_tree_rep;
alignment **a;
parameters *values;
{
     int grain,d_ntaxa,d_i,index;
     int i,j;
     int current_score;
     /*fprintf(stderr,"receieved and starting ");*/
     values->in_bandb_loop=1;

     pvm_upkint(&grain,1,1);
     pvm_upkint(&(values->current_bound),1,1);
     pvm_upkint(&d_ntaxa,1,1);
     pvm_upkint(&d_i,1,1);
     pvm_upkint(d_tree_rep,d_i+1,1);
     pvm_upkint(&index,1,1);
     nodes[0][0]=0;
     nodes[0][1]=d_ntaxa+1;
     nodes[1][0]=1;
     nodes[1][1]=2;

     d_nodes[0][0]=0;
     d_nodes[0][1]=d_ntaxa+1;
     d_nodes[1][0]=1;
     d_nodes[1][1]=2;

     /*PVM do align*/
     current_score=all_make_nodes(a,d_tree_rep,nodes,d_ntaxa,d_i+4,values);
     /*if only one then can simple make a[d_nodes[0][1]]=a[d_nodes[0][1]];*/
     pvm_initsend(PvmDataDefault);
     pvm_pkint(&grain,1,1);
     pvm_pkint(&index,1,1);
     pvm_pkint(&current_score,1,1);
     pvm_send(values->tids[0],PVM_SEND_SCORE_ONLY);
     values->in_bandb_loop=0;

}

void do_alignment_now_bandb(d_nodes,nodes,d_tree_rep,a,values)
int **nodes, **d_nodes;
int *d_tree_rep;
alignment **a;
parameters *values;
{
     int grain,d_ntaxa,d_i,index;
     int i,j;
     int current_score;

     values->in_bandb_loop=1;
     /*fprintf(stderr,"receieved and starting ");*/
     pvm_upkint(&grain,1,1);
     pvm_upkint(&(values->current_bound),1,1);
     pvm_upkint(&d_ntaxa,1,1);
     pvm_upkint(&d_i,1,1);
     pvm_upkint(d_tree_rep,d_i+1,1);
     pvm_upkint(&index,1,1);
     nodes[0][0]=0;
     nodes[0][1]=d_ntaxa+1;
     nodes[1][0]=1;
     nodes[1][1]=2;

     d_nodes[0][0]=0;
     d_nodes[0][1]=d_ntaxa+1;
     d_nodes[1][0]=1;
     d_nodes[1][1]=2;

     /*PVM do align*/
     current_score=all_make_nodes(a,d_tree_rep,nodes,d_ntaxa,d_i+4,values);
     /*if only one then can simple make a[d_nodes[0][1]]=a[d_nodes[0][1]];*/
     if (current_score==HUGE_COST) {
	  pvm_initsend(PvmDataDefault);
	  pvm_pkint(&grain,1,1);
	  pvm_pkint(&index,1,1);
	  pvm_pkint(&(current_score),1,1);
     }
     else {
	  all_other_make_nodes2(d_tree_rep,nodes,d_ntaxa,d_i+4,d_nodes);
	  /*if (current_alignment) dump_align(current_alignment);
			 current_alignment=make_align(a[d_nodes[0][1]]);*/
	  /*PVM send back*/
	  pvm_initsend(PvmDataDefault);
	  pvm_pkint(&grain,1,1);
	  pvm_pkint(&index,1,1);
	  pack_align(a[d_nodes[0][1]]);
     }
     pvm_send(values->tids[0],PVM_PACK_ALIGN_AND_SCORE);
     /*pvm dump align but no need if only assigned*/
     values->in_bandb_loop=0;

}

void do_new_build_parallel(a,ntax,values)
alignment **a;
int ntax;
parameters *values;
{
int time_to_die,tag,info,bufid,num_bytes,tid;
int i,j,ll,l,value,*temp_rep, cur_val,n,m;
int initial_length,d1,d2;
int *anc,check;
int *tree_rep,score_holder;
alignment ***a_up,*temp_align,*temp_align2;
alignment *cur_align,*cur_align2;
int d3,k,**nodes;
alignment *t1,*t2,*t3;
int ii;

/*tree_rep=(int *)malloc((ntax-3)*sizeof(int));
assert((int) tree_rep);*/
anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);
/*temp_rep=(int *)malloc((ntax-3)*sizeof(int));
assert((int) temp_rep);*/
a_up=(alignment ***)malloc(((2*ntax)-1)*sizeof(alignment **));
assert((int) a_up);
for (i=0;i<2*ntax-1;i++) {
	a_up[i]=(alignment **)malloc(2*sizeof(alignment *));
	assert((int) a_up[i]);
	a_up[i][0]=NULL;
	a_up[i][1]=NULL;
	}
cur_align=NULL;
cur_align2=NULL;
temp_align=NULL;
temp_align2=NULL;
t1=NULL;
t2=NULL;
t3=NULL;

nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) nodes);
for (i=0;i<ntax-1;i++) {
	nodes[i]=(int *)malloc(2*sizeof(int));
	assert((int) nodes);
	}

/*set so no automatic cladogram*/
score_holder=values->phylo_score;
values->phylo_score=0;
values->in_some=1;

/*get initial tree REMEMBER all < ntax are shifted one*/
nodes[0][0]=0;/*this should never be referenced*/
nodes[0][1]=ntax+1;
nodes[1][0]=1;
nodes[1][1]=2;

anc[0]=ntax;
anc[1]=ntax+1;
anc[2]=ntax+1;
anc[ntax+1]=ntax;

/* and optimized nodes
a[ntax+1]=nw(a[0],a[1],values);

initial_length=a[ntax+1]->score;
a_up[ntax+1][0]=make_align(a[1]);
a_up[ntax+1][1]=make_align(a[0]);
*/
time_to_die=0;
while (!time_to_die) {
	bufid=pvm_recv(values->tids[0],-1);
     info=pvm_bufinfo(bufid,&num_bytes, &tag,&tid);
	if (tag==PVM_NEW_BUILD_DONE) time_to_die=1; /*die*/
	else if (tag==PVM_NEW_BUILD_ADD)  { /*add_thang*/
		pvm_upkint(&j,1,1);
			pvm_upkint(&d1,1,1);
		pvm_upkint(&d2,1,1);
			pvm_upkint(&d3,1,1);
			pvm_upkint(&ii,1,1);
			/*if (t1) t1=dump_align(t1);
			t1=unpack_align_and_score(values);*/
			if (t2) t2=dump_align(t2);
			t2=unpack_align_and_score(values);
			if (d1>0) {
			if (t3) t3=dump_align(t3);
			t3=unpack_align_and_score(values);
			}
		temp_align=nw(a[ii+2],t2,values);/*would be ii+3 but the adjustment for jack-up*/
		if (d1>0) {
			if (temp_align2) temp_align2=dump_align(temp_align2);
			temp_align2=nw(temp_align,t3,values);
			if (score_holder>0) {
				if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
				else {
					values->phylo_score=score_holder;
					cur_val=temp_align2->score=cladogram(temp_align2,values);
					values->phylo_score=0;
					if (values->previous) values->previous=dump_align(values->previous);
					values->previous=make_align(temp_align2);
					}
				}
			else cur_val=temp_align2->score;
			}
		else {
			if (score_holder>0) {
				if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
				else {
					values->phylo_score=score_holder;
					cur_val=temp_align->score=cladogram(temp_align,values);
					values->phylo_score=0;
					if (values->previous) values->previous=dump_align(values->previous);
					values->previous=make_align(temp_align);
					}
				}
			else cur_val=temp_align->score;
			}
		pvm_initsend( PvmDataDefault );
	  pvm_pkint(&j,1,1);
	  pvm_pkint(&cur_val,1,1);
	  if (ii<ntax-4) pack_align(temp_align);
	  else {
		if (temp_align2) pack_align(temp_align2);
		else pack_align(temp_align);
		}
	  pvm_send(values->tids[0],PVM_NEW_BUILD_ALIGNMENTS_ONLY_SEND);
		if (temp_align) temp_align=dump_align(temp_align);
		if (temp_align2) temp_align2=dump_align(temp_align2);
		if (t1) t1=dump_align(t1);
	  if (t2) t2=dump_align(t2);
	  if (t3) t3=dump_align(t3);
	}
	else if (tag==PVM_NEW_BUILD_OPTIMIZE) {  /*optimize with best*/
		pvm_upkint(&i, 1, 1);
	  for (k=0;k<ntax-1;k++) pvm_upkint(nodes[k],2,1);
	  pvm_upkint(anc,(2*ntax)-1,1);
	  /*reoptimize internal*/
		check=all_diagnose_tree_local2(a,nodes,ntax,i+4,values); /*this gets the down pass need all they could be packed if necessary*/
		/*get up_pass alignments*/
		get_up_pass2(a,nodes,ntax,values,a_up,i+4,anc);
		}
	}

/*deallocate and free*/
end_shit:;
free(anc);
for (i=0;i<2*ntax-1;i++) {
	if (a_up[i][0]) a_up[i][0]=dump_align(a_up[i][0]);
	if (a_up[i][1]) a_up[i][1]=dump_align(a_up[i][1]);
	free(a_up[i]);
	}
free(a_up);
if (temp_align) temp_align=dump_align(temp_align);
if (temp_align2) temp_align2=dump_align(temp_align2);
if (cur_align) cur_align=dump_align(cur_align);
if (cur_align2) cur_align2=dump_align(cur_align2);
for (i=0;i<ntax-1;i++) free(nodes[i]);
free(nodes);
values->phylo_score=score_holder;
values->in_some=0;
}

void do_new_build_parallel2(a,ntax,values)
alignment **a;
int ntax;
parameters *values;
{
alignment *temp_align,*temp_align2;
int d1,d2,d3;
alignment *t1,*t2,*t3;
int j,ii,cur_val,best,score_holder;

temp_align=NULL;
temp_align2=NULL;
t2=NULL;
t3=NULL;

/*set so no automatic cladogram*/
score_holder=values->phylo_score;
values->phylo_score=0;
values->in_some=1;

pvm_upkint(&score_holder,1,1);
pvm_upkint(&j,1,1);
pvm_upkint(&d1,1,1);
pvm_upkint(&d2,1,1);
pvm_upkint(&d3,1,1);
pvm_upkint(&ii,1,1);
t2=unpack_align_and_score(values);
if (d1>0) t3=unpack_align_and_score(values);
if (!values->new_optimization) {
	temp_align=nw(a[ii+2],t2,values);/*would be ii+3 but the adjustment for jack-up*/
		if (d1>0) {
		if (temp_align2) temp_align2=dump_align(temp_align2);
		temp_align2=nw(temp_align,t3,values);
		if (score_holder>0) {
			if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
			else {
				values->phylo_score=score_holder;
				cur_val=temp_align2->score=cladogram(temp_align2,values);
				values->phylo_score=0;
				if (values->previous) values->previous=dump_align(values->previous);
				values->previous=make_align(temp_align2);
				}
			}
		cur_val=temp_align2->score;
		}
	else {
		if (score_holder>0) {
			if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
			else {
				values->phylo_score=score_holder;
				cur_val=temp_align->score=cladogram(temp_align,values);
				values->phylo_score=0;
				if (values->previous) values->previous=dump_align(values->previous);
				values->previous=make_align(temp_align);
				}
			}
		cur_val=temp_align->score;
		}
	}
else {
	if (d1>0) {
		values->in_optimize=1;
		temp_align2=nw(t2,t3,values);
		values->in_optimize=0;
		temp_align2=newer_make_ambig(temp_align2,values);
		temp_align=nw(temp_align2,a[ii+2],values);
		}
	else temp_align=nw(t2,a[ii+2],values);
	cur_val=temp_align->score;
	}
pvm_initsend( PvmDataDefault );
pvm_pkint(&j,1,1);
pvm_pkint(&cur_val,1,1);

if (!values->new_optimization) {
		if (ii<ntax-4)  pack_align(temp_align);
	else {
		if (temp_align2) pack_align(temp_align2);
		else pack_align(temp_align);
		}
	}
pvm_send(values->tids[0],PVM_NEW_BUILD_ALIGNMENTS_ONLY_SEND);
if (temp_align) temp_align=dump_align(temp_align);
if (temp_align2) temp_align2=dump_align(temp_align2);
if (t2) t2=dump_align(t2);
if (t3) t3=dump_align(t3);

values->phylo_score=score_holder;
values->in_some=0;
}

void do_new_parallel_swap(a,ntax,values)
alignment **a;
parameters *values;
int ntax;
{
int i,j,k,l,m,num_clipped,check,tot_best,best,d1,d2;
int l_clip, l_rest,cur_best,found_better,node_removed;
int **cur_nodes, *blocked, *is_it_or_desc;
int clipped_base,***reroot_array,cur_num_trees;
int old_length, tree_counter;
int ***best_nodes_new;
int new_start,all_others;
int where_it_was_anc,where_it_was_desc;
int whole_length;
int initial_length,d3,d4;
alignment *thing_to_add,*temp_align,*temp_align2;
int where_it_was_branch;
int is_unique,*anc,mult;
alignment ***a_up,**best_aligns;
int score_holder,**temp_nodes;
int will_be_OK,**cur_groups,ngroups2;
int ii,jj;
int ***corres,all_rest;
alignment **a_down;
int ***best_nodes_old,x,y;
alignment **a_down_buf, **a_buf, ***a_up_buf;
int branch_holder,down_node,up_node;
alignment *up_temp, **check_buf;
char *collapse_name;

best_nodes_old=NULL;
up_temp=NULL;
check_buf=NULL;
collapse_name=NULL;
values->support=NULL;
/*allocations*/
if (values->new_optimization) {
	a_down=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) a_down);
	for (i=0;i<2*ntax-1;i++) a_down[i]=NULL;
	corres=(int ***)malloc((2*(ntax))*sizeof(int **));
	assert((int) corres);
	for (i=0;i<2*ntax;i++) corres[i]=NULL;
	for (i=0;i<ntax-1;i++) a_down[i]=make_align(a[i]);
	check_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) check_buf);
	for (i=0;i<2*ntax-1;i++) check_buf[i]=NULL;
	for (i=0;i<ntax-1;i++) check_buf[i]=make_align(a[i]);
	}
a_down_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) a_down_buf);
for (i=0;i<2*ntax-1;i++) a_down_buf[i]=NULL;
for (i=0;i<ntax-1;i++) a_down_buf[i]=make_align(a[i]);

temp_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int)temp_nodes);
for (i=0;i<ntax-1;i++) {
	temp_nodes[i]=(int *)malloc(2*sizeof(int));
	assert((int) temp_nodes[i]);
	}

best_nodes_new=(int ***)malloc(values->keep_aligns*sizeof(int **));
assert((int) best_nodes_new);
for (i=0;i<values->keep_aligns;i++) {
	best_nodes_new[i]=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int) best_nodes_new[i]);
	for (j=0;j<ntax-1;j++) {
		best_nodes_new[i][j]=(int *)malloc(2*sizeof(int));
		assert((int) best_nodes_new[i][j]);
		}
	}
cur_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) cur_nodes);
for (k=0;k<ntax-1;k++) {
	cur_nodes[k]=(int *)malloc(2*sizeof(int));
	assert((int) cur_nodes[k]);
	}
blocked=(int *)malloc((2*(ntax-1))*sizeof(int));
assert((int) blocked);
is_it_or_desc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) is_it_or_desc);
if (values->atbr) {
	reroot_array=(int ***)malloc((2*ntax-3)*sizeof(int **));
	assert((int) reroot_array);
	for (i=0;i<2*ntax-3;i++) {
		reroot_array[i]=(int **)malloc((ntax-1)*sizeof(int *));
		assert((int) reroot_array[i]);
		for (j=0;j<ntax-1;j++) {
			reroot_array[i][j]=(int *)malloc(2*sizeof(int));
			assert((int) reroot_array[i][j]);
			}
		}
	}
if (values->tbr_align_shortcut) {
	a_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
	assert((int) a_buf);
	for (i=0;i<2*ntax-1;i++) a_buf[i]=NULL;
	for (i=0;i<ntax-1;i++) a_buf[i]=make_align(a[i]);
	a_up_buf=(alignment ***)malloc(((2*ntax)-1)*sizeof(alignment **));
	assert((int) a_up_buf);
	for (i=0;i<2*ntax-1;i++) {
		a_up_buf[i]=(alignment **)malloc(2*sizeof(alignment *));
		assert((int) a_up_buf[i]);
		a_up_buf[i][0]=NULL;
		a_up_buf[i][1]=NULL;
		}
	}

anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);
a_up=(alignment ***)malloc(((2*ntax)-1)*sizeof(alignment **));
	assert((int) a_up);
	for (i=0;i<2*ntax-1;i++) {
		a_up[i]=(alignment **)malloc(2*sizeof(alignment *));
		assert((int) a_up[i]);
		a_up[i][0]=NULL;
		a_up[i][1]=NULL;
		}
best_aligns=(alignment **)malloc(values->keep_aligns*sizeof(alignment *));
assert((int) best_aligns);
for (i=0;i<values->keep_aligns;i++) best_aligns[i]=NULL;
thing_to_add=NULL;
temp_align=NULL;
temp_align2=NULL;

will_be_OK=1;
/*receive the info*/
pvm_upkint(&num_clipped,1,1);
for (i=0;i<ntax-1;i++) pvm_upkint(cur_nodes[i],2,1);
pvm_upkint(anc,((2*ntax)-1),1);
pvm_upkint(&d1,1,1);
pvm_upkint(&d2,1,1);
pvm_upkint(&where_it_was_desc,1,1);
pvm_upkint(&where_it_was_anc,1,1);
pvm_upkint(&where_it_was_branch,1,1);
pvm_upkint(&score_holder,1,1);
pvm_upkint(&old_length,1,1);
pvm_upkint(blocked,(2*(ntax-1)),1);
pvm_upkint(&node_removed,1,1);
pvm_upkint(&mult,1,1);
pvm_upkint(is_it_or_desc,((2*ntax)-1),1);
for (ii=ntax+1;ii<2*ntax-1;ii++) a_down_buf[ii]=unpack_align_and_score(values);
if (values->tbr_align_shortcut) for (ii=ntax+1;ii<2*ntax-1;ii++) a_buf[ii]=unpack_align_and_score(values);
tree_counter=0;
tot_best=old_length;
				if ((num_clipped < 3) || (!values->atbr)) { /*if no alternative roots or sbr*/
					if (!values->new_optimization) {
						check=all_diagnose_tree_local2(a,cur_nodes,ntax,ntax,values);
						/*check=all_diagnose_tree_local3(a,cur_nodes,ntax,ntax,values,is_it_or_desc);
						for (k=ntax+1;k<2*ntax-1;k++) if (is_it_or_desc[k]) {
							if (a[k]) a[k]=dump_align(a[k]);
							a[k]=make_align(a_down_buf[k]);
							}*/
						get_up_pass2(a,cur_nodes,ntax,values,a_up,ntax,anc);
						}
					else {
						check=all_diagnose_tree_local_new_opt2(a_down,cur_nodes,ntax,ntax,values,node_removed,corres,is_it_or_desc);
						/*get up_pass alignments*/
						get_up_pass_new_opt2(a,a_down,a_up,cur_nodes,ntax,values,ntax,anc,corres,node_removed,is_it_or_desc);
						}
					/* get "best" set sequences[cur_nodes[d1][d2]] as XUY */
					d3=cur_nodes[d1][d2];
					if (d3<ntax) d3-=1;
					d4=where_it_was_desc;
					if (d4<ntax) d4-=1;
					if (thing_to_add) thing_to_add=dump_align(thing_to_add);
					thing_to_add=make_align(a_down_buf[d3]);
					if (!values->new_optimization) {
						if (temp_align) temp_align=dump_align(temp_align);
						if (where_it_was_anc>ntax) temp_align=nw(a_up[where_it_was_anc][where_it_was_branch],a[d4],values);
						else temp_align=make_align(a[d4]);
						if (score_holder>0) {
							values->phylo_score=score_holder;
							temp_align->score=cladogram(temp_align,values);
							values->phylo_score=0;
							}
						/*redo score with phyloshit*/
						best=old_length-temp_align->score;
						}
					else {
						all_rest=thing_to_add->score + a_down[cur_nodes[0][1]]->score;
						best=old_length - all_rest;
						}
					if ((best>0) || ((best==0) && (mult))) {/*to make sure that that sequence matters ie adds length*/
						/*loop through places to put back*/
						for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
							if (values->groups) {
								for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=cur_nodes[ii][jj];
								check=temp_nodes[j/2][j%2];
								temp_nodes[j/2][j%2]=node_removed;
								temp_nodes[d1][!d2]=check;
								cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
								will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
								for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
								free(cur_groups);
								}
							if (will_be_OK) {
								/* get_cost*/
								if (temp_align) temp_align=dump_align(temp_align);
								if (temp_align2) temp_align2=dump_align(temp_align2);
								d3=cur_nodes[j/2][j%2];
								if (d3<ntax) d3-=1;
								if (!values->new_optimization) {
									temp_align=nw(thing_to_add,a[d3],values);
									if (j>1) {
										/*fprintf(stderr,"{");*/
										d4=(j/2)+ntax;
										temp_align2=nw(temp_align,a_up[d4][j%2],values);
										if (score_holder>0) {
											if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
											else {
												values->phylo_score=score_holder;
												temp_align2->score=cladogram(temp_align2,values);
												values->phylo_score=0;
												if (values->previous) values->previous=dump_align(values->previous);
												values->previous=make_align(temp_align2);
												}
											}
										cur_best=temp_align2->score;
										/*fprintf(stderr,"}");*/
										}
									else {
										/*fprintf(stderr,"(");*/
										if (score_holder>0) {
											if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
											else {                                                                                  values->phylo_score=score_holder;
												temp_align->score=cladogram(temp_align,values);
												values->phylo_score=0;
												if (values->previous) values->previous=dump_align(values->previous);
												values->previous=make_align(temp_align);
												}
											}
										cur_best=temp_align->score;
										/*fprintf(stderr,")");*/
										}
									}
								else {
									d4=(j/2)+ntax;
									if (j>1) {
										values->in_optimize=1;
										temp_align2=nw(a[d3],a[d4],values);
										values->in_optimize=0;
										temp_align2=newer_make_ambig(temp_align2,values);
										temp_align=nw(temp_align2,thing_to_add,values);
										}
									else temp_align=nw(thing_to_add,a[d3],values);
									cur_best=temp_align->score+all_rest;
									}
								if (values->new_optimization) {
									if ((cur_best < tot_best) || ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns))) {
										for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=cur_nodes[ii][jj];
										check=temp_nodes[j/2][j%2];
										temp_nodes[j/2][j%2]=node_removed;
										temp_nodes[d1][!d2]=check;
										check=cur_best;
										cur_best=all_diagnose_tree_local_new_opt(check_buf,temp_nodes,ntax,ntax,values,0,corres);
										}
									}
								if (cur_best < tot_best) {
									tot_best=cur_best;
									found_better=1;
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[0][k][0]=cur_nodes[k][0];
										best_nodes_new[0][k][1]=cur_nodes[k][1];
										}
									check=best_nodes_new[0][j/2][j%2];
									best_nodes_new[0][j/2][j%2]=node_removed;
									best_nodes_new[0][d1][!d2]=check;
									/*add best_aligns stuf free first*/
									for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
									if (temp_align2) best_aligns[0]=make_align(temp_align2);
									else best_aligns[0]=make_align(temp_align);
									tree_counter=1;
									}
								else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
									/*fprintf(stderr,"F");*/
									for (k=0;k<ntax-1;k++) {
										best_nodes_new[tree_counter][k][0]=cur_nodes[k][0];
										best_nodes_new[tree_counter][k][1]=cur_nodes[k][1];
										}
									check=best_nodes_new[tree_counter][j/2][j%2];
									best_nodes_new[tree_counter][j/2][j%2]=node_removed;
									best_nodes_new[tree_counter][d1][!d2]=check;
									/*chec to see if novel*/
									is_unique=1;
									if (!values->new_optimization) {
										if (temp_align2) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align2,best_aligns[k]);
										else for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
										}
									else {
										if (1) /*(values->collapse==2) will get sorted out by master*/ is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
										else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*0,values,collapse_name);
										}
									if (is_unique) {
										found_better=1;
										if (temp_align2) best_aligns[tree_counter]=make_align(temp_align2);
										else best_aligns[tree_counter]=make_align(temp_align);
										++tree_counter;
										}
									/*fprintf(stderr,"O");*/
									}
								} /*SBR loop--regrafted*/
							}/*will it be OK*/
						} /*BEST split size ok*/
					}/*too few clipped or SBR*/
				else { /*TBR with at least three taxa*/
					check=get_root_array(reroot_array,num_clipped,node_removed,cur_nodes,d1,d2,ntax,is_it_or_desc,anc);
					if (!values->new_optimization) {
						check=all_diagnose_tree_local2(a,cur_nodes,ntax,ntax,values);
						/*check=all_diagnose_tree_local3(a,cur_nodes,ntax,ntax,values,is_it_or_desc);
						for (k=ntax+1;k<2*ntax-1;k++) if (is_it_or_desc[k]) {
							if (a[k]) a[k]=dump_align(a[k]);
							a[k]=make_align(a_down_buf[k]);
							}*/
						get_up_pass2(a,cur_nodes,ntax,values,a_up,ntax,anc);
						}
					else {
						/*only down pass if ! in picked*/
						check=all_diagnose_tree_local_new_opt2(a_down,cur_nodes,ntax,ntax,values,node_removed,corres,is_it_or_desc);
						/*get up_pass alignments*/
						get_up_pass_new_opt2(a,a_down,a_up,cur_nodes,ntax,values,ntax,anc,corres,node_removed,is_it_or_desc);
						}
					for (l=0;l<2*num_clipped-3;l++) {
						/*only down pass if in picked*/
						d3=reroot_array[l][d1][d2];
						if (d3<ntax) d3-=1;
						d4=where_it_was_desc;
						if (d4<ntax) d4-=1;
						if (thing_to_add) thing_to_add=dump_align(thing_to_add);
						/*if same rooting as original use the buffer*/
						if ((cur_nodes[cur_nodes[d1][d2]-ntax][0]==reroot_array[l][reroot_array[l][d1][d2]-ntax][0]) &&(cur_nodes[cur_nodes[d1][d2]-ntax][1]==reroot_array[l][reroot_array[l][d1][d2]-ntax][1])) thing_to_add=make_align(a_down_buf[d3]);
						else if (values->tbr_align_shortcut) { /*X U Y*/
							if (values->new_optimization) {
								x=reroot_array[l][reroot_array[l][d1][d2]-ntax][0];
								y=reroot_array[l][reroot_array[l][d1][d2]-ntax][1];
								if (x<ntax) x-=1;
								if (y<ntax) y-=1;
								thing_to_add=nw(a_buf[x],a_buf[y],values);
								}
							else { /*currently disabled*/
								/*get X and Y*/
								x=reroot_array[l][reroot_array[l][d1][d2]-ntax][0];
								y=reroot_array[l][reroot_array[l][d1][d2]-ntax][1];
								if (x<ntax) {
									down_node=x;
									up_node=y;
																			}
								else if (y<ntax) {
									down_node=y;
									up_node=x;
									}
								else {
									if (cur_nodes[y-ntax][0]==x) { down_node=x; up_node=y;}
									else if (cur_nodes[y-ntax][1]==x) { down_node=x; up_node=y; }
									else if (cur_nodes[x-ntax][0]==y) { down_node=y; up_node=x; }
									else { down_node=y; up_node=x; }

									}
								/*get the up*/
								if (up_temp) up_temp=dump_align(up_temp);
								up_temp=get_single_up(a_down_buf[down_node],a_down_buf[cur_nodes[d1][d2]],values);
								/*make it*/
								if (down_node<ntax) x-=1;
								thing_to_add=nw(a_down_buf[down_node],up_temp,values);
								if (up_temp) up_temp=dump_align(up_temp);

								}
							}
						else {
							if (!values->new_optimization) check=all_diagnose_tree_local3(a,reroot_array[l],ntax,ntax,values,is_it_or_desc); /*,is_it_or_desc); get down passes*/
							else check=all_diagnose_tree_local_new_opt3(a_down,reroot_array[l],ntax,ntax,values,node_removed,corres,is_it_or_desc);
							if (!values->new_optimization) thing_to_add=make_align(a[d3]);
							else thing_to_add=make_align(a_down[d3]);
							}
						if (!values->new_optimization) {
							if (temp_align) temp_align=dump_align(temp_align);
							if (where_it_was_anc>ntax) temp_align=nw(a_up[where_it_was_anc][where_it_was_branch],a[d4],values);
							else temp_align=make_align(a[d4]);
							if (score_holder>0) {
								values->phylo_score=score_holder;
								temp_align->score=cladogram(temp_align,values);
								values->phylo_score=0;
								}
							/*redo score with phyloshit*/
							best=old_length-temp_align->score;
							}
						else {
							all_rest=thing_to_add->score + a_down[reroot_array[l][0][1]]->score;
							best=old_length - all_rest;
							}
						if ((best>0) || ((best==0) && (mult))) {/*to make sure that that sequence matters ie adds length*/
							/*loop through places to put back*/
							for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
								if (values->groups) {
									for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=reroot_array[l][ii][jj];
									check=temp_nodes[j/2][j%2];
									temp_nodes[j/2][j%2]=node_removed;
									temp_nodes[d1][!d2]=check;
									cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
									will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
									for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
									free(cur_groups);
									}
								if (will_be_OK) {
									/* get_cost*/
									if (temp_align) temp_align=dump_align(temp_align);
									if (temp_align2) temp_align2=dump_align(temp_align2);
									d3=reroot_array[l][j/2][j%2];
									if (d3<ntax) d3-=1;
									if (!values->new_optimization) {
										temp_align=nw(thing_to_add,a[d3],values);
										if (j>1) {
											/*fprintf(stderr,"{");*/
											d4=(j/2)+ntax;
											temp_align2=nw(temp_align,a_up[d4][j%2],values);
											if (score_holder>0) {
												if (!compare_aligns(temp_align2,values->previous)) temp_align2->score=values->previous->score;
												else {
													values->phylo_score=score_holder;
													temp_align2->score=cladogram(temp_align2,values);
													values->phylo_score=0;
													if (values->previous) values->previous=dump_align(values->previous);
													values->previous=make_align(temp_align2);
													}
												}
											cur_best=temp_align2->score;
											/*fprintf(stderr,"}");*/
											}
										else {
											/*fprintf(stderr,"(");*/
											if (score_holder>0) {
												if (!compare_aligns(temp_align,values->previous)) temp_align->score=values->previous->score;
												else {
													values->phylo_score=score_holder;
													temp_align->score=cladogram(temp_align,values);
													values->phylo_score=0;
													if (values->previous) values->previous=dump_align(values->previous);
													values->previous=make_align(temp_align);
													}
												}
											cur_best=temp_align->score;
											/*fprintf(stderr,")");*/
											}
										}
									else {
										d4=(j/2)+ntax;
										if (j>1) {
											values->in_optimize=1;
											temp_align2=nw(a[d3],a[d4],values);
											values->in_optimize=0;
											temp_align2=newer_make_ambig(temp_align2,values);
											temp_align=nw(temp_align2,thing_to_add,values);
											}
										else temp_align=nw(thing_to_add,a[d3],values);
										cur_best=temp_align->score+all_rest;
										}
									if (values->new_optimization) {
										if ((cur_best < tot_best) || ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns))) {
											for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=reroot_array[l][ii][jj];
											check=temp_nodes[j/2][j%2];
											temp_nodes[j/2][j%2]=node_removed;
											temp_nodes[d1][!d2]=check;
											check=cur_best;
											cur_best=all_diagnose_tree_local_new_opt(check_buf,temp_nodes,ntax,ntax,values,0,corres);
											}
										}
									if (cur_best < tot_best) {
										tot_best=cur_best;
										found_better=1;
										for (k=0;k<ntax-1;k++) {
											best_nodes_new[0][k][0]=reroot_array[l][k][0];
											best_nodes_new[0][k][1]=reroot_array[l][k][1];
											}
										check=best_nodes_new[0][j/2][j%2];
										best_nodes_new[0][j/2][j%2]=node_removed;
										best_nodes_new[0][d1][!d2]=check;
										/*add best_aligns stuf free first*/
										for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
										if (temp_align2) best_aligns[0]=make_align(temp_align2);
										else best_aligns[0]=make_align(temp_align);
										tree_counter=1;
										}
									else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
										for (k=0;k<ntax-1;k++) {
											best_nodes_new[tree_counter][k][0]=reroot_array[l][k][0];
											best_nodes_new[tree_counter][k][1]=reroot_array[l][k][1];
											}
										check=best_nodes_new[tree_counter][j/2][j%2];
										best_nodes_new[tree_counter][j/2][j%2]=node_removed;
										best_nodes_new[tree_counter][d1][!d2]=check;
										/*chec to see if novel*/
										is_unique=1;
										if (!values->new_optimization) {
											if (temp_align2) for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align2,best_aligns[k]);
											else for (k=0;k<tree_counter;k++) is_unique*=compare_aligns(temp_align,best_aligns[k]);
											}
										else {
											if (1) /*(values->collapse==2) will get sorted out by master*/ is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
											else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*0,values,collapse_name);
											}
										if (is_unique) {
											found_better=1;
											if (temp_align2) best_aligns[tree_counter]=make_align(temp_align2);
											else best_aligns[tree_counter]=make_align(temp_align);
											++tree_counter;
											}
										}
									} /*SBR loop--regrafted*/
								}
							} /*BEST split size ok*/
						}
					} /*reroot TBR/SBR choice*/

/*pack and send the results*/
pvm_initsend( PvmDataDefault );
pvm_pkint(&tot_best,1,1); /*pack cost*/
pvm_pkint(&tree_counter,1,1);/*pack #solns*/
for (i=0;i<tree_counter;i++) { /*make sure OK when no solutions found*/
	/*loop nodes then align*/
	for (j=0;j<ntax-1;j++) pvm_pkint(best_nodes_new[i][j],2,1);
	best_aligns[i]->score=tot_best;
	pack_align(best_aligns[i]);
	}
pvm_send(values->tids[0],PVM_NEW_SWAP_DONE);

/*free everything*/
if (values->new_optimization) {
	values->support=NULL;
	for (i=0;i<2*ntax-1;i++) if (a_down[i]) a_down[i]=dump_align(a_down[i]);
	free(a_down);
	for (i=0;i<2*ntax;i++) {
		if (corres[i]) {
			if (corres[i][0]) free(corres[i][0]);
			if (corres[i][1]) free(corres[i][1]);
			free(corres[i]);
			}
		}
	free(corres);
	for (i=0;i<2*ntax-1;i++) if (check_buf[i]) check_buf[i]=dump_align(check_buf[i]);
	free(check_buf);
	}
for (i=0;i<ntax-1;i++) {
	free(cur_nodes[i]);
	}
free(cur_nodes);
free(is_it_or_desc);
free(blocked);
free(anc);
if (values->atbr) {
	for (i=0;i<2*ntax-3;i++) {
		for (j=0;j<ntax-1;j++) free(reroot_array[i][j]);
		free(reroot_array[i]);
		}
	free(reroot_array);
	}

for (i=0;i<values->keep_aligns;i++) {
	for (j=0;j<ntax-1;j++) {
		free(best_nodes_new[i][j]);
		}
	free(best_nodes_new[i]);
	}
free(best_nodes_new);
if (thing_to_add) thing_to_add=dump_align(thing_to_add);
if (temp_align) temp_align=dump_align(temp_align);
if (temp_align2) temp_align2=dump_align(temp_align2);
for (i=0;i<ntax-1;i++) free(temp_nodes[i]);
free(temp_nodes);

for (i=0;i<2*ntax-1;i++) if (a_down_buf[i]) a_down_buf[i]=dump_align(a_down_buf[i]);
free(a_down_buf);
if (values->tbr_align_shortcut) {
	for (i=0;i<2*ntax-1;i++) if (a_buf[i]) a_buf[i]=dump_align(a_buf[i]);
	free(a_buf);
	if (a_up_buf) {
		for (i=0;i<2*ntax-1;i++) {
			if (a_up_buf[i][0]) a_up_buf[i][0]=dump_align(a_up_buf[i][0]);
			if (a_up_buf[i][1]) a_up_buf[i][1]=dump_align(a_up_buf[i][1]);
			free(a_up_buf[i]);
			}
		free(a_up_buf);
		}
	}
if (a_up) {
	for (i=0;i<2*ntax-1;i++) {
		if (a_up[i][0]) a_up[i][0]=dump_align(a_up[i][0]);
		if (a_up[i][1]) a_up[i][1]=dump_align(a_up[i][1]);
		free(a_up[i]);
		}
	free(a_up);
	}
for (i=0;i<values->keep_aligns;i++) if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
free(best_aligns);
}


void redo_it_all(a,only_single_machine,values)
parameters *values;
alignment **a;
int only_single_machine;
{
int numt,info,i,time_holder,number_holder,time_holder2;
struct hostinfo *hostp;

/*first get times and num_trees*/
time_holder2=0;
pvm_initsend( PvmDataDefault );
pvm_pkint(&(values->phylo_time), 1, 1);
pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_REQUEST_TIME);
for (i=0;i<(values->num_hosts-1);i++) {
	pvm_recv(values->tids[i+1],PVM_SEND_TIME);
	pvm_upkint(&time_holder,1,1);
	time_holder2+=time_holder;
	}
values->phylo_time+=(time_holder2/(values->num_hosts-1));
pvm_initsend( PvmDataDefault );/*should this be modified for multiplatform?*/
pvm_pkint(&(values->number), 1, 1);
pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_REQUEST_NUM_TREES);
for (i=0;i<(values->num_hosts-1);i++) {
	pvm_recv(values->tids[i+1],PVM_SEND_NUM_TREES);
	pvm_upkint(&number_holder,1,1);
	values->number+=number_holder;
	}
/*kill processes*/
pvm_initsend( PvmDataDefault );
pvm_pkint(values->tids, values->num_hosts, 1);
pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_ALL_MUST_DIE);

/*spawn Processes*/
if (!only_single_machine) numt=pvm_spawn(EXE_FILE_NAME, (char**)0, PvmTaskHost|PvmHostCompl, ".", values->num_hosts-1, &(values->tids[1]));
else numt=pvm_spawn(EXE_FILE_NAME,(char **)0, 0,"",values->num_hosts-1,&(values->tids[1]));
if (numt!=(values->num_hosts-1)) {
	fprintf(stderr,"Problems spawning tasks only %d of %d actually spawned.\n",numt,values->num_hosts-1);
	for (i=0;i<values->num_hosts-1;i++) fprintf(stderr,"   host %s task id [%d]=%d\n",hostp->hi_name,i,values->tids[i]);
	fprintf(stderr,"Reset PVM.\n");
	pvm_exit();
	exit(-1);
	}

/*initialize*/
pvm_initsend( PvmDataDefault );
pvm_pkint(values->tids, values->num_hosts, 1);
info = pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_SEND_TIDS);
do_initialize_thang(values,a);
}

alignment *get_single_up(down,base,values)
alignment *base, *down;
parameters *values;
{
alignment *up;
int i,j,k,l;
int n_seqs,length;
char **names, **bases;
int found,seq_counter,*all_gaps;
int new_length,l_counter,d1;

	/*alloacte holders*/
	n_seqs=base->n_seqs;
	length=base->length;
	names=(char **)malloc(n_seqs*sizeof(char *));
	assert((int)names);
	bases=(char **)malloc(n_seqs*sizeof(char *));
	assert((int)bases);
	for (i=0;i<n_seqs;i++) {
		names[i]=(char *)malloc(100*sizeof(char));
		assert((int)names[i]);
		bases[i]=(char *)malloc((1+length)*sizeof(char));
		assert((int)bases[i]);
	}
	all_gaps=(int *)malloc(length*sizeof(int));
	assert((int)all_gaps);
		seq_counter=0;
		for (k=0;k<n_seqs;k++) {
			found=0;
			for (l=0;l<down->n_seqs;l++) if (!strcmp(base->taxon_name[k],down->taxon_name[l])) found=1;
			if (!found) { /*want these*/
				strcpy(names[seq_counter],base->taxon_name[k]);
				strcpy(bases[seq_counter],base->s[k]);
				++seq_counter;
			}
		}/*k*/
		/*check seqs copied are correct in number*/
		if (seq_counter==(n_seqs-down->n_seqs)) { /*else in sbr/tbr non-needed up-passes*/
			/*scan for all gaps in remainder*/
			for (k=0;k<length;k++) {
				all_gaps[k]=1;
				for (l=0;l<seq_counter;l++) if (bases[l][k]!='-') all_gaps[k]=0;
			}
			new_length=0;
			for (k=0;k<length;k++) if (!all_gaps[k]) ++new_length;

			/*allocate and copy to up*/
			up=(alignment *)malloc(sizeof(alignment));
			assert((int)up);
			up->n_seqs=seq_counter;
			up->score=0;
			up->length=new_length;
			up->type_weight=0;
			up->name=(char *)malloc(8*sizeof(char));
			assert((int)up->name);
			up->name[0]='U';
			up->name[1]='p';
			up->name[2]='-';
			up->name[3]='P';
			up->name[4]='a';
			up->name[5]='s';
			up->name[6]='s';
			up->name[7]='\0';
			up->s=(char **)malloc((up->n_seqs+up->type_weight)*sizeof(char *));
			assert((int)up->s);
			for (k=0;k<(up->n_seqs+up->type_weight);k++) {
				up->s[k]=(char *)malloc((up->length+1)*sizeof(char));
				assert((int)up->s[k]);
				l_counter=0;
				for (l=0;l<length;l++) if (!all_gaps[l]) up->s[k][l_counter++]=bases[k][l];
				up->s[k][new_length]='\0';
			}

			up->taxon_name=(char **)malloc((up->n_seqs)*sizeof(char *));
			assert((int)up->taxon_name);
			for (k=0;k<up->n_seqs;k++) {
				up->taxon_name[k]=(char *)malloc((strlen(names[k])+1)*sizeof(char));
				assert((int)up->taxon_name[k]);
				up->taxon_name[k]=(char *)strcpy(up->taxon_name[k],names[k]);
			}
		}


	/*freeing*/
	for (i=0;i<n_seqs;i++) {
		free(names[i]);
		free(bases[i]);
	}
	free(names);
	free(bases);
	free(all_gaps);


return up;
}

int do_dave_build(a,temp_nodes,ntax,i,values,order,value,a_best,best_k)
alignment **a, **a_best;
int **temp_nodes;
int ntax,i,*order,value;
int *best_k;
parameters *values;
{
int ii,k,kk,num_to_do;
int will_be_OK,done,same;
int **cur_groups,ngroups2;
int cur_val,temp_val;
alignment *temp_align;


will_be_OK=1;
if (values->groups) {
	cur_groups=get_cur_groups(temp_nodes,ntax,i+4,&ngroups2);
	will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,i+4-1);
	for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
	free(cur_groups);
	}
if (will_be_OK) {
	/*set up nodes to do and order*/
	for (k=0;k<ntax;k++) order[k]=(-1);
	order[0]=i+2;
	num_to_do=1;
	done=0;
	while (!done) {
		done=1;
		for (k=1;k<i+3;k++) {
			same=0;
			for (kk=0;kk<num_to_do;kk++) if (order[kk]==k) same=1;
			if (!same) {
				for (kk=0;kk<num_to_do;kk++) {
					if (temp_nodes[k][0]==(order[kk]+ntax)) {order[num_to_do++]=k;done=0;}
					else if (temp_nodes[k][1]==(order[kk]+ntax)) {order[num_to_do++]=k;done=0;}
						}
					}
			}
		}
	/*do_alignment*/
	/*could put loop in here for best adddiotn sequence*/
	*best_k=i+2;
	if (values->best_order) {
		temp_val=value;
	    	for (k=i+2;k<ntax-1;k++) {
 	   		if (k!=(i+2)) {
	    			temp_align=make_align(a[i+2]);
	        		a[i+2]=dump_align(a[i+2]);
	        		a[i+2]=make_align(a[k]);
	        		}
	        	cur_val=all_diagnose_tree_cut_off(a,temp_nodes,ntax,i+4,values,order,num_to_do,temp_val,a_best);
	        	/*fprintf(stderr,"O%d-%d-%d- ",k,cur_val,temp_val);*/
	        	if (cur_val < temp_val) {
	        	    temp_val=cur_val;
	        	    *best_k=k;
		             }
	    	    	if (k!=(i+2)) {
    				a[i+2]=dump_align(a[i+2]);
       				a[i+2]=make_align(temp_align);
 				temp_align=dump_align(temp_align);
                    	}
                     }
                     /*fprintf(stderr,"PVM%d ",(*best_k));*/
		}
	else temp_val=all_diagnose_tree_cut_off(a,temp_nodes,ntax,i+4,values,order,num_to_do,value,a_best);
	
	return temp_val;
	}
else return HUGE_COST;
}

void do_new_parallel_swap_dave(a,ntax,values)
alignment **a;
parameters *values;
int ntax;
{
int i,j,k,l,m,num_clipped,check,tot_best,best,d1,d2;
int l_clip, l_rest,cur_best,found_better,node_removed;
int **cur_nodes, *blocked, *is_it_or_desc;
int clipped_base,***reroot_array,cur_num_trees;
int old_length, tree_counter;
int ***best_nodes_new;
int new_start,all_others;
int where_it_was_anc,where_it_was_desc;
int initial_length,d3,d4;
int where_it_was_branch;
int is_unique,*anc,mult;
alignment **best_aligns,**check_buf;
int score_holder,**temp_nodes;
int will_be_OK,**cur_groups,ngroups2;
int ii,jj;
int all_rest;
int ***best_nodes_old,x,y;
int branch_holder,down_node,up_node;
char *collapse_name;
int *order,kk,num_to_do,same,done;
int collapse_holder;

collapse_holder=values->collapse;
if (mult) values->collapse=2;
order=NULL;
order=(int *)malloc(ntax*sizeof(int));
assert((int)order);
best_nodes_old=NULL;
check_buf=NULL;
collapse_name=NULL;
values->support=NULL;
/*allocations*/
check_buf=(alignment **)malloc((2*ntax-1)*sizeof(alignment *));
assert((int) check_buf);
for (i=0;i<2*ntax-1;i++) check_buf[i]=NULL;
for (i=0;i<ntax-1;i++) check_buf[i]=make_align(a[i]);


temp_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int)temp_nodes);
for (i=0;i<ntax-1;i++) {
	temp_nodes[i]=(int *)malloc(2*sizeof(int));
	assert((int) temp_nodes[i]);
	}

best_nodes_new=(int ***)malloc(values->keep_aligns*sizeof(int **));
assert((int) best_nodes_new);
for (i=0;i<values->keep_aligns;i++) {
	best_nodes_new[i]=(int **)malloc((ntax-1)*sizeof(int *));
	assert((int) best_nodes_new[i]);
	for (j=0;j<ntax-1;j++) {
		best_nodes_new[i][j]=(int *)malloc(2*sizeof(int));
		assert((int) best_nodes_new[i][j]);
		}
	}
cur_nodes=(int **)malloc((ntax-1)*sizeof(int *));
assert((int) cur_nodes);
for (k=0;k<ntax-1;k++) {
	cur_nodes[k]=(int *)malloc(2*sizeof(int));
	assert((int) cur_nodes[k]);
	}
blocked=(int *)malloc((2*(ntax-1))*sizeof(int));
assert((int) blocked);
is_it_or_desc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) is_it_or_desc);
if (values->atbr) {
	reroot_array=(int ***)malloc((2*ntax-3)*sizeof(int **));
	assert((int) reroot_array);
	for (i=0;i<2*ntax-3;i++) {
		reroot_array[i]=(int **)malloc((ntax-1)*sizeof(int *));
		assert((int) reroot_array[i]);
		for (j=0;j<ntax-1;j++) {
			reroot_array[i][j]=(int *)malloc(2*sizeof(int));
			assert((int) reroot_array[i][j]);
			}
		}
	}

anc=(int *)malloc(((2*ntax)-1)*sizeof(int));
assert((int) anc);

best_aligns=(alignment **)malloc(values->keep_aligns*sizeof(alignment *));
assert((int) best_aligns);
for (i=0;i<values->keep_aligns;i++) best_aligns[i]=NULL;

will_be_OK=1;
/*receive the info*/
pvm_upkint(&num_clipped,1,1);
for (i=0;i<ntax-1;i++) pvm_upkint(cur_nodes[i],2,1);
pvm_upkint(anc,((2*ntax)-1),1);
pvm_upkint(&d1,1,1);
pvm_upkint(&d2,1,1);
pvm_upkint(&old_length,1,1);
pvm_upkint(blocked,(2*(ntax-1)),1);
pvm_upkint(&node_removed,1,1);
pvm_upkint(&mult,1,1);
pvm_upkint(is_it_or_desc,((2*ntax)-1),1);
for (ii=ntax+1;ii<2*ntax-1;ii++) a[ii]=unpack_align_and_score(values);
tree_counter=0;
tot_best=old_length;

				/*when one or two taxa*/
				if ((num_clipped < 3) || (!values->atbr)) { /*if no alternative roots or sbr*/
					/*get nodes of clipped bit and rest*/
					all_diagnose_tree_here_pieces(a,cur_nodes,ntax,ntax,values,node_removed);
					/*then do the adding thing loop through places to put back*/
					/*check for inprovement*/
					best=0;
					for (j=1;j<ntax-1;j++) if ((j!=(node_removed-ntax)) && (a[ntax+j])) best+=a[ntax+j]->score;
					best=old_length-best;
					if ((best>0) || ((best==0) &&(mult))) {
					for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
						for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=cur_nodes[ii][jj];
						check=temp_nodes[j/2][j%2];
						temp_nodes[j/2][j%2]=node_removed;
						temp_nodes[d1][!d2]=check;
						if (values->groups) {
							cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
							will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
							for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
							free(cur_groups);
							}
						if (will_be_OK) {
							/*set up nodes to do and order*/
							for (k=0;k<ntax;k++) order[k]=(-1);
							order[0]=node_removed-ntax;
							num_to_do=1;
							/*if (j>1) {
								order[1]=(j/2);
								num_to_do=2;
								}*/
							done=0;
							while (!done) {
								done=1;
								for (k=1;k<ntax-1;k++) {
									same=0;
									for (kk=0;kk<num_to_do;kk++) if (order[kk]==k) same=1;
									if (!same) {
										for (kk=0;kk<num_to_do;kk++) {
											if (temp_nodes[k][0]==(order[kk]+ntax)) {order[num_to_do++]=k;done=0;}
											else if (temp_nodes[k][1]==(order[kk]+ntax)) {order[num_to_do++]=k;done=0;}
											}
										}
									}
								}
							/*do_alignment*/
							cur_best=all_diagnose_tree_cut_off_spr(a,temp_nodes,ntax,ntax,values,order,num_to_do,tot_best);
							/*if (cur_best <= tot_best) cur_best=all_diagnose_tree_here3(a,temp_nodes,ntax,ntax,values,HUGE_COST);
                            fprintf(stderr,"IF%d",cur_best);*/
							if (cur_best < tot_best) {
								/*fprintf(stderr," Found better at %d\n",cur_best);*/
								tot_best=cur_best;
								found_better=1;
								for (k=0;k<ntax-1;k++) {
									best_nodes_new[0][k][0]=temp_nodes[k][0];
									best_nodes_new[0][k][1]=temp_nodes[k][1];
									}
								/*add best_aligns stuf free first*/
								for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
								best_aligns[0]=make_align(a[best_nodes_new[0][0][1]]);
								best_aligns[0]->score=tot_best;
								/*if (values->new_optimization && (values->collapse != 2)) {
									if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*0,values);
									else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*0,values,values->support[0]);
									best_aligns[0]->score=tot_best;
									}*/
								tree_counter=1;
								}
							else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
								for (k=0;k<ntax-1;k++) {
									best_nodes_new[tree_counter][k][0]=temp_nodes[k][0];
									best_nodes_new[tree_counter][k][1]=temp_nodes[k][1];
									}
								/*chec to see if novel*/
								is_unique=1;
								if (1) /*(values->collapse==2) will get sorted out by master*/ is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
								else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*0,values,collapse_name);
								if (is_unique) {
									found_better=1;
									best_aligns[tree_counter]=make_align(a[best_nodes_new[tree_counter][0][1]]);
									best_aligns[tree_counter]->score=tot_best;
									if (values->new_optimization && (values->collapse!=2)) {
										free(best_aligns[tree_counter]->name);
										best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
										assert((int) best_aligns[tree_counter]->name);
										best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name);
										}
									++tree_counter;
									}
								/*fprintf(stderr,"O");*/
								}
							} /*OK*/
						} /*regraft 'j'*/
						} /*best*/
						else {values->number+=((2*(ntax-num_clipped))-3);}
					}/*too fewipped or SBR*/
					else { /*TBR*/
						check=get_root_array(reroot_array,num_clipped,node_removed,cur_nodes,d1,d2,ntax,is_it_or_desc,anc);
						for (l=0;l<2*num_clipped-3;l++) {
							/*get nodes of clipped bit and rest*/
							all_diagnose_tree_here_pieces(a,reroot_array[l],ntax,ntax,values,node_removed);
							/*then do the adding thing loop through places to put back*/
							best=0;
							for (j=1;j<ntax-1;j++) if ((j!=(node_removed-ntax)) && (a[ntax+j])) best+=a[ntax+j]->score;
							best=old_length-best;
							if ((best>0) || ((best==0) && (mult))) {
							for (j=1;j<(2*(ntax-1));j++) if (!blocked[j]) { /*regraft 'in principle' */
								for (ii=0;ii<ntax-1;ii++) for (jj=0;jj<2;jj++) temp_nodes[ii][jj]=reroot_array[l][ii][jj];
								check=temp_nodes[j/2][j%2];
								temp_nodes[j/2][j%2]=node_removed;
								temp_nodes[d1][!d2]=check;
								if (values->groups) {
									cur_groups=get_cur_groups(temp_nodes,ntax,ntax,&ngroups2);
									will_be_OK=compare_groups(values->groups,cur_groups,values->ngroups,ngroups2,ntax-1,ntax-1);
									for (ii=0;ii<ngroups2;ii++) free(cur_groups[ii]);
									free(cur_groups);
									}
								if (will_be_OK) {
									/*set up nodes to do and order*/
									for (k=0;k<ntax;k++) order[k]=(-1);
									order[0]=node_removed-ntax;
									num_to_do=1;
									/*if (j>1) {
										order[1]=(j/2);
										num_to_do=2;
										}*/
									done=0;
									while (!done) {
										done=1;
										for (k=1;k<ntax-1;k++) {
											same=0;
											for (kk=0;kk<num_to_do;kk++) if (order[kk]==k) same=1;
											if (!same) {
												for (kk=0;kk<num_to_do;kk++) {
													if (temp_nodes[k][0]==(order[kk]+ntax)) {order[num_to_do++]=k;done=0;}
													else if (temp_nodes[k][1]==(order[kk]+ntax)) {order[num_to_do++]=k;done=0;}
													}
												}
											}
										}
									/*do_alignment*/
									cur_best=all_diagnose_tree_cut_off_spr(a,temp_nodes,ntax,ntax,values,order,num_to_do,tot_best);
									/*if (cur_best <= tot_best) cur_best=all_diagnose_tree_here(a,temp_nodes,ntax,ntax,values,tot_best);
    						        fprintf(stderr,"IF%d",cur_best);*/
									if (cur_best < tot_best) {
										/*fprintf(stderr," Found better at %d\n",cur_best);*/
										tot_best=cur_best;
										found_better=1;
										for (k=0;k<ntax-1;k++) {
											best_nodes_new[0][k][0]=temp_nodes[k][0];
											best_nodes_new[0][k][1]=temp_nodes[k][1];
											}
										/*add best_aligns stuf free first*/
										for (k=0;k<tree_counter;k++) if (best_aligns[k]) best_aligns[k]=dump_align(best_aligns[k]);
										best_aligns[0]=make_align(a[best_nodes_new[0][0][1]]);
										best_aligns[0]->score=tot_best;
										/*if (values->new_optimization && (values->collapse != 2)) {

											if (!mult) check=get_collapsed_thang(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*0,values);
											else check=get_collapsed_thang_with_groups_old(best_aligns[0],best_nodes_new[0],ntax,check_buf,PARALLEL*0,values,values->support[0]);
											best_aligns[0]->score=tot_best;
											}*/
										tree_counter=1;
										}
									else if ((cur_best == tot_best) && mult && (tree_counter< values->keep_aligns)) {
										/*fprintf(stderr,"F");*/
										for (k=0;k<ntax-1;k++) {
											best_nodes_new[tree_counter][k][0]=temp_nodes[k][0];
											best_nodes_new[tree_counter][k][1]=temp_nodes[k][1];
											}
										/*chec to see if novel*/
										is_unique=1;
										if (1) /*(values->collapse==2) will get sorted out by master*/ is_unique=is_unique_tree(best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax);
										else is_unique=is_unique_tree_align(best_aligns,best_nodes_new[tree_counter],best_nodes_new,tree_counter,best_nodes_old,cur_num_trees,old_length+1,tot_best,ntax,check_buf,PARALLEL*0,values,collapse_name);
										if (is_unique) {
											/*fprintf(stderr," Found another for %d\n",tree_counter+1);*/
											found_better=1;
											best_aligns[tree_counter]=make_align(a[best_nodes_new[tree_counter][0][1]]);
											best_aligns[tree_counter]->score=tot_best;
											if (values->new_optimization && (values->collapse!=2)) {
												free(best_aligns[tree_counter]->name);
												best_aligns[tree_counter]->name=(char *)malloc((1+strlen(collapse_name))*sizeof(char));
												assert((int) best_aligns[tree_counter]->name);
												best_aligns[tree_counter]->name=(char *)strcpy(best_aligns[tree_counter]->name,collapse_name);
												}
											++tree_counter;
											}
										/*fprintf(stderr,"O");*/
										}
									} /*OK*/
								} /*regraft 'j'*/
							}/*best*/
							else {values->number+=((2*(ntax-num_clipped))-3);}
							} /*SBR*/
						}/*TBR*/


/*pack and send the results*/
pvm_initsend( PvmDataDefault );
pvm_pkint(&tot_best,1,1); /*pack cost*/
pvm_pkint(&tree_counter,1,1);/*pack #solns*/
for (i=0;i<tree_counter;i++) { /*make sure OK when no solutions found*/
	/*loop nodes then align*/
	for (j=0;j<ntax-1;j++) pvm_pkint(best_nodes_new[i][j],2,1);
	best_aligns[i]->score=tot_best;
	pack_align(best_aligns[i]);
	}
pvm_send(values->tids[0],PVM_NEW_SWAP_DONE_DAVE);

/*free everything*/
if (order) free(order);
values->support=NULL;
for (i=0;i<2*ntax-1;i++) if (check_buf[i]) check_buf[i]=dump_align(check_buf[i]);
free(check_buf);

for (i=0;i<ntax-1;i++) free(cur_nodes[i]);
free(cur_nodes);
free(is_it_or_desc);
free(blocked);
free(anc);
if (values->atbr) {
	for (i=0;i<2*ntax-3;i++) {
		for (j=0;j<ntax-1;j++) free(reroot_array[i][j]);
		free(reroot_array[i]);
		}
	free(reroot_array);
	}

for (i=0;i<values->keep_aligns;i++) {
	for (j=0;j<ntax-1;j++) {
		free(best_nodes_new[i][j]);
		}
	free(best_nodes_new[i]);
	}
free(best_nodes_new);
for (i=0;i<ntax-1;i++) free(temp_nodes[i]);
free(temp_nodes);

for (i=0;i<values->keep_aligns;i++) if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
free(best_aligns);
values->collapse=collapse_holder;
}





