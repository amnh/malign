/*Copyright 1992 Ward Wheeler all rights reserved*/
/*PVM parallel calls of both parent and child so only a single exe*/
/*
   need tp paralellize
   -pair
   -bandb
   -reorder
   -multi_swap

When all paralellized, disble cache for parent.
*/

#include "align3.h"

int ***best_nodes;

main (argc, argv)
int argc;
char *argv[];
{
     int i,j,k, bound,test,count,num_groups,kkk;
     int **advancement_scores;
     int ntaxa;
     FILE *file_in,*fopen(),*fp_input,*fp_output,*fp_groups;
     alignment **best_aligns;
     alignment **a;
     long int alpha, omega;
     parameters *values;
     char file_in_buffer[100],file_out_buffer[100],parameter_file[100],groups_file[100],check;
     int ran_rep;
     int mytid,pvm_checker,num_archs,me,i_am_the_parent;
     struct hostinfo *hostp;
     int killed, initialize,send_alignment,do_alignment,dump_message_for_align,info;
     int *d_tree_rep,d_sub_taxa,current_score, d_ntaxa;
     alignment *current_alignment;
     int **d_nodes,grain, **nodes,d_i, int_holder, numt, index;
     int num_bytes, tag, tid, bufid, number_holder, time_holder;
     int only_single_machine;
     alignment *temp1, *temp2, *temp3;
     int score_holder,node_received,time_holder2;
     int real_num_tot,thang;
     int *order,am_at_i,am_at_j,this_best,cur_best;
     alignment **a_best;
     char *temp_name;
     int num_seqs,bytes,source,gotten,seed,num_to_receive,type,best_k;
     int num_to_loop,initial_best,base,ii,jj,best_j,r_best_k;
     alignment **a_best_holder;

     a_best_holder=NULL;
     temp1=NULL;
     temp2=NULL;
     temp3=NULL;
     order=NULL;
     a_best=NULL;
     values=NULL;
     values=(parameters *)malloc(sizeof(parameters));
     assert((int)values);
     values->groups=NULL;
     get_initial_values(values);
     a=NULL;
     /*initialize*/
     fp_input=NULL;
     fp_output=NULL;
     fp_groups=NULL;
     advancement_scores=NULL;
     real_num_tot=0;

     i_am_the_parent=0;
     if (PARALLEL) {
	  /*PVM enroll in pvm*/
	  mytid = pvm_mytid();
	  /*PVM get hosts and numbers*/
	  pvm_checker=pvm_config(&(values->num_hosts), &num_archs, &hostp);/*check hostp for proper indirection*/
	  if (pvm_checker<0) {
	       fprintf(stderr,"Something bad happened.\n");
	       pvm_exit();
	       exit(-1);
	  }
	  /*PVM determine if parent or child*/
	  /*jack up processes to icrease number of jobs - performance tuning
	  values->num_hosts*=values->process_factor;*/
	  values->num_hosts*=1;
	  /*for test purposes jack up hosts to use argus alone*/
	  if (values->num_hosts==1) {
	    values->num_hosts+=10;
	    /*  fprintf(stderr,"NUM hosts=%d\n",values->single_host_addition);*/
	    only_single_machine=1;
	  }
	  else only_single_machine=0;
	  values->tids=(int *)malloc(values->num_hosts*sizeof(int));
	  assert((int) values->tids);
	  values->tids[0] = pvm_parent();
	  if(values->tids[0]==PvmNoParent) i_am_the_parent=1;
	  else i_am_the_parent=0;
	  /*PVM if parent*/
	  if (i_am_the_parent) {
	       values->tids[0] = mytid;
	       me = 0;
	       /*
	       fprintf(stderr,"I'm here with %d hosts/processes.\n",values->num_hosts);
	       */
	       /* start up copies of myself trying for only non-me hosts to be spawned*/
	       if (!only_single_machine) numt=pvm_spawn(EXE_FILE_NAME, (char**)0, PvmTaskHost|PvmHostCompl, ".", values->num_hosts-1, &(values->tids[1]));
	       else numt=pvm_spawn(EXE_FILE_NAME,(char **)0, 0,"",values->num_hosts-1,&(values->tids[1]));
	       if (numt!=(values->num_hosts-1)) {
		    fprintf(stderr,"Problems spawning tasks only %d of %d actually spawned.\n",numt,values->num_hosts-1);
		    for (i=0;i<values->num_hosts-1;i++) fprintf(stderr,"   host %s task id [%d]=%d\n",hostp->hi_name,i,values->tids[i]);
		    fprintf(stderr,"Reset PVM.\n");
		    pvm_exit();
		    exit(-1);
		    }
	       /* multicast values->tids array to children */
	       pvm_initsend( PvmDataDefault );
	       pvm_pkint(values->tids, values->num_hosts, 1);
	       info = pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_SEND_TIDS);
	       /*fprintf(stderr,"info = %d numnosts=%d tid[%d]=%d \n", info,values->num_hosts-1,1,values->tids[1]);
	       for (i=0;i<values->num_hosts;i++) {

		    }*/
	  }
     }

     if ((i_am_the_parent) || (!PARALLEL)) {
	  /*Print Banner*/
	  fprintf(stderr,"MALIGN version %s\n",malign_version);
	  if (PARALLEL) fprintf(stderr,"Parallel version %s\n",PARALLEL_VERSION);
	  fprintf(stderr,"by Ward Wheeler and David Gladstein\n");
	  fprintf(stderr,"Department of Invertebrates\n");
	  fprintf(stderr,"American Museum of Natural History\n");
	  fprintf(stderr,"Phone (212) 769-5754\n");
	  fprintf(stderr,"FAX   (212) 769-5233\n");
	  fprintf(stderr,"e-mail wheeler@amnh.org\n");
	  fprintf(stderr,"Copyright 1991-1995 all rights reserved\n\n");


	  /*Options to baby interface to redo stdin, and stdout*/
	  if (argc==1) {
	       fprintf(stderr,"What is the name of the sequence containing file?\n");
	       scanf("%s",file_in_buffer);
	       fprintf(stderr,"What is the name of output file?\n");
	       scanf("%s",file_out_buffer);
	       fprintf(stderr,"What is the name of the parameter file?\n");
	       scanf("%s",parameter_file);
	       get_file(parameter_file,values);
	       fprintf(stderr,"Will you specify alignment groups (y/n)?\n");
	       scanf("%s",&check);
	       if (isupper(check)) check=tolower(check);
	       if (check=='y') {
		    fprintf(stderr,"What is the name of this groups file?\n");
		    scanf("%s",groups_file);
	       }
	       fp_input=freopen(file_in_buffer,"r",stdin);
	       if (!fp_input) {
		    fprintf(stderr,"No such input file.\n");
		    exit(-1);
	       }
	       fp_output=freopen(file_out_buffer,"w",stdout);
	       if (!fp_output) {
		    fprintf(stderr,"Output file problems.\n");
		    exit(-1);
	       }
	  }
	  else if ((argc>1) && (argc<4)) {
	       get_file(argv[1],values);
	       check='n';
	  }
	  else if (argc>=4) {
	       get_command_line(argc,argv,values,groups_file,&check,file_in_buffer,file_out_buffer);
	       fp_input=freopen(file_in_buffer,"r",stdin);
	       if ((!fp_input) && (values->number_of_input_alignments<1)) {
		    fprintf(stderr,"No such input file %s.\n",file_in_buffer);
		    exit(-1);
	       }
	       fp_output=freopen(file_out_buffer,"w",stdout);
	       if (!fp_output) {
		    fprintf(stderr,"Output file problems.\n");
		    exit(-1);
	       }
	  }
	  /*have values now*/
	  /*use as if stdin*/
	  if (values->number_of_input_alignments>0) {
	    /*read input files*/
	    if (values->VERBOSE) fprintf(stderr,"Reading pre-aligned sequences.\n");
	    if ((values->cull) || (values->elision)) {
		    a=get_prealigned_sequences(values);
    		for (i=0;i<values->all_done;i++) for (j=0;j<a[i]->n_seqs;j++) for (k=0;k<a[i]->length;k++) if (islower(a[i]->s[j][k])) a[i]->s[j][k]=toupper(a[i]->s[j][k]);
    		for (i=0;i<values->all_done;i++) for (j=1;j<a[i]->n_seqs;j++) for (k=0;k<a[i]->length;k++) if (a[i]->s[j][k]=='.') a[i]->s[j][k]=a[i]->s[0][k];
	    	}
	    else {
		    if (!values->chop) values->chop=huge_length;
    		a=new_get_prealigned_sequences_chop(values);
    		if (values->chop==huge_length) values->chop=0;
       		}
    	  }
	  else {
		if (values->input_is_align) {
			if (values->VERBOSE) fprintf(stderr,"Reading pre-aligned sequences.\n");
			a=read_input_align(values,a,&ntaxa);
			}
		else a=read_old_file(a,values);
		}
	if (!values->cull && (!values->elision)) do_modifications(a,values);
	  /*for (i=0;i<values->all_done;i++) fprintf(stderr,"%s->%d seqs\n",a[i]->name,a[i]->n_seqs);*/
	  do_pre_stuff(values,a);
	 /*allocate cache*/
	 /* if (!(values->new_optimization && (values->number_of_input_alignments > 1))) if (!values->chop) for (i=0;i<values->all_done;i++) a[i]->type_weight=0;*/
	  /*could spawn here*/
	  if (!PARALLEL) if (values->cache_size>0) values->align_cache=allocate_align_cache(values->cache_size,values,a);

	  bound=HUGE_COST;
	  file_in=NULL;
	  if (argc==3) {
	       file_in=fopen(argv[2],"r");
	       if (!file_in) {
		    fprintf(stderr,"No such groups file\n");
		    exit(-1);
	       }
	  }
	  else if (check=='y') {
	       file_in=fopen(groups_file,"r");
	       if (!file_in) {
		    fprintf(stderr,"No such groups file\n");
		    exit(-1);
	       }
	  }
	  if (file_in) {
	       values->groups=get_groups2(file_in,values->all_done,a,values);
	       test=compare_groups(values->groups,values->groups,values->ngroups,values->ngroups,values->all_done,
		   values->all_done);
	       if (!test) {
		    fprintf(stderr,"Specified groups are inconsistent\n");
		    exit(1);
	       }
	       if (values->get_heur) {
		    fprintf(stderr,"Option (pair) does not operate with specified alignment groups.\n");
		    values->get_heur=0;
	       }
	  }

	  best_aligns = (alignment **) malloc((values->keep_aligns) * sizeof (alignment *));
	  assert((int)best_aligns);
	  for (j=0;j<(values->keep_aligns);j++) {
		best_aligns[j]=NULL;
		}
	  count=0;

	  values->input_names=(char **)malloc(values->all_done*sizeof(char *));
	  assert((int)values->input_names);
	  for (j=0;j<values->all_done;j++) {
	       values->input_names[j]=(char *)malloc((1+(strlen(a[j]->name)))*sizeof(char));
	       assert((int)values->input_names[j]);
	       strcpy(values->input_names[j],a[j]->name);
	  }
	  if ((values->phylo_score==2)||(values->phylo_score==3)||(values->phylo_score==4)||(values->how_many>0)) {
	       values->best_rep=(int **)malloc(values->keep_trees*sizeof(int *));
	       assert((int)values->best_rep);
	       for (j=0;j<values->keep_trees;j++) {
		    values->best_rep[j]=(int *)malloc((values->all_done-3)*sizeof(int));
		    assert((int)values->best_rep[j]);
	       }
	  }

	  if (values->iter) get_guess_max_gap(a,values->all_done,values);
	  alpha=(long int) time(NULL);
	  srand((int) alpha);


	  /*fixes for multiple param files*/
	  if (values->other_parm) {
		for (i=0;i<values->number_of_input_alignments;i++) if (values->other_parm[i]) {
			values->other_parm[i]->new_optimization=values->new_optimization;
			values->other_parm[i]->phylo_score=0;
			values->other_parm[i]->all_done=values->all_done;
			values->other_parm[i]->VERBOSE=values->VERBOSE;
			values->other_parm[i]->rep_error=values->rep_error;
			values->other_parm[i]->groups=NULL;
			values->other_parm[i]->ngroups=0;
			values->other_parm[i]->cache_size=0;
			values->other_parm[i]->align_cache=NULL;
			values->other_parm[i]->gap_must_cost=values->gap_must_cost;
			values->other_parm[i]->shortcut_tree=values->shortcut_tree;
			values->other_parm[i]->shortcut_tree2=values->shortcut_tree2;
			values->other_parm[i]->other_parm=NULL;
			}
		}
	  /*PVM broadcast initial values*/
	  if (PARALLEL) do_initialize_thang(values,a);

	  /*do something searches*/
	  if (values->get_heur) {
	       best_aligns[0]=(minimum_spanning_alignment (a, values->all_done, values));
	       if (values->VERBOSE) fprintf(stderr,"Found an alignment at cost %d\n",best_aligns[0]->score);
	       if (!values->groups) bound=best_aligns[0]->score;
	       count=1;
	       if (values->print_intermediates) printem(best_aligns,count,values,a);
	  }

	  if (values->rand_align) {
	       best_aligns=do_some(a,values->all_done,values->how_many,bound,best_aligns, values,&count);
	       if (values->print_intermediates==1) printem(best_aligns,count, values,a);
	       if (bound>best_aligns[0]->score) bound=best_aligns[0]->score;
		}

	  if ((values->get_heur2) || (values->best)) {
	       advancement_scores=(int **)malloc((values->all_done)*sizeof(int*));
	       assert((int)advancement_scores);
	       for (i=0;i<(values->all_done);i++) {
		    advancement_scores[i]=(int *)malloc((values->all_done)*sizeof(int));
		    assert((int)advancement_scores[i]);
		    for (j=0;j<(values->all_done);j++) *(advancement_scores[i]+j)=0;
	       }
	  }
/*        for (i=0;i<a[0]->length;i++) fprintf(stderr,"%d",a[0]->s[1][i]);
	  fprintf(stderr,"\n");
*/
if (!(PARALLEL*values->jackboot)) {
	  for (ran_rep=0;ran_rep<max(1,values->rand_order);ran_rep++) {
		if (values->jackboot>0) {
			bound=HUGE_COST;
			if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
				for (i=0;i<1000;i++) {
				values->jack_array[i]=(!(rand()/(632*(max_rand/1000))));
				/*fprintf(stderr,"%d",values->jack_array[i]);*/
				/*remember to resend if parallel*/
				}
			}
	       if (values->VERBOSE) if (values->rand_order>0) fprintf(stderr,"\nPerforming random rep %d of %d",ran_rep+1,values->rand_order);
		if (ran_rep==0) {
			 /*use minimum spanning name to reorder the a[i]'s*/
			if (values->reorder) {
				fprintf(stderr,"	");
				a=new_reorder(a,values->all_done,0,advancement_scores,values->groups,values);
				if (PARALLEL) {
					redo_it_all(a,only_single_machine,values);
					real_num_tot+=values->number;
					}
				}
			else if (values->rand_order > 0) {
				a=randomize_taxon_order(a,values,values->groups);
				if (PARALLEL) {
					redo_it_all(a,only_single_machine,values);
					real_num_tot+=values->number;
					}
				real_num_tot+=values->number;
				}
			}
		if ((ran_rep>0) || values->jackboot) {
			a=randomize_taxon_order(a,values,values->groups);
			if (PARALLEL) {
				redo_it_all(a,only_single_machine,values);
				real_num_tot+=values->number;
				}
			real_num_tot+=values->number;
			}
		else fprintf(stderr,"\n");
		if (values->get_heur2) best_aligns=make_heur_align(a,values->all_done,best_aligns,values,count,&count,1);
		if (values->print_intermediates) {
			printem(best_aligns,count,values,a);
			for (kkk=0;kkk<count;kkk++) fprintf(stderr,"%s\n",best_aligns[kkk]->name);
			}
		if (bound>best_aligns[0]->score) bound=best_aligns[0]->score;
		if (values->jackboot>0) {
			values->jack_name[ran_rep]=(char *)malloc((1+strlen(best_aligns[0]->name))*sizeof(char));
			assert((int) values->jack_name[ran_rep]);
			strcpy(values->jack_name[ran_rep],best_aligns[0]->name);
			/*fprintf(stderr,"%d %s\n",ran_rep,values->jack_name[ran_rep]);*/
			}
		}/*end ran rep*/
	}
else {
	gotten=0;
	fprintf(stderr,"Performing random replicate ");
	temp_name=(char *)malloc(FIFTEEN_BITS_MAX*sizeof(char));
	assert((int) temp_name);
	num_to_receive=0;
	for (ran_rep=0;((ran_rep<(values->num_hosts-1)) && (ran_rep<values->jackboot));ran_rep++) {
		/*send out*/
		pvm_initsend( PvmDataDefault );
		/*pack info*/
		seed=rand();
		pvm_pkint(&seed,1,1);
		pvm_send(values->tids[++num_to_receive],PVM_JACK_DO);
		fprintf(stderr,"%d ",ran_rep+1);
		}
	while (num_to_receive>0) {
		/*receive*/
		bufid=pvm_recv(-1,PVM_JACK_DONE);
		if (bufid) {
			info=pvm_bufinfo(bufid, &bytes, &type, &source);
			if (ran_rep<values->jackboot) {
				/*send out*/
				pvm_initsend( PvmDataDefault );
				/*pack info*/
				seed=rand();
				pvm_pkint(&seed,1,1);
				pvm_send(source,PVM_JACK_DO);
				ran_rep++;
				num_to_receive++;
				fprintf(stderr,"%d ",ran_rep);
				}
			pvm_upkstr(temp_name);
			/*fprintf(stderr,"%s\n",temp_name);*/
			values->jack_name[gotten]=(char *)malloc((1+strlen(temp_name))*sizeof(char));
			assert((int) values->jack_name[gotten]);
			strcpy(values->jack_name[gotten],temp_name);
			gotten++;
			--num_to_receive;
			values->number+=((2*(values->all_done+1))-3);
			}
		}
	fprintf(stderr,"Done\n");
	free(temp_name);
	if (!best_aligns[0]) best_aligns[0]=make_align(a[0]);
	}
	  if (values->jackboot>0) num_groups=get_jack_groups(best_aligns,a,values);
	  if (values->best) {
	       /*reorder does not seem to work ok due to leading and trailing gaps as missing*/
		a=new_reorder(a,values->all_done,1,advancement_scores,values->groups,values);
		if (PARALLEL) {
			redo_it_all(a,only_single_machine,values);
			real_num_tot+=values->number;
			}

		best_aligns=do_all (a, values->all_done,bound,best_aligns,values,count,&count);
		}
	  /*do cull and elision stuff*/
	  if (values->print_intermediates && values->cull) {
	       printf("Aligned alignments before 'culling':\n");
	       printem(best_aligns,count,values,a);
	  }
	  if (values->cull) best_aligns=do_cull(a,values,best_aligns,&count);
	  if (values->elision) best_aligns=do_elision(a,values,best_aligns,&count);
	  /*finish up*/
	  if (values->time) omega=(long int) time(NULL);
	  if (values->print_intermediates) printf("Best one(s)\n");
	  /*expand '-' => 'X' when contig*/
	  if (!values->new_optimization) for (i=0;i<count;i++) make_contig_X(best_aligns[i],values);
	  printem(best_aligns,count,values,a);
	  if (values->givem) {
	       for (i=0;i<values->all_done;i++) fprintf(stderr,"%s #\n",a[i]->name);
	       for (i=0;i<values->all_done;i++) print_alignment(a[i]);
	  }

	  if ((values->cost_only) || ((values->get_heur==0) && (!values->get_heur2) && (!values->rand_align) &&
	      (!values->get_best) && (!values->new_optimization)))    {
	       values->length_best_clad=cladogram(a[ntaxa],values);
	       fprintf(stderr,"Cost of alignment is %d\n",values->length_best_clad);
	       printf("Cost of alignment is %d\n",values->length_best_clad);
	  }

					/*optimization stuff for optali*/
					if (values->new_optimization && values->optimize_nodes) {
						fprintf(stderr,"	Optimizing nodes\n");
						printf("\n***********************\n");
						printf("   Optimizing nodes\n");
						printf("***********************\n");
						/*do_and_print_optimization(a,values,count,best_aligns);*/
						for (i=0;i<count;i++) {
						    fprintf(stderr,"Cladogram %d\n",i);
						    printf("Cladogram %d\n",i);
						    do_and_print_optimization2(a,values,best_aligns[i],best_nodes[i]);
						}
						}

	 /*free allocated stuff*/
	  if (advancement_scores) {
	       for (i=0;i<values->all_done;i++) if (advancement_scores[i]) free(advancement_scores[i]);
	       free(advancement_scores);
	  }

	  if (values->saw_farris+values->saw_goloboff+values->saw_wheeler>0) do_estimated_weights(best_aligns,values,
	      count);
	  if (fp_input)  fclose (fp_input);
	  if (fp_output) fclose (fp_output);
	  if (values->time==1) fprintf(stderr,"Alignment took %ld seconds\n",(long int) (omega-alpha));
	  if (PARALLEL) {
	       time_holder2=0;
		if (!values->new_optimization) {
			pvm_initsend( PvmDataDefault );
			pvm_pkint(&(values->phylo_time), 1, 1);
			pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_REQUEST_TIME);
			for (i=0;i<(values->num_hosts-1);i++) {
				pvm_recv(values->tids[i+1],PVM_SEND_TIME);
				pvm_upkint(&time_holder,1,1);
				time_holder2+=time_holder;
				}
			values->phylo_time+=(time_holder2/(values->num_hosts-1));
			}
			pvm_initsend( PvmDataDefault );/*should this be modified for multiplatform?*/
			pvm_pkint(&(values->number), 1, 1);
			pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_REQUEST_NUM_TREES);
			for (i=0;i<(values->num_hosts-1);i++) {
				pvm_recv(values->tids[i+1],PVM_SEND_NUM_TREES);
				pvm_upkint(&number_holder,1,1);
			values->number+=number_holder;
			     }
			 values->number+=real_num_tot;
		   }
	     else values->number+=real_num_tot;
	  if (values->make_groups_file) {
	  	fp_groups=(FILE *)fopen(values->output_groups_name,"w");
	  	if (!fp_groups) fprintf(stderr,"Error in making groups file\n");
	  	else {
	  		fprintf(fp_groups,"%s",best_aligns[0]->name);
	  		fclose(fp_groups);
	  		}
	  	}
	  if (!values->new_optimization && (values->phylo_time>-1)) fprintf(stderr,"  Phylogeny reconstructions took %ld seconds.\n",values->phylo_time);
	  if (CURRENT_PLATFORM==MAC) fprintf(stderr,"All done now.  You can exit.\n");
	  if (PARALLEL) {
	       /*PVM send kill to processes*/
	       pvm_initsend( PvmDataDefault );
	       pvm_pkint(values->tids, values->num_hosts, 1);
	       pvm_mcast(&(values->tids[1]), values->num_hosts-1, PVM_ALL_MUST_DIE);
	       /*PVM end parent*/
	       }
     }
     /*PVM if child */
     if ((!i_am_the_parent) && (PARALLEL)) {
       /*PVM receive values->tids array */
       pvm_recv(values->tids[0], PVM_SEND_TIDS);
       pvm_upkint(values->tids, values->num_hosts, 1);
       for( i=1; i<values->num_hosts ; i++ ) if( mytid == values->tids[i] ){
	 me = i;
	 break;
       }
       current_alignment=NULL;
       best_aligns = (alignment **) malloc((values->keep_aligns) * sizeof (alignment *));
       assert((int)best_aligns);
       for (j=0;j<(values->keep_aligns);j++) {
	 best_aligns[j]=NULL;
       }

       /*PVM recieve instruction-some logical loop*/
       killed=0;
       while (!killed) {
	 bufid=pvm_recv(values->tids[0],-1);
	 info=pvm_bufinfo(bufid,&num_bytes, &tag,&tid);
	 /*fprintf(stderr,"Tag%d",tag);*/
	 if (tag==PVM_ALL_MUST_DIE) killed=1;
	 else if (tag==PVM_UNPACK_INITIAL_SHIT) {
	   a=unpack_init_vals(a,values);
	   /*the values->all_done shit is because of the ++ntaxa fool thing*/
	   nodes=(int **)malloc(((values->all_done+1)-1)*sizeof(int *));
	   assert((int)nodes);
	   for (i=0;i<((values->all_done+1)-1);i++) {
	     nodes[i]=(int *)malloc(2*sizeof(int));
	     assert((int)nodes[i]);
	   }
	   d_nodes=(int **)malloc(((values->all_done+1)-1)*sizeof(int *));
	   assert((int) d_nodes);
	   for (i=0;i<((values->all_done+1)-1);i++) {
	     d_nodes[i]=(int *)malloc(2*sizeof(int));
	     assert((int) d_nodes[i]);
	   }
	   d_tree_rep=(int *)malloc((values->all_done+1)*sizeof(int));
	   assert((int) d_tree_rep);
	   current_alignment=NULL;
	   order=(int *)malloc((values->all_done+1)*sizeof(int));
	   assert((int) order);
	   a_best=(alignment **)malloc(((values->all_done+1)-1)*sizeof(alignment *));
	   assert((int) a_best);
	   for (i=0;i<values->all_done;i++) a_best[i]=NULL;

	   a_best_holder=(alignment **)malloc(((values->all_done+1)-1)*sizeof(alignment *));
	   assert((int) a_best_holder);
	   for (i=0;i<values->all_done;i++) a_best_holder[i]=NULL;

	   fprintf(stderr,"From %d of %d: I'm initialized with %d seqs (best %d)\n",me,values->num_hosts,values->all_done,values->best_order);
	   pvm_initsend( PvmDataDefault );
	   pvm_pkint(&me,1,1);
	   pvm_send(values->tids[0],PVM_INITIAL_SHIT_RECIEVED);
	   /*PVM all values->x shit*/
	 }
	 else if (tag==PVM_DO_ALIGNMENT_NOW) do_alignment_now(d_nodes,nodes,d_tree_rep,a,values);
	 else if (tag==PVM_DO_ALIGNMENT_SCORE_ONLY) do_alignment_score_only(d_nodes,nodes,d_tree_rep,a,values);
	 else if (tag==PVM_DO_ALIGNMENT_NOW_BANDB) do_alignment_now_bandb(d_nodes,nodes,d_tree_rep,a,values);
	 else if (tag==PVM_DO_ALIGNMENT_SCORE_ONLY_BANDB) do_alignment_score_only_bandb(d_nodes,nodes,d_tree_rep,a,values);
	 else if (tag==PVM_REQUEST_TIME) {
	   pvm_upkint(&time_holder,1,1);
	   pvm_initsend( PvmDataDefault );
	   pvm_pkint(&(values->phylo_time),1,1);
	   pvm_send(values->tids[0],PVM_SEND_TIME);
	 }
	 else if (tag==PVM_DUMP_INITIAL_GROUPS) {
	   /*pvm_upkint(&(values->groups_as_start),1,1);*/
	   for (i=0;i<values->ngroups;i++) free(values->groups[i]);
	   free(values->groups);
	   values->groups=NULL;
	   values->ngroups=0;
	 }
	 else if (tag==PVM_REQUEST_NUM_TREES) {
	   pvm_upkint(&number_holder,1,1);
	   pvm_initsend( PvmDataDefault );
	   pvm_pkint(&(values->number),1,1);
	   pvm_send(values->tids[0],PVM_SEND_NUM_TREES);
	 }
	 else if (tag==PVM_REQUEST_DIAGNOSE_ALIGNMENT) do_diagnose_alignment(nodes,a,values);
	 else if (tag==PVM_NEW_BUILD_START) do_new_build_parallel(a,values->all_done+1,values);
	 else if (tag==PVM_NEW_BUILD_ADD_FUNCTION) do_new_build_parallel2(a,values->all_done+1,values);
	 else if (tag==PVM_NEW_SWAP) do_new_parallel_swap(a,values->all_done+1,values);
	 else if (tag==PVM_NEW_SWAP_DAVE) do_new_parallel_swap_dave(a,values->all_done+1,values);
	 else if (tag==PVM_PARALLEL_NW_NEW_OPT_UP_DO) thang=up_node_remote(values);
	 else if (tag==PVM_PARALLEL_NW) {
	   if (temp1) temp1=dump_align(temp1);
	   if (temp2) temp2=dump_align(temp2);
	   if (temp3) temp3=dump_align(temp3);
	   score_holder=values->phylo_score;
	   pvm_upkint(&values->phylo_score,1,1);
	   pvm_upkint(&node_received,1,1);
	   temp1=unpack_align_and_score(temp1);
	   temp2=unpack_align_and_score(temp2);
	   temp3=nw(temp1,temp2,values);
	   pvm_initsend( PvmDataDefault );
	   pvm_pkint(&node_received,1,1);
	   pack_align(temp3);
	   pvm_send(values->tids[0],PVM_PARALLEL_NW_DONE);
	   if (temp1) temp1=dump_align(temp1);
	   if (temp2) temp2=dump_align(temp2);
	   if (temp3) temp3=dump_align(temp3);
	   values->phylo_score=score_holder;
	 }
	 else if (tag==PVM_PARALLEL_NW_NEW_OPT) {
	   if (temp1) temp1=dump_align(temp1);
	   if (temp2) temp2=dump_align(temp2);
	   if (temp3) temp3=dump_align(temp3);
	   score_holder=values->in_optimize;
	   pvm_upkint(&values->in_optimize,1,1);
	   pvm_upkint(&node_received,1,1);
	   temp1=unpack_align_and_score(temp1);
	   temp2=unpack_align_and_score(temp2);
	   /*cache stuff*/
	   if (values->align_cache && (!values->in_optimize)) {
	     if (strcmp (temp1->name,temp2->name) > 0) temp_name=other_pair_names(temp2->name,temp1->name);
	     else temp_name=other_pair_names(temp1->name,temp2->name);
	     num_seqs=1;
	     for (j=0;j<strlen(temp_name);j++) if (temp_name[j]=='(') ++num_seqs;
	     temp3=retrieve_align_cache(temp_name,values,values->align_cache[num_seqs-2]);
	     free(temp_name);
	   }
	   if (!temp3) {
	     temp3 = nw(temp1,temp2,values);
	     if (values->align_cache && (!values->in_optimize)) store_align_cache(temp3,values,values->align_cache[num_seqs-2]);
	   }
	   pvm_initsend( PvmDataDefault );
	   pvm_pkint(&node_received,1,1);
	   pack_align(temp3);
	   pvm_send(values->tids[0],PVM_PARALLEL_NW_NEW_OPT_DONE);
	   if (temp1) temp1=dump_align(temp1);
	   if (temp2) temp2=dump_align(temp2);
	   if (temp3) temp3=dump_align(temp3);
	   values->in_optimize=score_holder;
	 }
	 else if (tag==PVM_REDO_INITIAL_ALIGNS) {
	   for (j=0;j<(2*(values->all_done+1))-1;j++) {
	     /*for (j=0;j<values->all_done;j++) {*/
	     if (a[j]) a[j]=dump_align(a[j]);
	     a[j]=NULL;
	   }
	   free(a);
	   a=(alignment **)malloc(((2*(values->all_done+1))-1)*sizeof(alignment *));
	   assert((int) a);
	   for (j=0;j<values->all_done;j++) a[j]=(alignment *)unpack_align_and_score((parameters *) values);
	   if (values->ngroups>0) for (i=0;i<values->ngroups;i++) for (j=0;j<(values->all_done+1);j++) pvm_upkint(&(values->groups[i][j]),1,1);
	 }
	 else if (tag==PVM_SEND_ANODES) {
	   pvm_upkint(&am_at_i,1,1);
	   if (values->best_order)  {
	      for (j=0;j<values->all_done;j++) {
		if (a[j]) a[j]=dump_align(a[j]);
		a[j]=unpack_align_and_score(a[j]);
	     }
	   }
	   for (j=1;j<am_at_i+2;j++) {
	     if (a[(values->all_done+1)+j]) a[(values->all_done+1)+j]=dump_align(a[(values->all_done+1)+j]);
	     a[(values->all_done+1)+j]=(alignment *)unpack_align_and_score((parameters *) values);
	   }
	   for (ii=0;ii<am_at_i+3;ii++) for (jj=0;jj<2;jj++) pvm_upkint(&(nodes[ii][jj]),1,1);
	   /*fprintf(stderr,"GotNodes%d",am_at_i);*/
	 }
	 else if (tag==PVM_DAVE_DUMP_NODES) {for (j=1;j<am_at_i+3;j++) if (a[(values->all_done+1)+j]) a[(values->all_done+1)+j]=dump_align(a[(values->all_done+1)+j]);}
	 else if (tag==PVM_DAVE_BUILD)  {
	   pvm_upkint(&am_at_i,1,1);
	   /*if (values->best_order) {
	     if (am_at_i>0) {
	       for (j=am_at_i-1;j<values->all_done;j++) {
		 if (a[j]) a[j]=dump_align(a[j]);
		 a[j]=unpack_align_and_score(a[j]);
	       }
	     }

	   }*/
	   pvm_upkint(&am_at_j,1,1);
	   pvm_upkint(&cur_best,1,1);
	   initial_best=cur_best;
	   pvm_upkint(&num_to_loop,1,1);
	   pvm_upkint(&base,1,1); /*starting point==k in caller*/
	   best_j=base;
	   r_best_k=am_at_i+2;
	   for (ii=0;ii<am_at_i+2;ii++) for (jj=0;jj<2;jj++) d_nodes[ii][jj]=nodes[ii][jj];
	   d_nodes[am_at_i+2][0]=am_at_i+3;
	   for (k=0;k<num_to_loop;k++) {
	     d_nodes[am_at_i+2][1]=nodes[(base+k)/2][(base+k)%2];
	     d_nodes[(base+k)/2][(base+k)%2]=(values->all_done+1)+am_at_i+2;
	     this_best=do_dave_build(a,d_nodes,values->all_done+1,am_at_i,values,order,cur_best,a_best,&best_k);
	     if (this_best<cur_best) if (this_best<HUGE_COST) {
	       cur_best=this_best;
	       best_j=base+k;
	       r_best_k=best_k;
	       /*for (ii=0;ii<am_at_i+2;ii++) {
		 if (a_best_holder[ii]) a_best_holder[ii]=dump_align(a_best_holder[ii]);
		 a_best_holder[ii]=make_align(a_best[ii]);
		 }*/
	     }
	     d_nodes[am_at_i+2][1]=nodes[am_at_i+2][1];
	     d_nodes[(base+k)/2][(base+k)%2]=nodes[(base+k)/2][(base+k)%2];
	    /*fprintf(stderr,"(%d-%d-%d=>%d)",am_at_i,base+k,best_k,this_best);*/
	   }
	   pvm_initsend( PvmDataDefault );
	   pvm_pkint(&best_j,1,1);
	   pvm_pkint(&cur_best,1,1);
	   if (cur_best<initial_best) if (cur_best<HUGE_COST) {
	     pvm_pkint(&r_best_k,1,1);
	     for (i=0;i<am_at_i+2;i++) pack_align(a_best[i]);
	   }
	   pvm_send(values->tids[0],PVM_DAVE_BUILD_DONE);
	 }
	 else if (tag==PVM_JACK_DO) {
	   pvm_upkint(&seed,1,1);
	   srand(seed);
	   /*randomize taxon order and make array*/
	   /*fprintf(stderr,"Here");*/
	   if (!values->jack_array) {
	     values->jack_array=(int *)malloc(1000*sizeof(int));
	     assert((int) values->jack_array);
	   }
	   for (j=0;j<1000;j++) values->jack_array[j]=(!(rand()/(632*(max_rand/1000))));
	   a=randomize_taxon_order(a,values,values->groups);
	   bound=HUGE_COST;
	   if (best_aligns[0]) best_aligns[0]=dump_align(best_aligns[0]);
	   count=0;
	   if (values->get_heur2) best_aligns=make_heur_align(a,values->all_done,best_aligns,values,count,&count,0);/*do in parallel if 1, 0 if not*/
	   pvm_initsend( PvmDataDefault );
	   pvm_pkstr(best_aligns[0]->name);
	   pvm_send(values->tids[0],PVM_JACK_DONE);

	 }
	 else {
	   fprintf(stderr,"Bad tag %d received. Bye.\n",tag);
	   killed=1;
	 }

       }
       if (order) free(order);
       if (d_tree_rep) free(d_tree_rep);
       for (i=0;i<((values->all_done+1)-1);i++) {
	 if (nodes[i]) free(nodes[i]);
	 if (d_nodes[i]) free(d_nodes[i]);
	 if (a_best) if (a_best[i]) a_best[i]=dump_align(a_best[i]);
	 if (a_best_holder) if (a_best_holder[i]) a_best_holder[i]=dump_align(a_best_holder[i]);
       }
       if (a_best) free(a_best);
       if (a_best_holder) free(a_best_holder);
       if (nodes) free(nodes);
       if (d_nodes) free(d_nodes);
       /*remember to free up at end*/
       for (i=0;i<((2*values->all_done+1)-1);i++) if (a[i]) a[i]=dump_align(a[i]);
       if (a) free(a);
       end_free(values);
       free(values);
       /*PVM end child*/
       fprintf(stderr,"Done as child.\n");
     }
     /*PVM Leave PVM*/
     if (PARALLEL) {
       pvm_exit();
       if (values->tids) free(values->tids);

     }
     if (a) {
       for (i=0;i<((2*values->all_done+1)-1);i++) if (a[i]) a[i]=dump_align(a[i]);
       free(a);
     }
     if (values->VERBOSE) fprintf(stderr,"Performed %d multiple alignments.\n",values->number);
     for (i=0;i<values->keep_aligns;i++) if (best_aligns[i]) best_aligns[i]=dump_align(best_aligns[i]);
     free(best_aligns);
     end_free(values);
     free(values);

     /*make sure its not trying to free a cache not in parent*/
     return 0;
   }







