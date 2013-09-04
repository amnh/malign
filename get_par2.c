/*Copyright 1994 Ward Wheeler all rights reserved*/
#include "align3.h"

int charbit[256];
char bitchar[32];

void get_align_file(file_name,values)
char *file_name;
parameters *values;
{
     char s[100],cval[100],found,cval2[100],temp_option[100];
     int val,i,j,check;
     FILE *fp_in, *fopen(),*fp_param;

     fp_in=fopen(file_name,"r");
     temp_option[0]='\0';
     if (fp_in) {
	  while (!feof(fp_in)) {
	       if (temp_option[0]=='\0') fscanf(fp_in,"%s ",s);
	       else {
		for (i=0;i<strlen(temp_option);i++) s[i]=temp_option[i];
		s[strlen(temp_option)]='\0';
		}
		  temp_option[0]='\0';
	       found=0;
	       if (strnicmp("infile",s,6)==0) {
			 if (!values->number_of_input_alignments) {
			      values->data_sets=(int *)malloc((1+values->number_of_input_alignments)*sizeof(int));
			      assert((int) values->data_sets);
			      values->data_set_weights=(int *)malloc((1+values->number_of_input_alignments)*sizeof(int));
			      assert((int) values->data_set_weights);
			      }
			 else {
			      values->data_sets=(int *)realloc(values->data_sets,(1+values->number_of_input_alignments)*sizeof(int));
			      assert((int) values->data_sets);
			      values->data_set_weights=(int *)realloc(values->data_set_weights,(1+values->number_of_input_alignments)*sizeof(int));
			      assert((int) values->data_set_weights);
			      }
		    /*if (!values->new_optimization) {
			values->data_sets[values->number_of_input_alignments]=0;
			 values->data_set_weights[values->number_of_input_alignments]=1;
			}
		    else {*/
		    if (1) {
			 fscanf(fp_in,"%s ",cval2);
			 if (!strnicmp("a",cval2,1)) values->data_sets[values->number_of_input_alignments]=0;
			 else      if (!strnicmp("n",cval2,1)) values->data_sets[values->number_of_input_alignments]=1;
			 else      if (!strnicmp("m",cval2,1)) values->data_sets[values->number_of_input_alignments]=2;
			 else      if (!strnicmp("g",cval2,1)) values->data_sets[values->number_of_input_alignments]=3;
			 else {
			      fprintf(stderr,"Unrecognized option modifier (%s) after (infile) in file (%s).\n     Should be align, noalign, or morph when optalign specified.\n",cval2,file_name);
			      exit(-1);
			      }
			 fscanf(fp_in,"%d ",&val);
			 values->data_set_weights[values->number_of_input_alignments]=val;
			 }
		    ++values->number_of_input_alignments;
		    fscanf(fp_in,"%s ",cval);
		    values->input_file_name=(char **)realloc(values->input_file_name,values->number_of_input_alignments*sizeof(char *));
		    assert((int)values->input_file_name);
		    values->input_file_name[values->number_of_input_alignments-1]=(char *)malloc((1+strlen(cval))*sizeof(char));
		    assert((int)values->input_file_name[values->number_of_input_alignments-1]);
		    values->input_file_name[values->number_of_input_alignments-1]=(char *)strcpy(values->input_file_name[values->number_of_input_alignments-1],cval);
		    /*add parameter allocate and read*/
		    check=fscanf(fp_in,"%s ",cval);
		    if (check) {
		        if (!strnicmp(cval,"param",5)) {
            	     fscanf(fp_in,"%s ",cval);
        			if (values->other_parm) values->other_parm=(parameters **)realloc(values->other_parm,(values->number_of_input_alignments)*sizeof(parameters *));
        			else values->other_parm=(parameters **)malloc((values->number_of_input_alignments)*sizeof(parameters *));
    	    		assert((int) values->other_parm);
    		    	values->other_parm[values->number_of_input_alignments-1]=(parameters *)malloc(sizeof(parameters));
        			assert((int) values->other_parm[values->number_of_input_alignments-1]);
        			get_initial_values(values->other_parm[values->number_of_input_alignments-1]);
        			get_file(cval,values->other_parm[values->number_of_input_alignments-1]);
        			fix_all_after_values2(values->other_parm[values->number_of_input_alignments-1]);
        			}
    		    else {
    				for (i=0;i<strlen(cval);i++) temp_option[i]=cval[i];
    				temp_option[strlen(cval)]='\0';
    				if (values->other_parm) values->other_parm=(parameters **)realloc(values->other_parm,(values->number_of_input_alignments)*sizeof(parameters *));
    				else values->other_parm=(parameters **)malloc((values->number_of_input_alignments)*sizeof(parameters *));
    				assert((int) values->other_parm);
    				values->other_parm[values->number_of_input_alignments-1]=NULL;
    				/*fprintf(stderr,"Infile designations must be followed by 'param file_name' to specify parameter file.");
	    			exit(-1);*/
		    		}
		    }
		    found=1;
		    }
	       if (found==0)  {
		    fprintf(stderr,"Unrecognized option - (%s) in %s!\n",s,file_name);
		    exit(-1);
		    }
	       }
	  fclose(fp_in);
	  }
     else {
	  fprintf(stderr,"No such alignment file!\n");
	  exit(-1);
	  }
}

void get_file(file_name,values)
char *file_name;
parameters *values;
{
     char s[100],cval[100],found,cval2[100],temp_option[100];
     int val,i,j,check;
     FILE *fp_in, *fopen();

     fp_in=fopen(file_name,"r");
     temp_option[0]='\0';
     if (fp_in) {
	  while (!feof(fp_in)) {
		if (temp_option[0]=='\0') fscanf(fp_in,"%s ",s);
	       else {
		for (i=0;i<strlen(temp_option);i++) s[i]=temp_option[i];
		s[strlen(temp_option)]='\0';
		}
		  temp_option[0]='\0';
	       found=0;
	       if (strnicmp("intern",s,6)==0) {
		    fscanf(fp_in,"%d",&val);
		    values->gap_cost=(val);
		    found=1;
		    }
	       if (strnicmp("alignfi",s,7)==0) {
		    fscanf(fp_in,"%s ",cval);
		    values->align_file=(char *)malloc((1+strlen(cval))*sizeof(char));
		    assert((int) values->align_file);
		    values->align_file=(char *)strcpy(values->align_file,cval);
		    get_align_file(values->align_file,values);
		    found=1;
		    }
	       if (strnicmp("infile",s,6)==0) {
		if (!values->number_of_input_alignments) {
			values->data_sets=(int *)malloc(sizeof(int));
			 assert((int) values->data_sets);
			 values->data_set_weights=(int *)malloc(sizeof(int));
			 assert((int) values->data_set_weights);
			 }
		    else {
			 values->data_sets=(int *)realloc(values->data_sets,(1+values->number_of_input_alignments)*sizeof(int));
			 assert((int) values->data_sets);
			 values->data_set_weights=(int *)realloc(values->data_set_weights,(1+values->number_of_input_alignments)*sizeof(int));
			 assert((int) values->data_set_weights);
			 }
		    /*if (!values->new_optimization) {
			values->data_sets[values->number_of_input_alignments]=0;
			 values->data_set_weights[values->number_of_input_alignments]=1;
				}
		    else {*/
		    if(1) {
			 fscanf(fp_in,"%s ",cval2);
			 if (!strnicmp("a",cval2,1)) values->data_sets[values->number_of_input_alignments]=0;
			 else      if (!strnicmp("n",cval2,1)) values->data_sets[values->number_of_input_alignments]=1;
			 else      if (!strnicmp("m",cval2,1)) values->data_sets[values->number_of_input_alignments]=2;
			 else      if (!strnicmp("g",cval2,1)) values->data_sets[values->number_of_input_alignments]=3;
			 else {
			      fprintf(stderr,"Unrecognized option modifier (%s) after (infile) in file (%s).\n     Should be align, noalign, or morph when optalign specified.\n",cval2,file_name);
			      exit(-1);
			      }
			 fscanf(fp_in,"%d ",&val);
			 values->data_set_weights[values->number_of_input_alignments]=val;
			 }
		    fscanf(fp_in,"%s ",cval);
		    ++values->number_of_input_alignments;
		    values->input_file_name=(char **)realloc(values->input_file_name,values->number_of_input_alignments*sizeof(char *));
		    assert((int)values->input_file_name);
		    values->input_file_name[values->number_of_input_alignments-1]=(char *)malloc((1+strlen(cval))*sizeof(char));
		    assert((int)values->input_file_name[values->number_of_input_alignments-1]);
		    values->input_file_name[values->number_of_input_alignments-1]=(char *)strcpy(values->input_file_name[values->number_of_input_alignments-1],cval);
		    check=fscanf(fp_in,"%s ",cval);
		    if (check) {
    		    if (!strnicmp(cval,"param",5)) {
        		     fscanf(fp_in,"%s ",cval);
        			if (values->other_parm) values->other_parm=(parameters **)realloc(values->other_parm,(values->number_of_input_alignments)*sizeof(parameters *));
        			else values->other_parm=(parameters **)malloc((values->number_of_input_alignments)*sizeof(parameters *));
        			assert((int) values->other_parm);
        			values->other_parm[values->number_of_input_alignments-1]=(parameters *)malloc(sizeof(parameters));
        			assert((int) values->other_parm[values->number_of_input_alignments-1]);
        			get_initial_values(values->other_parm[values->number_of_input_alignments-1]);
        			get_file(cval,values->other_parm[values->number_of_input_alignments-1]);
        			fix_all_after_values2(values->other_parm[values->number_of_input_alignments-1]);
        			}
    		    else {
    				for (i=0;i<strlen(cval);i++) temp_option[i]=cval[i];
    				temp_option[strlen(cval)]='\0';
    				if (values->other_parm) values->other_parm=(parameters **)realloc(values->other_parm,(values->number_of_input_alignments)*sizeof(parameters *));
    				else values->other_parm=(parameters **)malloc((values->number_of_input_alignments)*sizeof(parameters *));
    				assert((int) values->other_parm);
    				values->other_parm[values->number_of_input_alignments-1]=NULL;
    				}
    			}
		found=1;
		    }
	       if (strnicmp("costo",s,5)==0) {
		    values->cost_only=1;
		    found=1;
		    }
	       /*
	       if (strnicmp("cladec",s,6)==0) {
		    values->check_cache_for_scores=1;
		    found=1;
		    }
		    */
	       if (strnicmp("reord",s,5)==0) {
		    values->reorder=1;
		    found=1;
		    }
	      if (strnicmp("poop",s,4)==0) {
		    values->poopsy=1;
		    found=1;
		    }
	        if (strnicmp("dalign",s,6)==0) {
		    values->do_dave_align=1;
		    found=1;
		    }
	       if (!strnicmp("noiter",s,6)) {
		    values->iter=0;
		    found=1;
		    }
	      if (!strnicmp("worst",s,5)) {
		    values->worst=1;
		    found=1;
		    }
	       if (!strnicmp("expand",s,6)) {
		    values->expand_X=1;
		    found=1;
		    }
	       if (!strnicmp("makeg",s,5)) {
		    values->make_groups_file=0;
		    fscanf(fp_in,"%s ",values->output_groups_name);
		    found=1;
		    }
	       if (!strnicmp("nolow",s,5)) {
		    values->low_mem=0;
		    found=1;
		    }
	       if (strnicmp("shortt",s,6)==0) {
		    values->shortcut_tree=1;
		    found=1;
		    }
	       if (strnicmp("short2",s,6)==0) {
		    values->shortcut_tree2=1;
		    found=1;
		    }

	       if (strnicmp("freqc",s,5)==0) {
		    values->freq_cost=1;
		    found=1;
		    }
	       if (strnicmp("hypancgap",s,9)==0) {
		    values->gap_in_hypanc=1;
		    found=1;
		    }
	       if (strnicmp("phylot",s,6)==0) {
		    values->phylo_time=0;
		    found=1;
		    }
	       if (strnicmp("optali",s,6)==0) {
		    values->new_optimization=1;
		    found=1;
		    }
	       if (strnicmp("opti",s,4)==0) {
		    values->optimize_nodes=1;
		    found=1;
		    }
	       if (strnicmp("apol",s,4)==0) {
		    values->apolist=1;
		    found=1;
		    }
	       if (strnicmp("spr",s,3)==0) {
		    values->sbr=1;
		    found=1;
		    }
	       if (strnicmp("tbr",s,3)==0) {
		    values->tbr=values->sbr=1;
		    found=1;
		    }
	       if (strnicmp("atbsh",s,5)==0) {
		    values->tbr_align_shortcut=1;
		    found=1;
		    }
	       if (strnicmp("optr",s,4)==0) {
		    values->rand_apo=1;
		    found=1;
		    }
		if (strnicmp("cainf",s,5)==0) {
		    values->cache_info=1;
		    found=1;
		    }
	       if (strnicmp("cull",s,4)==0) {
		    values->cull=1;
		    found=1;
		    }
	       if (strnicmp("elis",s,4)==0) {
		    values->elision=1;
		    found=1;
		    }
	       if (strnicmp("prefdi",s,6)==0) {
		    values->pref_direc=DIAGONAL;
		    found=1;
		    }
	       if (strnicmp("prefma",s,6)==0) {
		    values->pref_direc=DIAGONAL;
		    found=1;
		    }
	       if (strnicmp("prefdo",s,6)==0) {
		    values->pref_direc=DOWN;
		    found=1;
		    }
	       if (strnicmp("preflo",s,6)==0) {
		    values->pref_direc=DOWN;
		    found=1;
		    }
	       if (strnicmp("nolead",s,6)==0) {
		    values->no_leading_or_trailing_cost=1;
		    found=1;
		    }
	       if (strnicmp("prefri",s,6)==0) {
		    values->pref_direc=RIGHT;
		    found=1;
		    }
	       if (strnicmp("prefsh",s,6)==0) {
		    values->pref_direc=RIGHT;
		    found=1;
		    }
	       if (strnicmp("prefra",s,6)==0) {
		    values->pref_direc=RANDOM;
		    found=1;
		    }
	       if (strnicmp("inalign",s,7)==0) {
		    ++values->input_is_align;
		    found=1;
		    }
	       if (strnicmp("gapmustc",s,8)==0) {
		    values->gap_must_cost=1;
		    found=1;
		    }
	       if (strnicmp("maxgap",s,6)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->max_gap=(val);
		    found=1;
		    }
	       if (strnicmp("chop",s,4)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->chop=(val);
		    found=1;
		    }
	       if (strnicmp("jack",s,4)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->jackboot=(val);
		    found=1;
		    }
	       if (strnicmp("grains",s,6)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->grain_size=(val);
		    found=1;
		    }
	       if (strnicmp("proce",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->process_factor=(val);
		    found=1;
		    }
	       if (strnicmp("randor",s,6)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->rand_order=(val);
		    found=1;
		    }
	       if (strnicmp("treer",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->tree_rand_order_max=(val);
		    found=1;
		    }
	       if (strnicmp("cache",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
	values->cache_size=val;
		    found=1;
		    }
	       if (strnicmp("maxw",s,4)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->weight_range=(val);
		    found=1;
		    }
	       if (strnicmp("leadi",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->leading_gap_cost=(val);
		    found=1;
		    }
	       if (strnicmp("matq",s,4)==0) {
		    values->matquick=1;
		    found=1;
		    }
	       if (strnicmp("est",s,3)==0) {
		    values->get_weights_2=1;
	fscanf(fp_in,"%s ",cval);
	if (strnicmp("f",cval,1)==0) values->saw_farris=1;
	else if (strnicmp("g",cval,1)==0) {
	  values->saw_goloboff=1;
			 fscanf(fp_in,"%s ",cval);
	  if (isdigit(cval[0]))    values->saw_goloboff=atoi(cval);
	  else {
	       fprintf(stderr,"Goloboff must be followed by a value for (k).\n");
	    exit(-1);
	       }
	  }
	else if (strnicmp("w",cval,1)==0) values->saw_wheeler=1;
	else {
	  fprintf(stderr,"Option (est) must be followed by a valid type of estimation (Farris, Goloboff, or Wheeler).\n");
	  exit(-1);
	  }
		    found=1;
		    }
	       if (strnicmp("maxa",s,4)==0) {
		    values->max_out_aligns=1;
		    found=1;
		    }
	       if (strnicmp("alignn",s,6)==0) {
		    values->align_node_swap=1;
		    found=1;
		    }
	       if (strnicmp("aspr",s,4)==0) {
		    values->asbr=1;
		    found=1;
		    }
	       if (strnicmp("atbr",s,4)==0) {
		    values->asbr=values->atbr=1;
		    found=1;
		    }
	       if (strnicmp("ftbr",s,4)==0) {
		    values->tbr_first=1;
		    found=1;
		    }
	       if (strnicmp("arrt",s,4)==0) {
		    values->arrt=1;
		    found=1;
		    }
		 if (strnicmp("alignr",s,6)==0) {
		    values->align_root_swap=1;
		    found=1;
		    }
	       if (strnicmp("alignc",s,6)==0) {
		    values->align_complete_root_swap=1;
		    found=1;
		    }
	       if (strnicmp("alignp",s,6)==0) {
		    values->align_partial_root_swap=1;
		    found=1;
		    }
	       if (strnicmp("start",s,5)==0) {
		    values->groups_as_start=1;
		    found=1;
		    }
	       if (strnicmp("report",s,6)==0) {
		    values->dump_parameters=1;
		    found=1;
		    }
	       if (strnicmp("maxt",s,4)==0) {
		    values->max_out_trees=1;
		    found=1;
		    }
	       if (strnicmp("outo",s,4)==0) {
		    fscanf(fp_in,"%s ",s);
	       if (strnicmp("in",s,2)==0) {
			 values->output_order=1;
			 }
	       else if (strnicmp("al",s,2)==0) {
		    values->output_order=0;
			 }
		    found=1;
		    }
	       if (strnicmp("iter",s,4)==0) {
		    values->iter=1;
		    found=1;
		    }
	       if (strnicmp("contig",s,6)==0) {
		    values->min_num_gaps=1;
		    found=1;
		    }
	       if (strnicmp("discon",s,6)==0) {
		    values->min_num_gaps=2;
		    found=1;
		    }
	       if (strnicmp("lowm",s,4)==0) {
		    values->low_mem=1;
		    found=1;
		    }
	       if (strnicmp("trail",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->trailing_gap_cost=(val);
		    found=1;
		    }
	       if (strnicmp("coll",s,4)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->collapse=(val);
		    found=1;
		    }
	       if (strnicmp("givem",s,5)==0) {
		    values->givem=1;
		    found=1;
		    }
	       if (strnicmp("chang",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->change_cost=(val);
		    found=1;
		    }
	       if (strnicmp("ascii",s,5)==0) {
		    values->acgt=1;
		    found=1;
		    }
	       if (strnicmp("hen86",s,5)==0) {
		    values->farris=1;
		    values->max_name_length=10;
		    found=1;
		    }
	       if (strnicmp("nona",s,4)==0) {
		    values->nona=1;
		    values->max_name_length=10;
		    found=1;
		    }
	       if (strnicmp("clados",s,6)==0) {
		    values->clados=1;
		    found=1;
		    }
	       if (strnicmp("interl",s,6)==0) {
		    values->inter_dig=1;
		    found=1;
		    }
	       if (strnicmp("heng",s,4)==0)  {
		    values->hen_gap=1;
		    if (values->new_optimization) values->max_name_length=10;
		    found=1;
		    }
	       if (strnicmp("score",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->phylo_score=(val);
		    found=1;
		    }
	       if (strnicmp("randa",s,5)==0) {
		    values->rand_align=1;
		    fscanf(fp_in,"%d ",&val);
		    values->how_many=(val);
		    found=1;
		    }
	       if (strnicmp("keept",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->keep_trees=(val);
		    found=1;
		    }
	       if (strnicmp("length",s,6)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->length_at_end=(val);
		    found=1;
		    }
	       if (strnicmp("keepa",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->keep_aligns=(val);
		    found=1;
		    }
	       if (strnicmp("linel",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->line_length=(val);
		    found=1;
		    }
	       if (strnicmp("cutp",s,4)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->cutpoint=(val);
		    found=1;
		    }
	       if (strnicmp("randt",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->randtrees=(val);
		    found=1;
		    }
	       if (strnicmp("extrag",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->extra_adjustment=(val);
		    found=1;
		    }
	       if (strnicmp("coding",s,6)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->coding=(val);
		    found=1;
		    }
	       if (strnicmp("alignm",s,6)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->align_multi_swap=(val);
		    found=1;
		    }
	       if (strnicmp("treem",s,5)==0) {
		    fscanf(fp_in,"%d ",&val);
		    values->tree_multi_swap=(val);
		    found=1;
		    }
	       if (strnicmp("pair",s,4)==0)  {
		    values->get_heur=1;
		    found=1;
		    }
	       if (strnicmp("exact",s,5)==0) {
		    values->best=1;
		    found=1;
		    }
	       if (strnicmp("best",s,4)==0) {
		    values->best_order=1;
		    found=1;
		    }
	       if (strnicmp("build",s,5)==0) {
		    values->get_heur2=1;
		    found=1;
		    }
	       if (strnicmp("trees",s,5)==0) {
		    values->tree_swap=1;
		    found=1;
		    }
	       if (strnicmp("aligns",s,6)==0) {
		    values->align_swap=1;
		    found=1;
		    }
	       if (strnicmp("dot",s,3)==0)   {
		    values->dot=1;
		    found=1;
		    }
	       if (strnicmp("quick",s,5)==0) {
		    values->aquick=values->get_heur2=1;
		    found=1;
		    }
	       if (strnicmp("newq",s,4)==0) {
		    values->get_heur2=values->get_heur3=1;
		    found=1;
		    }
	       if (strnicmp("nogap",s,5)==0) {
		    values->phylo_gap=0;
		    found=1;
		    }
	       if (strnicmp("tacit",s,5)==0) {
		    values->VERBOSE=0;
		    found=1;
		    }
	       if (strnicmp("noerr",s,5)==0) {
		    values->rep_error=0;
		    found=1;
		    }
	       if (strnicmp("grainma",s,7)==0)    {
		    values->grain_size=HUGE_COST;
		    found=1;
		    }
	       if (strnicmp("grainmi",s,7)==0)    {
		    values->grain_size=1;
		    found=1;
		    }
	       if (strnicmp("silent",s,6)==0) {
		    values->VERBOSE=values->rep_error=0;
		    found=1;
		    }
	       if (strnicmp("showm",s,5)==0) {
		    values->show_mem=1;
		    found=1;
		    }
	       if (strnicmp("printi",s,6)==0) {
		    values->print_intermediates=1;
		    found=1;
		    }
	       if (strnicmp("paup",s,4)==0)  {
		    values->paup=1;
		    found=1;
		    }
	       if (strnicmp("pdot",s,4)==0)  {
		    values->paup_dot=1;
		    found=1;
		    }
	       if (strnicmp("time",s,4)==0)  {
		    values->time=1;
		    found=1;
		    }
	       if (strnicmp("aligna",s,6)==0)     {
		    values->swap_while_add=1;
		    found=1;
		    }
	       if (strnicmp("treea",s,5)==0) {
		    values->clade_swap_while_add=1;
		    found=1;
		    }
	       if (strnicmp("newc",s,4)==0)  {
		    if (values->new_codes) {
			 values->n_codes+=1;
			 values->new_codes=(char **)realloc(values->new_codes,values->n_codes*sizeof(char *));
			 assert((int)values->new_codes);
			 values->new_codes[values->n_codes-1]=(char *)malloc(4*sizeof(char));
			 assert((int)values->new_codes[values->n_codes-1]);
			 fscanf(fp_in,"%s ",cval);
			 values->new_codes[values->n_codes-1][0]=cval[0];
			 fscanf(fp_in,"%s ",cval);
			 values->new_codes[values->n_codes-1][1]=cval[0];
			 values->new_codes[values->n_codes-1][2]=cval[1];
			 values->new_codes[values->n_codes-1][3]=cval[2];
			 }
		    if (!values->new_codes) {
			 values->new_codes=(char **)malloc(1*sizeof(char *));
			 assert((int)values->new_codes);
			 values->new_codes[0]=(char *)malloc(4*sizeof(char));
			 assert((int)values->new_codes[0]);
			 values->n_codes=1;
			 fscanf(fp_in,"%s ",cval);
			 values->new_codes[0][0]=cval[0];
			 fscanf(fp_in,"%s ",cval);
			 values->new_codes[0][1]=cval[0];
			 values->new_codes[0][2]=cval[1];
			 values->new_codes[0][3]=cval[2];
			 }
		    found=1;
		    }

	       if (strnicmp("matr",s,4)==0) {
		    if (!values->delta) {
			 values->delta=(int **)malloc(4*sizeof(int*));
			 assert((int)values->delta);
			 for (i=0;i<4;i++) {
			      values->delta[i]=(int *)malloc(4*sizeof(int));
			      assert((int)values->delta[i]);
			      }
			 for (j=0;j<4;j++) {
			      for (i=0;i<4;i++) {
				   fscanf(fp_in,"%d ",&val);
				   values->delta[j][i]=(val);
				   }
			      }
			 }
		    found=1;
		    }
	       if (found==0)  {
		    fprintf(stderr,"Unrecognized option - (%s) in %s!\n",s,file_name);
		    exit(-1);
		    }
	       }
	  fclose(fp_in);
	  }
     else {
	  fprintf(stderr,"No such parameter file %s!\n",file_name);
	  exit(-1);
	  }
/*print_values(values);*/

}

void get_command_line(n_coms,commands,values,groups_name,check,fp_input,fp_output)
parameters *values;
int n_coms;
char *commands[],*groups_name,*check;
char *fp_input, *fp_output;
{
int i,j,k;

     char found;

     *check='n';

     for (i=1;i<n_coms;i++) {
	       found=0;
	       if (stricmp("input",commands[i])==0) {
		    fp_input=(char *)strcpy(fp_input,commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (stricmp("output",commands[i])==0) {
		    fp_output=(char *)strcpy(fp_output,commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("param",commands[i],5)==0) {
		    get_file(commands[i+1],values);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("group",commands[i],5)==0) {
		    groups_name=(char *)strcpy(groups_name,commands[i+1]);
		    *check='y';
		    found=1;
		    ++i;
		    }
	       if (strnicmp("cull",commands[i],4)==0) {
		    values->cull=1;
		    found=1;
		    }
	       if (strnicmp("poop",commands[i],4)==0) {
		    values->poopsy=1;
		    found=1;
		    }
	       if (strnicmp("elis",commands[i],4)==0) {
		    values->elision=1;
		    found=1;
		    }
	       if (strnicmp("intern",commands[i],6)==0) {
		    values->gap_cost=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("maxgap",commands[i],6)==0) {
		    values->max_gap=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("chop",commands[i],4)==0) {
		    values->chop=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("jack",commands[i],4)==0) {
		    values->jackboot=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("grains",commands[i],6)==0) {
		    values->grain_size=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("proce",commands[i],5)==0) {
		    values->process_factor=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("randor",commands[i],6)==0) {
		    values->rand_order=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("treer",commands[i],5)==0) {
		    values->tree_rand_order_max=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("maxw",commands[i],4)==0) {
		    values->weight_range=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("alignfi",commands[i],7)==0) {
		    values->align_file=(char *)malloc((1+strlen(commands[i+1]))*sizeof(char));
		    assert((int) values->align_file);
		    values->align_file=(char *)strcpy(values->align_file,commands[i+1]);
		    get_align_file(values->align_file,values);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("infile",commands[i],6)==0) {
		if (!values->number_of_input_alignments) {
			 values->data_sets=(int *)malloc((1+values->number_of_input_alignments)*sizeof(int));
			 assert((int) values->data_sets);
			 values->data_set_weights=(int *)malloc((1+values->number_of_input_alignments)*sizeof(int));
			 assert((int) values->data_set_weights);
			 }
		    else {
			 values->data_sets=(int *)realloc(values->data_sets,(1+values->number_of_input_alignments)*sizeof(int));
			 assert((int) values->data_sets);
			 values->data_set_weights=(int *)realloc(values->data_set_weights,(1+values->number_of_input_alignments)*sizeof(int));
			 assert((int) values->data_set_weights);
			 }
		  /*if (!values->new_optimization) {
			values->data_sets[values->number_of_input_alignments]=0;
			 values->data_set_weights[values->number_of_input_alignments]=1;
				}
		   else {*/
		   if (1) {
			 ++i;
			 if (!strnicmp("a",commands[i],1)) values->data_sets[values->number_of_input_alignments]=0;
			 else      if (!strnicmp("n",commands[i],1)) values->data_sets[values->number_of_input_alignments]=1;
			 else      if (!strnicmp("m",commands[i],1)) values->data_sets[values->number_of_input_alignments]=2;
			 else      if (!strnicmp("g",commands[i],1)) values->data_sets[values->number_of_input_alignments]=3;
			 else {
			      fprintf(stderr,"Unrecognized option modifier (%s) after (infile) in command line.\n  Should be align, noalign, or morph when optalign specified.\n",commands[i]);
			      exit(-1);
			      }
			 ++i;
			 values->data_set_weights[values->number_of_input_alignments]=atoi(commands[i]);
			 }
		    ++values->number_of_input_alignments;
		    values->input_file_name=(char **)realloc(values->input_file_name,values->number_of_input_alignments*sizeof(char *));
		    assert((int)values->input_file_name);
		    values->input_file_name[values->number_of_input_alignments-1]=(char *)malloc((1+strlen(commands[i+1]))*sizeof(char));
		    assert((int)values->input_file_name[values->number_of_input_alignments-1]);
		    values->input_file_name[values->number_of_input_alignments-1]=(char *)strcpy(values->input_file_name[values->number_of_input_alignments-1],commands[i+1]);
		    ++i;
		    if ((i+1)<n_coms) {
    		    if (!strnicmp(commands[++i],"param",5)) {
    		        if (values->other_parm) values->other_parm=(parameters **)realloc(values->other_parm,(values->number_of_input_alignments)*sizeof(parameters *));
        			else values->other_parm=(parameters **)malloc((values->number_of_input_alignments)*sizeof(parameters *));
        			assert((int) values->other_parm);
        			values->other_parm[values->number_of_input_alignments-1]=(parameters *)malloc(sizeof(parameters));
        			assert((int) values->other_parm[values->number_of_input_alignments-1]);
        			get_initial_values(values->other_parm[values->number_of_input_alignments-1]);
        			get_file(commands[++i],values->other_parm[values->number_of_input_alignments-1]);
        			fix_all_after_values2(values->other_parm[values->number_of_input_alignments-1]);
        			}
    		    else {
    	    		if (values->other_parm) values->other_parm=(parameters **)realloc(values->other_parm,(values->number_of_input_alignments)*sizeof(parameters *));
        			else values->other_parm=(parameters **)malloc((values->number_of_input_alignments)*sizeof(parameters *));
        			assert((int) values->other_parm);
        			values->other_parm[values->number_of_input_alignments-1]=NULL;
        			--i;
        			}
        	}
		   found=1;
		    }
	       if (strnicmp("leadi",commands[i],5)==0) {
		    values->leading_gap_cost=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
       if (strnicmp("cache",commands[i],5)==0) {
			values->cache_size=atoi(commands[i+1]);
		    found=1;
			++i;
		    }
       if (strnicmp("trail",commands[i],5)==0) {
		    values->trailing_gap_cost=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
       if (strnicmp("coll",commands[i],4)==0) {
		    values->collapse=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
       if (strnicmp("iter",commands[i],4)==0) {
		    values->iter=1;
		    found=1;
		    }
       if (strnicmp("expand",commands[i],6)==0) {
		    values->expand_X=1;
		    found=1;
		    }
       if (strnicmp("makeg",commands[i],5)==0) {
		    values->make_groups_file=1;
		    found=1;
		    strcpy(values->output_groups_name,commands[i+1]);
		    i++;
		    }
	if (!strnicmp("noiter",commands[i],6)) {
		    values->iter=0;
		    found=1;
		    }
	if (!strnicmp("worst",commands[i],5)) {
		    values->worst=1;
		    found=1;
		    }
	  if (!strnicmp("nolowm",commands[i],6)) {
		    values->low_mem=0;
		    found=1;
		    }
	  if (strnicmp("alignn",commands[i],6)==0) {
		    values->align_node_swap=1;
		    found=1;
		    }
	  if (strnicmp("alignr",commands[i],6)==0) {
		    values->align_root_swap=1;
		    found=1;
		    }
	  if (strnicmp("alignc",commands[i],6)==0) {
		    values->align_complete_root_swap=1;
		    found=1;
		    }
	  if (strnicmp("alignp",commands[i],6)==0) {
		    values->align_partial_root_swap=1;
		    found=1;
		    }
	  if (strnicmp("start",commands[i],5)==0) {
		    values->groups_as_start=1;
		    found=1;
		    }
	  if (strnicmp("report",commands[i],6)==0) {
		    values->dump_parameters=1;
		    found=1;
		    }
	  if (strnicmp("gapmustc",commands[i],8)==0) {
		    values->gap_must_cost=1;
		    found=1;
		    }
	       if (strnicmp("phylot",commands[i],6)==0) {
		    values->phylo_time=0;
		    found=1;
		    }
	       if (strnicmp("optali",commands[i],6)==0) {
		    values->new_optimization=1;
		    found=1;
		    }
	       if (strnicmp("opti",commands[i],4)==0) {
		    values->optimize_nodes=1;
		    found=1;
		    }
	       if (strnicmp("apol",commands[i],4)==0) {
		    values->apolist=1;
		    found=1;
		    }
	       if (strnicmp("spr",commands[i],3)==0) {
		    values->sbr=1;
		    found=1;
		    }
	       if (strnicmp("tbr",commands[i],3)==0) {
		    values->tbr=values->sbr=1;
		    found=1;
		    }
	       if (strnicmp("atbsh",commands[i],5)==0) {
		    values->tbr_align_shortcut=1;
		    found=1;
		    }
	       if (strnicmp("optr",commands[i],4)==0) {
		    values->rand_apo=1;
		    found=1;
		    }
	       if (strnicmp("inalign",commands[i],7)==0) {
		    ++values->input_is_align;
		    found=1;
		    }
	       if (strnicmp("cainf",commands[i],5)==0) {
		    values->cache_info=1;
		    found=1;
		    }
	       if (strnicmp("contig",commands[i],6)==0) {
		    values->min_num_gaps=1;
		    found=1;
		    }
	       if (strnicmp("discon",commands[i],6)==0) {
		    values->min_num_gaps=2;
		    found=1;
		    }
	       if (strnicmp("est",commands[i],3)==0) {
		    values->get_weights_2=1;
	if (strnicmp("f",commands[i+1],1)==0) values->saw_farris=1;
	else if (strnicmp("w",commands[i+1],1)==0) values->saw_wheeler=1;
	else if (strnicmp("g",commands[i+1],1)==0) {
	  values->saw_goloboff=1;
	  if (isdigit(commands[i+2][0])) {
	       values->saw_goloboff=atoi(commands[i+2]);
	    ++i;
	    }
	  else {
	       fprintf(stderr,"Goloboff must be followed by a value for (k).\n");
	    exit(-1);
	       }
	  }
	else {
	  fprintf(stderr,"Option (est) must be followed by a valid type of estimation (Farris, Goloboff, or Wheeler).\n");
	  exit(-1);
	  }
	++i;
		    found=1;
		    }
	       if (strnicmp("maxa",commands[i],4)==0) {
		    values->max_out_aligns=1;
		    found=1;
		    }
	       if (strnicmp("maxt",commands[i],4)==0) {
		    values->max_out_trees=1;
		    found=1;
		    }
	       /*
	       if (strnicmp("cladec",commands[i],6)==0) {
		    values->check_cache_for_scores=1;
		    found=1;
		    }
		    */
	       if (strnicmp("costo",commands[i],5)==0) {
		    values->cost_only=1;
		    found=1;
		    }
	       if (strnicmp("reord",commands[i],5)==0) {
		    values->reorder=1;
		    found=1;
		    }
	       if (strnicmp("dalign",commands[i],6)==0) {
		    values->do_dave_align=1;
		    found=1;
		    }
	       if (strnicmp("shortt",commands[i],6)==0) {
		    values->shortcut_tree=1;
		    found=1;
		    }
	       if (strnicmp("short2",commands[i],6)==0) {
		    values->shortcut_tree2=1;
		    found=1;
		    }

	       if (strnicmp("freqc",commands[i],5)==0) {
		    values->freq_cost=1;
		    found=1;
		    }
	       if (strnicmp("hypancgap",commands[i],9)==0) {
		    values->gap_in_hypanc=1;
		    found=1;
		    }
	       if (strnicmp("nolead",commands[i],6)==0) {
		    values->no_leading_or_trailing_cost=1;
		    found=1;
		    }
	       if (strnicmp("preflo",commands[i],6)==0) {
		    values->pref_direc=DOWN;
		    found=1;
		    }
	       if (strnicmp("prefma",commands[i],6)==0) {
		    values->pref_direc=DIAGONAL;
		    found=1;
		    }
			if (strnicmp("prefsh",commands[i],6)==0) {
		    values->pref_direc=RIGHT;
		    found=1;
		    }
	       if (strnicmp("prefdo",commands[i],6)==0) {
		    values->pref_direc=DOWN;
		    found=1;
		    }
	       if (strnicmp("prefdi",commands[i],6)==0) {
		    values->pref_direc=DIAGONAL;
		    found=1;
		    }
	       if (strnicmp("prefra",commands[i],6)==0) {
		    values->pref_direc=RANDOM;
		    found=1;
		    }
	       if (strnicmp("prefri",commands[i],6)==0) {
		    values->pref_direc=RIGHT;
		    found=1;
		    }
	       if (strnicmp("chang",commands[i],5)==0) {
		    values->change_cost=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("ascii",commands[i],5)==0) {
		    values->acgt=1;
		    found=1;
		    }
	       if (strnicmp("hen86",commands[i],5)==0) {
		    values->farris=1;
		    values->max_name_length=10;
		    found=1;
		    }
	       if (strnicmp("nona",commands[i],4)==0) {
		    values->nona=1;
		    values->max_name_length=10;
		    found=1;
		    }
	       if (strnicmp("clados",commands[i],6)==0) {
		    values->clados=1;
		    found=1;
		    }
	       if (strnicmp("interl",commands[i],6)==0) {
		    values->inter_dig=1;
		    found=1;
		    }
	       if (strnicmp("heng",commands[i],4)==0)  {
		    values->hen_gap=1;
		    if (values->new_optimization) values->max_name_length=10;
		    found=1;
		    }
	       if (strnicmp("lowm",commands[i],4)==0)  {
		    values->low_mem=1;
		    found=1;
		    }
	       if (strnicmp("showm",commands[i],5)==0) {
		    values->show_mem=1;
		    found=1;
		    }
	       if (strnicmp("score",commands[i],5)==0) {
		    values->phylo_score=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("outo",commands[i],4)==0) {
		    if (!(strnicmp("in",commands[i+1],2))) {
			 values->output_order=1;
			 found=1;
			 }
		    else if (!(strnicmp("al",commands[i+1],2))) {
			 values->output_order=0;
			 found=1;
			 }
		    ++i;
		    }
	       if (strnicmp("randa",commands[i],5)==0) {
		    values->rand_align=1;
		    values->how_many=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("keept",commands[i],5)==0) {
		    values->keep_trees=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("length",commands[i],6)==0) {
		    values->length_at_end=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("keepa",commands[i],5)==0) {
		    values->keep_aligns=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("linel",commands[i],5)==0) {
		    values->line_length=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("cutp",commands[i],4)==0) {
		    values->cutpoint=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("randt",commands[i],5)==0) {
		    values->randtrees=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("extrag",commands[i],5)==0) {
		    values->extra_adjustment=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("coding",commands[i],6)==0) {
		    values->coding=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("alignm",commands[i],6)==0) {
		    values->align_multi_swap=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("treem",commands[i],5)==0) {
		    values->tree_multi_swap=atoi(commands[i+1]);
		    found=1;
		    ++i;
		    }
	       if (strnicmp("pair",commands[i],4)==0)  {
		    values->get_heur=1;
		    found=1;
		    }
	       if (strnicmp("exact",commands[i],5)==0) {
		    values->best=1;
		    found=1;
		    }
	       if (strnicmp("best",commands[i],4)==0) {
		    values->best_order=1;
		    found=1;
		    }
	       if (strnicmp("build",commands[i],5)==0) {
		    values->get_heur2=1;
		    found=1;
		    }
	       if (strnicmp("trees",commands[i],5)==0) {
		    values->tree_swap=1;
		    found=1;
		    }
	       if (strnicmp("aligns",commands[i],6)==0) {
		    values->align_swap=1;
		    found=1;
		    }
	       if (strnicmp("aspr",commands[i],4)==0) {
		    values->asbr=1;
		    found=1;
		    }
	       if (strnicmp("atbr",commands[i],4)==0) {
		    values->asbr=values->atbr=1;
		    found=1;
		    }
	       if (strnicmp("ftbr",commands[i],4)==0) {
		    values->tbr_first=1;
		    found=1;
		    }
	       if (strnicmp("arrt",commands[i],4)==0) {
		    values->arrt=1;
		    found=1;
		    }
	       if (strnicmp("dot",commands[i],3)==0)   {
		    values->dot=1;
		    found=1;
		    }
	       if (strnicmp("quick",commands[i],5)==0) {
		    values->aquick=values->get_heur2=1;
		    found=1;
		    }
	       if (strnicmp("newq",commands[i],4)==0) {
		    values->get_heur2=values->get_heur3=1;
		    found=1;
		    }
	       if (strnicmp("noqui",commands[i],5)==0) {
		    values->aquick=values->get_heur2=0;
		    found=1;
		    }
	       if (strnicmp("nogap",commands[i],5)==0) {
		    values->phylo_gap=0;
		    found=1;
		    }
	       if (strnicmp("tacit",commands[i],5)==0) {
		    values->VERBOSE=0;
		    found=1;
		    }
	       if (strnicmp("noerr",commands[i],5)==0) {
		    values->rep_error=0;
		    found=1;
		    }
	       if (strnicmp("grainma",commands[i],7)==0)    {
		    values->grain_size=HUGE_COST;
		    found=1;
		    }
	       if (strnicmp("grainmi",commands[i],7)==0)    {
		    values->grain_size=1;
		    found=1;
		    }
	       if (strnicmp("silent",commands[i],6)==0) {
		    values->VERBOSE=values->rep_error=0;
		    found=1;
		    }
	       if (strnicmp("printi",commands[i],6)==0) {
		    values->print_intermediates=1;
		    found=1;
		    }
	       if (strnicmp("paup",commands[i],4)==0)  {
		    values->paup=1;
		    found=1;
		    }
	       if (strnicmp("pdot",commands[i],4)==0)  {
		    values->paup_dot=1;
		    found=1;
		    }
	       if (strnicmp("time",commands[i],4)==0)  {
		    values->time=1;
		    found=1;
		    }
	       if (strnicmp("aligna",commands[i],6)==0)     {
		    values->swap_while_add=1;
		    found=1;
		    }
	       if (strnicmp("treea",commands[i],5)==0) {
		    values->clade_swap_while_add=1;
		    found=1;
		    }
	       if (strnicmp("newc",commands[i],4)==0)  {
		    if (values->new_codes) {
			 values->n_codes+=1;
			 values->new_codes=(char **)realloc(values->new_codes,values->n_codes*sizeof(char *));
			 assert((int)values->new_codes);
			 values->new_codes[values->n_codes-1]=(char *)malloc(4*sizeof(char));
			 assert((int)values->new_codes[values->n_codes-1]);
			 values->new_codes[values->n_codes-1][0]=commands[i+1][0];
			 values->new_codes[values->n_codes-1][1]=commands[i+2][0];
			 values->new_codes[values->n_codes-1][2]=commands[i+2][1];
			 values->new_codes[values->n_codes-1][3]=commands[i+2][2];
			 }
		    if (!values->new_codes) {
			 values->new_codes=(char **)malloc(1*sizeof(char *));
			 assert((int)values->new_codes);
			 values->new_codes[0]=(char *)malloc(4*sizeof(char));
			 assert((int)values->new_codes[0]);
			 values->n_codes=1;
			 values->new_codes[0][0]=commands[i+1][0];
			 values->new_codes[0][1]=commands[i+2][0];
			 values->new_codes[0][2]=commands[i+2][1];
			 values->new_codes[0][3]=commands[i+2][2];
			 }
		    i+=2;
		    found=1;
		    }
	       if (strnicmp("matr",commands[i],4)==0) {
		    if (!values->delta) {
			 values->delta=(int **)malloc(4*sizeof(int*));
			 assert((int)values->delta);
			 for (k=0;k<4;k++) {
			      values->delta[k]=(int *)malloc(4*sizeof(int));
			      assert((int)values->delta[k]);
			      }
			 for (j=0;j<4;j++) {
			      for (k=0;k<4;k++) values->delta[j][k]=atoi(commands[i+1+(j*4)+k]);
			      }
			 i+=16;
			 }
		    found=1;
		    i+=1;
		    }
	       if (found==0)  {
		    fprintf(stderr,"Unrecognized option - (%s) in command line!\n",commands[i]);
		    exit(-1);
		    }
	       }
/*print_values(values);*/

}

int strnicmp_mine(a,b,n)
int n;
char *a,*b;
{
int i,to_return;
char *aa, *bb;

aa=(char *)malloc((strlen(a)+1)*sizeof(char));
bb=(char *)malloc((strlen(b)+1)*sizeof(char));
if ((!aa) || (!bb)) {
     fprintf(stderr,"Problems..bye\n");
     exit(-1);
     }
strcpy(aa,a);
strcpy(bb,b);

for (i=0;i<strlen(aa);i++) if (isupper(aa[i])) aa[i]=tolower(aa[i]);
for (i=0;i<strlen(bb);i++) if (isupper(bb[i])) bb[i]=tolower(bb[i]);

to_return=strncmp(aa,bb,n);
free(aa);
free(bb);
return to_return;
}

int stricmp_mine(a,b)
char *a,*b;
{
int i,to_return;
char *aa, *bb;

aa=(char *)malloc((strlen(a)+1)*sizeof(char));
bb=(char *)malloc((strlen(b)+1)*sizeof(char));
if ((!aa) || (!bb)) {
     fprintf(stderr,"Problems..bye\n");
     exit(-1);
     }

strcpy(aa,a);
strcpy(bb,b);

for (i=0;i<strlen(aa);i++) if (isupper(aa[i])) aa[i]=tolower(aa[i]);
for (i=0;i<strlen(bb);i++) if (isupper(bb[i])) bb[i]=tolower(bb[i]);

to_return=strcmp(aa,bb);
free(aa);
free(bb);
return to_return;
}

void get_initial_values(values)
parameters *values;
{
     bitchar[0]=0;
     bitchar[A_BIT]='A';
     bitchar[C_BIT]='C';
     bitchar[G_BIT]='G';
     bitchar[T_BIT]='T';
     bitchar[U_BIT]='T';
     bitchar[N_BIT]='N';
     bitchar[X_BIT]='X';
     bitchar[GAP_BIT]='-';
     bitchar[R_BIT]='R';
     bitchar[Y_BIT]='Y';
     bitchar[M_BIT]='M';
     bitchar[W_BIT]='W';
     bitchar[S_BIT]='S';
     bitchar[K_BIT]='K';
     bitchar[B_BIT]='B';
     bitchar[D_BIT]='D';
     bitchar[H_BIT]='H';
     bitchar[V_BIT]='V';
     bitchar[E_BIT]='E';
     bitchar[F_BIT]='F';
     bitchar[I_BIT]='I';
     bitchar[J_BIT]='J';
     bitchar[L_BIT]='L';
     bitchar[O_BIT]='O';
     bitchar[P_BIT]='P';
     bitchar[Q_BIT]='Q';
     bitchar[Z_BIT]='Z';
     bitchar[C1_BIT]='1';
     bitchar[C2_BIT]='2';
     bitchar[C3_BIT]='3';
     bitchar[C4_BIT]='4';
     bitchar[C5_BIT]='5';
     bitchar[C6_BIT]='6';
     bitchar[C7_BIT]='7';
     bitchar[C8_BIT]='8';

     charbit['A']=A_BIT;
     charbit['C']=C_BIT;
     charbit['G']=G_BIT;
     charbit['T']=T_BIT;
     charbit['U']=U_BIT;
     charbit['N']=N_BIT;
     charbit['X']=X_BIT;
     charbit['-']=GAP_BIT;
     charbit['R']=R_BIT;
     charbit['Y']=Y_BIT;
     charbit['M']=M_BIT;
     charbit['W']=W_BIT;
     charbit['S']=S_BIT;
     charbit['K']=K_BIT;
     charbit['B']=B_BIT;
     charbit['D']=D_BIT;
     charbit['H']=H_BIT;
     charbit['V']=V_BIT;

     charbit['E']=E_BIT;
     charbit['F']=F_BIT;
     charbit['I']=I_BIT;
     charbit['J']=J_BIT;
     charbit['L']=L_BIT;
     charbit['O']=O_BIT;
     charbit['P']=P_BIT;
     charbit['Q']=Q_BIT;
     charbit['Z']=Z_BIT;
     charbit['1']=C1_BIT;
     charbit['2']=C2_BIT;
     charbit['3']=C3_BIT;
     charbit['4']=C4_BIT;
     charbit['5']=C5_BIT;
     charbit['6']=C6_BIT;
     charbit['7']=C7_BIT;
     charbit['8']=C8_BIT;

     values->gap_cost = 4;
     values->leading_gap_cost = -32000;
     values->trailing_gap_cost = -32000;
     values->change_cost = 2;
     values->acgt = 0;
     values->farris = 0;
     values->inter_dig = 0;
     values->hen_gap = 0;
     values->phylo_score = 3;
     values->paup = 0;
     values->paup_dot = 0;
     values->how_many = 0;
     values->keep_trees = 100;
     values->length_at_end = 0;
     values->get_heur = 0;
     values->best= 0;
     values->get_heur2=0;
     values->rand_align=0;
     values->tree_swap=0;
     values->align_swap=0;
     values->dot=0;
     values->aquick=0;
     values->keep_aligns=100;
     values->line_length=60;
     values->phylo_gap=1;
     values->VERBOSE=1;
     values->rep_error=1;
     values->randtrees=0;
     values->print_intermediates=0;
     values->number=0;
     values->time=0;
     values->extra_adjustment=-32000;
     values->ngroups=0;
     values->groups=NULL;
     values->coding=-32000;
     values->in_some=0;
     values->align_multi_swap=0;
     values->tree_multi_swap=0;
     values->delta=NULL;
     values->all_made=NULL;
     values->tree_made=NULL;
     values->new_codes=NULL;
     values->n_codes=0;
     values->number_best_clad=0;
     values->length_best_clad=HUGE_COST;
     values->best_rep=NULL;
     values->ttr=0;
     values->transition=0;
     values->transversion=0;
     values->n_transv=0;
     values->givem=0;
     values->matquick=0;
     values->max_gap=MAX_SEQUENCE_SIZE;
     values->iter=1;
     values->low_mem=0;
     values->swap_while_add=0;
     values->clade_swap_while_add=0;
     values->output_order=0;
  values->show_mem=0;
  values->max_out_aligns=0;
  values->max_out_trees=0;
  values->cost_only=0;
  values->temp_changes=0;
  values->temp_gaps=0;
  values->overall_changes=0;
  values->temp_ti=0;
  values->overall_ti=0;
  values->temp_tv=0;
  values->overall_tv=0;
  values->overall_gaps=0;
  values->temp_changes_g=0;
  values->temp_gaps_g=0;
  values->overall_changes_g=0;
  values->temp_ti_g=0;
  values->overall_ti_g=0;
  values->temp_tv_g=0;
  values->overall_tv_g=0;
  values->overall_gaps_g=0;
  values->get_weights=0;
  values->get_weights_2=0;
  values->weight_range=10;
  values->cache_size=0;
  values->align_cache=NULL;
  values->cache_info=0;
  values->saw_farris=0;
  values->saw_goloboff=0;
  values->saw_wheeler=0;
  values->input_is_align=0;
     values->b_string1=NULL;
     values->b_string2=NULL;
  values->phylo_time=-1;
  values->pref_direc=0;
  values->min_num_gaps=0;
  values->gap_must_cost=0;
  values->reorder=0;
  values->rand_order=0;
  values->tree_rand_order_max=0;
  values->number_of_input_alignments=0;
  values->input_file_name=NULL;
  values->align_file=NULL;
  values->actual_num_sequences=0;
  values->cull=0;
  values->elision=0;
  values->ce_weights=NULL;
  values->previous=NULL;
  /*
  values->check_cache_for_scores=0;
  */
  values->grain_size=1;
  /*values->num_hosts=0;
  values->tids=NULL;
  values->process_factor=1;*/ /*should not be reset*/
  values->new_optimization=0;
  values->no_leading_or_trailing_cost=0;
  values->gap_in_hypanc=0;
  values->data_sets=NULL;
  values->data_set_weights=NULL;
  values->align_node_swap=0;
  values->align_root_swap=0;
  values->align_complete_root_swap=0;
  values->align_partial_root_swap=0;
  values->groups_as_start=0;
  values->dump_parameters=0;
  values->in_bandb_loop=0;
  values->current_bound=HUGE_COST;
  values->chop=0;
  values->optimize_nodes=0;
  values->in_optimize=0;
  values->apolist=0;
  values->rand_apo=0;
  values->freq_cost=0;
  values->sbr=0;
  values->tbr=0;
  values->shortcut_tree=0;
  values->shortcut_tree2=1;
  values->get_heur3=0;
  values->asbr=0;
  values->atbr=0;
  values->tbr_align_shortcut=0;
  values->collapse=0;
  values->arrt=0;
  values->max_name_length=1000;
  values->do_dave_align=0;
  values->other_parm=NULL;
  values->jackboot=0;
  values->jack_name=NULL;
  values->jack_array=NULL;
  values->clados=0;
  values->cutpoint=50;
  values->nona=0;
  values->make_groups_file=0;
  values->start=NULL;
  values->stop=NULL;
  values->expand_X=0;
  values->best_order=0;
  values->poopsy=0;
  values->input_names=NULL;
  values->tbr_first=0;
values->worst=0;
  }

void fix_all_after_values(a,values)
parameters *values;
alignment **a;
{
int i,j,k;

if (values->worst){
 for (i=0;i<values->number_of_input_alignments;i++) values->data_set_weights[i]*= (-1);

}
if (values->best_order && values->groups) {
        fprintf(stderr,"Cannot have both 'best order' and specified groups.  Cancelling best order.\n");
        values->best_order=0;
        }
if (values->jackboot>0) {
	values->jack_name=(char **)malloc(values->jackboot*sizeof(char *));
	assert((int) values->jack_name);
	for (i=0;i<values->jackboot;i++) values->jack_name[i]=NULL;
	values->rand_order=values->jackboot;
	values->aquick=1;
	values->jack_array=(int *)malloc(1000*sizeof(int));
	assert((int) values->jack_array);
	if (values->cutpoint <50 ) {fprintf(stderr,"Jackboot cutoff too low reset to 50%%\n"); values->cutpoint=50;}
	else if (values->cutpoint > 100) {fprintf(stderr,"Jackboot cutoff too high reset to 100%%\n"); values->cutpoint=100;}
	}
if (values->cull) values->get_heur3=0;
else if (!values->new_optimization) {
	if (values->aquick || values->get_heur2) values->get_heur2=values->get_heur3=1;
	}
if (values->new_optimization) {
       values->do_dave_align=1;
       values->expand_X=1;
}
if (values->phylo_score > 0) {
    if (values->phylo_score%2) values->phylo_score=7;
    else values->phylo_score=8;
}
if (values->tree_swap) values->sbr=1;
if (values->align_swap || values->swap_while_add) values->asbr=1;
if (values->align_node_swap || values->align_root_swap || values->align_partial_root_swap || values->align_complete_root_swap) values->asbr=1;
if (values->align_root_swap) values->atbr=1;
if (values->align_partial_root_swap || values->align_complete_root_swap) values->arrt=1;

if ((values->rand_align) && (values->how_many<1)) {
     fprintf(stderr,"The value for 'randaligns' must be greater than '0.'\n");
     exit(-1);
     }
if (values->phylo_score==1) {
     values->randtrees=0;
     values->keep_trees=1;
     }
if (!values->atbr && values->tbr_align_shortcut) values->tbr_align_shortcut=0;
if ((values->phylo_score==3) && (values->randtrees==0)) values->keep_trees=1;
if ((values->phylo_score==0) && (values->randtrees==0)) values->keep_trees=1;
if ((values->phylo_score==5) && (values->randtrees==0)) values->keep_trees=1;
if ((values->phylo_score==7) && (values->randtrees==0)) values->keep_trees=1;
if ((values->aquick==1) && (values->rand_align==0) && (values->best==0)) values->keep_aligns=1;
if ((values->get_heur==1) && (values->rand_align==0) && (values->best==0) && (values->get_heur2==0)) {
     values->keep_aligns=1;
  values->cache_size=0;
  values->cache_info=0;
  }
if (values->get_heur3) values->cache_size=values->cache_info=0;
if ((values->phylo_score==2) && (values->randtrees<1)) {
     fprintf(stderr,"A value for 'randtrees' must be specified with random cladogram generation.\n");
     exit(-1);
     }
if ((values->phylo_score==3) && (values->tree_multi_swap==1)) values->tree_swap=1;

if (values->leading_gap_cost==-32000) values->leading_gap_cost=1;
     else values->leading_gap_cost=values->gap_cost - values->leading_gap_cost;
if (values->trailing_gap_cost==-32000) values->trailing_gap_cost=1;
     else values->trailing_gap_cost=values->gap_cost - values->trailing_gap_cost;

if (values->extra_adjustment!=-32000) values->extra_adjustment=values->gap_cost - values->extra_adjustment;
if (values->coding!=-32000) values->coding=values->gap_cost - values->coding;
if ((values->coding!=-32000) && (values->extra_adjustment!=-32000)) {
     fprintf(stderr,"Cannot specify both 'coding' and 'extragaps.'\n");
     exit(-1);
     }

if (values->delta) {
     for (i=0;i<4;i++){
	  for (j=i+1;j<4;j++) {
	       if (values->delta[i][j]!=values->delta[j][i]) {
			 fprintf(stderr,"Character transformation matrix must be symmetrical.\n");
			 exit(-1);
			 }
	       }
	  }
     }
for (i=0;i<values->n_codes;i++) {
     fprintf(stderr,"Amino acid symbol (%c) recoded as %c%c%c in protein\n\n",values->new_codes[i][0],values->new_codes[i][1],values->new_codes[i][2],values->new_codes[i][3]);
     }

/*Check for transition transversion situation*/
/*check for metricity*/
if (values->delta) {
     for (i=0;i<4;i++){
	  for (j=0;j<4;j++) if (values->delta[i][j]<0) {
	  fprintf(stderr,"Character transformation costs cannot be negative.\n");
     exit(-1);
	  }
	  for (j=i+1;j<4;j++) {
	       if (values->delta[i][j]!=values->delta[j][i]) {
			 fprintf(stderr,"Character transformation matrix is asymmetrical\n");
			 exit(-1);
			 }
	       }
	  }
     }
     /*metricity*/
if (values->delta) {
 if (values->delta[0][1]>(values->delta[0][2]+values->delta[1][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][1]>(values->delta[0][3]+values->delta[1][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][2]>(values->delta[0][1]+values->delta[1][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][2]>(values->delta[0][3]+values->delta[2][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][3]>(values->delta[0][1]+values->delta[1][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][3]>(values->delta[0][2]+values->delta[2][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][2]>(values->delta[1][0]+values->delta[0][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][2]>(values->delta[1][3]+values->delta[3][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][3]>(values->delta[1][0]+values->delta[0][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][3]>(values->delta[1][2]+values->delta[2][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[2][3]>(values->delta[2][0]+values->delta[0][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[2][3]>(values->delta[2][1]+values->delta[1][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
}
/*Check for transition transversion situation*/
if (values->delta) {
     if ((values->delta[0][1]==values->delta[0][3]) && (values->delta[0][1]==values->delta[1][2]) && (values->delta[0][1]==values->delta[2][3])) {
	       if (values->delta[0][2]==values->delta[1][3]) {
		    values->ttr=1;
			 values->transition=values->delta[0][2];
	       values->transversion=values->delta[0][1]-values->transition;
			 /*for (i=0;i<4;i++) free(values->delta[i]);
			 free(values->delta);
			 values->delta=NULL;*/
	       fprintf(stderr,"Matrix specifies only transitions and transversions.\n");
			 fprintf(stderr,"    Transversion cost %d Transition cost %d\n\n",values->transversion+values->transition,values->transition);
	       }
	  }
     }

if ((!values->VERBOSE || !values->rep_error) && values->show_mem) values->show_mem=0;
if (values->iter) values->max_gap=100;
if (values->cost_only) if (!values->VERBOSE) fprintf(stderr,"Determining cost of preexisting alignment.\n");
if (values->cost_only) values->get_heur=values->get_heur2=values->best=values->rand_align=0;
if ((values->tree_rand_order_max>0) || (values->rand_order>0) || (values->pref_direc==RANDOM) || (values->how_many>1) || (values->randtrees>1)) {
     srand((int) time(NULL));
  }
if (values->rand_order > 0) {
     if (values->rep_error) {
	  if (values->get_heur) if (values->rep_error) fprintf(stderr,"\nRandom orders will not affect pair.\n");
	  if (values->best) if (values->rep_error) fprintf(stderr,"\nRandom orders will not affect an exact solution.\n");
    }
     }

/*until Sankoff fixed
if (values->delta) {
     fprintf(stderr,"This option is currently disabled.  Please specify only transitions and transversion costs in matrix.\n");
  exit(-1);
     }*/

/*print file names
for (i=0;i<values->number_of_input_alignments;i++) {
     fprintf(stderr,"input alignment %d from file %s\n",i,values->input_file_name[i]);
     }*/
/*disable cache for multiple input alignments*/
if (values->number_of_input_alignments>0) if ((values->elision) || (values->cull))  {
     values->cache_size=0;
     values->cache_info=0;
     fprintf(stderr,"Cache is disabled when input data are alignments.\n");
     }
/*checks for cull and elision*/
if (values->cull || values->elision) {
     values->output_order=0;
     if (values->number_of_input_alignments<2) {
	  fprintf(stderr,"Cannot use options 'cull' or 'elision' unless multiple alignments are input\n");
	  exit(-1);
	  }
     for (i=1;i<values->number_of_input_alignments;i++) {
	  if (a[0]->n_seqs!=a[i]->n_seqs) {
	       fprintf(stderr,"Input alignments have different numbers of sequences--cannot use 'cull' or 'elision'\n");
	       exit(-1);
	       }
	  }
     }
if (values->cull && values->elision) {
     fprintf(stderr,"Cannot use both 'cull' and 'ellision' in the same run--'cull' is disabled\n    Rerun with 'cull' only to get both results\n");
     values->cull=0;
     }
if (values->elision) {
     values->get_heur=0;
     values->best=0;
     values->get_heur2=0;
     values->rand_align=0;
     if (values->hen_gap) {
	  values->hen_gap=0;
	  values->farris=1;
	  }
     }
if (values->cull || values->elision) {
     /*Reorder 1->n-1 */
     if (values->VERBOSE) fprintf(stderr,"Reordering alignments...");
     for (i=1;i<values->number_of_input_alignments;i++) new_order_as_input(a[i],a[0]->taxon_name,values);
     fprintf(stderr,"done\n");
     values->chop=0;
     }
/*
if ((values->check_cache_for_scores) && (!values->cache_size)) {
     fprintf(stderr,"Cannot cache cladogram scores without alignment cache.\n");
     values->check_cache_for_scores=0;
     }
     */
if (values->grain_size==(-1)) values->grain_size=((2*values->all_done)-1);
if (values->new_optimization) values->phylo_score=0;
if ((values->new_optimization && values->number_of_input_alignments) || (values->chop)) {
     if (values->number_of_input_alignments) fprintf(stderr,"Input data sets : %d\n",values->number_of_input_alignments);
     else fprintf(stderr,"Input data sets : %d\n",values->number_of_input_alignments+1);
     for (i=0;i<values->number_of_input_alignments;i++) {
	  fprintf(stderr,"    Data set %d with weight %d will ",i+1,values->data_set_weights[i]);
	  if ((values->data_sets[i]==0) || (values->data_sets[i]==3)) fprintf(stderr,"be aligned.\n");
	  else {
	       fprintf(stderr,"not be aligned");
	       if (values->data_sets[i]==1) fprintf(stderr," since it is already aligned.\n");
	       if (values->data_sets[i]==2) fprintf(stderr," since it consists of character data.\n");
	       }
	  }

     }
if (values->get_heur && values->new_optimization) {
	fprintf(stderr,"Cannot perform 'pair' with optimization alignments.\n");
	values->get_heur=0;
	}
if (values->apolist) {
	values->optimize_nodes=1;
	if (!values->new_optimization) {
		fprintf(stderr,"Cannot optimize alignments.\n");
		values->optimize_nodes=values->apolist=0;
		}
	}
get_lookup2(values);
/*modify change cost to reflect titv*/
k=0;
if (values->delta) {
	for (i=0;i<4;i++) for (j=0;j<4;j++) if (values->delta[i][j]>k) k=values->delta[i][j];
	values->change_cost=k;
	}
}


void print_values(values)
parameters *values;
{
     fprintf(stderr,"cost change   %d\n",values->change_cost);
     fprintf(stderr,"cost gap %d\n",values->gap_cost);
     fprintf(stderr,"cost leading  %d\n",values->leading_gap_cost);
     fprintf(stderr,"cost trailig  %d\n",values->trailing_gap_cost);
     fprintf(stderr,"phyloscore    %d\n",values->phylo_score);
     fprintf(stderr,"cost extra    %d\n",values->extra_adjustment);
     fprintf(stderr,"cost coding   %d\n",values->coding);

}


void fix_all_after_values2(values) /* I donh't know what this does differently*/
parameters *values;
{
int i,j;

if (values->worst) {
  for (i=0;i<values->number_of_input_alignments;i++) values->data_set_weights[i]*= (-1);
  ;
}
if (values->phylo_score>2) values->phylo_score=7;
if (values->tree_swap) values->sbr=1;
if (values->align_swap || values->swap_while_add) values->asbr=1;
if (values->align_node_swap || values->align_root_swap || values->align_partial_root_swap || values->align_complete_root_swap) values->asbr=1;
if (values->align_root_swap) values->atbr=1;
if (values->align_partial_root_swap || values->align_complete_root_swap) values->arrt=1;

if ((values->rand_align) && (values->how_many<1)) {
     fprintf(stderr,"The value for 'randaligns' must be greater than '0.'\n");
     exit(-1);
     }
if (values->phylo_score==1) {
     values->randtrees=0;
     values->keep_trees=1;
     }
if (!values->atbr && values->tbr_align_shortcut) values->tbr_align_shortcut=0;
if ((values->phylo_score==3) && (values->randtrees==0)) values->keep_trees=1;
if ((values->phylo_score==0) && (values->randtrees==0)) values->keep_trees=1;
if ((values->phylo_score==5) && (values->randtrees==0)) values->keep_trees=1;
if ((values->phylo_score==7) && (values->randtrees==0)) values->keep_trees=1;
if ((values->aquick==1) && (values->rand_align==0) && (values->best==0)) values->keep_aligns=1;
if ((values->get_heur==1) && (values->rand_align==0) && (values->best==0) && (values->get_heur2==0)) {
     values->keep_aligns=1;
  values->cache_size=0;
  values->cache_info=0;
  }
if (values->get_heur3) values->cache_size=values->cache_info=0;
if ((values->phylo_score==2) && (values->randtrees<1)) {
     fprintf(stderr,"A value for 'randtrees' must be specified with random cladogram generation.\n");
     exit(-1);
     }
if ((values->phylo_score==3) && (values->tree_multi_swap==1)) values->tree_swap=1;

if (values->leading_gap_cost==-32000) values->leading_gap_cost=1;
     else values->leading_gap_cost=values->gap_cost - values->leading_gap_cost;
if (values->trailing_gap_cost==-32000) values->trailing_gap_cost=1;
     else values->trailing_gap_cost=values->gap_cost - values->trailing_gap_cost;

if (values->extra_adjustment!=-32000) values->extra_adjustment=values->gap_cost - values->extra_adjustment;
if (values->coding!=-32000) values->coding=values->gap_cost - values->coding;
if ((values->coding!=-32000) && (values->extra_adjustment!=-32000)) {
     fprintf(stderr,"Cannot specify both 'coding' and 'extragaps.'\n");
     exit(-1);
     }

if (values->delta) {
     for (i=0;i<4;i++){
	  for (j=i+1;j<4;j++) {
	       if (values->delta[i][j]!=values->delta[j][i]) {
			 fprintf(stderr,"Character transformation matrix must be symmetrical.\n");
			 exit(-1);
			 }
	       }
	  }
     }
for (i=0;i<values->n_codes;i++) {
     fprintf(stderr,"Amino acid symbol (%c) recoded as %c%c%c in protein\n\n",values->new_codes[i][0],values->new_codes[i][1],values->new_codes[i][2],values->new_codes[i][3]);
     }

/*Check for transition transversion situation*/
/*check for metricity*/
if (values->delta) {
     for (i=0;i<4;i++){
	  for (j=0;j<4;j++) if (values->delta[i][j]<0) {
	  fprintf(stderr,"Character transformation costs cannot be negative.\n");
     exit(-1);
	  }
	  for (j=i+1;j<4;j++) {
	       if (values->delta[i][j]!=values->delta[j][i]) {
			 fprintf(stderr,"Character transformation matrix is asymmetrical\n");
			 exit(-1);
			 }
	       }
	  }
     }
     /*metricity*/
if (values->delta) {
 if (values->delta[0][1]>(values->delta[0][2]+values->delta[1][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][1]>(values->delta[0][3]+values->delta[1][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][2]>(values->delta[0][1]+values->delta[1][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][2]>(values->delta[0][3]+values->delta[2][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][3]>(values->delta[0][1]+values->delta[1][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[0][3]>(values->delta[0][2]+values->delta[2][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][2]>(values->delta[1][0]+values->delta[0][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][2]>(values->delta[1][3]+values->delta[3][2])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][3]>(values->delta[1][0]+values->delta[0][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[1][3]>(values->delta[1][2]+values->delta[2][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[2][3]>(values->delta[2][0]+values->delta[0][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
     else if (values->delta[2][3]>(values->delta[2][1]+values->delta[1][3])) {
	  fprintf(stderr,"Transformation matrix is non-metric.\n");
     exit(-1);
	  }
}
/*Check for transition transversion situation*/
if (values->delta) {
     if ((values->delta[0][1]==values->delta[0][3]) && (values->delta[0][1]==values->delta[1][2]) && (values->delta[0][1]==values->delta[2][3])) {
	       if (values->delta[0][2]==values->delta[1][3]) {
		    values->ttr=1;
			 values->transition=values->delta[0][2];
	       values->transversion=values->delta[0][1]-values->transition;
			 /*for (i=0;i<4;i++) free(values->delta[i]);
			 free(values->delta);
			 values->delta=NULL;*/
	       fprintf(stderr,"Matrix specifies only transitions and transversions.\n");
			 fprintf(stderr,"    Transversion cost %d Transition cost %d\n\n",values->transversion+values->transition,values->transition);
	       }
	  }
     }

if ((!values->VERBOSE || !values->rep_error) && values->show_mem) values->show_mem=0;
if (values->iter) values->max_gap=100;
if (values->cost_only) if (!values->VERBOSE) fprintf(stderr,"Determining cost of preexisting alignment.\n");
if (values->cost_only) values->get_heur=values->get_heur2=values->best=values->rand_align=0;
if ((values->tree_rand_order_max>0) || (values->rand_order>0) || (values->pref_direc==RANDOM) || (values->how_many>1) || (values->randtrees>1)) {
     srand((int) time(NULL));
  }
if (values->rand_order > 0) {
     if (values->rep_error) {
	  if (values->get_heur) if (values->rep_error) fprintf(stderr,"\nRandom orders will not affect pair.\n");
	  if (values->best) if (values->rep_error) fprintf(stderr,"\nRandom orders will not affect an exact solution.\n");
    }
     }

/*until Sankoff fixed
if (values->delta) {
     fprintf(stderr,"This option is currently disabled.  Please specify only transitions and transversion costs in matrix.\n");
  exit(-1);
     }*/

/*print file names
for (i=0;i<values->number_of_input_alignments;i++) {
     fprintf(stderr,"input alignment %d from file %s\n",i,values->input_file_name[i]);
     }*/
/*disable cache for multiple input alignments*/
if (values->number_of_input_alignments>0) if ((values->elision) || (values->cull))  {
     values->cache_size=0;
     values->cache_info=0;
     fprintf(stderr,"Cache is disabled when input data are alignments.\n");
     }
/*checks for cull and elision*/
if (values->cull || values->elision) {
     values->output_order=0;
     if (values->number_of_input_alignments<2) {
	  fprintf(stderr,"Cannot use options 'cull' or 'elision' unless multiple alignments are input\n");
	  exit(-1);
	  }
     }
if (values->cull && values->elision) {
     fprintf(stderr,"Cannot use both 'cull' and 'ellision' in the same run--'cull' is disabled\n    Rerun with 'cull' only to get both results\n");
     values->cull=0;
     }
if (values->elision) {
     values->get_heur=0;
     values->best=0;
     values->get_heur2=0;
     values->rand_align=0;
     if (values->hen_gap) {
	  values->hen_gap=0;
	  values->farris=1;
	  }
     }
/*
if ((values->check_cache_for_scores) && (!values->cache_size)) {
     fprintf(stderr,"Cannot cache cladogram scores without alignment cache.\n");
     values->check_cache_for_scores=0;
     }
     */
if (values->grain_size==(-1)) values->grain_size=((2*values->all_done)-1);
if (values->new_optimization) values->phylo_score=0;
if ((values->new_optimization && values->number_of_input_alignments) || (values->chop)) {
     if (values->number_of_input_alignments) fprintf(stderr,"Input data sets : %d\n",values->number_of_input_alignments);
     else fprintf(stderr,"Input data sets : %d\n",values->number_of_input_alignments+1);
     for (i=0;i<values->number_of_input_alignments;i++) {
	  fprintf(stderr,"    Data set %d with weight %d will ",i+1,values->data_set_weights[i]);
	  if ((values->data_sets[i]==0) || (values->data_sets[i]==3)) fprintf(stderr,"be aligned.\n");
	  else {
	       fprintf(stderr,"not be aligned");
	       if (values->data_sets[i]==1) fprintf(stderr," since it is already aligned.\n");
	       if (values->data_sets[i]==2) fprintf(stderr," since it consists of character data.\n");
	       }
	  }

     }
if (values->get_heur && values->new_optimization) {
	fprintf(stderr,"Cannot perform 'pair' with optimization alignments.\n");
	values->get_heur=0;
	}
if (values->apolist) {
	values->optimize_nodes=1;
	if (!values->new_optimization) {
		fprintf(stderr,"Cannot optimize alignments.\n");
		values->optimize_nodes=values->apolist=0;
		}
	}
get_lookup2(values);
}

void get_lookup2(values)
parameters *values;
{
unsigned char i,j,t,h,ih,jh;
int MiNiMuM_cost;

for (i=0;i<32;i++) for (j=0;j<32;j++) values->lookup[i][j].cost=values->lookup[i][j].base=values->lookup[i][j].union_bit=0;
if (!values->ttr && (!values->delta)) {
    if ((values->gap_cost > values->change_cost)) {
    /*gaps > V=I */
    for (i=1;i<32;i++) {
            for (j=1;j<32;j++) {
                if (i & j) {
                    values->lookup[i][j].base=(i & j);
                    values->lookup[i][j].cost=0;
                    values->lookup[i][j].union_bit=0;
                    }
                else {
                   values->lookup[i][j].union_bit=1;
                   /*no intersection problems
                    so straight union but if a gap and base on one side
                    then only bases cause lower cost (A-) + (G) => R*/
                   t=(i | j);
                   if (i>GAP_BIT) t-=GAP_BIT;
                   else if (j>GAP_BIT) t-=GAP_BIT;
                   values->lookup[i][j].base=t;
                   /*cost is a base change when result is lower than gap*/
                   if (t<GAP_BIT) values->lookup[i][j].cost=values->change_cost;
                   else values->lookup[i][j].cost=values->gap_cost;
                }
            }
        }
    }
    else if ((values->gap_cost == values->change_cost)) {
    /*gaps > V=I */
    for (i=1;i<32;i++) {
            for (j=1;j<32;j++) {
                if (i & j) {
                    values->lookup[i][j].base=(i & j);
                    values->lookup[i][j].cost=0;
                    values->lookup[i][j].union_bit=0;
                    }
                else {
                   values->lookup[i][j].union_bit=1;
                   t=(i | j);
                   values->lookup[i][j].base=t;
                   values->lookup[i][j].cost=values->gap_cost;
                }
            }
        }
    }
    else {fprintf(stderr,"Not yet implemented\n");exit(-1);}
}
else  {
    for (i=1;i<32;i++) {
        for (j=1;j<32;j++) {
            if (i & j) {
               values->lookup[i][j].base=(i & j);
               values->lookup[i][j].cost=0;
               values->lookup[i][j].union_bit=0;
               }
          else {
                MiNiMuM_cost=HUGE_COST;
                values->lookup[i][j].union_bit=1;
                /*find cheapest cost possible*/
                if (i & A_BIT) {
                       if ((j & C_BIT) && (values->delta[0][1]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[0][1];
                       if ((j & G_BIT) && (values->delta[0][2]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[0][2];
                       if ((j & T_BIT) && (values->delta[0][3]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[0][3];
                       if ((j & GAP_BIT) && (values->gap_cost<MiNiMuM_cost)) MiNiMuM_cost=values->gap_cost;
                }
                if (i & C_BIT) {
                       if ((j & A_BIT) && (values->delta[1][0]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[1][0];
                       if ((j & G_BIT) && (values->delta[1][2]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[1][2];
                       if ((j & T_BIT) && (values->delta[1][3]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[1][3];
                       if ((j & GAP_BIT) && (values->gap_cost<MiNiMuM_cost)) MiNiMuM_cost=values->gap_cost;
                }
                if (i & G_BIT) {
                       if ((j & A_BIT) && (values->delta[2][0]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[2][0];
                       if ((j & C_BIT) && (values->delta[2][1]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[2][1];
                       if ((j & T_BIT) && (values->delta[2][3]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[2][3];
                       if ((j & GAP_BIT) && (values->gap_cost<MiNiMuM_cost)) MiNiMuM_cost=values->gap_cost;
                }
                if (i & T_BIT) {
                       if ((j & A_BIT) && (values->delta[3][0]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[3][0];
                       if ((j & C_BIT) && (values->delta[3][1]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[3][1];
                       if ((j & G_BIT) && (values->delta[3][2]<MiNiMuM_cost)) MiNiMuM_cost=values->delta[3][2];
                       if ((j & GAP_BIT) && (values->gap_cost<MiNiMuM_cost)) MiNiMuM_cost=values->gap_cost;
                }
                if ((i & GAP_BIT) && (values->gap_cost<MiNiMuM_cost)) MiNiMuM_cost=values->gap_cost;
                values->lookup[i][j].cost=MiNiMuM_cost;
                /*get ancestral*/
                if (i & A_BIT) {
                       if ((j & C_BIT) && (values->delta[0][1]==MiNiMuM_cost)) values->lookup[i][j].base |= (A_BIT | C_BIT);
                       if ((j & G_BIT) && (values->delta[0][2]==MiNiMuM_cost)) values->lookup[i][j].base |= (A_BIT | G_BIT);
                       if ((j & T_BIT) && (values->delta[0][3]==MiNiMuM_cost)) values->lookup[i][j].base |= (A_BIT | T_BIT);
                       if ((j & GAP_BIT) && (values->gap_cost==MiNiMuM_cost)) values->lookup[i][j].base |= (A_BIT | GAP_BIT);
                }
                if (i & C_BIT) {
                       if ((j & A_BIT) && (values->delta[1][0]==MiNiMuM_cost)) values->lookup[i][j].base |= (C_BIT | A_BIT);
                       if ((j & G_BIT) && (values->delta[1][2]==MiNiMuM_cost)) values->lookup[i][j].base |= (C_BIT | G_BIT);
                       if ((j & T_BIT) && (values->delta[1][3]==MiNiMuM_cost)) values->lookup[i][j].base |= (C_BIT | T_BIT);
                       if ((j & GAP_BIT) && (values->gap_cost==MiNiMuM_cost)) values->lookup[i][j].base |= (C_BIT | GAP_BIT);
                }
                if (i & G_BIT) {
                       if ((j & A_BIT) && (values->delta[2][0]==MiNiMuM_cost)) values->lookup[i][j].base |= (G_BIT | A_BIT);
                       if ((j & C_BIT) && (values->delta[2][1]==MiNiMuM_cost)) values->lookup[i][j].base |= (G_BIT | C_BIT);
                       if ((j & T_BIT) && (values->delta[2][3]==MiNiMuM_cost)) values->lookup[i][j].base |= (G_BIT | T_BIT);
                       if ((j & GAP_BIT) && (values->gap_cost==MiNiMuM_cost)) values->lookup[i][j].base |= (G_BIT | GAP_BIT);
                }
                if (i & T_BIT) {
                       if ((j & A_BIT) && (values->delta[3][0]==MiNiMuM_cost)) values->lookup[i][j].base |= (T_BIT | A_BIT);
                       if ((j & C_BIT) && (values->delta[3][1]==MiNiMuM_cost)) values->lookup[i][j].base |= (T_BIT | C_BIT);
                       if ((j & G_BIT) && (values->delta[3][2]==MiNiMuM_cost)) values->lookup[i][j].base |= (T_BIT | G_BIT);
                       if ((j & GAP_BIT) && (values->gap_cost==MiNiMuM_cost)) values->lookup[i][j].base |= (T_BIT | GAP_BIT);
                }
                if ((i & GAP_BIT) && (values->gap_cost==MiNiMuM_cost)) {
                       if (j & A_BIT)  values->lookup[i][j].base |= (GAP_BIT | A_BIT);
                       if (j & C_BIT)  values->lookup[i][j].base |= (GAP_BIT | C_BIT);
                       if (j & G_BIT)  values->lookup[i][j].base |= (GAP_BIT | G_BIT);
                       if (j & T_BIT) values->lookup[i][j].base |= (GAP_BIT | T_BIT);
                }
                assert((int) values->lookup[i][j].base);

          }
        }/*j*/
    }/*i*/
   }/*else*/
/*for (i=1;i<32;i++) {
    for (j=1;j<32;j++) printf("(%d,%d=>%d,%d,%d) ",i,j,values->lookup[i][j].base,values->lookup[i][j].cost,values->lookup[i][j].union_bit);
    printf("\n");
}*/
}

void get_lookup(values)
parameters *values;
{
unsigned char i,j,t,h,ih,jh;

for (i=0;i<32;i++) for (j=0;j<32;j++) values->lookup[i][j].cost=values->lookup[i][j].base=values->lookup[i][j].union_bit=0;
if (!values->ttr) {
    if ((values->gap_cost > values->change_cost)) {
    /*gaps > V=I */
    for (i=1;i<32;i++) {
            for (j=1;j<32;j++) {
                if (i & j) {
                    values->lookup[i][j].base=(i & j);
                    values->lookup[i][j].cost=0;
                    values->lookup[i][j].union_bit=0;
                    }
                else {
                   values->lookup[i][j].union_bit=1;
                   /*no intersection problems
                    so straight union but if a gap and base on one side
                    then only bases cause lower cost (A-) + (G) => R*/
                   t=(i | j);
                   if (i>GAP_BIT) t-=GAP_BIT;
                   else if (j>GAP_BIT) t-=GAP_BIT;
                   values->lookup[i][j].base=t;
                   /*cost is a base change when result is lower than gap*/
                   if (t<GAP_BIT) values->lookup[i][j].cost=values->change_cost;
                   else values->lookup[i][j].cost=values->gap_cost;
                }
            }
        }
    }
    else if ((values->gap_cost == values->change_cost)) {
    /*gaps > V=I */
    for (i=1;i<32;i++) {
            for (j=1;j<32;j++) {
                if (i & j) {
                    values->lookup[i][j].base=(i & j);
                    values->lookup[i][j].cost=0;
                    values->lookup[i][j].union_bit=0;
                    }
                else {
                   values->lookup[i][j].union_bit=1;
                   t=(i | j);
                   values->lookup[i][j].base=t;
                   values->lookup[i][j].cost=values->gap_cost;
                }
            }
        }
    }
    else {fprintf(stderr,"Not yet implemented--but really, indels > changes\n");exit(-1);}
}
else if (values->gap_cost > (values->transversion + values->transition)) {
    /*Gaps > V > I */
    fprintf(stderr,"%d %d %d,",values->gap_cost, values->transition, (values->transversion + values->transition));
    for (i=1;i<32;i++) {
    for (j=1;j<32;j++) {
         if (i & j) {
               values->lookup[i][j].base=(i & j);
               values->lookup[i][j].cost=0;
               values->lookup[i][j].union_bit=0;
               }
          else {
              values->lookup[i][j].union_bit=1;
         /*no intersection problems
           so straight union but if a gap and base on one side
           then only bases cause lower cost (A-) + (G) => R*/
          t=(i | j);
          if (i>GAP_BIT) t-=GAP_BIT;
          else if (j>GAP_BIT) t-=GAP_BIT;
          if (t>GAP_BIT) {/*indel*/
            values->lookup[i][j].base=t;
            values->lookup[i][j].cost=values->gap_cost;
            }
          else { /*must not be an indel--if indel done*/
              values->lookup[i][j].base=t;
              ih=i;
              jh=j;
              if (i>GAP_BIT) ih-=GAP_BIT;
              if (j>GAP_BIT) jh-=GAP_BIT;
              /*cost is a base change when result is lower than gap*/
              h=num_states(t);
              if (h==2) {
                if (t == R_BIT) values->lookup[i][j].cost=values->transition;
                else if (t == Y_BIT) values->lookup[i][j].cost=values->transition;
                else values->lookup[i][j].cost=(values->transversion + values->transition);
                }
              else if (h==3) {
                if (ih==R_BIT) values->lookup[i][j].cost=(values->transversion + values->transition);
                else if (ih==Y_BIT) values->lookup[i][j].cost=(values->transversion + values->transition);
                else if (jh==R_BIT) values->lookup[i][j].cost=(values->transversion + values->transition);
                else if (jh==Y_BIT) values->lookup[i][j].cost=(values->transversion + values->transition);
                else { /*must be a transition ans have the correct encestor*/
                    if (ih==A_BIT) values->lookup[i][j].base=R_BIT;
                    else if (ih==G_BIT) values->lookup[i][j].base=R_BIT;
                    else if (ih==C_BIT) values->lookup[i][j].base=Y_BIT;
                    else if (ih==T_BIT) values->lookup[i][j].base=Y_BIT;
                    else if (jh==A_BIT) values->lookup[i][j].base=R_BIT;
                    else if (jh==G_BIT) values->lookup[i][j].base=R_BIT;
                    else if (jh==C_BIT) values->lookup[i][j].base=Y_BIT;
                    else if (jh==T_BIT) values->lookup[i][j].base=Y_BIT;
                    values->lookup[i][j].cost=values->transition;
                    }
                }
              else {
                if ((ih==R_BIT) && (jh==Y_BIT)) values->lookup[i][j].cost=(values->transversion + values->transition);
                else {/*must be a transition ans have the correct encestor*/
                    if (ih==A_BIT) values->lookup[i][j].base=R_BIT;
                    else if (ih==G_BIT) values->lookup[i][j].base=R_BIT;
                    else if (ih==C_BIT) values->lookup[i][j].base=Y_BIT;
                    else if (ih==T_BIT) values->lookup[i][j].base=Y_BIT;
                    else if (jh==A_BIT) values->lookup[i][j].base=R_BIT;
                    else if (jh==G_BIT) values->lookup[i][j].base=R_BIT;
                    else if (jh==C_BIT) values->lookup[i][j].base=Y_BIT;
                    else if (jh==T_BIT) values->lookup[i][j].base=Y_BIT;
                    values->lookup[i][j].cost=values->transition;
                    }
                }
              }
            }
        }
    }
    }
   else {
    fprintf(stderr,"Lookup not yet implemented\n");
   }
/*for (i=1;i<32;i++) {
    for (j=1;j<32;j++) printf("(%d,%d=>%d,%d,%d) ",i,j,values->lookup[i][j].base,values->lookup[i][j].cost,values->lookup[i][j].union_bit);
    printf("\n");
}*/
fprintf(stderr,"On to 2\n");
get_lookup2(values);
}
