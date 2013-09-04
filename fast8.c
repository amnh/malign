#include "align3.h"
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

parameters *values_here;
extern int *cc_weights;

int NOTUs, NCharacters, NHTUs, NTaxa;
look_up_thang **Character;
unsigned char **Name;
unsigned char buf[1024];
NodeT *Node, **PlacedTaxa;
int NPlacedTaxa;
int NTrees = 0, NOptimizations = 0;
double Optimizations = 0;
int NDirtyNodes = 0;
NodeT **DirtyNode;
BitVectorT *SupportedClades;
int NSupportedClades;
unsigned int SupportedCladesHash;
int NConsensusClades;
BitVectorT *ConsensusClades;
BufferedTreeT *ScratchBufferedTree = NOBUFFEREDTREE;
int tasks_available = 0, max_tasks_available, n_workers = 0, next_task = 0;
int *PartialFinalMap;

struct hostinfo *hostp;
int nhost = 0, narch, tids[256];
struct {
  int tid;
  int treeid;
  int ntasks;
  char *name;
} task[256];
char **serverargv;
int treeserial = 0, buffertick = 0;

int EPSILON = 0;
int BACKUP = 1, VERBOSE = 0, PRUNE = 1, doSPR = 1, doTBR = 1,
LOCALREARRANGEMENT = 1,
PABLO = 0, INITIALPABLO = 0, STEEPEST = 0,
TREE = 0, INDENT = 0, STATS = 0,
MAXWORKERS = 0, MINWORKERS = 0, LOOPCHUNK = 10000,
DOUBLEBUFFER = 1,
MAXTREES, FINAL_MAXTREES = 1024, INITIAL_MAXTREES = 0, maxtrees_exceeded = 0,
CONSENSUS = 0, CONSENSUSBITS = 0, CONSENSUSTABLE = 0, NUMBERCONSENSUSNODES = 0,
RANDOMIZATIONS = 0, FILLMISSING = 0,
FILLINS = 0,
MISSINGREPORT = 0, MISSINGSTATS = 0, MISSINGCONSENSUS = 0, MISSINGCONSENSUSDATA = 0,
TAGS = 0;
char *POORHOST = NULL;
enum {CHARACTER, BITMAP} FORMAT = CHARACTER;
ComputeSupportT SupportComputer;
enum {STANDALONE, MASTER, SERVER} SERVERMODE = STANDALONE;

FILE *blather;

void NewReadData (sequences)
look_up_thang **sequences;
{
  int taxon, character, i, j;
/*  printf ("%d OTUs.\n", NOTUs);
  printf ("%d characters.\n", NCharacters);*/
  NHTUs = NOTUs - 1;
  NTaxa = NOTUs + NHTUs;
  NewMatrix (Character, NTaxa * 2, NCharacters, look_up_thang);
/*  NewMatrix (BitPosition, 256, NCharacters, int);
  NewArray (UsedPositions, NCharacters, int);*/
  NewArray (Name, NTaxa, unsigned char *);
  /*NewArray (MissingData, NCharacters, unsigned int);
  LoopBelow (character, NCharacters) {
    UsedPositions[character] = 0;
    MissingData[character] = 0;
    LoopBelow (i, 256)
      BitPosition[i][character] = -1;
  }*/
  /*LoopBelow (taxon, NOTUs) memcpy(Character[taxon],sequences[taxon],NCharacters*sizeof(look_up_thang));*/
  for (j=0;j<NOTUs;j++) {
  	for (i=0;i<NCharacters;i++) {
  		Character[j][i].base=sequences[j][i].base;
	   	/*if (!i) fprintf(stderr,"%2d vs. %2d ",sequences[j][i].base,Character[j][i].base);*/
	  	}
	  }
buf[0]='a';
buf[1]='\0';
LoopBelow (taxon, NOTUs) {
 	Name[taxon]=CopyString(buf);
 	}
}

int RandomCharacter (int universe) {
  int i, j, character;
  character = rand () % CountBits (universe);
  j = 0;
  for (i = 0;; i++)
    if ((1 << i) & universe) {
      if (j is character)
	return (1 << i);
      else
	j++;
    }
}

void WriteData (FILE *f) {
  int taxon, character;
  fprintf (f, "%d %d\n", NOTUs, NCharacters);
  LoopBelow (taxon, NOTUs) {
    fprintf (f, "%s", Node[taxon].name);
    LoopBelow (character, NCharacters) {
      fprintf (f, " %d", RandomCharacter (Character[taxon][character].base));
    }
    fprintf (f, "\n");
  }
}

unsigned char **RawData;


/* tree buffer stuff */


BufferedTreeT *NewBufferedTree () {
  int i;
  BufferedTreeT *bt;
  bt = New (1, BufferedTreeT);
  NewArray (bt->nodes, NTaxa, BufferedTreeNodeT);
  NewArray (bt->placed_taxa_indices, NTaxa, int);
  NewArray (bt->supported_clades, NHTUs - 1, BitVectorT);
  LoopBelow (i, NHTUs - 1)
    bt->supported_clades[i] = NewBitVector ();
  return (bt);
}

NodeT *FindRoot (NodeT *n)
{
  assert (n);
  while (n->parent) n = n->parent;
  return (n);
}

void CopyBufferedTree (BufferedTreeT *dest, BufferedTreeT *source) {
  int i;
  dest->n_placed_taxa = source->n_placed_taxa;
  LoopBelow (i, source->n_placed_taxa) {
    dest->nodes[i] = source->nodes[i];
    dest->placed_taxa_indices[i] = source->placed_taxa_indices[i];
  }
  dest->n_supported_clades = source->n_supported_clades;
  LoopBelow (i, source->n_supported_clades)
    CopyBitVector (dest->supported_clades[i], source->supported_clades[i]);
  dest->supported_clades_hash = source->supported_clades_hash;
  dest->cost = source->cost;
  dest->generation = source->generation;
}

void StoreTree (BufferedTreeT *bt, int taxa_in_use, int cost) {
  int i;
  NodeT *n, *p;
/*   fprintf (stderr, "store tree with %d taxa\n", taxa_in_use); */
  LoopBelow (i, taxa_in_use) {
    n = PlacedTaxa[i];
    if (n->left && !n->right) fail ("left but no right");
    if (n->right && !n->left) fail ("right but no left");
  }
  LoopBelow (i, taxa_in_use) {
    n = PlacedTaxa[i];
    p = n->parent;
    if (p is NONODE)
      bt->nodes[i].parent_index = -1;
    else {
      bt->nodes[i].parent_index = p->taxon;
      bt->nodes[i].on_left = (n is p->left);
    }
  }
  bt->n_placed_taxa = taxa_in_use;
  LoopBelow (i, taxa_in_use)
    bt->placed_taxa_indices[i] = PlacedTaxa[i]->taxon;
  bt->n_supported_clades = NSupportedClades;
  LoopBelow (i, NSupportedClades) {
    CopyBitVector (bt->supported_clades[i], SupportedClades[i]);
  }
  bt->supported_clades_hash = SupportedCladesHash;
  bt->cost = cost;
}

void RecomputeNode (NodeT *n) {
  int i,ii,jj;
  parameters *valuesI;

  if (n->left) {
    RecomputeNode (n->left);
    RecomputeNode (n->right);
    i=0;
    for (ii=0;ii<max(1,values_here->number_of_input_alignments);ii++) {
        if (!values_here->other_parm) valuesI=values_here;
        else if (values_here->other_parm[ii]) valuesI=values_here->other_parm[ii];
        else valuesI=values_here;
        for (jj=values_here->start[ii];jj<values_here->stop[ii];jj++) {
          n->character[i]=valuesI->lookup[n->left->character[i].base][n->right->character[i].base];
          n->localcost+=(n->character[i].cost*cc_weights[i]);
      	  i++;
    	} /*jj*/
      } /*end of characters loop*/
    n->totalcost = n->localcost + n->left->totalcost + n->right->totalcost;
  }
}

void LoadTree (BufferedTreeT *bt, int restore_supported_clades) {
  extern void ComputeSupportedClades ();
  int i;
  int taxa_in_use = bt->n_placed_taxa;
  NodeT *n;
  NPlacedTaxa = bt->n_placed_taxa;

/*   fprintf (stderr, "loading tree with %d taxa\n", NPlacedTaxa); */
  LoopBelow (i, NPlacedTaxa) {
    n = PlacedTaxa[i] = Node + bt->placed_taxa_indices[i];
    n->parent = n->left = n->right = NONODE;
  }
  LoopBelow (i, taxa_in_use) {
    n = PlacedTaxa[i];
    if (bt->nodes[i].parent_index < 0) {
      if (n->parent isnt NONODE) fail ("parent already set");
      n->parent = NONODE;
    } else {
      if (n->parent isnt NONODE) fail ("parent already set");
      n->parent = Node + bt->nodes[i].parent_index;
      if (bt->nodes[i].on_left) {
	if (n->parent->left isnt NONODE) fail ("left already set");
	n->parent->left = n;
      } else {
	if (n->parent->right isnt NONODE) fail ("right already set");
	n->parent->right = n;
      }
    }
    n->localcost = n->totalcost = n->dirty = 0;
  }
  LoopBelow (i, NPlacedTaxa) {
    n = PlacedTaxa[i];
    if (n->left && !n->right) fail ("left but no right");
    if (n->right && !n->left) fail ("right but no left");
  }
  RecomputeNode (FindRoot (PlacedTaxa[0]));
  if (restore_supported_clades) {
    ComputeSupportedClades ();
  }
}

/* int Combine (int x, int y) { */
/*   int z; */
/*   z = SETBITS & x & y; */
/*   if (! z) */
/*     z = SETBITS & (x | y); */
/*   return (z); */
/* } */

void ComputeFinal (NodeT *n) {
  int i;
  look_up_thang n_ch;
  NodeT *parent = n->parent;
  LoopBelow (i, NCharacters) {
    /*n_ch = n->character[i] & SETBITS;*/
    n_ch=n->character[i];
    n_ch.union_bit=0;
    /*if root*/
    if (parent is NONODE) n->final_character[i].base = n_ch.base;
    /*if parent is contained in node -- should be OK*/
    else if (((~n_ch.base) & parent->final_character[i].base) is 0) n->final_character[i].base = parent->final_character[i].base;
    /*if union on down add parent--needs to be fixed wrt matrix*/
    else if (n->character[i].union_bit) n->final_character[i].base = n_ch.base | parent->final_character[i].base;
    /*if was inrtersection add anything in parent and one of the descendents--need to be fixed wrt matrix*/
    else if (n->left) n->final_character[i].base = n_ch.base | (parent->final_character[i].base & (n->left->character[i].base | n->right->character[i].base));
    /*if terminal*/
    else n->final_character[i].base = n_ch.base;
    n->final_character[i].union_bit=0; /*cost should not matter*/
/*     if (parent && n->left) { */
/*       NodeT *left = n->left, *right = n->right; */
/*       int z = */
/* 	Combine (parent->final_character[i], Combine (left->character[i], right->character[i])) & */
/* 	  Combine (left->character[i], Combine (parent->final_character[i], right->character[i])) & */
/* 	    Combine (right->character[i], Combine (left->character[i], parent->final_character[i])); */
/*       if (z isnt n->final_character[i]) */
/* 	fprintf (stderr, "Discrepancy at taxon %d character %d: sm = %d, fitch = %d\n", */
/* 		 n->taxon, i, z, n->final_character[i]); */
/*     } */
  }
  if (n->left) {
    ComputeFinal (n->left);
    ComputeFinal (n->right);
  }
}


void ComputeSupportEverything (NodeT *n) {
  if (n->left is NONODE) {
    ClearBitVector (n->taxa);
    SetBit (n->taxa, n->taxon);
  } else {
    ComputeSupportEverything (n->left);
    ComputeSupportEverything (n->right);
    OrBitVectors (n->taxa, n->left->taxa, n->right->taxa);
    if (n->parent isnt NONODE)
      SupportedClades[NSupportedClades++] = n->taxa;
  }
}

void FindSupportMust (NodeT *n) {
  int supported = 0, i;
  NodeT *p = n->parent;
  if (p) {
    LoopBelow (i, NCharacters)
      if (!(n->final_character[i].base & p->final_character[i].base)) {
/* 	fprintf (stderr, "taxon %d character %d supported\n", n->taxon, i); */
	supported = 1;
	break;
      }
  }
  if (n->left is NONODE) {
    ClearBitVector (n->taxa);
    SetBit (n->taxa, n->taxon);
  } else {
    FindSupportMust (n->left);
    FindSupportMust (n->right);
    OrBitVectors (n->taxa, n->left->taxa, n->right->taxa);
    if (supported) {
      SupportedClades[NSupportedClades++] = n->taxa;
    }
  }
}

void ComputeSupportMust (NodeT *root) {
  ComputeFinal (root);
  FindSupportMust (root);
}

int IsSingleton (look_up_thang x) {
  int i, y;
  y = 1;
  LoopBelow (i, 26) /*wordsize in yapp but 26 bits in structure*/
    if (x.base is y)
      return (1);
    else
      y <<= 1;
  return (0);
}

void FindSupportMay (NodeT *n) {
  int supported = 0, i;
  NodeT *p = n->parent;
  if (p) {
    LoopBelow (i, NCharacters)
      if (n->final_character[i].base isnt p->final_character[i].base) {
	supported = 1;
	break;
      }
  }
  if (n->left is NONODE) {
    ClearBitVector (n->taxa);
    SetBit (n->taxa, n->taxon);
  } else {
    FindSupportMay (n->left);
    FindSupportMay (n->right);
    OrBitVectors (n->taxa, n->left->taxa, n->right->taxa);
    if (supported) {
      SupportedClades[NSupportedClades++] = n->taxa;
    }
  }
}

void ComputeSupportMay (NodeT *root) {
  ComputeFinal (root);
  FindSupportMay (root);
}

void FindSupportNotSameSingleton (NodeT *n) {
  int supported = 0, i;
  NodeT *p = n->parent;
  if (p) {
    LoopBelow (i, NCharacters)
      if (! (n->final_character[i].base is p->final_character[i].base &&
	     IsSingleton(n->final_character[i]))) {
	supported = 1;
	break;
      }
  }
  if (n->left is NONODE) {
    ClearBitVector (n->taxa);
    SetBit (n->taxa, n->taxon);
  } else {
    FindSupportNotSameSingleton (n->left);
    FindSupportNotSameSingleton (n->right);
    OrBitVectors (n->taxa, n->left->taxa, n->right->taxa);
    if (supported) {
      SupportedClades[NSupportedClades++] = n->taxa;
    }
  }
}

void ComputeSupportNotSameSingleton (NodeT *root) {
  ComputeFinal (root);
  FindSupportNotSameSingleton (root);
}

void ComputeSupportedClades () {
  int i;
  NodeT *root = FindRoot (PlacedTaxa[0]);

  NSupportedClades = 0;
  (*SupportComputer) (root);
/*   fprintf (stderr, "%d supported clades\n", NSupportedClades); */
  qsort (SupportedClades, NSupportedClades, sizeof (BitVectorT), QsortCompareBitVectors);
  SupportedCladesHash = 0;
  LoopBelow (i, NSupportedClades)
    SupportedCladesHash += BitVectorHash (SupportedClades[i]);
}

int CladeIsSupported (BitVectorT taxa) {
  return ((int)bsearch (&taxa, SupportedClades, NSupportedClades,
			sizeof (BitVectorT), QsortCompareBitVectors));
}

int IsSupported (NodeT *n) {
  return (CladeIsSupported (n->taxa));
}

void PrintConsensusTable (FILE *f) {
  int i, j;
  LoopBelow (i, NConsensusClades) {
    fprintf (f, "%d ", i+1);
    LoopBelow (j, NOTUs)
      if (GetBit (ConsensusClades[i], j))
	fprintf (f, "%s ", Node[j].name);
    fprintf (f, "\n");
  }
}

void PrintConsensusBits (FILE *f) {
  int i;
  char *tag = TAGS? "ConsensusBits ": "";
  fprintf (f, "%s%d\n", tag, NOTUs);
  LoopBelow (i, NConsensusClades) {
    fputs (tag, f);
    fprintfBitVector (f, ConsensusClades[i]);
    fprintf (f, "\n");
  }
}

void InitConsensus () {
  int i, n_clades = NSupportedClades + NOTUs + 1;
  NewArray (ConsensusClades, n_clades, BitVectorT);
  LoopBelow (i, n_clades)
    ConsensusClades[i] = NewBitVector ();
  LoopBelow (i, NSupportedClades)
    CopyBitVector (ConsensusClades[i], SupportedClades[i]);
  NConsensusClades = NSupportedClades;
}

void Consense () {
  int i, j;
  j = 0;
  LoopBelow (i, NConsensusClades) {
    if (CladeIsSupported (ConsensusClades[i])) {
      BitVectorT temp;
      temp = ConsensusClades[i];
      ConsensusClades[i] = ConsensusClades[j];
      ConsensusClades[j] = temp;
      j++;
    }
  }
  NConsensusClades = j;
}

BitVectorT *oldconsensusclades;
int noldconsensusclades;

void PrintConsensusNode (FILE *f, ConsensusNodeT *n) {
  int i, j;
  if (n->n_kids is 0) {
    LoopBelow (i, NTaxa)
      if (GetBit (ConsensusClades[n->clade], i)) {
	fprintf (f, "%s", Node[i].name);
	break;
      }
  } else {
    LoopBelow (i, n->n_kids) {
      if (i is 0) {
	LoopBelow (j, noldconsensusclades)
	  if (! CompareBitVectors (ConsensusClades[n->clade], oldconsensusclades[j]))
	    break;
	if (j < noldconsensusclades)
	  fprintf (f, "%d", j+1);
	fprintf (f, "(");
      } else
	fprintf (f, " ");
      PrintConsensusNode (f, n->kids[i]);
    }
    fprintf (f, ")");
  }
}

void PrintConsensus (FILE *f) {
  ConsensusNodeT *nodes;
  int i, j;
  if (NUMBERCONSENSUSNODES) {
    NewArray (oldconsensusclades, NConsensusClades, BitVectorT);
    LoopBelow (i, NConsensusClades) {
      oldconsensusclades[i] = NewBitVector ();
      CopyBitVector (oldconsensusclades[i], ConsensusClades[i]);
    }
    noldconsensusclades = NConsensusClades;
  } else
    noldconsensusclades = 0;
  ClearBitVector (ConsensusClades[NConsensusClades + NOTUs]);
  LoopBelow (i, NOTUs) {
    ClearBitVector (ConsensusClades[NConsensusClades + i]);
    SetBit (ConsensusClades[NConsensusClades + i], i);
    SetBit (ConsensusClades[NConsensusClades + NOTUs], i);
  }
  NConsensusClades += NOTUs + 1;
  qsort (ConsensusClades, NConsensusClades, sizeof (BitVectorT), QsortCompareBitVectors);
  NewArray (nodes, NConsensusClades, ConsensusNodeT);
  LoopBelow (i, NConsensusClades) {
    nodes[i].n_kids = 0;
    NewArray (nodes[i].kids, NConsensusClades, ConsensusNodeT *);
    nodes[i].has_parent = 0;
    nodes[i].clade = i;
  }
  LoopBelow (i, NConsensusClades) {
    LoopInterval (j, i + 1, NConsensusClades - 1) {
      if (BitVectorSubsetEq (ConsensusClades[i], ConsensusClades[j])) {
	nodes[j].kids[nodes[j].n_kids++] = &nodes[i];
	nodes[i].has_parent = 1;
	break;
      }
    }
  }
  LoopBelow (i, NConsensusClades) {
    if (! nodes[i].has_parent)
      PrintConsensusNode (f, &nodes[i]);
  }
}

void CompressTreeBuffer (TreeBufferT *tb) {
  int max_acceptable_cost = tb->min_tree_cost + Max (0, EPSILON);
  int i, j;
  j = 0;
  LoopBelow (i, tb->n_buffered_trees) {
    if (tb->buffered_tree[i]->cost <= max_acceptable_cost) {
      BufferedTreeT *temp = tb->buffered_tree[i];
      tb->buffered_tree[i] = tb->buffered_tree[j];
      tb->buffered_tree[j] = temp;
      j++;
    }
  }
  tb->n_buffered_trees = j;
  tb->compressed = 1;
}

#define INITIALMAXBUFFEREDTREES 1024
TreeBufferT *NewTreeBuffer () {
  int i;
  TreeBufferT *tb = New (1, TreeBufferT);
  tb->n_buffered_trees = 0;
  tb->max_buffered_trees = INITIALMAXBUFFEREDTREES;
  NewArray (tb->buffered_tree, tb->max_buffered_trees, BufferedTreeT *);
  LoopBelow (i, tb->max_buffered_trees)
    tb->buffered_tree[i] = NOBUFFEREDTREE;
  tb->min_tree_cost = INFINITY;
  tb->compressed = 0;
  return (tb);
}

BufferedTreeT *NextAvailableBufferedTree (TreeBufferT *b) {
  if (b->n_buffered_trees >= MAXTREES)
    return (NOBUFFEREDTREE);
  if (b->n_buffered_trees >= b->max_buffered_trees) {
    int i, m = b->max_buffered_trees;
    b->max_buffered_trees = Min (MAXTREES, b->max_buffered_trees * 2);
    fprintf (stderr, "\nIncreasing buffer size to %d\n", b->max_buffered_trees);
    b->buffered_tree = (BufferedTreeT **)Reallocate_D ((void *) b->buffered_tree,
						     b->max_buffered_trees * sizeof (BufferedTreeT *));
    LoopInterval (i, m, b->max_buffered_trees - 1)
      b->buffered_tree[i] = NOBUFFEREDTREE;
  }
  if (b->buffered_tree[b->n_buffered_trees] is NOBUFFEREDTREE) {
    b->buffered_tree[b->n_buffered_trees] = NewBufferedTree ();
  }
  return (b->buffered_tree[b->n_buffered_trees++]);
}

int NCandidateTrees = 0, NCandidatesRejected = 0, NCandidatesDiscarded = 0, NTreesReceived = 0;
int NServerDiscards = 0;

void BufferTree (TreeBufferT *b, int n_placed_taxa, int cost, BufferedTreeT *bt) {
  int i, j, old_cost;
  BufferedTreeT *t;
  if (bt is NOBUFFEREDTREE)
    NCandidateTrees++;
  if (bt) {

    SupportedCladesHash = bt->supported_clades_hash;
        NSupportedClades = bt->n_supported_clades;
    LoopBelow (i, NSupportedClades) {
      SupportedClades[i] = bt->supported_clades[i];
    }
  } else {
    ComputeSupportedClades ();
  }

  if (b->n_buffered_trees >= MAXTREES) {
    if (cost < b->min_tree_cost) {
      b->min_tree_cost = cost;
      /*fprintf (blather, "->%d", cost);*/
      CompressTreeBuffer (b);
    } else if (! b->compressed) {
      CompressTreeBuffer (b);
    }
    if (b->n_buffered_trees >= MAXTREES) {
      if (! maxtrees_exceeded) {
	/*fprintf (blather, "\nmaxtrees = %d exceeded\n", MAXTREES);*/
	maxtrees_exceeded = 1;
      }
      NCandidatesDiscarded++;
      return;
    }
  }

  LoopBelow (i, b->n_buffered_trees) {
    t = b->buffered_tree[i];
    if (t->supported_clades_hash is SupportedCladesHash &&
	t->n_supported_clades is NSupportedClades &&
	t->cost is cost) {
      LoopBelow (j, NSupportedClades) {
	if (CompareBitVectors (t->supported_clades[j], SupportedClades[j])) {
	  break;
	}
      }
      if (j >= NSupportedClades) {
	NCandidatesRejected++;
	return;
      }
    }
  }

/*   fprintf (stderr, "BufferTree storing tree with %d taxa\n", n_placed_taxa); */
  t = NextAvailableBufferedTree (b);
  if (bt) {
    assert (bt->n_placed_taxa is n_placed_taxa);
    CopyBufferedTree (t, bt);
  } else
    StoreTree (t, n_placed_taxa, cost);
  b->compressed = 0;
  t->generation = b->current_generation;
  old_cost = b->min_tree_cost;
  if (cost < old_cost) {
    b->min_tree_cost = cost;
/*    fprintf (blather, "->%d", cost);*/
  }
}

#ifdef _CLYDE_
void TransmitTreeToParent (TreeBufferT *b, int n_placed_taxa, int cost, BufferedTreeT *bt) {
  int me, parent;
  if (b->n_buffered_trees < MAXTREES ||
      cost < b->min_tree_cost) {
    /*   fprintf (stderr, "transmitting tree with %d taxa\n", n_placed_taxa); */
    assert (bt is NOBUFFEREDTREE);
    me = pvm_mytid ();
    parent = pvm_parent ();
    ComputeSupportedClades ();
    StoreTree (ScratchBufferedTree, n_placed_taxa, cost);
    pvm_initsend (PvmDataDefault);
    PackBufferedTree (ScratchBufferedTree, 1);
    assert (! pvm_send (pvm_parent (), RETURNTREE));
    if (ScratchBufferedTree->cost < b->min_tree_cost) {
      b->min_tree_cost = ScratchBufferedTree->cost;
    }
  } else {
    NServerDiscards++;
  }
}
#endif

#ifdef _CLYDE_
void TellParentDone () {
  pvm_initsend (PvmDataDefault);
  assert (! pvm_send (pvm_parent (), LOOPISDONE));
}
#endif

void DoNothing () {
}

void (*LOOPDONE) () = DoNothing;

void (*BUFFERTREE) (TreeBufferT *b, int n_placed_taxa, int cost, BufferedTreeT *bt) = BufferTree;

void ResetTreeBuffer (TreeBufferT *b) {
  b->n_buffered_trees = 0;
  b->min_tree_cost = INFINITY;
  b->compressed = 1;
}

void ResetGeneration (TreeBufferT *b) {
  int i;
  b->current_generation = 0;
  LoopBelow (i, b->n_buffered_trees)
    b->buffered_tree[i]->generation = 0;
}

void NextGeneration (TreeBufferT *b) {
 /* fprintf (blather, ";");*/
  b->current_generation++;
}

int NodeCost (NodeT *n)
{
  int character, c = 0;
  LoopBelow (character, NCharacters)
    if (n->character[character].union_bit) c++;
  return (c);
}

#define costof(c) ((c&UNIONBIT)? 1 : 0)
#define totalcostof(p) (p? p->totalcost: 0)

void SaveNode (NodeT *n)
{
  NodeT *b = n->backupnode;
  b->localcost = n->localcost;
  b->totalcost = n->totalcost;
  memcpy (b->character, n->character, NCharacters * sizeof (look_up_thang));
  n->dirty = 1;
  DirtyNode[NDirtyNodes++] = n;
/*   fprintf (stderr, "save node %s cost %d\n", n->name, n->totalcost); */
}

void RestoreNodes ()
{
  int i;
  NodeT *n, *b;
/*   fprintf (stderr, "restoring %d dirty nodes: ", NDirtyNodes); */
  LoopBelow (i, NDirtyNodes) {
    n = DirtyNode[i];
    b = n->backupnode;
/*     fprintf (stderr, "%s (cost %d->%d) ", n->name, n->totalcost, b->totalcost); */
    n->localcost = b->localcost;
    n->totalcost = b->totalcost;
    memcpy (n->character, b->character, NCharacters * sizeof (look_up_thang));
    n->dirty = 0;
  }
/*   fprintf (stderr, "\n"); */
  NDirtyNodes = 0;
}

#define BILLION 1000000000

int Repair (NodeT *node, int force, BackupT backup, int grandtotalcost, TreeBufferT *buffer)
{
  int character, count, highcount, newtotalcost, d, retval;
  look_up_thang newcharacter, oldcharacter;
/*  unsigned int newcharacter, oldcharacter;*/
  NodeT *n;
  static int lasthighcount = -1;
  int ii,jj;
  parameters *valuesI;

/*   fprintf (blather, "repair node %s\n", node?  node->name: "null"); */
  assert ((int) backup is NOBACKUP || backup is FROMBACKUP || backup is TOBACKUP);
  if (!node) {
/*     fprintf (stderr, "Repair NULL\n"); */
    if (BACKUP && backup is FROMBACKUP)
      RestoreNodes ();
    if (PRUNE && backup is FROMBACKUP)
      lasthighcount = -1;
    return (0);
  } else {
    if (BACKUP && backup is FROMBACKUP) {
      RestoreNodes ();
    } else {
      highcount = -1;
      if (PRUNE && backup is FROMBACKUP && lasthighcount > force)
	force = lasthighcount;
      character=0;
for (ii=0;ii<max(1,values_here->number_of_input_alignments);ii++) {
    if (!values_here->other_parm) valuesI=values_here;
    else if (values_here->other_parm[ii]) valuesI=values_here->other_parm[ii];
    else valuesI=values_here;
    /*fprintf(stderr,"%d->%d ",values_here->start[ii],values_here->stop[ii]);*/
    for (jj=values_here->start[ii];jj<values_here->stop[ii];jj++) {
/* LoopBelow (character, NCharacters) {*/
	for (n = node, count = 0; n; n = n->parent, count++) {
	  /*NOptimizations++;
	  newcharacter = SETBITS & n->left->character[character] & n->right->character[character];
	  if (! newcharacter) newcharacter = UNIONBIT | n->left->character[character] | n->right->character[character];
	  */
      newcharacter=valuesI->lookup[n->left->character[character].base][n->right->character[character].base];
	  oldcharacter = n->character[character];
	  	if (count < force || (*((unsigned int *) &newcharacter)) isnt *((unsigned int *) &oldcharacter)) {
        /*if (count < force || newcharacter isnt oldcharacter) {*/
	    if (BACKUP && backup is TOBACKUP && !n->dirty) SaveNode (n);
	      d=newcharacter.cost - oldcharacter.cost;
	      n->localcost += (cc_weights[character]*d);
	    /*d = (costof (newcharacter) - costof (oldcharacter));
	    n->localcost += d;*/
	    n->character[character] = newcharacter;
	    if (count > highcount) highcount = count;
	    if (backup is TOBACKUP) {
	      assert ((int) buffer);
	      grandtotalcost += (cc_weights[character]*d);
	    }
	  }
	else break;
	}
	if (backup is TOBACKUP && PRUNE && grandtotalcost > buffer->min_tree_cost + EPSILON) return (grandtotalcost);
  /*}character loop*/
  	  character ++;
	} /*jj*/
  } /*end of characters loop*/
    if (PRUNE && backup is TOBACKUP && highcount > lasthighcount)
	lasthighcount = highcount;
      for (n = node, count = 0;
	   n;
	   n = n->parent, count++) {
	newtotalcost = n->localcost + totalcostof (n->left) + totalcostof (n->right);
	if (count <= highcount ||
	    newtotalcost isnt n->totalcost) {
/* 	  if (BACKUP && backup is TOBACKUP) { */
/* 	    	  fprintf (stderr, "node %s is %sdirty\n", n->name, n->dirty? "": "not "); */
/* 		} */
	  if (BACKUP && backup is TOBACKUP && !n->dirty)
	    SaveNode (n);
	  n->totalcost = newtotalcost;
/* 	  fprintf (stderr, "%s->totalcost = %d\n", n->name, n->totalcost); */
	  if (PRUNE && backup is TOBACKUP && count > lasthighcount)
	    lasthighcount = count;
	} else {
	  break;
	}
      }
    }
  }
  retval = (FindRoot (node))->totalcost;
  if (backup is TOBACKUP)
    assert ((int) retval is grandtotalcost);
/*   fprintf (blather, "Repair returns %d\n", retval); */
  if (PRUNE && backup is FROMBACKUP)
    lasthighcount = -1;
  return (retval);
}

void ClearNode (NodeT *n)
{
  n->localcost = 0;
  memset (n->character, 0, NCharacters * sizeof (look_up_thang));
}

int correct_predictions = 0, incorrect_predictions = 0;

int AddSib (NodeT *old_node, NodeT *new_node, NodeT *HTU, BackupT backup, TreeBufferT *buffer)
{
  int grandtotalcost;
  NodeT *old_node_parent = old_node->parent;
  NodeT *n;

/*   fprintf (blather, "Add %s to %s via %s\n", new_node->name, old_node->name, HTU->name); */
  ClearNode (HTU);
  n = FindRoot (old_node);
  grandtotalcost = n->totalcost + new_node->totalcost + HTU->localcost;
/*   fprintf (stderr, "grandtotalcost %d = %s(%d) + %s(%d) + %s(%d)\n", */
/* 	   grandtotalcost, n->name, n->totalcost, */
/* 	   new_node->name, new_node->totalcost, */
/* 	   HTU->name, HTU->localcost); */
  HTU->parent = old_node->parent;
  HTU->left = old_node;
  HTU->right = new_node;
  old_node->parent = new_node->parent = HTU;
  if (old_node_parent) {
    if (old_node_parent->left is old_node) old_node_parent->left = HTU;
    else if (old_node_parent->right is old_node) old_node_parent->right = HTU;
    else fail ("Bogus adoption in AddSib");
  }
  return (Repair (HTU, 1, backup, grandtotalcost, buffer));
}

NodeT *RemoveAsSib (NodeT *node, BackupT backup)
{
  NodeT *parent, *sibling, *grandparent;
  parent = node->parent;
  assert (parent);
  if (node is parent->left)
    sibling = parent->right;
  else
    sibling = parent->left;
  grandparent = parent->parent;
/*   fprintf (blather, "Remove %s and %s from %s leaving %s\n", node->name, parent->name, grandparent? grandparent->name: "none", sibling->name); */
  parent->parent = parent->left = parent->right = NONODE;
  node->parent = NONODE;
  sibling->parent = grandparent;
  if (grandparent) {
    if (grandparent->left is parent) grandparent->left = sibling;
    else if (grandparent->right is parent) grandparent->right = sibling;
    else fail ("Bogus adoption in RemoveAsSib");
  }
  Repair (grandparent, 0, backup, 0, 0);
  return (parent);
}

#ifdef PARANOID
#define FindSlot(node,parent) (node is parent->left? &(parent->left): (node is parent->right? &(parent->right): (fail ("can't find slot"),&(parent->left))))
#define FindOtherSlot(node,parent) (node is parent->left? &(parent->right): (node is parent->right? &(parent->left): (fail ("can't find other slot"),&(parent->left))))
#define FindSib(node,parent) (node is parent->left? (parent->right): (node is parent->right? (parent->left): (fail ("can't find other slot"),(parent->left))))
#else
#define FindSlot(node,parent) (node is parent->left? &(parent->left):  &(parent->right))
#define FindOtherSlot(node,parent) (node is parent->left? &(parent->right): &(parent->left))
#define FindSib(node,parent) (node is parent->left? (parent->right): (parent->left))
#endif

#ifdef _CLYDE_
void ReceiveTree (TreeBufferT *tb) {
  int bytes, tag, tid, thiscost, i;
  assert (! pvm_bufinfo (pvm_recv (-1, -1), &bytes, &tag, &tid));
  switch (tag) {
  case RETURNTREE:
    NTreesReceived++;
    UnpackBufferedTree (ScratchBufferedTree, 1);
    thiscost = ScratchBufferedTree->cost;
/*      fprintf (stderr, "received tree with %d taxa\n", ScratchBufferedTree->n_placed_taxa); */
    if (thiscost < tb->min_tree_cost) {
      int i;
      pvm_initsend (PvmDataDefault);
      PackInt (buffertick);
      PackInt (thiscost);
      PackInt (tb->compressed? tb->n_buffered_trees: 0);
      LoopBelow (i, n_workers)
 	if (task[i].tid isnt tid)
	  assert (! pvm_send (task[i].tid, NEWMINTREECOST));
    }
    if (thiscost <= tb->min_tree_cost + EPSILON) {
/*        fprintf (stderr, "buffering it\n");  */
      (*BUFFERTREE) (tb, ScratchBufferedTree->n_placed_taxa, thiscost, ScratchBufferedTree);
    } else {
/*       fprintf (stderr, "cost of %d is too high compared to %d.\n", thiscost, ScratchBufferedTree->cost); */
    }
    break;
  case LOOPISDONE:
    LoopBelow (i, n_workers)
      if (task[i].tid is tid)
	break;
    assert (i < n_workers);
    task[i].ntasks--;
    tasks_available++;
    next_task = i;
    break;
  default:
    assert (0);
  }
}
#endif

#ifdef _CLYDE_
int GetAServer (TreeBufferT *tb) {
  int tid, i;
  while (tasks_available is 0 || pvm_probe (-1, -1) > 0) {
    ReceiveTree (tb);
  }
  for (i = 0; i < n_workers && task[next_task].ntasks is DOUBLEBUFFER; i++)
    next_task = (next_task + 1) % n_workers;
  assert (i < n_workers);
  task[next_task].ntasks++;
  tasks_available--;
  tid = task[next_task].tid;
  next_task = (next_task + 1) % n_workers;
  return (tid);
}
#endif

#ifdef _CLYDE_
void WaitForServers (TreeBufferT *tb) {
  int i;
  while (tasks_available isnt max_tasks_available)
    ReceiveTree (tb);
  buffertick++;
  pvm_initsend (PvmDataDefault);
  LoopBelow (i, n_workers)
    assert (! pvm_send (task[i].tid, BUFFERTICK));
}
#endif

void IgnoreTreeBuffer (TreeBufferT *tb) {
}

void (*WAITFORSERVERS) (TreeBufferT *tb) = IgnoreTreeBuffer;

int SERVERTASKSWITH = 0, SERVERTASKSWITHOUT = 0;

#ifdef _CLYDE_
void PackTreeForServer (TreeBufferT *tb, int tid) {
  int taskindex;
  LoopBelow (taskindex, n_workers)
    if (task[taskindex].tid is tid) break;
  assert (taskindex < n_workers);
  pvm_initsend (PvmDataDefault);
  PackInt (tb->min_tree_cost);
  PackInt (tb->compressed? tb->n_buffered_trees: 0);
  if (treeserial isnt task[taskindex].treeid) {
    PackInt (1);
    SERVERTASKSWITH++;
    StoreTree (ScratchBufferedTree, NPlacedTaxa, tb->min_tree_cost);
    PackBufferedTree (ScratchBufferedTree, 0);
    task[taskindex].treeid = treeserial;
  } else {
    SERVERTASKSWITHOUT++;
    PackInt (0);
  }
}
#endif

#ifdef _CLYDE_
void UpdateMinTreeCost (TreeBufferT *tb) {
  int i;
  while (pvm_nrecv (-1, NEWMINTREECOST) > 0)
    if (UnpackNextInt () is buffertick) {
      i = UnpackNextInt ();
      tb->min_tree_cost = Min (tb->min_tree_cost, i);
      tb->n_buffered_trees = UnpackNextInt ();
    }
}
#endif

void (*UPDATEMINTREECOST) (TreeBufferT *tb) = IgnoreTreeBuffer;

void LocalRearrangementLoop (TreeBufferT *tb, int firstplacedtaxon, int lastplacedtaxon)
{
  int thiscost, placedtaxon;
  NodeT *placednode, *parent,
  *t1, *t2, *t3, **t3slot;

  LoopInterval (placedtaxon, firstplacedtaxon, lastplacedtaxon) {
    placednode = PlacedTaxa[placedtaxon];
    parent = placednode->parent;
    t1 = placednode->left;
    t2 = placednode->right;
    if (parent && parent->parent && t1) {
      assert (t2);
      t3slot = FindOtherSlot (placednode, parent);
      t3 = *t3slot;
      *t3slot = t1;
      t1->parent = parent;
      placednode->left = t3;
      t3->parent = placednode;
      thiscost = Repair (placednode, 1, TOBACKUP, (FindRoot (parent))->totalcost, tb);
      if (thiscost <= tb->min_tree_cost + EPSILON) {
	(*BUFFERTREE) (tb, NPlacedTaxa, thiscost, NOBUFFEREDTREE);
      }
      NTrees++;
      *t3slot = t2;
      t2->parent = parent;
      placednode->left = t1;
      t1->parent = placednode;
      placednode->right = t3;
      t3->parent = placednode;
      thiscost = Repair (placednode, 1, TOBACKUP, thiscost, tb);
      (*UPDATEMINTREECOST) (tb);
      if (thiscost <= tb->min_tree_cost + EPSILON) {
	(*BUFFERTREE) (tb, NPlacedTaxa, thiscost, NOBUFFEREDTREE);
      }
      NTrees++;
      *t3slot = t3;
      t3->parent = parent;
      placednode->right = t2;
      t2->parent = placednode;
      Repair (placednode, 1, FROMBACKUP, 0, 0);
    }
  }
  (*LOOPDONE) ();
}

#ifdef _CLYDE_
void LocalRearrangementToServer (TreeBufferT *tb, int firstplacedtaxon, int lastplacedtaxon) {
  int tid = GetAServer (tb);
  PackTreeForServer (tb, tid);
  PackInt (firstplacedtaxon);
  PackInt (lastplacedtaxon);
  assert (! pvm_send (tid, LOCALREARRANGEMENTLOOP));
}
#endif

void (*LOCALREARRANGEMENTLOOPFUN) (TreeBufferT *tb, int firstplacedtaxon, int lastplacedtaxon) =
     LocalRearrangementLoop;


void LocalRearrangement (TreeBufferT *tb)
{
  int placedtaxon, lastplacedtaxon;
  int tree, previous_generation, oldprune;

  oldprune = PRUNE; PRUNE = 0;
  ResetGeneration (tb);
  while (tb->n_buffered_trees > 0 &&
	 tb->buffered_tree[tb->n_buffered_trees - 1]->generation is
	 tb->current_generation) {
    NextGeneration (tb);
    previous_generation = tb->current_generation - 1;
    for (tree = tb->n_buffered_trees - 1;
	 tree >= 0 &&
	 tb->buffered_tree[tree]->generation is previous_generation;
	 tree--) {
      LoadTree (tb->buffered_tree[tree], 0);
      treeserial++;
      for (placedtaxon = 2; placedtaxon <= NPlacedTaxa - 1; placedtaxon += LOOPCHUNK) {
	lastplacedtaxon = Min (NPlacedTaxa - 1, placedtaxon + LOOPCHUNK - 1);
	(*LOCALREARRANGEMENTLOOPFUN) (tb, placedtaxon, lastplacedtaxon);
      }
    }
    (*WAITFORSERVERS) (tb);
    CompressTreeBuffer (tb);
  }
  PRUNE = oldprune;
}

void ClearMarks ()
{
  int taxon;
  LoopBelow (taxon, NTaxa)
    Node[taxon].mark = 0;
}

void MarkTree (NodeT *n, int mark)
{
  if (n) {
    n->mark = mark;
    MarkTree (n->left, mark);
    MarkTree (n->right, mark);
  }
}



void ClearPartialFinal () {
  int i;
  LoopBelow (i, NCharacters)
    PartialFinalMap[i] = 0;
}

void ComputePartialFinal (int i, NodeT *n) {
  look_up_thang n_ch;
  NodeT *parent = n->parent;
  n_ch = n->character[i];
  n_ch.union_bit=0;
  if (parent is NONODE)
    n->partial_final_character[i].base = n_ch.base;
  else if (((~n_ch.base) & parent->partial_final_character[i].base) is 0) n->partial_final_character[i].base = parent->partial_final_character[i].base;
  else if (n->character[i].union_bit) n->partial_final_character[i].base = n_ch.base | parent->partial_final_character[i].base;
  else if (n->left)
    n->partial_final_character[i].base =
      n_ch.base |
	(parent->partial_final_character[i].base &
	 (n->left->character[i].base |
	  n->right->character[i].base));
  else
    n->partial_final_character[i].base = n_ch.base;
  n->partial_final_character[i].union_bit=0;
  if (n->left) {
    ComputePartialFinal (i, n->left);
    ComputePartialFinal (i, n->right);
  }
}

int GetAddCost (int basecost, NodeT *graftnode, NodeT *prunedclade, NodeT *graftroot, TreeBufferT *buffer) {
  int i,ii,jj;
  parameters *valuesI;
  NodeT *parent = graftnode->parent;
  i=0;
    for (ii=0;ii<max(1,values_here->number_of_input_alignments);ii++) {
        if (!values_here->other_parm) valuesI=values_here;
        else if (values_here->other_parm[ii]) valuesI=values_here->other_parm[ii];
        else valuesI=values_here;
        for (jj=values_here->start[ii];jj<values_here->stop[ii];jj++) {
            if (PRUNE && basecost > buffer->min_tree_cost + EPSILON) break;
            if (! PartialFinalMap[i]) {
              ComputePartialFinal (i, graftroot);
              PartialFinalMap[i] = 1;
            }
            if (valuesI->lookup[prunedclade->character[i].base][graftnode->partial_final_character[i].base].union_bit || valuesI->lookup[prunedclade->character[i].base][parent->partial_final_character[i].base].union_bit)
            /*if (((prunedclade->character[i].base & graftnode->partial_final_character[i].base) is 0) &&
            	((prunedclade->character[i].base & parent->partial_final_character[i]).base is 0))*/
            basecost+=(min(valuesI->lookup[prunedclade->character[i].base][graftnode->partial_final_character[i].base].cost,valuesI->lookup[prunedclade->character[i].base][parent->partial_final_character[i].base].cost)*cc_weights[i]);
      	    i++;
    	 } /*jj*/
      } /*end of characters loop*/
  return (basecost);
}

void SPRLoop (TreeBufferT *tb, int firstprunedtaxon, int lastprunedtaxon, int oldcost)
{
  int thiscost, grafttaxon, prunedtaxon, reallydidit, basecost, predcost;
  NodeT *prunedclade, *graftnode, *parent, *sibling, *HTU, *graftroot;
  LoopInterval (prunedtaxon, firstprunedtaxon, lastprunedtaxon) {
    prunedclade = PlacedTaxa[prunedtaxon];
    parent = prunedclade->parent;
    if (parent && parent->parent) {
      sibling = FindSib (prunedclade, parent);
      HTU = RemoveAsSib (prunedclade, NOBACKUP);
      ClearMarks ();
      MarkTree (prunedclade, 1);
      HTU->mark = 1;
      graftroot = FindRoot (sibling);
      basecost = prunedclade->totalcost + graftroot->totalcost;
      if (PABLO) {
	ClearPartialFinal ();
      }
      if (basecost <= tb->min_tree_cost + EPSILON)
      LoopInterval (grafttaxon, 2, NPlacedTaxa - 1) {
	graftnode = PlacedTaxa[grafttaxon];
	if (! graftnode->mark) {
	  if (PABLO) {
	    thiscost = GetAddCost (basecost, graftnode, prunedclade, graftroot, tb);
	    reallydidit = 0;
	  } else {
	    thiscost = AddSib (graftnode, prunedclade, HTU, TOBACKUP, tb);
	    reallydidit = 1;
	  }
	  NTrees++;
	  (*UPDATEMINTREECOST) (tb);
	  if (thiscost <= tb->min_tree_cost + EPSILON) {
	    if (! reallydidit) {
	      predcost = thiscost;
	      thiscost = AddSib (graftnode, prunedclade, HTU, TOBACKUP, tb);
	      if (PABLO && predcost isnt thiscost)
		fprintf (stderr, "discrepancy %d %d %d\n", grafttaxon, predcost, thiscost);
	      reallydidit = 1;
	    }
	    (*BUFFERTREE) (tb, NPlacedTaxa, thiscost, NOBUFFEREDTREE);
	  }
	  if (reallydidit) {
	    RemoveAsSib (prunedclade, FROMBACKUP);
	  }
	}
	if (STEEPEST && tb->min_tree_cost isnt oldcost)
	  break;
      }
      AddSib (sibling, prunedclade, HTU, NOBACKUP, 0);
    }
    if (STEEPEST && tb->min_tree_cost isnt oldcost)
      break;
  }  (*LOOPDONE) ();
}

#ifdef _CLYDE_
void SPRToServer (TreeBufferT *tb, int firstprunedtaxon, int lastprunedtaxon, int oldcost) {
  int tid = GetAServer (tb);
  PackTreeForServer (tb, tid);
  PackInt (firstprunedtaxon);
  PackInt (lastprunedtaxon);
  PackInt (oldcost);
  assert (! pvm_send (tid, SPRLOOP));
}
#endif

void (*SPRLOOPFUN) (TreeBufferT *tb, int firstprunedtaxon, int lastprunedtaxon, int oldcost) = SPRLoop;

void SPR (TreeBufferT *tb)
{
  int prunedtaxon, lastprunedtaxon, tree, previous_generation;
  int oldcost;
  ResetGeneration (tb);
  while (tb->n_buffered_trees > 0 &&
	 tb->buffered_tree[tb->n_buffered_trees - 1]->generation is
	 tb->current_generation) {
    NextGeneration (tb);
    previous_generation = tb->current_generation - 1;
    oldcost = tb->min_tree_cost;
    for (tree = tb->n_buffered_trees - 1;
	 tree >= 0 &&
	 tb->buffered_tree[tree]->generation is previous_generation;
	 tree--) {
      LoadTree (tb->buffered_tree[tree], 0);
      treeserial++;
      for (prunedtaxon = 2; prunedtaxon <= NPlacedTaxa - 1; prunedtaxon += LOOPCHUNK) {
	lastprunedtaxon = Min (NPlacedTaxa - 1, prunedtaxon + LOOPCHUNK - 1);
	(*SPRLOOPFUN) (tb, prunedtaxon, lastprunedtaxon, oldcost);
	if (STEEPEST && tb->min_tree_cost isnt oldcost)
	  break;
      }
      if (STEEPEST && tb->min_tree_cost isnt oldcost)
	break;
    }
    (*WAITFORSERVERS) (tb);
    CompressTreeBuffer (tb);
  }
}

NodeT *Reroot (NodeT *alpha)
{
  NodeT *alphaparent, *node, *prevnode, *parent, **prevslot, *root, *sib;
  assert (alpha);
  assert (alpha->parent);
  assert (alpha->parent->parent);
  alphaparent = alpha->parent;
  node = alphaparent;
  prevnode = alpha;
  while (1) {
    parent = node->parent;
    prevslot = FindSlot (prevnode, node);
    node->parent = prevnode;
    if (parent->parent) {
      *prevslot = parent;
      prevnode = node;
      node = parent;
    } else {
      root = parent;
      sib = FindSib (node, root);
      *prevslot = sib;
      sib->parent = node;
      root->left = alpha;
      root->right = alphaparent;
      alpha->parent = root;
      alphaparent->parent = root;
      Repair (node, INFINITY, NOBACKUP, 0, 0);
      return (*prevslot);
    }
  }
}

void TBRLoop (TreeBufferT *tb, int firstprunedtaxon, int lastprunedtaxon, int oldcost)
{
  int thiscost, grafttaxon, prunedtaxon, reallydidit, basecost, predcost;
  NodeT *prunedclade, *graftnode, *parent, *sibling, *HTU, *graftroot;
  NodeT *rerootnode, *unrerootnode, *currentroot;
  int reroottaxon;

  LoopInterval (prunedtaxon, firstprunedtaxon, lastprunedtaxon) {
    prunedclade = PlacedTaxa[prunedtaxon];
    parent = prunedclade->parent;
    if (parent && parent->parent) {
      sibling = FindSib (prunedclade, parent);
      HTU = RemoveAsSib (prunedclade, NOBACKUP);
      ClearMarks ();
      MarkTree (prunedclade, 1);
      HTU->mark = 1;
      graftroot = FindRoot (sibling);
      basecost = prunedclade->totalcost + graftroot->totalcost;
      if (PABLO) {
	ClearPartialFinal ();
      }
      if (basecost <= tb->min_tree_cost + EPSILON)
      LoopInterval (reroottaxon, 2, NPlacedTaxa - 1) {
	rerootnode = PlacedTaxa[reroottaxon];
	if (rerootnode->mark && rerootnode->parent && rerootnode->parent->parent) {
	  unrerootnode = Reroot (rerootnode);
	  currentroot = FindRoot (rerootnode);
	  LoopInterval (grafttaxon, 2, NPlacedTaxa - 1) {
	    graftnode = PlacedTaxa[grafttaxon];
	    if (! graftnode->mark) {
	      if (PABLO) {
		thiscost = GetAddCost (basecost, graftnode, currentroot, graftroot, tb);
		reallydidit = 0;
	      } else {
		thiscost = AddSib (graftnode, currentroot, HTU, TOBACKUP, tb);
		reallydidit = 1;
	      }
	      NTrees++;
	      (*UPDATEMINTREECOST) (tb);
	      if (thiscost <= tb->min_tree_cost + EPSILON) {
		if (! reallydidit) {
		  predcost = thiscost;
		  thiscost = AddSib (graftnode, currentroot, HTU, TOBACKUP, tb);
		  if (PABLO && predcost isnt thiscost)
		    fprintf (stderr, "discrepancy %d %d %d\n", grafttaxon, predcost, thiscost);
		  reallydidit = 1;
		}
		(*BUFFERTREE) (tb, NPlacedTaxa, thiscost, NOBUFFEREDTREE);
	      }
	      if (reallydidit) {
		RemoveAsSib (prunedclade, FROMBACKUP);
	      }
	    }
	    if (STEEPEST && tb->min_tree_cost isnt oldcost)
	      break;
	  }
	  Reroot (unrerootnode);
	}
	if (STEEPEST && tb->min_tree_cost isnt oldcost)
	  break;
      }
      AddSib (sibling, prunedclade, HTU, NOBACKUP, 0);
    }
    if (STEEPEST && tb->min_tree_cost isnt oldcost)
      break;
  }
  (*LOOPDONE) ();
}

#ifdef _CLYDE_
void TBRToServer (TreeBufferT *tb, int firstprunedtaxon, int lastprunedtaxon, int oldcost) {
  int tid = GetAServer (tb);
  PackTreeForServer (tb, tid);
  PackInt (firstprunedtaxon);
  PackInt (lastprunedtaxon);
  PackInt (oldcost);
  assert (! pvm_send (tid, TBRLOOP));
}
#endif

void (*TBRLOOPFUN) (TreeBufferT *tb, int firstprunedtaxon, int lastprunedtaxon, int oldcost) = TBRLoop;

void TBR (TreeBufferT *tb)
{
  int prunedtaxon, lastprunedtaxon, tree, previous_generation;
  int oldcost;
  ResetGeneration (tb);
  while (tb->n_buffered_trees > 0 &&
	 tb->buffered_tree[tb->n_buffered_trees - 1]->generation is
	 tb->current_generation) {
    NextGeneration (tb);
    previous_generation = tb->current_generation - 1;
    oldcost = tb->min_tree_cost;
    for (tree = tb->n_buffered_trees - 1;
	 tree >= 0 &&
	 tb->buffered_tree[tree]->generation is previous_generation;
	 tree--) {
      LoadTree (tb->buffered_tree[tree], 0);
      treeserial++;
      for (prunedtaxon = 2; prunedtaxon <= NPlacedTaxa - 1; prunedtaxon += LOOPCHUNK) {
	lastprunedtaxon = Min (NPlacedTaxa - 1, prunedtaxon + LOOPCHUNK - 1);
	(*TBRLOOPFUN) (tb, prunedtaxon, lastprunedtaxon, oldcost);
	if (STEEPEST && tb->min_tree_cost isnt oldcost)
	  break;
      }
      if (STEEPEST && tb->min_tree_cost isnt oldcost)
	break;
    }
    (*WAITFORSERVERS) (tb);
    CompressTreeBuffer (tb);
  }
}

void PrintTreeNode (FILE *f, NodeT *n)
{
  if (n->left is NONODE) {
    fprintf (f, "%s", n->name);
  } else {
    int supported = n->parent is NONODE || IsSupported (n);
    if (supported)
      fprintf (f, "(");
    PrintTreeNode (f, n->left);
    fprintf (f, " ");
    PrintTreeNode (f, n->right);
    if (supported)
      fprintf (f, ")");
  }
}

void PrintTree (FILE *f)
{
  NodeT *root = FindRoot (Node + 0);
  PrintTreeNode (f, root);
}

void SkipSpaces (FILE *f, int spaces)
{
  while (spaces--)
    fprintf (f, " ");
}

#define SPACESPERLEVEL 2

void IndentTreeNode (FILE *f, NodeT *n, int spaces)
{
  if (n->left is NONODE) {
    SkipSpaces (f, spaces);
    fprintf (f, "%s\n", n->name);
  } else {
    int newindent = IsSupported (n)? spaces + SPACESPERLEVEL: spaces;
    IndentTreeNode (f, n->left, newindent);
    IndentTreeNode (f, n->right, newindent);
  }
}

void IndentTree (FILE *f)
{
  NodeT *root = FindRoot (Node + 0);
  IndentTreeNode (f, root, 0);
}

void PrintSupportedClades (FILE *f, BufferedTreeT *t) {
  int i, j;

  fprintf (f, "%d supported clades:\n", t->n_supported_clades);
  LoopBelow (i, t->n_supported_clades) {
    LoopBelow (j, NTaxa) {
      if (GetBit (t->supported_clades[i], j)) {
	fprintf (f, "%s ", Node[j].name);
      }
    }
    fprintf (f, "\n");
  }
}

void PrintSupportBits (FILE *f, BufferedTreeT *t) {
  int i;

  fprintf (f, "%d supported clades, hash %d:\n", t->n_supported_clades, t->supported_clades_hash);
  LoopBelow (i, t->n_supported_clades) {
    fprintf (f, "%3d ", i);
    fprintfBitVector (f, t->supported_clades[i]);
    fprintf (f, "\n");
  }
}

void PrintLoadedSupportBits (FILE *f) {
  int i;

  fprintf (f, "%d loaded-supported clades, hash %d:\n", NSupportedClades, SupportedCladesHash);
  LoopBelow (i, NSupportedClades) {
    fprintf (f, "%3d ", i);
    fprintfBitVector (f, SupportedClades[i]);
    fprintf (f, "\n");
  }
}

void PhylogenyAllocations () {
  int taxon, character;
  assert (NOTUs >= 3);
  InitBitVectorSize (NOTUs);
  LoopBelow (taxon, NHTUs) {
    sprintf (buf, "HTU %d", taxon);
    Name[taxon+NOTUs] = CopyString (buf);
  }
  NewArray (Node, NTaxa * 2, NodeT);
  NewArray (DirtyNode, NTaxa, NodeT *);
  NewArray (PlacedTaxa, NTaxa, NodeT *);
  LoopBelow (taxon, NTaxa) {
    Node[taxon].taxon = taxon;
    Node[taxon].name = Name[taxon];
    Node[taxon].character = Character[taxon];
    Node[taxon].parent = Node[taxon].left = Node[taxon].right = NONODE;
    if (taxon >= NOTUs)
      LoopBelow (character, NCharacters)
	Character[taxon][character].base =Character[taxon][character].union_bit =Character[taxon][character].cost = 0;
    Node[taxon].localcost = Node[taxon].totalcost = 0;
    Node[taxon].backupnode = &Node[taxon+NTaxa];
    Node[taxon+NTaxa].character = Character[taxon+NTaxa];
    Node[taxon].dirty = 0;
    Node[taxon].taxa = NewBitVector ();
    NewArray (Node[taxon].final_character, NCharacters, look_up_thang);
    NewArray (Node[taxon].partial_final_character, NCharacters, look_up_thang);
  }
  NewArray (PartialFinalMap, NCharacters, int);
  NewArray (SupportedClades, NTaxa, BitVectorT);
  if (SERVERMODE isnt STANDALONE)
    ScratchBufferedTree = NewBufferedTree ();
}

#ifdef _CLYDE_
void AdditionToServer (TreeBufferT *tb, int firsttaxon, int lasttaxon, int OTU_index, int HTU_index) {
  int tid = GetAServer (tb);
  PackTreeForServer (tb, tid);
  PackInt (firsttaxon);
  PackInt (lasttaxon);
  PackInt (OTU_index);
  PackInt (HTU_index);
  assert (! pvm_send (tid, ADDITIONLOOP));
}
#endif

void AdditionLoop (TreeBufferT *DestBuffer, int firsttaxon, int lasttaxon, int OTU_index, int HTU_index)
{
  int thiscost, placedtaxon, basecost, oldcost, reallydidit;
  NodeT *HTUnode, *OTUnode, *placednode, *root;
  OTUnode = &Node[OTU_index];
  HTUnode = &Node[HTU_index];
  OTUnode->parent = NONODE;
  PlacedTaxa[NPlacedTaxa] = OTUnode;
  PlacedTaxa[NPlacedTaxa+1] = HTUnode;
  if (INITIALPABLO) {
    root = FindRoot (PlacedTaxa[firsttaxon]);
    basecost = root->totalcost;
    ClearPartialFinal ();
  }
  LoopInterval (placedtaxon, firsttaxon, lasttaxon) {
    placednode = PlacedTaxa[placedtaxon];
    if (INITIALPABLO) {
      thiscost = GetAddCost (basecost, placednode, OTUnode, root, DestBuffer);
      reallydidit = 0;
    } else {
      thiscost = AddSib (placednode, OTUnode, HTUnode, TOBACKUP, DestBuffer);
      reallydidit = 1;
    }
    (*UPDATEMINTREECOST) (DestBuffer);
    if (thiscost <= DestBuffer->min_tree_cost + EPSILON) {
      if (! reallydidit) {
	oldcost = thiscost;
	thiscost = AddSib (placednode, OTUnode, HTUnode, TOBACKUP, DestBuffer);
	if (INITIALPABLO && oldcost isnt thiscost)
	  fprintf (stderr, "discrepancy %d %d\n", oldcost, thiscost);
	reallydidit = 1;
      }
      (*BUFFERTREE) (DestBuffer, NPlacedTaxa+2, thiscost, NOBUFFEREDTREE);
    }
    if (reallydidit)
      RemoveAsSib (OTUnode, FROMBACKUP);
    NTrees++;
  }
  (*LOOPDONE) ();
}

void (*ADDITIONLOOPFUN) (TreeBufferT *tb, int firsttaxon, int lasttaxon, int OTU_index, int HTU_index) = AdditionLoop;

TreeBufferT *Phylogeny ()
{
  int newHTU, newOTU, thiscost, placedtaxon, lastplacedtaxon, randomization;
  NodeT *HTUnode, *OTUnode;
  TreeBufferT *TreeBuffer1 = NewTreeBuffer(), *TreeBuffer2 = NewTreeBuffer(),
  *TreeBuffer3 = NewTreeBuffer(), *DestBuffer;
  int tree;
  int *OTUindex;
  MAXTREES = INITIAL_MAXTREES;
  NewArray (OTUindex, NOTUs, int);
  LoopBelow (randomization, Max (1, RANDOMIZATIONS)) {
    int i, j;
    LoopBelow (i, NOTUs)
      OTUindex[i] = i;
    if (RANDOMIZATIONS > 0) {
      int temp;
      srand ((unsigned int) time (0));
      LoopInterval (i, 1, NOTUs - 2) {
	j = i + rand () % (NOTUs - i);
	temp = OTUindex[i]; OTUindex[i] = OTUindex[j]; OTUindex[j] = temp;
      }
    }
   /* fprintf (blather, "Outgroup is %s, starting with node %s", Node[0].name, Node[OTUindex[1]].name);*/
    PlacedTaxa[0] = &Node[OTUindex[0]];
    PlacedTaxa[1] = &Node[NOTUs];
    PlacedTaxa[2] = &Node[OTUindex[1]];
    PlacedTaxa[0]->parent = PlacedTaxa[1]->parent = NONODE;
    thiscost = AddSib (PlacedTaxa[0], PlacedTaxa[2], PlacedTaxa[1], NOBACKUP, TreeBuffer1);
    NPlacedTaxa = 3;
    ResetTreeBuffer (TreeBuffer1);
    (*BUFFERTREE) (TreeBuffer1, NPlacedTaxa, thiscost, NOBUFFEREDTREE);
    /*fprintf (blather, "\n");*/
    LoopInterval (newOTU, 2, NOTUs - 1) {
      ResetTreeBuffer (TreeBuffer2);
      DestBuffer = (newOTU is NOTUs - 1)? TreeBuffer3: TreeBuffer2;
      newHTU = newOTU + NOTUs - 1;
      OTUnode = &Node[OTUindex[newOTU]];
      HTUnode = &Node[newHTU];
      OTUnode->parent = NONODE;
      /*fprintf (blather, "add node %s", OTUnode->name);*/
      LoopBelow (tree, TreeBuffer1->n_buffered_trees) {
	LoadTree (TreeBuffer1->buffered_tree[tree], 0);
	treeserial++;
	for (placedtaxon = 2; placedtaxon <= NPlacedTaxa - 1; placedtaxon += LOOPCHUNK) {
	  lastplacedtaxon = Min (NPlacedTaxa - 1, placedtaxon + LOOPCHUNK - 1);
	  (*ADDITIONLOOPFUN) (DestBuffer, placedtaxon, lastplacedtaxon, OTUindex[newOTU], newHTU);
	}
      }
      (*WAITFORSERVERS) (DestBuffer);
      CompressTreeBuffer (DestBuffer);
      if (LOCALREARRANGEMENT) {
	LocalRearrangement (DestBuffer);
      }
      /*fprintf (blather, " %d trees\n", DestBuffer->n_buffered_trees);*/
      {TreeBufferT *temp = TreeBuffer1; TreeBuffer1 = TreeBuffer2; TreeBuffer2 = temp;}
    }
  }
  if (INITIAL_MAXTREES isnt FINAL_MAXTREES) {
    MAXTREES = FINAL_MAXTREES;
    maxtrees_exceeded = 0;
  }
  if (doSPR) {
    /*fprintf (blather, "SPR ");*/
    SPR (TreeBuffer3);
   /* fprintf (blather, " %d trees\n", TreeBuffer3->n_buffered_trees);*/
  }
  if (doTBR) {
    /*fprintf (blather, "TBR ");*/
    TBR (TreeBuffer3);
   /* fprintf (blather, " %d trees\n", TreeBuffer3->n_buffered_trees);*/
  }
  return (TreeBuffer3);
}

char *Setting (int x) {
  return (x? "on": "off");
}

void ProcessCommandLineArguments ()
{
      BACKUP = 1;
      PRUNE = 1;
      doTBR = values_here->tbr;
      doSPR = values_here->sbr;
      LOCALREARRANGEMENT = values_here->clade_swap_while_add;
      VERBOSE=0;
      SupportComputer = ComputeSupportMust;
      NTrees = NOptimizations = 0;
      NDirtyNodes = 0;
      if (values_here->phylo_score==7) {EPSILON= -1;MAXTREES=1;}
      else {EPSILON=0;MAXTREES=values_here->keep_trees;}
      INITIAL_MAXTREES=MAXTREES;
      RANDOMIZATIONS=values_here->tree_rand_order_max;
/*  FILE *devnull = fopen ("/dev/null", "w");
  int serverargc;
  blather = devnull;
  NewArray (serverargv, argc+2, char *);
  serverargv[0] = "-server";
  serverargc = 1;*/
  SupportComputer = ComputeSupportMust;
/*  while (argc--) {
    if (streq (argv[0], "-backup")) {
      BACKUP = 1;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-nobackup")) {
      serverargv[serverargc++] = argv[0];
      BACKUP = 0;
    } else if (streq (argv[0], "-prune")) {
      serverargv[serverargc++] = argv[0];
      PRUNE = 1;
    } else if (streq (argv[0], "-noprune")) {
      serverargv[serverargc++] = argv[0];
      PRUNE = 0;
    } else if (streq (argv[0], "-spr")) {
      serverargv[serverargc++] = argv[0];
      doSPR = 1;
    } else if (streq (argv[0], "-nospr")) {
      doSPR = 0;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-tbr")) {
      doTBR = 1;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-notbr")) {
      doTBR = 0;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-local")) {
      LOCALREARRANGEMENT = 1;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-nolocal")) {
      LOCALREARRANGEMENT = 0;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-pablo")) {
      PABLO = 1;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-nopablo")) {
      PABLO = 0;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-steepest")) {
      STEEPEST = 1;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-nosteepest")) {
      STEEPEST = 0;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-initialpablo")) {
      INITIALPABLO = 1;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-noinitialpablo")) {
      INITIALPABLO = 0;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-tree")) {
      TREE = 1;
    } else if (streq (argv[0], "-notree")) {
      TREE = 0;
    } else if (streq (argv[0], "-indent")) {
      INDENT = 1;
    } else if (streq (argv[0], "-noindent")) {
      INDENT = 0;
    } else if (streq (argv[0], "-consensus")) {
      CONSENSUS = 1;
    } else if (streq (argv[0], "-noconsensus")) {
      CONSENSUS = 0;
    } else if (streq (argv[0], "-consensusbits")) {
      CONSENSUSBITS = 1;
    } else if (streq (argv[0], "-noconsensusbits")) {
      CONSENSUSBITS = 0;
    } else if (streq (argv[0], "-consensustable")) {
      CONSENSUSTABLE = 1;
    } else if (streq (argv[0], "-noconsensustable")) {
      CONSENSUSTABLE = 0;
    } else if (streq (argv[0], "-numberconsensusnodes")) {
      NUMBERCONSENSUSNODES = 1;
    } else if (streq (argv[0], "-nonumberconsensusnodes")) {
      NUMBERCONSENSUSNODES = 0;
    } else if (streq (argv[0], "-missingreport")) {
      MISSINGREPORT = 1;
    } else if (streq (argv[0], "-nomissingreport")) {
      MISSINGREPORT = 0;
    } else if (streq (argv[0], "-missingconsensus")) {
      MISSINGCONSENSUS = 1;
    } else if (streq (argv[0], "-nomissingconsensus")) {
      MISSINGCONSENSUS = 0;
    } else if (streq (argv[0], "-missingconsensusdata")) {
      MISSINGCONSENSUSDATA = 1;
    } else if (streq (argv[0], "-nomissingconsensusdata")) {
      MISSINGCONSENSUSDATA = 0;
    } else if (streq (argv[0], "-missingstats")) {
      MISSINGSTATS = 1;
    } else if (streq (argv[0], "-nomissingstats")) {
      MISSINGSTATS = 0;
    } else if (streq (argv[0], "-stats")) {
      STATS = 1;
    } else if (streq (argv[0], "-nostats")) {
      STATS = 0;
    } else if (streq (argv[0], "-tags")) {
      TAGS = 1;
    } else if (streq (argv[0], "-notags")) {
      TAGS = 0;
    } else if (streq (argv[0], "-verbose")) {
      VERBOSE = 1;
      blather = stderr;
    } else if (streq (argv[0], "-noverbose")) {
      VERBOSE = 0;
      blather = devnull;
    } else if (streq (argv[0], "-standalone")) {
      SERVERMODE = STANDALONE;
    } else if (streq (argv[0], "-master")) {
      SERVERMODE = MASTER;
    } else if (streq (argv[0], "-server")) {
      SERVERMODE = SERVER;
    } else if (streq (argv[0], "-support")) {
      serverargv[serverargc++] = argv[0];
      argv++;
      argc--;
      if (streq (argv[0], "must"))
	SupportComputer = ComputeSupportMust;
      else if (streq (argv[0], "everything"))
	SupportComputer = ComputeSupportEverything;
      else if (streq (argv[0], "may"))
	SupportComputer = ComputeSupportMay;
      else if (streq (argv[0], "notsamesingleton"))
	SupportComputer = ComputeSupportNotSameSingleton;
      else {
	fprintf (stderr, "bad -support argument: %s\n", argv[0]);
	exit (-1);
      }
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-format")) {
      serverargv[serverargc++] = argv[0];
      argv++;
      argc--;
      if (streq (argv[0], "character"))
	FORMAT = CHARACTER;
      else if (streq (argv[0], "bitmap"))
	FORMAT = BITMAP;
      else {
	fprintf (stderr, "bad -format argument: %s\n", argv[0]);
	exit (-1);
      }
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-epsilon")) {
      serverargv[serverargc++] = argv[0];
      EPSILON = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-poorhost")) {
      serverargv[serverargc++] = argv[0];
      POORHOST = argv[1];
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-fillins")) {
      FILLINS = atoi (argv[1]);
      argv++;
      argc--;
    } else if (streq (argv[0], "-randomizations")) {
      RANDOMIZATIONS = atoi (argv[1]);
      argv++;
      argc--;
    } else if (streq (argv[0], "-fillmissing")) {
      FILLMISSING = atoi (argv[1]);
      argv++;
      argc--;
    } else if (streq (argv[0], "-maxtrees")) {
      serverargv[serverargc++] = argv[0];
      FINAL_MAXTREES = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-initialmaxtrees")) {
      serverargv[serverargc++] = argv[0];
      INITIAL_MAXTREES = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-maxworkers")) {
      serverargv[serverargc++] = argv[0];
      MAXWORKERS = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-minworkers")) {
      serverargv[serverargc++] = argv[0];
      MINWORKERS = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-loopchunk")) {
      serverargv[serverargc++] = argv[0];
      LOOPCHUNK = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else if (streq (argv[0], "-doublebuffer")) {
      serverargv[serverargc++] = argv[0];
      DOUBLEBUFFER = atoi (argv[1]);
      argv++;
      argc--;
      serverargv[serverargc++] = argv[0];
    } else {
      fprintf (stderr, "bad argument: %s\n", argv[0]);
      exit (-1);
    }
    argv++;
  }
  serverargv[serverargc] = NULL;
   if (SERVERMODE isnt MASTER) LOOPCHUNK = 1; */
/*  fprintf (blather, "Verbose %s.\n", Setting (VERBOSE));
  fprintf (blather, "Statistics printing %s.\n", Setting (STATS));
  fprintf (blather, "Node backup %s.\n", Setting (BACKUP));
  fprintf (blather, "Costly tree pruning %s.\n", Setting (PRUNE));
  fprintf (blather, "Local rearrangement %s.\n", Setting (LOCALREARRANGEMENT));
  fprintf (blather, "Steepest descent %s.\n", Setting (STEEPEST));
  fprintf (blather, "Pablo's hack %s.\n", Setting (PABLO));
  fprintf (blather, "Pablo's hack %s for initial tree construction.\n", Setting (INITIALPABLO));
  fprintf (blather, "SPR %s.\n", Setting (doSPR));
  fprintf (blather, "TBR %s.\n", Setting (doTBR));
  fprintf (blather, "Tree printing %s.\n", Setting (TREE));
  fprintf (blather, "Indented tree printing %s.\n", Setting (INDENT));
  fprintf (blather, "Consensus tree printing %s.\n", Setting (CONSENSUS));
  fprintf (blather, "Consensus bitmap printing %s.\n", Setting (CONSENSUSBITS));
  fprintf (blather, "Missing data assumption statistics %s.\n", Setting (MISSINGSTATS));
  fprintf (blather, "Missing data assumption report %s.\n", Setting (MISSINGREPORT));
  fprintf (blather, "Missing data assumption consensus report %s.\n", Setting (MISSINGCONSENSUS));
  fprintf (blather, "Missing data consensus dataset output %s.\n", Setting (MISSINGCONSENSUSDATA));
  fprintf (blather, "Output tags %s.\n", Setting (TAGS));
  fprintf (blather, "epsilon = %d.\n", EPSILON);*/
  if (INITIAL_MAXTREES is 0)
    INITIAL_MAXTREES = FINAL_MAXTREES;
/*  fprintf (blather, "initialmaxtrees = %d.\n", INITIAL_MAXTREES);
  fprintf (blather, "maxtrees = %d.\n", FINAL_MAXTREES);
  fprintf (blather, "randomizations = %d.\n", RANDOMIZATIONS);
  fprintf (blather, "server mode = %s.\n",
	   SERVERMODE is STANDALONE? "standalone":
	   SERVERMODE is MASTER? "master":
	   SERVERMODE is SERVER? "server":
	   "a bogus value!!!");
  fprintf (blather, "poorhost = %s\n", POORHOST? POORHOST: "NULL");
  fprintf (blather, "minworkers = %d\n", MINWORKERS);
  fprintf (blather, "maxworkers = %d\n", MAXWORKERS);
  fprintf (blather, "loopchunk = %d\n", LOOPCHUNK);
  fprintf (blather, "doublebuffer = %d\n", DOUBLEBUFFER);*/
}

#ifdef _CLYDE_
void ServerLoop () {
  BufferedTreeT *bt = NOBUFFEREDTREE;
  TreeBufferT *tb = NOTREEBUFFER;
  int tag = READDATA;
  int bytes, tid;

  while (tag isnt FINISH) {
    assert (! pvm_bufinfo (pvm_recv (-1, -1), &bytes, &tag, &tid));
    switch (tag) {
    case READDATA: {
      char fname[256];
      int fd;
      FILE *f;
      sprintf (fname, "/tmp/input%d", (int)getpid ());
      fd = open (fname, O_CREAT | O_WRONLY, 0644);
      assert (fd >= 0);
      UnpackFile (fd);
      close (fd);
      f = fopen (fname, "r");
      /*ReadData (f);*/
      PhylogenyAllocations ();
      fclose (f);
      unlink (fname);
      bt = NewBufferedTree ();
      tb = NewTreeBuffer ();
      break;
    }
    case ADDITIONLOOP:
    case LOCALREARRANGEMENTLOOP:
    case SPRLOOP:
    case TBRLOOP: {
      int mincost = UnpackNextInt ();
      int bufferedtrees = UnpackNextInt ();
      int treepresent = UnpackNextInt ();
      if (treepresent) {
	UnpackBufferedTree (bt, 0);
	LoadTree (bt, 0);
      }
      ResetTreeBuffer (tb);
      tb->min_tree_cost = mincost;
      tb->n_buffered_trees = bufferedtrees;
      switch (tag) {
      case ADDITIONLOOP: {
	int firsttaxon = UnpackNextInt(), lasttaxon = UnpackNextInt(),
	OTU_index = UnpackNextInt(), HTU_index = UnpackNextInt();
	AdditionLoop (tb, firsttaxon, lasttaxon, OTU_index, HTU_index);
	break;
      }
      case LOCALREARRANGEMENTLOOP: {
	int firstplacedtaxon = UnpackNextInt (), lastplacedtaxon = UnpackNextInt (), oldprune = PRUNE;
	PRUNE = 0;
	LocalRearrangementLoop (tb, firstplacedtaxon, lastplacedtaxon);
	PRUNE = oldprune;
	break;
      }
      case SPRLOOP: {
	int prunedtaxon = UnpackNextInt (), lastprunedtaxon = UnpackNextInt ();
	int oldcost = UnpackNextInt ();
	SPRLoop (tb, prunedtaxon, lastprunedtaxon, oldcost);
	break;
      }
      case TBRLOOP: {
	int firstprunedtaxon = UnpackNextInt (), lastprunedtaxon = UnpackNextInt ();
	int oldcost = UnpackNextInt ();
	TBRLoop (tb, firstprunedtaxon, lastprunedtaxon, oldcost);
	break;
      }
      default:
	assert (0);
	break;
      }
      break;
    case REPORTSTATS:
      pvm_initsend (PvmDataDefault);
      PackInt (NTrees);
      Optimizations += NOptimizations;
      NOptimizations = 0;
      PackDouble (Optimizations);
      PackInt (NServerDiscards);
      assert (! pvm_send (pvm_parent (), STATSREPORT));
      break;
    case NEWMINTREECOST:
      if (UnpackNextInt () is buffertick) {
	int i = UnpackNextInt ();
	if (i < tb->min_tree_cost) {
	  tb->min_tree_cost = i;
	}
	tb->n_buffered_trees = UnpackNextInt ();
      }
      break;
    case BUFFERTICK:
      buffertick++;
      break;
    case FINISH:
      break;
    }
    default:
      assert (0);
      break;
    }
  }
}
#endif

#define PVMDEBUG 0

#ifdef _CLYDE_
void SpawnServers () {
  int desired, i, host;
  char buf[99], *myname;
  assert (! pvm_config (&nhost, &narch, &hostp));
  desired = nhost;
  if (nhost > 1) {
    assert (! gethostname (buf, sizeof (buf)));
    myname = buf;
    desired--;
  } else
    myname = NULL;
  if (POORHOST isnt NULL)
    desired--;
  if (MAXWORKERS > 0)
    desired = Min (MAXWORKERS, desired);
  if (MINWORKERS > 0)
    desired = Max (MINWORKERS, desired);
  host = 0;
  fprintf (blather, "Using master and %d workers: ", desired);
  LoopBelow (i, desired) {
    while ((POORHOST && streq (hostp[host].hi_name, POORHOST)) ||
	   (myname && !strncmp (hostp[host].hi_name, myname, strlen (myname)))) {
      host++;
      assert (host < nhost);
    }
    assert (host < nhost);
    n_workers = pvm_spawn ("yapp",
			   serverargv,
			   PvmTaskHost | PVMDEBUG,
			   hostp[host].hi_name,
			   1,
			   tids + i);
    assert (n_workers is 1);
    task[i].tid = tids[i];
    task[i].treeid = 0;
    task[i].ntasks = 0;
    fprintf (blather, "%s ", hostp[host].hi_name);
    host++;
  }
  fprintf (blather, "\n");
  n_workers = desired;
  max_tasks_available = n_workers * DOUBLEBUFFER;
  tasks_available = max_tasks_available;
}
#endif

unsigned int **MissingConsensus = NULL;

void TabulateMissingConsensus (int i, int j, int isfirst) {
  if (! MissingConsensus)
    NewMatrix (MissingConsensus, NOTUs, NCharacters, unsigned int);
  if (isfirst)
    MissingConsensus[i][j] = Node[i].final_character[j].base;
  else
    MissingConsensus[i][j] &= Node[i].final_character[j].base;
}

void PrintMissingConsensusDataset (FILE *f) {
  int i, j;
  char *tag = "AssumptionConsensusData";
  fprintf (f, "%s %d %d\n", tag, NOTUs, NCharacters);
  LoopBelow (i, NOTUs) {
    fprintf (f, "%s %s\n%s ", tag, Node[i].name, tag);
    LoopBelow (j, NCharacters)
      fprintf (f, "%d ", MissingConsensus[i][j]);
    fprintf (f, "\n");
  }
}

int daves_cladogram(c_taxa,values,sequences,nbases,cc_weights)
parameters *values;
int nbases;
look_up_thang **sequences;
int *cc_weights;
{
  clock_t starttime, endtime, elapsedtime;
  double oct;
  int cost, tree, i, j, total_candidates;
  TreeBufferT *tb;

  values_here=values;
  NOTUs=c_taxa;
  NCharacters=nbases;

  time (&starttime);
  ProcessCommandLineArguments ();
#ifdef _CLYDE_
  if (SERVERMODE isnt STANDALONE)
    pvm_mytid ();
#endif
  if (SERVERMODE is SERVER) {
#ifdef _CLYDE_
    BUFFERTREE = TransmitTreeToParent;
    LOOPDONE = TellParentDone;
    UPDATEMINTREECOST = UpdateMinTreeCost;
    ServerLoop ();
#endif
  } else {
    if (SERVERMODE is STANDALONE)
        NewReadData (sequences);
      /*ReadData (stdin);*/
    else {
#ifdef _CLYDE_
      char fname[256];
      FILE *f;
      int fd;
      SpawnServers ();
      pvm_initsend (PvmDataDefault);
      PackFile (fileno (stdin));
      assert (! pvm_send (pvm_mytid (), READDATA));
      LoopBelow (i, n_workers)
	assert (! pvm_send (task[i].tid, READDATA));
      pvm_recv (-1, READDATA);
      sprintf (fname, "/tmp/input%d", (int)getpid ());
      fd = open (fname, O_CREAT | O_WRONLY, 0644);
      assert (fd >= 0);
      UnpackFile (fd);
      close (fd);
      f = fopen (fname, "r");
      ReadData (f);
      fclose (f);
      unlink (fname);
      ADDITIONLOOPFUN = AdditionToServer;
      TBRLOOPFUN = TBRToServer;
      SPRLOOPFUN = SPRToServer;
      LOCALREARRANGEMENTLOOPFUN = LocalRearrangementToServer;
      WAITFORSERVERS = WaitForServers;
#endif
    }
    PhylogenyAllocations ();
    tb = Phylogeny ();
    cost = tb->min_tree_cost;
    time (&endtime);
    if (SERVERMODE is MASTER) {
#ifdef _CLYDE_
      pvm_initsend (PvmDataDefault);
      LoopBelow (i, n_workers){
	pvm_send (task[i].tid, REPORTSTATS);
	pvm_recv (-1, STATSREPORT);
	NTrees += UnpackNextInt ();
	Optimizations += UnpackNextDouble ();
	NServerDiscards += UnpackNextInt ();
      }
#endif
    }
    NCandidateTrees += NServerDiscards;
    NCandidatesDiscarded += NServerDiscards;
    Optimizations += NOptimizations;
    NOptimizations = 0;
    if (NCharacters && NTrees)
      oct = Optimizations / ((double) NCharacters) / ((double) NTrees);
    else
      oct = 0;
    elapsedtime = endtime - starttime;
    /*fprintf (blather, "Examined %d trees and optimized %.0f characters in %ld seconds\n",
	     NTrees, Optimizations, elapsedtime);*/
    if (elapsedtime > 0) {
     /* fprintf (blather, "%ld trees/second, %.0f optimizations/second\n",
	       NTrees / elapsedtime, Optimizations / elapsedtime);*/
    }
   /* fprintf (blather, "%f optimizations/character/tree\n", oct);*/
    total_candidates = NCandidateTrees + NTreesReceived;
    if (total_candidates > 0)
     /* fprintf (blather, "%d trees considered for buffering, %d rejected (%d%%), %d discarded (%d%%)\n",
	       total_candidates,
	       NCandidatesRejected, NCandidatesRejected * 100 / total_candidates,
	       NCandidatesDiscarded, NCandidatesDiscarded * 100 / total_candidates
	       );*/

    LoopBelow (tree, tb->n_buffered_trees) {
      LoadTree (tb->buffered_tree[tree], 1);
      if (TAGS && tree is 0) {
	fprintf (stdout, "Multiplicity %d\n", tb->n_buffered_trees);
      }
      if (TREE) {
	fprintf (stdout, TAGS? "Tree %d ": "Tree%d = ", tree);
	PrintTree (stdout);
	fprintf (stdout, TAGS? "\n": "[%d];\n", tb->buffered_tree[tree]->cost);
	if (TAGS)
	  fprintf (stdout, "Length %d\n", tb->buffered_tree[tree]->cost);
      }
      if (INDENT) {
	if (tree > 0) fprintf (stdout, "\n");
	IndentTree (stdout);
      }
      if (CONSENSUS || CONSENSUSBITS || CONSENSUSTABLE) {
	if (tree is 0)
	  InitConsensus ();
	else
	  Consense ();
      }
      if (MISSINGREPORT || MISSINGSTATS) {
	int i, j, nmissing = 0, nchanged = 0, nsingletoned = 0;
	LoopBelow (i, NOTUs)
	  LoopBelow (j, NCharacters) {
	    if (MissingDataCode (RawData[i][j])) {
	      nmissing++;
	    }
	    if (Node[i].character[j].base isnt Node[i].final_character[j].base) {
	      if (MISSINGREPORT)
		fprintf (stdout, "Assuming %s %d %d %d\n", Node[i].name, i, j, Node[i].final_character[j]);
	      nchanged++;
	      if (CountBits ( Node[i].final_character[j].base) is 1)
		nsingletoned++;
	    }
	  }
	if (MISSINGSTATS && nmissing) {
	  fprintf (stdout, "MissingStats %d of %d (%f%%) are missing.\n",
		   nmissing, NOTUs * NCharacters, (nmissing * 100.0) / (NOTUs * NCharacters));
	  fprintf (stdout, "MissingStats %d of %d changed (%f%%).\n",
		   nchanged, nmissing, (nchanged * 100.0)/nmissing);
	  fprintf (stdout, "MissingStats %d of %d changed to singletons (%f%%).\n",
		   nsingletoned, nmissing, (nsingletoned * 100.0)/nmissing);
	}
      }
      if (MISSINGCONSENSUS || MISSINGCONSENSUSDATA)
	LoopBelow (i, NOTUs)
	  LoopBelow (j, NCharacters)
	    TabulateMissingConsensus (i, j, tree is 0);
    }
    if (MISSINGCONSENSUS)
      LoopBelow (i, NOTUs)
	LoopBelow (j, NCharacters)
	  if (MissingDataCode (RawData[i][j]))
	    fprintf (stdout, "AssumptionConsensus %s %d %d %d\n", Node[i].name, i, j, MissingConsensus[i][j]);
    if (MISSINGCONSENSUSDATA)
      PrintMissingConsensusDataset (stdout);
    if (CONSENSUSBITS) {
      PrintConsensusBits (stdout);
    }
    if (CONSENSUSTABLE) {
      PrintConsensusTable (stdout);
    }
    if (CONSENSUS) {
      fprintf (stdout, TAGS? "Consensus ": "Consensus_of_%d = ", tb->n_buffered_trees);
      PrintConsensus (stdout);
      fprintf (stdout, TAGS? "\n": ";\n");
    }
    if (STATS) {
      fprintf (stdout,
	       "%d OTUs %d characters %d trees %.0f optimizations %ld seconds %f o/c/t cost %d multiplicity %d\n",
	       NOTUs, NCharacters,
	       NTrees, Optimizations, endtime - starttime, oct, cost, tb->n_buffered_trees);
    }
    LoopBelow (i, FILLINS) {
      fprintf (stderr, " %d", i);
      fprintf (stdout, "yapp -notree -noconsensus -consensusbits -format bitmap << STOPSTOPSTOP | tally-bits consensus.tally \n");
      WriteData (stdout);
      fprintf (stdout, "STOPSTOPSTOP\n");
      fprintf (stdout, "echo %d\n", i);
    }
    if (FILLINS)
      fprintf (stderr, "\n");
    if (SERVERMODE is MASTER) {
#ifdef _CLYDE_
      pvm_initsend (PvmDataDefault);
      LoopBelow (i, n_workers)
	assert (! pvm_send (task[i].tid, FINISH));
#endif
    }
    /*fprintf (blather, "%d trees received, %d without support data, %d with.\n",
	     RECEIVEDWITH + RECEIVEDWITHOUT, RECEIVEDWITHOUT, RECEIVEDWITH);
    fprintf (blather, "%d trees sent, %d without support data, %d with.\n",
	     SENTWITH + SENTWITHOUT, SENTWITHOUT, SENTWITH);
    fprintf (blather, "%d server tasks dispatched, %d without trees, %d with.\n",
	     SERVERTASKSWITH + SERVERTASKSWITHOUT, SERVERTASKSWITHOUT, SERVERTASKSWITH);*/

  }
  DeallocateEverything ();
  /*exit (0);*/
  return cost;
}
