#define UNIONBIT (1 << (WORDSIZE-1))
#define SETBITS (~ UNIONBIT)
#define MAXCHARACTERPOSITIONS (WORDSIZE-1)
#define INFINITY (1 <<  20)

typedef struct t {
  struct t *parent, *left, *right;
  int taxon;
  char *name;
  int localcost, totalcost;
  int mark;
  int dirty;
  struct t *backupnode;
  look_up_thang  *final_character;
  look_up_thang  *partial_final_character;
  look_up_thang  *character;
  BitVectorT taxa;
} NodeT;
#define NONODE ((NodeT *)0)

typedef enum {NOBACKUP, TOBACKUP, FROMBACKUP} BackupT;

typedef struct {
  int parent_index: 31;
  int on_left: 1;
} BufferedTreeNodeT;

typedef struct {
  int n_placed_taxa;
  BufferedTreeNodeT *nodes;
  int *placed_taxa_indices;
  int n_supported_clades;
  BitVectorT *supported_clades;
  unsigned int supported_clades_hash;
  int cost;
  int generation;
} BufferedTreeT;

#define NOBUFFEREDTREE (BufferedTreeT *)NULL

typedef struct {
  int n_buffered_trees;
  int max_buffered_trees;
  BufferedTreeT **buffered_tree;
  int min_tree_cost;
  int current_generation;
  int compressed;
} TreeBufferT;

#define NOTREEBUFFER (TreeBufferT *)NULL

typedef struct t2 {
  int n_kids;
  struct t2 **kids;
  int has_parent;
  int clade;
} ConsensusNodeT;

typedef void (*ComputeSupportT) (NodeT *root);

typedef enum {READDATA, ADDITIONLOOP, LOCALREARRANGEMENTLOOP, SPRLOOP, TBRLOOP, FINISH, REPORTSTATS,
		NEWMINTREECOST, BUFFERTICK,
		RETURNTREE, LOOPISDONE, STATSREPORT} MSGTagT;

