#include "align3.h"
#include <sys/types.h>
#include <unistd.h>

int RECEIVEDWITH = 0, RECEIVEDWITHOUT = 0, SENTWITH = 0, SENTWITHOUT = 0;

#ifdef _CLYDE_
#define PACKFILEBUFFERSIZE 1024

void PackFile (int fd) {
  int n;
  char buf[PACKFILEBUFFERSIZE];
  while ((n = read (fd, buf, sizeof (buf))) > 0) {
    assert (! pvm_pkint (&n, 1, 1));
    assert (! pvm_pkbyte (buf, n, 1));
  }
  n = 0;
  assert (! pvm_pkint (&n, 1, 1));
}

int UnpackNextInt () {
  int i;
  assert (! pvm_upkint (&i, 1, 1));
  return (i);
}

void PackInt (int i) {
  assert (! pvm_pkint (&i, 1, 1));
}

double UnpackNextDouble () {
  double i;
  assert (! pvm_upkdouble (&i, 1, 1));
  return (i);
}

void PackDouble (double d) {
  assert (! pvm_pkdouble (&d, 1, 1));
}

void UnpackFile (int fd) {
  int n;
  char buf[PACKFILEBUFFERSIZE];
  while ((n = UnpackNextInt ()) > 0) {
    assert (n <= sizeof (buf));
    assert (! pvm_upkbyte (buf, n, 1));
    assert (write (fd, buf, n) is n);
  }
}

void PackBufferedTree (BufferedTreeT *bt, int support_stuff) {
  int i;
/*   { */
/*     int i; */
/*     fprintf (stderr, "packing tree:\n"); */
/*     LoopBelow (i, bt->n_placed_taxa) */
/*       fprintf (stderr, "%d %d %c\n", */
/* 	       bt->placed_taxa_indices[i], */
/* 	       bt->nodes[i].parent_index, */
/* 	       bt->nodes[i].on_left? 'l': 'r'); */
/*   } */
  assert (! pvm_pkint (&bt->n_placed_taxa, 1, 1));
  assert (! pvm_pkint ((int *)bt->nodes, bt->n_placed_taxa, 1));
  assert (! pvm_pkint (bt->placed_taxa_indices, bt->n_placed_taxa, 1));
  if (support_stuff) {
    SENTWITH++;
    assert (! pvm_pkint (&bt->n_supported_clades, 1, 1));
    LoopBelow (i, bt->n_supported_clades)
      assert (! pvm_pkint (bt->supported_clades[i], BitVectorWords, 1));
    assert (! pvm_pkint (&bt->supported_clades_hash, 1, 1));
  } else {
    SENTWITHOUT++;
    bt->n_supported_clades = bt->supported_clades_hash = 0;
  }
  assert (! pvm_pkint (&bt->cost, 1, 1));
  assert (! pvm_pkint (&bt->generation, 1, 1));
}

void UnpackBufferedTree (BufferedTreeT *bt, int support_stuff) {
  int i;
  assert (! pvm_upkint (&bt->n_placed_taxa, 1, 1));
  assert (! pvm_upkint ((int *)bt->nodes, bt->n_placed_taxa, 1));
  assert (! pvm_upkint (bt->placed_taxa_indices, bt->n_placed_taxa, 1));
  if (support_stuff) {
    RECEIVEDWITH++;
    assert (! pvm_upkint (&bt->n_supported_clades, 1, 1));
    LoopBelow (i, bt->n_supported_clades)
      assert (! pvm_upkint (bt->supported_clades[i], BitVectorWords, 1));
    assert (! pvm_upkint (&bt->supported_clades_hash, 1, 1));
  } else {
    RECEIVEDWITHOUT++;
    bt->n_supported_clades = bt->supported_clades_hash = 0;
  }
  assert (! pvm_upkint (&bt->cost, 1, 1));
  assert (! pvm_upkint (&bt->generation, 1, 1));
/*   { */
/*     int i; */
/*     fprintf (stderr, "unpacking tree:\n"); */
/*     LoopBelow (i, bt->n_placed_taxa) */
/*       fprintf (stderr, "%d %d %c\n", */
/* 	       bt->placed_taxa_indices[i], */
/* 	       bt->nodes[i].parent_index, */
/* 	       bt->nodes[i].on_left? 'l': 'r'); */
/*   } */
}

#endif
