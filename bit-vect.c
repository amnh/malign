#include "align3.h"

int BitVectorWords, BitVectorSize;

void InitBitVectorSize (unsigned n) {
  BitVectorSize = n;
  BitVectorWords = (n + BV_INDEX_MASK) >> BV_INDEX_SHIFT;
}

BitVectorT NewBitVector () {
  BitVectorT x;
  NewArray (x, BitVectorWords, unsigned int);
  memset (x, 0, BitVectorWords * sizeof (unsigned int));
  return (x);
}

int QsortCompareBitVectors (const void *xp, const void *yp) {
  BitVectorT x = *(BitVectorT *)xp, y = *(BitVectorT *)yp;
  int i;
  for (i = BitVectorWords - 1; i >= 0; i--)
    if (x[i] < y[i]) return (-1);
  else
    if (x[i] > y[i]) return (1);
  return (0);
}

int CompareBitVectors (BitVectorT x, BitVectorT y) {
  int i;
  for (i = BitVectorWords - 1; i >= 0; i--)
    if (x[i] < y[i]) return (-1);
  else
    if (x[i] > y[i]) return (1);
  return (0);
}

void SetBit(BitVectorT bv, unsigned int n) {
  bv[BitVectorWordIndex (n)] |= (1 << BitVectorBitIndex (n));
}

int GetBit(BitVectorT bv, unsigned int n) {
  return ((bv[BitVectorWordIndex (n)] & (1 << BitVectorBitIndex (n)))
	  >> BitVectorBitIndex (n));
}

int BitVectorSubsetEq (BitVectorT x, BitVectorT y) {
  int i;
  LoopBelow (i, BitVectorWords)
    if ((x[i] & ~y[i]) isnt 0)
      return (0);
  return (1);
}

void fprintfBitVector (FILE *f, BitVectorT x) {
  int i;
  for (i = BitVectorSize - 1; i >= 0; i--)
    fprintf (f, "%d", GetBit (x, i));
}

unsigned int BitVectorHash (BitVectorT x) {
  unsigned int h = 0;
  int i;

  LoopBelow (i, BitVectorWords)
    h += x[i];
  return (h);
}

void ClearBitVector (BitVectorT x) {
  int i;
  LoopBelow (i, BitVectorWords)
    x[i] = 0;
}

void CopyBitVector (BitVectorT to, BitVectorT from) {
  int i;
  LoopBelow (i, BitVectorWords)
    to[i] = from[i];
}

void OrBitVectors (BitVectorT to, BitVectorT from1, BitVectorT from2) {
  int i;
  LoopBelow (i, BitVectorWords)
    to[i] = from1[i] | from2[i];
}

int fscanfBitVector (FILE *f, BitVectorT x) {
  int i, bit;
  ClearBitVector (x);
  for (i = BitVectorSize - 1; i >= 0; i--) {
    if (fscanf (f, "%1d", &bit) is 1) {
      if (bit) SetBit (x, i);
    } else {
      assert (i is BitVectorSize - 1);
      return (0);
    }
  }
  return (1);
}
