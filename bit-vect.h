#define BV_INDEX_SHIFT 5
#define BV_INDEX_MASK 31
#define BitVectorWordIndex(i) ((i)>>BV_INDEX_SHIFT)
#define BitVectorBitIndex(i) ((i)&BV_INDEX_MASK)

typedef unsigned int *BitVectorT;
extern int BitVectorWords, BitVectorSize;

void InitBitVectorSize (unsigned n);
BitVectorT NewBitVector ();
int QsortCompareBitVectors (const void *xp, const void *yp);
int CompareBitVectors (BitVectorT x, BitVectorT y);
void SetBit(BitVectorT bv, unsigned int n);
int GetBit(BitVectorT bv, unsigned int n);
int BitVectorSubsetEq (BitVectorT x, BitVectorT y);
void fprintfBitVector (FILE *f, BitVectorT x);
unsigned int BitVectorHash (BitVectorT x);
void ClearBitVector (BitVectorT x);
void CopyBitVector (BitVectorT to, BitVectorT from);
void OrBitVectors (BitVectorT to, BitVectorT from1, BitVectorT from2);
int fscanfBitVector (FILE *f, BitVectorT x);
