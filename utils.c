#include "align3.h"

double UniformRand (void)
{
  double r = (rand () / (double) 0x7FFF);
  assert (0 <= r && r <= 1);
  return (r);
}

int BitRand (void)
{
  return (rand () % 2);
}

void **Allocations = (void **)NULL;
int NAllocations = 0;
int MaxAllocations = 0;

void *Allocate_D (unsigned u)
{
  void *p;
  if (NAllocations >= MaxAllocations) {
    MaxAllocations = Max (2048, 2 * MaxAllocations);
    Allocations = realloc (Allocations, MaxAllocations * sizeof (void *));
    assert (Allocations);
    assert (NAllocations < MaxAllocations);
  }
  p = (void *)malloc (u);
  assert (p != (void *)NULL);
  Allocations[NAllocations++] = p;
  return (p);
}

void *Reallocate_D (void *p, unsigned n)
{
  int i;
  for (i = NAllocations - 1; i >= 0; i--)
    if (p is Allocations[i])
      break;
  assert (i >= 0);
  p = realloc (p, n);
  assert (p);
  Allocations[i] = p;
  return (p);
}

void DeallocateEverything () {
  int i;
  LoopBelow (i, NAllocations)
    free (Allocations[i]);
  free (Allocations);
  Allocations = (void **) NULL;
  MaxAllocations = NAllocations = 0;
}

int streq (char *s1, char *s2)
{
  return (0 == strcmp (s1, s2));
}

int round (double x)
{
  return ((int) (0.5+x));
}

char *CopyString (char *s)
{
  char *t = Allocate_D (1 + strlen (s));
  strcpy (t, s);
  return (t);
}

void fail (char *s)
{
  fprintf (stderr, "Error: %s\n", s);
  strcpy (NULL, "Kaboom!");
  exit (-1);
}

int CountBits (int x) {
  int bits = 0, i;
  LoopBelow (i, WORDSIZE)
    if (x & (1 << i))
      bits++;
  return (bits);
}

int MissingDataCode (c)
char c;
{
  return (c is '-' ||
	  c is '.' ||
	  c is '?');
}

