#define WORDSIZE 32
#define is ==
#define isnt !=
#define LoopBelow(var,limit) for (var = 0; var < limit; var++)
#define LoopCount(var,count) for (var = count; var > 0; var--)
#define LoopInterval(var,start,limit) for (var = start; var <= limit; var++)
#define Min(x,y) (x<y? x: y)
#define Max(x,y) (x>y? x: y)
#define Infinity 0x7FFFFFF
#define New(n,what) ((what *)Allocate_D((n)*sizeof(what)))
#define NewArray(where,n,what) where=((what *)Allocate_D((n)*sizeof(what)))
#define NewMatrix(where,rows,columns,what) {int __i;where=New(rows,what *);LoopBelow(__i,rows)where[__i]=New(columns,what);}
#define scanstring(buf) assert (1 == scanf ("%s", buf))
#define fscanstring(f,buf) assert (1 == fscanf (f, "%s", buf))

double UniformRand (void);
int BitRand (void);
void *Allocate_D (unsigned int);
void *Reallocate_D (void *p, unsigned int);
int streq (char *s1, char *s2);
int round (double x);
char *CopyString (char *s);
void fail (char *s);
int CountBits (int x);
int MissingDataCode (char c);
void DeallocateEverything ();
