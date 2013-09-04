extern void PackFile (int fd);
extern int UnpackNextInt ();
extern void PackInt ();
extern void UnpackFile (int fd);
extern void PackDouble (double d);
extern double UnpackNextDouble ();
extern void PackBufferedTree (BufferedTreeT *b, int support_stuff);
extern void UnpackBufferedTree (BufferedTreeT *b, int support_stuff);
extern int RECEIVEDWITH, RECEIVEDWITHOUT, SENTWITH, SENTWITHOUT;
