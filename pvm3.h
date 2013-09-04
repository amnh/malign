
/*
 *           PVM 3.2:  Parallel Virtual Machine System 3.2
 *               University of Tennessee, Knoxville TN.
 *           Oak Ridge National Laboratory, Oak Ridge TN.
 *                   Emory University, Atlanta GA.
 *      Authors:  A. L. Beguelin, J. J. Dongarra, G. A. Geist,
 *    W. C. Jiang, R. J. Manchek, B. K. Moore, and V. S. Sunderam
 *                   (C) 1992 All Rights Reserved
 *
 *                              NOTICE
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted
 * provided that the above copyright notice appear in all copies and
 * that both the copyright notice and this permission notice appear in
 * supporting documentation.
 *
 * Neither the Institutions (Emory University, Oak Ridge National
 * Laboratory, and University of Tennessee) nor the Authors make any
 * representations about the suitability of this software for any
 * purpose.  This software is provided ``as is'' without express or
 * implied warranty.
 *
 * PVM 3.2 was funded in part by the U.S. Department of Energy, the
 * National Science Foundation and the State of Tennessee.
 */

/*
 *	pvm3.h
 *
 *	Libpvm3 includes.
 *
$Log$
 */

#ifndef	_PVM3_H_

#define	_PVM3_H_

/*
*	Data packing styles for pvm_initsend()
*/

#define	PvmDataDefault	0
#define	PvmDataRaw		1
#define	PvmDataInPlace	2
#define	PvmDataFoo		3

/*
*	pvm_spawn options
*/

#define	PvmTaskDefault	0
#define	PvmTaskHost		1	/* specify host */
#define	PvmTaskArch		2	/* specify architecture */
#define	PvmTaskDebug	4	/* start task in debugger */
#define	PvmTaskTrace	8	/* process generates trace data */
/* for MPP ports */
#define	PvmMppFront		16	/* spawn task on service node */
#define	PvmHostCompl	32	/* complement host set */

/*
*	pvm_notify types
*/

#define	PvmTaskExit		1	/* on task exit */
#define	PvmHostDelete	2	/* on host fail/delete */
#define	PvmHostAdd		3	/* on host startup */

/*
*	for pvm_setopt and pvm_getopt
*/

#define	PvmRoute		1	/* routing policy */
#define		PvmDontRoute	1	/* don't allow direct task-task links */
#define		PvmAllowDirect	2	/* allow direct links, but don't request */
#define		PvmRouteDirect	3	/* request direct links */
#define	PvmDebugMask	2	/* debugmask */
#define	PvmAutoErr		3	/* auto error reporting */
#define	PvmOutputTid	4	/* stdout device for children */
#define	PvmOutputCode	5
#define	PvmTraceTid		6	/* trace device for children */
#define	PvmTraceCode	7
#define	PvmFragSize		8	/* message fragment size */

/*
*	Libpvm error codes
*/

#define	PvmOk			0	/* okay */
							/* reserve -1 */
#define	PvmBadParam		-2	/* bad parameter (neg msg id, etc) */
#define	PvmMismatch		-3	/* barrier count mismatch */
#define	PvmNoData		-5	/* read past end of buffer */
#define	PvmNoHost		-6	/* no such host */
#define	PvmNoFile		-7	/* no such executable */
#define	PvmNoMem		-10	/* can't get memory */
#define	PvmBadMsg		-12	/* received msg can't be decoded */
#define	PvmSysErr		-14	/* can't contact our pvmd/some system error */
#define	PvmNoBuf		-15	/* no current buffer */
#define	PvmNoSuchBuf	-16	/* bad message id */
#define	PvmNullGroup	-17	/* null group name is illegal */
#define	PvmDupGroup		-18	/* already in group */
#define	PvmNoGroup		-19	/* no group with name */
#define	PvmNotInGroup	-20	/* task not in group */
#define	PvmNoInst		-21	/* no such instance in group */
#define	PvmHostFail		-22	/* host failed */
#define	PvmNoParent		-23	/* no parent task */
#define	PvmNotImpl		-24	/* function not implemented */
#define	PvmDSysErr		-25	/* pvmd system error */
#define	PvmBadVersion	-26	/* pvmd-pvmd protocol version mismatch */
#define	PvmOutOfRes		-27	/* out of resources */
#define	PvmDupHost		-28	/* host already configured */
#define	PvmCantStart	-29	/* failed to exec new slave pvmd */
#define	PvmAlready		-30	/* already doing operation */
#define	PvmNoTask		-31	/* no such task */
#define	PvmNoEntry		-32	/* no such name, index pair */
#define	PvmDupEntry		-33	/* name, index pair already exists */

/*
*	returned by pvm_config()
*/

struct hostinfo {
	int hi_tid;			/* pvmd tid */
	char *hi_name;		/* host name */
	char *hi_arch;		/* host arch */
	int hi_speed;		/* cpu relative speed */
};

/*
*	returned by pvm_tasks()
*/

struct taskinfo {
	int ti_tid;				/* task id */
	int ti_ptid;			/* parent tid */
	int ti_host;			/* pvmd tid */
	int ti_flag;			/* status flags */
	char *ti_a_out;			/* a.out name */
};


#ifdef __ProtoGlarp__
#undef __ProtoGlarp__
#endif
#if defined(__STDC__) || defined(__cplusplus)
#define __ProtoGlarp__(x) x
#else
#define __ProtoGlarp__(x) ()
#endif

#ifdef __cplusplus
extern "C" {
#endif

int	pvm_advise		__ProtoGlarp__(( int what ));	/* TB REMOVED */
int	pvm_addhosts	__ProtoGlarp__(( char **names, int count, int *svp ));
int	pvm_barrier		__ProtoGlarp__(( char *group, int count ));
int	pvm_bcast		__ProtoGlarp__(( char *group, int code ));
int	pvm_bufinfo		__ProtoGlarp__(( int mid, int *len, int *code, int *tid ));
int	pvm_config		__ProtoGlarp__(( int *nhostp, int *narchp,
										struct hostinfo **hostp ));
int	pvm_delete		__ProtoGlarp__(( char *name, int req ));
int	pvm_delhosts	__ProtoGlarp__(( char **names, int count, int *svp ));
int	pvm_exit		__ProtoGlarp__(( void ));
int	pvm_freebuf		__ProtoGlarp__(( int mid ));
int	pvm_getfds		__ProtoGlarp__(( int **fds ));
int	pvm_getinst		__ProtoGlarp__(( char *group, int tid ));
int	pvm_getopt		__ProtoGlarp__(( int what ));
int	pvm_getrbuf		__ProtoGlarp__(( void ));
int	pvm_getsbuf		__ProtoGlarp__(( void ));
int	pvm_gettid		__ProtoGlarp__(( char *group, int inst ));
int	pvm_gsize		__ProtoGlarp__(( char *group ));
int	pvm_halt		__ProtoGlarp__(( void ));
int	pvm_initsend	__ProtoGlarp__(( int encod ));
int	pvm_insert		__ProtoGlarp__(( char *name, int req, int data ));
int	pvm_joingroup	__ProtoGlarp__(( char *group ));
int	pvm_kill		__ProtoGlarp__(( int tid ));
int	pvm_lookup		__ProtoGlarp__(( char *name, int req, int *datap ));
int	pvm_lvgroup		__ProtoGlarp__(( char *group ));
int	pvm_mcast		__ProtoGlarp__(( int *tids, int count, int code ));
int	pvm_mkbuf		__ProtoGlarp__(( int encod ));
int	pvm_mstat		__ProtoGlarp__(( char *host ));
int	pvm_mytid		__ProtoGlarp__(( void ));
int	pvm_notify		__ProtoGlarp__(( int what, int code,
										int count, int *vals ));
int	pvm_nrecv		__ProtoGlarp__(( int tid, int code ));
int	pvm_packf		__ProtoGlarp__(( const char *fmt, ... ));
int	pvm_parent		__ProtoGlarp__(( void ));
int	pvm_perror		__ProtoGlarp__(( char *msg ));
int	pvm_pkbyte		__ProtoGlarp__(( char *cp, int cnt, int std ));
int	pvm_pkcplx		__ProtoGlarp__(( float *xp, int cnt, int std ));
int	pvm_pkdcplx		__ProtoGlarp__(( double *zp, int cnt, int std ));
int	pvm_pkdouble	__ProtoGlarp__(( double *dp, int cnt, int std ));
int	pvm_pkfloat		__ProtoGlarp__(( float *fp, int cnt, int std ));
int	pvm_pkint		__ProtoGlarp__(( int *np, int cnt, int std ));
int	pvm_pklong		__ProtoGlarp__(( long *np, int cnt, int std ));
int	pvm_pkshort		__ProtoGlarp__(( short *np, int cnt, int std ));
int	pvm_pkstr		__ProtoGlarp__(( char *cp ));
int	pvm_pkuint		__ProtoGlarp__(( unsigned int *np, int cnt, int std ));
int	pvm_pkulong		__ProtoGlarp__(( unsigned long *np, int cnt, int std ));
int	pvm_pkushort	__ProtoGlarp__(( unsigned short *np, int cnt, int std ));
int	pvm_probe		__ProtoGlarp__(( int tid, int code ));
int	pvm_pstat		__ProtoGlarp__(( int tid ));
int	pvm_recv		__ProtoGlarp__(( int tid, int code ));
int	(*pvm_recvf		__ProtoGlarp__(( int (*newf)(int, int, int) )) )();
int	pvm_send		__ProtoGlarp__(( int tid, int code ));
int	pvm_sendsig		__ProtoGlarp__(( int tid, int signum ));
int	pvm_serror		__ProtoGlarp__(( int how ));	/* TB REMOVED */
int	pvm_setdebug	__ProtoGlarp__(( int mask ));	/* TB REMOVED */
int	pvm_setopt		__ProtoGlarp__(( int what, int val ));
int	pvm_setrbuf		__ProtoGlarp__(( int mid ));
int	pvm_setsbuf		__ProtoGlarp__(( int mid ));
int	pvm_spawn		__ProtoGlarp__(( char *file, char **argv, int flags,
										char *where, int count, int *tids ));
int	pvm_start_pvmd	__ProtoGlarp__(( int argc, char **argv, int block ));
int	pvm_tasks		__ProtoGlarp__(( int where, int *ntaskp,
										struct taskinfo **taskp ));
int	pvm_tickle		__ProtoGlarp__(( int narg, int *argp,
										int *nresp, int *resp ));
int	pvm_tidtohost	__ProtoGlarp__(( int tid ));
int	pvm_unpackf		__ProtoGlarp__(( const char *fmt, ... ));
int	pvm_upkbyte		__ProtoGlarp__(( char *cp, int cnt, int std ));
int	pvm_upkcplx		__ProtoGlarp__(( float *xp, int cnt, int std ));
int	pvm_upkdcplx	__ProtoGlarp__(( double *zp, int cnt, int std ));
int	pvm_upkdouble	__ProtoGlarp__(( double *dp, int cnt, int std ));
int	pvm_upkfloat	__ProtoGlarp__(( float *fp, int cnt, int std ));
int	pvm_upkint		__ProtoGlarp__(( int *np, int cnt, int std ));
int	pvm_upklong		__ProtoGlarp__(( long *np, int cnt, int std ));
int	pvm_upkshort	__ProtoGlarp__(( short *np, int cnt, int std ));
int	pvm_upkstr		__ProtoGlarp__(( char *cp ));
int	pvm_upkuint		__ProtoGlarp__(( unsigned int *np, int cnt, int std ));
int	pvm_upkulong	__ProtoGlarp__(( unsigned long *np, int cnt, int std ));
int	pvm_upkushort	__ProtoGlarp__(( unsigned short *np, int cnt, int std ));
char *pvm_version	__ProtoGlarp__(( void ));

#ifdef __cplusplus
}
#endif

#endif	/*_PVM3_H_*/

