/*
 * CASINO shared memory
 * ====================
 *
 * Preprocessor flags used by this file:
 * -------------------------------------
 * SHM_SYSV      : activates System V shared memory (most widely used).
 * SHM_POSIX     : activates POSIX shared memory.
 * SHM_POSIX_BGQ : activates a specialised version of POSIX shared memory 
 *                 which tries to go around Blue Gene Q bugs (see below).
 * DBG           : activates lots of debug prints, especially if 
 *                 shm_debug=.true. in CASINO module shallocate_smp.f90.
 * 
 * Notes:
 * ------
 * - The preprocessor flags for alloc_shm.c can be set in the machine-specific
 *   *.arch files using make variable CFLAGS_SHM or by passing the
 *   preprocessor flags to the variable CFLAGS_SHM on the make command line.
 * 
 * - If neither SHM_SYSV, SHM_POSIX or SHM_POSIX_BGQ preprocessor flags are set
 *   and USE_SHM=yes the compilation will stop with an error.
 * 
 * - For preprocessor flags needed to generate the C function name expected
 *   by the Fortran compiler see the comments below.
 *
 * - Posix file name randomized to avoid possible name clashes, see
 *   get_randomTag function
 *
 * - Hack for Blue Gene systems: 
 *   Shared memory implementation on BG systems is buggy: 
 *   it does not removed unlinked files, ftruncate produces unexpected 
 *   results, mmap does not use the offset argument, and God knows what else.
 *   
 *   To avoid all this trouble the implementation opens only one shared
 *   memory file, the memory allocation needed by CASINO arrays is
 *   controlled in this file with a linked list of pointers to blocks of
 *   memory. The basic memory model of the algorithm is a stack of memory
 *   blocks which reuses freed blocks if the size of the new memory
 *   request fits into an already-existing block.
 *     
 *   It is not very smart but it should do for the current operation pattern 
 *   used by CASINO in shared memory (i.e. allocate 1-2 files to store data and
 *   used a third as a temporary buffer for IO).
 *
 *   On Blue Genes the size of the shared memory file must be known before the 
 *   run starts and it must be passed to the executable (in MB) using the 
 *   environment variable CASINO_MAXSHM.
 *
 *   In order to avoid segmentation of the shared memory file deshalloc the 
 *   last shalloc array whenever is possible in the application code.
 *
 *  LA 01.2011, 02.2014
 *
 * - Do not use the modern // comment designator; only the old-style one
 *   with asterisks is supported by some currently available C compilers.
 *
 *  MDT, 12.2012
 */

#ifdef SHM_POSIX_BGQ
#define SHM_POSIX
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#ifdef SHM_SYSV
  #include <sys/ipc.h>
  #include <sys/shm.h>
#else
  #ifdef SHM_POSIX
    #include <unistd.h>
    #include <sys/mman.h>
    #include <fcntl.h>
  #endif
#endif
#include <errno.h>
#include <mpi.h>
#ifdef _OPENMP
  #include <omp.h>
#endif

/*
 * The Fortran compiler may expect C functions to have:
 * - no underscores appended to end of name [use -DF90_NO_UNDERSCORE]
 * - one underscore appended to end of name [default]
 * - two underscores appended to end of name [use -DF90_DOUBLE_UNDERSCORE]
 * and:
 * - name in lower case [default]
 * - name in upper case [use -DF90_CAPITALS]
 */
#ifndef ALLOC_SHM
  #ifdef F90_DOUBLE_UNDERSCORE
    #ifdef F90_CAPITALS
      #define ALLOC_SHM ALLOC_SHM__
      #define DEALLOC_SHM DEALLOC_SHM__
      #define GET_SMP_LIST GET_SMP_LIST__
      #define SET_SHM_DEBUG SET_SHM_DEBUG__
      #define CLEAN_SHM CLEAN_SHM__
      #define SET_SHMSIZE SET_SHMSIZE__
      #define GET_SHMSIZE GET_SHMSIZE__
      #define GET_NNPSMP GET_NNPSMP__
    #else
      #define ALLOC_SHM alloc_shm__
      #define DEALLOC_SHM dealloc_shm__
      #define GET_SMP_LIST get_smp_list__
      #define SET_SHM_DEBUG set_shm_debug__
      #define CLEAN_SHM clean_shm__ 
      #define SET_SHMSIZE set_shmsize__
      #define GET_SHMSIZE get_shmsize__
      #define GET_NNPSMP get_nnpsmp__
    #endif
  #else
    #ifdef F90_NO_UNDERSCORE
      #ifdef F90_CAPITALS
        #define ALLOC_SHM ALLOC_SHM
        #define DEALLOC_SHM DEALLOC_SHM
        #define GET_SMP_LIST GET_SMP_LIST
        #define SET_SHM_DEBUG SET_SHM_DEBUG
        #define CLEAN_SHM CLEAN_SHM
        #define SET_SHMSIZE SET_SHMSIZE
        #define GET_SHMSIZE GET_SHMSIZE
        #define GET_NNPSMP GET_NNPSMP
      #else
        #define ALLOC_SHM alloc_shm
        #define DEALLOC_SHM dealloc_shm
        #define GET_SMP_LIST get_smp_list
        #define SET_SHM_DEBUG set_shm_debug
        #define CLEAN_SHM clean_shm
        #define SET_SHMSIZE set_shmsize
        #define GET_SHMSIZE get_shmsize
        #define GET_NNPSMP get_nnpsmp
      #endif
    #else
      #ifdef F90_CAPITALS
        #define ALLOC_SHM ALLOC_SHM_
        #define DEALLOC_SHM DEALLOC_SHM_
        #define GET_SMP_LIST GET_SMP_LIST_
        #define SET_SHM_DEBUG SET_SHM_DEBUG_
        #define CLEAN_SHM CLEAN_SHM_
        #define SET_SHMSIZE SET_SHMSIZE_
        #define GET_SHMSIZE GET_SHMSIZE_
        #define GET_NNPSMP GET_NNPSMP_
      #else
        #define ALLOC_SHM alloc_shm_
        #define DEALLOC_SHM dealloc_shm_
        #define GET_SMP_LIST get_smp_list_
        #define SET_SHM_DEBUG set_shm_debug_
        #define CLEAN_SHM clean_shm_
        #define SET_SHMSIZE set_shmsize_
        #define GET_SHMSIZE get_shmsize_
        #define GET_NNPSMP get_nnpsmp_
      #endif
    #endif
  #endif
#endif

static MPI_Comm node_comm;
static MPI_Group group_world, node_group;

static int shm_debug=0;
static int mype, npes, mast, err;
static const int errlen=256;
#ifdef DBG
  static int am_nodemaster;
#endif
static int numablk=0,nthreads=0;

void SET_SHM_DEBUG(void);
void GET_NNPSMP(int *nnpsmp);
void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes);
void ALLOC_SHM(void **ptr, int64_t *nelem, MPI_Fint *ftype, int *ret);
void DEALLOC_SHM(void **shm, int64_t *nelem, MPI_Fint *ftype);
/* These two functions are necessary only for BGQ hack */
void SET_SHMSIZE(int64_t *s);
void GET_SHMSIZE(int64_t *s);


#ifdef SHM_SYSV
  static void alloc_shm_sysv(void **ptr, size_t size, int *ret);
  static void dealloc_shm_sysv(void **shm);
#else
#ifdef SHM_POSIX
  static void alloc_shm_posix(void **ptr, size_t size, int *ret);
  static void dealloc_shm_posix(void **shm, int64_t *nelem, MPI_Datatype type);
  static int get_randomTag( MPI_Comm comm, char * tag, size_t len );
#ifdef SHM_POSIX_BGQ
  static void clean_shm_posix(void);
#endif
#endif
#endif

  static void set_shm_numablock(void);

/*
 The length of the string to write in NUMA group tag. The upper limit of
 99999 for the number of cores in a NUMA node should be enough for a while.
*/
#define NUMATAG_WIDTH 5

void SET_SHM_DEBUG(void) { shm_debug = 1; }

#ifdef SHM_SYSV

/*
  System V version using MPI_Get_procesor_name
*/

  void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
    int nlen_processor, nlen;
    #ifdef DBG
      int mynid;
    #endif
    int i, n, found, nn=0, np=0;
    char nodnam[MPI_MAX_PROCESSOR_NAME+NUMATAG_WIDTH], *cnidlst, *tc1, *tc2;
    char string[errlen];

    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    /*
     * We can't just call 'MPI_Get_processor_name(nodnam, &nlen_processor);'
     * here because there are machines whose node names have different string
     * lengths, producing 'Invalid count argument' or 'Message truncated'
     * error messages at runtime. Instead we get the node name and compute
     * the maximum length with MPI_Allreduce.  RQH 07.2009.
     */
    MPI_Get_processor_name(nodnam, &nlen_processor);

    #ifdef DBG
      printf(" mype %i procname %s \n", mype, nodnam);
    #endif
    MPI_Allreduce(&nlen_processor, &nlen, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    set_shm_numablock();

    if ( numablk > 0 ) {
      /*
       * We assume that the MPI ranks are assigned consecutively on the
       * node's cores, i.e. rank(i)->c0, rank(i+1)->c1 ...  Distribution of
       * ranks can usually be controlled with MPI environment variables.
       */
      if ( nthreads == 0 ) {
        /* Group id */
        n = mype/numablk;
      } else {
       /* Group id with OpenMP (assumes that threads are on adjacent cores) */
        n = mype/(numablk/nthreads);
      }
      nlen += NUMATAG_WIDTH; /* append group id to processor name */
      sprintf(string,"%-d",n);
      strncat(nodnam,string,strlen(string));
    }

    #ifdef DBG
      if (shm_debug) {
        printf("pe %d of %d, nodnam = %s (len = %d), node = %d\n", mype, npes,
           nodnam, nlen, n);
        MPI_Barrier(MPI_COMM_WORLD);
      }
    #endif

      if ( (cnidlst = (char *) malloc((npes*nlen+1)*sizeof(char))) == NULL){
        fprintf(stderr,"failed to allocate cnidlst pe %d\n",mype);
        MPI_Abort(MPI_COMM_WORLD,-1);
      }
    MPI_Allgather(&nodnam, nlen, MPI_CHAR, cnidlst, nlen, MPI_CHAR,
       MPI_COMM_WORLD);

    nn = np = 0;

    for (i=0; i<npes; i++){
      tc1=cnidlst+i*nlen;
      if( strncmp(nodnam,tc1,nlen) == 0 ) nodpes[np++]=i;
      found = 0;
      for (n=0; n<nn; n++) {
        tc2=cnidlst+nodlst[n]*nlen;
        if ( strncmp(tc1,tc2,nlen) == 0 )  found = 1;
      }
      if (!found) nodlst[nn++] = i;
    }

    *nnod = nn;
    *nppn = np;
    mast = nodpes[0];

    #ifdef DBG
      if (shm_debug) {
        am_nodemaster = (mype == nodpes[0]);
        printf("mype %d nn %d np %d ",mype, nn, np);
        for (i=0; i<np; i++) printf("%d ",nodpes[i]);
        printf(" ");
        for(i=0; i<nlen*npes; i++) printf("%c",cnidlst[i]); printf("\n");
        for (n=0; n<nn; n++) {
          if (am_nodemaster && nodlst[n]==mype) {
            printf("node %d (%s) of %d: master_pe = %d of (",n,nodnam,nn,
               nodpes[0]);
            for (i=0; i<np; i++) printf(" %d", nodpes[i]);
            printf(" )\n");
          }
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }
    #endif

    MPI_Comm_group(MPI_COMM_WORLD,&group_world);
    MPI_Group_incl(group_world,np,nodpes,&node_group);
    MPI_Comm_create(MPI_COMM_WORLD,node_group,&node_comm);

    free(cnidlst);

  }

#else /* #ifdef SHM_SYSV */

  #ifdef SHM_POSIX
    #ifdef SHM_POSIX_BGQ

   /* File name root for POSIX shared memory. */
      static char posix_shm_fn[256]="/casino.shm";
      static int shmfd; /* file descriptor */
      static off_t shmsize;
      static off_t shmoff = 0;
      static void *shmbase = NULL; /* base  address of shared file */
      static int is_shmsize_set = 0;

      void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
        #ifdef DBG
          int mynid;
        #endif
	  int i, n, found, nn=0, np=0, tag0=-1;
	  char string[errlen];
	  int *nd_shm, *tag;
	  size_t serrlen =  errlen;

	  /* Add random string to the share file name */
	  get_randomTag(MPI_COMM_WORLD, posix_shm_fn, sizeof(posix_shm_fn));

	  MPI_Comm_rank(MPI_COMM_WORLD, &mype);
	  MPI_Comm_size(MPI_COMM_WORLD, &npes);

      /* Get the amount of shared memory available */
	  {
	    if (!is_shmsize_set){
	      char *seg_size;
	      seg_size = getenv("CASINO_MAXSHM");
	      if ( seg_size == NULL){
		fprintf(stderr,"CASINO_MAXSHM not set, quitting ... \n");
		MPI_Abort(MPI_COMM_WORLD,-1);
	      }
	      else {
               shmsize = atoll(seg_size);
               shmsize = shmsize << 20;
              }
	    }
	    /* if (mype == 0) fprintf(stderr,"shm size %lld %s\n", (long long int) shmsize, posix_shm_fn); */
	    size_t sizeaux = npes * sizeof(int);
	    if (sizeaux > shmsize ) {
	      fprintf(stderr,"rank %d: not enough shared memory to find smp ranks: shmsize %lld needed %ld\n", mype, shmsize, sizeaux);
	      MPI_Abort(MPI_COMM_WORLD,-1);
	    }
	  
	  }
	  snprintf(string,serrlen,"pe %d : shm_open tmp",mype);
	  if ((shmfd=shm_open(posix_shm_fn,(O_CREAT|O_RDWR), 0666)) < 0) perror(string);
	  
	  snprintf(string,serrlen,"pe %d : ftruncate tmp",mype);
	  if (ftruncate(shmfd,shmsize) < 0) perror(string);
	  /* MPI_Barrier? */
	  MPI_Barrier(MPI_COMM_WORLD);
	  
	  snprintf(string, serrlen,"pe %d : mmap init",mype);
	  if ((shmbase = mmap(NULL, shmsize, (PROT_READ|PROT_WRITE),
			      MAP_SHARED, shmfd, 0))< 0) perror(string);
	  
	  /* Used to find the ranks on a node */
	  nd_shm = (int *) shmbase;
	    
      /*
       * According to the standard, ftruncate fills the file with 0, hence the 
       * loop from below is not needed.
       *
       * Set nd_shm to -1 an unnecessary number of times.
       * Probably this is faster than going to disk.
       * for(i=0;i<npes;i++) nd_shm[i]=-1;
       */

      #ifdef DBG
	  printf("pe %d : 1 nd_shm :",mype);
	  for (i=0;i<npes;i++) printf(" %d",nd_shm[i]);
	  printf("\n");
      #endif

	  if ( (tag = (int *)malloc(npes*sizeof(int))) == NULL ) {
	    fprintf(stderr,"failed to allocate tag pe %d\n",mype);
	    MPI_Abort(MPI_COMM_WORLD,-1);
	  }
	  
	for(i=0;i<npes;i++) tag[i]=-1;

      /* Each rank colours its position */
	int color=12345;
	nd_shm[mype]=color;
      
      #ifdef DBG
	printf("pe %d : 2 nd_shm :",mype);
	for (i=0;i<npes;i++) printf(" %d",nd_shm[i]);
	printf("\n");
      #endif

	MPI_Barrier(MPI_COMM_WORLD);
        /* Tagging with mype of the first rank found in nd_shm. */
	for(i=0; i<npes; i++) 
	  if(nd_shm[i] == color) { tag0=i; break; }
        
        /* Each rank gets a full version of the tag array. */
	  MPI_Allgather(&tag0,1,MPI_INT,tag,1,MPI_INT,MPI_COMM_WORLD);
	  
	  /* Find the list of ranks on the node and the list of nodes */
	  nn=np=0;
	  for(i=0;i<npes;i++){
	    if(tag[mype] == tag[i]) nodpes[np++]=i;
	    found=0;
	    for(n=0;n<nn;n++) {
		if(tag[i] == tag[nodlst[n]]) found = 1;
	    }
	    if (!found) nodlst[nn++]=i;
	  }
	    
	  *nnod = nn;
	  *nppn = np;
	  mast = nodpes[0];
	    
          #ifdef DBG
	    if (shm_debug) {
	      am_nodemaster = (mype == nodpes[0]);
	      printf("mype %d nn %d np %d ",mype, nn, np);
	      for (i=0; i<np; i++) printf("%d ",nodpes[i]);
	      printf(" ");
	      for(i=0; i<npes; i++) printf("%d ",tag[i]); printf("\n");
	      for (n=0; n<nn; n++) {
		if (am_nodemaster && nodlst[n]==mype) {
		  printf("node %d tag %d of %d: master_pe = %d of (", n, tag[mype],
			 nn, nodpes[0]);
		  for (i=0; i<np; i++) printf(" %d", nodpes[i]);
		  printf(" )"); printf(" nodlst");
		  for (i=0; i<nn; i++) printf(" %d", nodlst[i]); printf("\n");
		}
		MPI_Barrier(MPI_COMM_WORLD);
	      }
	    }
          #endif
	    
	    MPI_Comm_group(MPI_COMM_WORLD, &group_world);
	    MPI_Group_incl(group_world, np, nodpes, &node_group);
	    
	    MPI_Comm_create(MPI_COMM_WORLD, node_group, &node_comm);


      /* Free tag memory */
	    free(tag);
      }
									    
    #else /*ifdef SHM_POSIX_BGQ */

    /* POSIX version. */

    /* File name root for POSIX shared memory. */
    static char posix_shm_fn[256]="/casino.shm";

    void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
      #ifdef DBG
        int mynid;
      #endif
      int i, is, ie, n, fd, shift, found, nn=0, np=0, tag0=-1;
      char string[errlen];
      char fname[256];
      off_t size;
      int *nd_shm, *tag;
      
      /* Add random string to the share file name */
      get_randomTag(MPI_COMM_WORLD, posix_shm_fn, sizeof(posix_shm_fn));
      (void)strcpy(fname,posix_shm_fn);
      (void)strncat(fname,"_tmp",4);

      MPI_Comm_rank(MPI_COMM_WORLD, &mype);
      MPI_Comm_size(MPI_COMM_WORLD, &npes);

      snprintf(string, errlen, "pe %d : shm_open",mype);
      if ((fd=shm_open(fname,(O_CREAT|O_RDWR), 0666)) < 0) perror(string);

      size=(off_t)npes*(off_t)sizeof(int);
      snprintf(string, errlen, "pe %d : ftruncate",mype);
      if (ftruncate(fd,size) < 0) perror(string);

      snprintf(string, errlen, "pe %d : mmap",mype);
      if ((nd_shm = (int *)mmap(NULL, size, (PROT_READ|PROT_WRITE),
         MAP_SHARED, fd, 0))< 0) perror(string);

      /*
       * According to the standard, ftruncate fills the file with 0, hence the 
       * loop from below is not needed.
       *
       * Set nd_shm to -1 an unnecessary number of times.
       * Probably this is faster than going to disk.
       * for(i=0;i<npes;i++) nd_shm[i]=-1;
       */

      #ifdef DBG
        printf("pe %d : 1 nd_shm :",mype);
        for (i=0;i<npes;i++) printf(" %d",nd_shm[i]);
        printf("\n");
      #endif

      /* Each rank puts a 11111 to its position */
      int color=12345;
      nd_shm[mype]=color;

      set_shm_numablock();

      /* MPI_Barrier? */
      MPI_Barrier(MPI_COMM_WORLD);

      #ifdef DBG
        printf("pe %d : 2 nd_shm :",mype);
        for (i=0;i<npes;i++) printf(" %d",nd_shm[i]);
        printf("\n");
      #endif

      if(numablk == 0) {
        /* Tagging with mype of the first rank found in nd_shm. */
        for(i=0; i<npes; i++) {
          if(nd_shm[i] == color) { tag0=i; break; }
        }
        /* Unlink here, one core per node. */
        if (mype == tag0) {
          snprintf(string, errlen, "pe %d : shm_unlink",mype);
          if (shm_unlink(fname) < 0) perror(string);
        }
      } else {
        /* Test if the ranks are contiguous in nd_shm. */
        nn=0; found=1;
        for(i=0; i<npes; i++) {
          if(nd_shm[i] == color) {
            nn++;
            if(found) { found=0; is=i; }
            ie=i;
          }
          /* Unlink here, one core per node */
          if (mype == is) {
            snprintf(string, errlen, "pe %d : shm_unlink",mype);
            if (shm_unlink(fname) < 0) perror(string);
          }
        }
        #ifdef DBG
          printf("pe %d : is %d ie %d nn %d \n",mype,is,ie,nn);
        #endif

        if( (ie-is+1) != nn ) {
          /*
           * Problem: the ranks are not sequential on the node, so the
           * algorithm won't work. Quit.
           */
          fprintf(stderr,
             "pe %d : ranks are not sequentially placed on the NUMA node.\n"
             "In nd_shm we got start %d end %d total %d \n  "
             "Try to use MPI environment variable place sequential ranks\n"
             "on the same node. Exiting!", mype, is, ie, nn);
          MPI_Abort(MPI_COMM_WORLD,-1);
        } else {
          #ifdef DBG
            printf("pe %d : tag0 %d, is %d, ie %d \n ", mype,tag0,is,ie);
          #endif
          shift=numablk;
          if (nthreads>0) shift=numablk/nthreads;

          #ifdef DBG
            fprintf(stderr,"pe %d : shift %d \n", mype,shift);
          #endif

          tag0 = is+(mype-is)/shift;
        }
      }

      #ifdef DBG
        fprintf(stderr,"pe %d : tag0 %d \n", mype, tag0);
      #endif


      /* Each rank gets a full version of the tag array. */
      if ( (tag = (int *)malloc(npes*sizeof(int))) == NULL ) {
        fprintf(stderr,"failed to allocate tag pe %d\n",mype);
        MPI_Abort(MPI_COMM_WORLD,-1);
      }
      for(i=0;i<npes;i++) tag[i]=-1;
      MPI_Allgather(&tag0,1,MPI_INT,tag,1,MPI_INT,MPI_COMM_WORLD);

      /* Find the list of ranks on the node and the list of nodes */
      nn=np=0;
      for(i=0;i<npes;i++){
        if(tag[mype] == tag[i]) nodpes[np++]=i;
        found=0;
        for(n=0;n<nn;n++) {
          if(tag[i] == tag[nodlst[n]]) found = 1;
        }
        if (!found) nodlst[nn++]=i;
      }

      *nnod = nn;
      *nppn = np;
      mast = nodpes[0];

      #ifdef DBG
        if (shm_debug) {
          am_nodemaster = (mype == nodpes[0]);
          printf("mype %d nn %d np %d ",mype, nn, np);
          for (i=0; i<np; i++) printf("%d ",nodpes[i]);
          printf(" ");
          for(i=0; i<npes; i++) printf("%d ",tag[i]); printf("\n");
          for (n=0; n<nn; n++) {
            if (am_nodemaster && nodlst[n]==mype) {
              printf("node %d tag %d of %d: master_pe = %d of (", n, tag[mype],
                 nn, nodpes[0]);
              for (i=0; i<np; i++) printf(" %d", nodpes[i]);
              printf(" )"); printf(" nodlst");
              for (i=0; i<nn; i++) printf(" %d", nodlst[i]); printf("\n");
            }
            MPI_Barrier(MPI_COMM_WORLD);
          }
        }
      #endif

      MPI_Comm_group(MPI_COMM_WORLD, &group_world);
      MPI_Group_incl(group_world, np, nodpes, &node_group);
      MPI_Comm_create(MPI_COMM_WORLD, node_group, &node_comm);

      /* Add NUMATAG to the POSIX file name. */
      if (numablk>0) {
        /*
         * Note: the following formula generates different shm file names on
         * different nodes. Something like (mype%nppn)/shift would be enough
         * but this assumes that every node has the same number of cores and
         * is fully populated.
         */
        sprintf(string,"%-d",(mype/shift));
        strncat(posix_shm_fn,string,strlen(string));
      }

      /* Unmap nd_shm. */
      snprintf(string, errlen, "pe %d : munmap size %ld",mype,size);
      if (munmap((void *)nd_shm, size) < 0) perror(string);
      snprintf(string, errlen, "pe %d : close size %ld",mype,size);
      if (close(fd) < 0) perror(string);
      /* Free tag pointer */
      free(tag);
      MPI_Barrier(node_comm);

    }
    #endif /* ifdef SHM_POSIX_BGQ */
  #else /* #ifdef SHM_POSIX */

    void GET_SMP_LIST(int *nnod, int *nodlst, int *nppn, int *nodpes) {
      int mype;
      MPI_Comm_rank(MPI_COMM_WORLD, &mype);
      if (mype==0) printf("(%s\n%s\n%s\n",
         "The Makefile flag CFLAGS_SHM contains neither",
         "-DSHM_SYSV, -DSHM_POSIX, nor -DSHM_POSIX_BGQ.",
         "Quitting ...\n");
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD,-1);
    }

  #endif /* #ifdef SHM_POSIX */

#endif /* #ifdef SHM_SYSV */


void ALLOC_SHM(void **ptr, int64_t *nelem, MPI_Fint *ftype, int *ret) {
  int typsiz;
  size_t size;
  MPI_Datatype type;
  err = 0;
  type = MPI_Type_f2c( *ftype);
  MPI_Type_size(type, &typsiz);
  #ifdef DBG
    if (shm_debug) printf("pe %d : type size = %d\n", mype, typsiz);
  #endif
  size = (size_t)*nelem * (size_t)typsiz;
  #ifdef DBG
    if (shm_debug) printf("pe %d : full size = %ld\n", mype, size);
  #endif
  #ifdef SHM_SYSV
    #ifdef DBG
      printf("pe %d : using System V shm\n", mype);
    #endif
    alloc_shm_sysv(ptr, size, ret);
  #endif
  #ifdef SHM_POSIX
    #ifdef DBG
      printf("pe %d : using POSIX shm\n",mype);
    #endif
    alloc_shm_posix(ptr, size, ret);
  #endif
  if (*ret) *ret = err;
}


void DEALLOC_SHM(void **shm, int64_t *nelem, MPI_Fint *ftype) {
  MPI_Datatype type;
  type = MPI_Type_f2c( *ftype);
  #ifdef DBG
    if (shm_debug) printf("pe %d : dealloc_shm: address = %lx nelem %ld\n",
       mype, *shm, *nelem);
  #endif
  #ifdef SHM_SYSV
    dealloc_shm_sysv(shm);
  #else
    #ifdef SHM_POSIX
      dealloc_shm_posix(shm, nelem, type);
    #endif
  #endif
}

void CLEAN_SHM(void) {
  #ifdef SHM_POSIX_BGQ
      clean_shm_posix();
  #endif
}


  void SET_SHMSIZE(int64_t *s){
  #ifdef SHM_POSIX_BGQ
    /* warn if this function is called more than once */
    if (is_shmsize_set)
      fprintf(stderr,"Warning - shared memory already set, ignoring ... \n");
    else{
      shmsize = *s;
      is_shmsize_set = 1;
    }
  #endif
  }

  void GET_SHMSIZE(int64_t *s){
  #ifdef SHM_POSIX_BGQ
    /* Returns the amount of used shared memory, i.e. shmoff */
    *s = (int64_t) shmoff; 
  #endif
  }



#ifdef SHM_SYSV

  /* Implementation using System V shared memory [shmget & shmat]. */

  static void alloc_shm_sysv(void **ptr, size_t size, int *ret) {
    char string[errlen];
    int shmid;
    void *shm;
    struct shmid_ds ds;

    if (mype == mast) {
      if ((shmid = shmget(IPC_PRIVATE, size, 0666)) < 0) {
        sprintf(string,"pe %d: shmget",mype);
        perror(string);
        if(errno==EINVAL) printf(
           "EINVAL: requested size of shared memory segment too large (%lx).\n"
           "Ask administrator to increase system-wide limit.\n"
           "E.g.: on Linux use 'sysctl kernel.shmmax=$((1<<63))'\n", size);
        else if(errno==ENOMEM) printf(
           "ENOMEM: memory allocation for control structures failed.\n");
        else if(errno==ENOSPC) printf(
           "ENOSPC: either the system-wide limit for shared memory (SHMALL)\n"
           "or the maximum number of SHM ids (SHMMNI) has been reached.\n"
           "Maybe some other process is not cleaning up properly?\n"
           "Clean up, reboot or increase limits (e.g., on Linux use sysctl)\n");
        else printf(
           "Have hit an unidentified error when calling shmget.");
        if (*ret) err++; else exit(1);
      }
      MPI_Bcast(&shmid, 1, MPI_INT, 0, node_comm);
      shm = shmat(shmid, NULL, 0);
      if (shm == (void *) -1) {
        sprintf(string,"pe %d: shmat",mype);
        perror(string);
        if (*ret) err++; else exit(1);
      }

      /* Deallocate segment after last process has been detached. */
      if (shmctl(shmid, IPC_RMID, &ds) < 0) {
        sprintf(string,"pe %d: shmget",mype);
        perror(string);
        if (*ret) err++; else exit(1);
      }

    } else {

      MPI_Bcast(&shmid, 1, MPI_INT, 0, node_comm);
      shm = shmat(shmid, NULL, 0);
      if (shm == (void *) -1) {
        sprintf(string,"pe %d: shmat",mype);
        perror(string);
        if (*ret) err++; else exit(1);
      }
    }
    MPI_Barrier(node_comm);
    *ptr = shm;
  }


  static void dealloc_shm_sysv(void **ptr) {
    char string[errlen];
    snprintf(string,errlen, "pe %d: shmdt",mype);
    if (shmdt(*ptr) < 0) {
      perror(string);
      exit(1);
    }
  }

#else /* #ifdef SHM_SYSV */

  #ifdef SHM_POSIX
    #ifdef SHM_POSIX_BGQ

     
    struct list_t {
      off_t off;
      size_t size;
      void *p;
      struct list_t *next;
    };

    static struct list_t *head, *tail;
    
    static void alloc_shm_posix(void **ptr, size_t size, int *ret) {
      void *shm=NULL;
      int node_id;
      off_t offset;
      struct list_t *p;

      MPI_Comm_rank(node_comm, &node_id);
      int64_t psize = sysconf(_SC_PAGE_SIZE);
      int got_file = 0;
      /* Search for if there is already an available file */

      if (node_id == 0 && shmoff != 0){
        p = head;
        /* fprintf(stderr,"in alloc file check %d %d\n", shmid, got_file); */
        do {
          if ((p->p == NULL && size <= p->size) ){
            /* fprintf(stderr,"in alloc file check %d %ld %ld %ld\n", mype , p->size, size, p->off); */
            offset = p->off;
            /* fprintf(stderr,"alloc reuse %d %d \n", mype, shmid_reuse); */
            got_file = 1;
            break;
          }
          /* fprintf(stderr,"in alloc file check %p %p %p\n", p, p->p, p->next); */
          p = p->next;
        }while(p != NULL);
      }
      /* fprintf(stderr,"in alloc after file check %d %d %d \n", mype, shmid, got_file); */
      MPI_Bcast(&got_file, 1, MPI_INT, 0, node_comm);
      if (got_file){
        MPI_Bcast(&offset, 1, MPI_LONG_LONG_INT, 0, node_comm);
	
	/* if (mype == 0) 
	  fprintf(stderr,"rank %d open reuse offset %lld %ld %p %p %p\n", mype, offset, size, p, tail, head); */

	shm = (char *) shmbase + offset;  

        if (node_id == 0) {
          p->p = shm;
        }
        MPI_Barrier(node_comm);
        *ptr = shm;
        return ;
      }

      /* Existing segments are too small or no segment; get a new one at the tail*/

      shm = (char *) shmbase + shmoff; 
      if ( (char *) shm + size > (char *) shmbase + shmsize ) {
	fprintf(stderr, "rank %d: trying to get %llu segment, shared memory allocation %llu exceeded. Quitting ...\n", mype, size, shmsize);
	MPI_Abort(MPI_COMM_WORLD, -1);
      }

      /* Store file descriptor, index and pointer value */
     if ( node_id == 0){
        if (shmoff == 0){
          /* first element in the linked list */
          head = malloc(sizeof(struct list_t));
          tail = head;
        } else{
          tail->next = malloc(sizeof(struct list_t));
          tail = tail->next;
        }
        tail->p=shm;
        tail->size=size + (psize - size%psize)%psize; /* Rounding to the page size */
	tail->off=shmoff;
        tail->next=NULL;
       
        /* print fn and shmid on rank 0
        if (mype == 0)
	  fprintf(stderr,"rank %d open shmid %lld %lld %d %p %p %p\n",  mype, size, shmoff, psize, shm, shmbase, tail); */
     }
     /* Offsets are aligned to the page boundary */
     shmoff += size + (psize - size%psize)%psize;
     *ptr = shm;
    }

     
#ifdef DBG
#include <sys/stat.h>
#endif

    static void dealloc_shm_posix(void **shm, int64_t *nelem,
       MPI_Datatype type) {

      void *ps = *shm;
      *shm = NULL;

      /* Mark the corresponding element of the linked list as available*/
      /* if it is the last one remove it */
      int node_id, got_tail;
      off_t newoff;
      MPI_Comm_rank(node_comm, &node_id);
      got_tail = 0;
      if (node_id == 0){
	struct list_t *p, *pp;
	p = pp = head;
	do {

	  if ( p->p == ps ){
	    
	    /* if ( mype == 0)
	     fprintf(stderr,"rank %d unmap size %ld offset %ld ptr %p %p %p %p\n", mype, size, p->off, p->p, p, tail, head); */
	    
	    p->p = NULL;
	    
	    if ( p == tail ) {
	      got_tail=1;
	      newoff = p->off;
	      tail = pp;
              tail->next = NULL;
	      free(p);
	    }               
	    break;
	  }
	  if (p->next == NULL){
	    fprintf(stderr,"Linked list finished in shm_dealloc, no pointer found\n");
	    MPI_Abort(MPI_COMM_WORLD,-1);
	  }
	  pp = p; p = p->next;
	} while(1);
      }
      
      MPI_Bcast(&got_tail,1,MPI_INT,0,node_comm);
      if (got_tail) {
	MPI_Bcast(&newoff,1,MPI_LONG_LONG_INT,0,node_comm);
	shmoff = newoff;
      }
    }

    static void clean_shm_posix(void){
      char string[errlen];
      int node_id;
      size_t serrlen = errlen;

      MPI_Comm_rank(node_comm, &node_id);
      if (node_id == 0){
	/* struct list_t *p, *pp;
	p = pp = head; */
	do {
	  
	  snprintf(string,serrlen,"pe %d : clean munmap size %lld",mype,shmsize);
	  if (munmap(shmbase, shmsize) < 0) perror(string);
      
	  snprintf(string,serrlen,"pe %d: shm_unlink in clean shm",mype);
	  if (shm_unlink(posix_shm_fn) < 0 ){
	    perror(string);
	  }
	  snprintf(string,serrlen,"pe %d: close in clean shm",mype);
	  if (close(shmfd) < 0 ){
	    perror(string);
	  }
	  /* fprintf(stderr,"clean_shm  %s\n", posix_shm_fn);
	  pp=p; p=p->next; free(pp); */
	}while(0);
      }
    } 

    #else 
    /* Implementation using POSIX shared memory [shm_open & mmap]. */

    static void alloc_shm_posix(void **ptr, size_t size, int *ret) {
      void *shm=NULL;
      char string[errlen];
      int fd,node_id;
      static int shmid = 0;
      char fn[256];

      MPI_Comm_rank(node_comm, &node_id);
      MPI_Bcast(&shmid, 1, MPI_INT, 0, node_comm); /* shmid is increased at the bottom of this function */
      
	/* Generate a unique file name.
	   This algorithm assumes that a node runs only one CASINO process*/
      snprintf(fn, 255,"%s%s%-d",posix_shm_fn,"_",shmid); /* Should I use snprintf ? */
	
      /* Open the shm file in exclusive mode on nodes' 0 rank */
      
      if (node_id == 0){
	snprintf(string, errlen, "pe %d: shm_open excl",mype); 
	if ((fd = shm_open(fn, (O_CREAT|O_EXCL|O_RDWR), 0666)) < 0) {
	  perror(string);
	  if (*ret) err++; else exit(1);
	}
	snprintf(string, errlen, "pe %d: ftruncate",mype);
	if (ftruncate(fd, (off_t)size) < 0) {
	  perror(string);
	  if (*ret) err++; else exit(1);
	}
      }
      MPI_Barrier(node_comm);
      /* Attach the other ranks to the open file */
      if (node_id != 0)
	snprintf(string, errlen, "pe %d: shm_open attach",mype);
	if ((fd = shm_open(fn, (O_RDWR), 0666)) < 0) {
	  perror(string);
	  if (*ret) err++; else exit(1);
	} 
      #ifdef DBG
        fprintf(stderr,"pe %d: fd %d shared posix file name %s size %ld\n ",
           mype, fd, fn, size);
      #endif
     
      snprintf(string, errlen, "pe %d: mmap",mype);
      if ((shm = mmap(NULL, size, (PROT_READ|PROT_WRITE), MAP_SHARED, fd, 0))
         == (void *) -1) {
        perror(string);
        if (errno==ENOMEM) printf(
           "ENOMEM: requested shared memory size exceeds system limit\n");
        else if (errno==EBADF) printf(
           "EBADF: wrong value for fd argument. \n");
        else if (errno==EACCES) printf(
           "EACCES: wrong access flags.\n");
        else printf(
           "Have hit an unidentified error when calling mmap.");
        if (*ret) err++; else exit(1);
      }
      #ifdef DBG
        if (shm_debug) printf(" z: pe %d: mmap: shm = %lx\n",mype,shm);
      #endif

      snprintf(string, errlen, "pe %d: close in alloc",mype);
      if (close(fd) < 0 ){
         perror(string);
         if (*ret) err++; else exit(1);
      }

      MPI_Barrier(node_comm);

      /* Unlink, only on master of each node. */
      if (node_id == 0) {
        snprintf(string, errlen, "pe %d: shm_unlink",mype);
        if (shm_unlink(fn) < 0 ){
          perror(string);
          if (*ret) err++; else exit(1);
        }
	shmid++; /* next file id*/
      }

      *ptr = shm;
    }


    static void dealloc_shm_posix(void **shm, int64_t *nelem,
       MPI_Datatype type) {
      char string[errlen];
      int typsiz;
      size_t size;

      MPI_Type_size(type, &typsiz);
      size = (size_t)*nelem * (size_t)typsiz;

      #ifdef DBG
        printf("pe %d : munmap size %ld %d %ld shm %lx\n", mype, *nelem,
           typsiz, size, *shm);
      #endif

      snprintf(string, errlen, "pe %d : munmap size %ld", mype, size);
      if (munmap(*shm, size)<0) {
        perror(string); exit(1);
      }
    }
    #endif /* ifdef SHM_POSIX_BGQ */
 
  #endif /* #ifdef SHM_POSIX */

#endif /* #ifdef SHM_SYSV */


void set_shm_numablock(void) {
  /* Pick CASINO_NUMABLK and OMP_NUM_THREADS from environment. */
  char *nb;
  nb = getenv("CASINO_NUMABLK");
  if ( nb != NULL ) numablk=atoi(nb);
  #ifdef _OPENMP
    #pragma omp parallel
    {
      #pragma omp master
      nthreads=omp_get_num_threads();
    }
  #endif
  if (shm_debug) printf("numa_block %d : num OpenMP threads %i \n", numablk,
     nthreads);
}

void GET_NNPSMP(int *nnpsmp) {
 /* Get nnpsmp=CASINO_NUMABLK for CASINO allocations. */
  char *nb;
  nb = getenv("CASINO_NUMABLK");
  if ( nb != NULL ) 
   *nnpsmp=atoi(nb); 
  else {
   *nnpsmp=0;
  }
}

#ifdef SHM_POSIX
  static int get_randomTag( MPI_Comm comm, char * tag, size_t stag){
 /* appends a random string 32 characters long to tag but not more
  * than len(tag) -size(tag) -1 
  *
  * Note: on systems that reboot the compute node at the start of
  * each computation aux[1] could have the same value for different
  * runs (e.g. BGQ) but this is not a problem. PIDs will certainly be
  * different on a workstation if two or more users run CASINO at the
  * same time.
    */ 
    int rank;
    double t1, t2, tt;
    pid_t pid;
    long unsigned int tagProd[2], aux[2];
    char s[33];
    
    t1 = MPI_Wtime();
    MPI_Comm_rank(comm, &rank);
    pid = getpid();
    tt = MPI_Wtick();
    t2 = MPI_Wtime();


    /* The shift pushes some bits to the left in aux[0] for more entropy */
    int maxShift = 8 * sizeof(aux[0]) / 2;
    
    aux[1] = (long unsigned int ) pid << ( (3*rank) % maxShift);
    
    aux[0] = (long unsigned int ) ((t2-t1) / tt ) << ( (5*rank) % maxShift);
    
    MPI_Allreduce(aux, tagProd, 2, MPI_UNSIGNED_LONG, MPI_SUM, comm);
    
    sprintf(s,"%016lx%016lx", tagProd[0], tagProd[1]);
    size_t lentag = strnlen(tag, stag);
    int  diff = stag - lentag - 1;
    if ( diff > 0)
      strncat(tag, s, diff);
    
    return (0); /* some error  message when diff <= 0 */
  }
#endif
