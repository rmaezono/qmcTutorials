/*------------------------------------------------------------------------*
 * Emulate FORTRAN etime extension for compilers that don't have it (e.g. *
 * NAG). Avoids having to use mpi_wtime or NAG library to get CPU times.  *
 *------------------------------------------------------------------------*/

/*
 * The Fortran compiler may expect C functions to have:
 * - no underscores appended to end of name [use -DF90_NO_UNDERSCORE]
 * - one underscore appended to end of name [default]
 * - two underscores appended to end of name [use -DF90_DOUBLE_UNDERSCORE]
 * and:
 * - name in lower case [default]
 * - name in upper case [use -DF90_CAPITALS]
 */
#ifndef ETIME
  #ifdef F90_DOUBLE_UNDERSCORE
    #ifdef F90_CAPITALS
      #define ETIME ETIME__
      #define DTIME DTIME__
    #else
      #define ETIME etime__
      #define DTIME dtime__
    #endif
  #else
    #ifdef F90_NO_UNDERSCORE
      #ifdef F90_CAPITALS
        #define ETIME ETIME
        #define DTIME DTIME
      #else
        #define ETIME etime
        #define DTIME dtime
      #endif
    #else
      #ifdef F90_CAPITALS
        #define ETIME ETIME_
        #define DTIME DTIME_
      #else
        #define ETIME etime_
        #define DTIME dtime_
      #endif
    #endif
  #endif
#endif

#ifdef CASINO_NO_ETIME

  /* Use this if the include files are just not there */
  float ETIME(float*);
  float DTIME(float*);


  float ETIME(float *t) {
    t[0]=(float)0;
    t[1]=(float)0;
    return 0;
  }


  float DTIME(float *t) {
    t[0]=(float)0;
    t[1]=(float)0;
    return 0;
  }


#else

  #include<sys/times.h>
  #include<time.h>

  #ifdef CASINO_T3E_MODE

    /* Use this if the C compiler is weird in the T3E way */
    float ETIME(double [2]);
    float DTIME(double [2]);


    float ETIME(double t[2]) {
      struct tms buffer;
      times(&buffer);

      t[0]=(float)buffer.tms_utime/CLK_TCK;
      t[1]=(float)buffer.tms_stime/CLK_TCK;
      return t[0]+t[1];
    }


    float DTIME(double t[2]) {
      static clock_t old_values[2] = { 0, 0};
      struct tms buffer;
      times(&buffer);

      t[0]=(float)(buffer.tms_utime-old_values[0])/CLK_TCK;
      t[1]=(float)(buffer.tms_stime-old_values[1])/CLK_TCK;
      old_values[0]=buffer.tms_utime;
      old_values[1]=buffer.tms_stime;
      return t[0]+t[1];
    }


  #else

    /* This is what will be run under normal conditions */
    #ifndef CLK_TCK
      /* From bits/time.h; accounts for POSIX.1-1988 -> ISO C99 transition */
      #include <bits/types.h>
      extern long int __sysconf (int);
      #define CLK_TCK ((__clock_t) __sysconf (2))
    #endif
    float ETIME(float*);
    float DTIME(float*);


    float ETIME(float *t) {
      struct tms buffer;
      times(&buffer);

      t[0]=(float)buffer.tms_utime/CLK_TCK;
      t[1]=(float)buffer.tms_stime/CLK_TCK;
      return t[0]+t[1];
    }


    float DTIME(float *t) {
      static clock_t old_values[2] = { 0, 0};
      struct tms buffer;
      times(&buffer);

      t[0]=(float)(buffer.tms_utime-old_values[0])/CLK_TCK;
      t[1]=(float)(buffer.tms_stime-old_values[1])/CLK_TCK;
      old_values[0]=buffer.tms_utime;
      old_values[1]=buffer.tms_stime;
      return t[0]+t[1];
    }


  #endif

#endif
