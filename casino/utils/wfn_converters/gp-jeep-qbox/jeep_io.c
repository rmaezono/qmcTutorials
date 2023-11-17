/******************************************************************************
  Plug in routine to read in JEEP wavefunctions in binary format.
  First version AJW (2000), improved by BM (2001).
 *****************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/* swap integer and return new one */
int SwapInt(int n) {
  unsigned char *c;
  unsigned char tmp;
  c=(unsigned char*)(&n);
  tmp = c[3]; c[3] = c[0]; c[0] = tmp;
  tmp = c[2]; c[2] = c[1]; c[1] = tmp;
  return n;
}

/* swap the double value passed in, no return value */
void SwapDouble(double *d) {
  unsigned char *c;
  unsigned char tmp;
  c=(unsigned char *)d;
  tmp = c[7]; c[7] = c[0]; c[0] = tmp;
  tmp = c[6]; c[6] = c[1]; c[1] = tmp;
  tmp = c[5]; c[5] = c[2]; c[2] = tmp;
  tmp = c[4]; c[4] = c[3]; c[3] = tmp;
}

void readwf1_(int *nstr,int *ngwr)
{
  FILE *wfinfile;
  int ns,ng,ns2,ng2;
  long tsize,tsize2,size;
  int swapFlag;

  wfinfile = fopen( "jeep.wf", "r" );
  if ( !wfinfile )    {
    fprintf( stderr,"Could not open jeep.wf\n");
    return;}
  if (fread( &ns, sizeof(int), 1, wfinfile ) != 1) goto error1;
  if (fread( &ng, sizeof(int), 1, wfinfile ) != 1) goto error2;

  ns2=SwapInt(ns);
  ng2=SwapInt(ng);

  tsize = 2*sizeof(int) + /* size of nst, ngw */
    2*ns*ng*sizeof(double); /* size of 2 complex wfs */
  tsize2 = 2*sizeof(int) + /* size of nst, ngw */
    2*ns2*ng2*sizeof(double); /* size of 2 complex wfs */

  fseek(wfinfile, 0L, SEEK_END);
  size = ftell(wfinfile);

  if (tsize!=size && tsize2!=size && (tsize-8)*2+8!=size && (tsize2-8)*2+8!=size)
    {
    printf("File size error: Actual size: %ld Computed size %ld or swapped %ld\n",size,tsize,tsize2);
    exit(1);
  }
  swapFlag = (tsize2==size || (tsize2-8)*2+8==size);
  if (swapFlag) {
    printf("Converting wavefunction by swapping\n");
    ns=ns2;
    ng=ng2;
  }

  *nstr=ns; *ngwr=ng;
  fclose( wfinfile );

  return;

error1:
 fclose (wfinfile);
 return;

error2:
 fclose (wfinfile);
 return;

}

void readwf2_(double *cg,int *istate)
{
  FILE *wfinfile;
  int ns,ng,ns2,ng2,is,swapFlag,i;
  long ipos,tsize,tsize2,size;

  wfinfile = fopen( "jeep.wf", "r" );
  if ( !wfinfile )    {
    fprintf( stderr,"Could not open jeep.wf\n");
    return;}
  if (fread( &ns, sizeof(int), 1, wfinfile ) != 1) goto error1;
  if (fread( &ng, sizeof(int), 1, wfinfile ) != 1) goto error2;

  ns2=SwapInt(ns);
  ng2=SwapInt(ng);

  tsize = 2*sizeof(int) + /* size of nst, ngw */
    2*ns*ng*sizeof(double); /* size of 2 complex wfs */
  tsize2 = 2*sizeof(int) + /* size of nst, ngw */
    2*ns2*ng2*sizeof(double); /* size of 2 complex wfs */

  fseek(wfinfile, 0L, SEEK_END);
  size = ftell(wfinfile);

  if (tsize!=size && tsize2!=size && (tsize-8)*2+8!=size && (tsize2-8)*2+8!=size) {
    printf("File size error: Actual size: %ld Computed size %ld or swapped %ld\n",size,tsize,tsize2);
    exit(1);
  }
  swapFlag = (tsize2==size || (tsize2-8)*2+8==size);
  if (swapFlag) {
    ns=ns2;
    ng=ng2;
  }

  /* Jump to the start of the state istate */
  is=*istate;
  ipos = 2*sizeof(int) + /* size of nstr, ngwr */
    2*(is-1)*ng*sizeof(double); /* size of 2 complex wfs */
  /*  fprintf(stderr,"Jumping to point %d\n",ipos); */

  fseek(wfinfile,ipos, SEEK_SET);
  if (fread( cg, sizeof(double), 2*ng, wfinfile ) != 2*ng) goto error3;
  fclose( wfinfile );

  if (swapFlag)
    for(i=0; i<2*ng; i++) {
      SwapDouble(&(cg[i]));
    }

 return;

error1:
 fclose (wfinfile);
 return;

error2:
 fclose (wfinfile);
 return;

error3:
 fclose (wfinfile);
 return;

}

/* Dummy routines: AIX doesn't put the _ after subroutine names,
   and MS C compiler turns names into uppercase */
void readwf1(int *nstr,int *ngwr)
{
  readwf1_(nstr,ngwr);
}
void readwf2(double *cg,int *istate)
{
  readwf2_(cg,istate);
}
void READWF1(int *nstr,int *ngwr)
{
  readwf1_(nstr,ngwr);
}
void READWF2(double *cg,int *istate)
{
  readwf2_(cg,istate);
}
