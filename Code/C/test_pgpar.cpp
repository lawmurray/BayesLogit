
// #include "PolyaGamma.h"
#include "PolyaGammaPar.h"
#include <stdio.h>

const int N = 0x10000;
const int nthreads = 4;

int main(void) {
  
  int M = 500;

  double samp[N];
  int    df[N];
  double psi[N];

  for (int i = 0; i < N; i++) {
    df [i] = 1;
    psi[i] = 0.0;
  }

  //------------------------------------------------------------------------------

  PolyaGammaPar pgpar;

  time_t start, end;

  start = time(NULL);

  for (int i = 0; i < M; i++)
    pgpar.draw(samp, df, psi, N, nthreads);

  end = time(NULL);
  double diff = (double)(end - start);
  printf("Time: %f sec. using %i threads for %i.\n", diff, nthreads, N);

  // Write it.
  FILE *file; 
  file = fopen("par.txt","w");
  for (int i = 0; i < N; i++) {
    fprintf(file,"%g\n", samp[i]);
  }
  fclose(file);

  //------------------------------------------------------------------------------

  PolyaGamma pg;
  RNG r;

  start = time(NULL);

  for (int j = 0; j < M; j++)
    for(int i = 0; i < N; i++)
      samp[i] = pg.draw(df[i], psi[i], r);

  end = time(NULL);
  diff = (double)(end - start);
  printf("Time: %f sec. (serial) for %i.\n", diff, N);

  // Write it.
  file = fopen("ser.txt","w");
  for (int i = 0; i < N; i++) {
    fprintf(file,"%g\n", samp[i]);
  }
  fclose(file);

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

// When N = 0x10000 and M = 500 we have

// Serial: 13 or 14 s.

// Par: 
// 1 thread: 14 s.
// 2 thread: 8 s.
// 3 thread: 5 s.
// 4 thread: 4 s.
