// Following Nathan Reed's article:
// http://reedbeta.com/blog/quick-and-easy-gpu-random-numbers-in-d3d11/
// PRNGs are designed to go deep, i.e. have good distributions when sequentially updating state
// Hashes are designed to go wide, i.e. have good distributions across initial seeds
// Using thread index as seed, hashes map better to the GPU
// Hashed thread index can also be used as seed for the PRNGs
// Caution: do not initialize xorshift with 0 as the sequence stays 0
// Wang hash returns 0 for seed==61

// For normalized random number in [0, 1) use: (float)rng_state * RAND_NORM
// rng_state can be max 0xFFFFFFFF==4294967295 => rand in [0,1)
const float RAND_NORM = (1.0f / 4294967296.0f);

// Xorshift algorithm from George Marsaglia's paper
uint32_t rand_xorshift(uint32_t rng_state) {
  rng_state ^= (rng_state << 13);
  rng_state ^= (rng_state >> 17);
  rng_state ^= (rng_state << 5);
  return rng_state;
}

// LCG (Linear congruential generator) values from Numerical Recipes
uint32_t rand_lcg(uint32_t rng_state) {
  rng_state = 1664525 * rng_state + 1013904223;
  return rng_state;
}

// Hash function by Thomas Wang
// http://www.burtleburtle.net/bob/hash/integer.html
uint32_t wang_hash(uint32_t seed) {
  seed = (seed ^ 61) ^ (seed >> 16);
  seed *= 9;
  seed = seed ^ (seed >> 4);
  seed *= 0x27d4eb2d;
  seed = seed ^ (seed >> 15);
  return seed;
}

// Mersenne Twister
// https://us.fixstars.com/opencl/book/OpenCLProgrammingBook/mersenne-twister/

// Multiply With Carry
// https://www.ast.cam.ac.uk/~stg20/cuda/random/index.html
// used by CUDAMCML, they also have versions for [0,1) and (0,1]

// Quasirandom sequences (blue noise)
// http://extremelearning.com.au/unreasonable-effectiveness-of-quasirandom-sequences/

// Hash by Dave Hoskins
// https://www.shadertoy.com/view/4djSRW

// Random from original mcml (works only on CPU)
#ifdef CL2CPU
#define STANDARDTEST 1
  /* testing program using fixed rnd seed. */
/***********************************************************
 *	A random number generator from Numerical Recipes in C.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9
float ran3(int *idum) {
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  return mj*FAC;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC
/***********************************************************
 *	Generate a random number between 0 and 1.  Take a 
 *	number as seed the first time entering the function.  
 *	The seed is limited to 1<<15.  
 *	We found that when idum is too large, ran3 may return 
 *	numbers beyond 0 and 1.
 ****/
float RandomNum(void) {
  static int first_time=1;
  static int idum;	/* seed for ran3. */
  if(first_time) {
#if STANDARDTEST /* Use fixed seed to test the program. */
    idum = - 1;
#else
    idum = -(int)time(NULL)%(1<<15);
	  /* use 16-bit integer as the seed. */
#endif
    ran3(&idum);
    first_time = 0;
    idum = 1;
  }
  
  return( (float)ran3(&idum) );
}
#undef STANDARDTEST
#endif // CL2CPU