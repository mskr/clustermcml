/*********************************************************************************
*
* This file contains functions related to random number generation.

* Following Nathan Reed's article:
* http://reedbeta.com/blog/quick-and-easy-gpu-random-numbers-in-d3d11/
* PRNGs are designed to go deep, i.e. have good distributions when sequentially updating state
* Hashes are designed to go wide, i.e. have good distributions across initial seeds
* Using thread index as seed, hashes map better to the GPU
* Hashed thread index can also be used as seed for the PRNGs
*
* Caution: do not initialize xorshift with 0 as the sequence stays 0
* Wang hash returns 0 for seed==61
*
*********************************************************************************/

/**
* For normalized random number in [0, 1) use:
* (float)rng_state * RAND_NORM
*/
extern "C" const float RAND_NORM;

/**
* Xorshift algorithm from George Marsaglia's paper
*/
extern "C" uint32_t rand_xorshift(uint32_t rng_state);

/**
* LCG (Linear congruential generator) values from Numerical Recipes
*/
extern "C" uint32_t rand_lcg(uint32_t rng_state);

/**
* Hash function by Thomas Wang
* http://www.burtleburtle.net/bob/hash/integer.html
*/
extern "C" uint32_t wang_hash(uint32_t seed);

/**
* A random number generator from Numerical Recipes in C
*/
extern "C" float ran3(int *idum);

/*
* Random from original mcml (works only on CPU)
*/
extern "C" float RandomNum(void);