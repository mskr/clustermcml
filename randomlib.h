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