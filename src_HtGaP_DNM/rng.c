/*
00002 * Random number generators implemented by Mike Turmon between 1995-1997.
00003 * 
00004 * Normal RNG improved for better speed, 2003.
00005 * 
00006 * Some are LCG's, some are GFSR's.  Also includes a recent 'twisted'
00007 * GFSR, which is not mine.  The TGFSR may not be 64-bit ready.
00008 * 
00009 * Notes: 
00010 *  -- random is the standard unix shift register algorithm, with 256
00011 *  bytes of state (64 int's).
00012 *  -- drand48 is a linear congruential generator (LCG) with multiplier
00013 *  a = 25214903917, offset b = 11, and modulus M = 2^48.
00014 *  -- both minimal generators are LCG's with a = 16807, b = 0, M = 2^31-1.
00015 *  See S.K.Park and K.W.Miller, "Random number generators: good ones are hard
00016 *  to find", Comm. ACM, October 1988, vol 31, number 10, pp 1192-1201, and 
00017 *  many other refs.
00018 *  -- GFSR521 is the 521-bit shift-register generator described in:
00019 *  B.D. Ripley, "Thoughts on pseudorandom number generators," 
00020 *  J. Comput. Appl. Math., 31:153-163, 1990. 
00021 *  -- TGFSR is the very recent `twisted' shift-register of Matsumoto 
00022 *  and Nishimura
00023 * Not here:
00024 *  -- EICG, the particular modulus p=2^31-1 inversive generator described
00025 *  by P. Hellekalek in "Inversive Pseudorandom Number Generators:
00026 *  Concepts, Results, and Links,"  In: C. Alexopoulos et al.,
00027 *  editors, Proceedings of the 1995 Winter Simulation Conference, 
00028 *  pages 255-262, 1995.  The EICG iteration is defined by:
00029 *              y_n = inv(a*(n_0 + n) + b) (mod p)      n >= 0
00030 *  where y=inv(z) means y*z = 1 (mod p), which is unique for prime p.
00031 */

/* Timings on a sun ultrasparc170E, with the sparcworks v4 compiler and
00034 * option -fast, for 10^6 random numbers, are given below:
00035 *  <Generator>                                 <time>    <time/0.26>
00036 *  random from UNIX:                0.70 sec  (2.69)
00037 *  drand48 from UNIX:               2.26 sec  (8.69)
00038 *  LCG of Park/Miller (int version):        0.86 sec  (3.31)
00039 *  LCG of Park/Miller (real version):       2.17 sec  (8.35)
00040 *  GFSR of length-521 of Ripley 1990:       0.26 sec  (1.00)
00041 *  Matsumoto `twisted' GFSR:            0.36 sec  (1.38)
00042 */
 
/*LINTLIBRARY*/

00046 #include <sys/types.h>  /* for getpid() */
00047 #include <sys/time.h>   /* for gettimeofday() */
00048 #include <stdio.h>      /* for fprintf() */
00049 #include <stdlib.h>     /* for calloc() */
00050 #include <unistd.h>     /* for getpid() */
00051 #include <string.h>     /* for memcpy() */
00052 #include <math.h>       /* for M_PI, sin, log */
00053 #include <time.h>       /* clock(), time() */
00054 #include "rng.h"        /* for consistency */
 
00056 #define RANDOM       1 /* random()  random number generator, a SR */
00057 #define DRAND48      2 /* drand48() random number generator, an LCG */
00058 #define MINIMAL_INT  3 /* "minimal" LCG of Park and Miller, integer version */
00059 #define MINIMAL_REAL 4 /* "minimal" LCG of Park and Miller, real version */
00060 #define GFSR521      5 /* 521-bit SR of Ripley 1990 */
00061 #define TGFSR        6 /* `twisted' SR of Matsumoto and Nishimura 1997 */
00062 
00063 #ifndef GENERATOR
00064 #define GENERATOR    GFSR521
00065 #endif
 
00067 /* == GLOBAL DEFINITIONS ============================================== */
00068 
00069 #define GenIs(g) (GENERATOR == g)  /* are we using generator g? */
00070 
00071 #include <strings.h>
00072 #if GenIs(MINIMAL_REAL)
00073 #include <math.h> 
00074 #elif GenIs(RANDOM)
00075 long  random(void);
00076 int   srandom(unsigned);
00077 char *initstate(unsigned, char *, int);
00078 char *setstate(char *);
00079 #elif GenIs(DRAND48)
00080 #include <stdlib.h>
00081 #endif
00082 
00083 /* == UTILITIES ======================================================= */
00084 
00085 static int 
00086 rng_is_bigendian(void)
00087 {
00088   union {
00089     long l;
00090     char c[sizeof(long)];
00091   } u;
00092 
00093   u.l = 1L;
00094   return(u.c[sizeof(long) - 1] == 1);
00095 }
00096 
00097 
00098 /* 
00099  * generate a double nan.  routine is local to this file
00100  */
00101 static double 
00102 rng_getnand(void)
00103 {
00104   /* set the inf/nan buffer up only once */
00105   static int ieee_nan_inited = 0;
00106   /* the union ensures the underlying double is aligned properly */
00107   static union {
00108     double d;
00109     unsigned char c[8];
00110   } nan;
00111 
00112   if (!ieee_nan_inited) {
00113     int i;
00114     for (i = 0; i < 8; i++)
00115       nan.c[i] = 1;  /* have seen 0 here also */
00116     if (rng_is_bigendian()) {
00117       nan.c[0] = 0x7F;
00118       nan.c[1] = 0xF0; /* with 0 above goes 0xF8 here */
00119     } else {
00120       nan.c[7] = 0x7F;
00121       nan.c[6] = 0xF0; /* ditto */
00122     }
00123     ieee_nan_inited = 1; /* flag: buffer is set up now */
00124   }
00125   return(nan.d);
00126 }
00127 
00128 
00129 
00130 /*
00131  * short text description of the current generator
00132  */
00133 char *
00134 rng_name(void)
00135 {
00136 #if GenIs(RANDOM)
00137 return("random from UNIX");
00138 
00139 #elif GenIs(MINIMAL_REAL)
00140 return("LCG of Park/Miller (real version)");
00141 
00142 #elif GenIs(MINIMAL_INT)
00143   return("LCG of Park/Miller (int version)");
00144 
00145 #elif GenIs(DRAND48)
00146   return("drand48 from UNIX");
00147 
00148 #elif GenIs(GFSR521)
00149   return("GFSR of length-521 of Ripley 1990");
00150 
00151 #elif GenIs(TGFSR)
00152   return("Matsumoto `twisted' GFSR");
00153 
00154 #else
00155 #error "empty function"
00156 #endif
00157 }
00158 
00159 /*
00160  * short text description of the current generator
00161  */
00162 int
00163 rng_numeric_name(void)
00164 {
00165   return(GENERATOR);
00166 }
00167 
00168 /* == STATIC STORAGE FOR STATE ======================================== */
00169 
00170 #if GenIs(GFSR521)
00171 #define GFSR_P 521 /* how far back to look for first bit to XOR */
00172 #define GFSR_Q 32  /* distance from first bit to second bit  */
00173 #define GFSR_L 32  /* word length */
00174 static unsigned state[GFSR_P]; /* holds the shift register */
00175 static unsigned *state_dst; /* oldest value */
00176 static unsigned *state_src; /* GFSR_Q in front of oldest value */
00177 
00178 #elif GenIs(RANDOM)
00179 #define RANDOM_REGLENGTH 64
00180 static long state[RANDOM_REGLENGTH];
00181 
00182 #elif GenIs(MINIMAL_REAL)
00183 static double state; /* but this double only holds integers */
00184 
00185 #elif GenIs(MINIMAL_INT)
00186 static int state;
00187 
00188 #elif GenIs(TGFSR)
00189 #define TGFSR_P 624
00190 static unsigned long state[TGFSR_P]; /*  624 words */
00191 
00192 #endif
00193 /* == INITIALIZING THE STATE ============================================ */
00194 
00195 /* special for the GFSR521 */
00196 #if GenIs(GFSR521)
00197 #define MODULUS_I 2147483647  /* 2^31 - 1*/
00198 #define MULTIPLIER_I 16807    /* 7**5 */
00199 #define MODoverMULT_I 127773  /* MODULUS / MULTIPLIER */
00200 #define MODmodMULT_I 2836     /* MODULUS % MULTIPLIER */
00201 
00202 static
00203 void
00204 init_gfsr(int seed)
00205 {
00206   unsigned char x0[2*GFSR_P]; /* seed bits : 0/1 */
00207   int  x0_seed;
00208   int i, j;
00209 
00210   /* set up first GFSR_P entries via linear generator */
00211   for (x0_seed = seed, i = 0; i < GFSR_P; i++) {
00212     x0_seed = 
00213       MULTIPLIER_I * (x0_seed % MODoverMULT_I) - 
00214       MODmodMULT_I * (x0_seed / MODoverMULT_I);
00215     if (x0_seed <= 0)
00216       x0_seed += MODULUS_I;
00217     x0[i] = (x0_seed > (MODULUS_I/2)); /* convert to 0/1 for x0 */
00218   }
00219   /* set up rest of x0 */
00220   for (i = GFSR_P; i < 2*GFSR_P; i++) 
00221     x0[i] = x0[i - GFSR_P] ^ x0[i - GFSR_Q];
00222   /* finally set the state array */
00223   for (i = 0; i < GFSR_P; i++) 
00224     for (state[i] = j = 0; j < GFSR_L; j++)
00225       state[i] |= x0[i+16*(j+1)] << j;
00226   /* initialize the pointers */
00227   state_dst = state + (GFSR_P-GFSR_Q-1);
00228   state_src = state + (GFSR_P-1);
00229 }
00230 #endif
00231 
00232 /* special for the TGFSR */
00233 #if GenIs(TGFSR)
00234 static
00235 void
00236 init_tgfsr(unsigned long seed) /* seed should not be 0 */
00237 {
00238   int k;
00239   
00240   /* setting initial state using     */
00241   /* the generator Line 25 of Table 1 in          */
00242   /* [KNUTH 1981, The Art of Computer Programming */
00243   /*    Vol. 2 (2nd Ed.), pp102]                  */
00244 
00245   state[0]= seed & 0xffffffff;
00246   for (k=1; k<TGFSR_P; k++)
00247     state[k] = (69069 * state[k-1]) & 0xffffffff;
00248 }
00249 #endif
00250 
00251 /*
00252  * Initialize state according to a given seed
00253  */
00254 int
00255 rng_init_state_seeded(int seed)
00256 {
00257 #if GenIs(RANDOM)
00258   initstate((unsigned) seed, (char *) state, sizeof(state));
00259   setstate((char *) state);
00260 
00261 #elif GenIs(MINIMAL_REAL)
00262   state = seed; /* an assignment is enough */
00263 
00264 #elif GenIs(MINIMAL_INT)
00265   state = seed; /* an assignment is enough */
00266 
00267 #elif GenIs(DRAND48)
00268   srand48((long) seed);
00269 
00270 #elif GenIs(GFSR521)
00271   init_gfsr(seed);
00272 
00273 #elif GenIs(TGFSR)
00274   init_tgfsr((unsigned long) seed);
00275 
00276 #else
00277 #error "empty function"
00278 #endif
00279   return(seed);
00280 }
00281 
00282 /* 
00283  * form a positive int to seed the random number generator 
00284  */
00285 
00286 int 
00287 rng_new_seed(void)
00288 {
00289   unsigned int t = 0;
00290   FILE *fp;
00291   struct timeval tv;
00292   int ok = 0;
00293   int nread;
00294 
00295   /* try to read t from /dev/urandom */
00296   if (!ok) {
00297     fp = fopen("/dev/urandom", "r");
00298     if (fp != NULL) {
00299       nread = fread(&t, sizeof(t), (size_t) 1, fp);
00300       fclose(fp);
00301       if (nread == 1)
00302     ok = 1; // success
00303       else
00304     t = 0; // reset to uninitialized
00305     }
00306   }
00307   /* try to open a pipe to the system to get a seed */
00308   if (!ok) {
00309     // md5 is less portable, and no better for our purposes
00310     // plus, it returns a very long checksum
00311     fp = popen("ps auxww | cksum", "r");
00312     if (fp != NULL) {
00313       nread = fscanf(fp, "%u", &t);
00314       pclose(fp);
00315       if (nread == 1)
00316     ok = 1; // success
00317       else
00318     t = 0; // reset to uninitialized
00319     }
00320   }
00321   /* if not, use the old initialization method.
00322    * its seeds can repeat too often, although we have tried
00323    * to patch that up by using the "microseconds" from tv.
00324    */
00325   if (!ok) {
00326     gettimeofday(&tv, NULL);
00327     /* "logic" here:
00328      * tv.sec, tv.usec are time-of-day
00329      * clock() is time-since-invocation
00330      * Low-order bits of tv.sec change from run to run
00331      * Many bits of tv.usec change.
00332      * So, combine them so the changing bits do not overlap:
00333      * sec bits low, usec bits high.  Since 1e6 ~ 2^20, there
00334      * are 20 bits of tv.usec.
00335      * uid is not doing much.
00336      * For pids, make sure if the pid and ppid change by one each,
00337      * the changes do not cancel.
00338      */
00339     t = 
00340       ((tv.tv_sec  & 0xffffffff)  <<  0) ^
00341       ((tv.tv_usec & 0x000fffff)  << 11) ^
00342       ((clock()    & 0x00ffffff)  <<  7) ^
00343       ((getuid()   & 0x000fffff)  << 10) ^
00344       ((getpid()   & 0x000fffff)  <<  9) ^
00345       ((getppid()  & 0x0000ffff)  <<  1);
00346   }
00347   /* at this point, t is defined */
00348   t &= 0x7fffffff;   /* force t < 2^31 */
00349   if (t == 0) t = 1; /* force nonzero */
00350   return(t);
00351 }
00352 
00353 
00354 /*
00355  * Set the seed automatically.  Does not re-seed the rng if
00356  * it is already seeded.  This is useful if this is used as a 
00357  * dynamically loaded library: multiple uses of this routine will
00358  * not reset the seed.  Returns the seed used.
00359  *
00360  * assumes the generated seed is nonzero.
00361  */
00362 int
00363 rng_init_state(void)
00364 {
00365   static int seed = 0; /* 0 is also the sentinel for unseeded */
00366 
00367   /* do not re-seed */
00368   if (seed == 0) {
00369     /* returned seed is guaranteed nonzero */
00370     rng_init_state_seeded(seed = rng_new_seed());
00371   }
00372   return seed;
00373 }
00374 
00375 
00376 /* == SAVING/RESTORING THE STATE ======================================== */
00377 
00378 /*
00379  * number of bytes in the entire state of the generator
00380  */
00381 int
00382 rng_state_length(void)
00383 {
00384   int bufused; /* sizeof the space needed to store the state */
00385 
00386 #if GenIs(RANDOM)
00387   bufused = sizeof(state);
00388 
00389 #elif GenIs(MINIMAL_REAL)
00390   bufused = sizeof(state);
00391 
00392 #elif GenIs(MINIMAL_INT)
00393   bufused = sizeof(state);
00394 
00395 #elif GenIs(DRAND48)
00396   bufused = 3*sizeof(unsigned short);
00397 
00398 #elif GenIs(GFSR521)
00399   bufused = sizeof(state) + sizeof(state_dst) + sizeof(state_src);
00400 
00401 #elif GenIs(TGFSR)
00402   bufused = sizeof(state);
00403 
00404 #else
00405 #error "empty function"
00406 #endif
00407   return(bufused);
00408 }
00409 
00410 /*
00411  * get the entire state of the generator
00412  * inputs: buffer for the state, and its length
00413  * returns size of state, or 0 if failure (state too big for buf)
00414  */
00415 int
00416 rng_get_state(char *buf, int buflen)
00417 {
00418   /* sizeof the space needed to store the seed */
00419   int bufused = rng_state_length();
00420 
00421 #if GenIs(RANDOM)
00422   if (buflen >= bufused)
00423     bcopy((void *) state, (void *) buf, (size_t) bufused);
00424   else
00425     bufused = 0;
00426 
00427 #elif GenIs(MINIMAL_REAL)
00428   if (buflen >= bufused)
00429     bcopy((void *) &state, (void *) buf, (size_t) bufused);
00430   else
00431     bufused = 0;
00432 
00433 #elif GenIs(MINIMAL_INT)
00434   if (buflen >= bufused)
00435     bcopy((void *) &state, (void *) buf, (size_t) bufused);
00436   else
00437     bufused = 0;
00438 
00439 #elif GenIs(DRAND48)
00440   if (buflen >= bufused) {
00441     unsigned short *buf_internal;
00442 
00443     /* get address of internal buffer where state is located */
00444     /* (as side-effect: sets state to contents of buf!) */
00445     buf_internal = seed48((unsigned short *) buf); 
00446     /* copy the internal buffer to our own buffer */
00447     bcopy((void *) buf_internal, (void *) buf, bufused);
00448     /* restore the old state, undoing side-effect above */
00449     seed48((unsigned short *) buf); 
00450   }
00451   else
00452     bufused = 0;
00453 
00454 #elif GenIs(GFSR521)
00455   if (buflen >= bufused) {
00456     /*  save: state, state_dst, state_src */
00457     bcopy((void *) state, 
00458       (void *) buf, 
00459       (size_t) sizeof(state));
00460     bcopy((void *) &state_dst, 
00461       (void *) (buf + sizeof(state)), 
00462       (size_t) sizeof(state_dst));
00463     bcopy((void *) &state_src, 
00464       (void *) (buf + sizeof(state) + sizeof(state_dst)), 
00465       (size_t) sizeof(state_src));
00466   }
00467   else
00468     bufused = 0;
00469 
00470 #elif GenIs(TGFSR)
00471   if (buflen >= bufused)
00472     bcopy((void *) state, (void *) buf, (size_t) bufused);
00473   else
00474     bufused = 0;
00475 
00476 #else
00477 #error "empty function"
00478 #endif
00479   return(bufused);
00480 }
00481 
00482 /*
00483  * get the entire state to buf
00484  * inputs: buffer holding the state, and its length
00485  * returns size of state, or 0 if failure (buf too small for state)
00486  */
00487 int
00488 rng_set_state(char *buf, int buflen)
00489 {
00490   int bufneed = rng_state_length();
00491 
00492 #if GenIs(RANDOM)
00493   if (buflen >= bufneed) {
00494     bcopy((void *) buf, (void *) state, (size_t) bufneed);
00495     setstate((char *) state); /* tell random() to use this buffer */
00496   }
00497   else
00498     bufneed = 0;
00499 
00500 #elif GenIs(MINIMAL_REAL)
00501   if (buflen >= bufneed)
00502     bcopy((void *) buf, (void *) &state, (size_t) bufneed);
00503   else
00504     bufneed = 0;
00505 
00506 #elif GenIs(MINIMAL_INT)
00507   if (buflen >= bufneed)
00508     bcopy((void *) buf, (void *) &state, (size_t) bufneed);
00509   else
00510     bufneed = 0;
00511 
00512 #elif GenIs(DRAND48)
00513   if (buflen >= bufneed) {
00514     /* copy our own buffer to the internal buffer */
00515     seed48((unsigned short *) buf); 
00516   }
00517   else
00518     bufneed = 0;
00519 
00520 #elif GenIs(GFSR521)
00521   if (buflen >= bufneed) {
00522     /*  restore: state, state_dst, state_src */
00523     bcopy((void *) buf, (void *) state, (size_t) sizeof(state));
00524     bcopy((void *) (buf + sizeof(state)), 
00525       (void *) &state_dst, 
00526       (size_t) sizeof(state_dst));
00527     bcopy((void *) (buf + sizeof(state) + sizeof(state_dst)), 
00528       (void *) &state_src, 
00529       (size_t) sizeof(state_src));
00530   }
00531   else
00532     bufneed = 0;
00533 
00534 #elif GenIs(TGFSR)
00535   if (buflen >= bufneed)
00536     bcopy((void *) buf, (void *) state, (size_t) bufneed);
00537   else
00538     bufneed = 0;
00539 
00540 #else
00541 #error "empty function"
00542 #endif
00543   return(bufneed);
00544 }
00545 
00546 /* == GENERATING RANDOM NUMBERS ======================================== */
00547 
00548 #if GenIs(RANDOM)
00549 #define MODULUS 2147483648.0  /* 2^31 */
00550 #define MODULUS_INV 4.656612873077393e-10 /* 1/MODULUS */
00551 
00552 #elif GenIs(MINIMAL_REAL)
00553 #define MODULUS 2147483647.0  /* 2^31 - 1*/
00554 #define MODULUS_INV 4.656612875245797e-10 /* 1/MODULUS */
00555 #define MULTIPLIER 16807.0 /* 7**5 */
00556 
00557 #elif GenIs(MINIMAL_INT)
00558 #define MODULUS 2147483647  /* 2^31 - 1*/
00559 #define MODULUS_INV 4.656612875245797e-10 /* 1/MODULUS */
00560 #define MULTIPLIER 16807 /* 7**5 */
00561 #define MODoverMULT 127773 /* MODULUS / MULTIPLIER */
00562 #define MODmodMULT 2836 /* MODULUS % MULTIPLIER */
00563 
00564 #elif GenIs(GFSR521)
00565 #define MODULUS 4294967296 /* 2^32 */
00566 #define MODULUS_INV 2.3283064365387e-10 /* 1/MODULUS */
00567 
00568 #elif GenIs(TGFSR)
00569 /* Period parameters */  
00570 #define N_twist 624
00571 #define M_twist 397
00572 #define MATRIX_A 0x9908b0df   /* constant vector a */
00573 #define UPPER_MASK 0x80000000 /* most significant w-r bits */
00574 #define LOWER_MASK 0x7fffffff /* least significant r bits */
00575 /* for tempering */   
00576 #define TEMPERING_MASK_B 0x9d2c5680
00577 #define TEMPERING_MASK_C 0xefc60000
00578 #define TEMPERING_SHIFT_U(y)  (y >> 11)
00579 #define TEMPERING_SHIFT_S(y)  (y << 7)
00580 #define TEMPERING_SHIFT_T(y)  (y << 15)
00581 #define TEMPERING_SHIFT_L(y)  (y >> 18)
00582 
00583 #endif
00584 
00585 /*
00586  * uniform number in (0,1)
00587  */
00588 double 
00589 rng_uniform (void)
00590 {
00591 #if GenIs(RANDOM)
00592   return(((double) random()) * MODULUS_INV);
00593 
00594 #elif GenIs(MINIMAL_REAL)
00595   return(
00596      (state = fmod(MULTIPLIER * state, MODULUS)) * MODULUS_INV
00597      ); 
00598 #elif GenIs(MINIMAL_INT)
00599   state = 
00600     MULTIPLIER * (state % MODoverMULT) - 
00601     MODmodMULT * (state / MODoverMULT);
00602   if (state <= 0)
00603     state += MODULUS;
00604   return(state * MODULUS_INV); 
00605 
00606 #elif GenIs(DRAND48)
00607   return(drand48());
00608 
00609 #elif GenIs(GFSR521)
00610   double rval; 
00611 
00612   rval = *state_dst * MODULUS_INV; /* use this value */
00613   *state_dst ^= *state_src; /* update the shift register */
00614   state_dst = (state_dst > state) ? state_dst-1 : state+(GFSR_P-1);
00615   state_src = (state_src > state) ? state_src-1 : state+(GFSR_P-1);
00616   return(rval);
00617 
00618 #elif GenIs(TGFSR)
00619   unsigned long y;
00620   static int k = 1;
00621   static unsigned long mag01[2]={0x0, MATRIX_A};
00622   /* mag01[x] = x * MATRIX_A  for x=0,1 */
00623   
00624   if(k == N_twist){ /* generate N_twist words at one time */
00625     int kk;
00626     for (kk=0;kk<N_twist-M_twist;kk++) {
00627       y = (state[kk]&UPPER_MASK)|(state[kk+1]&LOWER_MASK);
00628       state[kk] = state[kk+M_twist] ^ (y >> 1) ^ mag01[y & 0x1];
00629     }
00630     for (;kk<N_twist-1;kk++) {
00631       y = (state[kk]&UPPER_MASK)|(state[kk+1]&LOWER_MASK);
00632       state[kk] = state[kk+(M_twist-N_twist)] ^ (y >> 1) ^ mag01[y & 0x1];
00633     }
00634     y = (state[N_twist-1]&UPPER_MASK)|(state[0]&LOWER_MASK);
00635     state[N_twist-1] = state[M_twist-1] ^ (y >> 1) ^ mag01[y & 0x1];
00636     
00637     k = 0;
00638   }
00639   
00640   y = state[k++];
00641   y ^= TEMPERING_SHIFT_U(y);
00642   y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
00643   y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
00644   y &= 0xffffffff; /* you may delete this line if word size = 32 */
00645   y ^= TEMPERING_SHIFT_L(y);
00646 
00647   return ( (double)y * (1.0/ (unsigned long)0xffffffff) );
00648 
00649 #else
00650 #error "empty function"
00651 #endif
00652 }
00653 
00654 /******************************************************************
00655  *
00656  * Normal random numbers
00657  *
00658  ******************************************************************/
00659 
00660 /*
00661  * normal, mean 0, unit variance
00662  * 
00663  * formerly just:
00664  *  return(sqrt(-2 * log(rng_uniform())) * sin(2 * M_PI * rng_uniform()));
00665  * as in rng_normal_old().
00666  * In this version, we retain the "orthogonal" random number as a double
00667  * called ace_in_the_hole, and indicate its presence (on alternate calls
00668  * to this routine) via have_ace.  
00669  * The method used eliminates the sin, cos calls of the "pure" bivariate 
00670  * transformation by starting from a uniform tuple (u1,u2) that lives 
00671  * on the unit disk.  The random angle that used to be generated by 
00672  * sin(2 pi random()) is now just encoded in the coordinates of 
00673  * the uniform tuple, scaled by the square root of its magnitude, 
00674  * i.e. u1/sqrt(norm) and u2/sqrt(norm).  
00675  * And,  conditioned to be < 1, norm := u1^2 + u2^2 is U(0,1).  So norm 
00676  * is used for the length of the bivariate normal, which is generated by:
00677  *   n1 = sqrt(-2 * log(norm)) * (u1/sqrt(norm))
00678  *   n2 = sqrt(-2 * log(norm)) * (u2/sqrt(norm))
00679  * Tucking the sqrt(norm) into the first factor yields the code below.
00680  * It runs about 2.5 times faster than the old version.
00681  */
00682 double
00683 rng_normal(void) {
00684   static double ace_in_the_hole; /* a N(0,1) side-effect of the prior call */
00685   static int have_ace = 0;       /* 1 if we have the above, 0 if not */
00686   double u1, u2;                 /* a pair of uniform variables */
00687   double norm;                   /* their norm taken as a tuple */
00688 
00689   if (have_ace) {
00690     have_ace = 0;                /* we're about to use it up */
00691     return ace_in_the_hole;
00692   } else {
00693     do {
00694       /* get a pair (u1,u2) that is uniform on the unit disk */
00695       u1 = 2.0*rng_uniform() - 1.0;  /* U([-1,1)) */
00696       u2 = 2.0*rng_uniform() - 1.0;  /* U([-1,1)) */
00697     } while (1.0 <= (norm = u1*u1 + u2*u2));
00698     norm = sqrt(-2.0 * log(norm)/norm); /* just a scale factor, see above */
00699     ace_in_the_hole = u1 * norm; /* first N(0,1): save for later */
00700     have_ace = 1;                /* remember we saved it */
00701     return u2 * norm;            /* second N(0,1): return it now */
00702   }
00703 }
00704 
00705 /*
00706  * normal, mean 0, unit variance
00707  */
00708 double
00709 rng_normal_old(void) {
00710     return (sqrt (-2.0 * log(rng_uniform())) *
00711         sin(2.0 * M_PI * rng_uniform()));
00712 }
00713 
00714 
00715 /*
00716  * normal, mean mu, standard deviation sigma
00717  */
00718 double
00719 rng_normal_sdev(double mu, double sigma) {
00720   return rng_normal()*sigma + mu;
00721 }
00722 
00723 /*
00724  * normal, mean mu, variance sigma2
00725  * if rng_normal were more clever, we could eliminate the sqrt()
00726  */
00727 double
00728 rng_normal_var(double mu, double sigma2) {
00729   return rng_normal()*sqrt(sigma2) + mu;
00730 }
00731 
00732 /******************************************************************
00733  *
00734  * Non-normal real-valued random numbers
00735  *
00736  ******************************************************************/
00737 
00738 /*
00739  * rng_gamma(alpha,theta): Returns a draw from gamma(alpha,theta).
00740  * 
00741  * The corresponding density is:
00742  * 
00743  *                       x^(alpha - 1) exp(-x/theta)
00744  * p(x; alpha, theta) =  ---------------------------
00745  *                       Gamma(alpha) theta^alpha
00746  * 
00747  * In this parameterization, theta is purely a scale parameter 
00748  * (similar to sigma in the Normal distribution) and alpha is 
00749  * a shape parameter.
00750  * 
00751  * For alpha <= 1, the distribution is one-sided and cusped like 
00752  * an exponential, and for alpha > 1, the distribution is (increasingly) 
00753  * two-sided and rounded like a gaussian.
00754  * 
00755  * If X ~ gamma(alpha, theta), then
00756  * 
00757  *   E X   = alpha theta
00758  *   var X = alpha theta^2
00759  * 
00760  * The algorithm below is primarily comprised of two rejection-based 
00761  * methods, recommended by Devroye, one for alpha < 1 and one for
00762  * alpha > 1.  (We use the exponential for alpha = 1.)  Both are 
00763  * crafted to have bounded (and small) rejection proportion for 
00764  * all values of alpha and of the variate x.  
00765  * I further separate out the alpha = 2 case because
00766  * it is easy and relatively common.
00767  * Regarding the alpha = 1 case, Devroye's book seems to indicate 
00768  * his version of the Best (1978) "XG" algorithm works for alpha = 1.
00769  * It has a problem there, as his web errata point out.  This is 
00770  * another reason to special-case the exponential.
00771  * 
00772  * Some of the acceptance crafting seems more useful for older hardware. 
00773  * For example, extra tests have been added for early rejection to avoid
00774  * computation of exponentials.  The exponentials give the exact rejection
00775  * criterion, but the extra tests are quicker to compute.  I have not 
00776  * checked the impact of tweaking the algorithms.
00777  * 
00778  * References: 
00779  *   Luc Devroye, Non-Uniform Random Variate Generation, Springer, 
00780  *     1986, chap. IX.3-IX.6.
00781  *   http://cg.scs.carleton.ca/~luc/rnbookindex.html
00782  *   http://cg.scs.carleton.ca/~luc/errors.pdf
00783  * 
00784  * which cites:
00785  * 
00786  *   Best 1978: D. J. Best, "A simple algorithm for the computer 
00787  *   generation of
00788  *   random samples from a Student's t or symmetric beta distribution," in
00789  *   Compstat 1978: Proceedings in Computational statistics, ed. L.C.A. 
00790  *   Corsten and J. Hermans, pp. 341-7, Physica Verlag, Vienna, 1978.
00791  *   => a.k.a. "XG" algorithm for alpha > 1.
00792  *   See Devroye, page 410.
00793  * 
00794  * and
00795  * 
00796  *   Best 1983: D. J. Best, "A note on gamma variate generators with shape
00797  *   parameter less than unity," Computing, vol. 30, pp. 185-8, 1983.
00798  *   => a.k.a. "XGS" algorithm for alpha < 1.
00799  *   See Devroye's exercise 6 on page 426.
00800  * 
00801  * I tested this algorithm for proper moments (1 thru 4), and 
00802  * distributional shape, for many values of alpha (excercising all 
00803  * cases below).
00804  *  
00805  */
00806 double 
00807 rng_gamma(double alpha, double theta)
00808 {
00809   // follow DeVroye's notation
00810   double X, U, V, W, Y, Z;  // random variates
00811   double b, c, t;           // parameters
00812 
00813   if (alpha <= 0 || theta < 0)
00814     return rng_getnand();
00815   if (alpha == 1.0) {
00816     // Use the exponential, e.g. Devroye p. 405
00817     X = -log(rng_uniform());
00818   } else if (alpha == 2.0) {
00819     // Use the sum of 2 exponentials (Devroye p. 405)
00820     X = -log(rng_uniform() * rng_uniform());
00821   } else if (alpha < 1.0) {
00822     // Use RGS algorithm (Best, 1983; Devroye, p. 426)
00823     c = 1.0 / alpha; 
00824     t = 0.07 + 0.75 * sqrt(1.0 - alpha); // breakpoint: min reject prob
00825     b = 1.0 + exp(-t) * alpha / t;
00826     // generate a variate x
00827     while (1) {
00828       U = rng_uniform(); 
00829       W = rng_uniform(); 
00830       V = b * U;
00831       if (V <= 1.0) {
00832     X = t * pow(V, c);
00833     /* the first step is an early bailout because
00834      * exp(-x) >= (2-x)/(2+x) for x>=0
00835      */
00836     // test is equivalent to w <= (2-x)/(2+x) since x>=0 
00837     if (W * (2.0 + X) <= (2.0 - X)) break;
00838     if (W <= exp(-X)) break;
00839       } else {
00840     X = -log(c * t * (b - V));
00841     Y = X / t;
00842     if ((W * (alpha + Y - alpha*Y)) <= 1.0) break;
00843     if (W <= pow(Y, alpha - 1.0)) break;
00844       }
00845     }
00846   } else {
00847     // Use rejection algorithm XG (Best, 1978; Devroye, p. 410)
00848     //   (This algorithm is OK for alpha > 1.  If alpha = 1,
00849     //   b = 0, and the last rejection step with log(X/b) is
00850     //   illegal.  This is noted in Devroye's errata on the web.)
00851     b = alpha - 1.0;
00852     c = 3.0 * alpha - 0.75;
00853     // generate a variate x
00854     while (1) {
00855       U = rng_uniform();  
00856       V = rng_uniform();
00857       W = U * (1.0 - U);  
00858       Y = sqrt(c / W) * (U - 0.5);
00859       X = b + Y;
00860       if (X >= 0.0) {
00861     Z = 64.0 * W*W*W * V*V;
00862     /* the first step is an early bailout which
00863      * Devroye remarks can be omitted at little cost.
00864      */
00865     // test is same as Z <= 1 - 2*Y*Y/X since X >= 0 here
00866     if ((Z * X) <= (X - 2.0 * Y * Y)) break;
00867     if (log(Z) <= (2.0 * (b * log(X/b) - Y))) break;
00868       }
00869     }
00870   }
00871   return X * theta;
00872 }
00873 
00874 /*
00875  * beta(a,b) sample: Just use the gamma generator above.
00876  * Based on the remarks of Devroye chap. IX.4, pp 432-433:
00877  * "Roughly speaking, we will be able to improve
00878  * over this generator by at most 50%."
00879  */
00880 double 
00881 rng_beta(double a, double b)
00882 {
00883   double Ga = rng_gamma(a, 1.0);
00884   double Gb = rng_gamma(b, 1.0);
00885   return Ga / (Ga + Gb);
00886   // note: Ga+Gb is independent of the returned value, and is
00887   // distributed as gamma(a+b).
00888 }
00889 
00890 /*
00891  * chi-square(r) sample: Use the gamma generator above.
00892  * See Devroye, chap. IX.3, pg 403.
00893  */
00894 double 
00895 rng_chi2(double r)
00896 {
00897   return rng_gamma(r * 0.5, 2.0);
00898 }
00899 
00900 
00901 /*
00902  * exponential(1) sample: Simple inversion method.
00903  * If you want a different rate, just multiply.
00904  */
00905 double 
00906 rng_exponential()
00907 {
00908   return -log(rng_uniform());
00909 }
00910 
00911 
00912 
00913 /******************************************************************
00914  *
00915  * Discrete-valued random numbers
00916  *
00917  ******************************************************************/
00918 
00919 /******************************************************************
00920  * Poisson
00921  */
00922 
00923 /*
00924  * return a Poisson r.v., using an inter-arrival time method
00925  * I.e., a Poisson process of unit intensity has exponential(1)
00926  * inter-arrival times, thus the number of arrivals between
00927  * time = 0 and time = lambda is Poisson(lambda), thus generate
00928  * Poisson R.V.'s by seeing how many exponential(1) r.v.'s it
00929  * takes until the sum exceeds lambda.  Since an exponential(1)
00930  * is just -log(uniform), this is the same as seeing when:
00931  *   exp(Sum_k [ -log(U(k)) ] ) > lambda, 
00932  * i.e. when 
00933  *   Prod_k(U(k)) < lambda.
00934  * See also Devroye, chapter X.3, Lemma 3.3.
00935  */
00936 static
00937 double
00938 rng_poisson_count(double lambda)
00939 {
00940   const double gate = exp(-lambda);
00941   double prod, x;
00942 
00943   if (gate == 1) {
00944     // e.g., lambda = 1e-20 => exp(-lambda) == 1 for IEEE doubles
00945     // The below line will always return 0 with most rng_uniform
00946     // implementations, because they don't output numbers less than
00947     // floating point epsilon.
00948     return (rng_uniform() < lambda) ? 1.0 : 0.0;
00949   }
00950   // it is guaranteed that rate < 1, so at exit, x is at least 1
00951   for (prod = 1.0, x = 0.0; prod > gate; prod *= rng_uniform(), x++)
00952     ;
00953   return x-1.0; // the previous "prod" was the last one > gate
00954 }
00955 
00956 /*
00957  * transfomed rejection algorithm, due to W. H?rmann.
00958  *
00959  * This is algorithm PTRS of the paper:
00960  *   W. H?rmann, 
00961  *   "The transformed rejection method for generating Poisson random 
00962  *   variables," Insurance: Mathematics and Economics 12, 39-45 (1993).
00963  *
00964  * See also:
00965  *   http://statmath.wu.ac.at/staff/hoermann/publications.html
00966  *
00967  * This generator is faster than that of Devroye, less complex,
00968  * and does not seem to suffer from unpleasant bugs!  Checked
00969  * 2/2010 against the pmf for several rates mu.
00970  *
00971  * Valid only for rate (mu) >= 10.
00972  */
00973 static
00974 double
00975 rng_poisson_trans_reject(double mu)
00976 {
00977   double b, a, v_r;
00978   double one_alpha;  // "(1/alpha)" in the paper
00979   double U, US, V;
00980   double k;
00981   
00982   // find parameters
00983   b = 0.931 + 2.53 * sqrt(mu);
00984   a = -0.059 + 0.02483 * b;
00985   v_r = 0.9277 - 3.6224 / (b - 2.0);
00986   // loop until accept
00987   while (1) {
00988     U = rng_uniform();
00989     V = rng_uniform();
00990     U = U - 0.5;
00991     US = 0.5 - fabs(U);
00992     k = floor((2.0*a/US + b) * U + mu + 0.43);
00993     if (k < 0)
00994       continue;
00995     if (US >= 0.07 && V <= v_r)
00996       break;
00997     if (US < 0.013 && V > US)
00998       continue;
00999     // further setup
01000     one_alpha = 1.1239 + 1.1328 / (b - 3.4);
01001     // acceptance-rejection test
01002     if (log(V * one_alpha / (a/(US*US) + b))
01003     <= (-mu + k*log(mu) - lgamma(k+1.0)))
01004       break;
01005   }
01006   return k;
01007 }
01008 
01009 
01010 double
01011 rng_poisson(double lambda)
01012 {
01013   // 30 is the approximate switchover point (intel core 2, 2010)
01014   const double mode_switch = 30;
01015 
01016   if (lambda <= 0.0)
01017     return rng_getnand(); // illegal
01018   else if (lambda < mode_switch)
01019     // fast for small lambda
01020     //   note, exp(-mode_switch) must comfortably fit in a double
01021     return rng_poisson_count(lambda);
01022   else
01023     // longer setup time, but better for large lambda
01024     return rng_poisson_trans_reject(lambda);
01025 }
01026 
01027 /******************************************************************
01028  * Geometric
01029  */
01030 
01031 /*
01032  * geometric(p) sample: Simple inversion method.
01033  * See Devroye chap. X.2.
01034  */
01035 double 
01036 rng_geometric(double p)
01037 {
01038   if (p < 0.25)
01039     // avoid loss of precision in computing log(1-p) 
01040     // when p is small by computing:
01041     //    log1p(-p) = log(1 + (-p)) = log(1-p)
01042     return ceil(log(rng_uniform()) / log1p(-p));
01043   else
01044     return ceil(log(rng_uniform()) / log(1.0 - p));
01045 }
01046 
01047 
01048 /******************************************************************
01049  * Binomial
01050  */
01051 
01052 /*
01053  * Waiting time algorithm.  Works for small n*p.
01054  * Idea is that a geometric(p) r.v. is the number of 
01055  * Bernoulli trials up to and including the first success.
01056  * After the total of trials exceeds n, we count the number 
01057  * of summands.  This is one more (because it exceeded n) than
01058  * the number of 1's in the corresponding Bernoulli sequence.
01059  * The latter is the binomial(p) draw.
01060  * See Devroye, "First waiting time algorithm", p. 525.
01061  */
01062 
01063 static
01064 double
01065 rng_binomial_wait(double n, double p)
01066 {
01067   double sum, x;
01068 
01069   // take care of endpoints now, for robustness
01070   if (n == 0.0)
01071     return 0.0;
01072   else if (p == 0.0)
01073     return 0.0;
01074   else if (p == 1.0)
01075     return n;
01076   // wait for geometric
01077   for (sum = 0, x = 0; sum <= n; sum += rng_geometric(p), x++)
01078     ;
01079   return x - 1.0;
01080 }
01081 
01082 /*
01083  * This is algorithm BTRS from H?rmann.
01084  * It requires n*p > 10.
01085  * This generator is faster than that of Devroye, less complex,
01086  * and (once again) does not seem to suffer from unpleasant bugs!  
01087  * Checked extensively against binomial probabilities from
01088  * Matlab, 2/2010.
01089  *
01090  * W. H?rmann, :The generation of binomial random variates,"
01091  * Journal of Statistical Computation and Simulation 46, 
01092  * 101-110 (1993), online in TR form at:
01093  * http://statmath.wu.ac.at/staff/hoermann/publications.html
01094  *
01095  * The tech report version has a missing log() in step
01096  * 3.1 of the algorithm; noted below.  The published paper
01097  * may correct this.  One can also determine this rather blatant
01098  * error by comparison with the more elaborate BTRD algorithm,
01099  * which uses the same v <- log(v) construction in step 3.2, 
01100  * before the final acceptance test in step 3.4 in which the 
01101  * new "v" is compared against the log-probability.  Thus, the
01102  * log needs to be there.
01103  */
01104 static
01105 double
01106 rng_binomial_trans_reject(double n, double p)
01107 {
01108   const double spq = sqrt(n*p*(1-p));
01109   double b, a, c, v_r;
01110   double alpha, lpq, m, h;
01111   double U, US, V;
01112   double k;
01113   
01114   // find parameters
01115   b = 1.15 + 2.53 * spq;
01116   a = -0.0873 + 0.0248 * b + 0.01 * p;
01117   c = n * p + 0.5;
01118   v_r = 0.92 - 4.2 / b;
01119   // loop until accept
01120   while (1) {
01121     U = rng_uniform();
01122     V = rng_uniform();
01123     U = U - 0.5;
01124     US = 0.5 - fabs(U);
01125     k = floor((2.0 * a/US + b) * U + c);
01126     if (k < 0 || k > n)
01127       continue;
01128     if (US >= 0.07 && V <= v_r)
01129       break;
01130     // further setup
01131     alpha = (2.83 + 5.1/b) * spq;
01132     lpq = log(p/(1-p));
01133     m = floor(n*p+p);
01134     h = lgamma(m+1) + lgamma(n-m+1);
01135     // acceptance-rejection test
01136     V = log(V * alpha/(a / (US*US) + b)); // log() missing from TR
01137     if (V <= h - lgamma(k+1) -lgamma(n-k+1) + (k-m)*lpq)
01138       break;
01139   }
01140   return k;
01141 }
01142 
01143 
01144 /*
01145  * final binomial random variable entry point.
01146  * 
01147  * Calls one of two RNG's, a simple one for small n*p, and
01148  * a more sophisticated one with larger setup time for
01149  * larger n*p.
01150  */
01151 
01152 double
01153 rng_binomial(double n, double p)
01154 {
01155   // 16 is the approximate switchover point (intel core 2, 2010)
01156   // (must be > 10 anyway to accommodate restrictions on fancy RNG)
01157   const double mode_switch = 16;
01158   int flip;
01159   double X;    // random draw
01160 
01161   // take care of errors now
01162   //   (tests written to yield error for nans and infinites)
01163   if (!(n >= 0.0) || !(floor(n) == n) || !finite(n))
01164     return rng_getnand();
01165   if (!(p >= 0.0 && p <= 1.0))
01166     return rng_getnand();
01167   // take care of endpoints to avoid trouble later
01168   if (n == 0.0)
01169     return 0.0;
01170   else if (p == 0.0)
01171     return 0.0;
01172   else if (p == 1.0)
01173     return n;
01174   // if p > 1/2, flip to: n - bin(n, 1-p) 
01175   if (p > 0.5) {
01176     p = 1.0 - p;
01177     flip = 1;
01178   } else {
01179     flip = 0;
01180   }
01181   // generate bin(n, p)
01182   if ((n*p) < mode_switch) {
01183     // fast for small rates
01184     X = rng_binomial_wait(n, p);
01185   } else {
01186     // np is large: transformed rejection method
01187     X = rng_binomial_trans_reject(n, p);
01188   }
01189   // account for flip if need be
01190   return flip ? (n - X) : X;
01191 }
01192 
01193 /******************************************************************
01194  *
01195  * Linear algebra tools
01196  *
01197  ******************************************************************/
01198 
01199 /* 
01200  * Cholesky decomposition (in-place).
01201  * No pivoting.  Works only on strictly positive-definite matrices 
01202  * for this reason.  
01203  * Checked against Matlab chol and found agreement to floating-point
01204  * precision, 4/2006.
01205  * We write  R = G' G, where G is the Cholesky decomposition.
01206  * References only the diagonal and the "upper triangle" of the input R, 
01207  * and over-writes it with the cholesky factor G.  The upper half of 
01208  * the input matrix is left alone.
01209  * For 0 <= j <= i <= N-1, mat[i*N+j] contains the Cholesky factor.  If we 
01210  * identify G(i,j) with mat[i*d+j], the "lower triangle" corresponds to 
01211  * G(i,j) for j <= i, and the other entries in G are zero.
01212  * About shape:
01213  *   The C Cholesky decomposition below produces a triangular
01214  *   Cholesky factor.  It is laid out in memory the same way
01215  *   as a Cholesky factor from matlab's "chol(R)".  
01216  *
01217  * tested versus matlab chol, turmon, 6/2006; again 2/2010.
01218  * (randomized inputs, size 2x2, 5x5, 12x12, 128x128, 256x256)
01219  *
01220  * 0 is returned on failure (ie, non pd matrix), 1 for success.
01221  */
01222 
01223 static
01224 int
01225 rng_cholesky(double *a, 
01226               int N)
01227 {
01228   double fac;
01229   int i, j, k;
01230 
01231   for (k = 0; k < N; k++) {
01232     if (a[k*N+k] <= 0.0) 
01233       return 0;
01234     fac = a[k*N+k] = sqrt(a[k*N+k]);
01235     fac = 1/fac;
01236     for (j = k+1; j < N; j++)
01237       a[j*N+k] *= fac;
01238     for (j = k+1; j < N; j++)
01239       for (i = j; i < N; i++)
01240         a[i*N+j] -= a[i*N+k] * a[j*N+k];
01241   }
01242   return 1;
01243 }
01244 
01245 /*
01246  * Find the inverse of a nonsingular upper-triangular nXn matrix a.
01247  * The computation is done in place, so the inverse replaces a.
01248  * The lower triangle is not referenced or used.
01249  * If a does not have non-zero diagonals, the inverse does not
01250  * exist, and 0 is returned.  Otherwise, 1 is returned.
01251  *
01252  * This is included here because it would help with an
01253  * inverse Wishart, if we ever need it.
01254  *
01255  * tested versus matlab inv, turmon, 6/2006 
01256  * (randomized inputs, size 2x2, 5x5, 12x12, 128x128, 256x256)
01257  *
01258  */
01259 
01260 // not currently used!
01261 #ifdef NOT_DEFINED
01262 
01263 static
01264 int
01265 rng_invert_upper(double *a,
01266          int N)
01267 {
01268   int i,j,k;
01269   double sum;
01270 
01271   for (j = 0; j < N; j++) {
01272     if (a[j*N+j] == 0.0)
01273       return 0; /* failure */
01274     a[j*N+j] = 1 / a[j*N+j];
01275     for (i = j+1; i < N; i++) {
01276       sum = 0.0;
01277       for (k = j; k < i; k++) {
01278         sum -= a[i*N+k] * a[k*N+j];
01279       }
01280       a[i*N+j] = sum / a[i*N+i];
01281     }
01282   }
01283   return 1;
01284 }
01285 
01286 #endif
01287 
01288 
01289 /* multiply c = a * b, with all square upper-triangular matrices,
01290  * and all elements stored in standard matlab ordering
01291  */
01292 static
01293 void 
01294 rng_mult_upper_triangles(double* c, double* a, double* b, int N)
01295 {
01296   int i, j, k;
01297 
01298   for (i = 0; i < N; i++)
01299     for (j = i; j < N; j++) {
01300       // zero out c(j,i) (lower triangle, or diagonal)
01301       c[j+i*N] = 0;
01302       // add up c(i,j) (j>=i)
01303       c[i+j*N] = 0;
01304       for (k = i; k <= j; k++)
01305     // Computing:
01306         //   c(i,j) += a(i,k) * b(k,j)
01307     // Note:
01308     //   a(i,k) = 0 if k < i, giving lower loop bound
01309     //   b(k,j) = 0 if k > j, giving upper loop bound
01310         c[i+j*N] += a[i+k*N] * b[k+j*N];
01311     }
01312 }
01313 
01314 
01315 /* multiply b = a * a', with a square and upper triangular,
01316  * and all elements stored in standard matlab ordering.
01317  * The result is a full matrix.
01318  */
01319 static
01320 void 
01321 rng_square_upper_triangle(double* b, double* a, int N)
01322 {
01323   int i, j, k;
01324 
01325   for (i = 0; i < N; i++)
01326     for (j = i; j < N; j++) {
01327       // sum up b(i,j) = b(j,i) (note j >= i)
01328       b[i+j*N] = 0;
01329       for (k = 0; k <= i; k++)
01330     // Note:
01331     //   a(k,i) = 0 if k > i, and
01332     //   a(k,j) = 0 if k > j, giving k-loop lower bound: min(i,j)
01333     // since j >= i, we start at i.
01334         b[i+j*N] += a[k+i*N] * a[k+j*N];
01335       // reflect into b(j,i)
01336       b[j+i*N] = b[i+j*N];
01337     }
01338 }
01339 
01340 /*************************************************************
01341  *
01342  * Normal random vectors
01343  *
01344  *************************************************************/
01345 
01346 /*
01347  * normal d-vector, mean mu (d), covariance sigma (dxd)
01348  *
01349  * If mu is supplied as NULL, assume mu is all zeros.
01350  * 
01351  * turmon 4/2006: tested for d=1 and d=6 with dense random covariances.
01352  * First, second (full covariance), and fourth (elementwise) moments check
01353  * out, and show correct convergence to their specified values as n grows.
01354  * 3D scatter plots (of the 6D vectors) look normal.  1-d histograms
01355  * also look normal.  The Kolmogorov-Smirnov distributional test along
01356  * each of 6 dimensions shows close adherence to normality (p-values
01357  * around 0.1-0.5 at n=1e6).
01358  *
01359  * Regarding strides.  stridex1 is the distance in doubles between
01360  * adjacent vector samples (typically equals d).  stridex2 is the
01361  * distance in doubles between adjacent elements within one vector
01362  * (typically equals 1).  stridemu is the same, between elements
01363  * of mu (typically equals 1).
01364  */
01365 int 
01366 rng_normal_vector(
01367           double* x,     /* samples */
01368           double* mu,    /* mean, taken to be 0.0 if NULL */
01369           double* sigma, /* covariance */
01370           int n,         /* number of samples */
01371           int d,         /* dims of mean, samples */
01372           int stridex1,  /* stride for each sample */
01373           int stridex2,  /* stride within each sample */
01374           int stridemu   /* stride for elements in mu */
01375           )
01376 {
01377   int i, j, k; /* counters */
01378   double *z;   /* temp space for N(0,1)'s */
01379   double *R;   /* temp space for cholesky factor */
01380 
01381   /* allocate temporary storage */
01382   R = calloc(d*d, sizeof(double));
01383   z = calloc(d  , sizeof(double));
01384   if (R==NULL || z==NULL) {
01385     (void)fprintf(stderr, "rng_normal_vector2: failed calloc\n");
01386     return 0;     
01387   }  
01388   /* copy sigma into R */
01389   memcpy((void *) R, (void *) sigma, sizeof(double)*d*d);
01390   /* find cholesky factor */
01391   if (!rng_cholesky(R, d)) {
01392     free(R); free(z); return(0);
01393   }
01394   /* now, R[i*d + j], j<=i, contains the Cholesky factor */
01395   /* sample each x(i), 1 <= i <= n */
01396   for (i = 0; i < n; i++) {
01397     /* make a fresh z for each x(i) */
01398     for (j = 0; j < d; j++) 
01399       z[j] = rng_normal(); /* z is d, iid, N(0,1) RV's */
01400     /* compute x(i) = mu + R * z  */
01401     /* note that R is lower triangular */
01402     for (j = 0; j < d; j++) {
01403       /* x(i,j) = mu(j) */
01404       x[i*stridex1+j*stridex2] = mu ? mu[j*stridemu] : 0.0;
01405       /* perform x(i,j) += R(j,1:j) <dot> z(1:j) */
01406       for (k = 0; k <= j; k++) 
01407     x[i*stridex1+j*stridex2] += R[j*d+k] * z[k];
01408     }
01409   }
01410   free((void *) R); free((void *) z); 
01411   return(1); /* OK */
01412 }
01413 
01414 
01415 /*************************************************************
01416  *
01417  * Wishart matrices
01418  *
01419  *************************************************************/
01420 
01421 /*
01422   Wishart Random Matrices: Background
01423   
01424   References:
01425     James E. Gentle, "Random Number Generation and Monte Carlo Methods,"
01426     Springer, 1998.
01427     See the description surrounding Algorithm 5.8.  We use this idea.
01428     But, the notation in step 3 of that algorithm does not make clear
01429     that simple matrix multiplication is what's happening.
01430 
01431     Mark E. Johnson, "Multivariate statistical simulation," 
01432     Wiley, 1987.
01433     The algorithm on p. 204 has clearer notation which we follow.
01434 
01435     T. W. Anderson, "An introduction to multivariate statistical
01436     analysis," Wiley, 1984.
01437     The ideas are on pages 249-251.
01438 
01439     W. B. Smith and R. R. Hocking, R. R. "Algorithm AS 53: 
01440     Wishart Variate Generator". JRSS C, 1972, 21 (3): 341-345. 
01441     This paper apparently is the first modern reference for the 
01442     simplified algorithm.
01443 
01444  
01445   turmon 2/2010: checked against Matlab wishrnd for agreement
01446   with wishrnd from the statistics toolbox.  Found good agreement
01447   for randomized sigma, both in expected value of ensembles of
01448   ~10^4 draws of W, and in standard deviation of the ensembles of W's.
01449 */
01450 
01451 /* 
01452  * Generate NxN Wishart matrix with dof degrees of freedom.
01453  * This is equivalent to an outer product of "dof" random vectors
01454  * from N(0,sigma), but done using only NxN matrix operations using
01455  * the "Bartlett decomposition" referred to above.
01456  * The idea is as follows.  Suppose:
01457  *    T = [upper triangular, as described below]
01458  * A Wishart matrix for sigma = I is, according to the Bartlett
01459  * decomposition,
01460  *    W = T T'
01461  * Suppose we have the upper-triangular cholesky factorization:
01462  *    sigma = G' G
01463  * A Wishart matrix for general sigma is
01464  *    W = G' T' T G = (T G)' * (T G) = Z' Z
01465  * where 
01466  *    Z = T G.
01467  * All operations are at at most O(N^3), which is better
01468  * than the outer product method of O(N^2 * dof) if dof >> N.
01469  */
01470 
01471 static
01472 int
01473 rng_wishart_chol(double *sigma, double dof, int N, double *w)
01474 {
01475   int i,j;
01476   double *G = NULL;  // cholesky factor
01477   double *T = NULL;  // triangular random matrix for Bartlett decomp.
01478   double *Z = NULL;  // workspace
01479   int ok = 0;        // suppose not OK
01480 
01481   /* create temporary storage */
01482   G = (double *) calloc(N*N, sizeof(double));
01483   T = (double *) calloc(N*N, sizeof(double));
01484   Z = (double *) calloc(N*N, sizeof(double));
01485   if (!G || !T || !Z)
01486     goto done; // ok=0, so will free and return failure
01487 
01488   /* get upper triangle and diagonal values:
01489    *   upper triangle: iid ~N(0,1)
01490    *   diagonal elements ~ Chi2(dof) ... Chi2(dof-N+1)
01491    */
01492   for (j = 0; j < N; j++)
01493     for (i = 0; i < N; i++)
01494       if (i == j)
01495     T[i+N*i] = sqrt(rng_chi2(dof - ((double) i)));
01496       else if (i < j)
01497     T[i+N*j] = rng_normal();
01498       else
01499     T[i+N*j] = 0.0;
01500 
01501   // G = chol(sigma)
01502   bcopy(sigma, G, N*N*sizeof(double)); // in-place cholesky
01503   if (!rng_cholesky(G, N))
01504     goto done;
01505 
01506   // Z = T * G.  Everything is upper-triangular.
01507   rng_mult_upper_triangles(Z, T, G, N);
01508 
01509   // w = Z * Z'
01510   rng_square_upper_triangle(w, Z, N);
01511 
01512   // mxt_put_matrix("S", -1, sigma, N, N);
01513   // mxt_put_matrix("G", -1, G, N, N);
01514   // mxt_put_matrix("T", -1, T, N, N);
01515   // mxt_put_matrix("Z", -1, Z, N, N);
01516 
01517   ok = 1; // successful exit
01518  done:
01519   free(G);
01520   free(T);
01521   free(Z);
01522   return ok;
01523 }
01524 
01525 /*
01526  * rng_wishart_outer: Wishart sample "w" using outer products
01527  * Generates the Wishart matrix directly, using outer products
01528  * of normal random vectors drawn from N(0,sigma).
01529  * Returns 0 on failure, 1 on success.
01530  * 
01531  * Best for dof not much more than N, but any dof is legal.
01532  */
01533 static
01534 int
01535 rng_wishart_outer(double *sigma, double dof, int N, double *w)
01536 {
01537   double *x, *x1;    /* indexes into random draw */
01538   int i, j, s;
01539   int ok;
01540 
01541   /* create temporary storage */
01542   if ((x = (double *) calloc(N*dof, sizeof(double))) == NULL)
01543     return 0;
01544   /* make normal RV's */
01545   ok = rng_normal_vector(x,     /* samples */
01546              NULL,    /* mean, taken to be 0.0 if NULL */
01547              sigma,   /* covariance */
01548              dof,     /* number of samples */
01549              N,       /* dims of mean, samples */
01550              N,       /* stride for each sample */
01551              1,       /* stride within each sample */
01552              1);      /* stride for elements in mu */
01553   if (!ok) {
01554     free(x);
01555     return 0; /* sigma not > 0, or alloc fail in rng */
01556   }
01557   /* find unnormalized outer product of dof vectors */
01558   bzero(w, N*N*sizeof(*w));
01559   for (s = 0; s < dof; s++) {
01560     x1 = x + s*N; /* the s'th random vector */
01561     for (i = 0; i < N; i++)
01562       for (j = 0; j < N; j++)
01563     w[i*N+j] += x1[i] * x1[j];
01564   }
01565   /* done */
01566   free(x);
01567   return 1;
01568 }
01569 
01570 /*
01571  * Draw one n-by-n matrix "w" from Wishart distribution with 
01572  * "dof" degrees of freedom, parameterized by covariance matrix 
01573  * sigma (symmetric positive definite, n-by-n).
01574  * Returns 1 for OK, 0 for not.
01575  *   (not OK: bad sigma, or allocation failure)
01576  * Any integer dof > 0 is OK.
01577  * For convenience, multiplies the sample by a constant "scale"
01578  *   (Use scale = 1 for plain outer product matrix, which is the
01579  *    true Wishart r.v.  Use scale = 1/dof for a sample covariance 
01580  *    matrix that approaches sigma as dof grows.)
01581  */
01582 
01583 int
01584 rng_wishart(double* sigma, double dof, double scale, int N, double* w)
01585 {
01586   int ok;
01587   int i;
01588 
01589   if (dof > N) {
01590     /* Wishart sample using Bartlett decomposition, good for large dof */
01591     ok = rng_wishart_chol(sigma, dof, N, w);
01592   } else {
01593     /* Wishart sample using naive outer products, good for small dof */
01594     ok = rng_wishart_outer(sigma, dof, N, w);
01595   }
01596   // rescale if needed
01597   if (ok && scale != 1.0)
01598     for (i = 0; i < N*N; i++)
01599       w[i] *= scale;
01600 
01601   return ok; /* 1 for OK, 0 for not */
01602 }
01603 