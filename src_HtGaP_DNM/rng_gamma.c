double rng_gamma(double alpha, double theta)
{
// follow DeVroye's notation
 double X, U, V, W, Y, Z;  // random variates
 double b, c, t;           // parameters

if (alpha <= 0 || theta < 0)
 return rng_getnand();
  if (alpha == 1.0) {
    // Use the exponential, e.g. Devroye p. 405
     X = -log(rng_uniform());
   } else if (alpha == 2.0) {
     // Use the sum of 2 exponentials (Devroye p. 405)
     X = -log(rng_uniform() * rng_uniform());
   } else if (alpha < 1.0) {
     // Use RGS algorithm (Best, 1983; Devroye, p. 426)
     c = 1.0 / alpha; 
     t = 0.07 + 0.75 * sqrt(1.0 - alpha); // breakpoint: min reject prob
     b = 1.0 + exp(-t) * alpha / t;
     // generate a variate x
    while (1) {
       U = rng_uniform(); 
      W = rng_uniform(); 
       V = b * U;
       if (V <= 1.0) {
    X = t * pow(V, c);
     /* the first step is an early bailout because
      * exp(-x) >= (2-x)/(2+x) for x>=0
     */
     // test is equivalent to w <= (2-x)/(2+x) since x>=0 
     if (W * (2.0 + X) <= (2.0 - X)) break;
    if (W <= exp(-X)) break;
       } else {
     X = -log(c * t * (b - V));
     Y = X / t;
    if ((W * (alpha + Y - alpha*Y)) <= 1.0) break;
     if (W <= pow(Y, alpha - 1.0)) break;
       }
    }
   } else {
     // Use rejection algorithm XG (Best, 1978; Devroye, p. 410)
     //   (This algorithm is OK for alpha > 1.  If alpha = 1,
    //   b = 0, and the last rejection step with log(X/b) is
     //   illegal.  This is noted in Devroye's errata on the web.)
    b = alpha - 1.0;
     c = 3.0 * alpha - 0.75;
     // generate a variate x
     while (1) {
       U = rng_uniform();  
       V = rng_uniform();
      W = U * (1.0 - U);  
       Y = sqrt(c / W) * (U - 0.5);
       X = b + Y;
       if (X >= 0.0) {
    Z = 64.0 * W*W*W * V*V;
     /* the first step is an early bailout which
     * Devroye remarks can be omitted at little cost.
      */
     // test is same as Z <= 1 - 2*Y*Y/X since X >= 0 here
    if ((Z * X) <= (X - 2.0 * Y * Y)) break;
     if (log(Z) <= (2.0 * (b * log(X/b) - Y))) break;
       }
    }
   }
   return X * theta;
}