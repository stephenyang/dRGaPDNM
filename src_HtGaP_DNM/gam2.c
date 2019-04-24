float rgama(float a) 
{
    float d,c,x,v,u;
    d = a-1.0/3.0;
    c = 1.0/sqrt(9.0*d);
        for(;;)
        {
            do {x5RNOR;
            v=1.+c*x;}
        while(v < 0.);
        v=v*v*v;
        u = UNI;
        if( u<1.-.0331*(x*x)*(x*x) )
            return (d*v);
        if( log(u) < 0.5*x*x + d*(1.0-v+log(v)))
            break;
        }
    return (d*v);
}