long long gcd_ext(long long a,long long b,long long &x,long long&y)
{
    if(b==0)
    {
        x=1, y=0;
        return a ;
    }
    else 
    {
        long long r = gcd_ext(b, a % b, x, y);
        long long t = x;
        x = y;
        y = t - a / b * y;
        return r ;
    }
}