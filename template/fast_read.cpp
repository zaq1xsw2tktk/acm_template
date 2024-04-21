void read(int &x)
{
    char s = getchar();
    x = 0;
    ll f = 1;
    while(s < '0' || s > '9'){if(s == '-')f = -1; s = getchar();}
    while(s >= '0' && s <= '9'){x = (x << 3) + (x << 1) + (s ^ 48); s = getchar();}
    x *= f;
}

void read(double &x)
{
    char s = getchar();
    x = 0;
    ll f = 1;
    while(s < '0' || s > '9'){if(s == '-')f = -1; s = getchar();}
    while(s >= '0' && s <= '9'){x = x * 10 + (s ^ 48); s = getchar();}
    int p = 1;
    if (s=='.') {
        s = getchar();
        while(s >= '0' && s <= '9'){x = x * 10 + (s ^ 48); s = getchar(); p *= 10;}
    }
    x *= f;
    x /= p;
}
