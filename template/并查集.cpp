const int N = 125 + 3;
int p[N];
int fnode(int x) {return p[x] == x ? x : (p[x] = fnode(p[x]));}
void connect(int x, int y) {p[fnode(x)] = fnode(y);}
void init() {for(int i = 0; i < N; i++)p[i] = i;}
