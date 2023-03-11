const int N = 200005;
int p[N];
int fnode(int x) {return p[x] == x ? x : (p[x] = fnode(p[x]));}
void connect(int x, int y) {p[fnode(x)] = fnode(y);}
void init(int n = N) {for(int i = 0; i < n; i++)p[i] = i;}