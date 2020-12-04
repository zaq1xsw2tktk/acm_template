const long long N = 1e5 + 4;
int t[4*N];

void upd(int pos, int val, int index = 1, int l = 1, int r = N - 1) {
    if (pos < l || pos > r) return;
    if (l == r) {
        t[index] = val;
        return;
    }
    int mid = (l + r) / 2;
    upd(pos, val, index * 2, l, mid);
    upd(pos, val, index * 2 + 1, mid + 1, r);
    t[index] = min(t[index * 2], t[index * 2 + 1]);
}

int get_min(int L, int R, int index = 1, int l = 1, int r = N - 1) {
    if (L > R) return 1e9;
    if (r < L || l > R) return 1e9;
    if (L <= l && r <= R) return t[index];
    int mid = (l + r) / 2;
    return std::min(get_min(L, R, index * 2, l, mid), get_min(L, R, index * 2 + 1, mid + 1, r));
}