template<class T>
class SegmentTree {
public:
    int N;
    vector<T> t;
    function<T(T, T)> fun;
    T base_val;
    SegmentTree<T>(int n, function<T(T, T)> fun=[](T a, T b){return a + b;}, T base_val=T()) {
        this->fun = fun;
        this->base_val = base_val;
        N = n + 4;
        t.resize(N * 4);
    }
    void upd(int pos, T val, int index = 1, int l = 1, int r = -1) {
        if (r == -1) r = N - 1;
        if (pos < l || pos > r) return;
        if (l == r) {
            t[index] = val;
            return;
        }
        long long mid = (l + r) / 2;
        upd(pos, val, index * 2, l, mid);
        upd(pos, val, index * 2 + 1, mid + 1, r);
        t[index] = fun(t[index * 2], t[index * 2 + 1]);
    }
    T get_val(int L, int R, int index = 1, int l = 1, int r = -1) {
        if (r == -1) r = N - 1;
        if (L > R) return base_val;
        if (r < L || l > R) return base_val;
        if (L <= l && r <= R) return t[index];
        long long mid = (l + r) / 2;
        return fun(get_val(L, R, index * 2, l, mid), get_val(L, R, index * 2 + 1, mid + 1, r));
    }
};
