template<class T>
class SegmentTree {
public:
    int N;
    vector<T> t;
    function<T(T, T)> fun;
    T base_val;
    SegmentTree(int n, function<T(T, T)> fun=[](T a, T b){return a + b;}, T base_val=T()) {
        this->fun = fun;
        this->base_val = base_val;
        N = n + 4;
        t.resize(N * 4);
    }
    void upd(int pos, T val) {
        upd(pos, val, 1, 1, N - 1);
    }
    void upd(int pos, T val, int index, int l, int r) {
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
    T get_val(int L, int R) {
        if (L > R) return base_val;
        return get_val(L, R, 1, 1, N - 1);
    }
    T get_val(int L, int R, int index, int l, int r) {
        if (r < L || l > R) return base_val;
        if (L <= l && r <= R) return t[index];
        long long mid = (l + r) / 2;
        return fun(get_val(L, R, index * 2, l, mid), get_val(L, R, index * 2 + 1, mid + 1, r));
    }
};

// seg update

// template<class T>
// class SegmentTree {
// public:
//     int N;
//     vector<T> t;
//     T base_val;
//     SegmentTree(int n, int base_val = -1) {
//         this->base_val = base_val;
//         N = n + 4;
//         t.resize(N * 4);
//     }
//     void upd(int l, int r, T val) {
//         upd(l, r, val, 1, 1, N - 1);
//     }
//     void upd(int L, int R, T val, int index, int l, int r) {
//         if (R < l || L > r) return;
//         if (L <= l && r <= R ) {
//             t[index] = val;
//             return;
//         }
//         long long mid = (l + r) / 2;
//         upd(L, R, val, index * 2, l, mid);
//         upd(L, R, val, index * 2 + 1, mid + 1, r);
//     }
//     T get_val(int pos) {
//         return get_val(pos, 1, 1, N - 1);
//     }
//     T get_val(int pos, int index, int l, int r) {
//         if (pos < l || pos > r) return base_val;
//         if (l == r) return t[index];
//         long long mid = (l + r) / 2;
//         return max(t[index], max(get_val(pos, index * 2, l, mid), get_val(pos, index * 2 + 1, mid + 1, r)));
//     }
// };

// 区间加，区间和
// SegmentTree<double> d1(n + 5), d2(n + 5);
// auto f = [&](int l, int r) {
//     double sr = r * d1.get_val(1, r) - d2.get_val(1, r);
//     double sl = (l - 1) * d1.get_val(1, l - 1) - d2.get_val(1, l - 1);
//     return sr - sl;
// };
// auto query_add = [&] (int l, int r, double val) {
//         d1.upd(l, d1.get_val(l, l) + val);
//         d1.upd(r + 1, d1.get_val(r + 1, r + 1) - val);
//         d2.upd(l, d2.get_val(l, l) + val * (l - 1));
//         d2.upd(r + 1, d2.get_val(r + 1, r + 1) - val * r);
// };

template<class T>
class SegmentTree {
public:
    int N;
    vector<T> d, p;
    function<T(T, T)> fun;
    T base_val;
    SegmentTree(int n, function<T(T, T)> fun=[](T a, T b){return max(a, b);}, T base_val=T()) {
        this->fun = fun;
        this->base_val = base_val;
        N = n + 4;
        d.resize(N * 4);
        p.resize(N * 4);
    }
    void segment_inc(int l, int r, int c, int s=1, int t=N-1, int p=1) {
        // [l, r] 为修改区间, c 为被修改的元素的变化量, [s, t] 为当前节点包含的区间, p
        // 为当前节点的编号
        if (l <= s && t <= r) {
            d[p] += (t - s + 1) * c, b[p] += c;
            return;
        }  // 当前区间为修改区间的子集时直接修改当前节点的值,然后打标记,结束修改
        int m = s + ((t - s) >> 1);
        if (b[p] && s != t) {
            // 如果当前节点的懒标记非空,则更新当前节点两个子节点的值和懒标记值
            d[p * 2] += b[p] * (m - s + 1), d[p * 2 + 1] += b[p] * (t - m);
            b[p * 2] += b[p], b[p * 2 + 1] += b[p];  // 将标记下传给子节点
            b[p] = 0;                                // 清空当前节点的标记
        }
        if (l <= m) segment_inc(l, r, c, s, m, p * 2);
        if (r > m) segment_inc(l, r, c, m + 1, t, p * 2 + 1);
        d[p] = d[p * 2] + d[p * 2 + 1];
    }
    int getsum(int l, int r, int s, int s=1, int t=N-1, int p=1) {
        // [l, r] 为查询区间, [s, t] 为当前节点包含的区间, p 为当前节点的编号
        if (l <= s && t <= r) return d[p];
        // 当前区间为询问区间的子集时直接返回当前区间的和
        int m = s + ((t - s) >> 1);
        if (b[p]) {
            // 如果当前节点的懒标记非空,则更新当前节点两个子节点的值和懒标记值
            d[p * 2] += b[p] * (m - s + 1), d[p * 2 + 1] += b[p] * (t - m);
            b[p * 2] += b[p], b[p * 2 + 1] += b[p];  // 将标记下传给子节点
            b[p] = 0;                                // 清空当前节点的标记
        }
        int sum = 0;
        if (l <= m) sum = getsum(l, r, s, m, p * 2);
        if (r > m) sum += getsum(l, r, m + 1, t, p * 2 + 1);
        return sum;
    }
};

// 区间add 区间max
template <typename T>
class SegTreeLazyRangeAdd {
  vector<T> tree, lazy;
  vector<T> *arr;
  int n, root, n4, end;

  void maintain(int cl, int cr, int p) {
    int cm = cl + (cr - cl) / 2;
    if (cl != cr && lazy[p]) {
      lazy[p * 2] += lazy[p];
      lazy[p * 2 + 1] += lazy[p];
      tree[p * 2] += lazy[p];// * (cm - cl + 1);
      tree[p * 2 + 1] += lazy[p];// * (cr - cm);
      lazy[p] = 0;
    }
  }

  T range_min(int l, int r, int cl, int cr, int p) {
    if (l <= cl && cr <= r) return tree[p];
    int m = cl + (cr - cl) / 2;
    T min_val = 2e18;
    maintain(cl, cr, p);
    if (l <= m) min_val = std::min(min_val, range_min(l, r, cl, m, p * 2));
    if (r > m) min_val = std::min(min_val, range_min(l, r, m + 1, cr, p * 2 + 1));
    return min_val;
  }

  void range_add(int l, int r, T val, int cl, int cr, int p) {
    if (l <= cl && cr <= r) {
      lazy[p] += val;
      tree[p] += val;
      return;
    }
    int m = cl + (cr - cl) / 2;
    maintain(cl, cr, p);
    if (l <= m) range_add(l, r, val, cl, m, p * 2);
    if (r > m) range_add(l, r, val, m + 1, cr, p * 2 + 1);
    tree[p] = std::min(tree[p * 2], tree[p * 2 + 1]);
  }

  void build(int s, int t, int p) {
    if (s == t) {
      tree[p] = (*arr)[s];
      return;
    }
    int m = s + (t - s) / 2;
    build(s, m, p * 2);
    build(m + 1, t, p * 2 + 1);
    tree[p] = min(tree[p * 2], tree[p * 2 + 1]);
  }

 public:
  explicit SegTreeLazyRangeAdd<T>(vector<T> v) {
    n = v.size();
    n4 = n * 4;
    tree = vector<T>(n4, 0);
    lazy = vector<T>(n4, 0);
    arr = &v;
    end = n - 1;
    root = 1;
    build(0, end, 1);
    arr = nullptr;
  }

  void upd(int l, T val) {
    T need_add_val = val - range_min(l, l);
    range_add(l, l, need_add_val);
  }

  void show(int p, int depth = 0) {
    if (p > n4 || tree[p] == 0) return;
    show(p * 2, depth + 1);
    for (int i = 0; i < depth; ++i) putchar('\t');
    printf("%d:%d\n", tree[p], lazy[p]);
    show(p * 2 + 1, depth + 1);
  }

  T range_min(int l, int r) { return range_min(l, r, 0, end, root); }

  void range_add(int l, int r, T val) { range_add(l, r, val, 0, end, root); }
};
