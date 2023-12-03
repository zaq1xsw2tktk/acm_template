#include <bits/stdc++.h>
using namespace std;
#define FOR(i, n) for (int (i) = 0; (i) < (n); (i)++)
#define FORI(i, a, b) for (int (i) = (a); (i) < (b); (i)++)
 
mt19937 rnd(time(0));
 
#define ll long long
#define vi vector<int>
#define vvi vector<vector<int> >
#define endl '\n'
 
#define mp(m, n) make_pair((m), (n))
 
template<typename T>
void read(vector<T> &t) {FOR(i, t.size()) {cin >> t[i];}}
template<typename T> string tostring(T a) { istringstream sin; sin >> a; return sin.str(); }
 
// #define DEBUG
 
#ifdef DEBUG
template<typename T>
void _debug(string s, T x) {
    cerr << s << ":";
    for (auto it = x.begin(); it != x.end(); ++it) {
        cerr << " " << *it;
    }
    cerr << endl;
}
 
template<typename T, typename K>
void _debug(string s, unordered_map<T, K> x) {
    cerr << s << ":";
    for (auto it = x.begin(); it != x.end(); ++it) {
        cerr << " " << it->first << ":" << it->second;
    }
    cerr << endl;
}
 
template<typename T, typename K>
void _debug(string s, map<T, K> x) {
    cerr << s << ":";
    for (auto it = x.begin(); it != x.end(); ++it) {
        cerr << " " << it->first << ":" << it->second;
    }
    cerr << endl;
}
 
template<typename T, typename K>
void _debug(string s, set<T, K> x) {
    cerr << s << ":";
    for (auto it = x.begin(); it != x.end(); ++it) {
        cerr << " " << *it;
    }
    cerr << endl;
}
 
template<typename T, typename K>
void _debug(string s, vector<pair<T, K> > x) {
    cerr << s << ":";
    for (auto it = x.begin(); it != x.end(); ++it) {
        cerr << " " << it->first << "," << it->second;
    }
    cerr << endl;
}
 
template<typename T, typename K>
void _debug(string s, set<pair<T, K> > x) {
    cerr << s << ":";
    for (auto it = x.begin(); it != x.end(); ++it) {
        cerr << " " << it->first << "," << it->second;
    }
    cerr << endl;
}
 
template<typename T, typename K>
void _debug(string s, pair<T, K> x) {
    cerr << s << ": " << x.first << "," << x.second << endl;
}
 
void _debug(string s, int x) {
    cerr << s << ": " << x << endl;
}
void _debug(string s, long long x) {
    cerr << s << ": " << x << endl;
}
void _debug(string s, double x) {
    cerr << s << ": " << x << endl;
}
void _debug(string s, string x) {
    cerr << s << ": " << x << endl;
}
void _debug(string s, char x) {
    cerr << s << ": " << x << endl;
}
void _debug(string s, size_t x) {
    cerr << s << ": " << x << endl;
}
void _debug(string s, bool x) {
    cerr << s << ": " << x << endl;
}
 
#define debug(x) _debug(#x, (x))
#else
#define debug(x)
#endif
 
#define db debug
 
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
 
void read(ll &x)
{
    char s = getchar();
    x = 0;
    ll f = 1;
    while(s < '0' || s > '9'){if(s == '-')f = -1; s = getchar();}
    while(s >= '0' && s <= '9'){x = (x << 3) + (x << 1) + (s ^ 48); s = getchar();}
    x *= f;
}
 
template <class T>
void print(vector<T> &v) {
    for (int i = 0; i < v.size(); i++) {
        cout << v[i];
        if (i == v.size() - 1) {
            cout << endl;
        } else {
            cout << " ";
        }
    }
}
 
void solve();
 
void pre_init();
 
void pre_init() {}
 
int N; // User number
int K; // Cell number
int T; // TTI number
int R; // RBG number
vector<vector<vector<vector<double> > > > s0_tkrn, s_tkrn;
vector<vector<vector<vector<double> > > > d_krmn;
vector<vector<vector<vector<int> > > > d_krmn_order;
int J;
vector<int> frameID, uID, TTI0, nTTI;
vector<double> TBS;
vector<vector<vector<vector<double> > > > p_tkrn;
vector<vector<vector<vector<int> > > > b_tkrn;
vector<vector<vector<double> > > s_tkn;
vector<vector<vector<int> > > s0_tkn_max_idxk;
vector<vector<int> > index_tn;
 
 
const double W = 192;
 
void solve();
 
void calc_b_tkrn() {
    FOR(t, T) FOR(k, K) FOR(r, R) FOR (n, N) {
        b_tkrn[t][k][r][n] = (p_tkrn[t][k][r][n] > 0);
    }
}
 
void calc_b_tkrn(int t) {
    FOR(k, K) FOR(r, R) FOR (n, N) {
        b_tkrn[t][k][r][n] = (p_tkrn[t][k][r][n] > 0);
    }
}
 
void calc_s_tkrn() {
    calc_b_tkrn();
    FOR(t, T) FOR(k, K) FOR(n, N) FOR(r, R) {
        // debug(t);
        if (p_tkrn[t][k][r][n] == 0) continue;
        double fz = s0_tkrn[t][k][r][n] * p_tkrn[t][k][r][n];
        // debug(t);
        FOR (m, N) {
            if (m == n) continue;
            fz *= exp(d_krmn[k][r][m][n] * b_tkrn[t][k][r][m]);
        }
        double fm = 1;
        FOR (k_, K) FOR (n_, N) {
            if (k_ == k || n_ == n) continue;
            fm += s0_tkrn[t][k_][r][n_] * p_tkrn[t][k_][r][n_] * exp(-d_krmn[k_][r][n_][n]);
        }
        // debug(s0_tkrn[t][k][r][n]);
        // debug(p_tkrn[t][k][r][n]);
        // debug(fz);
        // debug(fm);
        s_tkrn[t][k][r][n] = fz / fm;
        // debug(s_tkrn[t][k][r][n]);
    }
}
 
void calc_s_tkrn(int t) {
    calc_b_tkrn(t);
    FOR(k, K) FOR(n, N) FOR(r, R) {
        if (p_tkrn[t][k][r][n] == 0) continue;
        double fz = s0_tkrn[t][k][r][n] * p_tkrn[t][k][r][n];
        FOR (m, N) {
            if (m == n) continue;
            fz *= exp(d_krmn[k][r][m][n] * b_tkrn[t][k][r][m]);
        }
        double fm = 1;
        FOR (k_, K) FOR (n_, N) {
            if (k_ == k || n_ == n) continue;
            fm += s0_tkrn[t][k_][r][n_] * p_tkrn[t][k_][r][n_] * exp(-d_krmn[k_][r][n_][n]);
        }
        s_tkrn[t][k][r][n] = fz / fm;
        // debug(fm);
    }
}
 
void calc(int t) {
    calc_s_tkrn(t);
    FOR(k, K) FOR(n, N) {
        double val = 1;
        int cnt = 0;
        FOR(r, R) {
            if (b_tkrn[t][k][r][n] == 0) continue;
            val *= s_tkrn[t][k][r][n];
            cnt += 1;
        }
        s_tkn[t][k][n] = 0;
        if (cnt > 0) {
            s_tkn[t][k][n] = pow(val, 1.0 / cnt);
        }
    }
}
 
vector<double> g;
 
void calc(bool print_debug=true) {
    calc_s_tkrn();
    FOR(t, T) FOR(k, K) FOR(n, N) {
        double val = 1;
        int cnt = 0;
        FOR(r, R) {
            if (b_tkrn[t][k][r][n] == 0) continue;
            val *= s_tkrn[t][k][r][n];
            cnt += 1;
        }
        s_tkn[t][k][n] = 0;
        if (cnt > 0) {
            s_tkn[t][k][n] = pow(val, 1.0 / cnt);
        }
        // debug(s_tkn[t][k][n]);
    }
    double score = 0;
    for (int j = 0; j < J; j++) {
        g[j] = 0;
        int n = uID[j];
        for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
            FOR (k, K) FOR (r, R) {
                g[j] += W * b_tkrn[t][k][r][n] * log2(1 + s_tkn[t][k][n]);
            }
        }
        if (g[j] >= TBS[j]) score ++;
        if (g[j] > 0 && g[j] < TBS[j]) cerr << "bad ";
        cerr << j << " " << g[j] << " " << TBS[j] << endl;
        if (g[j] > 0 && g[j] < TBS[j]) {
            TBS[j] = 1e9;
        }
    }
    FOR(t, T) FOR(k, K) FOR(n, N) FOR(r, R) {
        score -= p_tkrn[t][k][r][n] * 1e-6;
    }
    if(print_debug) debug(g);
    if(print_debug) debug(score);
}
 
void clear(int t) {
    FOR(k, K) FOR (r, R) FOR (n, N) p_tkrn[t][k][r][n] = 0;
}
void clear(int t, int r) {
    FOR(k, K) FOR (n, N) p_tkrn[t][k][r][n] = 0;
}

pair<vector<double>, vector<vector<int> > > solve_one_t_list_v2_use_mask(int t, int r, const vector<pair<double, int> > &need_check, int N_, int k, double &save, double p_remain_k=4) {
    int N = need_check.size();
    vector<double> use_p(N + 1, 5.0);
    vector<vector<int> > filled_list(N + 1);
    if (N == 0) return {use_p, filled_list};
 
    vector<ll> filled_list_mask(N + 1);

    vector<int> idx_revers(N_, -1); // uID 映射到 need_check
    FOR (n, N) idx_revers[need_check[n].second] = n;
    vector<double> s0_n(N);
    vector<vector<double> > d_mn;
    vector<double> remain_TBS(N); FOR (n, N) remain_TBS[n] = need_check[n].first;
    FOR(n, N) s0_n[n] = s0_tkrn[t][k][r][need_check[n].second];
    d_mn.resize(N, vector<double>(N));
    FOR(m, N) FOR(n, N) d_mn[m][n] = d_krmn[k][r][need_check[m].second][need_check[n].second];
    int max_score = 0;
    set<ll> visited;

    vector<ll> masks;
    FOR (n, N) {
        masks.push_back(1LL << n);
    }
    while (true) {
        vector<pair<double, ll>> score_mask;
        for (auto mask : masks) {
            FOR (n, N) {
                if (mask & (1LL << n)) continue;
                ll new_u = mask | (1LL << n);
                if (visited.count(new_u)) continue;
                visited.insert(new_u);
                vector<double> need;
                for (int n = 0; n < N; n++) {
                    if (!(new_u & (1LL << n))) continue;
                    int now = n;
                    double p = 0;
                    FOR (m, N) if (m != n &&((1LL << m) & new_u)) p += d_mn[n][m];
                    p = exp(p);
                    double ne = pow(2.0, remain_TBS[now] / W) - 1;
                    ne /= p;
                    ne /= s0_n[now];
                    need.push_back(ne + 1e-6);
                    // need[n] = ne + 1e-6;
                }
                double up = std::accumulate(need.begin(), need.end(), 0.0);
                if (up < min(p_remain_k, 4.0) && up < use_p[__builtin_popcountll(new_u)]) {
                    use_p[__builtin_popcountll(new_u)] = up;
                    filled_list_mask[__builtin_popcountll(new_u)] = new_u;
                    score_mask.push_back({up, new_u});
                }
                if (up < min(p_remain_k, 4.0)) {
                    score_mask.push_back({up, new_u});
                }
            }
        }
        if (score_mask.size() == 0) break;
        sort(score_mask.begin(), score_mask.end());
        // debug(score_mask);
        vector<ll> new_mask;
        for (int i = 0; i < score_mask.size() && i < 100; i++) {
            new_mask.push_back(score_mask[i].second);
        }
        masks.swap(new_mask);
    }
    FOR (n, N) {
        FOR (n1, N) {
            if ((1LL << n1) & filled_list_mask[n]) {
                filled_list[n].push_back(n1);
            }
        }
        if (filled_list[n].size() == 0) continue;
        debug(use_p[n]);
        debug(filled_list[n]);
    }
    // if (t == 5) exit(0);
    return {use_p, filled_list};


    // // 贪心
    // FOR (n, N) {
    //     ll used = 0;
    //     used |= (1 << n);
    //     // uu.push_back(n); used[n] = true;
    //     while (true) {
    //         ll best_new_u = 0;
    //         double best_new_u_score = 1e5;
    //         FOR (n1, N) {
    //             if (used & (1LL << n1)) continue;
    //             ll new_u = used;
    //             new_u |= 1LL << n1;
    //             if (visited.count(new_u)) {
    //                 continue;
    //             }
    //             visited.insert(new_u);

    //             vector<double> need;
    //             for (int n = 0; n < N; n++) {
    //                 if (!(new_u & (1LL << n))) continue;
    //                 int now = n;
    //                 double p = 0;
    //                 FOR (m, N) if (m != n &&((1LL << m) & new_u)) p += d_mn[n][m];
    //                 p = exp(p);
    //                 double ne = pow(2.0, remain_TBS[now] / W) - 1;
    //                 ne /= p;
    //                 ne /= s0_n[now];
    //                 need.push_back(ne + 1e-6);
    //                 // need[n] = ne + 1e-6;
    //             }
    //             double up = std::accumulate(need.begin(), need.end(), 0.0);
    //             if (up < min(p_remain_k, 4.0) && up < use_p[__builtin_popcountll(new_u)]) {
    //                 use_p[__builtin_popcountll(new_u)] = up;
    //                 filled_list_mask[__builtin_popcountll(new_u)] = new_u;
    //             }
    //             if (up < min(p_remain_k, 4.0) && up < best_new_u_score) {
    //                 best_new_u_score = up;
    //                 best_new_u = new_u;
    //             }
    //         }
    //         if (__builtin_popcountll(best_new_u) == 0) break;
    //         used = best_new_u;
    //     }
    // }
    // debug(filled_list_mask);
    // FOR (n, N) {
    //     FOR (n1, N) {
    //         if ((1LL << n1) & filled_list_mask[n]) {
    //             filled_list[n].push_back(n1);
    //         }
    //     }
    //     if (filled_list[n].size() == 0) continue;
    //     debug(use_p[n]);
    //     debug(filled_list[n]);
    // }
    // // if (t == 5) exit(0);
    // return {use_p, filled_list};
}

pair<vector<double>, vector<vector<int> > > solve_one_t_list_v2(int t, int r, const vector<pair<double, int> > &need_check, int N_, int k, double &save, double p_remain_k=4) {
    int N = need_check.size();
    vector<double> use_p(N + 1, 5.0);
    vector<vector<int> > filled_list(N + 1);
    if (N == 0) return {use_p, filled_list};
    if (N <= 64) return solve_one_t_list_v2_use_mask(t, r, need_check, N_, k, save, p_remain_k);
 
    vector<int> idx_revers(N_, -1); // uID 映射到 need_check
    FOR (n, N) idx_revers[need_check[n].second] = n;
    vector<double> s0_n(N);
    vector<vector<double> > d_mn;
    vector<double> remain_TBS(N); FOR (n, N) remain_TBS[n] = need_check[n].first;
    FOR(n, N) s0_n[n] = s0_tkrn[t][k][r][need_check[n].second];
    d_mn.resize(N, vector<double>(N));
    FOR(m, N) FOR(n, N) d_mn[m][n] = d_krmn[k][r][need_check[m].second][need_check[n].second];
    int max_score = 0;
    // 贪心
    FOR (n, N) {
        vector<int> used(N);
        vector<int> uu;
        uu.push_back(n); used[n] = true;
        while (true) {
            vector<int> best_new_u;
            double best_new_u_score = 1e5;
            FOR (n1, N) {
                if (used[n1]) continue;
                vector<int> new_u = uu;
                new_u.push_back(n1);
                vector<double> need(new_u.size());
                for (int n = 0; n < new_u.size(); n++) {
                    int now = new_u[n];
                    double p = 0;
                    FOR (m, new_u.size()) if (m != n) p += d_mn[new_u[n]][new_u[m]];
                    p = exp(p);
                    double ne = pow(2.0, remain_TBS[now] / W) - 1;
                    ne /= p;
                    ne /= s0_n[now];
                    need[n] = ne + 1e-6;
                }
                double up = std::accumulate(need.begin(), need.end(), 0.0);
                if (up < min(p_remain_k, 4.0) && up < use_p[new_u.size()]) {
                    use_p[new_u.size()] = up;
                    filled_list[new_u.size()] = new_u;
                }
                if (up < min(p_remain_k, 4.0) && up < best_new_u_score) {
                    best_new_u_score = up;
                    best_new_u = new_u;
                }
            }
            if (best_new_u.size() == 0) break;
            uu = best_new_u;
            used[uu.back()] = 1;
        }
    }
    FOR (n, N) {
        if (filled_list[n].size() == 0) continue;
        debug(use_p[n]);
        debug(filled_list[n]);
    }
    // if (t == 5) exit(0);
    return {use_p, filled_list};
 
}
 
pair<vector<double>, vector<vector<int> > > solve_one_t_list(int t, int r, const vector<pair<double, int> > &need_check, int N_, int k, double &save, double p_remain_k=4) {
    int N = need_check.size();
    vector<double> use_p(N + 1, 5.0);
    vector<vector<int> > filled_list(N + 1);
    if (N == 0) return {use_p, filled_list};
 
    vector<int> idx_revers(N_, -1); // uID 映射到 need_check
    FOR (n, N) idx_revers[need_check[n].second] = n;
    vector<double> s0_n(N);
    vector<vector<double> > d_mn;
    vector<double> remain_TBS(N); FOR (n, N) remain_TBS[n] = need_check[n].first;
    FOR(n, N) s0_n[n] = s0_tkrn[t][k][r][need_check[n].second];
    d_mn.resize(N, vector<double>(N));
    FOR(m, N) FOR(n, N) d_mn[m][n] = d_krmn[k][r][need_check[m].second][need_check[n].second];
    int max_score = 0;
    // 贪心
    FOR (n, N) {
        vector<int> ori_order = d_krmn_order[k][r][need_check[n].second];
        vector<int> order;
        for (auto x : ori_order) {
            if (idx_revers[x] != -1) {
                order.push_back(idx_revers[x]);
            }
        }
 
        vector<int> u;
        for (int use = 0; use < order.size(); use++) {
            u.push_back(order[use]);
            if (use > 0 && u.size() == 1) {
                vector<pair<double, int> > n_idx;
                for (int kk = use + 1; kk < order.size(); kk++) {
                    int a1 = order[use];
                    int a2 = order[kk];
                    n_idx.push_back({d_mn[a1][a2], a2});
                }
                sort(n_idx.begin(), n_idx.end(), greater<pair<double, int>>());
                for (int k = use + 1; k < order.size(); k++) {
                    order[k] = n_idx[k - use - 1].second;
                }
            }
 
            vector<double> need(u.size());
            for (int n = 0; n < u.size(); n++) {
                int now = u[n];
                double p = 0;
                FOR (m, u.size()) if (m != n) p += d_mn[u[n]][u[m]];
                p = exp(p);
                double ne = pow(2.0, remain_TBS[now] / W) - 1;
                ne /= p;
                ne /= s0_n[now];
                need[n] = ne + 1e-6;
            }
            // debug(u);
            // debug(need);
            double up = std::accumulate(need.begin(), need.end(), 0.0);
            if (up < min(p_remain_k, 4.0) && up < use_p[u.size()]) {
                use_p[u.size()] = up;
                filled_list[u.size()] = u;
            }
            if (std::accumulate(need.begin(), need.end(), 0.0) >= min(2.5, p_remain_k - 0.5)) {
                u.pop_back();
            }
            if (u.size() == 0) {
                break;
            }
        }
        max_score = max<int>(max_score, u.size());
        // filled_list[n] = u;
        // debug(filled_list[n]);
    }
    FOR (n, N) {
        if (filled_list[n].size() == 0) continue;
        debug(use_p[n]);
        debug(filled_list[n]);
    }
    // if (t == 5) exit(0);
    return {use_p, filled_list};
 
}
 
 
pair<int, vector<int> > solve_one_t(int t, int r, const vector<pair<double, int> > &need_check, int N_, int k, double &save) {
    int N = need_check.size();
    vector<vector<int> > filled_list(N);
    if (N == 0) return {0, vector<int> ()};
    vector<int> idx_revers(N_, -1); // uID 映射到 need_check
    FOR (n, N) idx_revers[need_check[n].second] = n;
    vector<double> s0_n(N);
    vector<vector<double> > d_mn;
    vector<double> remain_TBS(N); FOR (n, N) remain_TBS[n] = need_check[n].first;
    FOR(n, N) s0_n[n] = s0_tkrn[t][k][r][need_check[n].second];
    d_mn.resize(N, vector<double>(N));
    FOR(m, N) FOR(n, N) d_mn[m][n] = d_krmn[k][r][need_check[m].second][need_check[n].second];
    int max_score = 0;
    // 贪心
    FOR (n, N) {
        vector<int> ori_order = d_krmn_order[k][r][need_check[n].second];
        vector<int> order;
        for (auto x : ori_order) {
            if (idx_revers[x] != -1) {
                order.push_back(idx_revers[x]);
            }
        }
 
        vector<int> u;
        for (int use = 0; use < order.size(); use++) {
            u.push_back(order[use]);
            if (use > 0 && u.size() == 1) {
                vector<pair<double, int> > n_idx;
                for (int kk = use + 1; kk < order.size(); kk++) {
                    int a1 = order[use];
                    int a2 = order[kk];
                    n_idx.push_back({d_mn[a1][a2], a2});
                }
                sort(n_idx.begin(), n_idx.end(), greater<pair<double, int>>());
                for (int k = use + 1; k < order.size(); k++) {
                    order[k] = n_idx[k - use - 1].second;
                }
            }
 
            vector<double> need(u.size());
            for (int n = 0; n < u.size(); n++) {
                int now = u[n];
                double p = 0;
                FOR (m, u.size()) if (m != n) p += d_mn[u[n]][u[m]];
                p = exp(p);
                double ne = pow(2.0, remain_TBS[now] / W) - 1;
                ne /= p;
                ne /= s0_n[now];
                need[n] = ne + 1e-6;
            }
            // debug(u);
            // debug(need);
            if (std::accumulate(need.begin(), need.end(), 0.0) > 2.0) {
                u.pop_back();
            }
            if (u.size() == 0) {
                break;
            }
        }
        max_score = max<int>(max_score, u.size());
        filled_list[n] = u;
        // debug(filled_list[n]);
    }
    if (max_score < 2) {
        return {0, vector<int> ()};
    }
    FOR (n, N) {
        vector<int> u = filled_list[n];
        if (u.size() == max_score) {
            for (int n = 0; n < u.size(); n++) {
                int now = u[n];
                double p = 0;
                FOR (m, u.size()) if (m != n) p += d_mn[u[n]][u[m]];
                p = exp(p);
                double ne = pow(2.0, remain_TBS[now] / W) - 1;
                ne /= p;
                ne /= s0_n[now];
                ne += 1e-6;
                p_tkrn[t][k][r][need_check[now].second] = ne;
            }
            return {max_score, u};
        }
    }
 
    return {0, vector<int> ()};
}
 
pair<int, vector<vector<int>> > solve_one_t(int t, const vector<pair<double, int> > &need_check, int N_, int k, double &save) {
    int N = need_check.size();
    vector<vector<int> > filled_list = vector<vector<int> >(R);
    if (N == 0) return {0, filled_list};
    // debug(need_check);
    vector<int> idx_revers(N_, -1); // uID 映射到 need_check
    FOR (n, N) idx_revers[need_check[n].second] = n;
    vector<vector<vector<double>>> s0_krn;
    s0_krn.resize(K, vector<vector<double>>(R, vector<double>(N)));
    vector<vector<vector<vector<double> > > > _d_krmn;
    _d_krmn.resize(K, vector<vector<vector<double>>>(R, vector<vector<double>>(N, vector<double>(N))));
    FOR(k, K) FOR(r, R) FOR(n, N) s0_krn[k][r][n] = s0_tkrn[t][k][r][need_check[n].second];
    FOR(k, K) FOR(r, R) FOR(m, N) FOR(n, N) _d_krmn[k][r][m][n] = d_krmn[k][r][need_check[m].second][need_check[n].second];
    vector<double> remain_TBS(N); FOR (n, N) remain_TBS[n] = need_check[n].first;
    vector<int> isFilled(N, false); int filled_cnt = 0;
    double saved_p = 0;
    FOR (r, R) {
        int begin_idx = -1;
        FOR (n, N) {
            if (!isFilled[n]) {
                begin_idx = n;
                break;
            }
        }
        if (begin_idx == -1) continue;
        vector<int> ori_order = d_krmn_order[k][r][need_check[begin_idx].second];
        vector<int> order;
        for (auto x : ori_order) {
            if (idx_revers[x] != -1 && !isFilled[idx_revers[x]]) {
                    order.push_back(idx_revers[x]);
            }
        }
        // debug(order);
        vector<int> u;
        for (int use = 0; use < order.size(); use++) {
            u.push_back(order[use]);
            if (use > 0 && u.size() == 1) {
                vector<pair<double, int> > n_idx;
                for (int kk = use + 1; kk < order.size(); kk++) {
                    int a1 = order[use];
                    int a2 = order[kk];
                    n_idx.push_back({_d_krmn[k][r][a1][a2], a2});
                }
                sort(n_idx.begin(), n_idx.end(), greater<pair<double, int>>());
                for (int k = use + 1; k < order.size(); k++) {
                    order[k] = n_idx[k - use - 1].second;
                }
            }
 
            vector<double> need(u.size());
            for (int n = 0; n < u.size(); n++) {
                int now = u[n];
                double p = 0;
                FOR (m, u.size()) if (m != n) p += _d_krmn[k][r][u[n]][u[m]];
                p = exp(p);
                double ne = pow(2.0, remain_TBS[now] / W) - 1;
                ne /= p;
                ne /= s0_krn[k][r][now];
                need[n] = ne + 1e-6;
            }
            // debug(u);
            // debug(need);
            if (std::accumulate(need.begin(), need.end(), 0.0) > 1.0) {
                u.pop_back();
            }
        }
 
        if (u.size() == 1) u.clear();
        for (auto u_: u) {
            isFilled[u_] = true;
            filled_cnt ++;
            filled_list[r].push_back(u_);
        }
        double max_all = 0;
        for (int n = 0; n < u.size(); n++) {
            int now = u[n];
            double p = 0;
            FOR (m, u.size()) if (m != n) p += _d_krmn[k][r][u[n]][u[m]];
            p = exp(p);
            double ne = pow(2.0, remain_TBS[now] / W) - 1;
            ne /= p;
            ne /= s0_krn[k][r][now];
            ne += 1e-6;
            max_all += ne;
            p_tkrn[t][k][r][need_check[now].second] = ne;
        }
        if (u.size() > 0) {
            saved_p += 1-max_all;
        }
        // debug(u);
    }
    // exit(0);
    // debug(t);
    // debug(filled_cnt);
    save = saved_p;
    return {filled_cnt, filled_list};
}
 
// int solve_one_t_old(int t, const vector<pair<double, int> > &need_check, int N_) {
//     int N = need_check.size();
//     if (N == 0) return 0;
//     debug(need_check);
//     vector<int> idx_revers(N_, -1); // uID 映射到 need_check
//     FOR (n, N) idx_revers[need_check[n].second] = n;
//     vector<vector<vector<double>>> s0_krn;
//     s0_krn.resize(K, vector<vector<double>>(R, vector<double>(N)));
//     vector<vector<vector<vector<double> > > > _d_krmn;
//     _d_krmn.resize(K, vector<vector<vector<double>>>(R, vector<vector<double>>(N, vector<double>(N))));
//     FOR(k, K) FOR(r, R) FOR(n, N) s0_krn[k][r][n] = s0_tkrn[t][k][r][need_check[n].second];
//     FOR(k, K) FOR(r, R) FOR(m, N) FOR(n, N) {
//         _d_krmn[k][r][m][n] = d_krmn[k][r][need_check[m].second][need_check[n].second];
//     }
//     vector<double> remain_TBS(N); // need_check
//     vector<int> isFilled(N, false); // need_check
//     int filled_cnt = 0;
//     FOR (n, N) remain_TBS[n] = need_check[n].first;
//     FOR (r, R) {
//         debug(r);
//         int k = 0;
//         if (filled_cnt == N) continue;
//         int begin_idx = -1;
//         FOR (n, N) {
//             if (!isFilled[n]) {
//                 begin_idx = n;
//                 break;
//             }
//         }
//         debug(begin_idx);
//         vector<int> ori_order = d_krmn_order[k][r][need_check[begin_idx].second];
//         vector<int> order;
//         for (auto x : ori_order) {
//             if (idx_revers[x] != -1 && !isFilled[idx_revers[x]]) {
//                 // double need_p = pow(2.0, (remain_TBS[idx_revers[x]] / W) - 1) / s0_krn[k][r][idx_revers[x]];
//                 // if (need_p < 1.0) {
//                     // debug(need_p);
//                     order.push_back(idx_revers[x]);
//                 // }
//             }
//         }
//         // debug(remain_TBS);
//         debug(order);
//         // debug(N);
//         vector<int> u;
//         for (int use = 0; use < order.size(); use++) {
//             debug(u);
//             u.push_back(order[use]);
//             vector<double> need(u.size());
//             for (int n = 0; n < u.size(); n++) {
//                 int now = u[n];
//                 double p = 0;
//                 FOR (m, u.size()) if (m != n) p += _d_krmn[k][r][order[u[m]]][order[u[n]]];
//                 p = exp(p);
//                 // debug(p);
//                 double ne = pow(2.0, remain_TBS[now] / W) - 1;
//                 // debug(ne);
//                 ne /= p;
//                 // debug(ne);
//                 ne /= s0_krn[k][r][now];
//                 // debug(ne);
//                 need[n] = ne + 1e-6;
//                 // debug(ne);
//             }
//             debug(u);
//             debug(need);
//             // debug(u);
//             // debug(std::accumulate(need.begin(), need.end(), 0.0));
//             if (std::accumulate(need.begin(), need.end(), 0.0) > 1.0) {
//                 u.pop_back();
//                 // debug(u);
//             }
//         }
 
//         debug(u);
//         vector<double> need(u.size());
//         for (int n = 0; n < u.size(); n++) {
//             int now = u[n];
//             double p = 0;
//             FOR (m, u.size()) if (m != n) p += _d_krmn[k][r][order[u[m]]][order[u[n]]];
//             p = exp(p);
//             double ne = pow(2.0, remain_TBS[now] / W) - 1;
//             ne /= p;
//             ne /= s0_krn[k][r][now];
//             need[n] = ne + 1e-6;
//             p_tkrn[t][k][r][need_check[now].second] = need[n];
//         }
//         debug(need);
//         debug(u);
//         for (auto u_ : u) {
//             isFilled[u_] = true;
//             filled_cnt ++;
//         }
//     }
//     debug(filled_cnt);
//     // exit(0);
//     // return 0;
//     return 0;
// }
 
vector<vector<int> > end_J;
 
void read_all_data() {
    read(N); read(K); read(T); read(R);
    s0_tkrn.resize(T, vector<vector<vector<double>>>(K, vector<vector<double>>(R, vector<double>(N))));
    s_tkrn.resize(T, vector<vector<vector<double>>>(K, vector<vector<double>>(R, vector<double>(N))));
    s_tkn.resize(T, vector<vector<double>>(K, vector<double>(N)));
    FOR(t, T) FOR(k, K) FOR(r, R) FOR(n, N) read(s0_tkrn[t][k][r][n]);
    d_krmn.resize(K, vector<vector<vector<double>>>(R, vector<vector<double>>(N, vector<double>(N))));
 
    FOR(k, K) FOR(r, R) FOR(m, N) FOR(n, N) read(d_krmn[k][r][m][n]);
    read(J);
    frameID.resize(J), TBS.resize(J), uID.resize(J), TTI0.resize(J), nTTI.resize(J);
    FOR (i, J) {
        read(frameID[i]);read(TBS[i]);read(uID[i]);read(TTI0[i]);read(nTTI[i]);
    }
    p_tkrn.resize(T, vector<vector<vector<double>>>(K, vector<vector<double>>(R, vector<double>(N, 0))));
    b_tkrn.resize(T, vector<vector<vector<int>>>(K, vector<vector<int>>(R, vector<int>(N, 0))));
    g.resize(J);
    index_tn.resize(T, vector<int>(N, -1));
    end_J.resize(T);
    for (int j = 0; j < J; j++) {
        for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
            index_tn[t][uID[j]] = j;
        }
        end_J[TTI0[j] + nTTI[j] - 1].push_back(j);
    }
    FOR (t, T) debug(index_tn[t]);
}
 
template <typename T>
struct hungarian {  // km
  int n;
  vector<int> matchx;  // 左集合对应的匹配点
  vector<int> matchy;  // 右集合对应的匹配点
  vector<int> pre;     // 连接右集合的左点
  vector<bool> visx;   // 拜访数组 左
  vector<bool> visy;   // 拜访数组 右
  vector<T> lx;
  vector<T> ly;
  vector<vector<T> > g;
  vector<T> slack;
  T inf;
  T res;
  queue<int> q;
  int org_n;
  int org_m;
 
  hungarian(int _n, int _m) {
    org_n = _n;
    org_m = _m;
    n = max(_n, _m);
    inf = numeric_limits<T>::max();
    res = 0;
    g = vector<vector<T> >(n, vector<T>(n));
    matchx = vector<int>(n, -1);
    matchy = vector<int>(n, -1);
    pre = vector<int>(n);
    visx = vector<bool>(n);
    visy = vector<bool>(n);
    lx = vector<T>(n, -inf);
    ly = vector<T>(n);
    slack = vector<T>(n);
  }
 
  void addEdge(int u, int v, T w) {
    g[u][v] = w;  // 负值还不如不匹配 因此设为0不影响
  }
 
  bool check(int v) {
    visy[v] = true;
    if (matchy[v] != -1) {
      q.push(matchy[v]);
      visx[matchy[v]] = true;  // in S
      return false;
    }
    // 找到新的未匹配点 更新匹配点 pre 数组记录着"非匹配边"上与之相连的点
    while (v != -1) {
      matchy[v] = pre[v];
      swap(v, matchx[pre[v]]);
    }
    return true;
  }
 
  void bfs(int i) {
    while (!q.empty()) {
      q.pop();
    }
    q.push(i);
    visx[i] = true;
    while (true) {
      while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v = 0; v < n; v++) {
          if (!visy[v]) {
            T delta = lx[u] + ly[v] - g[u][v];
            if (slack[v] >= delta) {
              pre[v] = u;
              if (delta) {
                slack[v] = delta;
              } else if (check(v)) {  // delta=0 代表有机会加入相等子图 找增广路
                                      // 找到就return 重建交错树
                return;
              }
            }
          }
        }
      }
      // 没有增广路 修改顶标
      T a = inf;
      for (int j = 0; j < n; j++) {
        if (!visy[j]) {
          a = min(a, slack[j]);
        }
      }
      for (int j = 0; j < n; j++) {
        if (visx[j]) {  // S
          lx[j] -= a;
        }
        if (visy[j]) {  // T
          ly[j] += a;
        } else {  // T'
          slack[j] -= a;
        }
      }
      for (int j = 0; j < n; j++) {
        if (!visy[j] && slack[j] == 0 && check(j)) {
          return;
        }
      }
    }
  }
 
  void solve(int j = -1) {
    if (j == -1) {
        j = n;
    }
    // 初始顶标
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        lx[i] = max(lx[i], g[i][j]);
      }
    }
 
    for (int i = 0; i < j; i++) {
      fill(slack.begin(), slack.end(), inf);
      fill(visx.begin(), visx.end(), false);
      fill(visy.begin(), visy.end(), false);
      bfs(i);
    }
    // custom
    for (int i = 0; i < n; i++) {
      if (g[i][matchx[i]] > 10000) {
        res += g[i][matchx[i]];
      } else {
        matchx[i] = -1;
      }
    }
 
  }
};
void solve_bipart_final2(vector<double> &remain_TBS_) {
    vector<double> remain_TBS = remain_TBS_;
    vector<vector<int> > use_tr(T, vector<int>(R));
    vector<int> remain_r_cnt_t(T, R);
    vector<vector<double>> remain_p_tk(T, vector<double>(K, R));
    FOR (t, T) FOR(k, K) FOR(r, R) FOR(n, N) {
        if (p_tkrn[t][k][r][n] > 0) {
            use_tr[t][r] = true;
            remain_p_tk[t][k] -= p_tkrn[t][k][r][n];
            // debug(remain_p_tk[t][k]);
        }
    }
    FOR (t, T) FOR(r, R) {
        if (use_tr[t][r]) {
            remain_r_cnt_t[t]--;
        }
    }
 
    vector<vector<int> > list_n_t(T);
    vector<vector<vector<int> > > trans_jtr(J);
    vector<vector<pair<double, pair<int, int> > > > vpp(T);
    FOR (j, J) {
        assert (remain_TBS[j] <= 0 || remain_TBS[j] == TBS[j]);
        if (remain_TBS[j] <= 0) continue;
        trans_jtr[j].resize(nTTI[j]);
        for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
            trans_jtr[j][t - TTI0[j]].resize(R);
            list_n_t[t].push_back(j);
            FOR (r, R) {
                if (use_tr[t][r]) continue;
                double trans_all = 0;
                FOR (k, K) {
                    double trans = log2(1 + s0_tkrn[t][k][r][uID[j]] * min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6))) * W;
                    trans_all += trans;
                }
                trans_jtr[j][t - TTI0[j]][r] = trans_all;
                vpp[t].push_back({trans_all / TBS[j], {j, r}});
                // debug(make_pair(j, t));
                // debug(make_pair(r, TBS[j]));
                // debug(trans_all);
            }
        }
    }
    vector<double> acc_score(J);
    vector<vector<vector<double> > > mul_jk(J);
    vector<vector<vector<int> > > cnt_jk(J);
 
    FOR (j, J) {
        mul_jk[j].resize(K, vector<double>(nTTI[j], 1));
        cnt_jk[j].resize(K, vector<int>(nTTI[j], 0));
    }
 
    FOR (t, T) {
        sort(vpp[t].begin(), vpp[t].end());
        reverse(vpp[t].begin(), vpp[t].end());
    }
 
    vector<int> cnt_t(T);
    vector<vector<int> > j_list_t(T);
    FOR (j, J) {
        if (remain_TBS[j] <= 0) continue;
        for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
            cnt_t[t]++;
            j_list_t[t].push_back(j);
        }
    }
    set<pair<int, int>> p_cnt_t;
    FOR (t, T) {
        p_cnt_t.insert({cnt_t[t], t});
    }
 
    vector<int> t_counter_j = nTTI;
 
    // reverse(vpp.begin(), vpp.end());
    while (p_cnt_t.size()) {
        pair<int, int> cnt_t1 = *p_cnt_t.begin();
        int t = cnt_t1.second;
        debug(t + 0);
        p_cnt_t.erase(cnt_t1);
        for (int i = 0; i < vpp[t].size(); i++) {
            auto [score, jr] = vpp[t][i];
            auto [j, r] = jr;
            if (acc_score[j] > TBS[j]) continue;
            if (use_tr[t][r]) continue;
            use_tr[t][r] = 1;
            double s_before = acc_score[j];
            double d_score = 0;
            FOR (k, K) {
                double ori_score = cnt_jk[j][k][t - TTI0[j]] == 0 ? 0 : (
                    cnt_jk[j][k][t - TTI0[j]] * log2(1 + pow(mul_jk[j][k][t - TTI0[j]], 1.0 / cnt_jk[j][k][t - TTI0[j]])) * W
                );
                cnt_jk[j][k][t - TTI0[j]] += 1;
                mul_jk[j][k][t - TTI0[j]] *= s0_tkrn[t][k][r][uID[j]] * min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6));
                double new_score = cnt_jk[j][k][t - TTI0[j]] == 0 ? 0 : (
                    cnt_jk[j][k][t - TTI0[j]] * log2(1 + pow(mul_jk[j][k][t - TTI0[j]], 1.0 / cnt_jk[j][k][t - TTI0[j]])) * W
                );
                // debug(new_score);
                // debug(ori_score);
                d_score += log2(1 + min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6)) * s0_tkrn[t][k][r][uID[j]]) * W;
                p_tkrn[t][k][r][uID[j]] = min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6));
                acc_score[j] += new_score - ori_score;
            }
            cerr << t << " " << j << " "; debug(TBS[j]); debug(d_score); debug(acc_score[j]);
            cerr << acc_score[j] - s_before << endl;
            if (acc_score[j] > TBS[j]) {
                for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
                    if (p_cnt_t.count({cnt_t[t], t})) {
                        p_cnt_t.erase({cnt_t[t], t});
                        cnt_t[t]--;
                        p_cnt_t.insert({cnt_t[t], t});
                    }
                }
            }
        }
        for (int j : j_list_t[t]) {
            t_counter_j[j]--;
        }
        for (int j : j_list_t[t]) {
            if (t_counter_j[j] == 0 && acc_score[j] < TBS[j]) {
                cerr << "failed " << j << " " << acc_score[j] << " " << TBS[j] << endl;
                acc_score[j] = 1e9;
                FOR (r, R) FOR (k, K) for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
                    if (p_tkrn[t][k][r][uID[j]] > 0) {
                        use_tr[t][r] = 0;
                        p_tkrn[t][k][r][uID[j]] = 0;
                    }
                }
                for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
                    cnt_t[t]--;
                    for (int j : j_list_t[t]) {
                        t_counter_j[j]++;
                    }
                    p_cnt_t.insert({cnt_t[t], t});
                }
 
            }
        }
    }
    // FOR (j, J) {
    //     cout << remain_TBS[j] - acc_score[j] << endl;
    // }
 
    // FOR (t, T) {
    //     debug(list_n_t[t]);
    //     debug(use_tr[t]);
    // }
    // FOR (j, J) {
    //     if (trans_jtr[j].size() == 0) continue;
    //     cerr << TBS[j] << " "; debug(j);
    //     FOR (tti, nTTI[j]) {
    //         cerr << tti + TTI0[j] << " "; debug(trans_jtr[j][tti]);
    //     }
    // }
 
}
 
void solve_bipart_final(vector<double> &remain_TBS) {
    vector<vector<int> > use_tr(T, vector<int>(R));
    vector<int> remain_r_cnt_t(T, R);
    vector<vector<double>> remain_p_tk(T, vector<double>(K, R));
    FOR (t, T) FOR(k, K) FOR(r, R) FOR(n, N) {
        if (p_tkrn[t][k][r][n] > 0) {
            use_tr[t][r] = true;
            remain_p_tk[t][k] -= p_tkrn[t][k][r][n];
            // debug(remain_p_tk[t][k]);
        }
    }
    FOR (t, T) FOR(r, R) {
        if (use_tr[t][r]) {
            remain_r_cnt_t[t]--;
        }
    }
 
    hungarian<double> h(J, T * R);
    FOR (j, J) {
        assert (remain_TBS[j] <= 0 || remain_TBS[j] == TBS[j]);
        if (remain_TBS[j] <= 0) continue;
        for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
            FOR (r, R) {
                if (use_tr[t][r]) continue;
                double trans_all = 0;
                FOR (k, K) {
                    // debug(remain_p_tk[t][k]);
                    // debug(((remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6)));
                    // debug(log2(1 + s0_tkrn[t][k][r][uID[J]] * ((remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6))));
                    // debug(trans);
                    double trans = log2(1 + s0_tkrn[t][k][r][uID[j]] * min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6))) * W;
                    trans_all += trans;
                }
                // debug(trans_all);
                // debug(TBS[j]);
                double score = 0;
                if (trans_all > TBS[j]) {
                    score = 10000 + TBS[j] / trans_all;
                } else {
                    score = trans_all / TBS[j];
                }
                h.addEdge(j, t * R + r, score);
 
                // if (trans_all > TBS[j]) {
                //     // cerr << j << " " << t << " " << r << endl;
                //     h.addEdge(j, t * R + r, 1);
                // }
            }
        }
    }
    h.solve(J);
    FOR (j, J) {
        if (h.matchx[j] >= 0) {
            int t = h.matchx[j] / R;
            int r = h.matchx[j] % R;
            double left = 1e-6, right = 1.0;
            while (right - left > 1e-3) {
                double mid = (right + left) / 2;
                double trans_all = 0;
                FOR (k, K) {
                    double trans = log2(1 + mid * min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6)) * s0_tkrn[t][k][r][uID[j]]) * W;
                    trans_all += trans;
                }
                if (trans_all > TBS[j]) {
                    right = mid;
                } else {
                    left = mid;
                }
            }
            double trans_all = 0;
            FOR (k, K) {
                double trans = log2(1 + right * min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6)) * s0_tkrn[t][k][r][uID[j]]) * W;
                trans_all += trans;
            }
            remain_TBS[j] -= trans_all;
            FOR (k, K) {
                p_tkrn[t][k][r][uID[j]] = right * min(4.0, (remain_p_tk[t][k]) / (remain_r_cnt_t[t] + 1e-6));
            }
        }
    }
    debug(vector<int>(h.matchx.begin(), h.matchx.begin() + J));
    // debug(TBS);
    // debug(remain_TBS);
}
 
 
void solve_bipart(int t, vector<pair<double, int> > need_check, vector<int> R_used, vector<double> remain_p_k) {
    if (need_check.size() == 0) return;
    // if (t == 4) {
    //     debug(R_used);
    //     debug(save);
    // }
    // return;
    int cnt_R_not_used = 0;
    FOR (r, R) {
        if (!R_used[r]) cnt_R_not_used++;
    }
 
    int NN = need_check.size();
    vector<double> remainTBS(NN);
    FOR (nn, NN) remainTBS[nn] = need_check[nn].first;
    vector<double> ori_remainTBS = remainTBS;
    vector<vector<double> > bipart(R, vector<double>(NN));
    FOR (r, R) {
        if (R_used[r]) continue;
        FOR (nn, NN) {
            double trans_all = 0;
            FOR (k, K) {
                double trans = log2(1 + s0_tkrn[t][k][r][need_check[nn].second] * ((remain_p_k[k]) / (cnt_R_not_used + 1e-6))) * W;
                trans_all += trans;
            }
            double score = 0;
            if (trans_all > need_check[nn].first) {
                score = 10000 + need_check[nn].first / trans_all;
            } else {
                score = trans_all / need_check[nn].first;
            }
            // trans_all / need_check[nn].first;
 
            bipart[r][nn] = score;
        }
    }
    // FOR (r, R) {
    //     debug(bipart[r]);
    // }
    hungarian<double> h(R, NN);
 
    FOR (r, R) FOR (nn, NN) {
        h.addEdge(r, nn, bipart[r][nn]);
    }
    h.solve();
    debug(h.matchx);
    vector<int> filled_need(NN);
    double use_r = 0;
    double use_r_cnt = 0;
    vector<double> remain_p_k_after = remain_p_k;
    FOR (r, R) {
        if (h.matchx[r] != -1) {
            R_used[r] = true;
            int nn = h.matchx[r];
            filled_need[nn] = true;
 
            double left = 1e-6, right = 1.0;
            while (right - left > 1e-3) {
                double mid = (right + left) / 2;
                double trans_all = 0;
                FOR (k, K) {
                    double trans = log2(1 + mid * ((remain_p_k[k]) / (cnt_R_not_used + 1e-6)) * s0_tkrn[t][k][r][need_check[nn].second]) * W;
                    trans_all += trans;
                }
                if (trans_all > remainTBS[nn]) {
                    right = mid;
                } else {
                    left = mid;
                }
            }
            use_r += right;
            use_r_cnt++;
            FOR (k, K) {
                p_tkrn[t][k][r][need_check[nn].second] = right * ((remain_p_k[k]) / (cnt_R_not_used + 1e-6));
                remain_p_k_after[k] -= right * ((remain_p_k[k]) / (cnt_R_not_used + 1e-6));
            }
        }
    }
    // double remain_p = 0;
    int remain_p_cnt = 0;
    for (int i = 0; i < R; i++) {
        if (!R_used[i]) {
            // remain_p += 1;
            remain_p_cnt++;
        }
    }
 
    // remain_p += use_r_cnt - use_r;
    // remain_p += save;
    // remain_p /= (1e-6 + remain_p_cnt);
    // remain_p = min(4.0, remain_p);
    // debug(remain_p);
 
    vector<pair<double, pair<int, int>>> vpp;
    FOR (r, R) {
        if (R_used[r]) continue;
        FOR (nn, NN) {
            if (filled_need[nn]) continue;
            vpp.push_back({bipart[r][nn], {r, nn}});
        }
    }
    vector<vector<double> > mul_kNN(K, vector<double>(NN, 1)), cnt_kNN(K, vector<double>(NN, 0));
    vector<vector<vector<double> >> mul_kNN_val(K, vector<vector<double>>(NN));
    sort(vpp.begin(), vpp.end(), greater<pair<double, pair<int, int>>>());
    for (int i = 0; i < vpp.size(); i++) {
        auto [score, r_nn] = vpp[i];
        auto [r, nn] = r_nn;
        if (R_used[r]) continue;
        if (filled_need[nn]) continue;
 
        FOR (k, K) {
            p_tkrn[t][k][r][need_check[nn].second] = min(4.0, remain_p_k_after[k] / (1e-6 + remain_p_cnt));
        }
        R_used[r] = true;
 
        // remainTBS[nn] -= ori_remainTBS[nn] * score;
        // if (remainTBS[nn] < 0) {
        //     filled_need[nn] = true;
        // }
 
        double s = 0;
        FOR (k, K) {
            mul_kNN[k][nn] *= s0_tkrn[t][k][r][need_check[nn].second] * min(4.0, remain_p_k_after[k] / (1e-6 + remain_p_cnt));
            mul_kNN_val[k][nn].push_back(s0_tkrn[t][k][r][need_check[nn].second] * min(4.0, remain_p_k_after[k] / (1e-6 + remain_p_cnt)));
            debug(k);
            debug(mul_kNN_val[k][nn]);
            cnt_kNN[k][nn]++;
            s += cnt_kNN[k][nn] * log2(1 + pow(mul_kNN[k][nn], 1.0 / cnt_kNN[k][nn])) * W;
        }
        debug(need_check[nn].second);
        debug(remainTBS[nn]);
        debug(s);
        if (remainTBS[nn] < s) {
            filled_need[nn] = true;
        }
    }
    // if (t == 4) { debug(save);
    //     exit(0);
    // }
}
 
void solve() {
    #ifdef DEBUG
    // freopen("/Users/zaq1xsw2tktk/Downloads/tests/06", "r", stdin);
    // freopen("/Users/zaq1xsw2tktk/Downloads/tests/51", "r", stdin);
    // freopen("/Users/zaq1xsw2tktk/Downloads/tests/21", "r", stdin);
    // freopen("21", "r", stdin);
    #endif
    read_all_data();
    // FOR (j, J) {
    //     double s = 0;
    //     for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
    //         clear(t);
    //         FOR (r, R) {
    //             FOR (k, K) {
    //                 p_tkrn[t][k][r][uID[j]] = s0_tkrn[t][k][r][uID[j]];
    //             }
    //         }
    //         calc(t);
    //         FOR (r, R) {
    //             FOR (k, K) {
    //                 s += log2(1 + s_tkn[t][k][uID[j]]) * W;
    //             }
    //         }
    //     }
    //     debug(TBS[j]);
    //     debug(s);
    // }
    // exit(0);
 
 
    FOR (t, T) clear(t);
    d_krmn_order.resize(K, vector<vector<vector<int>>>(R, vector<vector<int>>(N, vector<int>())));
    FOR(k, K) FOR(r, R) FOR(m, N) {
        // debug(k);
        // debug(r);
        vector<pair<double, int> > vp;
        FOR(n, N) {
            if (n == m) {
                vp.push_back({1, n});
            } else {
                vp.push_back({d_krmn[k][r][m][n], n});
            }
        }
        sort(vp.begin(), vp.end(), greater<pair<double, int>>());
        // debug(vp);
        // debug(vp.size());
        FOR (n, N) {
            // debug(n);
            d_krmn_order[k][r][m].push_back(vp[n].second);
        }
        // debug(d_krmn_order[k][r][m]);
    }
 
    vector<vector<double> > remain_TBS_by_t(T + 1);
    vector<double> remain_TBS = TBS;
    remain_TBS_by_t[0] = remain_TBS;
    // debug(TBS);
 
    vector<int> cnt_t(T);
    vector<vector<int> > j_list_t(T);
    FOR (j, J) {
        for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
            cnt_t[t]++;
            j_list_t[t].push_back(j);
        }
    }
    set<pair<int, int>> p_cnt_t;
    FOR (t, T) {
        p_cnt_t.insert({cnt_t[t], t});
    }
 
    vector<int> t_counter_j = nTTI;
 
    while (p_cnt_t.size()) {
        pair<int, int> cnt_t1 = *p_cnt_t.begin();
        int t = cnt_t1.second;
        p_cnt_t.erase(cnt_t1);
        debug(t);
        // remain_TBS = remain_TBS_by_t[t];
        clear(t);
        // debug(remain_TBS);
        vector<pair<double, int> > need_check;
        // vector<pair<pair<double, double>, int> > need_check2;
        FOR (n, N) {
            if (index_tn[t][n] == -1 || remain_TBS[index_tn[t][n]] <= 0) {
                continue;
            }
            double max_s0_tkrn = 0;
            FOR (r, R) FOR (k, K) {
                max_s0_tkrn = max(max_s0_tkrn, s0_tkrn[t][k][r][index_tn[t][n]]);
            }
            need_check.push_back({remain_TBS[index_tn[t][n]], n});
            // need_check2.push_back({{remain_TBS[index_tn[t][n]] / log2(1 + max_s0_tkrn), remain_TBS[index_tn[t][n]]}, n});
            // need_check2.push_back({{max_s0_tkrn, remain_TBS[index_tn[t][n]]}, n});
        }
        sort(need_check.begin(), need_check.end());
        
        if (need_check.size() == 0) continue;
        for (int i = 0; i < need_check.size(); i++) {
            cerr << need_check[i].first << "," << need_check[i].second << "," << index_tn[t][need_check[i].second] << (i != need_check.size() - 1 ? " " : "\n");
        }
 
        vector<double> mul_k = vector<double>(K, 1);
        vector<int> cnt_k = vector<int>(K, 0);
 
        // solve_one_t(t, s0_tkrn[t], d_krmn, need_check);
 
        int max_filled = 0;
        int max_filled_k = 0;
        double save = 0;
        double _ = 0;
        vector<int> used_R(R);
        vector<double> p_remain_k(K, R);
        while (true) {
            bool all_r_used = true;
            FOR (r, R) if (!used_R[r]) all_r_used = false;
            if (all_r_used) break;
            debug(need_check);
            double best_score = 0;
            int best_r = -1;
            vector<int> best_u;
            int k = std::max_element(p_remain_k.begin(), p_remain_k.end()) - p_remain_k.begin();
            int out_k = 0;
            FOR (r, R) {
                if (used_R[r]) continue;
                // debug(r);
                // double tt = clock();
                FOR (k, K) {
                    if (((double)(clock()) / CLOCKS_PER_SEC) < 1.5) {
                        auto [max_score_s, u_s] = solve_one_t_list_v2(t, r, need_check, N, k, _, p_remain_k[k]);
                        // debug((double)(clock() - tt) / CLOCKS_PER_SEC);
                        for (int i = need_check.size(); i >= 2; i--) {
                            if (u_s[i].size()) {
                                if (u_s[i].size() > best_u.size() || u_s[i].size() == best_u.size() && max_score_s[i] < best_score) {
                                    best_score = max_score_s[i];
                                    best_u = u_s[i];
                                    best_r = r;
                                    out_k = k;
                                }
                                debug(u_s[i]);
                                break;
                            }
                        }
                    } else {
                        auto [max_score_s, u_s] = solve_one_t_list(t, r, need_check, N, k, _, p_remain_k[k]);
                        // debug((double)(clock() - tt) / CLOCKS_PER_SEC);
                        for (int i = need_check.size(); i >= 2; i--) {
                            if (u_s[i].size()) {
                                if (u_s[i].size() > best_u.size() || u_s[i].size() == best_u.size() && max_score_s[i] < best_score) {
                                    best_score = max_score_s[i];
                                    best_u = u_s[i];
                                    best_r = r;
                                    out_k = k;
                                }
                                debug(u_s[i]);
                                break;
                            }
                        }
                    }
                }
            }
            k = out_k;
            int r = best_r;
 
            
            bool find_extra = false;
            if (best_u.size() <= 1 && false) {
                best_u = {0};
                FOR (r_, R) if (!used_R[r_]) r = r_;
                int remain_r_cnt = 0;
                FOR (r_, R) if (!used_R[r_]) remain_r_cnt += 1;
                debug(best_u);
                vector<int> used_idx(need_check.size());
                for (int i = 0; i < best_u.size(); i++) used_idx[best_u[i]] = 1;
                for (int i = 0; i < need_check.size(); i++) {
                    if (used_idx[i]) continue;
                    vector<double> s_use(K);
                    FOR(k_, K) {
                        if (k == k_) {continue;}
                        s_use[k_] = s0_tkrn[t][k_][r][need_check[i].second] * p_remain_k[k_] / (1 + remain_r_cnt + 1e-6);
                    }
                    vector<double> need_p(best_u.size());
                    for (int j = 0; j < best_u.size(); j++) {
                        int now = best_u[j];
                        double fm = 1;
                        FOR (k_, K) {
                            if (k == k_) {continue;}
                            fm += s_use[k_] * exp(-d_krmn[k_][r][need_check[i].second][need_check[best_u[j]].second]);
                        }
                        double p = 0;
                        FOR (m, best_u.size()) if (m != j) p += d_krmn[k][r][need_check[best_u[j]].second][need_check[best_u[m]].second];
                        p = exp(p);
                        double ne = pow(2.0, need_check[now].first / W) - 1;
                        ne /= p;
                        ne /= s0_tkrn[t][k][r][need_check[now].second];
                        ne *= fm;
                        ne += 1e-6;
                        need_p[j] = ne;
                        // debug(fm);
                    }
                    double new_score = 0;
                    FOR (k_, K) {
                        if (k == k_) continue;
                        double fm = 1;
                        for (int j = 0; j < best_u.size(); j++) {
                            fm += need_p[j] * s0_tkrn[t][k][r][need_check[best_u[j]].second] * exp(-d_krmn[k_][r][need_check[i].second][need_check[best_u[j]].second]);
                        }
                        double fz = s_use[k_];
                        new_score += W * log2(1 + fz / fm);
                    }

                    debug(need_p);
                    debug(new_score);
                    debug(need_check[i].first);
                    debug(t);
                    if (std::accumulate(need_p.begin(), need_p.end(), 0.0) < min(4.0, p_remain_k[k]) && new_score > need_check[i].first) {
                        find_extra = true;
                        p_remain_k[k] -= std::accumulate(need_p.begin(), need_p.end(), 0.0);
                        FOR (k_, K) {
                            if (k_ == k) continue;
                            p_remain_k[k_] -= p_remain_k[k_] / (1 + remain_r_cnt + 1e-6);
                        }
                        debug(need_p);
                        debug(new_score);
                        debug(need_check[i].first);
                        debug(t);
                        for (int n = 0; n < best_u.size(); n++) {
                            int now = best_u[n];
                            p_tkrn[t][k][r][need_check[now].second] = need_p[n];
                        }
                        FOR (k_, K) {
                            if (k == k_) continue;
                            p_tkrn[t][k_][r][need_check[i].second] = p_remain_k[k_] / (1 + remain_r_cnt + 1e-6);
                        }
                        vector<pair<double, int> > new_need_check;
                        for (int n = 0; n < need_check.size(); n++) {
                            if (used_idx[n] || n == i) continue;
                            new_need_check.push_back(need_check[n]);
                        }
                        new_need_check.swap(need_check);
                        break;
                    }
                }
            }
            used_R[r] = 1;
            if (find_extra) continue;
            // TRY
            if (best_u.size() <= 1) break;
            vector<pair<double, int> > new_need_check;
            vector<int> use_need_check(need_check.size(), false);
            vector<int> u = best_u;
            for (auto x : u) use_need_check[x] = true;
            for (int i = 0; i < need_check.size(); i++) {
                if (use_need_check[i]) {
                    continue;
                }
                new_need_check.push_back(need_check[i]);
            }
            debug(p_remain_k[k]);
            for (int n = 0; n < u.size(); n++) {
                int now = u[n];
                double p = 0;
                FOR (m, u.size()) if (m != n) p += d_krmn[k][r][need_check[u[n]].second][need_check[u[m]].second];
                p = exp(p);
                double ne = pow(2.0, need_check[now].first / W) - 1;
                ne /= p;
                ne /= s0_tkrn[t][k][r][need_check[now].second];
                ne += 1e-6;
                p_tkrn[t][k][r][need_check[now].second] = ne;
                debug(ne);
            }
            new_need_check.swap(need_check);
            p_remain_k[k] -= best_score;
            debug(best_u);
            debug(best_score);
        }
 
        auto before = remain_TBS;
        double s = 0;
        calc(t);
        FOR (n, N) {
            if (index_tn[t][n] == -1 || remain_TBS[index_tn[t][n]] <= 0) {
                continue;
            }
            int j = index_tn[t][n];
            double use = 0;
            if (t == 74) debug(remain_TBS[j] + 0);
            FOR(k, K) FOR(r, R) {
                s += W * b_tkrn[t][k][r][n] * log2(1 + s_tkn[t][k][n]);
                remain_TBS[j] -= W * b_tkrn[t][k][r][n] * log2(1 + s_tkn[t][k][n]);
                // if (W * b_tkrn[t][k][r][n] * log2(1 + s_tkn[t][k][n]) != 0) {
                //     debug(W * b_tkrn[t][k][r][n] * log2(1 + s_tkn[t][k][n]));
                // }
            }
            if (t == 74) debug(remain_TBS[j]);
        }
        // debug(s);
        FOR (j, J) before[j] -= remain_TBS[j];
        remain_TBS_by_t[t + 1] = remain_TBS;
        // continue;
        // if (((double)clock())/(double)CLOCKS_PER_SEC > 1.700) {
        //     continue;
        // }
        // bool has_end = false;
        // int min_t = 1e9;
        // for (auto j : end_J[t]) {
        //     if (remain_TBS[j] > 0 && remain_TBS[j] != TBS[j] && remain_TBS[j] < 1e9  ) {
        //         debug(j);
        //         debug(remain_TBS[j]);
        //         debug(TBS[j]);
        //         has_end = true;
        //         min_t = min(min_t, TTI0[j]);
        //         // cerr << "bad " << j << " " << remain_TBS[j] << " " << TBS[j] << endl;
        //     }
        // }
        // if (has_end) {
        //     int new_t = min_t - 1;
        //     for (auto j : end_J[t]) {
        //         if (remain_TBS[j] > 0 && remain_TBS[j] != TBS[j] && remain_TBS[j] < 1e9) {
        //             remain_TBS_by_t[new_t + 1][j] += 1e10;
        //             TBS[j] += 1e10;
        //             debug(remain_TBS_by_t[new_t + 1][j]);
        //             // min_t = min(min_t, TTI0[j]);
        //             // cerr << "bad " << j << " " << remain_TBS[j] << " " << TBS[j] << endl;
        //         }
        //     }
        //     t = new_t;
        // }
        
        // debug(before);
    }
    solve_bipart_final(remain_TBS);
    // solve_bipart_final(remain_TBS);
    FOR (t, T) {
        vector<int> v;
        vector<int> v2;
        FOR (n, N) {
            if (index_tn[t][n] == -1) {
                continue;
            }
            v2.push_back(index_tn[t][n]);
            if (remain_TBS[index_tn[t][n]] <= 0) {
                continue;
            }
            v.push_back(index_tn[t][n]);
        }
        cerr << "t: " << t << " "; debug(v);
        debug(v2);
    }
 
    solve_bipart_final2(remain_TBS);
 
    // FOR(t, T) FOR(k, K) FOR(r, R) {
    //     // if (k != 0) continue;
    //     FOR (n, N) {
    //         if (n == r) {
    //             p_tkrn[t][k][r][n] = 1.0;
    //         }
    //     }
    // }
 
    // FOR (t, T) FOR (k, K) {
    //     double t1 = 0;
    //     FOR (r, R) {
    //         double t2 = 0;
    //         FOR (n, N) {
    //             t1 += p_tkrn[t][k][r][n];
    //             t2 += p_tkrn[t][k][r][n];
    //         }
    //         if (t2) {
    //             // debug(t2);
    //             assert(t2 <= 4);
    //         }
    //     }
    //     if (t1) {
    //         debug(t);
    //         debug(R - t1);
    //         assert(t1 <= R);
    //         // debug(t1);
    //     }
    // }
 
    // calc();
    // int cc = 0;
    // FOR (j, J) if (TBS[j] == 2) cc++;
    // debug(cc);
    // FOR(j, J) {
    //     if (g[j] >= TBS[j] && cc > 0) {
    //         cc--;
    //         for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
    //             FOR (k, K) FOR (r, R) {
    //                 p_tkrn[t][k][r][uID[j]] = 0;
    //             }
    //         }
    //     }
    // }
    // bool has_gt2 = false;
    // FOR (j, J) {
    //     if (nTTI[j] > 1) has_gt2 = true;
    // }
    // calc();
    // int cc = (J) / 2;
    // if (has_gt2) FOR(j, J) {
    //         for (int t = TTI0[j]; t < TTI0[j] + nTTI[j]; t++) {
    //             FOR (k, K) FOR (r, R) {
    //                 p_tkrn[t][k][r][uID[j]] = 0;
    //             }
    //         }
    // }
}
 
int main() {
    ios_base::sync_with_stdio(false); std::cin.tie(0);
    pre_init();
    int t = 1;
    cout << setprecision(9);
    // cin >> t;
    solve();
 
    #ifndef DEBUG
    FOR(t, T) FOR(k, K) FOR(r, R) {
        print (p_tkrn[t][k][r]);
    }
    #endif
 
    #ifdef DEBUG
    calc();
    #endif
    debug(J);
 
    cerr<<"Time:"<<1000*((double)clock())/(double)CLOCKS_PER_SEC<<"ms\n";
 
    // FOR (k, K) FOR (r, R) {
    //     cerr << k << " " << r << " "; debug(p_tkrn[7][k][r]);
    //     FOR (n, N) {
    //         if (p_tkrn[7][k][r][n]) {
    //             cerr << n << " ";
    //         }
    //     }
    //     cerr << endl;
    // }
}
// 