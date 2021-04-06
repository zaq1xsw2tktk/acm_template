#include <bits/stdc++.h>
using namespace std;
#define FOR(i, n) for (int (i) = 0; (i) < (n); (i)++)
#define FORI(i, a, b) for (int (i) = (a); (i) < (b); (i)++)
 
#define ll long long
#define mp(m, n) make_pair((m), (n))
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

#define debug(x) _debug(#x, (x))
#else
#define debug(x)
#endif
 
#define db debug

void solve();

int main() {
    ios_base::sync_with_stdio(false); std::cin.tie(0);
    int t = 1;
    // cin >> t;
    while (t--) {
        solve();
    }
}

void solve() {

}

