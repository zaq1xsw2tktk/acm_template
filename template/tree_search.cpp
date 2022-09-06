#include <bits/stdc++.h>
using namespace std;
struct Tree {
    vector<vector<int> > edges;
    int n;
    Tree(const vector<vector<int>>& edges) : edges(edges), n(edges.size()) {}
    void bfs(vector<int> &distance, int node, int fnode=-1) {
        for (int i = 0; i < edges[node].size(); ++i) {
            if (edges[node][i] != fnode) {
                distance[edges[node][i]] = distance[node] + 1;
                bfs(distance, edges[node][i], node);
            }
        }
    }
    vector<int> get_distance(int fnode) {
        vector<int> distance(n);
        distance[fnode] = 0;
        bfs(distance, fnode);
        return distance;
    }
    int get_d(int* st = nullptr, int* ed = nullptr) {
        vector<int> dis = get_distance(0);
        int index = std::max_element(dis.begin(), dis.end()) - dis.begin();
        vector<int> dis2 = get_distance(index);
        int ed_index = std::max_element(dis2.begin(), dis2.end()) - dis2.begin();
        if (st != nullptr) {
            *st = index;
        }
        if (ed != nullptr) {
            *ed = ed_index;
        }
        return dis2[ed_index];
    }
};

// lca 模板
class Tree {
public:
    vector<vector<int> > e;
    vector<int> depth;
    vector<vector<int> > fnode;
    int root;
    int n;
    int get_fnode(int x, int t) {
        for (int i = 0; i < 20; i++) {
            if (t & (1 << i)) {
                x = fnode[i][x];
            }
        }
        return x;
    }
    int lca(int x, int y) {
        if (depth[x] < depth[y]) {
            swap(x, y);
        }
        int dpeth_diff = depth[x] - depth[y];
        x = get_fnode(x, dpeth_diff);
        if (x == y) return x;
        for (int i = 19; i >= 0; i--) {
            if (fnode[i][x] != fnode[i][y]) {
                x = fnode[i][x];
                y = fnode[i][y];
            }
        }
        return fnode[0][x];
    }
    void init() {
        auto dfs = [] (int x, int fnode, vector<int> &depth, vector<vector<int> > &e, vector<int> &f, auto dfs) -> void {
            for (int i = 0; i < e[x].size(); i++) {
                if (e[x][i] == fnode) {
                    continue;
                }
                f[e[x][i]] = x;
                depth[e[x][i]] = depth[x] + 1;
                dfs(e[x][i], x, depth, e, f, dfs);
            }
        };
        depth.resize(n);
        fnode.resize(20, vector<int>(n, -1));
        fnode[0][root] = 0;
        dfs(root, -1, depth, e, fnode[0], dfs);
        for (int i = 1; i < 20; i++) {
            for (int j = 0; j < n; j++) {
                fnode[i][j] = fnode[i - 1][fnode[i - 1][j]];
            }
        }
    }
    Tree(vector<pair<int, int> > edges, int root=0) {
        this->root = root;
        n = edges.size() + 1;
        e.resize(n);
        for (int i = 0; i < edges.size(); i++) {
            auto [x, y] = edges[i];
            e[x].push_back(y);
            e[y].push_back(x);
        }
        init();
    }
};