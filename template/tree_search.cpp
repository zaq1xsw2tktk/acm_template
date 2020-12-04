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