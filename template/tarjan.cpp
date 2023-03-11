void tarjan(int x, int fnode, vector<vector<int> > &edges, vector<int> &visited, vector<int> &depth, deque<int> &q, vector<vector<int> > &conn_group) {
    int val = depth[x] = q.size();
    q.push_back(x);
    visited[x] = 1;
    for (int i = 0; i < edges[x].size(); i++) {
        int t = edges[x][i];
        if (t == fnode) {
            continue;
        }
        if (visited[t]) {
            depth[x] = min(depth[x], depth[t]);
        } else {
            tarjan(t, x, edges, visited, depth, q, conn_group);
            depth[x] = min(depth[x], depth[t]);
        }
    }
    if (depth[x] == val) {
        vector<int> v;
        while(q.size()) {
            int t = q.back();
            v.push_back(t);
            q.pop_back();
            if (t == x) break;
        }
        debug(v);
        conn_group.push_back(v);
    }
}