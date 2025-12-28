#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <string>
#include <limits>
#include <utility>
#include <algorithm>
#include <iomanip>
#include <filesystem>

using namespace std;
using std::filesystem::u8path;

struct BFS {
    vector<int> dist;
    vector<int> parent;
    vector<int> order;
};

BFS bfs(const vector<vector<int>>& adj, int s) {
    const int V = (int)adj.size();
    const int INF = numeric_limits<int>::max();
    BFS out{vector<int>(V, INF), vector<int>(V,-1), {}};

    queue<int> q;
    out.dist[s] = 0;
    q.push(s);

    // Iterate over queue
    while(!q.empty()) {
        //Dequeue
        int u = q.front();
        q.pop();
        out.order.push_back(u);
        
        for(int x : adj[u]) {
            if(out.dist[x] == INF) {
                out.dist[x] = out.dist[u] + 1;
                out.parent[x] = u;
                q.push(x);
            }
        } 
    }
    return out;
}

int main() {
    ifstream file(u8path("ø5g1.txt"));
    if(!file) {
        cout << "Feil ved åpning av ø5g1.txt" << endl;
        return EXIT_FAILURE;
    }

    int N, M;
    if(!(file >> N >> M)) {
        cout << "Ugyldig filformat: forventet 'N M' på første linje.\n";
        return 1;
    }

    vector<vector<int>> adj(N);
    for(int i = 0; i < M; i++) {
        int u, v;
        if(!(file >> u >> v)) {
            cout << "Forventet" << M << "kanter, men fant færre.\n";
            return 1;
        }
        if(u < 0 || u >= N || v < 0 || v >= N) {
            cout << "Kant utenfor område på linje " << (i + 2) << ".\n";
            return 1;
        }
        adj[u].push_back(v);
    }

    for(auto& nbrs : adj) sort(nbrs.rbegin(), nbrs.rend());     // descending

    cout << "Startnode (default 0): ";
    string line;
    getline(cin >> ws, line);
    int s = 0;
    if(!line.empty()) s = stoi(line);
    if(s < 0 || s >= N) {
        cout << "Startnode utenfor område 0.." << (N-1) << "\n";
        return 1;
    }

    BFS res = bfs(adj, s);
    cout << left << setw(6) << "Node"
        << setw(8) << "Forgj"
        << setw(6) << "Dist" << "\n";
    for(int v = 0; v < N; ++v) {
    string parent = (res.parent[v] == -1) ? "-" : to_string(res.parent[v]);
    string dist   = (res.dist[v] == numeric_limits<int>::max()) ? "-" : to_string(res.dist[v]);
    
    cout << left << setw(6) << v
         << setw(8) << parent
         << setw(6) << dist << "\n";
    }

    return 0;
}