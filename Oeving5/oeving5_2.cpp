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
};

BFS bfs(const vector<int>& offsets, const vector<int>& to, int n, int s) {
    BFS out{vector<int>(n, -1), vector<int>(n,-1)};

    queue<int> q;
    out.dist[s] = 0;
    q.push(s);

    // Iterate over queue
    while(!q.empty()) {
        //Dequeue
        int u = q.front();
        q.pop();        
        for(int i = offsets[u]; i < offsets[u+1]; i++) {
            int v = to[i];
            if(out.dist[v] == -1) {
                out.dist[v] = out.dist[u] + 1;
                out.parent[v] = u;
                q.push(v);
            }
        } 
    }
    return out;
}

vector<int> topo_kahn(const vector<int>& offsets, const vector<int>& to, int N, vector<int> indeg) {
    queue<int> q;
    for(int u = 0; u < N; u++) if(indeg[u] == 0) q.push(u);
    vector<int> order; 
    order.reserve(N);
    while (!q.empty()) {
        int u = q.front(); q.pop();
        order.push_back(u);
        for (int i = offsets[u]; i < offsets[u+1]; ++i) {
            int v = to[i];
            if (--indeg[v] == 0) q.push(v);
        }
    }
    if((int)order.size() != N) return {};   // cycle 
    return order;
}

int main() {

    ifstream file(u8path("ø5Skandinavia.txt"));
    if(!file) {
        cout << "Feil ved åpning av ø5Skandinavia.txt" << endl;
        return EXIT_FAILURE;
    }

    int N, M;
    if(!(file >> N >> M)) {
        cout << "Ugyldig filformat: forventet 'N M' på første linje.\n";
        return 1;
    }

    vector<pair<int,int>> edges;
    edges.reserve(M);
    
    vector<int> outdeg(N, 0), indeg(N, 0);
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
        edges.emplace_back(u,v);
        outdeg[u]++;
        indeg[v]++;
    }

    //  CSR offsets
    vector<int> offsets(N+1, 0);
    for(int u = 0; u < N; u++) offsets[u+1] = offsets[u] + outdeg[u];

    vector<int> cursor = offsets;
    vector<int> to(M);
    for(const auto& e : edges) {
        int u = e.first; int v = e.second;
        to[cursor[u]++] = v;
    }
    edges.clear(); edges.shrink_to_fit();

    vector<int> topo = topo_kahn(offsets, to, N, indeg);
    if (topo.empty()) {
        cout << "Ikke DAG (inneholder en sykel)\n";
        return 0;
    }

    for (int i = 0; i < N; ++i) {
        if (i) cout << ' ';
        cout << topo[i];
    }
    cout << '\n';

    return 0;
}

