#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>
#include <cstdio>
#include <locale>
#include <clocale>
using namespace std;

// Types and constants
const int INF = INT32_MAX / 4;

struct Edge { 
    unsigned int to; 
    int w; 
};
struct Node { 
    double lat = 0.0, 
    lon = 0.0; 
};

struct Poi { 
    uint32_t node; 
    uint32_t code; 
    string name; 
};

struct Graph {
    vector<Node> nodes; // id -> coords
    vector<vector<Edge>> adj; // forward graph
    vector<vector<Edge>> radj; // reverse graph (for ALT preproc)
    vector<Poi> pois; // points of interest
};

// I/O helpers
static inline void trim_crlf(string &s){ 
    while(!s.empty() && (s.back() == '\r' || s.back() == '\n')) s.pop_back(); 
}

// Split ved ASCII
static inline int hsplit(const string &line, vector<string> &out, int maxFields){
    out.clear(); 
    out.reserve(maxFields);
    int n = (int)line.size(); 
    int j = 0;
    for(int k = 0; k < maxFields; k++) {
        while(j < n && line[j] <= ' ') j++; 
        if(j >= n) { 
            return (int)out.size(); 
        }
        int start = j; 
        while(j < n && line[j] > ' ') j++; 
        out.emplace_back(line.substr(start, j-start));
    }
    return (int)out.size();
}

// Loader
struct DataPaths { 
    string dir; 
    string nodes = "noder.txt"; 
    string edges = "kanter.txt"; 
    string pois = "interessepkt.txt"; 
};

static inline string join(const string &a, const string &b) { 
    if(a.empty()) return b; 
    if(a.back() == '/' || a.back() == '\\') return a+b; 
    return a + "/" + b; 
}

bool load_nodes(const string &file, Graph &G){
    ifstream in(file); 
    if(!in) { 
        cerr << "Could not open nodes: " << file << "\n"; 
        return false; 
    }
    string line; 
    getline(in,line); 
    trim_crlf(line); 
    if(line.empty()) return false; 
    unsigned long long N = stoull(line);
    G.nodes.assign(N, {}); 
    G.adj.assign(N, {}); 
    G.radj.assign(N, {});

    vector<string> f; 
    f.reserve(4);
    for(unsigned long long i = 0; i < N; i++){
        if(!getline(in,line)){ 
            cerr << "nodes: premature EOF" << "\n"; 
            return false; 
        }
        hsplit(line, f, 3);
        uint32_t id = (uint32_t)stoul(f[0]);
        double lat = stod(f[1]); double lon = stod(f[2]);
        if(id >= N){ 
            cerr << "Node id out of range: " << id << " vs N=" << N << "\n"; 
            return false; 
        }
        G.nodes[id].lat = lat; G.nodes[id].lon = lon;
    }
    cerr << "Loaded nodes: " << G.nodes.size() << "\n";
    return true;
}

bool load_edges(const string &file, Graph &G) {
    ifstream in(file); 
    if(!in) { 
        cerr << "Could not open edges: " << file << "\n";
        return false; 
    }
    string line; 
    getline(in,line); 
    trim_crlf(line); 
    if(line.empty()) return false; 
    unsigned long long M = stoull(line);
    vector<string> f; 
    f.reserve(6);
    unsigned long long cnt = 0; 
    while(getline(in,line)) {
        if(line.empty()) continue; 
        hsplit(line, f, 5);
        uint32_t u = (uint32_t)stoul(f[0]); 
        uint32_t v = (uint32_t)stoul(f[1]);
        int32_t t = (int32_t)stol(f[2]); // hundredths of seconds
        if(u >= G.adj.size() || v >= G.adj.size()) continue;
        G.adj[u].push_back({v,t});
        G.radj[v].push_back({u,t}); // reverse
        cnt++;
    }
    cerr << "Loaded edges: " << cnt << "\n";
    return true;
}

bool load_pois(const string &file, Graph &G){
    ifstream in(file);
    if(!in) {
        cerr << "Could not open POIs: " << file << "\n";
        return false;
    }

    string line;
    vector<string> f;
    f.reserve(3);

    while(getline(in, line)) {
        if(line.empty()) continue;
        trim_crlf(line);

        f.clear();
        hsplit(line, f, 3);
        if(f.size() < 2) continue;

        uint32_t node = (uint32_t)stoul(f[0]);
        uint32_t code = (uint32_t)stoul(f[1]);

        string name;
        size_t qpos = line.find('"');
        if (qpos != string::npos) {
            name = line.substr(qpos); // from first quote to end
            if(!name.empty() && name.front() == '"') name.erase(0, 1);
            if(!name.empty() && name.back()  == '"') name.pop_back();
        } else {
            name = ""; // no name found
        }

        G.pois.push_back({node, code, name});
    }

    cerr << "Loaded POIs: " << G.pois.size() << "\n";
    return true;
}

// Dijkstra / ALT core
struct SearchStats { 
    unsigned long long pq_pops = 0; 
    unsigned long long relaxations = 0; 
    long long ms = 0; 
};

struct ResultRoute {
    vector<uint32_t> path;
    int32_t total_time = INF; // hundredths of seconds
    SearchStats stats;
};

struct FrontierItem { 
    int32_t f; 
    uint32_t v; 
}; // for ALT (f = g+h). For Dijkstra, f==g.

struct PQCmp { 
    bool operator()(const FrontierItem &a, const FrontierItem &b) const { return a.f > b.f; } 
};

// Heuristic container for ALT
struct ALTData {
    vector<vector<int32_t>> fromL, toL; // sizes: L x N   and   N x L  (vector-of-vectors)
    vector<uint32_t> landmark_ids; // node ids used as landmarks
    bool loaded=false;
};

// Heuristic max over landmarks
static inline int32_t alt_heuristic(const ALTData &alt, uint32_t n, uint32_t goal){
    if(!alt.loaded) return 0; 
    int32_t best = 0;
    const size_t L = alt.landmark_ids.size();
    for(size_t l = 0; l < L; l++){
        int32_t a = alt.fromL[l][goal]; // d(L, t)
        int32_t b = alt.fromL[l][n]; // d(L, n)
        if(a < INF && b < INF) best = max(best, max(0, a - b)); // max(0, d(L,t) - d(L,n))
        int32_t c = alt.toL[n][l]; // d(n, L)
        int32_t d = alt.toL[goal][l]; // d(t, L)
        if(c < INF && d < INF) best = max(best, max(0, c - d)); // max(0, d(n,L) - d(t,L))
    }
    return best;
}

ResultRoute dijkstra_search(const Graph &G, uint32_t s, uint32_t t){
    const uint32_t N = (uint32_t)G.nodes.size();
    vector<int32_t> dist(N, INF); 
    vector<uint32_t> parent(N, UINT32_MAX); 
    vector<char> used(N, 0);
    priority_queue<FrontierItem, vector<FrontierItem>, PQCmp> pq;
    auto t0 = chrono::high_resolution_clock::now();
    dist[s]=0; 
    parent[s]=s; 
    pq.push({0, s});
    SearchStats st; 
    while(!pq.empty()){
        auto [f,u] = pq.top(); 
        pq.pop();
        if(used[u]) continue; 
        used[u] = 1; 
        st.pq_pops++;
        if(u == t) break;
        for(const auto &e: G.adj[u]){
            if(dist[u] != INF && dist[e.to] > dist[u] + e.w) {
                dist[e.to] = dist[u] + e.w; 
                parent[e.to] = u; 
                pq.push({dist[e.to], e.to}); 
                st.relaxations++;
            }
        }
    }
    auto t1 = chrono::high_resolution_clock::now();
    ResultRoute r; 
    r.stats = st; 
    r.stats.ms = chrono::duration_cast<chrono::milliseconds>(t1-t0).count(); 
    r.total_time = dist[t];
    if(dist[t]<INF){
        // rebuild
        vector<uint32_t> rev; 
        for(uint32_t v = t; v != s; v = parent[v]) { 
            rev.push_back(v); 
            if(parent[v] == UINT32_MAX) break; 
        }
        rev.push_back(s); 
        reverse(rev.begin(), rev.end()); 
        r.path = move(rev);
    }
    return r;
}

ResultRoute alt_search(const Graph &G, const ALTData &alt, uint32_t s, uint32_t t){
    const uint32_t N = (uint32_t)G.nodes.size();
    vector<int32_t> g(N, INF); 
    vector<int32_t> h(N, 0); 
    vector<uint32_t> parent(N, UINT32_MAX); 
    vector<char> used(N, 0);
    priority_queue<FrontierItem, vector<FrontierItem>, PQCmp> pq;
    auto t0 = chrono::high_resolution_clock::now();
    g[s] = 0; 
    h[s] = alt_heuristic(alt, s, t); 
    parent[s]=s; 
    pq.push({g[s]+h[s], s});
    SearchStats st; 
    while(!pq.empty()){
        auto [f,u] = pq.top(); 
        pq.pop();
        if(used[u]) continue; 
        used[u]=1; 
        st.pq_pops++;
        if(u==t) break;
        for(const auto &e: G.adj[u]) {
            int32_t ng = (g[u] == INF ? INF : g[u] + e.w);
            if(ng < g[e.to]) {
                g[e.to] = ng;
                if(h[e.to] == 0 && e.to!=s) { 
                    h[e.to] = alt_heuristic(alt, e.to, t); 
                }
                parent[e.to] = u;
                pq.push({g[e.to] + h[e.to], e.to}); 
                st.relaxations++;
            }
        }
    }
    auto t1 = chrono::high_resolution_clock::now();
    ResultRoute r; 
    r.stats = st; 
    r.stats.ms = chrono::duration_cast<chrono::milliseconds>(t1-t0).count();
    r.total_time = g[t];
    if(g[t]<INF) { 
        vector<uint32_t> rev; 
        for(uint32_t v = t; v != s; v = parent[v]) {
            rev.push_back(v); 
            if(parent[v] == UINT32_MAX) break; 
        } 
        rev.push_back(s); 
        reverse(rev.begin(), rev.end()); 
        r.path = move(rev); 
    }
    return r;
}

// Landmark selection and preprocessing
// Pick K landmarks based on extremes of lat/lon (min/max, corners)
vector<uint32_t> choose_landmarks(const Graph &G, int K){
    struct Cand { 
        double key; 
        uint32_t id;
    };

    vector<Cand> byLat, byLon; 
    byLat.reserve(G.nodes.size()); 
    byLon.reserve(G.nodes.size());
    for(uint32_t i = 0; i < G.nodes.size(); i++) { 
        byLat.push_back({G.nodes[i].lat, i}); 
        byLon.push_back({G.nodes[i].lon, i}); 
    }
    sort(byLat.begin(), byLat.end(), [](auto&a, auto&b){return a.key<b.key;});
    sort(byLon.begin(), byLon.end(), [](auto&a, auto&b){return a.key<b.key;});
    vector<uint32_t> L;
    auto push_unique=[&](uint32_t id){ 
        if(find(L.begin(), L.end(), id)==L.end()) L.push_back(id); 
    };
    // min/max lat/lon and rough corners
    push_unique(byLat.front().id); 
    push_unique(byLat.back().id);
    push_unique(byLon.front().id); 
    push_unique(byLon.back().id);

    if(byLat.size()>4) { 
        push_unique(byLat[byLat.size()/4].id); 
        push_unique(byLat[3*byLat.size()/4].id); 
    }
    if(byLon.size()>4) { 
        push_unique(byLon[byLon.size()/4].id); 
        push_unique(byLon[3*byLon.size()/4].id); 
    }
    while((int)L.size() < K) { 
        push_unique((uint32_t)(L.size()*2654435761u % G.nodes.size())); 
    }
    if((int)L.size() > K) L.resize(K);
    return L;
}

// One-to-all Dijkstra (preprocessing)
static vector<int32_t> dijkstra_all(const vector<vector<Edge>> &adj, uint32_t s){
    const uint32_t N = (uint32_t)adj.size(); 
    vector<int32_t> dist(N, INF); 
    vector<char> used(N,0);
    priority_queue<
        pair<int32_t,uint32_t>, vector<pair<int32_t,uint32_t>>,
        greater<pair<int32_t,uint32_t>>
    > pq;
    dist[s] = 0; pq.push({0,s});
    while(!pq.empty()){
        auto [d,u] = pq.top(); 
        pq.pop(); 
        if(used[u]) continue; 
        used[u]=1; 
        for(const auto &e: adj[u]){
            if(dist[e.to] > d + e.w) { 
                dist[e.to] = d+e.w; 
                pq.push({dist[e.to], e.to}); 
            }
        }
    }
    return dist;
}

// Save/load binary matrices
bool save_alt_cache(const string &base, const ALTData &alt){
    string f1 = base + ".fromL.bin"; 
    string f2 = base + ".toL.bin"; 
    string f3 = base + ".meta.txt";
    ofstream o1(f1, ios::binary), o2(f2, ios::binary), o3(f3);
    if(!o1 || !o2 || !o3) { 
        cerr << "Cannot write ALT cache" << "\n"; 
        return false; 
    }
    unsigned long long L = alt.landmark_ids.size(); 
    unsigned long long N = alt.fromL.empty() ? 0 : alt.fromL[0].size();

    o3 << N << " " << L << "\n"; 
    for(auto id: alt.landmark_ids) o3 << id << "\n";
    for(size_t l = 0; l < L; l++) { 
        o1.write((const char*)alt.fromL[l].data(), sizeof(int32_t)*N); 
    }
    for(size_t v = 0; v < N; v++) { 
        o2.write((const char*)alt.toL[v].data(), sizeof(int32_t)*L); 
    }
    return true;
}

bool load_alt_cache(const string &base, ALTData &alt) {
    string f1 = base + ".fromL.bin"; string f2 = base + ".toL.bin"; string f3 = base + ".meta.txt";
    ifstream i1(f1, ios::binary), i2(f2, ios::binary), i3(f3);
    if(!i1 || !i2 || !i3) return false;
    unsigned long long N,L; 
    i3>>N>>L; 
    vector<uint32_t> ids(L); 
    for(size_t i = 0; i < L; i++) i3 >> ids[i];
    alt.landmark_ids = move(ids); 
    alt.fromL.assign(L, vector<int32_t>(N));
    alt.toL.assign(N, vector<int32_t>(L));
    for(size_t l = 0; l < L; l++) { 
        i1.read((char*)alt.fromL[l].data(), sizeof(int32_t)*N); 
    }
    for(size_t v = 0; v < N; v++) { 
        i2.read((char*)alt.toL[v].data(), sizeof(int32_t)*L); 
    }
    alt.loaded = true; 
    return true;
}

ALTData preprocess_ALT(const Graph &G, int K, const string &cache_base){
    ALTData alt; 
    if(load_alt_cache(cache_base, alt)) { 
        cerr << "ALT cache loaded." << "\n"; 
        return alt; 
    }
    alt.landmark_ids = choose_landmarks(G, K);
    const uint32_t N = (uint32_t)G.nodes.size(); 
    const size_t L = alt.landmark_ids.size();
    alt.fromL.assign(L, vector<int32_t>(N, INF)); 
    alt.toL.assign(N, vector<int32_t>(L, INF));

    for(size_t l = 0; l < L; l++) {
        uint32_t s = alt.landmark_ids[l]; 
        cerr << "Preprocess L[" << l << "] node=" << s << " (forward)" << "\n";
        auto dF = dijkstra_all(G.adj, s); 
        alt.fromL[l] = move(dF);
    }
    for(size_t l = 0; l < L; l++) {
        uint32_t s = alt.landmark_ids[l]; 
        cerr << "Preprocess L[" << l << "] node=" << s << " (reverse)" << "\n";
        auto dR = dijkstra_all(G.radj, s); // distances on reverse => v->L
        for(uint32_t v = 0; v < N; v++) alt.toL[v][l] = dR[v];
    }
    alt.loaded = true; 
    save_alt_cache(cache_base, alt);
    return alt;
}

// Nearest-k POIs via Dijkstra
struct NearbyPOI { 
    uint32_t node; 
    uint32_t code; 
    string name; 
    int32_t time; 
};

vector<NearbyPOI> nearest_pois(const Graph &G, uint32_t src, uint32_t mask, int k){
    const uint32_t N = (uint32_t)G.nodes.size();
    vector<int32_t> dist(N, INF); 
    vector<char> used(N,0);
    vector<vector<pair<uint32_t,uint32_t>>> poiIndex(N); // node -> list of (code,name_idx)
    for(size_t i = 0; i < G.pois.size(); i++) { 
        const auto &p = G.pois[i]; 
        if(p.node < N) poiIndex[p.node].push_back({(uint32_t)i, p.code}); 
    }

    priority_queue<pair<int32_t,uint32_t>, vector<pair<int32_t,uint32_t>>, greater<pair<int32_t,uint32_t>>> pq;
    dist[src] = 0; 
    pq.push({0,src});
    vector<NearbyPOI> out; out.reserve(k);
    while(!pq.empty()) {
        auto [d,u] = pq.top(); 
        pq.pop(); 
        if(used[u]) continue; 
        used[u] = 1;
        for(auto [idx, code] : poiIndex[u]) {
            if((code & mask) == mask || (code & mask)) {
                const auto &P = G.pois[idx];
                out.push_back({u, code, P.name, d});
                if((int)out.size() >= k) return out;
            }
        }
        for(const auto &e: G.adj[u]) {
            if(dist[e.to] > d + e.w) { 
                dist[e.to] = d + e.w; 
                pq.push({dist[e.to], e.to}); 
            } 
        }
    }
    return out;
}

static string fmt_hms(int32_t hundredths){
    if(hundredths >= INF) return string("unreachable");
    long long total_ms = (long long)hundredths * 10; // millisekund
    long long s = total_ms/1000; 
    long long ms = total_ms%1000; 
    long long h=s/3600; s%=3600; 
    long long m=s/60; s%=60;
    char buf[64]; 
    snprintf(buf, sizeof(buf), "%lld:%02lld:%02lld.%03lld", (long long)h, (long long)m, (long long)s, (long long)ms);
    return string(buf);
}

bool write_csv_path(const string &file, const Graph &G, const vector<uint32_t> &path){
    ofstream o(file); 
    if(!o) return false; 
    o << "lat,lon\n";
    for(uint32_t v : path) { 
        o << fixed<<setprecision(7) << G.nodes[v].lat << "," << G.nodes[v].lon << "\n"; 
    } 
    return true;
}

bool write_csv_pois(const string &file,
                    const Graph &G,
                    const vector<NearbyPOI> &pois)
{
    ofstream o(file);
    if (!o) return false;

    o << "node,lat,lon,code,name,time_hms\n";

    for (const auto &p : pois) {
        o << p.node << ","
          << fixed << setprecision(7) << G.nodes[p.node].lat << ","
          << fixed << setprecision(7) << G.nodes[p.node].lon << ","
          << p.code << ","
          << "\"" << p.name << "\"" << ","
          << fmt_hms(p.time)
          << "\n";
    }
    return true;
}

// CLInterface
struct Args{
    string data_dir; 
    bool do_pre = false; 
    int L = 6; 
    string cache_base; 
    bool do_route = false; 
    uint32_t s = 0, t = 0; 
    string algo = "dijkstra"; 
    string out_csv = "";
    bool do_nearest = false; 
    uint32_t src = 0; uint32_t mask = 0; int k = 5;
};

Args parse_args(int argc, char**argv) {
    Args a; 
    for(int i = 1; i < argc; i++) { 
        string v = argv[i];
        auto need = [&](string flag) { 
            if(i + 1 >= argc) { 
                cerr << "Missing value for " << flag << "\n"; 
                exit(1);
            } 
            return string(argv[++i]); 
        };
        if(v == "--data") a.data_dir = need(v);
        else if(v == "--preprocess") a.do_pre = true;
        else if(v == "--landmarks") a.L = stoi(need(v));
        else if(v == "--route") { 
            a.do_route = true; 
            a.s=(uint32_t)stoul(need(v)); 
            a.t=(uint32_t)stoul(need(v)); 
        }
        else if(v == "--algo") a.algo = need(v);
        else if(v == "--kml" || v == "--csv") a.out_csv = need(v);
        else if(v == "--nearest") { 
            a.do_nearest = true; 
            a.src=(uint32_t)stoul(need(v)); 
        }
        else if(v == "--mask") a.mask = (uint32_t)stoul(need(v));
        else if(v == "--k") a.k = stoi(need(v));
        else { 
            cerr << "Unknown arg: " << v << "\n"; 
        }
    }
    if(a.data_dir.empty()) { 
        cerr << "Use --data <folder>" << "\n"; 
        exit(2);
    } 
    a.cache_base = join(a.data_dir, string("alt_cache_L")+to_string(a.L));
    return a;
}

int main(int argc, char**argv) {
    setlocale(LC_ALL, "");
    ios::sync_with_stdio(false); 
    cin.tie(nullptr);
    Args a = parse_args(argc, argv);
    DataPaths P; 
    P.dir = a.data_dir; 
    string nodes = join(P.dir, P.nodes), edges = join(P.dir, P.edges), pois = join(P.dir, P.pois);

    Graph G; 
    if(!load_nodes(nodes, G)) return 1; 
    if(!load_edges(edges, G)) return 2; 
    load_pois(pois, G);

    ALTData alt; 
    if(a.algo == "alt" || a.do_pre) { alt = preprocess_ALT(G, a.L, a.cache_base); 
    }

    if (a.do_route) {
        auto t_search_start = chrono::high_resolution_clock::now();

        ResultRoute r = (a.algo == "alt"
                        ? alt_search(G, alt, a.s, a.t)
                        : dijkstra_search(G, a.s, a.t));

        auto t_search_end = chrono::high_resolution_clock::now();
        auto search_ms = chrono::duration_cast<chrono::milliseconds>(
                            t_search_end - t_search_start
                        ).count();

        cout << "ALGO = " << ((a.algo == "alt") ? "ALT" : "Dijkstra");
        cout << "\nstart = " << a.s << " goal = " << a.t;
        cout << "\ntravel_time(h:mm:ss.mmm) = " << fmt_hms(r.total_time);
        cout << "\nqueue_pops = " << r.stats.pq_pops
            << " relaxations = " << r.stats.relaxations
            << " time_ms = " << r.stats.ms; 
        cout << "\n";

        if(!r.path.empty()) cout << "path_nodes = " << r.path.size() << "\n";
        if(!a.out_csv.empty() && !r.path.empty()) {
            if(write_csv_path(a.out_csv, G, r.path))
                cout << "wrote " << a.out_csv << " (lat, lon) for plotting\n";
            else
                cerr << "Failed to write " << a.out_csv << "\n";
        }
    }

    if (a.do_nearest) {
        auto t_near_start = chrono::high_resolution_clock::now();

        auto res = nearest_pois(G, a.src, a.mask, a.k);

        auto t_near_end = chrono::high_resolution_clock::now();
        auto near_ms = chrono::duration_cast<chrono::milliseconds>(
                        t_near_end - t_near_start
                    ).count();

        cout << "nearest k = " << a.k << " from node = " << a.src
            << " mask = " << a.mask << "\n";
        for(size_t i = 0; i < res.size(); i++) {
            auto &p = res[i];
            cout << i+1 << ". node = " << p.node
                << " time = " << fmt_hms(p.time)
                << " code = " << p.code
                << " name = \"" << p.name << "\""
                << " lat = " << G.nodes[p.node].lat
                << " lon = " << G.nodes[p.node].lon << "\n";
        }
        if(res.empty()) cout << "No POIs found for given mask.\n";

        cout << "nearest_pois_time_ms = " << near_ms << "\n";
        if (!a.out_csv.empty()) {
            if (write_csv_pois(a.out_csv, G, res)) {
                cout << "wrote " << a.out_csv << "\n";
            } else {
                cerr << "Failed to write " << a.out_csv << "\n";
            }
        }
    }
    return 0;
}
