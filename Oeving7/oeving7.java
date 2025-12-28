import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.*;
import java.util.concurrent.TimeUnit;

public class oeving7 {

    // Types and constants
    static final int INF = Integer.MAX_VALUE / 4;

    // ----- Basic data structures -----

    static class Edge {
        int to;
        int w; // hundredths of seconds

        Edge(int to, int w) {
            this.to = to;
            this.w = w;
        }
    }

    static class Node {
        double lat = 0.0;
        double lon = 0.0;
    }

    static class Poi {
        int node;
        int code;
        String name;

        Poi(int node, int code, String name) {
            this.node = node;
            this.code = code;
            this.name = name;
        }
    }

    static class Graph {
        Node[] nodes;           // id -> coords
        List<Edge>[] adj;       // forward graph
        List<Edge>[] radj;      // reverse graph (for ALT preproc)
        List<Poi> pois = new ArrayList<>(); // points of interest
    }

    // I/O helpers

    // split by ASCII whitespace, up to maxFields
    static int hsplit(String line, List<String> out, int maxFields) {
        out.clear();
        int n = line.length();
        int j = 0;
        for (int k = 0; k < maxFields; k++) {
            while (j < n && line.charAt(j) <= ' ') j++;
            if (j >= n) {
                return out.size();
            }
            int start = j;
            while (j < n && line.charAt(j) > ' ') j++;
            out.add(line.substring(start, j));
        }
        return out.size();
    }

    // Loader

    static class DataPaths {
        String dir;
        String nodes = "noder.txt";
        String edges = "kanter.txt";
        String pois = "interessepkt.txt";
    }

    static String join(String a, String b) {
        if (a == null || a.isEmpty()) return b;
        char last = a.charAt(a.length() - 1);
        if (last == '/' || last == '\\') return a + b;
        return a + "/" + b;
    }

    static boolean load_nodes(String file, Graph G) {
        try (BufferedReader in = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
            String line = in.readLine();
            if (line == null) return false;
            line = line.trim();
            if (line.isEmpty()) return false;
            long Nlong = Long.parseLong(line);
            if (Nlong > Integer.MAX_VALUE) {
                System.err.println("Too many nodes for Java int: " + Nlong);
                return false;
            }
            int N = (int) Nlong;

            G.nodes = new Node[N];
            G.adj = new ArrayList[N];
            G.radj = new ArrayList[N];
            for (int i = 0; i < N; i++) {
                G.nodes[i] = new Node();
                G.adj[i] = new ArrayList<>();
                G.radj[i] = new ArrayList<>();
            }

            List<String> f = new ArrayList<>(4);
            for (int i = 0; i < N; i++) {
                line = in.readLine();
                if (line == null) {
                    System.err.println("nodes: premature EOF");
                    return false;
                }
                hsplit(line, f, 3);
                if (f.size() < 3) {
                    System.err.println("nodes: bad line: " + line);
                    return false;
                }
                int id = (int) Long.parseLong(f.get(0));
                double lat = Double.parseDouble(f.get(1));
                double lon = Double.parseDouble(f.get(2));
                if (id < 0 || id >= N) {
                    System.err.println("Node id out of range: " + id + " vs N=" + N);
                    return false;
                }
                G.nodes[id].lat = lat;
                G.nodes[id].lon = lon;
            }
            System.err.println("Loaded nodes: " + G.nodes.length);
            return true;
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error loading nodes: " + e.getMessage());
            return false;
        }
    }

    static boolean load_edges(String file, Graph G) {
        try (BufferedReader in = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
            String line = in.readLine();
            if (line == null) return false;
            line = line.trim();
            if (line.isEmpty()) return false;
            long M = Long.parseLong(line);

            List<String> f = new ArrayList<>(6);
            long cnt = 0;
            while ((line = in.readLine()) != null) {
                if (line.isEmpty()) continue;
                hsplit(line, f, 5);
                if (f.size() < 3) continue;

                int u = (int) Long.parseLong(f.get(0));
                int v = (int) Long.parseLong(f.get(1));
                int t = (int) Long.parseLong(f.get(2)); // hundredths of seconds

                if (u < 0 || v < 0 || u >= G.adj.length || v >= G.adj.length) continue;

                G.adj[u].add(new Edge(v, t));
                G.radj[v].add(new Edge(u, t)); // reverse
                cnt++;
            }
            System.err.println("Loaded edges: " + cnt + " (header M=" + M + ")");
            return true;
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error loading edges: " + e.getMessage());
            return false;
        }
    }

    static boolean load_pois(String file, Graph G) {
        try (BufferedReader in = new BufferedReader(new FileReader(file, StandardCharsets.UTF_8))) {
            String line;
            List<String> f = new ArrayList<>(3);

            while ((line = in.readLine()) != null) {
                if (line.isEmpty()) continue;

                f.clear();
                hsplit(line, f, 3);
                if (f.size() < 2) continue;

                int node = (int) Long.parseLong(f.get(0));
                int code = (int) Long.parseLong(f.get(1));

                String name;
                int qpos = line.indexOf('"');
                if (qpos != -1) {
                    name = line.substring(qpos);
                    if (!name.isEmpty() && name.charAt(0) == '"') {
                        name = name.substring(1);
                    }
                    if (!name.isEmpty() && name.charAt(name.length() - 1) == '"') {
                        name = name.substring(0, name.length() - 1);
                    }
                } else {
                    name = "";
                }

                G.pois.add(new Poi(node, code, name));
            }

            System.err.println("Loaded POIs: " + G.pois.size());
            return true;
        } catch (IOException | NumberFormatException e) {
            System.err.println("Error loading POIs: " + e.getMessage());
            return false;
        }
    }

    // ----- Dijkstra / ALT core -----

    static class SearchStats {
        long pq_pops = 0;
        long relaxations = 0;
        long ms = 0;
    }

    static class ResultRoute {
        List<Integer> path = new ArrayList<>();
        int total_time = INF; // hundredths of seconds
        SearchStats stats = new SearchStats();
    }

    static class FrontierItem {
        int f; // f = g+h, or g for Dijkstra
        int v;

        FrontierItem(int f, int v) {
            this.f = f;
            this.v = v;
        }
    }

    // Heuristic container for ALT
    static class ALTData {
        int[][] fromL;       // L x N   (distances from landmark to nodes)
        int[][] toL;         // N x L   (distances from nodes to landmark)
        List<Integer> landmark_ids = new ArrayList<>(); // node ids used as landmarks
        boolean loaded = false;
    }

    // Heuristic max over landmarks
    static int alt_heuristic(ALTData alt, int n, int goal) {
        if (!alt.loaded) return 0;
        int best = 0;
        int L = alt.landmark_ids.size();
        for (int l = 0; l < L; l++) {
            int a = alt.fromL[l][goal]; // d(L, t)
            int b = alt.fromL[l][n];    // d(L, n)
            if (a < INF && b < INF) {
                best = Math.max(best, Math.max(0, a - b)); // max(0, d(L,t) - d(L,n))
            }

            int c = alt.toL[n][l];      // d(n, L)
            int d = alt.toL[goal][l];   // d(t, L)
            if (c < INF && d < INF) {
                best = Math.max(best, Math.max(0, c - d)); // max(0, d(n,L) - d(t,L))
            }
        }
        return best;
    }

    static ResultRoute dijkstra_search(Graph G, int s, int t) {
        int N = G.nodes.length;
        int[] dist = new int[N];
        Arrays.fill(dist, INF);
        int[] parent = new int[N];
        Arrays.fill(parent, -1);
        boolean[] used = new boolean[N];

        PriorityQueue<FrontierItem> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a.f));
        long t0 = System.nanoTime();

        dist[s] = 0;
        parent[s] = s;
        pq.add(new FrontierItem(0, s));
        SearchStats st = new SearchStats();

        while (!pq.isEmpty()) {
            FrontierItem item = pq.poll();
            int f = item.f;
            int u = item.v;

            if (used[u]) continue;
            used[u] = true;
            st.pq_pops++;

            if (u == t) break;

            for (Edge e : G.adj[u]) {
                if (dist[u] != INF && dist[e.to] > dist[u] + e.w) {
                    dist[e.to] = dist[u] + e.w;
                    parent[e.to] = u;
                    pq.add(new FrontierItem(dist[e.to], e.to));
                    st.relaxations++;
                }
            }
        }

        long t1 = System.nanoTime();
        ResultRoute r = new ResultRoute();
        r.stats = st;
        r.stats.ms = TimeUnit.NANOSECONDS.toMillis(t1 - t0);
        r.total_time = dist[t];

        if (dist[t] < INF) {
            List<Integer> rev = new ArrayList<>();
            for (int v = t; v != s; v = parent[v]) {
                rev.add(v);
                if (parent[v] == -1) break;
            }
            rev.add(s);
            Collections.reverse(rev);
            r.path = rev;
        }
        return r;
    }

    static ResultRoute alt_search(Graph G, ALTData alt, int s, int t) {
        int N = G.nodes.length;
        int[] g = new int[N];
        int[] h = new int[N];
        int[] parent = new int[N];
        boolean[] used = new boolean[N];

        Arrays.fill(g, INF);
        Arrays.fill(h, 0);
        Arrays.fill(parent, -1);

        PriorityQueue<FrontierItem> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a.f));
        long t0 = System.nanoTime();

        g[s] = 0;
        h[s] = alt_heuristic(alt, s, t);
        parent[s] = s;
        pq.add(new FrontierItem(g[s] + h[s], s));
        SearchStats st = new SearchStats();

        while (!pq.isEmpty()) {
            FrontierItem item = pq.poll();
            int f = item.f;
            int u = item.v;

            if (used[u]) continue;
            used[u] = true;
            st.pq_pops++;

            if (u == t) break;

            for (Edge e : G.adj[u]) {
                int ng = (g[u] == INF ? INF : g[u] + e.w);
                if (ng < g[e.to]) {
                    g[e.to] = ng;
                    if (h[e.to] == 0 && e.to != s) {
                        h[e.to] = alt_heuristic(alt, e.to, t);
                    }
                    parent[e.to] = u;
                    pq.add(new FrontierItem(g[e.to] + h[e.to], e.to));
                    st.relaxations++;
                }
            }
        }

        long t1 = System.nanoTime();
        ResultRoute r = new ResultRoute();
        r.stats = st;
        r.stats.ms = TimeUnit.NANOSECONDS.toMillis(t1 - t0);
        r.total_time = g[t];

        if (g[t] < INF) {
            List<Integer> rev = new ArrayList<>();
            for (int v = t; v != s; v = parent[v]) {
                rev.add(v);
                if (parent[v] == -1) break;
            }
            rev.add(s);
            Collections.reverse(rev);
            r.path = rev;
        }
        return r;
    }

    // Landmark selection and preprocessing

    static class Cand {
        double key;
        int id;

        Cand(double key, int id) {
            this.key = key;
            this.id = id;
        }
    }

    static void pushUnique(List<Integer> list, int id) {
        if (!list.contains(id)) {
            list.add(id);
        }
    }

    static List<Integer> choose_landmarks(Graph G, int K) {
        int N = G.nodes.length;
        List<Cand> byLat = new ArrayList<>(N);
        List<Cand> byLon = new ArrayList<>(N);
        for (int i = 0; i < N; i++) {
            byLat.add(new Cand(G.nodes[i].lat, i));
            byLon.add(new Cand(G.nodes[i].lon, i));
        }
        byLat.sort(Comparator.comparingDouble(a -> a.key));
        byLon.sort(Comparator.comparingDouble(a -> a.key));

        List<Integer> L = new ArrayList<>();

        if (!byLat.isEmpty()) {
            pushUnique(L, byLat.get(0).id);
            pushUnique(L, byLat.get(byLat.size() - 1).id);
        }
        if (!byLon.isEmpty()) {
            pushUnique(L, byLon.get(0).id);
            pushUnique(L, byLon.get(byLon.size() - 1).id);
        }

        if (byLat.size() > 4) {
            pushUnique(L, byLat.get(byLat.size() / 4).id);
            pushUnique(L, byLat.get(3 * byLat.size() / 4).id);
        }
        if (byLon.size() > 4) {
            pushUnique(L, byLon.get(byLon.size() / 4).id);
            pushUnique(L, byLon.get(3 * byLon.size() / 4).id);
        }

        while (L.size() < K && N > 0) {
            int id = (int) ((L.size() * 2654435761L) % N);
            pushUnique(L, id);
        }
        if (L.size() > K) {
            L = new ArrayList<>(L.subList(0, K));
        }
        return L;
    }

    // One-to-all Dijkstra (preprocessing)
    static int[] dijkstra_all(List<Edge>[] adj, int s) {
        int N = adj.length;
        int[] dist = new int[N];
        boolean[] used = new boolean[N];
        Arrays.fill(dist, INF);
        Arrays.fill(used, false);

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[0]));
        dist[s] = 0;
        pq.add(new int[]{0, s});

        while (!pq.isEmpty()) {
            int[] item = pq.poll();
            int d = item[0];
            int u = item[1];
            if (used[u]) continue;
            used[u] = true;

            for (Edge e : adj[u]) {
                if (dist[e.to] > d + e.w) {
                    dist[e.to] = d + e.w;
                    pq.add(new int[]{dist[e.to], e.to});
                }
            }
        }
        return dist;
    }

    // Save/load binary matrices
    static boolean save_alt_cache(String base, ALTData alt) {
        String f1 = base + ".fromL.bin";
        String f2 = base + ".toL.bin";
        String f3 = base + ".meta.txt";

        try (DataOutputStream o1 = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f1)));
             DataOutputStream o2 = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(f2)));
             PrintWriter o3 = new PrintWriter(new OutputStreamWriter(new FileOutputStream(f3), StandardCharsets.UTF_8))) {

            long L = alt.landmark_ids.size();
            long N = (alt.fromL == null || alt.fromL.length == 0) ? 0 : alt.fromL[0].length;

            o3.println(N + " " + L);
            for (int id : alt.landmark_ids) {
                o3.println(id);
            }

            // fromL: L x N
            for (int l = 0; l < L; l++) {
                int[] row = alt.fromL[l];
                for (int v = 0; v < N; v++) {
                    o1.writeInt(row[v]);
                }
            }

            // toL: N x L
            for (int v = 0; v < N; v++) {
                int[] col = alt.toL[v];
                for (int l = 0; l < L; l++) {
                    o2.writeInt(col[l]);
                }
            }

            return true;
        } catch (IOException e) {
            System.err.println("Cannot write ALT cache: " + e.getMessage());
            return false;
        }
    }

    static boolean load_alt_cache(String base, ALTData alt) {
        String f1 = base + ".fromL.bin";
        String f2 = base + ".toL.bin";
        String f3 = base + ".meta.txt";

        try (DataInputStream i1 = new DataInputStream(new BufferedInputStream(new FileInputStream(f1)));
             DataInputStream i2 = new DataInputStream(new BufferedInputStream(new FileInputStream(f2)));
             BufferedReader i3 = new BufferedReader(new InputStreamReader(new FileInputStream(f3), StandardCharsets.UTF_8))) {

            String header = i3.readLine();
            if (header == null) return false;
            String[] parts = header.trim().split("\\s+");
            if (parts.length < 2) return false;

            long Nlong = Long.parseLong(parts[0]);
            long Llong = Long.parseLong(parts[1]);
            if (Nlong > Integer.MAX_VALUE || Llong > Integer.MAX_VALUE) return false;

            int N = (int) Nlong;
            int L = (int) Llong;

            List<Integer> ids = new ArrayList<>(L);
            for (int i = 0; i < L; i++) {
                String line = i3.readLine();
                if (line == null) return false;
                ids.add(Integer.parseInt(line.trim()));
            }
            alt.landmark_ids = ids;

            alt.fromL = new int[L][N];
            alt.toL = new int[N][L];

            // fromL: L x N
            for (int l = 0; l < L; l++) {
                for (int v = 0; v < N; v++) {
                    alt.fromL[l][v] = i1.readInt();
                }
            }

            // toL: N x L
            for (int v = 0; v < N; v++) {
                for (int l = 0; l < L; l++) {
                    alt.toL[v][l] = i2.readInt();
                }
            }

            alt.loaded = true;
            return true;
        } catch (IOException | NumberFormatException e) {
            return false;
        }
    }

    static ALTData preprocess_ALT(Graph G, int K, String cache_base) {
        ALTData alt = new ALTData();
        if (load_alt_cache(cache_base, alt)) {
            System.err.println("ALT cache loaded.");
            return alt;
        }
        alt.landmark_ids = choose_landmarks(G, K);
        int N = G.nodes.length;
        int L = alt.landmark_ids.size();
        alt.fromL = new int[L][N];
        alt.toL = new int[N][L];

        for (int l = 0; l < L; l++) {
            Arrays.fill(alt.fromL[l], INF);
        }
        for (int v = 0; v < N; v++) {
            Arrays.fill(alt.toL[v], INF);
        }

        for (int l = 0; l < L; l++) {
            int s = alt.landmark_ids.get(l);
            System.err.println("Preprocess L[" + l + "] node=" + s + " (forward)");
            int[] dF = dijkstra_all(G.adj, s);
            alt.fromL[l] = dF;
        }

        for (int l = 0; l < L; l++) {
            int s = alt.landmark_ids.get(l);
            System.err.println("Preprocess L[" + l + "] node=" + s + " (reverse)");
            int[] dR = dijkstra_all(G.radj, s); // distances on reverse => v->L
            for (int v = 0; v < N; v++) {
                alt.toL[v][l] = dR[v];
            }
        }

        alt.loaded = true;
        save_alt_cache(cache_base, alt);
        return alt;
    }

    // ----- Nearest-k POIs via Dijkstra -----

    static class NearbyPOI {
        int node;
        int code;
        String name;
        int time; // hundredths of seconds

        NearbyPOI(int node, int code, String name, int time) {
            this.node = node;
            this.code = code;
            this.name = name;
            this.time = time;
        }
    }

    static class POIIndexEntry {
        int idx;
        int code;

        POIIndexEntry(int idx, int code) {
            this.idx = idx;
            this.code = code;
        }
    }

    static List<NearbyPOI> nearest_pois(Graph G, int src, int mask, int k) {
        int N = G.nodes.length;
        int[] dist = new int[N];
        boolean[] used = new boolean[N];
        Arrays.fill(dist, INF);
        Arrays.fill(used, false);

        // node -> list of (poiIndex, code)
        List<POIIndexEntry>[] poiIndex = new ArrayList[N];
        for (int i = 0; i < N; i++) poiIndex[i] = new ArrayList<>();

        for (int i = 0; i < G.pois.size(); i++) {
            Poi p = G.pois.get(i);
            if (p.node >= 0 && p.node < N) {
                poiIndex[p.node].add(new POIIndexEntry(i, p.code));
            }
        }

        PriorityQueue<int[]> pq = new PriorityQueue<>(Comparator.comparingInt(a -> a[0]));
        dist[src] = 0;
        pq.add(new int[]{0, src});

        List<NearbyPOI> out = new ArrayList<>(k);

        while (!pq.isEmpty()) {
            int[] item = pq.poll();
            int d = item[0];
            int u = item[1];

            if (used[u]) continue;
            used[u] = true;

            for (POIIndexEntry entry : poiIndex[u]) {
                int code = entry.code;
                if ((code & mask) == mask || (code & mask) != 0) {
                    Poi P = G.pois.get(entry.idx);
                    out.add(new NearbyPOI(u, code, P.name, d));
                    if (out.size() >= k) return out;
                }
            }

            for (Edge e : G.adj[u]) {
                if (dist[e.to] > d + e.w) {
                    dist[e.to] = d + e.w;
                    pq.add(new int[]{dist[e.to], e.to});
                }
            }
        }

        return out;
    }

    // ----- Formatting & CSV -----

    static String fmt_hms(int hundredths) {
        if (hundredths >= INF) return "unreachable";
        long total_ms = (long) hundredths * 10L;
        long s = total_ms / 1000L;
        long ms = total_ms % 1000L;
        long h = s / 3600L;
        s %= 3600L;
        long m = s / 60L;
        s %= 60L;
        return String.format("%d:%02d:%02d.%03d", h, m, s, ms);
    }

    static boolean write_csv_path(String file, Graph G, List<Integer> path) {
        try (PrintWriter o = new PrintWriter(new OutputStreamWriter(new FileOutputStream(file), StandardCharsets.UTF_8))) {
            o.println("lat,lon");
            for (int v : path) {
                String lat = String.format(Locale.US, "%.7f", G.nodes[v].lat);
                String lon = String.format(Locale.US, "%.7f", G.nodes[v].lon);
                o.println(lat + "," + lon);
            }
            return true;
        } catch (IOException e) {
            return false;
        }
    }

    static boolean write_csv_pois(String file, Graph G, List<NearbyPOI> pois) {
        try (PrintWriter o = new PrintWriter(new OutputStreamWriter(new FileOutputStream(file), StandardCharsets.UTF_8))) {
            o.println("node,lat,lon,code,name,time_hms");

            for (NearbyPOI p : pois) {
                String lat = String.format(Locale.US, "%.7f", G.nodes[p.node].lat);
                String lon = String.format(Locale.US, "%.7f", G.nodes[p.node].lon);
                o.println(p.node + "," +
                        lat + "," +
                        lon + "," +
                        p.code + "," +
                        "\"" + p.name + "\"," +
                        fmt_hms(p.time));
            }
            return true;
        } catch (IOException e) {
            return false;
        }
    }

    // ----- CLI Interface -----

    static class Args {
        String data_dir;
        boolean do_pre = false;
        int L = 6;
        String cache_base;
        boolean do_route = false;
        int s = 0, t = 0;
        String algo = "dijkstra";
        String out_csv = "";
        boolean do_nearest = false;
        int src = 0;
        int mask = 0;
        int k = 5;
    }

    static Args parse_args(String[] argv) {
        Args a = new Args();
        int i = 0;

        while (i < argv.length) {
            String v = argv[i];

            switch (v) {
                case "--data":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for --data");
                        System.exit(1);
                    }
                    a.data_dir = argv[i + 1];
                    i += 2;
                    break;

                case "--preprocess":
                    a.do_pre = true;
                    i += 1;
                    break;

                case "--landmarks":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for --landmarks");
                        System.exit(1);
                    }
                    a.L = Integer.parseInt(argv[i + 1]);
                    i += 2;
                    break;

                case "--route":
                    if (i + 2 >= argv.length) {
                        System.err.println("Missing values for --route");
                        System.exit(1);
                    }
                    a.do_route = true;
                    a.s = (int) Long.parseLong(argv[i + 1]);
                    a.t = (int) Long.parseLong(argv[i + 2]);
                    i += 3;
                    break;

                case "--algo":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for --algo");
                        System.exit(1);
                    }
                    a.algo = argv[i + 1];
                    i += 2;
                    break;

                case "--kml":
                case "--csv":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for " + v);
                        System.exit(1);
                    }
                    a.out_csv = argv[i + 1];
                    i += 2;
                    break;

                case "--nearest":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for --nearest");
                        System.exit(1);
                    }
                    a.do_nearest = true;
                    a.src = (int) Long.parseLong(argv[i + 1]);
                    i += 2;
                    break;

                case "--mask":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for --mask");
                        System.exit(1);
                    }
                    a.mask = (int) Long.parseLong(argv[i + 1]);
                    i += 2;
                    break;

                case "--k":
                    if (i + 1 >= argv.length) {
                        System.err.println("Missing value for --k");
                        System.exit(1);
                    }
                    a.k = Integer.parseInt(argv[i + 1]);
                    i += 2;
                    break;

                default:
                    System.err.println("Unknown arg: " + v);
                    i += 1;
                    break;
            }
        }

        if (a.data_dir == null || a.data_dir.isEmpty()) {
            System.err.println("Use --data <folder>");
            System.exit(2);
        }
        // IMPORTANT: Java-only cache name, so it doesn't clash with C++ cache
        a.cache_base = join(a.data_dir, "alt_cache_java_L" + a.L);
        return a;
    }

    public static void main(String[] args) {
        Args a = parse_args(args);
        DataPaths P = new DataPaths();
        P.dir = a.data_dir;
        String nodesPath = join(P.dir, P.nodes);
        String edgesPath = join(P.dir, P.edges);
        String poisPath = join(P.dir, P.pois);

        Graph G = new Graph();
        if (!load_nodes(nodesPath, G)) System.exit(1);
        if (!load_edges(edgesPath, G)) System.exit(2);
        load_pois(poisPath, G);

        ALTData alt = new ALTData();
        if ("alt".equals(a.algo) || a.do_pre) {
            alt = preprocess_ALT(G, a.L, a.cache_base);
        }

        if (a.do_route) {
            long t_search_start = System.nanoTime();

            ResultRoute r = "alt".equals(a.algo)
                    ? alt_search(G, alt, a.s, a.t)
                    : dijkstra_search(G, a.s, a.t);

            long t_search_end = System.nanoTime();
            long search_ms = TimeUnit.NANOSECONDS.toMillis(t_search_end - t_search_start);

            System.out.println("ALGO = " + ("alt".equals(a.algo) ? "ALT" : "Dijkstra"));
            System.out.println("start = " + a.s + " goal = " + a.t);
            System.out.println("travel_time(h:mm:ss.mmm) = " + fmt_hms(r.total_time));
            System.out.println("queue_pops = " + r.stats.pq_pops +
                    " relaxations = " + r.stats.relaxations +
                    " time_ms = " + r.stats.ms);
            if (!r.path.isEmpty()) {
                System.out.println("path_nodes = " + r.path.size());
            }
            if (!a.out_csv.isEmpty() && !r.path.isEmpty()) {
                if (write_csv_path(a.out_csv, G, r.path)) {
                    System.out.println("wrote " + a.out_csv + " (lat, lon) for plotting");
                } else {
                    System.err.println("Failed to write " + a.out_csv);
                }
            }
        }

        if (a.do_nearest) {
            long t_near_start = System.nanoTime();

            List<NearbyPOI> res = nearest_pois(G, a.src, a.mask, a.k);

            long t_near_end = System.nanoTime();
            long near_ms = TimeUnit.NANOSECONDS.toMillis(t_near_end - t_near_start);

            System.out.println("nearest k = " + a.k + " from node = " + a.src +
                    " mask = " + a.mask);
            for (int i = 0; i < res.size(); i++) {
                NearbyPOI p = res.get(i);
                System.out.println((i + 1) + ". node = " + p.node +
                        " time = " + fmt_hms(p.time) +
                        " code = " + p.code +
                        " name = \"" + p.name + "\"" +
                        " lat = " + G.nodes[p.node].lat +
                        " lon = " + G.nodes[p.node].lon);
            }
            if (res.isEmpty()) {
                System.out.println("No POIs found for given mask.");
            }

            System.out.println("nearest_pois_time_ms = " + near_ms);
            if (!a.out_csv.isEmpty()) {
                if (write_csv_pois(a.out_csv, G, res)) {
                    System.out.println("wrote " + a.out_csv);
                } else {
                    System.err.println("Failed to write " + a.out_csv);
                }
            }
        }
    }
}
