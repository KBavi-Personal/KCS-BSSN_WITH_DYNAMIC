#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <climits>
#include <string>

using namespace std;

/* ============================================================
   Graph data structure
   - Adjacency list
   - Fast edge existence check
   - Degree tracking
   ============================================================ */
struct Graph {
    int n = 0;
    vector<vector<int>> adj;
    vector<unordered_set<int>> has;
    vector<int> deg;

    Graph() {}

    int add_node() {
        adj.emplace_back();
        has.emplace_back();
        deg.push_back(0);
        return n++;
    }

    bool has_edge(int u, int v) const {
        return has[u].count(v);
    }

    void add_edge(int u, int v) {
        if (u == v || has_edge(u, v)) return;
        has[u].insert(v);
        has[v].insert(u);
        adj[u].push_back(v);
        adj[v].push_back(u);
        deg[u]++; deg[v]++;
    }

    long long num_edges() const {
        long long s = 0;
        for (int d : deg) s += d;
        return s / 2;
    }
};

/* ============================================================
   Deterministic hashing utilities
   - Used instead of RNG so expansion is stable
   ============================================================ */
static inline uint64_t mix64(uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

static inline uint64_t H(uint64_t seed, uint64_t a, uint64_t b, uint64_t c) {
    return mix64(seed ^ mix64(a) ^ mix64(b) ^ mix64(c));
}

static inline int rnd(uint64_t seed, uint64_t a, uint64_t b, uint64_t c, int lo, int hi) {
    return lo + (int)(H(seed, a, b, c) % (uint64_t)(hi - lo + 1));
}

static inline double rnd01(uint64_t seed, uint64_t a, uint64_t b, uint64_t c) {
    return (double)H(seed, a, b, c) / (double)ULLONG_MAX;
}

/* ============================================================
   Add clique edges (used for k-truss backbone)
   ============================================================ */
static void add_clique(Graph& G, const vector<int>& nodes) {
    for (int i = 0; i < (int)nodes.size(); i++)
        for (int j = i + 1; j < (int)nodes.size(); j++)
            G.add_edge(nodes[i], nodes[j]);
}

/* ============================================================
   Community backbone
   - Chain of overlapping k-cliques
   - Controls diameter
   ============================================================ */
struct Community {
    int d;                  // target diameter
    vector<int> backbone;   // nodes of backbone (size k + d)
};

/* ============================================================
   Expandable realistic social network generator
   ============================================================ */
class SocialNetwork {
public:
    Graph G;
    uint64_t seed;
    int k;                         // GLOBAL k-truss parameter
    vector<Community> comms;
    int min_degree;
    int max_degree;

    /* ---- realism parameters ---- */
    double p_core = 0.4;//core: attach into truss cores (triangle-rich)
    double p_triadic = 0.4;//triadic: friend-of-friend edges
    double p_pref = 0.15;//pref: connect to higher-degree nodes (hubs)
    double p_random = 0.05;//random: weak ties

    SocialNetwork(uint64_t seed_, int k_, int min_degree_, int max_degree_) : seed(seed_), k(k_), min_degree(min_degree_), max_degree(max_degree_) {}

    /* --------------------------------------------------------
       Initialize communities
       - random diameters
       - overlapping k-clique backbones
       -------------------------------------------------------- */
    void init_communities(int num_communities, int d_min, int d_max) {
        for (int c = 0; c < num_communities; c++) {
            Community C;
            C.d = rnd(seed, c, 11, 22, d_min, d_max);

            int len = k + C.d;
            for (int i = 0; i < len; i++) C.backbone.push_back(G.add_node());

            // sliding k-cliques → k-truss safe
            for (int i = 0; i <= C.d; i++) {
                vector<int> win;
                for (int j = 0; j < k; j++) win.push_back(C.backbone[i + j]);
                add_clique(G, win);
            }
            comms.push_back(C);
        }
    }

    /* --------------------------------------------------------
       Expand graph to N nodes
       -------------------------------------------------------- */
    void expand_to(int N) {
        while (G.n < N) {
            int u = G.add_node();
            add_user(u);
        }
    }

    /* --------------------------------------------------------
       Write edge list
       -------------------------------------------------------- */
    void save(const string& file) {
        ofstream out(file);
        for (int u = 0; u < G.n; u++)
            for (int v : G.adj[u])
                if (u < v) out << u << " " << v << "\n";
    }

private:
    /* --------------------------------------------------------
       Add edges for one new user
       -------------------------------------------------------- */
    void add_user(int u) {
        int target_deg = rnd(seed, u, 1, 2, min_degree, max_degree);

        int e_core = target_deg * p_core;
        int e_tri = target_deg * p_triadic;
        int e_pref = target_deg * p_pref;
        int e_rand = target_deg - (e_core + e_tri + e_pref);

        /* ---- Core k-truss attachments ---- */
        for (int i = 0; i < e_core; i++) {
            int c = rnd(seed, u, i, 10, 0, (int)comms.size() - 1);
            attach_core(u, c, i);
        }

        /* ---- Triadic closure ---- */
        for (int i = 0; i < e_tri; i++) {
            int w = triadic(u, i);
            if (w != -1) G.add_edge(u, w);
        }

        /* ---- Preferential attachment ---- */
        for (int i = 0; i < e_pref; i++) {
            int v = preferential(u, i);
            if (v != -1) G.add_edge(u, v);
        }

        /* ---- Weak random ties ---- */
        for (int i = 0; i < e_rand; i++) {
            int v = rnd(seed, u, i, 99, 0, G.n - 2);
            if (v != u) G.add_edge(u, v);
        }
    }

    /* --------------------------------------------------------
       k-truss safe attachment
       -------------------------------------------------------- */
    void attach_core(int u, int c, int tag) {
        auto& C = comms[c];
        int win = rnd(seed, u, c, tag, 0, C.d);
        int drop = rnd(seed, u, c, tag + 100, 0, k - 1);

        for (int j = 0; j < k; j++) {
            if (j == drop) continue;
            G.add_edge(u, C.backbone[win + j]);
        }
    }

    /* --------------------------------------------------------
       Triadic closure
       -------------------------------------------------------- */
    int triadic(int u, int tag) {
        if (G.adj[u].empty()) return -1;
        int v = G.adj[u][rnd(seed, u, tag, 30, 0, G.adj[u].size() - 1)];
        if (G.adj[v].empty()) return -1;
        int w = G.adj[v][rnd(seed, u, v, tag, 0, G.adj[v].size() - 1)];
        if (w != u && !G.has_edge(u, w)) return w;
        return -1;
    }

    /* --------------------------------------------------------
       Preferential attachment (degree-biased)
       -------------------------------------------------------- */
    int preferential(int u, int tag) {
        int best = -1;
        for (int t = 0; t < 30; t++) {
            int v = rnd(seed, u, tag, t, 0, G.n - 2);
            if (v == u || G.has_edge(u, v)) continue;
            if (best == -1 || G.deg[v] > G.deg[best]) best = v;
        }
        return best;
    }
};


//int run() {
//    uint64_t seed = 42;
//    int k = 6;
//
//    SocialNetwork net(seed, k);
//
//    // random community diameters in [2,12]
//    net.init_communities(6, 2, 12);
//
//    net.expand_to(10000);
//    cout << "10k nodes, edges=" << net.G.num_edges() << "\n";
//    net.save("social_10k.edgelist");
//
//    net.expand_to(20000);
//    cout << "20k nodes, edges=" << net.G.num_edges() << "\n";
//    net.save("social_20k.edgelist");
//
//    net.expand_to(30000);
//    cout << "30k nodes, edges=" << net.G.num_edges() << "\n";
//    net.save("social_30k.edgelist");
//}
