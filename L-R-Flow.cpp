const int N = 3e5 + 9;
const long long inf = 1LL << 61;
struct Dinic {
  struct edge {
    int to, rev;
    long long flow, w;
    int id; };
  int n, s, t, mxid;
  vector<int> d, flow_through;
  vector<int> done;
  vector<vector<edge>> g;
  Dinic() {}
  Dinic(int _n) {
    n = _n + 10; mxid = 0; g.resize(n);}
  void add_edge(int u, int v, long long w, int id = -1) {