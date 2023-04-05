// if we move from b to c, around a
//negative if right,positive if left,zero if collinear
T orient(pt a, pt b, pt c) {return cross(b-a,c-a);}
// angle made at a while we turn from b to c  [0,2*PI]
ld orientedAngle(pt a,pt b,pt c) {
    if (orient(a, b, c) >= 0) return angle(b - a, c - a);
    else return 2 * PI - angle(b - a, c - a);
}
// amplitude travelled at A, from P to Q [-PI,PI]
ld angleTravelled(pt a,pt p,pt q) {
    ld ampli = angle(p - a, q - a);
    if (orient(a, p, q) > 0) return ampli;return -ampli;}
// returns if p lies inside the angle with a at center
bool isInAngle(pt a,pt b,pt c,pt p) {
    assert(orient(a, b, c) != 0);
    if (orient(a, b, c) < 0) swap(b, c);
    return orient(a, b, p) >= 0 && orient(a, c, p) <= 0;
}
//check if a polygon is convex,vertices in order
bool isConvex(vector<pt>&p) {
    bool pos = false, neg = false;
    for (int i = 0, n = sz(p); i < n; i++) {
        int o = orient(p[(i + 1) % n], p[i], p[(i + 2)%n]);
        if (o > 0) pos = true;if (o < 0) neg = true;
    }return !(pos && neg);}
bool half(pt p) {
    assert(p.x != 0 || p.y != 0); // (0,0) is not defined
    return p.y > 0 || (p.y == 0 && p.x < 0);
}
// sort in order of (-PI,PI]
void polarSort(vector<pt>&v) {
    sort(v.begin(), v.end(), [](pt a, pt b) {
        return make_tuple(half(a), 0) < make_tuple(half(b), cross(a, b));
    });}
struct line {
    pt v;T c; //parallel to line  // from ax+by=x
    //so our equation of line is cross(v,xi+yj)=c
    line(pt v, T c) : v(v), c(c) {}
    // from equation ax+by=c
    line(T a, T b, T c) : v({b, -a}), c(c) {}
    line(pt p, pt q) : v(q - p), c(cross(v, p)) {}
    T side(pt p) { return cross(v, p) - c; }
    ld dist(pt p) { return abs(side(p)) / abs(v); }
    ld sqDist(pt p){return(side(p)*side(p))/(double)sq(v);}
    // perpendicular line that passes through p
    line perpThrough(pt p) { return {p, p + prep(v)}; }
    // return in direction of line which point comes first
    bool cmpProj(pt p,pt q){return dot(v, p) < dot(v, q);}
    // return the translated line by vector t
    line translate(pt t) { return {v, c + cross(v, t)};}
    line shiftLeft(ld dist){return {v,T(c+dist * abs(v))};}
    pt proj(pt p) { return p - prep(v) * side(v) / sq(v); }
    pt ref(pt p) { return p - prep(v) * 2 * side(v)/sq(v);}
};
bool inter(line l1,line l2,pt &out){
    T d=cross(l1.v,l2.v);if(d==0) return false;
    out=(l2.v*l1.c-l1.v*l2.c)/d; // requires floating point
    return true;
}
line bisector(line l1, line l2, bool interior) {
    assert(cross(l1.v, l2.v) != 0); // l1 and l2 should not be parallel
    double sign = interior ? 1 : -1;
    return {l2.v / abs(l2.v) + l1.v / abs(l1.v) * sign,
            T(l2.c / abs(l2.v) + l1.c / abs(l1.v) * sign)};}
// p in disk with a,b as diameter
bool inDisk(pt a,pt b,pt p) {return dot(a - p, b - p) <=0;}
// to check if point p lie on segment [a,b]
// orient(a, b, p) == 0 && inDisk(a, b, p);
//if segment [a,b], [c,d] intersection properly(~endpoints)
bool properInter(pt a,pt b,pt c,pt d,pt &out) {
    ld oa = orient(c, d, a), ob = orient(c, d, b);
    ld oc = orient(a, b, c), od = orient(a, b, d);
    if (oa * ob < 0 && oc * od < 0) {
        out = (a * ob - b * oa) / (ob - oa);
        return true;}return false;}
//To create sets of points
struct cmp {
    bool operator()(pt a, pt b) {
        return make_pair(a.x, a.y) < make_pair(b.x, b.y);}};
set<pt,cmp> inter(pt a, pt b,pt c, pt d) {
    pt out;
    if (properInter(a, b, c, d, out)) return {out};
    set<pt, cmp> s;
    if (onSegment(c, d, a)) s.insert(a);
    if (onSegment(c, d, b)) s.insert(b);
    if (onSegment(a, b, c)) s.insert(c);
    if (onSegment(a, b, d)) s.insert(d);return s;}
ld segPointDist(pt a,pt b,pt p) {
    if (a != b) {line l(a, b);
        if (l.cmpProj(a, p) && l.cmpProj(p, b)) return l.dist(p);
    }return min(abs(p - a), abs(p - b));}
ld segSegDist(pt a,pt b,pt c,pt d) {
    pt dummy; if (properInter(a, b, c, d, dummy)) return 0;
    return min({segPointDist(a, b, c),segPointDist(a,b,d),
        segPointDist(c, d, a),segPointDist(c, d, b)});}
// counterclockwise order,if in clockwise order
// area is negative for both concave and convex
int areaPolygonTw(vector<pt>&p) {
    int area = 0;for (int i = 0, n = p.size(); i < n; i++){
        area+=cross(p[i],p[(i + 1) % n]);}returnabs(area);}
// point in order of convex polygon
// if area is greater than zero then points are in counterClockwise
// otherwise clockwise
bool counterClockwise(vector<pt>&p) {
    int area = 0;for (int i = 0, n = p.size(); i < n; i++){
    area += cross(p[i], p[(i + 1) % n]);}return area > 0;}
bool above(pt a,pt p) {return p.y >= a.y;}
//check if [PQ] crosses ray from A(going right)
bool crossesRay(pt a,pt p,pt q) {
    return (above(a, q) - above(a, p)) * orient(a, p, q) > 0;}
//check if a point lie inside polygon(can be convex or concave)
bool inPolygon(vector<pt>&v,pt a,bool strict=true) {
    int numCrossings = 0;
    for (int i = 0, n = v.size(); i < n; i++) {
        if (onSegment(v[i],v[(i+1)%n],a)) {return !strict;}
        numCrossings+= crossesRay(a, v[i],v[(i + 1) % n]);}
    return numCrossings&1;}
// reorder by first element being left bottom point
void reorder(vector<pt>&a){
    int pos=0;for(int i=1;i<a.size();i++)
        if(a[i].x<a[pos].x || (a[i].x==a[pos].x && a[i].y<a[pos].y))pos=i;
    rotate(a.begin(),a.begin()+pos,a.end());}
// lattice points in a polygon strictly inside, s=i+b/2-1
pair<int,int> latticePoints(vector<pt>&p) {
    int area = areaPolygonTw(p),bdryPoints =0,n = p.size();
    for (int i = 0; i < p.size(); i++) {
        bdryPoints += latticePoints(p[i], p[(i + 1) % n]);}
    bdryPoints += sz(p);
    return {(area - bdryPoints + 2) / 2, bdryPoints};
}
bool ok(pt &a,pt &b,pt &c,bool includeCollinear){
    int o = orient(a, b, c);
return o > 0 || (includeCollinear && o == 0);
}
vector<pt> convexHullGC(vector<pt>a,bool includeCollinear=false) {
    if (a.size() <= 2) return a; int n = sz(a);
    for (int i = 1; i < n; i++) {
        if (a[i].y < a[0].y || (a[i].y == a[0].y && a[i].x < a[0].x)) {
            swap(a[0], a[i]);}}pt o = a[0];
    sort(a.begin() + 1, a.end(), [&o](pt a, pt b) {
        int val = orient(o, a, b);if (val == 0) {
            return sq(o - a) < sq(o - b);}return val> 0;});
    vector<pt> stk;stk.push_back(a[0]);stk.push_back(a[1]);
    for (int i = 2; i < n; i++) {
        while (stk.size() >= 2 && !ok(stk[sz(stk) - 2], stk[sz(stk) - 1], a[i], includeCollinear)) {
            stk.pop_back();}stk.push_back(a[i]);}
    if (includeCollinear) {
        for (int j = n - 2; j > 0; j--) {
            if (orient(a[0], a[j], a.back()) == 0) {
                stk.push_back(a[j]);}}}return stk;}
// minimum euclid distance using line sweep
int minEuclidDistance(vector<pt>&v) {
    set<pair<int, int>> s;int dist = 8e18;
    sort(v.begin(), v.end());int j = 0;
    for (auto &c: v) {
        int d = ceil(sqrt(dist));
        while (!s.empty() && c.x - v[j].x >= d) {
            s.erase({v[j].y, v[j].x});j++;}
        auto it1 = s.lower_bound({c.y - d, c.x});
        auto it2 = s.upper_bound({c.y + d, c.x});
        while (it1 != it2) {
            dist = min((T) dist, (T) sq(c - pt{(T) it1->second, (T) it1->first}));
            ++it1;}s.insert({c.y, c.x});}return dist;}
int diameter(vector<pt>a) {
    a = convexHullGC(a);int n = a.size(),ptr_a = 0,ptr_b=1;
    while (cross(a[(ptr_a + 1) % n] - a[ptr_a], a[(ptr_b + 1) % n] - a[ptr_b]) > 0) {
        ptr_b++;ptr_b %= n;}
    int ans = dot(a[ptr_b] - a[ptr_a], a[ptr_b] - a[ptr_a]);
    int begin_a = ptr_a,begin_b = ptr_b;
    do {
        if (cross(a[(ptr_a + 1) % n] - a[ptr_a], a[(ptr_b + 1) % n] - a[ptr_b]) > 0) {
            ptr_b++,ptr_b %= n;} else {ptr_a++;ptr_a %= n;}
        ans = max(ans, dot(a[ptr_b] - a[ptr_a], a[ptr_b] - a[ptr_a]));
    } while (begin_a !=ptr_a||begin_b!=ptr_b);return ans;}