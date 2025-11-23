#include <iostream>
#include <sstream>
#include <fstream>
#include <queue>
#include <stack>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>
#include <utility>
#include <string>
#include <initializer_list>

#include <assert.h>
#include <errno.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <random>
#include <chrono>
#include <sys/stat.h>
#include "bucket.h"

using namespace std;

#define PHOTON_THRESHOLD 0.6 // threshold for PHOTON 

#define ELECTRON_THRESHOLD 0.1 // threshold for ELECTRON

#define RESIDUALS 1 // 1 if the difference between parent and its child subgraph is reported. 0 otherwise.

#define DEPTH 2 // only the subgraphs with <= DEPTH are reported. Depth of a leaf is 0, depth of a leaf's parent is 1 and so on.
#define LOWERBOUND 10 // lower bound for the size of subgraphs that are shown in *Hierarchy file (and command line output)
#define UPPERBOUND 1000 // compute densities of subgraphs with at most this size, set to INT_MAX to compute all -- takes a lot of time

typedef long long lol;
typedef int vertex; // vertices are 32 bytes
typedef int edge; // edges are 32 bytes
typedef pair<vertex, vertex> vp;
typedef tuple<vertex, vertex, vertex> vt;
typedef vector<vector<vertex> > Graph;

typedef chrono::duration<double> tms;
typedef tuple<vertex, vertex> couple;

class subcore {
	public:
		int levelFromLeaf; // 0 is leaf, 1 is leaf's parent
		bool visible;
		vertex rank;
		vertex K;
		vertex parent;
		vertex root;
		vector<vertex> children;
		vertex size;
		edge nEdge;
		double ed;
		double balance;
		double metrics[4];
		int ne; // number of negative edges
		vertex vl; // size of left set
		vertex vr; // size of right set

		subcore (vertex k) {
			K = k;
			rank = 0;
			parent = -1;
			root = -1;
			visible = true;
			size = 0;
			nEdge = 0;
			ed = -1;
			balance = -1;
			metrics[0] = metrics[1] = metrics[2] = metrics[3] = -1;
			ne = -1;
			vl = vr = -1;
			levelFromLeaf = -1;
		}

		void print (int a, int index = -1, FILE* fp = NULL) {
			string s = "%d\t%d\t%d\t%d\t%.2lf\t%.2lf\t%d\t%.2lf\t%d\t%d\t%.2lf\t%.2lf\t%.2lf\n";
			string t = "%d  %d  %d  %d  %d  %.2lf  %.2lf  %d  %d  %.2lf  %d  %d  %.2lf  %.2lf  %.2lf\t";
			if (fp == NULL) // for command line output
				printf (s.c_str(),
						K, size, nEdge, ne, ed, balance, levelFromLeaf,
						metrics[0],	vl,	vr, 	metrics[1],	metrics[2],	metrics[3]);
			else
				if (a == 1) // for Hierarchy file
					fprintf (fp, "id: %d  "
							"K: %d  "
							"|V|: %d  "
							"|E|: %d  "
							"|E-|: %d  "
							"ed: %.2lf  "
							"balance: %.2lf  "
							"LEAF?: %d  "
							"parent id: %d  "
							"polarity: %.2lf  "
							"|V_l|: %d  "
							"|V_r|: %d  "
							"hLeft: %.2lf  "
							"hRight: %.2lf  "
							"tension: %.2lf\n",
							index, K, size,	nEdge, ne, ed, balance,	levelFromLeaf, parent,
							metrics[0],	vl, vr,	metrics[1],	metrics[2],	metrics[3]);
				else if (a == 2) // for NUCLEI file
					fprintf (fp, t.c_str(),
							index, K, size, nEdge, ne, ed, balance, levelFromLeaf, parent,
							metrics[0],	vl,	vr,	metrics[1],	metrics[2],	metrics[3]);
		}
};

struct helpers {
	helpers (vector<vp>* ael) {
		el = ael;
	}
	helpers (vector<vt>* atris) {
		tris = atris;
	}
	helpers () {}

	vector<vp>* el;
	vector<vt>* tris;
};


template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
	std::hash<T> hasher;
	seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
template<typename S, typename T> struct hash<pair<S, T>>
{
	inline size_t operator()(const pair<S, T> & v) const
	{
		size_t seed = 0;
		::hash_combine(seed, v.first);
		::hash_combine(seed, v.second);
		return seed;
	}
};
}

inline void clearQueue (queue<int> &q )
{
   queue<int> empty;
   swap( q, empty );
}

inline bool hashUniquify (vector<vertex>& vertices) {
	unordered_map<vertex, bool> hermap;
	for (size_t i = 0; i < vertices.size(); i++) {
		int t = vertices[i];
		if (hermap.find (t) == hermap.end())
			hermap[t] = true;
		else {
			vertices.erase (vertices.begin() + i);
			i--;
		}
	}
	if (vertices.size() > UPPERBOUND)
		return false;

	sort (vertices.begin(), vertices.end());
	return true;
}

inline int ind (int val, vector<int>& v) {
	for (size_t i = 0; i < v.size(); i++)
		if (v[i] == val)
			return i;
	return -1;
}

inline void print_time (FILE* fp, const string& str, tms t) {
	fprintf (fp, "%s %.6lf\n", str.c_str(), t.count());
	fflush(fp);
}

inline bool less_than (vertex u, vertex v, Graph& graph) {
	return (graph[u].size() < graph[v].size() || (graph[u].size() == graph[v].size() && u < v));
}

inline bool orderedConnected (Graph& graph, Graph& orderedGraph, vertex u, vertex v) {
	vertex a = u, b = v;
	if (less_than (v, u, graph))
		swap (a, b);
	for (auto c : orderedGraph[a])
		if (c == b)
			return true;
	return false;
}

inline void createOrderedIndexEdges (Graph& graph, vector<vp>& el, vector<vertex>& xel, Graph& orderedGraph) {
	xel.push_back(0);
	orderedGraph.resize(graph.size());
	for (size_t u = 0; u < graph.size(); u++) {
		for (size_t j = 0; j < graph[u].size(); j++) {
			vertex v = graph[u][j];
			if (less_than (u, v, graph)) {
				orderedGraph[u].push_back(v);
				vp c (u, v);
				el.push_back(c);
			}
		}
		xel.push_back(el.size());
	}
}

inline void createOrdered (Graph& orderedGraph, Graph& graph) {
	orderedGraph.resize (graph.size());
	for (vertex i = 0; i < graph.size(); ++i)
		for (auto u : graph[i])
			if (less_than (i, u, graph))
				orderedGraph[i].push_back(u);
}

inline vertex getEdgeId (vertex u, vertex v, vector<vertex>& xel, vector<vp>& el, Graph& graph) {

	vertex a = u, b = v;
	if (less_than (b, a, graph))
		swap (a, b);

	for (vertex i = xel[a]; i < xel[a+1]; i++)
		if (el[i].second == b)
			return i;

	printf ("getEdgeId returns -1\n");
	exit(1);
}

inline void intersection (vector<vertex>& a, vector<vertex>& b, vector<vertex>& c) {
	size_t i = 0, j = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			c.push_back(a[i]);
			i++;
			j++;
		}
	}
}

inline int checkConnectedness (Graph& graph, Graph& orderedGraph, vertex u, vertex v, vector<vertex>* xel = NULL) {

	int ret = 0;
	vertex a = u, b = v;
	if (less_than (b, a, graph)) {
		swap (a, b);
		ret -= INT_MAX;
	}
	// ret is subtracted INT_MAX if b is less than a (useful when setting TC value later)
	for (size_t k = 0; k < orderedGraph[a].size(); k++)
		if (orderedGraph[a][k] == b)
			return ret + k;

	return INT_MAX;
}

inline void get_mode (int mode, bool* modes) {
	int m = mode;
	modes[0] = modes[1] = modes[2] = modes[3] = false; // +++ +-- ++- ---
	while (m > 0) {
		if (m >= 8) {
			modes[3] = true;
			m-=8;
		}
		else if (m >= 4) {
			modes[2] = true;
			m-=4;
		}
		else if (m >= 2) {
			modes[1] = true;
			m-=2;
		}
		else if (m >= 1) {
			modes[0] = true;
			m-=1;
		}
	}
}
template <typename VtxType, typename EdgeType>
void readGraph (char *filename, vector<vector<VtxType>>& graph, vector<vector<VtxType>>& signs, EdgeType* nEdge);

int base_ktruss (Graph& graph, Graph& signs, bool hierarchy, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp, vector<vertex>* badK = NULL, int* threshold = NULL);
int base_ktruss_neutron (Graph& graph, Graph& signs, bool hierarchy, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp, vector<vertex>* badK = NULL, int* threshold = NULL);
void base_ktruss_zhao (Graph& graph, Graph& signs, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile);

int filterBad (Graph& graph, Graph& signs, bool hierarchy, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp);
int filterFromPolarizeds (Graph& graph, Graph& signs, bool hierarchy, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp);

bool computeMetrics (FILE* f, Graph& graph, Graph& signs, edge nEdge, FILE* fp);
bool computeMetricsZhao (Graph& graph, Graph& signs, edge nEdge, FILE* fp, set<vertex>& subG);

void createSkeleton (vertex u, initializer_list<vertex> neighbors, vertex* nSubcores, vector<vertex>& K, vector<subcore>& skeleton,	vector<vertex>& component, vector<vertex>& unassigned, vector<vp>& relations);
void updateUnassigned (vertex t, vector<vertex>& component, vertex* cid, vector<vp>& relations, vector<vertex>& unassigned);
void buildHierarchy (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex nVtx);
void buildHierarchyNeutron (vertex cn, vector<vp>& relations, vector<subcore>& skeleton, vertex* nSubcores, edge nEdge, vertex nVtx);
int reportSubgraph (vector<vertex>& vset, int variant, vertex index, unordered_map<vertex, vertex>& orderInFile, vector<vertex>& component, helpers& ax, vector<subcore>& skeleton, Graph& graph, Graph& signs, edge nEdge, FILE* fp, FILE* gp);
int presentNuclei (int variant, vector<subcore>& skeleton, vector<vertex>& component, Graph& graph, Graph& signs, edge* nEdge, helpers& ax, string vfile, FILE* gp);

int stats (vector<vertex>& vset, int variant, vertex index, vector<subcore>& skeleton, Graph& graph, Graph& signs, edge nEdge, FILE* fp, FILE* gp, unordered_map<vertex, int>& lr);

void past (Graph& newg, Graph& newsigns, vector<vertex>& leftRight);
