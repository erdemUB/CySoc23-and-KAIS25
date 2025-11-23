#include "main.h"

#define MAXLINE 1000000

template <typename VtxType, typename EdgeType>
void readEdgeList (bool mm, char *filename, Graph& graph, Graph& signs,EdgeType* nEdge) {
	char* line = (char*) malloc (sizeof (char) * MAXLINE);
	FILE* matfp = fopen(filename, "r");

	// skip comments
	do
		fgets(line, MAXLINE, matfp);
	while (line[0] == '%' || line[0] == '#');

	VtxType nVtx;
	vector<couple> coords;
	VtxType u, v, nv = 0;
	stringstream ss (line);
	if (mm) {
		ss >> nVtx >> *nEdge;
		printf ("|V|: %d   |E|: %d\n", nVtx, *nEdge);
	}
	else {
		ss >> u >> v;
		nv = max (nv, (max (u, v)));
		if (u != v) {
			coords.push_back (make_tuple (u, v));
			coords.push_back (make_tuple (v, u));
		}
	}

	while (fgets(line, MAXLINE, matfp)) {
		stringstream ss (line);
		ss >> u >> v;
		nv = max (nv, (max (u, v)));
		if (u != v) {
			coords.push_back (make_tuple (u, v));
			coords.push_back (make_tuple (v, u));
		}
	}
	fclose (matfp);

	if (mm) {
		if (nVtx != nv + 1) {
			printf ("nVtx in header (%d) is wrong, must be %d\n", nVtx, nv + 1);
			nVtx = nv + 1;
		}
	}
	else
		nVtx = nv + 1;

	sort (coords.begin(), coords.end());

	// begin constructing graph
	graph.resize (nVtx);
	EdgeType i = 0;
	graph[get<0>(coords[i])].push_back(get<1>(coords[i]));
	for (i = 1; i < coords.size(); i++)
		if (coords[i] != coords[i-1])
			graph[get<0>(coords[i])].push_back(get<1>(coords[i]));

	// sort each neighbor list
	EdgeType ne = 0;
	for (auto v : graph) {
		sort (v.begin(), v.end());
		ne += v.size();
	}
	*nEdge = ne / 2;

	if (mm) {
		if (*nEdge != ne / 2) {
			printf ("nEdge in header (%d) is wrong, must be %d\n", *nEdge, ne/2);
			*nEdge = ne / 2;
		}
	}
	else
		*nEdge = ne / 2;


	// Read Signs
	matfp = fopen(filename, "r");

	signs.resize(graph.size());
	for (vertex i = 0; i < graph.size(); i++)
		signs[i].resize (graph[i].size(), 0);

	// skip comments
	do
		fgets(line, MAXLINE, matfp);
	while (line[0] == '%' || line[0] == '#');

	int sign;
	if (mm) {
		ss >> nVtx >> *nEdge;
		printf ("|V|: %d   |E|: %d\n", nVtx, *nEdge);
	}
	else {
		stringstream ss (line);
		ss >> u >> v >> sign;
		signs[u][ind (v, graph[u])] = sign;
		signs[v][ind (u, graph[v])] = sign;
	}
	while (fgets(line, MAXLINE, matfp)) {
		stringstream ss (line);
		ss >> u >> v >> sign;
		signs[u][ind (v, graph[u])] = sign;
		signs[v][ind (u, graph[v])] = sign;
	}
	fclose (matfp);

	if (mm) {
		if (nVtx != nv + 1) {
			printf ("nVtx in header (%d) is wrong, must be %d\n", nVtx, nv + 1);
			nVtx = nv + 1;
		}
	}
	else
		nVtx = nv + 1;
}

template <typename VtxType, typename EdgeType>
void readGraph (char *filename, vector<vector<VtxType>>& graph, vector<vector<VtxType>>& signs, EdgeType* nEdge) {
	string st (filename);
	string gname = st.substr (st.find_last_of("/") + 1);
	int idx = gname.find_last_of(".");
	string ext = gname.substr(idx);
	readEdgeList<VtxType, EdgeType> (false, filename, graph, signs, nEdge);
	return;
}

template void readGraph (char *filename, vector<vector<vertex>>& graph, vector<vector<vertex>>& signs, edge* nEdge);
