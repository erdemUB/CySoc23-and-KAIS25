#include "main.h"

// per edge
lol countTriangles (Graph& graph, Graph& signs, bool* modes, Graph& orderedGraph, Graph& TC) {

	lol tc = 0;
	int ppp = 0, ppn = 0, pnn = 0, nnn = 0;
	int s1, s2, s3, sum, ab;
	for (size_t i = 0; i < orderedGraph.size(); i++) {
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orderedGraph[i].size(); k++) {
				vertex a = orderedGraph[i][j];
				vertex b = orderedGraph[i][k];
				int c = checkConnectedness (graph, orderedGraph, a, b);

				if (c == INT_MAX)
					continue;

				if (c < 0) {
					c += INT_MAX;
					ab = b;
				}
				else {
					ab = a;
				}

				s1 = signs[i][ind (a, graph[i])]; // i-a
				s2 = signs[i][ind (b, graph[i])]; // i-b
				s3 = signs[a][ind (b, graph[a])]; // a-b
				sum = s1 + s2 + s3;

				if ((modes[0] && sum == 3) || (modes[1] && sum == -1) || (modes[2] && sum == 1) || (modes[3] && sum == -3)) {
					TC[i][j]++;
					TC[i][k]++;
					TC[ab][c]++;
					tc++;
				}

				if (sum == 3)
					ppp++;
				else if (sum == 1)
					ppn++;
				else if (sum == -1)
					pnn++;
				else if (sum == -3)
					nnn++;

			}
		}
	}

	int a = ppp + ppn + pnn + nnn;
	printf ("total # triangles: %d\n",a );

	return tc;
}

int base_ktruss (Graph& graph, Graph& signs, bool hierarchy, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp, vector<vertex>* badK, int* threshold) {

	bool modes[4];
	get_mode (mode, modes); // set the modes flags: 0: +++ 1: +-- 2: ++- 3: ---

	const auto t1 = chrono::steady_clock::now();
	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

	lol tric;
	// Compute
	tric = countTriangles (graph, signs, modes, orderedGraph, TC);

	fprintf (fp, "# triangles: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Triangle counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	printf ("nEdge: %d\n", *nEdge);
	K.resize (*nEdge, -1);

	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, *nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			if (TC[i][j] > 0)
				nBucket.Insert (id++, TC[i][j]);
			else
				K[id++] = 0;
		}

	vertex tc_e = 0;

	// required for hierarchy
	vertex cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<vertex> component; // subcore ids for each vertex
	vector<vp> relations;
	vector<vertex> unassigned;
	vertex nSubcores;

	if (hierarchy) {
		cid = 0;
		nSubcores = 0;
		component.resize (*nEdge, -1);
	}

	vertex monitor = 0, s1, s2, s3, sum;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;
#ifdef MONITOR
		if (monitor % 10000 == 0)
			printf ("e: %d    val: %d    counter: %d  nEdge: %d\n", e, val, monitor, *nEdge);
		monitor++;
#endif
		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		tc_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		s1 = signs[u][ind (v, graph[u])]; // sign of u-v
		vector<vertex> commonNeighbors;
		intersection (graph[u], graph[v], commonNeighbors);
		for (auto w : commonNeighbors) { // decrease the TC of the neighbor edges with greater TC
			edge f = getEdgeId (u, w, xel, el, graph);
			edge g = getEdgeId (v, w, xel, el, graph);
			s2 = signs[u][ind (w, graph[u])]; // u-w
			s3 = signs[v][ind (w, graph[v])]; // v-w
			sum = s1 + s2 + s3;

			if ((modes[0] && sum == 3) || (modes[1] && sum == -1) || (modes[2] && sum == 1) || (modes[3] && sum == -3)) { // change these to accept multiple types of signed triangles
				if (K[f] == -1 && K[g] == -1) {
					if (nBucket.CurrentValue(f) > tc_e)
						nBucket.DecVal(f);
					if (nBucket.CurrentValue(g) > tc_e)
						nBucket.DecVal(g);
				}
				else if (hierarchy)
					createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxtruss = tc_e;
	printf ("maxt: %d\n", *maxtruss);

	const auto p2 = chrono::steady_clock::now();

	int numOfSubgraphs = -1;
	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}
	else {
	
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchy (*maxtruss, relations, skeleton, &nSubcores, *nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 2,3 nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);
		numOfSubgraphs = presentNuclei (23, skeleton, component, graph, signs, nEdge, hp, vfile, fp);

		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 2,3 nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
		
	}

	return numOfSubgraphs;
}

// per edge
lol countTrianglesNeutron (Graph& graph, Graph& signs, Graph& orderedGraph, Graph& TC) {

	lol tc = 0;
	int ppp = 0, ppn = 0, pnn = 0, nnn = 0;
	int s1, s2, s3, sum, ab;
	for (size_t i = 0; i < orderedGraph.size(); i++) {
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orderedGraph[i].size(); k++) {
				vertex a = orderedGraph[i][j];
				vertex b = orderedGraph[i][k];
				int c = checkConnectedness (graph, orderedGraph, a, b);

				if (c == INT_MAX)
					continue;

				if (c < 0) {
					c += INT_MAX;
					ab = b;
				}
				else {
					ab = a;
				}

				s1 = signs[i][ind (a, graph[i])]; // i-a
				s2 = signs[i][ind (b, graph[i])]; // i-b
				s3 = signs[a][ind (b, graph[a])]; // a-b
				sum = s1 + s2 + s3;

				if (sum == 3 || sum == -1) {
					TC[i][j]++;
					TC[i][k]++;
					TC[ab][c]++;
					tc++;
				}
				else if (sum == 1 || sum == -3) {
					TC[i][j]--;
					TC[i][k]--;
					TC[ab][c]--;
					tc++;
				}

				if (sum == 3)
					ppp++;
				else if (sum == 1)
					ppn++;
				else if (sum == -1)
					pnn++;
				else if (sum == -3)
					nnn++;

			}
		}
	}

	int a = ppp + ppn + pnn + nnn;
	printf ("total # triangles: %d\n",a );
	return tc;
}

int base_ktruss_neutron (Graph& graph, Graph& signs, bool hierarchy, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp, vector<vertex>* badK, int* threshold) {


	const auto t1 = chrono::steady_clock::now();
	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

	lol tric;
	// Compute
	tric = countTrianglesNeutron (graph, signs, orderedGraph, TC);

	fprintf (fp, "# triangles: %lld\n", tric);
	const auto t2 = chrono::steady_clock::now();
	print_time (fp, "Triangle counting: ", t2 - t1);

	// Peeling
	const auto p1 = chrono::steady_clock::now();
	printf ("nEdge: %d\n", *nEdge);
	K.resize (*nEdge, -1);

	int minTC = 0;

	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++)
			if (TC[i][j] < minTC)
				minTC = TC[i][j];

	if (minTC > 0)
		minTC = 0;

	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, *nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			if (TC[i][j] - minTC > 0)
				nBucket.Insert (id++, TC[i][j] - minTC);
			else
				K[id++] = 0;
		}

	vertex tc_e = 0;

	// required for hierarchy
	vertex cid; // subcore id number
	vector<subcore> skeleton; // equal K valued cores
	vector<vertex> component; // subcore ids for each vertex
	vector<vp> relations;
	vector<vertex> unassigned;
	vertex nSubcores;

	if (hierarchy) {
		cid = 0;
		nSubcores = 0;
		component.resize (*nEdge, -1);
	}

	vertex monitor = 0, s1, s2, s3, sum;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;
#ifdef MONITOR
		if (monitor % 10000 == 0)
			printf ("e: %d    val: %d    counter: %d  nEdge: %d\n", e, val, monitor, *nEdge);
		monitor++;
#endif
		if (hierarchy) {
			unassigned.clear();
			subcore sc (val);
			skeleton.push_back (sc);
		}

		tc_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		s1 = signs[u][ind (v, graph[u])]; // sign of u-v
		vector<vertex> commonNeighbors;
		intersection (graph[u], graph[v], commonNeighbors);
		for (auto w : commonNeighbors) { // update the TC of the neighbor edges with greater TC
			edge f = getEdgeId (u, w, xel, el, graph);
			edge g = getEdgeId (v, w, xel, el, graph);
			s2 = signs[u][ind (w, graph[u])]; // u-w
			s3 = signs[v][ind (w, graph[v])]; // v-w
			sum = s1 + s2 + s3;
			if (sum == 3 || sum == -1) { // change these to accept multiple types of signed triangles
				if (K[f] == -1 && K[g] == -1) {
					if (nBucket.CurrentValue(f) > tc_e)
						nBucket.DecVal(f);
					if (nBucket.CurrentValue(g) > tc_e)
						nBucket.DecVal(g);
				}
				else if (hierarchy)
					createSkeleton (e, {f, g}, &nSubcores, K, skeleton, component, unassigned, relations);
			}
		}

		if (hierarchy)
			updateUnassigned (e, component, &cid, relations, unassigned);
	}

	nBucket.Free();
	*maxtruss = tc_e;
	printf ("maxt: %d\n", *maxtruss);

	const auto p2 = chrono::steady_clock::now();

	int numOfSubgraphs = -1;
	if (!hierarchy) {
		print_time (fp, "Only peeling time: ", p2 - p1);
		print_time (fp, "Total time: ", (p2 - p1) + (t2 - t1));
	}
	else {
		print_time (fp, "Only peeling + on-the-fly hierarchy construction time: ", p2 - p1);
		const auto b1 = chrono::steady_clock::now();
		buildHierarchyNeutron (*maxtruss, relations, skeleton, &nSubcores, *nEdge, nVtx);
		const auto b2 = chrono::steady_clock::now();

		print_time (fp, "Building hierarchy time: ", b2 - b1);
		print_time (fp, "Total 2,3 nucleus decomposition time (excluding density computation): ", (p2 - p1) + (t2 - t1) + (b2 - b1));

		fprintf (fp, "# subcores: %d\t\t # subsubcores: %d\t\t |V|: %d\n", nSubcores, skeleton.size(), graph.size());

		const auto d1 = chrono::steady_clock::now();
		helpers hp (&el);

		for (subcore& s : skeleton)
			s.K = s.K + minTC;
		numOfSubgraphs = presentNuclei (23, skeleton, component, graph, signs, nEdge, hp, vfile, fp);

		const auto d2 = chrono::steady_clock::now();

		print_time (fp, "Total 2,3 nucleus decomposition time: ", (p2 - p1) + (t2 - t1) + (b2 - b1) + (d2 - d1));
	}
	return minTC;
}

// for avoiding bad triangles
void filterBadPreProcessing (Graph& graph, Graph& signs, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss) {

	bool modes[4];
	get_mode (mode, modes); // set the modes flags: 0: +++ 1: +-- 2: ++- 3: ---

	vertex nVtx = graph.size();

	// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
	vector<vp> el;
	vector<vertex> xel;
	Graph orderedGraph;
	createOrderedIndexEdges (graph, el, xel, orderedGraph);

	// Triangle counting for each edge
	vector<vector<vertex> > TC (nVtx);
	for (vertex i = 0; i < nVtx; i++)
		TC[i].resize (orderedGraph[i].size(), 0);

	lol tric;
	// Compute
	tric = countTriangles (graph, signs, modes, orderedGraph, TC);

	// Peeling
	K.resize (*nEdge, -1);

	Naive_Bucket nBucket;
	nBucket.Initialize (nVtx, *nEdge); // maximum triangle count of an edge is nVtx
	vertex id = 0;
	for (size_t i = 0; i < orderedGraph.size(); i++)
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			if (TC[i][j] > 0)
				nBucket.Insert (id++, TC[i][j]);
			else
				K[id++] = 0;
		}

	vertex tc_e = 0;

	vertex s1, s2, s3, sum;
	while (true) {
		edge e;
		vertex val;
		if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
			break;
		tc_e = K[e] = val;

		vertex u = el[e].first;
		vertex v = el[e].second;
		s1 = signs[u][ind (v, graph[u])]; // sign of u-v
		vector<vertex> commonNeighbors;
		intersection (graph[u], graph[v], commonNeighbors);
		for (auto w : commonNeighbors) { // decrease the TC of the neighbor edges with greater TC
			edge f = getEdgeId (u, w, xel, el, graph);
			edge g = getEdgeId (v, w, xel, el, graph);
			s2 = signs[u][ind (w, graph[u])]; // u-w
			s3 = signs[v][ind (w, graph[v])]; // v-w
			sum = s1 + s2 + s3;

			if ((modes[0] && sum == 3) || (modes[1] && sum == -1) || (modes[2] && sum == 1) || (modes[3] && sum == -3)) { // change these to accept multiple types of signed triangles
				if (K[f] == -1 && K[g] == -1) {
					if (nBucket.CurrentValue(f) > tc_e)
						nBucket.DecVal(f);
					if (nBucket.CurrentValue(g) > tc_e)
						nBucket.DecVal(g);
				}
			}
		}
	}

	nBucket.Free();
	*maxtruss = tc_e;

	float threshold = *maxtruss * PHOTON_THRESHOLD;

	// Filter the endpoints of edges based on their unbalanced atom numbers and the threshold (PHOTON_THRESHOLD)
	for (edge e = 0; e < *nEdge; ++e) {
		if (K[e] >= threshold) {
			vertex u = el[e].first;
			vertex v = el[e].second;

			for (vertex w : graph[u]) {
				vector<int>::iterator pos = std::find(graph[w].begin(), graph[w].end(), u);
				int idx = pos - graph[w].begin();
				graph[w].erase(pos);
				signs[w].erase(signs[w].begin() + idx);
			}

			for (vertex w : graph[v]) {
				vector<int>::iterator pos = std::find(graph[w].begin(), graph[w].end(), v);
				int idx = pos - graph[w].begin();
				graph[w].erase(pos);
				signs[w].erase(signs[w].begin() + idx);
			}

			graph[u].clear();
			signs[u].clear();
			graph[v].clear();
			signs[v].clear();
		}
	}
}

// for avoiding bad triangles
int filterBad (Graph& graph, Graph& signs, bool hierarchy, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	vector<vertex> badK;
	vertex max_badK;

	filterBadPreProcessing (graph, signs, 12, nEdge, badK, &max_badK);
	return base_ktruss (graph, signs, hierarchy, mode, nEdge, K, maxtruss, vfile, fp);
}

// for tension enhancement
void filterPolarizedsPreProcessing (Graph& graph, Graph& signs) {

	vertex nVtx = graph.size();

	unordered_map<vertex, double> tri;
    unordered_map<vertex, bool> removed;
	for (vertex v = 0; v < nVtx; ++v) {
		tri.emplace(v, 0);
		removed.emplace(v, false);
	}

	// counts number of +-- minus unbalanced triangles for each node
	for (int u = 0; u < nVtx; ++u) {
		for (int j = 0; j < graph[u].size(); ++j) {
			int v = graph[u][j];
			if (u < v) {
				vector<vertex> commonNeighbors;
				int s1 = signs[u][ind (v, graph[u])]; // sign of u-v
				intersection (graph[u], graph[v], commonNeighbors);
				for (int w : commonNeighbors) {
					if (v < w) {
						int s2 = signs[u][ind (w, graph[u])]; // u-w
						int s3 = signs[v][ind (w, graph[v])]; // v-w
						int sum = s1 + s2 + s3;
						if (sum == -1) {
							tri[u]++;
							tri[v]++;
							tri[w]++;
						}
						else if (sum != 3){
							tri[u]--;
							tri[v]--;
							tri[w]--;
						}
					}
				}
			}
		}
	}
	
	double n = graph.size();
	
	// peeling process - iteratively peels node with minimum friction until minimum friction >= ELECTRON_THRESHOLD
	while (n > 0) {
	
    	double minVal = 2;
    	std::vector<int> minNodes;
    	
    	double maxTriDeg = (n - 1) * (n - 2) / 2;
    	
    	for (auto it : tri) {
    	    if (!removed[it.first]) {
				// friction
        	    double val = (it.second) / maxTriDeg;
        	    if (val < minVal) {
        	        minVal = val;
					minNodes.clear();
					minNodes.push_back(it.first);
        	    }
				else if (val == minVal) {
					minNodes.push_back(it.first);
				}
    	    }
    	}

    	if (minVal < ELECTRON_THRESHOLD) {
			for (int u : minNodes) {

				n--;
				removed[u] = true;
				
				for (vertex v : graph[u]) {
					if (!removed[v]) {
						vector<vertex> commonNeighbors;
						int s1 = signs[u][ind (v, graph[u])]; // sign of u-v
						intersection (graph[u], graph[v], commonNeighbors);
						for (int w : commonNeighbors) {
							if (v < w && !removed[w]) {
								int s2 = signs[u][ind (w, graph[u])]; // u-w
								int s3 = signs[v][ind (w, graph[v])]; // v-w
								int sum = s1 + s2 + s3;
								if (sum == -1) {
									tri[v]--;
									tri[w]--;
								}
								else if (sum != 3){
									tri[v]++;
									tri[w]++;
								}
							}
						}
					}
				}

				for (vertex w : graph[u]) {
					vector<int>::iterator pos = std::find(graph[w].begin(), graph[w].end(), u);
					int idx = pos - graph[w].begin();
					graph[w].erase(pos);
					signs[w].erase(signs[w].begin() + idx);
				}

				graph[u].clear();
				signs[u].clear();
			}
    	}
    	else {
    	    break;
    	}
	}
}

// for tension enhancement
int filterFromPolarizeds (Graph& graph, Graph& signs, bool hierarchy, int mode, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile, FILE* fp) {

	filterPolarizedsPreProcessing (graph, signs);
	return base_ktruss (graph, signs, hierarchy, mode, nEdge, K, maxtruss, vfile, fp);
}

// per edge
lol countTrianglesZhao (Graph& graph, Graph& signs, Graph& orderedGraph, Graph& TC) {

	lol tc = 0;
	int ppp = 0, ppn = 0, pnn = 0, nnn = 0;
	int s1, s2, s3, sum, ab;
	for (size_t i = 0; i < orderedGraph.size(); i++) {
		for (size_t j = 0; j < orderedGraph[i].size(); j++) {
			for (size_t k = j + 1; k < orderedGraph[i].size(); k++) {
				vertex a = orderedGraph[i][j];
				vertex b = orderedGraph[i][k];
				int c = checkConnectedness (graph, orderedGraph, a, b);

				if (c == INT_MAX)
					continue;

				if (c < 0) {
					c += INT_MAX;
					ab = b;
				}
				else {
					ab = a;
				}

				s1 = signs[i][ind (a, graph[i])]; // i-a
				s2 = signs[i][ind (b, graph[i])]; // i-b
				s3 = signs[a][ind (b, graph[a])]; // a-b
				sum = s1 + s2 + s3;

				if (sum == 3 || sum == -1) {
					TC[i][j]++;
					TC[i][k]++;
					TC[ab][c]++;
					tc++;
				}

				if (sum == 3)
					ppp++;
				else if (sum == 1)
					ppn++;
				else if (sum == -1)
					pnn++;
				else if (sum == -3)
					nnn++;

			}
		}
	}
	return tc;
}

void base_ktruss_zhao (Graph& graph, Graph& signs, edge* nEdge, vector<vertex>& K, vertex* maxtruss, string vfile) {

	int numSubGs = 10;

	Graph graphCopy = graph;
	Graph signsCopy = signs;
	edge nEdgeCopy = *nEdge;

	string nFile = vfile + "_NUCLEI";
	FILE* fp = fopen (nFile.c_str(), "w");

	while (nEdge > 0 && numSubGs > 0) {
		vertex nVtx = graph.size();
		K.clear();

		// Create directed graph from low degree vertices to higher degree vertices AND prepare a CSR-like structure to index the edges
		vector<vp> el;
		vector<vertex> xel;
		Graph orderedGraph;
		createOrderedIndexEdges (graph, el, xel, orderedGraph);

		// Triangle counting for each edge
		vector<vector<vertex> > TC (nVtx);
		for (vertex i = 0; i < nVtx; i++)
			TC[i].resize (orderedGraph[i].size(), 0);

		lol tric;
		// Compute
		tric = countTrianglesZhao (graph, signs, orderedGraph, TC);

		// Peeling
		K.resize (*nEdge, -1);

		Naive_Bucket nBucket;
		nBucket.Initialize (nVtx, *nEdge); // maximum triangle count of an edge is nVtx
		vertex id = 0;
		for (size_t i = 0; i < orderedGraph.size(); i++)
			for (size_t j = 0; j < orderedGraph[i].size(); j++) {
				if (TC[i][j] > 0)
					nBucket.Insert (id++, TC[i][j]);
				else
					K[id++] = 0;
			}

		vertex tc_e = 0;

		vertex monitor = 0, s1, s2, s3, sum;
		while (true) {
			edge e;
			vertex val;
			if (nBucket.PopMin(&e, &val) == -1) // if the bucket is empty
				break;

			tc_e = K[e] = val;

			vertex u = el[e].first;
			vertex v = el[e].second;
			s1 = signs[u][ind (v, graph[u])]; // sign of u-v
			vector<vertex> commonNeighbors;
			intersection (graph[u], graph[v], commonNeighbors);
			for (auto w : commonNeighbors) { // update the TC of the neighbor edges with greater TC
				edge f = getEdgeId (u, w, xel, el, graph);
				edge g = getEdgeId (v, w, xel, el, graph);
				s2 = signs[u][ind (w, graph[u])]; // u-w
				s3 = signs[v][ind (w, graph[v])]; // v-w
				sum = s1 + s2 + s3;
				if (sum == 3 || sum == -1) { // change these to accept multiple types of signed triangles
					if (K[f] == -1 && K[g] == -1) {
						if (nBucket.CurrentValue(f) > tc_e)
							nBucket.DecVal(f);
						if (nBucket.CurrentValue(g) > tc_e)
							nBucket.DecVal(g);
					}
				}
			}
		}

		nBucket.Free();

		// ZHAO on maximum balanced atom number
		vector<bool> edges(*nEdge);
		for (int i = 0; i < nVtx; ++i) {
			for (int o : graph[i]) {
				if (o > i) {
					edge g = getEdgeId (i, o, xel, el, graph);
					if (K[g] >= tc_e)
						edges[g] = true;
				}
			}
		}
		while (1) {
			vector<bool> unbalEdges(*nEdge);
			vector<bool> unvisited(*nEdge);
			for (edge e = 0; e < *nEdge; ++e) {
				if (edges[e]) {
					vertex u = el[e].first;
					vertex v = el[e].second;
					s1 = signs[u][ind (v, graph[u])]; // sign of u-v
					vector<vertex> commonNeighbors;
					intersection (graph[u], graph[v], commonNeighbors);
					for (auto w : commonNeighbors) {
						edge f = getEdgeId (u, w, xel, el, graph);
						edge g = getEdgeId (v, w, xel, el, graph);
						if (edges[f] && edges[g]) {
							s2 = signs[u][ind (w, graph[u])]; // u-w
							s3 = signs[v][ind (w, graph[v])]; // v-w
							sum = s1 + s2 + s3;
							if (sum == -3 || sum == 1) {
								unbalEdges[f] = true;
								unbalEdges[g] = true;
								unbalEdges[e] = true;
								unvisited[f] = true;
								unvisited[g] = true;
								unvisited[e] = true;
							}
						}
					}
				}
			}

			// Support group construction
			vector<vector<int>> G;
			vector<int> GID(*nEdge);
			int id = 0;
			for (edge e = 0; e < *nEdge; ++e) {
				if (unbalEdges[e] && unvisited[e] && K[e] == tc_e) {
					unvisited[e] = false;
					vector<int> supp;
					queue<int> Q;
					Q.push(e);
					while (!Q.empty()) {
						edge e2 = Q.front();
						Q.pop();
						supp.push_back(e2);
						GID[e2] = id;

						vertex u = el[e2].first;
						vertex v = el[e2].second;
						s1 = signs[u][ind (v, graph[u])]; // sign of u-v
						vector<vertex> commonNeighbors;
						intersection (graph[u], graph[v], commonNeighbors);
						for (auto w : commonNeighbors) {
							edge f = getEdgeId (u, w, xel, el, graph);
							edge g = getEdgeId (v, w, xel, el, graph);
							if (edges[f] && edges[g]) {
								s2 = signs[u][ind (w, graph[u])]; // u-w
								s3 = signs[v][ind (w, graph[v])]; // v-w
								sum = s1 + s2 + s3;
								if (sum == 3 || sum == -1) {
									if (unvisited[f] && K[f] == tc_e) {
										unvisited[f] = false;
										Q.push(f);
									}
									if (unvisited[g] && K[g] == tc_e) {
										unvisited[g] = false;
										Q.push(g);
									}
								}
							}
						}

					}
					G.push_back(supp);
					id++;
				}
			}


			int minSuppSz = 0;
			int minID = 0;
			int gSz = G.size();
			for (int id = 0; id < gSz; ++id) {
				int suppSz = G[id].size();
				if (minSuppSz == 0 || suppSz < minSuppSz) {
					minSuppSz = suppSz;
					minID = id;
				}
			}
			if (minSuppSz == 0)
				break;

			for (edge e : G[minID]) {
				edges[e] = false;
				vertex u = el[e].first;
				vertex v = el[e].second;
				s1 = signs[u][ind (v, graph[u])]; // sign of u-v
				vector<vertex> commonNeighbors;
				intersection (graph[u], graph[v], commonNeighbors);
				for (auto w : commonNeighbors) {
					edge f = getEdgeId (u, w, xel, el, graph);
					edge g = getEdgeId (v, w, xel, el, graph);
					if (edges[f] && edges[g]) {
						s2 = signs[u][ind (w, graph[u])]; // u-w
						s3 = signs[v][ind (w, graph[v])]; // v-w
						sum = s1 + s2 + s3;
						if (sum == 3 || sum == -1) {
							K[f] -= 1;
							K[g] -= 1;
						}
					}
				}
			}
		}

		set<vertex> subG;
		int numE = 0;
		for (edge e = 0; e < *nEdge; ++e) {
			if (edges[e]) {
				numE++;
				vertex u = el[e].first;
				vertex v = el[e].second;
				subG.insert(u);
				subG.insert(v);
				int u_v = ind (v, graph[u]);
				int v_u = ind (u, graph[v]);
				graph[u].erase(graph[u].begin() + u_v);
				signs[u].erase(signs[u].begin() + u_v);
				graph[v].erase(graph[v].begin() + v_u);
				signs[v].erase(signs[v].begin() + v_u);
			}
		}

		computeMetricsZhao (graphCopy, signsCopy, nEdgeCopy, fp, subG);
		
		*nEdge -= numE;

		numSubGs -= 1;
	}
	fclose (fp);
}
