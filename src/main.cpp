#include "main.h"

int main (int argc, char *argv[]) {

	const auto t1 = chrono::steady_clock::now();
	if (argc < 2) {
		// 1:+++,  2:+--, 4:++-, 8:---; any combination is possible and expressed as the sum; e.g., 13 is 8+4+1, so ---, ++-, +++ are all accepted.
		fprintf(stderr, "usage: %s "
				"\n <filename>"
				"\n <nucleus type: 23, -23, 239, 234>"
				"\n <hierarchy?: YES or NO>"
				"\n <triangle for 23: 1(+++) 2(+--) 4(++-) 8(---) OR four-clique for 34: 6 -6 4 -4 21 22 -21 -22 01 02 03 >\n", argv[0]);
		exit(1);
	}

	char *filename = argv[1];
	string tmp (argv[1]);
	string gname = tmp.substr (tmp.find_last_of("/") + 1);

	string nd (argv[2]);
	
	if (!(nd == "234" || nd == "239" || nd == "23" || nd == "-1" || nd == "-23" || nd == "zhao")) { // -1 is for metric computations
		printf ("Invalid algorithm, options are 23, -23, 239, and 234\n");
		exit(1);
	}

	// make sure the input graph is simple: no duplicated edges; no duplicated edge with a different sign
	// read the graph, give sorted edges in graph and signs in signs
	edge nEdge = 0;
	Graph graph, signs;
	readGraph<vertex, edge> (filename, graph, signs, &nEdge);


	if (nd == "-1") {
		FILE* f = fopen (argv[3], "r");
		string t (argv[3]);
		string nFile = gname + "_" + t + "_METRICS";
		FILE* fp = fopen (nFile.c_str(), "w");

		while (1) {
			bool b = computeMetrics (f, graph, signs, nEdge, fp);
			if (!b)
				break;
		}

		fclose (f);
		fclose (fp);
		return 1;
	}

	string hrc;

	int mode;
	string vfile, out_file;
	
	if (nd == "-23" || nd == "zhao") {
		vfile = gname + "_" + nd;
	}
	else {
		string sm (argv[4]);
		mode = atoi (argv[4]);
		vfile = gname + "_" + nd + "_" + sm;
	}

	vertex maxK; // maximum K value in the graph
	vector<vertex> K;

	if (nd == "zhao") {
		base_ktruss_zhao (graph, signs, &nEdge, K, &maxK, vfile);
		return 0;
	}

	hrc = argv[3];
	bool hierarchy = (hrc == "YES" ? true : false);
	if (hierarchy)
		out_file = vfile + "_Hierarchy";
	else
		out_file = vfile + "_K";

	FILE* fp = fopen (out_file.c_str(), "w");

	if (nd == "239") {
		filterBad (graph, signs, hierarchy, mode, &nEdge, K, &maxK, vfile, fp);
	}
	else if (nd == "234") {
		filterFromPolarizeds (graph, signs, hierarchy, mode, &nEdge, K, &maxK, vfile, fp);
	}
	else if (nd == "23") {
		base_ktruss (graph, signs, hierarchy, mode, &nEdge, K, &maxK, vfile, fp);
	}
	else if (nd == "-23") {
		int minTC = base_ktruss_neutron (graph, signs, hierarchy, &nEdge, K, &maxK, vfile, fp);
		#ifdef DUMP_K
			string kfile = vfile + "_K_values";
			FILE* kf = fopen (kfile.c_str(), "w");
			for (vertex i = 0; i < K.size(); i++)
				fprintf (kf, "%d\n", K[i] + minTC);
			fclose (kf);
		#endif
		const auto t2 = chrono::steady_clock::now();
		printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s-nucleus: %d\n", gname.c_str(), graph.size(), nEdge, nd.c_str(), maxK);
		print_time (fp, "End-to-end Time: ", t2 - t1);
		fclose (fp);
		return 0;
	}

#ifdef DUMP_K
	string kfile = vfile + "_K_values";
	FILE* kf = fopen (kfile.c_str(), "w");
	for (vertex i = 0; i < K.size(); i++)
		fprintf (kf, "%lld\n", K[i]);
	fclose (kf);
#endif

	const auto t2 = chrono::steady_clock::now();
	printf ("%s\t|V|: %d\t|E|: %d\tmaxK for %s-nucleus: %d\n", gname.c_str(), graph.size(), nEdge, nd.c_str(), maxK);
	print_time (fp, "End-to-end Time: ", t2 - t1);
	fclose (fp);

	return 0;
}
