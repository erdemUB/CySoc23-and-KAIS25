#include "main.h"

double color (double ed) {
        double a = (1 - ed) / 0.25;
        int X = floor(a);
        int Y = floor(255 * (a-X));
        if (X < 4)
                return (3 - X) + (255 - (double) Y) / 255;
        else
                return 0;
}

void print_nested_circle (vector<subcore>& hrc, int ind, FILE* fp, string cfl) {
	if (hrc[ind].size < LOWERBOUND)
            return;
        double parent_color = color(hrc[ind].ed);
        fprintf(fp, "{\"color\": %lf, \"fl\": \"%s\", \"index\": \"%d\", \"name\": \"%ld %.2lf (%d)\", \"size\": %ld",
                        parent_color, cfl.c_str(), ind, hrc[ind].size, hrc[ind].ed, hrc[ind].K, hrc[ind].size);
        if (hrc[ind].children.size() == 1) {
                fprintf(fp, ", \"children\": [\n");
                int ch = hrc[ind].children[0];
                // ghost child
                fprintf(fp, "{\"color\": %lf, \"fl\": \"\", \"name\": \"\", \"size\": %ld}, ",
                                parent_color, hrc[hrc[ch].parent].size - hrc[ch].size);
                // real child
                print_nested_circle (hrc, ch, fp, cfl);
                fprintf(fp, "\n]\n");
        }
        else if (hrc[ind].children.size() > 1) {
                fprintf(fp, ", \n\"children\": [\n");
                size_t i;
                for (i = 0; i < hrc[ind].children.size() - 1; i++) {
                        print_nested_circle (hrc, hrc[ind].children[i], fp, cfl);
                        fprintf(fp, ",\n");
                }
                print_nested_circle (hrc, hrc[ind].children[i], fp, cfl);
                fprintf(fp, "\n]");
        }
        fprintf(fp, "}\n");
}

inline vertex commons (vector<vertex>& a, vector<vertex>& b) {
	vertex i = 0, j = 0;
	vertex count = 0;
	while (i < a.size() && j < b.size()) {
		if (a[i] < b[j])
			i++;
		else if (b[j] < a[i])
			j++;
		else {
			count++;
			i++;
			j++;
		}
	}
	return count;
}

bool pullChildrenSets (int variant, FILE* fp, vertex index, unordered_map<vertex, vertex>& orderInFile, vector<vertex>& vset, vector<subcore>& skeleton) {

	int limit = UPPERBOUND;
	char c;
	for (vertex eda : skeleton[index].children) {
		if (skeleton[eda].size == -1)
			return false;

		if (orderInFile.find (eda) == orderInFile.end()) {
			printf ("PROBLEM: %d has -1 as order\n", eda);
			exit(1);
		}

		vertex sc = orderInFile[eda];
		fseek (fp, 0, SEEK_SET);
		vertex ln = 0;
		if (sc != 0) {
			do {
				c = fgetc (fp);
				if (c == '\n') {
					ln++;
					if (ln == sc)
						break;
				}
			} while (c != EOF);
		}

		// now you are at the correct line of the file and can get the vertices of eda
		int d;
		double f;

		//           id, K, |V|, |E|,|E-| ed, balance,LEAF?,parentId, polarity, |V_l|, |V_r|, harmonyLeft, harmonyRight, tension
		fscanf (fp, "%d  %d  %d  %d  %d  %lf  %lf  %d  %d  %lf  %d  %d  %lf  %lf  %lf",
					 &d, &d, &d, &d, &d, &f,   &f, &d, &d, &f,  &d, &d, &f,  &f,  &f);
		bool flag = false; // set to true when the first -1 is seen
		while (fscanf (fp, "%d", &d) != EOF) {
			if (d != -1) {
				vset.push_back (d);
				if (vset.size() > limit) {
					fseek (fp, 0, SEEK_END);
					return false;
				}
			}
			else if (!flag)
				flag = true;
			else if (flag)
				break;
		}
		fseek (fp, 0, SEEK_END);
	}
	return true;
}

inline void dummyLine (subcore* sc, FILE* fp, vertex index) {
	sc->size = -1;
	sc->ed = -1;
	sc->nEdge = -1;
	sc->parent = -1;
	sc->print (2, index, fp);
	fprintf (fp, "-1\n");
}

inline void removeChild (vertex i, vector<subcore>& backup) {
	if (backup[i].parent != -1) {
		vertex pr = backup[i].parent;
		for (vertex j = 0; j < backup[pr].children.size(); j++)
			if (backup[pr].children[j] == i) {
				backup[pr].children.erase (backup[pr].children.begin() + j);
				break;
			}
	}
}

void rearrange (vector<subcore>& skeleton) { // rearrange children and parents based on visibility
	for (size_t i = 0; i < skeleton.size(); i++)
		skeleton[i].children.clear();
	for (size_t i = 0; i < skeleton.size() - 1; i++) {
		if (skeleton[i].visible) {
			int pr = skeleton[i].parent;
			while (!skeleton[pr].visible)
				pr = skeleton[pr].parent;
			skeleton[i].parent = pr;
			skeleton[pr].children.push_back (i);
		}
	}
}



void polarityMetrics (Graph& newg, Graph& newsigns, vector<vertex>& leftRight, subcore* sc) {
	double score = 0;
	int leftSame = 0, leftDiff = 0, rightSame = 0, rightDiff = 0, acrossDiff = 0, acrossSame = 0;
	int vLeft = 0, vRight = 0;
	for (auto u = 0; u < newg.size(); u++) {
		for (auto j = 0; j < newg[u].size(); j++) {
			vertex v = newg[u][j];
			if (leftRight[u] == leftRight[v]) {
				if (leftRight[u] == -1) {
					if (newsigns[u][j] == 1)
						leftSame++;
					else
						leftDiff++;
				}
				else {
					if (newsigns[u][j] == 1)
						rightSame++;
					else
						rightDiff++;
				}
			}
			else {
				if (newsigns[u][j] == 1)
					acrossSame++;
				else
					acrossDiff++;
			}
		}
		if (leftRight[u] == -1)
			vLeft++;
	}

	leftSame /= 2;
	rightSame /= 2;
	leftDiff /= 2;
	rightDiff /= 2;
	acrossDiff /= 2;
	acrossSame /= 2;
	int e = leftSame + rightSame + leftDiff + rightDiff + acrossDiff + acrossSame;
	int eLeft = leftSame - leftDiff;
	int eRight = rightSame - rightDiff;
	vRight = newg.size() - vLeft;

	double polarity = (double) ((2 * (eLeft + eRight)) + (2 * (acrossDiff - acrossSame))) / (newg.size()); // this is cc^- from Bonchi CIKM'19; larger the better

	double harmonyLeft, harmonyRight, tension;  // larger the better
	if (vLeft > 1)
		harmonyLeft = (double) eLeft / (vLeft * (vLeft - 1) / 2);
	else if (vLeft == 1)
		harmonyLeft = 1;
	else if (vLeft == 0)
		harmonyLeft = -1;

	if (vRight > 1)
		harmonyRight = (double) eRight / (vRight * (vRight - 1) / 2);
	else if (vRight == 1)
		harmonyRight = 1;
	else if (vRight == 0)
		harmonyRight = -1;

//	printf ("all: %d\n", eLeft + eRight + acrossSame - acrossDiff);
	if (vLeft > 0 && vRight > 0)
		tension = (double) (acrossDiff - acrossSame) / (vLeft * vRight);
	else
		tension = -1;
	sc->vl = vLeft;
	sc->vr = vRight;
	sc->metrics[0] = polarity;
	sc->metrics[1] = harmonyLeft;
	sc->metrics[2] = harmonyRight;
	sc->metrics[3] = tension;
	return;
}

// smaller the better
int agreementScore (Graph& newg, Graph& newsigns, vector<vertex>& leftRight) {
	int score = 0; // number of edges that do not obey the bipartition
	for (auto u = 0; u < newg.size(); u++)
		for (auto j = 0; j < newg[u].size(); j++) {
			vertex v = newg[u][j];
			if (u < v) {
				if (newsigns[u][j] == 1 && leftRight[u] != leftRight[v])
					score++;
				if (newsigns[u][j] == -1 && leftRight[u] == leftRight[v])
					score++;
			}
		}
	return score;
}

// from "A Local-Search 2-Approximation for 2-Correlation-Clustering" by Coleman et al.
void past (Graph& newg, Graph& newsigns, vector<vertex>& leftRight) {
	int minVertex = -1; // vertex id for which the agreement score is minimum
	int globalMin = INT_MAX;
	queue<vertex> q;
	vector<vertex> best;
	for (auto u = 0; u < newg.size(); u++) {
		vector<vertex> visited (newg.size(), false);
		clearQueue(q);
		q.push (u);
		visited[u] = true;
		while (!q.empty()) {
			vertex s = q.front();
			q.pop();
			leftRight[s] = -1;
			for (auto j = 0; j < newg[s].size(); j++) {
				vertex v = newg[s][j];
				if (newsigns[s][j] == 1 && !visited[v]) {
					q.push (v);
					visited[v] = true;
				}
			}
		}
		for (auto i = 0; i < leftRight.size(); i++) {
			if (leftRight[i] == 0)
				leftRight[i] = 1;
		}
		int s = agreementScore (newg, newsigns, leftRight);
		if (s < globalMin) {
			globalMin = s;
			minVertex = u;
			best = leftRight;
			if (s == 0)
				break; // means subgraph is perfectly balanced
		}
	}
	leftRight = best;
	int s = 0;
	for (auto i = 0; i < leftRight.size(); i++)
		s += leftRight[i];
	if (s > 0) { // means Right is larger
		for (auto i = 0; i < leftRight.size(); i++)
			leftRight[i] *= -1;
	}
	return;
}

// compute the balance of the subgraph induced by vset
double computeBalance (bool metricComputation, vector<vertex>& vset, Graph& graph, Graph& signs, vector<vertex>& leftRight, subcore* sc, unordered_map<vertex, int>* lr = NULL) {

	// find the induced subgraph by vset
	unordered_map<vertex, int> subg;
	int id = 0;
	for (auto v: vset) {
		subg.emplace (v, id++);
	}

	// obtain the induced subgraph by vset
	Graph newg (vset.size());
	Graph newsigns (vset.size());
	int pe = 0, ne = 0; // # positive and negative edges
	for (int i = 0; i < vset.size(); i++) {
		int u = vset[i];
		for (int j = 0; j < graph[u].size(); j++) {
			int w = graph[u][j];
			if (subg.find (w) != subg.end()) {
				newg[i].push_back (subg[w]);
				int sw = signs[u][j];
				newsigns[i].push_back (sw);
				(sw == 1) ? pe++ : ne++;
			}
		}
	}

	pe /= 2;
	ne /= 2;
	sc->ne = ne; // number of negative edges

	if (metricComputation) { // metricComputation is to just read the polarized communities (competitors) and output the metrics in a *NUCLEI format
		for (int i = 0; i < vset.size(); i++) {
			vertex u = vset[i];
			leftRight[i] = (*lr)[u];
		}
	}
	else {
		past (newg, newsigns, leftRight);
	}

	// leftRight[i] is -1 if i is in the left; left set is larger than or equal to right set
	polarityMetrics (newg, newsigns, leftRight, sc);

	// count the balanced and unbalanced triangles
	int tc = 0, btc = 0;
	int s1, s2, s3, sum, ab;
	for (size_t i = 0; i < newg.size(); i++) {
		for (size_t j = 0; j < newg[i].size(); j++) {
			vertex a = newg[i][j];
			if (newg[i].size() <= newg[a].size()) {
				for (size_t k = j + 1; k < newg[i].size(); k++) {
					vertex b = newg[i][k];
					vertex idx = ind (b, newg[a]);
					if (idx != -1) {
						tc++;

						s1 = newsigns[i][j]; // i-a
						s2 = newsigns[i][k]; // i-b
						s3 = newsigns[a][ind (b, newg[a])]; // a-b
						sum = s1 + s2 + s3;

						if (sum == 3 || sum == -1)
							btc++;
					}
				}
			}
		}
	}

	if (tc > 0)
		return ((double) btc) / tc;
	return 0;
}

// read the left and right sets (in each line) from a file and compute metrics
bool computeMetrics (FILE* f, Graph& graph, Graph& signs, edge nEdge, FILE* fp) {

	subcore sc (-1); // K, parent and other unnecessary things are set to -1
	vector<vertex> vset;
	unordered_map<vertex, int> lr;
	bool flag = false; // set to true when the first -1 is seen
	int d;
	while (fscanf (f, "%d", &d) != EOF) {
		if (d != -1) {
			vset.push_back (d);
			if (!flag)
				lr[d] = -1;
			else
				lr[d] = 1;
		}
		else if (!flag)
			flag = true;
		else if (flag)
			break;
	}
	if (!flag) // nothing is read since EOF is reached
		return false;

	sort (vset.begin(), vset.end());

	// edge density
	edge edge_count = 0;
	if (vset.size() <= UPPERBOUND)
		for (size_t i = 0; i < vset.size(); i++)
			edge_count += commons (vset, graph[vset[i]]);

	edge_count /= 2;
	sc.nEdge = edge_count;
	sc.size = vset.size();
	if (vset.size() > 1)
		sc.ed = (double) edge_count / (sc.size * (sc.size - 1) / 2);

	// compute the balance of the subgraph: D_3 metric: #balanced triangles/# all triangles
	vector<vertex> leftRight (vset.size(), 0); // two partitions in the subgraph, left is -1; right is 1. Wlog, left is the larger one.
	sc.balance = computeBalance (true, vset, graph, signs, leftRight, &sc, &lr);

	sc.print (2, -1, fp); // output to *METRICS file

	// print the vertices in the left set
	for (size_t i = 0; i < leftRight.size(); i++)
		if (leftRight[i] == -1)
			fprintf(fp, "%d ", vset[i]);
	fprintf(fp, "-1 ");

	// print the vertices in the right set
	for (size_t i = 0; i < leftRight.size(); i++)
		if (leftRight[i] == 1)
			fprintf(fp, "%d ", vset[i]);
	fprintf(fp, "-1\n");

	return true;
}

// read the left and right sets (in each line) from a file and compute metrics
bool computeMetricsZhao (Graph& graph, Graph& signs, edge nEdge, FILE* fp, set<vertex>& subG) {

	subcore sc (-1); // K, parent and other unnecessary things are set to -1
	vector<vertex> vset;
	unordered_map<vertex, int> lr;

	for (vertex u : subG) {
		vset.push_back(u);
		lr[u] = -1;
	}

	sort (vset.begin(), vset.end());

	// edge density
	edge edge_count = 0;
	if (vset.size() <= UPPERBOUND)
		for (size_t i = 0; i < vset.size(); i++)
			edge_count += commons (vset, graph[vset[i]]);

	edge_count /= 2;
	sc.nEdge = edge_count;
	sc.size = vset.size();
	if (vset.size() > 1)
		sc.ed = (double) edge_count / (sc.size * (sc.size - 1) / 2);

	// compute the balance of the subgraph: D_3 metric: #balanced triangles/# all triangles
	vector<vertex> leftRight (vset.size(), 0); // two partitions in the subgraph, left is -1; right is 1. Wlog, left is the larger one.
	sc.balance = computeBalance (false, vset, graph, signs, leftRight, &sc, &lr);

	sc.print (2, -1, fp); // output to *METRICS file

	// print the vertices in the left set
	for (size_t i = 0; i < leftRight.size(); i++)
		if (leftRight[i] == -1)
			fprintf(fp, "%d ", vset[i]);
	fprintf(fp, "-1 ");

	// print the vertices in the right set
	for (size_t i = 0; i < leftRight.size(); i++)
		if (leftRight[i] == 1)
			fprintf(fp, "%d ", vset[i]);
	fprintf(fp, "-1\n");

	return true;
}


// find the r-cliques whose component is index, append the vertices in those cliques to the vertices of its all children, sort, compute its density
int reportSubgraph (vector<vertex>& vset, int variant, vertex index, unordered_map<vertex, vertex>& orderInFile, vector<vertex>& component, helpers& ax, vector<subcore>& skeleton, Graph& graph, Graph& signs, edge nEdge, FILE* fp, FILE* gp) {

	bool pass = true;
	if (variant > 0) {
		if (skeleton[index].parent == -1) {

			skeleton[index].size = graph.size();
			skeleton[index].nEdge = nEdge;
			skeleton[index].ed = 0;
			skeleton[index].print (2, index, fp);
			fprintf(fp, "-1\n");
			return -1;
		}

		vertex i;

		if (variant == 12 || variant == 13 || variant == 14) {
			for (i = 0; i < component.size(); i++) {
				if (component[i] == index)
					vset.push_back (i);
			}
		}
		else if (variant == 230 || variant == 23 || variant == 24) {
			for (vertex i = 0; i < component.size(); i++)
				if (component[i] == index) {
					vset.push_back (get<0>((*ax.el)[i]));
					vset.push_back (get<1>((*ax.el)[i]));
				}
		}
		else if (variant == 34) {
			for (vertex i = 0; i < component.size(); i++)
				if (component[i] == index) {
					vset.push_back (get<0>((*ax.tris)[i]));
					vset.push_back (get<1>((*ax.tris)[i]));
					vset.push_back (get<2>((*ax.tris)[i]));
				}
		}

		pass = hashUniquify (vset);
		if (!pass) {
			dummyLine (&skeleton[index], fp, index);
			return -1;
		}

		pass = true;
		pass = pullChildrenSets (variant, fp, index, orderInFile, vset, skeleton);
		if (!pass) {
			dummyLine (&skeleton[index], fp, index);
			return -1;
		}

		pass = hashUniquify (vset);
		if (!pass) {
			dummyLine (&skeleton[index], fp, index);
			return -1;
		}

	}

	unordered_map<vertex, int> lr;
	if (variant == -19) {
		int leftsize = vset[0];
		vset.erase (vset.begin());
		for (int i = 0; i < leftsize; i++)
			lr[vset[i]] = -1;
		for (int i = leftsize; i < vset.size(); i++)
			lr[vset[i]] = 1;
		sort (vset.begin(), vset.end());
	}

	int r = stats (vset, variant, index, skeleton, graph, signs, nEdge, fp, gp, lr);

	return r;
}

int stats (vector<vertex>& vset, int variant, vertex index, vector<subcore>& skeleton, Graph& graph, Graph& signs, edge nEdge, FILE* fp, FILE* gp, unordered_map<vertex, int>& lr) {
	// edge density
	edge edge_count = 0;
	if (vset.size() <= UPPERBOUND)
		for (size_t i = 0; i < vset.size(); i++)
			edge_count += commons (vset, graph[vset[i]]);

	edge_count /= 2;
	skeleton[index].nEdge = edge_count;
	skeleton[index].size = vset.size();
	if (vset.size() > 1)
		skeleton[index].ed = (double) edge_count / (skeleton[index].size * (skeleton[index].size - 1) / 2);


	// compute the balance of the subgraph: D_3 metric: #balanced triangles/# all triangles
	vector<vertex> leftRight (vset.size(), 0); // two partitions in the subgraph, left is -1; right is 1. Wlog, left is the larger one.
	if (variant == -19)
		skeleton[index].balance = computeBalance (true, vset, graph, signs, leftRight, &(skeleton[index]), &lr);
	else
		skeleton[index].balance = computeBalance (false, vset, graph, signs, leftRight, &(skeleton[index]));

	bool highlight = (skeleton[index].size >= LOWERBOUND); // we didn't know the size before reportSubgraph call; that's why this is here

	if (highlight) {
		skeleton[index].print (1, index, gp); // output to *Hierarchy file
		skeleton[index].print (0); // output to command line
	}

	skeleton[index].print (2, index, fp); // output to *NUCLEI file

	// print the vertices in the left set
	for (size_t i = 0; i < leftRight.size(); i++)
		if (leftRight[i] == -1)
			fprintf(fp, "%d ", vset[i]);
	fprintf(fp, "-1 ");

	// print the vertices in the right set
	for (size_t i = 0; i < leftRight.size(); i++)
		if (leftRight[i] == 1)
			fprintf(fp, "%d ", vset[i]);
	fprintf(fp, "-1\n");


	if (skeleton[index].size >= 10)
		return 1;
	else
		return 0;
}

void bfsHierarchy (vector<subcore>& skeleton, stack<vertex>& scs) {
	rearrange (skeleton);
	queue<vertex> bfsorder; // we are doing bfs on the hierarchy tree and push the dequeued nodes to the stack
	bfsorder.push(skeleton.size() - 1);
	while (!bfsorder.empty()) {
		vertex s = bfsorder.front();
		bfsorder.pop();
		if (skeleton[s].children.empty()) {
			skeleton[s].levelFromLeaf = 0;
			int p = s;
			int l = 1;
			while (skeleton[p].parent != -1) {
				int pp = skeleton[p].parent;
				if (skeleton[pp].visible) {
					int cl = skeleton[pp].levelFromLeaf;
					if (cl > l || cl == -1)
						skeleton[pp].levelFromLeaf = l++;
				}
				p = pp;
			}
		}
		scs.push (s);
		for (vertex r : skeleton[s].children)
			bfsorder.push (r);
	}
}

inline void findRepresentative (vertex* child, vector<subcore>& skeleton) {
	vertex u = *child;
	if (skeleton[u].parent != -1) {
		vertex pr = skeleton[u].parent;
		while (skeleton[u].K == skeleton[pr].K) {
			u = pr;
			if (skeleton[u].parent != -1)
				pr = skeleton[u].parent;
			else
				break;
		}
	}
	*child = u;
}

int presentNuclei (int variant, vector<subcore>& skeleton, vector<vertex>& component, Graph& graph, Graph& signs, edge* nEdge, helpers& ax, string vfile, FILE* gp) {

	// assign unassigned items to top subcore
	for (vertex i = 0; i < component.size(); i++)
		if (component[i] == -1)
			component[i] = skeleton.size() - 1;

	// match each component with its representative
	unordered_map<vertex, vertex> replace;
	for (vertex i = 0; i < skeleton.size(); i++) {
		vertex sc = i;
		vertex original = sc;
		findRepresentative (&sc, skeleton);
		if (original != sc)
			skeleton[original].visible = false;
		replace[original] = sc;
	}

	// each component takes its representative's component number
	for (vertex i = 0; i < component.size(); i++)
		if (replace.find (component[i]) != replace.end())
			component[i] = replace[component[i]];

	stack<vertex> subcoreStack;
	bfsHierarchy (skeleton, subcoreStack);

	string nFile = vfile + "_NUCLEI";
	FILE* fp = fopen (nFile.c_str(), "w+");
	vector<subcore> backup (skeleton);

	vector<vertex> allReportedNodes;

	unordered_map<vertex, vertex> orderInFile; // key is the skeleton index, value is the order
	vertex o = 0; // order of subcores in file
	int numOfSubgraphs = 0; // number of subgraphs with at least 10 nodes
	printf("TrussN\t|V|\t|E|\t|E-|\ted\tbalance\tleaf?\tpol\tLsize\tRsize\thLeft\thRight\ttension\n");
	while (!subcoreStack.empty()) {
		vertex i = subcoreStack.top();
		subcoreStack.pop();
		if (backup[i].visible && backup[i].children.empty() && backup[i].levelFromLeaf <= DEPTH) {
			orderInFile[i] = o++;
			vector<vertex> reportedNodes;
			int ret = reportSubgraph (reportedNodes, variant, i, orderInFile, component, ax, skeleton, graph, signs, *nEdge, fp, gp);
			if (variant == 230) {
				if (ret == 1)
					numOfSubgraphs++;
				allReportedNodes.insert (allReportedNodes.end(), reportedNodes.begin(), reportedNodes.end());
			}
			removeChild (i, backup);
		}

	}
	fclose (fp);

	// remove the found subgraphs from the graph
	if (variant == 230) {
		printf ("allReportedNodes.size: %d\n", allReportedNodes.size());

		for (auto u: allReportedNodes) {
			for (int i = 0; i < graph[u].size(); i++) {
				vertex v = graph[u][i];
				for (int j = 0; j < graph[v].size(); j++) {
					if (graph[v][j] == u) {
						graph[v].erase (graph[v].begin() + j);
						signs[v].erase (signs[v].begin() + j);
						break;
					}
				}
			}
		}

		for (auto u: allReportedNodes) {
			graph[u].clear();
			signs[u].clear();
		}

		int s = 0;
		for (auto g: graph)
			s += g.size();

		*nEdge = s/2;
		printf ("new nEdge: %d\n", *nEdge);

	}
	return numOfSubgraphs;
}
