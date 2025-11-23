"""
After running the signed nucleus code, run the following to get several additional metrics for the corresponding subgraphs:

python3 metrics.py [input graph file] [_23_3_NUCLEI file] [output text file] [Optional: Metrics (1,2,3,...)]

- Input Graph Format
    # node  neighbor  sign
    0   1   1
    1   2   -1

Example Command:
python3 metrics.py example.txt example.txt_23_3_NUCLEI example_metrics.txt 2,3,7

Default metrics: 1,2,6,7

Metrics:
1 - Frustration index
2 - Normalized frustration index
3 - Algebraic conflict
4 - Normalized algebraic conflict
5 - Walk-based measure of balance
6 - Cohesiveness
7 - Divisiveness

NOTE: You need to follow the Gurobi instructions here: https://github.com/saref/frustration-index-XOR.
However, Anaconda is not required, so you may skip those steps if desired.
Make sure you have all other required modules installed (see requirements.txt). This is tested on Python3.8.6.
"""

from copy import deepcopy
import math
import sys
import numpy as np
import concurrent.futures as futures

# For frustration index calculation
import networkx as nx
import multiprocessing
from gurobipy import *

# Returns frustration index (Adapted from https://github.com/saref/frustration-index-XOR)
def get_frustration_index(G, S, graph):
    signedMatrix = nx.to_numpy_matrix(graph)

    weighted_edges = nx.get_edge_attributes(graph, "weight")
    sorted_weighted_edges = {}
    for (u, v) in weighted_edges:
        if u < v:
            sorted_weighted_edges[(u, v)] = weighted_edges[(u, v)]
        if u > v:
            sorted_weighted_edges[(v, u)] = weighted_edges[(u, v)]

    # lazyParam=int(input("What is the lazy parameter for unbalanced triangle lazy cuts? (0/1/2/3)"))
    # See "lazy" as a tunable parameter in linear constraint attributes in Gurobi optimizer reference manual below:
    # https://www.gurobi.com/documentation/8.1/refman/lazy.html
    lazyParam = int(1)

    # methodParam=int(input("What method do you want to use?"))
    # (-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent)
    methodParam = int(-1)

    # speedupParam=int(input("Do you want to use the speedups? (0=No, 1=Yes)"))
    speedupParam = int(1)

    order = len(signedMatrix)
    NumberOfNegative = ((-1 == signedMatrix)).sum() / 2
    size = int(np.count_nonzero(signedMatrix) / 2)

    neighbors = {}
    Degree = []
    for u in sorted(graph.nodes()):
        neighbors[u] = list(graph[u])
        Degree.append(len(neighbors[u]))
    unsignedDegree = Degree

    # Finding the node with the highest unsigned degree
    maximum_degree = max(unsignedDegree)
    [node_to_fix] = [
        ([i for i, j in enumerate(unsignedDegree) if j == maximum_degree]).pop()
    ]

    # Model parameters
    model = Model("XOR model for computing frustration index")

    # The method you have selected above will be used:
    # (-1=automatic, 0=primal simplex, 1=dual simplex, 2=barrier, 3=concurrent, 4=deterministic concurrent)
    model.setParam(GRB.param.Method, methodParam)

    # What is the time limit in second?
    model.setParam("TimeLimit", 10 * 3600)

    # Do you want details of branching to be reported? (0=No, 1=Yes)
    model.setParam(GRB.param.OutputFlag, 1)

    # How many threads to be used for exploring the feasible space in parallel?
    model.setParam(GRB.Param.Threads, multiprocessing.cpu_count())

    # This chunk of code lists the graph triangles
    GraphTriangles = []
    for n1 in sorted((graph).nodes()):
        neighbors1 = set((graph)[n1])
        for n2 in filter(lambda x: x > n1, neighbors1):
            neighbors2 = set((graph)[n2])
            common = neighbors1 & neighbors2
            for n3 in filter(lambda x: x > n2, common):
                GraphTriangles.append([n1, n2, n3])

    # This chunk of code lists the balanced and unbalanced triangles
    w = nx.get_edge_attributes(graph, "weight")
    unbalanced_triangles = []
    balanced_triangles = []
    for triad in GraphTriangles:
        if (
            sorted_weighted_edges[(triad[0], triad[1])]
            * sorted_weighted_edges[(triad[0], triad[2])]
            * sorted_weighted_edges[(triad[1], triad[2])]
            == -1
        ):
            unbalanced_triangles.append(triad)
        elif (
            sorted_weighted_edges[(triad[0], triad[1])]
            * sorted_weighted_edges[(triad[0], triad[2])]
            * sorted_weighted_edges[(triad[1], triad[2])]
            == 1
        ):
            balanced_triangles.append(triad)

    # Create decision variables and update model to integrate new variables
    x = []
    for i in range(0, order):
        x.append(model.addVar(vtype=GRB.BINARY, name="x" + str(i)))  # arguments by name
    model.update()

    f = {}
    for (i, j) in sorted_weighted_edges:
        f[(i, j)] = model.addVar(
            lb=0.0, ub=1, vtype=GRB.BINARY, name="f" + str(i) + "," + str(j)
        )
    model.update()

    # Set the objective function
    OFV = 0
    for (i, j) in sorted_weighted_edges:
        OFV = OFV + f[(i, j)]
    model.setObjective(OFV, GRB.MINIMIZE)

    # Add constraints to the model and update model to integrate new constraints

    ## ADD CORE CONSTRAINTS ##

    for (i, j) in sorted_weighted_edges:
        model.addConstr(
            f[(i, j)]
            >= x[i]
            - (sorted_weighted_edges[(i, j)]) * x[j]
            - (1 - sorted_weighted_edges[(i, j)]) / 2,
            "1st Edge" + "," + str(i) + "," + str(j),
        )
        model.addConstr(
            f[(i, j)]
            >= -x[i]
            + (sorted_weighted_edges[(i, j)]) * x[j]
            + (1 - sorted_weighted_edges[(i, j)]) / 2,
            "2nd Edge" + "," + str(i) + "," + str(j),
        )
    model.update()

    ## ADD ADDITIONAL CONSTRAINTS (speed-ups) ##

    if speedupParam == 1:

        # Triangle valid inequalities
        # lazyTriangles=int(input("What triangles do you want to implement the inequalities on? (1:all, 2:unbalanced)"))
        lazyTriangles = int(2)

        if lazyTriangles == 1:
            triangleInequalityCount = len(GraphTriangles) * 4
            for triangle in GraphTriangles:
                [i, j, k] = triangle
                model.addConstr(
                    (
                        (f[(i, j)] + (sorted_weighted_edges[(i, j)] - 1) / 2)
                        / sorted_weighted_edges[(i, j)]
                    )
                    + (
                        (f[(i, k)] + (sorted_weighted_edges[(i, k)] - 1) / 2)
                        / sorted_weighted_edges[(i, k)]
                    )
                    >= (
                        (f[(j, k)] + (sorted_weighted_edges[(j, k)] - 1) / 2)
                        / sorted_weighted_edges[(j, k)]
                    ),
                    "triangle1" + "," + str(i) + "," + str(j) + "," + str(k),
                )
                model.addConstr(
                    (
                        (f[(i, j)] + (sorted_weighted_edges[(i, j)] - 1) / 2)
                        / sorted_weighted_edges[(i, j)]
                    )
                    + (
                        (f[(j, k)] + (sorted_weighted_edges[(j, k)] - 1) / 2)
                        / sorted_weighted_edges[(j, k)]
                    )
                    >= (
                        (f[(i, k)] + (sorted_weighted_edges[(i, k)] - 1) / 2)
                        / sorted_weighted_edges[(i, k)]
                    ),
                    "triangle2" + "," + str(i) + "," + str(j) + "," + str(k),
                )
                model.addConstr(
                    (
                        (f[(j, k)] + (sorted_weighted_edges[(j, k)] - 1) / 2)
                        / sorted_weighted_edges[(j, k)]
                    )
                    + (
                        (f[(i, k)] + (sorted_weighted_edges[(i, k)] - 1) / 2)
                        / sorted_weighted_edges[(i, k)]
                    )
                    >= (
                        (f[(i, j)] + (sorted_weighted_edges[(i, j)] - 1) / 2)
                        / sorted_weighted_edges[(i, j)]
                    ),
                    "triangle3" + "," + str(i) + "," + str(j) + "," + str(k),
                )
                model.addConstr(
                    (
                        (f[(i, j)] + (sorted_weighted_edges[(i, j)] - 1) / 2)
                        / sorted_weighted_edges[(i, j)]
                    )
                    + (
                        (f[(i, k)] + (sorted_weighted_edges[(i, k)] - 1) / 2)
                        / sorted_weighted_edges[(i, k)]
                    )
                    + (
                        (f[(j, k)] + (sorted_weighted_edges[(j, k)] - 1) / 2)
                        / sorted_weighted_edges[(j, k)]
                    )
                    <= 2,
                    "triangle4" + "," + str(i) + "," + str(j) + "," + str(k),
                )

        elif lazyTriangles == 2:
            triangleInequalityCount = len(unbalanced_triangles)
            for triangle in unbalanced_triangles:
                [i, j, k] = triangle
                model.addConstr(
                    f[(i, j)] + f[(i, k)] + f[(j, k)] >= 1,
                    "UnbalancedTriangle" + "," + str(i) + "," + str(j) + "," + str(k),
                )

        elif lazyTriangles == 3:
            triangleInequalityCount = 0

        model.update()
        model.setAttr(
            "Lazy",
            model.getConstrs()[2 * size : 2 * size + triangleInequalityCount],
            [lazyParam] * triangleInequalityCount,
        )
        model.update()

        # Colour the node with the highest degree as 1
        model.addConstr(x[node_to_fix] == 1, "1stnodecolour")
        model.update()

        # branching priority is based on unsigned degree
        model.setAttr("BranchPriority", model.getVars()[:order], unsignedDegree)
        model.update()

    # Calculate Frustration Index
    model.optimize()

    # Save optimal objective function value
    obj = model.getObjective()
    objectivevalue = obj.getValue()

    return objectivevalue

# Returns metrics given G (adjacency list), S (edge signs), and metrics (dictionary)
def get_metrics(G, S, metrics):

    graph = nx.Graph()
    for u in G:
        for v in G[u]:
            e = (min(int(u), int(v)), max(int(u), int(v)))
            graph.add_edge(e[0], e[1], weight=S[e])

    if "frustration index" in metrics or "normalized frustration index" in metrics:
        # Prevents printing to terminal
        old_stdout = sys.stdout
        sys.stdout = open(os.devnull, "w")

        F = get_frustration_index(G, S, graph)

        # Allows printing to terminal
        sys.stdout = old_stdout

        if "frustration index" in metrics:
            metrics["frustration index"] = int(F)

    numNodes = len(G)

    # Signed Adjacency Matrix (V x V)
    A = np.zeros((numNodes, numNodes))

    # Diagonal Matrix (V x V)
    D = np.zeros((numNodes, numNodes))

    dmax = 0
    m = 0
    for u, neighbors in G.items():
        deg_u = len(neighbors)
        D[u, u] = deg_u
        for v in neighbors:
            e = (min(u, v), max(u, v))
            A[u, v] = S[e]
            d = (deg_u + len(G[v])) / 2
            if d > dmax:
                dmax = d
            m += 1
    m = m / 2

    if "normalized frustration index" in metrics:
        # Normalized Frustration Index
        if m == 0:
            normalizedF = 1
        else:
            normalizedF = 1 - (int(F) / (m / 2))
        metrics["normalized frustration index"] = round(normalizedF, 2)

    if "algebraic conflict" in metrics or "normalized algebraic conflict" in metrics:
        # Signed Laplacian Matrix
        L = np.subtract(D, A)
        # Algebraic Conflict
        eigvals = np.linalg.eigvalsh(L)
        # Numpy eigenvalues are not 100% accurate (ex. values close to 0 should be 0)
        algebraicConflict = round(eigvals[0], 12)
        # Converts 0.0 and -0.0 to 0
        if algebraicConflict == 0:
            algebraicConflict = 0

        if "algebraic conflict" in metrics:
            metrics["algebraic conflict"] = round(algebraicConflict, 2)

        if "normalized algebraic conflict" in metrics:
            # Normalized Algebraic Conflict
            if dmax == 1:
                normalizedAlgebraicConflict = 1
            else:
                normalizedAlgebraicConflict = 1 - (algebraicConflict / (dmax - 1))
            metrics["normalized algebraic conflict"] = round(
                normalizedAlgebraicConflict, 2
            )

    if "walk-based measure of balance" in metrics:
        # Walk-based measure of balance
        eigvalsA = np.linalg.eigvalsh(A)
        absA = np.absolute(A)
        eigvalsabsA = np.linalg.eigvalsh(absA)
        W_num = 0
        W_den = 0
        for val in eigvalsA:
            W_num += math.exp(val)
        for val in eigvalsabsA:
            W_den += math.exp(val)
        W = W_num / W_den
        metrics["walk-based measure of balance"] = round(W, 2)

    return metrics

metrics = {}
if len(sys.argv) > 4:
    inputMetrics = sys.argv[4].split(",")
    if "1" in inputMetrics:
        metrics["frustration index"] = None
    if "2" in inputMetrics:
        metrics["normalized frustration index"] = None
    if "3" in inputMetrics:
        metrics["algebraic conflict"] = None
    if "4" in inputMetrics:
        metrics["normalized algebraic conflict"] = None
    if "5" in inputMetrics:
        metrics["walk-based measure of balance"] = None
    if "6" in inputMetrics:
        metrics["cohesiveness"] = None
    if "7" in inputMetrics:
        metrics["divisiveness"] = None
else:
    # Default Metrics
    metrics = {
        "frustration index": None,
        "normalized frustration index": None,
        "cohesiveness": None,
        "divisiveness": None,
    }

# Adjacency list
G = {}

# Edge signs
S = {}

numNodes = None

# Reads dataset
f = open(sys.argv[1], "r")

for line in f:
    line = line.strip()
    if "#" in line:
        values = line.split()
        numNodes = int(values[1])

        for u in range(numNodes):
            G[u] = []

    delim = ""
    if line[:1].isdigit():
        if delim == "":
            for char in line:
                if not char.isdigit():
                    delim = char
                    break
        values = line.split(delim)
        values = [i for i in values if i]
        # No self loops
        if values[0] != values[1]:
            # Populate edge signs
            e = (
                min(int(values[0]), int(values[1])),
                max(int(values[0]), int(values[1])),
            )
            if int(values[2]) == 1:
                S[e] = 1
            else:
                S[e] = -1
            values = list(map(int, values[:2]))
            # Populate adjacency list
            if values[1] not in G[values[0]]:
                G[values[0]].append(values[1])
            if values[0] not in G[values[1]]:
                G[values[1]].append(values[0])
f.close()

subgraphs = {}
signs = {}
otherInfo = {}

def obtain_rank(data):
    sort_data = sorted(range(len(data)), key=lambda i: data[i])
    result = [0] * len(data)
    for i, idx in enumerate(sort_data, 1):
        result[idx] = i - 1
    return result

cohesiveness = {}
divisiveness = {}

# Reads nuclei file
count = 0
prev_lines = []
f = open(sys.argv[2], "r")
for line in f:
    line = line.strip()
    if "\t" in line:
        line = line.split("\t")
        data = line[0].split()
        subG_id = int(data[0])
        if subG_id == -1:
            subG_id = count
            count += 1
        
        left = []
        right = []
        is_right = False

        N = [False] * numNodes
        nodes = []
        edges = []
        subG = {}
        for i in line[1].split()[:-1]:
            u = int(i)
            if u == -1:
                is_right = True
            else:
                if is_right:
                    right.append(u)
                else:
                    left.append(u)
                if u not in subG:
                    nodes.append(u)
                    N[u] = True
                    subG[u] = []
                    for v in G[u]:
                        if N[v]:
                            subG[u].append(v)
                            subG[v].append(u)
                            edges.append((min(u, v), max(u, v)))
        subG_copy = deepcopy(subG)
        for key, val in subG_copy.items():
            if len(val) == 0:
                del subG[key]
                nodes.remove(key)
        if nodes not in prev_lines:
            prev_lines.append(nodes)
        else:
            continue
        if len(nodes) >= 10 and len(nodes) <= 1000:
            internal_edges = 0
            pos_internal_edges = 0
            for v1 in left:
                for v2 in left:
                    e = (
                        min(v1, v2),
                        max(v1, v2),
                    )
                    if e in S:
                        if S[e] == 1:
                            pos_internal_edges += 1
                        internal_edges += 1
            for v1 in right:
                for v2 in right:
                    e = (
                        min(v1, v2),
                        max(v1, v2),
                    )
                    if e in S:
                        if S[e] == 1:
                            pos_internal_edges += 1
                        internal_edges += 1
            
            # partitions are unconnected within
            if internal_edges == 0:
                continue

            if "cohesiveness" in metrics:
                cohesiveness[subG_id] = round(pos_internal_edges / internal_edges, 2)
            
            if "divisiveness" in metrics:
                external_edges = 0
                neg_external_edges = 0
                for v1 in left:
                    for v2 in right:
                        e = (
                            min(v1, v2),
                            max(v1, v2),
                        )
                        if e in S:
                            if S[e] == -1:
                                neg_external_edges += 1
                            external_edges += 1
                if external_edges == 0:
                    divisiveness[subG_id] = 1
                else:
                    divisiveness[subG_id] = round(neg_external_edges / external_edges, 2)
                
            # Changes node labels to match frustration index parameters
            ranking = obtain_rank(nodes)
            node_to_rank = {}
            for i in range(len(nodes)):
                node_to_rank[nodes[i]] = ranking[i]
            new_subG = {}
            for u, neighbors in subG.items():
                new_subG[node_to_rank[u]] = []
                for v in subG[u]:
                    new_subG[node_to_rank[u]].append(node_to_rank[v])
            newS = {}
            for e in edges:
                u, v = e
                newE = (node_to_rank[u], node_to_rank[v])
                newS[newE] = S[e]
            subgraphs[subG_id] = new_subG
            signs[subG_id] = newS

            numV = data[2]
            numE = data[3]
            numNegE = data[4]
            density = data[5]
            balance = data[6]
            v_l = data[10]
            v_r = data[11]
            harmonyLeft = data[12]
            harmonyRight = data[13]
            tension = data[14]
            polarity = data[9]
            otherInfo[subG_id] = (numV, numE, numNegE, density, balance, v_l, v_r, harmonyLeft, harmonyRight, tension, polarity)
f.close()

if len(subgraphs) == 0:
    print('NO SUBGRAPHS ERROR')
    print(sys.argv[2])
    sys.exit()

del G
del S

# Parallel Processing
executor = futures.ProcessPoolExecutor()
futures = {executor.submit(get_metrics, subG, signs[subG_id], metrics): subG_id for subG_id, subG in subgraphs.items()}
executor.shutdown()

f = open(sys.argv[3], "w")

for future in futures:
    metrics = future.result()
    subG_id = futures[future]

    numV, numE, numNegE, density, balance, v_l, v_r, harmonyLeft, harmonyRight, tension, polarity = otherInfo[subG_id]

    f.write("id: " + str(subG_id) + " |V|: " + str(numV) + " |E|: " + str(numE) + " |E-|: " + str(numNegE) + " density: " + str(density) + " |V_l|: " + str(v_l) + " |V_r|: " + str(v_r) + " harmony left: " + str(harmonyLeft) + " harmony right: " + str(harmonyRight) + " tension: " + str(tension) + " polarity: " + str(polarity) + " relative 3-balance: " + str(balance))

    if "cohesiveness" in metrics:
        metrics["cohesiveness"] = cohesiveness[subG_id]
    if "divisiveness" in metrics:
        metrics["divisiveness"] = divisiveness[subG_id]
        
    space = False
    for metric, val in metrics.items():
        f.write(" " + metric + ": " + str(val))
    f.write("\n")

f.close()
