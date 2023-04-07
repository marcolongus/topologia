# Cargamos en matrices los datos para construir el grafo recursivamente.
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from networkx.drawing.nx_pydot import write_dot


# =====================================================================================================================#
# Data structure of topologia.txt an how to handle it:

# The first column of topologia.txt correspond to an infected particle. The non-childs agents are written first.
# The following columns are the childs of the first column particle.
# In order to build the epidemic graph we need to get access to every agent and their correspondent childs.
# The np.array particle_index[] get wich particle is in one of the n rows of the archive.
# With particle_inverse[particle] we get acces to the posistion of the "particle" in particle_index in order to get their childs.
# So in particle_index[particle_inverse[particle]] we get acces to the row of a "particle" and their childs.
# =====================================================================================================================#

# Deep First Search, tree building, breadth first search, hiarchy_walk
def particle_tree(n, truncate=False):
    particle = particle_inverse[n]
    if particle > N:
        return [-1]

    tree = np.loadtxt(data, skiprows=particle + shift[n_simul], max_rows=1, dtype=int)
    if tree.size == 1:
        return []
    if truncate and (tree.size >= 2):
        tree = tree[1 : tree.size]

    return tree.tolist()


def dfs(visited, node, graph):
    graph[node] = particle_tree(node, True)
    #print(graph[node], node)
    if node not in visited:
        visited.add(node)
        for neighbour in graph.get(node):
            prev_node = neighbour
            dfs(visited, neighbour, graph)



def hiarchy_walk(graph, nodes):
    branch_out_degree = 0
    descendientes = set()
    for element in nodes:
        branch_out_degree += graph.out_degree(element)
        for next_branch_node in list(graph.successors(element)):
            descendientes.add(next_branch_node)
    return [branch_out_degree, descendientes]


# ============================================================================#
# Data, Source node, simulation paramaters:
# ============================================================================#
data = "main/data/topologia.txt"  # The data of the n-simulation where continuously written on this archive.
N = 1000
start_node = [i for i in range(10)]
simulation = np.loadtxt("main/data/evolution.txt", usecols=0, ndmin=1, dtype=int)  # size of every epidemic in data.
print("simulation: ", simulation)
n_simulation = 1

generation_ditribution = []
total_degree = []
deepest_generation = 0


for n_simul in range(n_simulation):
    # shift vector tells where to search for the n-epidemic in data.
    shift = simulation.copy()

    for i in range(1, shift.size):
        shift[i] = simulation[0:i].sum()
    shift[0] = 0

    args = {
        "usecols": 0,
        "max_rows": simulation[n_simul],
        "skiprows": shift[n_simul],
        "ndmin": 1,
        "dtype": int,
    }

    particle_index = np.loadtxt(data, **args)
    particle_inverse = np.full(N, -1, dtype=int)

    for i in range(N):
        if i < particle_index.size:
            particle_inverse[particle_index[i]] = i
    np.savetxt("main/data/topologia_inverso.txt", particle_inverse, fmt="%1.f")
    
    # ============================================================================#
    # Nodes analysis
    # ============================================================================#
    for node in start_node:
        graph = {}
        max_degree = 0  # Guardamos nodo con más conexiones.
        visited = set()
        dfs(visited, node, graph) # visualization
        print()

        # ============================================================================#
        # Define nx.graph from dfs
        # ============================================================================#
        DG = nx.DiGraph(graph)  # directed graph
        G = nx.Graph(graph)     # Undirected

        # Hiarchy_walk driver code and graph analysis:
        walk_set = set()
        walk_set.add(node)

        # Np array for childs per generation [[childs, generation]]:
        p = 0
        p_distribution = np.array([hiarchy_walk(DG, walk_set)[0], p])
        print("Start Node:", node)
        if hiarchy_walk(DG, walk_set)[0] == 0: 
            generation_ditribution.append([0])
            continue

        while hiarchy_walk(DG, walk_set)[0] > 0:
            print(hiarchy_walk(DG, walk_set), p)
            p += 1
            walk_set = hiarchy_walk(DG, walk_set)[1]
            p_distribution = np.vstack((
                p_distribution, 
                np.array([hiarchy_walk(DG, walk_set)[0], p])
                )
            )

        print("\n [Childs, Generatoin]:\n ", p_distribution)
        # Degree and degree distribution of the tree:
        try:
            childs     = p_distribution[:, 1]
            generation = p_distribution[:, 0]
        except Exception as error:
            raise(error)

        try:
            deep_x = childs.max() + 1
            if deep_x > deepest_generation:
                deepest_geneneration = deep_x

            generation_ditribution.append(generation.tolist())
            degree_array = np.array(G.degree())[:, 1]
            total_degree.append(degree_array)

            max_sim = degree_array.max()

            if max_sim > max_degree:
                max_degree = max_sim
        
        except Exception as error:
            raise error


        # ============================================================================#
        # GRAFICOS
        # ============================================================================#
        Grafico = True
        if Grafico:
            #============================================#
            # GRAFICO DE GENERACIONES
            #============================================#
            plt.title("Infected per generations")
            plt.xlabel("Generation")
            plt.ylabel("Generation childs")
            plt.xlim(0, deep_x)
            plt.xticks([i * 2 for i in range(deep_x // 2 + 1)])
            plt.ylim(0, generation.max() + 5)
            plt.plot(childs, generation)
            plt.savefig("figuras/epidemic%i.png" % n_simul)
            plt.show()
            plt.cla()

            #============================================#
            # GRAFICO DE DEGREE DISTRIBUCION
            #============================================#
            bin_limit = max_degree
            bins = np.linspace(1, bin_limit, bin_limit)
            plt.xticks([i * 2 for i in range(bin_limit // 2 + 1)])
            plt.hist(degree_array, bins, label="Degree Dist.")
            plt.legend()
            plt.savefig("figuras/degree_histogram%i.png" % n_simul)
            plt.show()
            plt.cla()
            
            #============================================#
            # GRAFICO DEL GRAFO
            #============================================#
            info = False
            if info and n_simul > -1:
                # print("nodos:",DG.nodes)
                # print("vertices:",DG.edges),print()
                for element in DG.nodes:
                    pass
                    # print(element, DG.degree(element),np.array(list((DG.successors(element))),dtype=int), DG.out_degree(element))
                # print()
                # print(nx.info(DG))
                pos = nx.kamada_kawai_layout(G, scale=1000)
                # pos = nx.spring_layout(G,iterations=1000)
                args_g = {
                    "with_labels": False,
                    "pos": pos,
                    "node_size": 15,
                    "node_color": "blue"
                    # "node_color"  : range(simulation[n_simul]),
                    # "cmap"        : None
                }

                nx.draw(DG, **args_g)
                nx.draw_networkx_nodes(
                    G, pos, nodelist=[node], node_color="r", node_size=25
                )
                plt.savefig("figuras/grafo%i.png" % n_simul)
                plt.show()
                plt.cla()

        # ============================================================================#
        # ============================================================================#


# ============================================================================#
# ESCRITURA DE DATOS Y ORGANIZACIÓN DE ARCHIVOS
# ============================================================================#
G_l = []
for gen in generation_ditribution:
    print(gen)

for i, element in enumerate(generation_ditribution):
    G_l.append(len(generation_ditribution[i]))

G_l.sort()
print("Profundidades:", G_l)

with open("estadistica.txt", "w") as f:
    for i in range(n_simulation, -1, -1):
        for gen in generation_ditribution:
            #local = np.array(generation_ditribution[i])
            local = np.array(gen)
            np.savetxt(f, local.reshape(1, local.size), fmt="%3i")


# Plot de promedio de histogramas.
TotalH = []
for element in total_degree:
    for conections in element:
        TotalH.append(conections)

bin_limit = max_degree
bins = np.linspace(1, bin_limit, bin_limit)

plt.xlabel("Degree")
plt.ylabel("Simulations")
plt.xticks([i * 5 for i in range(max_degree // 5 + 1)])
plt.hist(TotalH, bins, alpha=0.5, color="darkgreen")
plt.savefig("figuras/Histograma_promedio.png")
plt.show()
plt.clf()

print()
print("Deepest Generation:", deepest_generation)


def plot_loghist(x, bins):
    hist, bins = np.histogram(x, bins=bins)
    #logbins = np.logspace(np.log10(bins[0]), np.log10(bins[-1]), len(bins))
    plt.hist(x, bins=bins)
    #plt.loglog()
    #plt.xscale("log")


plot_loghist(TotalH, bins)
plt.show()
