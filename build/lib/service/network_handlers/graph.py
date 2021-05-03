def is_in_tuple_list(tl, val):
    res = False
    for (x, y) in tl:
        if val == x:    return True
    return res


class Graph:

    def __init__(self, g={}):
        self.graph = g

    def delete_node(self,node):
        predecessors = self.get_predecessors(node)
        for pre in predecessors:
            self.graph[pre].remove(node)
        del self.graph[node]

    def print_graph(self):
        for v in self.graph.keys():
            print(v, "->", self.graph[v])

    def get_nodes(self):
        return list(self.graph.keys())

    def get_edges(self):
        edges = []
        for v in self.graph.keys():
            for d in self.graph[v]:
                edges.append((v, d))
        return edges

    def add_vertex(self, v):
        if v not in self.graph:
            self.graph[v] = []

    def add_edge(self, o, d):
        if o not in self.graph:
            self.add_vertex(o)

        if d not in self.graph:
            self.add_vertex(d)

        if d not in self.graph[o]:
            self.graph[o].append(d)

    def get_successors(self, v):
        return self.graph[v]

    def get_predecessors(self, v):
        lst = []
        for pre in self.graph:
            if v in self.graph[pre]:
                lst.append(pre)
        return sorted(lst)

    def get_adjacents(self, v):
        lst_suc = self.get_successors(v)
        lst_pre = self.get_predecessors(v)
        res = lst_pre
        for v1 in lst_suc:
            if v1 not in res:
                res.append(v1)
        return res

    def out_degree(self, v):
        return len(self.graph[v])

    def in_degree(self, v):
        return len(self.get_predecessors(v))

    def degree(self, v):
        return len(self.get_adjacents(v))

    def reachable_bfs(self, v):
        l = [v]
        res = [v]
        dist = 0
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l and elem != node:
                    l.append(elem)
        return res

    def reachable_dfs(self, v):
        l = [v]
        res = []
        while len(l) > 0:
            node = l.pop(0)
            if node != v: res.append(node)
            for elem in self.graph[node]:
                if elem not in res and elem not in l:
                    l.insert(0, elem)
        return res

    def distance(self, s, d):
        if s == d:
            return 0
        l = [(s, 0)]
        visited = [s]
        while len(l) > 0:
            node, dist = l.pop(0)
            for elem in self.graph[node]:
                if elem == d:
                    return dist + 1
                elif elem not in visited:
                    l.append((elem, dist + 1))
                    visited.append(elem)
        return None

    def reachable_with_dist(self, s):
        res = []
        l = [(s, 0)]
        while len(l) > 0:
            node, dist = l.pop(0)
            if node != s:    res.append((node, dist))
            for elem in self.graph[node]:
                if not is_in_tuple_list(l, elem) and not is_in_tuple_list(res, elem):
                    l.append((elem, dist + 1))
        return res

    def shortest_path(self, s, d):
        if s == d:
            return []
        l = [(s, [])]
        visited = [s]
        while len(l) > 0:
            node, preds = l.pop(0)
            for elem in self.graph[node]:
                if elem == d:
                    return preds + [node, elem]
                elif elem not in visited:
                    l.append((elem, preds + [node]))
                    visited.append(elem)
        return None

    def node_has_cycle(self, v):
        l = [v]
        res = False
        visited = [v]
        while len(l) > 0:
            node = l.pop(0)
            for elem in self.graph[node]:
                if elem == v:
                    return True
                elif elem not in visited:
                    l.append(elem)
                    visited.append(elem)
        return res

    def has_cycle(self):
        res = False
        for v in self.graph.keys():
            if self.node_has_cycle(v):    return True
        return res

    def all_degrees(self, deg_type="inout"):
        degs = {}
        for v in self.graph.keys():
            if deg_type == "out" or deg_type == "inout":
                degs[v] = len(self.graph[v])
            else:
                degs[v] = 0
        if deg_type == "in" or deg_type == "inout":
            for v in self.graph.keys():
                for d in self.graph[v]:
                    if deg_type == "in" or v not in self.graph[d]:
                        degs[d] = degs[d] + 1
        return degs

    def mean_degree(self, deg_type="inout"):
        degs = self.all_degrees(deg_type)
        return sum(degs.values()) / float(len(degs))

    def prob_degree(self, deg_type="inout"):
        degs = self.all_degrees(deg_type)
        res = {}
        for k in degs.keys():
            if degs[k] in res.keys():
                res[degs[k]] += 1
            else:
                res[degs[k]] = 1
        for k in res.keys():
            res[k] /= float(len(degs))
        return res

    def mean_distance(self):
        tot = 0
        num_reachable = 0
        for k in self.graph.keys():
            distk = self.reachable_with_dist(k)
            for _, dist in distk:
                tot += dist
            num_reachable += len(distk)
        meandist = float(tot) / num_reachable
        n = len(self.getNodes())
        return meandist, float(num_reachable) / ((n - 1) * n)

    def clustering_coef(self, v):
        adjs = self.get_adjacents(v)
        if len(adjs) <= 1:    return 0.0
        ligs = 0
        for i in adjs:
            for j in adjs:
                if i != j:
                    if j in self.graph[i] or i in self.graph[j]:
                        ligs = ligs + 1
        return float(ligs) / (len(adjs) * (len(adjs) - 1))

    def all_clustering_coefs(self):
        ccs = {}
        for k in self.graph.keys():
            ccs[k] = self.clustering_coef(k)
        return ccs

    def mean_clustering_coef(self):
        ccs = self.all_clustering_coefs()
        return sum(ccs.values()) / float(len(ccs))

    def mean_clustering_perdegree(self, deg_type="inout"):
        degs = self.all_degrees(deg_type)
        ccs = self.all_clustering_coefs()
        degs_k = {}
        for k in degs.keys():
            if degs[k] in degs_k.keys():
                degs_k[degs[k]].append(k)
            else:
                degs_k[degs[k]] = [k]
        ck = {}
        for k in degs_k.keys():
            tot = 0
            for v in degs_k[k]:    tot += ccs[v]
            ck[k] = float(tot) / len(degs_k[k])
        return ck

    def check_if_valid_path(self, p):
        if p[0] not in self.graph.keys():    return False
        for i in range(1, len(p)):
            if p[i] not in self.graph.keys() or p[i] not in self.graph[p[i - 1]]:
                return False
        return True

    def check_if_hamiltonian_path(self, p):
        if not self.check_if_valid_path(p):
            return False
        to_visit = list(self.get_nodes())
        if len(to_visit) != len(to_visit):
            return False

        for i in range(len(p)):
            if p[i] in to_visit:
                to_visit.remove(p[i])
            else:
                return False

        if not to_visit:
            return True
        else:
            return False

    def check_balanced_node(self, node):
        return self.in_degree(node) == self.out_degree(node)

    def check_balanced_graph(self):
        for n in self.graph.keys():
            if not self.check_balanced_node(n):    return False
        return True

    def eulerian_cycle(self):
        if not self.check_balanced_graph():
            return None
        edges_visit = list(self.get_edges())
        res = []
        while edges_visit:
            pair = edges_visit[0]
            i = 1
            if res != []:
                while pair[0] not in res:
                    pair = edges_visit[i]
                    i = i + 1
            edges_visit.remove(pair)
            start, nxt = pair
            cycle = [start, nxt]
            while nxt != start:
                for suc in self.graph[nxt]:
                    if (nxt, suc) in edges_visit:
                        print(pair)
                        pair = (nxt, suc)
                        nxt = suc
                        cycle.append(nxt)
                        print(cycle)
                        edges_visit.remove(pair)
            if not res:
                res = cycle
            else:
                pos = res.index(cycle[0])
                print("res: ", res, pos)
                for i in range(len(cycle) - 1): res.insert(pos + i + 1, cycle[i + 1])
                print("res: ", res)
        return res

    def check_nearly_balanced_graph(self):
        res = None, None
        for n in self.graph.keys():
            indeg = self.in_degree(n)
            outdeg = self.out_degree(n)
            if indeg - outdeg == 1 and res[1] is None:
                res = res[0], n
            elif indeg - outdeg == -1 and res[0] is None:
                res = n, res[1]
            elif indeg == outdeg:
                pass
            else:
                return None, None
        return res

    def eulerian_path(self):
        unb = self.check_nearly_balanced_graph()
        if unb[0] is None or unb[1] is None: return None
        self.graph[unb[1]].append(unb[0])
        cycle = self.eulerian_cycle()
        for i in range(len(cycle) - 1):
            if cycle[i] == unb[1] and cycle[i + 1] == unb[0]:
                break
        path = cycle[i + 1:] + cycle[1:i + 1]
        return path

    def closeness_centrality(self, node):
        dist = self.reachable_with_dist(node)
        if len(dist) == 0: return 0.0
        s = 0.0
        for d in dist: s += d[1]
        return len(dist) / s

    def highest_closeness(self, top=10):
        cc = {}
        for k in self.graph.keys(): cc[k] = self.closeness_centrality(k)
        ord_cl = sorted(list(cc.items()), key=lambda x: x[1], reverse=True)
        return list(map(lambda x: x[0], ord_cl[:top]))

    def betweenness_centrality(self, node):
        total_sp = 0
        sps_with_node = 0
        for s in self.graph.keys():
            for t in self.graph.keys():
                if s != t and s != node and t != node:
                    sp = self.shortest_path(s, t)
                    if sp is not None:
                        total_sp += 1
                        if node in sp: sps_with_node += 1
        return sps_with_node / total_sp
























