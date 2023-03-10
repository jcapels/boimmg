from boimmgpy.service.network_handlers.graph import Graph


class MetabolicNetwork(Graph):

    def __init__(self, network_type="metabolite-reaction", split_rev=False):
        Graph.__init__(self, {})
        self.net_type = network_type
        self.node_types = {}
        if network_type == "metabolite-reaction":
            self.node_types["metabolite"] = []
            self.node_types["reaction"] = []
        elif network_type == "metabolite-metabolite":
            self.node_types["metabolite"] = []
        self.split_rev = split_rev

    def convert_metabolite_net(self, gmr):
        for m in gmr.node_types["metabolite"]:
            self.add_vertex(m)
            sucs = gmr.get_successors(m)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if m != s2:
                        self.add_edge(m, s2)

    def convert_reaction_graph(self, gmr):
        for r in gmr.node_types["reaction"]:
            self.add_vertex(r)
            sucs = gmr.get_successors(r)
            for s in sucs:
                sucs_r = gmr.get_successors(s)
                for s2 in sucs_r:
                    if r != s2:
                        self.add_edge(r, s2)

    def add_vertex_type(self, v, nodetype):
        self.add_vertex(v)
        self.node_types[nodetype].append(v)

    def get_nodes_type(self, node_type):
        if node_type in self.node_types:
            return self.node_types[node_type]
        else:
            return None

    def get_pathway(self):
        longer_path = []
        for node in self.get_nodes():
            path = self.reachable_bfs(node)
            if len(path) > len(longer_path):
                longer_path = path
        return longer_path

    def get_longest_pathway(self, node):
        longer_path = []
        path = self.reachable_bfs(node)
        if len(path) > len(longer_path):
            longer_path = path
        return longer_path

    def get_complete_pathway(self):
        pathway = self.get_pathway()
        node = pathway[0]
        return self.reachable_bfs(node)

    def get_all_targets(self):
        res = []
        for token in self.graph:
            successors = self.get_successors(token)
            predecessors = self.get_predecessors(token)
            if not successors:
                res.append(token)

            elif successors == predecessors:
                res.append(token)

        return res

    def prune_redundant_cycles(self):
        old_graph = self.graph.copy()

        if len(old_graph) > 1:
            for token in old_graph:
                successors = self.get_successors(token)
                predecessors = self.get_predecessors(token)

                if successors == predecessors:
                    self.delete_node(token)

    def get_starting_point(self):
        res = []
        for token in self.graph:
            predecessors = self.get_predecessors(token)
            if not predecessors:
                res.append(token)
        return res

    def get_all_target_pathways(self):
        starting_points = self.get_starting_point()
        targets = self.get_all_targets()

        res = []
        if len(self.graph) > 1:
            for point in starting_points:
                for target in targets:
                    sh_path = self.shortest_path(point, target)

                    if sh_path:
                        res.append(sh_path)
        else:
            res = [starting_points]
        return res
