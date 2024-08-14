import bisect
from math import inf


class BaseNode:
    def neighbors_edges(self):
        # return 2 lists so that it can be zipped or one list omitted
        # and so that it matches the intuition of the name of the func
        raise NotImplementedError

    def node_search(self, goal_condition, cost_func, heuristic=None):
        h = heuristic or (lambda x: 0)

        costs = {self: cost_func(self)}
        frontier_queue = [(costs[self] + h(self), self)]
        to_explore = {self}

        while frontier_queue:
            _, current = frontier_queue.pop()
            if current in to_explore:
                to_explore.remove(current)
                if goal_condition(current):
                    return current, costs[current]
                for neighbor, _ in zip(*current.neighbors_edges()):
                    cost = cost_func(neighbor)
                    if cost < costs.get(neighbor, inf):
                        costs[neighbor] = cost
                        to_explore.add(neighbor)
                        h_cost = cost + h(neighbor)
                        bisect.insort_right(frontier_queue, (h_cost, neighbor),
                            key=lambda n: -n[0])
        return None, None

    def path_search(self, goal_condition, cost_func, heuristic=None):
        h = heuristic or (lambda x: 0)

        costs = {self: cost_func(self)}
        frontier_queue = [(costs[self] + h(self), self)]
        to_explore = {self}
        prev_node = dict()
        prev_edge = dict()

        def instruction_path(node):
            edges = []
            nodes = []
            while node in prev_node:
                nodes.insert(0, node)
                edges.insert(0, prev_edge[node])
                node = prev_node[node]
            return nodes, edges

        while frontier_queue:
            _, current = frontier_queue.pop()
            if current in to_explore:
                to_explore.remove(current)
                if goal_condition(current):
                    return *instruction_path(current), costs[current]
                for neighbor, edge in zip(*current.neighbors_edges()):
                    cost = cost_func(neighbor, edge, costs[current])
                    if cost < costs.get(neighbor, inf):
                        costs[neighbor] = cost
                        to_explore.add(neighbor)
                        prev_node[neighbor] = current
                        prev_edge[neighbor] = edge
                        h_cost = cost + h(neighbor)
                        bisect.insort_right(frontier_queue, (h_cost, neighbor),
                            key=lambda n: -n[0])
        return None, None, None
