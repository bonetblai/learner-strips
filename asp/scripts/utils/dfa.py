from pathlib import Path
from random import choices, choice
import re

class DFA:
    def __init__(self, fname: Path = None, logger = None):
        self.fname = fname
        self.num_nodes = 0
        self.num_edges = 0
        self.root = -1
        self.labels = dict()
        self.goals = []
        self.successors = []
        self.logger = logger

        if self.fname is not None:
            if self.fname.suffix == '.dfa':
                if self.logger: self.logger.info(f"Reading '{self.fname}'")
                with self.fname.open('r') as fd:
                    line = fd.readline().split()
                    assert line[0] == 'dfa'
                    self.num_nodes = int(line[1])
                    self.root = int(line[2])
                    self.labels = fd.readline().split()[1:]
                    self.labels = dict(zip(self.labels, range(1, 1 + len(self.labels))))
                    self.goals = [ int(g) for g in fd.readline().split()[1::2] ]
                    if self.logger: self.logger.info(f'num_nodes={self.num_nodes}, root={self.root}, labels={self.labels}, goals={self.goals}')
                    self.successors = []
                    self.num_edges = 0
                    for node in range(self.num_nodes):
                        raw = fd.readline().split()[1:]
                        node_successors = list(zip(raw[0::2], map(int, raw[1::2])))
                        if self.logger: self.logger.debug(f'successors of node {node}: {node_successors}')
                        self.successors.append(node_successors)
                        self.num_edges += len(node_successors)
            elif self.fname.suffix == '.lp':
                if self.logger: self.logger.info(f"Reading '{self.fname}'")
                with self.fname.open('r') as fd:
                    lines = [ line.rstrip('\n') for line in fd.readlines() if line[0] != '%' ]
                    nodes = [ int(re.search('node\((\d*)\).', line).group(1)) for line in lines if line[:5] == 'node(' ]
                    labels = [ re.search('labelname\((\d*),"(.*)"\).', line).groups() for line in lines if line[:10] == 'labelname(' ]
                    labels = sorted([ (label, int(index)) for (index, label) in labels ], key=lambda item: item[1])
                    edges = [ tuple(map(int, re.search('tlabel\(\((\d*),(\d*)\),(\d*)\).', line).groups())) for line in lines if line[:7] == 'tlabel(' ]
                    self.num_nodes = len(nodes)
                    self.num_edges = len(edges)
                    self.labels = dict(labels)
                    if self.logger: self.logger.info(f'num_nodes={self.num_nodes}, root={self.root}, labels={self.labels}, goals={self.goals}')
                    self.successors = [ [] for i in range(self.num_nodes) ]
                    for (src, dst, label) in edges:
                        edge = (labels[label-1][0], dst)
                        self.successors[src].append(edge)
                    for node in range(self.num_nodes):
                        if self.logger: self.logger.debug(f'successors of node {node}: {self.successors[node]}')
            else:
                raise RuntimeError(f'unknown suffix {self.fname.suffix} for DFA file')
            assert len(self.successors) == self.num_nodes

    def dump_as_lp(self, fd):
        fd.write(f'% {self.num_nodes} nodes, {self.num_edges} edges\n')
        for node in range(self.num_nodes):
            fd.write(f'node({node}).\n')
        for label in self.labels:
            fd.write(f'labelname({self.labels[label]},"{label}").\n')
        for node, node_successors in enumerate(self.successors):
            for (label, next) in node_successors:
                fd.write(f'edge(({node},{next})).\n')
                fd.write(f'tlabel(({node},{next}),{self.labels[label]}).\n')

    def sample_nodes(self, n, repeat=True, avoid=None):
        nodes = list(range(self.num_nodes))
        if avoid is not None:
            eligible = [ node for node in nodes if node not in avoid ]
        else:
            eligible = nodes
        return choices(eligible, k=n)

    def sample_path(self, src, max_length, repeat=True):
        visited = set()
        path = [ src ]
        visited.add(src)
        while len(path) < max_length:
            successors = [ node for (label, node) in self.successors[path[-1]] if node not in visited ]
            if len(successors) == 0: break
            current = choice(successors)
            path.append(current)
            visited.add(current)
        return path

