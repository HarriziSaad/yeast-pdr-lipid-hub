import re

class Tree:

    def __init__(self, newick):
        self.adj = {}
        self.name = {}
        self.edge = {}
        self._parse(newick)

    def _parse(self, s):
        s = s.strip().rstrip(';')
        pos, n_nodes, n_edges = (0, 0, 0)
        stack = []
        cur_children = None
        i = 0
        L = len(s)

        def new_node():
            nonlocal n_nodes
            n_nodes += 1
            self.adj[n_nodes] = []
            return n_nodes

        def parse_tail(i):
            j = i
            while j < L and s[j] not in ':,)':
                j += 1
            label = s[i:j]
            length = None
            if j < L and s[j] == ':':
                k = j + 1
                while k < L and s[k] not in ',)':
                    k += 1
                length = float(s[j + 1:k])
                j = k
            return (label, length, j)
        stack = []
        children_stack = []
        root = None
        while i < L:
            ch = s[i]
            if ch == '(':
                children_stack.append([])
                i += 1
            elif ch == ',':
                i += 1
            elif ch == ')':
                kids = children_stack.pop()
                label, length, i = parse_tail(i + 1)
                nid = new_node()
                for cid, clen, clab in kids:
                    self._link(nid, cid, clen, clab)
                if children_stack:
                    children_stack[-1].append((nid, length, label))
                else:
                    root = nid
            else:
                j = i
                while j < L and s[j] not in ':,)':
                    j += 1
                nm = s[i:j]
                length = None
                if j < L and s[j] == ':':
                    k = j + 1
                    while k < L and s[k] not in ',)':
                        k += 1
                    length = float(s[j + 1:k])
                    j = k
                nid = new_node()
                self.name[nid] = nm
                children_stack[-1].append((nid, length, ''))
                i = j
        self.root0 = root
        if len(self.adj[root]) == 2:
            (a, ea), (b, eb) = self.adj[root]
            e1, e2 = (self.edge[ea], self.edge[eb])
            self.adj[a] = [(x, e) for x, e in self.adj[a] if x != root]
            self.adj[b] = [(x, e) for x, e in self.adj[b] if x != root]
            length = (e1['length'] or 0) + (e2['length'] or 0)
            sh = e1['sh'] if e1['sh'] is not None else e2['sh']
            ufb = e1['ufb'] if e1['ufb'] is not None else e2['ufb']
            eid = ea
            self.edge[eid] = dict(a=a, b=b, length=length, sh=sh, ufb=ufb)
            self.adj[a].append((b, eid))
            self.adj[b].append((a, eid))
            del self.edge[eb]
            del self.adj[root]

    def _link(self, p, c, length, label):
        eid = len(self.edge) + 1
        while eid in self.edge:
            eid += 1
        sh = ufb = None
        if label:
            parts = label.split('/')
            try:
                sh = float(parts[0])
                if len(parts) > 1:
                    ufb = float(parts[1])
            except ValueError:
                pass
        self.edge[eid] = dict(a=p, b=c, length=length, sh=sh, ufb=ufb)
        self.adj[p].append((c, eid))
        self.adj[c].append((p, eid))

    def leaves(self):
        return [n for n in self.adj if n in self.name]

    def root_on_edge(self, eid):
        return Rooted(self, eid)

class Rooted:

    def __init__(self, tree, eid):
        self.t = tree
        self.root_edge = eid
        e = tree.edge[eid]
        self.children = {}
        self.parent = {}
        self.stem = {}
        self.order = []
        for top in (e['a'], e['b']):
            other = e['b'] if top == e['a'] else e['a']
            self.parent[top] = 'ROOT'
            self.stem[top] = eid
            st = [top]
            seen = {top, other}
            while st:
                n = st.pop()
                self.order.append(n)
                kids = []
                for m, ee in tree.adj[n]:
                    if m in seen:
                        continue
                    seen.add(m)
                    kids.append(m)
                    self.parent[m] = n
                    self.stem[m] = ee
                    st.append(m)
                self.children[n] = kids
        self.top = (e['a'], e['b'])
        self._leafset = {}
        for n in reversed(self.order):
            if n in tree.name:
                self._leafset[n] = frozenset([n])
            else:
                u = set()
                for k in self.children[n]:
                    u |= self._leafset[k]
                self._leafset[n] = frozenset(u)

    def leafnodes(self, n):
        return self._leafset[n]

    def leafnames(self, n):
        return {self.t.name[x] for x in self._leafset[n]}

    def mrca(self, nodes):
        nodes = list(nodes)
        s = set(nodes)
        best, bs = (None, 10 ** 9)
        n = nodes[0]
        while True:
            if s <= self._leafset[n]:
                return n
            p = self.parent[n]
            if p == 'ROOT':
                return None
            n = p

    def support(self, n):
        e = self.t.edge[self.stem[n]]
        return (e['sh'], e['ufb'])

    def ancestors(self, n):
        out = []
        while n != 'ROOT':
            out.append(n)
            n = self.parent[n]
        return out
