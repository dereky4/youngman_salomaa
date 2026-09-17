from collections import defaultdict, deque
from typing import Optional
import random
import math
import copy
from sep import seperator
import itertools
import networkx as nx
import icdfa


def rand(p: float = 0.5) -> bool:
    return random.random() < p


class FA:

    def __init__(self, n, edges, require_trim=True) -> None:
        self.n = n
        self.outgoing = defaultdict(dict)
        self.incoming = defaultdict(dict)
        for q0, q1, w in edges:
            assert 0 <= q0 < n
            assert 0 <= q1 < n
            self.outgoing[q0][q1] = w
            self.incoming[q1][q0] = w
        if require_trim:
            assert self.trim()
        else:
            assert self.reachable()

    def reachable(self) -> bool:
        visited = {0}
        queue = deque([0])
        while queue:
            q0 = queue.popleft()
            for q1 in self.outgoing[q0]:
                if q1 == self.n - 1:
                    return True
                if q1 in visited:
                    continue
                visited.add(q1)
                queue.append(q1)
        return False

    def trim(self) -> bool:
        return self.allReachable() and self.reachableFromAll()

    def reachableStates(self):
        visited = {0}
        queue = deque([0])
        while queue:
            q0 = queue.popleft()
            for q1 in self.outgoing[q0]:
                if q1 in visited:
                    continue
                visited.add(q1)
                queue.append(q1)
        return visited

    def allReachable(self) -> bool:
        return len(self.reachableStates()) == self.n

    def usefulStates(self):
        visited = {self.n - 1}
        queue = deque([self.n - 1])
        while queue:
            q0 = queue.popleft()
            for q1 in self.incoming[q0]:
                if q1 in visited:
                    continue
                visited.add(q1)
                queue.append(q1)
        return visited

    def reachableFromAll(self) -> bool:
        return len(self.usefulStates()) == self.n

    def reachableAndUsefulStates(self):
        return self.reachableStates() & self.usefulStates()


class WFA(FA):

    def __init__(
        self, n: int, edges: list[tuple[int, int, int]], require_trim=True
    ) -> None:
        super().__init__(n, edges, require_trim)

    def eliminateState(self, q, states, outgoing, incoming):
        assert 0 < q < self.n - 1
        states.remove(q)
        loop = outgoing[q].get(q, 0)
        for q0, w0 in incoming[q].items():
            assert incoming[q][q0] == outgoing[q0][q]
            if q0 not in states:
                continue
            for q1, w1 in outgoing[q].items():
                assert outgoing[q][q1] == incoming[q1][q]
                if q1 not in states:
                    continue
                cur = outgoing[q0].get(q1, 0)
                updated = cur + w0 + loop + w1
                outgoing[q0][q1] = updated
                incoming[q1][q0] = updated

    def weight(self, q, states, outgoing, incoming):
        inq = sum(q0 != q and q0 in states for q0 in incoming[q])
        outq = sum(q1 != q and q1 in states for q1 in outgoing[q])
        loop = outgoing[q].get(q, 0)
        w = inq * loop * outq
        w += outq * sum(
            w0 for q0, w0 in incoming[q].items() if q0 != q and q0 in states
        )
        w += inq * sum(w1 for q1, w1 in outgoing[q].items() if q1 != q and q1 in states)
        return w

    def gates(self, q, states, outgoing, incoming):
        res = 0
        for q0 in states:
            if q == q0:
                continue
            for q1 in states:
                if q == q1:
                    continue
                if q1 not in outgoing[q0] and q in outgoing[q0] and q in incoming[q1]:
                    res += outgoing[q0][q] + outgoing[q].get(q, 0) + incoming[q1][q]
        return res

    def gatesWeight(self, q, states, outgoing, incoming):
        return (
            self.gates(q, states, outgoing, incoming),
            self.weight(q, states, outgoing, incoming),
        )

    def optEliminate(self, block, states, outgoing, incoming):
        originalOutgoing = copy.deepcopy(outgoing)
        originalIncoming = copy.deepcopy(incoming)
        remaining = [q for q in states if q not in block]
        for q0 in remaining:
            for q1 in remaining:
                if q1 in outgoing[q0]:
                    outgoing[q0][q1] = math.inf
                if q1 in incoming[q0]:
                    incoming[q0][q1] = math.inf
        orders = list(itertools.permutations(block))
        for _, order in enumerate(orders):
            originalOutgoingCopy = copy.deepcopy(originalOutgoing)
            originalIncomingCopy = copy.deepcopy(originalIncoming)
            for q in order:
                self.eliminateState(
                    q, states, originalOutgoingCopy, originalIncomingCopy
                )
            for q0 in remaining:
                for q1 in remaining:
                    if q1 in originalOutgoingCopy[q0]:
                        outgoing[q0][q1] = min(
                            outgoing[q0].get(q1, math.inf), originalOutgoingCopy[q0][q1]
                        )
                    if q1 in originalIncomingCopy[q0]:
                        incoming[q0][q1] = min(
                            incoming[q0].get(q1, math.inf), originalIncomingCopy[q0][q1]
                        )
            for q in block:
                states.add(q)
        for q in block:
            states.remove(q)

    def getBblock(self, states, outgoing, incoming, b):
        block = [None] * b
        weights = [math.inf] * b
        i = 0
        for q in states:
            if q == 0 or q == self.n - 1:
                continue
            i += 1
            w = self.weight(q, states, outgoing, incoming)
            for j in range(b):
                if w >= weights[j]:
                    continue
                for k in range(b - 1, j, -1):
                    block[k], weights[k] = block[k - 1], weights[k - 1]
                block[j], weights[j] = q, w
                break
        assert i >= b
        assert None not in block
        assert len(set(block)) == len(block)
        assert all(weights[idx] <= weights[idx + 1] for idx in range(b - 1))
        return block

    def totalWidth(
        self,
        order: Optional[list[int]] = None,
        w=False,
        wopt=False,
        bs=None,
        gate=False,
    ) -> int:
        if not order:
            order = [i for i in range(1, self.n - 1)]
        # assert set(order) == {i for i in range(1, self.n - 1)}
        states = {i for i in range(self.n)}
        # assert 0 in states and self.n - 1 in states
        outgoing = copy.deepcopy(self.outgoing)
        incoming = copy.deepcopy(self.incoming)
        if not w and not gate:
            for q in order:
                self.eliminateState(q, states, outgoing, incoming)
        elif w:
            if wopt:
                for b in bs:
                    block = self.getBblock(states, outgoing, incoming, b)
                    self.optEliminate(block, states, outgoing, incoming)
            while len(states) > 2:
                qm = min(
                    (q for q in states if 0 < q < self.n - 1),
                    key=lambda q: self.weight(q, states, outgoing, incoming),
                )
                self.eliminateState(qm, states, outgoing, incoming)
        else:
            while len(states) > 2:
                qm = min(
                    (q for q in states if 0 < q < self.n - 1),
                    key=lambda q: self.gatesWeight(q, states, outgoing, incoming),
                )
                self.eliminateState(qm, states, outgoing, incoming)
        assert len(states) == 2
        f = self.n - 1
        width = outgoing[0].get(0, 0) + outgoing[0].get(f, 0)
        if 0 in outgoing[f]:
            width += outgoing[f][0] + width
        width += outgoing[f].get(f, 0)
        return width

    def eliminateStates(self, order):
        states = {i for i in range(self.n)}
        outgoing = copy.deepcopy(self.outgoing)
        incoming = copy.deepcopy(self.incoming)
        for q in order:
            self.eliminateState(q, states, outgoing, incoming)
        return states, outgoing, incoming

    def sepH(self, blocks, w=True, opt=False):
        states = {i for i in range(self.n)}
        outgoing = copy.deepcopy(self.outgoing)
        incoming = copy.deepcopy(self.incoming)
        for b, block in enumerate(blocks):
            if opt and b != len(blocks) - 1:
                self.optEliminate(block, states, outgoing, incoming)
            elif w:
                for _ in range(len(block)):
                    qm = min(
                        (q for q in block if q in states),
                        key=lambda q: self.weight(q, states, outgoing, incoming),
                    )
                    self.eliminateState(qm, states, outgoing, incoming)
            else:
                random.shuffle(block)
                for q in block:
                    self.eliminateState(q, states, outgoing, incoming)
        assert len(states) == 2
        f = self.n - 1
        width = outgoing[0].get(0, 0) + outgoing[0].get(f, 0)
        if 0 in outgoing[f]:
            width += outgoing[f][0] + width
        width += outgoing[f].get(f, 0)
        return width


class RE:
    OPS = ("+", "", "*")
    BIN = ("+", "")
    EPSILON = "3"

    class Node:
        def __init__(self, c, left=None, right=None) -> None:
            self.c = c
            self.left = left
            self.right = right
            self.size = (
                int(c not in RE.OPS and c != RE.EPSILON)
                + (left.size if left else 0)
                + (right.size if right else 0)
            )

    def __init__(self, c=None, left=None, right=None):
        if c is None:
            self.root = None
            return
        self.root = RE.Node(
            c, left.root if left else None, right.root if right else None
        )

    def isEpsilon(self):
        return self.root and self.root.c == RE.EPSILON

    @staticmethod
    def union(left, right):
        if left.root is None:
            return right
        if right.root is None:
            return left
        return RE("+", left, right)

    @staticmethod
    def concat(left, right):
        if left.root is None or right.root is None:
            return RE()
        if left.root.c == RE.EPSILON:
            return right
        if right.root.c == RE.EPSILON:
            return left
        return RE("", left, right)

    @staticmethod
    def star(left):
        if left.root is None or left.root.c == RE.EPSILON:
            return RE(RE.EPSILON)
        return RE("*", left)

    def __str__(self) -> str:
        return self.__repr__()

    def __repr__(self) -> str:
        if self.root is None:
            return "emptyset"
        buffer = []

        def dfs(node):
            if node.c not in RE.OPS:
                buffer.append(node.c)
                return
            bracket = (node.c == "*" and node.left.c in RE.BIN) or node.left.c == "+"
            if bracket:
                buffer.append("(")
            dfs(node.left)
            if bracket:
                buffer.append(")")
            if node.c == "*":
                buffer.append("*")
                return
            if node.c == "+":
                buffer.append(" + ")
            if node.right.c in RE.BIN:
                buffer.append("(")
            dfs(node.right)
            if node.right.c in RE.BIN:
                buffer.append(")")

        dfs(self.root)
        return "".join(buffer)

    def count_nodes(self):
        def dfs(cur):
            if cur is None:
                return 0
            return 1 + dfs(cur.left) + dfs(cur.right)

        return dfs(self.root) if self.root else None

    def size(self):
        return self.root.size if self.root else 0


class NFA(FA):

    def __init__(self, n: int, edges: list[tuple[int, int, str]]):
        edgesMap = defaultdict(lambda: defaultdict(RE))
        for q0, q1, c in edges:
            edgesMap[q0][q1] = RE.union(edgesMap[q0][q1], RE(c))
        newEdges = [
            (q0, q1, edgesMap[q0][q1]) for q0 in edgesMap for q1 in edgesMap[q0]
        ]
        super().__init__(n, newEdges)

    def eliminateState(self, q, states, outgoing, incoming):
        assert q != 0 and q != self.n - 1
        states.remove(q)
        loop = outgoing[q].get(q, RE(RE.EPSILON))
        loop = RE.star(loop)
        for q0, t1 in incoming[q].items():
            if q0 not in states:
                continue
            for q1, t2 in outgoing[q].items():
                if q1 not in states:
                    continue
                cur = outgoing[q0].get(q1, RE())
                newPath = RE.concat(RE.concat(t1, loop), t2)
                updated = RE.union(cur, newPath)
                outgoing[q0][q1] = updated
                incoming[q1][q0] = updated

    def weight(self, q, states, outgoing, incoming):
        inq = sum(q0 != q and q0 in states for q0 in incoming[q])
        outq = sum(q1 != q and q1 in states for q1 in outgoing[q])
        loop = outgoing[q].get(q, RE(RE.EPSILON)).size()
        w = inq * loop * outq
        w += outq * sum(
            w0.size() for q0, w0 in incoming[q].items() if q0 != q and q0 in states
        )
        w += inq * sum(
            w1.size() for q1, w1 in outgoing[q].items() if q1 != q and q1 in states
        )
        return w

    def getRE(self, order: Optional[list[int]] = None, w=False) -> RE:
        if not order:
            order = [i for i in range(1, self.n - 1)]
        states = {i for i in range(self.n)}
        outgoing = copy.deepcopy(self.outgoing)
        incoming = copy.deepcopy(self.incoming)
        if not w:
            for q in order:
                self.eliminateState(q, states, outgoing, incoming)
        else:
            worder = []
            for _ in range(self.n - 2):
                qm = min(
                    (q for q in range(1, self.n - 1) if q in states),
                    key=lambda q: self.weight(q, states, outgoing, incoming),
                )
                worder.append(qm)
                self.eliminateState(qm, states, outgoing, incoming)
            print("weight", worder)
        f = self.n - 1
        startLoop = outgoing[0].get(0, RE(RE.EPSILON))
        startLoop = RE.star(startLoop)
        res = RE.concat(startLoop, outgoing[0][f])
        if 0 in outgoing[f]:
            bigLoop = RE.concat(RE.concat(outgoing[f][0], startLoop), outgoing[0][f])
        else:
            bigLoop = RE()
        endLoop = outgoing[f].get(f, RE(RE.EPSILON))
        if not endLoop.isEpsilon():
            bigLoop = RE.union(endLoop, bigLoop)
        return RE.concat(res, RE.star(bigLoop))

    def sepH(self, blocks) -> RE:
        states = {i for i in range(self.n)}
        outgoing = copy.deepcopy(self.outgoing)
        incoming = copy.deepcopy(self.incoming)
        worder = []
        for block in blocks:
            for _ in range(len(block)):
                qm = min(
                    (q for q in block if q in states),
                    key=lambda q: self.weight(q, states, outgoing, incoming),
                )
                worder.append(qm)
                self.eliminateState(qm, states, outgoing, incoming)
        print("blocks", worder)
        f = self.n - 1
        startLoop = outgoing[0].get(0, RE(RE.EPSILON))
        startLoop = RE.star(startLoop)
        res = RE.concat(startLoop, outgoing[0][f])
        if 0 in outgoing[f]:
            bigLoop = RE.concat(RE.concat(outgoing[f][0], startLoop), outgoing[0][f])
        else:
            bigLoop = RE()
        endLoop = outgoing[f].get(f, RE(RE.EPSILON))
        if not endLoop.isEpsilon():
            bigLoop = RE.union(endLoop, bigLoop)
        return RE.concat(res, RE.star(bigLoop))


class DFA(NFA):
    def __init__(self, n: int, edges: list[tuple[int, int, str]]):
        self.delta = defaultdict(dict)
        self.sigma = set()
        for q0, q1, c in edges:
            assert c not in self.delta[q0] or self.delta[q0][c] == q1
            self.delta[q0][c] = q1
            self.sigma.add(c)
        super().__init__(n, edges)


def generate_A_n(n):
    edges = [(i, j, f"({i}, {j})") for i in range(n) for j in range(n)]
    return DFA(n, edges)


def generate_wA_n(n):
    edges = [(i, j, 1) for i in range(n) for j in range(n)]
    return WFA(n, edges)


def verify_prop6():
    for n in range(3, 8):
        A_n = generate_A_n(n)
        re = A_n.getRE()
        s = str(re)
        w = re.size()
        for j in range(1, n - 2 + 1):
            for k in range(1, n - 2 + 1):
                assert s.count(f"({j}, {k})") == w * 2 ** (-(j + k))
            assert 3 * s.count(f"(0, {j})") == w * 2 ** (-(n + j - 3))
            assert s.count(f"({j}, 0)") == s.count(f"({j}, {n - 1})")
            assert s.count(f"({j}, 0)") == w * 2 ** (-(n + j - 1))
            assert 3 * s.count(f"({n - 1}, {j})") == w * 2 ** (-(n + j - 2))
        assert s.count("(0, 0)") == s.count(f"(0, {n - 1})")
        assert 3 * s.count("(0, 0)") == w * 4 ** (-(n - 2))
        assert s.count(f"({n - 1}, 0)") == s.count(f"({n - 1}, {n - 1})")
        assert 6 * s.count(f"({n - 1}, 0)") == w * 4 ** (-(n - 2))


def add_edge(edges, Sigma, j, k):
    p = 0.5
    for a in Sigma:
        if rand(p):
            edges.append((j, k, a))
        if j != k and rand(p):
            edges.append((k, j, a))


def generate_artificial_seperators(n, c, w, Sigma):
    edges = []
    blocks = []

    def helper(n, end):
        m1 = math.floor(c * n ** (1 - w))
        if m1 >= n:
            blocks.append([i for i in range(end - n, end)])
            for j in range(end - n, end):
                for k in range(end - n, end):
                    add_edge(edges, Sigma, j, k)
            return
        m2 = (n - m1 + 1) // 2
        m3 = n - m1 - m2
        helper(m3, end - m1 - m2)
        helper(m2, end - m1)
        blocks.append([i for i in range(end - m1, end)])
        for j in range(end - m1, end):
            for k in range(end - n, end):
                add_edge(edges, Sigma, j, k)

    helper(n, n + 1)
    edges.append((0, n, Sigma[0]))
    edges.append((n, n + 1, Sigma[0]))
    return NFA(n + 2, edges), blocks, edges


def grid_wfa(rows, cols):
    n = rows * cols + 2
    edges = []
    for j in range(1, n):
        edges.append((0, j, 1))
        edges.append((j, 0, 1))
    for j in range(1, n - 1):
        edges.append((n - 1, j, 1))
        edges.append((j, n - 1, 1))
    for r in range(rows):
        for c in range(cols):
            j = cols * r + c + 1
            assert 1 <= j < n - 1
            if r < rows - 1:
                k = cols * (r + 1) + c + 1
                assert 1 <= k < n - 1
                edges.append((j, k, 1))
                edges.append((k, j, 1))
            if c < cols - 1:
                k = cols * r + (c + 1) + 1
                assert 1 <= k < n - 1
                edges.append((j, k, 1))
                edges.append((k, j, 1))

    blocks = []

    def gen_blocks(top, bottom, left, right):
        if top > bottom or left > right:
            return
        if bottom - top <= 1 and right - left <= 1:
            blocks.append(
                [
                    cols * r + c + 1
                    for r in range(top, bottom + 1)
                    for c in range(left, right + 1)
                ]
            )
            return
        if bottom - top >= right - left:
            mid = top + (bottom - top) // 2
            gen_blocks(top, mid - 1, left, right)
            gen_blocks(mid + 1, bottom, left, right)
            blocks.append([cols * mid + c + 1 for c in range(left, right + 1)])
        else:
            mid = left + (right - left) // 2
            gen_blocks(top, bottom, left, mid - 1)
            gen_blocks(top, bottom, mid + 1, right)
            blocks.append([cols * r + mid + 1 for r in range(top, bottom + 1)])

    gen_blocks(0, rows - 1, 0, cols - 1)
    return WFA(n, edges), blocks


def bridge_wfa(n):
    mid = n // 2
    edges = []
    for j in range(mid + 1):
        for k in range(mid + 1):
            edges.append((j, k, 1))
    for j in range(mid, n):
        for k in range(mid, n):
            edges.append((j, k, 1))
    blocks = [[j for j in range(1, mid)], [j for j in range(mid + 1, n - 1)], [mid]]
    return WFA(n, edges), blocks


def circular3regular_wfa(c):
    edges = []
    n = 2 * c + 2
    for j in range(1, c + 1):
        for k in (j + c, (j % c) + 1):
            edges.append((j, k, 1))
            edges.append((k, j, 1))
    sep = [1, 1 + c, (c // 2) + 1, (c // 2) + 1 + c]
    nonsep = [j for j in range(2, n - 1) if j not in sep]
    for j in (2,):
        edges.append((0, j, 1))
        # edges.append((j, 0, 1))
        # edges.append((n - 1, j, 1))
        edges.append((j, n - 1, 1))
    return WFA(n, edges), [nonsep, sep]


def assert_seperator(B, C, edges):
    edge_set = {(j, k) for j, k, _ in edges}
    for j in B:
        for k in C:
            assert (j, k) not in edge_set
            assert (k, j) not in edge_set


def is_planar(edges, exclude):
    graph = nx.Graph()
    for j, k, _ in edges:
        if j < k and j not in exclude and k not in exclude:
            graph.add_edge(j, k)
    res = nx.check_planarity(graph)[0]
    return res


def test_eliminate_most(n, p, m):
    while True:
        edges = [(i, j, 1) for i in range(n) for j in range(n) if rand(p)]
        try:
            wfa = WFA(n, edges)
            break
        except AssertionError:
            continue
    elim = random.sample([i for i in range(1, n - 1)], m)
    remaining, outgoing, _ = wfa.eliminateStates(elim)
    total_degree = sum(
        sum(outgoing[q0].get(q1, 0) > 0 for q1 in remaining) for q0 in remaining
    )
    return total_degree / (len(remaining) * len(remaining))


def eliminate_most_main(ns):
    for n in ns:
        for m in (5, 10, 15, 20):
            if (n == 500 and m != 20) or m >= n - 1:
                continue
            for p in (0.5,):
                TRIALS = 10000
                total_rho = perfect = 0
                for trial in range(TRIALS):
                    rho = test_eliminate_most(n, p, m)
                    total_rho += rho
                    perfect += rho == 1.0
                    print(trial / TRIALS)
                print("DONE!", total_rho / TRIALS, perfect)
                with open("results.txt", "a") as f:
                    f.write("----------eliminate-most---------\n")
                    f.write(f"n={n}, m={m}, p={p}\n")
                    f.write(f"{total_rho / TRIALS}, {perfect}\n")


def testGateMotivation(n, sep):
    s = (n - 2 - sep) // 2
    B = {j for j in range(sep + 1, sep + 1 + s)}
    S = {j for j in range(sep + 1)}
    S.add(n - 1)
    C = {j for j in range(sep + 1 + s, n - 1)}

    print(f"B = {min(B)}..{max(B)}")
    print(f"S = {sorted(S)}")
    print(f"C = {min(C)}..{max(C)}")
    print(f"|B| = {len(B)}, |S| = {len(S)}, |C| = {len(C)}")

    edges = []
    for j in range(n):
        for k in range(n):
            if (j in B and k in B) or (j in C and k in C):
                edges.append((j, k, 5))
            elif j in S or k in S:
                edges.append((j, k, 1))

    wfa = WFA(n, edges)
    weight = wfa.totalWidth(w=True)
    gate = wfa.totalWidth(gate=True)
    order = list(B) + list(C) + [j for j in S if j not in (0, n - 1)]
    separator = wfa.totalWidth(order)

    print(f"Weight heuristic: {weight}")
    print(f"Gate score heuristic: {gate}")
    print(f"Separator order: {separator}")
    print(f"Weight / separator: {weight / separator}")
    

def testFig():
    B = {1, 2}
    C = {3, 4}
    S = {5, 6}
    edges = [(0, 6, 1), (6, 7, 1)]
    for leaf in (B, C):
        for j in leaf:
            for k in leaf:
                if j != k:
                    edges.append((j, k, 7))
    for j in S:
        for k in range(1, 7):
            if j == k or (k in S and j > k):
                continue
            edges.append((j, k, 1))
            edges.append((k, j, 1))
    wfa = WFA(8, edges)
    print(wfa.totalWidth(w=True))
    print(wfa.totalWidth(gate=True))
    print(wfa.totalWidth(order=[1, 2, 3, 4, 5, 6]))


def icdfa_test(gen):
    edges, n = next(gen)
    wfa = WFA(n, edges, require_trim=False)
    B, S, C = seperator(n, edges)
    # print(len(B), len(C))
    blocks = [B, C, S]
    order = [j for j in range(1, n - 1)]
    random.shuffle(order)
    M0 = wfa.totalWidth(order)
    M1 = wfa.totalWidth(w=True)
    M2 = wfa.sepH(blocks, w=False, opt=False)
    M3 = wfa.sepH(blocks, w=True, opt=False)
    M4 = wfa.totalWidth(gate=True)
    M5 = wfa.sepH(blocks, w=True, opt=True)
    low, high = min(len(B), len(C)), max(len(B), len(C))
    M6 = wfa.totalWidth(w=True, wopt=True, bs=[high, low])
    return M0, M1, M2, M3, M4, M5, M6


def icdfa_main(ns):
    k = 7
    TRIALS = 10000
    for n in ns:
        for alph in (5, n // 2, n):
            gen = icdfa.generator(n, alph)
            total_raw = [0] * k
            sole_wins = [0] * k
            wins = [0] * k
            total_low = 0
            for trial in range(1, TRIALS + 1):
                res = icdfa_test(gen)
                low = min(res)
                total_low += low
                low_count = sum(size == low for size in res)
                for i in range(k):
                    total_raw[i] += res[i]
                    if res[i] == low:
                        sole_wins[i] += low_count == 1
                        wins[i] += 1
                print(f"trial {trial}/{TRIALS}")
                if not trial % 1000:
                    with open("results.txt", "a") as f:
                        f.write("-----------dense------------\n")
                        f.write(f"n={n}, alph={alph}\n")
                        f.write(f"{sole_wins}\n")
                        total_raw.append(total_low)
                        f.write(f"{total_raw}\n")
                        total_raw.pop()
            with open("results.txt", "a") as f:
                f.write("-----------dense------------\n")
                f.write(f"n={n}, alph={alph}\n")
                f.write(f"{sole_wins}\n")
                total_raw.append(total_low)
                f.write(f"{[total / total_raw[1] for total in total_raw]}\n")


def icdfa_sparse_test(gen):
    edges, n = next(gen)
    wfa = WFA(n, edges, require_trim=False)
    order = [j for j in range(1, n - 1)]
    blocks = []

    def helper(sub):
        if len(sub) <= 3:
            blocks.append(sub)
            return
        B, S, C = seperator(n, edges, sub)
        helper(B)
        helper(C)
        if S:
            blocks.append(S)

    helper([i for i in range(1, n - 1)])
    # assert all(any(i in block for block in blocks) for i in range(1, n - 1))
    # print([len(block) for block in blocks])
    random.shuffle(order)
    M0 = wfa.totalWidth(order)
    M1 = wfa.totalWidth(w=True)
    M2 = wfa.sepH(blocks, w=False, opt=False)
    M3 = wfa.sepH(blocks, w=True, opt=False)
    M4 = wfa.totalWidth(gate=True)
    M5 = wfa.sepH(blocks, w=True, opt=True)
    return M0, M1, M2, M3, M4, M5


def icdfa_sparse_main(ns):
    k = 6
    TRIALS = 10000
    for alph in (2, 5):
        for n in ns:
            if n == 20 and alph == 5:
                continue
            gen = icdfa.generator(n, alph)
            total_raw = [0] * k
            sole_wins = [0] * k
            wins = [0] * k
            total_low = 0
            for trial in range(1, TRIALS + 1):
                while True:
                    try:
                        res = icdfa_sparse_test(gen)
                        break
                    except KeyError:
                        continue
                low = min(res)
                total_low += low
                low_count = sum(size == low for size in res)
                for i in range(k):
                    total_raw[i] += res[i]
                    if res[i] == low:
                        sole_wins[i] += low_count == 1
                        wins[i] += 1
                print(trial)
                if not trial % 1000:
                    with open("results.txt", "a") as f:
                        f.write("----------sparse-------------\n")
                        f.write(f"trials={trial}, n={n}, alph={alph}\n")
                        f.write(f"{sole_wins}\n")
                        total_raw.append(total_low)
                        f.write(f"{total_raw}\n")
                        total_raw.pop()
            with open("results.txt", "a") as f:
                f.write("----------sparse-------------\n")
                f.write(f"trials={trial}, n={n}, alph={alph}\n")
                f.write(f"{sole_wins}\n")
                total_raw.append(total_low)
                f.write(f"{[total / total_raw[1] for total in total_raw]}\n")


def icdfa_planar_test():
    gen = icdfa.generator(20, 2)
    r = 0
    TRIALS = 8737
    for i in range(1, TRIALS + 1):
        edges, n = next(gen)
        r += is_planar(edges, {0, n - 1})
        print(i / TRIALS)
    print("DONE!")
    print(r / TRIALS)


def weight_vs_gate_test(gen):
    edges, n = next(gen)
    wfa = WFA(n, edges, require_trim=False)
    weight = wfa.totalWidth(w=True)
    gate = wfa.totalWidth(gate=True)
    return weight, gate


def weight_vs_gate_main(ns):
    TRIALS = 10000
    for n in ns:
        res = []
        for alph in range(2, 51):
            gen = icdfa.generator(n, alph)
            lowest, highest = math.inf, -math.inf
            total_weight = total_gate = 0
            gate_wins = ties = 0
            total_best = 0
            for trial in range(1, TRIALS + 1):
                weight, gate = weight_vs_gate_test(gen)
                lowest = min(lowest, gate / weight)
                highest = max(highest, gate / weight)
                total_weight += weight
                total_gate += gate
                total_best += min(weight, gate)
                if weight > gate:
                    gate_wins += 1
                elif weight == gate:
                    ties += 1
                print(f"{alph}.{trial}")
            print(f"DONE alph={alph}")
            print(lowest, highest, total_gate / total_weight, gate_wins, ties)
            res.append(
                # (lowest, highest, total_gate / total_weight, total_gate, total_weight, total_best, gate_wins, ties)
                total_gate / total_weight
            )
            with open("results.txt", "a") as f:
                f.write("----------weight-vs-gate-v2----------\n")
                f.write(f"alph max = {alph}\n")
                f.write(f"n={n}\n")
                f.write(f"{res}\n")
        with open("results.txt", "a") as f:
            f.write("----------weight-vs-gate-v2----------\n")
            f.write(f"n={n}\n")
            f.write(f"{res}\n")


def sep_test(ns):
    for n in ns:
        for alph in (2, 5, n // 2, n):
            gen = icdfa.generator(n, alph)
            TRIALS = 10000
            totalS = 0
            totalAbs = 0
            for i in range(1, TRIALS + 1):
                edges, n = next(gen)
                B, S, C = seperator(n, edges)
                totalS += len(S)
                totalAbs += abs(len(B) - len(C))
                print(f"trial {i}/{TRIALS}")
            with open("results.txt", "a") as f:
                f.write("----------sep_test----------\n")
                f.write(f"n={n}, alph={alph}\n")
                f.write(f"{totalS / TRIALS}, {totalAbs / TRIALS}\n")


if __name__ == "__main__":
    # Section 5 gate score motivation example
    # testGateMotivation(50, 10)

    # Table 1
    # sep_test((20, 40, 50))

    # # Tables 2 and 3
    icdfa_main([20, 40, 50])
    icdfa_sparse_main([20, 40, 50])

    # # Additional weight vs. gate score experiments
    # weight_vs_gate_main([10, 20, 25, 30, 35, 40, 50])

    # # Table 4
    # eliminate_most_main((20, 30, 50, 100, 200, 500))
