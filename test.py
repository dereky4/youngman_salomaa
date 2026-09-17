from main import WFA, NFA, RE


n = 6
edges = [
    (5, 5, 1),
    (0, 1, 0),
    (1, 2, 1),
    (2, 1, 1),
    (1, 3, 1),
    (2, 3, 0),
    (3, 3, 1),
    (3, 4, 1),
    (3, 5, 0),
    (4, 3, 1),
    (4, 2, 1),
    (4, 5, 0),
]
edgesc = [(j, k, 'a' if w else RE.EPSILON) for j, k, w in edges]
wfa = WFA(n, edges)
nfa = NFA(n, edgesc)
print(wfa.totalWidth())
print(nfa.getRE().size())
