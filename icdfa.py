import FAdo
import FAdo.rndfap
import collections
import random


def generator(n, k):
    # FAdo's default seed=0 triggers time.clock_gettime_ns, which is
    # Unix-only and crashes on Windows. Pass an explicit non-zero seed
    # to skip that code path.
    seed = random.SystemRandom().randrange(1, 2 ** 63)
    gen = FAdo.rndfap.ICDFArgen(n=n, k=k, seed=seed)
    for dfa in gen:
        edgeMaps = [collections.Counter() for _ in range(n)]
        assert len(dfa.States) == n
        assert len(dfa.Sigma) == k
        assert dfa.Initial == 0
        for q0 in dfa.States:
            for a in dfa.Sigma:
                q1 = dfa.Delta(q0, a)
                if q1 is not None:
                    edgeMaps[q0][q1] += 1
        edges = [(q0, q1, w) for q0 in dfa.States for q1, w in edgeMaps[q0].items()]
        yield edges, n


if __name__ == "__main__":
    i = 0
    m = 999 * 999 * 999
    M = -1
    for edges, _ in generator(10, 10):
        if i == 100:
            break
        i += 1
        print(i / 1000)
        edgeCount = sum(w > 0 for _, _, w in edges)
        m = min(m, edgeCount)
        M = max(M, edgeCount)
    print("DONE")
    print(m, M)
