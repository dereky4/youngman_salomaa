from gurobipy import Model, GRB, quicksum
import random
from main import WFA
import itertools


def get_order(n, edgeMat):
    model = Model("order")
    x = [[None] * (n - 1) for _ in range(n - 1)]
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            x[i][j] = model.addVar(vtype=GRB.BINARY, name=f"x_{(i,j)}")
    for i in range(1, n - 1):
        model.addConstr(quicksum(x[i][j] for j in range(1, n - 1)) == 1)
    for j in range(1, n - 1):
        model.addConstr(quicksum(x[i][j] for i in range(1, n - 1)) == 1)
    model.setObjective(
        quicksum(
            quicksum(
                quicksum(
                    quicksum(
                        x[je][j] * x[ke][k] * edgeMat[j][k] * 2 ** (-(je + ke))
                        for ke in range(1, n - 1)
                    )
                    for k in range(1, n - 1)
                )
                + x[je][j] * (edgeMat[0][j] / 3) * 2 ** (-(n + je - 3))
                + x[je][j] * (edgeMat[j][0] + edgeMat[j][n - 1]) * 2 ** (-(n + je - 1))
                + x[je][j] * (edgeMat[n - 1][j] / 3) * 2 ** (-(n + je - 2))
                for je in range(1, n - 1)
            )
            for j in range(1, n - 1)
        ),
        GRB.MINIMIZE,
    )
    model.setParam("OutputFlag", 0)
    model.optimize()
    print(model.Status == GRB.OPTIMAL)
    order = [None] * (n - 2)
    for i in range(1, n - 1):
        for j in range(1, n - 1):
            if x[i][j].x:
                order[i - 1] = j
    return order


if __name__ == "__main__":
    n = 15
    edgeMat = [[random.randint(0, 1000) for _ in range(n)] for _ in range(n)]
    edges = [(j, k, edgeMat[j][k]) for j in range(n) for k in range(n)]
    order = get_order(n, edgeMat)
    wfa = WFA(n, edges)
    w = wfa.totalWidth(w=True)
    opt = wfa.totalWidth(order)
    orders = list(itertools.permutations(order))
    for order in orders:
        assert wfa.totalWidth(order) >= opt
    print(w / opt)
