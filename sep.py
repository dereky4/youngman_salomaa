from gurobipy import Model, GRB, quicksum
import random


def seperator(n, edges, sub=None):
    # linear program for finding small seperators
    if sub is None:
        sub = {i for i in range(1, n - 1)}
    else:
        sub = set(sub)
    model = Model("Vertex Seperator")
    B = {}
    C = {}
    for i in sub:
        B[i] = model.addVar(vtype=GRB.BINARY, name=f"B_{i}")
        C[i] = model.addVar(vtype=GRB.BINARY, name=f"C_{i}")
    for i, j, _ in edges:
        if i != j and i in sub and j in sub:
            model.addConstr(B[i] + C[j] <= 1)
            model.addConstr(C[i] + B[j] <= 1)
    for i in sub:
        model.addConstr(B[i] + C[i] <= 1)
    z = model.addVar(vtype=GRB.CONTINUOUS, name="z")
    model.addConstr(z <= quicksum(B[i] for i in sub))
    model.addConstr(z <= quicksum(C[i] for i in sub))
    model.setObjective(
        z, GRB.MAXIMIZE
    )
    model.setParam("OutputFlag", 0)
    model.optimize()
    return (
        [i for i in sub if B[i].x],
        [i for i in sub if not B[i].x and not C[i].x],
        [i for i in sub if C[i].x],
    )


def main():
    n = 50
    m = n // 2
    edges = {
        (i, j, 1)
        for i in range(n)
        for j in range(n)
        if random.random() < 0.1
    }

    def hasEdge(i, j, edges):
        return i != j and ((i, j) in edges or (j, i) in edges)

    B, S, C = seperator(n, edges)
    for i in B:
        for j in C:
            assert not hasEdge(i, j, edges)
    print(len(B), len(S), len(C))


if __name__ == "__main__":
    num = 3.14159
    rounded_num = round(num, 2)
    print(rounded_num)  # Output: 3.14

