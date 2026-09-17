# Experiments

The code reproduces the experimental results reported in Sections 5 and 6 of the paper.

The separator computations use Gurobi, so a Gurobi license is required for the relevant experiments.

## Section 5 example

```python
testGateMotivation(50, 10)
```

Output is printed to the terminal.

## Table 1

```python
sep_test((20, 40, 50))
```

Results are appended to `results.txt` as

```text
----------sep_test----------
n=<n>, alph=<alphabet size>
<average |S|>, <average ||B| - |C||>
```

## Tables 2 and 3

```python
icdfa_main([20, 40, 50])
icdfa_sparse_main([20, 40, 50])
```

Results are appended to `results.txt`.

`icdfa_main` outputs

```text
-----------dense------------
n=<n>, alph=<alphabet size>
[<M0 wins>, <M1 wins>, ..., <M6 wins>]
[<M0/M1>, <M1/M1>, ..., <M6/M1>, <MB/M1>]
```

`icdfa_sparse_main` outputs

```text
----------sparse-------------
trials=10000, n=<n>, alph=<alphabet size>
[<M0 wins>, <M1 wins>, ..., <M5 wins>]
[<M0/M1>, <M1/M1>, ..., <M5/M1>, <MB/M1>]
```

Intermediate results are also written every 1000 trials.

## Weight vs. gate score

```python
weight_vs_gate_main([40, 50])
```

Results are appended to `results3.txt` as

```text
----------weight-vs-gate-v2----------
alph max = <alphabet size>
n=<n>
[<gate/weight ratios up to this alphabet size>]
```

## Table 4

```python
eliminate_most_main((20, 30, 50, 100, 200, 500))
```

Results are appended to `results.txt` as

```text
----------eliminate-most---------
n=<n>, m=<number of eliminated states>, p=0.5
<average rho>, <number of trials with rho = 1>
```