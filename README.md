# Incremental cut procedure for Densest Subgraph Problem

This repository contains an implementation of the incremental cut procedure for the Densest Subgraph Problem

# Compiling

A makefile is provided. To compile the `incremental` executable, use:

```
mkdir bin
make
```

# Usage

The executable reads the input file from stdin. To execute one of the examples, use:
```
./bin/incremental < examples/example_unweighted.txt
```

# Options

`--help` : Print a help message

`--accuracy N` : Set the accuracy of the density to N decimal places. Note that a value of N too high may produce errors, especially on big graphs or graphs with big weights. By default, 4 decimal places are used.

`--injectLambda L` : Start the incremental procedure with density L. Note that a value of L larger than the optimal may produce errors.

`--weightedEdges` : Use this option if the edges have non-unitary benefits.

`--weightedNodes` : Use this option if the nodes have non-unitary costs.

`--dumpSourceSet FILE` : Dump the list of nodes that are in the optimal solution to FILE.

# Input format

The input graph, with N nodes and M edges, is to be provided in the following text format:

```
N M
[WEIGHT_NODE_1]
[WEIGHT_NODE_2]
...
[WEIGHT_NODE_N]
FROM_EDGE_1 TO_EDGE_1 [WEIGHT_EDGE_1]
FROM_EDGE_2 TO_EDGE_2 [WEIGHT_EDGE_2]
...
FROM_EDGE_M TO_EDGE_M [WEIGHT_EDGE_M]
```

Remember that if the graph is weighted, option `--weightedEdges` must be used if the edges have weights (non-unitary benefits) and `--weightedNodes` if the nodes have weights (non-unitary costs).

See directory [examples](examples) for sample input files.

# License

