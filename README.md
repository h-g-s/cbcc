# cbcc
COIN-OR Branch&amp;Cut Cbc examples using the C API

## Usage

After installing the [COIN-OR CBC MIP Solver](https://github.com/coin-or/Cbc), to build  just run:

```console
$ make
```

## tsp-compact

This example solves the Traveling Salesman Problem using a compact (weak) MIP
formulation.

Example of usage:

```console
$ ./tsp-compact data/ulysses22.tsp
```

