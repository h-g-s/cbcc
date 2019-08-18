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

## queens

Solves the n-queens problem using a binary programming formulation.

Example of usage:

```console
$ ./queens 10
```

## queens-lazy

Solves the n-queens problem using a binary programming formulation. Constraints
to forbid conflicts in the two diagonals are added as lazy conntraints in the 
cut generator callback.

Example of usage:

```console
$ ./queens-lazy 10
```


