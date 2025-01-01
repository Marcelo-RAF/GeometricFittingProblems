---
layout: default
title: How to use
nav_order: 3
---

# How to use

This section introduces three fitting algorithms: `AGCGA`, proposed in ..., `AACGA`, proposed in ..., and the algorithm we propose in this work, `ICGA`. Additionally, we present a combined fitting and classification algorithm, `fittingclass`, and a detection algorithm, `LOVO-CGA`.

**Fitting**

To use any of the fitting algorithms, you need a matrix containing a set of points, where each row represents a point in Euclidean space, and the target object to fit specified as a string. Supported objects include: `sphere`, `plane`, `circle`, `line`.

**Example**

Consider a matrix $D$ representing a set of points from a sphere. You can fit the object as follows:

```julia
julia> ICGA(D, "sphere")
```

**Fitting and Classification**

To use the `fittingclass`,  algorithm, provide a matrix containing a set of points corresponding to the object set $Q$ = {hyperspheres, hyperplanes, hypercircles, lines} and two threshold values, $\varepsilon_1$ and $\varepsilon_2$. Para ruídos criados por uma distribuição normal, os testes feitos com $\varepsilon_1 = 1.0$ e $\varepsilon_2 = 10e-3$ se mostraram com bons resultados.

**Example**

Consider a matrix D representing a set of points. You can apply the algorithm as follows:

```julia
julia> fittingclass(D, 1.0, 10e-3)
```

**Detection**

LOVOCGA(data::Matrix, nout::Int, θ::Vector, method::string, object::string, ε=1.0e-4)

Para a detecção de objetos do conjunto $Q$ = {hyperspheres, hyperplanes, hypercircles, lines}, utilizaremos a função `LOVOCGA`. 
Os parâmetros de entrada são: matriz de pontos, quantidade de outliers do problema, chute inicial, o método a ser utilizado e o objeto a ser detectado.

**Example**

Considere uma matriz $D$ representando um conjunto de 20 pontos de um círculo, em que 4 deles são outliers. Podemos então fazer:

```julia
julia> LOVOCGA(D, 4, [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], "AGCGA", "circle", ε=1.0e-4)
```