---
layout: default
title: Generate Examples
---

To generate examples, we will use the `build_problem` function, which can create artificial examples of 2D spheres, 3D spheres, 2D lines, 3D lines, planes, and circles. The input parameters for this function depend on the object to be fitted and can optionally include noise in the points.

## Esfera 2D

To generate an example of a two-dimensional sphere, you should call the `build_problem` function with the string `"sphere2D"` and a vector containing the center, radius, the number of points, the number of outliers in the point set, and a boolean parameter indicating whether the example should include noise or not.

For example, to create an example of a two-dimensional sphere with a center **c = (2.0, 3.0)**, radius **9**, containing **23** points, **2** outliers, and including noise, you would use the following:

```julia
build_problem("sphere2D", [2.0, 3.0, 9.0, 23.0, 2.0], true)
```
