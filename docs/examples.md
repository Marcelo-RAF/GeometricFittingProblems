---
layout: default
title: Generate Examples
nav_order: 4
---

To generate examples, we will use the `build_problem` function, which can create artificial examples of 2D spheres, 3D spheres, 2D lines, 3D lines, planes, and circles. The input parameters for this function depend on the object to be fitted and can optionally include noise in the points. The function will generate a file with the `csv` extension, which will contain a matrix with the problem's points, the number of points, the number of outliers, and whether or not it includes noise.

# 2D-Sphere

To create a 2D sphere example, use the `build_problem` function with the string `"sphere2D"` and a vector containing:

- The center coordinates
- The radius
- The number of points
- The number of outliers
- A boolean indicating whether to include noise

**Example:**
For a sphere with center **c = (2.0, 3.0)**, radius **9**, containing **23** points, **2** outliers, and including noise:

```julia
build_problem("sphere2D", [2.0, 3.0, 9.0, 23.0, 2.0], true)
```

# 3D-Sphere

For 3D spheres, the process is similar to 2D spheres, with an additional center coordinate.

**Example:**
To generate a sphere with center **c = (-3.0, 2.0, 1.0)**, radius **31**, containing **80** points, **0** outliers, and no noise:

```julia
build_problem("sphere3D", [-3.0, 2.0, 1.0, 31.0, 80.0, 0.0], false)
```

# Planes

To create examples of planes, provide:

- A point on the plane
- Two non-collinear direction vectors
- The number of points
- The number of outliers
- A boolean indicating whether to include noise

**Example:**
For a plane passing through point **p = (-5.0, 2.0, -1.0)** with directions **v₁ = (2.0, 3.0, -1.0)** and **v₂ = (1.0, -3.0, 1.0)**, containing **120** points, **12** outliers, and including noise:

```julia
build_problem("plane", [-5.0, 2.0, -1.0, 2.0, 3.0, -1.0, 1.0, -3.0, 1.0, 120, 12], true)
```

# 2D Lines

For 2D lines, specify:

- The slope
- The intercept
- The number of points
- The number of outliers
- A boolean indicating whether to include noise

**Example:**
For a line with slope **3**, intercept **-8**, containing **14** points, **1** outlier, and no noise:

```julia
build_problem("line2D", [3.0, -8.0, 14.0, 1.0], false)
```

# 3D Lines

To generate a 3D line, provide:

- A point on the line
- A direction vector
- The number of points
- The number of outliers
- A boolean indicating whether to include noise

**Example:**
For a line passing through point **p = (-2.0, 0.0, 1.0)** with direction vector **v = (1.0, 2.0, 1.0)**, containing **50** points, **12** outliers, and including noise:

```julia
build_problem("line3D", [-2.0, 0.0, 1.0, 1.0, 2.0, 1.0, 50.0, 12], true)
```

# Circles

To create circles, specify:

- The center
- The radius
- Two directions for the plane containing the circle
- The number of points
- The number of outliers
- A boolean indicating whether to include noise

**Example:**
For a circle centered at **c = (-5.0, -2.0, 15.0)**, radius **29**, directions **v₁ = (2.0, 1.0, 3.0)** and **v₂ = (-1.0, 2.0, 2.0)**, containing **25** points, **2** outliers, and no noise:

```julia
build_problem("circle", [-5.0, -2.0, 15.0, 29.0, 2.0, 1.0, 3.0, -1.0, 2.0, 2.0, 25.0, 2.0], false)
```

Para carregar os exemplos criados, basta digitar

# Learn More

For more details on how to use the `build_problem` function and how to fit geometric objects using Geometric Algebra algorithms, visit the [How to Use](how_to_use.md) guide.
