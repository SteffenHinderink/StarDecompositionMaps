# Star Decomposition Maps

<img src=sdm-image.jpeg width=375 align=right>

This is an implementation of the volumetric mapping method presented in the SIGGRAPH Asia 2024 paper
"[**Bijective Volumetric Mapping via Star Decomposition**](http://graphics.cs.uos.de/papers/sdm-preprint.pdf)".
Previous restrictions to star-shaped or convex shapes are lifted.
The central idea is a systematic decomposition into simpler subproblems.

The implementation covers the steps:

1. Target decomposition into star-shaped parts
2. Compatible source decomposition
3. Boundary map extension onto interior part boundaries

It outputs the subproblems which can subsequently be handled by bijective mapping methods that require star-shaped targets, e.g.
[GalaxyMaps](https://github.com/SteffenHinderink/GalaxyMaps).

## Citation

If you find this code useful, please consider citing our paper:

```bibtex
@article{StarDecompositionMaps,
  author  = {Hinderink, Steffen and Br{\"u}ckler, Hendrik and Campen, Marcel},
  title   = {Bijective Volumetric Mapping via Star Decomposition},
  year    = {2024},
  journal = {ACM Trans. Graph.},
  volume  = {43},
  number  = {6},
  pages   = {x:1--x:11}
}
```

## Prerequisites

The code can be used as either a **library** or an **executable**.

In either case, the following libraries need to be installed:

- [Eigen](https://eigen.tuxfamily.org)
- [CGAL](https://www.cgal.org)
- [GMP](https://gmplib.org)

The following libraries are included in the [ext](ext) directory:

- [OpenVolumeMesh](https://gitlab.vci.rwth-aachen.de:9000/OpenVolumeMesh/OpenVolumeMesh.git)
- [OpenMesh](https://gitlab.vci.rwth-aachen.de:9000/OpenMesh/OpenMesh.git)
- [Exact predicates](https://www.cs.cmu.edu/~quake/robust.html)
- [TetGen](http://www.tetgen.org)

For the viewer:

- [Dear ImGui](https://github.com/ocornut/imgui.git)
  - [Glad](https://glad.dav1d.de)
  - [GLFW](https://github.com/glfw/glfw.git)
- [stb](https://github.com/nothings/stb)
- [tinyfiledialogs](https://sourceforge.net/projects/tinyfiledialogs/)

Remember to clone this repository recursively:\
```git clone --recursive https://github.com/SteffenHinderink/StarDecompositionMaps.git```

## Library

### CMake

The library can be included via CMake like this:

```cmake
add_subdirectory(path/to/StarDecompositionMaps)
target_link_libraries(yourTarget StarDecompositionMaps)
```

### Usage

The function ```sdm``` implements the method of the paper.
Input are the tetrahedral meshes $M$ and $N$ represented using
[OpenVolumeMesh](https://www.graphics.rwth-aachen.de/software/openvolumemesh/).
They are required to have no degenerate or inverted tetrahedra.
The function ```retetrahedrize``` using
[TetGen](http://www.tetgen.org)
can be used to achieve this.
To represent boundary map $\psi$, corresponding boundary vertices of $M$ and $N$ must have the same indices.
Output are the compatible pairs of meshes $M_i$ and $N_i$.
Their vertices are available in rational numbers in the property ```Q``` of type ```Vector3q```.
The meshes $N_i$ are star-shaped such that a method like
[GalaxyMaps](https://github.com/SteffenHinderink/GalaxyMaps)
can be used subsequently.

The algorithms for the different steps of the methods can also be used independently.
- The function ```decompose``` implements the decomposition of a mesh into star-shaped parts.
- The function ```cut``` finds a surface in a mesh given a boundary loop.
  It includes the surface shift algorithm.
- The function ```align``` implements the recursive triangulation alignment algorithm.

```cpp
#include <retet.h>
#include <sdm.h>

Mesh M = ...;
Mesh N = ...;

retetrahedrize(N); // optional, in case N has degenerate or inverted tetrahedra

std::vector<std::pair<Mesh, Mesh>> components = sdm(M, N);

for (int i = 0; i < components.size(); i++) {
    Mesh& Mi = components[i].first;
    Mesh& Ni = components[i].second; // star-shaped

    auto Q = Ni.property<Vertex, Vector3q>("Q");
    for (auto v : Ni.vertices()) {
        std::cout << "(" << Q[v].transpose() << ")" << std::endl;
    }
}
```

## Executable

### Building

- ```mkdir build```
- ```cd build```
- ```cmake [-DGUI=OFF] ..``` (The option ```GUI``` controls if the program is built with the viewer)
- ```make -j```

### Usage

```./sdm <M> <N> [-o <out_dir>]```

- ```<M>```:
Tetrahedral mesh $M$.
- ```<N>```:
Tetrahedral mesh $N$.
This may contain degenerate or inverted tetrahedra.
In that case,
[TetGen](http://www.tetgen.org)
is used to create a new tetrahedrization.
Corresponding boundary vertices of $M$ and $N$ must have the same indices.
- ```-o <out_dir>```:
Optional argument to output the mesh pairs into the directory ```out_dir```.
Meshes are output in ```.ovm``` format
and have the property ```Q_string``` of type ```std::string```
containing the rational vertices as strings.

e.g. ```./sdm ../meshes/open.vtk ../meshes/thumb.vtk```

The viewer shows the components in different colors.
Pressing <kbd>P</kbd> switches between $M$ and $N$.
