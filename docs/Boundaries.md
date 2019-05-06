Layer boundaries can be represented as heightfields.

This should be an optional feature since it will cause (much?) longer simulation times.

### 1D vs. 2D heightfields

The heightfields can be 1D and assumed to be radially symmetric. Like the original MCML, I record photon distribution in symmetric fashion, i.e. only as function of radius and depth resp. transmission angle. Nevertheless, 2D heightfields that are not symmetric can produce different results. Imagine there are 3 boundaries, then the exact scattering location on the middle boundary can determine the radius of transmission on the top or bottom boundary.

Therefore 2D heightfields would be nice to have, but then also a real 2D output is desired, otherwise you could not tell which brightness a particular point on the heightfield has. I start with implementing 1D heightfields.

### Input format

As of the latest update, both types of boundaries (flat and heightfield) are allowed, choosable on per-boundary basis. Explicit heightfield boundaries can also be turned on or off globally for performance reasons.

**Update 1:** The original mci file format must be adapted because it assumed boundaries to be implicitely given from the layer thickness d.

Instead of

```
1                                               # Number of layers
#n	mua	mus	g	d               # One line for each layer
1                                               # n for medium above
1.4	0.1	90	0	100             # layer 1
1                                               # n for medium below
```

I need

```
1                                               # Number of layers
#n	mua	mus	g	d               # One line for each layer
1                                               # n for medium above
0.0 0.5 0.0 0.5 0.0 0.5 0.0 0.5 0.0 0.5         # Upper boundary as heighfield with 10 samples
1.4	0.1	90	0	100             # layer 1
0.0 0.5 0.0 0.5 0.0 0.5 0.0 0.5 0.0 0.5         # Lower boundary as heighfield with 10 samples
1                                               # n for medium below
```
or something similar.

To simplify this, I will use the same number of samples for all heightfields and a regular radial grid with its center at xy==(0,0), represented as an 1D array.

**Update 2:** To enable irregular sampling I use explicit spacing values after each height value (here with 0.1 spacing everywhere):

```
1                                               # Number of layers
#n	mua	mus	g	d               # One line for each layer
1                                               # n for medium above
0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1         # Upper boundary as heighfield with 10 samples
1.4	0.1	90	0	100             # layer 1
0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1 0.0 0.1 0.5 0.1         # Lower boundary as heighfield with 10 samples
1                                               # n for medium below
```


Furthermore, lines containing boundary data should be marked by a special character so that they can be skipped in case the user wants to run a simulation with implicit flat boundaries.

**Update 3:** To enable mixing implicit and explicit boundaries in the same simulation, I have extended the format further. Lines containing explicit boundary data start with 'b', then number of samples and then the list of sample values (alternating heights and spacings). Example:

```
2                        	# Number of layers
#n	mua	mus	g	d         	# One line for each layer
1                         	# n for medium above
b 50 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 # explicit heightfield boundary with 50 samples
1.4	0.1	90	0	0.5    	# layer 1
# implicit flat boundary
1.4	0.1	90	0	0.5    	# layer 2
b 10 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 0.0 0.1 # explicit heightfield boundary with 10 samples
1                        	# n for medium below
```

To be able to mix implicit and explicit boundaries, the internal Boundary object should have a Boolean indicating if it is a heightfield and a heightfield that is only filled with data, if the Boolean is true. Otherwise a single depth value indicates the position of the implicit flat boundary.

When the kernel define IGNORE_HEIGHTFIELDS is set, any explicit boundary data from the file is 

To support different numbers of samples per boundary, I use dynamic memory allocation during reading the file, but the resulting non-contiguous data structure cannot be used uploaded and used as is in GPU memory. Therefore there is another Boundary structure, that points to index ranges in separate buffers containg the actual height samples and spacing values for all boundaries.

**TODO:**
- Support file paths instead of reading sample data directly, to prevent unreadable input files (for humans). Our custom MCI format can even be backwards compatible when we write the boundary info into (special) comments. Unfortunately CUDAMCMLio does not skip comments as it should.

### Output data (RAT detection arrays)

The arrays are already clamped, so that all weights that fall outside the array bounds are accumulated at the edges. In original MCML this could only happen in xy directions.

Now also the z-index of the A detection array must be clamped, as we can have for example weight drops at z < 0. Note that heights are interpreted to go from heightfield origin into -z direction, because +z goes "down" into the material.

### Layer model validity check

One must check if heightfields overlap each other before starting the simulation. Overlapping heightfields would constitute an invalid layer model.

### Extrapolating the heightfield

Because layers are assumed to stretch infinitely in xy directions, there must be an edge case handling for heightfields. For an outside point, one could simply take the nearest heightfield value, thus continuing it infinitely.

### Interpolating the heightfield and finding normals

To obtain a height value at arbitrary floating point coordinates, the heightfield must be interpolated. Therefore it would make sense to store them in textures and use built-in access operators on the GPU. Currently, the program uses raw buffers and does interpolation in code.

To calculate reflections and transmissions at heightfield boundaries, a normal will be needed. The normal can be calculated from the gradient using forward-, backward- or central difference methods.

When interpolating linearly, there will be sharp local minima/maxima where the gradient cannot be defined. For these cases I treat the boundary as a xy-plane with normal along z-axis.

The direction of the normal is always chosen to point into the half space from where an intersecting ray comes in.

### Intersection testing heightfields

Some methods for intersecting a ray with a heightfield need to be investigated.

First approach: The problem is usually not analytically solved (in real-time apps). Instead, a ray-marching with some step size is performed and as soon as the heightfield is crossed, a binary search refines the intersection point. This can be seen in this [figure from a 2005 graphics paper](https://www.researchgate.net/figure/Ray-intersection-with-a-height-field-surface-using-binary-search-Starting-with-A-and-B_fig5_220792017), in the [2006 Parallax Occlusion Mapping paper by Tatarchuk](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.115.3066&rep=rep1&type=pdf) and in a [documentation of the vterrain software](http://vterrain.org/Implementation/Libs/ray-heighhtfield.html). In the case of a monte carlo simulation, one must not step further than the photon step length. For points along the ray, the corresponding point on the heightfield is evaluated. Finally, of course, as a stop criterion for the binary search a distance epsilon must be defined.

Second approach: The above methods are not suited for MCML. Instead an exact analytical intersection with conic-section will be implemented. The cone shape arises when radial 1D heightfields are used. Before the actual intersection test, one must find all potential cones in the path of the ray (fortunately we have a relatively small step length s), which can be done by looking only at the xy plane. Then we apply the [analytic method from the Geometric Tools Website](https://www.geometrictools.com/Source/Intersection3D.html#LinearVolumetric) to find the intersection point. [A debug visualization has been done in this Shadertoy](https://www.shadertoy.com/view/ltdfDs).


Intersection algorithm outline:
- Iterate over all cones in the path of the line and take the closest hit.
- To make sure not to miss any cones, get cone under the closest point on the line to heightfield center.
- By assigning an ascending index from inner to outer cones, all in-between cones can easily be iterated.

To test the implementation, following ideas came up:
- Shoot many random rays at the heightfield
- Visualize results in a pointcloud renderer to see the radial heightfield shape
- This can also be tested for multiple layers by saving points but following the rays further

### Regular vs. irregular grids

The heightfield can be regularly or irregularly sampled. This will determine the memory format to choose. Regular heightfields only need one spacing value (or alternatively one width value) and an array of heights. To support irregular heightfields, we have a spacing value for each height value (see file format above). For processing, the following data structure is used:

```
/**
* Radial heightfield
* Holds center, but not actual sample values,
* because they are dynamic. Instead it points
* to ranges in arrays of heights and spacings.
*/
struct RHeightfield {
    Real3 center;

    uint i_heights;
    uint n_heights;

    uint i_spacings;
    uint n_spacings;
};
```

### Heightfields represented as functions

A Heightfield can also be represented by a function which takes an x (and y) value and calculates some height based on some formula. This would also be a nice to have feature.

### Future work

Finally it should be possible to use arbitrary 3D meshes as boundaries.