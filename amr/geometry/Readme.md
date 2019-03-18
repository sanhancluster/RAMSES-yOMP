# Fortran-geometry

This is a very specific and very limited module in Fortran to easily volumetric properties of geometrical objects.
For now it supports:
* [x] cube
* [x] cylinder
* [ ] box (easy to implement, just not done)
* [ ] sphere (should be easy too)

## How to use it?

The soft is presented as a module that you can `use` (`use CylinderGeometry`). You can then create geometrical objects. Properties of these objects are accessed using the Fortran 2008 Object-Oriented syntax, e.g.:
```fortran
program test
  use CylinderGeometry

  type(Cylinder_t) :: cyl
  type(Cube_t) :: cube

  real(8) :: volume
  logical :: intersect

  ! Define a cylinder
  cyl%center = (/1., 1., 1./)
  cyl%r = 1.0  ! This is the radius
  cyl%h = 1.0  ! This is the height
  cyl%normal1 = (/1., 1., 1./)

  ! Define a cube
  cube%center = (/10., 8., 0./)
  cube%a = 1.5 ! This is the width of the cube

  ! Let's get the volume
  print*, 'Volume of the cube:', cube%volume()

  ! Let's get the volume that intersect between the cube and the cylinder
  call cyl%intersectionVolume(cube, volume, intersect)

  print*, 'Cube and volume intersect?', intersect
  print*, 'Volume of intersection   :', volume
end program test
```

## A note on the "algorithm"

There are two algorithms implemented to compute the intersection volume. The first one is using a Monte-Carlo based approach, where random points are drawn from within the smallest volume. The ratio of points also falling in the other one to the total number of points gives the relative volume of the intersection. The convergence is in `1/sqrt(N)`, where `N` is the number of draws. To use it, use the method `intersectionMethodMC`.

The other algorithm is a recursive one. At each step, it divides the box in two equal-sized boxes along one axis and calls itself recursively for the two subboxes. At each depth, the axis is changed (first split along the x axis, then along the y axis, then along the z axis). The algorithm stops at the maximal depth or when all the points of the box are within the cylinder. Once the maximal depth is reached, we approximate the volume of the box in the cylinder by the number of corners in the cylinders / the number of corners (8).

## Need some more examples?

You can read the test or send me an email to get more information.

## Licensing
The code follows the MIT license.
