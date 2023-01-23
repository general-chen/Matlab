## This code plots a non-slender delta wing in matlab, and plots the pressure contour on the surface of the delta wing.
  - first, we need the airfoil interpolation function curve for cross section, in this demo, the airfoil is NACA0012.
  - for different cross-sections, we need to scale the airfoil appropriately to ensure that the model is thick in the middle and thin at the wingtips.
  - once we plot each airfoil of cross sections (depend on the accuracy of the model you need), we can combine them together in matlab as a surface:
  - ![Aaron Swartz](Geometry_contour_animation/delta_wing.png)
