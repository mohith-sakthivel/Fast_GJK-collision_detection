# Collision Detection using GJK Algorithm

## Description:
This is a collision detection software based on the GJK Algorithm to detect collision between convex polytopes. It is based on a fast and robust implementation of the GJK algorithm by M Montanari et al. ([Paper](https://dl.acm.org/doi/10.1145/3072959.3083724)) 

This algorithm can be used for applications like
- collision prediction for robotic systems
- collision checking in simulations
- Video games

This repo includes pre-built classes only for convex polytopes. The same code can be extended to other shapes. Classes for new shapes should include variables to store the shape, and fuctions for updating the shape position in 3D space and the corresponding support function for the shape. Note that the algorithm is valid only for convex shapes.

## Example:
This repo includes an example program in which the collision detection algorithm is used for detecting collision between two ABB IRB6620 robots. The collision is moitored between the robot end effectors and the adjacent walls. All the objects in this environment are modeled as convex polytopes.

## Additional:
The module rotation_ops.py contains many useful helper functions for performing 3D transformations.
