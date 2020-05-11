"""
    This module contains the implemntations of the geometries
    which can be used to model the objects for monitoring collisions

    Contains class for convex polytopes

    Other geometries can be added as well and be used with the GJK algorithm
"""
import numpy as np
from rotation_ops import RotatePos
# np.random.seed(0)


class Polytope():
    def __init__(self, PointList):
        """
        Create and intialize a convex ploytope at the robot
        end effector
        """
        self.Vertices = np.array(PointList, dtype='float64')
        assert len(self.Vertices.shape) == 2, \
            "Atleast 1 point must be specified.\n" +\
            "All points must be in 3D space"
        assert self.Vertices.shape[1] == 3, "All points must be in 3D space"
        self.NoOfVertices = len(PointList)
        self.CurrentVertices = np.copy(self.Vertices)

    def UpdatePosition(self, pos, quat):
        """
        Updates the vertices of the polytope based on the
        provided position vector and quaternion
        """
        pos = np.array(pos, dtype="float64")
        quat = np.array(quat, dtype="float64")
        for i in range(len(self.Vertices)):
            self.CurrentVertices[i] = pos+RotatePos(self.Vertices[i], quat)

    def SupportFunc(self, dVector):
        """
        Provides the farthest point in the polytope along the
        provided direction dVector

        Inputs:
            dVector      - direction vector in 3D space
        Output:
            Farthest point in polytope along the direction dVector,
            Projection of the farthest point on the direction vector
        """
        MaxDist = np.dot(dVector, self.CurrentVertices[0])
        VertexIndex = 0
        for i in range(1, len(self.CurrentVertices)):
            dist = np.dot(dVector, self.CurrentVertices[i])
            if dist > MaxDist:
                MaxDist = dist
                VertexIndex = i
        return np.array([self.CurrentVertices[VertexIndex],
                         MaxDist, VertexIndex])
