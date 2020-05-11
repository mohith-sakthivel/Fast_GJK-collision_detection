"""
Contains the class implementataion of the GJK Algorithm
"""

import numpy as np
# np.random.seed(0)


class ComputeDist():
    def __init__(self, ObjectList=[], tolerance=0.00001):
        """
        Initialize the cobject environment
        """
        self.ObjCount = 0
        self.ObjID = {}
        self.ConvObject = []
        self.tolerance = tolerance
        for shape in ObjectList:
            assert len(shape) == 2, \
                "Shape name and a shape class object should be provided."
            assert type(shape[0]) == str, "Shape name must be a string"
            self.ObjID[shape[0]] = self.ObjCount
            self.ConvObject.append(shape[1])
            self.ObjCount += 1
        self.WitnessPts = {}

    def SupportMinkowskiDiff(self, Shape1, Shape2, Dir):
        """
        Calculate the support of the Minkowski difference Shape1-Shape2
        """
        Pt1, Dist1, Index1 = Shape1.SupportFunc(Dir)
        Pt2, Dist2, Index2 = Shape2.SupportFunc(-Dir)
        return np.array([Pt1-Pt2, Dist1+Dist2, (Index1, Index2)])

    def S1D(self, Simplex):
        """
        Sub-routine for 1-simplex.
        Searches the voronoi regions of a line.
        """
        if all(Simplex[1] == Simplex[0]):
            return Simplex[1:], np.ones(1, dtype='float64')
        t = Simplex[1]-Simplex[0]
        po = -(np.dot(Simplex[1], t)/np.dot(t, t))*t+Simplex[1]
        u_max = 0
        # Find on which axis to project the simplex to avoid degeneracy
        for i in range(3):
            u = Simplex[0][i]-Simplex[1][i]
            if abs(u) > abs(u_max):
                u_max = u
                Index = i
        # Only i-th co-ordinate is retained
        k = 1
        C2 = np.zeros(2, dtype='float64')
        # Check if the origin is on the line segment
        for j in range(2):
            C2[j] = ((-1)**(j+1))*(Simplex[k][Index]-po[Index])
            k = j
        if (u_max > 0 and all(C2 > 0)) or (u_max < 0 and all(C2 < 0)):
            return Simplex, C2/u_max
        else:
            # Find which end point of the line segment is closest to the origin
            if (np.linalg.norm(Simplex[0]) < np.linalg.norm(Simplex[1])):
                return Simplex[0:1], np.ones(1, dtype='float64')
            else:
                return Simplex[1:], np.ones(1, dtype='float64')

    def S2D(self, Simplex):
        """
        Sub-routine for 2-simplex.
        Searches the voronoi regions of a planar triangle.
        """
        n = np.cross(Simplex[1]-Simplex[0], Simplex[2]-Simplex[0])
        po = (np.dot(Simplex[0], n)/np.dot(n, n))*n

        # Find on which plane to project the simplex to avoid degeneracy
        u_max = 0
        k = 1
        l = 2
        for i in range(0, 3):
            u = ((-1)**(i))*(
                             Simplex[0][k]*Simplex[1][l]
                             + Simplex[1][k]*Simplex[2][l]
                             + Simplex[2][k]*Simplex[0][l]
                             - Simplex[1][k]*Simplex[0][l]
                             - Simplex[2][k]*Simplex[1][l]
                             - Simplex[0][k]*Simplex[2][l])
            if abs(u) > abs(u_max):
                u_max = u
                J = i
            k = l
            l = i
        # co-ordinate J is discarded
        x, y = np.delete(np.arange(0, 3), J)
        k = 1
        l = 2
        C3 = np.zeros(3, dtype='float64')
        # Check if the origin is within the triangle
        for j in range(3):
            C3[j] = (
                     po[x]*Simplex[k][y]+po[y]*Simplex[l][x]
                     + Simplex[k][x]*Simplex[l][y]
                     - po[x]*Simplex[l][y]-po[y]*Simplex[k][x]
                     - Simplex[l][x]*Simplex[k][y])
            k = l
            l = j
        if (u_max > 0 and all(C3 > 0)) or (u_max < 0 and all(C3 < 0)):
            return Simplex, C3/u_max
        d = np.Infinity
        # Find which side of the triangle is closest to the origin
        for j in range(0, 3):
            if (u_max >= 0 and -C3[j] >= 0) or (u_max <= 0 and -C3[j] <= 0):
                Simplex1D = np.delete(Simplex, j, axis=0)
                W_astrix, Lambda_astrix = self.S1D(Simplex1D)
                d_astrix = np.linalg.norm(np.matmul(Lambda_astrix, W_astrix))
                if d_astrix < d:
                    W = W_astrix
                    Lamda = Lambda_astrix
                    d = d_astrix
        return W, Lamda

    def S3D(self, Simplex):
        """
        Sub-routine for 3-simplex.
        Searches the voronoi regions of a tetrahedron.
        """
        M = np.vstack([Simplex.T, np.ones((1, 4), dtype='float64')])
        detM = 0
        C4 = np.zeros(4, dtype='float64')
        # Check if the origin is within the tetrahedron
        for j in range(4):
            C4[j] = ((-1)**((j+1)+4))*np.linalg.det(
                            np.hstack([M[0:3, 0:j], M[0:3, j+1:]]))
            detM += C4[j]
        if (detM > 0 and all(C4 > 0)) or (detM < 0 and all(C4 < 0)):
            return Simplex, C4/detM
        d = np.Infinity
        # Find which face of the tetrahedron is closest to the origin
        for j in range(0, 4):
            if (detM >= 0 and -C4[j] >= 0) or (detM <= 0 and -C4[j] <= 0):
                Simplex2D = np.delete(Simplex, j, axis=0)
                W_astrix, Lambda_astrix = self.S2D(Simplex2D)
                d_astrix = np.linalg.norm(np.matmul(Lambda_astrix, W_astrix))
                if d_astrix < d:
                    W = W_astrix
                    Lamda = Lambda_astrix
                    d = d_astrix
        return W, Lamda

    def SignedVolumes(self, Simplex):
        """
        Performs the signed volumes distance algorithm
        on the simplex specified
        """
        # Call routine based on simplex size
        if len(Simplex) == 4:
            return self.S3D(Simplex)
        elif len(Simplex) == 3:
            return self.S2D(Simplex)
        elif len(Simplex) == 2:
            return self.S1D(Simplex)
        elif len(Simplex) == 1:
            return Simplex, np.ones((1), dtype='float64')
        else:
            print("error")

    def GetDist(self, ShapeID1, ShapeID2):
        """
        Computes the distance between 2 convex bodies
        """
        k = 0
        Shape1 = self.ObjID[ShapeID1]
        Shape2 = self.ObjID[ShapeID2]
        # Initialize the simplex
        # Select a random initial point for the simplex
        index_1 = np.random.randint(0, self.ConvObject[Shape1].NoOfVertices)
        index_2 = np.random.randint(0, self.ConvObject[Shape2].NoOfVertices)
        InitialPt = (
                     self.ConvObject[Shape1].CurrentVertices[index_1] -
                     self.ConvObject[Shape2].CurrentVertices[index_2])

        Simplex = InitialPt.reshape((1, 3))
        # Set initial direction opposite to the selected point to maximize area
        NewPt, _, _ = self.SupportMinkowskiDiff(
                                                self.ConvObject[Shape1],
                                                self.ConvObject[Shape2],
                                                -InitialPt)

        Simplex = np.vstack([NewPt, Simplex])
        while True:
            k += 1
            # Use signed volumes to get the supporting points and their weights
            Simplex, Lambda = self.SignedVolumes(Simplex)
            vk = np.matmul(Lambda, Simplex)
            NewPt, hk, _ = self.SupportMinkowskiDiff(
                                                     self.ConvObject[Shape1],
                                                     self.ConvObject[Shape2],
                                                     -vk)
            vk_square = np.dot(vk, vk)
            gk = vk_square+hk
            if (gk < self.tolerance or len(Simplex) == 4):
                return np.sqrt(vk_square)
            Simplex = np.vstack([NewPt, Simplex])
