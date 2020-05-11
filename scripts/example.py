# -*- coding: utf-8 -*-
"""
This module contains an example implementation of the  GJK collision
detecton algorithm by modelling theobjects as two cuboids.

This example program was used to monitor the collision between
two PUMA robots each mounted on a linear track
"""
from geometry import Polytope
from rotation_ops import RotateQuat, RotatePos
from compute_distance import ComputeDist
import time
import socket
import numpy as np
from robot_comm import GetData, RobotStop


# Safety Threshold
Track_Limit = 15300
EE_Limit = 15000
Wall_Limit = 100

# Station Positions and Shifts
quat_AB = np.array([0, 0, 0, 1])
Track_Contact_Dist = 15455
X_Shift = 17079
X_Shift_Site_Determined = 17088
X_Shift_RobotStudio = 17079
Rob_Wrist = np.array([0, 0, 200])


# Communication Parameters
Rob1_Addr = ('192.168.131.1', 8888)
Rob2_Addr = ('192.168.131.2', 8889)


# Setup the Ende Effector Polytope
Flange1 = Polytope([[250, 250, 0], [250, -250, 0],
                    [-250, -250, 0], [-250, 250, 0],
                    [250, 250, 10], [250, -250, 10],
                    [-250, -250, 10], [-250, 250, 10]])
Flange2 = Polytope([[250, 250, 0], [250, -250, 0],
                    [-250, -250, 0], [-250, 250, 0],
                    [250, 250, 10], [250, -250, 10],
                    [-250, -250, 10], [-250, 250, 10]])
RobA = ("Rob1", Flange1)
RobB = ("Rob2", Flange2)


# Setup the rearside walls
Wall1 = Polytope([[-1500, 5000, 5000], [-1500, 5000, -5000],
                  [-1500, -5000, -5000], [-1500, -5000, 5000],
                  [-1510, 5000, 5000], [-1510, 5000, -5000],
                  [-1510, -5000, -5000], [-1510, -5000, 5000]])
Wall2 = Polytope([[18230, 5000, 5000], [18230, 5000, -5000],
                  [18230, -5000, -5000], [18230, -5000, 5000],
                  [18250, 5000, 5000], [18250, 5000, -5000],
                  [18250, -5000, -5000], [18250, -5000, 5000]])
Wall_Metal_Door = ("Metal Door Wall", Wall1)
Wall_Concrete = ("Concrete Wall", Wall2)

# Setup the station
Station = ComputeDist([RobA, RobB, Wall_Metal_Door, Wall_Concrete])

# Program  Start
global Rob1_sock
global Rob2_sock
try:
    Rob1_sock = socket.create_connection(Rob1_Addr, 2)
    Rob2_sock = socket.create_connection(Rob2_Addr, 2)
except TimeoutError:
    print('Connnection to server timed out')

count = 0
while True:
    Rob1_Pose = []
    Rob2_Pose = []
    start_time = time.time()

    # Get Data from the robots
    if not GetData(Rob1_Addr, Rob1_Pose, Rob1_sock):
        break
    if not GetData(Rob2_Addr, Rob2_Pose, Rob2_sock):
        break

    Rob1_W_Pose = np.array(Rob1_Pose)
    Rob2_W_Pose = np.array(Rob2_Pose)

    # Transform the robot pose into the world frame
    Rob2_W_Pose[1] = RotateQuat(quat_AB, Rob2_W_Pose[1])
    Rob2_W_Pose[0] = RotatePos(Rob2_W_Pose[0], quat_AB)

    Rob2_W_Pose[0][0] = Rob2_W_Pose[0][0]+X_Shift-Rob2_W_Pose[2]
    Rob1_W_Pose[0][0] = Rob1_W_Pose[0][0]+Rob1_W_Pose[2]

    # Find the position of the Endeffector polytope in world 3D space
    Flange1.UpdatePosition(Rob1_W_Pose[0], Rob1_W_Pose[1])
    Flange2.UpdatePosition(Rob2_W_Pose[0], Rob2_W_Pose[1])

    # Calculate the concerned distances
    EE_Dist = Station.GetDist("Rob1", "Rob2")
    Rob1_Wall_Dist = Station.GetDist("Rob1", "Metal Door Wall")
    Rob2_Wall_Dist = Station.GetDist("Rob2", "Concrete Wall")

    # Check Wrist Distance from Walls
    Rob1_Wrist_Pose = Rob1_W_Pose[0] - RotatePos(Rob_Wrist, Rob1_W_Pose[1])
    Rob2_Wrist_Pose = Rob2_W_Pose[0] - RotatePos(Rob_Wrist, Rob2_W_Pose[1])

    if ((EE_Dist < EE_Limit) or (Rob1_Wall_Dist < Wall_Limit) or
            (Rob2_Wall_Dist < Wall_Limit) or
            (Rob1_W_Pose[2]+Rob2_W_Pose[2] > Track_Limit)):

        Rob1_Stopped = RobotStop(Rob1_Addr, Rob1_sock)
        Rob2_Stopped = RobotStop(Rob2_Addr, Rob2_sock)
        print('Collision Predicted')
        if not (Rob1_Stopped or Rob2_Stopped):
            print('Both Robots did not stop')
        elif not Rob1_Stopped:
            print('Robot 1 did not stop')
        if not Rob2_Stopped:
            print('Robot 2 did not stop')
        else:
            print('Robots stopped')
        break

    end_time = time.time()
    if count % 5 == 0:
        print('Distance: ' + str(EE_Dist))
        print('Robo1 1 and Wall Distance: ' + str(Rob1_Wall_Dist))
        print('Robo1 2 and Wall Distance: ' + str(Rob2_Wall_Dist))
        print('Robot 1 Wrist: '+str(np.round(Rob1_Wrist_Pose, 2)))
        print('Robot 1 Wrist: '+str(np.round(Rob1_Wrist_Pose, 2)))
        print('Time: ' + str(end_time-start_time), end='\n\n')

    count += 1

    while (time.time()-start_time) < 0.25:
        pass

print('Connection Closed')
