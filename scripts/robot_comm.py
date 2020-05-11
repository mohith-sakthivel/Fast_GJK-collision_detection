# -*- coding: utf-8 -*-
"""
    This module contains functions for communication with
    ABB robots over TCP/IP
"""

import numpy as np
import socket


def GetData(RobAddr, Pose, sock, close_comm=False):
    """
    Recieves position data from the robot as a string
    """
    try:
        sock.send("Give Pose".encode('utf-8'))
    except OSError:
        sock = socket.socket()
        sock.connect((RobAddr[0], RobAddr[1]))
        sock.send("Give Pose".encode('utf-8'))
    except ConnectionRefusedError:
        print('Connection Refused by server')
        return False
    Data = []
    Data.append((sock.recv(1024)).decode())
    if close_comm:
        sock.close()
    index = np.char.find(Data[0], ']')
    Pos = np.fromstring(Data[0][1:index], dtype='float64', sep=',')
    Orient = np.fromstring(Data[0][(index+3):-1], dtype='float64', sep=',')
    Track = float(Data[0][(np.char.find(Data[0], '&')+1):])
    Pose.append(Pos)
    Pose.append(Orient)
    Pose.append(Track)
    return True


def RobotStop(RobAddr, sock, close_comm=False):
    """
    Sends a stop message to the robot
    """
    try:
        sock.send("Stop".encode('utf-8'))
    except OSError:
        sock = socket.socket()
        sock.connect((RobAddr[0], RobAddr[1]))
        sock.send("Stop".encode('utf-8'))
    except ConnectionRefusedError:
        return False
    if close_comm:
        sock.close()
    return True
