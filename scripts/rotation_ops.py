# -*- coding: utf-8 -*-
# This module contains implementaions of popular 3D rotation operations
import numpy as np
from math import cos, sin, pi, sqrt


def Transf_Matrix(pos, quat):
    """
    Converts a position and quaternion value into
    a (4x4) Transformation matrix

    Arguments:
        pos - array (x, y, z)
        quat - array (q0, qx, qy, qz)
    """
    Result = np.zeros((4, 4))
    Result[3][3] = 1
    Result[0:3, 3] = pos
    Result[0:3, 0:3] = quat_to_homo(quat)
    return Result


def quat_to_homo(quat):
    """
    Converts a quaternion into a homogeneous matrix

    Arguments:
        quat - array (q0, qx, qy, qz)
    """
    ResultMatrix = np.ones((3, 3), dtype="float64")
    ResultMatrix[0, 0] = 1-2*(quat[2]**2+quat[3]**2)
    ResultMatrix[0, 1] = 2*(quat[1]*quat[2]-quat[3]*quat[0])
    ResultMatrix[0, 2] = 2*(quat[1]*quat[3]+quat[2]*quat[0])
    ResultMatrix[1, 0] = 2*(quat[1]*quat[2]+quat[3]*quat[0])
    ResultMatrix[1, 1] = 1-2*(quat[1]**2+quat[3]**2)
    ResultMatrix[1, 2] = 2*(quat[2]*quat[3]-quat[1]*quat[0])
    ResultMatrix[2, 0] = 2*(quat[1]*quat[3]-quat[2]*quat[0])
    ResultMatrix[2, 1] = 2*(quat[2]*quat[3]+quat[1]*quat[0])
    ResultMatrix[2, 2] = 1-2*(quat[1]**2+quat[2]**2)
    return ResultMatrix


def eulerzxz_to_quat(euler_z_x_z):
    """
    Converts a Euler ZXZ angles into quaternion values

    Arguments:
        euler_z_x_z - array (eu_z1, eu_x, eu_z2)
    """
    euz = euler_z_x_z[0]*pi/180
    eux_d = euler_z_x_z[1]*pi/180
    euz_dd = euler_z_x_z[2]*pi/180
    cos_euz = cos(euz)
    sin_euz = sin(euz)
    cos_eux_d = cos(eux_d)
    sin_eux_d = sin(eux_d)
    cos_euz_dd = cos(euz_dd)
    sin_euz_dd = sin(euz_dd)

    m_11 = cos_euz*cos_euz_dd-cos_eux_d*sin_euz*sin_euz_dd
    m_12 = cos_euz_dd*sin_euz+cos_eux_d*cos_euz*sin_euz_dd
    m_13 = sin_eux_d*sin_euz_dd
    m_21 = -cos_euz*sin_euz_dd-cos_eux_d*cos_euz_dd*sin_euz
    m_22 = -sin_euz*sin_euz_dd+cos_eux_d*cos_euz*cos_euz_dd
    m_23 = sin_eux_d*cos_euz_dd
    m_31 = sin_eux_d*sin_euz
    m_32 = -cos_euz*sin_eux_d
    m_33 = cos_eux_d

    rot = np.array([[m_11, m_12, m_13],
                    [m_21, m_22, m_23],
                    [m_31, m_32, m_33]], dtype='float64')

    rot = rot.transpose()

    if (rot[2][1] - rot[1][2]) < 0:
        q2m = -1
    else:
        q2m = 1

    if (rot[0][2] - rot[2][0]) < 0:
        q3m = -1
    else:
        q3m = 1

    if (rot[1][0] - rot[0][1]) < 0:
        q4m = -1
    else:
        q4m = 1

    q = np.zeros(4, dtype='float64')
    q[0] = sqrt(rot[0][0] + rot[1][1] + rot[2][2] + 1)/2
    q[1] = q2m * sqrt(rot[0][0] - rot[1][1] - rot[2][2] + 1)/2
    q[2] = q3m * sqrt(rot[1][1] - rot[0][0] - rot[2][2] + 1)/2
    q[3] = q4m * sqrt(rot[2][2] - rot[0][0] - rot[1][1] + 1)/2
    return q


def eulerzyx_to_quat(euler_z_y_x):
    """
    Converts a Euler ZYX angles into quaternion values

    Arguments:
        euler_z_y_x - array (eu_z, eu_y, eu_x)
    """
    eux = euler_z_y_x[2]*pi/180
    euy = euler_z_y_x[1]*pi/180
    euz = euler_z_y_x[0]*pi/180
    cos_eux = cos(eux)
    sin_eux = sin(eux)
    cos_euy = cos(euy)
    sin_euy = sin(euy)
    cos_euz = cos(euz)
    sin_euz = sin(euz)

    m_11 = cos_euy*cos_euz
    m_12 = sin_eux*sin_euy*cos_euz - sin_euz*cos_eux
    m_13 = sin_eux*sin_euz + sin_euy*cos_eux*cos_euz
    m_21 = sin_euz*cos_euy
    m_22 = sin_eux*sin_euy*sin_euz + cos_eux*cos_euz
    m_23 = -sin_eux*cos_euz + sin_euy*sin_euz*cos_eux
    m_31 = -sin_euy
    m_32 = sin_eux*cos_euy
    m_33 = cos_eux*cos_euy

    rot = np.array([[m_11, m_12, m_13],
                    [m_21, m_22, m_23],
                    [m_31, m_32, m_33]], dtype='float64')

    if (rot[2][1] - rot[1][2]) < 0:
        q2m = -1
    else:
        q2m = 1
    if (rot[0][2] - rot[2][0]) < 0:
        q3m = -1
    else:
        q3m = 1
    if (rot[1][0] - rot[0][1]) < 0:
        q4m = -1
    else:
        q4m = 1

    q = np.zeros(4, dtype='float64')
    q[0] = sqrt(rot[0][0] + rot[1][1] + rot[2][2] + 1)/2
    q[1] = q2m * sqrt(rot[0][0] - rot[1][1] - rot[2][2] + 1)/2
    q[2] = q3m * sqrt(rot[1][1] - rot[0][0] - rot[2][2] + 1)/2
    q[3] = q4m * sqrt(rot[2][2] - rot[0][0] - rot[1][1] + 1)/2
    return q


def rot_to_quat(rot):
    """
    Converts a rotation matrix into a quaternion

    Arguments:
        rot - array (3x3)
    """
    if (rot[2][1] - rot[1][2]) < 0:
        q2m = -1
    else:
        q2m = 1

    if (rot[0][2] - rot[2][0]) < 0:
        q3m = -1
    else:
        q3m = 1

    if (rot[1][0] - rot[0][1]) < 0:
        q4m = -1
    else:
        q4m = 1

    q = np.zeros(4, dtype='float64')

    q[0] = sqrt(rot[0][0] + rot[1][1] + rot[2][2] + 1)/2
    q[1] = q2m * sqrt(rot[0][0] - rot[1][1] - rot[2][2] + 1)/2
    q[2] = q3m * sqrt(rot[1][1] - rot[0][0] - rot[2][2] + 1)/2
    q[3] = q4m * sqrt(rot[2][2] - rot[0][0] - rot[1][1] + 1)/2
    return q


def QuatMult(quat1, quat2):
    """
    Multiplies two quaternion

    Arguments:
        quat1 - array (q0, qx, qy, qz)
        quat2 - array (q0, qx, qy, qz)
    """
    result = np.zeros(4, dtype='float64')
    result[0] = quat1[0]*quat2[0]-np.dot(quat1[1:], quat2[1:])
    result[1:] = quat1[0]*quat2[1:]+quat2[0]*quat1[1:] +\
        np.cross(quat1[1:], quat2[1:])
    return result


def RotatePos(pos, quat):
    """
    Rotates the provided pos by the specified quaternion

    Arguments:
        pos - array (x, y, z)
        quat - array (q0, qx, qy, qz)
    """
    assert len(pos) == 3, "Specify coordinates in 3D space"
    assert len(quat) == 4, "Provide unit quaternion values"
    assert abs(np.linalg.norm(quat)-1) < 0.00001, \
        "Quaternion does not have unit norm"
    temp = np.zeros(4, dtype='float64')
    result = np.zeros(4, dtype='float64')
    temp[0] = np.dot(quat[1:], pos)
    temp[1:] = np.cross(quat[1:], pos)+quat[0]*pos
    # result[0] = np.dot(temp[1:],-quat[1:])+temp[0]*quat[0]
    result[1:] = (
                  np.cross(temp[1:], -quat[1:]) +
                  temp[0]*quat[1:]+quat[0]*temp[1:])
    return result[1:]


def RotateQuat(quat1, quat2):
    """
    Rotates quat 2 by quat1

    Arguments:
        quat1 - array (q0, qx, qy, qz)
        quat2 - array (q0, qx, qy, qz)
    """
    quat1_conj = np.zeros(4, dtype='float64')
    quat1_conj[0] = quat1[0]
    quat1_conj[1:] = -quat1[1:]
    return QuatMult(QuatMult(quat1, quat2), quat1_conj)
