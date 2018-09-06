#!/usr/bin/env python

# Copyright (C) 2017 Udacity Inc.
#
# This file is part of Robotic Arm: Pick and Place project for Udacity
# Robotics nano-degree program
#
# All Rights Reserved.

# Author: Harsh Pandya

# import modules
import numpy as np
import rospy
import tf
from kuka_arm.srv import *
from trajectory_msgs.msg import JointTrajectory, JointTrajectoryPoint
from geometry_msgs.msg import Pose
from mpmath import *
from sympy import *


def get_cosine(a, b, c):
    # Use atan2 function to avoid sin/cos errors

    cos_angle_c = (c * c - a * a - b * b) / (-2 * a * b)
    sin_angle_c = sqrt(1 - cos_angle_c * cos_angle_c)
    angle_c_rad = atan2(sin_angle_c, cos_angle_c)

    return angle_c_rad  # return angle in radians


def t_matrix(alpha, a, d, q):
    trans_matrix = Matrix([
        [cos(q), -sin(q), 0, a],
        [sin(q) * cos(alpha), cos(q) * cos(alpha), -sin(alpha), -sin(alpha) * d],
        [sin(q) * sin(alpha), cos(q) * sin(alpha), cos(alpha), cos(alpha) * d],
        [0, 0, 0, 1]
    ])


    return trans_matrix


def handle_calculate_IK(req):
    rospy.loginfo("Received %s eef-poses from the plan" % len(req.poses))
    if len(req.poses) < 1:
        print
        "No valid poses received"
        return -1
    else:

        ### Your FK code here
        # Create symbols
        q1, q2, q3, q4, q5, q6, q7 = symbols('q1:8')  # theta
        a0, a1, a2, a3, a4, a5, a6 = symbols('a0:7')
        d1, d2, d3, d4, d5, d6, d7 = symbols('d1:8')
        alpha0, alpha1, alpha2, alpha3, alpha4, alpha5, alpha6 = symbols('alpha0:7')
        # Create Modified DH parameters

        s = {alpha0: 0, a0: 0, d1: 0.75, q1: q1,
             alpha1: -pi / 2, a1: 0.35, d2: 0, q2: -pi / 2 + q2,
             alpha2: 0, a2: 1.25, d3: 0, q3: q3,
             alpha3: -pi / 2, a3: -0.054, d4: 1.5, q4: q4,
             alpha4: pi / 2, a4: 0, d5: 0, q5: q5,
             alpha5: -pi / 2, a5: 0, d6: 0, q6: q6,
             alpha6: 0, a6: 0, d7: 0.303, q7: 0}

        # Define Modified DH Transformation matrix
        T0_1 = t_matrix(alpha0, a0, d1, q1)
        T0_1 = T0_1.subs(s)

        T1_2 = t_matrix(alpha1, a1, d2, q2)
        T1_2 = T1_2.subs(s)

        T2_3 = t_matrix(alpha2, a2, d3, q3)
        T2_3 = T2_3.subs(s)

        T3_4 = t_matrix(alpha3, a3, d4, q4)
        T3_4 = T3_4.subs(s)

        T4_5 = t_matrix(alpha4, a4, d5, q5)
        T4_5 = T4_5.subs(s)

        T5_6 = t_matrix(alpha5, a5, d6, q6)
        T5_6 = T5_6.subs(s)

        T6_G = t_matrix(alpha6, a6, d7, q7)
        T6_G = T6_G.subs(s)

        # Create individual transformation matrices
        T0_3 = simplify(T0_1 * T1_2 * T2_3)  # base to link_3
        T3_6 = simplify(T3_4 * T4_5 * T5_6)  # Link3 to link_6
        T0_6 = simplify(T0_3 * T3_6)  # base to link_6
        T0_G = simplify(T0_6 * T6_G)  # base to Gripper

        roll, pitch, yaw = symbols('roll pitch yaw')
        R_z = Matrix([[cos(yaw), -sin(yaw), 0],
                      [sin(yaw), cos(yaw), 0],
                      [0, 0, 1]])

        R_y = Matrix([[cos(pitch), 0, sin(pitch)],
                      [0, 1, 0],
                      [-sin(pitch), 0, cos(pitch)]])

        R_x = Matrix([[1, 0, 0],
                      [0, cos(roll), -sin(roll)],
                      [0, sin(roll), cos(roll)]])

        R_corr = R_z.subs(yaw, pi) * R_y.subs(pitch, -pi / 2)
        R_G = R_z * R_y * R_x
        # Extract rotation matrices from the transformation matrices
        R0_3 = T0_3[0:3, 0:3]

        EE_length = 0.303

        # Initialize service response
        joint_trajectory_list = []
        for x in xrange(0, len(req.poses)):
            # IK code starts here
            joint_trajectory_point = JointTrajectoryPoint()

            # Extract end-effector position and orientation from request
            # px,py,pz = end-effector position
            # roll, pitch, yaw = end-effector orientation
            px = req.poses[x].position.x
            py = req.poses[x].position.y
            pz = req.poses[x].position.z

            (roll, pitch, yaw) = tf.transformations.euler_from_quaternion(
                [req.poses[x].orientation.x, req.poses[x].orientation.y,
                 req.poses[x].orientation.z, req.poses[x].orientation.w])

            # IK code here
            # Compensate for rotation discrepancy between DH parameters and Gazebo
            RB_G = R_G.subs({'roll': roll, 'pitch': pitch, 'yaw': yaw}) * R_corr

            # Calculate joint angles using Geometric IK method

            wx = px - EE_length * RB_G[0, 2]
            wy = py - EE_length * RB_G[1, 2]
            wz = pz - EE_length * RB_G[2, 2]

            # Calculate Theta 1 to 3
            theta1 = atan2(wy, wx)
            xy_hypo = sqrt(wx * wx + wy * wy)

            length_j2_3 = sqrt(s[a2] * s[a2] + s[d3] * s[d3])
            length_j3_WC = sqrt(s[a3] * s[a3] + s[d4] * s[d4])
            vectorB = sqrt((wz - s[d1]) * (wz - s[d1]) + (xy_hypo - s[a1]) * (xy_hypo - s[a1]))
            angleVectorB = atan2(wz - s[d1], xy_hypo - s[a1])

            angle_a = get_cosine(length_j2_3, vectorB, length_j3_WC)
            angle_b = get_cosine(length_j3_WC, length_j2_3, vectorB)

            theta2 = pi / 2. - ((angleVectorB + angle_a))
            theta3 = pi / 2. - angle_b - 0.036

            # Calculate Theta 4 to 6

            R0_3_inv = R0_3.transpose()
            RHS = R0_3_inv * RB_G
            RHS_num = RHS.evalf(subs={q1: theta1, q2: theta2, q3: theta3})

            theta4 = atan2(RHS_num[2, 2], -RHS_num[0, 2])
            theta5 = atan2(sqrt(RHS_num[0, 2] * RHS_num[0, 2] + RHS_num[2, 2] * RHS_num[2, 2]), RHS_num[1, 2])
            theta6 = atan2(-RHS_num[1, 1], RHS_num[1, 0])

            # In the next line replace theta1,theta2...,theta6 by your joint angle variables
            joint_trajectory_point.positions = [theta1, theta2, theta3, theta4, theta5, theta6]
            joint_trajectory_list.append(joint_trajectory_point)
        rospy.loginfo("length of Joint Trajectory List: %s" % len(joint_trajectory_list))
        return CalculateIKResponse(joint_trajectory_list)


def IK_server():
    # initialize node and declare calculate_ik service
    rospy.init_node('IK_server')
    s = rospy.Service('calculate_ik', CalculateIK, handle_calculate_IK)
    print
    "Ready to receive an IK request"
    rospy.spin()


if __name__ == "__main__":
    IK_server()
