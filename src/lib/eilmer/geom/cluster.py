# cluster.py
"""
Cluster functors for grid generation.

We use these so that we can pass around the parameters and cluster function
for generating an arrays of parametric ordinates [0.0 ... 1.0].

PJ, 2020-07-06
"""

from abc import ABC, abstractmethod
import numpy as np
from eilmer.roberts import distribute_points_1


class ClusterFunction(ABC):
    @abstractmethod
    def __repr__():
        pass

    @abstractmethod
    def distribute_parameter_values(self, nv):
        pass


class LinearFunction(ClusterFunction):
    def __init__(self):
        return

    def __repr__(self):
        return "LinearFunction()"

    def distribute_parameter_values(self, nv):
        """
        Returns an array of uniformly-distributed sample points in range [0.0 ... 1.0].
        """
        return np.linspace(0.0, 1.0, nv)


class RobertsFunction(ClusterFunction):
    """
    Roberts' boundary-layer-like clustering function.

    Note that, compared with the old definition that we delegate the actual work to,
    the ends are renamed 0, 1 to align with the Eilmer4 notation.
    """
    _slots_ = ['end0', 'end1', 'beta']

    def __init__(self, end0, end1, beta):
        """
        Store the cluster parameters for later use.
        """
        self.end0 = end0
        self.end1 = end1
        self.beta = beta
        return

    def __repr__(self):
        return f"RobertsFunction(end0={self.end0}, end1={self.end1}, beta={self.beta})"

    def distribute_parameter_values(self, nv):
        """
        Returns an array of clustered sample points in range [0.0 ... 1.0].
        """
        return distribute_points_1(0.0, 1.0, nv-1, self.end0, self.end1, self.beta)
