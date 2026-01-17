"""
Test module for the busemann.py module.

See busemann.py for references.

.. Author: RJG

.. Version: 2024-07-04
"""

from busemann import BusemannDiffuser
from math import asin
import pytest

# Use values from example in Section 5 of Moelder (2003)
M2 = 3.0
k = 1.4
theta_23 = asin(k/M2)
bd = BusemannDiffuser(M2, theta_23)
    
def test_init():
    props = bd.properties()
    assert props.M1 == None
    assert props.M2 == pytest.approx(3.0)
    assert props.M3 == pytest.approx(2.48155)
    assert props.Pi == pytest.approx(0.958194)

def test_field_integration():
    r = 1.0
    dtheta = 0.001
    bd.generate_contour(r, dtheta)
    props = bd.properties()
    assert props.M1 == pytest.approx(5.76788)
    assert props.M2 == pytest.approx(3.0)
    assert props.M3 == pytest.approx(2.48155)
    assert props.Pi == pytest.approx(0.958194)

