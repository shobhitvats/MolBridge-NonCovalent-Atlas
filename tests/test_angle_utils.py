import numpy as np
from utils.angle_utils import angles_between, filter_by_angle

def test_angles_between_basic():
    v1 = np.array([[1,0,0],[0,1,0]])
    v2 = np.array([[1,0,0],[1,0,0]])
    ang = angles_between(v1, v2)
    assert ang.shape == (2,)
    assert ang[0] == 0.0
    assert 89.0 < ang[1] < 91.0

def test_filter_by_angle():
    v1 = np.array([[1,0,0],[0,1,0]])
    v2 = np.array([[1,0,0],[1,0,0]])
    mask = filter_by_angle(v1, v2, 10)
    assert mask.tolist() == [True, False]