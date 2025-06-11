import math

def must_be_greater_than_1_by_pi(k):
    if k <= 1 / math.pi:
        raise ValueError("k must be greater than 1/pi.")