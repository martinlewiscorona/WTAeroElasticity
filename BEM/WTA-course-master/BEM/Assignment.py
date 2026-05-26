import matlab.engine
import pandas as pd
import numpy as np
from dataclasses import dataclass

@dataclass
class Element:
    name: str()
    rho: float()
    Area: float()
    EI: float()
    xi: float()

# class StructuralModel:
#     def __init__():
        
#     def M():
        
#     def K():
        
#     def C():
        
#     def F():
        
        
eng = matlab.engine.start_matlab()


# example input
V = 12.0
omega = 1.267
pitch = 0.0

# call MATLAB
Rx, FN, FT, P = eng.BEM(V, omega, pitch, nargout=4)

# convert to numpy
Rx = np.array(Rx).flatten()
FN = np.array(FN).flatten()
FT = np.array(FT).flatten()

