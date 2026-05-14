import pandas as pd
import numpy as np
from combat.pycombat import pycombat

np.random.seed(42)
data = pd.DataFrame(np.random.randn(5, 10)) # 5 genes, 10 samples
batch = [1,1,1,1,1, 2,2,2,2,2]
mod = [[0,1,0,1,0, 1,0,1,0,1], [25, 30, 40, 50, 60, 25, 30, 40, 50, 60]] # two covariates
try:
    res = pycombat(data, batch, mod=mod)
    print("Success with mod!")
except Exception as e:
    print("Error:", e)
