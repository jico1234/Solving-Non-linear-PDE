import numpy as np
import pandas as pd
from scipy.sparse import diags, identity
from scipy.sparse.linalg import spsolve


def solve_nade(vel, dis, sou, rr, ic, l_bc, r_bc, x_start, x_end, t_end, n_points, n_step):
  x=np.linespace(x_start, x_end, n_points)
  dx=x[1]-x[0]
  dt=t_end/n_steps
  x_mid=(x[:-1]+ x[1:])/2
  def build_trans(t):
    v_mid=vel(x_mid, t)
    d_mid=dis(x_mid, t)
    coef_i=v_mid/2+d_mid/dx
    coef_next=v_mid/2-d_mid/dx
