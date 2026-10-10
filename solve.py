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
    diag=np.zeros(n_points)
    upper=np.zeros(n_points-1)
    lower=np.zeros(n_points-1)
    boundary_term=np.zeros(n_points)
    diag[1:-1]=-(coef_i[1:]-coef_next[:-1])/dx
    upper[1:]=-coef_next[1:]/dx
    lower[:-1]=coef_i[:-1]/dx
    vel_ends = vel(x[[0, -1]], t)
    disp_ends = dis(x[[0, -1]], t)
    left_gradient = left_bc[1](t) if left_bc[0] == "neumann" else 0
    right_gradient = right_bc[1](t) if right_bc[0] == "neumann" else 0
    diag[0] = (-2 * coef_i[0] + 2 * vel_ends[0]) / dx
    upper[0] = -2 * coef_next[0] / dx
    boundary_term[0] = -2 * disp_ends[0] * left_gradient / dx
    diag[-1] = (-2 * vel_ends[1] + 2 * coef_next[-1]) / dx
    lower[-1] = 2 * coef_i[-1] / dx
    boundary_term[-1] = 2 * disp_ends[1] * right_gradient / dx
    matrix = diags([lower, diag, upper], [-1, 0, 1], format="csr")
    return matrix, boundary_term


  def rhs(conc, t, transport_matrix, boundary_term):
    return (transport_matrix @ conc + boundary_term - reaction_rate * conc**2 + source(x, t))
  dirichlet_nodes = []
  if left_bc[0] == "dirichlet":
      dirichlet_nodes.append((0, left_bc[1]))
  if right_bc[0] == "dirichlet":
      dirichlet_nodes.append((n_points - 1, right_bc[1]))

  conc = initial_conc(x).astype(float)
  for node, value_at in dirichlet_nodes:
      conc[node] = value_at(0.0)

  matrix, boundary_term = build_transport(0.0)
  rhs_old = rhs(conc, 0.0, matrix, boundary_term)

  times = [0.0]
  conc_history = [conc.copy()]
  
