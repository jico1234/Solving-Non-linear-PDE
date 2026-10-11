import numpy as np
import pandas as pd
from scipy.sparse import diags, identity
from scipy.sparse.linalg import spsolve
from pysr import PySRRegressor
import sympy as sp

def solve_nade(vel, dis, sou, rr, ic, l_bc, r_bc, x_start, x_end, t_end, n_points, n_steps):
    x = np.linspace(x_start, x_end, n_points)
    dx = x[1] - x[0]
    dt = t_end / n_steps
    x_mid = (x[:-1] + x[1:]) / 2
    def build_transport(t):
        v_mid = vel(x_mid, t)
        d_mid = dis(x_mid, t)
        coef_i = v_mid / 2 + d_mid / dx
        coef_next = v_mid / 2 - d_mid / dx
        diag = np.zeros(n_points)
        upper = np.zeros(n_points - 1)
        lower = np.zeros(n_points - 1)
        boundary_term = np.zeros(n_points)
        diag[1:-1] = -(coef_i[1:] - coef_next[:-1]) / dx
        upper[1:] = -coef_next[1:] / dx
        lower[:-1] = coef_i[:-1] / dx
      
        vel_ends = vel(x[[0, -1]], t)
        disp_ends = dis(x[[0, -1]], t)

        left_gradient = l_bc[1](t) if l_bc[0] == "neumann" else 0
        right_gradient = r_bc[1](t) if r_bc[0] == "neumann" else 0

        diag[0] = (-2 * coef_i[0] + 2 * vel_ends[0]) / dx
        upper[0] = -2 * coef_next[0] / dx
        boundary_term[0] = (
            -2 * disp_ends[0] * left_gradient / dx
        )

        diag[-1] = (-2 * vel_ends[1] + 2 * coef_next[-1]) / dx
        lower[-1] = 2 * coef_i[-1] / dx
        boundary_term[-1] = (
            2 * disp_ends[1] * right_gradient / dx
        )

        matrix = diags(
            [lower, diag, upper],
            [-1, 0, 1],
            format="csr"
        )

        return matrix, boundary_term

    def rhs(conc, t, transport_matrix, boundary_term):
        return (
            transport_matrix @ conc
            + boundary_term
            - rr * conc**2
            + sou(x, t)
        )

    dirichlet_nodes = []

    if l_bc[0] == "dirichlet":
        dirichlet_nodes.append((0, l_bc[1]))

    if r_bc[0] == "dirichlet":
        dirichlet_nodes.append((n_points - 1, r_bc[1]))

    conc = ic(x).astype(float)

    for node, value_at in dirichlet_nodes:
        conc[node] = value_at(0.0)

    matrix, boundary_term = build_transport(0.0)
    rhs_old = rhs(conc, 0.0, matrix, boundary_term)

    times = [0.0]
    conc_history = [conc.copy()]

    for step in range(n_steps):
        t = (step + 1) * dt

        matrix, boundary_term = build_transport(t)
        conc_new = conc.copy()

        for iteration in range(20):
            rhs_new = rhs(conc_new, t, matrix, boundary_term)

            residual = (
                conc_new - conc
                - 0.5 * dt * (rhs_new + rhs_old)
            )

            jacobian = (
                identity(n_points)
                - 0.5 * dt
                * (matrix - 2 * rr * diags(conc_new))
            ).tolil()

            for node, value_at in dirichlet_nodes:
                residual[node] = conc_new[node] - value_at(t)
                jacobian[node, :] = 0
                jacobian[node, node] = 1

            change = spsolve(jacobian.tocsc(), -residual)
            conc_new = conc_new + change

            if np.max(np.abs(change)) < 1e-10:
                break

        conc = conc_new
        rhs_old = rhs(conc, t, matrix, boundary_term)

        times.append(t)
        conc_history.append(conc.copy())

    return x, np.array(times), np.array(conc_history)

def make_dataframe(x, times, conc_history, vel, dis, sou):
  tables = []
  for t, conc in zip(times, conc_history):
      tables.append(pd.DataFrame({
          "t": t, "x": x, "C": conc,
          "v": vel(x, t), "D": dis(x, t), "Sou": sou(x, t),
      }))
  return pd.concat(tables, ignore_index=True)
#This is input for the solution that can be changed anytime, and different data sets may be generated.
if __name__ == "__main__":
    length = 10.0
    total_time = 10.0
    vel = lambda x, t: 0.5 * (1 + 0.2 * np.sin(2 * np.pi * t / total_time)) + 0 * x
    dis = lambda x, t: 0.05 * (1 + 0.5 * x / length)
    sou = lambda x, t: 0.02 * np.exp(-(x - 3) ** 2 / 0.5) * (t < 5)
    rr = 0.05
    ic=lambda x: 0*x
    l_bc = ("dirichlet", lambda t: 1.0)
    r_bc = ("neumenn", lambda t: 0)

    x, times, conc_history=solve_nade(vel, dis, sou, rr, ic, l_bc, r_bc, x_start=0, x_end=length, t_end=total_time, n_points=201, n_steps=400)

    df = make_dataframe(x, times, conc_history, vel, dis, sou)
    df.to_csv("ade_data.csv", index=False)

#Generated the data set, now we will find the pattern using pysr, verify the residual and r2 from the pysr_equations.csv to find out the best fit solution

data=pd.read_csv("ade_data.csv")
data = data[data["t"] > 0]
data = data.sample(frac=1, random_state=0)
train = data.iloc[:3000]
test = data.iloc[3000:5000]
inputs = ["x", "t", "v", "D", "Sou"]
max_res=0.001
min_r2=0.98
r2_g=0.98
x, t, v, D, So = sp.symbols("x t v D Sou")
sy_map = {"x": x, "t": t, "v": v, "D": D, "Sou": So}
length, total_time, rr = 10.0, 10.0, 0.05
v_expr = 0.5 * (1 + 0.2 * sp.sin(2 * sp.pi * t / total_time))
D_expr = 0.05 * (1 + 0.5 * x / length)
S_expr = sp.Piecewise((0.02 * sp.exp(-(x - 3) ** 2 / 0.5), t < 5), (0, True))
def r2_score(actual, predicted):
    return 1 - np.sum((actual - predicted) ** 2) / np.sum((actual - actual.mean()) ** 2)
def pde_residual_rms(C_expr, points):
    C = C_expr.subs({v: v_expr, D: D_expr, So: S_expr})
    residual = (sp.diff(C, t) + sp.diff(v_expr * C, x) - sp.diff(D_expr * sp.diff(C, x), x) + rr * C ** 2 - S_expr)
    f = sp.lambdify((x, t), residual, modules=["scipy", "numpy"])
    values = np.broadcast_to(f(points["x"].values, points["t"].values), (len(points),))
    return np.sqrt(np.mean(values ** 2))

def check_equations(csv_file, data, inputs, n_points=5000):
    equations = pd.read_csv(csv_file)
    actual = data["C"].values
    points = data.sample(n=min(n_points, len(data)), random_state=1)
    results = []
    for _, row in equations.iterrows():
        r2, residual = np.nan, np.nan
        try:
            C_expr = sp.sympify(row["equation"], locals=symbol_map)
            f = sp.lambdify((x, t, v, D, So), C_expr, modules=["scipy", "numpy"])
            predicted = np.broadcast_to(f(*[data[c].values for c in inputs]), actual.shape)
            r2 = r2_score(actual, predicted)
            if r2 >= r2_gate:
                residual = pde_residual_rms(C_expr, points)
        except Exception as err:
            print(f"  skipped complexity {row['complexity']}: {err}")

        accepted = bool(np.isfinite(r2) and np.isfinite(residual)
                        and r2 >= r2_gate and r2 > min_r2 and residual < max_res)
        results.append({"complexity": row["complexity"], "equation": row["equation"],
                        "r2": r2, "residual": residual, "accepted": accepted})
    report = pd.DataFrame(results)
    checked = report[report["r2"] >= r2_gate]
    if checked.empty:
        print(f"  No equation reached R2 >= {r2_gate}, so no residual was checked.")
    else:
        print(checked[["complexity", "r2", "residual", "accepted"]].to_string(index=False))

    good = report[report["accepted"]]
    return None if good.empty else good.sort_values("complexity").iloc[0]
unary_sets = [
    ["exp", "erfc", "sqrt"],
    ["exp", "erfc", "tanh"],
    ["exp", "erfc", "sin", "cos"],
    ["exp", "erfc", "log", "tanh"],
    ["exp", "erfc", "square", "tanh"],
]
best = None
run=0
for unary in unary_sets:
    run = run + 1
    print("Run", run, "- unary operators:", unary)
    model = PySRRegressor(niterations=100, maxsize=25, binary_operators=["+", "-", "*", "/"], unary_operators=unary, model_selection="best",)
    model.fit(train[inputs], train["C"])

csv_file = "pysr_equations_run" + str(run) + ".csv"
model.equations_.to_csv(csv_file, index=False)
best = check_equations(csv_file, data, inputs)
if best is not None:
    break

print("No equation")
if best is None:
    print("No solution")
else:
    print("solution")
    print("solution is :", best["equation"])
    
