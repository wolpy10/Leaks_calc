import math
from data_process import Pump
from units import m3_hr

Pump.load_params()


def calc_v1u(Q: float):

    return NotImplemented


def calc_leak(holes_exist=False, inducer_exist=False) -> float:
    g = 9.807

    [Q] = Pump.get_params(["Q"], "flow")
    [n, D2, b2, z, betta2, sigma2, D_seal, L_seal, delta, nu] = Pump.get_params(
        ["n", "D2", "b2", "z", "betta2", "sigma2", "D_seal", "L_seal", "delta", "nu"],
        "impeller",
    )

    R2 = D2 / 2

    y = 1 - 3.625 * math.sin(betta2) / (z + 3.625 * math.sin(betta2))
    psi_2 = 1 - z * sigma2 / (2 * math.pi * R2 * math.sin(betta2))

    u2 = n * R2

    H_theory = (
        n
        / g
        * (
            pow(R2, 2) * y * n
            - R2 * Q / (2 * math.pi * b2 * psi_2 * R2 * math.tan(betta2))
        )
    )
    if inducer_exist:
        [D1, d_sleeve, n_ind, betta2_ind] = Pump.get_params(
            ["D1", "d_sleeve", "n_ind", "betta2_ind"], "inducer"
        )
        v1u = NotImplemented

    v2u = H_theory * g / (n * R2)
    H_seal = (
        H_theory
        - pow(v2u, 2) / (2 * g)
        - pow(u2, 2) / (8 * g) * (1 - pow(D_seal / D2, 2))
    )

    mu0, mu1 = 1, 0.01
    while abs(mu1 - mu0) / mu0 * 100 > 0.01:
        mu0 = mu1
        vu = D_seal * n / 2
        Q_leak = mu0 * math.pi * D_seal * delta * math.sqrt(2 * g * H_seal)
        v0 = Q_leak / (math.pi * D_seal * delta)
        Re = 2 * delta * math.sqrt(pow(v0, 2) + pow(vu, 2)) / nu

        lambda_smooth = 0.316 / pow(Re, 0.25)
        lambda_seal = lambda_smooth * math.sqrt(
            1
            + 1
            / (4 * pow(1 + 1.3 * math.sqrt(lambda_smooth), 2))
            * pow(D_seal * n / v0, 2)
        )

        mu1 = 1 / math.sqrt(1.3 + lambda_seal * L_seal / (2 * delta))
    mu = mu1
    Q_leak = mu * math.pi * D_seal * delta * math.sqrt(2 * g * H_seal)
    return Q_leak


Pump.write_leaks(
    {
        "Q_leak_shroud": round(calc_leak(Pump.get_params(["Q"], "flow")) / m3_hr, 3),
        "Q_leak_hub": None,
        "Q_leak_shaft": None,
    }
)
