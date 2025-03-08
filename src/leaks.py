import math
from data_process import Pump
from units import m3_hr

Pump.load_params()


def calc_v1u() -> float:
    [D1_ind, d_sleeve, n_ind, betta2_ind] = Pump.get_params(
        ["D1_ind", "d_sleeve", "n_ind", "betta2_ind"], "inducer"
    )
    [Q] = Pump.get_params(["Q"], "flow")

    D_calc = (D1_ind + d_sleeve) / 2
    u_ind = D_calc / 2 * n_ind
    vm2_ind = 4 * Q / (math.pi * (pow(D1_ind, 2) - pow(d_sleeve, 2)))
    vu2_ind = u_ind - vm2_ind / math.tan(betta2_ind)

    v1u = vu2_ind * D_calc / D1_ind
    return v1u


def calc_leak(seal: str, holes_exist=False, inducer_exist=False) -> float:
    g = 9.807

    [Q, nu] = Pump.get_params(["Q", "nu"], "flow")
    [
        n,
        D2,
        b2,
        D1,
        d_sleeve,
        z,
        betta2,
        sigma2,
    ] = Pump.get_params(
        [
            "n",
            "D2",
            "b2",
            "D1",
            "d_sleeve",
            "z",
            "betta2",
            "sigma2",
        ],
        "impeller",
    )
    
    [D_seal, L_seal, delta] = Pump.get_params(["D_seal", "L_seal", "delta"], seal)

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
        v1u = calc_v1u()
    else:
        v1u = 0

    R1 = (D1 + d_sleeve) / 2
    v2u = (H_theory * g / n + v1u * R1) / R2

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
        "Q_leak_shroud": round(
            calc_leak(inducer_exist=True, seal="shroud") / m3_hr, 3
        ),
        "Q_leak_hub": None,
        "Q_leak_shaft": None,
    }
)
