import math
import json

def process_data(data: dict, part: str) -> dict:
    if part == "impeller":
        data_impeller = {
            "n": float(data["impeller"]["n"]) * math.pi / 30,
            "D2": float(data["impeller"]["D2"]) * 1e-3,
            "b2": float(data["impeller"]["b2"]) * 1e-3,
            "z": int(data["impeller"]["z"]),
            "betta2": float(data["impeller"]["betta2"]) * math.pi / 180,
            "sigma2": float(data["impeller"]["sigma2"]) * 1e-3,
            "D_seal": float(data["impeller"]["D_seal"]) * 1e-3,
            "L_seal": float(data["impeller"]["L_seal"]) * 1e-3,
            "delta": float(data["impeller"]["delta"]) * 1e-3,
            "nu": float(data["impeller"]["nu"]) * 1e-6,
        }
        data_inducer = {}
    else:
        data_impeller = {}
        data_inducer = {
            "D1": float(data["inducer"]["D1"]) * 1e-3,
            "d_sleeve": float(data["inducer"]["d_sleeve"]) * 1e-3,
            "n": float(data["inducer"]["n"]) * math.pi / 30,
            "betta2_ind": float(data["inducer"]["betta2_ind"]),
        }

    return {"impeller": data_impeller, "inducer": data_inducer}


def get_data(part: str) -> list:
    if part == "impeller":
        data = process_data(data, part)["impeller"]
        return [data["n"], data["D2"], data["b2"], data["z"], data["betta2"], data["sigma2"],
                data["D_seal"], data["L_seal"], data["delta"], data["nu"]]
    data = process_data(data, part)["inducer"]
    return [data["D1"], data["d_sleeve"], data["n"], data["betta2_ind"]]


def calc_v1u(Q: float):
    
    return NotImplemented


def calc_leak(Q: float, holes_exist=False, inducer_exist=False) -> float:
    g = 9.807

    Q /= 3600
    [n, D2, b2, z, betta2, sigma2, D_seal, L_seal,
        delta, nu] = get_data(part="impeller")
    R2 = D2 / 2

    y = 1 - 3.625 * math.sin(betta2) / (z + 3.625 * math.sin(betta2))
    psi_2 = 1 - z * sigma2 / (2 * math.pi * R2 * math.sin(betta2))

    u2 = n * R2

    H_theory = n / g * (pow(R2, 2) * y * n - R2 * Q /
                        (2 * math.pi * b2 * psi_2 * R2 * math.tan(betta2)))
    if inducer_exist:
        v1u = NotImplemented
        
    v2u = H_theory * g / (n * R2)
    H_seal = H_theory - pow(v2u, 2) / (2 * g) - \
        pow(u2, 2) / (8 * g) * (1 - pow(D_seal / D2, 2))

    mu0, mu1 = 1, 0.01
    while abs(mu1 - mu0) / mu0 * 100 > 0.01:
        mu0 = mu1
        vu = D_seal * n / 2
        Q_leak = mu0 * math.pi * D_seal * delta * math.sqrt(2 * g * H_seal)
        v0 = Q_leak / (math.pi * D_seal * delta)
        Re = 2 * delta * math.sqrt(pow(v0, 2) + pow(vu, 2)) / nu

        lambda_smooth = 0.316 / pow(Re, 0.25)
        lambda_seal = lambda_smooth * math.sqrt(1 + 1 / (4 * pow(1 + 1.3 * math.sqrt(lambda_smooth), 2))
                                                * pow(D_seal * n / v0,  2))

        mu1 = 1 / math.sqrt(1.3 + lambda_seal * L_seal / (2 * delta))
    mu = mu1
    Q_leak = mu * math.pi * D_seal * delta * math.sqrt(2 * g * H_seal)
    return Q_leak


print(f"{(calc_leak(Q=450) * 3600):.3f}")
