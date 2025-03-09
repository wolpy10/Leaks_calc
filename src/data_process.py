from pathlib import Path
import json
import math
from units import *


class Pump:
    __data_path = Path(__file__).resolve().parent.parent / "data"
    __data = {}
    __leaks = {}

    def load_params():
        try:
            with open(
                Pump.__data_path / "input_data.json", "r", encoding="utf-8"
            ) as file:
                Pump.__data = json.load(file)
            Pump.__process_data()
        except FileNotFoundError:
            print("File not found")
            raise FileNotFoundError

    def __process_data():
        try:
            data_impeller = {
                "n": float(Pump.__data["impeller"]["n"]) * rpm,
                "D2": float(Pump.__data["impeller"]["D2"]) * mm,
                "D1": float(Pump.__data["impeller"]["D1"]) * mm,
                "d_sleeve": float(Pump.__data["impeller"]["d_sleeve"]) * mm,
                "b2": float(Pump.__data["impeller"]["b2"]) * mm,
                "z": int(Pump.__data["impeller"]["z"]),
                "betta2": float(Pump.__data["impeller"]["betta2"]) * deg,
                "sigma2": float(Pump.__data["impeller"]["sigma2"]) * mm,
            }
            data_inducer = {
                "D1_ind": float(Pump.__data["inducer"]["D1_ind"]) * mm,
                "d_sleeve": float(Pump.__data["inducer"]["d_sleeve"]) * mm,
                "n_ind": float(Pump.__data["inducer"]["n_ind"]) * rpm,
                "betta2_ind": float(Pump.__data["inducer"]["betta2_ind"]) * deg,
            }
            data_flow = {
                "Q": Pump.__data["flow"]["Q"] * m3_hr,
                "nu": float(Pump.__data["flow"]["nu"]) * sSt,
                "H_stage": float(Pump.__data["flow"]["H_stage"]) * m,
            }
            data_shroud_seal = {
                "D_seal": float(Pump.__data["shroud_seal"]["D_seal"]) * mm,
                "L_seal": float(Pump.__data["shroud_seal"]["L_seal"]) * mm,
                "delta": float(Pump.__data["shroud_seal"]["delta"]) * mm,
            }
            data_hub_seal = {
                "D_seal": float(Pump.__data["hub_seal"]["D_seal"]) * mm,
                "L_seal": float(Pump.__data["hub_seal"]["L_seal"]) * mm,
                "delta": float(Pump.__data["hub_seal"]["delta"]) * mm,
            }
            Pump.__data = {
                "impeller": data_impeller,
                "inducer": data_inducer,
                "flow": data_flow,
                "shroud": data_shroud_seal,
                "hub": data_hub_seal,
                "shaft": data_hub_seal,
            }
        except KeyError:
            print("JSON doesn't match the template")
            raise KeyError

    def __get_params(variables: list[str], part: str) -> list:
        try:
            extracted_variables = []
            for var in variables:
                extracted_variables.append(Pump.__data[part.lower()][var])
            return extracted_variables
        except TypeError:
            print("The data hasn't been defined yet")
            raise TypeError

    def __write_leaks(results: dict):
        Pump.__leaks = results
        try:
            with open(
                Pump.__data_path / "output_data.json", "w", encoding="utf-8"
            ) as file:
                file.writelines(json.dumps(Pump.__leaks))
        except FileNotFoundError:
            print("File not found")
            raise FileNotFoundError

    def __calc_v1u() -> float:
        [D1_ind, d_sleeve, n_ind, betta2_ind] = Pump.__get_params(
            ["D1_ind", "d_sleeve", "n_ind", "betta2_ind"], "inducer"
        )
        [Q] = Pump.__get_params(["Q"], "flow")

        D_calc = (D1_ind + d_sleeve) / 2
        u_ind = D_calc / 2 * n_ind
        vm2_ind = 4 * Q / (math.pi * (pow(D1_ind, 2) - pow(d_sleeve, 2)))
        vu2_ind = u_ind - vm2_ind / math.tan(betta2_ind)

        v1u = vu2_ind * D_calc / D1_ind
        return v1u

    def __calc_leak(seal: str, inducer_exist=False, holes_exist=False) -> float:
        g = 9.807

        [Q, nu] = Pump.__get_params(["Q", "nu"], "flow")
        [
            n,
            D2,
            b2,
            D1,
            d_sleeve,
            z,
            betta2,
            sigma2,
        ] = Pump.__get_params(
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

        [D_seal, L_seal, delta] = Pump.__get_params(["D_seal", "L_seal", "delta"], seal)

        R2 = D2 / 2

        y = 1 - 3.625 * math.sin(betta2) / (z + 3.625 * math.sin(betta2))
        psi_2 = 1 - z * sigma2 / (2 * math.pi * R2 * math.sin(betta2))

        u2 = n * R2

        if holes_exist and seal == "shaft":
                [H_seal] = Pump.__get_params(["H_stage"], "flow")
        else:
            H_theory = (
                n
                / g
                * (
                    pow(R2, 2) * y * n
                    - R2 * Q / (2 * math.pi * b2 * psi_2 * R2 * math.tan(betta2))
                )
            )

            if inducer_exist:
                v1u = Pump.__calc_v1u()
            else:
                v1u = 0

            R1 = (D1 + d_sleeve) / 2
            v2u = (H_theory * g / n + v1u * R1) / R2

            H_seal = (
                H_theory
                - pow(v2u, 2) / (2 * g)
                - pow(u2, 2) / (8 * g) * (1 - pow(D_seal / D2, 2))
            )

            if seal == "shaft":
                [H_stage] = Pump.__get_params(["H_stage"], "flow")
                H_seal = H_stage - H_seal
                print(H_seal)
                if H_seal < 0:
                    print("The direction of shaft flow is incorrect")
                    raise ValueError

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

    def calc_leaks(inducer_exist: bool, account_shaft_leak: bool, holes_exist: bool):
        results = {}
        leak_shroud = Pump.__calc_leak(
            seal="shroud", inducer_exist=inducer_exist, holes_exist=holes_exist
        )
        results.update({"Q_leak_shroud": round(leak_shroud / m3_hr, 3)})
        if holes_exist:
            leak_hub = Pump.__calc_leak(
                seal="hub", inducer_exist=inducer_exist, holes_exist=holes_exist
            )
            results.update({"Q_leak_hub": round(leak_hub / m3_hr, 3)})

        if account_shaft_leak:
            leak_shaft = Pump.__calc_leak(
                seal="shaft", inducer_exist=inducer_exist, holes_exist=holes_exist
            )
            results.update({"Q_leak_shaft": round(leak_shaft / m3_hr, 3)})

        Pump.__write_leaks(results)

    def get_leaks():
        return Pump.__leaks
