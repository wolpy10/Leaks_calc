from pathlib import Path
import json
import math
from units import *


class Pump:
    __data_path = Path(__file__).resolve().parent / "data"
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
                "b2": float(Pump.__data["impeller"]["b2"]) * mm,
                "z": int(Pump.__data["impeller"]["z"]),
                "betta2": float(Pump.__data["impeller"]["betta2"]) * deg,
                "sigma2": float(Pump.__data["impeller"]["sigma2"]) * mm,
                "D_seal": float(Pump.__data["impeller"]["D_seal"]) * mm,
                "L_seal": float(Pump.__data["impeller"]["L_seal"]) * mm,
                "delta": float(Pump.__data["impeller"]["delta"]) * mm,
                "nu": float(Pump.__data["impeller"]["nu"]) * sSt,
            }
            data_inducer = {
                "D1": float(Pump.__data["inducer"]["D1"]) * mm,
                "d_sleeve": float(Pump.__data["inducer"]["d_sleeve"]) * mm,
                "n_ind": float(Pump.__data["inducer"]["n_ind"]) * rpm,
                "betta2_ind": float(Pump.__data["inducer"]["betta2_ind"]) * deg,
            }
            data_flow = {"Q": Pump.__data["flow"]["Q"] * m3_hr}
            Pump.__data = {
                "impeller": data_impeller,
                "inducer": data_inducer,
                "flow": data_flow,
            }
        except KeyError:
            print("JSON doesn't match the template")
            raise KeyError

    def get_params(variables: list[str], part: str) -> list:
        try:
            extracted_variables = []
            for var in variables:
                extracted_variables.append(Pump.__data[part.lower()][var])
            return extracted_variables
        except KeyError:
            print("The data hasn't been defined yet")

    def get_leaks():
        return Pump.__leaks

    def write_leaks(results: dict):
        Pump.__leaks = results
        print(results)
        try:
            with open(
                Pump.__data_path / "output_data.json", "w", encoding="utf-8"
            ) as file:
                file.writelines(json.dumps(Pump.__leaks))
        except FileNotFoundError:
            print("File not found")
            raise FileNotFoundError
