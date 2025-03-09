import math
from data_process import Pump

Pump.load_params()

Pump.calc_leaks(inducer_exist=True, account_shaft_leak=True, holes_exist=True)
