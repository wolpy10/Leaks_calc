import sys
import importlib.util
import sys

module_path = "D:\\Files\\Work\\Hydrodynamic\\leaks\\src"
module_name = "data_process"

spec = importlib.util.spec_from_file_location(
    module_name, f"{module_path}\\{module_name}.py"
)
module = importlib.util.module_from_spec(spec)
sys.modules[module_name] = module
spec.loader.exec_module(module)

from data_process import Pump

Pump.load_params()

Pump.calc_leaks()
