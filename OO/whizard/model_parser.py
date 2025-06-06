import re
from math import sqrt, pi


parameter_pattern = re.compile(r"parameter\s+(\w+)\s+=\s+(\d+\.?\d*(?:[eE]-?\d+)?)")
derived_pattern = re.compile(r"derived\s+(\w+)\s+=\s+([\w\s\(\)\-\+\*\/\.]*[\w\)])")


class ModelParser():

    _parameters: dict[str, float]
    _raw_derived_parameters: dict[str, str]

    def __init__(self, path: str):
        self._parameters = {}
        self._raw_derived_parameters = {}
        with open(path) as file:
            text = file.read()
        for match in re.finditer(parameter_pattern, text):
            name, value = match.group(1, 2)
            self._parameters[name] = float(value)
        for match in re.finditer(derived_pattern, text):
            name, value = match.group(1, 2)
            self._raw_derived_parameters[name] = value
        self._update_derived()


    def _update_derived(self):
        for name, value in self._raw_derived_parameters.items():
            derived_value = eval(value, (self._parameters | {"sqrt": sqrt, "pi": pi}))
            self._parameters[name] = derived_value


    def get_parameters_dict(self) -> dict[str, float]:
        return self._parameters.copy()


    def get_parameters_list(self) -> list[float]:
        return list(self._parameters.values())


    def set_parameter(self, name: str, value: float, update: bool = True):
        if name in self._parameters:
            if name in self._raw_derived_parameters:
                print(f"ERROR: {name} is a derived parameter, calculated as {self._raw_derived_parameters[name]}")
                raise KeyError
            self._parameters[name] = value
        else:
            raise KeyError
        if update:
            self._update_derived()


    def set_parameters(self, parameter_dict: dict[str, float]):
        for name, value in parameter_dict.items():
            self.set_parameter(name, value, update=False)
        # we call update manually after all other parameters are set
        self._update_derived()


    def add_derived_parameter(self, name: str, expression: str):
        self._raw_derived_parameters[name] = expression
        self._update_derived()


if __name__  == "__main__":
    model_parser = ModelParser("SM_ac.mdl")
    print(model_parser.get_parameters_dict())
    model_parser.set_parameter("mW", 55.1337)
    print(model_parser.get_parameters_dict())
