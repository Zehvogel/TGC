import json


class AltSetupHandler():
    _sm_ref: dict[str, float]
    _variations: dict[str, list[float]]
    _alt_setup: dict[str, dict[str, float]]
    _mirror: bool


    def __init__(self, path: str, mirror: bool = True):
        with open(path) as f:
            conf = json.load(f)
        self._sm_ref = conf["SM"]
        self._variations = conf["variations"]
        self._mirror = mirror
        self._alt_setup = self._make_alt_configs()


    def make_name(self, parameter: str, var: float) -> str:
        # need to get rid of - to use as c++ identifier later
        alt_name = f"{parameter}_{'pos' if var > 0. else 'neg'}_{abs(var):.0e}".replace("-", "m")
        return alt_name


    def _make_variation(self, parameter: str, var: float) -> dict[str, float]:
        alt_conf = self._sm_ref.copy()
        alt_conf[parameter] += var
        return alt_conf


    def _get_parameters_variations_it(self):
        for parameter, variations in self._variations.items():
            for var in variations:
                yield parameter, var
                if self._mirror:
                    yield parameter, -var


    def _make_alt_configs(self) -> dict[str, dict[str, float]]:
        alt_configs = {}
        for parameter, var in self._get_parameters_variations_it():
            alt_name = self.make_name(parameter, var)
            alt_conf = self._make_variation(parameter, var)
            alt_configs[alt_name] = alt_conf
        return alt_configs


    def get_alt_setup(self):
        # would probably be cleaner to return a deepcopy :/
        return self._alt_setup


    def get_variations(self):
        return self._variations


    def get_named_variations_it(self):
        for parameter, var in self._get_parameters_variations_it():
            alt_name = self.make_name(parameter, var)
            yield alt_name, var
