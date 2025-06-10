import json


def make_variation(sm_cf: dict[str, float], parameter: str, var: float, negate: bool = False) -> tuple[str, dict[str, float]]:
    alt_name = f"{parameter}_{'pos' if not negate else 'neg'}_{var:.0e}".replace("-", "m") # need to get rid of - to use as c++ identifier later
    alt_conf = sm_cf.copy()
    alt_conf[parameter] += var if not negate else -1.*var
    return alt_name, alt_conf


def make_alt_configs(sm_cf: dict[str, float], var_cf: dict[str, list[float]], mirror: bool = True):
    alt_configs = {}
    for parameter, variations in var_cf.items():
        for var in variations:
            alt_name, alt_conf = make_variation(sm_cf, parameter, var)
            alt_configs[alt_name] = alt_conf
            if mirror:
                alt_name, alt_conf = make_variation(sm_cf, parameter, var, negate=True)
                alt_configs[alt_name] = alt_conf
    return alt_configs


def load_config(path: str, mirror: bool = True):
    with open(path) as f:
        conf = json.load(f)
    return make_alt_configs(conf["SM"], conf["variations"])
