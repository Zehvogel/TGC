from .Dataset import Dataset
import json
import ROOT
from typing import Any
from itertools import chain, filterfalse
from multiprocessing import Pool

def id_genrange():
    # the generator overlords have chosen these numbers...
    higgs_inc_range = range(402001, 402014 + 1)
    # some are missing:
    missing = [500000 + i for i in [
        5, 7, 9, 11, 23, 24, 31, 32, 41, 42, 49, 50, 61, 63, 65, 67, 69,
        71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 109,
        111, 121, 123, 129, 133, 137
        ]]
    all2fto5f_range = filterfalse(lambda i: i in missing, range(500001, 500248 + 1))
    the6f_range     = range(402301, 402324 + 1)
    for id in chain(higgs_inc_range, all2fto5f_range, the6f_range):
        yield str(id)


class ILDMC2020Dataset(Dataset):
    """Implementation to create a Dataset from the ILD 2020 mass production

    Relies on the filename scheme of the production.
    """

    _meta_dict: dict[str, Any] = {}
    _name2id: dict[str, str] = {}

    def __init__(self, files: str | list[str], genmeta: str):
        super().__init__()
        with open(genmeta) as meta_file:
            self._meta_dict = json.load(meta_file)
            for id in id_genrange():
                meta = self._meta_dict[id]
                process_name = meta["process_names"]
                e_pol = meta["polarization1"]
                p_pol = meta["polarization2"]
                # print(f"adding {process_name}_e{e_pol}p{p_pol} to name2id")
                self._name2id[f"{process_name}_e{e_pol}p{p_pol}"] = id

        # helper
        def process_lines(content):
            for line in content:
                path = line.strip()
                self._add_file(path)

        if isinstance(files, str):
            with open(files) as file:
                process_lines(file)
        elif isinstance(files, list):
            process_lines(files)
        else:
            return TypeError


    @staticmethod
    def _parse_path(path: str) -> str:
        dir, _, fname = path.rpartition("/")
        # the mc-2020 filenames follow a certain fixed naming scheme, example:
        # rv02-02.sv02-02.mILD_l5_o1_v02.E250-SetA.I500102.P4f_sze_sl.eL.pR.n024.d_dstm_15180_122_mini-DST.edm4hep.root
        # <reco v>.<sim v>.<det model>.<machine setting>.<process id>.<process name>.<e pol>.<p pol>.<... other stuff
        parts = fname.split(".")
        process_name = parts[5].lstrip("P")
        e_pol = parts[6]
        p_pol = parts[7]
        return f"{process_name}_{e_pol}{p_pol}"


    def _get_meta_for_process(self, process_name: str) -> dict[str, Any]:
        id = self._name2id[process_name]
        meta = self._meta_dict[id]
        # now decide what parts we actually want
        process_meta = {}
        process_meta["xsec_fb"] = float(meta["cross_section_in_fb"])
        process_meta["n_events"] = 0
        process_e_pol = -1.0 if "eL" in process_name else 1.0 if "eR" in process_name else 0.0
        process_p_pol = -1.0 if "pL" in process_name else 1.0 if "pR" in process_name else 0.0
        process_meta["e_pol"] = process_e_pol
        process_meta["p_pol"] = process_p_pol
        return process_meta


    def _add_file(self, path):
        # does all the heavy lifting
        super()._add_file(path)
        # still need to add our events to the count
        with ROOT.TFile(path) as file:
            tree = file["events"]
            process_name = self._parse_path(path)
            dataset = self._dataset[process_name]
            meta = dataset.metadata
            meta["n_events"] += tree.GetEntries()
            # print(process_name, vars(dataset))

    # TODO: what about the completeness check?
    # Implement this as an ILD specific method?