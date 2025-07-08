from .Dataset import Dataset
import ROOT
import re
from typing import Any
from .Backports import kP10
from collections import deque


class Analysis:
    """Holds a Dataset, all the data frames and possibly also all the histograms"""

    _dataset: Dataset
    _df: dict[str, Any] = {}
    _histograms = {}
    _varied_histograms = {}
    _sums = {}
    _stacks = {}
    _canvases = {}
    # need a place to "park" them somewhere as the THStacks do not take ownership :(
    _scaled_histograms = {}
    _legends = {}
    _reports = {}
    _categories: dict[str, list[str]] = {}
    _snapshots = {}
    _booked_objects: list[Any] = []
    _arrows = {}
    _lines = {}


    def __init__(self, dataset: Dataset):
        self._dataset = dataset
        for name, tree_name, files in dataset.get_samples():
            # filter out meta only
            if tree_name == "" and files == [""]:
                continue
            df = ROOT.RDataFrame(tree_name, files)
            self._df[name] = df


    def _get_frames(self, categories: list[str]|None):
        if not categories:
            for k, df in self._df.items():
                yield k, df
        else:
            for category_name in categories:
                category = self._categories[category_name]
                for df_name in category:
                    df = self._df[df_name]
                    yield df_name, df


    def _get_frame_names(self, categories: list[str]|None):
        if not categories:
            for k in self._df:
                yield k
        else:
            for category_name in categories:
                category = self._categories[category_name]
                for df_name in category:
                    yield df_name


    def _call_on_frames(self, method_name: str, args, categories: list[str]|None = None):
        for k, df in self._get_frames(categories):
            res = getattr(df, method_name)(*args)
            yield k, res


    def _call_on_frames_typed(self, method_name: str, typename: str, args, categories: list[str]|None = None):
        for k, df in self._get_frames(categories):
            res = getattr(df, method_name)[typename](*args)
            yield k, res


    def _apply_to_frames(self, method_name: str, args, categories: list[str]|None = None):
        for k, res_df in self._call_on_frames(method_name, args, categories):
            self._df[k] = res_df


    def _define(self, args, categories: list[str]|None = None):
        self._apply_to_frames("Define", args, categories)


    def _vary(self, args, categories: list[str]|None = None):
        self._apply_to_frames("Vary", args, categories)


    # refactor to use internal implementation above
    def Define(self, *args):
        for k, df in self._df.items():
            self._df[k] = df.Define(*args)


    # refactor to use internal implementation above
    def define_only_on(self, categories: list[str], *args):
        for category_name in categories:
            category = self._categories[category_name]
            for df_name in category:
                df = self._df[df_name]
                self._df[df_name] = df.Define(*args)


    def col_exists(self, col_name: str) -> bool:
        # take first df as all have the same
        cols = list(self._df.values())[0].GetColumnNames()
        return col_name in cols


    def init_categories(self):
        for k in self._df:
            category_name = self._dataset.get_category(k)
            # TODO: error handling if there is no category?
            try:
                self._categories[category_name].append(k)
            except KeyError:
                self._categories[category_name] = [k]


    def init_parameters(self, params: list[tuple[str, str, str]]):
        """Inits the podio generic parameters supplied as a list of (name, c++ typename, alternative value). If parameters were already initialised in a previous step, only the missing ones are added."""
        ROOT.gInterpreter.Declare("#include <podio/GenericParameters.h>")
        if not self.col_exists("Parameters"):
            self.Define("Parameters", "podio::GenericParameters par; par.loadFrom(GPDoubleKeys, GPDoubleValues); par.loadFrom(GPFloatKeys, GPFloatValues); par.loadFrom(GPIntKeys, GPIntValues); par.loadFrom(GPStringKeys, GPStringValues); return par;")
        for p_name, p_type, p_alt in params:
            par_col_name = f"params_{p_name.replace('.', '_')}"
            if not self.col_exists(par_col_name):
                self.Define(par_col_name, f"Parameters.get<{p_type}>(\"{p_name}\").value_or({p_alt})")


    def book_some_method(self, method_name: str, args, categories: list[str]|None = None, typename: str|None = None):
        results = {}
        if typename:
            it = self._call_on_frames_typed(method_name, typename, args, categories)
        else:
            it = self._call_on_frames(method_name, args, categories)
        for k, res in it:
            results[k] = res
            self._booked_objects.append(res)
        return results


    def book_sum(self, name: str, column: str, categories: list[str]|None = None, typename: str|None = None):
        self._sums[name] = self.book_some_method("Sum", (column,), categories, typename)


    def book_histogram_1D(self, name: str, column: str, config: tuple, categories: list[str]|None = None):
        self._histograms[name] = self.book_some_method("Histo1D", (config, column), categories)

    # TODO refactor
    def book_weighted_histogram_1D(self, name: str, column: str, weight: str, config: tuple, categories: list[str]|None = None):
        self._histograms[name] = self.book_some_method("Histo1D", (config, column, weight), categories)

    def book_histogram_2D(self, name: str, column1: str, column2: str, config: tuple, categories: list[str]|None = None):
        self._histograms[name] = self.book_some_method("Histo2D", (config, column1, column2), categories)


    def book_variations(self, name: str, categories: list[str]|None = None):
        res = {}
        histograms = self._histograms[name]
        for k in self._get_frame_names(categories):
            res[k] = ROOT.RDF.Experimental.VariationsFor(histograms[k])
        self._varied_histograms[name] = res


    def book_templates(self, configs: dict[str, tuple], weights: list[str], categories: list[str]|None = None):
        """Book template creation for the given named histogram configs and weight names. The first weight will be used as the nominal one. Expects weight names to start with 'weight_'"""
        nominal_name = weights[0]
        other_names = weights[1:]
        short_names = [name.removeprefix("weight_") for name in other_names]
        vary_args = (nominal_name, f"ROOT::RVecD{{{','.join(other_names)}}}", short_names, "weight")
        self._vary(vary_args, categories)
        for hist_name, hist_conf in configs.items():
            name = f"tmpl_{hist_name}"
            self.book_weighted_histogram_1D(name, hist_name, nominal_name, hist_conf, categories)
            self.book_variations(name, categories)


    def run(self):
        """Execute all booked computations"""
        ROOT.RDF.RunGraphs(self._booked_objects)


    # TODO: extend to allow to select categories
    # TODO: reconsider if this is needed
    def get_sum(self, name: str, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0) -> float:
        # ideally I could also make a functor that does this, but there is this
        # additional issue that histograms need to be cloned and numbers not etc.
        result = 0.
        sums = self._sums[name]
        for i, (category_name, dataframes) in enumerate(self._categories.items()):
            for k in dataframes:
                _sum = sums[k].GetValue()
                weight_key = k.removesuffix("_signal").removesuffix("_bkg")
                weight = self._dataset.get_weight(weight_key, int_lumi, e_pol, p_pol)
                result += _sum * weight
        return result


    # TODO: extend to allow to select categories
    # TODO: reconsider if this is needed
    def get_mean(self, name: str, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0) -> float:
        weighted_counts, errors2 = self._calc_cutflow(int_lumi, e_pol, p_pol)
        last_filter = list(weighted_counts.keys())[-1]
        counts = weighted_counts[last_filter]
        count = sum(counts)
        _sum = self.get_sum(name, int_lumi, e_pol, p_pol)
        return _sum / count


    def store_raw_histograms(self, variable_names: list[str], output_path: str):
        """Stores the unscaled histograms into a ROOT file, together with the meta information on the original integrated luminosity and beam polarisation"""
        with ROOT.TFile(f"{output_path}/raw_histograms.root", "recreate") as output_file:
            for category_name, dataframes in self._categories.items():
                cat_dir = output_file.mkdir(category_name)
                for k in dataframes:
                    dir = cat_dir.mkdir(k)
                    # add meta information
                    meta_dir = dir.mkdir("meta")
                    meta_dir.cd()
                    # FIXME
                    # weight_key = k.removesuffix("_signal").removesuffix("_bkg")
                    weight_key = k
                    lumi, e_pol, p_pol = self._dataset.get_lumi_and_pol(weight_key)
                    # FIXME: great now instead of a stupid error this even crashes the kernel!
                    lumi_par = ROOT.TParameter["float"]("lumi", lumi)
                    lumi_par.Write()
                    e_pol_par = ROOT.TParameter["float"]("e_pol", e_pol)
                    e_pol_par.Write()
                    p_pol_par = ROOT.TParameter["float"]("p_pol", p_pol)
                    p_pol_par.Write()
                    dir.cd()
                    # write histos
                    for var_name in variable_names:
                        h = self._histograms[var_name][k]
                        h.Write(var_name)


    # TODO: add also storage of raw histograms and a way to read them back in for mixing
    def store_histograms(self, variable_names: list[str], output_path: str, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0):
        with ROOT.TFile(f"{output_path}/histograms_{int_lumi}_{e_pol}_{p_pol}.root", "recreate") as output_file:
            for name in variable_names:
                params = (name, int_lumi, e_pol, p_pol)
                for cat_name, h in self._scaled_histograms[params].items():
                    h.Write(f"h_{name}_{cat_name}")



    def draw_histogram(self, name: str, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0, draw_opt: str = "hist", categories: list[str]|None = None, logY: bool = False, plot_dir: str|None = None, x_arrowl: float|None = None, x_arrowr: float|None = None):
        histograms = self._histograms[name]
        stack = ROOT.THStack()
        params = (name, int_lumi, e_pol, p_pol)
        self._scaled_histograms[params] = {}
        legend = ROOT.TLegend(0.6, 0.7, 1., 1,)
        # FIXME: put the calculation of the histograms into a separate method and only calculate them here if needed
        for i, (category_name, dataframes) in enumerate(self._categories.items()):
            if categories and category_name not in categories:
                # skip
                continue
            if len(dataframes) == 0:
                # nothing to do for empty categories
                continue
            scaled_cat_histograms = []
            for k in dataframes:
                # get histogram, clone it, scale it, put it in a list
                h = histograms[k].Clone()
                # our crossection table only knows about the full processes
                weight = self._dataset.get_weight(k, int_lumi, e_pol, p_pol)
                h.Scale(weight)
                scaled_cat_histograms.append(h)
            # sum up scaled category histograms, starting with the first one
            h = sum(scaled_cat_histograms[1:], start=scaled_cat_histograms[0])
            h.SetFillColor(kP10[i].GetNumber())
            # TODO: parse k to nice process name
            legend.AddEntry(h, category_name, "f")
            stack.Add(h)
            # store h so that it does not get deleted
            self._scaled_histograms[params][category_name] = h
        legend.SetNColumns(2)
        self._legends[params] = legend
        self._stacks[params] = stack
        canvas = ROOT.TCanvas()
        self._canvases[params] = canvas
        stack.SetTitle(f";{name};events")
        stack.Draw(draw_opt)
        legend.Draw()
        y_max = stack.GetMaximum()
        x_length = abs(stack.GetXaxis().GetXmax() - stack.GetXaxis().GetXmin())
        if x_arrowr:
            x1 = x_arrowr
            x2 = x1 - 0.05 * x_length
            y = 0.9 * y_max
            arrowr = ROOT.TArrow(x2, y, x1, y, 0.025, "<")
            arrowr.Draw()
            self._arrows[(params, "r")] = arrowr
            liner = ROOT.TLine(x1, 0, x1, y)
            liner.Draw()
            self._lines[(params, "r")] = liner
        if x_arrowl:
            x1 = x_arrowl
            x2 = x1 + 0.05 * x_length
            y = 0.9 * y_max
            arrowl = ROOT.TArrow(x1, y, x2, y, 0.025, ">")
            arrowl.Draw()
            self._arrows[(params, "l")] = arrowl
            linel = ROOT.TLine(x1, 0, x1, y)
            linel.Draw()
            self._lines[(params, "l")] = linel
        if logY:
            canvas.SetLogy()
        canvas.Draw()
        if plot_dir:
            canvas.SaveAs(f"{plot_dir}/{params}.pdf")
        return stack, legend, canvas


    def add_filter(self, expression: str, name: str):
        for k, df in self._df.items():
            self._df[k] = df.Filter(expression, name)


    def book_reports(self):
        for k, df in self._df.items():
            report = df.Report()
            self._reports[k] = report
            self._booked_objects.append(report)


    def _calc_cutflow(self, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0) -> tuple[dict[str, list[float]], dict[str, list[float]]]:
        # all processes should have the same *named* filters applied:
        names = list(list(self._df.values())[0].GetFilterNames())
        n_filters = len(names)
        # aggregate all the information
        numbers = {}
        errors2 = {}
        for category_name, frames in self._categories.items():
            nums = [0.0] * (n_filters + 1)
            errs2 = [0.0] * (n_filters + 1)
            # print(frames)
            for k in frames:
                reports = self._reports[k]
                weight_key = k
                # FIXME only do this step is needed!
                weight_key  = weight_key.removesuffix("_signal").removesuffix("_bkg")
                weight = self._dataset.get_weight(weight_key, int_lumi, e_pol, p_pol)
                for i, cut_info in enumerate(reports):
                    # print(f"processing cutinfo for category: {category_name}, frame: {k}, cut_name: {cut_info.GetName()}, pass: {cut_info.GetPass()}")
                    if i == 0:
                        n = cut_info.GetAll() * weight
                        nums[0] += n
                        errs2[0] += n * weight # need to square weight
                    n = cut_info.GetPass() * weight
                    nums[i+1] += n
                    errs2[i+1] += n * weight # need to square weight
            # print(nums)
            numbers[category_name] = nums
            errors2[category_name] = errs2
        return numbers, errors2


    def draw_cutflow(self, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0, plot_dir: str|None = None):
        numbers, errors2 = self._calc_cutflow(int_lumi, e_pol, p_pol)
        names = list(list(self._df.values())[0].GetFilterNames())
        n_filters = len(names)
        stack = ROOT.THStack()
        legend = ROOT.TLegend(0.6, 0.7, 1., 1,)
        for i, category_name in enumerate(self._categories):
            h = ROOT.TH1D("", "", n_filters+1, 0, n_filters+1)
            nums = numbers[category_name]
            for j, count in enumerate(nums):
                h.Fill(j, count)
            h.SetFillColor(kP10[i].GetNumber())
            legend.AddEntry(h, category_name, "f")
            stack.Add(h)
        name = "cut flow"
        params = (name, int_lumi, e_pol, p_pol)
        legend.SetNColumns(2)
        self._legends[params] = legend
        self._stacks[params] = stack
        canvas = ROOT.TCanvas()
        self._canvases[params] = canvas
        stack.SetTitle(";;events")
        stack.Draw("hist")
        x_axis = stack.GetXaxis()
        for i, name in enumerate(["All"] + names):
            x_axis.SetBinLabel(i+1, str(name))
        legend.Draw()
        canvas.SetLogy()
        canvas.Draw()
        if plot_dir:
            canvas.SaveAs(f"{plot_dir}/{params}.pdf")
        # canvas.SetPadBottomMargin(0.3)


    # TODO: adapt to use the new categories
    def print_reports(self, int_lumi: float = 5000, e_pol: float = 0.0, p_pol: float = 0.0):
        """Print cut-flow report for each category/prefix"""
        numbers, errors2 = self._calc_cutflow(int_lumi, e_pol, p_pol)
        largest = max([num for nums in numbers.values() for num in nums])
        names = list(list(self._df.values())[0].GetFilterNames())
        # do the actual printing
        # TODO: just also give a fixed format for the yields like Teresa
        from math import log10, ceil, sqrt
        max_digits = ceil(log10(largest+1)) #+ len(" 0e-00")
        max_name = max((len(name) for name in self._categories))
        offset = max(max_digits, max_name) + 1
        line = ""
        for category_name in self._categories:
            line_piece = f"{category_name: >{offset + len(' 0e-00')}}"
            line += line_piece
        print(line)
        for i, name in enumerate(["All"] + names):
            line = ""
            for category_name in self._categories:
                num = numbers[category_name][i]
                err = sqrt(errors2[category_name][i])
                line_piece = f"{round(num): >{offset}} ({err:.0e})"
                # print(f"::{line_piece}::")
                line += line_piece
            print(line, name)
        line = ""
        for category_name in self._categories:
            eff = numbers[category_name][-1] / numbers[category_name][0] if numbers[category_name][0] > 0 else 0
            if eff >= 0.01:
                line_piece = f"{eff: >{offset}.2f}"
            else:
                line_piece = f"{eff: >{offset}.0e}"
            line += line_piece
        print(line, "efficiency")


    def set_categories(self, input: dict[str, dict[str, str]]):
        """Defines categories that can be used that group multiple samples together or just use parts of a sample. The categories are expected to be disjoint. Two signal categories can not share the same processes."""
        # example_input = {
        # "4f_sw_sl_signal": {"pattern": "4f_sw_sl", "cut": "pre_defined_signal_cut"},
        # "4f_sl_bkg": {"pattern": r"4f\w+sl", "cut": None}, # inverse signal cut is applied automatically where necessary
        # "3f": {"pattern": "ea_3f|ae_3f", "cut": None} # FIXME: cutting not even supported at the moment, add a signal category if needed :(
        #}
        # TODO: The second part could also just be a named tuple
        # TODO: at the moment signals are supposed to end in '_signal' and backgrounds that need the inverse signal cut need to end in '_bkg' but cleaner would be to set a bool
        categories: dict[str, list[str]] = {}
        signals: list[str] = []
        for category_name, category_def in input.items():
            dataframes = []
            pattern = re.compile(category_def["pattern"])
            for df_name in filter(pattern.search, self._df.keys()):
                # print(category_name, df_name)
                dataframes.append(df_name)
            if category_name.endswith("_signal"):
                signals.append(category_name)
            categories[category_name] = dataframes
            # now transfer metadata
            # need to get metadata from before the signal definition cut for correct weight
            # FIXME: add the new samples here and also the category
            # old_frame = frame.removesuffix("_bkg").removesuffix("_signal")
            # *_, old_meta = self._dataset.get_sample(old_frame)
            # meta = old_meta.copy()
            # meta["category"] = category_name
            # new_sample = ([tree_name], [file_name], meta)
            # dataset.add_sample(frame, *new_sample)
        # now apply cuts and clean up the backgrounds
        for signal_name in signals:
            signal_def = input[signal_name]
            # need to apply cut to each dataframe of the signal
            cut = signal_def.get("cut")
            new_df_names = []

            # maybe a bit ugly, sorry :(
            def rename_and_cut(df_name, background=False):
                new_df = self._df[df_name]
                if cut:
                    if not background:
                        new_df = new_df.Filter(f"return {cut};")
                    else:
                        new_df = new_df.Filter(f"return !({cut});")
                new_df_name = f"{df_name}_signal" if not background else f"{df_name}_bkg"
                # print(f"creating {new_df_name} from {df_name} with background == {background}")
                # FIXME: also add explicitly to dataset to be able to look up metadata!
                self._df[new_df_name] = new_df
                return new_df_name

            for df_name in categories[signal_name]:
                new_df_names.append(rename_and_cut(df_name))
            # overwrite old df name list
            categories[signal_name] = new_df_names

            # now look for any backgrounds to rename and apply inverse cut
            pattern = re.compile(signal_def["pattern"])
            for category_name, dataframes in categories.items():
                if not category_name.endswith("_bkg"):
                    continue
                new_df_names = dataframes.copy()
                for df_name in filter(pattern.search, dataframes):
                    # print(f"signal {signal_name} has background {df_name}")
                    new_df_names.remove(df_name)
                    new_df_names.append(rename_and_cut(df_name, True))
                # overwrite old df name list
                categories[category_name] = new_df_names

        self._categories = categories


    def is_complete_categorisation(self):
        df_names = list(self._df.keys())
        categories = self._categories
        for c, frames in categories.items():
            for frame in frames:
                short_frame = frame.removesuffix("_bkg").removesuffix("_signal")
                if short_frame in df_names:
                    df_names.remove(short_frame)
                if frame in df_names:
                    df_names.remove(frame)
        return len(df_names) == 0


    def book_snapshots(self, tree_name: str, out_dir: str, meta_outname: str, column_list = None, no_rvec = False):
        """Will be written out when the event loop runs"""
        # Need to avoid double counting, write out all non-suffixed dataframes or start from the categories
        dataset = Dataset()
        categories = self._categories
        for category_name, frames in categories.items():
            for frame in frames:
                df = self._df[frame]
                file_name = f"{out_dir}/{category_name}_{frame}.snapshot.root"
                args = [tree_name, file_name]
                if column_list:
                    args.append(column_list)
                else:
                    # append empty column selection regex -> all columns
                    selection = ""
                    args.append(selection)
                snapshot_options = ROOT.RDF.RSnapshotOptions()
                snapshot_options.fLazy = True
                snapshot_options.fCompressionLevel = 9
                snapshot_options.fCompressionAlgorithm = ROOT.RCompressionSetting.EAlgorithm.kZSTD
                if no_rvec:
                    snapshot_options.fVector2RVec = False
                args.append(snapshot_options)
                snapshot = df.Snapshot(*args)
                self._snapshots[frame] = snapshot
                self._booked_objects.append(snapshot)
                # now transfer metadata
                # need to get metadata from before the signal definition cut for correct weight
                # FIXME: this has to be done as a change to the dataset at category creation time instead!!!
                # old_frame = frame.removesuffix("_bkg").removesuffix("_signal"
                old_frame = frame
                *_, old_meta = self._dataset.get_sample(old_frame)
                meta = old_meta.copy()
                meta["category"] = category_name
                new_sample = ([tree_name], [file_name], meta)
                dataset.add_sample(frame, *new_sample)
        with open(meta_outname, "w") as out_file:
            out_file.write(dataset.to_json(indent=2))


    def check_snapshots(self, tree_name: str, out_dir: str, meta_outname: str):
        # repeat some of the logic of book snapshots: create a dataset but put only working files inside
        # warn on broken files where the count after filtering is non-zero!
        dataset = Dataset()
        categories = self._categories
        for category_name, frames in categories.items():
            for frame in frames:
                file_name = f"{out_dir}/{category_name}_{frame}.snapshot.root"

                # TODO: check if file exists and contains the tree and has as many entries as expected
                # first: check expected count
                # we can take them from the cut reports
                report = self._reports[frame]
                # get the last cut
                last_cut = deque(report, maxlen=1).pop()
                expected_count = last_cut.GetPass()
                try:
                    with ROOT.TFile.Open(file_name) as file:
                        try:
                            tree = file[tree_name]
                            written_count = tree.GetEntries()
                            if written_count != expected_count:
                                print(f"ERROR: unequal counts! written: {written_count} expected: {expected_count}")
                                continue
                        except KeyError:
                            if expected_count > 0:
                                print(f"ERROR: missing tree {tree_name} in {file_name}")
                            else:
                                print(f"INFO: missing tree {tree_name} but expected_count is zero in {file_name}")
                            continue
                except OSError:
                    if expected_count > 0:
                        print(f"ERROR: missing file {file_name} with {expected_count} expected events")
                    else:
                        print(f"INFO: missing empty file {file_name}")
                    continue

                # now transfer metadata
                # need to get metadata from before the signal definition cut for correct weight
                old_frame = frame.removesuffix("_bkg").removesuffix("_signal")
                *_, old_meta = self._dataset.get_sample(old_frame)
                meta = old_meta.copy()
                meta["category"] = category_name
                new_sample = ([tree_name], [file_name], meta)
                dataset.add_sample(frame, *new_sample)
        with open(meta_outname, "w") as out_file:
            out_file.write(dataset.to_json(indent=2))
