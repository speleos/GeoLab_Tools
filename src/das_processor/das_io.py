"""HDF5 / DAS file I/O: DictFile, DASDataFrame, file discovery and loading."""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import datetime
import logging
import os
import re
import warnings
from collections.abc import Iterable

# --- Third party -------------------------------------------------------------
import numpy as np
import pandas as pd
from h5py import Dataset, File, Group
from sympy import symbols, sympify

# --- Local -------------------------------------------------------------------
from ._log import _ensure_console_handler

logger = logging.getLogger(__name__)


class DictFile(File):
    """Read and write HDF5 files as nested Python dictionaries."""

    def __init__(self, filename, mode, verbose=False, **kwargs):
        self.verbose = verbose
        if verbose:
            _ensure_console_handler()
            logger.info("Open file: " + str(filename))
        super().__init__(filename, mode, **kwargs)

    def save_dict(self, data):
        """Save a nested dictionary tree to this HDF5 file."""
        if self.verbose:
            logger.info("Save dictionary tree to: " + str(self.filename))
        self._save_tree(data, self, 0)

    def load_dict(self, field=None, skipFields=None, getData=True):
        """Load a dictionary tree, optionally excluding fields or datasets."""
        if self.verbose:
            line = "Load dictionary from "
            if field:
                line += f'field "{field}"'
            line += " within " + str(self.filename) + ":"
            logger.info("%s", line)

        self.skipFields = [] if skipFields is None else skipFields
        if isinstance(field, str):
            field = self[field]
        if isinstance(field, Dataset):
            return field[()]

        data = {}
        self._load_tree(self if field is None else field, data, 0, getData)
        return data

    def load_keys(self):
        """Load the HDF5 tree with dataset shapes in place of values."""
        if self.verbose:
            logger.info("Loading tree structure: " + str(self.filename))
        data = {}
        self.skipFields = []
        self._load_tree(self, data, 0, False)
        return data

    def _save_tree(self, tree, parent_group, depth):
        if tree is None or len(tree) == 0:
            if self.verbose:
                logger.info("%s %s", "\t" * depth, "-")
            return

        for key, value in tree.items():
            if isinstance(value, dict):
                if self.verbose:
                    logger.info("\t" * depth + str(key))
                child_group = parent_group.create_group(str(key))
                self._save_tree(value, child_group, depth + 1)
            else:
                if not isinstance(value, np.ndarray):
                    value = np.array(value)
                try:
                    parent_group.create_dataset(str(key), data=value)
                except TypeError:
                    if not isinstance(value, np.ndarray):
                        raise
                    parent_group.create_dataset(
                        str(key), data=np.array(value, dtype="S")
                    )

            if self.verbose:
                self._print_info(str(key), value, depth)

    def _load_tree(self, tree, parent_dict, depth, get_data):
        if tree is None or not isinstance(tree, Group):
            if self.verbose:
                logger.info("%s %s", "\t" * depth, "-")
            return

        for key, value in tree.items():
            if key in self.skipFields:
                continue
            key = int(key) if isinstance(key, str) and key.isdigit() else str(key)
            if isinstance(value, Group):
                if self.verbose:
                    logger.info("\t" * depth + key)
                parent_dict[key] = {}
                self._load_tree(value, parent_dict[key], depth + 1, get_data)
                continue

            if get_data is True or (
                isinstance(get_data, list) and key in get_data
            ):
                if value.dtype.char in "SO":
                    if value.shape == ():
                        raw_value = value[()]
                        if raw_value == b"None":
                            parent_dict[key] = None
                        else:
                            try:
                                parent_dict[key] = raw_value.decode("utf-8")
                            except UnicodeDecodeError:
                                parent_dict[key] = raw_value
                    else:
                        parent_dict[key] = np.array(value).astype(str)
                elif value.shape == ():
                    parent_dict[key] = value[()].item()
                else:
                    parent_dict[key] = np.array(value[()])
            else:
                parent_dict[key] = value.shape

            if self.verbose:
                self._print_info(key, value, depth)

    def _print_info(self, key, value, depth, maxlength=16):
        prefix = "\t" * depth + ("{} ({}):".format(key, str(value.dtype))).ljust(30)
        if len(value.shape) == 0 or len(value) <= maxlength:
            logger.info(prefix + str(value[()]))
        else:
            logger.info(prefix + str(value.shape))


def save(filename, data, **kwargs):
    """Save a nested dictionary to an HDF5 file."""
    with DictFile(filename, "w", **kwargs) as hdf_file:
        hdf_file.save_dict(data)


def load(filename, field=None, skipFields=None, **kwargs):
    """Load a nested dictionary from an HDF5 file."""
    with DictFile(filename, "r", **kwargs) as hdf_file:
        return hdf_file.load_dict(field=field, skipFields=skipFields)


def load_keys(filename, **kwargs):
    """Load an HDF5 tree with shapes in place of dataset values."""
    with DictFile(filename, "r", **kwargs) as hdf_file:
        return hdf_file.load_keys()


class DASDataFrame(pd.DataFrame):
    """Pandas DataFrame carrying DAS metadata."""

    _metadata = ["meta"]

    @property
    def _constructor(self):
        return DASDataFrame

    def __init__(
        self,
        data: np.ndarray | Iterable | dict | pd.DataFrame = None,
        index=None,
        columns=None,
        dtype=None,
        copy=None,
        meta: dict | None = None,
    ):
        super().__init__(data, index, columns, dtype, copy)
        self.meta = {} if meta is None else meta.copy()
        if self.meta:
            self._set_dimensionRange()

    def _set_dimensionRange(self):
        """Compute dimension ranges and sizes from the dataframe axes."""
        if np.issubdtype(self.columns.dtype, np.str_) or np.issubdtype(
            self.columns.dtype, object
        ):
            chs = np.arange(len(self.columns))
        else:
            chs = np.array(self.columns)

        if len(chs) > 2:
            dchs = np.diff(chs)
            sdchs = dchs[:-1] != dchs[1:]
            sdchs &= np.r_[sdchs[1:], False]
            stopind = np.r_[np.where(sdchs)[0] + 1, len(chs) - 1]
            startind = np.r_[0, stopind[:-1] + 1]
            size = np.array(
                [len(chs[start:end]) + 1 for start, end in zip(startind, stopind)]
            )
        else:
            startind = 0
            stopind = -1
            size = len(chs)

        self.meta["dimensionRanges"] = {
            "dimension0": {
                "max": self.shape[0] - 1,
                "min": 0,
                "size": self.shape[0],
                "unitScale": self.meta["dt"],
            },
            "dimension1": {
                "max": chs[stopind],
                "min": chs[startind],
                "size": size,
                "unitScale": self.meta["dx"],
            },
        }
        self.meta.update(
            dimensionSizes=np.r_[self.shape[0], np.sum(size)],
            dimensionUnits=np.array([b"s", b"m"]),
        )


def _combine_units(units: list[str], operator: str = "*") -> str:
    """Combine unit strings using multiplication or division."""
    known_units = (
        "rad",
        "dt",
        "strain",
        "W",
        "NormCts",
        "m",
        "K",
        "s",
        "Hz",
        "N",
        "Pa",
        "J",
    )
    unit_symbols = symbols(known_units)
    unit_dict = dict(zip(known_units, unit_symbols))
    units = [unit.replace("·", "*").replace("ε", "strain") for unit in units]
    combined_unit = sympify(units[0], unit_dict)
    for unit in units[1:]:
        if unit:
            if operator == "/":
                combined_unit /= sympify(unit, unit_dict)
            elif operator == "*":
                combined_unit *= sympify(unit, unit_dict)
    return str(combined_unit)


def _set_sensitivity(
    header: dict, sensitivitySelect: int, userSensitivity: None | dict = None
) -> tuple[float, str, float, str]:
    if sensitivitySelect >= 0:
        try:
            sensitivity = header["sensitivities"][sensitivitySelect, :]
            sensitivity_unit = header["sensitivityUnits"][sensitivitySelect]
        except KeyError:
            sensitivity = header["sensitivity"]
            sensitivity_unit = header["sensitivityUnit"]
        except IndexError:
            logger.warning("SensitivitySelect index %d not found in file." % sensitivitySelect)
            raise

        unit_out = _combine_units([header["unit"], sensitivity_unit], "/")
        sensitivities_out = np.ones((1, 1), dtype=np.float32)
        sensitivity_units_out = [""]
    elif sensitivitySelect == -1:
        unit_out = header["unit"]
        sensitivity = 1.0
        sensitivities_out = header["sensitivities"]
        sensitivity_units_out = header["sensitivityUnits"]
    elif sensitivitySelect == -2:
        unit_new = _combine_units([header["unit"], "m"])
        if not any(unit == "m" for unit in re.findall(r"[\w']+", unit_new)):
            unit_out = unit_new
            sensitivity = 1.0 / header["gaugeLength"]
            sensitivities_out = header["sensitivities"] * header["gaugeLength"]
            sensitivity_units_out = [
                _combine_units([unit, "m"]) for unit in header["sensitivityUnits"]
            ]
        else:
            unit_out = header["unit"]
            sensitivity = 1.0
            sensitivities_out = header["sensitivities"]
            sensitivity_units_out = header["sensitivityUnits"]
    elif sensitivitySelect == -3 and userSensitivity is not None:
        unit_out = _combine_units(
            [header["unit"], userSensitivity["sensitivityUnit"]], "/"
        )
        sensitivity = userSensitivity["sensitivity"]
        sensitivities_out = np.ones((1, 1), dtype=np.float32)
        sensitivity_units_out = [""]
    else:
        raise ValueError("Undefined sensitivitySelect")

    return sensitivity, unit_out, sensitivities_out, sensitivity_units_out


def _fix_meta_back_compability(meta: dict, nSamples: int, nChannels: int):
    """Fill metadata fields missing from older DAS file versions."""
    speed_of_light = 2.99792458e8
    header = meta["header"]

    if "experiment" not in header and "exp" in header:
        header["experiment"] = header.pop("exp")

    if "cableSpec" not in meta:
        meta["cableSpec"] = {
            "fiberOverLength": 1.0,
            "refractiveIndex": 1.4677,
            "zeta": 0.78,
        }

    if "dx" not in header or "gaugeLength" not in header:
        dx_fiber = (
            meta["demodSpec"]["dTau"]
            * speed_of_light
            / (2 * meta["cableSpec"]["refractiveIndex"])
        )
        if "dx" not in header:
            header["dx"] = dx_fiber / meta["cableSpec"]["fiberOverLength"]
        if "gaugeLength" not in header:
            header["gaugeLength"] = meta["demodSpec"]["nDiffTau"] * dx_fiber

    if "dataScale" not in header:
        header["dataScale"] = (
            np.pi / 2**29 / header["dt"] / header["gaugeLength"]
        )
        header["unit"] = "rad/m/s"
    if "spatialUnwrRange" not in header:
        header["spatialUnwrRange"] = (
            8 * np.pi / header["dt"] / header["gaugeLength"]
        )

    if "dimensionRanges" not in header:
        dimension_ranges = {}
        dimension_ranges["dimension0"] = {
            "min": 0,
            "max": nSamples - 1,
            "size": nSamples,
            "unitScale": header["dt"],
        }
        try:
            processing_chain = (
                meta["processingChain"]
                if "processingChain" in meta
                else meta["dimChain"]
            )
            step_number = 0
            while "step-%d" % step_number in processing_chain:
                step_number += 1
            columns = processing_chain["step-%d" % (step_number - 1)]["cols"]
        except (KeyError, TypeError):
            columns = meta["demodSpec"]

        dimension_ranges["dimension1"] = {
            "min": columns["roiStart"],
            "max": columns["roiEnd"],
            "size": (
                columns["roiEnd"] - columns["roiStart"] + 1
            )
            // columns["roiDec"],
            "unitScale": header["dx"],
        }
        header["dimensionRanges"] = dimension_ranges
        header.update(
            dimensionSizes=np.r_[nSamples, nChannels],
            dimensionNames=["time", "distance"],
            dimensionUnits=["s", "m"],
            name="Phase rate per distance",
        )

    dimension = header["dimensionRanges"]["dimension1"]
    if not isinstance(dimension["min"], Iterable):
        dimension.update({key: np.r_[value] for key, value in dimension.items()})

    if "channels" not in header:
        roi_channels = []
        for start, stop, size in zip(
            dimension["min"], dimension["max"], dimension["size"]
        ):
            roi_channels.append(
                np.arange(start, stop + 1, (stop - start + 1) // size)
            )
        header["channels"] = np.sort(np.unique(np.concatenate(roi_channels)))

    try:
        itu = int(meta["monitoring"]["Laser"]["itu"])
    except (KeyError, TypeError, ValueError):
        itu = 51
    wavelength = speed_of_light / (190e12 + itu * 1e11)

    if "sensitivity" in header:
        header["sensitivities"] = np.atleast_2d(header.pop("sensitivity"))
        header["sensitivityUnits"] = np.atleast_1d(header.pop("sensitivityUnit"))

    if "sensitivities" not in header or "NA" in header["sensitivityUnits"]:
        header["sensitivities"] = np.atleast_2d(
            4
            * np.pi
            * meta["cableSpec"]["zeta"]
            * meta["cableSpec"]["refractiveIndex"]
            / wavelength
        )
        header["sensitivityUnits"] = np.atleast_1d("rad/m/strain")
    elif len(header["sensitivityUnits"]) and "unit" in header:
        normalized_units = []
        for sensitivity_unit in header["sensitivityUnits"].copy():
            output_unit = _combine_units([header["unit"], sensitivity_unit], "/")
            if "rad" in re.findall(r"[\w']+", output_unit):
                sensitivity_unit = output_unit
            normalized_units.append(
                sensitivity_unit.replace("·", "*").replace("ε", "strain")
            )
        header["sensitivityUnits"] = np.array(normalized_units)

    if "unit" in header:
        header["unit"] = header["unit"].replace("·", "*").replace("ε", "strain")
    if "sensorType" not in header:
        header["sensorType"] = ["D"]


def get_data_indexes(
    ffidpath: str | list[str],
    trange: None | list[float] | np.ndarray = None,
    channels: None | slice | list[int] | np.ndarray = None,
    show_header_info: bool = True,
):
    """Compute DAS sample and channel indexes from time and channel selections."""
    if isinstance(ffidpath, list):
        ffidpath = ffidpath[0]
    with DictFile(ffidpath, "r") as das_file:
        metadata = das_file.load_dict(skipFields=["data"])
        n_samples, n_channels = das_file["data"].shape
        _fix_meta_back_compability(metadata, n_samples, n_channels)

    header = metadata["header"]
    header_time = datetime.datetime.utcfromtimestamp(float(header["time"]))
    if trange is None:
        samples = slice(None)
    else:
        if isinstance(trange, (int, float)):
            trange = (0.0, trange)
        elif len(trange) == 1:
            trange = (0.0, trange[0])
        samples_to_read = int(np.diff(trange)[0] / header["dt"] + 0.5)
        sample_start = int(trange[0] / header["dt"] + 0.5)
        samples = slice(sample_start, sample_start + samples_to_read)

    if show_header_info:
        _ensure_console_handler()
        logger.info("-- Header info file: %s --" % os.path.basename(ffidpath))
        logger.info("\tExperiment:            %s" % header["experiment"])
        logger.info("\tFile timestamp:        %s" % header_time.strftime("%Y-%m-%d %H:%M:%S"))
        logger.info("\tType of data:          {}, unit: {}".format(
                header["name"], header["unit"]
            ))
        logger.info("\tSampling frequency:    %.2f Hz" % (1.0 / header["dt"]))
        logger.info("\tData shape:            %d samples x %d channels"
            % (n_samples, n_channels))
        logger.info("\tGauge length:          %.1f m" % header["gaugeLength"])
        logger.info("\tSensitivities:         %s"
            % ",".join(
                "{:.2e} {}".format(sensitivity, unit)
                for sensitivity, unit in zip(
                    header["sensitivities"][:, 0], header["sensitivityUnits"]
                )
            ))
        try:
            dimension = header["dimensionRanges"]["dimension1"]
            regions = ",".join(
                "%d:%d:%d"
                % (
                    start,
                    stop,
                    1 if stop <= start else (stop - start) // (size - 1),
                )
                for start, stop, size in zip(
                    dimension["min"], dimension["max"], dimension["size"]
                )
            )
            logger.info("\tRegions of interest:   %s" % regions)
        except (KeyError, TypeError, ZeroDivisionError):
            pass

    if isinstance(channels, slice):
        channels = np.arange(header["channels"][-1])[channels]

    header_channels = header["channels"]
    if channels is not None:
        channels_found, ch_index, _ = np.intersect1d(
            header_channels, channels, return_indices=True
        )
        if show_header_info:
            with np.printoptions(threshold=20):
                logger.info("\t%d channels requested:%s" % (len(channels), channels))
                logger.info("\t%d channels found:    %s" % (len(channels_found), channels_found))
    else:
        ch_index = np.arange(len(header_channels))
    if show_header_info:
        logger.info("-----------------------------------")

    if len(ch_index) == 0:
        raise ValueError("Selected channels are not in found in file!")
    return ch_index, samples


def get_filemeta(filepath: str, metaDetail: int = 1):
    """Read file metadata, optionally limiting output to interpretation fields."""
    with DictFile(filepath, "r") as das_file:
        metadata = das_file.load_dict(skipFields=["data"])
    if metaDetail != 1:
        return metadata

    monitoring = metadata["monitoring"]
    return {
        "fileVersion": metadata["fileVersion"],
        "header": metadata["header"],
        "timing": metadata["timing"],
        "cableSpec": metadata["cableSpec"],
        "monitoring": {
            "Gps": monitoring["Gps"],
            "Laser": {"itu": monitoring["Laser"]["itu"]},
        },
    }


def save_to_DAS_file(
    dfdas: DASDataFrame,
    filepath: None | str = None,
    datatype: str = "processed",
    **auxmeta,
) -> tuple[str, dict]:
    """Save a DAS dataframe in a format readable by :func:`load_DAS_files`."""
    if filepath is None:
        timestr = dfdas.meta["time"].strftime("%H%M%S")
        filename = timestr + ".hdf5"
        original_filepath = dfdas.meta["filepaths"][0]
        original_directory, _ = os.path.split(original_filepath)
        experiment_path, _ = os.path.split(original_directory)
        filepath = os.path.join(experiment_path, datatype, filename)

    directory, _ = os.path.split(filepath)
    if not os.path.exists(directory):
        os.makedirs(directory)
    logger.info("Saving to file: " + filepath)

    timestamp = dfdas.meta["time"].replace(
        tzinfo=datetime.timezone.utc
    ).timestamp()
    dfdas._set_dimensionRange()
    header = {
        "channels": dfdas.columns.values,
        "dataScale": np.float32(1.0),
        "spatialUnwrRange": np.float32(0.0),
    }
    header.update(**dfdas.meta)
    header.update(time=timestamp)

    file_version = header.pop("fileVersion")
    ffid_dict = {"header": {}, "fileVersion": file_version}
    ffid_dict.update(data=np.asarray(dfdas, dtype=np.float32), **auxmeta)
    ffid_dict["header"].update(**header)

    save(filepath, ffid_dict)
    return filepath, ffid_dict


def str2datetime(sdt: str, dateonly: bool = False) -> datetime.datetime:
    """Convert DAS time strings to datetime objects."""
    if isinstance(sdt, datetime.datetime):
        return sdt
    if dateonly:
        return datetime.datetime.strptime(sdt, "%Y%m%d")
    try:
        return datetime.datetime.strptime(sdt, "%H%M%S")
    except ValueError:
        return datetime.datetime.strptime(sdt, "%Y%m%d %H%M%S")


def convert_index_to_rel_time(dfdas: DASDataFrame):
    """Replace a datetime index with seconds relative to its first sample."""
    dfdas["t"] = (dfdas.index - dfdas.index[0]).total_seconds()
    dfdas.set_index("t", inplace=True)


def create_time_axis(
    tstart: datetime.datetime, nSamples: int, dt: float
) -> pd.DatetimeIndex:
    """Create a datetime index for DAS samples at interval ``dt`` seconds."""
    return tstart + pd.timedelta_range(
        start="0s", periods=nSamples, freq="%dns" % int(dt * 1e9 + 0.5)
    )


def unwrap(
    phi: np.ndarray, wrapStep: float = 2 * np.pi, axis: int = -1
) -> np.ndarray:
    """Unwrap phase values with the specified wrap interval."""
    scale = 2 * np.pi / wrapStep
    return (np.unwrap(phi * scale, axis=axis) / scale).astype(phi.dtype)


def wrap(x: np.ndarray, wrapStep: float = 2 * np.pi) -> np.ndarray:
    """Wrap values to ``[-wrapStep/2, wrapStep/2)``."""
    if wrapStep > 0:
        return (x + wrapStep / 2) % wrapStep - wrapStep / 2
    return x


def find_DAS_files(
    experiment_path: str,
    start: str | datetime.datetime,
    duration: float | datetime.timedelta,
    channels: None | slice | list[int] | np.ndarray = None,
    datatype: str = "dphi",
    show_header_info: bool = True,
    load_file_from_start: bool = True,
) -> tuple[list[str], np.ndarray | None, slice | None]:
    """Find DAS files overlapping a time interval and select data indexes."""

    def checkdate(date_string):
        return datetime.datetime.strptime(date_string, "%Y%m%d").strftime("%Y%m%d")

    exper_path = experiment_path
    dates_in_exp = []
    exppath_has_dates = False
    for date_entry in sorted(os.listdir(exper_path)):
        if os.path.isdir(os.path.join(exper_path, date_entry)):
            try:
                dates_in_exp.append(checkdate(date_entry))
                exppath_has_dates = True
            except ValueError:
                pass

    if not exppath_has_dates:
        path = exper_path
        previous_component = ""
        for _ in range(2):
            path, component = os.path.split(path)
            try:
                dates_in_exp.append(checkdate(component))
                exppath_has_dates = True
                exper_path = path
                if previous_component:
                    datatype = previous_component
                break
            except ValueError:
                previous_component = component
        if not dates_in_exp:
            dates_in_exp.append("19000101")

    if isinstance(start, str):
        datetime_start = str2datetime(start)
    else:
        datetime_start = start
    if datetime_start.year == 1900 and exppath_has_dates:
        datetime_start += str2datetime(dates_in_exp[0], True) - datetime.datetime(
            1900, 1, 1
        )

    if not isinstance(duration, datetime.timedelta):
        duration = datetime.timedelta(seconds=duration)
    datetime_stop = datetime_start + duration

    filepaths = []
    header_times = []
    for date_in_exp in dates_in_exp:
        date = str2datetime(date_in_exp, True)
        if datetime_start.date() <= date.date() <= datetime_stop.date():
            data_directory = (
                os.path.join(exper_path, date_in_exp, datatype)
                if exppath_has_dates
                else exper_path
            )
            for filename in sorted(os.listdir(data_directory)):
                try:
                    file_time = str2datetime(
                        date_in_exp + " " + os.path.splitext(filename)[0]
                    )
                except ValueError:
                    continue

                if not (
                    datetime_start - datetime.timedelta(seconds=11)
                    < file_time
                    <= datetime_stop + datetime.timedelta(seconds=1)
                ):
                    continue

                filepath = os.path.join(data_directory, filename)
                with DictFile(filepath, "r") as das_file:
                    header_time = datetime.datetime.utcfromtimestamp(
                        float(das_file["header"]["time"][()])
                    )
                if datetime_start.year == 1900:
                    header_time += datetime.datetime(1900, 1, 1) - datetime.datetime.combine(
                        header_time.date(), datetime.time(0)
                    )
                if (
                    datetime_start - datetime.timedelta(seconds=10)
                    < header_time
                    <= datetime_stop
                ):
                    filepaths.append(filepath)
                    header_times.append(header_time)

    if filepaths:
        start_time = (
            0.0
            if load_file_from_start
            else max(0.0, (datetime_start - header_times[0]).total_seconds())
        )
        ch_index, samples = get_data_indexes(
            filepaths,
            (start_time, start_time + duration.total_seconds()),
            channels,
            show_header_info,
        )
    else:
        ch_index = None
        samples = None
        warnings.warn(
            "DAS_processor.find_DAS_files did not find any files", UserWarning
        )

    return filepaths, ch_index, samples


def load_DAS_files(
    filepaths: str | list[str],
    chIndex: None | slice | list[int] | np.ndarray = None,
    samples: None | slice | int = None,
    sensitivitySelect: int = 0,
    integrate: bool = True,
    unwr: bool = False,
    spikeThr: None | float = None,
    userSensitivity: None | dict = None,
    useLabeledColumns: None | bool = None,
    verbose: bool = False,
) -> DASDataFrame:
    """Load OptoDAS HDF5 files into a metadata-aware dataframe."""
    loadstart = datetime.datetime.now()

    if isinstance(filepaths, str):
        filepaths = [filepaths]
    if filepaths is None or len(filepaths) == 0:
        raise ValueError("No filepaths provided")

    if isinstance(samples, range):
        samples = slice(samples.start, samples.stop, samples.step)
    if isinstance(samples, int):
        samples = slice(0, samples)
    elif samples is None:
        samples = slice(None)
    if samples.step is not None and samples.step > 1:
        raise ValueError("Sample step>1 in time is not supported")

    if chIndex is None:
        chIndex = slice(None)
    elif not isinstance(chIndex, slice):
        chIndex = np.atleast_1d(chIndex)
        if len(chIndex) == 0:
            raise ValueError("Selected channels are not in found in file!")
        if len(chIndex) == 1:
            chIndex = slice(chIndex[0], chIndex[0] + 1, 1)
        else:
            channel_steps = np.diff(chIndex)
            if np.all(channel_steps == channel_steps[0]):
                chIndex = slice(chIndex[0], chIndex[-1] + 1, channel_steps[0])

    samples_read = 0
    if verbose:
        _ensure_console_handler()
        logger.info("Loading %d file(s)" % len(filepaths))

    for file_number, filepath in enumerate(filepaths):
        if not os.path.exists(filepath):
            raise FileNotFoundError(f"{filepath} not found")
        if verbose:
            logger.info("  " + os.path.basename(filepath))

        with DictFile(filepath, "r") as das_file:
            if das_file["data"].ndim > 2:
                raise IndexError(
                    "simpleDAS does not handle data with more than 2 dimensions"
                )
            samples_in_file, channels_in_file = das_file["data"].shape

            if file_number == 0:
                channels_out = np.arange(channels_in_file)[chIndex].size
                samples_max = (samples_in_file + 1) * len(filepaths)
                samples_out = np.arange(samples_max)[samples].size
                samples_remaining = samples_out
                data = np.zeros((samples_out, channels_out), dtype=np.float32)

                metadata = das_file.load_dict(skipFields=["data"])
                _fix_meta_back_compability(
                    metadata, samples_in_file, channels_in_file
                )

                if metadata["header"]["channels"][chIndex].size == 0:
                    raise IndexError(
                        f"chIndexes {chIndex} not in data "
                        f"(nChannels = {channels_in_file})"
                    )

                (
                    sensitivity,
                    unit_out,
                    sensitivities_out,
                    sensitivity_units_out,
                ) = _set_sensitivity(
                    metadata["header"], sensitivitySelect, userSensitivity
                )
                if not isinstance(sensitivity, np.ndarray) or sensitivity.size == 1:
                    if sensitivity == 0.0:
                        raise ValueError(
                            "Sensitivity value set to zero, use "
                            "sensitivitySelect=-1 or -2 to avoid error"
                        )
                    scale = np.float32(
                        metadata["header"]["dataScale"] / sensitivity
                    )
                else:
                    if channels_in_file != len(sensitivity):
                        raise ValueError(
                            "Length of sensitivity vector does not match "
                            "number of channels"
                        )
                    if not np.any(sensitivity != 0.0):
                        raise ValueError(
                            "Sensitivity value set to zero, use "
                            "sensitivitySelect=-1 or -2 to avoid error"
                        )
                    sensitivity = np.atleast_2d(sensitivity[chIndex])
                    scale = np.array(
                        metadata["header"]["dataScale"] / sensitivity,
                        dtype=np.float32,
                    )

                if samples.stop is None or samples.stop > samples_in_file:
                    samples_to_read = slice(samples.start, None)
                else:
                    samples_to_read = slice(samples.start, samples.stop)
            elif samples_remaining >= samples_in_file:
                samples_to_read = slice(None)
            else:
                samples_to_read = slice(0, samples_remaining)

            try:
                sub_data = das_file["data"][samples_to_read, chIndex]
            except TypeError:
                warnings.warn(
                    "h5py does not accept complex slices. Slow read...", UserWarning
                )
                sub_data = das_file["data"][:, chIndex][samples_to_read, :]

            samples_read_in_file = sub_data.shape[0]
            data[samples_read : samples_read + samples_read_in_file, :] = (
                np.atleast_2d(sub_data) * scale
            )
            samples_read += samples_read_in_file
            samples_remaining = samples_out - samples_read

    if verbose:
        logger.info("Read in %.1f s; conditioning..." % (datetime.datetime.now() - loadstart).total_seconds())
    data = data[:samples_read, :]

    if (unwr or spikeThr or integrate) and metadata["header"]["dataType"] == 2:
        unwr, spikeThr, integrate = (False,) * 3
        warnings.warn(
            "Options unwr, spikeThr or integrate can only be used with "
            "time differentiated phase data",
            UserWarning,
        )

    if unwr and metadata["header"]["spatialUnwrRange"]:
        data = unwrap(
            data, metadata["header"]["spatialUnwrRange"] / sensitivity, axis=1
        )
    if spikeThr is not None:
        data[np.abs(data) > spikeThr / sensitivity] = 0

    if integrate:
        unit_new = _combine_units([unit_out, "s"])
        if not any(unit == "s" for unit in re.findall(r"[\w']+", unit_new)):
            data = np.cumsum(data, axis=0) * metadata["header"]["dt"]
            unit_out = unit_new
        else:
            warnings.warn(
                "Data unit %s is not differentiated. Integration skipped." % unit_out,
                UserWarning,
            )
        metadata["header"]["name"] = metadata["header"]["name"].replace(
            " rate", ""
        )

    sample_start = 0 if samples.start is None else samples.start
    time_start = datetime.datetime.utcfromtimestamp(
        metadata["header"]["time"] + sample_start * metadata["header"]["dt"]
    )
    n_samples, _ = data.shape
    time_axis = create_time_axis(time_start, n_samples, metadata["header"]["dt"])

    meta = {
        key: metadata["header"][key]
        for key in (
            "dt",
            "dx",
            "gaugeLength",
            "experiment",
            "dataType",
            "dimensionRanges",
            "dimensionUnits",
            "dimensionNames",
            "name",
            "sensorType",
        )
    }
    meta.update(
        fileVersion=metadata["fileVersion"],
        time=time_start,
        unit=unit_out,
        sensitivities=sensitivities_out,
        sensitivityUnits=sensitivity_units_out,
        filepaths=filepaths,
    )

    if useLabeledColumns is not False:
        if len(metadata["header"]["sensorType"]) != channels_in_file:
            if useLabeledColumns is None:
                columns = metadata["header"]["channels"][chIndex]
            else:
                raise ValueError("Cannot use labeled columns")
        else:
            columns = metadata["header"]["sensorType"][chIndex]
    else:
        columns = metadata["header"]["channels"][chIndex]

    dfdas = DASDataFrame(data, columns=columns, index=time_axis, meta=meta)
    if verbose:
        logger.info("Total load time: %.1f s"
            % (datetime.datetime.now() - loadstart).total_seconds())
    return dfdas
