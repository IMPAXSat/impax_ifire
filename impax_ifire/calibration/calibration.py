"""
A module for all things calibration.
"""

import os
import random
import sys
from copy import deepcopy
from encodings.punycode import T
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.time import Time
from astropy.timeseries import TimeSeries
from spacepy import pycdf
from swxsoc.swxdata import SWXData

from impax_ifire import log
from impax_ifire.util.telemetry_config import telemetry_config

l1_depot = os.path.join(os.path.dirname(__file__), "..", "data", "log")

__all__ = [
    "calibrate_file",
    "get_calibration_file",
    "process_file",
    "read_calibration_file",
]


def process_file(data_filename: Path) -> list[Path | None]:
    """
    This is the entry point for the pipeline processing.
    It runs all of the various processing steps required.

    Parameters
    ----------
    data_filename: Path
        Fully specificied filename of an input file

    Returns
    -------
    output_files: List[Union[Path, None]]
        Fully specificied filenames for the output files.

    Notes
    -----
    If you want to have the pipeline mark files as processes, but don't want to return an output file,
    return a list with None element.
    If you don't want to mark the file as processed, return an empty list.
    """
    log.info(f"Processing file {data_filename}.")
    output_files = []

    calibrated_file = calibrate_file(data_filename)
    # calibrated_file = None
    output_files.append(calibrated_file)
    #  data_plot_files = plot_file(data_filename)
    #  calib_plot_files = plot_file(calibrated_file)

    # add other tasks below
    return output_files


def calibrate_file(data_filename: Path, output_level=2) -> Path:
    """
    Given an input file, calibrate it and return a new file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)
    output_level: int
        The requested data level of the output file.

    Returns
    -------
    output_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)

    Examples
    --------
    """

    # print(f"{cdfin.attrs['Data_type']}:{type(cdfin.attrs['Data_type'])}")
    if "l0" in str(data_filename):
        # calibrate L0 to L1
        cdfout = do_L0L1(data_filename)
        return cdfout

    elif "l1" in str(data_filename) or "l2" in str(data_filename):
        # calibrate L1 to L2
        raise NotImplementedError
    else:
        raise ValueError

    # below: template from SWxSOC. Eventually implement this.
    """
    log.info(f"Despiking removing {random.randint(0, 10)} spikes")
    log.warning(f"Despiking could not remove {random.randint(1, 5)}")

    calib_file = get_calibration_file(data_filename)
    if calib_file is None:
        raise ValueError(f"Calibration file for {data_filename} not found.")
    else:
        calib_data = read_calibration_file(calib_file)  # NOQA

    # example log messages

    return None
    """


def get_calibration_file(time: Time) -> Path:
    """
    Given a time, return the appropriate calibration file.

    Parameters
    ----------
    data_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)
    time: ~astropy.time.Time

    Returns
    -------
    calib_filename: str
        Fully specificied filename for the appropriate calibration file.

    Examples
    --------
    """
    return None


def read_calibration_file(calib_filename: Path):
    """
    Given a calibration, return the calibration structure.

    Parameters
    ----------
    calib_filename: str
        Fully specificied filename of the non-calibrated file (data level < 2)

    Returns
    -------
    output_filename: str
        Fully specificied filename of the appropriate calibration file.

    Examples
    --------
    """

    # if can't read the file

    return


def fillval(T):
    typename = str(T)
    if "float" in typename or "double" in typename:
        fillval = -1.0e31
    elif "uint" in typename:
        fillval = np.iinfo(T).max
    elif "int" in typename:
        fillval = np.iinfo(T).min
    else:
        raise ValueError

    return fillval


changes_L1 = [
    "Threshold1",
    "Threshold2",
    "SEU_mask",
    "Auto_accum",
    "Offset1",
    "Offset2",
    "Edge1",
    "Edge2",
    "Edge3",
    "Edge4",
    "Edge5",
    "Edge6",
    "Bus5V_voltage",
    "BusHV_current",
    "Detector_temperature",
    "CPU_temperature",
    "Bus3V3_voltage",
    "Calib_pulse1",
    "Calib_pulse2",
    "Accum_rate1",
    "Accum_rate2",
    # note: ebins counted as support
]


def do_L0L1(cdf_L0, save_dir=l1_depot):
    # bins to energy space
    # HK to SI units

    tlm_conf = telemetry_config()

    indata = SWXData.load(cdf_L0)
    meta_L1 = SWXData.global_attribute_template(
        instr_name="ifire", data_level="l1", version="1.0.0"
    )

    outdict = {k: indata.timeseries[k] for k in indata.timeseries.keys() if k != "time"}

    # TH1 and TH2
    outdict["Threshold1"] = outdict["Threshold1"].value / 2.5 * u.keV

    outdict["Threshold2"] = outdict["Threshold2"].value / 2.5 * u.keV

    # SEU_mask
    outdict["SEU_mask"] = (
        outdict["SEU_mask"].value * outdict["Accumulations"].value * 60 / 2.3 * u.s
    )

    # Auto_accum
    outdict["Auto_accum"] = outdict["Auto_accum"].value / 0.05 * u.Hz

    # Offset1
    outdict["Offset1"] = (outdict["Offset1"].value - 16000) / 32.7 * u.keV

    # Offset2
    outdict["Offset2"] = (outdict["Offset2"].value - 16000) / 32.7 * u.keV

    # Edge1-Edge6
    outdict["Edge1"] = outdict["Edge1"].value / 32.7 * u.keV
    outdict["Edge2"] = outdict["Edge2"].value / 32.7 * u.keV
    outdict["Edge3"] = outdict["Edge3"].value / 32.7 * u.keV
    outdict["Edge4"] = outdict["Edge4"].value / 32.7 * u.keV
    outdict["Edge5"] = outdict["Edge5"].value / 32.7 * u.keV
    outdict["Edge6"] = outdict["Edge6"].value / 32.7 * u.keV

    # Bus5V_voltage
    outdict["Bus5V_voltage"] = outdict["Bus5V_voltage"].value * 0.001615 * u.Volt
    # BusHV_current
    outdict["BusHV_current"] = outdict["BusHV_current"].value * 7.3 * 806e-12 * u.amp
    # Bus3V3_voltage
    outdict["Bus3V3_voltage"] = 5000 / (outdict["Bus3V3_voltage"].value + 1) * u.Volt

    # Detector_temperature
    outdict["Detector_temperature"] = (
        (outdict["Detector_temperature"].value * 806e-6 - 0.5) / 1e-2 * u.deg_C
    )
    # CPU_temperature
    outdict["CPU_temperature"] = (
        (outdict["CPU_temperature"].value * 806e-6 - 0.76) / 2.5e-3 + 31
    ) * u.deg_C

    # Calib_pulse1
    outdict["Calib_pulse1"] = (outdict["Calib_pulse1"].value - 16000) / 32.7 * u.keV
    # Calib_pulse2
    outdict["Calib_pulse2"] = (outdict["Calib_pulse2"].value - 16000) / 32.7 * u.keV

    # Accum_rate1
    outdict["Accum_rate1"] = (
        outdict["Accum_rate1"].value * 10 / outdict["Accumulations"].value * u.Hz
    )
    # Accum_rate2
    outdict["Accum_rate2"] = (
        outdict["Accum_rate2"].value * 10 / outdict["Accumulations"].value * u.Hz
    )
    # ebins

    ts = TimeSeries(time=indata.timeseries["time"], data=outdict)

    outdata = SWXData(
        timeseries=ts,
        support=indata.support,
        spectra=indata.spectra,
        meta=meta_L1,
    )

    # apply CATDESC to each TimeSeries column
    for k in outdata.timeseries.keys():
        tlm_names = list(filter(lambda x: tlm_conf[x].cdf_fieldnam == k, tlm_conf))
        if not tlm_names:
            # raise IndexError(f"No CATDESC def found for {k}!")
            print(f"No CATDESC def found for {k}!")
            continue
        tlm_obj = tlm_conf[tlm_names[0]]
        print(k)
        outdata.timeseries[k].meta["CATDESC"] = tlm_obj.description

    # manually set CATDESC for the following:
    # ts["time"].meta["CATDESC"] = "UTC time since some epoch"

    # convert all fillvals to their float equivalents.
    #   ACTUALLY, SWXData loads /everything/ into a float64, with NaNs where int FILLVALs were found. So just check for nan.

    # or maybe do nothing?
    # for k in changes_L1:
    #     outdata.timeseries[k][np.isnan(indata.timeseries[k])] = (
    #         fillval(np.float32) * outdata.timeseries[k].unit
    #     )

    # todo: do something to this as well...
    #
    # a curiosity of SWXData: I put this (ebins) into the source object as a support var, not as spectra. But when I read the file indata it shows up as a spectra variable. Not sure where it moved (in readback vs. write).

    # outdata.spectra["ebins"]

    cdf = outdata.save(save_dir)
    return cdf


if __name__ == "__main__":
    process_file(Path(sys.argv[1]))
