"""StationXML construction for DAS streams."""
from __future__ import annotations

# --- Standard library --------------------------------------------------------
import logging
import math
import os
from typing import Any

# --- Third party -------------------------------------------------------------
import numpy as np
from obspy import Stream, UTCDateTime
from obspy.core.inventory import (
    Channel, CoefficientsTypeResponseStage, Equipment, InstrumentSensitivity, Inventory,
    Network, Operator, Person, PolesZerosResponseStage, Response, Site, Station,
)

# --- Local -------------------------------------------------------------------
from ._log import TRACE
from .coords import get_azimuth_dip
from .scaling import _INT32_STAGE_TEXT, INT32_SCALES

logger = logging.getLogger(__name__)


def station_xml_creation(st: Stream, global_params: dict[str, Any]) -> None:
    """
    Create a StationXML file from an ObsPy Stream `st` using configuration in `global_params`.

    Parameters:
    - st: ObsPy Stream containing station traces with coordinates and metadata.
    - global_params: dictionary with acquisition, path, exec, freq, trace header and stationxml params.

    The decimation stages are built from sp["decimation_steps"], i.e. from what
    was actually applied to the data, not from freq_params["decimation_factor"].

    Side effects:
    - Writes StationXML to `pp.get("MSEED_path")`.

    Raises:
    - KeyError if a trace has no latitude/longitude/elevation.
    - ValueError if the recorded decimation history does not reproduce the
      stream's sampling rate, or if int32_resolution is not a known label.
    """
    ap, pp, ep = global_params['acquisition_params'], global_params['path_params'], global_params['exec_params']
    th, sp = global_params['trace_header'], global_params['stationxml_params']

    for tr in st:
        for coord in ["latitude", "longitude", "elevation"]:
            if coord not in tr.stats.coordinates:
                logger.warning(f'Station {tr.stats.station} has no {coord} field!\nStationXML file not created.')
                raise KeyError(f"Missing '{coord}' in tr.stats.coordinates")

    # --- Decimation history -------------------------------------------------
    steps = sp.get("decimation_steps", [])
    total_decimation = math.prod(s["factor"] for s in steps)   # 1 if no steps

    # The recorded history must explain the data's sampling rate. If steps were
    # applied before this data was loaded (or state drifted), the response
    # would be wrong, so stop here instead of writing it.
    expected_rate = ap["original_sample_rate"] / total_decimation
    actual_rate = st[0].stats.sampling_rate
    if not np.isclose(expected_rate, actual_rate, rtol=1e-6):
        raise ValueError(
            f"Decimation history is inconsistent with the data: "
            f"original_sample_rate / total decimation = {expected_rate} Hz, "
            f"but the stream is at {actual_rate} Hz. Steps applied before this "
            f"data was loaded are not recorded, so the StationXML response "
            f"would be wrong.")

    # --- Response stages ----------------------------------------------------
    response_stages = [
        PolesZerosResponseStage(
            stage_sequence_number=1,
            name="Gain",
            stage_gain=1,
            stage_gain_frequency=1.0,
            input_units="counts",
            output_units="counts",
            input_units_description=ap["sensitivityUnits"],
            output_units_description=ap["sensitivityUnits"],
            normalization_factor=1,
            normalization_frequency=0,
            pz_transfer_function_type="LAPLACE (RADIANS/SECOND)",
            zeros=[],
            poles=[],
        )
    ]

    for i, step in enumerate(steps):
        fir = step.get("fir")
        if fir is not None:
            # Linear-phase FIR applied by resample_poly (decimate_dfdas).
            # The filter delay, (ntaps-1)/2 samples, is removed by resample_poly
            # (output is centred), so delay == correction.
            fir = np.asarray(fir, dtype=float)
            delay = ((len(fir) - 1) / 2) / step["input_sampling_rate"]
            stage_name = (f"Linear-phase FIR (Kaiser) low-pass + Decimation step {i + 1}. "
                          f"Decimation factor: {step['factor']}")
            numerator, denominator = fir.tolist(), []
        else:
            # Zero-phase Chebyshev II applied by decimate_st (filtfilt)
            delay = 0
            stage_name = (f"Low-pass Chebyshev filter stage + Decimation step {i + 1}. "
                          f"Decimation factor: {step['factor']}")
            numerator = np.asarray(step["b"]).tolist()
            denominator = np.asarray(step["a"]).tolist()

        response_stages.append(CoefficientsTypeResponseStage(
            stage_sequence_number=i + 2,
            name=stage_name,
            cf_transfer_function_type="DIGITAL",
            stage_gain_frequency=1.0,
            stage_gain=1,
            input_units="counts",
            output_units="counts",
            input_units_description=ap["sensitivityUnits"],
            output_units_description=ap["sensitivityUnits"],
            decimation_factor=step["factor"],
            decimation_input_sample_rate=step["input_sampling_rate"],
            decimation_offset=0,
            decimation_delay=delay,
            decimation_correction=delay,
            denominator=denominator,
            numerator=numerator,
        ))

    last_stage = len(response_stages) + 1
    output_rate = ap["original_sample_rate"] / total_decimation

    response_stages.append(CoefficientsTypeResponseStage(
        stage_sequence_number=last_stage,
        name="Sensitivity",
        cf_transfer_function_type="DIGITAL",
        numerator=[],
        denominator=[],
        decimation_factor=1,
        decimation_correction=0,
        decimation_delay=0,
        decimation_offset=0,
        decimation_input_sample_rate=output_rate,
        stage_gain=ap["sensitivity"],
        stage_gain_frequency=1.0,
        input_units="counts",
        output_units="counts",
        input_units_description=ap["sensitivityUnits"],
        output_units_description=ap["sensitivityUnits"],
    ))

    # Counts -> physical units stage, only for int32 output. stage_gain comes
    # from the same table dfdas_to_mseed multiplies by, so they cannot differ.
    resolution = sp.get("int32_resolution")
    logger.debug(f"int32_resolution: {resolution}")
    if resolution is not None:
        if resolution not in INT32_SCALES or resolution not in _INT32_STAGE_TEXT:
            raise ValueError(
                f"Unknown int32_resolution {resolution!r}; "
                f"expected one of {list(INT32_SCALES)}.")
        stage_name, out_units, out_desc = _INT32_STAGE_TEXT[resolution]
        response_stages.append(CoefficientsTypeResponseStage(
            stage_sequence_number=last_stage + 1,
            name=stage_name,
            cf_transfer_function_type="DIGITAL",
            stage_gain=INT32_SCALES[resolution],
            stage_gain_frequency=1.0,
            input_units="counts",
            output_units=out_units,
            input_units_description="Digital counts mapped to int32",
            output_units_description=out_desc,
            numerator=[],
            denominator=[],
            decimation_factor=1,
            decimation_offset=0,
            decimation_delay=0,
            decimation_correction=0,
            decimation_input_sample_rate=output_rate,
        ))

    response = Response(
        resource_id="OptoDAS header",
        instrument_sensitivity=InstrumentSensitivity(
            value=ap["sensitivity"],
            frequency=1.0,
            input_units="counts",
            output_units="counts",
            output_units_description="strain per meter",
            input_units_description=ap["sensitivityUnits"] + " phase difference expressed in radians per strain per meter",
            frequency_range_start=0.001,
            frequency_range_end=2000.0,
            frequency_range_db_variation=0,
        ),
        response_stages=response_stages)

    # --- Equipment / operator / network ------------------------------------
    install_date = UTCDateTime(*th.get("net_starttime"))

    def _das_equipment():
        return Equipment(
            type="DAS",
            description="OptoDAS",
            manufacturer="ASN",
            serial_number="44",
            installation_date=install_date,
            resource_id="OptoDAS_44",
        )

    sensor, data_logger, pre_amplifier = _das_equipment(), _das_equipment(), _das_equipment()

    operator = [
        Operator(
            agency="ARDITI",
            contacts=[
                Person(
                    names=["Afonso Loureiro"],
                    agencies=["ARDITI, Funchal, Portugal",
                              "Instituto Dom Luiz, Lisboa, Portugal"],
                    emails=["maloureiro@fc.ul.pt"])],
            website="www.arditi.pt",
        )
    ]

    inv = Inventory(networks=[], source=th.get("inventory_source"))

    net = Network(
        code=st[0].stats.network,
        stations=[],
        description=th.get("net_description"),
        start_date=UTCDateTime(ap.get("start_date")),
        end_date=UTCDateTime(ap.get("end_date")),
        total_number_of_stations=th.get("net_no_stations"),
        selected_number_of_stations=len(st),
    )

    for tr in st:
        station_code = tr.stats.station
        if ep.get("encode_station_name"):
            das_channel = tr.stats['decoded_station']
        else:
            das_channel = tr.stats.station

        sta = Station(
            code=station_code,
            latitude=tr.stats.coordinates["latitude"],
            longitude=tr.stats.coordinates["longitude"],
            elevation=tr.stats.coordinates["elevation"],
            creation_date=install_date,
            start_date=UTCDateTime(ap.get("start_date")),
            end_date=UTCDateTime(ap.get("end_date")),
            site=Site(name=f'{th.get("sta_name_prefix")} {tr.stats.channel} channel no.{str(das_channel).zfill(5)}'),
            total_number_of_channels=1,
        )
        sta.operators = operator

        azimuth, dip = get_azimuth_dip(int(tr.stats.distance / ap.get("dx") + 1), global_params, roi_channel=False)
        cha = Channel(
            code="HSF",
            location_code=th.get("location"),
            latitude=tr.stats.coordinates["latitude"],
            longitude=tr.stats.coordinates["longitude"],
            elevation=tr.stats.coordinates["elevation"],
            start_date=UTCDateTime(ap.get("start_date")),
            end_date=UTCDateTime(ap.get("end_date")),
            sample_rate=tr.stats.sampling_rate,
            depth=0,
            azimuth=azimuth,
            dip=dip,
        )
        cha.response = response
        cha.equipments = [sensor, data_logger, pre_amplifier]
        cha.sensor = sensor
        cha.data_logger = data_logger
        cha.pre_amplifier = pre_amplifier
        sta.channels.append(cha)
        net.stations.append(sta)

    inv.networks.append(net)

    file = os.path.join(pp.get("MSEED_path"), f'RESP.{th.get("network")}.HSF.xml')
    os.makedirs(pp.get("MSEED_path"), exist_ok=True)
    inv.write(file, format="stationxml", validate=True)

    logger.debug(f'StationXML file written to: {file}')

    if ep.get("verbose", False) == 2:
        logger.log(TRACE, "%s", inv)
        inv[0][-1][0].response.plot(min_freq=0.01, unwrap_phase=True)

    return
