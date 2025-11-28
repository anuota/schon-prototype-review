# SPDX-License-Identifier: EUPL-1.2
import os
from pathlib import Path
from xml.dom import minidom
from datetime import datetime

BASE_DIR = Path("/Users/anya/Coding/ML_FT-ICR-MS/Measurement_d-files")
TEST_DIR = Path("/Users/anya/Coding/CoreMS/tests/tests_data/ftms")


def find_method_file(d_folder: Path):
    """Locate apexAcquisition.method inside a .d folder."""
    for root, dirs, files in os.walk(d_folder):
        for f in files:
            if f == "apexAcquisition.method":
                return Path(root) / f
        # Ensure we descend into .m subdirectories (Bruker sometimes nests method files there)
        for sub in dirs:
            if sub.endswith(".m"):
                sub_path = Path(root) / sub
                for r2, d2, f2s in os.walk(sub_path):
                    for f2 in f2s:
                        if f2 == "apexAcquisition.method":
                            return Path(r2) / f2
    return None



# Helper: parse apexAcquisition.method into a dict[name] = value, as in the comparison script
def parse_apex_method(path: Path) -> dict:
    """Parse apexAcquisition.method into a dict[name] = value.

    This matches the logic used in the comparison script, so downstream
    interpretation of fields like acquisition_time and EXC_Freq_* is consistent.
    """
    if not path.is_file():
        raise FileNotFoundError(path)

    xmldoc = minidom.parse(path.open())
    root = xmldoc.documentElement
    params: dict[str, str] = {}

    # acquisition_time from <methodmetadata><primarykey><date>
    for child in root.childNodes:
        if child.nodeName == "methodmetadata":
            for section in child.childNodes:
                if section.nodeName == "primarykey":
                    for el in section.childNodes:
                        if el.nodeName == "date" and el.childNodes:
                            params["acquisition_time"] = el.childNodes[0].nodeValue

    # Polarity from <reportinfo> section
    for child in root.childNodes:
        if child.nodeName == "reportinfo":
            for section in child.childNodes:
                if section.nodeName == "section" and section.getAttribute("title") == "Main":
                    for el in section.childNodes:
                        if el.nodeName == "section" and el.getAttribute("title") == "Polarity":
                            try:
                                val = el.childNodes[1].getAttribute("value")
                            except Exception:
                                val = None
                            if val:
                                params["Polarity"] = val

    # Generic <paramlist><param name="..."><value>...</value></param>
    for child in root.childNodes:
        if child.nodeName == "paramlist":
            for param in child.childNodes:
                if param.nodeName != "param":
                    continue
                name = param.getAttribute("name")
                value = None
                for el in param.childNodes:
                    if el.nodeName == "value" and el.firstChild is not None:
                        value = el.firstChild.nodeValue
                params[name] = value

    return params


# def extract_acquisition_time(method_file: Path):
#     """Parse apexAcquisition.method and extract acquisition_time + EXC_Freq_*.

#     Uses the same parsing strategy as in the comparison script (parse_apex_method)
#     so that all fields are interpreted consistently.
#     """
#     try:
#         params = parse_apex_method(method_file)
#     except Exception:
#         return None, None, None

#     acq_time = None
#     raw_acq = params.get("acquisition_time")

#     if raw_acq:
#         # Same format as used in the comparison code:
#         # e.g. "Feb_5_2019 13:44:41.642"
#         try:
#             acq_time = datetime.strptime(raw_acq, "%b_%d_%Y %H:%M:%S.%f")
#         except Exception:
#             # If parsing fails, just leave it as None so sorting still works
#             acq_time = None

#     exc_high = params.get("EXC_Freq_High")
#     exc_low = params.get("EXC_Freq_Low")

#     return acq_time, exc_high, exc_low
def extract_acquisition_time(method_file: Path):
    """Parse apexAcquisition.method and extract acquisition_time and selected parameters."""
    try:
        params = parse_apex_method(method_file)
    except Exception:
        return None

    # Time
    raw_acq = params.get("acquisition_time")
    if raw_acq:
        try:
            acq_time = datetime.strptime(raw_acq, "%b_%d_%Y %H:%M:%S.%f")
        except Exception:
            acq_time = None
    else:
        acq_time = None

    # Parameters to extract (easy to extend later)
    EXTRA_KEYS = [
        "EXC_Freq_High",
        "EXC_Freq_Low",
        # "SW",
        "SW_h",
        "SW_h_Broadband",
        # "O1",
        "FR_low",
        "P_30",
    ]

    extracted = {key: params.get(key) for key in EXTRA_KEYS}

    return {
        "folder": str(method_file.parent),
        "time": acq_time,
        **extracted,
    }

# def scan_all_d_folders(base_dir: Path):
#     results = []

#     for root, dirs, files in os.walk(base_dir):
#         for d in dirs:
#             if d.endswith(".d"):
#                 d_folder = Path(root) / d
#                 method_file = find_method_file(d_folder)
#                 if not method_file:
#                     results.append((str(d_folder), None, None, None))
#                     continue

#                 acq_time, exc_high, exc_low = extract_acquisition_time(method_file)
#                 results.append((str(d_folder), acq_time, exc_high, exc_low))

#     return results
def scan_all_d_folders(base_dir: Path):
    """Scan all .d folders under base_dir and collect parsed method info.

    Returns a list of dicts with at least the keys:
    - folder
    - time
    - EXC_Freq_High
    - EXC_Freq_Low
    - SW
    - SW_h
    - SW_h_Broadband
    - O1
    - FR_low
    - P_30
    """
    results = []

    # default skeleton for entries where we cannot parse anything
    def empty_entry(d_folder: Path) -> dict:
        return {
            "folder": str(d_folder),
            "time": None,
            "EXC_Freq_High": None,     # high freq of excitation sweep (Hz)
            "EXC_Freq_Low": None,      # low freq of excitation sweep (Hz)
            # "SW": None,               # width of the frequency band that is turned into your spectrum. frequency window ~ [FR_low, FR_low + SW]
            "SW_h": None,              # spectral width of main FT-ICR acquisition (Hz)
            "SW_h_Broadband": None,    # spectral width for broadband acquisition (Hz)
            # "O1": None,               # reference / center-ish frequency inside that window
            "FR_low": None,            # low frequency cutoff / boundary (Hz)
            "P_30": None,              # calibration low-frequency point ~ FR_low
        }

    for root, dirs, files in os.walk(base_dir):
        for d in dirs:
            if not d.endswith(".d"):
                continue

            d_folder = Path(root) / d
            method_file = find_method_file(d_folder)

            if not method_file:
                # No method file found – keep an empty entry so we still see the folder
                results.append(empty_entry(d_folder))
                continue

            info = extract_acquisition_time(method_file)
            if not info:
                info = empty_entry(d_folder)

            # Ensure all expected keys are present even if extract_acquisition_time
            # returns only a subset
            full_info = empty_entry(d_folder)
            full_info.update(info)

            results.append(full_info)

    return results

if __name__ == "__main__":
    all_results = scan_all_d_folders(BASE_DIR)
    all_results += scan_all_d_folders(TEST_DIR)

    # print("\n=== Acquisition Times in .d folders ===\n")
    # for folder, acq_time, exc_high, exc_low in sorted(
    #         all_results, key=lambda x: (x[1] if x[1] else datetime.min)
    #     ):
    #     if acq_time:
    #         print(f"{acq_time}  EXC_Freq_High={exc_high}  EXC_Freq_Low={exc_low}   —   {folder}")
    #     else:
    #         print(f"(no acquisition_time)  EXC_Freq_High={exc_high}  EXC_Freq_Low={exc_low}   —   {folder}")
    print("\n=== Acquisition Times in .d folders ===\n")

    # unified sorted list
    sorted_results = sorted(
        all_results,
        key=lambda x: (x["time"] if x["time"] else datetime.min)
    )

    for r in sorted_results:
        t = r["time"]
        if t:
            t_str = t.strftime("%Y-%m-%d %H:%M:%S.%f")
        else:
            t_str = "(no acquisition_time)"

        print(
            f"{t_str}  "
            f"EXC_H={r.get('EXC_Freq_High')}  "
            # f"SW={r.get('SW')}  "
            f"SW_h={r.get('SW_h')}  "
            f"SW_h_BB={r.get('SW_h_Broadband')}  "
            # f"O1={r.get('O1')}  "
            f"EXC_L={r.get('EXC_Freq_Low')}  "
            f"FR_low={r.get('FR_low')}  "
            f"P_30={r.get('P_30')}  "
            f"—  {r['folder']}"
        )