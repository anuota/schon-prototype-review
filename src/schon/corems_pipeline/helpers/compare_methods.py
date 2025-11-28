#!/usr/bin/env python3
# SPDX-License-Identifier: EUPL-1.2
from pathlib import Path
from xml.dom import minidom


def parse_apex_method(path: Path) -> dict:
    """Parse apexAcquisition.method into a dict[name] = value."""
    if not path.is_file():
        raise FileNotFoundError(path)

    xmldoc = minidom.parse(path.open())
    root = xmldoc.documentElement
    params = {}

    # acquisition_time from <methodmetadata><date>
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
                            # usually in attribute "value"
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


def compare_param_dicts(a_params: dict, b_params: dict, label_a: str, label_b: str):
    keys_a = set(a_params.keys())
    keys_b = set(b_params.keys())

    only_a = sorted(keys_a - keys_b)
    only_b = sorted(keys_b - keys_a)
    common = sorted(keys_a & keys_b)

    print(f"Total params in {label_a}: {len(keys_a)}")
    print(f"Total params in {label_b}: {len(keys_b)}")
    print(f"Common params:           {len(common)}")
    print(f"Only in {label_a}:       {len(only_a)}")
    print(f"Only in {label_b}:       {len(only_b)}")
    print()

    if only_a:
        print(f"=== Parameters only in {label_a} ===")
        for k in only_a:
            print(f"  {k} = {a_params[k]}")
        print()

    if only_b:
        print(f"=== Parameters only in {label_b} ===")
        for k in only_b:
            print(f"  {k} = {b_params[k]}")
        print()

    diff = []
    for k in common:
        va = a_params[k]
        vb = b_params[k]
        if va != vb:
            diff.append((k, va, vb))

    if diff:
        print(f"=== Parameters with different values ({label_a} vs {label_b}) ===")
        for k, va, vb in diff:
            print(f"  {k}: {label_a}='{va}'  |  {label_b}='{vb}'")


if __name__ == "__main__":
    # Paths you gave:
    path_a = Path(
        "/Users/anya/Coding/CoreMS/tests/tests_data/ftms/ESI_NEG_SRFA.d/"   #test data works
        "1901_Neg_CompMix_1.4s_150-1000.m/apexAcquisition.method"
    )
    # path_a = Path(
    #     "/Users/anya/Coding/ML_FT-ICR-MS/Measurement_d-files/calibrated_d-files/"
    #     "ESI_neg_G017736-0_100_200sc_000001.d/"                                             #works
    #     "Broad_150-1000_crudeoil_20180118_IA100ms.m/apexAcquisition.method"
    # )
    path_b = Path(
        "/Users/anya/Coding/ML_FT-ICR-MS/Measurement_d-files/calibrated_d-files/"
        "ESI_neg_G015548-0_100_200sc_000001.d/"                                    #doesnt work
        "ESI_neg_150_750_4M_200sc_cid.m/apexAcquisition.method"
    )

    # path_b = Path(
    #     "/Users/anya/Coding/ML_FT-ICR-MS/Measurement_d-files/calibrated_d-files/"
    #     "ESI_neg_G017736-0_100_200sc_000001.d/"                                             #works
    #     "Broad_150-1000_crudeoil_20180118_IA100ms.m/apexAcquisition.method"
    # )

    params_a = parse_apex_method(path_a)
    params_b = parse_apex_method(path_b)

    compare_param_dicts(params_a, params_b, "ESI_NEG_SRFA", "ESI_neg_G015548")
    #compare_param_dicts(params_a, params_b, "ESI_NEG_G017736", "ESI_neg_G015548")

    