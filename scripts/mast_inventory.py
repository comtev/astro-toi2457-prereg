"""Record what MAST holds for TIC 117799904, to document that no Sector 109 or 110 data were public.

Usage: uv run python scripts/mast_inventory.py OUT_JSON

Queries (all public, no login): MAST observations by target name, MAST observations within 0.01 deg
of the target, and the TESScut sector list. Writes every returned observation with its sector,
provenance, observation window and public release date, plus the query time and parameters.
Added on 2026-10-11 with the second revision of the pre-registration; not part of the original freeze.
"""

import json
import sys
from datetime import datetime, timezone

from astropy.time import Time
from astroquery.mast import Observations, Tesscut

TIC, RA, DEC = 117799904, 70.16725, -12.89081


def rows(table):
    out = []
    for r in table:
        def iso(v):
            try:
                return Time(float(v), format="mjd").utc.isot
            except Exception:
                return None
        out.append({"obs_collection": str(r["obs_collection"]), "provenance_name": str(r["provenance_name"]),
                    "target_name": str(r["target_name"]), "sequence_number": str(r["sequence_number"]),
                    "obs_id": str(r["obs_id"]), "t_min_utc": iso(r["t_min"]), "t_max_utc": iso(r["t_max"]),
                    "t_obs_release_utc": iso(r["t_obs_release"]), "dataproduct_type": str(r["dataproduct_type"])})
    return sorted(out, key=lambda x: (x["provenance_name"], x["sequence_number"], x["obs_id"]))


res = {"queried_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
       "target": {"tic": TIC, "ra_deg": RA, "dec_deg": DEC},
       "queries": {
           "by_target_name": f"Observations.query_criteria(target_name='{TIC}')",
           "by_position": f"Observations.query_region('{RA} {DEC}', radius='0.01 deg') filtered to obs_collection TESS or HLSP",
           "tesscut": f"Tesscut.get_sectors(objectname='TIC {TIC}')"}}
res["by_target_name"] = rows(Observations.query_criteria(target_name=str(TIC)))
reg = Observations.query_region(f"{RA} {DEC}", radius="0.01 deg")
reg = reg[[str(c) in ("TESS", "HLSP") for c in reg["obs_collection"]]]
res["by_position"] = rows(reg)
res["tesscut_sectors"] = [int(s) for s in Tesscut.get_sectors(objectname=f"TIC {TIC}")["sector"]]
seq = {x["sequence_number"] for x in res["by_target_name"] + res["by_position"]}
res["any_sector_109_or_110"] = bool({"109", "110"} & seq) or bool({109, 110} & set(res["tesscut_sectors"]))
with open(sys.argv[1], "w") as f:
    json.dump(res, f, indent=1)
print(json.dumps({k: res[k] for k in ("queried_utc", "tesscut_sectors", "any_sector_109_or_110")}, indent=1))
print(sorted(seq))
