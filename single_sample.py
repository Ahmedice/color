# single_sample.py
import math
from typing import Dict, Any

import numpy as np

from utils import compute_concentration_from_a260, calculate_ratios, calculate_dilution, assess_purity


def process_single_sample(inputs: Dict[str, Any], protocol_settings: Dict[str, Any]) -> Dict[str, Any]:
    """Process single sample inputs and return results dict."""
    # Extract inputs
    a260 = inputs.get("a260")
    a280 = inputs.get("a280")
    ratio_260_230 = inputs.get("ratio_260_230")  # might be None
    factor = inputs.get("factor", 50.0)
    sample_type = inputs.get("sample_type", "DNA")
    target_conc = inputs.get("target_conc", protocol_settings.get("target_conc"))
    final_vol = inputs.get("final_vol", protocol_settings.get("final_vol"))

    # Input validation
    if final_vol <= 0:
        return {"Notes": "Final volume must be positive."}
    if target_conc < 0:
        return {"Notes": "Target concentration must be non-negative."}

    # Compute concentration
    concentration = compute_concentration_from_a260(a260, factor) if not math.isnan(a260) else np.nan

    # Ratios
    r260_280, r260_230 = calculate_ratios(inputs)
    ratio_note = ""

    # Dilution
    v1, v2, dilution_note = calculate_dilution(concentration, target_conc, final_vol)

    # Purity
    purity_verdict, purity_reco = assess_purity(sample_type, r260_280, r260_230)

    # Prepare notes
    notes = []
    if ratio_note:
        notes.append(ratio_note)
    if dilution_note:
        notes.append(dilution_note)
    if math.isnan(a260):
        notes.append("A260 is missing or invalid")
    if math.isnan(concentration):
        notes.append("Concentration invalid (NaN or missing A260)")

    return {
        "Concentration_ng_per_ul": round(concentration, 2) if not math.isnan(concentration) else np.nan,
        "260/280": round(r260_280, 2) if not np.isnan(r260_280) else np.nan,
        "260/230": round(r260_230, 2) if not np.isnan(r260_230) else np.nan,
        "Purity": purity_verdict,
        "Purity_notes": purity_reco,
        "Volume Sample (µl)": v1,
        "Volume Diluent (µl)": v2,
        "Protocol_target_ng_per_ul": target_conc,
        "Protocol_final_vol_ul": final_vol,
        "Notes": "; ".join(notes) if notes else ""
    }