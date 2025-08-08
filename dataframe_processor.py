# dataframe_processor.py
import math
from typing import Dict, Any

import pandas as pd
import numpy as np

from utils import safe_float, compute_concentration_from_a260, calculate_ratios, calculate_dilution, assess_purity, safe_ratio


def process_dataframe(df: pd.DataFrame, detected_map: Dict[str, str], protocol_settings: Dict[str, Any], default_factor: float = 50.0) -> pd.DataFrame:
    """Process uploaded dataframe row-wise and return augmented dataframe."""
    results = []
    for idx, row in df.iterrows():
        # read values using detected_map, robustly
        raw_a260 = row[detected_map["a260"]] if detected_map.get("a260") in row.index else np.nan
        raw_a280 = row[detected_map["a280"]] if detected_map.get("a280") in row.index else np.nan
        raw_ratio_260_230 = row[detected_map["ratio_260_230"]] if detected_map.get("ratio_260_230") in row.index else np.nan
        raw_nucleic = row[detected_map["nucleic_acid"]] if detected_map.get("nucleic_acid") in row.index else np.nan
        raw_factor = row[detected_map["factor"]] if detected_map.get("factor") in row.index else np.nan
        raw_sample_type = row[detected_map["sample_type"]] if detected_map.get("sample_type") in row.index else None
        raw_sample_id = row[detected_map["sample_id"]] if detected_map.get("sample_id") in row.index else f"row_{idx+1}"

        # safe conversions
        a260, note_a260 = safe_float(raw_a260)
        a280, note_a280 = safe_float(raw_a280)
        ratio_260_230, note_ratio_260_230 = safe_ratio(raw_ratio_260_230) if not pd.isna(raw_ratio_260_230) else (np.nan, "missing")
        nucleic_val, note_nuc = safe_float(raw_nucleic)
        factor_val, note_factor = safe_float(raw_factor)

        # Determine factor to use
        factor_to_use = factor_val if not math.isnan(factor_val) else default_factor

        # Determine concentration: use device reported if present and numeric else compute
        if not math.isnan(nucleic_val):
            conc = nucleic_val
            conc_note = "Concentration from device"
        else:
            if math.isnan(a260):
                conc = np.nan
                conc_note = "No A260 and no device concentration"
            else:
                conc = compute_concentration_from_a260(a260, factor_to_use)
                conc_note = f"Computed from A260 with factor {factor_to_use}"

        # Ratios
        r260_280, r260_230, ratio_note = calculate_ratios(a260, a280, ratio_260_230)

        # Dilution
        v1, v2, dilution_note = calculate_dilution(conc, protocol_settings.get("target_conc"), protocol_settings.get("final_vol"))

        # Purity
        purity_verdict, purity_reco = assess_purity(r260_280, r260_230, (raw_sample_type or "").strip())

        # Notes aggregation
        notes_list = []
        if note_a260 == "missing":
            notes_list.append("A260: القيمة مفقودة")
        elif note_a260 == "non-numeric":
            notes_list.append("A260: القيمة غير رقمية")
        if note_a280 == "missing":
            notes_list.append("A280: القيمة مفقودة")
        elif note_a280 == "non-numeric":
            notes_list.append("A280: القيمة غير رقمية")
        if note_ratio_260_230 == "missing":
            notes_list.append("نسبة 260/230: القيمة مفقودة")
        if note_nuc == "non-numeric":
            notes_list.append("تركيز الجهاز: القيمة غير رقمية")
        if math.isnan(conc):
            notes_list.append("التركيز غير صالح (NaN)")
        if dilution_note:
            notes_list.append(dilution_note)
        if ratio_note:
            # ratio_note is English fallback; translate minimal keywords
            if "A230" in ratio_note:
                notes_list.append("A230 مفقود/صفر أو غير صالح")
            if "A280" in ratio_note:
                notes_list.append("A280 مفقود/صفر أو غير صالح")

        # Build result row
        result_row = dict(row)  # keep original columns
        result_row.update({
            "Conc_ng_per_ul": round(conc, 2) if not math.isnan(conc) else np.nan,
            "260/280": round(r260_280, 2) if not math.isnan(r260_280) else np.nan,
            "260/230": round(r260_230, 2) if not math.isnan(r260_230) else np.nan,
            "Purity": purity_verdict,
            "Purity_notes": purity_reco,
            "Volume Sample (µl)": v1,
            "Volume Diluent (µl)": v2,
            "Notes": "; ".join(notes_list) if notes_list else ""
        })
        results.append(result_row)
    result_df = pd.DataFrame(results)
    # Reorder columns: original columns first, then added ones
    added_cols = ["Conc_ng_per_ul", "260/280", "260/230", "Purity", "Purity_notes", "µl_DNA", "µl_H2O", "Notes"]
    existing = [c for c in result_df.columns if c not in added_cols]
    result_df = result_df[existing + [c for c in added_cols if c in result_df.columns]]
    return result_df