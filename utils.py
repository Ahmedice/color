# utils.py
import json
import math
from typing import Dict, Tuple, Any

import pandas as pd
import numpy as np


def load_config(path: str = "config.json") -> Dict[str, Dict[str, Any]]:
    """Load JSON config file with default protocols."""
    with open(path, "r", encoding="utf-8") as f:
        config = json.load(f)
    return config


def map_columns(columns: pd.Index) -> Dict[str, str]:
    """Map various CSV header names to canonical names.
    Input: pandas Index of columns. Output: dict mapping canonical -> actual column name or None.
    """
    # Known name variants
    mapping_variants = {
        "sample_id": ["sample id", "id", "sample", "sample_name", "sample name"],
        "a260": ["a260", "a260 (abs)", "a260nm", "abs260", "a260 (abs)"],
        "a280": ["a280", "a280 (abs)", "a280nm", "abs280", "a280 (abs)"],
        "a230": ["a230", "a230 (abs)", "a230nm", "abs230", "a230 (abs)"],
        "nucleic_acid": ["nucleic acid", "nucleicacid", "conc", "concentration", "concentration (ng/µl)", "concentration (ng/ul)", "nucleic acid (ng/µl)"],
        "sample_type": ["sample type", "type", "sample_type", "nucleic_type", "nucleic acid type"],
        "factor": ["factor", "conversion factor", "factor50", "cfactor"]
    }
    # canonical -> detected column or None
    detected = {k: None for k in mapping_variants.keys()}
    lower_map = {col.lower(): col for col in columns}

    for canon, variants in mapping_variants.items():
        for var in variants:
            # try exact match
            if var in lower_map:
                detected[canon] = lower_map[var]
                break
        if detected[canon] is None:
            # try contains match
            for col_lower, original in lower_map.items():
                for var in variants:
                    if var in col_lower:
                        detected[canon] = original
                        break
                if detected[canon]:
                    break
    return detected


def safe_float(x) -> Tuple[float, str]:
    """Try to convert to float. Return (value, note)."""
    if pd.isna(x):
        return (np.nan, "missing")
    try:
        val = float(x)
        return (val, "")
    except Exception:
        return (np.nan, "non-numeric")


def compute_concentration_from_a260(a260: float, factor: float) -> float:
    """Compute concentration from A260 and factor."""
    return a260 * factor


def calculate_ratios(a260: float, a280: float, ratio_260_230: float = None) -> Tuple[float, float, str]:
    """Calculate 260/280 and 260/230 ratios. Returns (r260_280, r260_230, note)."""
    note = ""
    # 260/280
    try:
        if a280 == 0 or math.isnan(a280):
            r260_280 = np.nan
            note += "A280 missing or zero; "
        else:
            r260_280 = a260 / a280
    except Exception:
        r260_280 = np.nan
        note += "Error computing 260/280; "

    # 260/230
    if ratio_260_230 is None or math.isnan(ratio_260_230):
        r260_230 = a260 / 0.5 if not math.isnan(a260) else np.nan
        note += "نسبة 260/230 مفقودة — استخدم افتراضي 0.5 لحسابها; "
    else:
        r260_230 = ratio_260_230
    return (r260_280, r260_230, note.strip())


def calculate_dilution(concentration: float, target_conc: float, final_vol: float) -> Tuple[float, float, str]:
    """Calculate V1 and V2 for dilution. Return rounded V1 and V2 and note if issues."""
    note = ""
    if concentration == 0 or math.isnan(concentration):
        return (np.nan, np.nan, "zero concentration or invalid concentration")
    # V1 = (C2 * V2) / C1  where V2 is final_vol
    v1 = (target_conc * final_vol) / concentration
    if v1 > final_vol:
        note = "sample too dilute to reach target — need to concentrate or use more sample"
    v2 = final_vol - v1
    v1_rounded = round(v1, 2)
    v2_rounded = round(v2, 2)
    return (v1_rounded, v2_rounded, note)


def assess_purity(r260_280: float, r260_230: float, sample_type: str) -> Tuple[str, str]:
    """Assess purity based on ratios and sample type.
    Returns (verdict, recommendation) both in Arabic short texts.
    """
    sample_type_up = (sample_type or "").strip().lower()
    verdict = "حسن"  # OK by default
    recommendations = []

    # DNA rules
    if sample_type_up == "dna":
        if np.isnan(r260_280):
            recommendations.append("غير ممكن حساب نسبة 260/280 بسبب نقص البيانات.")
        else:
            if r260_280 < 1.7:
                recommendations.append("نسبة 260/280 منخفضة (<1.7): قد تشير إلى وجود بقايا بروتين أو ملوثات أخرى قد تؤثر على دقة القياس.")
        if np.isnan(r260_230):
            recommendations.append("غير ممكن حساب نسبة 260/230 بسبب نقص البيانات.")
        else:
            if r260_230 < 1.8:
                recommendations.append("نسبة 260/230 منخفضة (<1.8): قد تشير إلى تلوث بالملح أو مواد حاملة تؤثر على نقاوة العينة.")

    # RNA rules
    elif sample_type_up == "rna":
        if np.isnan(r260_280):
            recommendations.append("غير ممكن حساب نسبة 260/280 بسبب نقص البيانات.")
        else:
            if r260_280 < 1.9:
                recommendations.append("نسبة 260/280 منخفضة (<1.9): قد تشير إلى تحلل الحمض النووي أو تلوث في العينة.")
        if np.isnan(r260_230):
            recommendations.append("غير ممكن حساب نسبة 260/230 بسبب نقص البيانات.")
        else:
            if r260_230 < 1.8:
                recommendations.append("نسبة 260/230 منخفضة (<1.8): قد تدل على وجود ملوثات مثل الأملاح أو المواد الحاملة.")

    # Protein rules
    elif sample_type_up in ["protein", "prot"]:
        if not np.isnan(r260_280) and r260_280 > 1.2:
            recommendations.append("نسبة 260/280 مرتفعة (>1.2): مما قد يدل على وجود تلوث بالحمض النووي داخل العينة البروتينية.")

    # Unknown or other types
    else:
        if not np.isnan(r260_280) and r260_280 < 1.7:
            recommendations.append("نسبة 260/280 منخفضة: قد تشير إلى وجود ملوثات أو شوائب في العينة.")
        if not np.isnan(r260_230) and r260_230 < 1.8:
            recommendations.append("نسبة 260/230 منخفضة: قد تعكس وجود ملوثات مثل الأملاح أو المواد الحاملة.")

    # تحديد النتيجة النهائية
    if recommendations:
        verdict = "تحذير"
    else:
        verdict = "حسن"
        recommendations = ["نقاوة العينة ضمن الحدود المقبولة والمعايير العلمية المتعارف عليها"]

    recommendation_text = "؛ ".join(recommendations)
    return (verdict, recommendation_text)


def safe_ratio(x) -> Tuple[float, str]:
   """Try to convert to float. Return (value, note)."""
   if pd.isna(x):
       return (np.nan, "missing")
   try:
       val = float(x)
       return (val, "")
   except Exception:
       return (np.nan, "non-numeric")