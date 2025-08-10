import json
import numpy as np
import pandas as pd
from typing import Dict, Tuple, Any

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
    missing = [k for k, v in detected.items() if v is None]
    return detected, missing

def safe_float(val):
    try:
        return float(val), ""
    except:
        return np.nan, "قيمة غير صالحة"

def compute_concentration_from_a260(a260: float, factor: float) -> float:
    """Compute concentration from A260 and factor."""
    return a260 * factor

def calculate_ratios(row):
    # حساب 260/280
    a260 = row.get("a260")
    a280 = row.get("a280")
    a230 = row.get("a230")

    ratio_260_280 = np.nan
    if pd.notna(a260) and pd.notna(a280) and a280 != 0:
        ratio_260_280 = a260 / a280

    # حساب 260/230
    ratio_260_230 = np.nan
    if pd.notna(row.get("ratio_260_230")):
        ratio_260_230 = row.get("ratio_260_230")
    elif pd.notna(a260) and pd.notna(a230) and a230 != 0:
        ratio_260_230 = a260 / a230

    return ratio_260_280, ratio_260_230

def compute_concentration(row, factor):
    # إذا كانت قيمة التركيز موجودة ومستعملة لا تغيرها
    if "Conc_ng_per_ul" in row and pd.notna(row["Conc_ng_per_ul"]):
        return row["Conc_ng_per_ul"]
    # غير ذلك حسب القانون: A260 * factor
    a260 = row.get("a260")
    if pd.notna(a260):
        return a260 * factor
    return np.nan

def assess_purity(sample_type, ratio_260_280, ratio_260_230):
    """
    تقييم نقاوة العينة بناء على نوع العينة و النسب.
    القواعد مبنية على قواعد عامة:
    DNA: 260/280 ~ 1.8, 260/230 ~ 2.0-2.2
    RNA: 260/280 ~ 2.0, 260/230 ~ 2.0-2.2
    Protein: 260/280 ~ 0.6-0.7 (عادي)
    """
    notes = []
    verdict = "حسن"

    if sample_type == "DNA":
        if np.isnan(ratio_260_280) or ratio_260_280 < 1.7:
            verdict = "تحذير"
            notes.append("نسبة 260/280 منخفضة (<1.7): قد تشير إلى وجود بقايا بروتين أو ملوثات أخرى.")
        elif ratio_260_280 > 2.0:
            verdict = "تحذير"
            notes.append("نسبة 260/280 مرتفعة (>2.0): قد تشير إلى تلوث بمواد أخرى.")

        if np.isnan(ratio_260_230) or ratio_260_230 < 1.8:
            verdict = "تحذير"
            notes.append("نسبة 260/230 منخفضة (<1.8): قد تشير إلى ملوثات عضوية أو أوساخ.")

    elif sample_type == "RNA":
        if np.isnan(ratio_260_280) or ratio_260_280 < 1.8:
            verdict = "تحذير"
            notes.append("نسبة 260/280 منخفضة (<1.8): قد تشير إلى تحلل الحمض النووي أو تلوث.")
        elif ratio_260_280 > 2.2:
            verdict = "تحذير"
            notes.append("نسبة 260/280 مرتفعة (>2.2): قد تشير إلى تلوث.")
        
        if np.isnan(ratio_260_230) or ratio_260_230 < 1.8:
            verdict = "تحذير"
            notes.append("نسبة 260/230 منخفضة (<1.8): قد تشير إلى ملوثات.")

    elif sample_type == "Protein":
        if np.isnan(ratio_260_280) or ratio_260_280 < 0.5:
            verdict = "تحذير"
            notes.append("نسبة 260/280 منخفضة (<0.5): قد تدل على نقاوة منخفضة للبروتين.")
        elif ratio_260_280 > 0.8:
            verdict = "تحذير"
            notes.append("نسبة 260/280 مرتفعة (>0.8): قد تدل على تلوث بالحمض النووي.")
        
        # 260/230 ليس مقياسًا دقيقًا للبروتين عادة، يمكن تخطيه

    # تحقق إذا كانت النسب مفقودة وأضف ملاحظات
    if np.isnan(ratio_260_280):
        notes.append("نسبة 260/280 غير متوفرة أو غير صالحة.")
    if np.isnan(ratio_260_230):
        notes.append("نسبة 260/230 غير متوفرة أو غير صالحة.")

    purity_notes = "؛ ".join(notes)
    return verdict, purity_notes

def calculate_dilution(conc_sample, target_conc, final_vol):
    note = ""
    if conc_sample is None or conc_sample <= 0 or np.isnan(conc_sample):
        return np.nan, np.nan, "تركيز العينة غير صالح للحساب"
    
    volume_sample = (target_conc * final_vol) / conc_sample
    
    if volume_sample > final_vol:
        volume_sample = final_vol
        volume_diluent = 0
        note = "العينة مخففة جدًا للوصول للتركيز الهدف — يلزم تركيز العينة أو استخدام كمية أكبر."
    else:
        volume_diluent = final_vol - volume_sample
    
    return round(volume_sample, 2), round(volume_diluent, 2), note
