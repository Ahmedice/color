import numpy as np
import pandas as pd
from utils import safe_float, calculate_ratios, compute_concentration, assess_purity, calculate_dilution

def process_dataframe(df: pd.DataFrame, detected_columns: dict, protocol_settings: dict, default_factor: float):
    # أعد تسمية الأعمدة حسب المعيار الداخلي
    col_map = {
        detected_columns.get("sample_id", "Sample ID"): "Sample ID",
        detected_columns.get("sample_type", "Sample Type"): "Sample Type",
        detected_columns.get("a260", "A260 (Abs)"): "a260",
        detected_columns.get("a280", "A280 (Abs)"): "a280",
        detected_columns.get("a230", "A230 (Abs)"): "a230",
        detected_columns.get("ratio_260_230", "ratio_260_230"): "ratio_260_230",
        detected_columns.get("concentration", "Conc_ng_per_ul"): "Conc_ng_per_ul",
    }
    df = df.rename(columns=col_map)

    # تأكد من وجود الأعمدة المطلوبة
    for c in ["Sample ID", "Sample Type", "a260", "a280", "a230", "ratio_260_230"]:
        if c not in df.columns:
            df[c] = np.nan

    results = []
    for idx, row in df.iterrows():
        sample_id = row.get("Sample ID", "")
        sample_type = row.get("Sample Type", "")

        # استخرج القيم الرقمية الآمنة
        a260, _ = safe_float(row.get("a260"))
        a280, _ = safe_float(row.get("a280"))
        a230, _ = safe_float(row.get("a230"))
        ratio_260_230_input, _ = safe_float(row.get("ratio_260_230"))
        
        # حساب النسب
        ratio_260_280, ratio_260_230 = calculate_ratios({
            "a260": a260,
            "a280": a280,
            "a230": a230,
            "ratio_260_230": ratio_260_230_input
        })

        # تركيز العينة
        factor = default_factor
        target_conc = protocol_settings.get("target_conc", 10)
        final_vol = protocol_settings.get("final_vol", 20)

        concentration = compute_concentration(row, factor)

        # تقييم نقاوة العينة
        purity, purity_notes = assess_purity(sample_type, ratio_260_280, ratio_260_230)

        # حساب التخفيف المطلوب
        vol_sample, vol_diluent, dilution_note = calculate_dilution(concentration, target_conc, final_vol)

        # تجميع النتائج للصف الحالي
        result_row = {
            "Sample ID": sample_id,
            "Sample Type": sample_type,
            "A260 (Abs)": a260,
            "A280 (Abs)": a280,
            "A230 (Abs)": a230,
            "Conc_ng_per_ul": round(concentration, 2) if not np.isnan(concentration) else np.nan,
            "260/280": round(ratio_260_280, 2) if not np.isnan(ratio_260_280) else np.nan,
            "260/230": round(ratio_260_230, 2) if not np.isnan(ratio_260_230) else np.nan,
            "Purity": purity,
            "Purity_notes": purity_notes,
            "Volume Sample (µl)": vol_sample,
            "Volume Diluent (µl)": vol_diluent,
            "Notes": dilution_note
        }
        results.append(result_row)

    return pd.DataFrame(results)