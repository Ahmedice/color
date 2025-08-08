# app.py
"""
NanoDrop Analyzer - Streamlit App
How to run:
1. pip install -r requirements.txt
2. streamlit run app.py

Notes:
- All internal variable names and comments are in English.
- All user-facing labels/messages are in Arabic.
- The app loads protocols from config.json at startup.
"""

import json
import io
from typing import Dict, Tuple, Any
import math

import pandas as pd
import numpy as np
import streamlit as st

# -------------------------
# Utility & Core Functions
# -------------------------

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
    r260_280, r260_230, ratio_note = calculate_ratios(a260, a280, ratio_260_230)

    # Dilution
    v1, v2, dilution_note = calculate_dilution(concentration, target_conc, final_vol)

    # Purity
    purity_verdict, purity_reco = assess_purity(r260_280, r260_230, sample_type)

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


def safe_ratio(x) -> Tuple[float, str]:
   """Try to convert to float. Return (value, note)."""
   if pd.isna(x):
       return (np.nan, "missing")
   try:
       val = float(x)
       return (val, "")
   except Exception:
       return (np.nan, "non-numeric")

# -------------------------
# Streamlit App UI
# -------------------------

st.set_page_config(page_title="NanoDrop Analyzer", layout="wide")

# Load protocols
try:
    protocols = load_config("config.json")
except Exception as e:
    st.error("فشل تحميل ملف الإعدادات config.json. تأكد من وجود الملف. Error: " + str(e))
    protocols = {
        "PCR": {"target_conc": 10, "final_vol": 20},
        "qPCR": {"target_conc": 5, "final_vol": 10},
        "Sanger": {"target_conc": 15, "final_vol": 12}
    }

st.title("محلل NanoDrop")
st.markdown("<style>body {direction: rtl;}</style>", unsafe_allow_html=True)
st.sidebar.header("شرح التطبيق")
st.sidebar.markdown(
    """
    شرح التطبيق:
    تطبيق "NanoDrop Analyzer" هو أداة لتحليل بيانات جهاز NanoDrop لقياس تركيز ونقاء العينات الحيوية (DNA, RNA, Protein).
    التطبيق يسمح بتحليل عينة واحدة يدوياً أو تحميل دفعات عينات من ملف CSV لتحليلها دفعة واحدة.
    يعتمد التطبيق على حساب تركيز العينة استناداً إلى قراءة الامتصاص عند 260 نانومتر مع معامل تحويل (factor)، وحساب نسب النقاء 260/280 و260/230.
    بالإضافة إلى ذلك، يحسب التطبيق التخفيف المطلوب للوصول إلى تركيز الهدف وحجم المحلول النهائي حسب البروتوكولات المختارة (PCR, qPCR, Sanger).
    يتم تقييم نقاوة العينة وفق معايير علمية مع عرض توصيات واضحة لتوجيه المستخدم بخصوص جودة العينات.
    التطبيق يعمل محلياً وبواجهة عربية سهلة الاستخدام، مع إمكانية تحميل النتائج النهائية للعينات المحللة.
    
    **ملاحظات مهمة:**
    - جهاز NanoDrop عادة لا يعطي قراءة A230 مباشرة، لذلك التطبيق يدعم إدخال نسبة 260/230 أو حسابها تقديرياً.
    - تأكد من إدخال القيم بدقة لأن النتائج تعتمد عليها بشكل كبير.
    - التطبيق يعمل محلياً بالكامل ولا يحتاج اتصال بالإنترنت.
    - الكود مكتوب بطريقة منظمة لتسهيل التعديل وإضافة بروتوكولات أو تحسينات مستقبلية.
    
    **طريقة التشغيل:**
    1. ثبت المتطلبات عبر:
       `pip install -r requirements.txt`
    2. شغل التطبيق:
       `streamlit run app.py`
    
    لأي استفسار أو اقتراح، يرجى التواصل مع المطور.
    """
)

st.sidebar.header("Config.json")
st.sidebar.json(protocols, expanded=False)

tab1, tab2 = st.tabs(["عينة واحدة", "ملف CSV كامل"])

# -------------------------
# TAB 1: Single Sample
# -------------------------
with tab1:
    st.header("تحليل عينة واحدة")
    with st.form("single_sample_form"):
        st.subheader("نوع الاختبار")
        col1_test, col2_test, col3_test = st.columns([1, 1, 1])
        with col1_test:
            sample_id = st.text_input("Sample ID", value="")
        with col2_test:
            sample_type = st.selectbox("Sample Type", options=["DNA", "RNA", "Protein"])
        with col3_test:
            protocol_choice = st.selectbox("Protocol", options=list(protocols.keys()))

        st.subheader("القرائات")
        col1_readings, col2_readings, col3_readings = st.columns([1, 1, 1])
        with col1_readings:
            a260_input = st.text_input("A260", value="", help="مثال: 0.12")
        with col2_readings:
            a280_input = st.text_input("A280", value="")
        with col3_readings:
            ratio_260_230_input = st.text_input(
                "260/230",
                value="", help="راجع جهازك أو نتائجك التجريبية."
            )

        st.subheader("التخفيف")
        col1_dilution, col2_dilution, col3_dilution = st.columns([1, 1, 1])
        with col1_dilution:
            factor_input = st.text_input("Factor", value="", help="A260 لحساب التركيز.")
        with col2_dilution:
            final_vol_input = st.text_input("Final Vol (µl)", value="")
        with col2_dilution:
            final_vol_input = st.text_input("Final Vol (µl)", value="", help=f"القيمة الافتراضية حسب البروتوكول: {protocols.get(protocol_choice, {}).get('final_vol', 20)} µl")
        with col3_dilution:
            target_conc_input = st.text_input("Target Conc", value="", help=f"القيمة الافتراضية حسب البروتوكول: {protocols.get(protocol_choice, {}).get('target_conc', 10)} ng/µl")

        # load defaults from protocol
        default_final_vol = protocols.get(protocol_choice, {}).get("final_vol", 20)
        default_target_conc = protocols.get(protocol_choice, {}).get("target_conc", 10)

        submitted = st.form_submit_button("حساب")

    if submitted:
        # Convert inputs safely
        a260_val, note_a260 = safe_float(a260_input)
        a280_val, note_a280 = safe_float(a280_input)
        ratio_260_230_val, note_260_230 = safe_float(ratio_260_230_input)

        # Compute
        factor_val, note_factor = safe_float(factor_input)
        final_vol_val, note_final_vol = safe_float(final_vol_input)
        target_conc_val, note_target_conc = safe_float(target_conc_input)

        inputs_dict = {
            "a260": a260_val,
            "a280": a280_val,
            "ratio_260_230": ratio_260_230_val,
            "factor": factor_val,
            "sample_type": sample_type,
            "target_conc": target_conc_val,
            "final_vol": final_vol_val
        }

        # load defaults from protocol
        target_conc_protocol, _ = safe_float(protocols.get(protocol_choice, {}).get("target_conc", 10))
        final_vol_protocol, _ = safe_float(protocols.get(protocol_choice, {}).get("final_vol", 20))
        
        protocol_settings = {"target_conc": target_conc_protocol, "final_vol": final_vol_protocol}
        results = process_single_sample(inputs_dict, protocol_settings)

        # Display results in two columns
        left, right = st.columns(2)
        with left:
            st.subheader("النتائج العددية")
            st.markdown(f"- **تركيز (ng/µl):** `{results['Concentration_ng_per_ul']}`")
            st.markdown(f"- **نسبة 260/280:** `{results['260/280']}`")
            st.markdown(f"- **نسبة 260/230:** `{results['260/230']}`")
            st.markdown(f"- **حكم النقاء:** `{results['Purity']}`")
            st.markdown(f"- **ملاحظات النقاء:** `{results['Purity_notes']}`")
            if results["Notes"]:
                st.warning("ملاحظات: " + results["Notes"])
        with right:
            st.subheader("خطوات التخفيف / التحضير")
            if isinstance(results["V1_ul"], float) and not math.isnan(results["V1_ul"]):
                if results["V1_ul"] > results["Protocol_final_vol_ul"]:
                    st.error("العينة مخففة جدًا للوصول للتركيز الهدف — يلزم تركيز العينة أو استخدام مزيد من العينة.")
                else:
                    st.markdown(f"خذ **{results['Volume Sample (µl)']} µl** من العينة + **{results['Volume Diluent (µl)']} µl** ماء معقم للوصول إلى **{results['Protocol_target_ng_per_ul']} ng/µl** في **{results['Protocol_final_vol_ul']} µl**.")
            else:
                st.error("لا يمكن حساب V1/V2 بسبب تركيز غير صالح.")
        with st.expander("تفاصيل الإدخال"):
            st.write({
                "Sample ID": sample_id or "(لم يُعطَ)",
                "A260": a260_input,
                "A280": a280_input,
                "260/230 Ratio": ratio_260_230_input or "(مفقود)",
                "Factor": factor_input,
                "Protocol": protocol_choice,
                "Final Volume (µl)": final_vol_input,
                "Target Concentration (ng/µl)": target_conc_input
            })

# -------------------------
# TAB 2: CSV Batch
# -------------------------
with tab2:
    st.header("تحميل ملف CSV للمعالجة الدفعيّة")
    st.markdown("**ملاحظة:** تنسيق الأعمدة المفضّل: `Sample ID, Sample Type, A260 (Abs), A280 (Abs), A230 (Abs) (اختياري), Nucleic Acid (concentration)`")
    with st.expander("مثال ملف CSV سريع (يمكن نسخه ولصقه):"):
        st.code("""Sample ID,Sample Type,A260 (Abs),A280 (Abs),A230 (Abs)
DNA_1,DNA,0.12,0.07,0.06
RNA_1,RNA,0.08,0.04,
PROT_1,Protein,0.02,0.03,0.01
""", language="csv")

    uploaded_file = st.file_uploader("اختر ملف CSV لتحميله", type=["csv"])
    protocol_choice_batch = st.selectbox("اختر البروتوكول لكل العينات", options=list(protocols.keys()), key="protocol_batch")
    protocol_settings_batch = protocols.get(protocol_choice_batch, {"target_conc": 10, "final_vol": 20})
    default_factor_batch = 50.0

    if uploaded_file is not None:
        try:
            df = pd.read_csv(uploaded_file)
        except Exception as e:
            st.error("فشل قراءة ملف CSV: " + str(e))
            df = None

        if df is not None:
            st.success(f"تم تحميل الملف (الصفوف: {len(df)})")
            detected = map_columns(df.columns)
            st.write("اكتشاف أعمدة (مطابقة تقريبية):")
            st.json(detected)
            
            # look for 260/230 column
            possible_260_230_cols = [col for col in df.columns if "260/230" in col.lower()]
            if possible_260_230_cols:
               ratio_col = possible_260_230_cols[0]
               df['ratio_260_230'] = pd.to_numeric(df[ratio_col], errors='coerce')
               detected["ratio_260_230"] = ratio_col
            else:
               df['ratio_260_230'] = np.nan
               detected["ratio_260_230"] = None

            if detected["ratio_260_230"] is not None:
                result_df = process_dataframe(df, detected, protocol_settings_batch, default_factor_batch)
            else:
                st.warning("Column 'ratio_260_230' not found in CSV file. Please include this column for accurate results.")
                result_df = process_dataframe(df, detected, protocol_settings_batch, default_factor_batch)
            st.subheader("النتائج")
            st.dataframe(result_df, use_container_width=True)

            # Download button
            csv_bytes = result_df.to_csv(index=False).encode("utf-8")
            st.download_button(label="تحميل النتائج كـ CSV", data=csv_bytes, file_name="nanodrop_results.csv", mime="text/csv")
    else:
        # Show sample CSV content and a quick button to load sample into app for demo
        st.info("يمكنك تحميل ملف CSV أو استخدام المثال أسفل لتحميل سريع للاختبار.")
        if st.button("تحميل المثال إلى المعالج"):
            sample_csv = """Sample ID,Sample Type,A260 (Abs),A280 (Abs),A230 (Abs)
DNA_1,DNA,0.12,0.07,0.06
RNA_1,RNA,0.08,0.04,
PROT_1,Protein,0.02,0.03,0.01"""
            df_sample = pd.read_csv(io.StringIO(sample_csv))
            detected = map_columns(df_sample.columns)
            result_df = process_dataframe(df_sample, detected, protocol_settings_batch, default_factor_batch)
            st.dataframe(result_df, use_container_width=True)
            csv_bytes = result_df.to_csv(index=False).encode("utf-8")
            st.download_button(label="تحميل نتائج المثال كـ CSV", data=csv_bytes, file_name="nanodrop_example_results.csv", mime="text/csv")

# -------------------------
# Footer / Small Help
# -------------------------
st.markdown("---")
st.markdown("**ملاحظات:** التطبيق يعمل محليًا. لتعديل بروتوكول افتراضي حرّر ملف `config.json` ثم أعد تشغيل التطبيق.")
