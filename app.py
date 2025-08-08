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

import io
from typing import Dict, Tuple, Any
import math

import pandas as pd
import numpy as np
import streamlit as st
from utils import load_config, map_columns, safe_float, safe_ratio
from single_sample import process_single_sample
from dataframe_processor import process_dataframe

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

# تهيئة session_state للقيم الافتراضية عند أول تحميل
if "protocol_choice" not in st.session_state:
    st.session_state.protocol_choice = "PCR"
    st.session_state.target_conc = protocols["PCR"]["target_conc"]
    st.session_state.final_vol = protocols["PCR"]["final_vol"]
    st.session_state.factor = 50

# دالة لتحديث القيم عند تغيير البروتوكول
def on_protocol_change():
    st.session_state.target_conc = protocols[st.session_state.protocol_choice]["target_conc"]
    st.session_state.final_vol = protocols[st.session_state.protocol_choice]["final_vol"]
    st.session_state.factor = 50

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
    col1_test, col2_test, col3_test = st.columns([1, 1, 1])
    with col1_test:
        sample_id = st.text_input("Sample ID", value="")
    with col2_test:
        sample_type = st.selectbox("Sample Type", options=["DNA", "RNA", "Protein"])
    with col3_test:
        protocol_choice = st.selectbox(
            "Protocol",
            options=list(protocols.keys()),
            key="protocol_choice",
            on_change=on_protocol_change
        )
    with st.form("single_sample_form"):
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
            final_vol_input = st.text_input("Final Vol (µl)", value=str(st.session_state.final_vol), key="final_vol_input")
        with col3_dilution:
            target_conc_input = st.text_input("Target Conc", value=str(st.session_state.target_conc), key="target_conc_input")
        with col1_dilution:
            factor_input = st.text_input("Factor", value=str(st.session_state.factor), help="A260 لحساب التركيز.", key="factor_input")

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
            if isinstance(results["Volume Sample (µl)"], float) and not math.isnan(results["Volume Sample (µl)"]):
                if results["Volume Sample (µl)"] > results["Protocol_final_vol_ul"]:
                    st.error("العينة مخففة جدًا للوصول للتركيز الهدف — يلزم تركيز العينة أو استخدام مزيد من العينة.")
                else:
                    st.markdown(f"خذ **{results['Volume Sample (µl)']} µl** من العينة + **{results['Volume Diluent (µl)']} µl** ماء معقم للوصول إلى **{results['Protocol_target_ng_per_ul']} ng/µl** في **{results['Protocol_final_vol_ul']} µl**.")
            else:
                st.error("لا يمكن حساب V1/V2 بسبب تركيز غير صالح.")
        with st.expander("تفاصيل الإدخال"):
            st.write({
                "Sample ID": sample_id,
                "Sample Type": sample_type,
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
PROT_1,Protein,0.02,0.03,0.01""", language="csv")

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
