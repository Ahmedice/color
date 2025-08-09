import io
import pandas as pd
import numpy as np
import streamlit as st
from utils import map_columns, safe_float, calculate_ratios, compute_concentration, assess_purity, calculate_dilution
from dataframe_processor import process_dataframe

def batch_processing_tab(protocols):
    st.header("تحميل ملف CSV للمعالجة ")
    st.markdown("**ملاحظة:** تنسيق الأعمدة المفضّل: `Sample ID, Sample Type, A260 (Abs), A280 (Abs), A230 (Abs) (اختياري), Nucleic Acid (concentration)`")
    with st.expander("مثال ملف CSV سريع (يمكن نسخه ولصقه):"):
        st.code("""Sample ID,Sample Type,A260 (Abs),A280 (Abs),A230 (Abs)
DNA_1,DNA,0.28,0.16,0.12
DNA_2,DNA,0.10,0.07,0.01
DNA_3,DNA,0.45,0.50,0.15
RNA_1,RNA,0.22,0.13,0.08
RNA_2,RNA,0.15,0.05,0.03
PROT_1,Protein,0.05,0.32,0.03
PROT_2,Protein,0.02,0.10,0.01
PROT_3,Protein,0.04,0.03,0.00""", language="csv")

    uploaded_file = st.file_uploader("اختر ملف CSV أو Excel لتحميله", type=["csv", "xlsx"])
    protocol_choice_batch = st.selectbox("اختر البروتوكول لكل العينات", options=list(protocols.keys()), key="protocol_batch")
    protocol_settings_batch = protocols.get(protocol_choice_batch, {"target_conc": 10, "final_vol": 20})
    default_factor_batch = 50.0

    if uploaded_file is not None:
        try:
            if uploaded_file.name.endswith('.csv'):
                df = pd.read_csv(uploaded_file)
            else:
                df = pd.read_excel(uploaded_file)
        except Exception as e:
            st.error("فشل قراءة ملف CSV: " + str(e))
            return

        st.success(f"تم تحميل الملف (الصفوف: {len(df)})")

        detected, missing = map_columns(df.columns)
        st.write("اكتشاف أعمدة (مطابقة تقريبية):")
        st.json(detected)
        if "a260" in missing or "sample_type" in missing:
            st.warning("أعمدة A260 أو Sample Type مفقودة. يرجى التأكد من وجودها في الملف.")

        # التعامل مع عمود 260/230
        possible_260_230_cols = [col for col in df.columns if "260/230" in col.lower()]
        if possible_260_230_cols:
            ratio_col = possible_260_230_cols[0]
            df['ratio_260_230'] = pd.to_numeric(df[ratio_col], errors='coerce')
            detected["ratio_260_230"] = ratio_col
        else:
            df['ratio_260_230'] = np.nan
            detected["ratio_260_230"] = None
            st.info("لم يتم العثور على عمود '260/230' في الملف. سيتم التعامل مع هذه القيمة كغير موجودة.")

        # معالجة الداتا بإضافة التقييم والتخفيف
        result_df = process_dataframe(df, detected, protocol_settings_batch, default_factor_batch)

        # عرض الجدول النهائي مع كل القيم المضافة
        st.subheader("النتائج")
        st.dataframe(result_df)

        # عرض ملخصات تحذيرية للمستخدم (مثلاً عدد العينات التي بها تحذيرات نقاوة)
        warnings_count = result_df[result_df['Purity'] == 'تحذير'].shape[0]
        if warnings_count > 0:
            st.warning(f"هناك {warnings_count} عينة بحاجة إلى مراجعة نقاوتها حسب التقييم.")

        # زر تحميل النتائج مع كل الأعمدة الجديدة
        if uploaded_file.name.endswith('.csv'):
            csv_bytes = result_df.to_csv(index=False).encode("utf-8")
            st.download_button(label="تحميل النتائج كـ CSV", data=csv_bytes, file_name="nanodrop_results.csv", mime="text/csv")
        else:
            excel_bytes = io.BytesIO()
            result_df.to_excel(excel_bytes, index=False)
            excel_bytes = excel_bytes.getvalue()
            st.download_button(label="تحميل النتائج كـ Excel", data=excel_bytes, file_name="nanodrop_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    else:
        # عرض مثال للتحميل السريع
        st.info("يمكنك تحميل ملف CSV أو استخدام المثال أسفل لتحميل سريع للاختبار.")
        if st.button("تحميل المثال إلى المعالج"):
            sample_csv = """Sample ID,Sample Type,A260 (Abs),A280 (Abs),A230 (Abs)
DNA_1,DNA,0.28,0.16,0.12
DNA_2,DNA,0.10,0.07,0.01
DNA_3,DNA,0.45,0.50,0.15
RNA_1,RNA,0.22,0.13,0.08
RNA_2,RNA,0.15,0.05,0.03
PROT_1,Protein,0.05,0.32,0.03
PROT_2,Protein,0.02,0.10,0.01
PROT_3,Protein,0.04,0.03,0.00"""
            df_sample = pd.read_csv(io.StringIO(sample_csv))
            detected, _ = map_columns(df_sample.columns)
            result_df = process_dataframe(df_sample, detected, protocol_settings_batch, default_factor_batch)
            st.dataframe(result_df)
            csv_bytes = result_df.to_csv(index=False).encode("utf-8")
            st.download_button(label="تحميل نتائج المثال كـ CSV", data=csv_bytes, file_name="nanodrop_example_results.csv", mime="text/csv")