import streamlit as st
import plotly.graph_objects as go
import numpy as np
from scipy.interpolate import interp1d
from utils import assess_purity

class TheoryTab:
    def __init__(self):
        self.render()

    def render(self):
        st.header("النظرية: تحليل نقاوة العينات حسب النوع")

        # تقسيم الشاشة لعمودين: اليمين للسلايدر، اليسار للرسم البياني
        col_sliders, col_chart = st.columns([1, 2])

        with col_sliders:
            # اختيار نوع العينة
            sample_type = st.selectbox("اختر نوع العينة", ["DNA", "RNA", "Protein"])

            # Sliders
            salt = st.slider("Salt (230 nm Absorbance)", 0.0, 4.0, 0.5)
            dna_rna = st.slider("DNA/RNA (260 nm Absorbance)", 0.0, 8.0, 3.0)
            protein = st.slider("Protein (280 nm Absorbance)", 0.0, 4.0, 1.0)
            factor = st.slider("Factor (معامل الحساب)", 0.0, 100.0, 50.0)

            # حساب القيم
            ratio_260_280 = dna_rna / protein if protein != 0 else 0
            ratio_260_230 = dna_rna / salt if salt != 0 else 0
            concentration = dna_rna * factor

            purity_verdict, purity_reco = assess_purity(ratio_260_280, ratio_260_230, sample_type)

            # تنسيق القيم في صفين، كل صف به 3 أعمدة (نص عربي + القيمة)
            st.markdown("---")

            # الصف الأول: نوع العينة، نسبة 260/280، نسبة 260/230
            row1_col1, row1_col2, row1_col3 = st.columns(3)
            row1_col1.markdown(f"**نوع العينة:** {sample_type}")
            row1_col2.markdown(f"**نسبة 260/280:** {ratio_260_280:.2f}")
            row1_col3.markdown(f"**نسبة 260/230:** {ratio_260_230:.2f}")

            # الصف الثاني: Factor، التركيز، النقاء (تحذير مع أيقونة)
            row2_col1, row2_col2, row2_col3 = st.columns(3)
            row2_col1.markdown(f"**Factor:** {factor:.2f}")
            row2_col2.markdown(f"**التركيز (ng/µL):** {concentration:.2f}")
            row2_col3.markdown(f"**النقاء:** ℹ️ {purity_verdict}")

            # تعليق تحت الصفوف بشكل كامل
            st.markdown(f"**التعليق:** {purity_reco}")

        with col_chart:
            # بيانات للرسم البياني (المثال مشابه للصورة)
            x = [230, 260, 280]
            y = [salt, dna_rna, protein]
            x = [220] + x + [340]
            y = [0] + y + [0]
            f = interp1d(x, y, kind='quadratic')
            xnew = np.linspace(min(x), max(x), num=100, endpoint=True)
            ynew = f(xnew)

            fig = go.Figure()
            fig.add_trace(go.Scatter(x=xnew, y=ynew, mode='lines', line=dict(color='red')))
            fig.add_trace(go.Scatter(
                x=[230, 260, 280], y=[salt, dna_rna, protein],
                mode='markers+text',
                text=["Salt", "DNA/RNA", "Protein"],
                textposition="top center",
                textfont=dict(size=16)
            ))
            fig.update_layout(
                xaxis_title="Wavelength (nm)",
                yaxis_title="10mm Absorbance",
                xaxis_range=[220, 340],
                yaxis_range=[0, 7],
                margin=dict(t=30, b=20, l=50, r=30)
            )
            st.plotly_chart(fig, use_container_width=True)
