import streamlit as st
import plotly.graph_objects as go
import numpy as np
from scipy.interpolate import interp1d

class TheoryTab:
    def __init__(self):
        self.render()

    def render(self):
        st.header("Theory Tab Content")

        # Center the image
        st.markdown(
            """
            <style>
            .center {
                display: block;
                margin-left: auto;
                margin-right: auto;
                width: 50%;
            }
            </style>
            """,
            unsafe_allow_html=True
        )
        st.image("The-graph-depicting-the-RNA-purity-measured-by-NanoDrop-1000-vs-52.png", caption="Typical nucleic acid spectrum (DOI: 10.1371/journal.pone.0291949)", width=400, use_container_width=True, output_format='PNG', clamp=False)

        # Add a horizontal line
        st.markdown("---")

        # Explanation
        st.markdown("This tab demonstrates the relationship between absorbance values at different wavelengths (230nm, 260nm, 280nm) and the purity of a nucleic acid sample. Adjust the sliders to see how changes in these values affect the overall spectrum.")

        # Create columns for sliders and chart
        col1, col2 = st.columns([1, 2])

        with col1:
            st.markdown("#### Sliders")
            st.markdown(
                """
                <style>
                .stSlider > div > div > div > div {
                    display: none;
                }
                .stSlider > div > div::before {
                    content: attr(aria-valuenow);
                    position: absolute;
                    top: -20px;
                    left: 0;
                    font-size: 14px;
                    color: #888;
                }
                </style>
                """,
                unsafe_allow_html=True,
            )
            salt = st.slider("Salt", 0.0, 2.0, 0.5)
            dna_rna = st.slider("DNA/RNA", 0.0, 8.0, 3.0)
            protein = st.slider("Protein", 0.0, 4.0, 1.0)
            factor = st.slider("Factor", 0.0, 100.0, 50.0)

            ratio_260_280 = dna_rna / protein if protein != 0 else 0
            ratio_260_230 = dna_rna / salt if salt != 0 else 0
            concentration = dna_rna * factor

            st.markdown(f"260/280: <span style='background-color: rgba(0, 0, 0, 0.1); padding: 2px;'>{ratio_260_280:.2f}</span>", unsafe_allow_html=True)
            st.markdown(f"260/230: <span style='background-color: rgba(0, 0, 0, 0.1); padding: 2px;'>{ratio_260_230:.2f}</span>", unsafe_allow_html=True)
            st.markdown(f"Factor: <span style='background-color: rgba(0, 0, 0, 0.1); padding: 2px;'>{factor:.2f}</span>", unsafe_allow_html=True)
            st.markdown(f"التركيز: <span style='background-color: rgba(0, 0, 0, 0.1); padding: 2px;'>{concentration:.2f} ng/µL</span>", unsafe_allow_html=True)

        # Sample data based on the image
        x = [230, 260, 280]
        y = [salt, dna_rna, protein]

        # Add data points to connect to zero and extend to 340
        x = [220] + x + [340]
        y = [0] + y + [0]

        # Create a spline interpolation
        f = interp1d(x, y, kind='quadratic')

        # Generate x values for the curve
        xnew = np.linspace(min(x), max(x), num=100, endpoint=True)

        # Generate y values for the curve
        ynew = f(xnew)

        # Create the plot
        fig = go.Figure()

        # Add the scatter trace
        fig.add_trace(go.Scatter(x=xnew, y=ynew, mode='lines', line=dict(color='red')))

        # Add labels for the wavelengths
        fig.add_trace(go.Scatter(x=[230, 260, 280], y=[salt, dna_rna, protein], mode='markers+text', text=["Salt", "DNA/RNA", "Protein"], textposition="top center",
                                 textfont=dict(size=16)))


        # Set the axis labels
        fig.update_layout(
            xaxis_title="Wavelength (nm)",
            yaxis_title="10mm Absorbance",
            xaxis_range=[220, 340],
            yaxis_range=[0, 7]
        )

        # Display the plot in Streamlit
        with col2:
            st.plotly_chart(fig, key="updated_chart")