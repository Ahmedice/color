# NanoDrop Analyzer

This is a Streamlit application for analyzing NanoDrop data. It allows users to analyze single samples or batch samples from a CSV or Excel file.

## How to run

1.  Install the requirements:
    ```bash
    pip install -r requirements.txt
    ```
2.  Run the application:
    ```bash
    streamlit run app.py
    ```

## Files

*   `app.py`: Main Streamlit application file.
*   `batch_processor.py`: Contains the batch processing tab logic.
*   `dataframe_processor.py`: Contains functions for processing dataframes.
*   `utils.py`: Contains utility functions such as loading the configuration file and assessing purity.
*   `config.json`: Configuration file containing default protocols.
*   `requirements.txt`: List of required Python packages.
*   `single_sample.py`: Contains the single sample processing logic.
*   `theory_tab.py`: Contains the theory tab logic and plot.

## Contact

For any questions or suggestions, please contact [Black.ice.onet@gmail.com](mailto:Black.ice.onet@gmail.com).

## شرح إضافي (Arabic Explanation)

تطبيق NanoDrop Analyzer هو أداة لتحليل بيانات NanoDrop. يمكنك من خلال هذا التطبيق تحليل عينة واحدة أو مجموعة من العينات من خلال ملف CSV أو Excel.

## كيفية التشغيل (How to Run - Arabic)

1.  تثبيت المتطلبات:
    ```bash
    pip install -r requirements.txt
    ```
2.  تشغيل التطبيق:
    ```bash
    streamlit run app.py
    ```