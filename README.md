# Bioninformatics Unified Platform
## Installation Guide
1. python3 -m venv env/
2. source env/bin/activate
3. python3 -m pip install --upgrade pip
4. pip3 install wheel
5. pip3 install absl-py
6. pip install "apache-airflow[celery]==2.10.2" --constraint "https://raw.githubusercontent.com/apache/airflow/constraints-2.10.2/constraints-3.8.txt"

On Windows:
1. python3 -m venv env/
2. env\Scripts\activate
3. python3 -m pip install --upgrade pip
4. pip3 install wheel
5. pip3 install absl-py
6. pip install "apache-airflow[celery]==2.10.2" --constraint "https://raw.githubusercontent.com/apache/airflow/constraints-2.10.2/constraints-3.8.txt"
