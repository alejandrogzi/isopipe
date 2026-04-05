# syntax=docker/dockerfile:1.7
# version 0.0.2

# ---------- Runtime Stage ----------
FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    procps \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
    "intronIC==2.1.4" \
    "numpy >= 2" \
    "pandas >= 2" \
    "scipy >= 1.2.1" \
    "tensorflow >= 1.13.1" \
    "keras >= 2.2.4" \
    "scikit-learn>=0.22" \
    biogl \
    matplotlib

# APARENT
COPY assets/py/aparent/aparent.py /usr/local/lib/aparent/aparent.py
RUN ln -s /usr/local/lib/aparent/aparent.py /usr/local/bin/aparent \
    && chmod +x /usr/local/lib/aparent/aparent.py

# VEREDICT
COPY assets/py/veredict/veredict.py /usr/local/lib/veredict/veredict.py
RUN ln -s /usr/local/lib/veredict/veredict.py /usr/local/bin/veredict \
    && chmod +x /usr/local/lib/veredict/veredict.py 

RUN python -c "import pandas; import numpy; import scipy; print('All dependencies OK')"
RUN intronIC --version
RUN aparent --version
RUN veredict --version
