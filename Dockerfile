FROM python:3.10-slim

# System libraries required for scientific stack and CoreMS
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    gfortran \
    git \
    libpq-dev \
    libhdf5-dev \
    libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Install CoreMS with a fixed version + required base libraries
# Important: pin the CoreMS version to avoid unexpected changes
RUN pip install --no-cache-dir \
    "git+https://github.com/EMSL-Computing/CoreMS.git@v3.9.3" \
    "psycopg2-binary==2.9.10" \
    "matplotlib==3.7.5" \
    "pandas==1.5.3" \
    "streamlit==1.51.0"



# Copy the entire project into /app
COPY . .

EXPOSE 8501

# ENV PYTHONPATH=/app/src
ENV PYTHONPATH=/app/src
# By default, run generation of _peaks.csv from .d folders
CMD ["python", "-m", "SCHON.corems_pipeline.generate_peaks_from_fid"]