FROM continuumio/miniconda3:latest

WORKDIR /app

# ── 1. Create the conda environment ──────────────────────────────────────────
COPY conda_env.yaml .
RUN conda env create -f conda_env.yaml && \
    conda clean --all --yes

# ── 2. Copy source and install the mitoedit package ──────────────────────────
COPY . .
RUN conda run --no-capture-output -n mitoedit pip install --no-deps .

# ── 3. Runtime directories ───────────────────────────────────────────────────
# final_output and running are written to at request time
RUN mkdir -p /app/certs /app/final_output /app/running

# ── 4. Ports ─────────────────────────────────────────────────────────────────
# Port 80  — plain HTTP (default, or when PORT env var is 80)
# Port 443 — HTTPS (set PORT=443 and mount certs at /app/certs)
EXPOSE 80
EXPOSE 443

# ── 5. Environment variables ─────────────────────────────────────────────────
# PORT: which port uvicorn listens on (default 80; use 8000 for local dev)
ENV PORT=80
# ALLOWED_HOST: validated host for HTTPS redirect; unused in plain-HTTP mode
ENV ALLOWED_HOST=mitoedit.stjude.org
# SSL: paths are fixed; mount real certs at runtime for HTTPS:
#   docker run -v /path/to/certs:/app/certs ...
# If the files are absent the app falls back to plain HTTP automatically.
ENV SSL_CERTFILE=/app/certs/mitoedit.pem
ENV SSL_KEYFILE=/app/certs/mitoedit.key
# MITOEDIT_PASSWORD must be supplied at runtime — never bake it in:
#   docker run -e MITOEDIT_PASSWORD=your_secure_password ...

# ── 6. Entrypoint ─────────────────────────────────────────────────────────────
CMD ["conda", "run", "--no-capture-output", "-n", "mitoedit", \
     "python", "-m", "mitoedit.web.main"]
