FROM python:3.11-slim
WORKDIR /app
RUN apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*
COPY requirements.txt pyproject.toml ./
RUN pip install --no-cache-dir -r requirements.txt
# The engine imports all three modules; local_db/ carries the pinned
# reference data (TargetScan zip, miRBase FASTAs, miRTarBase GMT, GMTs).
COPY kenryu_engine.py kenryu_core.py kenryu_updater.py kenryu_cli.py ./
COPY data/ ./data/
COPY local_db/ ./local_db/
COPY static/ ./static/
RUN mkdir -p local_db && chmod 777 local_db
RUN useradd -m kenryu && chown -R kenryu:kenryu /app
USER kenryu
EXPOSE 7860
HEALTHCHECK --interval=60s --timeout=10s --start-period=30s \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:7860/health', timeout=5)"
CMD ["uvicorn", "kenryu_engine:app", "--host", "0.0.0.0", "--port", "7860"]
