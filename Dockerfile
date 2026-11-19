# H-22 (Ronda 3, v1.49): base image pinned by digest (reproducible
# builds: the tag 3.11-slim alone can move); build-essential is purged
# after pip install so no compiler ships in the final image; local_db
# is 755, not 777 (the app user owns it, nobody else needs write).
FROM python:3.11-slim@sha256:07d4824036af339fa1b6ef769cc99d9c10a141769538b52b3b597c84f7a7b535
WORKDIR /app
RUN apt-get update && apt-get install -y --no-install-recommends build-essential \
  && rm -rf /var/lib/apt/lists/*
COPY requirements.txt pyproject.toml ./
RUN pip install --no-cache-dir -r requirements.txt \
  && apt-get purge -y build-essential && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/*
# The engine imports all three modules; local_db/ carries the pinned
# reference data (TargetScan zip, miRBase FASTAs, miRTarBase GMT, GMTs).
COPY kenryu_engine.py kenryu_core.py kenryu_updater.py kenryu_cli.py ./
COPY data/ ./data/
COPY local_db/ ./local_db/
COPY static/ ./static/
RUN mkdir -p local_db && chmod 755 local_db
RUN useradd -m kenryu && chown -R kenryu:kenryu /app
USER kenryu
EXPOSE 7860
HEALTHCHECK --interval=60s --timeout=10s --start-period=30s \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:7860/health', timeout=5)"
CMD ["uvicorn", "kenryu_engine:app", "--host", "0.0.0.0", "--port", "7860"]
