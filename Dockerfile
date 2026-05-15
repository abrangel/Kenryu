FROM python:3.11-slim
WORKDIR /app
RUN apt-get update && apt-get install -y build-essential && rm -rf /var/lib/apt/lists/*
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY kenryu_engine.py .
COPY targetscan_full.json.zip .
COPY static/ ./static/
RUN mkdir -p local_db && chmod 777 local_db
EXPOSE 7860
CMD ["uvicorn", "kenryu_engine:app", "--host", "0.0.0.0", "--port", "7860"]
