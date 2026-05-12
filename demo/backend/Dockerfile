# Usar una imagen base ligera de Python
FROM python:3.11-slim

# Establecer directorio de trabajo
WORKDIR /app

# Instalar dependencias del sistema necesarias para Matplotlib
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copiar requerimientos e instalar dependencias de Python
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copiar el código del backend y los datos core
COPY main.py .
COPY targetscan_full.json.zip .

# Crear directorio para la caché de datos
RUN mkdir -p local_db

# Exponer el puerto de FastAPI
EXPOSE 8000

# Comando para iniciar el servidor (optimizado para producción)
CMD ["uvicorn", "main.py:app", "--host", "0.0.0.0", "--port", "8000", "--workers", "4"]
