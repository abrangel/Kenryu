# Usar una imagen base de Python con soporte para librerías científicas
FROM python:3.11-slim

# Establecer directorio de trabajo
WORKDIR /app

# Instalar dependencias del sistema
RUN apt-get update && apt-get install -y \
    build-essential \
    libpq-dev \
    && rm -rf /var/lib/apt/lists/*

# Copiar requerimientos e instalar dependencias de Python
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# --- ESTRUCTURA DE ARCHIVOS DE LA TESIS ---
# Copiar el backend
COPY main.py .

# Copiar los datos RAW (Dianas reales de los microRNAs de la tesis)
RUN mkdir -p repo_tesis/data/raw
COPY repo_tesis/data/raw/ repo_tesis/data/raw/

# Crear directorio para la caché de datos
RUN mkdir -p local_db

# Hugging Face requiere el puerto 7860
EXPOSE 7860

# Comando para iniciar el servidor
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "7860"]
