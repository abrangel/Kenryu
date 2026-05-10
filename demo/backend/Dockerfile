# Usar una imagen base ligera de Python
FROM python:3.11-slim

# Establecer directorio de trabajo
WORKDIR /app

# Instalar dependencias del sistema necesarias
RUN apt-get update && apt-get install -y \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

# Copiar requerimientos e instalar dependencias de Python
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# COPIAR TODO EL CONTENIDO (Esto incluye el ZIP que subiste por la web)
COPY . .

# Crear directorio para la caché de datos
RUN mkdir -p local_db

# Hugging Face requiere el puerto 7860
EXPOSE 7860

# Comando para iniciar el servidor
CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "7860"]
