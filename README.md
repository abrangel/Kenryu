# Kenryu - Motor Científico para Análisis Bioinformático

**Desarrollado por: Cesar Manzo**

Kenryu es una plataforma bioinformática avanzada diseñada para la identificación y análisis de biomarcadores core a partir de paneles de microRNAs (miRNAs). El sistema integra bases de datos genómicas de alta confianza (TargetScanHuman, miRTarBase) y evidencia científica en tiempo real de PubMed para proporcionar reportes diagnósticos de precisión.

## 🚀 Características Principales

*   **Identificación de Genes Core:** Algoritmos de convergencia molecular mediante intersección de dianas (Consenso Estricto, N-1, N-2).
*   **Contexto Biológico Real:** Integración con **MyGene.info** para obtener descripciones patológicas, funciones proteicas y rutas metabólicas (KEGG, Reactome).
*   **Evidencia Científica Automática:** Búsqueda y vinculación de artículos en **PubMed** relacionados con los genes identificados.
*   **Editor de Informe Pro:** Interfaz interactiva para la creación de reportes académicos con soporte para exportación a PDF, Markdown y texto plano.
*   **Visualización Científica:** Generación dinámica de Diagramas de Venn y gráficos de convergencia molecular.

## 🛠️ Arquitectura del Sistema

El proyecto sigue una arquitectura desacoplada para garantizar escalabilidad y rendimiento:

*   **Backend (Nube):** API construida con **FastAPI (Python)**, desplegada en **Hugging Face Spaces**. Incluye un sistema de caché binaria `.pkl` para procesamiento instantáneo de millones de registros.
*   **Frontend:** Interfaz moderna en **HTML5, CSS3 y JavaScript vanilla**, desplegada en **GitHub Pages**.
*   **Datos:** Base de datos biográfica comprimida y optimizada para despliegues ligeros.

## 📦 Instalación y Ejecución Local

### Backend
1. Navega a `demo/backend`.
2. Instala las dependencias:
   ```bash
   pip install -r requirements.txt
   ```
3. Ejecuta el servidor:
   ```bash
   python main.py
   ```
   *El servidor iniciará en `http://localhost:8000`.*

### Frontend
Simplemente abre el archivo `demo/frontend/index.html` en cualquier navegador moderno.

## 🌐 Despliegue en la Nube

### Backend (Hugging Face)
El código está optimizado para funcionar en Hugging Face Spaces usando Docker.
1. Crea un Space de tipo **Docker**.
2. Sube los archivos `main.py`, `Dockerfile`, `requirements.txt` y `targetscan_full.json.zip`.
3. El puerto de escucha es el **7860**.

### Frontend (GitHub Pages)
1. Sube la carpeta `demo/frontend` a un repositorio de GitHub.
2. Activa **GitHub Pages** desde la configuración del repositorio.
3. Asegúrate de actualizar la constante `API_URL` en `script.js` con la URL de tu Space de Hugging Face.

## 📄 Licencia
Este proyecto fue desarrollado bajo rigor académico para fines de investigación bioinformática y medicina de precisión.

---
*Bioinformatics Engine v1.0 — Cesar Manzo*
