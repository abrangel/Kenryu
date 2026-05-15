---
title: KENRYU Bioinformatics
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 7860
pinned: false
---

# KENRYU: Plataforma Avanzada de Convergencia Genómica y Regulación por microARNs

KENRYU representa un avance significativo en la interpretación bioinformática de los mecanismos de silenciamiento génico post-transcripcional. Diseñado para investigadores y clínicos, este motor permite identificar con precisión quirúrgica cómo grupos específicos de microARNs (miRNAs) convergen sobre redes de genes diana, revelando nodos críticos de regulación en diversas patologías.

## Fundamentos Científicos

La plataforma no se limita a una simple búsqueda de bases de datos; aplica un flujo de análisis multidimensional:

### 1. Integración de Evidencia Híbrida
El sistema consolida datos de dos vertientes esenciales:
*   **Predicción Termodinámica:** Utiliza algoritmos de TargetScanHuman v8.0 para identificar sitios de unión basados en la afinidad de secuencia y conservación evolutiva.
*   **Validación Experimental:** Consulta en tiempo real la base de datos miRTarBase para extraer interacciones confirmadas mediante ensayos de Luciferasa, Western Blot o CLIP-seq.

### 2. Algoritmos de Intersección y Consenso
KENRYU permite ajustar el rigor del análisis mediante tres modos de convergencia:
*   **Modo Estricto:** Identifica únicamente los genes regulados por la totalidad del panel de miRNAs ingresado.
*   **Modos Flexibles (N-1 / N-2):** Detecta genes que son dianas de la mayoría de los miRNAs, permitiendo descubrir redes regulatorias más amplias pero estadísticamente robustas.

### 3. Análisis de Enriquecimiento Funcional Balanceado
A diferencia de herramientas convencionales que pueden presentar sesgos hacia procesos biológicos generales, KENRYU implementa un algoritmo de balanceo que extrae los hallazgos más significativos de cuatro fuentes globales: **KEGG**, **Reactome**, **WikiPathways** y **Gene Ontology (GO)**. Esto asegura que el investigador obtenga una visión equitativa y profunda de las rutas metabólicas y de señalización afectadas.

## Innovaciones Tecnológicas

*   **Motor de Traducción Clínica:** Integración con APIs de procesamiento de lenguaje natural para traducir automáticamente descripciones técnicas y rutas biológicas del inglés al español, facilitando la interpretación inmediata de resultados.
*   **Paginación Inteligente de Reportes:** El sistema genera informes académicos en formato A4 con una maquetación profesional. Utiliza un motor de medición de altura para distribuir el contenido, gráficos y bibliografía en hojas independientes, garantizando un documento listo para su publicación o presentación clínica.
*   **Visualización Dinámica:** Generación automática de diagramas de convergencia molecular (tipo pétalos para paneles complejos), paisajes de significancia biológica (Volcano Plots) y redes de interacción proteína-proteína (PPI).

## Guía de Instalación y Uso

### Requisitos Previos
*   Python 3.9 o superior
*   Docker (opcional, para despliegue en contenedores)

### Instalación Local
1.  Clonar el repositorio:
    ```bash
    git clone https://github.com/abrangel/Kenryu.git
    cd Kenryu
    ```
2.  Instalar dependencias:
    ```bash
    pip install -r requirements.txt
    ```
3.  Ejecutar el motor:
    ```bash
    python kenryu_engine.py
    ```
4.  Acceder a la interfaz: `http://localhost:3001`

### Despliegue en Hugging Face
El archivo `Dockerfile` y los metadatos en el `README.md` están preconfigurados para un despliegue inmediato en Hugging Face Spaces. Solo es necesario subir el código y el archivo de datos `targetscan_full.json.zip`.

## Seguridad y Arquitectura
El motor utiliza **FastAPI** para una respuesta de alta velocidad y cuenta con políticas de seguridad **CORS restrictivas**, protegiendo la integridad de la API y restringiendo su uso a dominios autorizados y entornos de desarrollo seguros.

---
**Desarrollado por:** Cesar Manzo
**Especialidad:** Análisis Genómico y Bioinformática Clínica
