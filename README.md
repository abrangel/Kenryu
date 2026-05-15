---

# KENRYU — Plataforma de Análisis Bioinformático de microARNs

**Desarrollado por: Cesar Manzo**  
**Especialidad:** Análisis Genómico y Bioinformática Clínica | **Estado:** En producción

> Kenryu es una plataforma bioinformática avanzada diseñada para la identificación de dianas moleculares convergentes, el análisis de co-regulación mediada por microARNs (miRNAs) y la recuperación automatizada de evidencia científica en PubMed.

🔗 **Aplicación en línea:** [KENRYU en Hugging Face Spaces](https://huggingface.co/spaces/Kenryu007/Bioinformatica)
🔗 **Presentación Profesional:** [KENRYU en GitHub Pages](https://abrangel.github.io/Kenryu)

---

## ¿Qué es este proyecto?

Los microARNs (miRNAs) son reguladores maestros que controlan la expresión de redes génicas complejas. En la investigación genómica, identificar qué genes son regulados simultáneamente por múltiples miRNAs es crucial para descubrir nodos críticos en patologías como el cáncer o enfermedades neurodegenerativas.

Kenryu automatiza este flujo de trabajo: integra predicciones computacionales con validación experimental, aplica algoritmos de intersección avanzados, identifica rutas biológicas significativas (KEGG, Reactome, GO) y genera informes académicos con traducción automática al español y referencias científicas reales.

---

## Cómo funciona — Flujo de Análisis Bioinformático

### 1. Integración de Evidencia Híbrida
El sistema consolida datos de dos vertientes para maximizar la precisión:
*   **Predicción Termodinámica:** Utiliza el motor de TargetScanHuman v8.0 para identificar sitios de unión basados en afinidad de secuencia.
*   **Validación Experimental:** Consulta la base de datos miRTarBase (vía Harmonizome) para extraer interacciones confirmadas en laboratorio mediante CLIP-seq, Luciferasa o Western Blot.

### 2. Algoritmos de Convergencia Molecular
El investigador puede ajustar el rigor del consenso biológico:
*   **Modo Estricto (Strict):** Identifica genes regulados por la totalidad del panel de miRNAs.
*   **Modos N-1 / N-2:** Detecta genes que son dianas de la mayoría de los miRNAs, permitiendo descubrir redes regulatorias robustas incluso ante datos incompletos en bases de datos.

### 3. Enriquecimiento Funcional Balanceado
Consulta simultáneamente cuatro fuentes globales: **KEGG**, **Reactome**, **WikiPathways** y **GO (Biological Process)**. El sistema aplica un algoritmo de balanceo para asegurar una representación equitativa de todas las bases de datos en el informe final.

### 4. Motor de Traducción y Paginación
*   **Traducción Automática Real:** Traduce dinámicamente descripciones técnicas y rutas biológicas del inglés al español mediante la API MyMemory.
*   **Maquetación Académica:** Genera reportes en formato A4 con paginación inteligente por medición de altura, asegurando que el contenido, los gráficos y la bibliografía se distribuyan profesionalmente sin solapamientos.

---

## Arquitectura y APIs Integradas

| API / Base de Datos | Función en el sistema |
|---|---|
| **TargetScan v8.0** | Motor local de predicciones termodinámicas |
| **Harmonizome (miRTarBase)** | Validación experimental remota |
| **Enrichr (KEGG/Reactome)** | Enriquecimiento funcional dinámico |
| **MyGene.info** | Anotación clínica y patológica de genes |
| **NCBI PubMed** | Búsqueda activa de evidencia científica y PMIDs |
| **MyMemory API** | Traducción técnica en tiempo real |

---

## Estructura del Repositorio

```
Kenryu/
├── kenryu_engine.py         # Motor principal (FastAPI + Bioinformática)
├── requirements.txt         # Dependencias del entorno
├── Dockerfile               # Configuración de despliegue en contenedores
├── targetscan_full.json.zip # Base de datos TargetScan indexada
├── index.html               # Interfaz de usuario (Página principal)
├── script.js                # Lógica del frontend y paginación inteligente
├── style.css                # Diseño visual profesional
├── local_db/                # Caché persistente y bases de datos locales
└── README.md                # Documentación técnica
```

---

## Instalación y Ejecución

### Requisitos
*   Python 3.9 o superior
*   pip

### Instalación Local
1.  Clonar repositorio: `git clone https://github.com/abrangel/Kenryu.git`
2.  Instalar dependencias: `pip install -r requirements.txt`
3.  Iniciar motor: `python kenryu_engine.py`
4.  Acceso: `http://localhost:3001`

### Despliegue Híbrido
Este proyecto está optimizado para funcionar con el frontend en **GitHub Pages** y el procesamiento pesado (backend) en **Hugging Face Spaces**, garantizando una experiencia de usuario rápida y segura.

---

## Contexto Académico

Este proyecto ha sido desarrollado como una herramienta de grado clínico para el soporte de investigaciones en genómica y medicina traslacional. El enfoque modular y la integración de múltiples fuentes de evidencia permiten distinguir entre dianas predichas y aquellas con respaldo de literatura funcional publicada.

---
*KENRYU Bioinformatics Engine — Cesar Manzo*
