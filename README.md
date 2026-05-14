---
title: KENRYU Bioinformatica
emoji: 🧬
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 7860
pinned: false
---

# KENRYU: Motor de Análisis de Convergencia Genómica

KENRYU es un motor bioinformático avanzado diseñado para la identificación de dianas moleculares y el análisis de co-regulación mediada por microARNs. El sistema integra predicciones termodinámicas con validación experimental para ofrecer un panorama clínico de alta precisión sobre el silenciamiento génico post-transcripcional.

## Características Principales

*   **Identificación Multi-fuente:** Integración de TargetScanHuman v8.0 y miRTarBase para una cobertura total de dianas predichas y validadas experimentalmente.
*   **Análisis de Convergencia Molecular:** Algoritmos de intersección avanzada (Strict, N-1, N-2) para identificar genes regulados por múltiples microARNs simultáneamente.
*   **Enriquecimiento Funcional Balanceado:** Consultas simultáneas a KEGG, Reactome, WikiPathways y GO con algoritmos de equilibrio de fuentes para evitar sesgos estadísticos.
*   **Motor de Traducción Técnica:** Traducción automática en tiempo real de términos biológicos y resúmenes clínicos mediante la API MyMemory.
*   **Informe Académico Profesional:** Generación de reportes clínicos paginados con portadas independientes, perfiles de genes de élite y referencias bibliográficas automatizadas.
*   **Visualización Científica Dinámica:** Generación automática de diagramas de Venn, gráficos de convergencia molecular (tipo flor), Volcano Plots e Interactomas PPI (Redes STRING).

## Tecnologías Utilizadas

*   **Backend:** FastAPI (Python 3.9+)
*   **Bioinformática:** GSEApy, Pandas, NumPy, Biopython
*   **Frontend:** JavaScript (Vanilla), HTML5, CSS3
*   **Bases de Datos:** TargetScanHuman, miRTarBase, MyGene.info, PubMed (NCBI)

## Despliegue en Hugging Face Spaces

El proyecto está configurado para ejecutarse en entornos Docker. Para desplegar en Hugging Face:

1.  Crea un nuevo Space de tipo **Docker**.
2.  Sube todos los archivos de este repositorio.
3.  El motor utilizará automáticamente el puerto configurado en las variables de entorno o el puerto por defecto 7860.

## Seguridad

El sistema incluye una configuración de CORS (Cross-Origin Resource Sharing) endurecida para restringir el acceso a orígenes de confianza (Hugging Face y entornos locales), protegiendo la integridad del motor de análisis.

---
**Desarrollado por:** Cesar Manzo
**Especialidad:** Bioinformática y Genómica Clínica
