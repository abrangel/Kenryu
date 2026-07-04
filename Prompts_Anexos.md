# Prompts Anexos — TFM KENRYU

Este documento contiene la recopilación exhaustiva e ingeniería de *prompts* utilizada para el desarrollo, diseño y generación de los recursos visuales e infografías científicas incluidas en el Trabajo Fin de Máster titulado **"Plataforma integrada para la automatización de flujos de trabajo transcriptómicos"**.

La separación de estos contenidos del manuscrito principal responde a las directrices académicas de optimización de texto y reproducibilidad técnica independiente.

---

## Anexo A. Prompt para la Imagen del Dogma Central (Figura 1)

**Modelo de IA:** GPT-2

### Prompt Estructurado
A highly detailed, professional medical and biological infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra-sharp vector art style, precise corporate color-coding. The layout is divided into two main sections: a left vertical side-panel and a main right grid with 4 distinct columns separated by clean, subtle vertical lines.

*   **LEFT PANEL: "EL DOGMA CENTRAL"**
    A vertical flowchart tracking molecular biology's central dogma.
    *   Top: Text "ADN" with a blue double helix icon.
    *   Down arrow with the label "Transcripción".
    *   Middle: Text "ARN" with a green single-stranded wavy line icon.
    *   Down arrow with the label "Traducción".
    *   Bottom: Text "Proteína" with a purple 3D globular protein complex icon.

*   **MAIN GRID: "TIPOS DE ARN"** (4 Columns from left to right with white background cards)
    *   **Column 1 (Green theme):** Title header "ARN mensajero (ARNm)"
        *   Graphic: A single-stranded mRNA wavy line showing nucleotide bases, with labels "5'" and "3'" at the ends.
        *   Subtitle header: "FUNCIÓN"
        *   Body text: "Lleva la información genética del ADN al ribosoma para ser traducida en proteína."
        *   Subtitle header: "CARACTERÍSTICAS"
        *   Bullet points text: "- Lineal", "- Codificante", "- Contiene codones que especifican aminoácidos"
    *   **Column 2 (Blue theme):** Title header "ARN ribosomal (ARNr)"
        *   Graphic: Detailed ribosome molecular structure showing a large subunit labeled "Subunidad mayor" and a small subunit labeled "Subunidad menor" interacting with an RNA strand.
        *   Subtitle header: "FUNCIÓN"
        *   Body text: "Es el componente estructural y catalítico de los ribosomas, donde ocurre la síntesis de proteínas."
        *   Subtitle header: "CARACTERÍSTICAS"
        *   Bullet points text: "- Muy abundante", "- Forma parte del ribosoma", "- Tiene actividad catalítica (peptidil transferasa)"
    *   **Column 3 (Purple theme):** Title header "ARN de transferencia (ARNt)"
        *   Graphic: Cloverleaf secondary structure of tRNA, showing the 3' end labeled "3'" and "CCA", and the bottom loop labeled "Anticodón".
        *   Subtitle header: "FUNCIÓN"
        *   Body text: "Transporta aminoácidos al ribosoma and reconoce los codones del ARNm mediante su anticodón."
        *   Subtitle header: "CARACTERÍSTICAS"
        *   Bullet points text: "- Estructura en forma de trébol", "- Posee anticodón", "- Especifica para cada aminoácido"
    *   **Column 4 (Orange theme):** Title header "ARN no codificantes (ARNnc)"
        *   Graphic: Four small diverse RNA stem-loop structures labeled individually as "miARN", "lncARN", "snARN", "piARN", with a small text below saying "y otros...".
        *   Subtitle header: "FUNCIÓN"
        *   Body text: "Regulan la expresión génica a diferentes niveles (epigenético, transcripcional, postranscripcional y de estabilidad)."
        *   Subtitle header: "CARACTERÍSTICAS"
        *   Bullet points text: "- Diversos tamaños y estructuras", "- No codifican proteínas", "- Reguladores clave de procesos celulares"

**Aesthetic Specs:** Flat design, minimalist boxes, matching color fonts for headers, perfect alignment, crisp sans-serif typography, high contrast, 8k resolution, textbook layout, perfect graphic design symmetry. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted ribosomes, missing columns, asymmetric boxes, dark background, borders around the canvas, signatures, watermarks, bad anatomy of molecules.

---

## Anexo B. Prompt para Tipos de Datos Ómicos y sus Tecnologías (Figura 2)

**Modelo de IA:** Nano Banana Pro

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra sharp vector art style, precise corporate color coding. The layout is divided into two main sections: a left vertical side panel and a main right grid with 6 distinct columns separated by clean, subtle vertical lines.

*   **LEFT PANEL: "TECNOLOGÍAS ÓMICAS"**
    A vertical conceptual flowchart showing the hierarchy of omics technologies.
    *   Top: Title "ÓMICAS" with a multicolor molecular network icon.
    *   Down arrow labeled "Niveles biológicos".
    *   Middle: Icons representing DNA, RNA, proteins, metabolites, epigenetic marks, and microbiota.
    *   Down arrow labeled "Tecnologías de alto rendimiento".
    *   Bottom: A stylized sequencing machine + mass spectrometer + chromatograph icons.
    *   Color theme: grayscale + accent colors for each omic layer.

*   **MAIN GRID: "TIPOS DE DATOS ÓMICOS Y SUS TECNOLOGÍAS"**
    Six columns, each with a white card, flat design, color coded header, crisp sans serif typography.
    *   **Column 1 — Genómica (ADN) — Blue theme:** Graphic: DNA double helix + sequencing machine icon. Subtitle: "Tecnologías" Bullet points: Secuenciación de Nueva Generación (NGS), Whole Genome Sequencing (WGS), Whole Exome Sequencing (WES), SNP arrays. Subtitle: "Descripción" Body text: "Analiza la secuencia completa o parcial del ADN para identificar variantes genéticas."
    *   **Column 2 — Transcriptómica (ARN) — Green theme:** Graphic: Single stranded RNA ribbon + expression heatmap. Tecnologías: RNA Seq, scRNA Seq (single cell RNA sequencing), Microarrays de expresión génica. Descripción: "Cuantifica y caracteriza los niveles de ARN para estudiar la expresión génica."
    *   **Column 3 — Proteómica (Proteínas) — Purple theme:** Graphic: 3D protein complex + mass spectrometry peaks. Tecnologías: Espectrometría de masas (MS/MS), LC MS, MALDI TOF, Protein microarrays. Descripción: "Analiza la estructura, abundancia y modificaciones de las proteínas."
    *   **Column 4 — Metabolómica (Metabolitos) — Orange theme:** Graphic: Small metabolite molecules + chromatogram peaks. Tecnologías: LC MS, GC MS, RMN (Resonancia Magnética Nuclear). Descripción: "Detecta y cuantifica metabolitos para estudiar rutas bioquímicas."
    *   **Column 5 — Epigenómica — Red theme:** Graphic: DNA with methylation marks + chromatin structure. Tecnologías: Bisulfite sequencing, ChIP Seq, ATAC Seq. Descripción: "Estudia modificaciones químicas del ADN y la cromatina que regulan la expresión génica."
    *   **Column 6 — Microbiómica — Teal theme:** Graphic: Microbial community icons + 16S rRNA phylogenetic tree. Tecnologías: 16S rRNA sequencing, Metagenómica Shotgun. Descripción: "Caracteriza comunidades microbianas y su composición funcional."

**Aesthetic Specs:** Flat design, minimalist boxes, matching color fonts for headers, perfect alignment, crisp sans serif typography, high contrast, 8k resolution, textbook layout, perfect graphic design symmetry. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted molecules, missing columns, asymmetric boxes, dark background, borders around the canvas, signatures, watermarks.

---

## Anexo C. Prompt para Algoritmos en la Predicción de Dianas Génicas (Figura 3)

**Modelo de IA:** Nano Banana Pro

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra sharp vector art style, precise corporate color coding. The layout is divided into a left vertical side panel and a main right grid with 3 large columns separated by subtle vertical lines.

*   **LEFT PANEL: "ALGORITMOS EN LA PREDICCIÓN DE DIANAS GÉNICAS"**
    Vertical conceptual flowchart showing the evolution of prediction algorithms.
    *   Top icon: microRNA hairpin + mRNA strand.
    *   Arrow downward labeled "Evolución de los métodos".
    *   Middle icon: phylogenetic tree (conservación evolutiva).
    *   Arrow downward labeled "Integración de características".
    *   Bottom icon: machine learning neural network + ΔG thermodynamic symbol.
    *   Color theme: grayscale + blue accents.

*   **MAIN GRID: 3 COLUMN LAYOUT** (white cards, flat design)
    *   **Column 1 — Primera generación: Complementariedad de secuencia (Green theme):** Header: "Algoritmos de primera generación" Subheader: "Complementariedad de secuencia". Graphic: Single stranded miRNA aligning to a 3’UTR region, highlighting the seed region (nt 2–7) in bright green. Text blocks: Criterios: Coincidencia exacta en la región semilla (nt 2–7), Sitios complementarios en la 3’UTR. Limitaciones: Alta tasa de falsos positivos, No considera estructura ni evolución. Ejemplos (icons): miRanda (versión inicial), RNAhybrid (primeras versiones).
    *   **Column 2 — Segunda generación: Conservación evolutiva (Blue theme):** Header: "Algoritmos de segunda generación" Subheader: "Conservación evolutiva". Graphic: Phylogenetic tree with conserved binding sites highlighted across species (human, mouse, chicken, zebrafish). Criterios: Comparación entre especies, Sitios conservados en mamíferos y vertebrados, Penalización de sitios no conservados. Razonamiento: “Un sitio conservado es más probable que sea funcional.” Ejemplos: TargetScan (primeras versiones), PicTar.
    *   **Column 3 — Tercera generación: Modelos integrativos avanzados (Purple theme):** Header: "Algoritmos de tercera generación" Subheader: "Modelos integrativos avanzados". Graphic: Multilayer diagram combining: seed match, ΔG thermodynamic stability, RNA secondary structure, machine learning neural network, Context++ score bar chart. Criterios utilizados: Región semilla (6mer, 7mer, 8mer), Conservación evolutiva, Estabilidad termodinámica (ΔG), Accesibilidad estructural del ARNm, Context++ score (AU-richness, posición, distancia al stop codon), Modelos estadísticos y machine learning. Ejemplos: TargetScan 7/8 (Context++), miRanda mirSVR, DIANA microT CDS, PITA (accesibilidad estructural).

**Aesthetic Specs:** Flat design, minimalist boxes, matching color fonts for headers, perfect alignment, crisp sans serif typography, high contrast, 8k resolution, textbook layout, perfect graphic design symmetry. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted molecules, missing columns, asymmetric boxes, dark background, borders around the canvas, signatures, watermarks.

---

## Anexo D. Prompt para Tipos de Evidencia Experimental en miRTarBase (Figura 4)

**Modelo de IA:** Nano Banana Pro

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra‑sharp vector art style, precise corporate color‑coding. The layout is divided into a left vertical side‑panel and a main right grid with 2 large columns separated by subtle vertical lines.

*   **LEFT PANEL: "EVIDENCIA EXPERIMENTAL EN miRTarBase"**
    Vertical conceptual diagram showing the hierarchy of evidence.
    *   Top icon: microRNA hairpin + mRNA strand.
    *   Arrow downward labeled "Validación experimental".
    *   Middle icon: test tube + fluorescence reporter.
    *   Arrow downward labeled "Clasificación por nivel de evidencia".
    *   Bottom icon: database symbol (stacked cylinders) labeled "miRTarBase".
    *   Color theme: grayscale + teal accents.

*   **MAIN GRID: 2 COLUMN LAYOUT** (white cards, flat design)
    *   **Column 1 — Evidencia fuerte (Green theme):** Header: "Evidencia fuerte" Subheader: "Interacciones validadas experimentalmente". Graphic: Icons representing each technique: luciferase reporter assay (glowing reporter gene), Western blot bands, qPCR amplification curve, DNA mutation symbol (mutagénesis dirigida), CLIP‑Seq schematic (crosslink + immunoprecipitation + sequencing). Text blocks - Incluye: Ensayos de reportero luciferasa, Western blot, qPCR, Mutagénesis dirigida, CLIP‑Seq. Description: “Evidencia directa y robusta que confirma la interacción miARN–gen.”
    *   **Column 2 — Evidencia funcional moderada (Blue theme):** Header: "Evidencia funcional moderada" Subheader: "Soporte indirecto o correlacional". Graphic: Icons for: microarray heatmap, RNA‑Seq read alignment, correlation graph miRNA–gene expression. Incluye: Microarrays, RNA‑Seq, Estudios de correlación expresión miARN–gen. Description: “Evidencia funcional que sugiere interacción, pero no la demuestra directamente.”

**Aesthetic Specs:** Flat design, minimalist boxes, matching color fonts for headers, perfect alignment, crisp sans‑serif typography, high contrast, 8k resolution, textbook layout, perfect graphic design symmetry. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted molecular icons, missing columns, asymmetric boxes, dark background, borders around the canvas, signatures, watermarks.

---

## Anexo E. Prompt para Principios FAIR en Bioinformática (Figura 5)

**Modelo de IA:** Nano Banana Pro

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra sharp vector art style, precise corporate color coding. The layout is divided into a left vertical side panel and a main right grid with 4 horizontal sections, each representing one FAIR principle.

*   **LEFT PANEL: "PRINCIPIOS FAIR EN BIOINFORMÁTICA"**
    Vertical conceptual flowchart illustrating how FAIR improves data quality and integration.
    *   Top icon: database symbol + magnifying glass (Findable).
    *   Arrow downward labeled "Estandarización".
    *   Middle icon: cloud download symbol (Accessible).
    *   Arrow downward labeled "Interoperabilidad".
    *   Lower icon: puzzle piece connectors (Interoperable).
    *   Arrow downward labeled "Reutilización".
    *   Bottom icon: circular arrows + document symbol (Reusable).
    *   Color theme: grayscale + blue accents.

*   **MAIN GRID: 4 SECTION LAYOUT** (white cards, flat design)
    *   **Section 1 — F: Findable (Green theme):** Header: "Findable" Subheader: "Fáciles de encontrar". Graphic: Database icon with unique identifiers (DOI, ENSG ID) highlighted in bright green. Text blocks: Identificadores únicos y estables (DOI, Ensembl, NCBI), Metadatos claros que describen el contenido, Repositorios con motores de búsqueda eficientes.
    *   **Section 2 — A: Accessible (Blue theme):** Header: "Accessible" Subheader: "Fáciles de acceder". Graphic: Cloud storage icon with open format files (FASTA, GFF, CSV). Text blocks: Acceso mediante API o descarga directa, Formatos abiertos y documentados, Disponibilidad incluso si la plataforma cambia.
    *   **Section 3 — I: Interoperable (Purple theme):** Header: "Interoperable" Subheader: "Compatibles entre sí". Graphic: Diagram showing multiple databases (TargetScan, miRTarBase, Ensembl, NCBI Gene) connected by standardized identifiers. Text blocks: Compatibilidad entre versiones del genoma, Uso de vocabularios controlados, Integración sin conflictos entre plataformas.
    *   **Section 4 — R: Reusable (Orange theme):** Header: "Reusable" Subheader: "Reutilizables sin perder calidad". Graphic: Circular arrows around a document with metadata tags. Text blocks: Metadatos completos (versión del genoma, parámetros, métodos), Estándares internacionales (FAIR), Datos listos para reproducir, comparar y ampliar.

**Aesthetic Specs:** Flat design, minimalist boxes, matching color fonts for headers, perfect alignment, crisp sans serif typography, high contrast, 8k resolution, textbook layout, perfect graphic design symmetry. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted molecules, missing columns, asymmetric boxes, dark background, borders around the canvas, signatures, watermarks.

---

## Anexo F. Prompt para Relación entre Ontologías, Rutas y Análisis Funcional (Figura 6)

**Modelo de IA:** GPT-2

### Prompt Estructurado
A highly detailed professional biomedical infographic designed for a master's thesis or scientific textbook. Ultra-clean white background, publication-quality vector graphics, sharp typography, perfect symmetry, scientific educational style, horizontal 16:9 layout. The infographic explains the conceptual progression from genes to biological interpretation through ontology databases, pathway databases, and enrichment analysis.

*   **CENTRAL VISUAL STRUCTURE:** A large horizontal scientific workflow running from left to right across the entire canvas.
*   **STEP 1 — GENES (Gray):** Large DNA helix icon. Gene identifiers surrounding the DNA: ENSG000001, ENSG000002, TP53, EGFR, BRCA1. Header: "Genes". Subheader: "Punto de partida del análisis biológico". Small text: Lista de genes obtenidos mediante análisis genómicos, transcriptómicos o proteómicos. Arrow pointing right.
*   **STEP 2 — GENE ONTOLOGY (Green):** Large hierarchical DAG diagram. Three interconnected branches: BP (Biological Process), MF (Molecular Function), CC (Cellular Component). Header: "Gene Ontology (GO)". Subtitle: "¿Qué hace el gen?". Visual emphasis on hierarchical relationships and controlled vocabulary. Callout: Anotación olitológica funcional. Small explanatory text: Describe funciones, procesos biológicos y localizaciones celulares mediante una ontología estructurada. Arrow pointing right.
*   **STEP 3 — PATHWAY DATABASES (Blue):** Large section divided into three equal panels:
    *   *KEGG:* Metabolic pathway map. Enzymes, metabolites and signaling modules. Label: "KEGG". Subtitle: "Modelado de rutas biológicas". Question: "¿En qué ruta participa?".
    *   *Reactome:* Detailed molecular interaction pathway. Protein complexes, binding events, phosphorylation, signal transduction. Label: "Reactome". Subtitle: "Representación de eventos moleculares". Question: "¿Cómo ocurre el mecanismo?".
    *   *WikiPathways:* Collaborative pathway network. Community-curated pathway diagrams. Multiple contributors connected to pathway maps. Label: "WikiPathways". Subtitle: "Curación colaborativa de vías biológicas". Question: "¿Qué conocimiento adicional existe?".
    *   Arrow pointing right.
*   **STEP 4 — ENRICHMENT ANALYSIS (Purple):** Large enrichment analysis panel. Barplot, dotplot, adjusted p-values, FDR correction symbols, highlighted enriched pathways. Header: "Análisis de enriquecimiento funcional". Subtitle: "Integración computacional de conocimiento biológico". Question: "¿Qué procesos están sobrerrepresentados?". Small text: Integra información procedente de GO, KEGG, Reactome y WikiPathways para identificar mecanismos biológicos relevantes. Arrow pointing right.
*   **STEP 5 — BIOLOGICAL INTERPRETATION (Orange):** Central biological system illustration. Human silhouette, cell, signaling pathways highlighted, disease and phenotype icons. Header: "Interpretación biológica". Subtitle: "Comprensión de mecanismos celulares". Question: "¿Qué significado biológico tienen los resultados?". Small text: Transforma datos ómicos en conocimiento biológico interpretable.
*   **BOTTOM SUMMARY STRIP:** A horizontal comparison table: GO → ¿Qué hace el gen?; KEGG → ¿En qué ruta participa?; Reactome → ¿Cómo ocurre el proceso?; WikiPathways → ¿Qué conocimiento adicional existe?; Enrichment Analysis → ¿Qué mecanismos están sobrerrepresentados?

**Design Specifications:** Scientific textbook figure, Nature Reviews style, Cell Press style, publication-quality biomedical infographic, vector graphics only, clean white background, flat design, perfect alignment, balanced spacing, professional color coding (green for GO, blue for pathway databases, purple for enrichment analysis, orange for biological interpretation), high readability, large typography, 8K resolution, master thesis quality, academic publication style, landscape composition, 16:9 aspect ratio.

**Negative Prompt:** photorealistic, 3D render, gradients, shadows, glossy effects, dark background, handwritten text, low resolution, blurry text, overlapping labels, distorted diagrams, decorative elements, artistic style, watercolor, sketch, random molecules, cluttered layout, asymmetrical composition, watermark, logo, signatures, border frame, gibberish text, spelling errors.

---

## Anexo G. Prompt para Gene Ontology (GO): Conocimiento, Estructura e Interacción (Figura 7)

**Modelo de IA:** GPT-2

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra sharp vector art style, precise corporate color coding. The layout is divided into a left vertical side panel and a main right grid with 3 horizontal sections, each representing a core dimension of Gene Ontology: conocimiento, estructura e interacción.

*   **LEFT PANEL: “GENE ONTOLOGY (GO)”**
    Vertical conceptual flowchart illustrating GO as a computable knowledge system.
    *   Top icon: open knowledge graph + DNA helix (ontología biológica).
    *   Arrow downward labeled "Vocabulario controlado".
    *   Middle icon: DAG hierarchy tree (estructura formal).
    *   Arrow downward labeled "Relaciones lógicas".
    *   Bottom icon: interconnected nodes (interoperabilidad entre bases de datos).
    *   Color theme: grayscale + blue accents.

*   **MAIN GRID: 3 SECTION LAYOUT** (white cards, flat design)
    *   **Section 1 — Conocimiento (Green theme):** Header: “Conocimiento biológico estandarizado” Subheader: “Qué representa GO”. Graphic: Knowledge graph with labeled nodes (BP, MF, CC) connected to gene IDs. Text blocks: GO es un vocabulario controlado que describe funciones génicas de forma uniforme. Permite anotar genes con términos estandarizados y comparables entre organismos. Facilita la interpretación funcional en análisis ómicos y bioinformáticos.
    *   **Section 2 — Estructura (Blue theme):** Header: “Estructura jerárquica (DAG)” Subheader: “Cómo se organiza la ontología”. Graphic: Directed acyclic graph (DAG) with parent–child relationships, showing increasing specificity. Text blocks: GO está organizado como un grafo dirigido acíclico. Cada término puede tener múltiples padres e hijos. Las relaciones incluyen: “is a”, “part of”, “regulates”, “positively regulates”, “negatively regulates”. La jerarquía permite niveles de detalle desde general hasta altamente específico.
    *   **Section 3 — Interacción (Purple theme):** Header: “Interacción entre categorías GO” Subheader: “BP, MF y CC como dimensiones complementarias”. Graphic: Three layer diagram showing BP, MF, CC interacting around a gene node. Text blocks: BP (Biological Process): describe procesos celulares y fisiológicos. MF (Molecular Function): describe actividades bioquímicas específicas. CC (Cellular Component): describe la localización subcelular. Las tres dimensiones interactúan para describir qué hace un gen, cómo lo hace y dónde lo hace. Esta integración permite inferencia biológica computacional.

**Aesthetic Specs:** Flat design, minimalist boxes, matching color fonts for headers, perfect alignment, crisp sans serif typography, high contrast, 8k resolution, textbook layout, perfect graphic design symmetry. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted molecules, missing columns, asymmetric boxes, dark background, borders around the canvas, signatures, watermarks.

---

## Anexo H. Prompt para Redes PPI e Identificación de Hubs Terapéuticos (Figura 8)

**Modelo de IA:** GPT-2

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra‑sharp vector art style, precise corporate color coding. The layout is divided into a top conceptual header and a main vertical pipeline with six sections illustrating the construction and análisis de redes proteína‑proteína (PPI) y la identificación de hubs como posibles blancos terapéuticos.

*   **SECTION 1 — Moléculas biológicas (Proteínas de interés):** Graphic: protein icons labeled A–H, derived from omics studies (genomics, transcriptomics, proteomics). Text: “Proteínas codificadas por genes diferencialmente expresados.”
*   **SECTION 2 — Red proteína‑proteína (PPI):** Graphic: network graph with nodes (proteins) and edges (interactions). Text: “Interacciones físicas o funcionales representadas como un grafo.”
*   **SECTION 3 — Análisis de centralidad:** Graphic: three metrics displayed with icons: Degree (número de conexiones), Betweenness (control de rutas cortas), Closeness (proximidad global en la red). Text: “Medidas para identificar nodos clave.”
*   **SECTION 4 — Identificación de hubs:** Graphic: nodes highlighted in red/orange showing high centrality; peripheral nodes in gray. Text: “Hubs: proteínas altamente conectadas o estratégicas.”
*   **SECTION 5 — Funciones biológicas de los hubs:** Graphic: modules and pathways connected through hub nodes. Text: “Coordinan procesos celulares, integran módulos y participan en rutas esenciales.”
*   **SECTION 6 — Blancos terapéuticos potenciales:** Graphic: drug‑target icon pointing to hub nodes. Text: “Los hubs pueden ser dianas críticas para intervención terapéutica.”

**Aesthetic Specs:** flat design, crisp sans‑serif typography, perfect alignment, high contrast, 8k resolution, textbook layout, symmetric composition, no gradients, no shadows. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted molecules, missing sections, dark background, borders, watermarks.

---

## Anexo I. Prompt para el Diagrama de Venn (Figura 9)

**Modelo de IA:** GPT-2

### Prompt Estructurado
A highly detailed, professional biomedical infographic canvas layout, formatted for a scientific textbook. Clean solid white background, ultra‑sharp vector art style, precise academic color coding. The layout presents a two‑set Venn diagram comparing gene predictions from TargetScan and miRTarBase.

*   **TOP HEADER — Conceptual Overview:** Graphic: minimalist title banner “Diagrama de Venn Bioinformático”. Text: “Intersección entre predicciones teóricas y evidencias experimentales para identificar genes consenso.”
*   **SECTION 1 — TargetScan (Predicciones teóricas):** Graphic: left circle in blue tones labeled “TargetScan (Predicciones teóricas)”. Inside label: “2,342 genes predichos únicamente por TargetScan”. Text: “Predicciones computacionales basadas en complementariedad de secuencia.”
*   **SECTION 2 — miRTarBase (Evidencias experimentales):** Graphic: right circle in green tones labeled “miRTarBase (Evidencias experimentales)”. Inside label: “1,126 genes validados únicamente por miRTarBase”. Text: “Interacciones miRNA‑gen confirmadas experimentalmente.”
*   **SECTION 3 — Intersección (Genes consenso):** Graphic: overlapping region in teal. Inside label: “1,289 Genes consenso”. Text: “Genes respaldados por predicción y validación experimental.”
*   **SECTION 4 — Interpretación funcional:** Graphic: DNA icon below the Venn diagram. Text: “Genes consenso: dianas con mayor respaldo biológico y relevancia funcional.”
*   **SECTION 5 — Legend:** Graphic: color legend on the right side: Blue = TargetScan only, Green = miRTarBase only, Teal = Consensus genes. Text: “Codificación cromática para interpretación rápida.”

**Aesthetic Specs:** flat design, crisp sans‑serif typography, perfect alignment, high contrast, 8k resolution, textbook layout, symmetric composition, no gradients, no shadows. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted circles, incorrect labels, dark background, borders, watermarks.

---

## Anexo J. Prompt para el Volcano Plot (Figura 10)

**Modelo de IA:** GPT-2

### Prompt Estructurado
A highly detailed, professional biomedical data visualization canvas layout, formatted for a scientific textbook. Clean solid white background, ultra sharp vector art style, precise academic color coding. The layout presents a full volcano plot for functional enrichment analysis, showing statistical significance versus biological relevance.

*   **TOP HEADER — Conceptual Overview:** Graphic: minimalist title banner “Volcano Plot de Enriquecimiento”. Text: “Visualización de significancia estadística (–log10 p value) frente a magnitud del cambio (log2 Fold Change).”
*   **SECTION 1 — Ejes y estructura del gráfico:** Graphic: X axis labeled “log2(Fold Change)” with centered zero line; Y axis labeled “–log10(p value)”. Graphic: horizontal dashed threshold line marking “p value = 0.05”. Text: “El eje X representa la magnitud del cambio; el eje Y representa la significancia estadística.”
*   **SECTION 2 — Codificación de colores (tres grupos):** Graphic: Red points on the right side = “Aumentado / Enriquecido en condición B”. Blue points on the left side = “Disminuido / Enriquecido en condición A”. Gray points in the center = “No significativo”. Text: “Clasificación visual de rutas según relevancia y significancia.”
*   **SECTION 3 — Etiquetas de procesos biológicos (lado izquierdo):** Graphic: labeled points on the left cluster. Labels: “Respuesta inmune adaptativa”, “Señalización de citoquinas”, “Regulación del ciclo celular”. Text: “Procesos disminuidos o menos representados.”
*   **SECTION 4 — Etiquetas de procesos biológicos (lado derecho):** Graphic: labeled points on the right cluster. Labels: “Metabolismo de aminoácidos”, “Ruta de oxidación de ácidos grasos”, “Biogénesis mitocondrial”. Text: “Procesos aumentados o más representados.”
*   **SECTION 5 — Anotaciones laterales:** Graphic: left annotation bar “Enriquecido en condición A (Disminuido)”. Graphic: right annotation bar “Enriquecido en condición B (Aumentado)”. Text: “Indicadores laterales para interpretación rápida.”

**Aesthetic Specs:** Flat design, crisp sans serif typography, perfect alignment, high contrast, 8k resolution, textbook layout, symmetric composition, no gradients, no shadows, no artistic distortion. --ar 16:9

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted axes, incorrect labels, missing threshold line, dark background, borders, watermarks.

---

## Anexo K. Prompt para la Función de APIs en Bioinformática (Figura 11)

**Modelo de IA:** GPT-2

### Prompt Estructurado
Título del canvas: Función de las APIs en Bioinformática. Estilo: Infografía técnica, diseño limpio, fondo blanco, tipografía sans serif, colores azul/gris, estilo vectorial, 16:9, alta resolución.

*   **SECCIÓN SUPERIOR — Encabezado**
    *   Título grande: “Función de las APIs en Bioinformática”
    *   Subtítulo: Conectividad, automatización y acceso dinámico a datos biomédicos
*   **SECCIÓN IZQUIERDA — Concepto general**
    *   *Bloque 1: ¿Qué es una API?* Ícono: engrane + nube. Texto breve: Una API es un puente que permite que dos sistemas intercambien información de forma estandarizada, sin necesidad de mover archivos completos o realizar procesos manuales.
    *   *Bloque 2: ¿Por qué son esenciales en bioinformática?* Ícono: ADN + servidor. Texto breve: Permiten acceder a datos genómicos, clínicos o transcriptómicos en tiempo real, integrando múltiples fuentes de evidencia dentro de un mismo flujo de trabajo.
*   **SECCIÓN DERECHA — Funciones principales (tres columnas)**
    *   *Columna 1 — Integración de datos externos:* Ícono: base de datos conectada. Texto: Las APIs permiten consultar bases de datos públicas (NCBI, Ensembl, GDC, UniProt) sin descargar datasets completos.
    *   *Columna 2 — Automatización de procesos:* Ícono: robot + script. Texto: Facilitan pipelines que ejecutan consultas, análisis y validaciones sin intervención manual.
    *   *Columna 3 — Interoperabilidad entre herramientas:* Ícono: nodos conectados. Texto: Hacen posible que plataformas clínicas, herramientas bioinformáticas y sistemas hospitalarios “hablen el mismo idioma”.
*   **SECCIÓN INFERIOR — Ejemplos reales**
    *   Lista con íconos: FHIR Genomics API → integración con EHR; GDC API → acceso a datos de cáncer; Ensembl REST API → anotación genética; NCBI E-utilities → búsquedas programáticas.

**Especificaciones Estéticas:** Diseño plano (flat design), líneas delgadas, cajas alineadas, simetría visual. Colores: azul (#1E88E5), gris oscuro (#424242), acentos en cian. Nada de sombras, fotos, ni texturas. Sin marcas de agua, sin firmas, sin bordes gruesos.

**Negative Prompt:** Fotos, 3D, sombras, gradientes, fondos oscuros, texto borroso, íconos deformados, estilos infantiles, manuscritos, ruido visual, desalineación, baja resolución, marcas de agua, elementos realistas.

---

## Anexo L. Prompt para Uso de Recursos Bioinformáticos Externos (Figura 12)

**Modelo de IA:** Nano Banana Pro

### Prompt Estructurado
Imagen horizontal dividida en tres módulos conectados por flechas. Estilo minimalista, fondo blanco, líneas finas, iconos planos. Cada módulo incluye un icono grande, un título claro y un subtítulo breve en español. Las flechas deben mostrar el flujo de información entre los módulos.

*   **Módulo 1 — “Acceso mediante APIs”**
    *   Icono: nube + base de datos + flechas bidireccionales.
    *   Título visible: “Acceso mediante APIs”
    *   Subtítulo: “Consulta de datos en tiempo real”
    *   Flecha hacia el módulo 2.
*   **Módulo 2 — “Integración de recursos externos”**
    *   Icono central: ADN.
    *   Alrededor, cuatro iconos simples con nombres debajo: KEGG (círculo azul con símbolo molecular), PubMed (libro abierto estilizado), Enrichr (cuadro con barras simples), Harmonizome (red de puntos).
    *   Título visible: “Integración de recursos externos”
    *   Subtítulo: “Bases de conocimiento conectadas”
    *   Flecha hacia el módulo 3.
*   **Módulo 3 — “Interoperabilidad”**
    *   Icono: red de nodos conectados.
    *   Título visible: “Interoperabilidad”
    *   Subtítulo: “Sistemas que se comunican entre sí”

**Estilo General:** Colores suaves azul y gris, tipografía sans serif, texto corto y legible. Nada de bloques largos, nada de estilo infografía tradicional, sin sombras, sin 3D, sin logos reales. Diseño limpio, profesional y fácil de entender.

**Negative Prompt:** no infografías complejas, no texto largo, no párrafos, no 3D, no sombras, no colores fuertes, no caricaturas, no desorden visual, no logos oficiales.

---

## Anexo M. Prompt para Fases de Creación de la Plataforma KENRYU (Figura 13)

**Modelo de IA:** Nano Banana Pro

### Prompt Estructurado
Infografía científica horizontal, estilo profesional, fondo blanco, líneas delgadas, colores azul y gris. Diseño limpio, modular y altamente legible. La imagen muestra un diagrama de flujo dividido en seis fases numeradas, cada una con título, subtítulo, iconos y cajas de contenido. Flechas conectan cada fase en orden secuencial. Todo el texto está en español, con tipografía sans serif moderna.

*   **Fase 1 — Input y Preprocesamiento:** Caja con título: “Paso 1: Input y Parámetros”. Subtítulo: “Sanitización y estandarización de miRNAs”. Lista de ejemplos: hsa-miR-25-5p, hsa-miR-22-3p, hsa-miR-21. Parámetro visible: Ventana temporal PubMed 2015–2024. Iconos: documento, engrane, símbolo de RegEx. Flecha hacia la Fase 2. Texto breve: “Salida: Lista de miRNAs normalizados y validados”.
*   **Fase 2 — Recolección Multi‑Fuente:** Título: “Paso 2: Recolección Multi‑Fuente”. Iconos de herramientas: TargetScanHuman v8.0, miRTarBase, Harmonizome API. Flecha hacia la Fase 3. Texto breve: “Salida: Genes predichos y validados por miRNA”.
*   **Fase 3 — Convergencia Genómica (Core Engine):** Título: “Paso 3: Convergencia Genómica”. Subtítulo: “Modos de consenso: Strict, N‑1, N‑2”. Tabla simplificada con columnas: GENE, miR‑21, miR‑22, miR‑25, miR‑29, Conteo. Icono: engrane central. Flecha hacia Fase 4. Texto breve: “Salida: Lista de genes convergentes”.
*   **Fase 4 — Enriquecimiento Funcional y Visualización:** Título: “Paso 4: Enriquecimiento Funcional y Visualización”. Iconos de bases de datos: KEGG, Reactome, WikiPathways. Iconos de gráficos: diagrama de Venn, red STRING‑DB, volcano plot, bubble plot. Flecha hacia Fase 5. Texto breve: “Salida: Visualizaciones generadas”.
*   **Fase 5 — Investigación Multi‑Fuente (Evidencia Clínica):** Título: “Paso 5: Módulo de Evidencia”. Iconos de fuentes: PubMed, OMIM, ClinVar, ClinicalTrials.gov. Tabla con columnas: GENE, PubMed, OMIM, ClinVar, ClinicalTrials, Score. Flecha hacia Fase 6. Texto breve: “Salida: Evidencia integrada por gen”.
*   **Fase 6 — Reporte Académico:** Título: “Paso 6: Reporte Académico (Salida Final)”. Iconos: PDF, Markdown, documento con texto. Subtítulo: “Síntesis narrativa + Bibliografía estilo Vancouver”. Vista previa de reporte: “Reporte de Análisis Multi‑Fuente de Convergencia de miRNAs”. Texto breve: “Salida: Reporte listo para publicación”.

*   **Sección Inferior — Stack Tecnológico:** Título: “Stack Tecnológico Actualizado”. Lista con iconos:
    *   Backend: Python, FastAPI
    *   Frontend: HTML5, CSS3, JS
    *   Comunicación: API REST
    *   Motor: Procesamiento asíncrono
    *   Almacenamiento: Base de datos relacional
    *   Deployment: Docker, GitHub Pages, Hugging Face Spaces
*   **Beneficios Clave:** Análisis multi‑fuente robusto y reproducible, Evidencia clínica y genómica integrada, Reportes académicos listos para publicación.

**Negative Prompt:** photo, realistic, 3D render, shadows, gradients, messy layout, overlapping text, spelling mistakes, gibberish characters, blur, grainy texture, organic backgrounds, artistic chaos, handwritten fonts, low resolution, distorted graphs, misaligned text layers.
