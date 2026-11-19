"""KENRYU core: pure, deterministic helpers (no network, no FastAPI).

Extracted from kenryu_engine.py so the logic can be unit-tested in
isolation. kenryu_engine.py re-exports these names, so existing imports
(`kenryu_engine.deep_clean_mirna`, ...) keep working.
"""
import base64
import io
import logging
import random
import re
from collections import Counter

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns
from matplotlib_venn import venn2, venn3

logger = logging.getLogger(__name__)


def deep_clean_mirna(m: str) -> str:
    m = m.replace('\u2011', '-').replace('\u2013', '-').replace('\u2014', '-').strip()
    # Strip species prefixes (hsa-, mmu-, rno-, ...), repeated if needed —
    # only when a miR/let/mir name follows, so "mir-33a-5p" is never touched.
    m = re.sub(r'^([a-z]{3}-)+(?=miR|let|mir)', '', m, flags=re.IGNORECASE)
    return m.lower()

def clean_mirna_name_cloud(m: str, species: str = "hsa") -> str:
    c = deep_clean_mirna(m)
    parts = c.split('-')
    # The let-7 family keeps its own name (hsa-let-7a-5p / mmu-let-7a-5p).
    # The generic rule below would drop "let" and wrongly return
    # hsa-miR-7a-5p, so the remote miRTarBase lookup queried a name
    # that does not exist.
    if "let" in parts[:2]:
        i = parts.index("let")
        rest = parts[i + 1:]
        r = f"{species}-let"
        if rest:
            r += "-" + "-".join(rest)
        return r
    if len(parts) > 1:
        num = parts[1]; suf = '-'.join(parts[2:])
        r = f"{species}-miR-{num}"
        if suf: r += f"-{suf}"
        return r
    return f"{species}-miR-{c}"

# ── CARGA DE DATOS ────────────────────────────────────────────────────────────
def consensus_genes(gene_sets: list, mode: str) -> set:
    """Convergent genes of a miRNA panel, by consensus mode (pure function).

    strict: intersection of all sets. n-1: genes in at least n-1 sets.
    n-2: genes in at least n-2 sets. Unknown modes behave like n-2.
    Extracted verbatim from analyze() so the logic is unit-testable.
    """
    if mode == "strict" or len(gene_sets) == 1:
        return set.intersection(*gene_sets)
    n = len(gene_sets)
    threshold = max(1, n - 1) if mode == "n-1" else max(1, n - 2)
    all_genes = [g for s in gene_sets for g in s]
    return {g for g, c in Counter(all_genes).items() if c >= threshold}


def create_visuals(gene_sets: list, mirnas: list, common_all: list,
                   pubmed_counts: dict | None = None,
                   ppi_edges: list | None = None):
    plt.rcParams.update(plt.rcParamsDefault)
    plt.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'text.color': '#1a1a2e', 'axes.labelcolor': '#1a1a2e',
        'xtick.color': '#444444', 'ytick.color': '#444444',
        'axes.edgecolor': '#cccccc', 'font.family': 'sans-serif',
        'axes.spines.top': False, 'axes.spines.right': False,
    })

    # ── VENN / CO-REGULACIÓN ─────────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(8, 8), facecolor='white')
    ax1.set_facecolor('white'); ax1.axis('off')
    n = len(gene_sets)
    PETAL_COLORS = ['#a8d8ea', '#aa96da', '#fcbad3', '#ffffd2', '#b5ead7',
                    '#ffdac1', '#e2f0cb', '#c7ceea', '#f8b195', '#f67280']
    try:
        if 2 <= n <= 3:
            if n == 2:
                v = venn2(gene_sets[:2], set_labels=mirnas[:2], ax=ax1)
            else:
                v = venn3(gene_sets[:3], set_labels=mirnas[:3], ax=ax1)
            for text in ax1.texts:
                text.set_color('#1a1a2e'); text.set_fontsize(10)
        else:
            from matplotlib.patches import Ellipse
            ax1.set_xlim(0, 1); ax1.set_ylim(0, 1)
            for i in range(n):
                angle_deg = (360 / n) * i
                angle_rad = np.radians(angle_deg)
                offset = 0.18
                cx = 0.5 + offset * np.cos(angle_rad)
                cy = 0.5 + offset * np.sin(angle_rad)
                color = PETAL_COLORS[i % len(PETAL_COLORS)]
                ell = Ellipse(xy=(cx, cy), width=0.38, height=0.55,
                              angle=angle_deg, facecolor=color, alpha=0.55,
                              edgecolor='white', linewidth=1.5, zorder=2)
                ax1.add_patch(ell)
                lx = 0.5 + 0.44 * np.cos(angle_rad)
                ly = 0.5 + 0.44 * np.sin(angle_rad)
                short = mirnas[i].replace('hsa-', '')
                ax1.text(lx, ly, short, ha='center', va='center',
                         color='#1a1a2e', fontsize=8.5, fontweight='bold', zorder=5,
                         bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                                   edgecolor='none', alpha=0.7))
            center_circ = patches.Circle((0.5, 0.5), 0.13, facecolor='white',
                                         edgecolor='#aaaaaa', linewidth=1.5, zorder=6)
            ax1.add_patch(center_circ)
            ax1.text(0.5, 0.52, "CO-REGULACIÓN", ha='center', va='center',
                     color='#1a1a2e', fontsize=7.5, fontweight='bold', zorder=7)
            ax1.text(0.5, 0.46, f"{len(common_all)}", ha='center', va='center',
                     color='#1a3a6b', fontsize=16, fontweight='bold', zorder=7)
            ax1.text(0.5, 0.40, "GENES CORE", ha='center', va='center',
                     color='#1a1a2e', fontsize=7.5, fontweight='bold', zorder=7)
    except Exception as e:
        logger.warning(f"Venn: {e}")
        ax1.text(0.5, 0.5, f"{len(common_all)} genes\ncomunes",
                 ha='center', va='center', transform=ax1.transAxes,
                 color='#1a1a2e', fontsize=16, fontweight='bold')
    ax1.set_title("Diagrama de Co-regulación miRNA",
                  color='#1a1a2e', fontsize=13, pad=14, fontweight='bold')
    buf1 = io.BytesIO()
    plt.savefig(buf1, format='png', facecolor='white', bbox_inches='tight', dpi=150)
    plt.close(fig1); buf1.seek(0)

    # ── VOLCANO (DATOS REALES) ───────────────────────────────────────────────
    # Adaptación propia de la tesis (Fig. 15): X = frecuencia en PubMed,
    # Y = nº de miRNAs del panel que regulan el gen. Cada punto es un gen
    # real de la unión (top 25 por nº de miRNAs); los genes core
    # (regulados por TODOS los miRNAs del panel) van en rojo y con nombre.
    # Sin conteos reales de PubMed (modo offline) el gráfico NO se dibuja:
    # no se inventa nada. (La versión anterior dibujaba 200 puntos
    # aleatorios con np.random: era una ilustración, no datos.)
    buf2 = None
    if pubmed_counts:
        union_counts = {}
        for s in gene_sets:
            for g in s:
                union_counts[g] = union_counts.get(g, 0) + 1
        n_mir = len(mirnas)
        ranked = sorted(union_counts,
                        key=lambda g: (-union_counts[g],
                                       -(pubmed_counts.get(g, 0)), g))[:25]
        xs = [np.log10(pubmed_counts.get(g, 0) + 1) for g in ranked]
        ys = [union_counts[g] for g in ranked]
        is_core = [union_counts[g] == n_mir for g in ranked]
        fig2, ax2 = plt.subplots(figsize=(10, 6), facecolor='white')
        ax2.set_facecolor('white')
        xs_bg = [x for x, c in zip(xs, is_core) if not c]
        ys_bg = [y for y, c in zip(ys, is_core) if not c]
        xs_co = [x for x, c in zip(xs, is_core) if c]
        ys_co = [y for y, c in zip(ys, is_core) if c]
        core_genes = [g for g, c in zip(ranked, is_core) if c]
        if xs_bg:
            ax2.scatter(xs_bg, ys_bg, c='#6baed6', s=30, alpha=0.7,
                        label='Genes de la unión')
        if xs_co:
            ax2.scatter(xs_co, ys_co, c='#e05c5c', s=70, zorder=5,
                        edgecolors='#8b1a1a', linewidth=0.6,
                        label='Genes core')
            for g, x, y in zip(core_genes, xs_co, ys_co):
                ax2.annotate(g, (x, y), textcoords="offset points",
                             xytext=(6, 4), fontsize=9, color='#1a3a6b',
                             fontweight='bold')
        ax2.set_xlabel("Artículos en PubMed (log10)", color='#444', fontsize=11)
        ax2.set_ylabel("miRNAs del panel que regulan el gen", color='#444', fontsize=11)
        ax2.set_title("Volcán: evidencia en PubMed vs. regulación por el panel",
                      color='#1a1a2e', fontsize=13, pad=14, fontweight='bold')
        ax2.set_ylim(0.5, n_mir + 0.8)
        ax2.set_xlim(left=0)
        ax2.legend(fontsize=9, framealpha=0.9, edgecolor='#cccccc')
        buf2 = io.BytesIO()
        plt.savefig(buf2, format='png', facecolor='white', bbox_inches='tight', dpi=150)
        plt.close(fig2); buf2.seek(0)

    # ── INTERACTOMA PPI (DATOS REALES: STRING-DB) ────────────────────────────
    # Cada arista es una interacción real de STRING-DB (confianza >= 0,4);
    # el grosor es la confianza. Sin aristas reales (modo offline) el
    # gráfico NO se dibuja. (La versión anterior dibujaba un anillo con
    # aristas inventadas por posición: no eran interacciones.)
    buf3 = None
    if ppi_edges is not None and common_all:
        # Todos los genes core son nodos (los aislados también se ven:
        # que STRING no registre interacción entre dos genes a >= 0,4
        # también es información). Las aristas son solo las reales, y los
        # nodos grises son proteínas vecinas de primera corona (las que
        # más interactúan con los genes core según STRING).
        core_nodes = list(dict.fromkeys(common_all))[:12]
        degree = {}
        for a, b, _s in ppi_edges:
            for g in (a, b):
                degree[g] = degree.get(g, 0) + 1
        neighbor_nodes = sorted(
            {g for a, b, _s in ppi_edges for g in (a, b)
             if g not in core_nodes},
            key=lambda g: (-degree[g], g))
        nodes = (sorted(core_nodes, key=lambda g: (-degree.get(g, 0), g))
                 + neighbor_nodes)[:16]
        edges = [(a, b, s) for a, b, s in ppi_edges
                 if a in nodes and b in nodes]
        fig3, ax3 = plt.subplots(figsize=(9, 9), facecolor='white')
        ax3.set_facecolor('white'); ax3.axis('off')
        n_nodes = len(nodes)
        r = 0.33
        pos = {}
        for i, gene in enumerate(nodes):
            angle = 2 * np.pi * i / max(n_nodes, 1) - np.pi / 2
            pos[gene] = (0.5 + r * np.cos(angle), 0.5 + r * np.sin(angle))
        for a, b, s in edges:
            x1, y1 = pos[a]
            x2, y2 = pos[b]
            ax3.plot([x1, x2], [y1, y2], c='#8aa8c8',
                     linewidth=0.6 + 2.4 * s, zorder=1, alpha=0.75)
        for gene in nodes:
            gx, gy = pos[gene]
            color = '#4e9af1' if gene in core_nodes else '#9aa5b1'
            ax3.scatter([gx], [gy], s=220 + 90 * degree.get(gene, 0), c=color,
                        zorder=5, edgecolors='white', linewidth=1.8)
            ax3.text(gx, gy + 0.06, gene, ha='center', va='center',
                     fontsize=9, fontweight='bold', color='#1a1a2e', zorder=7,
                     bbox=dict(boxstyle='round,pad=0.2', facecolor='white',
                               edgecolor='none', alpha=0.8))
        n_neighbors = len([g for g in nodes if g not in core_nodes])
        ax3.set_title(f"Interactoma PPI — STRING-DB (confianza ≥ 0,4), "
                      f"{len(edges)} interacciones"
                      + (f" — {n_neighbors} proteínas vecinas (grises)"
                         if n_neighbors else ""),
                      color='#1a1a2e', fontsize=13, pad=14, fontweight='bold')
        ax3.set_xlim(0.05, 0.95); ax3.set_ylim(0.05, 0.95)
        buf3 = io.BytesIO()
        plt.savefig(buf3, format='png', facecolor='white', bbox_inches='tight', dpi=150)
        plt.close(fig3); buf3.seek(0)

    def _b64(buf):
        return base64.b64encode(buf.getvalue()).decode() if buf is not None else None
    return (_b64(buf1), _b64(buf2), _b64(buf3))

# ── SÍNTESIS ACADÉMICA ───────────────────────────────────────────────────────��────────
def build_synthesis(mirnas: list, common: list, gene_details: dict, enrichment: list, references: list, gene_research: dict = None) -> dict:
    """Construye las 3 secciones narrativas del informe.

    gene_research: dict opcional con datos multi-fuente por gen
      {gene: {pubmed: [...], omim: [...], clinvar: {count, url}, trials: [...]}}
    Cuando está disponible, enriquece la narrativa de cada gen core con
    menciones contextualizadas a OMIM, ClinVar y ClinicalTrials.
    """
    gene_research = gene_research or {}
    n = len(common)
    core_str = ", ".join(common[:6]) + ("..." if n > 6 else "")
    
    # 1. Introducción
    p1 = (f"El presente estudio ha identificado un núcleo genómico de alta convergencia. A diferencia de un análisis "
          f"estadístico convencional, la integración de datos de PubMed y bases de datos biológicas revela una "
          f"interconexión funcional profunda entre los {n} genes identificados ({core_str}). "
          f"El análisis multiplataforma ejecutado sobre {len(mirnas)} microARN(s) ({', '.join(mirnas)}) "
          f"identificó estos biomarcadores mediante intersección de TargetScan y miRTarBase.")

    # 2. Detalles por Gen Core (Top 5) — enriquecido con evidencia multi-fuente
    gene_paragraphs = []
    for i, g in enumerate(common[:5]):
        d = gene_details.get(g, {})
        full_name = d.get('full_name', g)
        pathology = d.get('pathology', "Nodo regulador identificado.")
        para = (f"En relación al gen {g} ({full_name}), la literatura científica lo describe como {pathology.lower().replace('.', '')}. "
                f"Diversas investigaciones destacan su papel como un nodo de control esencial.")

        # La investigación en vivo de este gen (PubMed, OMIM, ClinVar, ensayos).
        gr = gene_research.get(g, {}) or {}

        # Separación texto / cita (v1.44): la descripción de cada ruta es una
        # NOTA EDITORIAL del equipo — contexto, sin cita — y el PMID que se
        # imprimía junto a ella era solo el primer resultado de una búsqueda
        # por gen + término, sin verificar que respalde la frase. Van en
        # bloques distintos y se dice así.
        evidence_list = d.get('clinical_evidence', [])
        if evidence_list:
            route_lines = [f"• {ev['source']} {ev['term']}: {ev['desc']}" for ev in evidence_list]
            para += "\n\nNotas editoriales sobre las rutas (contexto del equipo, sin cita):\n" + "\n".join(route_lines)
            live_refs = gr.get("pubmed", []) or []
            if live_refs:
                live_lines = [
                    f"• {r['term']}: primer resultado en PubMed para {g} + «{r['term']}» — PMID: {r['pmid']}"
                    for r in live_refs
                ]
                para += ("\n\nPrimeros resultados en PubMed (búsqueda por gen + término; "
                         "no se verificó que respalden la nota anterior):\n" + "\n".join(live_lines))
        else:
            para += f" Su regulación por el panel de miRNAs analizado tiene consecuencias directas en la estabilidad tisular [{i+1}].\n"

        # Evidencia genómica multi-fuente (OMIM / ClinVar / ClinicalTrials)
        # Se construye un párrafo narrativo continuo, omitiendo silenciosamente
        # las fuentes sin resultados. Si NINGUNA fuente devolvió datos, no se
        # añade el bloque.
        omim_list = gr.get("omim", []) or []
        clinvar = gr.get("clinvar", {}) or {}
        trials_list = gr.get("trials", []) or []

        evidence_phrases = []

        # OMIM: enfermedades hereditarias asociadas
        valid_omim = [o for o in omim_list if o.get("title")]
        if valid_omim:
            if len(valid_omim) == 1:
                o = valid_omim[0]
                evidence_phrases.append(
                    f"el catálogo OMIM (#{o['omim_id']}) registra la entrada \"{o['title']}\" "
                    f"como referencia mendeliana del gen"
                )
            else:
                ids_str = ", ".join(f"#{o['omim_id']}" for o in valid_omim[:3])
                evidence_phrases.append(
                    f"OMIM registra {len(valid_omim)} entradas asociadas al gen ({ids_str}), "
                    f"reforzando su relevancia en enfermedades hereditarias"
                )

        # ClinVar: variantes patogénicas
        cv_count = clinvar.get("count", 0)
        if cv_count > 0:
            if cv_count == 1:
                evidence_phrases.append("ClinVar reporta 1 variante patogénica documentada en este locus")
            elif cv_count < 10:
                evidence_phrases.append(
                    f"ClinVar reporta {cv_count} variantes patogénicas/probablemente patogénicas "
                    f"documentadas en este locus"
                )
            else:
                evidence_phrases.append(
                    f"ClinVar registra {cv_count} variantes patogénicas/probablemente patogénicas "
                    f"clasificadas para este gen — un volumen que respalda su impacto clínico"
                )

        # ClinicalTrials: ensayos activos con el gen como protagonista
        if trials_list:
            if len(trials_list) == 1:
                t = trials_list[0]
                evidence_phrases.append(
                    f"actualmente existe 1 ensayo clínico ({t['nct_id']}) centrado en este gen: "
                    f"\"{t['title']}\""
                )
            else:
                ncts = ", ".join(t['nct_id'] for t in trials_list[:3])
                evidence_phrases.append(
                    f"actualmente hay {len(trials_list)} ensayos clínicos en ClinicalTrials.gov "
                    f"centrados en este gen ({ncts})"
                )

        # Componer el párrafo narrativo solo si hay al menos una fuente con datos
        if evidence_phrases:
            if len(evidence_phrases) == 1:
                evidence_narrative = f"Como evidencia genómica complementaria, {evidence_phrases[0]}."
            elif len(evidence_phrases) == 2:
                evidence_narrative = (
                    f"Como evidencia genómica complementaria, {evidence_phrases[0]}; "
                    f"asimismo, {evidence_phrases[1]}."
                )
            else:
                head = "; ".join(evidence_phrases[:-1])
                evidence_narrative = (
                    f"Como evidencia genómica complementaria, {head}; "
                    f"finalmente, {evidence_phrases[-1]}."
                )
            para += f"\n\n{evidence_narrative}"

        if d.get('conclusion'):
            para += f"\n\nConclusión clínica: {d['conclusion']}"

        gene_paragraphs.append(para)

    academic_text = p1 + "\n\n" + "\n\n".join(gene_paragraphs)

    # 3. Contexto Funcional Global — cada ruta menciona su evidencia PubMed
    # cuando está disponible, integrando la cita al narrativa en lugar de
    # dejarla huérfana en la bibliografía.
    route_details = []
    for item in enrichment[:15]:
        term = item.get('Term', 'Ruta biológica')
        source = item.get('Source', 'Base de datos')
        pval = item.get('Pval', 1.0)
        evidence = item.get('Evidence')

        # Cualificador de significancia según p-adj
        if pval < 0.001:
            significance = "alta significancia estadística"
        elif pval < 0.01:
            significance = "significancia robusta"
        else:
            significance = "significancia funcional"

        para = (
            f"La vía de {term} ({source}, p-adj={pval:.4f}) destaca por su {significance}. "
            f"Las investigaciones asocian esta ruta con la respuesta adaptativa celular, "
            f"integrando las señales de los genes core para mantener el equilibrio fisiológico."
        )
        # Si el enriquecimiento devolvió evidencia PubMed, anclar la cita aquí
        if evidence and evidence.get("id"):
            para += f" La evidencia bibliográfica de respaldo está catalogada en PubMed (PMID: {evidence['id']})."
        route_details.append(para)

    functional_text = "Contexto Funcional y Rutas Biológicas Globales:\n\n" + "\n\n".join(route_details)

    # 4. Referencias (Desglosadas para paginación)
    # Cita estilo Vancouver adaptada al tipo de fuente (PubMed / OMIM / ClinVar /
    # ClinicalTrials). El identificador correcto va al final junto con la URL.
    ref_lines = []
    for ref in references:
        parts = []
        authors = ref.get("authors", "")
        year = ref.get("year", "")
        journal = ref.get("journal", "")
        title = ref.get("title", "").rstrip(".")
        ref_type = ref.get("ref_type", "pubmed")
        source = ref.get("source", "")  # contiene contexto (gen / ruta)

        if authors:
            parts.append(f"{authors}.")
        if title:
            parts.append(f"{title}.")
        if journal and year:
            parts.append(f"{journal}. {year}.")
        elif journal:
            parts.append(f"{journal}.")
        elif year:
            parts.append(f"{year}.")

        # Identificador específico por tipo de fuente
        if ref_type == "pubmed" and ref.get("pmid"):
            parts.append(f"PMID: {ref['pmid']}.")
        elif ref_type == "omim":
            parts.append(f"OMIM: {ref.get('source','').replace('OMIM:','').split(' ')[0]}.")
        elif ref_type == "clinvar":
            parts.append(f"ClinVar.")
        elif ref_type == "clinicaltrials":
            nct = ref.get('source','').replace('NCT:','').split(' ')[0]
            parts.append(f"NCT: {nct}.")

        # Contexto: por qué esta referencia está aquí (gen core / ruta)
        if source and ref_type == "pubmed":
            parts.append(f"[{source}]")

        parts.append(ref["url"])
        citation = " ".join(parts)
        ref_lines.append(f"[{ref['id']}] {citation}")

    references_text = "Bibliografía:\n\n" + "\n\n".join(ref_lines) if ref_lines else ""

    return {
        "academic": academic_text,
        "functional": functional_text,
        "references": references_text
    }

def load_mirtarbase_gmt(path):
    """Load the local miRTarBase gene-set library (Harmonizome GMT:
    one row per miRNA: ``miRNA \\t MIRT-id \\t gene1 \\t gene2 ...``).
    Returns {miRNA name: set(genes)}. Pure, no network."""
    import gzip
    db = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0].strip()
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if name and genes:
                db.setdefault(name, set()).update(genes)
    return db


def get_validated_local(mirna, db, species="hsa"):
    """Validated targets for one miRNA from the local miRTarBase library.
    Tries the cleaned name and the name without the -5p/-3p suffix.
    The library is human (hsa-); other species return empty here and
    fall back to the live Harmonizome API."""
    if not db:
        return set()
    name = clean_mirna_name_cloud(mirna, species)
    if name in db:
        return set(db[name])
    base = re.sub(r"-(5p|3p)$", "", name, flags=re.IGNORECASE)
    if base in db:
        return set(db[base])
    return set()


def load_gmt_library(path):
    """Parse one GMT gene-set library (Enrichr format: one row per term:
    ``Term \\t description \\t gene1 \\t gene2 ...``).
    Returns {term: set(genes)}."""
    lib = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0].strip()
            genes = {g.strip() for g in parts[2:] if g.strip()}
            if term and genes:
                lib[term] = genes
    return lib


def benjamini_hochberg(pvalues):
    """Benjamini-Hochberg adjusted p-values. Returns a list aligned with
    the input order. Pure stdlib."""
    n = len(pvalues)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvalues[i])
    adj = [0.0] * n
    prev = 1.0
    for rank in range(n - 1, -1, -1):
        i = order[rank]
        prev = min(prev, pvalues[i] * n / (rank + 1))
        adj[i] = min(1.0, prev)
    return adj


def fisher_enrichment(genes, library, background=None, top_n=8, p_adj_max=0.05):
    """Over-representation analysis of `genes` against one GMT library.
    Fisher's exact test (scipy) + Benjamini-Hochberg. Returns up to
    `top_n` terms with adjusted p < `p_adj_max`, sorted by adjusted p.
    Each entry: {Term, Pval (adjusted), RawPval, Overlap, TermSize}.

    CORREGIDO (auditoría externa, 2026-11-19): antes se descartaban los
    términos con overlap 0 ANTES de BH, así que m era el nº de términos
    con solape y no el total de la biblioteca (subestimación de ~28× en
    GO, con falsos positivos en 40/40 listas nulas a p_adj < 0.3); y
    `genes` no se intersecaba con el universo, así que genes ausentes de
    la GMT inflaban b y podían hacer d negativo. Ahora: genes se recorta
    al universo, TODOS los términos no vacíos de la biblioteca entran en
    BH (los de overlap 0 con p = 1.0) y el corte es p_adj < 0.05."""
    from scipy.stats import fisher_exact
    genes = set(genes)
    if not genes or not library:
        return []
    universe = background or set().union(*library.values())
    n_bg = len(universe)
    if n_bg == 0:
        return []
    genes = genes & universe
    if not genes:
        return []
    rows = []
    for term, term_genes in library.items():
        term_genes = term_genes & universe
        if not term_genes:
            continue
        overlap = len(genes & term_genes)
        a = overlap
        b = len(genes - term_genes)
        c = len(term_genes - genes)
        d = n_bg - a - b - c
        if d < 0:
            continue
        if overlap == 0:
            # Sin solape: p = 1.0. Cuenta para m en BH pero nunca es
            # significativo; se descarta de la salida, no del ajuste.
            rows.append({"Term": term, "Pval": 1.0, "RawPval": 1.0,
                         "Overlap": 0, "TermSize": len(term_genes)})
            continue
        _, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({"Term": term, "Pval": p, "RawPval": p,
                     "Overlap": overlap, "TermSize": len(term_genes)})
    if not rows:
        return []
    adj = benjamini_hochberg([r["Pval"] for r in rows])
    for r, p_adj in zip(rows, adj):
        r["Pval"] = p_adj
    rows = [r for r in rows if r["Pval"] < p_adj_max and r["Overlap"] > 0]
    rows.sort(key=lambda r: r["Pval"])
    return rows[:top_n]


def local_enrichment(genes, libraries, top_total=24):
    """Enrichment against several local GMT libraries (offline).
    `libraries` maps a label (KEGG/Reactome/WikiPathways/GO) to a
    {term: set(genes)} library. Returns the combined top terms in the
    same shape the live Enrichr path produces."""
    out = []
    for label, library in libraries.items():
        for row in fisher_enrichment(genes, library):
            out.append({"Term": row["Term"], "Source": label,
                        "Pval": row["Pval"], "RawTerm": row["Term"]})
    out.sort(key=lambda r: r["Pval"])
    return out[:top_total]


def build_export(mirnas, common_genes, gene_details, enrichment,
                 validation, gene_sets=None, schema_version="1.0"):
    """Structured export of one analysis: a stable, documented schema
    (genes table, miRNA->gene edges, enrichment table, validation table)
    for downstream tools, Cytoscape and notebooks. Pure."""
    import datetime
    genes = []
    for g in common_genes:
        d = gene_details.get(g, {}) if gene_details else {}
        v = (validation or {}).get("genes", {}).get(g, {})
        genes.append({
            "gene": g,
            "full_name": d.get("full_name", g),
            "validation_score": v.get("score"),
            "validation_level": v.get("level"),
            "validation_layers": v.get("layers", []),
        })
    edges = []
    if gene_sets and mirnas:
        for m, targets in zip(mirnas, gene_sets):
            for g in sorted(set(targets) & set(common_genes)):
                edges.append({"mirna": m, "gene": g, "evidence": "TargetScan"})
    return {
        "schema_version": schema_version,
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "panel": list(mirnas),
        "data_sources": (validation or {}).get("data_sources", {}),
        "genes": genes,
        "edges": edges,
        "enrichment": [
            {"term": e.get("Term"), "source": e.get("Source"),
             "p_adj": e.get("Pval")} for e in (enrichment or [])
        ],
        "validation_summary": (validation or {}).get("summary", {}),
    }


def heatmap_matrix(mirnas, common_genes, target_sets):
    """Genes x miRNAs binary matrix: 1 when the gene is a target of the
    miRNA. Returns (genes, mirnas, matrix) with genes sorted."""
    genes = sorted(common_genes)
    matrix = [[1 if g in s else 0 for s in target_sets] for g in genes]
    return genes, list(mirnas), matrix


def heatmap_matrix_union(mirnas, target_sets, top_n=20):
    """Genes x miRNAs binary matrix over the UNION of target sets, ranked
    by how many panel miRNAs target each gene (desc). Core genes
    (targeted by all) sit on top; genes targeted by a single miRNA at
    the bottom — so the figure actually varies instead of being a solid
    block (the old core-only matrix was all 1s by construction).
    Returns (genes, mirnas, matrix, union_size)."""
    union = set()
    for s in target_sets:
        union |= set(s)
    counts = {g: sum(1 for s in target_sets if g in s) for g in union}
    ranked = sorted(union, key=lambda g: (-counts[g], g))
    genes = ranked[:top_n]
    matrix = [[1 if g in s else 0 for s in target_sets] for g in genes]
    return genes, list(mirnas), matrix, len(union)


def create_heatmap(mirnas, common_genes, target_sets, top_n=20):
    """Heatmap PNG (base64) of the genes x miRNAs target matrix, over the
    UNION of targets (top `top_n` by miRNA count) so the pattern varies.
    `common_genes` is kept for signature compatibility: the core genes
    always sit on top, because every panel miRNA targets them."""
    import io as _io
    genes, mirnas, matrix, union_size = heatmap_matrix_union(
        mirnas, target_sets, top_n=top_n)
    if not genes or not mirnas:
        return None
    arr = np.array(matrix, dtype=float)
    fig, ax = plt.subplots(figsize=(max(4, len(mirnas) * 1.2),
                                    max(3, len(genes) * 0.45)))
    sns.heatmap(arr, ax=ax, cmap="YlGnBu", cbar=True, linewidths=0.5,
                xticklabels=mirnas, yticklabels=genes)
    ax.set_title(f"Genes diana por miRNA — top {len(genes)} de {union_size} (1 = diana)")
    plt.tight_layout()
    buf = _io.BytesIO()
    fig.savefig(buf, format="png", dpi=100)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode()


def reverse_lookup(gene, targetscan_db):
    """miRNAs (families) whose local TargetScan set contains `gene`.
    Pure scan of the local database. Returns a sorted list of miRNA names."""
    gene = (gene or "").strip().upper()
    if not gene:
        return []
    hits = []
    for mirna, targets in (targetscan_db or {}).items():
        if gene in targets:
            hits.append(mirna)
    return sorted(hits)


def normalize_isomir(name, species="hsa"):
    """Map an isomiR notation to its canonical miRNA name.
    Strips 3' trimming/extension marks (+N/-N) ONLY when they follow the
    -5p/-3p arm suffix (e.g. hsa-miR-21-5p+2), shift annotations (|...),
    and sequence suffixes, then applies the canonical cleaner.
    Bare numbers are part of the name (hsa-miR-375 stays hsa-miR-375).
    Quantifying isomiRs from FASTQ is out of scope; this only normalizes
    the name so an isomiR label resolves to its canonical miRNA."""
    n = (name or "").strip()
    n = n.split("|")[0]
    n = re.sub(r"(-[53]p)[+-]\d+$", r"\1", n)
    n = re.sub(r"[^A-Za-z0-9\-]+$", "", n)
    return clean_mirna_name_cloud(n, species)


# ── VALIDATION OF EVERYTHING ──────────────────────────────────────────────
# Every gene in the core gets a validation score from independent evidence
# layers, each one verifiable on its own: predicted (TargetScan), experimentally
# validated (miRTarBase), literature (PubMed PMIDs), disease (OMIM),
# variants (ClinVar), trials (ClinicalTrials.gov). The score is a summary;
# the layers are the proof. Nothing here claims what the data does not show:
# "validated" means the gene is in miRTarBase's validated set for at least
# one panel miRNA, checked independently of the intersection itself (the
# analysis uses TargetScan as primary and miRTarBase only as fallback, so
# this layer is a separate, second verification, not a repetition of the
# intersection).

_VALIDATION_WEIGHTS = {
    "predicted": 25,   # TargetScan prediction
    "validated": 30,   # miRTarBase experimental validation
    "pubmed_each": 5,  # per PubMed reference, capped
    "pubmed_cap": 15,
    "omim": 10,
    "clinvar": 10,
    "trials": 10,
}


def validation_score(gene, research, predicted, validated):
    """Score one gene 0-100 from its evidence layers. Pure and
    deterministic. `research` is the gene_research entry for the gene
    ({pubmed: [...], omim: [...], clinvar: {count}, trials: [...]})."""
    research = research or {}
    layers, missing = [], []
    score = 0
    if predicted:
        layers.append("TargetScan (predicho)")
        score += _VALIDATION_WEIGHTS["predicted"]
    else:
        missing.append("TargetScan (predicho)")
    if validated:
        layers.append("miRTarBase (validado experimentalmente)")
        score += _VALIDATION_WEIGHTS["validated"]
    else:
        missing.append("miRTarBase (validado experimentalmente)")
    n_pubmed = len(research.get("pubmed") or [])
    if n_pubmed:
        layers.append(f"PubMed ({n_pubmed} refs)")
        score += min(_VALIDATION_WEIGHTS["pubmed_cap"],
                     n_pubmed * _VALIDATION_WEIGHTS["pubmed_each"])
    else:
        missing.append("PubMed")
    if research.get("omim"):
        layers.append(f"OMIM ({len(research['omim'])} entradas)")
        score += _VALIDATION_WEIGHTS["omim"]
    else:
        missing.append("OMIM")
    clinvar_count = (research.get("clinvar") or {}).get("count", 0)
    if clinvar_count:
        layers.append(f"ClinVar ({clinvar_count} variantes)")
        score += _VALIDATION_WEIGHTS["clinvar"]
    else:
        missing.append("ClinVar")
    if research.get("trials"):
        layers.append(f"ClinicalTrials ({len(research['trials'])} ensayos)")
        score += _VALIDATION_WEIGHTS["trials"]
    else:
        missing.append("ClinicalTrials")
    score = min(100, score)
    # "Alta" exige validación experimental de la INTERACCIÓN (miRTarBase):
    # PubMed, OMIM, ClinVar y ClinicalTrials describen lo importante que es
    # el gen, no que el microARN lo regule. Sin esa capa, un gen muy estudiado
    # llegaría a "alta" por popularidad, no por evidencia de la interacción.
    gated = score >= 70 and not validated
    level = "media" if gated else ("alta" if score >= 70 else ("media" if score >= 40 else "baja"))
    return {
        "gene": gene,
        "score": score,
        "level": level,
        "layers": layers,
        "missing": missing,
        "explanation": (
            f"{gene}: validación {level} ({score}/100) — "
            + (", ".join(layers) if layers else "sin capas de evidencia")
            + (f"; falta: {', '.join(missing)}" if missing else "")
            + ("; nivel limitado a media: sin validación experimental de la interacción" if gated else "")
        ),
    }


def build_validation_report(mirnas, common_genes, gene_research,
                            local_by_mirna, remote_by_mirna, data_versions):
    """Validation dossier for one analysis: per-gene scores plus the exact
    database versions the result rests on. `local_by_mirna` / `remote_by_mirna`
    map each panel miRNA to its TargetScan / miRTarBase target sets."""
    import datetime
    predicted_set = set().union(*local_by_mirna.values()) if local_by_mirna else set()
    validated_set = set().union(*remote_by_mirna.values()) if remote_by_mirna else set()
    gene_research = gene_research or {}
    genes = {}
    for gene in common_genes:
        genes[gene] = validation_score(
            gene, gene_research.get(gene),
            predicted=gene in predicted_set,
            validated=gene in validated_set,
        )
    scores = [g["score"] for g in genes.values()]
    levels = [g["level"] for g in genes.values()]
    return {
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "mirnas": list(mirnas),
        "data_sources": dict(data_versions or {}),
        "method": (
            "Cada gen se valida por capas independientes y verificables: "
            "predicción TargetScan, validación experimental miRTarBase "
            "(comprobada por separado de la intersección), PubMed (PMID), "
            "OMIM, ClinVar y ClinicalTrials.gov. La puntuación resume las "
            "capas; las capas son la prueba. El nivel alta exige validación "
            "experimental de la interacción; sin ella, el máximo es media. "
            "Las referencias de PubMed se deduplican por PMID antes de puntuar. "
            "La validación experimental se comprueba contra el snapshot local "
            "de miRTarBase (Harmonizome): un microARN sin entrada en ese "
            "snapshot no aporta validación aunque la literatura reporte dianas "
            "validadas, y su gen queda en media como máximo. Las notas "
            "editoriales sobre las rutas son contexto del equipo y no llevan "
            "cita; los PMIDs de PubMed son los primeros resultados de la "
            "búsqueda por gen + término y no se verificó que respalden la nota."
        ),
        "genes": genes,
        "summary": {
            "n_genes": len(genes),
            "high": levels.count("alta"),
            "medium": levels.count("media"),
            "low": levels.count("baja"),
            "mean_score": round(sum(scores) / len(scores), 1) if scores else None,
            "not_validated": [g for g, s in genes.items() if s["level"] == "baja"],
        },
    }
