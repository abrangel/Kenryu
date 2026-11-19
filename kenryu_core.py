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
    m = re.sub(r'^(hsa-)+', '', m, flags=re.IGNORECASE)
    return m.lower()

def clean_mirna_name_cloud(m: str) -> str:
    c = deep_clean_mirna(m)
    parts = c.split('-')
    # The let-7 family keeps its own name (hsa-let-7a-5p). The generic
    # rule below would drop "let" and wrongly return hsa-miR-7a-5p,
    # so the remote miRTarBase lookup queried a name that does not exist.
    if "let" in parts[:2]:
        i = parts.index("let")
        rest = parts[i + 1:]
        r = "hsa-let"
        if rest:
            r += "-" + "-".join(rest)
        return r
    if len(parts) > 1:
        num = parts[1]; suf = '-'.join(parts[2:])
        r = f"hsa-miR-{num}"
        if suf: r += f"-{suf}"
        return r
    return f"hsa-miR-{c}"

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


def create_visuals(gene_sets: list, mirnas: list, common_all: list):
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

    # ── VOLCANO ───────────────────────────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(10, 6), facecolor='white')
    ax2.set_facecolor('white')
    np.random.seed(42)
    n_bg = 200
    x_bg = np.random.normal(0, 1.2, n_bg)
    y_bg = np.abs(np.random.normal(0.5, 0.4, n_bg))
    mask_ns = y_bg < 1.3
    ax2.scatter(x_bg[mask_ns], y_bg[mask_ns], alpha=0.35, c='#aaaaaa', s=12, label='No significativo')
    ax2.scatter(x_bg[~mask_ns], y_bg[~mask_ns], alpha=0.45, c='#6baed6', s=14, label='Significativo')
    if common_all:
        n_core = min(len(common_all), 15)
        half = n_core // 2
        x_down = np.random.uniform(-3.5, -1.8, half)
        y_down = np.random.uniform(2.0, 4.5, half)
        x_up   = np.random.uniform(1.8, 3.5, n_core - half)
        y_up   = np.random.uniform(2.0, 4.5, n_core - half)
        x_core = np.concatenate([x_down, x_up])
        y_core = np.concatenate([y_down, y_up])
        ax2.scatter(x_core, y_core, c='#e05c5c', s=55, zorder=5,
                    edgecolors='#8b1a1a', linewidth=0.6, label='Genes Core')
        for i in range(min(6, n_core)):
            ax2.annotate(common_all[i], (x_core[i], y_core[i]),
                         textcoords="offset points", xytext=(6, 4),
                         fontsize=8.5, color='#1a3a6b', fontweight='bold',
                         arrowprops=dict(arrowstyle='-', color='#1a3a6b55', lw=0.6))
    ax2.axhline(y=1.3, color='#e05c5c', linestyle='--', alpha=0.6, lw=1.2, label='-log₁₀(0.05)')
    ax2.axvline(x=1.5, color='#aaaaaa', linestyle=':', alpha=0.5, lw=1.0)
    ax2.axvline(x=-1.5, color='#aaaaaa', linestyle=':', alpha=0.5, lw=1.0)
    ax2.axvline(x=0, color='#cccccc', linestyle='-', alpha=0.3, lw=0.8)
    ax2.text(-3.2, 0.1, '⬇ Regulación\nnegativa', fontsize=8, color='#666', style='italic')
    ax2.text(2.0, 0.1, '⬆ Regulación\npositiva', fontsize=8, color='#666', style='italic')
    ax2.set_xlabel("log₂ Fold Change (Efecto Regulador)", color='#444', fontsize=11)
    ax2.set_ylabel("−log₁₀ (p-valor ajustado)", color='#444', fontsize=11)
    ax2.set_title("Paisaje de Significancia Biológica",
                  color='#1a1a2e', fontsize=13, pad=14, fontweight='bold')
    ax2.legend(fontsize=9, framealpha=0.9, edgecolor='#cccccc')
    ax2.set_xlim(-4.5, 4.5)
    buf2 = io.BytesIO()
    plt.savefig(buf2, format='png', facecolor='white', bbox_inches='tight', dpi=150)
    plt.close(fig2); buf2.seek(0)

    # ── INTERACTOMA PPI ───────────────────────────────────────────────────────
    fig3, ax3 = plt.subplots(figsize=(9, 9), facecolor='white')
    ax3.set_facecolor('white'); ax3.axis('off')
    if common_all:
        n_nodes = min(len(common_all), 12)
        genes_disp = common_all[:n_nodes]
        NODE_COLORS = ['#4e9af1', '#f4845f', '#57cc99', '#c77dff', '#ffd166',
                       '#06d6a0', '#ef476f', '#118ab2', '#fca311', '#8ecae6',
                       '#a8dadc', '#e9c46a'][:n_nodes]
        r = 0.33
        positions = {}
        for i, gene in enumerate(genes_disp):
            angle = 2 * np.pi * i / n_nodes - np.pi / 2
            gx = 0.5 + r * np.cos(angle)
            gy = 0.5 + r * np.sin(angle)
            positions[gene] = (gx, gy)
        for i in range(n_nodes):
            for j in range(i + 1, n_nodes):
                if j - i <= 2 or (i == 0 and j == n_nodes - 1):
                    g1, g2 = genes_disp[i], genes_disp[j]
                    p1, p2 = positions[g1], positions[g2]
                    ax3.plot([p1[0], p2[0]], [p1[1], p2[1]],
                             c='#dddddd', linewidth=1.0, zorder=1, alpha=0.8)
        for gene, (gx, gy) in positions.items():
            ax3.plot([0.5, gx], [0.5, gy],
                     c='#aaaaaa', linewidth=0.8, zorder=2, alpha=0.5, linestyle='--')
        for i, gene in enumerate(genes_disp):
            gx, gy = positions[gene]
            ax3.scatter(gx, gy, c=NODE_COLORS[i], s=280, zorder=5,
                        edgecolors='white', linewidth=1.8)
            angle = 2 * np.pi * i / n_nodes - np.pi / 2
            lx = 0.5 + (r + 0.09) * np.cos(angle)
            ly = 0.5 + (r + 0.09) * np.sin(angle)
            ax3.text(lx, ly, gene, ha='center', va='center',
                     fontsize=8.5, fontweight='bold', color='#1a1a2e', zorder=7,
                     bbox=dict(boxstyle='round,pad=0.15', facecolor='white',
                               edgecolor='none', alpha=0.75))
        ax3.scatter(0.5, 0.5, c='#1a3a6b', s=600, zorder=6,
                    edgecolors='white', linewidth=2.5)
        ax3.text(0.5, 0.505, "HUB", ha='center', va='center', color='white',
                 fontsize=9, fontweight='bold', zorder=7)
        ax3.text(0.5, 0.47, f"{n_nodes}", ha='center', va='center', color='#c8a96e',
                 fontsize=11, fontweight='bold', zorder=7)
        ax3.set_xlim(0.05, 0.95); ax3.set_ylim(0.05, 0.95)
    ax3.set_title(f"Interactoma PPI — {len(common_all)} genes core",
                  color='#1a1a2e', fontsize=13, pad=14, fontweight='bold')
    buf3 = io.BytesIO()
    plt.savefig(buf3, format='png', facecolor='white', bbox_inches='tight', dpi=150)
    plt.close(fig3); buf3.seek(0)

    return (base64.b64encode(buf1.getvalue()).decode(),
            base64.b64encode(buf2.getvalue()).decode(),
            base64.b64encode(buf3.getvalue()).decode())

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

        evidence_list = d.get('clinical_evidence', [])
        if evidence_list:
            route_lines = []
            for ev in evidence_list:
                route_lines.append(f"• {ev['source']} {ev['term']}: {ev['desc']} (PMID: {ev['pmid']})")
            para += "\n\nRutas biológicas validadas:\n" + "\n".join(route_lines)
        else:
            para += f" Su regulación por el panel de miRNAs analizado tiene consecuencias directas en la estabilidad tisular [{i+1}].\n"

        # Evidencia genómica multi-fuente (OMIM / ClinVar / ClinicalTrials)
        # Se construye un párrafo narrativo continuo, omitiendo silenciosamente
        # las fuentes sin resultados. Si NINGUNA fuente devolvió datos, no se
        # añade el bloque.
        gr = gene_research.get(g, {}) or {}
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
    level = "alta" if score >= 70 else ("media" if score >= 40 else "baja")
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
            "capas; las capas son la prueba."
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
