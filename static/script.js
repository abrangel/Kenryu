// Detección automática del Backend para Local, GitHub Pages y HuggingFace
const API_URL = (window.location.hostname === 'localhost' || window.location.hostname === '127.0.0.1')
  ? `http://localhost:${window.location.port || 3001}`
  : (window.location.hostname.includes('github.io'))
    ? 'https://kenryu007-bioinformatica.hf.space' // URL directa del backend en Hugging Face
    : ''; // En HuggingFace: rutas relativas (el mismo servidor sirve el frontend)

let lastData = null;

// H-08 (v1.47): escapar datos externos y de usuario antes de insertarlos
// con innerHTML. Sin esto, un gen, un término de Enrichr o un mensaje de
// error con HTML se ejecuta como código en la página.
function escapeHtml(s) {
  return String(s ?? '')
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;')
    .replace(/'/g, '&#39;');
}
// Unreleased (2026-11-19): especificidad exacta (cola hipergeométrica: p de
// Fisher sin ajustar, con el universo de TargetScan como población). Se
// muestra "<1%" por debajo del 1% y un decimal entre 1% y 10%, para no
// borrar el umbral del 5% (un 4,7% no debe leerse como 5%).
const SPEC_TITLE = 'Probabilidad exacta de que listas aleatorias de genes del mismo tamaño, extraídas del universo de TargetScan, solapen con este término al menos tanto como la lista observada (cola hipergeométrica, sin ajustar por pruebas múltiples: guarda contra términos genéricos, no es una prueba de significancia). 100% = genérico; <1% = específico.';
function fmtSpec(x) {
  if (typeof x !== 'number' || !isFinite(x)) return '—';
  if (x < 0.01) return '<1%';
  if (x < 0.10) return `${(x * 100).toFixed(1)}%`;
  return `${Math.round(x * 100)}%`;
}
let pageCount = 1;
let sectionCount = 4;
let isPreview = false;

// Poblar años dinámicamente desde 2001 hasta el actual
function populateYearSelect() {
  const select = document.getElementById('pubmed-years');
  if (!select) return;
  const currentYear = new Date().getFullYear();
  let html = '';
  for (let y = currentYear; y >= 2001; y--) {
    html += `<option value="${y}">Año ${y}</option>`;
  }
  select.innerHTML = html;
  // Por defecto, ventana de 10 años (año de inicio = hace 9 años).
  select.value = String(currentYear - 9);
}

function switchView(v) {
  document.querySelectorAll('.view').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
  
  document.getElementById('view-' + v).classList.add('active');
  const navId = v === 'dashboard' ? 'nav-dash' : 'nav-report';
  const navElem = document.getElementById(navId);
  if (navElem) navElem.classList.add('active');
  
  document.getElementById('view-title').textContent = v === 'dashboard' ? 'Panel de Análisis' : 'Editor de Informe Pro';
  
  if (v === 'report' && lastData && document.getElementById('report-canvas-content').innerHTML === '') initReportWithData(lastData);
}

function addLog(msg, type = '') {
  const log = document.getElementById('process-logs');
  if (!log) return;
  const d = document.createElement('div');
  d.className = 'log-line';
  const ts = new Date().toLocaleTimeString('es-ES', {hour:'2-digit',minute:'2-digit',second:'2-digit'});
  const cls = type === 'ok' ? 'ok' : (type === 'err' ? 'err' : '');
  d.innerHTML = `<span class="log-ts">${ts}</span><span class="log-msg ${cls}">${escapeHtml(msg)}</span>`;
  log.appendChild(d);
  log.scrollTop = log.scrollHeight;
}

async function analyzePro() {
  const btn = document.querySelector('.btn-run');
  const input = document.getElementById('mirna-input').value;
  const mirnas = input.split(',').map(m => m.trim()).filter(Boolean);
  const mode = document.getElementById('consensus-mode').value;
  const startYear = document.getElementById('pubmed-years').value;
  // La UI elige un AÑO DE INICIO; la API espera una EDAD en años (v1.45:
  // antes se enviaba el año calendario 2026 como `years`, lo que anulaba
  // el filtro temporal de PubMed).
  const years = Math.max(1, new Date().getFullYear() - parseInt(startYear) + 1);
  
  if (!mirnas.length) return alert("Ingrese al menos un miRNA.");

  btn.disabled = true;
  btn.innerHTML = '<i class="fas fa-spinner fa-spin"></i> Analizando...';
  addLog(`Sincronizando motor para búsqueda desde ${startYear}...`, 'ok');

  try {
    const res = await fetch(`${API_URL}/api/v1/analyze`, {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({mirnas, years: years, mode: mode,
        require_experimental: document.getElementById('require-experimental')?.checked || false})
    });
    if (!res.ok) throw new Error(`HTTP ${res.status}: ${await res.text()}`);
    const data = await res.json();
    lastData = data;
    renderDashboard(data);
    
    addLog(`✓ Análisis completo — ${data.common_genes.length} biomarcadores.`, 'ok');
    
    initReportWithData(data);
  } catch (e) {
    addLog(`✕ ERROR: ${e.message}`, 'err');
  } finally {
    btn.disabled = false;
    btn.innerHTML = '<i class="fas fa-play"></i> Ejecutar';
  }
}

function renderDashboard(data) {
  try {
    document.getElementById('results-grid').style.display = 'grid';
    document.getElementById('gene-count-badge').textContent = `${data.common_genes ? data.common_genes.length : 0} genes`;

    const sysColor = s => (s && (s.includes('Cardio') || s.includes('Metab'))) ? 'var(--red)'
      : (s && s.includes('Neuro')) ? 'var(--blue)'
      : (s && s.includes('Onco')) ? 'var(--gold)'
      : 'var(--teal)';
      
    const genesEl = document.getElementById('core-genes-list');
    if (genesEl && data.common_genes) {
      genesEl.innerHTML = data.common_genes.map(g => {
        const d = data.gene_details && data.gene_details[g];
        const dot = sysColor(d && d.system);
        // H-08: el gen va en data-gene (escapado), no en un onclick con
        // interpolación — una comilla en el símbolo rompería el handler.
        return `<span class="gene-pill" data-gene="${escapeHtml(g)}" onclick="openGenePanel(this.dataset.gene)"><i class="fas fa-circle" style="font-size:5px;color:${dot};margin-right:5px;"></i>${escapeHtml(g)}</span>`;
      }).join('');
    }

    // Prueba de convergencia (Unreleased): ¿es inusual el número de
    // genes en común? p como rango entre tolerancias, nunca un solo valor.
    const cvEl = document.getElementById('convergence-container');
    if (cvEl) {
      const cv = data.convergence;
      if (cv) {
        let html = '';
        if (cv.skipped) {
          html = `<div style="color:var(--text-dim); font-size:12px; line-height:1.6;">Prueba de convergencia no calculada (${escapeHtml(cv.skipped)}).</div>`;
        } else {
          const tolText = Object.entries(cv.p_by_tolerance).map(([t, p]) => `±${Math.round(parseFloat(t) * 100)}%: ${p.toFixed(3)}`).join(' · ');
          html = `<div style="font-size:13px; line-height:1.6;">${cv.k_names} nombres de miRNA, ${cv.k_effective} familias efectivas, ${cv.genes_observed} genes en común.<br>p entre ${cv.p_min.toFixed(3)} y ${cv.p_max.toFixed(3)} (${tolText}; resolución ${cv.resolution.toFixed(4)}, ${cv.perms} permutaciones por tolerancia).<br><i>${escapeHtml(cv.label)}</i><br><span style="color:var(--text-dim); font-size:11px;">${escapeHtml(cv.limit)}</span></div>`;
          if (cv.same_family_warning) {
            html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">⚠️ Hay miRNAs de la misma familia: el panel tiene ${cv.k_names} nombres pero ${cv.k_effective} semillas independientes.</div>`;
          }
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">Manifiesto: KENRYU ${escapeHtml(data.engine_version || '')} (commit ${escapeHtml(data.commit || '')}), semilla de convergencia ${cv.seed}.</div>`;
        }
        cvEl.innerHTML = html;
      } else {
        cvEl.innerHTML = '<div style="color:var(--text-faint); padding:20px; text-align:center;">Sin prueba de convergencia en esta respuesta (versión anterior del motor).</div>';
      }
    }

    const setImg = (id, b64) => {
      const el = document.getElementById(id);
      if (el) {
        // H-08: solo base64 válido llega al src; cualquier otra cosa es
        // "Imagen no disponible".
        if (b64 && /^[A-Za-z0-9+/=]+$/.test(b64)) el.innerHTML = `<img src="data:image/png;base64,${b64}">`;
        else el.innerHTML = '<div style="color:var(--text-faint); font-size:11px; padding:40px;">Imagen no disponible.</div>';
      }
    };
    setImg('venn-container', data.venn_plot);
    setImg('volcano-container', data.volcano_plot);
    setImg('ppi-container', data.ppi_plot);
    setImg('heatmap-container', data.heatmap_plot);

    // Validación de cada gen: puntaje 0-100 por capas de evidencia
    // independientes + versiones exactas de las bases usadas.
    const vEl = document.getElementById('validation-container');
    if (vEl) {
      const v = data.validation;
      if (v && v.genes && Object.keys(v.genes).length) {
        const s = v.summary || {};
        const badge = document.getElementById('validation-badge');
        if (badge) badge.textContent = `${s.high || 0} alta · ${s.medium || 0} media · ${s.low || 0} baja` + (s.n_not_evaluated ? ` · ${s.n_not_evaluated} sin evaluar` : '') + (s.n_strong !== undefined ? ` · ${s.n_strong} fuerte · ${s.n_weak} débil · ${s.n_only_predicted} solo predicho` : '');
        // v1.54: aviso visible del modo legacy (sin índice oficial).
        let html = (data.warnings || []).map(w => `<div style="border:1px solid var(--gold); border-radius:8px; padding:10px 12px; margin-bottom:12px; font-size:12px; line-height:1.5;">⚠️ ${escapeHtml(w)}</div>`).join('');
        html += `<table class="sci-table"><thead><tr><th>Gen</th><th>Puntaje</th><th>Nivel</th><th>Evidencia miRTarBase</th><th>Dossier (4 niveles)</th><th>Capas de evidencia</th><th>Qué falta</th></tr></thead><tbody>`;
        for (const [gene, gv] of Object.entries(v.genes)) {
          const levelColor = gv.level === 'alta' ? 'var(--teal)' : gv.level === 'media' ? 'var(--gold)' : gv.level === 'not_evaluated' ? 'var(--text-dim)' : 'var(--red)';
          const levelLabel = gv.level === 'not_evaluated' ? 'no evaluado' : gv.level;
          const scoreText = gv.score === null || gv.score === undefined ? '—' : `${gv.score}/100`;
          // Ronda 3 (v1.48): nivel de evidencia miRTarBase por gen
          // (fuerte / apoyo de alto rendimiento / solo predicho) con sus
          // PMIDs. Matiz 2 (v1.51): la evidencia débil no se llama
          // "validado": es apoyo de alto rendimiento.
          const evLabel = {strong: 'fuerte', weak: 'apoyo de alto rendimiento', clip_only: 'solo CLIP sin PMID', legacy: 'validado (snapshot)', none: 'solo predicho'}[gv.evidence_level] || '—';
          const evPmids = (gv.mirtarbase_pmids || []).filter(p => Number.isInteger(p)).map(p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${p}/" target="_blank" rel="noopener">${p}</a>`).join(', ');
          const evCell = `${escapeHtml(evLabel)}${evPmids ? `<br><span style="font-family:var(--font-mono); font-size:10px;">${evPmids}</span>` : ''}`;
          // Ronda 3 (v1.49): dossier de 4 niveles + enfermedades asociadas
          // (Open Targets, o DisGeNET cuando está configurado).
          const tiers = gv.evidence_tiers || {};
          const tierOrder = ['prediction', 'experimental', 'disease', 'clinical_literature'];
          const tierNames = {prediction: 'Predicción', experimental: 'Experimental', disease: 'Enfermedad', clinical_literature: 'Literatura'};
          const tierMarks = tierOrder.map(k => {
            const t = tiers[k];
            const mark = t && t.present ? '●' : '○';
            return `<span title="${escapeHtml(t ? t.summary : '')}">${mark} ${tierNames[k]}</span>`;
          }).join('<br>');
          const disNames = (gv.diseases || []).map(dd => `${dd.name} (${dd.score})`).join('; ');
          const dossierCell = `<b>${gv.tiers_present ?? '—'}/4</b><br><span style="font-size:10px;">${tierMarks}</span>${disNames ? `<br><span style="font-size:10px; color:var(--text-dim);">${escapeHtml(disNames)}</span>` : ''}`;
          html += `<tr><td><b>${escapeHtml(gene)}</b></td><td><span style="font-family:var(--font-mono);">${scoreText}</span></td><td><span style="color:${levelColor}; font-weight:600;">${escapeHtml(levelLabel)}</span></td><td style="font-size:11px;">${evCell}</td><td style="font-size:11px;">${dossierCell}</td><td style="font-size:11px;">${(gv.layers || []).map(escapeHtml).join('<br>') || '—'}</td><td style="font-size:11px; color:var(--text-dim);">${(gv.missing || []).map(escapeHtml).join('<br>') || '—'}</td></tr>`;
        }
        html += '</tbody></table>';
        const ds = v.data_sources || {};
        const dsText = Object.entries(ds).map(([k, val]) => `${k} ${val}`).join(' · ');
        html += `<div style="margin-top:12px; color:var(--text-dim); font-size:11px; font-family:var(--font-mono);">Bases de datos: ${dsText || '—'}</div>`;
        if (data.research_scope && data.research_scope.note) {
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">${data.research_scope.note}</div>`;
        }
        if (s.mean_score !== null && s.mean_score !== undefined) {
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">${s.mean_score_label || 'Media de los genes evaluados'}: <b style="color:var(--text);">${s.mean_score}</b> — es la media solo de los evaluados, no de todo el panel.</div>`;
        }
        if (data.source_by_mirna) {
          const srcText = Object.entries(data.source_by_mirna).map(([m, src]) => `${m} — ${src}`).join('<br>');
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">Fuente por miRNA (un consenso puede mezclar predicción y validación experimental):<br>${srcText}</div>`;
        }
        // H-17 (v1.47): avisar en la interfaz de los miRNAs sin datos y
        // del panel de un solo miRNA (sin convergencia multi-miRNA).
        if (data.missing_mirnas && data.missing_mirnas.length) {
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">miRNAs sin datos en ninguna fuente (no entraron al consenso): ${data.missing_mirnas.map(escapeHtml).join(', ')}</div>`;
        }
        // v1.58: brazos que TargetScan no modela (su semilla no está en
        // la secuencia del brazo pedido): sus dianas se bloquearon, no
        // se sustituyeron por las del otro brazo.
        if (data.arm_unverified && data.arm_unverified.length) {
          const arms = data.arm_unverified.map(e => `${e.mirna} (${escapeHtml(e.reason)})`).join('<br>');
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">Brazos no modelados por TargetScan (sin dianas predichas):<br>${arms}</div>`;
        }
        if (data.arm_verification === 'unavailable') {
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">Verificación de brazos no disponible (faltan los archivos de semillas): los brazos no se comprobaron.</div>`;
        }
        if (data.single_mirna_warning) {
          html += `<div style="margin-top:8px; color:var(--text-dim); font-size:11px;">${escapeHtml(data.single_mirna_warning)}</div>`;
        }
        if (data.enrichment_note) {
          html += `<div style="margin-top:4px; color:var(--text-dim); font-size:11px;">${data.enrichment_note}</div>`;
        }
        vEl.innerHTML = html;
      } else {
        vEl.innerHTML = '<div style="color:var(--text-faint); padding:20px; text-align:center;">Sin dossier de validación en esta respuesta (versión anterior del motor).</div>';
      }
    }

    // Relevancia de las proteínas (UniProt): qué hace cada proteína,
    // dónde trabaja y con qué enfermedad se asocia.
    const prEl = document.getElementById('proteins-container');
    if (prEl) {
      const gd = data.gene_details || {};
      const rows = (data.common_genes || []).map(g => {
        const d = gd[g] || {};
        const p = d.protein;
        if (!p) return `<tr><td><b>${escapeHtml(g)}</b></td><td colspan="3" style="color:var(--text-faint); font-size:11px;">Sin datos de UniProt (sin internet o sin entrada revisada).</td></tr>`;
        const func = (p.function || '—').slice(0, 220);
        const dis = (p.diseases || []).join('; ') || '—';
        const loc = (p.locations || []).join(', ') || '—';
        // H-08: el enlace solo se renderiza si es una URL https de UniProt.
        const link = (p.url && p.url.startsWith('https://rest.uniprot.org/'))
          ? `<br><a href="${escapeHtml(p.url)}" target="_blank" style="color:var(--gold); font-size:10px;">UniProt ${escapeHtml(p.uniprot_id || '')}</a>` : '';
        return `<tr><td><b>${escapeHtml(g)}</b><br><span style="color:var(--text-dim); font-size:11px;">${escapeHtml(p.protein_name || '')}</span>${link}</td><td style="font-size:11px;">${escapeHtml(func)}</td><td style="font-size:11px;">${escapeHtml(loc)}</td><td style="font-size:11px;">${escapeHtml(dis)}</td></tr>`;
      }).join('');
      prEl.innerHTML = rows
        ? `<table class="sci-table"><thead><tr><th>Gen / Proteína</th><th>Qué hace</th><th>Dónde trabaja</th><th>Enfermedad asociada</th></tr></thead><tbody>${rows}</tbody></table>`
        : '<div style="color:var(--text-faint); padding:20px; text-align:center;">Sin genes núcleo.</div>';
    }

    const enEl = document.getElementById('enrich-container');
    if (enEl) {
      if (data.enrichment && data.enrichment.length) {
        let html = `<table class="sci-table"><thead><tr><th>Ruta Biológica</th><th>Fuente</th><th>p-valor</th><th title="${SPEC_TITLE}">Especificidad</th><th>PubMed</th></tr></thead><tbody>`;
        data.enrichment.filter(it => it.is_group_representative !== false).forEach(item => {
          const pvalVal = (item.Pval !== undefined && item.Pval !== null) 
            ? (typeof item.Pval === 'number' ? item.Pval.toExponential(2) : item.Pval)
            : '—';
          
          // MÁXIMA ROBUSTEZ: Detectar ID en objeto, array o string directo
          let pmid = null;
          if (item.Evidence) {
            if (Array.isArray(item.Evidence) && item.Evidence.length > 0) pmid = item.Evidence[0].id || item.Evidence[0];
            else pmid = item.Evidence.id || item.Evidence;
          }
          // H-08: el PMID debe ser numérico; cualquier otra cosa no enlaza.
          const pmidOk = pmid !== null && pmid !== undefined && /^[0-9]+$/.test(String(pmid));
          const pLink = pmidOk ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${pmid}" target="_blank" style="color:var(--gold); font-family:var(--font-mono); font-size:10px; font-weight:600; text-decoration:none; border-bottom:1px solid var(--gold-dim);">PMID: ${pmid}</a>` : '<span style="color:var(--text-faint); font-size:10px;">Sin evidencia</span>';

          // v1.54: especificidad — fracción de listas aleatorias de genes
          // del mismo tamaño, extraídas del universo de TargetScan, en las
          // que el término aparece por azar.
          // Es un número del servidor, no texto libre: no necesita escape.
          // Unreleased: con 100 muestras la resolución es del 1%: un 0
          // observado se muestra como "<1%", nunca como "0%".
          const gm = Array.isArray(item.group_members) ? item.group_members : [];
          const groupBadge = gm.length
            ? `<br><small style="color:var(--text-faint);" title="${escapeHtml(gm.join('; '))}">+${gm.length} término${gm.length === 1 ? '' : 's'} similar${gm.length === 1 ? '' : 'es'} (mismos genes del panel)</small>`
            : '';
          const spec = fmtSpec(item.specificity);
          const specTitle = SPEC_TITLE;

          html += `<tr><td><span style="font-weight:500;">${escapeHtml(item.Term || '—')}</span><br><small style="color:var(--text-dim);">${escapeHtml(item.ScientificDesc || '')}</small>${groupBadge}</td>
            <td><span class="source-tag">${escapeHtml(item.Source || '—')}</span></td>
            <td><span style="font-family:var(--font-mono); font-size:11px;">${pvalVal}</span></td>
            <td><span style="font-family:var(--font-mono); font-size:11px;" title="${specTitle}">${spec}</span></td>
            <td>${pLink}</td></tr>`;
        });
        html += '</tbody></table>';
        enEl.innerHTML = html;
      } else {
        enEl.innerHTML = '<div style="color:var(--text-faint); padding:20px; text-align:center;">Sin resultados de enriquecimiento funcional significativo en la ruta local.</div>';
      }
    }

    // Unreleased (2026-11-19): la ruta viva de Enrichr es una
    // comprobación exploratoria, en su propia sección. Su p ajustado no
    // es un FDR estándar con listas pequeñas, así que la columna se
    // rotula como exploratoria y nunca como significancia.
    const enrEl = document.getElementById('enrich-container-enrichr');
    if (enrEl) {
      if (data.enrichment_enrichr && data.enrichment_enrichr.length) {
        let html = `<table class="sci-table"><thead><tr><th>Ruta Biológica</th><th>Fuente</th><th>p ajustado de Enrichr (exploratorio)</th><th title="${SPEC_TITLE}">Especificidad</th></tr></thead><tbody>`;
        data.enrichment_enrichr.filter(it => it.is_group_representative !== false).forEach(item => {
          const pvalVal = (item.Pval !== undefined && item.Pval !== null)
            ? (typeof item.Pval === 'number' ? item.Pval.toExponential(2) : item.Pval)
            : '—';
          const gm = Array.isArray(item.group_members) ? item.group_members : [];
          const groupBadge = gm.length
            ? `<br><small style="color:var(--text-faint);" title="${escapeHtml(gm.join('; '))}">+${gm.length} término${gm.length === 1 ? '' : 's'} similar${gm.length === 1 ? '' : 'es'} (mismos genes del panel)</small>`
            : '';
          const spec = fmtSpec(item.specificity);
          const specTitle = SPEC_TITLE;
          html += `<tr><td><span style="font-weight:500;">${escapeHtml(item.Term || '—')}</span><br><small style="color:var(--text-dim);">${escapeHtml(item.ScientificDesc || '')}</small>${groupBadge}</td>
            <td><span class="source-tag">${escapeHtml(item.Source || '—')}</span></td>
            <td><span style="font-family:var(--font-mono); font-size:11px;">${pvalVal}</span></td>
            <td><span style="font-family:var(--font-mono); font-size:11px;" title="${specTitle}">${spec}</span></td></tr>`;
        });
        html += '</tbody></table>';
        enrEl.innerHTML = html;
      } else {
        enrEl.innerHTML = '<div style="color:var(--text-faint); padding:20px; text-align:center;">Sin términos de Enrichr en esta respuesta.</div>';
      }
    }
  } catch (err) {
    console.error("Error renderizando dashboard:", err);
  }
}

function panelFamily(m) {
  // hsa-miR-33a-5p -> mir-33-5p ; hsa-miR-106b-5p -> mir-106-5p ;
  // hsa-let-7a-5p -> let-7-5p (las familias TargetScan se nombran por semilla).
  let f = (m || '').trim().toLowerCase().replace(/^hsa-/, '');
  f = f.replace(/^((?:mir|let)-)(\d+)[a-z](?=-[53]p)/, '$1$2');
  return f;
}

function currentPanelFamilies() {
  const input = document.getElementById('mirna-input').value || '';
  return input.split(',').map(panelFamily).filter(Boolean);
}

const EXAMPLES = {
  lipidos: 'hsa-miR-33a-5p, hsa-miR-144-3p, hsa-miR-758-3p, hsa-miR-106b-5p, hsa-miR-33b-5p',
  cancer: 'hsa-miR-21-5p, hsa-miR-155-5p, hsa-let-7a-5p',
  control: 'hsa-miR-375, hsa-miR-29a-3p, hsa-miR-21-5p'
};

function loadExample(key) {
  const panel = EXAMPLES[key];
  if (!panel) return;
  document.getElementById('mirna-input').value = panel;
  analyzePro();
}

async function reverseLookup() {
  const gene = (document.getElementById('reverse-gene').value || '').trim();
  const out = document.getElementById('reverse-results');
  if (!gene) {
    out.innerHTML = '<div style="color:var(--text-dim); font-size:12px;">Escribe un gen, por ejemplo ABCA1.</div>';
    return;
  }
  out.innerHTML = '<div style="color:var(--text-dim); font-size:12px;">Buscando…</div>';
  try {
    const res = await fetch(`${API_URL}/api/v1/reverse?gene=${encodeURIComponent(gene)}`);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    const data = await res.json();
    const panelFams = currentPanelFamilies();
    const pills = data.mirnas.map(m => {
      const inPanel = panelFams.some(f => m.split('/').includes(f));
      const style = inPanel ? 'background:var(--gold-dim); border-color:var(--gold); color:var(--gold);' : '';
      return `<span class="gene-pill" style="${style}">${escapeHtml(m)}</span>`;
    }).join('');
    out.innerHTML = `<div style="margin-bottom:8px; font-size:13px;"><b>${escapeHtml(data.gene)}</b> es diana de <b>${data.n}</b> familias de miRNAs.</div>` +
      (data.n
        ? `<div style="color:var(--text-dim); font-size:11px; margin-bottom:6px;">Las doradas son de tu panel actual.</div><div class="genes-wrap">${pills}</div>`
        : '<div style="color:var(--text-dim); font-size:12px;">Ningún miRNA de la base lo regula (o el gen no está en TargetScan).</div>');
  } catch (e) {
    out.innerHTML = `<div style="color:var(--red); font-size:12px;">Error: ${escapeHtml(e.message)}</div>`;
  }
}

async function openGenePanel(gene) {
  document.getElementById('gp-name').textContent = gene;
  document.getElementById('gene-panel').classList.add('open');
  const d = lastData?.gene_details?.[gene] || {full_name: gene, system:'Multisistémico', pathology:'Diana de alta confianza.', associated_routes:[]};
  // Puntaje REAL del dossier de validación (no un número fijo).
  // H-18 (v1.47): sin puntaje no hay barra — "Sin evaluar", nunca un 85 fijo.
  const vScore = lastData?.validation?.genes?.[gene]?.score;
  const score = (vScore !== undefined && vScore !== null) ? vScore : null;
  // Dossier de 4 niveles (v1.49): predicción, experimental, enfermedad,
  // literatura/datos clínicos + enfermedades asociadas (Open Targets).
  const gv = lastData?.validation?.genes?.[gene] || {};
  const tiers = gv.evidence_tiers || {};
  const tierOrder = [['prediction', 'Predicción'], ['experimental', 'Experimental'], ['disease', 'Enfermedad'], ['clinical_literature', 'Literatura y datos clínicos']];
  const tierRows = tierOrder.map(([k, name]) => {
    const t = tiers[k];
    const mark = t && t.present ? '●' : '○';
    return `<div style="font-size:12px; color:var(--text-dim); line-height:1.7;">${mark} <b style="color:var(--text);">${name}</b> — ${escapeHtml(t ? t.summary : 'Sin datos')}</div>`;
  }).join('');
  const disRows = (gv.diseases || []).map(dd => `<div style="font-size:12px; color:var(--text-dim); line-height:1.7;">${escapeHtml(dd.name)} <span style="font-family:var(--font-mono); font-size:11px;">(${escapeHtml(String(dd.score))})</span></div>`).join('');
  // Tabla miRNA x gen (v1.51, idea 1): qué miRNA aporta cada nivel de
  // evidencia para este gen, con sus PMIDs y experimentos. "Fuerte" es
  // por pareja: no dice que todos los miRNAs del panel lo regulen.
  const evTable = gv.evidence_by_mirna || {};
  const panelMirnas = (lastData?.validation?.mirnas || Object.keys(evTable));
  const evLevelLabel = {strong: 'fuerte', weak: 'apoyo de alto rendimiento', clip_only: 'solo CLIP sin PMID (sin puntos)', legacy: 'validado (snapshot)', none: 'solo predicho'};
  const evRows = panelMirnas.map(m => {
    const e = evTable[m] || {level: 'none', pmids: [], experiments: []};
    const pmids = (e.pmids || []).filter(p => Number.isInteger(p)).map(p => `<a href="https://pubmed.ncbi.nlm.nih.gov/${p}/" target="_blank" rel="noopener">${p}</a>`).join(', ');
    const exps = (e.experiments || []).map(escapeHtml).join('; ');
    const bare = e.bare_name ? ' <span style="color:var(--text-dim);">(brazo no especificado)</span>' : '';
    return `<tr><td style="font-family:var(--font-mono); font-size:11px;">${escapeHtml(m)}${bare}</td><td style="font-size:11px;">${escapeHtml(evLevelLabel[e.level] || '—')}</td><td style="font-family:var(--font-mono); font-size:10px;">${pmids || '—'}</td><td style="font-size:11px;">${exps || '—'}</td></tr>`;
  }).join('');
  const evBlock = evRows ? `<div style="margin-bottom:20px; padding:15px; background:rgba(79,195,161,0.03); border-radius:8px; border-left:4px solid var(--teal);">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase; display:block; margin-bottom:8px;">Evidencia por miRNA — qué miRNA valida este gen</small>
      <table class="sci-table" style="width:100%;"><thead><tr><th style="text-align:left;">miRNA</th><th style="text-align:left;">Nivel</th><th style="text-align:left;">PMIDs</th><th style="text-align:left;">Experimentos</th></tr></thead><tbody>${evRows}</tbody></table>
      <div style="font-size:11px; color:var(--text-dim); margin-top:8px;">«Fuerte» significa que ese miRNA tiene evidencia fuerte para este gen; no dice que todos los miRNAs del panel lo regulen. La evidencia débil es apoyo de alto rendimiento y no demuestra unión directa. «Solo CLIP sin PMID» es evidencia débil sin referencia citable: se muestra, no suma puntos y nunca desbloquea «alta». «Brazo no especificado» significa que la evidencia viene de una entrada de miRTarBase sin -5p/-3p, así que no dice de qué brazo es.</div>
    </div>` : '';
  
  document.getElementById('gene-panel-body').innerHTML = `
    <div style="margin-bottom:20px;">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase;">Identificación Genómica</small>
      <div style="font-family:var(--font-serif); font-size:16px; color:var(--text-main); margin-top:4px;">${escapeHtml(d.full_name)}</div>
    </div>
    <div style="margin-bottom:20px;">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase;">Sistema Fisiológico</small>
      <div style="color:var(--gold); font-weight:600; font-size:13.5px; margin-top:4px;">${escapeHtml(d.system)}</div>
    </div>
    <div style="margin-bottom:20px;">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase;">Rigor Científico</small>
      ${score === null
        ? `<div style="margin-top:8px; font-size:12px; color:var(--text-dim);">Sin evaluar</div>`
        : `<div style="display:flex; align-items:center; gap:12px; margin-top:8px;">
        <div style="flex:1; height:7px; background:rgba(255,255,255,0.06); border-radius:10px; overflow:hidden; border:1px solid rgba(255,255,255,0.04);">
          <div style="width:${score}%; height:100%; background:linear-gradient(90deg, #c8a96e, #f2d295); border-radius:10px; box-shadow:0 0 10px rgba(200,169,110,0.3);"></div>
        </div>
        <span style="font-family:var(--font-mono); font-size:11px; color:var(--gold); font-weight:800;">${score}%</span>
      </div>`}
    </div>
    <div style="margin-bottom:20px; padding:15px; background:rgba(200,169,110,0.03); border-radius:8px; border-left:4px solid var(--gold);">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase; display:block; margin-bottom:8px;">Relevancia Patológica Detallada</small>
      <div style="font-family:'Spectral', serif; font-size:13.5px; text-align:justify; color:var(--text); line-height:1.7;">${escapeHtml(d.pathology)}</div>
    </div>
    ${Object.keys(tiers).length ? `<div style="margin-bottom:20px; padding:15px; background:rgba(79,195,161,0.03); border-radius:8px; border-left:4px solid var(--teal);">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase; display:block; margin-bottom:8px;">Dossier de evidencia — ${gv.tiers_present ?? '—'}/4 niveles</small>
      ${tierRows}
      ${disRows ? `<div style="margin-top:8px; font-size:11px; color:var(--text-dim);">Enfermedades asociadas (Open Targets):</div>${disRows}` : ''}
    </div>` : ''}
    ${evBlock}
    ${d.protein ? `<div style="margin-bottom:20px; padding:15px; background:rgba(79,195,161,0.05); border-radius:8px; border-left:4px solid var(--teal);">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase; display:block; margin-bottom:8px;">Proteína (UniProt)</small>
      <div style="font-size:13.5px; font-weight:600; color:var(--text);">${escapeHtml(d.protein.protein_name)} ${(d.protein.url && d.protein.url.startsWith('https://rest.uniprot.org/')) ? `<a href="${escapeHtml(d.protein.url)}" target="_blank" style="color:var(--gold); font-size:11px; font-family:var(--font-mono);">${escapeHtml(d.protein.uniprot_id || '')}</a>` : escapeHtml(d.protein.uniprot_id || '')}</div>
      <div style="font-size:12px; color:var(--text-dim); line-height:1.6; margin-top:6px;">${escapeHtml(d.protein.function || 'Sin descripción en UniProt.')}</div>
      ${d.protein.locations && d.protein.locations.length ? `<div style="font-size:11px; color:var(--text-dim); margin-top:6px;">Dónde: ${d.protein.locations.map(escapeHtml).join(', ')}</div>` : ''}
      ${d.protein.diseases && d.protein.diseases.length ? `<div style="font-size:11px; color:var(--text-dim); margin-top:4px;">Enfermedad: ${d.protein.diseases.map(escapeHtml).join('; ')}</div>` : ''}
    </div>` : ''}
    <div style="margin-bottom:10px;">
      <small style="color:var(--text-faint); font-family:var(--font-mono); font-size:9px; letter-spacing:1.5px; text-transform:uppercase;">Rutas Metabólicas de Consenso</small>
      <div style="font-size:11.5px; color:var(--text-dim); line-height:1.7; margin-top:8px; font-style:italic;">
        ${(d.associated_routes||[]).map(r => `<i class="fas fa-microscope" style="font-size:10px; margin-right:8px; color:var(--teal);"></i>${escapeHtml(r)}`).join('<br>')}
      </div>
    </div>
    ${d.pmid && /^[0-9]+$/.test(String(d.pmid)) ? `<div style="margin-top:25px; padding-top:15px; border-top:1px solid var(--border);"><a href="https://pubmed.ncbi.nlm.nih.gov/${d.pmid}" target="_blank" style="color:var(--teal); font-family:var(--font-mono); font-size:11px; text-decoration:none; display:flex; align-items:center; gap:8px; font-weight:600;"><i class="fas fa-book-medical"></i> VERIFICAR EVIDENCIA (PMID: ${d.pmid})</a></div>` : ''}
  `;
}
function closeGenePanel() { document.getElementById('gene-panel').classList.remove('open'); }

/* ========== LÓGICA DEL EDITOR PROFESIONAL ========== */

function romanize(n) {
  const map = [[1000,'M'],[900,'CM'],[500,'D'],[400,'CD'],[100,'C'],[90,'XC'],[50,'L'],[40,'XL'],[10,'X'],[9,'IX'],[5,'V'],[4,'IV'],[1,'I']];
  return map.reduce((acc,[v,r]) => { while (n >= v) { acc += r; n -= v; } return acc; }, '');
}

function fmt(cmd) { document.execCommand(cmd, false, null); }
function setFontSize(size) {
  document.execCommand('fontSize', false, '7');
  document.querySelectorAll('font[size="7"]').forEach(el => { el.removeAttribute('size'); el.style.fontSize = size; });
}

function setAccent(color, light, ev) {
  document.documentElement.style.setProperty('--accent', color);
  document.documentElement.style.setProperty('--accent-light', light);
  document.querySelectorAll('.color-chip').forEach(c => c.classList.remove('selected'));
  if (ev && ev.target) ev.target.classList.add('selected');
}

function moveSection(id, dir) {
  const el = document.getElementById(id);
  if (!el) return;
  if (dir === 'up' && el.previousElementSibling && el.previousElementSibling.classList.contains('report-section')) {
    el.parentNode.insertBefore(el, el.previousElementSibling);
  } else if (dir === 'down' && el.nextElementSibling && el.nextElementSibling.classList.contains('report-section')) {
    el.parentNode.insertBefore(el.nextElementSibling, el);
  }
  updateOutline();
}

window.addEventListener('beforeunload', (e) => {
  if (lastData) { e.preventDefault(); e.returnValue = ''; }
});

function setPaperColor(color) { document.documentElement.style.setProperty('--paper', color); }
function setInkColor(color) { document.documentElement.style.setProperty('--ink', color); }

function switchPanel(name) {
  ['outline','insert','props'].forEach(p => {
    document.getElementById('panel-'+p).style.display = (p === name ? 'block' : 'none');
    document.getElementById('tab-'+p).classList.toggle('active', p === name);
  });
  if (name === 'outline') updateOutline();
}

function updateMeta() {
  const inst = document.getElementById('prop-inst').value;
  const ver = document.getElementById('prop-ver').value;
  const dateElem = document.getElementById('prop-date').value;
  document.querySelectorAll('.meta-inst-val').forEach(el => el.innerText = inst);
  document.querySelectorAll('.meta-ver-val').forEach(el => el.innerText = 'v' + ver);
  if (dateElem) {
    const d = new Date(dateElem + 'T12:00:00');
    const dateStr = d.toLocaleDateString('es-ES', { month: 'long', year: 'numeric' });
    document.querySelectorAll('.meta-date-val').forEach(el => el.innerText = dateStr);
  }
}

function updateOutline() {
  const list = document.getElementById('outline-list');
  const sections = document.querySelectorAll('.report-section');
  if (!list) return;
  list.innerHTML = '';
  sections.forEach((sec, i) => {
    const heading = sec.querySelector('.section-heading');
    let text = heading ? heading.innerText.replace(/^\w+\.\s*/, '').trim() : ('Sección '+(i+1));
    const div = document.createElement('div');
    div.className = 'outline-item';
    div.innerHTML = `<span style="margin-right:6px;color:#888;">${romanize(i+1)}</span> ${escapeHtml(text)}`;
    div.onclick = () => sec.scrollIntoView({ behavior: 'smooth', block: 'start' });
    list.appendChild(div);
  });
  const allText = document.getElementById('report-canvas-content').innerText || '';
  const words = allText.trim().split(/\s+/).filter(w=>w.length>0).length;
  const chars = allText.replace(/\s/g,'').length;
  document.getElementById('doc-stats').innerHTML = `Páginas: ${pageCount}<br>Palabras: ${words}<br>Caracteres: ${chars}`;
}

function insertBlock(type) {
  const activePage = document.querySelector('.a4-page:last-of-type .page-inner');
  if (!activePage) return;
  sectionCount++;
  const id = `sec-${sectionCount}`;
  let html = '';
  if (type === 'note') {
    html = `<div class="report-section note-block" id="${id}"><div class="editable-block" contenteditable="true" data-placeholder="Escriba su nota aquí..."></div><div class="section-actions"><button class="section-act-btn" onclick="moveSection('${id}','up')">↑</button><button class="section-act-btn" onclick="moveSection('${id}','down')">↓</button><button class="section-act-btn danger" onclick="removeSection('${id}')">✕</button></div></div>`;
  } else if (type === 'quote') {
    html = `<div class="report-section" id="${id}"><div class="quote-block editable-block" contenteditable="true">Texto de la cita...</div><div class="section-actions"><button class="section-act-btn danger" onclick="removeSection('${id}')">✕</button></div></div>`;
  } else if (type === 'gene') {
    html = `<div class="report-section gene-block" id="${id}"><div class="gene-block-name" contenteditable="true">NOMBRE_GEN</div><div class="gene-block-text" contenteditable="true">Descripción...</div><div class="section-actions"><button class="section-act-btn danger" onclick="removeSection('${id}')">✕</button></div></div>`;
  } else if (type === 'divider') {
    html = `<div class="report-section" id="${id}" style="margin:10px 0;"><hr style="border-top:1px solid var(--border-sci);"><div class="section-actions"><button class="section-act-btn danger" onclick="removeSection('${id}')">✕</button></div></div>`;
  }
  const addBtn = activePage.querySelector('.add-section-btn');
  if (addBtn) addBtn.insertAdjacentHTML('beforebegin', html);
  updateOutline();
}

function addTableRow() {
  const tbody = document.getElementById('bio-tbody');
  if (tbody) tbody.insertAdjacentHTML('beforeend', `<tr><td style="font-weight:600;">NUEVO</td><td>Sistema</td><td style="font-style:italic;">Ruta</td><td style="text-align:justify;">Relevancia</td></tr>`);
}

function addPage() {
  pageCount++;
  const canvas = document.getElementById('report-canvas-content');
  const div = document.createElement('div');
  div.className = 'a4-page';
  div.id = `page-${pageCount}`;
  div.innerHTML = `<div class="page-inner"><div class="report-header"><div><div class="report-title">INFORME</div><div class="report-subtitle">CONTINUACIÓN</div></div><div class="report-meta"><div class="meta-date-val"></div><div class="meta-inst-val"></div><div class="meta-ver-val"></div></div></div><button class="add-section-btn" onclick="addNewSection()">+ Agregar nueva sección</button></div><div class="page-footer"><span>KENRYU Bioinformatics Engine</span><span class="page-num">${pageCount}</span></div>`;
  canvas.appendChild(div);
  updateMeta();
  updateOutline();
}

function addNewSection() {
  sectionCount++;
  const id = `sec-${sectionCount}`;
  const html = `<div class="report-section" id="${id}"><div class="section-heading" contenteditable="true"><span class="s-num">${romanize(sectionCount)}.</span> Nueva sección</div><div class="editable-block" contenteditable="true">Escriba el contenido...</div><div class="section-actions"><button class="section-act-btn" onclick="moveSection('${id}','up')">↑</button><button class="section-act-btn" onclick="moveSection('${id}','down')">↓</button><button class="section-act-btn danger" onclick="removeSection('${id}')">✕</button></div></div>`;
  event.target.insertAdjacentHTML('beforebegin', html);
  updateOutline();
}

function removeSection(id) { if (confirm('¿Eliminar?')) { document.getElementById(id).remove(); updateOutline(); } }
function togglePreview() { isPreview = !isPreview; document.body.classList.toggle('preview-mode', isPreview); document.getElementById('preview-btn').textContent = isPreview ? '✏️ Editar' : '👁 Vista previa'; }
function exportText() { const blob = new Blob([document.getElementById('report-canvas-content').innerText], {type:'text/plain'}); const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download='Reporte.txt'; a.click(); }

// ── EXPORTACIÓN MARKDOWN PROFESIONAL ──────────────────────────────────────────
// Genera un ZIP con:
//   /Reporte.md      ← markdown estructurado (encabezados, tablas, enlaces)
//   /assets/         ← imágenes PNG extraídas (Venn, Volcano, PPI)
// Estructura compatible con VSCode, Obsidian, Pandoc, GitHub. Lo que un
// bioinformático espera de un reporte reproducible.

function _b64ToBlob(b64Data, contentType) {
  const byteCharacters = atob(b64Data);
  const byteArrays = [];
  for (let offset = 0; offset < byteCharacters.length; offset += 512) {
    const slice = byteCharacters.slice(offset, offset + 512);
    const byteNumbers = new Array(slice.length);
    for (let i = 0; i < slice.length; i++) byteNumbers[i] = slice.charCodeAt(i);
    byteArrays.push(new Uint8Array(byteNumbers));
  }
  return new Blob(byteArrays, {type: contentType});
}

function _htmlToMarkdown(root, assetCounter) {
  // Walker recursivo que convierte HTML del reporte a Markdown.
  // assetCounter es un objeto compartido {n, images: [{name, blob}]} para
  // poder extraer imágenes y referenciarlas por nombre relativo.
  if (!root) return '';
  const out = [];

  const walk = (node) => {
    if (node.nodeType === Node.TEXT_NODE) {
      return node.textContent;
    }
    if (node.nodeType !== Node.ELEMENT_NODE) return '';

    const tag = node.tagName.toLowerCase();
    const cls = node.className || '';

    // Saltar elementos ocultos, decorativos o de UI
    if (node.style && (node.style.display === 'none' || node.style.visibility === 'hidden')) return '';
    if (cls.includes('section-actions') || cls.includes('add-section-btn') ||
        cls.includes('page-footer') || cls.includes('report-header') ||
        tag === 'script' || tag === 'style' || tag === 'button') return '';

    // Encabezados del editor: .section-heading = ## (sección II nivel)
    if (cls.includes('section-heading')) {
      // Normalizar whitespace: el span .s-num puede tener saltos internos
      const txt = node.innerText.replace(/\s+/g, ' ').trim();
      return '\n## ' + txt + '\n\n';
    }

    // Imágenes: extraer base64 y crear referencia relativa al ZIP
    if (tag === 'img') {
      const src = node.src || '';
      if (src.startsWith('data:image/')) {
        const match = src.match(/^data:image\/(png|jpe?g|svg\+xml|webp);base64,(.+)$/);
        if (match) {
          const ext = match[1] === 'jpeg' ? 'jpg' : (match[1] === 'svg+xml' ? 'svg' : match[1]);
          assetCounter.n++;
          // Nombre del asset: prefer alt si existe, sino 'figura-N'
          const baseName = (node.alt && node.alt.trim()) ? node.alt.trim().replace(/[^a-zA-Z0-9_-]/g, '-') : ('figura-' + assetCounter.n);
          const filename = baseName + '.' + ext;
          const contentType = match[1] === 'svg+xml' ? 'image/svg+xml' : 'image/' + match[1];
          assetCounter.images.push({
            name: filename,
            blob: _b64ToBlob(match[2], contentType),
          });
          const altText = node.alt || ('Figura ' + assetCounter.n);
          return '\n![' + altText + '](assets/' + filename + ')\n\n';
        }
      } else if (src) {
        // URL externa: dejar referencia directa
        return '\n![' + (node.alt || 'imagen') + '](' + src + ')\n\n';
      }
      return '';
    }

    // Tablas: conversión completa a markdown tipo GFM
    if (tag === 'table') {
      const rows = Array.from(node.querySelectorAll('tr'));
      if (rows.length === 0) return '';
      const md = [];
      let isFirst = true;
      for (const row of rows) {
        const cells = Array.from(row.children).map(c => c.innerText.trim().replace(/\|/g, '\\|').replace(/\n+/g, ' '));
        if (cells.length === 0) continue;
        md.push('| ' + cells.join(' | ') + ' |');
        if (isFirst) {
          md.push('| ' + cells.map(() => '---').join(' | ') + ' |');
          isFirst = false;
        }
      }
      return '\n' + md.join('\n') + '\n\n';
    }

    // Enlaces: [texto](url)
    if (tag === 'a' && node.href) {
      const inner = Array.from(node.childNodes).map(walk).join('').trim();
      return '[' + (inner || node.href) + '](' + node.href + ')';
    }

    // Bold / italic / inline code
    if (tag === 'strong' || tag === 'b') return '**' + Array.from(node.childNodes).map(walk).join('') + '**';
    if (tag === 'em' || tag === 'i') {
      // i sin clase fa- es itálica real, ignoramos los íconos FontAwesome
      if (cls.includes('fa') || cls.includes('fab') || cls.includes('fas') || cls.includes('far')) return '';
      return '*' + Array.from(node.childNodes).map(walk).join('') + '*';
    }
    if (tag === 'code') return '`' + node.innerText + '`';

    // Encabezados nativos h1-h6
    if (/^h[1-6]$/.test(tag)) {
      const level = parseInt(tag[1]);
      return '\n' + '#'.repeat(level) + ' ' + node.innerText.trim() + '\n\n';
    }

    // Listas
    if (tag === 'ul' || tag === 'ol') {
      const items = Array.from(node.children).filter(c => c.tagName.toLowerCase() === 'li');
      const lines = items.map((li, i) => {
        const prefix = tag === 'ol' ? (i + 1) + '. ' : '- ';
        const content = Array.from(li.childNodes).map(walk).join('').trim();
        return prefix + content;
      });
      return '\n' + lines.join('\n') + '\n\n';
    }

    // Saltos de página (separadores)
    if (cls.includes('a4-page')) {
      const content = Array.from(node.childNodes).map(walk).join('');
      return content + '\n\n---\n\n';
    }

    // Párrafos y divs: contenido + doble salto
    if (tag === 'p' || tag === 'div') {
      const content = Array.from(node.childNodes).map(walk).join('');
      // Si el contenido es solo whitespace, ignorar
      if (!content.trim()) return '';
      // Mantener saltos visuales como párrafos
      const isBlock = (tag === 'p') || (cls.includes('editable-block') || cls.includes('report-section'));
      return content + (isBlock ? '\n\n' : '');
    }

    // br → salto de línea explícito
    if (tag === 'br') return '\n';

    // Span e inline: solo contenido
    return Array.from(node.childNodes).map(walk).join('');
  };

  return walk(root);
}

async function exportMarkdown() {
  const rc = document.getElementById('report-canvas-content');
  if (!rc) { alert('No hay reporte para exportar. Ejecute primero un análisis.'); return; }
  if (typeof JSZip === 'undefined') { alert('Librería de empaquetado no cargada. Recargue la página.'); return; }

  const assetCounter = { n: 0, images: [] };
  let md = _htmlToMarkdown(rc, assetCounter);

  // Limpieza: colapsar saltos de línea excesivos
  md = md.replace(/\n{3,}/g, '\n\n').trim();

  // Frontmatter YAML estilo Pandoc/Jupyter para metadatos del reporte
  const today = new Date().toISOString().split('T')[0];
  const mirnasInput = document.getElementById('mirna-input');
  const mirnasList = mirnasInput ? mirnasInput.value.split(',').map(s => s.trim()).filter(Boolean) : [];
  const frontmatter = [
    '---',
    'title: "Reporte Bioinformático KENRYU"',
    'subtitle: "Convergencia Molecular y Silenciamiento Génico"',
    'date: "' + today + '"',
    // H-21 (v1.47): la versión viene del backend (lastData.engine_version),
    // nunca de un literal fijo; los miRNAs se serializan con JSON.stringify
    // para que una comilla o un salto de línea no rompa el YAML.
    'generator: "KENRYU Bioinformatics Engine v' + (lastData?.engine_version || '') + '"',
    mirnasList.length ? 'mirnas:\n' + mirnasList.map(m => '  - ' + JSON.stringify(m)).join('\n') : '',
    'lang: es',
    '---',
    '',
  ].filter(Boolean).join('\n');

  const finalMd = frontmatter + '\n' + md + '\n';

  // Si no hay imágenes, descargar .md directo
  if (assetCounter.images.length === 0) {
    const blob = new Blob([finalMd], { type: 'text/markdown;charset=utf-8' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = 'Reporte-KENRYU.md';
    a.click();
    URL.revokeObjectURL(a.href);
    return;
  }

  // Con imágenes: crear ZIP con estructura /Reporte.md + /assets/
  const zip = new JSZip();
  zip.file('Reporte.md', finalMd);
  const assetsFolder = zip.folder('assets');
  for (const img of assetCounter.images) {
    assetsFolder.file(img.name, img.blob);
  }
  // README explicativo
  zip.file('README.md',
    '# Reporte KENRYU — Paquete de exportación Markdown\n\n' +
    'Este ZIP contiene:\n\n' +
    '- **Reporte.md** — Reporte bioinformático completo en formato Markdown estándar (compatible con VSCode, Obsidian, GitHub, Pandoc).\n' +
    '- **assets/** — Carpeta con las figuras del reporte en formato PNG (' + assetCounter.images.length + ' imagen(es)).\n\n' +
    '## Uso\n\n' +
    '1. Para visualizar: abra `Reporte.md` con cualquier editor markdown.\n' +
    '2. Para convertir a PDF/Word: ejecute `pandoc Reporte.md -o Reporte.pdf` (requiere LaTeX).\n' +
    '3. Para publicar en GitHub: copie el contenido tal cual; las imágenes se renderizan automáticamente.\n\n' +
    'Generado por KENRYU Bioinformatics Engine — ' + today + '\n'
  );

  const zipBlob = await zip.generateAsync({ type: 'blob', compression: 'DEFLATE', compressionOptions: { level: 6 } });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(zipBlob);
  a.download = 'Reporte-KENRYU.zip';
  a.click();
  URL.revokeObjectURL(a.href);
}

function createNewPage() {
  pageCount++;
  const page = document.createElement('div');
  page.className = 'a4-page';
  page.id = `page-${pageCount}`;
  const inner = document.createElement('div');
  inner.className = 'page-inner';
  page.appendChild(inner);
  page.insertAdjacentHTML('beforeend', `<div class="page-footer"><span>KENRYU Bioinformatics Engine</span><span class="page-num">${pageCount}</span></div>`);
  return page;
}

function paginateReport(sourceElement) {
  const canvas = document.getElementById('report-canvas-content');
  canvas.innerHTML = '';
  pageCount = 0;

  // Página temporal para medir
  const testPage = document.createElement('div');
  testPage.className = 'a4-page';
  testPage.style.position = 'absolute';
  testPage.style.visibility = 'hidden';
  testPage.style.top = '-9999px';
  document.body.appendChild(testPage);
  const testInner = document.createElement('div');
  testInner.className = 'page-inner';
  testPage.appendChild(testInner);

  const maxHeight = 820; // Altura útil reducida para evitar solapamiento con el pie de página

  const children = Array.from(sourceElement.children);
  let currentPage = createNewPage();
  let currentInner = currentPage.querySelector('.page-inner');
  canvas.appendChild(currentPage);
  let currentHeight = 0;

  children.forEach(child => {
    // SOPORTE PARA SALTOS DE PÁGINA FORZADOS
    if (child.classList.contains('force-page-break')) {
      if (currentHeight > 0) {
        currentPage = createNewPage();
        currentInner = currentPage.querySelector('.page-inner');
        canvas.appendChild(currentPage);
        currentHeight = 0;
      }
      if (child.innerHTML.trim() === "") return; // Si es solo un marcador de salto
    }

    const clone = child.cloneNode(true);
    testInner.appendChild(clone);
    const childHeight = clone.offsetHeight + 20;

    if (currentHeight + childHeight > maxHeight && currentHeight > 0) {
      currentPage = createNewPage();
      currentInner = currentPage.querySelector('.page-inner');
      canvas.appendChild(currentPage);
      currentHeight = 0;
    }

    currentInner.appendChild(child.cloneNode(true));
    currentHeight += childHeight;
    testInner.innerHTML = '';
  });

  document.body.removeChild(testPage);
  updateMeta();
  updateOutline();
}

function initReportWithData(data) {
  document.getElementById('editor-empty').style.display = 'none';
  const canvas = document.getElementById('report-canvas-content');
  canvas.style.display = 'block';
  pageCount = 0; sectionCount = 0;

  const dateStr = new Date().toLocaleDateString('es-ES', {day:'numeric', month:'long', year:'numeric'});
  const repId = 'KR-' + Math.random().toString(36).substr(2,6).toUpperCase();

  const tempContainer = document.createElement('div');

  // 1. PORTADA (CON SALTO FORZADO)
  const cover = document.createElement('div');
  cover.className = 'force-page-break';
  cover.innerHTML = `
    <div class="report-header" style="margin-top:40mm; border-bottom:3px solid #1a3a6b; padding-bottom:30px; margin-bottom: 20px;">
      <div>
        <div class="report-title" style="font-size:48px; color:#1a3a6b;">INFORME</div>
        <div class="report-title" style="font-size:32px; font-weight:300; letter-spacing:8px;">BIOINFORMÁTICO</div>
      </div>
      <div class="report-meta">
        <div class="meta-id" style="font-size:18px;">${repId}</div>
        <div class="meta-date-val">${dateStr}</div>
        <div class="meta-inst-val">Análisis Genómico Avanzado</div>
        <div class="meta-ver-val">v${escapeHtml(data.engine_version || '—')} (Stable)</div>
      </div>
    </div>
    <div style="margin-top:60mm; text-align:center;">
      <div style="font-family:'IBM Plex Mono',monospace; font-size:12px; color:#666; letter-spacing:2px; text-transform:uppercase;">Preparado para la interpretación clínica de microARNs</div>
      <div style="margin-top:20px; font-size:14px; font-style:italic; color:#1a3a6b;">Convergencia Molecular y Silenciamiento Génico</div>
    </div>`;
  tempContainer.appendChild(cover);

  // 2. SÍNTESIS ACADÉMICA (CON SALTO FORZADO AL INICIO)
  const synthHead = document.createElement('div');
  synthHead.className = 'report-section force-page-break';
  synthHead.innerHTML = `<div class="section-heading"><span class="s-num">I.</span> Síntesis de investigación académica</div>`;
  tempContainer.appendChild(synthHead);

  const synthesisParagraphs = (data.scientific_synthesis || '').split('\n\n');
  synthesisParagraphs.forEach(p => {
    const pTag = document.createElement('p');
    pTag.className = 'editable-block';
    pTag.contentEditable = 'true';
    pTag.style.cssText = "width:100%; font-family:'Spectral', serif; line-height:1.75; font-size:14px; text-align:justify; margin-bottom:18px;";
    pTag.innerHTML = escapeHtml(p).replace(/\n/g, '<br>');
    tempContainer.appendChild(pTag);
  });

  // 2.1 CONTEXTO FUNCIONAL (DENTRO DE SECCIÓN I - SIN TÍTULO REDUNDANTE)
  if (data.functional_context) {
    const funcParagraphs = data.functional_context.split('\n\n');
    funcParagraphs.forEach(p => {
      if (p.trim() === "" || p.toLowerCase().includes("contexto funcional")) return;
      const pTag = document.createElement('p');
      pTag.className = 'editable-block';
      pTag.contentEditable = 'true';
      pTag.style.cssText = "width:100%; font-family:'Spectral', serif; line-height:1.75; font-size:13px; text-align:justify; margin-bottom:15px;";
      pTag.innerHTML = escapeHtml(p).replace(/\n/g, '<br>');
      tempContainer.appendChild(pTag);
    });
  }

  // 2.2 BIBLIOGRAFÍA (SECCIÓN FINAL DE INVESTIGACIÓN)
  if (data.references_text) {
    const refHead = document.createElement('div');
    refHead.className = 'report-section force-page-break';
    refHead.innerHTML = `<div class="section-heading"><span class="s-num">I.1</span> Bibliografía</div>`;
    tempContainer.appendChild(refHead);

    const refParagraphs = data.references_text.split('\n\n');
    refParagraphs.forEach(p => {
      if (p.trim() === "" || p.toLowerCase().startsWith("bibliografía") || p.toLowerCase().startsWith("referencias")) return;
      const pTag = document.createElement('p');
      pTag.className = 'editable-block';
      pTag.contentEditable = 'true';
      pTag.style.cssText = "width:100%; font-family:'Spectral', serif; line-height:1.45; font-size:11.5px; text-align:justify; margin-bottom:12px;";
      
      // Asegurar que las URLs sean hipervínculos reales (H-08: primero se
      // escapa el texto y luego se enlazan las URLs ya escapadas).
      const linkedText = escapeHtml(p).replace(/(https?:\/\/[^\s<]+)/g, '<a href="$1" target="_blank" style="color:#1a3a6b; text-decoration:underline;">$1</a>');
      pTag.innerHTML = linkedText.replace(/\n/g, '<br>');
      tempContainer.appendChild(pTag);
    });
  }

  // 3. TABLA (SALTO FORZADO)
  const tableSec = document.createElement('div');
  tableSec.className = 'report-section force-page-break';
  let tableRows = '';
  data.common_genes.slice(0, 40).forEach(g => {
    const d = (data.gene_details && data.gene_details[g]) || {system:'—', associated_routes:[], pathology:'—'};
    const pathText = (d.pathology||'').length > 180 ? d.pathology.substring(0, 180) + '...' : (d.pathology||'—');
    tableRows += `<tr><td style="font-weight:600; color:#1a3a6b; padding:8px;">${escapeHtml(g)}</td><td style="padding:8px;">${escapeHtml(d.system||'—')}</td><td style="font-style:italic; padding:8px;">${escapeHtml((d.associated_routes||[]).slice(0,2).join('; ')||'—')}</td><td style="text-align:justify; padding:8px; font-size:11px;">${escapeHtml(pathText)}</td></tr>`;
  });
  tableSec.innerHTML = `
    <div class="section-heading"><span class="s-num">II.</span> Panel de Biomarcadores Core Identificados</div>
    <table class="rep-table" style="width:100%; border-collapse:collapse; margin-top:10px;"><thead style="background:#1a3a6b; color:white;"><tr><th style="padding:10px;">Gen Core</th><th style="padding:10px;">Sistema</th><th style="padding:10px;">Rutas Asociadas</th><th style="padding:10px;">Relevancia Patológica</th></tr></thead><tbody id="bio-tbody">${tableRows}</tbody></table>`;
  tempContainer.appendChild(tableSec);

  // 4. DETALLES POR GEN (SALTO FORZADO)
  const detailHead = document.createElement('div');
  detailHead.className = 'report-section force-page-break';
  detailHead.innerHTML = `<div class="section-heading"><span class="s-num">III.</span> Traducción Patológica Detallada</div>`;
  tempContainer.appendChild(detailHead);

  data.common_genes.slice(0, 12).forEach(g => {
    const d = (data.gene_details && data.gene_details[g]) || {system:'—', associated_routes:[], pathology:'—'};
    const block = document.createElement('div');
    block.className = 'report-section gene-block';
    block.style.cssText = "margin-bottom:15px; padding:10px; border-left:3px solid #1a3a6b; background:#f9fbfc;";
    block.innerHTML = `<div class="gene-block-name" style="font-weight:bold; color:#1a3a6b; border-bottom:1px solid #eee; padding-bottom:4px; margin-bottom:6px;">${escapeHtml(g)} — ${escapeHtml(d.system||'—')}</div><div class="gene-block-text" contenteditable="true" style="font-size:12px; line-height:1.4; text-align:justify;">${escapeHtml(d.pathology||'—').replace(/\n/g, '<br>')}</div><div style="font-family:'IBM Plex Mono',monospace; font-size:9px; color:#1a3a6b; margin-top:6px;">Rutas: <i>${(d.associated_routes||[]).map(escapeHtml).join(', ')||'—'}</i></div>`;
    tempContainer.appendChild(block);
  });

  // 5. GRÁFICOS (SALTO FORZADO CON TAMAÑOS FIJOS PARA PDF)
  const vizSec = document.createElement('div');
  vizSec.className = 'report-section force-page-break';
  // v1.52 (bloque 2): mismo guard de base64 que setImg — solo base64
  // válido llega al src de una <img>; cualquier otra cosa es "no disponible".
  const b64Img = (b64, style, fallback) => (b64 && /^[A-Za-z0-9+/=]+$/.test(b64))
    ? `<img src="data:image/png;base64,${b64}" style="${style}">`
    : `<div style="padding:40px; color:#999; border:1px dashed #ccc;">${fallback}</div>`;
  const vennImg = b64Img(data.venn_plot, 'width:380px; height:300px; object-fit:contain; border:1px solid #eee; background:white; padding:10px; margin-bottom:10px;', 'Gráfico no disponible');
  const volcanoImg = b64Img(data.volcano_plot, 'width:480px; height:320px; object-fit:contain; border:1px solid #eee; background:white; padding:10px;', 'Gráfico no disponible');
  
  vizSec.innerHTML = `
    <div class="section-heading"><span class="s-num">IV.</span> Evidencia Gráfica de Convergencia</div>
    <div style="display:flex; flex-direction:column; align-items:center; gap:20px; margin-top:10px;">
      <div class="viz-print" style="text-align:center;">
        <div style="font-size:10px; font-weight:bold; color:#333; margin-bottom:5px;">IV.1 Diagrama de Co-regulación</div>
        ${vennImg}
      </div>
      <div class="viz-print" style="text-align:center;">
        <div style="font-size:10px; font-weight:bold; color:#333; margin-bottom:5px;">IV.2 Paisaje de Significancia Biológica</div>
        ${volcanoImg}
      </div>
    </div>`;
  tempContainer.appendChild(vizSec);

  // 6. PPI (CON SALTO FORZADO)
  const ppiSec = document.createElement('div');
  ppiSec.className = 'report-section force-page-break';
  const ppiImg = b64Img(data.ppi_plot, 'width:85%; max-height:400px; object-fit:contain; border:1px solid #eee; box-shadow:0 2px 4px rgba(0,0,0,0.05);', 'Interactoma no disponible');
  ppiSec.innerHTML = `
    <div class="section-heading"><span class="s-num">V.</span> Interactoma Proteico (STRING-DB)</div>
    <div class="viz-print" style="margin-top:10px; text-align:center;">${ppiImg}</div>`;
  tempContainer.appendChild(ppiSec);

  // Ejecutar paginación real
  setTimeout(() => paginateReport(tempContainer), 100);
}

function init() {
  console.log("KENRYU: Iniciando componentes...");
  populateYearSelect();
  const canvas = document.getElementById('report-canvas-content');
  if (canvas) {
    const observer = new MutationObserver(() => updateOutline());
    observer.observe(canvas, { childList: true, subtree: true, characterData: true });
  }
}

if (document.readyState === 'loading') {
  document.addEventListener('DOMContentLoaded', init);
} else {
  init();
}




