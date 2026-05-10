const API_URL = "http://localhost:8080";
let lastAnalysisData = null;

function switchTab(tab) {
    const dashView = document.getElementById('dashboard-view');
    const repView = document.getElementById('report-view');
    const title = document.getElementById('view-title');
    const menuDash = document.getElementById('menu-dashboard');
    const menuRep = document.getElementById('menu-report');

    if (tab === 'dashboard') {
        dashView.style.display = 'grid';
        repView.style.display = 'none';
        title.innerText = "Panel Principal";
        menuDash.classList.add('active');
        menuRep.classList.remove('active');
    } else {
        dashView.style.display = 'none';
        repView.style.display = 'block';
        title.innerText = "Editor de Informe";
        menuDash.classList.remove('active');
        menuRep.classList.add('active');
        if (lastAnalysisData) populateReport();
    }
}

async function analyzePro() {
    const elInput = document.getElementById('mirna-input');
    const elConsensus = document.getElementById('consensus-mode');
    if (!elInput) return;
    const mirnas = elInput.value.split(',').map(m => m.trim()).filter(m => m);
    const consensusMode = elConsensus ? elConsensus.value : "strict";

    document.getElementById('loading-overlay').classList.remove('hidden');
    document.getElementById('process-logs').innerHTML = "";
    addLog(`Iniciando orquestación bioinformática integral (Modo: ${consensusMode})...`);

    try {
        const response = await fetch(`${API_URL}/api/v1/analyze`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ mirnas, years: 10, consensus_mode: consensusMode })
        });

        if (!response.ok) throw new Error("Fallo en la comunicación con el motor.");
        const data = await response.json();
        lastAnalysisData = data;

        renderDashboardUI(data);
        addLog("✅ Análisis integral finalizado. Evidencia científica consolidada.");
    } catch (error) {
        addLog(`❌ ERROR: ${error.message}`, true);
    } finally {
        document.getElementById('loading-overlay').classList.add('hidden');
    }
}

function renderDashboardUI(data) {
    const list = document.getElementById('core-genes-list');
    if (list) {
        list.innerHTML = `<h3 style="width:100%; font-size:0.85rem; margin-bottom:10px; color:#1a3a6b;">Genes Core Identificados</h3>` + 
            data.common_genes.map(g => `<span class="gene-pill" onclick="showGeneSidebar('${g}')">${g}</span>`).join("");
    }

    const render = (id, b64) => {
        const el = document.getElementById(id);
        if (el && b64) el.innerHTML = `<img src="data:image/png;base64,${b64}" style="width:100%; border-radius:8px; border: 1px solid #ddd;">`;
    };
    render('venn-container', data.venn_plot);
    render('volcano-container', data.volcano_plot);
    render('ppi-container', data.ppi_plot);

    const enrich = document.getElementById('enrich-container');
    if (data.enrichment && enrich) {
        let html = `<table class="data-table"><thead><tr><th>Ruta Biológica</th><th>Fuente</th><th>p-valor</th><th>PubMed</th></tr></thead><tbody>`;
        data.enrichment.forEach(item => {
            let pLink = item.Evidence ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${item.Evidence.id}" target="_blank">PMID: ${item.Evidence.id}</a>` : '—';
            html += `<tr><td><strong>${item.Term}</strong><br><small>${item.ScientificDesc}</small></td><td>${item.Source}</td><td>${item.Pval.toExponential(2)}</td><td>${pLink}</td></tr>`;
        });
        enrich.innerHTML = html + "</tbody></table>";
    }
}

async function showGeneSidebar(gene) {
    let sidebar = document.getElementById('gene-detail-panel');
    if (!sidebar) {
        sidebar = document.createElement('div');
        sidebar.id = 'gene-detail-panel';
        sidebar.style.position = 'fixed'; sidebar.style.right = '0'; sidebar.style.top = '0';
        sidebar.style.height = '100vh'; sidebar.style.width = '340px';
        sidebar.style.background = 'white'; sidebar.style.zIndex = '2000';
        sidebar.style.padding = '25px'; sidebar.style.borderLeft = '1px solid #ddd';
        sidebar.style.boxShadow = '-5px 0 20px rgba(0,0,0,0.1)'; sidebar.style.overflowY = 'auto';
        document.body.appendChild(sidebar);
    }
    sidebar.innerHTML = `<div style="padding:40px; text-align:center;">Analizando evidencia clínica...</div>`;
    sidebar.style.display = 'block';

    try {
        const resp = await fetch(`${API_URL}/api/v1/gene-info/${gene}`);
        const info = await resp.json();
        const detail = lastAnalysisData.gene_details[gene];
        sidebar.innerHTML = `
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:20px;">
                <h2 style="margin:0; color:#1a3a6b; font-size:22px;">${gene}</h2>
                <button onclick="document.getElementById('gene-detail-panel').style.display='none'" style="cursor:pointer; border:none; background:none; font-size:28px; color:#999;">&times;</button>
            </div>
            <p style="font-weight:700; color:#333; margin-bottom:15px;">${info.fullName || 'Biomarcador Core'}</p>
            <div style="background:#e8eef7; padding:12px; border-radius:6px; margin-bottom:20px; border:1px solid #d0ddec;">
                <span style="font-size:11px; font-weight:800; color:#1a3a6b; text-transform:uppercase;">CONFIANZA: ${detail.confidence.score}%</span>
            </div>
            <h4 style="font-size:14px; color:#1a3a6b; border-bottom:1px solid #eee; padding-bottom:5px;">Impacto Patológico</h4>
            <p style="font-size:13px; line-height:1.6; text-align:justify; color:#444;">${detail.pathology}</p>
            <h4 style="font-size:14px; color:#1a3a6b; margin-top:20px;">Referencia PubMed</h4>
            ${detail.pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${detail.pmid}" target="_blank" style="font-size:12px; color:#1a3a6b; font-weight:700; text-decoration:underline;">Ver Evidencia Original (PMID: ${detail.pmid})</a>` : '<p style="font-size:12px; color:#888;">Validación bioinformática interna.</p>'}
        `;
    } catch (e) { sidebar.innerHTML = `<div style="padding:20px; color:red;">Error al cargar información clínica.</div>`; }
}

function populateReport() {
    if (!lastAnalysisData) return;
    const data = lastAnalysisData;

    document.getElementById('rep-id-tag').innerText = "INFORME-" + Math.random().toString(36).substr(2, 6).toUpperCase();
    document.getElementById('rep-date-tag').innerText = new Date().toLocaleDateString('es-ES', { day: 'numeric', month: 'long', year: 'numeric' });

    const pathways = data.enrichment.slice(0, 8).map(e => `<li><b>${e.Term}:</b> ${e.ScientificDesc}</li>`).join("");

    // SÍNTESIS ACADÉMICA TÉCNICA
    document.getElementById('rep-synthesis').innerHTML = `
        <p>El análisis genómico integral revela una convergencia molecular significativa mediada por el panel de microRNAs analizado. Mediante el empleo de algoritmos de intersección multivariante, se ha identificado un núcleo regulatorio compuesto por <b>${data.common_genes.length} biomarcadores core</b> sometidos a represión post-transcripcional simultánea.</p>
        <p style="margin-top:10px;">La desregulación coordinada de estos genes impacta de manera directa en la estabilidad de procesos celulares y rutas metabólicas críticas para la homeostasis celular, destacando los siguientes hallazgos prioritarios identificados en la red molecular:</p>
        <ul style="margin:10px 0; padding-left:20px;">${pathways}</ul>
        <p style="margin-top:10px;">Este paisaje epigenético, validado mediante la integración de predictores bioinformáticos y evidencia experimental recuperada de la base de datos PubMed (NCBI), sugiere que el panel analizado actúa como un orquestador central de la respuesta sistémica en los fenotipos identificados, proporcionando un fundamento sólido para futuras validaciones funcionales.</p>
    `;

    // TABLA PROFESIONAL 4 COLUMNAS
    let tableHtml = `<table style="width:100%; border-collapse:collapse; border:2.5px solid #000; font-size:9.5px; color:black;">
        <thead style="background:#eee;"><tr>
            <th style="border:1px solid #000; padding:8px; text-align:left; width:15%;">Gen Core</th>
            <th style="border:1px solid #000; padding:8px; text-align:left; width:20%;">Sistemas Afectados</th>
            <th style="border:1px solid #000; padding:8px; text-align:left; width:25%;">Rutas Biológicas Asociadas</th>
            <th style="border:1px solid #000; padding:8px; text-align:left; width:40%;">Relevancia Patológica Consolidada</th>
        </tr></thead><tbody>`;
    
    data.common_genes.slice(0, 35).forEach(gene => {
        const d = data.gene_details[gene] || { system: "Sistémico", pathology: "Pendiente de validación." };
        const associated = data.enrichment.filter(e => e.Genes.split(';').map(x => x.trim()).includes(gene)).slice(0, 3).map(e => e.Term).join("; ");
        tableHtml += `<tr>
            <td style="border:1px solid #000; padding:6px; font-weight:bold; font-size:10.5px;">${gene}</td>
            <td style="border:1px solid #000; padding:6px;">${d.system}</td>
            <td style="border:1px solid #000; padding:6px; font-style:italic;">${associated || 'Interacción directa miRNA-mRNA'}</td>
            <td style="border:1px solid #000; padding:6px; text-align:justify;">${d.pathology}</td>
        </tr>`;
    });
    tableHtml += `</tbody></table>`;
    document.getElementById('rep-table-container').innerHTML = tableHtml;

    // TRADUCCIÓN PATOLÓGICA DETALLADA
    let detailHtml = "";
    data.common_genes.slice(0, 12).forEach(gene => {
        const d = data.gene_details[gene];
        detailHtml += `
            <div class="evidence-block" style="border-left: 5px solid #1a3a6b; background: #f8f9fa; padding: 15px; margin-bottom: 20px; border-radius: 0 4px 4px 0;">
                <div style="font-weight:800; color:#1a3a6b; font-size:15px; margin-bottom:8px; display:flex; justify-content:space-between; border-bottom: 1px solid #d0ddec; padding-bottom:5px;">
                    <span>${gene}</span>
                    <span style="font-size:10px; font-weight:400; color:#666;">REFERENCIA CIENTÍFICA</span>
                </div>
                <div style="font-size:13px; color:#333; line-height:1.6; text-align:justify;">
                    Este biomarcador constituye un nodo regulatorio crítico del sistema <b>${d.system.toLowerCase()}</b>. Su asociación patológica consolidada incluye: ${d.pathology}. Presenta un índice de confianza del <b>${d.confidence.score}%</b> basado en la convergencia de bases de datos multi-ómicas.
                </div>
                <div style="margin-top:12px; font-size:11px; color:#1a3a6b; font-weight:700;">
                    ${d.pmid ? `<i class="fas fa-external-link-alt"></i> Fuente PubMed: <a href="https://pubmed.ncbi.nlm.nih.gov/${d.pmid}" target="_blank" style="color:#1a3a6b; text-decoration:underline;">PMID ${d.pmid} (Verificado)</a>` : 'Sustento bioinformático integral.'}
                </div>
            </div>
        `;
    });
    document.getElementById('rep-detail-container').innerHTML = detailHtml;

    // RENDER IMÁGENES AL REPORTE
    const setImg = (id, b64) => {
        const el = document.getElementById(id);
        if (el && b64) el.innerHTML = `<img src="data:image/png;base64,${b64}" style="width:100%; border:1.5px solid #000; margin:15px 0; background:white;">`;
    };
    setImg('rep-venn-img', data.venn_plot);
    setImg('rep-volcano-img', data.volcano_plot);
    setImg('rep-ppi-img', data.ppi_plot);
}

function execCmd(cmd) { document.execCommand(cmd, false, null); }
function addCustomPage() {
    const canvas = document.getElementById('report-canvas');
    const newPage = document.createElement('div');
    newPage.className = "a4-page";
    newPage.style.pageBreakBefore = "always";
    newPage.innerHTML = `<div class="academic-text" contenteditable="true" data-placeholder="Escriba aquí contenido adicional de investigación..."></div>`;
    canvas.appendChild(newPage);
}

function exportToPDF() {
    if (!lastAnalysisData) return alert("Realice un análisis primero.");
    const element = document.getElementById('report-canvas');
    const opt = { 
        margin: 0, 
        filename: 'Informe_Investigacion_KENRYU.pdf', 
        html2canvas: { scale: 2, useCORS: true }, 
        jsPDF: { unit: 'mm', format: 'a4', orientation: 'portrait' }
    };
    html2pdf().set(opt).from(element).save();
}

function exportToMarkdown() {
    if (!lastAnalysisData) return alert("Sin datos disponibles.");
    const data = lastAnalysisData;
    let md = `# INFORME DE INVESTIGACIÓN CIENTÍFICA - KENRYU\n\n**Fecha:** ${new Date().toLocaleDateString()}\n\n`;
    md += `## I. SÍNTESIS DE INVESTIGACIÓN\n${document.getElementById('rep-synthesis').innerText}\n\n`;
    md += `## II. PANEL DE BIOMARCADORES\n| Gen Core | Sistema Afectado | Rutas Asociadas | Patología Consolidada |\n| :--- | :--- | :--- | :--- |\n`;
    data.common_genes.slice(0, 30).forEach(g => {
        const d = data.gene_details[g];
        md += `| **${g}** | ${d.system} | - | ${d.pathology} |\n`;
    });
    md += `\n## III. EVIDENCIA PUBMED\n`;
    data.common_genes.slice(0, 15).forEach(g => {
        const d = data.gene_details[g];
        if (d.pmid) md += `- **${g}**: Evidencia en [PMID ${d.pmid}](https://pubmed.ncbi.nlm.nih.gov/${d.pmid})\n`;
    });
    const blob = new Blob([md], {type: 'text/markdown;charset=utf-8'});
    const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download = "Reporte_KENRYU.md"; a.click();
}

function addLog(msg, error = false) {
    const log = document.getElementById('process-logs');
    const div = document.createElement('div');
    if (error) div.style.color = '#ff7b72';
    div.innerHTML = `<span style="color:#888; font-size: 11px;">[${new Date().toLocaleTimeString()}]</span> <span style="font-size: 12px; color:#e1e4e8;">${msg}</span>`;
    log.appendChild(div);
    log.scrollTop = log.scrollHeight;
}
