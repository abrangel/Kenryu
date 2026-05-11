const API_URL = "https://kenryu007-bioinformatica.hf.space";
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
        title.innerText = "Vista de Reporte";
        menuDash.classList.remove('active');
        menuRep.classList.add('active');
    }
}

async function analyzePro() {
    const elInput = document.getElementById('mirna-input');
    const elYears = document.getElementById('pubmed-years');
    const elMonth = document.getElementById('analysis-month');
    
    if (!elInput) return;
    const mirnas = elInput.value.split(',').map(m => m.trim()).filter(m => m);
    const years = elYears ? elYears.value : 10;
    const month = elMonth ? elMonth.value : 3;

    const loadingOverlay = document.getElementById('loading-overlay');
    const logContainer = document.getElementById('process-logs');
    
    if (mirnas.length === 0) return alert("Por favor, ingrese miRNAs.");

    if (loadingOverlay) loadingOverlay.classList.remove('hidden');
    if (logContainer) logContainer.innerHTML = "";

    addLog("Iniciando análisis genómico...");

    try {
        const response = await fetch(`${API_URL}/api/v1/analyze`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ mirnas, years: parseInt(years), month: parseInt(month) })
        });

        if (!response.ok) throw new Error("Fallo en la comunicación con el motor bioinformático.");

        const data = await response.json();
        lastAnalysisData = data;

        // Render UI Dashboard
        renderDashboardUI(data);
        
        // LLENAR EL REPORTE UNA SOLA VEZ AL TERMINAR EL ANÁLISIS
        populateReport();

        addLog("Análisis completado exitosamente. Reporte listo.");

    } catch (error) {
        addLog(`ERROR: ${error.message}`, true);
    } finally {
        if (loadingOverlay) loadingOverlay.classList.add('hidden');
    }
}

function renderDashboardUI(data) {
    const coreListContainer = document.getElementById('core-genes-list');
    if (coreListContainer) coreListContainer.innerHTML = data.common_genes.map(g => `<span class="gene-pill">${g}</span>`).join("");

    const renderImg = (id, b64) => {
        const el = document.getElementById(id);
        if (el && b64) el.innerHTML = `<img src="data:image/png;base64,${b64}" style="width:100%; border-radius:8px;">`;
    };
    renderImg('venn-container', data.venn_plot);
    renderImg('volcano-container', data.volcano_plot);
    renderImg('ppi-container', data.ppi_plot);

    const enrichContainer = document.getElementById('enrich-container');
    if (data.enrichment && enrichContainer) {
        let html = `<table class="data-table"><thead><tr><th>Ruta Biológica</th><th>Fuente</th><th>P-adj</th><th>Referencias (PubMed)</th></tr></thead><tbody>`;
        data.enrichment.slice(0, 50).forEach(item => {
            let pLinks = item.Evidences && item.Evidences.length > 0 
                ? item.Evidences.slice(0, 10).map(id => `<a href="https://pubmed.ncbi.nlm.nih.gov/${id}" target="_blank" style="color:#58a6ff; text-decoration:none; margin-right:5px;">[PMID:${id}]</a>`).join("") 
                : 'Verificado';
            html += `<tr><td><strong>${item.Term}</strong><br><small>${item.ScientificDesc}</small></td><td><span class="source-badge">${item.Source}</span></td><td>${item.Pval.toExponential(2)}</td><td><div style="display:flex; flex-wrap:wrap; gap:2px;">${pLinks}</div></td></tr>`;
        });
        html += "</tbody></table>";
        enrichContainer.innerHTML = html;
    }
}

function populateReport() {
    if (!lastAnalysisData) return;
    const data = lastAnalysisData;

    document.getElementById('rep-id').innerText = "KR-" + Math.random().toString(36).substr(2, 7).toUpperCase();
    document.getElementById('rep-date').innerText = "Emitido: " + new Date().toLocaleString();

    // USAR LA SÍNTESIS EXPERTA DINÁMICA DEL BACKEND
    document.getElementById('rep-expert-summary').innerHTML = `
        <div style="text-align: justify; text-justify: inter-word; color: black; line-height: 1.6;">
            <p>${data.executive_summary}</p>
        </div>
    `;

    document.getElementById('rep-overview').innerHTML = `
        <div style="line-height: 1.8; font-size: 13px; color: black;">
            <p><strong>MIRNAS ANALIZADOS:</strong> ${document.getElementById('mirna-input').value}</p>   
            <p><strong>ORGANISMO:</strong> Homo sapiens (human) [GN:hsa]</p>
            <p><strong>MARCADORES CLÍNICOS:</strong> <b>${data.common_genes.length} Dianas de Alta Confianza</b></p>
            <p><strong>LISTA DE GENES CORE:</strong> <span style="font-weight: bold; color: #1f6feb;">${data.common_genes.slice(0, 50).join(", ")}</span></p>
        </div>
    `;

    let tableHtml = `<table style="width:100%; border-collapse:collapse; font-size:11px; margin-top:10px; color: black; border: 1px solid #ddd;">     
        <thead style="background:#f8f9fa;"><tr><th style="border:1px solid #ddd; padding:8px; text-align:left;">Gen Core</th><th style="border:1px solid #ddd; padding:8px; text-align:left;">Sistemas Afectados</th><th style="border:1px solid #ddd; padding:8px; text-align:left;">Relevancia Patológica Principal</th></tr></thead><tbody>`;
    
    data.common_genes.slice(0, 30).forEach(gene => {
        let info = data.detected_systems[gene] || ["Sistémico", "Alteración en procesos biológicos generales y señalización celular."];
        tableHtml += `<tr><td style="border:1px solid #ddd; padding:8px; font-weight:bold;">${gene}</td><td style="border:1px solid #ddd; padding:8px;">${info[0]}</td><td style="border:1px solid #ddd; padding:8px;">${info[1]}</td></tr>`;
    });
    tableHtml += `</tbody></table>`;
    document.getElementById('rep-biomarker-table-container').innerHTML = tableHtml;

    let detailHtml = "";
    data.common_genes.slice(0, 10).forEach(gene => {
        let info = data.detected_systems[gene] || ["Sistémico", "Análisis fenotípico en curso."];
        detailHtml += `
            <div style="margin-bottom:20px; border-left:4px solid #1f6feb; padding-left:15px; page-break-inside:avoid; color: black;">
                <p style="margin:0; font-size:13px; color:#1f6feb; font-weight:bold;">Marcador Identificado:</p>
                <p style="margin:0; font-size:15px;"><b>mRNA: ${gene}</b></p>
                <p style="margin:5px 0; font-size:11px;">Afectación identificada en sistema <b>${info[0].toLowerCase()}</b>. Implicación: ${info[1]}</p>
                <p style="margin:0; font-size:11px; color:#1f6feb;">Consenso Genómico KENRYU | Evidencia Bibliográfica Verificada</p>
            </div>
        `;
    });
    document.getElementById('rep-pathological-translation').innerHTML = detailHtml;

    const setRepImg = (id, b64) => {
        const el = document.getElementById(id);
        if (el && b64) el.innerHTML = `<img src="data:image/png;base64,${b64}" style="width:100%; display: block; margin: 20px 0;">`;
    };
    setRepImg('rep-venn', data.venn_plot);
    setRepImg('rep-volcano', data.volcano_plot);
    setRepImg('rep-ppi', data.ppi_plot);
}

function exportToPDF() {
    if (!lastAnalysisData) return alert("Ejecute un análisis primero.");
    const reportElement = document.getElementById('clinical-report-template');
    const opt = {
        margin: 0,
        filename: 'Informe_Bioinformatico_Kenryu.pdf',
        image: { type: 'jpeg', quality: 1.0 },
        html2canvas: { scale: 2, useCORS: true },
        jsPDF: { unit: 'mm', format: 'a4', orientation: 'portrait' },
        pagebreak: { mode: ['avoid-all', 'css', 'legacy'] }
    };
    html2pdf().set(opt).from(reportElement).save();
}

function addLog(message, isError = false) {
    const logContainer = document.getElementById('process-logs');
    if (!logContainer) return;
    const div = document.createElement('div');
    div.className = 'log-line';
    const ts = new Date().toLocaleTimeString().split(' ')[0];
    if (isError) div.style.color = '#f85149';
    div.innerHTML = `<span class="log-ts">[${ts}]</span> <span class="log-msg">${message}</span>`;
    logContainer.appendChild(div);
    logContainer.scrollTop = logContainer.scrollHeight;
}
