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
        title.innerText = "Panel de Control";
        menuDash.classList.add('active');
        menuRep.classList.remove('active');
    } else {
        dashView.style.display = 'none';
        repView.style.display = 'block';
        title.innerText = "Vista del Informe";
        menuDash.classList.remove('active');
        menuRep.classList.add('active');
    }
}

async function updateLocalDB() {
    addLog("Sincronizando fuentes genómicas...");
    try {
        const response = await fetch(`${API_URL}/api/v1/update_db`, { method: 'POST' });
        const data = await response.json();
        addLog("✅ " + data.message);
        alert("Sincronización completada con éxito.");
    } catch (error) {
        addLog("❌ Error: " + error.message, true);
    }
}

async function analyzePro() {
    const elInput = document.getElementById('mirna-input');
    const mirnas = elInput.value.split(',').map(m => m.trim()).filter(m => m);
    if (mirnas.length === 0) return alert("Ingrese miRNAs para analizar.");

    const consensusMode = document.getElementById('consensus-mode')?.value || 'strict';

    document.getElementById('loading-overlay').classList.remove('hidden');
    document.getElementById('process-logs').innerHTML = "";
    addLog("Iniciando análisis integral...");

    try {
        const response = await fetch(`${API_URL}/api/v1/analyze`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ mirnas, years: 10, consensus_mode: consensusMode })
        });

        if (!response.ok) throw new Error("Fallo en la comunicación con el servidor.");
        const data = await response.json();
        lastAnalysisData = data;
        localStorage.setItem('last_kenryu_analysis', JSON.stringify(data));

        renderDashboardUI(data);
        await populateReport();
        addLog("✅ Análisis finalizado. Datos reales cargados.");
    } catch (error) {
        addLog(`❌ Error Crítico: ${error.message}`, true);
    } finally {
        document.getElementById('loading-overlay').classList.add('hidden');
    }
}

function renderDashboardUI(data) {
    const list = document.getElementById('core-genes-list');
    
    // Crear mapa de genes enriquecidos
    const enrichedGenes = new Set();
    data.enrichment?.forEach(item => {
        item.Genes.split(';').forEach(g => enrichedGenes.add(g.trim()));
    });

    list.innerHTML = `<h3 style="width:100%; font-size:0.9rem; margin-bottom:10px;">Biomarcadores Core (${data.common_genes.length})</h3>` + 
                     data.common_genes.map(gene => {
                         const isEnriched = enrichedGenes.has(gene);
                         return `<span class="gene-pill ${isEnriched ? 'gene-pill--enriched' : ''}" 
                                       onclick="showGeneSidebar('${gene}')"
                                       title="${isEnriched ? 'Gen con enriquecimiento funcional' : 'Gen diana común'}">
                                    ${gene}${isEnriched ? ' <span class="pill-dot"></span>' : ''}
                                 </span>`;
                     }).join("");

    const render = (id, b64) => {
        const el = document.getElementById(id);
        if (el && b64) el.innerHTML = `<img src="data:image/png;base64,${b64}" style="width:100%; border-radius:4px;">`;
    };
    render('venn-container', data.venn_plot);
    render('volcano-container', data.volcano_plot);

    const enrich = document.getElementById('enrich-container');
    let html = `<table class="data-table"><thead><tr><th>Término Biológico</th><th>p-valor</th><th>Fuente</th></tr></thead><tbody>`;
    data.enrichment.slice(0, 30).forEach(item => {
        const pval = item.PvalAdj || item.Pval;
        html += `<tr><td><b>${item.Term}</b></td><td>${pval.toExponential(2)}</td><td><span class="status-badge" style="background:#333;">${item.Source}</span></td></tr>`;
    });
    enrich.innerHTML = html + "</tbody></table>";
}

async function showGeneSidebar(geneSymbol) {
    const panel = document.getElementById('gene-detail-panel');
    panel.innerHTML = `
        <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:16px;">
            <h3 style="margin:0; color:#58a6ff; font-size:18px;">${geneSymbol}</h3>
            <button onclick="document.getElementById('gene-detail-panel').style.transform='translateX(100%)'"
                    style="background:none; border:none; color:#666; font-size:18px; cursor:pointer;">✕</button>
        </div>
        <div style="color:#666; font-size:12px;">⏳ Consultando bases de datos clínicas...</div>`;
    panel.style.transform = 'translateX(0)';

    try {
        const [infoResp, diseaseResp] = await Promise.all([
            fetch(`${API_URL}/api/v1/gene-info/${geneSymbol}`),
            fetch(`${API_URL}/api/v1/diseases/${geneSymbol}`)
        ]);

        const info = infoResp.ok ? await infoResp.json() : {};
        const diseases = diseaseResp.ok ? await diseaseResp.json() : { diseases: [] };

        const diseasesHtml = diseases.diseases?.length > 0
            ? diseases.diseases.slice(0, 6).map(d => `
                <div style="padding:6px 8px; background:#252525; border-radius:4px; margin-bottom:4px; font-size:11px;">
                    <span style="color:#e8eef7;">${d.name}</span>
                    <span style="float:right; color:#555; font-size:10px;">${d.source}</span>
                </div>`).join('')
            : '<p style="color:#555; font-size:11px;">Sin asociaciones curadas disponibles</p>';

        const confidence = lastAnalysisData?.gene_details[geneSymbol]?.confidence || { score: 0, label: 'N/A' };

        panel.innerHTML = `
            <div style="display:flex; justify-content:space-between; align-items:center; margin-bottom:16px;">
                <h3 style="margin:0; color:#58a6ff; font-size:18px;">${geneSymbol}</h3>
                <button onclick="document.getElementById('gene-detail-panel').style.transform='translateX(100%)'"
                        style="background:none; border:none; color:#666; font-size:18px; cursor:pointer;">✕</button>
            </div>
            ${info.fullName ? `<p style="color:#aaa; font-size:12px; margin:0 0 12px;">${info.fullName}</p>` : ''}
            <div style="margin-bottom:16px; padding: 10px; background: #252525; border-radius: 4px;">
                <div style="font-size:10px; color:#555; text-transform:uppercase; margin-bottom:4px;">Score de Confianza</div>
                <div style="display:flex; align-items:center; gap:10px;">
                    <div style="flex:1; height:6px; background:#444; border-radius:3px; overflow:hidden;">
                        <div style="width:${confidence.score}%; height:100%; background:#56d364;"></div>
                    </div>
                    <span style="font-size:12px; font-weight:bold; color:#56d364;">${confidence.score}</span>
                </div>
                <div style="font-size:10px; color:#888; margin-top:4px;">${confidence.label}</div>
            </div>
            ${info.summary ? `<div style="margin-bottom:16px;"><div style="font-size:10px; color:#555; text-transform:uppercase; margin-bottom:6px;">Función</div><p style="font-size:11px; color:#ccc; line-height:1.5; margin:0;">${info.summary}</p></div>` : ''}
            <div style="margin-bottom:16px;"><div style="font-size:10px; color:#555; text-transform:uppercase; margin-bottom:8px;">Enfermedades asociadas</div>${diseasesHtml}</div>
            <div style="font-size:10px; color:#444; margin-top:16px; padding-top:12px; border-top:1px solid #2a2a2a;">Fuentes: MyGene.info · DisGeNET · OpenTargets</div>`;
    } catch (err) {
        panel.innerHTML += `<p style="color:#ff7b72; font-size:11px;">Error al cargar: ${err.message}</p>`;
    }
}

async function populateReport() {
    if (!lastAnalysisData) return;
    const data = lastAnalysisData;

    document.getElementById('rep-id').innerText = "INFORME-" + Math.random().toString(36).substr(2, 6).toUpperCase();
    document.getElementById('rep-date').innerText = new Date().toLocaleDateString();

    document.getElementById('rep-expert-summary').innerHTML = `
        <p style="text-align:justify; line-height:1.6;">Este documento presenta los resultados del análisis bioinformático sobre el panel de microRNAs. Se ha identificado una convergencia regulatoria sobre <b>${data.common_genes.length} genes comunes</b> mediante la intersección de fuentes locales y validadas utilizando el modo de consenso <b>${document.getElementById('consensus-mode').options[document.getElementById('consensus-mode').selectedIndex].text}</b>.</p>
    `;

    document.getElementById('rep-overview').innerHTML = `
        <div style="font-size:12px; line-height:1.5;">
            <p><b>miRNAs Analizados:</b> ${document.getElementById('mirna-input').value}</p>
            <p><b>Modo Consenso:</b> ${document.getElementById('consensus-mode').value}</p>
            <p><b>Número de Genes Core:</b> ${data.common_genes.length}</p>
        </div>
    `;

    // TABLA DE BIOMARCADORES DINÁMICA
    const genes = data.common_genes.slice(0, 15);
    let table = `<table style="width:100%; border-collapse:collapse; border:1.5px solid #000; font-size:10px;">
        <thead><tr style="background:#eee;">
            <th style="border:1px solid #000; padding:8px; text-align:left;">Gen Core</th>
            <th style="border:1px solid #000; padding:8px; text-align:left;">Nombre completo</th>
            <th style="border:1px solid #000; padding:8px; text-align:left;">Enfermedades asociadas</th>
            <th style="border:1px solid #000; padding:8px; text-align:left;">Confianza</th>
        </tr></thead><tbody id="dynamic-bio-tbody">
            <tr><td colspan="4" style="padding:12px; text-align:center; color:#888;">⏳ Consultando bases de datos clínicas...</td></tr>
        </tbody></table>`;
    document.getElementById('rep-biomarker-table-container').innerHTML = table;

    const tbody = document.getElementById('dynamic-bio-tbody');
    tbody.innerHTML = '';

    const promises = genes.map(gene =>
        fetch(`${API_URL}/api/v1/gene-info/${gene}`).then(r => r.ok ? r.json() : null).catch(() => null)
    );
    const results = await Promise.all(promises);

    results.forEach((info, idx) => {
        const gene = genes[idx];
        const fullName = info?.fullName || '—';
        const diseases = info?.diseases?.length > 0
            ? info.diseases.map(d => `<span style="background:#e8eef7; color:#1a3a6b; padding:1px 4px; border-radius:4px; margin:1px; display:inline-block; font-size:9px;">${d.name}</span>`).join(' ')
            : '<span style="color:#888; font-size:9px;">Sin asociaciones significativas</span>';
        const confidence = data.gene_details[gene]?.confidence?.score || 0;

        const row = document.createElement('tr');
        row.innerHTML = `
            <td style="border:1px solid #000; padding:6px; font-weight:bold;">${gene}</td>
            <td style="border:1px solid #000; padding:6px; font-size:9px;">${fullName}</td>
            <td style="border:1px solid #000; padding:6px;">${diseases}</td>
            <td style="border:1px solid #000; padding:6px; text-align:center;">${confidence}%</td>`;
        tbody.appendChild(row);
    });

    // DETALLES DESCRIPTIVOS USANDO PMIDs REALES
    let details = "";
    Object.keys(data.gene_details).slice(0, 8).forEach(gene => {
        const info = data.gene_details[gene];
        const pmid = info.pmid;
        details += `
            <div style="margin-top:15px; border-left:4px solid #1a3a6b; padding-left:12px; page-break-inside:avoid;">
                <p style="margin:0; font-size:13px; font-weight:bold;">Biomarcador: ${gene}</p>
                <p style="margin:2px 0; font-size:10px; color:#666;">Sistema: ${info.system}</p>
                <p style="margin:4px 0; font-size:11px; text-align:justify;">${info.desc}</p>
                <p style="margin:0; font-size:9px; color:#1a3a6b;">
                    ${pmid ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${pmid}" target="_blank" style="color:#1a3a6b; text-decoration:none;">📄 PMID: ${pmid}</a>` : '<span style="color:#888;">Sin evidencia PubMed directa</span>'}
                </p>
            </div>
        `;
    });
    document.getElementById('rep-pathological-translation').innerHTML = details;

    // BIBLIOGRAFÍA
    let bib = '<ol style="padding-left: 20px; font-size: 8.5px; line-height:1.4;">';
    Object.keys(data.gene_details).forEach(gene => {
        const p = data.gene_details[gene].pmid;
        if (p) bib += `<li style="margin-bottom:4px;">Evidencia molecular para <b>${gene}</b>. [PMID: ${p}](https://pubmed.ncbi.nlm.nih.gov/${p})</li>`;
    });
    data.enrichment.slice(0, 15).forEach(item => {
        if (item.Evidence) bib += `<li style="margin-bottom:4px;">Ruta: <b>${item.Term}</b>. [PMID: ${item.Evidence.id}](https://pubmed.ncbi.nlm.nih.gov/${item.Evidence.id})</li>`;
    });
    bib += '</ol>';
    document.getElementById('rep-bibliography').innerHTML = bib;

    const setImg = (id, b64) => {
        const el = document.getElementById(id);
        if (el && b64) el.innerHTML = `<img src="data:image/png;base64,${b64}" style="max-width:100%; display:block; margin:20px auto; border:1px solid #ddd;">`;
    };
    setImg('rep-venn', data.venn_plot);
    setImg('rep-volcano', data.volcano_plot);
}


async function exportToPDF() {
    if (!lastAnalysisData) return alert("Ejecute un análisis primero.");
    
    addLog("Iniciando generación de PDF profesional...");
    
    // Clonación del reporte para asegurar renderizado perfecto e independiente de la vista
    const original = document.getElementById('clinical-report-template');
    const clone = original.cloneNode(true);
    
    // Configuración de visibilidad forzada para la captura
    clone.style.position = 'fixed';
    clone.style.left = '-10000px';
    clone.style.top = '0';
    clone.style.display = 'block';
    clone.style.width = '210mm';
    document.body.appendChild(clone);

    // Pausa técnica para renderizado de imágenes Base64
    await new Promise(r => setTimeout(r, 2500));

    const opt = {
        margin: 10,
        filename: 'Reporte_Bioinformatico_Final.pdf',
        image: { type: 'jpeg', quality: 1.0 },
        html2canvas: { 
            scale: 2, 
            useCORS: true, 
            backgroundColor: '#ffffff',
            logging: false
        },
        jsPDF: { unit: 'mm', format: 'a4', orientation: 'portrait' },
        pagebreak: { mode: ['avoid-all', 'css', 'legacy'] }
    };

    try {
        await html2pdf().set(opt).from(clone).save();
        addLog("✅ PDF Descargado con éxito.");
    } catch (e) {
        addLog("❌ Error en exportación: " + e.message, true);
    } finally {
        document.body.removeChild(clone);
    }
}

function exportToMarkdown() {
    if (!lastAnalysisData) return alert("Sin datos.");
    const data = lastAnalysisData;
    let md = `# Informe Bioinformático Profesional\n\n**Fecha:** ${new Date().toLocaleDateString()}\n\n`;
    md += `## 1. Resumen\nSe identificaron **${data.common_genes.length} genes comunes**.\n\n`;
    md += `## 2. Biomarcadores Core\n| Gen | Sistema | Referencia |\n| :--- | :--- | :--- |\n`;
    Object.keys(data.gene_details).forEach(gene => {
        const info = data.gene_details[gene];
        md += `| **${gene}** | ${info.system} | PMID ${info.pmid} |\n`;
    });
    md += `\n## 3. Bibliografía Detallada\n`;
    Object.keys(data.gene_details).forEach(gene => {
        const p = data.gene_details[gene].pmid;
        if(p !== "N/A") md += `- Evidencia Gen ${gene}: [PubMed ${p}](https://pubmed.ncbi.nlm.nih.gov/${p})\n`;
    });

    const blob = new Blob([md], {type: 'text/markdown'});
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = "Informe_Final.md";
    a.click();
}

function openEditor() {
    window.open('../../kenryu_report_editor.html', '_blank');
}

function addLog(msg, error = false) {
    const logs = document.getElementById('process-logs');
    const div = document.createElement('div');
    div.style.padding = "2px 0";
    if (error) div.style.color = "#ff7b72";
    div.innerText = `[${new Date().toLocaleTimeString()}] ${msg}`;
    logs.appendChild(div);
    logs.scrollTop = logs.scrollHeight;
}
