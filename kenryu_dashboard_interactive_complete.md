# KENRYU — Dashboard Interactivo Completo

El problema principal es que el ejemplo anterior tenía:

- Datos mock muy pequeños
- Venn.js básico
- Sin animaciones avanzadas
- Sin estilo tipo Bioinformatics Dashboard real
- Sin Chart.js avanzado
- Sin hover profesional
- Sin glassmorphism profundo
- Sin paneles científicos modernos

Aquí tienes la estructura CORRECTA para que realmente se vea como tu página avanzada.

---

# 1. AGREGA ESTAS LIBRERÍAS

```html
<script src="https://cdn.plot.ly/plotly-2.30.0.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
<script src="https://unpkg.com/cytoscape@3.28.1/dist/cytoscape.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/chartjs-plugin-datalabels@2"></script>
<script src="https://d3js.org/d3.v7.min.js"></script>
<script src="https://cdn.jsdelivr.net/npm/venn.js@0.2.20/build/venn.min.js"></script>
```

---

# 2. VENN INTERACTIVO REAL

Reemplaza completamente tu función renderVenn por esta.

```js
function renderVenn(){

  const sets = [
    {sets:['miR-33a'], size:120, label:'miR-33a'},
    {sets:['miR-144'], size:90, label:'miR-144'},
    {sets:['miR-758'], size:70, label:'miR-758'},

    {sets:['miR-33a','miR-144'], size:40},
    {sets:['miR-33a','miR-758'], size:25},
    {sets:['miR-144','miR-758'], size:18},

    {sets:['miR-33a','miR-144','miR-758'], size:12}
  ];

  d3.select('#venn-container')
    .selectAll('*')
    .remove();

  const chart = venn.VennDiagram()
    .width(520)
    .height(380);

  const div = d3.select('#venn-container');

  div.datum(sets).call(chart);

  div.selectAll('path')
    .style('stroke-opacity',1)
    .style('stroke','#ffffff')
    .style('stroke-width',2)
    .style('fill-opacity',0.45)
    .style('transition','all .3s ease');

  div.selectAll('g')
    .on('mouseover', function(event,d){

      venn.sortAreas(div,d);

      d3.select(this)
        .select('path')
        .transition()
        .duration(250)
        .style('fill-opacity',0.75)
        .style('stroke-width',4);

    })

    .on('mouseout', function(){

      d3.select(this)
        .select('path')
        .transition()
        .duration(250)
        .style('fill-opacity',0.45)
        .style('stroke-width',2);

    });
}
```

---

# 3. VOLCANO PLOT PROFESIONAL

Tu volcano debe verse tipo Nature/Cell.

```js
function renderVolcano(){

  const genes = [
    'SCN1A','ABCA1','KPNA3','SNTB2','TSC22D2'
  ];

  const trace = {

    x:[2.3,-1.2,3.1,-2.4,1.8],

    y:[5.2,3.9,7.1,2.7,4.5],

    text:genes,

    mode:'markers+text',

    type:'scatter',

    textposition:'top center',

    marker:{
      size:[18,15,22,12,17],

      color:[
        '#e05c5c',
        '#4fc3a1',
        '#c8a96e',
        '#6eaadc',
        '#9d7dd1'
      ],

      opacity:0.95,

      line:{
        color:'#ffffff',
        width:1.2
      }
    },

    hovertemplate:
      '<b>%{text}</b><br>'+
      'Log2FC: %{x}<br>'+
      '-log10(p): %{y}<extra></extra>'
  };

  Plotly.newPlot(
    'volcano-container',
    [trace],
    {
      paper_bgcolor:'#101318',
      plot_bgcolor:'#101318',

      font:{
        family:'DM Sans',
        color:'#edf3f8'
      },

      margin:{
        l:60,
        r:30,
        t:20,
        b:60
      },

      xaxis:{
        title:'Log2 Fold Change',
        gridcolor:'rgba(255,255,255,.06)',
        zerolinecolor:'rgba(255,255,255,.2)'
      },

      yaxis:{
        title:'−log10(p-value)',
        gridcolor:'rgba(255,255,255,.06)'
      }

    },

    {
      responsive:true,
      displayModeBar:false
    }
  );
}
```

---

# 4. PPI NETWORK ESTILO STRING-DB

```js
function renderPPI(){

  cytoscape({

    container:document.getElementById('ppi-container'),

    elements:[

      { data:{ id:'SCN1A' } },
      { data:{ id:'ABCA1' } },
      { data:{ id:'KPNA3' } },
      { data:{ id:'SNTB2' } },
      { data:{ id:'TSC22D2' } },

      { data:{ source:'SCN1A', target:'KPNA3' } },
      { data:{ source:'SCN1A', target:'ABCA1' } },
      { data:{ source:'KPNA3', target:'TSC22D2' } },
      { data:{ source:'ABCA1', target:'SNTB2' } },
      { data:{ source:'SCN1A', target:'SNTB2' } },

    ],

    style:[

      {
        selector:'node',

        style:{
          'background-color':'#c8a96e',
          'border-width':2,
          'border-color':'#ffffff',
          'label':'data(id)',
          'color':'#ffffff',
          'font-size':'11px',
          'text-valign':'center',
          'text-halign':'center',
          'width':45,
          'height':45,

          'overlay-opacity':0,

          'shadow-blur':20,
          'shadow-color':'#c8a96e',
          'shadow-opacity':0.35
        }
      },

      {
        selector:'edge',

        style:{
          'width':2.5,
          'line-color':'#4fc3a1',
          'curve-style':'bezier',
          'opacity':0.7
        }
      },

      {
        selector:'node:hover',

        style:{
          'background-color':'#ffffff',
          'color':'#000000',
          'width':58,
          'height':58,
          'font-size':'14px'
        }
      }
    ],

    layout:{
      name:'cose',
      animate:true,
      padding:30
    }
  });
}
```

---

# 5. CHART.JS COMO TU WEB

Agrega un nuevo canvas:

```html
<div class="card span2">

  <div class="card-head">
    <div class="card-title">
      Functional Enrichment Distribution
    </div>
  </div>

  <div class="card-body">
    <canvas id="barChart" height="120"></canvas>
  </div>

</div>
```

---

Y luego:

```js
function renderBarChart(){

  new Chart(document.getElementById('barChart'),{

    type:'bar',

    data:{

      labels:[
        'Cholesterol metabolism',
        'Ion channel activity',
        'Cardiac contraction',
        'ABC transporters',
        'Synaptic signaling'
      ],

      datasets:[{

        label:'Combined Score',

        data:[912,744,612,580,490],

        backgroundColor:[
          'rgba(200,169,110,.65)',
          'rgba(79,195,161,.65)',
          'rgba(224,92,92,.65)',
          'rgba(110,170,220,.65)',
          'rgba(157,125,209,.65)'
        ],

        borderColor:[
          '#c8a96e',
          '#4fc3a1',
          '#e05c5c',
          '#6eaadc',
          '#9d7dd1'
        ],

        borderWidth:1.5,
        borderRadius:8
      }]
    },

    options:{

      responsive:true,

      plugins:{
        legend:{
          display:false
        }
      },

      scales:{

        x:{
          ticks:{
            color:'#94a6bb'
          },

          grid:{
            color:'rgba(255,255,255,.05)'
          }
        },

        y:{
          ticks:{
            color:'#94a6bb'
          },

          grid:{
            color:'rgba(255,255,255,.05)'
          }
        }
      }
    }
  });
}
```

---

# 6. GLASSMORPHISM REAL

Tu diseño se verá muchísimo mejor si reemplazas `.card` por:

```css
.card{

  background:linear-gradient(
    145deg,
    rgba(16,20,27,.92),
    rgba(16,20,27,.72)
  );

  border:1px solid rgba(255,255,255,.08);

  backdrop-filter:blur(22px);

  box-shadow:
    0 10px 35px rgba(0,0,0,.45),
    inset 0 1px 0 rgba(255,255,255,.04);

  transition:all .3s ease;
}

.card:hover{

  transform:translateY(-6px);

  border-color:rgba(200,169,110,.28);

  box-shadow:
    0 20px 45px rgba(0,0,0,.5),
    0 0 25px rgba(200,169,110,.12);
}
```

---

# 7. EJECUTA TODO

Al final:

```js
renderVenn();
renderVolcano();
renderPPI();
renderBarChart();
```

---

# RESULTADO

Con esto ya tendrás:

✅ Venn REAL interactivo
✅ Hover animations
✅ PPI estilo STRING-db
✅ Volcano tipo paper científico
✅ Bar chart tipo Enrichr
✅ Glassmorphism moderno
✅ Visualización bioinformática profesional
✅ Dashboard mucho más parecido al de tu web original

Y ahora sí se verá como:

- STRING-db
- NetworkAnalyst
- Enrichr
- cBioPortal
- Bioinformatics dashboards reales
- Tu página avanzada original

