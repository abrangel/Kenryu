import gseapy as gp
import pandas as pd

def test_enrich():
    genes = ["ABCA1", "KPNA3", "SCN1A", "TSC22D2", "SNTB2"]
    print(f"Probando enriquecimiento para: {genes}")
    try:
        enr = gp.enrichr(gene_list=genes,
                         gene_sets=['KEGG_2021_Human', 'GO_Biological_Process_2021'],
                         organism='Human',
                         outdir=None)
        results = enr.results
        print("Columnas encontradas:", results.columns.tolist())
        
        if results.empty:
            print("No se encontraron resultados.")
            return
            
        # Verificar si existe la columna de p-value ajustado
        col = 'Adjusted P-value'
        if col not in results.columns:
            # En algunas versiones de gseapy o bases de datos el nombre puede variar
            print("Buscando columna alternativa para P-value...")
            alternate = [c for c in results.columns if 'P-value' in c]
            print(f"Columnas candidatas: {alternate}")
        else:
            significant = results[results[col] < 0.05]
            print(f"Encontrados {len(significant)} resultados significativos.")
            
    except Exception as e:
        print(f"ERROR DETECTADO: {str(e)}")

if __name__ == "__main__":
    test_enrich()
