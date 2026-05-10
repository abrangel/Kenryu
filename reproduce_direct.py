
import asyncio
from demo.backend.main import analyze, PipelineRequest, TARGETSCAN_DB
import os

async def reproduce():
    print("Verificando estabilidad con 1 solo miRNA...")
    from demo.backend import main
    main.TARGETSCAN_DB = {"mir-33": {"ABCA1", "SCN1A"}}
    
    request = PipelineRequest(mirnas=["hsa-miR-33a-5p"], years=10, mode="strict")
    
    try:
        print("Llamando a analyze con 1 miRNA...")
        result = await analyze(request)
        print("ÉXITO: El sistema manejó correctamente el caso de 1 miRNA.")
        if "venn_plot" in result:
            print("Gráfico de fallback generado correctamente.")
    except Exception as e:
        print(f"FALLO DETECTADO: {type(e).__name__} - {e}")

if __name__ == "__main__":
    asyncio.run(reproduce())
