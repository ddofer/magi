import json
from playwright.sync_api import sync_playwright
SD="/tmp/claude-1000/-mnt-d-Research-OpenTargetsTransfer/8787aebc-ee4a-4493-9615-88e540dcfa09/scratchpad"
E="https://www.ensembl.org"; G="ENSG00000151224"
def log(*a): print(*a, flush=True)
OFF=["contig","seq","chr_band_core","age_of_base","alignment_compara_9593_constrained",
     "dna_align_core_alt_seq_mapping","gencode_primary","mane_select","nstd166",
     "regulatory_build","variation_set_ph_variants"]
REGIONS=[("full","10:80269780-80290126"),("zoom","10:80274100-80277000")]

def xhr(pg,u):
    try:
        return pg.evaluate("""async (u)=>{const r=await fetch(u,{headers:{'X-Requested-With':'XMLHttpRequest'},credentials:'same-origin'});return r.status;}""",u)
    except Exception as e:
        return repr(e)

def load(pg,url,tag):
    pg.goto(url,wait_until="domcontentloaded",timeout=300000)
    pg.wait_for_timeout(18000)
    try: pg.wait_for_selector("img.imagemap",timeout=200000)
    except Exception as e: log("  imgwait",tag,str(e)[:120])
    pg.wait_for_timeout(8000)

def export(pg,r,fmt,name):
    f={"filename":name,"format":"custom","image_format":fmt,"resize":"",
       "scale":"5" if fmt=="png" else "",
       "r":r,"data_type":"Location","extra":'{"highlightedTracks":[]}',
       "component":"ViewBottom","data_action":"View","strain":"0",
       "decodeURL":"1","db":"core","submit":"Download"}
    resp=pg.request.post(f"{E}/Homo_sapiens/ImageExport/ImageOutput",form=f,timeout=250000)
    bd=resp.body(); log("   ",fmt,resp.status,len(bd)); return bd

with sync_playwright() as p:
    b=p.chromium.launch(headless=True,args=["--no-sandbox","--disable-dev-shm-usage"])
    ctx=b.new_context(viewport={"width":1900,"height":1300})
    pg=ctx.new_page(); pg.set_default_timeout(300000)
    pg.goto(f"{E}/Homo_sapiens/Info/Index",wait_until="domcontentloaded",timeout=300000)
    pg.wait_for_timeout(4000)
    try: pg.click("text=I Agree",timeout=8000)
    except Exception: pass
    for tag,r in REGIONS:
        log("===",tag,r)
        load(pg,f"{E}/Homo_sapiens/Location/View?r={r};g={G}","pre")
        t=pg.evaluate("()=>document.body.innerText")
        log("  after g-load Comprehensive?", "Comprehensive" in t)
        base=f"{E}/Homo_sapiens/Config/Location/ViewBottom?db=core;g={G};r={r};submit=1;"
        log("  cfg transcript:", xhr(pg, base+"transcript_core_ensembl=transcript_label"))
        for k in OFF: xhr(pg, base+k+"=off")
        load(pg,f"{E}/Homo_sapiens/Location/View?r={r};g={G}","post")
        t=pg.evaluate("()=>document.body.innerText")
        log("  Comprehensive?","Comprehensive" in t,"| MANE?","MANE Select" in t,
            "| gnomAD?","gnomAD variants" in t,"| GERP?","GERP" in t)
        pg.screenshot(path=SD+"/shots/b4_"+tag+".png",full_page=True)
        for fmt in ["svg","pdf","png"]:
            try:
                bd=export(pg,r,fmt,"MAT1A_"+tag+"."+fmt)
                if len(bd)>2000: open(SD+"/b4_"+tag+"."+fmt,"wb").write(bd)
            except Exception as e: log("   FAIL",fmt,repr(e))
    b.close()
