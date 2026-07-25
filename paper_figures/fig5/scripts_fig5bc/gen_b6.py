from playwright.sync_api import sync_playwright
SD = "/tmp/claude-1000/-mnt-d-Research-OpenTargetsTransfer/8787aebc-ee4a-4493-9615-88e540dcfa09/scratchpad"
E = "https://www.ensembl.org"; G = "ENSG00000151224"; R = "10:80273890-80276800"
def log(*a): print(*a, flush=True)
OFF = ["contig","seq","chr_band_core","age_of_base","alignment_compara_9593_constrained",
       "dna_align_core_alt_seq_mapping","gencode_primary","mane_select","nstd166",
       "regulatory_build","variation_set_ph_variants"]
def xhr(pg,u):
    try: return pg.evaluate("""async (u)=>{const r=await fetch(u,{headers:{'X-Requested-With':'XMLHttpRequest'},credentials:'same-origin'});return r.status;}""",u)
    except Exception as e: return repr(e)
with sync_playwright() as p:
    b=p.chromium.launch(headless=True,args=["--no-sandbox","--disable-dev-shm-usage"])
    ctx=b.new_context(viewport={"width":1900,"height":1300}); pg=ctx.new_page()
    pg.set_default_timeout(300000)
    pg.goto(E+"/Homo_sapiens/Info/Index",wait_until="domcontentloaded",timeout=300000)
    pg.wait_for_timeout(4000)
    try: pg.click("text=I Agree",timeout=8000)
    except Exception: pass
    for phase in (1,2):
        pg.goto(E+"/Homo_sapiens/Location/View?r="+R+";g="+G,wait_until="domcontentloaded",timeout=300000)
        pg.wait_for_timeout(16000)
        try: pg.wait_for_selector("img.imagemap",timeout=250000)
        except Exception as e: log("imgwait",phase,str(e)[:90])
        pg.wait_for_timeout(6000)
        if phase==1:
            base=E+"/Homo_sapiens/Config/Location/ViewBottom?db=core;g="+G+";r="+R+";submit=1;"
            log("cfg",xhr(pg,base+"transcript_core_ensembl=transcript_label"))
            for k in OFF: xhr(pg,base+k+"=off")
    f={"filename":"MAT1A_zoom3.svg","format":"custom","image_format":"svg","resize":"","scale":"",
       "r":R,"data_type":"Location","extra":'{"highlightedTracks":[]}',"component":"ViewBottom",
       "data_action":"View","strain":"0","decodeURL":"1","db":"core","submit":"Download"}
    resp=pg.request.post(E+"/Homo_sapiens/ImageExport/ImageOutput",form=f,timeout=250000)
    bd=resp.body(); log("svg",resp.status,len(bd))
    if len(bd)>2000: open(SD+"/b6_zoom.svg","wb").write(bd)
    b.close()
