import json, time
from playwright.sync_api import sync_playwright
SD = "/tmp/claude-1000/-mnt-d-Research-OpenTargetsTransfer/8787aebc-ee4a-4493-9615-88e540dcfa09/scratchpad"
E = "https://www.ensembl.org"; G = "ENSG00000151224"
def log(*a): print(*a, flush=True)
TC_URL = (E + "/Homo_sapiens/Gene/TranscriptComparison?db=core;g=" + G +
          ";r=10:80269780-80290150;t1=ENST00000372213;t2=ENST00000871619;"
          "t3=ENST00000871620;t4=ENST00000871621;t5=ENST00000871622")

MK_RANGE = r"""
() => {
  const pres=[...document.querySelectorAll('pre.text_sequence')];
  const pre=pres.find(p=>p.textContent.includes('MAT1A-201  15481'));
  if(!pre) return null;
  const full=pre.textContent;
  const s=full.indexOf('MAT1A-201  15481');
  const e=full.indexOf('\n', full.indexOf('MAT1A-208  15481'));
  const w=document.createTreeWalker(pre,NodeFilter.SHOW_TEXT);
  let acc=0,sN=null,sO=0,eN=null,eO=0,n;
  while((n=w.nextNode())){const L=n.nodeValue.length;
    if(sN===null&&acc+L>s){sN=n;sO=s-acc;}
    if(eN===null&&acc+L>=e){eN=n;eO=e-acc;} acc+=L;}
  window.__r=document.createRange(); window.__r.setStart(sN,sO); window.__r.setEnd(eN,eO);
  const d=document.createElement('div'); d.appendChild(window.__r.cloneContents());
  const rows=d.innerHTML.split(/(?=<a href="\/Homo_sapiens\/Transcript\/Summary)/)
      .filter(x=>/>MAT1A-\d+<\/a>/.test(x))
      .map(x=>[(x.match(/>(MAT1A-\d+)<\/a>/)||[])[1], (x.match(/sequence_info/g)||[]).length]);
  const rect=window.__r.getBoundingClientRect();
  return {rows, absTop: rect.top + window.scrollY,
          styled: document.querySelectorAll('pre.text_sequence span.adorn [style*=background]').length};
}
"""

BUILD_JS = open(SD + "/build_js.txt").read() if False else r"""
() => {
  const d=document.createElement('div'); d.appendChild(window.__r.cloneContents());
  const key=[...document.querySelectorAll('._adornment_key')]
     .find(k=>[...k.querySelectorAll('dt')].some(t=>t.textContent.trim()==='Variants'));
  const dl=key.querySelector('dl'); const kids=[...dl.children];
  const keep=document.createElement('dl');
  for(let i=0;i<kids.length;i++){
    if(kids[i].tagName==='DT'&&kids[i].textContent.trim()==='Variants'){
      keep.appendChild(kids[i].cloneNode(true)); if(kids[i+1]) keep.appendChild(kids[i+1].cloneNode(true)); } }
  const keyDiv=document.createElement('div'); keyDiv.className='fig5c-key'; keyDiv.appendChild(keep);
  const st=document.createElement('style');
  st.textContent=`
    html,body{background:#fff !important;margin:0;padding:0;}
    #FIG5C{background:#fff;padding:12px 16px;display:inline-block;}
    #FIG5C .fig5c-key{background:#f0f0f0;padding:6px 8px;margin:0 0 12px 0;display:block;
                      font-family:Helvetica,Arial,sans-serif;font-size:11px;}
    #FIG5C .fig5c-key dl{margin:0;overflow:hidden;}
    #FIG5C .fig5c-key dt{float:left;font-weight:bold;color:#666;margin:3px 10px 0 0;}
    #FIG5C .fig5c-key dd{margin:0;overflow:hidden;}
    #FIG5C .fig5c-key ul{list-style:none;margin:0;padding:0;overflow:hidden;}
    #FIG5C .fig5c-key li{float:left;margin:0 4px 4px 0;}
    #FIG5C .fig5c-key .adorn-key-entry{padding:2px 6px;display:inline-block;}
    #FIG5C pre.text_sequence{margin:0;background:#fff;}`;
  const box=document.createElement('div'); box.id='FIG5C';
  const npre=document.createElement('pre'); npre.className='text_sequence';
  while(d.firstChild) npre.appendChild(d.firstChild);
  box.appendChild(keyDiv); box.appendChild(npre);
  document.head.appendChild(st);
  document.body.insertBefore(box, document.body.firstChild);
  [...document.body.children].forEach(c=>{if(c.id!=='FIG5C') c.style.display='none';});
  window.scrollTo(0,0);
  const b=box.getBoundingClientRect();
  return {w:b.width,h:b.height,rows:npre.textContent.split('\n').map(l=>l.slice(0,22))};
}
"""

with sync_playwright() as p:
    b = p.chromium.launch(headless=True, args=["--no-sandbox", "--disable-dev-shm-usage"])
    ctx = b.new_context(viewport={"width": 2400, "height": 1600}, device_scale_factor=6)
    pg = ctx.new_page(); pg.set_default_timeout(300000)
    pg.goto(E + "/Homo_sapiens/Info/Index", wait_until="domcontentloaded", timeout=300000)
    pg.wait_for_timeout(4000)
    try: pg.click("text=I Agree", timeout=8000)
    except Exception: pass
    pg.goto(TC_URL, wait_until="domcontentloaded", timeout=300000)
    pg.wait_for_timeout(25000)
    try: pg.wait_for_selector("pre.text_sequence", timeout=280000)
    except Exception as ex: log("wait", str(ex)[:100])
    st = pg.evaluate(MK_RANGE); log("initial", json.dumps(st)[:400])
    target = st["absTop"]
    docH = pg.evaluate("()=>document.body.scrollHeight")
    log("docH", docH, "targetY", target)

    # walk the viewport down to the target to trigger viewport-driven adornment
    y = 0
    while y < target:
        y += 1000
        pg.evaluate("(y)=>window.scrollTo(0,y)", min(y, target))
        pg.wait_for_timeout(250)
    pg.evaluate("(y)=>window.scrollTo(0,y-400)", target)
    log("scrolled to target")

    t0 = time.time(); best = None
    while time.time() - t0 < 600:
        pg.wait_for_timeout(15000)
        st = pg.evaluate(MK_RANGE)
        log("  t=%4.0f" % (time.time()-t0), json.dumps(st["rows"]), "styled", st["styled"])
        ok = st["rows"] and all(c > 5 for _, c in st["rows"])
        best = st
        if ok:
            log("ALL ROWS ADORNED"); break

    res = pg.evaluate(BUILD_JS); log("BUILD", json.dumps(res)[:600])
    pg.wait_for_timeout(4000)
    el = pg.query_selector("#FIG5C"); bb = el.bounding_box()
    el.screenshot(path=SD + "/fig5c_x6.png")
    pg.pdf(path=SD + "/fig5c.pdf", print_background=True,
           width=str(int(bb["width"]) + 6) + "px", height=str(int(bb["height"]) + 6) + "px",
           margin={"top": "0", "bottom": "0", "left": "0", "right": "0"})
    open(SD + "/fig5c_block.html", "w").write(pg.evaluate("()=>document.getElementById('FIG5C').outerHTML"))
    log("done", bb)
    b.close()
