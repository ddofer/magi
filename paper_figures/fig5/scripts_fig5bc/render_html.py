import sys
from playwright.sync_api import sync_playwright
src, png, pdf, dsf = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
with sync_playwright() as p:
    b=p.chromium.launch(headless=True,args=["--no-sandbox","--disable-dev-shm-usage"])
    ctx=b.new_context(viewport={"width":2600,"height":1000}, device_scale_factor=dsf)
    pg=ctx.new_page()
    pg.goto("file://"+src, wait_until="load"); pg.wait_for_timeout(1500)
    el=pg.query_selector("#FIG5C"); bb=el.bounding_box(); print("bbox", bb)
    el.screenshot(path=png)
    pg.pdf(path=pdf, print_background=True,
           width=str(int(bb['width'])+6)+"px", height=str(int(bb['height'])+6)+"px",
           margin={"top":"0","bottom":"0","left":"0","right":"0"})
    b.close()
