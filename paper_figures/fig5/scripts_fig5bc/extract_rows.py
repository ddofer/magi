#!/usr/bin/env python3
"""Extract selected transcript rows from an Ensembl 'Region in detail' SVG export
and restack them into a compact, print-ready vector panel.

Elements are moved with <g transform="translate(0,dy)"> so path data
(which mixes absolute M and relative l commands) is never rewritten.
"""
import re, sys

NUM = r'-?\d+(?:\.\d+)?'
PAIR = re.compile(r'([MmLl])\s*(%s),(%s)((?:\s+%s,%s)*)' % (NUM, NUM, NUM, NUM))


def parse_svg(path):
    s = open(path, encoding='latin-1').read()
    m = re.search(r'<svg\s+width="([\d.]+)"\s+height="([\d.]+)"', s)
    W, H = float(m.group(1)), float(m.group(2))
    style = re.search(r'<defs>.*?</defs>', s, re.S).group(0)
    body = s[s.index('</defs>') + len('</defs>'):s.rindex('</svg>')]
    els = re.findall(r'<(?:rect|text|path|line|poly|polygon)\b.*?(?:/>|</(?:text|path|poly|polygon|line|rect)>)',
                     body, re.S)
    return W, H, style, els


def path_ys(d):
    """Return absolute y coordinates traced by a path 'd' string."""
    ys, cy = [], 0.0
    for m in re.finditer(r'([MmLlHhVv])((?:\s*%s\s*,\s*%s)+|\s*%s)' % (NUM, NUM, NUM), d):
        cmd, args = m.group(1), m.group(2)
        pts = re.findall(r'(%s)\s*,\s*(%s)' % (NUM, NUM), args)
        for _, y in pts:
            y = float(y)
            cy = y if cmd in 'ML' else cy + y
            ys.append(cy)
    return ys


def yrange(e):
    t = re.match(r'<(\w+)', e).group(1)
    if t == 'rect':
        y = float(re.search(r'\by="(%s)"' % NUM, e).group(1))
        h = float(re.search(r'\bheight="(%s)"' % NUM, e).group(1))
        return y, y + h
    if t == 'text':
        y = float(re.search(r'\by="(%s)"' % NUM, e).group(1))
        return y - 7, y + 4
    d = re.search(r'\bd="([^"]*)"', e)
    if d:
        ys = path_ys(d.group(1))
        return (min(ys), max(ys)) if ys else None
    p = re.search(r'\bpoints="([^"]*)"', e)
    if p:
        ys = [float(b) for _, b in re.findall(r'(%s)\s*,\s*(%s)' % (NUM, NUM), p.group(1))]
        return (min(ys), max(ys)) if ys else None
    return None


def rect_wh(e):
    return (float(re.search(r'\bwidth="(%s)"' % NUM, e).group(1)),
            float(re.search(r'\bheight="(%s)"' % NUM, e).group(1)))


def split_multiline_text(e, line_height=13.0):
    m = re.match(r'(<text[^>]*>)(.*?)(</text>)', e, re.S)
    if not m:
        return [e]
    head, txt, tail = m.groups()
    parts = txt.split('\n')
    if len(parts) == 1:
        return [e]
    ty = float(re.search(r'\by="(%s)"' % NUM, head).group(1))
    out = []
    for j, part in enumerate(parts):
        h2 = re.sub(r'\by="(%s)"' % NUM, 'y="%.3f"' % (ty + j * line_height), head, count=1)
        out.append(h2 + part + tail)
    return out


def build(src, targets, out, pitch=44.0, top_pad=8.0, bot_pad=10.0,
          keep_gridlines=True, drop_highlight=True, label_line_height=13.0,
          x0=148.0):
    W, H, style, els = parse_svg(src)

    rows = []
    for name in targets:
        hit = next((e for e in els if e.startswith('<text') and name in e), None)
        if hit is None:
            raise SystemExit('label not found: %s' % name)
        ly = float(re.search(r'\by="(%s)"' % NUM, hit).group(1))
        rows.append((name, ly - 14.0, ly))
        print('  row %-34s glyph_y=%7.2f' % (name, ly - 14.0))

    grid = [e for e in els if e.startswith('<rect') and rect_wh(e)[1] > H * 0.5 and rect_wh(e)[0] <= 2]
    newH = top_pad + pitch * len(rows) + bot_pad
    parts = ['<rect x="0" y="0" width="%.0f" height="%.1f" fill="rgb(255,255,255)" stroke="none"/>' % (W + 10, newH)]

    if keep_gridlines:
        for e in grid:
            e2 = re.sub(r'\by="(%s)"' % NUM, 'y="0"', e, count=1)
            e2 = re.sub(r'\bheight="(%s)"' % NUM, 'height="%.1f"' % newH, e2, count=1)
            parts.append(e2)

    for i, (name, gy, ly) in enumerate(rows):
        lo, hi = gy - 3.0, gy + pitch - 3.0
        dy = (top_pad + i * pitch) - gy
        grp = []
        for e in els:
            if e in grid:
                continue
            r = yrange(e)
            if r is None:
                continue
            c = (r[0] + r[1]) / 2.0
            if not (lo <= c < hi):
                continue
            if e.startswith('<rect'):
                w, h = rect_wh(e)
                if drop_highlight and w > 400 and h > 20:
                    continue
            if e.startswith('<text'):
                grp.extend(split_multiline_text(e, label_line_height))
            else:
                grp.append(e)
        parts.append('<g id="%s" transform="translate(0,%.3f)">' %
                     (re.search(r'(MAT1A-\d+)', name).group(1), dy))
        parts.extend(grp)
        parts.append('</g>')

    Wc = W - x0
    hdr = ('<?xml version="1.0" encoding="UTF-8"?>\n'
           '<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" '
           'width="%.1f" height="%.1f" viewBox="%.1f 0 %.1f %.1f">\n' % (Wc, newH, x0, Wc, newH))
    open(out, 'w', encoding='latin-1').write(hdr + style + '\n' + '\n'.join(parts) + '\n</svg>\n')
    print('  wrote %s  %.1f x %.1f (viewBox x0=%.0f)' % (out, Wc, newH, x0))
    return Wc, newH


if __name__ == '__main__':
    src, out = sys.argv[1], sys.argv[2]
    targets = ['MAT1A-205 - ENST00000871619',
               'MAT1A-208 - ENST00000871622',
               'MAT1A-207 - ENST00000871621',
               'MAT1A-201 - ENST00000372213']
    build(src, targets, out)
