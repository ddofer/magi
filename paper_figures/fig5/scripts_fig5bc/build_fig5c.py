#!/usr/bin/env python3
"""Rebuild Fig. 5c (Ensembl 'Transcript comparison' style alignment) from
Ensembl REST data, using Ensembl's own variant-consequence colour scheme.

Needed because Ensembl release 116 no longer renders per-transcript variant
markup for the newer MAT1A transcripts on the live Transcript-comparison page.
"""
import json, html

SD = "/tmp/claude-1000/-mnt-d-Research-OpenTargetsTransfer/8787aebc-ee4a-4493-9615-88e540dcfa09/scratchpad/"

WIN_HI, WIN_LO = 80274670, 80274551          # genomic, inclusive; row = HI -> LO (minus strand)
ROW_START = 15481                            # Ensembl gene-relative coordinate of WIN_HI
GENE5 = 80290150                             # gene 5' end on the - strand (release 116)
VARIANT = 80274599

ORDER = [("MAT1A-201", "ENST00000372213"),
         ("MAT1A-205", "ENST00000871619"),
         ("MAT1A-206", "ENST00000871620"),
         ("MAT1A-207", "ENST00000871621"),
         ("MAT1A-208", "ENST00000871622")]

# Ensembl variant-consequence key: term -> (background, text colour)
CONS = {
 'splice_acceptor_variant':              ('#FF581A', '#ffffff'),
 'splice_donor_variant':                 ('#FF581A', '#ffffff'),
 'stop_gained':                          ('#ff0000', '#ffffff'),
 'frameshift_variant':                   ('#9400D3', '#ffffff'),
 'stop_lost':                            ('#ff0000', '#ffffff'),
 'start_lost':                           ('#ffd700', '#000000'),
 'inframe_insertion':                    ('#ff69b4', '#ffffff'),
 'inframe_deletion':                     ('#ff69b4', '#ffffff'),
 'missense_variant':                     ('#ffd700', '#000000'),
 'protein_altering_variant':             ('#ffd700', '#000000'),
 'splice_donor_5th_base_variant':        ('#ff7f50', '#ffffff'),
 'splice_region_variant':                ('#ff7f50', '#ffffff'),
 'splice_donor_region_variant':          ('#ff7f50', '#ffffff'),
 'splice_polypyrimidine_tract_variant':  ('#ff7f50', '#ffffff'),
 'incomplete_terminal_codon_variant':    ('#76ee00', '#000000'),
 'start_retained_variant':               ('#76ee00', '#000000'),
 'stop_retained_variant':                ('#76ee00', '#000000'),
 'synonymous_variant':                   ('#76ee00', '#000000'),
 'coding_sequence_variant':              ('#458b00', '#ffffff'),
 '5_prime_UTR_variant':                  ('#7ac5cd', '#000000'),
 '3_prime_UTR_variant':                  ('#7ac5cd', '#000000'),
 'non_coding_transcript_exon_variant':   ('#458b00', '#ffffff'),
 'intron_variant':                       ('#02599c', '#ffffff'),
}
# Ensembl consequence severity, most severe first (subset that can occur here)
RANK = ['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained',
        'frameshift_variant','stop_lost','start_lost','transcript_amplification',
        'inframe_insertion','inframe_deletion','missense_variant','protein_altering_variant',
        'splice_donor_5th_base_variant','splice_region_variant','splice_donor_region_variant',
        'splice_polypyrimidine_tract_variant','incomplete_terminal_codon_variant',
        'start_retained_variant','stop_retained_variant','synonymous_variant',
        'coding_sequence_variant','mature_miRNA_variant','5_prime_UTR_variant',
        '3_prime_UTR_variant','non_coding_transcript_exon_variant','intron_variant',
        'NMD_transcript_variant','non_coding_transcript_variant',
        'upstream_gene_variant','downstream_gene_variant','intergenic_variant']
RANKI = {t: i for i, t in enumerate(RANK)}

# Ensembl exon/intron key
FG_TRANSLATED, FG_UTR, FG_NONCODING, FG_INTRON = '#1044ee', '#cd6839', '#333333', '#aaaaaa'

gene = json.load(open(SD + 'mat1a.json'))
tx = {t['id']: t for t in gene['Transcript']}
seq = open(SD + 'win_seq.txt').read().strip()          # WIN_HI -> WIN_LO, minus strand
assert len(seq) == WIN_HI - WIN_LO + 1, (len(seq),)

overlap = json.load(open(SD + 'win_vars.json'))
vep = {v['id']: v for v in json.load(open(SD + 'vep_win.json'))}


def base_state(tid, pos):
    """Ensembl exon/intron colour for one genomic base of one transcript."""
    t = tx[tid]
    inexon = any(e['start'] <= pos <= e['end'] for e in t['Exon'])
    if not inexon:
        return FG_INTRON
    tr = t.get('Translation')
    if tr and tr['start'] <= pos <= tr['end']:
        return FG_TRANSLATED
    return FG_UTR if tr else FG_NONCODING


# per transcript, per genomic position -> most severe consequence term
mark = {tid: {} for _, tid in ORDER}
used, skipped = 0, []
for v in overlap:
    vid = v['id']
    rec = vep.get(vid)
    if not rec or 'transcript_consequences' not in rec:
        # HGMD_MUTATION features: VEP returns no per-transcript block. Fall back to the
        # region-overlap consequence, resolved per transcript by exon/intron state
        # (this is what Ensembl's own page shows for these features).
        skipped.append(vid)
        ct = v.get('consequence_type')
        if ct in RANKI:
            for _, tid in ORDER:
                for pos in range(max(v['start'], WIN_LO), min(v['end'], WIN_HI) + 1):
                    term = ct if base_state(tid, pos) != FG_INTRON else 'intron_variant'
                    cur = mark[tid].get(pos)
                    if cur is None or RANKI[term] < RANKI[cur]:
                        mark[tid][pos] = term
        continue
    used += 1
    for tc in rec['transcript_consequences']:
        tid = tc.get('transcript_id')
        if tid not in mark:
            continue
        terms = [t for t in tc.get('consequence_terms', []) if t in RANKI]
        if not terms:
            continue
        best = min(terms, key=lambda t: RANKI[t])
        for pos in range(max(v['start'], WIN_LO), min(v['end'], WIN_HI) + 1):
            cur = mark[tid].get(pos)
            if cur is None or RANKI[best] < RANKI[cur]:
                mark[tid][pos] = best

print('variants in window: %d, with transcript consequences: %d, skipped(no VEP tx data): %s'
      % (len(overlap), used, skipped))

# ---------------------------------------------------------------- render
rows_html = []
for name, tid in ORDER:
    cells = []
    for i, ch in enumerate(seq):
        pos = WIN_HI - i
        term = mark[tid].get(pos)
        if term and term in CONS:
            bg, fg = CONS[term]
            cells.append('<span style="background-color:%s;color:%s;text-decoration:underline" '
                         'title="%d %s">%s</span>' % (bg, fg, pos, term, ch))
        else:
            cells.append('<span style="color:%s">%s</span>' % (base_state(tid, pos), ch))
    rows_html.append('<span class="tn">%-9s</span> <span class="num">%5d</span> %s <span class="num">%5d</span>'
                     % (name, ROW_START, ''.join(cells), ROW_START + len(seq) - 1))

KEY = [('3 prime UTR', '3_prime_UTR_variant'), ('5 prime UTR', '5_prime_UTR_variant'),
       ('Coding sequence', 'coding_sequence_variant'), ('Frameshift', 'frameshift_variant'),
       ('Inframe deletion', 'inframe_deletion'), ('Intronic', 'intron_variant'),
       ('Missense', 'missense_variant'), ('Splice acceptor', 'splice_acceptor_variant'),
       ('Splice donor', 'splice_donor_variant'), ('Splice region', 'splice_region_variant'),
       ('Start lost', 'start_lost'), ('Stop gained', 'stop_gained'),
       ('Stop lost', 'stop_lost'), ('Synonymous', 'synonymous_variant')]
key_html = ''.join('<li><span class="chip" style="background-color:%s;color:%s">%s</span></li>'
                   % (CONS[t][0], CONS[t][1], html.escape(lab)) for lab, t in KEY)

doc = """<!doctype html><html><head><meta charset="utf-8"><style>
html,body{margin:0;padding:0;background:#fff;}
#FIG5C{display:inline-block;background:#fff;padding:14px 18px;}
.key{background:#f0f0f0;padding:7px 9px;margin:0 0 14px 0;font:11px Helvetica,Arial,sans-serif;
     display:block;overflow:hidden;max-width:700px;}
.key .lbl{float:left;font-weight:bold;color:#666;margin:3px 12px 0 0;}
.key ul{list-style:none;margin:0;padding:0;overflow:hidden;}
.key li{float:left;margin:0 4px 4px 0;}
.key .chip{display:inline-block;padding:2px 6px;}
pre.seq{margin:0;font:14px/1.45 "DejaVu Sans Mono","Courier New",monospace;white-space:pre;}
pre.seq .tn{color:#3b74b3;text-decoration:underline;}
pre.seq .num{color:#666;}
</style></head><body>
<div id="FIG5C">
  <div class="key"><span class="lbl">Variants</span><ul>%s</ul></div>
  <pre class="seq">%s</pre>
</div></body></html>""" % (key_html, '\n'.join(rows_html))

open(SD + 'fig5c_rebuilt.html', 'w').write(doc)
col = GENE5 - VARIANT + 1 - ROW_START + 1
print('wrote fig5c_rebuilt.html; variant column = %d; seq[col-1] = %s' % (col, seq[col - 1]))

# ---------------------------------------------------------------- SVG output
CW, LH, FS = 8.4, 20.0, 14.0            # char advance, line height, font size
PADX, PADY = 18.0, 14.0
LBL = 17                                # "MAT1A-201 15481 " prefix, in characters
KEY_H = 51.0

nrow, ncol = len(ORDER), len(seq)
seq_w = (LBL + ncol + 6) * CW
W = PADX * 2 + max(seq_w, 900)
H = PADY * 2 + KEY_H + 14 + nrow * LH

out = ['<?xml version="1.0" encoding="UTF-8"?>',
       '<svg xmlns="http://www.w3.org/2000/svg" width="%.1f" height="%.1f" viewBox="0 0 %.1f %.1f">'
       % (W, H, W, H),
       '<style>.m{font-family:"DejaVu Sans Mono","Courier New",monospace;font-size:%.1fpx;}'
       '.s{font-family:Helvetica,Arial,sans-serif;font-size:11px;}</style>' % FS,
       '<rect x="0" y="0" width="%.1f" height="%.1f" fill="#ffffff"/>' % (W, H)]

# legend
KEY_W = 0
for _lab, _t in KEY[:7]:
    KEY_W += 6.3 * len(_lab) + 16
KEY_W += 64
out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="#f0f0f0"/>'
           % (PADX, PADY, KEY_W, KEY_H))
out.append('<text class="s" x="%.1f" y="%.1f" font-weight="bold" fill="#666">Variants</text>'
           % (PADX + 8, PADY + 16))
for ri, chunk in enumerate((KEY[:7], KEY[7:])):
    x = PADX + 62 if ri == 0 else PADX + 8
    yb = PADY + 7 + ri * 21
    for lab, term in chunk:
        bg, fg = CONS[term]
        w = 6.3 * len(lab) + 12
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="17" fill="%s"/>' % (x, yb, w, bg))
        out.append('<text class="s" x="%.1f" y="%.1f" fill="%s">%s</text>'
                   % (x + 6, yb + 12.5, fg, html.escape(lab)))
        x += w + 4

y0 = PADY + KEY_H + 14
for r, (name, tid) in enumerate(ORDER):
    ybase = y0 + r * LH + FS
    out.append('<text class="m" x="%.1f" y="%.1f" fill="#3b74b3" '
               'text-decoration="underline">%s</text>' % (PADX, ybase, name))
    out.append('<text class="m" x="%.1f" y="%.1f" fill="#666">%d</text>'
               % (PADX + 10 * CW, ybase, ROW_START))
    glyphs = []
    for i, ch in enumerate(seq):
        pos = WIN_HI - i
        cx = PADX + (LBL + i) * CW
        term = mark[tid].get(pos)
        if term and term in CONS:
            bg, fg = CONS[term]
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.1f" fill="%s"/>'
                       % (cx, ybase - FS + 3.5, CW + 0.05, LH - 2, bg))
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="0.9" fill="%s"/>'
                       % (cx, ybase + 1.8, CW + 0.05, fg))
            glyphs.append((cx, fg, ch))
        else:
            glyphs.append((cx, base_state(tid, pos), ch))
    for cx, fg, ch in glyphs:
        out.append('<text class="m" x="%.2f" y="%.2f" fill="%s">%s</text>' % (cx, ybase, fg, ch))
    out.append('<text class="m" x="%.1f" y="%.1f" fill="#666">%d</text>'
               % (PADX + (LBL + ncol + 1) * CW, ybase, ROW_START + ncol - 1))
out.append('</svg>')
open(SD + 'fig5c_rebuilt.svg', 'w').write('\n'.join(out))
print('wrote fig5c_rebuilt.svg  %.0f x %.0f' % (W, H))
