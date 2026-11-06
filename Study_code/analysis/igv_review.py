"""Stratified random sampling and descriptive summaries of manual IGV review."""
import argparse
import csv
import math
from collections import Counter
from pathlib import Path
import numpy as np
GROUPS = ['Mutect2-only','Strelka2-only','MuSE-only','FreeBayes-only','Multi-caller consensus']
OUTCOMES = ['Supportive','Indeterminate','Not supportive']
SEED = 20261027

def load(path):
    with path.open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(path, rows):
    with path.open('x', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def sample(events):
    assert len(events) == 4541 and len({r['event_id'] for r in events}) == 4541
    assert len({r['pair_id'] for r in events}) == 83
    frame = sorted(events, key=lambda r:(r['pair_id'],r['chrom'],int(r['pos']),r['ref'],r['alt']))
    rng=np.random.default_rng(SEED);selected=[]
    for g in GROUPS:
        pool=[r for r in frame if r['group']==g]
        selected.extend(dict(pool[int(i)]) for i in rng.choice(len(pool),size=20,replace=False))
    selected=[selected[int(i)] for i in rng.permutation(100)]
    return [{'review_id': f'IGV{i:03d}', **r} for i, r in enumerate(selected, 1)]

def number(x):
 if x is None or str(x).strip() in ('','NA','nan','None'):return None
 v=float(x);return v if math.isfinite(v) else None


def fraction(n,d):return n/d if d else None


def base_summary(rows):
 c=Counter(r['overall_IGV'] for r in rows);n=len(rows)
 return {'n':n,'patients':len({r['pair_id'] for r in rows}),
  'unique_genomic_SNVs':len({(r['chrom'],r['pos'],r['ref'],r['alt']) for r in rows}),
  **{s:c[s] for s in OUTCOMES},**{s+'_fraction':fraction(c[s],n) for s in OUTCOMES}}


def qstats(values):
 vals=[number(v) for v in values];vals=[v for v in vals if v is not None]
 if not vals:return {'valid_n':0,'median':None,'q1':None,'q3':None,'min':None,'max':None}
 q=np.quantile(vals,[.25,.5,.75],method='linear')
 return {'valid_n':len(vals),'median':float(q[1]),'q1':float(q[0]),'q3':float(q[2]),'min':min(vals),'max':max(vals)}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('mode', choices=['sample', 'summarize'])
    p.add_argument('input', type=Path)
    p.add_argument('output', type=Path)
    args = p.parse_args()
    rows = load(args.input)
    if args.mode == 'sample':
        write(args.output, sample(rows))
        return
    assert len(rows) == 100
    assert all(r['overall_IGV'] in OUTCOMES for r in rows)
    assert all(sum(r['group'] == g for r in rows) == 20 for g in GROUPS)
    results = [{'group': g, **base_summary([r for r in rows if r['group'] == g])} for g in GROUPS]
    write(args.output, results)

if __name__ == '__main__':
    main()
