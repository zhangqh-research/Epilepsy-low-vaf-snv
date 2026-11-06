"""Caller and strategy performance against successful spike-in SNVs."""
import argparse
import csv
import gzip
import statistics
from collections import Counter
from decimal import Decimal
from pathlib import Path
CALLERS=['Mutect2','Strelka2','MuSE','FreeBayes']
MS='Mutect2 + Strelka2 Union'
MSE='Mutect2 + Strelka2 + MuSE Union'
MSF='Mutect2 + Strelka2 + FreeBayes Union'
FOUR='Four-caller Union'
CON='Four-caller Consensus (>=2)'
NAMES=CALLERS+[MS,MSE,MSF,FOUR,CON]
BINS=[('0%-<2%',Decimal('0'),Decimal('0.02')),('2%-<4%',Decimal('0.02'),Decimal('0.04'))]
def read(p):
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def key(r):return (r['pair_id'],r['chrom'],int(r['pos']),r['ref'],r['alt'])


def summary(values):
 vals=[x for x in values if x is not None]
 assert vals
 q1,_,q3=statistics.quantiles(vals,n=4,method='inclusive')
 return dict(evaluable_pairs=len(vals),median=statistics.median(vals),q1=q1,q3=q3,min=min(vals),max=max(vals))


def stratify(sets):
 out=dict(sets)
 out[MS]=sets['Mutect2']|sets['Strelka2']
 out[MSE]=out[MS]|sets['MuSE']
 out[MSF]=out[MS]|sets['FreeBayes']
 out[FOUR]=set.union(*sets.values())
 counts=Counter(k for s in sets.values() for k in s)
 out[CON]={k for k,n in counts.items() if n>=2}
 return out


def metrics(calls,truth,base):
 tp=len(calls&truth);fn=len(truth-calls);fp=len(calls-(truth|base));excluded=len((calls&base)-truth)
 assert len(calls)==tp+fp+excluded
 return dict(successful_truth=len(truth),retained_calls=len(calls),strategy_baseline=len(base),
             TP=tp,FP_baseline_adjusted=fp,FN=fn,excluded_baseline_calls=excluded,
             recall=tp/len(truth),miss_rate=fn/len(truth),
             baseline_adjusted_precision=tp/(tp+fp) if tp+fp else None,
             baseline_adjusted_F1=2*tp/(2*tp+fp+fn) if 2*tp+fp+fn else None)

def write(path, rows):
    with path.open('x', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader()
        w.writerows(rows)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--truth', type=Path, required=True)
    p.add_argument('--calls', type=Path, required=True)
    p.add_argument('--baseline', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(exist_ok=False)
    truth = {key(r): r for r in read(args.truth)}
    assert len(truth) == 48995
    pairs = sorted({k[0] for k in truth})
    assert len(pairs) == 10
    cs = {c: set() for c in CALLERS}
    bs = {c: set() for c in CALLERS}
    for src, sets in [(args.calls, cs), (args.baseline, bs)]:
        for r in read(src):
            k, caller = key(r), r['caller']
            assert k[0] in pairs and k not in sets[caller]
            sets[caller].add(k)
    calls, baseline = stratify(cs), stratify(bs)
    truthset = set(truth)
    pooled, paired, low, complement = [], [], [], []
    for name in NAMES:
        pooled.append({'strategy': name, **metrics(calls[name], truthset, baseline[name])})
        for pair in pairs:
            pt = {k for k in truthset if k[0] == pair}
            pc = {k for k in calls[name] if k[0] == pair}
            pb = {k for k in baseline[name] if k[0] == pair}
            paired.append({'pair_id': pair, 'strategy': name, **metrics(pc, pt, pb)})
    for label, lo, hi in BINS:
        bt = {k for k, r in truth.items() if lo <= Decimal(r['observed_spike_vaf']) < hi}
        for name in NAMES:
            tp = len(calls[name] & bt)
            low.append({'VAF_bin': label, 'strategy': name, 'successful_truth': len(bt),
                        'TP': tp, 'FN': len(bt)-tp, 'recall': tp/len(bt)})
        m, st = cs['Mutect2'] & bt, cs['Strelka2'] & bt
        complement.append({'VAF_bin': label, 'Mutect2_not_Strelka2': len(m-st),
                           'Strelka2_not_Mutect2': len(st-m), 'shared': len(m & st),
                           'union_TP': len(m | st)})
    write(args.output/'pooled_metrics.tsv', pooled)
    write(args.output/'pair_metrics.tsv', paired)
    write(args.output/'low_VAF_recall.tsv', low)
    write(args.output/'caller_complementarity.tsv', complement)

if __name__ == '__main__':
    main()
