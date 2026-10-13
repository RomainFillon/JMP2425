"""releve_pa1.py -- releve du point d'arret 1 (lot 51, course du 28/09/2026).

Compare l'etat courant de output/ et intermediate_outputs/ a l'archive d'avant la
course (JMP/_ARCHIVE_avant_course_20260928) et ecrit le releve dans un fichier.
Il n'interprete rien et ne corrige rien : il mesure et evalue les conditions 2-4
fixees par l'auteur (author's decision ; la condition 1, les neuf controles, est lancee a part par
tools/reference_phase.py et son resultat est joint au releve par --controles).

  macros \\EMP*        : avant / apres, ecart relatif, changement de signe
  paires Lo/Hi         : celles qui EXCLUAIENT zero dans l'archive gardent-elles
                         le signe ? (macros, puis IC des CSV : ci_lo/ci_hi,
                         ci95_lo/ci95_hi, diff_lo/diff_hi)
  tests joints         : verdict a 5 % par ligne, colonnes p_* / pval_*
  effectifs            : nobs, n_obs*, Nobs des macros
  pays                 : entrees / sorties dans ERA_tas_annual.csv
  jours changeant de bin : 1/2 somme |count_apres - count_avant| / somme count,
                         par pays-annee, dans ERA_tas_annual.csv

Usage : python tools/releve_pa1.py SORTIE.txt [--controles FICHIER]
Code  : 0 si les conditions 2, 3 et 4 passent, 1 sinon.
"""
import glob
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ARCH = os.path.join(os.path.dirname(os.path.dirname(ROOT)), '_ARCHIVE_avant_course_20260928')
SORTIE = sys.argv[1]
CTRL = sys.argv[sys.argv.index('--controles') + 1] if '--controles' in sys.argv else None

L = []
def w(s=''):
    L.append(s)

RE_DEF = re.compile(r'\\newcommand\{\\([A-Za-z]+)\}\{(.*?)\}\s*(?:%|$)', re.M)


def macros(p):
    return dict(RE_DEF.findall(open(p, encoding='utf-8').read()))


def num(s):
    t = s.replace('{,}', '').replace('\\,', '').replace('$', '').replace('\\%', '')
    t = t.replace('\\text{-}', '-').replace('--', '-').replace('−', '-').strip()
    try:
        return float(t)
    except ValueError:
        return None


def signe(x):
    return 0 if x is None or x == 0 else (1 if x > 0 else -1)


echecs = []

# ------------------------------- EN TETE (author's decision, lot 57) : les gardes qui ont refuse
# Une garde qui a refuse ne se lit pas enterree au milieu des ecarts de macros.
w('=' * 78)
w('0. GUARDS THAT REFUSED AND WERE REPORTED (declared resumption)')
w('   SENTENCE guards (emit_values_empirics.garde_phrase, emit_values) and frozen guards')
w('   (1.0052, POP_ATTENDUE); source: output/ecarts_reprise_gardes.tsv')
w('=' * 78)
pg = os.path.join(ROOT, 'output', 'ecarts_reprise_gardes.tsv')
if os.path.exists(pg):
    vues = set()
    for l in open(pg, encoding='utf-8'):
        c = l.rstrip('\n').split('\t')
        if c[0] == 'date' or len(c) < 5:
            continue
        cle = tuple(c[1:])          # une meme garde rapportee a chaque passage : une ligne
        if cle in vues:
            continue
        vues.add(cle)
        w('   %s | %s' % (c[1], c[2]))
        w('        expected: %s' % c[3])
        w('        obtained: %s' % c[4])
    w('   %d distinct guard(s)' % len(vues))
else:
    w('   (none)')
pe = os.path.join(ROOT, 'output', 'ecarts_reprise_empirics.tsv')
if os.path.exists(pe):
    n_e = sum(1 for l in open(pe, encoding='utf-8') if l.strip() and not l.startswith(('#', 'macro\t')))
    w('   (separately: %d \\EMP* macro differences, expected value written -> emitted, detail in %s)'
      % (n_e, os.path.relpath(pe, ROOT)))
w()

# ---------------------------------------------------------------- macros EMP
A = macros(os.path.join(ARCH, 'output', 'values_empirics.tex'))
N = macros(os.path.join(ROOT, 'output', 'values_empirics.tex'))
w('=' * 78)
w('1. \\EMP* MACROS: %d in the archive, %d after' % (len(A), len(N)))
w('=' * 78)
w('   gone: %s' % (sorted(set(A) - set(N)) or 'none'))
w('   new : %s' % (sorted(set(N) - set(A)) or 'none'))
n_id = n_mv = n_sg = 0
w('%-34s %-18s %-18s %10s  %s' % ('macro', 'before', 'after', 'rel. diff', 'sign'))
for k in sorted(set(A) & set(N)):
    a, n = A[k], N[k]
    xa, xn = num(a), num(n)
    if a == n:
        n_id += 1
        tag, er = '', ''
    else:
        n_mv += 1
        er = ('%+.2e' % ((xn - xa) / abs(xa))) if xa not in (None, 0) and xn is not None else 'n.d.'
        tag = ''
        if xa is not None and xn is not None and signe(xa) != signe(xn):
            tag = 'FLIPPED'
            n_sg += 1
    w('%-34s %-18s %-18s %10s  %s' % (k, a[:18], n[:18], er, tag))
w('   identical %d; changed %d; sign changes %d' % (n_id, n_mv, n_sg))

# ------------------------------------------------ condition 2 : paires Lo/Hi
w()
w('=' * 78)
w('2. Lo/Hi PAIRS OF THE MACROS: those that excluded zero in the archive')
w('=' * 78)
pairs = sorted(k[:-2] for k in A if k.endswith('Lo') and k[:-2] + 'Hi' in A)
exc_av = exc_ap = 0
for b in pairs:
    lo, hi = num(A[b + 'Lo']), num(A[b + 'Hi'])
    if lo is None or hi is None:
        w('   %-30s not numeric, set aside' % b)
        continue
    ex_a = signe(lo) == signe(hi) != 0
    lo2 = num(N.get(b + 'Lo', '')) if b + 'Lo' in N else None
    hi2 = num(N.get(b + 'Hi', '')) if b + 'Hi' in N else None
    ex_n = lo2 is not None and hi2 is not None and signe(lo2) == signe(hi2) != 0
    exc_av += ex_a
    exc_ap += ex_n
    if not ex_a:
        w('   %-30s [%s ; %s] did not exclude zero -- not kept' % (b, A[b + 'Lo'], A[b + 'Hi']))
        continue
    garde = ex_n and signe(lo2) == signe(lo)
    if not garde:
        echecs.append('condition 2: %s' % b)
    w('   %-30s KEPT     [%s ; %s] -> [%s ; %s]  %s' % (
        b, A[b + 'Lo'], A[b + 'Hi'], N.get(b + 'Lo', '?'), N.get(b + 'Hi', '?'),
        'sign kept' if garde else 'FAIL'))
w('   pairs excluding zero (macros): %d before, %d after' % (exc_av, exc_ap))


def csvs(base):
    out = {}
    for p in glob.glob(os.path.join(base, 'output', '**', '*.csv'), recursive=True):
        r = os.path.relpath(p, os.path.join(base, 'output')).replace('\\', '/')
        if 'past_outputs' in r or '_wide_' in r:
            continue
        out[r] = p
    return out


CA, CN = csvs(ARCH), csvs(ROOT)
IDS = ('degree', 'df', 'winsor_lo', 'winsor_hi', 'wval_lower', 'wval_upper')


def aligne(a, b):
    cles = [c for c in a.columns if c in b.columns and
            (a[c].dtype == object or c in IDS)]
    if not cles:
        return None, None
    a = a.copy(); b = b.copy()
    a['_k'] = a.groupby(cles, dropna=False).cumcount()
    b['_k'] = b.groupby(cles, dropna=False).cumcount()
    m = a.merge(b, on=cles + ['_k'], how='outer', suffixes=('_av', '_ap'), indicator=True)
    return m, cles


PAIRES_CSV = (('ci_lo', 'ci_hi'), ('ci95_lo', 'ci95_hi'), ('diff_lo', 'diff_hi'))
w()
w('   CIs OF THE CSV FILES (quantities emitted in tables):')
exc_csv_av = exc_csv_ap = 0
for r in sorted(set(CA) & set(CN)):
    try:
        a = pd.read_csv(CA[r]); b = pd.read_csv(CN[r])
    except Exception:
        continue
    for lo, hi in PAIRES_CSV:
        if not {lo, hi} <= set(a.columns) & set(b.columns):
            continue
        m, cles = aligne(a, b)
        if m is None:
            continue
        ea = (np.sign(m[lo + '_av']) == np.sign(m[hi + '_av'])) & (m[lo + '_av'] != 0)
        en = (np.sign(m[lo + '_ap']) == np.sign(m[hi + '_ap'])) & (m[lo + '_ap'] != 0)
        garde = ea & en & (np.sign(m[lo + '_av']) == np.sign(m[lo + '_ap']))
        exc_csv_av += int(ea.sum()); exc_csv_ap += int(en.sum())
        perdus = m[ea & ~garde]
        w('   %-44s %-16s excluding zero %3d -> %3d; sign lost %d' % (
            r, lo + '/' + hi, int(ea.sum()), int(en.sum()), len(perdus)))
        for _, row in perdus.iterrows():
            nom = ' '.join(str(row[c]) for c in cles)
            echecs.append('condition 2: %s %s' % (r, nom))
            w('      FAIL %s: [%.4g ; %.4g] -> [%.4g ; %.4g]' % (
                nom, row[lo + '_av'], row[hi + '_av'], row[lo + '_ap'], row[hi + '_ap']))
        orph = int((m['_merge'] != 'both').sum())
        if orph:
            w('      unmatched rows: %d' % orph)

# --------------------------------------------------- tests joints : verdicts
w()
w('   JOINT TESTS: verdict at 5 %% per row (columns p_* / pval_*)')
for r in sorted(set(CA) & set(CN)):
    try:
        a = pd.read_csv(CA[r]); b = pd.read_csv(CN[r])
    except Exception:
        continue
    pc = [c for c in a.columns if c in b.columns and (c.startswith('p_') or c.startswith('pval'))
          and pd.api.types.is_numeric_dtype(a[c])]
    if not pc:
        continue
    m, cles = aligne(a, b)
    if m is None:
        continue
    for c in pc:
        va = m[c + '_av'] < 0.05; vn = m[c + '_ap'] < 0.05
        ok = m[c + '_av'].notna() & m[c + '_ap'].notna()
        bascule = m[ok & (va != vn)]
        w('   %-44s %-12s rejections at 5 %%: %3d -> %3d out of %3d; flips %d' % (
            r, c, int(va[ok].sum()), int(vn[ok].sum()), int(ok.sum()), len(bascule)))
        for _, row in bascule.iterrows():
            nom = ' '.join(str(row[k]) for k in cles)
            echecs.append('condition 2 (joint test): %s %s %s' % (r, c, nom))
            w('      FAIL %s: p %.4g -> %.4g' % (nom, row[c + '_av'], row[c + '_ap']))

# ------------------------------------------------------------- condition 3
w()
w('=' * 78)
w('3. INTERVALS EXCLUDING ZERO: macros %d -> %d; CSV %d -> %d' % (
    exc_av, exc_ap, exc_csv_av, exc_csv_ap))
w('=' * 78)
tot_av, tot_ap = exc_av + exc_csv_av, exc_ap + exc_csv_ap
if tot_ap < tot_av / 2.0:
    echecs.append('condition 3: %d -> %d (below half)' % (tot_av, tot_ap))
w('   total %d -> %d: %s' % (tot_av, tot_ap, 'PASS' if tot_ap >= tot_av / 2.0 else 'FAIL'))

# ------------------------------------------------------------- condition 4
w()
w('=' * 78)
w('4. SAMPLE SIZES AND COUNTRIES')
w('=' * 78)
for k in sorted(k for k in A if 'Nobs' in k or k.endswith('Ncells') or 'Ncountr' in k):
    xa, xn = num(A[k]), num(N.get(k, ''))
    rap = xn / xa if xa and xn is not None else None
    w('   %-34s %-14s -> %-14s %s' % (k, A[k], N.get(k, '?'),
                                    'ratio %.4f' % rap if rap else ''))
    if rap is None or not 0.5 <= rap <= 2.0:
        echecs.append('condition 4: %s' % k)
for r in sorted(set(CA) & set(CN)):
    try:
        a = pd.read_csv(CA[r]); b = pd.read_csv(CN[r])
    except Exception:
        continue
    for c in [c for c in a.columns if c in b.columns and (c == 'nobs' or c.startswith('n_obs'))]:
        ua, ub = sorted(a[c].dropna().unique()), sorted(b[c].dropna().unique())
        w('   %-44s %-14s %s -> %s' % (r, c, ua[:4], ub[:4]))
E0 = pd.read_csv(os.path.join(ARCH, 'intermediate_outputs', 'ERA_tas_annual.csv'), keep_default_na=False)
E1 = pd.read_csv(os.path.join(ROOT, 'intermediate_outputs', 'ERA_tas_annual.csv'), keep_default_na=False)
p0, p1 = set(E0.GID_0), set(E1.GID_0)
w('   countries (ERA_tas_annual.csv): %d -> %d' % (len(p0), len(p1)))
w('   entered: %s' % (' '.join(sorted(p1 - p0)) or 'none'))
w('   left   : %s' % (' '.join(sorted(p0 - p1)) or 'none'))
if not 0.8 <= len(p1) / max(len(p0), 1) <= 1.25:
    echecs.append('condition 4: number of countries %d -> %d' % (len(p0), len(p1)))

# ------------------------------------------------ jours qui changent de bin
w()
w('=' * 78)
w('5. DAYS THAT CHANGE BIN (ERA_tas_annual.csv, common countries)')
w('=' * 78)
m = E0.merge(E1, on=['GID_0', 'year', 'value'], how='outer', suffixes=('_av', '_ap'))
m = m[m.GID_0.isin(p0 & p1)]
m[['count_av', 'count_ap']] = m[['count_av', 'count_ap']].fillna(0)
g = m.groupby(['GID_0', 'year']).apply(
    lambda d: pd.Series({'tv': 0.5 * (d.count_ap - d.count_av).abs().sum(),
                         'n': d.count_av.sum()}), include_groups=False)
g['f'] = g.tv / g.n
w('   overall (weighted days): %.4f' % (g.tv.sum() / g.n.sum()))
w('   by country-year: median %.4f, p90 %.4f, max %.4f' % (
    g.f.median(), g.f.quantile(.9), g.f.max()))
pc = g.groupby(level=0).apply(lambda d: d.tv.sum() / d.n.sum()).sort_values(ascending=False)
w('   by country, the first 15:')
for iso, v in pc.head(15).items():
    w('      %s %.4f' % (iso, v))
w('   countries with fraction 0: %d out of %d' % (int((pc <= 1e-12).sum()), len(pc)))

# ------------------------------------------- delta* d'Oster (lot 57, author's decision)
w()
w('=' * 78)
w('6. OSTER DELTA* BY SPECIFICATION (oster_delta_full.csv, delta_term1), before / after')
w('   (needed to rewrite the sentence of §2.3; reference = 00_parameters.R)')
w('=' * 78)
try:
    oa = pd.read_csv(os.path.join(ARCH, 'output', 'oster_delta_full.csv')).set_index('tag')
    on = pd.read_csv(os.path.join(ROOT, 'output', 'oster_delta_full.csv')).set_index('tag')
    tags = [t for t in oa.index if t in on.index] + [t for t in on.index if t not in oa.index]
    rang = lambda d: d[d.delta_term1 > 0].delta_term1.rank(method='min').astype('Int64')
    ra, rn = rang(oa), rang(on)
    w('%-26s %12s %12s  %9s %9s' % ('specification', 'before', 'after', 'rank+ bef', 'rank+ aft'))
    for t in tags:
        a = oa.delta_term1.get(t, float('nan')); n = on.delta_term1.get(t, float('nan'))
        w('%-26s %12.4f %12.4f  %9s %9s' % (t, a, n, ra.get(t, '-'), rn.get(t, '-')))
    w('   (rank+: rank among the STRICTLY positive delta*, 1 = the smallest)')
except Exception as ex:   # un morceau du releve qui echoue s'ecrit, il ne tue pas le releve
    w('   INCONCLUSIVE: %s' % ex)

# ------------------------------------------------------------- controles
w()
w('=' * 78)
w('CONDITION 1: nine checks')
w('   NB (author\'s decision): in the MIDDLE of the chain, known_bad_emission and check_fraicheur')
w('   fail BY CONSTRUCTION -- values.tex and the GE artefacts are refreshed only')
w('   in 02d. These failures are EXPECTED and do not count as defects.')
w('=' * 78)
if CTRL and os.path.exists(CTRL):
    for l in open(CTRL, encoding='utf-8'):
        w('   ' + l.rstrip())
else:
    w('   (not attached)')

w()
w('=' * 78)
w('VERDICT conditions 2-4: %s' % ('PASS' if not echecs else 'FAIL'))
for e in echecs:
    w('   ' + e)
open(SORTIE, 'w', encoding='utf-8').write('\n'.join(L) + '\n')
print('\n'.join(L[-(len(echecs) + 3):]))
sys.exit(0 if not echecs else 1)
