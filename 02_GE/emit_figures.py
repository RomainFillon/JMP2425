#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""
emit_figures.py -- the two general-equilibrium figures of the paper.

Writes, into output/ (where \graphicspath points):
  - ge_welfare_trajectories.pdf   (fig:ge_traj)   three margins by horizon,
                                  population-weighted fan
  - map_welfare_2100_country.pdf  (fig:ge_maps)   country-aggregated map,
                                  three margins
  - fig_ge_tiers.pdf              (fig:ge_tiers)  the three margins by quartile
                                  of 1990 GDP per capita, with the arm B band
                                  on the daily total

fig_ge_tiers.pdf RECALCULE RIEN. Il trace output/fig_ge_tiers.csv, ecrit par
emit_values.py dans le meme passage que les macros \GEtier*. C'est ce qui rend
impossible qu'une barre et la phrase qui la commente divergent : il n'y a pas
deux chemins de calcul a tenir d'accord, il n'y en a qu'un. Ce script tourne
APRES emit_values.py, et 02d_emit_latex.R les lance dans cet ordre.

Reads:  output/welfare_sc{1,2,3}_q50.csv (the solver grid),
        intermediate_outputs/{indices_drop,identifier_xy_final,population0,
        poids_par_cellule}.csv
Needs:  the solver grid must have run, and 02a must have produced the inputs.

Shares its conventions with emit_values.py: same margins, same horizon
convention (welfare is (T, N), the last transition is labelled 2100), and the
same RUNS_DIR and CAP_TAG overrides.

This script reads the welfare files only, never the hat-paths, so it is
unaffected by their orientation.

The map aggregates to the country before plotting rather than colouring each
pixel. A per-pixel map is not invariant to the column-mass cap on the trade
matrix, whereas the country aggregate is, so a pixel map would show a pattern
that the sensitivity analysis does not support. Each pixel therefore takes the
aggregate value of its dominant country, which also removes any dependence on
an external shapefile.
"""
import os
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

# LOT 74-D (01/10/2026) : PDF SANS DATE. Le backend pdf ecrivait /CreationDate a
# chaque passage : contenu identique, octets differents, et reference_phase.py (qui
# epingle les octets des artefacts) en echec apres TOUT rejeu de 02d. None retire le
# champ (matplotlib >= 3.x) ; SOURCE_DATE_EPOCH n'est alors pas lu. Les autres champs
# d'en-tete -- Creator « Matplotlib vX », Producer « Matplotlib pdf backend vX » --
# ne dependent que de la version, et le fichier ne porte pas de /ID.
PDF_META = {'CreationDate': None, 'ModDate': None}

HERE = os.path.dirname(os.path.abspath(__file__))
PKG  = os.path.dirname(HERE)

## The package directory is the root; there is nothing above it. Pointing REPO
## one level higher would send the figures outside the tree the manuscript reads,
## and nothing would report it: the file would simply never be seen.
REPO = PKG

## Horizon labels. Single source, shared with emit_values.py and 00_parameters.R.
## The terminal tick sits at the END of the terminal block, which is where the
## label belongs; the plotted data still stops at the block's own index year.
sys.path.insert(0, os.path.join(PKG, '01_estimation'))
from importlib import import_module as _import_module
_PARAMS = _import_module('00_parameters')
HORIZON_LABEL = _PARAMS.HORIZON_LABEL

## The solver inputs live in intermediate_outputs/.
INPUTS = os.path.join(PKG, 'intermediate_outputs')
OUTPUT = os.path.join(PKG, 'output')

## RUNS_DIR arrives in three forms: an absolute path, "output" (the directory
## itself, which is what the orchestrator exports), or the name of a subdirectory
## of output/. Joining it naively to PKG/output would give output/output.
_rd = os.environ.get('RUNS_DIR', '')
if os.path.isabs(_rd):
    RUNS_DIR = _rd
elif _rd in ('', 'output'):
    RUNS_DIR = os.path.join(PKG, 'output')
else:
    RUNS_DIR = os.path.join(PKG, 'output', _rd)
RUNS_DIR = os.path.normpath(RUNS_DIR)

## Quantile of the arm plotted. q50 is the headline; the other arms carry their
## quantile in the file name.
QUANTILE = os.environ.get('QUANTILE', 'q50')

for _n, _p in (('RUNS_DIR', RUNS_DIR), ('INPUTS', INPUTS)):
    if not os.path.isdir(_p):
        raise SystemExit(
            f"[STOP] {_n} missing: {_p}\n"
            f"  RUNS_DIR is produced by the solver's run grid;\n"
            f"  INPUTS by 02a_calibration_at_1deg.r.")
os.makedirs(OUTPUT, exist_ok=True)
CAP_TAG = os.environ.get('CAP_TAG', '')

## L'ETIQUETTE D'AMENITE. Le corps du papier tourne a la borne basse f_a = 0,
## donc les figures qu'il montre en sortent aussi. Le suffixe est celui que
## main.py compose ; sa presence est VERIFIEE sur les trois scenarios, parce
## qu'un bras a moitie present donnerait trois panneaux de deux modeles.
## AMEN_TAG='' force explicitement le bras interieur (sensibilite de S5.5).
_amen_env = os.environ.get('AMEN_TAG')
_NULL = '_NULLamen'
if _amen_env is None:
    AMEN_TAG = _NULL if all(
        os.path.exists(os.path.join(
            RUNS_DIR, f'welfare_sc{_s}_{QUANTILE}{CAP_TAG}{_NULL}.csv'))
        for _s in (1, 2, 3)) else ''
else:
    AMEN_TAG = _amen_env
print("[bras] GE figures: amenity %s (suffix %r)"
      % ('f_a = 0 (lower bound, BODY ARM)' if AMEN_TAG == _NULL
         else 'interior f_a^q50', AMEN_TAG))

YEAR0, DT = 2015, 5
NORTH_LAT = 45.0
MARGINS = {'daily': (2, 1), 'annual': (3, 1), 'extra': (2, 3)}
MARGIN_TITLE = {
    'daily':  'Daily total',
    'annual': 'Annual total',
    'extra':  'Daily extra',
}


def load_welfare(sc):
    df = pd.read_csv(os.path.join(
        RUNS_DIR, f'welfare_sc{sc}_{QUANTILE}{CAP_TAG}{AMEN_TAG}.csv'))
    return df.loc[:, ~df.columns.str.contains('^Unnamed')].values   # (T, N)


W = {sc: load_welfare(sc) for sc in (1, 2, 3)}
T, N = W[1].shape
years = np.array([YEAR0 + (t + 1) * DT for t in range(T)])
LABELS = {y: (HORIZON_LABEL if i == T - 1 else str(y))
          for i, y in enumerate(years)}

CUM = {m: np.cumsum(np.log(W[a] / W[b]), axis=0) for m, (a, b) in MARGINS.items()}

# ---- pixel metadata: id, lat/lon, population, dominant country ----
kept = pd.read_csv(os.path.join(INPUTS, 'indices_drop.csv')).iloc[:, 0].astype(int).values
assert len(kept) == N, f"indices_drop {len(kept)} != N {N}"
xy = pd.read_csv(os.path.join(INPUTS, 'identifier_xy_final.csv'))
xy.columns = [c.strip().strip('"') for c in xy.columns]
id_xy = dict(zip(xy['id'].astype(int), zip(xy['y'], xy['x'])))
pop = pd.read_csv(os.path.join(INPUTS, 'population0.csv'))
id_pop = dict(zip(pop['id'].astype(int), pop['value']))

pix = pd.DataFrame({'id': kept})
pix['y'] = [round(id_xy[i][0], 1) for i in kept]
pix['x'] = [round(id_xy[i][1], 1) for i in kept]
pix['pop'] = [id_pop[i] for i in kept]
# Lot 53 : meme garde que emit_values.py l.~540, meme traitement (rapport en reprise).
sys.path.insert(0, os.path.join(REPO, 'tools'))
from reprise import reprise_declaree, rapporte   # noqa: E402
if abs(pix['pop'].min() - 1.0052) >= 1e-3:
    if reprise_declaree():
        rapporte('emit_figures.py', 'min population0 over the retained cells',
                 '1.0052', '%.4f' % pix['pop'].min())
    else:
        raise AssertionError("unexpected population0 (see MANIFEST)")
w = pix['pop'].values

## Lot 51 (28/09/2026). Ce bloc prenait l'argmax de `share_population` de
## lengths_country.csv (1995), une part normalisee PAR PAYS : il designait le pays
## dont la cellule porte la plus grande part de SA propre population (RECAP 3.229),
## et son commentaire appelait cela « a correct use ». Desormais :
##   - la VALEUR d'un pays agrege ses cellules en proportion de leur population
##     (poids_par_cellule.csv, part_de_la_cellule, periode 2015) ;
##   - pour DESSINER, un pixel prend la couleur du pays qui porte la plus grande
##     part de SA population a lui (argmax de part_de_la_cellule, normalisee PAR
##     CELLULE) -- ce que l'ancien commentaire voulait dire.
sys.path.insert(0, os.path.join(REPO, 'tools'))
from temoins_poids import temoin_part_de_la_cellule   # noqa: E402
PARTS = pd.read_csv(os.path.join(INPUTS, 'poids_par_cellule.csv'), keep_default_na=False)
PARTS = PARTS[PARTS['GID_0'] != ''].copy()
temoin_part_de_la_cellule(PARTS)
PARTS['x'] = PARTS['x'].round(1)
PARTS['y'] = PARTS['y'].round(1)
PARTS = pix[['id', 'x', 'y']].merge(PARTS[['x', 'y', 'GID_0', 'part_de_la_cellule']],
                                    on=['x', 'y'])[['id', 'GID_0', 'part_de_la_cellule']]
dom = (PARTS.sort_values('part_de_la_cellule', ascending=False)
            .drop_duplicates('id')[['id', 'GID_0']])
pix = pix.merge(dom, on='id', how='left')        # pays AFFICHE du pixel (dessin seulement)


# ======================================================================
# Figure 1: welfare trajectories (fig:ge_traj)
# ======================================================================
def wq(vals, ww, q):
    o = np.argsort(vals)
    v, cw = vals[o], np.cumsum(ww[o]) / ww[o].sum()
    return np.interp(q, cw, v)


def emit_trajectories():
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), constrained_layout=True, sharey=True)
    for ax, m in zip(axes, ('daily', 'annual', 'extra')):
        lvl = CUM[m]   # (T, N)
        med = np.array([wq(lvl[t], w, .50) for t in range(T)])
        q25 = np.array([wq(lvl[t], w, .25) for t in range(T)])
        q75 = np.array([wq(lvl[t], w, .75) for t in range(T)])
        q05 = np.array([wq(lvl[t], w, .05) for t in range(T)])
        q95 = np.array([wq(lvl[t], w, .95) for t in range(T)])
        mean = (lvl * (w / w.sum())[None, :]).sum(axis=1)
        ax.fill_between(years, q05, q95, color='tab:blue', alpha=0.12, lw=0)
        ax.fill_between(years, q25, q75, color='tab:blue', alpha=0.25, lw=0)
        ax.plot(years, med, color='tab:blue', lw=2, label='median (pop-w)')
        ax.plot(years, mean, color='tab:red', lw=1.5, ls='--', label='mean (pop-w)')
        ax.axhline(0, color='k', lw=0.6, ls=':')
        ax.set_title(MARGIN_TITLE[m].replace('\\,', ' '), fontsize=10)
        ax.set_xlabel('Year')
        ax.set_xticks([2015, 2050, int(HORIZON_LABEL)])
    axes[0].set_ylabel(r'Cumulative $\Delta\log U$')
    axes[0].legend(fontsize=8, loc='lower left')
    out = os.path.join(OUTPUT, 'ge_welfare_trajectories.pdf')
    fig.savefig(out, bbox_inches='tight', metadata=PDF_META)
    plt.close(fig)
    print(f"[fig] {os.path.relpath(out, REPO)}")


# ======================================================================
# Figure 2: country-aggregated map (fig:ge_maps)
# ======================================================================
def emit_country_map():
    # Population-weighted aggregate by dominant country, per margin, at the
    # terminal horizon
    for m in MARGINS:
        pix[m] = CUM[m][-1, :]
    cty_val = {}
    pc = pix.drop(columns=['GID_0']).merge(PARTS, on='id', how='inner')
    pc['pop'] = pc['pop'] * pc['part_de_la_cellule']       # agregation PROPORTIONNELLE
    for m in MARGINS:
        g = pc.groupby('GID_0')
        cty_val[m] = {k: np.average(v[m], weights=v['pop']) if v['pop'].sum() > 0 else np.nan
                      for k, v in g}

    # One-degree raster canvas: each pixel takes the aggregate value of its own
    # country, so the map shows the country-level object, not a pixel one
    nlon, nlat = 360, 180
    def grid(m):
        Z = np.full((nlat, nlon), np.nan)
        for _, r in pix.dropna(subset=['GID_0']).iterrows():
            val = cty_val[m].get(r['GID_0'], np.nan)
            row = int(min(max(np.floor(90 - r['y']), 0), nlat - 1))
            col = int(min(max(np.floor(r['x'] + 180), 0), nlon - 1))
            Z[row, col] = val
        return Z

    vmax = max(abs(np.nanmin(list(cty_val['daily'].values()) + list(cty_val['annual'].values()))),
               abs(np.nanmax(list(cty_val['annual'].values()))))
    vmax = float(np.ceil(vmax * 20) / 20)   # round up to a clean tick
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    fig, axes = plt.subplots(3, 1, figsize=(9, 11), constrained_layout=True)
    for ax, m in zip(axes, ('daily', 'annual', 'extra')):
        im = ax.imshow(grid(m), cmap='RdBu', norm=norm,
                       extent=[-180, 180, -90, 90], aspect='auto', interpolation='nearest')
        ax.set_title(MARGIN_TITLE[m].replace('\\,', ' '), fontsize=10)
        ax.set_xticks([]); ax.set_yticks([])
        ax.axhline(NORTH_LAT, color='k', lw=0.5, ls='--', alpha=0.5)
    cb = fig.colorbar(
        im, ax=axes, shrink=0.6,
        label=r'Country-aggregate cumulative $\Delta\log U$ by ' + HORIZON_LABEL)
    out = os.path.join(OUTPUT, f'map_welfare_{HORIZON_LABEL}_country.pdf')
    fig.savefig(out, bbox_inches='tight', metadata=PDF_META)
    plt.close(fig)
    print(f"[fig] {os.path.relpath(out, REPO)}")


# ======================================================================
# Figure 3: welfare by 1990 GDP per capita quartile (fig:ge_tiers)
# ======================================================================
# CE QUE LA FIGURE MONTRE, ET CE QU'ELLE NE MONTRE PAS
# ----------------------------------------------------
# Trois barres par quartile -- daily total, annual total, daily extra -- a la
# borne basse d'amenite. Le SEGMENT ne porte que sur la barre daily total, et
# ce n'est pas une economie de trait : la bande vient de soixante solves du
# bras JOURNALIER, donc elle ne dit rien du bras annuel, et le daily extra,
# qui est une difference entre les deux bras, n'est pas une quantite tiree.
# Un segment sur les trois barres affirmerait trois incertitudes la ou une
# seule a ete mesuree.
#
# LES DEUX OBJETS SUPERPOSES NE SONT PAS LE MEME. La barre est le solve sur la
# courbe MEDIANE du Monte Carlo ; le segment est la distribution des SORTIES.
# Aucun des deux n'approxime l'autre, et emit_values.py asserte seulement que
# la barre tombe dans l'etendue des tirages -- ce qui est la condition pour
# que les deux puissent figurer sur le meme axe sans se contredire.
TIER_LABEL = {'Q1': 'Q1\n(poorest)', 'Q2': 'Q2', 'Q3': 'Q3',
              'Q4': 'Q4\n(richest)'}
TIER_ORDER = ('Q1', 'Q2', 'Q3', 'Q4')

## Trois tons qui se distinguent aussi en niveaux de gris : la revue imprime en
## noir et blanc. Le daily extra reprend le rouge des figures CRH du paquet.
BAR_STYLE = [
    ('daily',  'Daily total',  '#1f4e79'),
    ('annual', 'Annual total', '#a6bddb'),
    ('extra',  'Daily extra',  '#d73027'),
]


def emit_tiers():
    src = os.path.join(OUTPUT, 'fig_ge_tiers.csv')
    if not os.path.exists(src):
        raise SystemExit(
            "[STOP] %s missing: the quartile figure cannot be produced.\n"
            "  Its producer is emit_values.py, which writes it only if the\n"
            "  f_a = 0 arm AND the arm B grid are present. A gap here stops the\n"
            "  compilation on a missing includegraphics, which is the intended\n"
            "  polarity: a silent figure reads as a result."
            % os.path.relpath(src, REPO))
    d = pd.read_csv(src).set_index('tier').loc[list(TIER_ORDER)]

    ## LE CONTROLE QUE LA FIGURE PEUT FAIRE ELLE-MEME : une barre hors de son
    ## segment. emit_values.py l'asserte deja sur l'etendue des tirages ; ici on
    ## le refait sur ce qui est REELLEMENT TRACE, c'est-a-dire sur les colonnes
    ## du CSV, et non sur les objets Python en amont. Deux fichiers qui se
    ## desynchroniseraient passeraient le premier controle et pas celui-ci.
    for t in TIER_ORDER:
        r = d.loc[t]
        if not (r['band_lo'] <= r['band_hi']):
            raise SystemExit("[STOP] %s: band bounds inverted." % t)
        if not (min(r['band_lo'], r['daily']) <= r['daily']
                <= max(r['band_hi'], r['daily'])):
            raise SystemExit("[STOP] %s: bar outside the band." % t)

    plt.rcParams.update({'font.size': 9})
    CM = 1 / 2.54
    fig, ax = plt.subplots(figsize=(12 * CM, 8 * CM), constrained_layout=True)
    x = np.arange(len(TIER_ORDER), dtype=float)
    w = 0.26
    for k, (col, lab, colr) in enumerate(BAR_STYLE):
        ax.bar(x + (k - 1) * w, d[col].values, width=w, label=lab,
               color=colr, edgecolor='black', linewidth=0.4, zorder=2)
    ## Le segment est pose sur la barre daily total, donc au meme decalage.
    lo = d['daily'].values - d['band_lo'].values
    hi = d['band_hi'].values - d['daily'].values
    ax.errorbar(x - w, d['daily'].values, yerr=np.vstack([lo, hi]),
                fmt='none', ecolor='black', elinewidth=0.9, capsize=2.5,
                zorder=3)
    ax.axhline(0, color='black', lw=0.7, zorder=1)
    ax.set_xticks(x)
    ax.set_xticklabels([TIER_LABEL[t] for t in TIER_ORDER])
    ax.set_xlabel('Quartile of 1990 GDP per capita (population-weighted)')
    ax.set_ylabel(r'Cumulative $\Delta\log U$ by ' + HORIZON_LABEL)
    ## EN BAS A DROITE, et la place est MESUREE : le segment de Q1 descend
    ## a la borne basse de la bande, qui est la plus profonde des quatre,
    ## et il occupe le coin gauche. Une legende posee la recouvrait la
    ## moustache -- verifie a l'image, pas au code.
    ax.legend(frameon=False, fontsize=8, loc='lower right')
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    out = os.path.join(OUTPUT, 'fig_ge_tiers.pdf')
    fig.savefig(out, metadata=PDF_META)
    plt.close(fig)
    plt.rcParams.update({'font.size': plt.rcParamsDefault['font.size']})
    print("[fig] %s  (%d draws per quartile)"
          % (os.path.relpath(out, REPO), int(d['n_draws'].iloc[0])))
    for t in TIER_ORDER:
        r = d.loc[t]
        print("      %s  daily %+.4f [%+.4f, %+.4f]  annual %+.4f  "
              "extra %+.4f  neg %.0f%%"
              % (t, r['daily'], r['band_lo'], r['band_hi'], r['annual'],
                 r['extra'], r['share_neg']))


if __name__ == '__main__':
    os.makedirs(OUTPUT, exist_ok=True)
    emit_trajectories()
    emit_country_map()
    emit_tiers()
    print("[note] pixel-level map (not cap-invariant) NOT regenerated; "
          "it is not used by the paper.")
