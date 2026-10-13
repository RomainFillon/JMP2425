"""garde_oster.py -- l'emission Oster de §2.3 et sa garde D3 (lot 63, 29/09/2026).

CE QUE LA PHRASE AFFIRME. §2.3 : « the smallest positive selection ratio over the
\\EMPgridNspec specifications ..., \\delta^* = \\EMPosterDelta (\\EMPosterDeltaSpec) ...
against a threshold of one. Among the specifications, \\EMPosterDeltaNegN return a
negative \\delta^* ...; the ratio is read off the others. » Supplement l.1167-1169 :
« the informative entry is the least favourable entry on the grid, \\delta^*=\\EMPosterDelta,
whose coefficient moves by \\EMPosterMoveRel\\% ». Aucune des deux ne dit QUELLE
specification gagne : elles disent que la moins favorable depasse un.

AUTHOR'S DECISION (lot 63). Les quatre macros lisent la cellule la MOINS
FAVORABLE (plus petit delta* strictement positif), et non plus la reference de
00_parameters.R. L'ancienne D3 surveillait une coincidence de grille (la reference
porte-t-elle le minimum ?), qui s'est defaite le 29/09. La nouvelle garde surveille :
  1. \\EMPosterDelta     = le minimum des delta* positifs de la table ;
  2. \\EMPosterDeltaSpec = la specification qui l'atteint, en mots ;
  3. \\EMPosterDeltaNegN = le compte des delta* strictement negatifs ;
  4. ce minimum depasse UN -- l'affirmation sur le monde, la seule qui puisse
     refuser sur des donnees correctes ;
  5. \\EMPosterMoveRel   decrit la meme specification que \\EMPosterDelta.
1-3 et 5 sont des coherences entre macros et table ; 4 est une phrase du papier.
Les verifications RECALCULENT depuis la table, independamment de emissions() : une
emission fausse (test connu-faux) est attrapee.
"""
import pandas as pd

DEG_MOT = {2: 'quadratic', 3: 'cubic', 4: 'quartic'}


def _pct(x):
    """1.0 -> '1', 0.5 -> '0.5'.  Pas de zero terminal dans une prose."""
    return '%g' % float(x)


def spec_en_mots(row):
    """Nom en clair d'une ligne de oster_delta_full.csv. Seuils RELATIFS (%), pas
    des degres : voir la correction du 22/09 dans emit_values_empirics.py."""
    deg, fam = int(row['degree']), str(row['spec'])
    if fam == 'polynomial':
        if deg not in DEG_MOT:
            raise SystemExit("[STOP] degree %d has no name in DEG_MOT; add the word "
                             "rather than print a bare digit." % deg)
        forme = 'the %s polynomial' % DEG_MOT[deg]
    elif fam == 'spline':
        forme = 'the natural spline with %d degrees of freedom' % deg
    else:
        raise SystemExit("[STOP] unknown specification family: %r" % fam)
    wl, wu = _pct(row['wval_lower']), _pct(row['wval_upper'])
    wins = ('winsorised at $%s\\%%$' % wl) if wl == wu else \
           ('winsorised at $%s\\%%$ and $%s\\%%$' % (wl, wu))
    return '%s %s' % (forme, wins)


def move_rel(row):
    return 100 * abs((row['b_dot'] - row['b_tilde']) / row['b_tilde'])


def cellule_la_moins_favorable(od):
    pos = od[od['delta_term1'] > 0]
    if len(pos) == 0:
        raise SystemExit("[STOP] oster_delta_full.csv: no positive delta*.")
    return pos['delta_term1'].idxmin()


def emissions(od):
    """Les quatre macros, depuis la cellule la moins favorable. `od` indexe par tag."""
    tag = cellule_la_moins_favorable(od)
    row = od.loc[tag]
    return {'tag': tag, 'tag_moverel': tag,
            'EMPosterDelta': '%.1f' % float(row['delta_term1']),
            'EMPosterDeltaSpec': spec_en_mots(row),
            'EMPosterDeltaNegN': str(int((od['delta_term1'] < 0).sum())),
            'EMPosterMoveRel': '%.1f' % move_rel(row)}


def verifie(od, emis):
    """Rend (coherence, monde) : la liste des coherences violees (points 1, 2, 3, 5 ;
    bloquantes) et None ou le message du point 4 (phrase du papier)."""
    v = od['delta_term1'].astype(float)
    pos = v[v > 0]
    err = []
    if len(pos) == 0:
        return ['no positive delta* in the table'], None
    vmin, tmin = float(pos.min()), pos.idxmin()
    if emis['EMPosterDelta'] != '%.1f' % vmin:
        err.append("1. \\EMPosterDelta = %s, minimum of the positives = %.4f (%s)"
                   % (emis['EMPosterDelta'], vmin, tmin))
    if emis['tag'] != tmin or emis['EMPosterDeltaSpec'] != spec_en_mots(od.loc[tmin]):
        err.append("2. \\EMPosterDeltaSpec names %s (%r), the minimum is reached by %s"
                   % (emis['tag'], emis['EMPosterDeltaSpec'], tmin))
    n_neg = int((v < 0).sum())
    if emis['EMPosterDeltaNegN'] != str(n_neg):
        err.append("3. \\EMPosterDeltaNegN = %s, the table carries %d negative delta*"
                   % (emis['EMPosterDeltaNegN'], n_neg))
    if emis['tag_moverel'] != emis['tag'] or \
            emis['EMPosterMoveRel'] != '%.1f' % move_rel(od.loc[emis['tag']]):
        err.append("5. \\EMPosterMoveRel describes %s (%s), \\EMPosterDelta describes %s"
                   % (emis['tag_moverel'], emis['EMPosterMoveRel'], emis['tag']))
    monde = None if vmin > 1 else (
        "4. the smallest positive delta* is %.4f (%s): it does not exceed the threshold of one"
        % (vmin, tmin))
    return err, monde
