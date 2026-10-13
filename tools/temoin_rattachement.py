"""temoin_rattachement.py -- la regle de rattachement des pixels hors polygone (lot 47).

LA REGLE. Un pixel GHSL dont le centre tombe hors de tout polygone GADM est rattache au
polygone le plus proche SI ET SEULEMENT SI la distance de son centre a ce polygone est
inferieure au PAS DE LA TRAME. Elle dit quelque chose de falsifiable : « ceci est un
artefact de discretisation, a la precision de la grille pres ». Un pixel HABITE dont le
plus proche polygone est au-dela du pas n'est pas un artefact : c'est une anomalie, le
temoin REFUSE et la nomme au lieu de l'absorber.

Le temoin compte les pixels rattaches et leur population (absolu, part du monde),
imprime le seuil et la liste des rattaches, et refuse toute anomalie habitee.
"""
import numpy as np


def temoin_rattachement(x, y, pop, dist, pas, monde, sortie=None, cles=None, connues=None):
    """x, y : centres des pixels hors polygone et habites ; pop : population ;
    dist : distance (degres) au polygone le plus proche ; pas : pas de la trame.
    cles, connues (lot 49, decision (a)) : la position (ligne, colonne) de chaque pixel, et
    l'ensemble des positions des ANOMALIES CONNUES arretees par l'auteur (author's decision). Un pixel au-dela
    du pas qui figure dans la liste est rapporte, pas refuse ; un pixel au-dela du pas qui
    n'y figure PAS est refuse : la liste ne peut que rester stable ou raccourcir."""
    x, y, pop, dist = map(np.asarray, (x, y, pop, dist))
    if connues is not None:
        loin0 = dist >= pas
        dans = np.array([c in connues for c in cles]) & loin0
        print('  known anomalies (list fixed by the author): %d in the list, %d present '
              'in this pass, population %.4g (%.3f %% of the world)'
              % (len(connues), dans.sum(), pop[dans].sum(), 100 * pop[dans].sum() / monde))
        x, y, pop, dist = x[~dans], y[~dans], pop[~dans], dist[~dans]
    ok = dist < pas
    print('[temoin rattachement] threshold = grid step = %.8f deg (~%.0f m at the equator)'
          % (pas, pas * 111320))
    print('  inhabited pixels outside any polygon: %d, population %.4g (%.3f %% of the world)'
          % (len(x), pop.sum(), 100 * pop.sum() / monde))
    print('  ATTACHED (distance < step): %d, population %.4g (%.3f %% of the world)'
          % (ok.sum(), pop[ok].sum(), 100 * pop[ok].sum() / monde))
    if sortie is not None:
        np.savetxt(sortie, np.column_stack([x[ok], y[ok], pop[ok], dist[ok]]), delimiter=',',
                   header='x,y,pop,dist_deg', comments='', fmt='%.8f')
        print('  list of attached pixels: %s' % sortie)
    loin = ~ok
    if loin.any():
        k = np.argsort(-pop[loin])[:10]
        noms = ', '.join('(%.4f ; %.4f) pop %.1f at %.4f deg' % (x[loin][i], y[loin][i], pop[loin][i],
                                                                 dist[loin][i]) for i in k)
        raise ValueError('[STOP] temoin rattachement: %d inhabited pixel(s) outside any polygon beyond '
                         'the step (%.4g people) -- an anomaly, not an artefact: %s'
                         % (loin.sum(), pop[loin].sum(), noms))
    return int(ok.sum()), float(pop[ok].sum())
