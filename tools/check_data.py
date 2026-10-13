#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""check_data.py -- presence et md5 de chaque intrant de DATA.md (lot 84, 03/10/2026).

NE TELECHARGE RIEN ET NE MODIFIE RIEN. Lance en premier par run_all.sh.

Deux manifestes, une seule verite chacun :
  00_fetch/MD5_data_inputs.tsv       donnees publiques, chemins relatifs a data_inputs/
  00_fetch/MD5_intrants_derives.tsv  sorties des etapes SERVEUR et intrants sans producteur,
                                     chemins relatifs a la racine du paquet

Controle de coherence : chaque ligne de la table de DATA.md porte un md5. Pour une
ligne qui designe un fichier, c'est le md5 du fichier ; pour une ligne qui designe
un repertoire, c'est le md5 de l'ensemble, calcule sur les md5 COMPLETS du manifeste
(lignes « md5  chemin » triees, une par fichier, terminees par \n). Une table qui
diverge du manifeste arrete ici : la documentation ne peut pas promettre un autre
fichier que celui que la chaine verifie.

Usage :  python tools/check_data.py              presence + md5 (lit environ 30 Go)
         python tools/check_data.py --presence   presence seule (rapide)
Code  :  0 si tout concorde ; 1 si un intrant manque ou differe (liste rendue) ;
         2 si DATA.md et les manifestes divergent.
"""
import hashlib
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
MAN_PUB = os.path.join(ROOT, '00_fetch', 'MD5_data_inputs.tsv')
MAN_DER = os.path.join(ROOT, '00_fetch', 'MD5_intrants_derives.tsv')
DATA_MD = os.path.join(ROOT, 'DATA.md')


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def lire_manifeste(path, prefixe):
    out = []
    for l in open(path, encoding='utf-8'):
        l = l.rstrip('\n')
        if not l.strip() or l.startswith('#'):
            continue
        h, rel = l.split('\t')[:2]
        out.append((h.strip(), (prefixe + rel.strip()).replace('\\', '/')))
    return out


def md5_ensemble(entrees):
    s = ''.join('%s  %s\n' % (h, p) for h, p in sorted(entrees, key=lambda e: e[1]))
    return hashlib.md5(s.encode('utf-8')).hexdigest()


def coherence(attendus):
    """Chaque ligne de la table de DATA.md : `chemin` | ... | `md5`."""
    if not os.path.exists(DATA_MD):
        print('[STOP] DATA.md missing: %s' % DATA_MD)
        return 1
    # Le chemin retenu est le DERNIER chemin entre accents graves avant le md5 : une autre
    # colonne de la ligne peut citer un fichier (la licence renvoie a PROVENANCE.md).
    rx_md5 = re.compile(r'`([0-9a-f]{32})`')
    rx_chemin = re.compile(r'`((?:data_inputs|intermediate_outputs|output)/[^`]*)`')
    n = defauts = 0
    for l in open(DATA_MD, encoding='utf-8'):
        m = rx_md5.search(l) if l.startswith('|') else None
        chemins = rx_chemin.findall(l[:m.start()]) if m else []
        if not chemins:
            continue
        chemin, h = chemins[-1].rstrip('/'), m.group(1)
        n += 1
        sous = [(a, p) for a, p in attendus if p == chemin or p.startswith(chemin + '/')]
        if not sous:
            print('  DIVERGE  DATA.md cites %s, which is in no manifest' % chemin)
            defauts += 1
            continue
        attendu = sous[0][0] if (len(sous) == 1 and sous[0][1] == chemin) else md5_ensemble(sous)
        if attendu != h:
            print('  DIVERGE  %s: DATA.md %s, manifest %s' % (chemin, h, attendu))
            defauts += 1
    print('[coherence] %d line(s) of DATA.md, %d divergence(s)' % (n, defauts))
    return defauts


def main(argv):
    presence = '--presence' in argv
    attendus = lire_manifeste(MAN_PUB, 'data_inputs/') + lire_manifeste(MAN_DER, '')
    if coherence(attendus):
        print('[STOP] DATA.md and the manifests diverge. Code 2.')
        return 2
    manquants, differents, ok = [], [], 0
    for h, rel in attendus:
        p = os.path.join(ROOT, rel)
        if not os.path.isfile(p):
            manquants.append(rel)
            continue
        if not presence and md5(p) != h:
            differents.append(rel)
            continue
        ok += 1
    for rel in manquants:
        print('  MISSING  %s' % rel)
    for rel in differents:
        print('  DIFFERS  %s' % rel)
    mode = 'presence' if presence else 'presence + md5'
    print('[check_data] %s: %d match, %d missing, %d differ (out of %d)'
          % (mode, ok, len(manquants), len(differents), len(attendus)))
    if manquants or differents:
        print('  -> see DATA.md, section "Obtaining the data".')
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
