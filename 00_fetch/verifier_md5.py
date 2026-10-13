#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""verifier_md5.py -- verifie data_inputs/ contre MD5_data_inputs.tsv (lot 34, 25/09/2026).

NE TELECHARGE RIEN ET NE MODIFIE RIEN. C'est la moitie « verification a
l'arrivee » des scripts de recuperation : le reproducteur depose les sources
dans data_inputs/, puis lance ce script, qui dit fichier par fichier si ce qu'il
a obtenu est octet pour octet ce qui a produit les resultats.

MD5_data_inputs.tsv porte le md5 COMPLET de chaque fichier (chemin relatif a
data_inputs/). data_inputs/PROVENANCE.md n'en porte que les douze premiers
caracteres, et son md5 d'ensemble hache ces prefixes : une empreinte affaiblie
et une recette qui a deja fait conclure a tort a une discordance (registre,
lot 33).

Usage :  python 00_fetch/verifier_md5.py [source ...]
         sans argument : tout ; sinon les sous-repertoires nommes (GLORIA, BACI_CEPII_trade...)
Code  :  0 si tout concorde, 1 si un fichier manque ou differe.
"""
import hashlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DI = os.path.join(ROOT, 'data_inputs')
MAN = os.path.join(HERE, 'MD5_data_inputs.tsv')


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def main(sources):
    attendu = []
    for l in open(MAN, encoding='utf-8'):
        h, rel = l.rstrip('\n').split('\t')
        if not sources or rel.split('/')[0] in sources:
            attendu.append((rel, h))
    if not attendu:
        sys.exit('[STOP] aucune entree du manifeste pour %s' % ', '.join(sources))
    ok = manq = diff = 0
    for rel, h in attendu:
        p = os.path.join(DI, rel)
        if not os.path.exists(p):
            print('  ABSENT   %s' % rel)
            manq += 1
        elif md5(p) != h:
            print('  DIFFERE  %s' % rel)
            diff += 1
        else:
            ok += 1
    print('%d concordent, %d absents, %d differents (sur %d)' % (ok, manq, diff, len(attendu)))
    return 0 if manq == diff == 0 else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
