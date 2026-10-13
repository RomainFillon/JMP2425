#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""known_bad_emission.py -- le test qui a trouve le defaut, devenu permanent.

CE QU'IL VERIFIE
----------------
Un refus en cours d'emission doit laisser l'artefact PRECEDENT intact.

LE DEFAUT QU'IL A TROUVE (11/09/2026)
-------------------------------------
    with open(VALUES_TEX, 'w') as fh:
        fh.write(emit_values_tex())

Python evalue `open(..., 'w')` AVANT l'appel : le fichier est tronque a zero
avant que le contenu n'existe. Toute garde qui s'arrete dans emit_values_tex()
-- le refus md5 de decompose_metric_channels, par exemple -- detruisait donc
values.tex. Mesure : exit 2, values.tex a 0 octet, 104 macros perdues. Un
latexmk qui suit rend un PDF ou chaque valeur GE est indefinie.

Correctif en place : contenu evalue d'abord, ecriture dans un .tmp, os.replace.

POURQUOI IL TOURNE A CHAQUE EMISSION
------------------------------------
Un test qu'on ne lance que sur demande ne protege que le jour ou on y pense.
Il garde une propriete dont la violation est SILENCIEUSE -- l'emetteur sort en
erreur, ce qui a l'air correct, pendant que l'artefact est detruit.

LE BAC A SABLE (lot 85o, 07/10/2026)
------------------------------------
Jusqu'au lot 85n, ce controle perturbait output/welfare_sc2_q50.csv EN PLACE et
le restaurait dans un `finally`. Deux fois (lots 84 et 85l) la source est restee
perturbee (une cellule x 0,98) : un `finally` ne s'execute pas quand le
processus est TUE, et un 02d interrompu laissait l'arbre vivant faux.
Desormais le controle n'ecrit JAMAIS dans output/ :
  1. output/ est COPIE dans un repertoire temporaire (<tmp>/output) ;
  2. l'emetteur y tourne avec EMIT_SANDBOX_DIR=<tmp>/output et
     RUNS_DIR=<tmp>/output : il y lit et y ecrit (emit_values.py, reprise.py) ;
  3. TEMOIN POSITIF : sur la copie intacte l'emetteur doit rendre 0 et
     reproduire values.tex et tab_ge_north.tex A L'OCTET -- sans lui, un bac a
     sable casse (fichier manquant, chemin faux) ferait refuser l'emetteur et
     le test passerait pour une mauvaise raison ;
  4. perturbation de la COPIE de welfare_sc2_q50.csv ; l'emetteur doit sortir
     NON NUL et les artefacts de la copie doivent etre INCHANGES au bit ;
  5. output/ : md5 ET horodatage de chaque fichier identiques avant/apres.
     Un ecart fait echouer le controle.
Le repertoire temporaire est retire a la fin (il ne contient que des copies).
Un processus tue laisse au pire une copie orpheline dans le repertoire
temporaire du systeme, jamais un output/ perturbe.

Usage:  python tools/known_bad_emission.py
Sortie: rc=0 si la propriete tient, rc=1 sinon.
"""
import hashlib
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'output')
EMITTER = os.path.join(ROOT, '02_GE', 'emit_values.py')
SOURCE_NOM = 'welfare_sc2_q50.csv'
WATCHED_NOMS = ['values.tex', 'tab_ge_north.tex']


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as fh:
        for c in iter(lambda: fh.read(1 << 20), b''):
            h.update(c)
    return h.hexdigest()


def etat(racine):
    """{chemin relatif: (mtime_ns, md5)} de tous les fichiers sous racine."""
    e = {}
    for d, _, fs in os.walk(racine):
        for f in fs:
            q = os.path.join(d, f)
            e[os.path.relpath(q, racine)] = (os.stat(q).st_mtime_ns, md5(q))
    return e


def run_emitter(bac):
    env = dict(os.environ, RUNS_DIR=bac, EMIT_SANDBOX_DIR=bac)
    p = subprocess.run([sys.executable, '-B', 'emit_values.py'],
                       cwd=os.path.join(ROOT, '02_GE'), env=env,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    return p.returncode, p.stdout.decode('utf-8', 'replace')


def derniere_ligne_utile(txt):
    l = [x for x in txt.splitlines() if x.strip()]
    return l[-1][:200] if l else '(empty output)'


def main():
    for p in [os.path.join(OUT, SOURCE_NOM), EMITTER] + \
             [os.path.join(OUT, n) for n in WATCHED_NOMS]:
        if not os.path.exists(p):
            print('[SKIP] missing: %s' % os.path.relpath(p, ROOT))
            return 0

    avant = etat(OUT)
    failures = []
    tmp = tempfile.mkdtemp(prefix='known_bad_')
    bac = os.path.join(tmp, 'output')
    try:
        shutil.copytree(OUT, bac)
        print('[known-bad] sandbox: copy of output/ (%d files) in %s'
              % (len(avant), bac))

        # --- temoin positif : copie intacte -> rc 0 et memes artefacts
        rc, txt = run_emitter(bac)
        print('[known-bad] emitter on the intact copy: code %d' % rc)
        if rc != 0:
            failures.append('positive control: the emitter fails on the INTACT copy '
                            '(code %d: %s). The sandbox does not reproduce output/.'
                            % (rc, derniere_ligne_utile(txt)))
        for n in WATCHED_NOMS:
            if md5(os.path.join(bac, n)) != avant[n][1]:
                failures.append('positive control: %s of the copy differs from output/%s.'
                                % (n, n))
        if failures:
            raise RuntimeError('positive control')

        # --- perturbation de la COPIE
        import pandas as pd
        src = os.path.join(bac, SOURCE_NOM)
        d = pd.read_csv(src)
        cols = [c for c in d.columns if not c.startswith('Unnamed')]
        d.loc[d.index[-1], cols[0]] *= 0.98
        d.to_csv(src, index=False)
        print('[known-bad] copy of %s perturbed (md5 %s -> %s)'
              % (SOURCE_NOM, avant[SOURCE_NOM][1][:12], md5(src)[:12]))
        temoins = {n: md5(os.path.join(bac, n)) for n in WATCHED_NOMS}

        rc, txt = run_emitter(bac)
        print('[known-bad] emitter under perturbation: code %d (%s)'
              % (rc, derniere_ligne_utile(txt)))
        if rc == 0:
            failures.append(
                "the emitter returned 0 despite a perturbed source: no check "
                "bites.")
        for n in WATCHED_NOMS:
            if md5(os.path.join(bac, n)) != temoins[n]:
                failures.append(
                    '%s WAS MODIFIED by an emission that failed (%d bytes). '
                    "This is the truncation defect: a refusal must leave "
                    "the previous artefact intact."
                    % (n, os.path.getsize(os.path.join(bac, n))))
            else:
                print('[known-bad] %-28s INTACT' % n)
    except RuntimeError:
        pass
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    # --- l'arbre vivant : rien n'a bouge, ni contenu ni horodatage
    apres = etat(OUT)
    bouges = sorted(k for k in set(avant) | set(apres) if avant.get(k) != apres.get(k))
    if bouges:
        failures.append('output/ CHANGED during the check (%d file(s)): %s'
                        % (len(bouges), ', '.join(bouges[:10])))
    else:
        print('[known-bad] output/: %d files, md5 and timestamps identical '
              'before/after' % len(apres))

    if failures:
        print('\n[known-bad] FAILURE:')
        for f in failures:
            print('   - ' + f)
        return 1
    print('[known-bad] OK: a refusal leaves the artefacts intact, output/ untouched.')
    return 0


if __name__ == '__main__':
    sys.exit(main())
