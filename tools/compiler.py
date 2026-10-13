#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""compiler.py -- compile JMP_V11.tex et ECRIT ce que la compilation a rendu.

POURQUOI IL EXISTE (23/09/2026). Le 22/09 a 21:06, un edit a introduit un
appel a une macro qu'aucun producteur n'emet. Le PDF imprimait « for
1515specifications » et « and so specifications: ». La reference de phase est
restee a 0 ecart pendant douze heures : elle epingle des VALEURS, et un appel
sans valeur ne deplace aucune valeur. « Rien de publie n'a bouge » et « le
document compile proprement » etaient deux affirmations distinctes, traitees
comme une seule.

Le code de retour de la compilation n'est ecrit NULLE PART par pdflatex : le
log ne le porte pas. Ce script le capture, compte ce que le log dit, et ecrit
tout dans tools/COMPILATION.tsv, avec le md5 du .tex compile. reference_phase.py
le lit, et refuse de conclure si le .tex a change depuis.

Sequence : celle du registre -- .aux .bbl .blg .out .toc .brf effaces, puis
pdflatex x1 -> bibtex -> pdflatex x3, en -interaction=nonstopmode.

Usage:  python tools/compiler.py
Sortie: rc=0 si tout est propre (rc nuls, 0 indefinie, 0 renvoi ni citation
        non resolus), rc=1 sinon. Le TSV est ecrit dans les deux cas.
"""
import hashlib
import io
import os
import re
import subprocess
import sys
import datetime

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = 'JMP_V11'
SORTIE = os.path.join(ROOT, 'tools', 'COMPILATION.tsv')


def lance(cmd):
    r = subprocess.run(cmd, cwd=ROOT, stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
    return r.returncode


def comptes(log, blg):
    L = io.open(log, encoding='latin-1').read() if os.path.exists(log) else ''
    B = io.open(blg, encoding='latin-1').read() if os.path.exists(blg) else ''
    pages = re.findall(r'Output written on .*?\((\d+) pages', L)
    err_bib = re.search(r'\(There (?:was|were) (\d+|one) error message', B)
    return {
        'undefined_cs': len(re.findall(r'^! Undefined control sequence', L, re.M)),
        'erreurs_latex': len(re.findall(r'^! ', L, re.M)),
        'refs_non_resolus': len(re.findall(
            r"LaTeX Warning: Reference `[^']+' on page \d+\s+undefined", L)),
        'citations_non_resolues': len(re.findall(
            # Meme forme que les renvois. Une premiere version portait
            # `[^\n]*\n?[^\n]*undefined` : gourmande, elle enjambait jusqu'au
            # « undefined » de la ligne SUIVANTE et avalait une seconde
            # citation -- 1 compte pour 2 sur le cas synthetique du 23/09.
            r"(?:LaTeX|Package natbib) Warning: Citation `[^']+' on page \d+\s+undefined", L)),
        'rerun': len(re.findall(r'Rerun to get', L)),
        'erreurs_bibtex': (0 if not err_bib else
                           (1 if err_bib.group(1) == 'one' else int(err_bib.group(1)))),
        'pages': pages[-1] if pages else '',
    }


def main():
    tex = os.path.join(ROOT, BASE + '.tex')
    md5 = hashlib.md5(open(tex, 'rb').read()).hexdigest()
    for ext in ('.aux', '.bbl', '.blg', '.out', '.toc', '.brf'):
        p = os.path.join(ROOT, BASE + ext)
        if os.path.exists(p):
            os.remove(p)
    pdf = ['pdflatex', '-interaction=nonstopmode', BASE + '.tex']
    rc1 = lance(pdf)
    rcb = lance(['bibtex', BASE])
    rcs = [lance(pdf) for _ in range(3)]
    c = comptes(os.path.join(ROOT, BASE + '.log'), os.path.join(ROOT, BASE + '.blg'))
    if hashlib.md5(open(tex, 'rb').read()).hexdigest() != md5:
        print('[STOP] le .tex a change PENDANT la compilation. NON CONCLUANT.')
        return 2
    lignes = [('tex_md5', md5),
              ('date', datetime.datetime.now().isoformat(timespec='seconds')),
              ('rc_pdflatex_passe1', str(rc1)), ('rc_bibtex', str(rcb)),
              ('rc_pdflatex_finale', str(rcs[-1]))] + \
             [(k, str(v)) for k, v in c.items()]
    with io.open(SORTIE, 'w', encoding='utf-8', newline='') as fh:
        fh.write('# ecrit par tools/compiler.py -- NE PAS EDITER A LA MAIN.\n')
        for k, v in lignes:
            fh.write('%s\t%s\n' % (k, v))
    for k, v in lignes:
        print('   %-24s %s' % (k, v))
    propre = (rcs[-1] == 0 and rcb == 0 and c['undefined_cs'] == 0
              and c['refs_non_resolus'] == 0 and c['citations_non_resolues'] == 0
              and c['erreurs_bibtex'] == 0)
    print('VERDICT : %s' % ('compilation PROPRE' if propre else 'compilation NON PROPRE'))
    return 0 if propre else 1


if __name__ == '__main__':
    sys.exit(main())
