#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""check_labels.py -- aucun numero de scenario dans un artefact genere.

CE QU'IL REFUSE
---------------
Le manuscrit ecrit Sc.1 = annuel et Sc.2 = journalier. Le CODE numerote
autrement : sc3 est le bras annuel, sc2 le bras journalier, sc1 la baseline
sans rechauffement. La correspondance entre les deux est la plus fragile du
paquet, et elle n'est ecrite nulle part qu'un controle puisse relire.

Un titre de panneau qui epelle « Sc.2 vs. Baseline » fige donc cette
correspondance DANS UNE IMAGE, ou plus rien ne la verifie : ni la compilation,
ni check_macros, ni une relecture du .tex, qui ne voit qu'un \includegraphics.
Trois figures l'ont portee jusqu'au 17/09/2026.

Les noms de marge -- « Daily total », « Annual total », « Daily extra » --
n'ont pas ce probleme : ils nomment l'objet, pas son rang.

DEUX PASSES, ET LA SECONDE NE CONCLUT PAS TOUJOURS
--------------------------------------------------
  1. LES PRODUCTEURS. Un littera Sc.<chiffre> dans une chaine de .R/.py qui
     produit une figure. C'est la passe qui fait foi : elle lit la source,
     donc elle voit ce qui SERA trace, pas seulement ce qui l'a ete.
  2. LES PDF EMIS. Le texte des flux de contenu, quand il est lisible. Les
     devices qui convertissent le texte en glyphes (cairo_pdf) ne le rendent
     pas : ces fichiers sont comptes NON CONCLUANTS et non « propres ». Un
     fichier illisible n'est pas un fichier sans defaut.

rc=0 vert ; rc=1 une occurrence trouvee ; rc=2 la passe 1 n'a pas pu tourner.
"""
import os
import re
import sys
import zlib

# La console Windows est en cp1252 : un caractere hors table ferait tomber le
# controle sur un UnicodeEncodeError, c'est-a-dire echouer pour une raison qui
# n'a rien a voir avec ce qu'il mesure.
try:
    sys.stdout.reconfigure(encoding='utf-8', errors='replace')
except (AttributeError, ValueError):
    pass

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'output')

# CE QUI EST INTERDIT, ET CE QUI NE L'EST PAS.
#
#   INTERDIT   « Sc.2 », « Sc. 2 », « Scenario 2 » -- des LIBELLES rendus, ou
#              le numero du scenario est presente au lecteur.
#   AUTORISE   « sc2 », « welfare_sc1_q50.csv », « sc2 - sc3 » -- les noms de
#              RUNS du solveur. C'est la numerotation du CODE, elle est
#              partout, et l'interdire rendrait le controle inutilisable.
#
# La distinction tient au S majuscule suivi d'un separateur, ou au mot entier
# suivi d'un espace. C'est etroit, et c'est voulu : un controle qui hurle sur
# trois cents noms de fichiers se fait desarmer au bout d'une semaine.
RE_SCEN = re.compile(r'\bSc\.\s?[0-9]|\b[Ss]cenarios?[\s~]+[0-9]')

# Les repertoires de producteurs. Les sauvegardes .bak_* sont hors flux vivant :
# les inclure ferait echouer le controle sur du code que plus rien n'execute.
SRC_DIRS = ['01_estimation', '02_GE', '03_inference', '04_contrast']
SRC_EXT = ('.R', '.r', '.py')

# Un producteur ne trace un libelle que dans une CHAINE. Le reste du fichier
# parle legitimement de sc1/sc2/sc3 : ce sont les noms des runs.
RE_STR = re.compile(r"'([^'\n]*)'|\"([^\"\n]*)\"")


def is_live(name):
    return '.bak' in name or name.endswith('~')


def scan_sources():
    hits = []
    n_files = 0
    for d in SRC_DIRS:
        p = os.path.join(ROOT, d)
        if not os.path.isdir(p):
            continue
        for fn in sorted(os.listdir(p)):
            if not fn.endswith(SRC_EXT) or is_live(fn):
                continue
            q = os.path.join(p, fn)
            n_files += 1
            with open(q, encoding='utf-8', errors='replace') as fh:
                for i, line in enumerate(fh, 1):
                    st = line.lstrip()
                    if st.startswith('#') or st.startswith('##'):
                        continue
                    for m in RE_STR.finditer(line):
                        val = m.group(1) if m.group(1) is not None else m.group(2)
                        if val and RE_SCEN.search(val):
                            hits.append((os.path.join(d, fn), i, val.strip()))
    return hits, n_files


# Le contenu d'une chaine PDF, sans l'operateur qui la montre.
RE_SHOW = re.compile(rb'\(((?:\\.|[^()\\])*)\)\s*(?:Tj|TJ|\]|\')')


def pdf_text(path):
    r"""Le texte des flux de contenu, ou None si le PDF ne le rend pas.

    DEUX PIEGES, MESURES LE 17/09/2026 SUR L'ANCIENNE FIGURE 9.

    1. matplotlib emet UN GLYPHE PAR `Tj` : `(S) Tj (c) Tj (.) Tj (2) Tj`.
       Joindre les chaines par une espace -- ce que faisait la premiere
       version -- donne « S c . 2 », que le motif ne reconnait pas. Le
       controle rendait donc 0 occurrence sur un fichier qui en portait une.
       Les chaines se concatenent SANS separateur.

    2. Un PDF ou le texte est converti en glyphes (cairo_pdf, Type 3) rend
       quelques chaines parasites et rien de lisible. Trouver « un Tj » ne
       prouve donc pas qu'on sait lire le fichier. On exige un TEMOIN : au
       moins un mot alphabetique de quatre lettres. Sans lui, la fonction
       rend None -- le fichier est NON CONCLUANT, pas propre.
    """
    try:
        with open(path, 'rb') as fh:
            raw = fh.read()
    except OSError:
        return None
    out = []
    for m in re.finditer(rb'stream\r?\n(.*?)endstream', raw, re.S):
        blob = m.group(1)
        try:
            blob = zlib.decompress(blob)
        except zlib.error:
            pass
        for t in RE_SHOW.finditer(blob):
            out.append(t.group(1))
    if not out:
        return None
    txt = b''.join(out).decode('latin-1', 'replace')
    txt = re.sub(r'\\([()\\])', r'\1', txt)
    # LE TEMOIN DE LISIBILITE.
    if not re.search(r'[A-Za-z]{4}', txt):
        return None
    return txt


def scan_pdfs():
    hits, unread, n = [], [], 0
    if not os.path.isdir(OUT):
        return hits, unread, 0
    for fn in sorted(os.listdir(OUT)):
        if not fn.lower().endswith('.pdf'):
            continue
        n += 1
        txt = pdf_text(os.path.join(OUT, fn))
        if txt is None:
            unread.append(fn)
            continue
        for m in RE_SCEN.finditer(txt):
            hits.append((fn, txt[max(0, m.start() - 30):m.end() + 30]))
    return hits, unread, n


def main():
    print('check_labels.py -- scenario numbers in the generated artefacts')
    src_hits, n_src = scan_sources()
    if n_src == 0:
        print('[STOP] no producer read: the source directories are '
              'not found. The check is INCONCLUSIVE.')
        return 2

    print('\n1. PRODUCER STRINGS: %d file(s) read, %d occurrence(s)'
          % (n_src, len(src_hits)))
    for f, i, v in src_hits:
        print('   FORBIDDEN %s:%d  %r' % (f, i, v))
    if not src_hits:
        print('   no producer label carries a scenario number')

    pdf_hits, unread, n_pdf = scan_pdfs()
    print('\n2. TEXT OF THE EMITTED PDFS: %d file(s), %d readable, '
          '%d occurrence(s)' % (n_pdf, n_pdf - len(unread), len(pdf_hits)))
    for f, ctx in pdf_hits:
        print('   FORBIDDEN %s  ...%s...' % (f, ctx.strip()))
    if unread:
        print('   INCONCLUSIVE on %d file(s) (text converted to glyphs):'
              % len(unread))
        for f in unread[:8]:
            print('      %s' % f)
        if len(unread) > 8:
            print('      ... and %d more' % (len(unread) - 8))
        print('   These files are NOT declared clean: pass 1 is '
              'what counts for them.')

    bad = src_hits or pdf_hits
    print('\nVERDICT: %s'
          % ('%d occurrence(s) of a scenario number' % len(bad) if bad
             else 'no scenario number in a generated artefact'))
    return 1 if bad else 0


if __name__ == '__main__':
    sys.exit(main())
