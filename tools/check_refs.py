#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""check_refs.py -- renvois et citations orphelins du manuscrit.

Trois controles, sur le FLUX VIVANT (lignes de commentaire et blocs
\begin{comment} retires) :

  1. \ref / \eqref vers un \label qui n'existe pas      -> « ?? » au rendu
  2. \label defini et jamais appele                      -> information seule
  3. \cite* vers une cle absente du .bib                 -> « [?] » au rendu

Aucun de ces trois ne fait echouer latexmk : le PDF sort, avec des « ?? ».
C'est pour cela qu'ils se mesurent ici.
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BS = chr(92)

RE_COMMENT = re.compile(r'(?m)^[ \t]*%.*$')
# Dans un motif, un antislash LITTERAL s'ecrit avec DEUX caracteres : BS seul
# echappe la parenthese qui suit et casse le lookbehind.
RE_INLINE = re.compile(r'(?<!' + BS * 2 + r')%.*$', re.M)
RE_COMMENT_ENV = re.compile(r'\\begin\{comment\}.*?\\end\{comment\}', re.S)
RE_LABEL = re.compile(r'\\label\{([^{}]+)\}')
RE_REF = re.compile(r'\\(?:page)?(?:eq)?ref\{([^{}]+)\}')
RE_CITE = re.compile(r'\\cite[a-zA-Z]*\s*(?:\[[^\]]*\])*\s*\{([^{}]+)\}')
RE_BIBKEY = re.compile(r'(?m)^\s*@[A-Za-z]+\s*\{\s*([^,\s]+)\s*,')


def strip(t):
    t = RE_COMMENT_ENV.sub(' ', t)
    t = RE_COMMENT.sub('', t)
    return RE_INLINE.sub('', t)


def main():
    tex_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'JMP_V11.tex')
    bib_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, 'References_JMP.bib')
    raw = io.open(tex_path, encoding='utf-8', errors='replace').read()
    tex = strip(raw)
    lines = raw.split('\n')

    # LES \\label DES FICHIERS INCLUS COMPTENT. Les tables generees
    # portent leur propre \\label dans output/tab_*.tex ; ne lire que le
    # manuscrit les fait passer pour orphelines alors qu'elles sont definies
    # a l'inclusion. Un verificateur qui sur-rapporte se fait ignorer aussi
    # surement qu'un qui sous-rapporte.
    # Motif bati depuis chr(92) : la valeur de la chaine est exactement
    # le motif voulu, sans compter les niveaux d'echappement.
    _re_input = re.compile(chr(92) * 2 + 'input' + chr(92) + '{([^{}]+)'
                           + chr(92) + '}')
    inclus = []
    base = os.path.dirname(os.path.abspath(tex_path))
    for rel in _re_input.findall(tex):
        f = os.path.join(base, rel if rel.endswith('.tex') else rel + '.tex')
        if os.path.exists(f):
            inclus.append(io.open(f, encoding='utf-8', errors='replace').read())
    print('[refs] %d included file(s) read for their labels' % len(inclus))
    labels = set(RE_LABEL.findall(tex + chr(10).join(inclus)))

    refs = {}
    for m in RE_REF.finditer(tex):
        refs.setdefault(m.group(1), 0)
        refs[m.group(1)] += 1

    bib = set()
    if os.path.exists(bib_path):
        bib = set(RE_BIBKEY.findall(
            io.open(bib_path, encoding='utf-8', errors='replace').read()))

    cites = set()
    for m in RE_CITE.finditer(tex):
        for k in m.group(1).split(','):
            if k.strip():
                cites.add(k.strip())

    def where(pat):
        return [i for i, l in enumerate(lines, 1)
                if pat in l and not l.lstrip().startswith('%')]

    orph = sorted(k for k in refs if k not in labels)
    print("[refs] %d labels, %d ref keys, %d bib keys" % (len(labels), len(refs), len(bib)))
    print("\n1. REF WITHOUT LABEL: %d" % len(orph))
    for k in orph:
        print("   ORPHAN    %-28s lines %s"
              % (k, ', '.join(map(str, where('ref{' + k + '}')))))

    unused = sorted(k for k in labels if k not in refs)
    print("\n2. LABEL NEVER CALLED: %d" % len(unused))
    for k in unused:
        print("   UNUSED    %s" % k)

    badcite = sorted(k for k in cites if bib and k not in bib)
    print("\n3. CITATION WITHOUT BIB ENTRY: %d" % len(badcite))
    for k in badcite:
        print("   MISSING   %-34s lines %s"
              % (k, ', '.join(map(str, where(k)))))

    ##########################################################################
    # 4. ENTREE BIB JAMAIS CITEE -- pose le 22/09/2026.
    #
    # LE DEFAUT QU'IL REND VISIBLE. Le controle 3 regardait dans UN SEUL SENS :
    # une citation sans entree. Le sens inverse -- une entree que plus rien ne
    # cite -- n'etait vu par personne, et c'est exactement ce qu'un lot de
    # RETRAIT produit. Mesure du 22/09 : la clause de l'equation (5) a retire
    # les deux dernieres citations de `caliendo2015estimates`, dont l'une
    # attribuait a tort theta = 6,5. L'entree est devenue orpheline dans le
    # meme geste, et rien ne l'aurait dit.
    #
    # ⚠️ IL RAPPORTE, IL N'ARRETE PAS -- meme polarite que le controle 2. Une
    # entree non citee ne degrade aucun rendu : bibtex ne l'imprime pas. Ce
    # qu'elle signale est editorial (« ce papier ne s'appuie plus sur cette
    # reference »), et c'est une decision d'auteur, pas un defaut.
    ##########################################################################
    # ⛔ PREMIERE VERSION, JETEE, ET LA RAISON VAUT D'ETRE ECRITE. Elle listait
    # TOUTE entree du .bib jamais citee : mesure du 22/09, **80 entrees**. Le
    # `.bib` est une BIBLIOTHEQUE de travail, bien plus large que le jeu de
    # citations du papier ; « jamais citee » y est l'etat NORMAL, et
    # `caliendo2015estimates` s'y serait perdue parmi 79 autres. Un controle qui
    # rend 80 evenements dont un seul compte est un controle qu'on cesse de
    # lire. Ce qu'on veut n'est pas l'etat, c'est la TRANSITION.
    nocite = sorted(k for k in bib if k not in cites)
    base_p = os.path.join(ROOT, 'output', 'cited_keys.txt')
    base = set()
    if os.path.exists(base_p):
        base = {l.strip() for l in io.open(base_p, encoding='utf-8')
                if l.strip() and not l.startswith('#')}
    partis = sorted(base - cites)
    arrives = sorted(cites - base) if base else []
    print("\n4. SET OF CITATIONS: %d cited, %d bib entries not cited"
          % (len(cites), len(nocite)))
    if not base:
        print("   [!] no reference set: this check is INCONCLUSIVE "
              "this time. Base written -> output/cited_keys.txt")
    else:
        for k in partis:
            print("   NO LONGER CITED   %-30s -- its .bib entry becomes "
                  "an orphan" % k)
        for k in arrives:
            print("   NEW               %s" % k)
        if not partis and not arrives:
            print("   set of citations unchanged since the last reference")
    with io.open(base_p, 'w', encoding='utf-8') as fh:
        fh.write('# cles citees par le manuscrit -- reference du controle 4 de\n'
                 '# check_refs.py. Regeneree a chaque passage : le controle\n'
                 '# rapporte la TRANSITION, pas l etat.\n')
        for k in sorted(cites):
            fh.write(k + '\n')

    ##########################################################################
    # 5. UNE COMMANDE AMPUTEE DE SON ANTISLASH -- pose le 22/09/2026.
    #
    # LE DEFAUT, ET IL ETAIT SUR DISQUE QUAND CE CONTROLE A ETE ECRIT. Un
    # emetteur Python ecrivait "Table~\ref{tab:ge}" dans une chaine NORMALE :
    # Python y lit `\r` comme un RETOUR CHARIOT, et l'artefact publie portait
    # « Table~ef{tab:ge} ».
    #
    # ⛔ POURQUOI AUCUN AUTRE CONTROLE NE PEUT LE VOIR. Un renvoi ampute ne
    # DEVIENT PAS un renvoi casse : il cesse d'etre un renvoi. LaTeX imprime
    # « ef{tab:ge} » comme du texte ordinaire, sans avertissement ; le controle 1
    # ne le compte pas, puisqu'il cherche `\ref{`. Le manuscrit compile avec
    # zero renvoi non resolu et une reference litterale dans une legende.
    #
    # LE MOTIF : un `~` suivi de lettres puis d'une accolade. En LaTeX c'est
    # l'insecable d'un « Table~\ref{...} », « Figure~\ref{...} »,
    # « Section~\ref{...} » -- il n'y a pas d'autre raison d'ecrire cela.
    # Balaye les .tex GENERES, ou l'erreur naît, en plus du manuscrit.
    ##########################################################################
    # Le `~` et la commande amputee peuvent etre separes par un BLANC : le
    # retour chariot ecrit a la place de l'antislash devient un saut de ligne,
    # puis une espace apres normalisation. Sans `\s*` ce controle rate
    # exactement le defaut qui l'a fait naitre -- mesure du 22/09.
    RE_AMPUTE = re.compile(r'~\s*([a-zA-Z]+)\{')
    ampute = []
    cibles = [(os.path.basename(tex_path), tex)]
    outdir = os.path.join(ROOT, 'output')
    if os.path.isdir(outdir):
        for f in sorted(os.listdir(outdir)):
            if f.endswith('.tex') and not f.startswith('values'):
                cibles.append((f, strip(io.open(os.path.join(outdir, f),
                               encoding='utf-8', errors='replace').read())))
    # ⛔ BLANCS NORMALISES D'ABORD, ET C'EST LE TEST QUI L'A IMPOSE.
    # L'amputation observee le 22/09 a produit un SAUT DE LIGNE la ou le
    # retour chariot a ete ecrit : le fichier portait « Table~ » en fin de
    # ligne et « ef{tab:ge}. » au debut de la suivante. Un motif a fenetre
    # bornee ne franchit pas un saut de ligne, et le controle passait au VERT
    # sur le defaut meme qui l'a fait naitre -- il etait sur disque pendant
    # que je le declarais propre.
    for nom, corps in cibles:
        corps = re.sub(r'\s+', ' ', corps)
        for m in RE_AMPUTE.finditer(corps):
            ampute.append('%s: "~%s{" -- command stripped of its backslash. '
                          'It breaks nothing: it DISAPPEARS, and the text is '
                          'printed as is.' % (nom, m.group(1)))
    print("\n5. COMMAND STRIPPED OF ITS BACKSLASH: %d" % len(ampute))
    for a in ampute:
        print("   STRIPPED  " + a)

    return 1 if (orph or badcite or ampute) else 0


if __name__ == '__main__':
    sys.exit(main())
