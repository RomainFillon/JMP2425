#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""scan_octets_controle.py -- les caracteres de controle qu'une chaine Python
non brute laisse dans un .tex genere.

LA FAMILLE, ET ELLE EST ENUMERABLE EXACTEMENT. En Python, SEPT sequences
`\` + lettre se transforment silencieusement en caractere de controle :

    \a  cloche          0x07        \n  saut de ligne   0x0A
    \b  retour arriere  0x08        \r  retour chariot  0x0D
    \f  saut de page    0x0C        \t  tabulation      0x09
    \v  tab. verticale  0x0B

Toutes les autres lettres laissent l'antislash intact. `\x` et `\u` levent une
SyntaxError -- c'est le comportement SUR : ils plantent au lieu de corrompre.

⇒ TOUTE MACRO LATEX DONT LE NOM COMMENCE PAR a, b, f, n, r, t OU v EST DETRUITE
  PAR UNE CHAINE PYTHON ORDINAIRE. Dans ce depot : \theta, \times, \textbf,
  \tau, \to (t) ; \begin, \bar, \bottomrule (b) ; \frac, \footnotesize (f) ;
  \newcommand, \num (n) ; \ref, \rho, \right (r) ; \alpha, \arraystretch (a) ;
  \vspace, \vdots (v). `\begin` ouvre chaque table generee.

CE QUE CE CONTROLE FAIT. Cinq des sept -- 0x07, 0x08, 0x0B, 0x0C, 0x0D -- se
balaient sur les octets, RETROACTIVEMENT : si \alpha ou \vspace a ete detruit
ailleurs dans l'arbre, il est dans le PDF en ce moment.

🚨 AVEC UNE EXCEPTION SUR 0x0D, ET SANS ELLE LE CONTROLE EST DU BRUIT. La
premisse « ces cinq octets ne sont jamais legitimes dans un .tex » est vraie des
quatre autres et FAUSSE de celui-la sur cette plateforme : ce depot est sous
Windows, la moitie de ses artefacts sont en CRLF, et 0x0D y precede CHAQUE 0x0A.
Mesure du 22/09/2026, premiere version : **53 fichiers en defaut, zero defaut
reel**. Seul un 0x0D ISOLE -- non suivi de 0x0A -- signe un `\r` ecrit a la
place d'un antislash.

✅ 26/09/2026 (lot 35) : \t EST DESORMAIS COUVERT -- 0x09 est interdit dans les
   .tex GENERES (voir INTERDITS). Six octets sur sept se balaient donc sur les
   octets ; seul \n reste au volet faible.

⛔ CE QU'IL NE FAIT PAS, ET C'EST ECRIT AVEC LUI. Il ne couvre pas \n,
   legitime comme saut de ligne. (Texte d'origine, pour \t aussi :) Leur residu
   a une signature FAIBLE -- une ligne qui s'ouvre sur ce qui reste d'un nom de
   macro : `ewcommand`, `extbf`, `imes`, `heta`, `ext`, `au`, `o`, `abular`. Ce
   volet-la RAPPORTE, il n'arrete pas : la liste est finie mais heuristique.

⛔ ET IL N'EST PAS LA VRAIE REPARATION. Un balayage d'artefacts attrape ce qui
   est deja casse. Ce qui ferme la famille est que les producteurs n'ecrivent
   jamais de LaTeX dans une chaine non brute -- verifiable par AST sur les
   emetteurs Python, et R a le meme piege avec "\\ref". Liste post-soumission 18.

🔒 27/09/2026 (lot 37) : LE CHAMP NE DEPEND PLUS DE L'ENVIRONNEMENT. Les
   producteurs R ecrivent dans DIR_OUT, qui vaut <racine>/output SAUF si
   JMP_OUTPUT_DIR est posee (00_parameters.R l.98-106) -- et ce sont TOUS leurs
   .tex, pas seulement ceux des re-emetteurs. Le scanner ne balayait que
   output/ en dur : sous redirection, les .tex generes sortaient de son champ en
   silence. Il balaie desormais output/ (ce que le manuscrit lit) ET le DIR_OUT
   que la meme regle resout, quand il differe. Temoin :
   tools/test_scan_octets_dirout.py, qui verifie aussi que la resolution Python
   rend le meme chemin que 00_parameters.R.

Usage:  python tools/scan_octets_controle.py            (output/ + DIR_OUT)
        python tools/scan_octets_controle.py <racine>   (cette racine seule)
Sortie: rc=0 si aucun octet de controle, rc=1 sinon, rc=2 si rien a lire.
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def dir_out():
    """DIR_OUT par la regle de 00_parameters.R : relatif => lu depuis la racine."""
    v = os.environ.get('JMP_OUTPUT_DIR', '')
    if not v:
        return os.path.join(ROOT, 'output')
    return v if re.match(r'^([A-Za-z]:|/|\\)', v) else os.path.join(ROOT, v)


if len(sys.argv) > 1:
    RACINES = [sys.argv[1]]
else:
    RACINES = [os.path.join(ROOT, 'output')]
    if os.path.normcase(os.path.abspath(dir_out())) != \
            os.path.normcase(os.path.abspath(RACINES[0])):
        RACINES.append(dir_out())

# Les six qu'aucun .tex GENERE ne contient. 0x0A seul en est EXCLU.
#
# 0x09 AJOUTE LE 26/09/2026 (lot 35). Il etait exclu comme « indentation
# legitime ». Or ce balayage ne porte QUE sur output/, c'est-a-dire sur des
# fichiers GENERES, et aucun producteur n'y indente par tabulation : mesure du
# 26/09, 0 tabulation dans tous les .tex de output/. L'exclusion ne protegeait
# donc rien, et elle laissait passer le cas reel : 01c_robustness.R portait
# "importer$\times$year" a UN antislash, soit « importer$<TAB>imes$year », et le
# volet faible ne cherche ses residus qu'en DEBUT de ligne. Temoin :
# tools/test_scan_octets_tab.py. La regle ne vaut que pour les fichiers
# generes : un .tex ecrit a la main peut indenter par tabulation.
INTERDITS = {0x07: r'\a (bell)', 0x08: r'\b (backspace)',
             0x09: r'\t (tab)',
             0x0B: r'\v (vertical tab)', 0x0C: r'\f (form feed)',
             0x0D: r'\r (carriage return)'}

# Le volet FAIBLE : ce qui reste d'un nom de macro apres \n ou \t.
RESIDUS = ('ewcommand', 'ewline', 'ewpage', 'um{', 'ote{',
           'extbf', 'extit', 'extsc', 'ext{', 'imes', 'heta', 'abular',
           'oprule', 'au{', 'o{')
RE_RESIDU = re.compile(r'(?m)^\s*(' + '|'.join(re.escape(r) for r in RESIDUS) + r')')


def fichiers(racine):
    for d, _, fs in os.walk(racine):
        for f in sorted(fs):
            if f.endswith('.tex'):
                yield os.path.join(d, f)


def main():
    cibles = []
    for racine in RACINES:
        # Une racine absente ou vide est NON CONCLUANTE, redirigee ou non :
        # une variable posee promet des ecritures a cet endroit.
        if not os.path.isdir(racine):
            print('[STOP] %s not found: check INCONCLUSIVE' % racine)
            return 2
        trouves = list(fichiers(racine))
        if not trouves:
            print('[STOP] no .tex under %s: check INCONCLUSIVE' % racine)
            return 2
        cibles += [(racine, p) for p in trouves]

    durs, faibles, crlf = [], [], 0
    for OUT, p in cibles:
        brut = open(p, 'rb').read()
        for code, nom in sorted(INTERDITS.items()):
            ##################################################################
            # 🚨 0x0D EXIGE UNE EXCEPTION, ET SANS ELLE CE CONTROLE EST DU
            #    BRUIT. Mesure du 22/09/2026 : la premiere version a rendu
            #    **53 fichiers** en defaut, tous sur 0x0D. Aucun n'etait un
            #    defaut : ce depot est sous Windows et la moitie de ses
            #    artefacts sont en CRLF, ou 0x0D precede CHAQUE 0x0A.
            #
            # ⇒ UN RETOUR CHARIOT SUIVI D'UN SAUT DE LIGNE EST UNE FIN DE LIGNE
            #   WINDOWS. Seul un 0x0D ISOLE -- non suivi de 0x0A -- est la
            #   signature d'un `\r` ecrit a la place d'un antislash.
            #
            # La premisse « ces cinq octets ne sont jamais legitimes dans un
            # .tex » est vraie des quatre autres et FAUSSE de celui-la, sur
            # cette plateforme. Un controle qui rend 53 faux positifs est un
            # controle qu'on desactive le lendemain.
            ##################################################################
            if code == 0x0D:
                pos = [i for i in range(len(brut))
                       if brut[i] == 0x0D
                       and not (i + 1 < len(brut) and brut[i + 1] == 0x0A)]
                crlf += brut.count(b'\r\n')
            else:
                pos = [i for i in range(len(brut)) if brut[i] == code]
            if pos:
                ligne = brut[:pos[0]].count(b'\n') + 1
                durs.append('%s: %d occurrence(s) of 0x%02X, %s -- '
                            'first at l.%d. A LaTeX macro lost its '
                            'backslash there in a non-raw Python string.'
                            % (os.path.relpath(p, OUT), len(pos), code, nom,
                               ligne))
        txt = io.open(p, encoding='utf-8', errors='replace').read()
        for m in RE_RESIDU.finditer(txt):
            faibles.append('%s l.%d: a line opens on %r -- possible residue '
                           'of a macro name cut by \\n or \\t.'
                           % (os.path.relpath(p, OUT),
                              txt[:m.start()].count('\n') + 1, m.group(1)))

    print('[octets] %d .tex file(s) scanned under %s; %d CRLF line ending(s) '
          'seen and IGNORED -- a 0x0D followed by 0x0A is a Windows line '
          'ending, not a cut'
          % (len(cibles), ' + '.join(os.path.relpath(r, ROOT) for r in RACINES),
             crlf))
    print('\n1. FORBIDDEN CONTROL BYTES: %d' % len(durs))
    for d in durs:
        print('   FORBIDDEN ' + d)
    if not durs:
        print('   none of the six bytes that \\a \\b \\t \\v \\f \\r produce')
    print('\n2. WEAK RESIDUES OF \\n AND \\t (report, not a stop): %d'
          % len(faibles))
    for f in faibles:
        print('   RESIDUE?  ' + f)
    print('\nVERDICT: %s' % ('%d file(s) carry a control byte'
                             % len({d.split(': ')[0] for d in durs})
                             if durs else
                             'no control byte in the generated .tex files'))
    return 1 if durs else 0


if __name__ == '__main__':
    sys.exit(main())
