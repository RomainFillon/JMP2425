#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""scan_litteraux.py -- les nombres ecrits en dur dans le corps du manuscrit.

CE QU'IL FAIT, ET CE QU'IL NE FAIT PAS
--------------------------------------
Il liste les litteraux numeriques du FLUX VIVANT -- hors lignes de commentaire,
hors blocs \begin{comment} -- et les classe en deux tas :

  CONSTANTE  un parametre du modele ou un repere de lecture : theta, mu, une
             annee, un seuil de latitude, une division ISIC, un degre Celsius.
             Ils restent en dur, c'est leur place.
  RESULTAT   tout le reste : une magnitude, une part, un intervalle. Aucun ne
             doit rester en dur -- c'est la regle du paquet.

Il ne DECIDE pas : il propose un classement et montre la ligne, pour que le
tri se fasse sur le contexte. Un classement automatique qui se tromperait en
silence serait pire que pas de scan du tout.

⚠️ Le detecteur est teste par MUTATION : --selftest injecte un litteral qui
DOIT etre attrape, et un qui doit etre ignore. Un detecteur qu'on n'a pas vu
echouer ne mesure rien -- ce paquet a deja publie un « 0 Overfull » contre 14.

Usage:  python tools/scan_litteraux.py [manuscrit.tex] [--selftest]
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Nombres a ignorer : reperes de lecture et parametres declares du modele.
CONSTANTES = {
    '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '12', '100',
    '6.5', '0.6',                      # theta, mu (caliendo2015estimates)
    '45', '55', '60', '66.5', '35', '30', '25', '20', '15',   # latitudes / degres
    '1995', '2010', '2014', '2015', '2019', '2030', '2050', '2095', '2100',
    '1990', '2005', '2020',
    '0.05', '0.95', '0.01', '0.001', '1.96',                  # seuils statistiques
    '02', '05', '13', '14', '18', '33', '40',                 # divisions ISIC
}

RE_COMMENT = re.compile(r'(?m)^[ \t]*%.*$')
RE_INLINE_COMMENT = re.compile(r'(?<!\\)%.*$', re.M)
RE_COMMENT_ENV = re.compile(r'\\begin\{comment\}.*?\\end\{comment\}', re.S)
# Un nombre : signe optionnel (ASCII ou moins typographique), chiffres, decimales.
RE_NUM = re.compile(r'(?<![\w.])[-\u2212+]?\d+(?:[.,]\d+)?(?![\w.])')
RE_MACRO_ARG = re.compile(r'\\(?:label|ref|cite[a-z]*|input|includegraphics|'
                          r'newcommand|providecommand|usepackage|documentclass|'
                          r'begin|end|hspace|vspace|width|textwidth|linewidth)'
                          r'\s*(?:\[[^\]]*\])?\s*(?:\{[^{}]*\})?')


def strip(tex):
    tex = RE_COMMENT_ENV.sub(' ', tex)
    tex = RE_COMMENT.sub('', tex)
    tex = RE_INLINE_COMMENT.sub('', tex)
    tex = RE_MACRO_ARG.sub(' ', tex)
    return tex


def scan(tex):
    out = []
    for i, line in enumerate(strip(tex).split('\n'), 1):
        if not line.strip():
            continue
        for m in RE_NUM.finditer(line):
            raw = m.group(0)
            key = raw.lstrip('+-\u2212')
            if key in CONSTANTES:
                continue
            out.append((i, raw, line.strip()))
    return out




##################################################################################
# UNE MAGNITUDE ECRITE EN TOUTES LETTRES EST UN LITTERAL. Ajoute le 20/09/2026.
#
# LE TROU, ET IL EST DE NATURE. Ce detecteur cherchait des CHIFFRES. Il ne pouvait
# donc pas voir la phrase de §2.5 l.475 : « The two functions differ by MORE THAN
# A LOG POINT over the cold range ». C'est un litteral de resultat -- il enonce une
# magnitude, dans une unite du papier, sans passer par une macro -- et apres la
# division par theta du 19/09 la valeur vaut environ 0,25 point de log. Le chiffre
# est faux et aucun garde ne pouvait le dire.
#
# CE QUI DOIT SE DECLENCHER : une quantite en toutes lettres COLLEE a une unite du
# papier -- « a log point », « half a log point », « two degrees », « a one-degree
# rise ». C'est une magnitude, elle depend de l'echelle, elle doit venir d'une macro.
#
# CE QUI NE DOIT PAS SE DECLENCHER, ET C'EST LE CAS NEGATIF OBLIGATOIRE :
#   - « an order of magnitude », « orders of magnitude » : un RAPPORT, invariant a
#     l'echelle, et la phrase de §2.5 en contient un dans sa seconde moitie, qui
#     lui reste vraie. Un detecteur qui attrape les deux moities ne distingue plus
#     ce qui tombe de ce qui tient.
#   - « by a factor of ... » : un rapport aussi.
#   - une unite SANS magnitude : « reported in log points », « per degree ».
##################################################################################
_QUANT = (r'(?:an?|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve|'
          r'half|quarter|third|tenth|several|dozens?|couple|few)')
_UNITE = (r'(?:log[ -]points?|percentage[ -]points?|degrees?(?:[ -]Celsius)?|'
          r'log[ -]units?)')
# La quantite precede l'unite, eventuellement via « of a » ou un trait d'union.
RE_MOT = re.compile(
    r'(?<![A-Za-z])'
    r'((?:more|less|fewer|no)\s+than\s+)?'          # « more than a log point »
    r'(?:' + _QUANT + r')'
    r'(?:\s+(?:of\s+a|and\s+a)\b)?'
    r'[ -]'
    r'(?:' + _UNITE + r')'
    r'(?![A-Za-z])'
    # « three or four degrees OF FREEDOM » n'est pas une temperature : c'est le
    # parametre d'une spline. Mesure du 20/09 : sans cette exclusion le detecteur
    # rendait deux faux positifs (l.368 et l.1236) pour un vrai. Une garde qui se
    # declenche a tort echange un faux negatif contre un faux positif.
    r'(?!\s+of\s+freedom)', re.I)
# Les tournures de RAPPORT, retirees avant l'examen : elles sont invariantes a
# l'echelle et ne sont donc pas des litteraux de resultat.
RE_RAPPORT = re.compile(
    r'(?<![A-Za-z])(?:an?\s+order\s+of\s+magnitude|orders\s+of\s+magnitude|'
    r'(?:by\s+)?a\s+factor\s+of)(?![A-Za-z])', re.I)


def scan_mots(tex):
    """Les magnitudes ecrites en toutes lettres. Meme sortie que scan()."""
    out = []
    for i, line in enumerate(strip(tex).split('\n'), 1):
        if not line.strip():
            continue
        masque = RE_RAPPORT.sub(lambda m: ' ' * len(m.group(0)), line)
        for m in RE_MOT.finditer(masque):
            out.append((i, m.group(0).strip(), line.strip()))
    return out


def selftest_mots():
    """MUTATION : quatre cas qui DOIVENT matcher, six qui NE DOIVENT PAS."""
    doivent = [
        ("the real defect of 20/09",
         r"The two functions differ by more than a log point over the cold range."),
        ("a half unit", r"The gap is half a log point at the margin."),
        ("a plural", r"It falls by two log points between the two bounds."),
        ("a degree in words", r"the effect of a one-degree rise in the mean"),
    ]
    doivent_pas = [
        ("a RATIO, invariant",
         r"and by an order of magnitude less over $20$--$30^\circ$C."),
        ("a plural ratio", r"two orders of magnitude smaller"),
        ("a factor", r"exceeds it by a factor of \EMPmdeRatio."),
        ("a unit WITHOUT magnitude", r"The three margins are reported in log points."),
        ("degrees of freedom", r"natural splines with three or four degrees of freedom"),
        ("in a comment", "%% The two functions differ by more than a log point."),
    ]
    ok = True
    for nom, s in doivent:
        if scan_mots(s):
            print("  OK    %-28s caught: %r" % (nom, scan_mots(s)[0][1]))
        else:
            print("  FAIL  %-28s NOT caught: %s" % (nom, s[:60])); ok = False
    for nom, s in doivent_pas:
        h = scan_mots(s)
        if h:
            print("  FAIL  %-28s wrongly caught: %r" % (nom, h[0][1])); ok = False
        else:
            print("  OK    %-28s ignored" % nom)
    return ok


def selftest():
    """MUTATION : un cas qui DOIT matcher, un qui NE DOIT PAS."""
    doit = r"The effect is $-17.6$ log points at the bound."
    doit_pas_commentaire = "%% The effect is $-17.6$ log points."
    doit_pas_constante = r"with $\theta=6.5$ and $\mu=0.6$."
    ok = True
    a = scan(doit)
    if not any(x[1] in ('-17.6', '\u221217.6') for x in a):
        print("  FAIL: -17.6 in live text is NOT caught."); ok = False
    else:
        print("  OK: -17.6 caught in the live text")
    b = scan(doit_pas_commentaire)
    if b:
        print("  FAIL: a literal in a comment is caught: %r" % (b,)); ok = False
    else:
        print("  OK: literal in a comment ignored")
    c = scan(doit_pas_constante)
    if c:
        print("  FAIL: a declared parameter is caught: %r" % (c,)); ok = False
    else:
        print("  OK: theta/mu ignored")
    return ok


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    if '--selftest' in sys.argv:
        print("--- mutation: NUMERIC literals ---")
        a = selftest()
        print("--- mutation: magnitudes IN WORDS ---")
        b = selftest_mots()
        return 0 if (a and b) else 1
    tex = args[0] if args else os.path.join(ROOT, 'JMP_V11.tex')
    body = io.open(tex, encoding='utf-8', errors='replace').read()
    hits = scan(body)
    print("[scan] %s: %d NUMERIC literal(s) outside the declared constants\n"
          % (os.path.basename(tex), len(hits)))
    for ln, raw, ctx in hits:
        print("  l.%-5d %-10s %s" % (ln, raw, ctx[:104]))
    mots = scan_mots(body)
    print("\n[scan] %s: %d magnitude(s) IN WORDS\n"
          % (os.path.basename(tex), len(mots)))
    for ln, raw, ctx in mots:
        print("  l.%-5d %-24s %s" % (ln, raw, ctx[:90]))
    return 0


if __name__ == '__main__':
    sys.exit(main())
