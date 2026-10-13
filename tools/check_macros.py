#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""check_macros.py -- test d'acceptation de la couche d'emission.

Neuf controles. Les HUIT premiers ne cherchent aucun litteral dans le code ; le
neuvieme ne peut faire que cela, et c'est son objet -- il tient la regle du
18/09/2026, « la grille de specifications est definie en un seul endroit » :

  1. APPELEES MAIS NON DEFINIES -> rendront "??" ou "Undefined control sequence".
  2. DEFINIES MAIS JAMAIS APPELEES -> un exhibit genere que rien ne consomme.
  3. DEFINIES A \PROVISOIRE ET APPELEES -> le PDF imprime "??". REFUS (rc=1), macro
     nommee. Lot 85x (11/10/2026) : ce controle listait sans refuser (« voulu ») ;
     le clone du 10/10 a compile PROPRE en imprimant « ?? million » pour
     \GEmigInternalCrossM et \GEmigInternalCrossShare, dont le producteur etait hors
     chaine. Une macro \PROVISOIRE que le .tex n'appelle pas reste permise.
     Connu-faux : tools/test_check_macros_provisoire.py.

Les lignes de commentaire LaTeX entieres sont retirees avant le comptage : un
%% NEW-MACRO : \EMPfoo est une note de travail, pas un appel.

Usage:  python tools/check_macros.py [manuscrit.tex]
Sortie: rc=0 si aucune macro appelee n'est indefinie ni \PROVISOIRE, rc=1 sinon.
        CHECK_MACROS_SOURCES (chemins separes par os.pathsep) remplace les deux
        fichiers de macros : pour les tests seulement.
"""
import io
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'JMP_V11.tex')
SOURCES = [os.path.join(ROOT, 'output', 'values.tex'),
           os.path.join(ROOT, 'output', 'values_empirics.tex')]
if os.environ.get('CHECK_MACROS_SOURCES'):
    SOURCES = os.environ['CHECK_MACROS_SOURCES'].split(os.pathsep)

PREFIXES = ('GE', 'EMP', 'CALIB')
BS = chr(92)

RE_DEF = re.compile(r'newcommand\{' + BS + BS + r'([A-Za-z]+)\}')
RE_CALL = re.compile(BS + BS + r'((?:' + '|'.join(PREFIXES) + r')[A-Za-z]+)')
RE_PROV = re.compile(r'newcommand\{' + BS + BS + r'([A-Za-z]+)\}\{'
                     + BS + BS + r'PROVISOIRE\}')
RE_COMMENT = re.compile(r'(?m)^[ \t]*%.*$')
RE_INPUT = re.compile(r'input\{([^{}]+\.tex)\}')


def read(p):
    if not os.path.exists(p):
        return ''
    return io.open(p, encoding='utf-8', errors='replace').read()




# ======================================================================
# CONTROLE 10 -- UNE PAIRE DE BORNES EST RENDUE DANS LA FORME DE SA NATURE
# ======================================================================
# CE QU'IL FERME. Le suffixe `Lo`/`Hi` porte DEUX conventions dans ce paquet :
# une borne de confiance (`ci95_lo/hi`) et un min/max sur un ensemble. Releve du
# 20/09/2026 : 39 paires emises, 14 citees dans le flux vivant. Le `.tex` les
# rend identiquement, et ce qui les distingue est la TOURNURE, pas le nom --
# `$[a,b]$` pour un intervalle, « A to B » pour une etendue. Une discipline qui
# ne vit que dans la prose est une discipline qu'un futur auteur ne voit pas.
#
# ⛔ ET ELLE AVAIT DEJA MORDU : `\EMPohSlopeLo/Hi` sont le min et le max de
# `slope_q50` sur les schemas identifies -- une ETENDUE -- et le manuscrit les
# rendait `$[\EMPohSlopeLo,\EMPohSlopeHi]$`, c'est-a-dire comme un intervalle de
# confiance sur la pente chaude. Ce n'en est pas un.
#
# POURQUOI LA GARDE ET PAS LE RENOMMAGE. Renommer neuf paires citees dans une
# prose actuellement correcte ne corrige rien aujourd'hui et ouvre neuf surfaces
# d'edition. La garde, elle, protege AUSSI la paire qui n'existe pas encore.
# Avec elle un mauvais nom devient inoffensif ; sans elle, le renommage laisse
# la prochaine paire reproduire le piege.
#
# L'EMETTEUR DECLARE, LE CONTROLE VERIFIE -- meme dispositif que
# values_provenance.tsv et le controle 5. Une paire CITEE et NON DECLAREE
# echoue ; une paire declaree et non citee est ignoree, parce qu'elle ne peut
# tromper personne.
PAIR_SUF = (('Lo', 'Hi'), ('Min', 'Max'))


# Les deux gabarits de message. Ecrits SANS apostrophe et hors de la boucle :
# une concatenation d apostrophes par chr(39) casse la concatenation implicite
# de litteraux, et le % ne se lie plus qu au dernier fragment (mesure, 20/09).
_MSG_NODECL = ('pair CITED and NOT DECLARED: %s%s/%s%s -- to be added to '
               'PAIR_NATURE in emit_values_empirics.py')
_MSG_FORME = ('%s%s/%s%s is declared %r and rendered %r at occurrence %d/%d '
              '(link %r). An interval is rendered in brackets with a '
              'comma; a range is rendered with to or --.')


def _pair_form(txt, lo, hi):
    """La forme dans laquelle le .tex relie deux bornes : 'interval', 'range',
    ou None si elles n'apparaissent pas ensemble.

    `interval` : un crochet ouvrant precede la premiere borne et une virgule les
                 relie -- $[\\Lo,\\Hi]$.
    `range`    : tout autre lien court -- « to », « -- », « and » apres
                 « between ».
    """
    B = chr(92)
    pat = re.compile(re.escape(B + lo) + r'\{?\}?(.{0,24}?)' + re.escape(B + hi))
    ## (!) TOUTES LES OCCURRENCES, PAS LA PREMIERE. Mesure du 20/09/2026 :
    ## EMPohSlopeLo/Hi est citee DEUX fois -- l.449 et l.1183 -- et une paire
    ## peut etre juste a un endroit et fausse a l'autre. Un controle cale sur
    ## la premiere occurrence absout la seconde.
    formes = []
    for m in pat.finditer(txt):
        lien = m.group(1)
        avant = txt[max(0, m.start() - 12):m.start()]
        crochet = '[' in avant.split(B + lo)[0][-6:]
        if crochet and re.match(r'^\s*[,;]\s*\$?\s*$', lien.replace('$', ' ')):
            formes.append(('interval', lien))
        else:
            formes.append(('range', lien))
    return formes


def controle_paires(tex_vivant, pairs_path):
    """Rend (lignes de rapport, liste d'erreurs).

    (!) LES BLANCS SONT NORMALISES D'ABORD. Le `.tex` est enroule : deux bornes
    d'une meme paire tombent souvent de part et d'autre d'un saut de ligne, et
    un motif a fenetre bornee ne franchit pas un `
`. Sans cette ligne le
    controle passait au VERT sur deux paires qu'il n'avait pas regardees --
    exactement le mode de defaillance que tout ce dossier ferme.
    """
    B = chr(92)
    tex_vivant = re.sub(r'\s+', ' ', tex_vivant)
    out, bad = [], []
    if not os.path.exists(pairs_path):
        return ([], ['values_pairs.tsv missing -- re-run emit_values_empirics.py'])
    decl = {}
    for i, l in enumerate(io.open(pairs_path, encoding='utf-8').read().split('\n')):
        if i == 0 or not l.strip():
            continue
        c = l.split('\t')
        if len(c) >= 2:
            decl[c[0]] = c[1].strip()
    appelees = set(re.findall(re.escape(B) + r'((?:EMP|GE)[A-Za-z]+)', tex_vivant))
    vues = 0
    for lo_s, hi_s in PAIR_SUF:
        for n in sorted(appelees):
            if not n.endswith(lo_s) or len(n) <= len(lo_s) + 2:
                continue
            base = n[:-len(lo_s)]
            if base + hi_s not in appelees:
                continue
            _f = _pair_form(tex_vivant, base + lo_s, base + hi_s)
            if not _f:
                continue
            vues += 1
            if base not in decl:
                bad.append(_MSG_NODECL % (base, lo_s, base, hi_s))
                continue
            attendu = decl[base]
            mauvaises = [(k, f, g) for k, (f, g) in enumerate(_f, 1)
                         if f != attendu]
            for k, f, g in mauvaises:
                bad.append(_MSG_FORME
                           % (base, lo_s, base, hi_s, attendu, f, k, len(_f), g))
            if not mauvaises:
                out.append('   %-28s %-8s %d occurrence(s), all conforming'
                           % (base, attendu, len(_f)))
    out.insert(0, "   %d pair(s) cited in the live text, %d declared"
               % (vues, len(decl)))
    return out, bad


##############################################################################
# CONTROLE 11 -- UNE CONSTANTE IMPOSEE CITE SA SOURCE DANS LA MEME PHRASE.
# Pose le 22/09/2026.
#
# LE DEFAUT QU'IL FERME, ET IL S'ETAIT PRODUIT. `\GEtheta` remplace le litteral
# 6,5 a deux endroits du manuscrit (§2.3 et le supplement A). La VALEUR ne peut
# plus diverger : elle est lue dans specification_study.csv. L'ATTRIBUTION, elle,
# reste du texte libre -- et c'est exactement elle qui avait diverge : le
# supplement ecrivait \citep{caliendo2015estimates} pour une valeur que le
# fichier de configuration annote « Desmet et al., 2018 ».
#
# ⇒ LA MACRO FERME LA VALEUR, PAS LA PHRASE. Le controle exige donc que chaque
#   phrase qui appelle la macro porte AUSSI la cle de citation attendue.
#
# ⛔ LIMITE, ECRITE AVEC LE CONTROLE -- ET CORRIGEE PAR LA MUTATION.
#   Premiere redaction : « il verifie une co-occurrence, pas une attribution ;
#   changer la cle aux DEUX sites passerait au vert ». C'EST FAUX, et le test
#   P3 l'a montre : la cle attendue est EPINGLEE dans CITATIONS_ATTENDUES, donc
#   changer les deux sites echoue aussi. Le controle est plus fort que la
#   description que j'en avais faite.
#
#   LA VRAIE LIMITE : il ne sait pas si la cle EPINGLEE est la bonne. Si
#   `desmet2018geography` etait elle-meme une mauvaise attribution, le controle
#   l'imposerait partout avec la meme energie. La justesse de l'epingle se
#   verifie EN SOURCE, une fois -- ce qui a ete fait le 22/09 sur
#   specification_study.csv -- et pas par cet outil.
#
#   ⚠️ Et une garde decrite plus FAIBLE qu'elle n'est reste une garde mal
#   decrite : la prochaine personne lui fera moins confiance qu'elle ne merite,
#   ou reecrira un controle qui existe deja.
##############################################################################
# ✅ LA JUSTESSE DES DEUX CLES A ETE ETABLIE EN SOURCE LE 22/09/2026, et c'est
#    inscrit ICI parce que ce controle EPINGLE une cle sans pouvoir la juger.
#    Sans cette note, l'epingle reposerait sur une justification restee dans une
#    conversation -- le defaut que le registre porte deja sous son nom.
#
#    theta = 6,5 -> desmet2018geography. Desmet, Nagy et Rossi-Hansberg (2018),
#    calibration : « We choose a trade elasticity, nu, equal to 6.5, somewhere in
#    the middle between the 8.3 value estimated by Eaton and Kortum (2002) and
#    the 4.6 value estimated by Simonovska and Waugh (2014). » Les auteurs la
#    notent nu ; c'est le meme objet -- le parametre de forme de Frechet -- et
#    §2.1 definit deja theta comme tel.
#
#    alpha = 0,6 -> desmet2015spatial (et NON caliendo2015estimates). Desmet et
#    Rossi-Hansberg (2015), calibration : « Estimates of labor shares vary
#    between 0.45 and 0.7 [...] we choose labor shares to be the same in both
#    sectors and set them equal to 0.6. »
#    ⛔ ET ON NE PEUT PAS RABATTRE alpha SUR LA SOURCE DE THETA : DNRH 2018 fixe
#    mu = 0,8 et precise que ce n'est PAS une part du travail (« should be
#    interpreted as the nonland share »). La calibration de ce papier est mixte,
#    et elle a raison de l'etre. alpha n'est pas encore une macro, donc pas
#    encore dans CITATIONS_ATTENDUES ; l'entree l'attend.
CITATIONS_ATTENDUES = {'GEtheta': 'desmet2018geography'}

RE_PHRASE = re.compile(r'(?<=[.!?])\s+(?=[A-Z' + BS + BS + r'])')


def controle_citations(tex_vivant):
    bad, rap = [], []
    plat = re.sub(r'\s+', ' ', tex_vivant)
    phrases = RE_PHRASE.split(plat)
    for macro, cle in sorted(CITATIONS_ATTENDUES.items()):
        appel = re.compile(BS + BS + macro + r'(?![A-Za-z])')
        sites = [p for p in phrases if appel.search(p)]
        if not sites:
            rap.append('   %-12s not called: nothing to check'
                       % ('\\' + macro))
            continue
        sans = [p for p in sites if cle not in p]
        if sans:
            for p in sans:
                bad.append('%s is called in a sentence that does not cite '
                           '%r. An IMPOSED constant names its source where '
                           'it is read, otherwise two sites diverge on the '
                           'attribution -- which happened. Sentence: %s'
                           % ('\\' + macro, cle, '...' + p.strip()[:120]))
        else:
            rap.append('   %-12s %d site(s), all citing %s'
                       % ('\\' + macro, len(sites), cle))
    return rap, bad


def main():
    tex_raw = read(TEX)
    if not tex_raw:
        print('[STOP] manuscript not found: %s' % TEX)
        return 2
    # Un appel dans une ligne de commentaire n'est pas un appel.
    tex = RE_COMMENT.sub('', tex_raw)

    defined, provisoire = set(), set()
    for p in SOURCES:
        body = read(p)
        if not body:
            print('[!] source missing or empty: %s' % os.path.basename(p))
            continue
        defined |= set(RE_DEF.findall(body))
        provisoire |= set(RE_PROV.findall(body))

    called = set(RE_CALL.findall(tex))
    ours = {d for d in defined if d.startswith(PREFIXES)}

    missing = sorted(called - defined)
    orphan = sorted(ours - called)
    prov_called = sorted(provisoire & called)

    print('[macros] defined %d (of which %d with prefix GE/EMP/CALIB) | called %d'
          % (len(defined), len(ours), len(called)))

    print('\n1. CALLED BUT NOT DEFINED: %d' % len(missing))
    for m in missing:
        print('   MISSING  ' + BS + m)

    print('\n2. DEFINED BUT NEVER CALLED: %d' % len(orphan))
    for m in orphan:
        print('   ORPHAN   ' + BS + m)

    print('\n3. DEFINED AS PROVISOIRE AND CALLED: %d' % len(prov_called))
    for m in prov_called:
        print('   PROVISOIRE ' + BS + m + '   (called by the .tex: the PDF prints'
              ' "??" -- REFUSED)')

    # ---- 4. \\input en double ---------------------------------
    # Une table \\input deux fois est composee deux fois, avec un
    # \\label en double : le \\ref se resout alors vers la derniere
    # occurrence. Le PDF compile sans rien dire.
    inputs = RE_INPUT.findall(tex)
    dups = sorted({p for p in inputs if inputs.count(p) > 1})
    print('\n4. FILES \\\\input MORE THAN ONCE: %d' % len(dups))
    for p in dups:
        lignes = [i for i, l in enumerate(tex_raw.split('\n'), 1)
                  if ('input{' + p + '}') in l
                  and not l.lstrip().startswith('%')]
        print('   DOUBLE   %s   (lines %s)'
              % (p, ', '.join(str(x) for x in lignes)))

    # ---- 5. LE BRAS D'AMENITE DES MACROS DU CORPS --------------
    # Le corps du papier tourne a la borne basse f_a = 0. Une macro qu'il
    # consomme doit donc sortir d'un run `amenity = none`. Ce controle
    # n'existait pas, et trois macros du Nord ont porte les valeurs du bras
    # interieur sous un nom que le corps presentait comme la borne basse --
    # au point que S5 imprimait « -0.153 at the lower bound and -0.153 at the
    # interior estimate ». Une relecture ne l'attrape pas : les deux valeurs
    # sont plausibles, et rien dans le .tex ne dit d'ou elles viennent.
    #
    # FRONTIERE : le corps s'arrete a la premiere \label{supp:...}. Elle est
    # MESUREE dans le fichier, pas supposee : si le marqueur disparait, le
    # controle s'arrete au lieu de verifier un corps vide et de passer au vert.
    prov_path = os.path.join(ROOT, 'output', 'values_provenance.tsv')
    bad = []
    n_body = 0
    if not os.path.exists(prov_path):
        print('\n5. ARM OF THE BODY MACROS: NOT CHECKED'
              ' (values_provenance.tsv missing -- re-run emit_values.py)')
        return 2
    m_sup = re.search(r'(?m)^.*' + BS + BS + r'label\{supp:', tex)
    if m_sup is None:
        print('\n5. ARM OF THE BODY MACROS: NOT CHECKED'
              ' (no ' + BS + 'label{supp:...}: body boundary'
              ' not found)')
        return 2
    corps = tex[:m_sup.start()]
    appel_corps = set(RE_CALL.findall(corps))
    arm = {}
    with io.open(prov_path, encoding='utf-8') as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            ch = line.rstrip('\n').split('\t')
            if len(ch) >= 4:
                arm[ch[0]] = (ch[3], ch[1])          # amenity, valeur
    for nm in sorted(appel_corps):
        if not nm.startswith('GE') or nm not in arm:
            continue                                  # hors perimetre GE
        amen, val = arm[nm]
        if val.strip() in (BS + 'PROVISOIRE', ''):
            continue                                  # trou assume, pas un bras
        # `invariant` : une constante de grille ou d'etiquette, qui ne descend
        # d'aucun run. La qualifier d'un bras serait faux dans les deux sens.
        if amen == 'invariant':
            continue
        # UN NOM QUI DIT SON BRAS EST AUTORISE PARTOUT.
        # `...Qfifty` designe explicitement l'estimation interieure : la
        # sous-section « Limits » du corps la cite a dessein, pour opposer les
        # deux bornes. Ce que la regle interdit, c'est une valeur interieure
        # sous un nom NU -- celui que le lecteur lit comme la borne basse.
        if nm.endswith('Qfifty'):
            continue
        n_body += 1
        # L'etiquette du manifeste pour la borne basse est « null » ; le brief
        # dit « none ». Les deux nomment le meme bras, et c'est le MANIFESTE
        # qui fait foi (regle 3), donc on accepte les deux plutot que de
        # renommer une etiquette dans l'emetteur.
        if amen not in ('none', 'null'):
            bad.append((nm, amen, val))

    print('\n5. ARM OF THE BODY MACROS: %d \\GE* macro(s) checked'
          % n_body)
    for nm, amen, val in bad:
        print('   ARM      ' + BS + nm + ' = %s  <- amenity=%s'
              ' (the body runs at f_a = 0)' % (val, amen))
    if not bad:
        print('   all from an amenity=none run')

    # ---- 6. DEUX MACROS DIFFERENTES, MEME VALEUR, MEME PHRASE --------
    # « bracketed between \GEnorthNet at the lower bound and
    #   \GEnorthNetQfifty at the interior estimate »
    # a imprime « -0.153 ... and -0.153 ... » pendant des semaines. Les deux
    # macros existaient, etaient definies, etaient appelees : les controles 1
    # a 3 les voyaient toutes vertes. Ce qui manquait, c'est qu'une phrase qui
    # OPPOSE deux quantites et en imprime une seule est fausse, et ca ne se
    # voit qu'en lisant la valeur RENDUE.
    # Le controle est un RAPPORT, pas un blocage : deux macros peuvent
    # legitimement coincider (un arrondi commun), et un blocage forcerait a
    # contourner. Il rend la coincidence visible, ce qui suffit.
    valeurs = {nm: v for nm, (a, v) in arm.items()} if arm else {}
    for p in SOURCES:
        for m in re.finditer(r'newcommand\{' + BS + BS + r'([A-Za-z]+)\}\{([^}]*)\}',
                             read(p)):
            valeurs.setdefault(m.group(1), m.group(2))
    doublons = []
    for phrase in re.split(r'(?<=[.;:])\s', corps if 'corps' in dir() else tex):
        noms = [n for n in dict.fromkeys(RE_CALL.findall(phrase))
                if n in valeurs and valeurs[n].strip()
                and valeurs[n].strip() != BS + 'PROVISOIRE']
        vus = {}
        for n in noms:
            vus.setdefault(valeurs[n], []).append(n)
        for v, ns in vus.items():
            if len(ns) > 1:
                doublons.append((v, ns, ' '.join(phrase.split())[:90]))
    print('\n6. SAME VALUE, SAME SENTENCE, TWO MACROS: %d' % len(doublons))
    for v, ns, ph in doublons:
        print('   COINCIDE %s = %s' % (', '.join(BS + n for n in ns), v))
        print('            ... %s ...' % ph)
    if not doublons:
        print('   no sentence opposes two macros of the same value')

    # ---- 7. LES TIRAGES SELECTIONNES N'ONT PAS LEUR PLACE DANS LE CORPS --
    # \GEnorthDailyRqFive / RqFifty / RqNinetyFive et leurs voisins nomment
    # TROIS TIRAGES SOLVES, choisis aux 5e/50e/95e centiles d'un classement
    # en equilibre PARTIEL. Ce classement ne se transmet que faiblement au GE
    # (13 inversions sur 29 paires, correlation 0,66), donc les lire comme des
    # percentiles de sortie est faux. Le corps cite la distribution des
    # sorties ; ces macros restent emises pour S5.6b, qui les presente comme
    # des tirages et non comme des percentiles.
    RE_SELECTED = re.compile(
        r'^GE(?:north|world)[A-Za-z]*Rq(?:Five|Fifty|NinetyFive)$')
    selected_in_body = sorted(n for n in appel_corps if RE_SELECTED.match(n))
    print('\n7. SELECTED DRAWS CALLED IN THE BODY: %d'
          % len(selected_in_body))
    for n in selected_in_body:
        print('   FORBIDDEN ' + BS + n + '  (three solved draws, not'
              ' output percentiles; the body cites the distribution)')
    if not selected_in_body:
        print('   none: the body cites only the distribution of the outputs')

    # ---- 8. LE BRAS DES TABLES, CONTRE LE BRAS DU TEXTE QUI LES ENTOURE --
    # Le controle 5 verifie les MACROS du corps. Il ne voyait pas les TABLES :
    # jusqu'au 17/09/2026, S5.1-S5.3 citait des macros a f_a = 0 en regard de
    # quatre tables emises au bras interieur. Sur la meme page, la phrase
    # disait -0.010 et la cellule -0.031. Les deux nombres etaient corrects,
    # chacun dans son bras, et aucun des sept controles ne pouvait le dire.
    #
    # LA REGLE : on decoupe le manuscrit en blocs (une \subsection = un bloc),
    # on releve dans chaque bloc les macros \GE* dont values_provenance.tsv
    # donne le bras, et si TOUTES disent `none`, alors chaque table \input
    # dans ce bloc doit porter `% amenity = none` en tete.
    #
    # POLARITE : une table \input dans un bloc a f_a = 0 et SANS etiquette
    # est signalee, pas ignoree. Un marqueur disparu ne doit pas se lire
    # comme une absence de defaut -- c'est la regle du README de tools/.
    RE_SUB = re.compile(r'(?m)^' + BS + BS + r'(?:sub)*section' + BS + BS + '*')
    cuts = [m.start() for m in RE_SUB.finditer(tex)] or [0]
    if cuts[0] != 0:
        cuts = [0] + cuts
    cuts.append(len(tex))
    blocs = [tex[cuts[i]:cuts[i + 1]] for i in range(len(cuts) - 1)]

    tab_bad, tab_nolabel, n_tab = [], [], 0
    for bloc in blocs:
        bras_bloc = set()
        for nm in set(RE_CALL.findall(bloc)):
            if not nm.startswith('GE') or nm not in arm:
                continue
            amen, val = arm[nm]
            if amen == 'invariant' or val.strip() in (BS + 'PROVISOIRE', ''):
                continue
            bras_bloc.add(amen)
        # Le bloc n'est concerne que s'il est UNANIMEMENT a f_a = 0. Un bloc
        # muet (aucune macro datee) ou mixte -- S5.4 compare les deux bras --
        # n'impose rien a ses tables.
        # Meme convention qu'au controle 5 : le MANIFESTE ecrit « null »,
        # le contrat dit « none ». Les deux nomment la borne basse.
        if not bras_bloc or bras_bloc - {'none', 'null'}:
            continue
        for rel in RE_INPUT.findall(bloc):
            if 'tab_' not in rel:
                continue
            q = os.path.join(ROOT, rel)
            body = read(q)
            if not body:
                continue
            n_tab += 1
            lab = None
            for line in body.split(chr(10))[:6]:
                t = line.strip()
                if t.startswith('%') and 'amenity' in t and '=' in t:
                    lab = t.split('=', 1)[1].strip()
                    break
            if lab is None:
                tab_nolabel.append(rel)
            elif lab not in ('none', 'null'):
                tab_bad.append((rel, lab))

    print('\n8. ARM OF THE TABLES IN A BLOCK AT f_a = 0: %d table(s) checked'
          % n_tab)
    for rel, lab in tab_bad:
        print('   ARM      %s  emitted from amenity = %s, in a block whose'
              ' macros are at f_a = 0' % (rel, lab))
    for rel in tab_nolabel:
        print('   NO TAG   %s  no `%% amenity = ...` label at the top:'
              ' the arm cannot be established' % rel)
    if not tab_bad and not tab_nolabel:
        print('   all labelled amenity = none')

    # ---- 9. LA GRILLE DE SPECIFICATIONS, DEFINIE EN UN SEUL ENDROIT ------
    # Les huit controles precedents ne cherchent aucun litteral dans le code.
    # Celui-ci ne peut faire que cela, et c'est son objet : il tient la regle
    # posee le 18/09/2026, « la grille est definie en un seul endroit ».
    #
    # POURQUOI ELLE EXISTE. Jusqu'a cette date, ct1 et ct2 portaient CHACUN
    # leur copie de la liste des formes, et ct4 comme ct6 ne comptaient les
    # specifications que sur le CSV recu. Changer la grille demandait deux
    # editions coordonnees dans deux fichiers ; rien n'aurait signale la
    # seconde oubliee. Trois producteurs portaient en plus, chacun dans son
    # langage, une copie en clair de la specification de reference.
    #
    # DEUX MOITIES :
    #   (a) les deux declarations, R et Python, portent la MEME grille ;
    #   (b) aucun producteur ne porte de liste de specifications, ni
    #       d'etiquette de specification, en dur.
    # Les commentaires sont BLANCHIS avant la recherche -- ligne a ligne, pour
    # que les numeros de ligne restent ceux du fichier. Un commentaire qui
    # explique une specification n'est pas une liste de specifications.
    grid_bad = []

    def _g(x):
        """Le %g de R sur un seuil de winsorisation : 1.0 -> "1", 0.5 -> "0.5".
        C'est ce format qui est grave dans toutes les etiquettes deja ecrites."""
        return '%g' % float(x)

    def _blanchir(text, lang):
        """Blanchit commentaires et docstrings LIGNE A LIGNE : les numeros de
        ligne rapportes restent ceux du fichier.

        Les triple-quotes sont suivies par un AUTOMATE, pas par une expression
        reguliere. Le motif de check_horizon.py part en retour arriere
        exponentiel sur un fichier de deux mille lignes portant plusieurs
        docstrings -- mesure sur 02_GE/emit_values_empirics.py, ou l'appel ne
        rendait plus la main."""
        DELIMS = ('"' * 3, "'" * 3)
        lignes, ouvert = [], None
        for line in text.split(chr(10)):
            if lang == 'py':
                reste, sortie = line, []
                while reste:
                    if ouvert is None:
                        pos = [reste.find(d) for d in DELIMS]
                        pos = [x for x in pos if x >= 0]
                        if not pos:
                            sortie.append(reste)
                            break
                        i = min(pos)
                        ouvert = reste[i:i + 3]
                        sortie.append(reste[:i])
                        reste = reste[i + 3:]
                    else:
                        j = reste.find(ouvert)
                        if j < 0:
                            break
                        ouvert = None
                        reste = reste[j + 3:]
                line = ''.join(sortie)
            cut = line.find('#')
            lignes.append(line if cut < 0 else line[:cut])
        return chr(10).join(lignes)

    P_R = os.path.join(ROOT, '01_estimation', '00_parameters.R')
    P_PY = os.path.join(ROOT, '01_estimation', '00_parameters.py')

    # (a) les deux declarations. Cote R on LIT le fichier (pas de session R a
    # ouvrir depuis Python) : les formes, les seuils et la reference y sont
    # ecrits en clair, et l'etiquette se reconstruit par la meme formule que
    # grid_spec_label. Cote Python on importe le module.
    r_src = read(P_R)
    m_forms = re.search(r'(?ms)^GRID_FORMS\s*<-\s*list\((.*?)^\)', r_src)
    m_ref = re.search(r'(?m)^GRID_REF_SPEC\s*<-\s*grid_spec_label\('
                      r'\s*"([a-z]+)"\s*,\s*(\d+)L?\s*,\s*([0-9.]+)\s*,\s*([0-9.]+)\s*\)',
                      r_src)
    r_labels, r_ref = [], None
    if m_forms:
        formes = re.findall(r'spec\s*=\s*"([a-z]+)"\s*,\s*degree\s*=\s*(\d+)L?',
                            m_forms.group(1))
        seuils = re.findall(r'lower\s*=\s*([0-9.]+)\s*,\s*upper\s*=\s*([0-9.]+)',
                            re.search(r'(?ms)^WINSOR_CONFIGS\s*<-\s*list\((.*?)^\)',
                                      r_src).group(1))
        r_labels = ['%s%d_wl%swu%s' % (sp, int(dg), _g(lo), _g(hi))
                    for sp, dg in formes for lo, hi in seuils]
    if m_ref:
        r_ref = '%s%d_wl%swu%s' % (m_ref.group(1), int(m_ref.group(2)),
                                   _g(m_ref.group(3)), _g(m_ref.group(4)))

    sys.path.insert(0, os.path.join(ROOT, '01_estimation'))
    try:
        from importlib import import_module
        _P = import_module('00_parameters')
        py_labels, py_ref = list(_P.GRID_LABELS), _P.GRID_REF_SPEC
    except Exception as exc:                      # pragma: no cover
        py_labels, py_ref = [], None
        grid_bad.append('00_parameters.py unreadable: %s' % exc)

    print('\n9. THE SPECIFICATION GRID')
    if not r_labels:
        grid_bad.append('00_parameters.R does not declare GRID_FORMS')
    elif r_labels != py_labels:
        grid_bad.append('R and Python do not carry the same grid: '
                        'R %d / Python %d; difference %s'
                        % (len(r_labels), len(py_labels),
                           ', '.join(sorted(set(r_labels) ^ set(py_labels)))))
    elif r_ref != py_ref:
        grid_bad.append('reference specification: R="%s" vs Python="%s"'
                        % (r_ref, py_ref))
    else:
        print('   the two declarations agree: %d specifications, '
              'reference %s' % (len(r_labels), r_ref))

    # (a-bis) ARMB_N : LE NOMBRE DE TIRAGES DE LA BANDE, MEME REGLE QUE LA GRILLE.
    # Ajoute le 20/09/2026. Ce nombre vivait dans une variable d'environnement
    # avec defaut (ARMB_EXTRA) et etait recopie en dur dans quatre endroits qui
    # ne disaient pas tous la meme chose ; la nuit du 19 y a laisse neuf heures.
    # Il est desormais declare dans 00_parameters.R, repete dans son jumeau
    # Python, et une divergence entre les deux est un defaut au meme titre
    # qu'une divergence de grille.
    m_armb = re.search(r'(?m)^ARMB_N\s*<-\s*(\d+)L?', r_src)
    m_armb_base = re.search(r'(?m)^ARMB_BASE_N\s*<-\s*(\d+)L?', r_src)
    try:
        py_armb, py_armb_base = int(_P.ARMB_N), int(_P.ARMB_BASE_N)
    except Exception as exc:                      # pragma: no cover
        py_armb, py_armb_base = None, None
        grid_bad.append('00_parameters.py does not declare ARMB_N: %s' % exc)
    if not m_armb or not m_armb_base:
        grid_bad.append('00_parameters.R does not declare ARMB_N / ARMB_BASE_N')
    elif py_armb is not None and (int(m_armb.group(1)) != py_armb
                                  or int(m_armb_base.group(1)) != py_armb_base):
        grid_bad.append(
            'arm B band: R declares ARMB_N=%s ARMB_BASE_N=%s, Python %s / %s. '
            'ct5 would write one grid and the emitter would expect another.'
            % (m_armb.group(1), m_armb_base.group(1), py_armb, py_armb_base))
    elif py_armb is not None and py_armb not in (py_armb_base, 2 * py_armb_base):
        grid_bad.append('arm B band: ARMB_N=%d is neither the base (%d) nor its '
                        'double; ct5 cannot build this grid.'
                        % (py_armb, py_armb_base))
    elif py_armb is not None:
        print('   arm B band: the two declarations agree, ARMB_N = %d '
              '(base %d)' % (py_armb, py_armb_base))

    # (a-ter) OMEGA ET CMIP6_MODELS (lot 38, 27/09/2026). Omega pilote l'echelle
    # du bras amenite (R) et la MDE emise (Python) ; CMIP6_MODELS est la liste
    # que le script de projection impose au repertoire (R) et que \GEcmipModels
    # imprime (Python). Deux copies libres de diverger : meme regle que la grille.
    m_om = re.search(r'(?m)^Omega\s*<-\s*([0-9.]+)', r_src)
    m_cm = re.search(r'(?ms)^CMIP6_MODELS\s*<-\s*c\((.*?)\)', r_src)
    try:
        py_om, py_cm = float(_P.OMEGA), tuple(_P.CMIP6_MODELS)
    except Exception as exc:                      # pragma: no cover
        py_om, py_cm = None, None
        grid_bad.append('00_parameters.py does not declare OMEGA / CMIP6_MODELS: %s' % exc)
    if not m_om or not m_cm:
        grid_bad.append('00_parameters.R does not declare Omega / CMIP6_MODELS')
    elif py_om is not None:
        r_om = float(m_om.group(1))
        r_cm = tuple(re.findall(r'"([^"]+)"', m_cm.group(1)))
        if r_om != py_om:
            grid_bad.append('Omega: R %s, Python %s -- the MDE and the amenity arm '
                            'no longer carry the same scale.' % (r_om, py_om))
        if r_cm != py_cm:
            grid_bad.append('CMIP6_MODELS: R %s, Python %s -- the text would no longer '
                            'name the models the projection read.'
                            % (', '.join(r_cm), ', '.join(py_cm)))
        if r_om == py_om and r_cm == py_cm:
            print('   Omega and CMIP6_MODELS: the two declarations agree '
                  '(Omega = %g, %d GCM)' % (r_om, len(r_cm)))

    # (b) aucun producteur n'en porte une en dur.
    PROD_DIRS = ('01_estimation', '02_GE', '03_inference', '04_contrast', 'tools')
    EXEMPTS = {os.path.normcase(P_R), os.path.normcase(P_PY)}
    # UNE ENTREE DE GRILLE ECRITE EN CLAIR, ET LE MOTIF DEPEND DU LANGAGE.
    #   R      : list(spec = "spline", degree = 4L) -- c'est l'idiome du paquet,
    #            et le seul par lequel une grille y a jamais ete ecrite.
    #   Python : ("spline", 4) -- le couple positionnel de GRID_FORMS.
    # ⚠️ LE MOTIF POSITIONNEL NE S'APPLIQUE PAS AU R, ET C'EST MESURE. Applique
    # aux deux, il attrapait daily_basis(clim_p, "spline", 3L, ...) dans ct1 --
    # un APPEL de fonction, pas une grille. Le remede n'a pas ete d'exempter la
    # ligne : elle portait bien une specification en dur, elle a ete retiree
    # (ct1 §A2 lit desormais GRID_FORMS). Le motif, lui, reste specifique au
    # langage, sans quoi tout appel a trois arguments serait une grille.
    PAT_ENTREE = {
        'r': re.compile(r'spec\s*=\s*["' + "'" + r'](polynomial|spline)'),
        'py': re.compile(r'\(\s*["' + "'" + r'](polynomial|spline)["' + "'"
                         + r']\s*,\s*\d'),
    }
    # UNE ETIQUETTE DE SPECIFICATION ECRITE EN CLAIR : "spline4_wl0.5wu0.5".
    # Meme motif dans les deux langages : c'est une chaine, pas une syntaxe.
    PAT_LABEL = re.compile(r'["' + "'" + r'](?:polynomial|spline)'
                           r'\d+_wl[0-9.]+wu[0-9.]+["' + "'" + r']')
    n_scan = 0
    for d in PROD_DIRS:
        dd = os.path.join(ROOT, d)
        if not os.path.isdir(dd):
            continue
        for f in sorted(os.listdir(dd)):
            if not (f.endswith('.R') or f.endswith('.r') or f.endswith('.py')):
                continue
            q = os.path.join(dd, f)
            if os.path.normcase(q) in EXEMPTS or not os.path.isfile(q):
                continue
            n_scan += 1
            src = _blanchir(read(q), 'py' if f.endswith('.py') else 'r')
            for i, line in enumerate(src.split('\n'), 1):
                for pat, quoi in ((PAT_ENTREE['py' if f.endswith('.py') else 'r'], 'grid entry'),
                                  (PAT_LABEL, 'specification label')):
                    if pat.search(line):
                        grid_bad.append('%s/%s:%d hard-coded %s: %s'
                                        % (d, f, i, quoi, line.strip()[:70]))
    print('   %d producer(s) scanned' % n_scan)

    ##################################################################
    # (c) LES ARTEFACTS PORTENT-ILS LA GRILLE ENTIERE ? (lot 10, 18/09)
    #
    # POURQUOI CE VOLET EXISTE, ET CE QUE LES DEUX AUTRES NE VOYAIENT PAS.
    # Au lot 9, la table « Robustness: productivity dose-response across
    # specifications » portait DOUZE lignes alors que la grille en declarait
    # QUINZE, et le controle 9 etait vert. Il l'etait a juste titre sur ses
    # propres termes : (a) les deux jumeaux concordaient, et (b) aucun
    # producteur ne portait de liste en dur -- 01a lisait SPEC_CONFIGS, qui
    # est bien la grille du CORPS. Le defaut n'etait ni dans les declarations
    # ni dans le code : il etait dans le RESULTAT. Un producteur peut lire la
    # bonne constante et emettre le mauvais sous-ensemble, et rien de ce qui
    # regarde le code ne l'apercevra.
    #
    # ⛔ ET LE COMPTE NE SUFFIT PAS. Quinze lignes dont une serait NS 5 au lieu
    # de Poly 4 donneraient le bon nombre et la mauvaise grille. On compare
    # donc des ENSEMBLES D'ETIQUETTES partout ou le fichier en porte, et on ne
    # se rabat sur un comptage que la ou il n'y en a pas -- la table rendue,
    # qui n'imprime aucune etiquette machine.
    #
    # POLARITE : un artefact absent ne passe PAS au vert. Un controle qui ne
    # peut pas conclure et se tait se lit comme une absence de defaut, ce qui
    # est exactement le mecanisme que ce fichier existe pour interdire.
    ##################################################################
    def _labels_of(rel, col, groupby=None):
        """{cle -> ensemble d'etiquettes} d'un CSV producteur, sans pandas :
        ce controle doit tourner sur une machine qui n'a que la bibliotheque
        standard."""
        q = os.path.join(ROOT, rel)
        if not os.path.exists(q):
            return None
        ## write.csv() de R GUILLEMETTE ses en-tetes et ses champs texte ;
        ## fwrite() de data.table ne le fait pas. Les deux formes arrivent ici
        ## -- prod_joint_tests_grid_*.csv vient de write.csv, test1_delta_north
        ## de fwrite -- donc on deshabille systematiquement. Sans cela le
        ## controle rendait « ne porte pas de colonne model_name » sur un
        ## fichier qui la porte, c'est-a-dire un faux rouge qui apprend a
        ## ignorer le controle.
        ## Le decoupage naif sur la virgule tient parce qu'aucune etiquette de
        ## specification n'en contient : elles sont de la forme
        ## <forme><degre>_wl<seuil>wu<seuil>, posee par grid_spec_label.
        _dq = lambda t: t.strip().strip('"')
        with io.open(q, encoding='utf-8', errors='replace') as fh:
            head = [_dq(h) for h in fh.readline().strip().split(',')]
            if col not in head:
                return 'nocol'
            ic = head.index(col)
            ig = head.index(groupby) if groupby and groupby in head else None
            out = {}
            for line in fh:
                f = line.rstrip('\r\n').split(',')
                if len(f) <= max(ic, ig if ig is not None else 0):
                    continue
                out.setdefault(_dq(f[ig]) if ig is not None else '',
                               set()).add(_dq(f[ic]))
        return out

    GRID_ART = [
        ('output/prod_joint_tests_grid_daily.csv',  'model_name', None),
        ('output/prod_joint_tests_grid_annual.csv', 'model_name', None),
        ('output/contribA_grid.csv',                'spec_label', None),
        ('output/contribA/test1_delta_north.csv',   'spec_label', 'arm'),
    ]
    want = set(py_labels or r_labels)
    for rel, col, grp in GRID_ART:
        got = _labels_of(rel, col, grp)
        if got is None:
            grid_bad.append('%s missing: check INCONCLUSIVE on this '
                            'artefact' % rel)
            continue
        if got == 'nocol':
            grid_bad.append('%s carries no column %s' % (rel, col))
            continue
        for key, labs in sorted(got.items()):
            if labs != want:
                grid_bad.append(
                    '%s%s: %d label(s) against %d; difference %s'
                    % (rel, (' [%s=%s]' % (grp, key)) if grp else '',
                       len(labs), len(want),
                       ', '.join(sorted(labs ^ want))))

    ## LA TABLE RENDUE, c'est-a-dire ce que le lecteur compte. Elle n'imprime
    ## aucune etiquette machine : on compte donc ses lignes de donnees, celles
    ## du corps du tabular qui portent au moins quatre separateurs.
    _tab = os.path.join(ROOT, 'output', 'tab_robustness_prod_joint.tex')
    if not os.path.exists(_tab):
        grid_bad.append('output/tab_robustness_prod_joint.tex missing: '
                        'check INCONCLUSIVE')
    else:
        _t = read(_tab)
        _body = _t.split(BS + 'midrule')[-1].split(BS + 'bottomrule')[0]
        _n = sum(1 for ln in _body.split(chr(10))
                 if ln.count('&') >= 4 and ln.rstrip().endswith(BS + BS))
        if _n != len(want):
            grid_bad.append('tab_robustness_prod_joint.tex: %d data row(s) '
                            'for %d specifications' % (_n, len(want)))
        else:
            print('   tab_robustness_prod_joint.tex: %d data rows = '
                  'the grid' % _n)

    ## EMPgridNcells doit valoir DEUX fois la grille : specifications x deux
    ## metriques d'exposition. C'est la macro que les phrases du supplement
    ## citent, donc celle qu'un lecteur verifie.
    ## ⚠️ EMPamenGridNcells n'est PAS controlee ici : le canal amenite reste sur
    ## SPEC_CONFIGS (douze), et l'y confondre reintroduirait la collision de nom
    ## que le lot 10 vient de separer.
    _ve = os.path.join(ROOT, 'output', 'values_empirics.tex')
    if not os.path.exists(_ve):
        grid_bad.append('output/values_empirics.tex missing: check '
                        'INCONCLUSIVE on EMPgridNcells')
    else:
        _src = read(_ve)
        for _nm, _exp in (('EMPgridProdNspec', len(want)),
                          ('EMPgridNcells', 2 * len(want))):
            _m = re.search(r'newcommand\{' + BS + BS + _nm + r'\}\{([^}]*)\}',
                           _src)
            if not _m:
                grid_bad.append('%s not emitted: check INCONCLUSIVE' % _nm)
            elif _m.group(1).strip() != str(_exp):
                grid_bad.append('%s = %s, expected %d'
                                % (_nm, _m.group(1).strip(), _exp))

    for b in grid_bad:
        print('   GRID     ' + b)
    if not grid_bad:
        print('   no producer' + chr(39)
              + 's code carries a hard-coded list or label;')
        print('   the artefacts and the rendered table carry the whole grid')

    # ---- 10. LA FORME D'UNE PAIRE DE BORNES SUIT SA NATURE ----------------
    print('\n10. PAIRS OF BOUNDS: INTERVAL OR RANGE')
    _pl, pair_bad = controle_paires(tex, os.path.join(ROOT, 'output',
                                                      'values_pairs.tsv'))
    for _l in _pl:
        print(_l)
    for b in pair_bad:
        print('   PAIR     ' + b)

    # ---- 11. UNE CONSTANTE IMPOSEE CITE SA SOURCE DANS LA MEME PHRASE -----
    print('\n11. IMPOSED CONSTANTS: THE VALUE COMES FROM THE MACRO, '
          'THE ATTRIBUTION MUST FOLLOW')
    _cl, cite_bad = controle_citations(tex)
    for _l in _cl:
        print(_l)
    for b in cite_bad:
        print('   CITATION ' + b)

    verdicts = []
    if missing:
        verdicts.append('%d called macro(s) undefined' % len(missing))
    if prov_called:
        verdicts.append('%d called macro(s) at PROVISOIRE: %s'
                        % (len(prov_called), ', '.join(BS + m for m in prov_called)))
    if selected_in_body:
        verdicts.append('%d selected-draw macro(s) in the body'
                        % len(selected_in_body))
    if bad:
        verdicts.append('%d body macro(s) from the interior arm'
                        % len(bad))
    if tab_bad:
        verdicts.append('%d interior-arm table(s) under a text at f_a = 0'
                        % len(tab_bad))
    if tab_nolabel:
        verdicts.append('%d table(s) without an arm label' % len(tab_nolabel))
    if grid_bad:
        verdicts.append('%d defect(s) in the specification grid'
                        % len(grid_bad))
    if pair_bad:
        verdicts.append('%d pair(s) of bounds rendered against their nature'
                        % len(pair_bad))
    if cite_bad:
        verdicts.append('%d call(s) of an imposed constant without its citation'
                        % len(cite_bad))
    print('\nVERDICT: %s' % ('; '.join(verdicts) if verdicts
                             else 'NO called macro undefined, '
                                  'no body macro outside f_a = 0'))
    return 0 if not (missing or prov_called or bad or selected_in_body or tab_bad
                     or tab_nolabel or grid_bad or pair_bad or cite_bad) else 1


if __name__ == '__main__':
    sys.exit(main())
