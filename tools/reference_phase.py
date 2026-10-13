#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""reference_phase.py -- l'invariant de la phase PAQUET.

    ⛔ LE MANUSCRIT EST CLOS. Tout lot de cette phase touche des producteurs,
       des scripts et l'arbre que les controles lisent. AUCUN n'a le droit de
       changer une valeur publiee.

Ce script PREND une reference (`--poser`) et la VERIFIE (par defaut). Il ne
modifie rien d'autre que son propre fichier de reference.

Ce qu'il epingle :
  - le md5 de chacune des macros emises (values.tex + values_empirics.tex),
    macro par macro et pas en bloc : un diff global dit « ca a bouge », un diff
    par macro dit LAQUELLE ;
  - le md5 de chaque artefact inclus par JMP_V11.tex, plus les tables generees ;
  - le nombre de pages du PDF et le compte de mots de l'abstract ;
  - LA COMPILATION (23/09/2026) : codes de retour, « Undefined control
    sequence », renvois et citations non resolus -- lus dans
    tools/COMPILATION.tsv, ecrit par tools/compiler.py, et REFUSES si le .tex
    a change depuis cette compilation ;
  - LES CONTROLES, NOMMES UN PAR UN, avec leur code de retour.

🚨 POURQUOI LES DEUX DERNIERES CLASSES (23/09/2026). Le 22/09 a 21:06, un
   edit a appele une macro qu'aucun producteur n'emet ; le PDF imprimait
   « for 1515specifications ». Cette reference est restee a 0 ecart : elle
   epinglait des VALEURS, et un appel SANS valeur n'en deplace aucune. C'est
   la regle du 22/09 sous une autre forme -- un controle qui inspecte des
   objets ne voit pas la destruction d'un objet. « Rien de publie n'a bouge »
   et « le document compile proprement » sont deux affirmations ; la
   reference porte desormais les deux.
   Et « sept controles au vert » circulait depuis le 21/09 SANS LISTE : un
   compte annonce sans enumeration ne dit ni lesquels, ni qu'il y en a un
   de plus le lendemain. Les controles sont donc NOMMES ici.

⚠️ ET IL COMPTE. Pour chaque classe -- macros, artefacts -- il epingle aussi le
   NOMBRE. C'est la regle du 22/09 : un controle qui inspecte des objets ne voit
   pas la disparition d'un objet. Une macro retiree ne se verrait pas si on ne
   comparait que les macros presentes des deux cotes.

Usage:  python tools/reference_phase.py --declarer-reprise MOTIF
                 (lot 53 : ouvre le mode rapport des emetteurs, tools/reprise.py ;
                  --poser le referme)
        python tools/reference_phase.py --poser     (une fois, en tete de phase)
        python tools/reference_phase.py             (apres chaque lot)
Sortie: rc=0 identique, rc=1 une valeur a bouge, rc=2 non concluant.
"""
import hashlib
import io
import os
import re
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TEX = os.path.join(ROOT, 'JMP_V11.tex')
REF = os.path.join(ROOT, 'tools', 'REFERENCE_PHASE.tsv')
SRC_MACROS = ('output/values.tex', 'output/values_empirics.tex')

COMPILATION = os.path.join(ROOT, 'tools', 'COMPILATION.tsv')
COMPIL_CLES = ('rc_pdflatex_passe1', 'rc_bibtex', 'rc_pdflatex_finale',
               'undefined_cs', 'erreurs_latex', 'refs_non_resolus',
               'citations_non_resolues', 'erreurs_bibtex')

# LES CONTROLES, NOMMES. Ceux de tools/README.md qui REFUSENT quelque chose
# (emit_manifest_md5 et audit_paquet n'y sont pas : ils enregistrent ou
# rapportent, ils ne refusent rien), plus check_edits_differes et
# check_labels. Ajouter un controle ici le fait APPARAITRE a la verification
# suivante : la liste est elle-meme epinglee, par son compte.
CONTROLES = ('check_macros', 'check_horizon', 'scan_litteraux',
             'known_bad_emission', 'check_refs', 'scan_octets_controle',
             'check_fraicheur', 'check_edits_differes', 'check_labels')
# Lot 86f: launched and printed, but NOT pinned. check_fraicheur compares file
# timestamps, so its code depends on the order in which a machine wrote the
# files (rc=1 in the workshop, rc=0 in a fresh clone). CONTROLES keeps all nine,
# because tools/etape_releve_pa1.py runs that list.
CONTROLES_INFORMATIFS = ('check_fraicheur',)


def compilation():
    """[(cle, valeur)] de la derniere compilation -- ou NON CONCLUANT."""
    tex_md5 = hashlib.md5(open(TEX, 'rb').read()).hexdigest()
    if not os.path.exists(COMPILATION):
        return [(k, '<INCONCLUSIVE: tools/compiler.py never run>')
                for k in COMPIL_CLES]
    d = {}
    for l in io.open(COMPILATION, encoding='utf-8'):
        if l.startswith('#') or '\t' not in l:
            continue
        k, v = l.rstrip('\n').split('\t', 1)
        d[k] = v
    if d.get('tex_md5') != tex_md5:
        return [(k, '<INCONCLUSIVE: .tex modified since the compilation of %s>'
                 % d.get('date', '?')) for k in COMPIL_CLES]
    return [(k, d.get(k, '<missing>')) for k in COMPIL_CLES]


def controles():
    """[(nom, 'rc=N')] -- chaque controle lance, son code rapporte."""
    out = []
    env = dict(os.environ, PYTHONIOENCODING='utf-8')
    for c in CONTROLES:
        r = subprocess.run([sys.executable, os.path.join(ROOT, 'tools', c + '.py')],
                           cwd=ROOT, env=env, stdout=subprocess.DEVNULL,
                           stderr=subprocess.DEVNULL)
        out.append((c, 'rc=%d' % r.returncode))
        print('   [check] %-25s rc=%d%s' % (c, r.returncode,
              '  (informative, not pinned)' if c in CONTROLES_INFORMATIFS else ''))
    return out


RE_DEF = re.compile(r'newcommand\{\\([A-Za-z]+)\}\{(.*?)\}\s*(?:%|$)', re.M)
RE_GRAPH = re.compile(r'\\includegraphics(?:\[[^\]]*\])?\{([^{}]+)\}')
RE_INPUT = re.compile(r'\\input\{([^{}]+)\}')


def flux_vivant(t):
    t = re.sub(r'(?m)^[ \t]*%.*$', '', t)
    return re.sub(r'(?<!\\)%.*', '', t)


def md5f(p):
    q = os.path.join(ROOT, p)
    if not os.path.exists(q):
        return None
    return hashlib.md5(open(q, 'rb').read()).hexdigest()


def macros():
    """{nom: valeur} -- la VALEUR, pas l'empreinte du fichier.

    Deux raisons : un commentaire de provenance qui change ne doit pas lever
    l'alarme, et une valeur qui change doit la lever en se nommant.
    """
    out = {}
    for f in SRC_MACROS:
        p = os.path.join(ROOT, f)
        if not os.path.exists(p):
            continue
        for m in RE_DEF.finditer(io.open(p, encoding='utf-8',
                                         errors='replace').read()):
            out[m.group(1)] = m.group(2).strip()
    return out


def artefacts():
    t = flux_vivant(io.open(TEX, encoding='utf-8', errors='replace').read())
    noms = set()
    for m in RE_GRAPH.finditer(t):
        n = m.group(1)
        for c in (n, n + '.pdf', n + '.png'):
            if os.path.exists(os.path.join(ROOT, 'output', c)):
                noms.add('output/' + c)
                break
    for m in RE_INPUT.finditer(t):
        n = m.group(1)
        n = n if n.endswith('.tex') else n + '.tex'
        if os.path.exists(os.path.join(ROOT, n)):
            noms.add(n)
    return sorted(noms)


def pages_et_mots():
    log = os.path.join(ROOT, 'JMP_V11.log')
    pages = ''
    if os.path.exists(log):
        m = re.findall(r'\((\d+) pages',
                       io.open(log, encoding='utf-8', errors='replace').read())
        pages = m[-1] if m else ''
    ##########################################################################
    # 🚨 ON NE FIGE PAS UNE CONVENTION, ON COMPTE LE RENDU.
    #
    # Premiere version : un depouillement du SOURCE. Deux conventions ont
    # circule le 22/09 -- retirer `\mot` avec `*` ou avec `+` -- et elles
    # rendaient 98 et 100 sur le MEME texte. Le defaut n'etait pas de choisir
    # entre elles : c'est que LES DEUX sont des approximations d'une quantite
    # qu'on peut lire directement.
    #
    # Une revue compte le resume TEL QU'IL EST COMPOSE. Dans le PDF,
    # `\GEpopGainAnnual\%` est le mot « 11.0% » et `\EMPgapColdRatioFloor{}`
    # est le mot « 2.2 ». Une convention qui retire l'antislash et ce qui suit
    # EFFACE des mots qui existent dans le texte imprime -- et l'ecart de deux
    # est exactement la.
    #
    # ⇒ ON EXTRAIT LE RESUME DE LA COUCHE TEXTE DU PDF, entre « Abstract » et
    #   « Keywords », et on compte la. Aucune convention a geler.
    #
    # ⛔ MESURE DU 22/09 : le rendu fait **100 mots**. Le budget est ATTEINT,
    #   pas depasse, et IL N'Y A AUCUNE MARGE. Le « 98 » de la premiere
    #   reference etait un artefact de depouillement ; ecrit sans qualificatif,
    #   il se serait lu dans trois semaines comme deux mots de marge.
    ##########################################################################
    mots = ''
    pdf = os.path.join(ROOT, 'JMP_V11.pdf')
    if os.path.exists(pdf):
        try:
            r = subprocess.run(['pdftotext', '-f', '1', '-l', '1', pdf, '-'],
                               capture_output=True)
            s = ' '.join(r.stdout.decode('utf-8', 'replace').split())
            i, j = s.find('Abstract'), s.find('Keywords')
            if i >= 0 and j > i:
                mots = str(len(s[i + len('Abstract'):j].split()))
        except (OSError, ValueError):
            mots = ''
    if not mots:
        mots = '<unreadable -- pdftotext missing or abstract not found>'
    return pages, mots


def etat():
    lignes = []
    # Les controles d'ABORD : si l'un d'eux mutait l'arbre sans le restaurer,
    # les macros et artefacts lus ensuite le montreraient.
    cc = [(n, v) for n, v in controles() if n not in CONTROLES_INFORMATIFS]
    lignes.append(('COMPTE', 'controles', str(len(cc))))
    for n, v in cc:
        lignes.append(('CONTROLE', n, v))
    for k, v in compilation():
        lignes.append(('COMPIL', k, v))
    mm = macros()
    lignes.append(('COMPTE', 'macros', str(len(mm))))
    for n, v in sorted(mm.items()):
        lignes.append(('MACRO', n, v))
    aa = artefacts()
    lignes.append(('COMPTE', 'artefacts', str(len(aa))))
    for a in aa:
        lignes.append(('ARTEFACT', a, md5f(a) or '<missing>'))
    p, w = pages_et_mots()
    lignes.append(('PAPIER', 'pages', p))
    lignes.append(('PAPIER', 'mots_abstract', w))
    return lignes


def declarer_reprise():
    """Lot 53 : ecrit tools/REPRISE_DELIBEREE.tsv, seule porte du mode rapport des
    emetteurs (tools/reprise.py). La declaration porte le md5 de la reference en
    vigueur : elle ne vaut que contre la reference qu'elle suspend."""
    import datetime
    import reprise
    i = sys.argv.index('--declarer-reprise')
    motif = ' '.join(sys.argv[i + 1:]).strip()
    if not motif:
        print('[STOP] --declarer-reprise requires a reason.')
        return 2
    if not os.path.exists(REF):
        print('[STOP] no reference set: there is nothing to suspend.')
        return 2
    with io.open(reprise.DECLARATION, 'w', encoding='utf-8') as fh:
        fh.write('# REPRISE DELIBEREE -- ecrite par tools/reference_phase.py --declarer-reprise.\n'
                 '# Tant que ce fichier existe et porte le md5 de la reference en vigueur,\n'
                 '# les gardes de non-regression des emetteurs RAPPORTENT au lieu de bloquer.\n'
                 '# `--poser` le supprime.\n')
        fh.write('date\t%s\n' % datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))
        fh.write('motif\t%s\n' % motif)
        fh.write('reference_md5\t%s\n' % reprise.md5f(REF))
    d = reprise.lire()
    print('[reprise] declared: %s | reference suspended %s -> %s'
          % (d['motif'], d['reference_md5'][:8], os.path.relpath(reprise.DECLARATION, ROOT)))
    return 0


def main():
    if '--declarer-reprise' in sys.argv:
        return declarer_reprise()
    poser = '--poser' in sys.argv
    cur = etat()
    if poser:
        with io.open(REF, 'w', encoding='utf-8') as fh:
            fh.write('# reference de la phase PAQUET -- posee le 22/09/2026.\n'
                     '# Toute valeur publiee doit rester identique jusqu a la\n'
                     '# soumission du paquet. Re-verifier APRES CHAQUE LOT.\n'
                     '#\n'
                     '# RE-POSEE le 22/09/2026 : mots_abstract 98 -> 100.\n'
                     '# CE N EST PAS UN CHANGEMENT DU PAPIER, c est un changement\n'
                     '# de METRIQUE : le compte porte desormais sur le RENDU et\n'
                     '# non sur un depouillement du source. Preuve que le papier\n'
                     '# n a pas bouge : au moment de la re-pose, 0 macro et 0\n'
                     '# artefact en ecart, pages identiques.\n'
                     '# ⛔ LE BUDGET EST ATTEINT A 100 MOTS, SANS AUCUNE MARGE.\n'
                     '#\n'
                     '# RE-POSEE le 23/09/2026 (passe zero, avant le lot 20) :\n'
                     '# ajout des classes COMPIL et CONTROLE. L ancienne reference\n'
                     '# epinglait des valeurs et est restee a 0 ecart pendant que\n'
                     '# le PDF imprimait une macro indefinie (edit du 22/09 21:06).\n'
                     '# Re-pose prise APRES reparation, sur une compilation propre.\n')
            for c in cur:
                fh.write('\t'.join(c) + '\n')
        print('[reference] set: %d lines -> %s'
              % (len(cur), os.path.relpath(REF, ROOT)))
        for k, n, v in cur:
            if k in ('COMPTE', 'PAPIER', 'COMPIL', 'CONTROLE'):
                print('   %-9s %-14s %s' % (k, n, v))
        ## Lot 53 : une nouvelle reference REARME les gardes des emetteurs.
        import reprise
        if os.path.exists(reprise.DECLARATION):
            os.remove(reprise.DECLARATION)
            print('[reprise] declaration removed: the non-regression guards block again.')
        return 0

    if not os.path.exists(REF):
        print('[STOP] no reference set. Run `--poser` first. INCONCLUSIVE.')
        return 2
    ref = {}
    for l in io.open(REF, encoding='utf-8'):
        if l.startswith('#') or not l.strip():
            continue
        c = l.rstrip('\n').split('\t')
        if len(c) == 3:
            ref[(c[0], c[1])] = c[2]
    cur_d = {(k, n): v for k, n, v in cur}

    bouge = [(k, n, ref[(k, n)], cur_d[(k, n)]) for (k, n) in ref
             if (k, n) in cur_d and ref[(k, n)] != cur_d[(k, n)]]
    partis = sorted((k, n) for (k, n) in ref if (k, n) not in cur_d)
    venus = sorted((k, n) for (k, n) in cur_d if (k, n) not in ref)

    print('[reference] %d entry(ies) pinned, %d read back' % (len(ref), len(cur_d)))
    print('\n1. VALUES THAT CHANGED: %d' % len(bouge))
    for k, n, a, b in sorted(bouge):
        print('   CHANGED   %-9s %-26s %r -> %r' % (k, n, a, b))
    print('\n2. GONE: %d' % len(partis))
    for k, n in partis:
        print('   GONE      %-9s %s' % (k, n))
    print('\n3. NEW: %d' % len(venus))
    for k, n in venus:
        print('   NEW       %-9s %s' % (k, n))

    n = len(bouge) + len(partis) + len(venus)
    print('\nVERDICT: %s'
          % ('%d difference(s) -- STOP and ask before going on' % n if n else
             '0 differences: nothing has changed since the reference'))
    return 1 if n else 0


if __name__ == '__main__':
    sys.exit(main())
