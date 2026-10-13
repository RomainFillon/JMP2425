# -*- coding: utf-8 -*-
"""Controle des edits differes : etat de pose, ancre, temoin.

CE QU'IL FERME. Le 20/09/2026, la liste fermee D contenait DEUX edits deja poses
-- D4 et D7. Reposer D7 aurait casse le paragraphe : le manuscrit portait deja
son etat cible, et la substitution se serait appliquee a un texte qui n'existait
plus sous la forme attendue. Une liste d'edits differes qui ne porte pas l'etat
de pose fait poser deux fois.

TROIS CHAMPS PAR ENTREE, ET TROIS VERIFICATIONS :
  Etat    'a poser' | 'posee' | 'retiree'
  Ancre   le texte AVANT, qui doit encore exister dans le .tex
  Temoin  le texte APRES, qui ne doit pas encore exister

  a poser  -> ancre PRESENTE et temoin ABSENT
  posee    -> temoin PRESENT
  retiree  -> rien a verifier

Un edit 'a poser' dont l'ancre a disparu ARRETE : le sol a bouge sous lui, et la
substitution ne s'appliquerait pas -- silencieusement, comme toujours.

(!) LES BLANCS SONT NORMALISES. Le .tex est enroule et une ancre tombe souvent
de part et d'autre d'un saut de ligne. Sans cette normalisation le controle
rendrait 'ancre disparue' sur des ancres parfaitement presentes -- mesure du
20/09 : trois ancres sur quinze.
"""
import io
import os
import re
import sys
import unicodedata

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
## The register is a workshop file, not shipped with the package: when it is
## absent, the check is skipped (rc=0) and writes nothing. JMP_EDITS_DIFFERES_MD
## overrides the path (tools/test_edits_differes.py, case "absent").
MD = os.environ.get('JMP_EDITS_DIFFERES_MD') or os.path.join(ROOT, 'tools', 'EDITS_DIFFERES.md')
TEX = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, 'JMP_V11.tex')

RE_ENTREE = re.compile(r'^## (D\d+)\b(.*)$')
RE_CHAMP = re.compile(r'^- \*\*(Etat|État|Ancre|Témoin|Temoin)\*\*\s*:\s*(.*)$')


def flux_vivant(path):
    """Le .tex sans commentaires ni blocs comment, blancs normalises."""
    B = chr(92)
    E = re.escape(B)
    src = io.open(path, encoding='utf-8', errors='replace').read().split('\n')
    out, dans = [], False
    for l in src:
        if re.search(E + r'begin' + re.escape('{comment}'), l):
            dans = True
        if dans:
            if re.search(E + r'end' + re.escape('{comment}'), l):
                dans = False
            continue
        if l.lstrip().startswith('%'):
            continue
        out.append(re.sub(r'(?<!' + E + r')%.*$', '', l))
    return re.sub(r'\s+', ' ', ' '.join(out))


def lire_entrees(path):
    entrees, cur = [], None
    for l in io.open(path, encoding='utf-8').read().split('\n'):
        m = RE_ENTREE.match(l)
        if m:
            cur = {'id': m.group(1), 'titre': m.group(2).strip(' —-'),
                   'etat': None, 'ancre': None, 'temoin': None}
            entrees.append(cur)
            continue
        if cur is None:
            continue
        c = RE_CHAMP.match(l)
        if not c:
            continue
        cle = c.group(1).replace('É', 'E').replace('é', 'e').lower()
        val = c.group(2).strip()
        b = re.search(r'`(.+?)`', val)
        if cle == 'etat':
            # ⛔ UN JETON MACHINE, JAMAIS DE LA PROSE. Le champ portait un
            # libelle francais et le controle le normalisait : `'a poser' not
            # in 'a poser'` accentue est VRAI, et dix entrees sur treize sont
            # parties dans la mauvaise branche a cause d une seule lettre.
            # La cause n etait pas l accent, c etait qu un meme champ servait
            # a la machine et au lecteur. Le jeton est en ASCII entre
            # accents graves ; le libelle francais vit a cote, et le controle
            # ne le lit pas.
            cur['etat'] = b.group(1).strip().upper() if b else None
            if cur['etat'] not in ('TODO', 'DONE', 'DROPPED'):
                cur['etat_bad'] = val
        elif cle == 'ancre':
            cur['ancre'] = b.group(1) if b else None
        else:
            cur['temoin'] = b.group(1) if b else None
    return entrees


def controle(entrees, tex):
    bad, rap = [], []
    for e in entrees:
        i, et = e['id'], (e['etat'] or '')
        if e.get('etat_bad') is not None or not et:
            bad.append('%s: state %r is not a TODO / DONE / DROPPED token'
                       % (i, e.get('etat_bad', et)))
            continue
        if et == 'DROPPED':
            rap.append('   %-4s dropped' % i)
            continue
        if not e['ancre'] and not e['temoin']:
            bad.append('%s: neither anchor nor witness declared' % i)
            continue
        # ######################################################################
        # UNICITE DE L'ANCRE -- posee le 22/09/2026, a l'occasion de D3.
        #
        # Le controle ne verifiait que la PRESENCE. Une ancre qui correspond a
        # plusieurs endroits ne designe aucun site : l'edit s'appliquerait au
        # premier venu. Le controle passait au vert dessus.
        #
        # ⛔ CE CONTROLE N'AURAIT PAS ATTRAPE D3, ET IL FAUT LE DIRE. Mesure du
        # 22/09 : l'ancre declaree par D3 -- `the reference specification,
        # $\delta^*=\EMPosterDelta$` -- correspondait UNE SEULE FOIS dans le
        # flux vivant, l.1133, dans le supplement. Elle etait donc unique et
        # MAL PLACEE : le texte de remplacement visait le corps, l.445-448, ou
        # la meme idee s'ecrit dans l'ordre inverse (« gives ... for the
        # reference specification ») et ne correspond pas a l'ancre. Ce que ce
        # controle ferme est une famille VOISINE -- l'ancre ambigue -- pas
        # celle de D3, qui demanderait de verifier que le site de l'ancre est
        # coherent avec le texte de remplacement. A verser a la liste.
        # ######################################################################
        a_n = tex.count(e['ancre']) if e['ancre'] else 0
        t_n = tex.count(e['temoin']) if e['temoin'] else 0
        a_pres = a_n > 0
        t_pres = t_n > 0
        if et == 'DONE':                                    # etat : posee
            if not t_pres:
                bad.append('%s declared DONE but its witness is missing from the '
                           '.tex: %r' % (i, e['temoin']))
            else:
                # Le compte du temoin est RAPPORTE, pas bloquant : un temoin
                # multiple signale une pose repetee ou un fragment trop court,
                # mais il ne peut pas faire appliquer l edit au mauvais endroit
                # -- l edit est deja pose. Seule l ANCRE bloque.
                rap.append('   %-4s done, witness present%s'
                           % (i, '' if t_n == 1 else ' (%d occurrences)' % t_n))
        else:                                               # « a poser »
            if t_pres:
                bad.append('%s declared TODO but its witness is ALREADY in '
                           'the .tex: %r -- stale state, the edit is applied'
                           % (i, e['temoin']))
            elif not a_pres:
                bad.append('%s declared TODO but its anchor has DISAPPEARED from the '
                           '.tex: %r -- the ground has moved, the substitution would '
                           'not apply' % (i, e['ancre']))
            elif a_n > 1:
                bad.append('%s: its anchor matches %d TIMES in the live '
                           'text: %r -- an anchor present several times '
                           'designates no site. Lengthen the anchor until it is '
                           'unique.' % (i, a_n, e['ancre']))
            else:
                rap.append('   %-4s todo, unique anchor, witness missing' % i)
    return rap, bad


RE_LIGNE_TAB = re.compile(r'^\|\s*(D\d+)\s*\|(.*?)\|(.*?)\|\s*$', re.M)


def classe_entree(e):
    """closed / dropped / open, depuis le champ Etat d'une entree."""
    et = (e['etat'] or '').upper()
    if 'DROPPED' in et:
        return 'dropped'
    return 'closed' if 'DONE' in et else 'open'


def classe_ligne(cellule):
    """closed / dropped / open, depuis la cellule `etat` du tableau de tete.

    Les trois marqueurs sont EXCLUSIFS : une cellule qui n'en porte aucun, ou
    qui en porte deux, ne peut pas etre classee -- et un controle qui ne peut
    pas conclure NE PASSE PAS AU VERT.
    """
    vus = [n for n, m in (('closed', '✅'), ('open', '⏸'),
                          ('dropped', '❌')) if m in cellule]
    return vus[0] if len(vus) == 1 else None


def controle_tableau(entrees, md_txt):
    """################################################################
    LE TABLEAU DE TETE EST UN RESUME, ET IL DOIT DIRE LA MEME CHOSE QUE LES
    ENTREES. Pose le 22/09/2026.

    LE DEFAUT QU'IL FERME, ET IL S'EST PRODUIT. Ce controle lisait les CHAMPS
    des entrees et jamais l'en-tete. Releve du 22/09 : quatre lignes mentaient
    -- D9, D10 et D13 annoncees « a poser » avec un bloqueur alors que leurs
    entrees portaient DONE depuis la veille, D5 avec un bloqueur inexistant --
    et DIX entrees, D17 a D26, n'avaient aucune ligne. Le controle etait au
    VERT pendant ce temps. Consequence mesuree le matin meme : D9 a ete reprise
    comme bloquee, et il a fallu aller en source pour etablir qu'elle ne l'etait
    pas.

    ⇒ `EDITS_DIFFERES.md` est la memoire qui survit a un `/clear`. Si son
      en-tete ment, le prochain qui le lit repart faux. C'est la regle « un etat
      lu a une seule source n'est pas une mesure », appliquee au fichier qui
      porte cette regle.
    ################################################################"""
    bad = []
    lignes = {m.group(1): m.group(2) for m in RE_LIGNE_TAB.finditer(md_txt)}
    par_id = {e['id']: e for e in entrees}

    for i, e in par_id.items():
        if i not in lignes:
            bad.append('%s has an entry but NO row in the header '
                       'table -- the summary announces fewer edits than there are'
                       % i)
    for i in lignes:
        if i not in par_id:
            bad.append('%s has a row in the header table but NO '
                       'entry -- the summary announces an edit that does not exist' % i)
    for i, e in sorted(par_id.items(), key=lambda kv: int(kv[0][1:])):
        if i not in lignes:
            continue
        ce, cl = classe_entree(e), classe_ligne(lignes[i])
        if cl is None:
            bad.append('%s: the state cell of the table carries neither exactly one '
                       '%s, nor exactly one %s, nor exactly one %s -- it cannot be '
                       'classified, and a check that cannot conclude does '
                       'not pass. Cell: %r'
                       % (i, '✅', '⏸', '❌', lignes[i].strip()))
        elif ce != cl:
            bad.append('%s: the header table says %r and the entry says %r '
                       '(State: %r). The summary is wrong, and it is what is read '
                       'first.' % (i, cl, ce, (e['etat'] or '').strip()))
    return bad


def main():
    if not os.path.exists(MD):
        print('[edits] %s: register not shipped; check skipped' % os.path.basename(MD))
        return 0
    tex = flux_vivant(TEX)
    entrees = lire_entrees(MD)
    rap, bad = controle(entrees, tex)
    bad_tab = controle_tableau(
        entrees, io.open(MD, encoding='utf-8', errors='replace').read())
    print('[edits] %s: %d entry(ies)' % (os.path.basename(MD), len(entrees)))
    for l in rap:
        print(l)
    for b in bad:
        print('   EDIT     ' + b)
    for b in bad_tab:
        print('   TABLE    ' + b)
    if not bad_tab:
        print('   header table: %d row(s), all consistent'
              % len(entrees))
    n = len(bad) + len(bad_tab)
    print('\nVERDICT: %s' % ('%d defect(s)' % n if n
                             else 'states, anchors, witnesses and header '
                                  'table agree'))
    return 1 if n else 0


if __name__ == '__main__':
    sys.exit(main())
