"""reprise.py -- la REPRISE DELIBEREE, seule porte du mode rapport (lot 53, 28/09/2026).

Les gardes de non-regression des emetteurs (la valeur attendue d'une macro est le
resultat precedent) ont raison quand rien ne doit changer, c'est-a-dire sous la
regle de phase. Pendant une re-estimation deliberee, elles doivent RAPPORTER
l'ecart au lieu de bloquer : la nouvelle valeur est emise et l'ecart est ecrit.

Le mode n'est PAS un drapeau : aucun argument de ligne de commande, aucune variable
d'environnement ne l'ouvre. Il n'existe que si tools/REPRISE_DELIBEREE.tsv a ete
ecrit par `python tools/reference_phase.py --declarer-reprise "motif"`, et que ce
fichier porte le md5 de la reference de phase EN VIGUEUR. `--poser` (la nouvelle
reference, apres la course) supprime la declaration : les gardes se rearment.

    from reprise import reprise_declaree
    R = reprise_declaree()      # None hors reprise ; sinon dict(date, motif, ...)
"""
import hashlib
import io
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DECLARATION = os.path.join(ROOT, 'tools', 'REPRISE_DELIBEREE.tsv')
REFERENCE = os.path.join(ROOT, 'tools', 'REFERENCE_PHASE.tsv')
CLES = ('date', 'motif', 'reference_md5')


def md5f(p):
    return hashlib.md5(open(p, 'rb').read()).hexdigest()


def lire(chemin=DECLARATION, reference=REFERENCE):
    """dict de la declaration, None si absente ; refuse une declaration invalide."""
    if not os.path.exists(chemin):
        return None
    d = {}
    for l in io.open(chemin, encoding='utf-8'):
        if l.startswith('#') or '\t' not in l:
            continue
        k, v = l.rstrip('\n').split('\t', 1)
        d[k] = v
    manque = [k for k in CLES if not d.get(k)]
    if manque:
        raise SystemExit('[STOP] %s incomplete: %s missing.' % (chemin, ', '.join(manque)))
    if not os.path.exists(reference) or md5f(reference) != d['reference_md5']:
        raise SystemExit('[STOP] %s is STALE: it was declared against reference '
                         '%s, the reference in force is %s. A resumption is declared against '
                         'the reference it suspends.' % (
                             chemin, d['reference_md5'][:8],
                             md5f(reference)[:8] if os.path.exists(reference) else '<missing>'))
    return d


def reprise_declaree():
    return lire()


## LOT 85o : sous le bac a sable de known_bad_emission.py (EMIT_SANDBOX_DIR), le
## journal des gardes va dans la copie, jamais dans output/.
GARDES = os.path.join(os.environ.get('EMIT_SANDBOX_DIR') or os.path.join(ROOT, 'output'),
                      'ecarts_reprise_gardes.tsv')


def rapporte(qui, garde, attendu, obtenu):
    """Une garde gelee qui, en reprise, rapporte au lieu de bloquer : une ligne datee
    dans output/ecarts_reprise_gardes.tsv, et la meme ligne dans le journal."""
    import datetime
    neuf = not os.path.exists(GARDES)
    with io.open(GARDES, 'a', encoding='utf-8') as fh:
        if neuf:
            fh.write('date\tscript\tgarde\tattendu\tobtenu\n')
        fh.write('%s\t%s\t%s\t%s\t%s\n' % (datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S'),
                                           qui, garde, attendu, obtenu))
    print('[REPORT MODE] %s: guard "%s" -- expected %s, obtained %s (reported, non-blocking)'
          % (qui, garde, attendu, obtenu))


def bandeau(d, qui):
    """La ligne que chaque script en mode rapport ecrit dans le journal de la course."""
    return ('[REPORT MODE -- DECLARED RESUMPTION] %s: declared on %s, reason: %s '
            '(reference suspended %s). Differences from the expected values are REPORTED, '
            'not blocking.' % (qui, d['date'], d['motif'], d['reference_md5'][:8]))
