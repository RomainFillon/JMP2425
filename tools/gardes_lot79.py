"""gardes_lot79.py -- la garde de phrase du lot 79 (supplement E), en fonctions pures.

  G8  S4 : « The design has no cell fixed effect: identification is across cells within
      a sub-region and period ». Deux moitiés, et les deux doivent tenir :

  (a) LE CODE. Chaque estimation CRH du paquet (01d_crh_triangulation.R,
      03_inference/step*.R) ecrit ses effets fixes « | fe_block + fe », ou fe_block est
      passe au site d'appel ("grid_2c" pour l'amenite, "land_type" pour la productivite)
      et fe = sous-unite x periode, (year - 1985) / 5 * max(subnal) + subnal, subnal
      etant "subnal_ctry" ou "subnal_reg". Toute autre partie effets fixes, tout autre
      fe_block, toute autre definition de fe refuse : la phrase ne serait plus verifiee.

  (b) LES DONNEES. Si un effet fixe de cellule (ou plus fin) etait absorbe, le residu de
      chaque regresseur aurait une moyenne nulle dans chaque cellule, et sa part
      intra-cellule vaudrait EXACTEMENT 1. Les parts mesurees par ct9
      (crh_within_cell.csv, intra_resid_max) doivent donc etre nettement sous 1.
"""
import re

FE_BLOCKS_PERMIS = {'grid_2c', 'land_type'}
SUBNAL_PERMIS = {'subnal_ctry', 'subnal_reg'}
## Deux formes permises : la regression d'un bras, et la regression EMPILEE des deux bras
## (step3_vcov.R l.97, step7 l.146 : la covariance entre bras), ou les memes effets fixes
## sont interagis avec l'indicateur de bras -- toujours aucun effet fixe de cellule.
FE_PARTS_PERMISES = {'fe_block+fe', 'fe_block^arm+fe^arm'}
_FE_DEF = re.compile(r'fe\s*:=\s*\(\s*year\s*-\s*1985\s*\)\s*/\s*5\s*\*\s*max\(\s*subnal\s*\)\s*\+\s*subnal\b')
_EMPILE = re.compile(r'paste0\(\s*fe_block\s*,\s*"\^arm"\s*\)')
_FE_PART = re.compile(r'"\|"\s*,\s*([^)]*?)\)')
## le fe_block est l'argument qui PRECEDE le subnal au site d'appel ; une liste de
## colonnes ("subnal_reg", "subnal_ctry") n'en est pas un.
_APPEL = re.compile(r'"(?!subnal_)([A-Za-z0-9_.^]+)"\s*,\s*"(subnal_[A-Za-z0-9_]+)"')


def garde_g9(inc_sub, inc_arc, inc_mid, share_sub, share_arc, seuil_part=90.0):
    """G9 (lot 79, R7 ; re-ancree au lot 85o) -- §3.4 : « What the daily distribution adds
    is an income loss that grows with latitude. Net of the amenity term, ..., the daily
    extra is ‹Mid› in the temperate band, ‹Sub› between 45 and 60N and ‹Arc› above 60N
    ... Nearly all of it is produced by warming in the North itself ».
      - une perte qui grandit avec la latitude : inc_arc < inc_sub < inc_mid < 0 ;
      - « nearly all » : les deux parts « warming confined to the North » > 90 %.
    LOT 85o : la condition « close to zero » (|inc_mid| < 0,0005) est remplacee par
    inc_mid < 0 et inc_sub < inc_mid ; la bande temperee est desormais une macro
    (\\GEchanIncomeMid) dans le texte et dans tab:ge_decomp_lat."""
    v = (inc_sub, inc_arc, inc_mid, share_sub, share_arc)
    if not all(isinstance(x, float) and x == x for x in v):
        return False, 'non-finite value %r' % (v,)
    fautes = []
    if not inc_sub < 0:
        fautes.append('income 45-60N %+.4f >= 0: not a loss' % inc_sub)
    if not inc_arc < inc_sub:
        fautes.append('income >60N %+.4f not more negative than 45-60N %+.4f' % (inc_arc, inc_sub))
    if not inc_mid < 0:
        fautes.append('temperate band %+.4f >= 0: not a loss' % inc_mid)
    if not inc_sub < inc_mid:
        fautes.append('income 45-60N %+.4f not more negative than the temperate band %+.4f'
                      % (inc_sub, inc_mid))
    if not (share_sub > seuil_part and share_arc > seuil_part):
        fautes.append('North-only shares %.1f / %.1f %% (threshold %g)' % (share_sub, share_arc, seuil_part))
    return (not fautes), ('; '.join(fautes) if fautes else
                          'income %+.4f (45-60N) > %+.4f (>60N); temperate %+.4f; shares %.1f / %.1f %%'
                          % (inc_sub, inc_arc, inc_mid, share_sub, share_arc))


def _sans_commentaires(txt):
    return '\n'.join(l.split('#', 1)[0] for l in txt.split('\n'))


def analyse_code(textes):
    """textes : {nom: source R}. Rend (parties effets fixes, fe_blocks, subnals,
    fichiers sans definition de fe conforme parmi ceux qui ont des parties effets fixes)."""
    parts, blocks, subnals, sans_def = [], set(), set(), []
    for nom, t in textes.items():
        s = _EMPILE.sub('fe_block^arm', re.sub(r'\s+', ' ', _sans_commentaires(t)))
        p = [re.sub(r'[\s"]', '', m.group(1)).replace(',', '') for m in _FE_PART.finditer(s)]
        parts += [(nom, x) for x in p]
        for b, sn in _APPEL.findall(s):
            blocks.add(b); subnals.add(sn)
        if p and not _FE_DEF.search(s):
            sans_def.append(nom)
    return parts, blocks, subnals, sans_def


def garde_g8(textes, intra_resid_max, seuil=0.99):
    """textes : {nom: source R} des estimateurs CRH ; intra_resid_max : liste des parts
    intra-cellule maximales par canal x bras (crh_within_cell.csv)."""
    fautes = []
    parts, blocks, subnals, sans_def = analyse_code(textes)
    if not parts:
        fautes.append('no fixed-effects part found in the CRH code')
    autres = sorted({x for _, x in parts if x not in FE_PARTS_PERMISES})
    if autres:
        fautes.append('fixed-effects part(s) outside the allowed forms %s: %s'
                      % (sorted(FE_PARTS_PERMISES), autres))
    if not blocks or blocks - FE_BLOCKS_PERMIS:
        fautes.append('fe_block passed: %s (allowed %s)' % (sorted(blocks), sorted(FE_BLOCKS_PERMIS)))
    if not subnals or subnals - SUBNAL_PERMIS:
        fautes.append('subnal passed: %s (allowed %s)' % (sorted(subnals), sorted(SUBNAL_PERMIS)))
    if sans_def:
        fautes.append('definition of fe (sub-unit x period) missing: %s' % sans_def)
    if not intra_resid_max or not all(isinstance(v, float) and v == v for v in intra_resid_max):
        fautes.append('within-cell shares missing or non-finite')
    elif max(intra_resid_max) >= seuil:
        fautes.append('maximal within-cell share %.4f >= %.2f: a cell fixed effect is '
                      'absorbed' % (max(intra_resid_max), seuil))
    return (not fautes), ('; '.join(fautes) if fautes else
                          '%d formula(s) in the allowed forms, fe_block %s, subnal %s; max '
                          'within-cell share %.3f < %.2f'
                          % (len(parts), sorted(blocks), sorted(subnals), max(intra_resid_max), seuil))
