"""gardes_lot78.py -- les gardes de phrase du lot 78 (biais de binning), en
fonctions pures. Chacune rend (ok, detail) ; emit_values_empirics.py les appelle
sur les sorties de ct8/ct9 et tools/test_gardes_lot78.py les exerce contre leur
connu-faux.

Chaque garde porte la phrase qu'elle tient. Une garde qui refuse sur les donnees
reelles dit que la phrase est FAUSSE : elle se reecrit, la garde ne se relache pas.

  G2  « the design does not use them » / « draws on neither source » (B1, B2, B3)
  G3  « positive and significant », « reverses the contrast », « turns negative
      under both arms » (B2, B3)
  G4  « includes zero », « the cold arm falls to » (B2, B3)
  G5  « ‹min›--‹max›% of the variance ... is within-pair », « correlates at » :
      ingredients petits (B2, B3)
  G6  supplement E : la part intra-cellule est plus faible en productivite qu'en
      amenite, et celle de la temperature moyenne est petite (B4)
  G7  B3 : « the population-weighted warming trend is faster in cold countries,
      the trend of each bin count correlates positively with baseline temperature
      in both tails and not in the centre »
  B6  §2.2 : la forme ajoutee « carries no Akaike weight in either arm »
"""
import math


def _fin(*xs):
    return all(x is not None and isinstance(x, (int, float)) and math.isfinite(x) for x in xs)


def garde_g2(shift_d, shift_cold, seuil=0.5):
    """La variante between deplace D et le bras froid de moins de `seuil` ecart-type."""
    if not _fin(shift_d, shift_cold):
        return False, 'non-finite value (%r, %r)' % (shift_d, shift_cold)
    ok = abs(shift_d) < seuil and abs(shift_cold) < seuil
    return ok, 'D shift %+.3f se, cold arm %+.3f se (threshold %.2f)' % (shift_d, shift_cold, seuil)


def garde_g3(cold_pp_lo, d_pp_lo, d_pp_hi, d_base, north_pp_daily, north_pp_annual,
             north_base_daily, north_base_annual):
    """pp : bras froid positif et significatif ; D exclut zero, signe oppose a base ;
    Delta_north negatif dans les deux bras, positif dans les deux sous base."""
    if not _fin(cold_pp_lo, d_pp_lo, d_pp_hi, d_base, north_pp_daily, north_pp_annual,
                north_base_daily, north_base_annual):
        return False, 'non-finite value'
    fautes = []
    if not cold_pp_lo > 0:
        fautes.append('lower bound of the pp cold arm %+.4f <= 0' % cold_pp_lo)
    excl = d_pp_lo * d_pp_hi > 0
    if not excl:
        fautes.append('CI of D pp [%+.4f, %+.4f] contains zero' % (d_pp_lo, d_pp_hi))
    elif not (math.copysign(1, d_pp_lo) != math.copysign(1, d_base)):
        fautes.append('D pp and D base have the same sign (%+.4f, %+.4f)' % (d_pp_lo, d_base))
    if not (north_pp_daily < 0 and north_pp_annual < 0):
        fautes.append('Delta_north pp not negative in both arms (%+.4f, %+.4f)'
                      % (north_pp_daily, north_pp_annual))
    if not (north_base_daily > 0 and north_base_annual > 0):
        fautes.append('Delta_north base not positive in both arms (%+.4f, %+.4f): '
                      '"turns negative" requires a change of sign'
                      % (north_base_daily, north_base_annual))
    return (not fautes), ('; '.join(fautes) if fautes else
                          'pp cold arm > 0 significant; D pp excludes 0, opposite sign; '
                          'Delta_north pp < 0 in both arms')


def garde_g4(d_pp5_lo, d_pp5_hi, cold_pp5, cold_pp):
    """pp5 : IC de D contient zero ; |bras froid pp5| < |bras froid pp|."""
    if not _fin(d_pp5_lo, d_pp5_hi, cold_pp5, cold_pp):
        return False, 'non-finite value'
    fautes = []
    if not (d_pp5_lo <= 0 <= d_pp5_hi):
        fautes.append('CI of D pp5 [%+.4f, %+.4f] excludes zero' % (d_pp5_lo, d_pp5_hi))
    if not abs(cold_pp5) < abs(cold_pp):
        fautes.append('|pp5 cold arm| %.4f >= |pp cold arm| %.4f' % (abs(cold_pp5), abs(cold_pp)))
    return (not fautes), ('; '.join(fautes) if fautes else
                          'CI of D pp5 contains 0; cold arm %.3f -> %.3f' % (cold_pp, cold_pp5))


def garde_g5(part_intra_max_pct, corr_y, seuil_pct=5.0, seuil_corr=0.3):
    """Ingredients petits : part intra-paire max < 5 % ; correlation de tendance < 0,3."""
    if not _fin(part_intra_max_pct, corr_y):
        return False, 'non-finite value'
    ok = part_intra_max_pct < seuil_pct and corr_y < seuil_corr
    return ok, 'max within share %.2f %% (threshold %g); y trend corr %.3f (threshold %g)' % (
        part_intra_max_pct, seuil_pct, corr_y, seuil_corr)


def garde_g6(intra_prod_pct, intra_amen_pct, tmean_prod_pct, tmean_amen_pct, seuil_pct=5.0):
    """CRH : intra prod < intra amen ; les deux parts de la temperature moyenne < 5 %."""
    if not _fin(intra_prod_pct, intra_amen_pct, tmean_prod_pct, tmean_amen_pct):
        return False, 'non-finite value'
    fautes = []
    if not intra_prod_pct < intra_amen_pct:
        fautes.append('within prod %.2f %% >= within amen %.2f %%' % (intra_prod_pct, intra_amen_pct))
    if not (tmean_prod_pct < seuil_pct and tmean_amen_pct < seuil_pct):
        fautes.append('mean temperature %.2f / %.2f %% (threshold %g)'
                      % (tmean_prod_pct, tmean_amen_pct, seuil_pct))
    return (not fautes), ('; '.join(fautes) if fautes else
                          'within prod %.1f < amen %.1f %%; tmean %.2f / %.2f %%'
                          % (intra_prod_pct, intra_amen_pct, tmean_prod_pct, tmean_amen_pct))


def garde_g7(cor_trend_popw, bins):
    """bins : liste de (borne_basse, correlation). Queues = borne basse < 0 ou >= 25 :
    correlation > 0 partout ; centre = [0, 25) : aucune correlation > 0. Et la
    tendance ponderee est plus rapide dans les pays froids (correlation < 0)."""
    if not _fin(cor_trend_popw) or not bins or not all(_fin(lo, c) for lo, c in bins
                                                       if lo != float('-inf')):
        return False, 'non-finite value or missing bins'
    queues = [(lo, c) for lo, c in bins if lo < 0 or lo >= 25]
    centre = [(lo, c) for lo, c in bins if 0 <= lo < 25]
    fautes = []
    if not cor_trend_popw < 0:
        fautes.append('weighted cor(trend, Tbase) %+.3f >= 0' % cor_trend_popw)
    if not queues or not centre:
        fautes.append('empty tails or centre')
    mauvaises_q = [(lo, c) for lo, c in queues if not (_fin(c) and c > 0)]
    mauvais_c = [(lo, c) for lo, c in centre if not (_fin(c) and c <= 0)]
    if mauvaises_q:
        fautes.append('tail(s) without positive correlation: %s' % mauvaises_q)
    if mauvais_c:
        fautes.append('centre with positive correlation: %s' % mauvais_c)
    return (not fautes), ('; '.join(fautes) if fautes else
                          'weighted trend %+.3f; %d tail bins > 0, %d centre bins <= 0'
                          % (cor_trend_popw, len(queues), len(centre)))


def garde_b6(poids_par_bras, seuil=1e-6, bras=('daily', 'annual')):
    """poids_par_bras : {bras: [poids d'Akaike des specifications de la forme ajoutee]}.
    Max < seuil dans chacun des deux bras ; un bras absent ou vide refuse."""
    fautes = []
    for b in bras:
        w = poids_par_bras.get(b)
        if not w:
            fautes.append('arm %s: no specification of the added form' % b)
            continue
        if not all(_fin(x) for x in w):
            fautes.append('arm %s: non-finite weight' % b)
            continue
        if max(w) >= seuil:
            fautes.append('arm %s: maximal weight %.3e >= %.0e' % (b, max(w), seuil))
    return (not fautes), ('; '.join(fautes) if fautes else
                          'maximal weight of the added form < %.0e in both arms' % seuil)
