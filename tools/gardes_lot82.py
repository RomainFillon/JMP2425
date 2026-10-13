r"""gardes_lot82.py -- la garde d'egalite des effectifs (lot 82), en fonction pure.

  §2.3 et B.1 : « Before the exclusions the sample has \EMPtradeNobs observations;
  \EMPtradeNoClimObs of them belong to pairs involving one of \EMPtradeNoClimCountries
  countries without an ERA5 record and are dropped, and the sample on which the
  dose-responses are estimated has \EMPtradeNobsStageThree. »

  Trois artefacts, deux calculs : le nombre de lignes de BACI_raw_pairs_annual.rds doit etre
  l'effectif de la premiere etape (stage1_trade_diagnostics.csv), et cet effectif moins les
  lignes sans climat doit etre l'effectif d'estimation (prod_joint_tests_*.csv). Coherence
  entre deux calculs : BLOQUANTE meme en reprise.
"""


def garde_effectifs(n_baci, n_stage1, n_sans_climat, n_estimation, n_pays):
    v = (n_baci, n_stage1, n_sans_climat, n_estimation, n_pays)
    if not all(isinstance(x, int) and x >= 0 for x in v):
        return False, 'non-integer or negative counts %r' % (v,)
    fautes = []
    if n_baci != n_stage1:
        fautes.append('BACI %d rows != first stage %d' % (n_baci, n_stage1))
    if n_stage1 - n_sans_climat != n_estimation:
        fautes.append('%d - %d = %d != estimation %d'
                      % (n_stage1, n_sans_climat, n_stage1 - n_sans_climat, n_estimation))
    if (n_sans_climat == 0) != (n_pays == 0):
        fautes.append('%d rows without climate for %d countries' % (n_sans_climat, n_pays))
    return (not fautes), ('; '.join(fautes) if fautes else
                          '%d - %d (%d countries without ERA5) = %d'
                          % (n_stage1, n_sans_climat, n_pays, n_estimation))
