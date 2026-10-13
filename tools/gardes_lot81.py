"""gardes_lot81.py -- la garde de phrase du lot 81, en fonction pure.

  G10  §2.5 : « Under the annual mean the projected Northern productivity change is
       significantly positive in ‹six› cases out of ‹six›; under the daily distribution
       it is significantly positive in ‹one› case out of ‹six› ».
       - chaque compte est entre 0 et le nombre d'unites ;
       - le compte annuel est au moins le compte journalier (l'opposition de la phrase) ;
       - l'accord du nom : « cases » (pluriel) suit le compte annuel, qui ne doit donc pas
         valoir un ; « case » (singulier) suit le compte journalier, qui doit valoir un.
       La consigne du lot disait « somme <= 6 » : sur les donnees, 6 + 1 = 7, la somme de
       deux comptes sur les memes six unites n'a pas de borne a six. Remplacee par la borne
       de chaque compte.
"""


def garde_g10(n_annual, n_daily, n_units):
    v = (n_annual, n_daily, n_units)
    if not all(isinstance(x, int) for x in v) or n_units <= 0:
        return False, 'non-integer counts or zero units %r' % (v,)
    fautes = []
    if not (0 <= n_annual <= n_units and 0 <= n_daily <= n_units):
        fautes.append('count outside [0, %d]: annual %d, daily %d' % (n_units, n_annual, n_daily))
    if not n_annual >= n_daily:
        fautes.append('annual %d < daily %d' % (n_annual, n_daily))
    if n_annual == 1:
        fautes.append('plural "cases" with an annual count of 1')
    if n_daily != 1:
        fautes.append('singular "case" with a daily count of %d' % n_daily)
    return (not fautes), ('; '.join(fautes) if fautes else
                          'annual %d/%d, daily %d/%d' % (n_annual, n_units, n_daily, n_units))


MOTS = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six',
        7: 'seven', 8: 'eight', 9: 'nine', 10: 'ten'}
