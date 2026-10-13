"""temoins_poids.py -- les deux normalisations des poids de population, chacune gardee
par son propre temoin (lot 46, 27/09/2026).

LE DEFAUT QU'ILS FERMENT. Une meme colonne, `share_population`, designait deux objets :
normalisee PAR PAYS dans lengths_country.csv, PAR CELLULE dans population_weight_all.csv.
Les consommateurs l'ont lue a l'envers (argmax d'une part par pays pour designer le pays
dominant d'une cellule). Desormais chaque part porte son denominateur dans son nom, et chaque
temoin NOMME dans son message celle qu'il verifie :

  part_de_la_cellule : population du pays dans la cellule / population de la cellule
                       -> somme 1 PAR CELLULE (et par periode)
  part_du_pays       : population du pays dans la cellule / population du pays
                       -> somme 1 PAR PAYS (et par periode)

Un message ambigu sur une garde de normalisation ne vaudrait pas mieux que le nom partage.
"""
TOL = 1e-9


def temoin_part_non_triviale(serie, nom, tol=1e-12):
    """Refuse une part qui vaut 1 PARTOUT (lot 51) : la cle de groupement contient la
    dimension sur laquelle la part se normalise. Portee : ce defaut-la seulement (les
    comptes de polygones et la population entiere comptee par pays ne valaient pas 1)."""
    x = serie[serie.notna()]
    if len(x) and ((x - 1).abs() <= tol).all():
        raise ValueError('[STOP] part non triviale (%s) : les %d valeurs valent 1. La cle de '
                         'groupement contient la dimension sur laquelle la part se normalise.'
                         % (nom, len(x)))
    return int(((x - 1).abs() > tol).sum())


def temoin_part_de_la_cellule(df, col='part_de_la_cellule', par=('x', 'y', 'year')):
    """Refuse si `col` ne somme pas a 1 PAR CELLULE."""
    cles = [c for c in par if c in df.columns]
    s = df.groupby(cles)[col].sum()
    mauvais = s[(s - 1).abs() > TOL]
    if len(mauvais):
        raise ValueError(
            "[STOP] temoin PAR CELLULE (%s) : %d cellule(s) sur %d ne somment pas a 1 "
            "(min %.4g, max %.4g). Une part normalisee PAR PAYS fournie a la place d'une "
            "part PAR CELLULE produit exactement ceci." % (col, len(mauvais), len(s),
                                                         s.min(), s.max()))
    temoin_part_non_triviale(df[col], col)
    return len(s)


def temoin_part_du_pays(df, col='part_du_pays', par=('GID_0', 'year')):
    """Refuse si `col` ne somme pas a 1 PAR PAYS."""
    cles = [c for c in par if c in df.columns]
    s = df.groupby(cles)[col].sum()
    mauvais = s[(s - 1).abs() > TOL]
    if len(mauvais):
        raise ValueError(
            "[STOP] temoin PAR PAYS (%s) : %d pays-periode(s) sur %d ne somment pas a 1 "
            "(min %.4g, max %.4g). Une part normalisee PAR CELLULE fournie a la place "
            "d'une part PAR PAYS produit exactement ceci." % (col, len(mauvais), len(s),
                                                            s.min(), s.max()))
    temoin_part_non_triviale(df[col], col)
    return len(s)
