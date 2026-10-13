#!/usr/bin/env python
# -*- coding: utf-8 -*-
r"""cmp_rerun.py -- comparer un rejeu a l'etat publie, categorie par categorie.

CE QU'IL FAIT
-------------
Compare deux repertoires fichier par fichier et CLASSE chaque ecart, au lieu de
rendre un oui/non. La distinction qui compte :

  identique        md5 egal, rien a dire
  horodatage       PDF dont seuls /CreationDate, /ModDate et /ID different
  arrondi          CSV numerique dont tous les ecarts sont < tol
  DIVERGENT        tout le reste : valeur changee, colonne changee, forme changee
  ABSENT_REF       present dans le rejeu, absent des archives
  ABSENT_NEW       present dans les archives, absent du rejeu
  ILLISIBLE        n'a pas pu etre compare -> NE PASSE PAS AU VERT

⚠️ CE FICHIER MANQUAIT. `tools/README.md` et `tools/AMORCE.md` le prescrivent
depuis le 17/09 (« rejeux sous JMP_OUTPUT_DIR puis cmp_rerun.py avant toute
copie ») et la section reproductibilite du README decrit meme ses categories.
Il n'a jamais ete ecrit : la discipline citait un instrument inexistant, et
chaque rejeu etait donc compare a la main ou pas du tout. Ecrit le 19/09/2026.

POURQUOI « ILLISIBLE » NE PASSE PAS AU VERT
-------------------------------------------
Meme regle que les autres outils du dossier : rc=0 vert, rc=1 ecart trouve,
rc=2 impossible de conclure. Un fichier qu'on ne sait pas lire n'est pas un
fichier sans ecart -- c'est la distinction qui fait qu'un marqueur disparu ne
se lit pas comme une absence de defaut.

Usage:
  python tools/cmp_rerun.py --ref output --new output_01a_ctrl
  python tools/cmp_rerun.py --ref output --new output_01a_ctrl --only "MC_damage_prod*"
  python tools/cmp_rerun.py --ref A --new B --tol 1e-9
"""
import argparse
import fnmatch
import hashlib
import os
import re
import sys

PDF_VOLATILE = re.compile(
    rb"/(CreationDate|ModDate|ID)\s*(\([^)]*\)|\[[^\]]*\]|<[^>]*>)")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pdf_canonical(path):
    """Le PDF prive de ses champs volatils."""
    with open(path, "rb") as fh:
        return PDF_VOLATILE.sub(b"", fh.read())


def compare_csv(a, b, tol):
    """(categorie, detail). Compare valeur par valeur, colonnes comprises."""
    try:
        import pandas as pd
    except ImportError:
        return "ILLISIBLE", "pandas absent"
    try:
        da = pd.read_csv(a)
        db = pd.read_csv(b)
    except Exception as exc:                       # noqa: BLE001
        return "ILLISIBLE", f"lecture impossible : {exc}"

    if list(da.columns) != list(db.columns):
        only_a = [c for c in da.columns if c not in db.columns]
        only_b = [c for c in db.columns if c not in da.columns]
        return "DIVERGENT", (f"colonnes differentes "
                             f"(ref seul: {only_a or '-'}, "
                             f"new seul: {only_b or '-'})")
    if da.shape != db.shape:
        return "DIVERGENT", f"forme {da.shape} vs {db.shape}"

    worst, worst_col = 0.0, None
    for c in da.columns:
        sa, sb = da[c], db[c]
        if sa.dtype.kind in "biufc" and sb.dtype.kind in "biufc":
            import numpy as np
            va, vb = sa.to_numpy(dtype=float), sb.to_numpy(dtype=float)
            both_nan = np.isnan(va) & np.isnan(vb)
            d = np.abs(va - vb)
            d[both_nan] = 0.0
            if np.isnan(d).any():
                return "DIVERGENT", f"colonne '{c}' : NaN d'un seul cote"
            m = float(d.max()) if d.size else 0.0
            if m > worst:
                worst, worst_col = m, c
        else:
            if not sa.astype(str).equals(sb.astype(str)):
                return "DIVERGENT", f"colonne texte '{c}' differente"
    if worst == 0.0:
        return "identique", "valeurs egales au bit (md5 differe : en-tete ou fin de ligne)"
    if worst < tol:
        return "arrondi", f"ecart max {worst:.3e} sur '{worst_col}' (< {tol:g})"
    return "DIVERGENT", f"ecart max {worst:.3e} sur '{worst_col}'"


def compare_text(a, b):
    try:
        la = open(a, encoding="utf-8", errors="replace").read().splitlines()
        lb = open(b, encoding="utf-8", errors="replace").read().splitlines()
    except Exception as exc:                       # noqa: BLE001
        return "ILLISIBLE", f"lecture impossible : {exc}"
    if la == lb:
        return "identique", "lignes egales (md5 differe : fin de ligne)"
    n = sum(1 for x, y in zip(la, lb) if x != y) + abs(len(la) - len(lb))
    first = next((i + 1 for i, (x, y) in enumerate(zip(la, lb)) if x != y), None)
    return "DIVERGENT", f"{n} ligne(s) differente(s), la premiere en l.{first}"


def compare_one(ref, new, tol):
    if md5(ref) == md5(new):
        return "identique", "md5 egal"
    ext = os.path.splitext(ref)[1].lower()
    if ext == ".pdf":
        if pdf_canonical(ref) == pdf_canonical(new):
            return "horodatage", "identique hors /CreationDate /ModDate /ID"
        return "DIVERGENT", "contenu PDF modifie"
    if ext in (".csv", ".tsv"):
        return compare_csv(ref, new, tol)
    if ext in (".tex", ".txt", ".md", ".log"):
        return compare_text(ref, new)
    return "DIVERGENT", f"md5 differe (type {ext or 'sans extension'} non interprete)"


def lister(racine, motif):
    """Tous les fichiers SOUS `racine`, en chemins relatifs a separateur '/'.

    ⚠️ RECURSIF DEPUIS LE 20/09/2026. Cette fonction etait `os.listdir`, donc le
    balayage s'arretait au premier niveau. Sur output/, cela laissait CINQ
    sous-repertoires hors du bilan -- contribA/, crh_inference/, drf_prod_specs/,
    past_outputs/, _retires_20260911/ -- sans les declarer ni ABSENT ni
    ILLISIBLE : ils n'apparaissaient nulle part, et le total se lisait comme
    complet. Le 20/09, contribA/ portait cinq divergences reelles, dont
    test1_delta_north.csv, et drf_prod_specs/ douze : dix-sept ecarts invisibles
    sur un outil qui affichait « 27 ecart(s) reel(s) ».

    C'est la faute que ce fichier reproche aux autres : un controle qui ne
    regarde pas ne doit pas se lire comme un controle qui n'a rien trouve.

    Le motif s'applique au nom de base OU au chemin relatif, de sorte que
    --only "MC_damage_prod*" continue de designer un fichier ou qu'il soit, et
    que --only "contribA/*" designe un sous-arbre.
    """
    out = set()
    for dossier, _sous, fichiers in os.walk(racine):
        for nom in fichiers:
            plein = os.path.join(dossier, nom)
            if not os.path.isfile(plein):
                continue
            rel = os.path.relpath(plein, racine).replace(os.sep, "/")
            if fnmatch.fnmatch(nom, motif) or fnmatch.fnmatch(rel, motif):
                out.add(rel)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--ref", required=True, help="repertoire de reference (publie)")
    p.add_argument("--new", required=True, help="repertoire du rejeu")
    p.add_argument("--tol", type=float, default=1e-6)
    p.add_argument("--only", default="*", help="motif fnmatch sur le nom de fichier")
    p.add_argument("--quiet-identical", action="store_true")
    a = p.parse_args()

    for d in (a.ref, a.new):
        if not os.path.isdir(d):
            print(f"[STOP] repertoire absent : {d}")
            return 2

    fa = lister(a.ref, a.only)
    fb = lister(a.new, a.only)

    rows = []
    for f in sorted(fb - fa):
        rows.append((f, "ABSENT_REF", "present dans le rejeu seulement"))
    for f in sorted(fa - fb):
        rows.append((f, "ABSENT_NEW", "present dans les archives seulement"))
    for f in sorted(fa & fb):
        cat, det = compare_one(os.path.join(a.ref, f), os.path.join(a.new, f), a.tol)
        rows.append((f, cat, det))

    order = ["DIVERGENT", "ILLISIBLE", "ABSENT_NEW", "ABSENT_REF",
             "arrondi", "horodatage", "identique"]
    counts = {k: 0 for k in order}
    print(f"=== cmp_rerun : ref={a.ref}  new={a.new}  tol={a.tol:g}  motif={a.only} ===")
    for cat in order:
        sel = [r for r in rows if r[1] == cat]
        counts[cat] = len(sel)
        if not sel:
            continue
        if cat == "identique" and a.quiet_identical:
            continue
        print(f"\n-- {cat} ({len(sel)}) --")
        for f, _, det in sel:
            print(f"   {f:58s} {det}")

    print("\n=== BILAN ===")
    for cat in order:
        print(f"   {cat:12s} {counts[cat]:4d}")

    # `ABSENT_REF` n'est pas un defaut : un rejeu peut produire un artefact neuf.
    # `ABSENT_NEW` en est un : une cible declaree que le rejeu n'a pas produite.
    bad = counts["DIVERGENT"] + counts["ABSENT_NEW"]
    if counts["ILLISIBLE"]:
        print("\nVERDICT : rc=2 -- des fichiers n'ont pas pu etre compares.")
        return 2
    if bad:
        print(f"\nVERDICT : rc=1 -- {bad} ecart(s) reel(s).")
        return 1
    print("\nVERDICT : rc=0 -- le rejeu reproduit l'etat de reference.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
