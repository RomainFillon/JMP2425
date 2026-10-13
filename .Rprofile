source("renv/activate.R")
# =============================================================================
# .Rprofile -- racine du depot de replication
#
# R lit ce fichier AUTOMATIQUEMENT au demarrage, quand le repertoire de travail
# est celui qui le contient. Il n'y a donc rien a sourcer a la main : ouvrir R
# depuis la racine du depot suffit.
#
# CE QU'IL FAIT, ET POURQUOI
#   Il pose JMP_ROOT une fois, deduite de l'emplacement reel du depot. Aucun
#   script n'a plus besoin d'ecrire un chemin absolu : la racine ne peut pas
#   diverger du depot puisqu'elle n'est ecrite nulle part.
#
#   C'est ce qui remplace les douze redefinitions de .jmp_root qui portaient
#   chacune le meme chemin en dur -- douze endroits ou une machine differente
#   aurait fait echouer le paquet.
#
# ⚠ UNE VARIABLE DEJA POSEE N'EST PAS ECRASEE. L'orchestrateur .sh pose
# JMP_ROOT explicitement avant de lancer R ; s'il le fait, c'est lui qui decide.
# =============================================================================

local({
  if (nzchar(Sys.getenv("JMP_ROOT"))) {
    message("[.Rprofile] JMP_ROOT deja pose : ", Sys.getenv("JMP_ROOT"))
    return(invisible(NULL))
  }

  root <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)

  ## Controle de vraisemblance : on verifie qu'on est bien a la racine du
  ## depot, et pas dans un sous-dossier ou ailleurs. Sans lui, ouvrir R depuis
  ## 01_estimation/ poserait JMP_ROOT une couche trop bas, et tous les chemins
  ## seraient faux d'un cran -- sans qu'aucune erreur ne se leve avant le
  ## premier fread.
  attendus <- c("01_estimation", "data_inputs")
  manquants <- attendus[!dir.exists(file.path(root, attendus))]

  if (length(manquants)) {
    warning("[.Rprofile] ", root, "\n",
            "  ne ressemble pas a la racine du depot : ",
            paste(manquants, collapse = ", "), " absent(s).\n",
            "  JMP_ROOT N'A PAS ETE POSE. Ouvrir R depuis la racine, ou poser\n",
            "  JMP_ROOT a la main avant de lancer un script.",
            call. = FALSE)
    return(invisible(NULL))
  }

  Sys.setenv(JMP_ROOT = root)
  message("[.Rprofile] JMP_ROOT = ", root)
})

# =============================================================================
# Reglages de session
#
# Ils ne changent aucun resultat : ils rendent les sorties lisibles et evitent
# deux comportements par defaut qui surprennent sur un poste different.
# =============================================================================

options(
  ## ⚠ stringsAsFactors est FALSE par defaut depuis R 4.0, mais l'ecrire rend
  ## le paquet explicite sur une convention dont dependent plusieurs lectures.
  stringsAsFactors = FALSE,
  ## Sur des telechargements de plusieurs Go, le defaut de 60 s coupe -- et
  ## download.file rend alors un fichier TRONQUE sans lever d'erreur.
  timeout = 7200,
  scipen = 10,          # eviter la notation scientifique dans les CSV ecrits
  width = 100,
  warn = 1              # afficher les avertissements a mesure, pas a la fin
)
