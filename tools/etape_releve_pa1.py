"""etape_releve_pa1.py -- le point d'arret 1 devenu RELEVE (lot 57, 28/09/2026).

Author's decision (28/09/2026) : la course va au bout ; apres 02d (lot 85x, 11/10/2026 ;
apres ct7 jusque-la), run_all.sh appelle cette etape, NON BLOQUANTE, qui ecrit le
releve et rend la main a la chaine (inventaire des cibles). Elle ne corrige rien,
n'interprete rien et sort toujours en 0 : un echec de l'un de ses morceaux est ECRIT
dans le releve, pas propage. Elle passe apres 02d parce que emit_values_empirics.py
lit output/config_params.csv et que les controles lisent output/values.tex, deux
sorties de 02d : apres ct7, un clone neuf rendait « config_params.csv absent ».

  1. emit_values_empirics.py (mode rapport de la reprise declaree) : reecrit
     output/values_empirics.tex et output/ecarts_reprise_empirics.tsv, que 02d
     vient d'ecrire a partir des memes entrees ;
  2. les neuf controles de tools/reference_phase.py, chacun avec son code ;
  3. tools/releve_pa1.py : \\EMP* avant/apres, paires Lo/Hi, tests joints,
     effectifs, pays, jours qui changent de bin, et le verdict des conditions 2-4.

Sortie : logs/releve_pa1_<AAAAMMJJ_HHMM>.txt (+ .controles.txt).
"""
import datetime
import os
import subprocess
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'tools'))
from reference_phase import CONTROLES   # noqa: E402  la liste nommee, pas une copie

stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M')
os.makedirs(os.path.join(ROOT, 'logs'), exist_ok=True)
SORTIE = os.path.join(ROOT, 'logs', 'releve_pa1_%s.txt' % stamp)
CTRL = SORTIE.replace('.txt', '.controles.txt')
env = dict(os.environ, PYTHONIOENCODING='utf-8')


def lance(args, nom):
    r = subprocess.run([sys.executable] + args, cwd=ROOT, env=env, capture_output=True,
                       text=True, encoding='utf-8', errors='replace')
    print('[releve_pa1] %-28s rc=%d' % (nom, r.returncode), flush=True)
    return r


e = lance([os.path.join(ROOT, '02_GE', 'emit_values_empirics.py')], 'emit_values_empirics.py')
queue = [l for l in e.stdout.splitlines() if 'REPORT MODE' in l or '[write]' in l][-6:]
lignes = ['emit_values_empirics.py rc=%d' % e.returncode] + ['   ' + l for l in queue]
if e.returncode != 0:
    lignes += ['   ' + l for l in (e.stdout + e.stderr).splitlines()[-15:]]
for c in CONTROLES:
    r = lance([os.path.join(ROOT, 'tools', c + '.py')], c)
    lignes.append('%-24s rc=%d' % (c, r.returncode))
    if r.returncode != 0:
        lignes += ['      ' + l for l in (r.stdout + r.stderr).strip().splitlines()[-6:]]
with open(CTRL, 'w', encoding='utf-8') as fh:
    fh.write('\n'.join(lignes) + '\n')
r = lance([os.path.join(ROOT, 'tools', 'releve_pa1.py'), SORTIE, '--controles', CTRL],
          'releve_pa1.py')
print(r.stdout[-3000:])
if r.returncode not in (0, 1):
    print(r.stderr[-3000:])
print('[releve_pa1] record written: %s (conditions 2-4: %s) -- NON-BLOCKING, the chain goes on'
      % (os.path.relpath(SORTIE, ROOT), {0: 'PASS', 1: 'FAIL'}.get(r.returncode, 'INCONCLUSIVE')))
sys.exit(0)
