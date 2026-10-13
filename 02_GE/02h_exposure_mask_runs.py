# -*- coding: utf-8 -*-
r"""02h_exposure_mask_runs.py -- the four GE runs of the source decomposition.

WHAT RUNS
    For each mask built by 02g_exposure_masks.py (north, south), the daily arm
    and the annual arm of the headline configuration:
        f_a = 0 (AMENITY_SPEC=null), central draw (QUANTILE=q50), CAP_COLMASS=none,
        annual arm read off the annual dose-response (ANNUAL_DRF_SOURCE=annual).
    Four runs, one draw. There is no band on the masks.

    The no-warming run is not repeated: with no warming there is nothing to mask,
    so the masked no-warming run is the headline one, welfare_sc1_q50_NULLamen.csv
    in output/. 02i_exposure_mask_decomposition.py reads it from there.

WHERE IT WRITES
    intermediate_outputs/exposure_masks/<mask>/  (JMP_DATA and JMP_OUT of main.py)
    intermediate_outputs/exposure_masks/runs_manifest.csv
    logs/exposure_masks_<timestamp>/<mask>_<arm>.log

    Not output/ge_runs_manifest.csv. That manifest is indexed by file name, and the
    masked runs carry the same file names as the headline runs: a line there would
    make every consumer of the headline grid see two rows for one file.

THE INPUT FINGERPRINT
    Same recipe as 02_GE/run_one_ge.sh: the md5 of the nine files' md5, one per
    line, followed by the scale convention read from param_theta. Two values are
    recorded per run:
      inputs_md5       the recipe on the files the run actually read (the masked
                       climate files included). It differs from the headline grid
                       by construction, since two of the nine files are masked.
      inputs_md5_base  the same recipe with the two climate files taken from
                       intermediate_outputs/. It MUST equal the inputs_md5 of the
                       headline runs in output/ge_runs_manifest.csv; the run
                       refuses to start otherwise. The masked set then differs from
                       the headline set in the two climate files and nothing else.
    wage0.csv is read by the solver and is not in the recipe; it is checked here
    by md5 against intermediate_outputs/ instead.

THE DOSE-RESPONSES ARE FINGERPRINTED TOO. The recipe above hashes calibration
inputs, not the dose-response tables the solver reads. The masked folders hold
COPIES of output/MC_damage_*.csv; a re-estimation after 02g would leave them
stale. Before any run, every copy must equal its source in output/, and each
manifest row records drf_md5 (md5 of the sorted 'name:md5' lines of the
copies), which 02i recomputes from output/ and refuses on a mismatch.

EVERY FILE A RUN WRITES IS RECORDED. Besides the welfare matrix, each row
carries the md5 of the two hat-paths: the amenity one, which 02i reads, and
the productivity one, which nothing reads and which is recorded so that no file
in these folders is unknown to a check.

ONE SOLVER AT A TIME; 14 GB of free memory required before each run (a solver
takes about 7 GB). STOPS AT THE FIRST FAILURE.

Usage:  python 02_GE/02h_exposure_mask_runs.py            # checks, then the four runs
        python 02_GE/02h_exposure_mask_runs.py --check    # checks only, no run
"""
import os, sys, hashlib, subprocess, time, datetime
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get('JMP_ROOT') or os.path.dirname(HERE)
INP = os.path.join(ROOT, 'intermediate_outputs')
DEST = os.path.join(INP, 'exposure_masks')
MAIN = os.path.join(HERE, '02c_model', 'main.py')
HEADLINE_MANIFEST = os.path.join(ROOT, 'output', 'ge_runs_manifest.csv')
MANIFEST = os.path.join(DEST, 'runs_manifest.csv')
PYTHON = os.environ.get('PYTHON') or sys.executable

SETTINGS = {'AMENITY_SPEC': 'null', 'QUANTILE': 'q50', 'CAP_COLMASS': 'none',
            'ANNUAL_DRF_SOURCE': 'annual', 'DUMP_OUTER': '0'}
ATAG = '_NULLamen'
HEADLINE_FILES = ['welfare_sc%d_q50%s.csv' % (s, ATAG) for s in (1, 2, 3)]
# Order and names of the nine fingerprinted files: run_one_ge.sh, INPUTS_MD5.
RECIPE = ['lambda.h5', 'mig_flows.h5', 'indices_drop.csv', 'identifier_xy_final.csv',
          'population0.csv', 'tas_ssp2_daily.csv', 'tas_ssp2_annual.csv',
          'data_final_population.csv', 'data_final_productivity.csv']
CLIMATE = ('tas_ssp2_daily.csv', 'tas_ssp2_annual.csv')
MEM_KB = 14_000_000


def md5(p):
    h = hashlib.md5()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def scale_convention():
    """param_theta at its single source, same rule as run_one_ge.sh: one
    occurrence, no fallback."""
    p = os.path.join(HERE, '02c_model', 'specification_study.csv')
    hits = [l for l in open(p, encoding='utf-8')
            if l.strip().startswith('param_theta') and l.split(';')[0].strip() == 'param_theta']
    if len(hits) != 1:
        sys.exit('[FATAL] cle param_theta absente ou multiple dans %s' % p)
    theta = hits[0].split(';')[1].strip()
    try:
        float(theta)
    except ValueError:
        sys.exit('[FATAL] param_theta illisible : %r' % hits[0])
    return 'SCALE_PROD=1/THETA_GE;THETA_GE=%s' % theta


def fingerprint(folder_of):
    """folder_of(file) -> folder holding it. Reproduces
    { md5sum f | cut ; ... ; printf '%s\n' "$SCALE_CONVENTION" ; } | md5sum."""
    lines = [md5(os.path.join(folder_of(f), f)) for f in RECIPE] + [scale_convention()]
    return hashlib.md5(('\n'.join(lines) + '\n').encode('ascii')).hexdigest()


def free_kb():
    r = subprocess.run(['powershell.exe', '-NoProfile', '-Command',
                        '(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory'],
                       capture_output=True, text=True)
    digits = ''.join(ch for ch in r.stdout if ch.isdigit())
    if not digits:
        sys.exit('[STOP] memoire libre illisible : %r' % r.stdout)
    return int(digits)


def wait_memory():
    for _ in range(60):
        kb = free_kb()
        if kb >= MEM_KB:
            return kb
        print('  [memoire] %d ko libres < 14 Go, attente 60 s' % kb, flush=True)
        time.sleep(60)
    sys.exit('[STOP] moins de 14 Go libres depuis une heure')


# ---- 1. the base fingerprint must be the headline grid's ---------------------
# dtype=str, keep_default_na=False : the manifest writes amenity_spec=null, which
# pandas would otherwise read as NaN.
man = pd.read_csv(HEADLINE_MANIFEST, dtype=str, keep_default_na=False).set_index('file')
missing = [f for f in HEADLINE_FILES if f not in man.index]
if missing:
    sys.exit('[STOP] absents de output/ge_runs_manifest.csv : %s' % missing)
head_fp = sorted(set(man.loc[HEADLINE_FILES, 'inputs_md5'].astype(str)))
if len(head_fp) != 1:
    sys.exit('[STOP] les runs de tete portent %d empreintes : %s' % (len(head_fp), head_fp))
for f in HEADLINE_FILES:
    if md5(os.path.join(ROOT, 'output', f)) != man.loc[f, 'md5']:
        sys.exit('[STOP] %s : md5 disque != manifeste' % f)
base_live = fingerprint(lambda f: INP)
if base_live != head_fp[0]:
    sys.exit('[STOP] empreinte des intrants vivants %s != empreinte des runs de tete %s. '
             'Quelque chose differe entre les intrants et la grille publiee.'
             % (base_live[:8], head_fp[0][:8]))
print('[fingerprint] intrants vivants = runs de tete : %s' % base_live[:8])

# ---- 2. each masked set differs from the live one in the climate files only --
fp_mask = {}
for mask in ('north', 'south'):
    dd = os.path.join(DEST, mask)
    if not os.path.isdir(dd):
        sys.exit('[STOP] %s absent : lancer 02g_exposure_masks.py' % dd)
    base_m = fingerprint(lambda f, dd=dd: INP if f in CLIMATE or f == 'identifier_xy_final.csv' else dd)
    if base_m != base_live:
        sys.exit('[STOP] %s : hors climat, les intrants masques different des intrants vivants.' % mask)
    if md5(os.path.join(dd, 'wage0.csv')) != md5(os.path.join(INP, 'wage0.csv')):
        sys.exit('[STOP] %s : wage0.csv differe' % mask)
    for f in CLIMATE:
        if md5(os.path.join(dd, f)) == md5(os.path.join(INP, f)):
            sys.exit('[STOP] %s : %s identique au vivant, le masque n\'a pas ete applique' % (mask, f))
    fp_mask[mask] = fingerprint(lambda f, dd=dd: INP if f == 'identifier_xy_final.csv' else dd)
    print('[fingerprint] %s : inputs_md5 %s | base %s' % (mask, fp_mask[mask][:8], base_m[:8]))

# ---- 2bis. the dose-response copies must equal output/ ----------------------
DRF = sorted(f for f in os.listdir(os.path.join(ROOT, 'output'))
             if f.startswith('MC_damage_') and f.endswith('.csv'))
drf_md5 = {}
for mask in ('north', 'south'):
    dd = os.path.join(DEST, mask)
    lines = []
    for f in DRF:
        q = os.path.join(dd, f)
        h = md5(os.path.join(ROOT, 'output', f))
        if not os.path.exists(q) or md5(q) != h:
            sys.exit('[STOP] %s : copie de %s absente ou perimee. Relancer '
                     '02g_exposure_masks.py.' % (mask, f))
        lines.append('%s:%s' % (f, h))
    drf_md5[mask] = hashlib.md5(('\n'.join(lines) + '\n').encode('ascii')).hexdigest()
    print('[drf] %s : %d tables identiques a output/ ; drf_md5 %s'
          % (mask, len(DRF), drf_md5[mask][:8]))

if '--check' in sys.argv:
    print('[check] empreintes verifiees ; aucun run lance (--check).')
    sys.exit(0)

# ---- 3. the four runs ---------------------------------------------------------
stamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
logdir = os.path.join(ROOT, 'logs', 'exposure_masks_' + stamp)
os.makedirs(logdir, exist_ok=True)
HDR = ['file', 'mask', 'md5', 'run_id', 'scenario', 'quantile', 'cap_colmass',
       'amenity_spec', 'annual_drf_source', 'inputs_md5', 'inputs_md5_base',
       'drf_md5', 'amenity_md5', 'productivity_md5', 'md5_provenance',
       'seconds', 'rc', 'timestamp']
rows = (pd.read_csv(MANIFEST, dtype=str, keep_default_na=False) if os.path.exists(MANIFEST)
        else pd.DataFrame(columns=HDR))

for mask in ('north', 'south'):
    dd = os.path.join(DEST, mask)
    for arm, sc in (('annual', 3), ('daily', 2)):
        kb = wait_memory()
        env = dict(os.environ, JMP_ROOT=ROOT, JMP_DATA=dd, JMP_OUT=dd,
                   SCENARIO=str(sc), **SETTINGS)
        log = os.path.join(logdir, '%s_%s.log' % (mask, arm))
        print('[run] %s %s arm : debut %s (%d ko libres)'
              % (mask, arm, datetime.datetime.now().strftime('%H:%M:%S'), kb), flush=True)
        t0 = time.time()
        with open(log, 'w', encoding='utf-8') as fh:
            rc = subprocess.run([PYTHON, '-u', MAIN], env=env, stdout=fh,
                                stderr=subprocess.STDOUT).returncode
        sec = int(time.time() - t0)
        print('[run] %s %s arm : rc=%d, %d min %02d s' % (mask, arm, rc, sec // 60, sec % 60), flush=True)
        if rc != 0:
            sys.exit('[STOP] %s %s arm en echec (rc=%d). Log : %s' % (mask, arm, rc, log))
        wf = 'welfare_sc%d_q50%s.csv' % (sc, ATAG)
        p = os.path.join(dd, wf)
        if not os.path.exists(p):
            sys.exit('[STOP] %s absent apres un run a rc=0' % p)
        row = {'file': wf, 'mask': mask, 'md5': md5(p), 'run_id': 'exposure_masks_' + stamp,
               'scenario': sc, 'quantile': SETTINGS['QUANTILE'],
               'cap_colmass': SETTINGS['CAP_COLMASS'], 'amenity_spec': SETTINGS['AMENITY_SPEC'],
               'annual_drf_source': SETTINGS['ANNUAL_DRF_SOURCE'],
               'inputs_md5': fp_mask[mask], 'inputs_md5_base': base_live,
               'drf_md5': drf_md5[mask],
               'amenity_md5': md5(os.path.join(dd, 'amenity_hat_path_sc%d_q50%s.csv' % (sc, ATAG))),
               'productivity_md5': md5(os.path.join(dd, 'productivity_hat_path_sc%d_q50%s.csv' % (sc, ATAG))),
               # Says WHO wrote the md5 columns. 'run' : this script, right after
               # the run. The four rows of 25/09/2026 carry a backfill note instead.
               'md5_provenance': 'run',
               'seconds': sec, 'rc': rc,
               'timestamp': datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S')}
        rows = rows[~((rows['file'] == wf) & (rows['mask'] == mask))]
        rows = pd.concat([rows, pd.DataFrame([row])], ignore_index=True)[HDR]
        # Literal file name on the writing line, so tools/audit_paquet.py sees the producer.
        rows.to_csv(os.path.join(DEST, 'runs_manifest.csv'), index=False)   # == MANIFEST

print('[runs] quatre runs termines ; manifeste : %s' % os.path.relpath(MANIFEST, ROOT))
