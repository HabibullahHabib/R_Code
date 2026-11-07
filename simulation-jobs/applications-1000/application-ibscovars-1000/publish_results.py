from pathlib import Path
import subprocess,hashlib,json,shutil,zipfile,tempfile,fcntl
ROOT=Path(__file__).resolve().parent; FULL=ROOT/'full'; REPO=Path('/workspace/R_Code'); REL=Path('results')/ROOT.name
URL='https://github.com/HabibullahHabib/R_Code.git'
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def state(**d): (ROOT/'github_publication_status.json').write_text(json.dumps(d,indent=2)+'\n')
try:
 assert (FULL/'COMPLETED.txt').exists()
 lock=Path('/workspace/paper-feedback-publication.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX)
 assert not git('status','--porcelain'),'Repository has uncommitted work'
 assert git('branch','--show-current')=='main' and git('remote','get-url','origin')==URL
 git('fetch','origin','main');git('merge','--ff-only','origin/main')
 dest=REPO/REL;assert not dest.exists();dest.mkdir(parents=True)
 archive=dest/(ROOT.name+'_Tables_Plots_and_Logs.zip')
 files=[p for p in FULL.rglob('*') if p.is_file() and p.suffix in ('.csv','.pdf','.html','.txt','.log')]
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,str(p.relative_to(FULL)))
  z.write(ROOT/'original_script.R','original_script.R');z.write(ROOT/'run_all.R','run_all.R');z.write(ROOT/'run.sh','run.sh');z.write(ROOT/'export_report.R','export_report.R');z.write(ROOT/'full_console.log','console_output.log')
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 for p in files:
  target=dest/p.relative_to(FULL);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 shutil.copy2(ROOT/'original_script.R',dest/'original_script.R')
 for p in FULL.rglob('raw_*.csv.gz'):
  if p.stat().st_size<90*1024**2:shutil.copy2(p,dest/p.name)
  else:
   with p.open('rb') as f:
    i=1
    while block:=f.read(90*1024**2):
     (dest/(p.name+f'.part{i:03d}')).write_bytes(block);i+=1
 checks={str(p.relative_to(dest)):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in dest.rglob('*') if p.is_file()}
 (dest/'SHA256SUMS.json').write_text(json.dumps(checks,indent=2)+'\n')
 (dest/'README.md').write_text(f"""# {ROOT.name}

Completed **1,000 repeated-deletion experiments total**. ACTG175 uses 1,000 experiments at 30% MCAR on the observed-outcome benchmark subset; its original natural-missingness primary analysis is also retained. IBScovars uses 250 experiments at each of 10%, 20%, 30% and 40% missingness. Every experiment evaluates all three MI methods with 40 imputations and both shrinkage and non-shrinkage intervals.

[Results ZIP]({archive.name}) · [Full PDF report](simulation_report.pdf) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf)

CSV tables and raw validation results are included, together with source, logs, warnings/errors and session information. Benchmark recovery is measured against fixed complete-sample contrasts, not population parameters. Both variance versions use the same imputations. All original statistical functions and seeds are preserved; only validation_reps is overridden. No increase beyond 1,000 experiments is enabled.
""")
 readme=REPO/'README.md';t=readme.read_text();t+=f'\n[Professor feedback: {ROOT.name}]({REL}/README.md).\n';readme.write_text(t)
 git('add','--','README.md',str(REL));git('-c','user.name=Codex','-c','user.email=codex@openai.com','commit','-m',f'Add completed {ROOT.name} simulation results')
 commit=git('rev-parse','HEAD');state(state='pushing',commit=commit);git('push','origin','main')
 remote=git('ls-remote','origin','refs/heads/main').split()[0];assert remote==commit
 # Fetch back and verify every published result against its local checksum.
 git('fetch','origin','main')
 for name,sha in checks.items():
  proc=subprocess.Popen(['git','-C',str(REPO),'show',f'FETCH_HEAD:{REL}/{name}'],stdout=subprocess.PIPE)
  assert hashlib.file_digest(proc.stdout,'sha256').hexdigest()==sha;assert proc.wait()==0
 state(state='published_and_verified',commit=commit,url=f'https://github.com/HabibullahHabib/R_Code/tree/main/{REL}')
except Exception as e:
 state(state='publication_failed',error=str(e));raise
