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
  z.write(ROOT/'original_script.R','original_script.R');z.write(ROOT/'full_console.log','console_output.log')
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
 (dest/'README.md').write_text(f'''# {ROOT.name}\n\nSample-size sensitivity analysis completed. See completion_summary.csv and run_summary.txt for exact replicate and success/failure counts. Forty imputations. Original statistical code and seeds preserved; the execution wrapper records the requested replicate count. Seven sample sizes (20,50,100,150,200,250,300 per arm), two mean scenarios and five missing rates. Shrinkage and non-shrinkage results are included.\n\n[Tables, plots, logs and source ZIP]({archive.name}) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf)\n\nOriginal primary tables are provided individually as CSV and HTML. Raw replicate CSVs are gzip compressed and provided separately. Files ending .part001, .part002, etc. are consecutive binary segments of a single gzip file; concatenate in numeric order before decompressing. SHA256SUMS.json verifies downloaded files. No reduced-replicate trial results are included. Supplementary plots show coverage and mean interval length across sample sizes; CSVs retain the original detailed tables.\n''')
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
