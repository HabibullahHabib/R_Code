from pathlib import Path
import json,hashlib,shutil,subprocess,zipfile,fcntl,traceback
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';REPO=Path('/workspace/R_Code');REL=Path('results')/ROOT.name;STATUS=ROOT/'github_publication_status.json'
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def state(**data):STATUS.write_text(json.dumps(data,indent=2)+'\n')
try:
 lock=Path('/workspace/paired-study-publication.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX)
 assert (FULL/'COMPLETED.txt').exists()
 meta=json.loads((FULL/'run_metadata.json').read_text());assert meta['patients']==369 and meta['imputations']==40 and meta['imputation_method']=='regression'
 assert not git('status','--porcelain'),'Uncommitted repository work'
 git('fetch','origin','refs/heads/main');git('merge','--ff-only','FETCH_HEAD')
 dest=REPO/REL;assert not dest.exists();dest.mkdir(parents=True)
 archive=dest/'IBScovars_ReferenceStyle_Regression_Results.zip'
 files=[p for p in FULL.rglob('*') if p.is_file() and p.suffix in ('.csv','.pdf','.png','.html','.txt','.log','.json','.rds','.tex')]
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,str(p.relative_to(FULL)))
  for name in ['original_script.R','executed_script.R','run_all.R','run.sh','make_report.py','render_highlighted_plot.R','prepare_latex.py','README.md']:z.write(ROOT/name,name)
  z.write(ROOT/'full_console.log','console_output.log')
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 assert archive.stat().st_size<95*1024**2
 for p in files:
  target=dest/p.relative_to(FULL);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 for name in ['original_script.R','executed_script.R','run_all.R','run.sh','make_report.py','render_highlighted_plot.R','prepare_latex.py']:shutil.copy2(ROOT/name,dest/name)
 links=f'\n[Full PDF report](simulation_report.pdf) · [Primary table PDF](main_tables.pdf) · [Forest plot](comparison_graph.png) · [Results ZIP]({archive.name})\n\n[Full-precision pooled results](computed_pooled_results_full_precision.csv) · [Compiled LaTeX example](latex-ready/Compile_Application.pdf).\n'
 (dest/'README.md').write_text((ROOT/'README.md').read_text()+links)
 checks={str(p.relative_to(dest)):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in dest.rglob('*') if p.is_file()};(dest/'SHA256SUMS.json').write_text(json.dumps(checks,indent=2)+'\n')
 p=REPO/'README.md';p.write_text(p.read_text()+f'\n[Corrected IBScovars application, 40 regression imputations]({REL}/README.md).\n')
 older=REPO/'results/ibscovars-reference-style-application/README.md'
 if older.exists():
  previous=older.read_text();notice='**Superseded for the user-requested regression analysis.** Use [the corrected regression result](../ibscovars-reference-style-regression/README.md). This earlier propensity result is retained for audit.\n\n'
  if not previous.startswith('**Superseded'):older.write_text(notice+previous)
  oldchecks=older.parent/'SHA256SUMS.json'
  if oldchecks.exists():
   c=json.loads(oldchecks.read_text());c['README.md']=hashlib.file_digest(older.open('rb'),'sha256').hexdigest();oldchecks.write_text(json.dumps(c,indent=2)+'\n')
  git('add','--',str(older.relative_to(REPO)),str(oldchecks.relative_to(REPO)))
 git('add','-f','--',str(REL));git('add','--','README.md');git('-c','user.name=Codex','-c','user.email=codex@openai.com','commit','-m','Correct IBScovars application to regression and label proposed procedure results')
 commit=git('rev-parse','HEAD');state(state='pushing',commit=commit);git('push','origin','main');git('fetch','origin','refs/heads/main');assert git('rev-parse','FETCH_HEAD')==commit
 for name,sha in checks.items():
  p=subprocess.Popen(['git','-C',str(REPO),'show',f'{commit}:{REL}/{name}'],stdout=subprocess.PIPE)
  got=hashlib.file_digest(p.stdout,'sha256').hexdigest();assert p.wait()==0 and got==sha,name
 state(state='published_and_verified',commit=commit,verified_files=len(checks),url=f'https://github.com/HabibullahHabib/R_Code/tree/main/{REL}')
except Exception:
 state(state='publication_failed',error=traceback.format_exc());raise
