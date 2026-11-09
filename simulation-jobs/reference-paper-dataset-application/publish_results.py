from pathlib import Path
import json,hashlib,shutil,subprocess,zipfile,fcntl,traceback
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';REPO=Path('/workspace/R_Code');REL=Path('results')/ROOT.name;STATUS=ROOT/'github_publication_status.json'
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def state(**data):STATUS.write_text(json.dumps(data,indent=2)+'\n')
try:
 lock=Path('/workspace/paired-study-publication.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX)
 assert (FULL/'COMPLETED.txt').exists()
 meta=json.loads((FULL/'run_metadata.json').read_text());assert meta['patients']==100 and meta['imputations']==40 and meta['original_missing_responses']==32
 assert not git('status','--porcelain'),'Uncommitted repository work'
 git('fetch','origin','refs/heads/main');git('merge','--ff-only','FETCH_HEAD')
 dest=REPO/REL;assert not dest.exists();dest.mkdir(parents=True)
 archive=dest/'ReferencePaper_Dataset_Application_Results.zip'
 files=[p for p in FULL.rglob('*') if p.is_file() and p.suffix in ('.csv','.pdf','.png','.html','.txt','.log','.json','.rds')]
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,str(p.relative_to(FULL)))
  for name in ['original_script.R','run_all.R','run.sh','make_report.py','DATA_VERIFICATION.json','README.md']:z.write(ROOT/name,name)
  z.write(ROOT/'full_console.log','console_output.log')
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 assert archive.stat().st_size<95*1024**2
 for p in files:
  target=dest/p.relative_to(FULL);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 for name in ['original_script.R','run_all.R','run.sh','make_report.py']:shutil.copy2(ROOT/name,dest/name)
 links=f'\n[Full PDF report](simulation_report.pdf) · [Primary table PDF](main_tables.pdf) · [Forest plot](comparison_graph.png) · [Results ZIP]({archive.name})\n\n[Full-precision computed results](computed_results_full_precision.csv) · [Dataset verification](DATA_VERIFICATION.json).\n'
 (dest/'README.md').write_text((ROOT/'README.md').read_text()+links)
 checks={str(p.relative_to(dest)):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in dest.rglob('*') if p.is_file()};(dest/'SHA256SUMS.json').write_text(json.dumps(checks,indent=2)+'\n')
 p=REPO/'README.md';p.write_text(p.read_text()+f'\n[Reference-paper dataset application, 40 ABB imputations]({REL}/README.md).\n')
 git('add','-f','--',str(REL));git('add','--','README.md');git('-c','user.name=Codex','-c','user.email=codex@openai.com','commit','-m','Add reference paper dataset application results with 40 ABB imputations')
 commit=git('rev-parse','HEAD');state(state='pushing',commit=commit);git('push','origin','main');git('fetch','origin','refs/heads/main');assert git('rev-parse','FETCH_HEAD')==commit
 for name,sha in checks.items():
  p=subprocess.Popen(['git','-C',str(REPO),'show',f'{commit}:{REL}/{name}'],stdout=subprocess.PIPE)
  got=hashlib.file_digest(p.stdout,'sha256').hexdigest();assert p.wait()==0 and got==sha,name
 state(state='published_and_verified',commit=commit,verified_files=len(checks),url=f'https://github.com/HabibullahHabib/R_Code/tree/main/{REL}')
except Exception:
 state(state='publication_failed',error=traceback.format_exc());raise
