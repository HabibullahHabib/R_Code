from pathlib import Path
import subprocess,json,hashlib,shutil,zipfile,csv,fcntl,traceback
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';REPO=Path('/workspace/R_Code');REL=Path('results')/ROOT.name
STATUS=ROOT/'github_publication_status.json';URL='https://github.com/HabibullahHabib/R_Code.git'
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def state(**d):STATUS.write_text(json.dumps(d,indent=2)+'\n')
try:
 lock=Path('/workspace/paired-study-publication.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX)
 assert (FULL/'COMPLETED.txt').exists()
 row=next(csv.DictReader((FULL/'completion_summary.csv').open()));assert int(row['attempted'])==10000 and int(row['imputations_per_repetition'])==40 and int(row['successful'])+int(row['failed'])==10000
 assert not git('status','--porcelain'),'Repository has uncommitted work'
 assert git('branch','--show-current')=='main' and git('remote','get-url','origin')==URL
 git('fetch','origin','refs/heads/main');git('merge','--ff-only','FETCH_HEAD')
 dest=REPO/REL;assert not dest.exists(),'Existing result folder requires review';dest.mkdir(parents=True)
 archive=dest/(ROOT.name+'_Results.zip')
 files=[p for p in FULL.iterdir() if p.is_file() and p.suffix in ('.csv','.pdf','.png','.txt','.log','.json')]
 with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,p.name)
  for name in ['original_application.R','run_study.R','run.sh','export_report.R','README.md']:z.write(ROOT/name,name)
  z.write(ROOT/'full_console.log','console_output.log')
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 assert archive.stat().st_size<95*1024**2
 for p in files:shutil.copy2(p,dest/p.name)
 for name in ['original_application.R','run_study.R','run.sh','export_report.R']:shutil.copy2(ROOT/name,dest/name)
 (dest/'README.md').write_text((ROOT/'README.md').read_text()+f'\nCompleted {row["attempted"]} attempted repetitions: {row["successful"]} successful and {row["failed"]} failed.\n\n[PDF report](simulation_report.pdf) · [Plain-language graph](comparison_graph.png) · [Graph PDF](comparison_graph.pdf) · [Results ZIP]({archive.name})\n\n[Summary CSV](comparison_summary.csv) · [Paired comparison CSV](paired_comparison.csv) · [Raw paired repetitions](paired_repetitions.csv) · [Raw intervals](raw_intervals.csv). Protocol and source files document the exact settings.\n')
 checks={p.name:hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in dest.iterdir() if p.is_file()};(dest/'SHA256SUMS.json').write_text(json.dumps(checks,indent=2)+'\n')
 readme=REPO/'README.md';text=readme.read_text();text+=f'\n[Paired regression-only study: {ROOT.name}]({REL}/README.md).\n';readme.write_text(text)
 git('add','-f','--',str(REL));git('add','--','README.md');git('-c','user.name=Codex','-c','user.email=codex@openai.com','commit','-m',f'Add completed paired shrinkage study {ROOT.name}')
 commit=git('rev-parse','HEAD');state(state='pushing',commit=commit);git('push','origin','main')
 remote=git('ls-remote','origin','refs/heads/main').split()[0];assert remote==commit
 git('fetch','origin','refs/heads/main');fetched=git('rev-parse','FETCH_HEAD');assert fetched==commit
 for name,sha in checks.items():
  p=subprocess.Popen(['git','-C',str(REPO),'show',f'{fetched}:{REL}/{name}'],stdout=subprocess.PIPE)
  got=hashlib.file_digest(p.stdout,'sha256').hexdigest();assert p.wait()==0 and got==sha,name
 state(state='published_and_verified',commit=commit,verified_files=len(checks),url=f'https://github.com/HabibullahHabib/R_Code/tree/main/{REL}')
except Exception:
 state(state='publication_failed',error=traceback.format_exc());raise
