from pathlib import Path
import subprocess,os,time,json,hashlib,tarfile,shutil,fcntl,traceback,datetime,csv
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';REPO=Path('/workspace/R_Code');REL=Path('results/litter-regression-10000-checkpoints');STATUS=ROOT/'controller_status.json'
def write_status(**x):
 x['updated_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();p=STATUS.with_suffix('.tmp');p.write_text(json.dumps(x,indent=2));p.replace(STATUS)
def git(*a):return subprocess.check_output(['git','-C',str(REPO),*a],text=True,stderr=subprocess.STDOUT).strip()
def checkpoint_files():return sorted(p for p in (FULL/'checkpoints').glob('*.rds') if p.with_suffix('.count').exists())
def count():return sum(int(p.with_suffix('.count').read_text()) for p in checkpoint_files())
def sync():
 files=checkpoint_files();saved=sum(int(p.with_suffix('.count').read_text()) for p in files)
 tmp=ROOT/'checkpoint_snapshot.tar.gz.tmp';manifest={'saved_repetitions':saved,'requested_repetitions':10000,'imputations_per_repetition':40,'sha256':{},'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 selected=files+[p.with_suffix('.count') for p in files]
 selected += [p for p in ROOT.iterdir() if p.is_file() and p.suffix in ('.R','.py','.sh','.md','.json','.rds') and p.name not in ('controller_status.json','remote_checkpoint_status.json')]
 selected += [p for p in [FULL/'protocol.rds',FULL/'progress.csv'] if p.exists()]
 with tarfile.open(tmp,'w:gz') as tar:
  for p in selected:
   key=str(p.relative_to(ROOT));manifest['sha256'][key]=hashlib.sha256(p.read_bytes()).hexdigest();tar.add(p,arcname=key,recursive=False)
 with tarfile.open(tmp,'r:gz') as tar:
  for key,sha in manifest['sha256'].items():assert hashlib.sha256(tar.extractfile(key).read()).hexdigest()==sha,key
 lock=Path('/workspace/paired-study-publication.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX)
 try:
  assert not git('status','--porcelain'),'Repository has other uncommitted changes'
  git('fetch','origin','refs/heads/main');git('merge','--ff-only','FETCH_HEAD')
  dest=REPO/REL;dest.mkdir(exist_ok=True,parents=True);shutil.copyfile(tmp,dest/'checkpoint_snapshot.tar.gz')
  (dest/'manifest.json').write_text(json.dumps(manifest,indent=2));(dest/'README.md').write_text(f'# Litter 10000-repetition run: checkpoint backup\n\n{saved}/10000 repetitions saved; 40 regression imputations per repetition. This is a resumable checkpoint snapshot, not the final result. Restore the archive under /workspace/litter-regression-10000 and run controller.py with the R runtime available.\n')
  git('add','-f','--',str(REL));git('-c','user.name=Codex','-c','user.email=codex@openai.com','commit','-m',f'Save Litter validation checkpoints {saved}/10000');commit=git('rev-parse','HEAD');git('push','origin','main');git('fetch','origin','refs/heads/main');assert git('rev-parse','FETCH_HEAD')==commit
  expected=hashlib.sha256((dest/'checkpoint_snapshot.tar.gz').read_bytes()).hexdigest();p=subprocess.Popen(['git','-C',str(REPO),'show',f'{commit}:{REL}/checkpoint_snapshot.tar.gz'],stdout=subprocess.PIPE);assert hashlib.file_digest(p.stdout,'sha256').hexdigest()==expected;assert p.wait()==0
  (ROOT/'remote_checkpoint_status.json').write_text(json.dumps({'state':'published_and_verified','saved_repetitions':saved,'commit':commit,'url':f'https://github.com/HabibullahHabib/R_Code/tree/main/{REL}'},indent=2));return saved
 finally:fcntl.flock(lock,fcntl.LOCK_UN);lock.close()
lock=(ROOT/'controller.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
try:
 remote=0
 try:remote=sync()
 except Exception:print('Initial remote snapshot failed:',traceback.format_exc(),flush=True)
 if not (FULL/'COMPLETED.txt').exists():
  with (ROOT/'full_console.log').open('ab') as log:
   p=subprocess.Popen(['bash',str(ROOT/'run.sh')],stdout=log,stderr=subprocess.STDOUT)
   (ROOT/'r_process.pid').write_text(str(p.pid));write_status(state='running',pid=p.pid,saved_repetitions=count(),requested_repetitions=10000,remote_saved=remote)
   while p.poll() is None:
    time.sleep(30);saved=count()
    if saved>=remote+500:
     try:remote=sync()
     except Exception:print('Remote snapshot failed:',traceback.format_exc(),flush=True)
    write_status(state='running',pid=p.pid,saved_repetitions=saved,requested_repetitions=10000,remote_saved=remote)
   if p.returncode:raise RuntimeError(f'R process exited {p.returncode}; see full_console.log')
 write_status(state='building_reports',saved_repetitions=count(),requested_repetitions=10000)
 with (ROOT/'report_generation.log').open('ab') as log:subprocess.run(['python',str(ROOT/'finish_outputs.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
 with (ROOT/'github_publication.log').open('ab') as log:subprocess.run(['python',str(ROOT/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT,check=True)
 write_status(state='completed_and_published',saved_repetitions=count(),requested_repetitions=10000)
except Exception:
 write_status(state='failed',saved_repetitions=count(),requested_repetitions=10000,error=traceback.format_exc());print(traceback.format_exc(),flush=True);raise
