from pathlib import Path
import os,subprocess,time,json,fcntl
BASE=Path('/workspace');NAMES=['paired-actg175-natural-regression','paired-ibscovars-mar30-regression'];lock=(BASE/'paired-study-controller.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
children={};published=set();failed=set()
def state(**d):(BASE/'paired-study-status.json').write_text(json.dumps(d,indent=2)+'\n')
for name in NAMES:
 root=BASE/name;full=root/'full';assert (root/'smoke/COMPLETED.txt').exists() and (root/'smoke/FIDELITY_PASSED.txt').exists()
 if (full/'COMPLETED.txt').exists():continue
 active=False
 try:
  pid=int((root/'full_run.pid').read_text());active=str(root).encode() in Path(f'/proc/{pid}/cmdline').read_bytes() and Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0] not in ('Z','T','t')
 except (OSError,ValueError):pass
 if active:children[name]=(pid,None);continue
 (full/'FAILED.txt').unlink(missing_ok=True)
 with (root/'full_console.log').open('ab') as log:
  p=subprocess.Popen([str(root/'run.sh'),'full'],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
 (root/'full_run.pid').write_text(str(p.pid));children[name]=(p.pid,p)
while True:
 running=[]
 for name in NAMES:
  root=BASE/name;full=root/'full'
  if (full/'COMPLETED.txt').exists():
   if name not in published:
    pub=root/'github_publication_status.json'
    if not pub.exists() or json.loads(pub.read_text()).get('state')!='published_and_verified':
     with (root/'github_publication.log').open('ab') as log:subprocess.run(['python',str(root/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT)
    published.add(name)
   continue
  if name in failed:continue
  pid,child=children[name]
  try:alive=Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]!='Z'
  except OSError:alive=False
  if alive:running.append(name)
  else:failed.add(name)
 state(state='running_concurrently' if running else ('some_runs_failed' if failed else 'both_completed_waiting_for_user_review'),running=running,failed=sorted(failed),workers_per_study=2,repetitions_per_dataset=10000,regression_imputations_per_repetition=40)
 if not running:break
 time.sleep(15)
