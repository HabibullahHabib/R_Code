from pathlib import Path
import os,time,subprocess,json,fcntl
BASE=Path('/workspace');lock=(BASE/'corrected-application-queue.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
apps=['ibscovars','actg175'];children={};published=set();failed=set()
def state(**d):(BASE/'corrected-application-queue-status.json').write_text(json.dumps(d,indent=2)+'\n')
for app in apps:
 root=BASE/f'application-{app}-corrected';full=root/'full';assert (root/'as-written/COMPLETED.txt').exists()
 if (full/'COMPLETED.txt').exists():continue
 active=False
 try:
  pid=int((root/'full_run.pid').read_text());active=str(root).encode() in Path(f'/proc/{pid}/cmdline').read_bytes() and Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0] not in ('Z','T','t')
 except (OSError,ValueError):pass
 if active:children[app]=(pid,None);continue
 (full/'FAILED.txt').unlink(missing_ok=True)
 env=os.environ.copy();env['SIM_WORKERS']='2'
 with (root/'full_console.log').open('ab') as log:
  p=subprocess.Popen([str(root/'run_10000.sh'),'full'],env=env,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
 (root/'full_run.pid').write_text(str(p.pid));children[app]=(p.pid,p)
while True:
 running=[]
 for app in apps:
  root=BASE/f'application-{app}-corrected';full=root/'full'
  if (full/'COMPLETED.txt').exists():
   if app not in published:
    pub=root/'github_publication_status.json'
    if not pub.exists() or json.loads(pub.read_text()).get('state')!='published_and_verified':
     with (root/'github_publication.log').open('ab') as log:r=subprocess.run(['python',str(root/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT)
    published.add(app)
   continue
  if app in failed:continue
  pid,child=children[app]
  try:alive=Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]!='Z'
  except OSError:alive=False
  if alive:running.append(app)
  else:failed.add(app)
 state(state='running_concurrently' if running else ('run_failed' if failed else 'both_applications_completed_waiting_for_user_review'),running=running,failed=sorted(failed),workers_per_application=2,imputations_per_method_per_analysis=10000)
 if not running:break
 time.sleep(15)
