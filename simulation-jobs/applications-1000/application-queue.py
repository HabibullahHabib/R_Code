from pathlib import Path
import os,time,subprocess,json,fcntl
BASE=Path('/workspace');lock=(BASE/'application-queue.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(**d):(BASE/'application-queue-status.json').write_text(json.dumps(d,indent=2)+'\n')
for app in ['ibscovars','actg175']:
 root=BASE/f'application-{app}-1000';full=root/'full'
 while not (root/'smoke/COMPLETED.txt').exists():
  if (root/'smoke/FAILED.txt').exists():state(state='trial_failed',application=app);raise SystemExit(1)
  state(state='waiting_for_trial',application=app);time.sleep(15)
 if not (full/'COMPLETED.txt').exists():
  state(state='running',application=app,total_experiments=1000)
  active=False
  try:
   pid=int((root/'full_run.pid').read_text());active=str(root).encode() in Path(f'/proc/{pid}/cmdline').read_bytes() and Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]!='Z'
  except (OSError,ValueError):pass
  if active:
   while Path(f'/proc/{pid}/stat').exists() and not (full/'COMPLETED.txt').exists():
    if Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]=='Z':break
    time.sleep(30)
   code=0 if (full/'COMPLETED.txt').exists() else 1
  else:
   (full/'FAILED.txt').unlink(missing_ok=True)
   with (root/'full_console.log').open('ab') as log:
    p=subprocess.Popen([str(root/'run.sh'),'full'],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True);(root/'full_run.pid').write_text(str(p.pid));code=p.wait()
  if code:state(state='run_failed',application=app);raise SystemExit(code)
 pub=root/'github_publication_status.json'
 if not pub.exists() or json.loads(pub.read_text()).get('state')!='published_and_verified':
  state(state='publishing',application=app)
  with (root/'github_publication.log').open('ab') as log:p=subprocess.run(['python',str(root/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT)
  if p.returncode:state(state='publication_failed_compute_continues',application=app)
state(state='both_applications_completed_waiting_for_user_review')
