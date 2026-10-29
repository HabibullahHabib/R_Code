from pathlib import Path
import subprocess,time,json,fcntl,os,signal
BASE=Path('/workspace'); state=BASE/'paper-feedback-queue-status.json'
lock=(BASE/'paper-feedback-queue.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def status(**d): state.write_text(json.dumps(d,indent=2)+'\n')
for k in (5,):
 root=BASE/f'paper-feedback-{k}-treatments';full=root/'full'
 while not (root/'smoke/COMPLETED.txt').exists():
  if (root/'smoke/FAILED.txt').exists():status(state='blocked_smoke_validation',treatments=k);raise SystemExit(1)
  status(state='waiting_for_smoke_validation',treatments=k);time.sleep(30)
 if not (root/'FIDELITY_PASSED.txt').exists():
  env=os.environ.copy();env.update(PATH='/workspace/.r-runtime/bin:'+env['PATH'],R_LIBS_SITE='/workspace/.r-runtime/usr/lib/R/site-library',R_LIBS_USER='/workspace/.r-runtime/library',LD_LIBRARY_PATH='/workspace/.r-runtime/usr/lib/x86_64-linux-gnu',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
  with (root/'fidelity_check.log').open('ab') as log: result=subprocess.run(['/workspace/.r-runtime/bin/Rscript','--vanilla',str(root/'verify_fidelity.R'),str(root)],env=env,stdout=log,stderr=subprocess.STDOUT)
  if result.returncode:status(state='blocked_fidelity_validation',treatments=k);raise SystemExit(1)
 if not (full/'COMPLETED.txt').exists():
  status(state='running',treatments=k)
  active=False
  try:
   pid=int((root/'full_run.pid').read_text())
   cmd=Path(f'/proc/{pid}/cmdline').read_bytes()
   active=str(root).encode() in cmd and Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]!='Z'
  except (OSError,ValueError):pass
  if active:
   while Path(f'/proc/{pid}/stat').exists():
    if Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]=='Z':break
    if (full/'COMPLETED.txt').exists():break
    time.sleep(30)
   code=0 if (full/'COMPLETED.txt').exists() else 1
  else:
   (full/'FAILED.txt').unlink(missing_ok=True)
   with (root/'full_console.log').open('ab') as log:
    p=subprocess.Popen([str(root/'run.sh'),'full'],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
    (root/'full_run.pid').write_text(str(p.pid));code=p.wait()
  if code:status(state='run_failed',treatments=k,exit_code=code);raise SystemExit(code)
 with (root/'github_publication.log').open('ab') as log:
  p=subprocess.run(['python',str(root/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT)
 if p.returncode:status(state='publication_failed_compute_continues',treatments=k)
status(state='five_treatment_completed_awaiting_user_review',ten_treatment_autostart=False,equal_variance_remains_paused=True)
