from pathlib import Path
import subprocess,time,json,fcntl,os
BASE=Path('/workspace'); state=BASE/'paper-feedback-queue-status.json'
lock=(BASE/'paper-feedback-queue.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def status(**d): state.write_text(json.dumps(d,indent=2)+'\n')
while not (BASE/'sim1-k5-equal-variance/full/COMPLETED.txt').exists():
 status(state='waiting_for_equal_variance_run',order=[5,10]);time.sleep(30)
for k in (5,10):
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
  (full/'FAILED.txt').unlink(missing_ok=True)
  with (root/'full_console.log').open('ab') as log:
   p=subprocess.Popen([str(root/'run.sh'),'full'],stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
   (root/'full_run.pid').write_text(str(p.pid));code=p.wait()
  if code:status(state='run_failed',treatments=k,exit_code=code);raise SystemExit(code)
 # Wait for earlier publication so concurrent git changes do not overlap.
 if k==5:
  pub=BASE/'sim1-k5-equal-variance/github_publication_status.json'
  while pub.exists() and json.loads(pub.read_text()).get('state') in ('waiting_for_completed_results','pushing'):time.sleep(30)
 with (root/'github_publication.log').open('ab') as log:
  p=subprocess.run(['python',str(root/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT)
 if p.returncode:status(state='publication_failed',treatments=k);raise SystemExit(p.returncode)
status(state='all_completed_and_published')
