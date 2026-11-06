from pathlib import Path
import os,fcntl,subprocess,json,time
BASE=Path('/workspace');lock=(BASE/'sample-size-queue.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
def state(**d): (BASE/'sample-size-queue-status.json').write_text(json.dumps(d,indent=2)+'\n')
for k in (5,10):
 root=BASE/f'sample-size-{k}-treatments-2000';full=root/'full'
 assert (root/'smoke/COMPLETED.txt').exists()
 if not (root/'FIDELITY_PASSED.txt').exists():
  state(state='verifying_original_calculations',treatments=k)
  env=os.environ.copy();env.update(PATH='/workspace/.r-runtime/bin:'+env['PATH'],R_LIBS_SITE='/workspace/.r-runtime/usr/lib/R/site-library',R_LIBS_USER='/workspace/.r-runtime/library',LD_LIBRARY_PATH='/workspace/.r-runtime/usr/lib/x86_64-linux-gnu',OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
  with (root/'fidelity_check.log').open('ab') as log:r=subprocess.run(['/workspace/.r-runtime/bin/Rscript','--vanilla',str(root/'verify_fidelity.R'),str(root)],env=env,stdout=log,stderr=subprocess.STDOUT)
  if r.returncode:state(state='fidelity_failed',treatments=k);raise SystemExit(1)
 if not (full/'COMPLETED.txt').exists():
  state(state='running',treatments=k,total_jobs=2000,next_treatments=10 if k==5 else None)
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
  if code:state(state='run_failed',treatments=k);raise SystemExit(code)
 pub=root/'github_publication_status.json'
 if not pub.exists() or json.loads(pub.read_text()).get('state')!='published_and_verified':
  state(state='publishing',treatments=k)
  with (root/'github_publication.log').open('ab') as log:r=subprocess.run(['python',str(root/'publish_results.py')],stdout=log,stderr=subprocess.STDOUT)
  if r.returncode:state(state='publication_failed_computation_continues',treatments=k)
state(state='both_sample_size_runs_completed')
