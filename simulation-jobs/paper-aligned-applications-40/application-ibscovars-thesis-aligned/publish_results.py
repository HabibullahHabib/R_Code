from pathlib import Path
import subprocess,hashlib,json,shutil,zipfile,tempfile,fcntl
ROOT=Path(__file__).resolve().parent; FULL=ROOT/'full'; REPO=Path('/workspace/R_Code'); REL=Path('results')/ROOT.name
URL='https://github.com/HabibullahHabib/R_Code.git'
def git(*args):return subprocess.check_output(['git','-C',str(REPO),*args],text=True).strip()
def state(**d): (ROOT/'github_publication_status.json').write_text(json.dumps(d,indent=2)+'\n')
try:
 assert (FULL/'COMPLETED.txt').exists()
 lock=Path('/workspace/paper-feedback-publication.lock').open('w');fcntl.flock(lock,fcntl.LOCK_EX)
 assert not git('status','--porcelain'),'Repository has uncommitted work'
 assert git('branch','--show-current')=='main' and git('remote','get-url','origin')==URL
 git('fetch','origin','main');git('merge','--ff-only','origin/main')
 dest=REPO/REL;assert not dest.exists();dest.mkdir(parents=True)
 archive=dest/(ROOT.name+'_Tables_Plots_and_Logs.zip')
 files=[p for p in FULL.rglob('*') if p.is_file() and p.suffix in ('.csv','.pdf','.html','.txt','.log','.md')]
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,str(p.relative_to(FULL)))
  z.write(ROOT/'original_script.R','original_script.R');z.write(ROOT/'run_all.R','run_all.R');z.write(ROOT/'run.sh','run.sh');z.write(ROOT/'export_report.R','export_report.R');z.write(ROOT/'alignment_helpers.R','alignment_helpers.R');z.write(ROOT/'render_report.R','render_report.R');z.write(ROOT/'METHOD_REVIEW.md','METHOD_REVIEW.md');z.write(ROOT/'full_console.log','console_output.log')
 with zipfile.ZipFile(archive) as z:assert z.testzip() is None
 for p in files:
  target=dest/p.relative_to(FULL);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
 shutil.copy2(ROOT/'original_script.R',dest/'original_script.R')
 for p in FULL.rglob('raw_*.csv.gz'):
  if p.stat().st_size<90*1024**2:shutil.copy2(p,dest/p.name)
  else:
   with p.open('rb') as f:
    i=1
    while block:=f.read(90*1024**2):
     (dest/(p.name+f'.part{i:03d}')).write_bytes(block);i+=1
 checks={str(p.relative_to(dest)):hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in dest.rglob('*') if p.is_file()}
 (dest/'SHA256SUMS.json').write_text(json.dumps(checks,indent=2)+'\n')
 (dest/'README.md').write_text(f"""# {ROOT.name}: 40-imputation real-data adaptation

The user requested 40 imputations after review of the uploaded thesis and reference article. Both documents use 20 for their original real-data application; this count difference is explicit. Original dataset applications differ from ACTG175 and IBScovars, so these results are adaptations, not exact replications.

[Full PDF](simulation_report.pdf) · [Tables PDF](main_tables.pdf) · [Plots PDF](simulation_plots.pdf) · [Results ZIP]({archive.name})

[Thesis-primary Regression MI comparison](THESIS_PRIMARY_regression_equal_vs_shrinkage.csv): equal-variance versus unequal-variance shrinkage intervals. [Reference-method Propensity MI comparison](REFERENCE_METHOD_propensity_equal_vs_shrinkage.csv) is reported separately. PMM and the original unshrunk unequal-variance comparison remain supplementary. The unshrunk unequal-variance procedure is not the common-variance reference baseline.

See METHOD_REVIEW.md and METHOD_ALIGNMENT.txt for source sections, formulas, dataset differences and user overrides. No repeated-deletion simulations were added. Both variance constructions share the same 40 imputations and point estimates. Existing verified seeded imputation prefixes were reused when available, retaining original data, seeds and imputation functions. Session information, logs, warnings/errors, detailed CSVs and completion counts are included.
""")
 readme=REPO/'README.md';t=readme.read_text();t+=f'\n[Professor feedback: {ROOT.name}]({REL}/README.md).\n';readme.write_text(t)
 git('add','--','README.md',str(REL));git('-c','user.name=Codex','-c','user.email=codex@openai.com','commit','-m',f'Add completed {ROOT.name} simulation results')
 commit=git('rev-parse','HEAD');state(state='pushing',commit=commit);git('push','origin','main')
 remote=git('ls-remote','origin','refs/heads/main').split()[0];assert remote==commit
 # Fetch back and verify every published result against its local checksum.
 git('fetch','origin','main')
 for name,sha in checks.items():
  proc=subprocess.Popen(['git','-C',str(REPO),'show',f'FETCH_HEAD:{REL}/{name}'],stdout=subprocess.PIPE)
  assert hashlib.file_digest(proc.stdout,'sha256').hexdigest()==sha;assert proc.wait()==0
 state(state='published_and_verified',commit=commit,url=f'https://github.com/HabibullahHabib/R_Code/tree/main/{REL}')
except Exception as e:
 state(state='publication_failed',error=str(e));raise
