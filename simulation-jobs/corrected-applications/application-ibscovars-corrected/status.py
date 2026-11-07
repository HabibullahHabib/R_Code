from pathlib import Path
import csv,json,os
root=Path('/workspace/application-ibscovars-corrected');full=root/'full'
s={'completed':(full/'COMPLETED.txt').exists(),'failed':(full/'FAILED.txt').exists(),'total_jobs':120000}
try:
 p=int((root/'full_run.pid').read_text());os.kill(p,0)
 st=Path(f'/proc/{p}/stat').read_text().split(') ',1)[1].split()[0]
 s['process_running']=st not in ('Z','T','t');s['paused']=st in ('T','t');s['pid']=p
except (OSError,ValueError):s['process_running']=False
files=list((full/'checkpoints').glob('*/batch_*.rds'))
s['saved_jobs']=sum(int(p.read_text().strip()) for p in (full/'checkpoints').glob('*/batch_*.count'))
if (full/'progress.csv').exists():s['progress']=next(csv.DictReader((full/'progress.csv').open()))
for kind in ('warnings','errors'):
 p=full/(kind+'.csv');s[kind+'_recorded']=sum(1 for _ in csv.DictReader(p.open())) if p.exists() else 0
p=root/'github_publication_status.json'
if p.exists():
 s['github_publication']=json.loads(p.read_text()).get('state')
 try:
  pid=int((root/'github_publication.pid').read_text())
  st=Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0]
  cmd=Path(f'/proc/{pid}/cmdline').read_bytes()
  s['publication_watcher_active']=st!='Z' and b'publish_results.py' in cmd
 except (OSError,ValueError):s['publication_watcher_active']=False
print(json.dumps(s,indent=2))
