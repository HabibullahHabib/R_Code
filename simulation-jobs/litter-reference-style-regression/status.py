from pathlib import Path
import json,csv
r=Path(__file__).resolve().parent;f=r/'full';s={'completed':(f/'COMPLETED.txt').exists(),'failed':(f/'FAILED.txt').exists(),'imputation_method':'REG','imputations':40,'validation_requested':500}
s['saved_validation_repetitions']=sum(int(p.read_text()) for p in (f/'checkpoints').glob('batch_*.count'))
try:
 pid=int((r/'full_run.pid').read_text());s['process_running']=Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0] not in ('Z','T','t')
except (OSError,ValueError):s['process_running']=False
p=f/'progress.csv'
if p.exists():s['progress']=next(csv.DictReader(p.open()))
p=r/'github_publication_status.json'
if p.exists():s['publication']=json.loads(p.read_text())
print(json.dumps(s,indent=2))
