from pathlib import Path
import csv,json
root=Path('/workspace/paired-actg175-natural-regression');f=root/'full';s={'completed':(f/'COMPLETED.txt').exists(),'failed':(f/'FAILED.txt').exists(),'total_repetitions':10000,'imputations_per_repetition':40}
s['saved_repetitions']=sum(int(p.read_text()) for p in (f/'checkpoints').glob('batch_*.count'))
try:
 pid=int((root/'full_run.pid').read_text());s['process_running']=Path(f'/proc/{pid}/stat').read_text().split(') ',1)[1].split()[0] not in ('Z','T','t');s['pid']=pid
except (OSError,ValueError):s['process_running']=False
p=f/'progress.csv'
if p.exists():s['progress']=next(csv.DictReader(p.open()))
p=root/'github_publication_status.json'
if p.exists():s['publication']=json.loads(p.read_text())
print(json.dumps(s,indent=2))
