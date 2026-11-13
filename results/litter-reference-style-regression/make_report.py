from pathlib import Path
import csv,json,hashlib,subprocess
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from xml.sax.saxutils import escape
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';OUT=FULL/'Litter_Application_Final';assert (FULL/'COMPLETED.txt').exists()
pdfmetrics.registerFont(TTFont('ThesisSerif','/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf'));pdfmetrics.registerFont(TTFont('ThesisSerifBold','/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf'))
s=getSampleStyleSheet();s.add(ParagraphStyle(name='CaptionLitter',fontName='ThesisSerif',fontSize=12,leading=16,alignment=TA_CENTER,spaceAfter=12));s.add(ParagraphStyle(name='BodyLitter',fontName='ThesisSerif',fontSize=10,leading=14,spaceAfter=9));s.add(ParagraphStyle(name='HeadLitter',fontName='ThesisSerifBold',fontSize=10,leading=13,alignment=TA_CENTER))
def P(x,sty='BodyLitter'):return Paragraph(escape(str(x)).replace('\n','<br/>'),s[sty])
rows=list(csv.DictReader((OUT/'litter_PRIMARY_REG_APPLICATION_TABLE.csv').open()));assert len(rows)==9
methods=['Proposed method with shrinkage','Proposed method without shrinkage','Dunnett-type procedure with naive case deletion']
data=[[P(x,'HeadLitter') for x in ['Dose comparison','Point estimate','95% simultaneous\nconfidence interval','Interval length']]];bands=[]
for method in methods:
 bands.append(len(data));data.append([P(method),'','',''])
 for r in rows:
  if r['method']==method:data.append([r['contrast'],f'{float(r["estimate"]):.3f}',r['CI'],f'{float(r["length"]):.3f}'])
t=Table(data,colWidths=[111,104,173,95]);cmd=[('FONTNAME',(0,0),(-1,-1),'ThesisSerif'),('FONTSIZE',(0,0),(-1,-1),10),('ALIGN',(0,0),(-1,-1),'CENTER'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEABOVE',(0,0),(-1,0),.8,colors.black),('LINEBELOW',(0,0),(-1,0),.6,colors.black),('LINEBELOW',(0,-1),(-1,-1),.8,colors.black)]
for i in bands:cmd += [('SPAN',(0,i),(-1,i)),('ALIGN',(0,i),(-1,i),'LEFT'),('TOPPADDING',(0,i),(-1,i),12)]
t.setStyle(TableStyle(cmd))
inputs=list(csv.DictReader((FULL/'analysis_input_incomplete.csv').open()));missing=sum(x['y'] in ('','NA') for x in inputs)
val=list(csv.DictReader((OUT/'litter_VALIDATION_SUMMARY.csv').open()));assert {x['MI'] for x in val}=={'REG'}
story=[P('Litter: regression-imputation application','CaptionLitter'),P('Simultaneous statistical inferences of treatment-control mean differences','CaptionLitter'),t,Spacer(1,12),P(f'Primary analysis: 40 regression imputations shared by both proposed procedures. Target expected MAR missingness is 30%; realized missingness is {missing}/{len(inputs)} = {100*missing/len(inputs):.2f}%. The covariate is standardized litter size.'),P('Only regression was run. The uploaded script also listed PMM/PS; those comparisons were disabled to honour the user\'s regression-only instruction. The original seeds and REG calculations are retained.'),P('Validation check specified by the script: 500 new MAR repetitions requested. Benchmark recovery is against the fixed complete-sample treatment-control contrasts, not population-parameter coverage.')]
for v in val:story.append(P(f"{v['version']}: {v['successful_replicates']} successful repetitions; joint benchmark recovery {float(v['joint_benchmark_recovery_pct']):.2f}%; mean interval length {float(v['mean_interval_length']):.3f}."))
SimpleDocTemplate(str(FULL/'main_tables.pdf'),pagesize=A4,leftMargin=50,rightMargin=50,topMargin=45,bottomMargin=45).build(story)
for k in ['warnings','errors']:
 p=FULL/(k+'.csv')
 if not p.exists():p.write_text(k[:-1]+'\n')
meta={'dataset':'multcomp::litter','observations':len(inputs),'target_missing_pct':30,'realized_missing':missing,'realized_missing_pct':100*missing/len(inputs),'imputations':40,'imputation_method':'REG','validation_requested':500,'validation_successful':{v['version']:int(v['successful_replicates']) for v in val},'original_script_sha256':hashlib.file_digest((ROOT/'original_script.R').open('rb'),'sha256').hexdigest(),'executed_script_sha256':hashlib.file_digest((ROOT/'executed_script.R').open('rb'),'sha256').hexdigest()}
(FULL/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps(meta,indent=2))
