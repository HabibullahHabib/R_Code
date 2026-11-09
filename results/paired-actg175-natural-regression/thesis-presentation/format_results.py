from pathlib import Path
import csv,hashlib,json,statistics,zipfile
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle,PageBreak
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
pdfmetrics.registerFont(TTFont('ThesisSerif','/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf'))
pdfmetrics.registerFont(TTFont('ThesisSerifBold','/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf'))
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from xml.sax.saxutils import escape
BASE=Path('/workspace');REPO=BASE/'R_Code'
NAMES=[('ACTG175','paired-actg175-natural-regression'),('IBScovars','paired-ibscovars-mar30-regression')]
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='CaptionThesis',fontName='ThesisSerif',fontSize=12,leading=16,alignment=TA_CENTER,spaceAfter=12));styles.add(ParagraphStyle(name='NoteThesis',fontName='ThesisSerif',fontSize=10,leading=14,spaceAfter=9));styles.add(ParagraphStyle(name='HeaderThesis',fontName='ThesisSerifBold',fontSize=10,leading=13,alignment=TA_CENTER))
P=lambda s,style='NoteThesis':Paragraph(escape(str(s)).replace('\n','<br/>'),styles[style])
for app,name in NAMES:
 source=BASE/name/'full/raw_intervals.csv';source_sha=hashlib.file_digest(source.open('rb'),'sha256').hexdigest()
 rows=list(csv.DictReader(source.open()));groups={}
 for row in rows:groups.setdefault((row['variance_method'],row['contrast']),[]).append(row)
 assert {k[0] for k in groups}=={'Shrinkage','Non-shrinkage'}
 contrasts=sorted({r['contrast'] for r in rows});averages=[];example=[]
 for method in ['Non-shrinkage','Shrinkage']:
  for contrast in contrasts:
   g=groups[(method,contrast)];assert len(g)==10000
   avg={'procedure':method,'comparison':contrast,'repetitions':len(g)}
   for k in ['estimate','lower','upper','width']:avg['mean_'+k]=statistics.mean(float(x[k]) for x in g)
   assert abs(avg['mean_upper']-avg['mean_lower']-avg['mean_width'])<1e-9
   averages.append(avg)
   a=[x for x in g if int(x['sim'])==1];assert len(a)==1
   x=a[0];example.append({'procedure':method,'comparison':contrast,'repetition':1,**{k:float(x[k]) for k in ['estimate','lower','upper','width']}})
 dest=REPO/'results'/name/'thesis-presentation';dest.mkdir(parents=True,exist_ok=True)
 for filename,data in [('thesis_style_10000_averages.csv',averages),('thesis_style_40_imputation_example.csv',example)]:
  with (dest/filename).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=data[0]);w.writeheader();w.writerows(data)
 def make_table(data,mean=False):
  headers=['Treatment–control\ncomparison','Mean point estimate' if mean else 'Point estimate','Mean lower / upper\n95% interval limits*' if mean else '95% simultaneous\nconfidence interval','Mean interval\nlength' if mean else 'Interval\nlength']
  out=[[P(h,'HeaderThesis') for h in headers]];bands=[]
  for method,label in [('Non-shrinkage','Proposed procedure without shrinkage'),('Shrinkage','Proposed procedure with shrinkage')]:
   bands.append(len(out));out.append([P(label), '', '', ''])
   for x in data:
    if x['procedure']!=method:continue
    pre='mean_' if mean else ''
    out.append([x['comparison'],f'{x[pre+"estimate"]:.3f}',f'({x[pre+"lower"]:.3f}, {x[pre+"upper"]:.3f})',f'{x[pre+"width"]:.3f}'])
  t=Table(out,colWidths=[126,104,164,89],hAlign='CENTER')
  commands=[('FONTNAME',(0,0),(-1,-1),'ThesisSerif'),('FONTSIZE',(0,0),(-1,-1),11),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,0),10),('BOTTOMPADDING',(0,0),(-1,0),12),('TOPPADDING',(0,1),(-1,-1),7),('BOTTOMPADDING',(0,1),(-1,-1),7),('LINEABOVE',(0,0),(-1,0),.8,colors.black),('LINEBELOW',(0,0),(-1,0),.6,colors.black),('LINEBELOW',(0,-1),(-1,-1),.8,colors.black)]
  for r in bands:commands += [('SPAN',(0,r),(-1,r)),('ALIGN',(0,r),(-1,r),'LEFT'),('TOPPADDING',(0,r),(-1,r),15)]
  t.setStyle(TableStyle(commands));return t
 doc=SimpleDocTemplate(str(dest/'Thesis_Style_Results.pdf'),pagesize=A4,rightMargin=50,leftMargin=50,topMargin=55,bottomMargin=50)
 design='Natural missing values retained; repetitions vary only regression-imputation draws.' if app=='ACTG175' else 'A new 30% expected MAR pattern is imposed in each repetition on the complete dataset.'
 story=[P(f'{app}: results in the layout of thesis Table 4.2','CaptionThesis'),P('Table 1. Descriptive averages of simultaneous inference results over 10,000 repetitions.','CaptionThesis'),make_table(averages,True),Spacer(1,18),P('* The limits above are averages of the 10,000 separately computed interval endpoints. They are descriptive summaries, not a single 95% confidence interval.'),P('Each repetition used 40 regression imputations shared by the shrinkage and non-shrinkage procedures. '+design),P('The table layout follows thesis Table 4.2; only the two procedures requested by the user are shown. No equal-variance or imputation-method comparison is added.'),PageBreak(),P(f'{app}: one actual 40-imputation application analysis','CaptionThesis'),P('Table 2. Simultaneous statistical inferences from the first prespecified repetition.','CaptionThesis'),make_table(example,False),Spacer(1,18),P('This page follows the application meaning of thesis Table 4.2: each displayed interval is an actual simultaneous 95% interval from one 40-imputation analysis. It uses repetition 1, selected by its prespecified seed—not by its results.'),P('Table 2 illustrates the application table; it does not replace the full 10,000-repetition results on page 1. These two types of results must not be interpreted interchangeably.')]
 doc.build(story)
 (dest/'README.md').write_text(f'''# {app}: thesis-style presentation

Reformats the existing completed results without rerunning or changing any analysis.

- [PDF tables](Thesis_Style_Results.pdf): page 1 gives descriptive averages across all 10,000 repetitions; page 2 gives actual simultaneous 95% intervals from the first prespecified 40-imputation repetition, matching the application interpretation of thesis Table 4.2.
- [10,000-repetition averages CSV](thesis_style_10000_averages.csv). Mean lower/upper bounds are averages of separate intervals; they are **not** a single 95% confidence interval.
- [Actual 40-imputation example CSV](thesis_style_40_imputation_example.csv), repetition 1, not selected for its results.

Only shrinkage and non-shrinkage are shown. Regression imputation is fixed throughout. {design}

Original source: ../raw_intervals.csv; SHA256 {source_sha}. All original analyses, figures and logs remain unchanged.
''')
 with zipfile.ZipFile(dest/'Thesis_Style_Results.zip','w',zipfile.ZIP_DEFLATED) as z:
  for f in ['Thesis_Style_Results.pdf','thesis_style_10000_averages.csv','thesis_style_40_imputation_example.csv','README.md']:z.write(dest/f,f)
 manifest={p.name:hashlib.file_digest(p.open('rb'),'sha256').hexdigest() for p in dest.iterdir() if p.is_file() and p.name!='SHA256SUMS.json'};(dest/'SHA256SUMS.json').write_text(json.dumps(manifest,indent=2)+'\n')
 assert hashlib.file_digest(source.open('rb'),'sha256').hexdigest()==source_sha
 print(app,json.dumps(averages))
