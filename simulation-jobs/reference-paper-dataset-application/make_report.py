from pathlib import Path
import csv,json,hashlib,subprocess,shutil
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from xml.sax.saxutils import escape
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';OUT=FULL/'ReferencePaper_Data_FINAL_Application'
pdfmetrics.registerFont(TTFont('ThesisSerif','/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf'));pdfmetrics.registerFont(TTFont('ThesisSerifBold','/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf'))
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='CaptionRef',fontName='ThesisSerif',fontSize=12,leading=16,alignment=TA_CENTER,spaceAfter=13));styles.add(ParagraphStyle(name='TextRef',fontName='ThesisSerif',fontSize=10,leading=14,spaceAfter=9));styles.add(ParagraphStyle(name='HeaderRef',fontName='ThesisSerifBold',fontSize=10,leading=13,alignment=TA_CENTER))
def P(s,style='TextRef'):return Paragraph(escape(str(s)).replace('\n','<br/>'),styles[style])
rows=list(csv.DictReader((FULL/'computed_results_full_precision.csv').open()));assert len(rows)==8
assert {r['Procedure'] for r in rows}=={'With shrinkage','Without shrinkage'}
data=[[P('Treatment–control\ncomparison','HeaderRef'),P('Point estimate\n(treatment − control)','HeaderRef'),P('95% simultaneous\nconfidence interval','HeaderRef'),P('Interval\nlength','HeaderRef')]];bands=[]
for method,label in [('Without shrinkage','Proposed procedure without shrinkage'),('With shrinkage','Proposed procedure with shrinkage')]:
 bands.append(len(data));data.append([P(label),'','',''])
 for r in rows:
  if r['Procedure']==method:data.append([r['Comparison'],f'{float(r["Estimate"]):.3f}',f'({float(r["Lower"]):.3f}, {float(r["Upper"]):.3f})',f'{float(r["Length"]):.3f}'])
t=Table(data,colWidths=[122,112,164,85],hAlign='CENTER');cmd=[('FONTNAME',(0,0),(-1,-1),'ThesisSerif'),('FONTSIZE',(0,0),(-1,-1),11),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),8),('BOTTOMPADDING',(0,0),(-1,-1),8),('LINEABOVE',(0,0),(-1,0),.8,colors.black),('LINEBELOW',(0,0),(-1,0),.6,colors.black),('LINEBELOW',(0,-1),(-1,-1),.8,colors.black)]
for r in bands:cmd += [('SPAN',(0,r),(-1,r)),('ALIGN',(0,r),(-1,r),'LEFT'),('TOPPADDING',(0,r),(-1,r),15)]
t.setStyle(TableStyle(cmd))
story=[P('Reference-paper dose-response dataset','CaptionRef'),P('Simultaneous statistical inferences of treatment–control mean differences','CaptionRef'),t,Spacer(1,18),P('Original 100-patient dataset; five arms of 20 patients; 32 original missing responses retained. All responses and missing-value locations were verified against published Table 5.'),P('One application analysis with 40 treatment-restricted approximate Bayesian bootstrap (ABB) imputations. Shrinkage and non-shrinkage share the same imputations and point estimates. These are actual simultaneous 95% intervals, not averages over repeated simulations.'),P('This script uses within-dose ABB rather than regression because the published dataset contains no auxiliary covariate. The original reference paper reported 20 propensity-score imputations. Its reported Table 6 values are stored separately and are not results computed by this script.'),P('Conclusion: dose 0.05 is not significant; doses 0.20, 0.60 and 1.00 are significant under both procedures. Shrinkage shortens the interval only for dose 0.20; the other intervals become longer.')]
SimpleDocTemplate(str(FULL/'main_tables.pdf'),pagesize=A4,leftMargin=50,rightMargin=50,topMargin=50,bottomMargin=45).build(story)
shutil.copy2(OUT/'ReferencePaper_PRIMARY_FOREST_PLOT.pdf',FULL/'comparison_graph.pdf');shutil.copy2(OUT/'ReferencePaper_PRIMARY_FOREST_PLOT.png',FULL/'comparison_graph.png')
subprocess.run(['pdfunite',str(FULL/'main_tables.pdf'),str(FULL/'comparison_graph.pdf'),str(FULL/'simulation_report.pdf')],check=True)
for kind,header in [('warnings','warning'),('errors','error')]:
 p=FULL/(kind+'.csv')
 if not p.exists():p.write_text(header+'\n')
shutil.copy2(ROOT/'DATA_VERIFICATION.json',FULL/'DATA_VERIFICATION.json')
metadata={'dataset':'Jiang, Zhou and Liang (2023), Table 5 original incomplete data','patients':100,'arms':5,'patients_per_arm':20,'original_missing_responses':32,'imputations':40,'imputation_method':'Treatment-restricted approximate Bayesian bootstrap','repeated_simulation_experiments':0,'comparison':['Without shrinkage','With shrinkage'],'original_script_sha256':hashlib.file_digest((ROOT/'original_script.R').open('rb'),'sha256').hexdigest()}
(FULL/'run_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
print('Report generated:',FULL/'simulation_report.pdf')
