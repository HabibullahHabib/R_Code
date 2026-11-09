from pathlib import Path
import csv,json,hashlib,shutil,subprocess
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,Table,TableStyle
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from reportlab.lib import colors
from xml.sax.saxutils import escape
ROOT=Path(__file__).resolve().parent;FULL=ROOT/'full';OUT=FULL/'IBScovars_ReferenceStyle_Output'
assert (FULL/'COMPLETED.txt').exists()
pdfmetrics.registerFont(TTFont('ThesisSerif','/usr/share/fonts/truetype/liberation/LiberationSerif-Regular.ttf'));pdfmetrics.registerFont(TTFont('ThesisSerifBold','/usr/share/fonts/truetype/liberation/LiberationSerif-Bold.ttf'))
styles=getSampleStyleSheet();styles.add(ParagraphStyle(name='CaptionIBS',fontName='ThesisSerif',fontSize=12,leading=16,alignment=TA_CENTER,spaceAfter=12));styles.add(ParagraphStyle(name='TextIBS',fontName='ThesisSerif',fontSize=10,leading=14,spaceAfter=9));styles.add(ParagraphStyle(name='HeaderIBS',fontName='ThesisSerifBold',fontSize=10,leading=13,alignment=TA_CENTER))
def P(s,style='TextIBS'):return Paragraph(escape(str(s)).replace('\n','<br/>'),styles[style])
rows=list(csv.DictReader((OUT/'IBScovars_ReferenceStyle_Table.csv').open()));assert len(rows)==12
methods=['Proposed method with shrinkage','Proposed method without shrinkage',"Dunnett's procedure with naive case deletion"]
data=[[P('Methods / dose comparison','HeaderIBS'),P('Point estimate\n(treatment − control)','HeaderIBS'),P('95% simultaneous\nconfidence interval','HeaderIBS')]];bands=[]
for method in methods:
 bands.append(len(data));data.append([P(method),'',''])
 selected=[x for x in rows if x['method']==method];assert len(selected)==4
 for x in selected:data.append([x['dose_comparison'],f'{float(x["estimate"]):.3f}',x['CI']])
t=Table(data,colWidths=[175,130,178],hAlign='CENTER');cmd=[('FONTNAME',(0,0),(-1,-1),'ThesisSerif'),('FONTSIZE',(0,0),(-1,-1),11),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),6),('BOTTOMPADDING',(0,0),(-1,-1),6),('LINEABOVE',(0,0),(-1,0),.8,colors.black),('LINEBELOW',(0,0),(-1,0),.6,colors.black),('LINEBELOW',(0,-1),(-1,-1),.8,colors.black)]
for r in bands:cmd += [('SPAN',(0,r),(-1,r)),('ALIGN',(0,r),(-1,r),'LEFT'),('TOPPADDING',(0,r),(-1,r),13)]
t.setStyle(TableStyle(cmd))
inputs=list(csv.DictReader((FULL/'analysis_input_incomplete.csv').open()));missing=sum(x['y'] in ('NA','') for x in inputs);achieved=100*missing/len(inputs)
story=[P('IBScovars: reference-paper-style application','CaptionIBS'),P('Simultaneous statistical inferences of treatment–control mean differences','CaptionIBS'),t,Spacer(1,15),P(f'One application analysis; 40 regression imputations as specified by the uploaded script. The target missingness is 30% MAR; the realized missingness is {missing}/{len(inputs)} = {achieved:.2f}%.'),P('Shrinkage and non-shrinkage share the same 40 imputations and point estimates. Naive case deletion is the comparator explicitly included in the supplied script. No additional simulation loop or imputation-method comparison was added.'),P('The accompanying forest plot reproduces the script’s reference-paper layout: 40 individual completed-dataset results, both proposed procedures after Rubin combination of the 40 imputations, and naive case deletion, separately for each dose contrast.')]
SimpleDocTemplate(str(FULL/'main_tables.pdf'),pagesize=A4,leftMargin=50,rightMargin=50,topMargin=45,bottomMargin=45).build(story)
shutil.copy2(OUT/'IBScovars_ReferenceStyle_ForestPlot.pdf',FULL/'comparison_graph.pdf');shutil.copy2(OUT/'IBScovars_ReferenceStyle_ForestPlot.png',FULL/'comparison_graph.png')
subprocess.run(['pdfunite',str(FULL/'main_tables.pdf'),str(FULL/'comparison_graph.pdf'),str(FULL/'simulation_report.pdf')],check=True)
for kind,header in [('warnings','warning'),('errors','error')]:
 p=FULL/(kind+'.csv')
 if not p.exists():p.write_text(header+'\n')
meta={'dataset':'DoseFinding::IBScovars','patients':len(inputs),'target_missing_pct':30,'realized_missing':missing,'realized_missing_pct':achieved,'imputations':40,'imputation_method':'regression','repeated_simulation_experiments':0,'table_methods':methods,'individual_imputation_intervals':160,'forest_plot_rows':172,'executed_script_sha256':hashlib.file_digest((ROOT/'executed_script.R').open('rb'),'sha256').hexdigest(),'original_script_sha256':hashlib.file_digest((ROOT/'original_script.R').open('rb'),'sha256').hexdigest()}
(FULL/'run_metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
print(json.dumps(meta,indent=2))
