import pandas as pd
import csv

MKAN_PlasGene= pd.read_csv('Mycobacterium_kansasii_ATCC_12478_Plas_Annot.gff_topchopped.tsv', sep='\t',names=['seqname','source','feature','start','end','score','strand','frame','attribute'])
#MKAN_PlasGene[0:3]
  
len(MKAN_PlasGene)
#305

#MKAN_PlasGene.dtypes

MKAN_PlasGene.feature.unique()
#array(['region', 'gene', 'CDS', 'pseudogene'], dtype=object)

#MKAN_PlasGene_Select= MKAN_PlasGene[(MKAN_PlasGene['feature'] != 'CDS') & (MKAN_PlasGene['feature'] != 'region')]
#MKAN_PlasGene_Select[['ID','Additional']] = MKAN_PlasGene_Select.attribute.str.split(";",1,expand=True,)

#MKAN_PlasGene_Select[0:5]
 
#MKAN_PlasGene_Select= MKAN_PlasGene[(MKAN_PlasGene['feature'] != 'CDS') & (MKAN_PlasGene['feature'] != 'region')]
#MKAN_PlasGene_Select[['ID','Dbxref','Name','old_locus_tag','Additional']] = MKAN_PlasGene_Select.attribute.str.split(";",4,expand=True,)
#MKAN_PlasGene_Select[0:5]  

MKAN_PlasGene_Select= MKAN_PlasGene[(MKAN_PlasGene['feature'] != 'CDS') & (MKAN_PlasGene['feature'] != 'region')]
#MKAN_PlasGene_Select[['ID','Dbxref','Name','gbkey','gene_biotype','locus_tag','old_locus_tag','Additional']] = MKAN_PlasGene_Select.attribute.str.split(";",8,expand=True,)

MKAN_PlasGene_Select[['ID','Dbxref','Name','gbkey','gene_biotype','locus_tag','old_locus_tag','Additional']] = MKAN_PlasGene_Select.attribute.str.split(";",7,expand=True,)
#MKAN_PlasGene_Select[0:5]
       
#MKAN_PlasGene_Select.Additional.unique()

MKAN_Plas_TA = pd.read_csv('190313_Plas_Tnhits_Zaid.csv',index_col=False, skiprows=1, names=['Plas_Tnhits'])

#MKAN_Plas_TA_list = MKAN_Plas_TA.values.tolist()
MKAN_Plas_TA_list = MKAN_Plas_TA['Plas_Tnhits'].tolist()
MKAN_Plas_TA_list[0:4]
#[720, 946, 1075, 1226]
Inter_TAs_Plas = []
Intra_TAs_Plas = []

def TA_select_2(entry_a, Genes_Gff):
    MKAN_annot_notCDS_Chr_Ins =Genes_Gff[(Genes_Gff['start'] <= entry_a) & (Genes_Gff['end'] >= entry_a)]
    MKAN_annot_notCDS_Chr_Ins_list = MKAN_annot_notCDS_Chr_Ins.values.tolist()
    flat_list_1 = [item for sublist in MKAN_annot_notCDS_Chr_Ins_list for item in sublist]
    flat_list_1.append(entry_a)
    return(flat_list_1)

def TA_select_intergenic_2(entry_a, Genes_Gff):
    Genes_Gff['endDiff'] = entry_a - Genes_Gff['end']
    Last_End_diff_grtr_0 = Genes_Gff[Genes_Gff['endDiff'] >= 0]
    Last_End_row = Last_End_diff_grtr_0[Last_End_diff_grtr_0.endDiff == Last_End_diff_grtr_0.endDiff.min()]
    Genes_Gff['startDiff'] = Genes_Gff['start'] - entry_a
    Next_Start_diff_grtr_0 = Genes_Gff[Genes_Gff['startDiff'] >= 0]
    Next_Start_row = Next_Start_diff_grtr_0[Next_Start_diff_grtr_0.startDiff == Next_Start_diff_grtr_0.startDiff.min()]
    TA_interGen_rows = pd.merge(Last_End_row,Next_Start_row, on='feature')
    TA_interGen_rows_list = TA_interGen_rows.values.tolist()
    flat_list_1 = [item for sublist in TA_interGen_rows_list for item in sublist]
    flat_list_1.append(entry_a)
    return(flat_list_1)


Intra_TAs_Plas = []
for line in MKAN_Plas_TA_list:
    Intra_TAs_Plas.append(TA_select_2(line, MKAN_PlasGene_Select))


Intra_TAs_Plas[0:19]

len(Intra_TAs_Plas)
#463
Inter_TAs_Plas = []
for line in MKAN_Plas_TA_list:
    Inter_TAs_Plas.append(TA_select_intergenic_2(line, MKAN_PlasGene_Select))


len(Inter_TAs_Plas)
#463

Inter_TAs_Plas[0:10]
#[[720], [946], [1075], [1226], [1783], [2566], [2590], ['NC_022654.1', 'RefSeq', 'gene', 5488, 5736, '.', '+', '.', 'ID=gene-MKAN_RS27780;Dbxref=GeneID:29696832;Name=MKAN_RS27780;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27780;old_locus_tag=MKAN_29095', 'ID=gene-MKAN_RS27780', 'Dbxref=GeneID:29696832', 'Name=MKAN_RS27780', 'gbkey=Gene', 'gene_biotype=protein_coding', 'locus_tag=MKAN_RS27780', 'old_locus_tag=MKAN_29095', None, 102, 2898, 'NC_022654.1', 'RefSeq', 5840, 6730, '.', '+', '.', 'ID=gene-MKAN_RS27785;Dbxref=GeneID:29696805;Name=MKAN_RS27785;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27785;old_locus_tag=MKAN_29100', 'ID=gene-MKAN_RS27785', 'Dbxref=GeneID:29696805', 'Name=MKAN_RS27785', 'gbkey=Gene', 'gene_biotype=protein_coding', 'locus_tag=MKAN_RS27785', 'old_locus_tag=MKAN_29100', None, -892, 2, 5838], ['NC_022654.1', 'RefSeq', 'gene', 5488, 5736, '.', '+', '.', 'ID=gene-MKAN_RS27780;Dbxref=GeneID:29696832;Name=MKAN_RS27780;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27780;old_locus_tag=MKAN_29095', 'ID=gene-MKAN_RS27780', 'Dbxref=GeneID:29696832', 'Name=MKAN_RS27780', 'gbkey=Gene', 'gene_biotype=protein_coding', 'locus_tag=MKAN_RS27780', 'old_locus_tag=MKAN_29095', None, 126, -350, 'NC_022654.1', 'RefSeq', 6957, 7430, '.', '+', '.', 'ID=gene-MKAN_RS27790;Dbxref=GeneID:29696757;Name=MKAN_RS27790;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27790;old_locus_tag=MKAN_29105', 'ID=gene-MKAN_RS27790', 'Dbxref=GeneID:29696757', 'Name=MKAN_RS27790', 'gbkey=Gene', 'gene_biotype=protein_coding', 'locus_tag=MKAN_RS27790', 'old_locus_tag=MKAN_29105', None, -1568, 1095, 5862], ['NC_022654.1', 'RefSeq', 'gene', 5488, 5736, '.', '+', '.', 'ID=gene-MKAN_RS27780;Dbxref=GeneID:29696832;Name=MKAN_RS27780;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27780;old_locus_tag=MKAN_29095', 'ID=gene-MKAN_RS27780', 'Dbxref=GeneID:29696832', 'Name=MKAN_RS27780', 'gbkey=Gene', 'gene_biotype=protein_coding', 'locus_tag=MKAN_RS27780', 'old_locus_tag=MKAN_29095', None, 271, -374, 'NC_022654.1', 'RefSeq', 6957, 7430, '.', '+', '.', 'ID=gene-MKAN_RS27790;Dbxref=GeneID:29696757;Name=MKAN_RS27790;gbkey=Gene;gene_biotype=protein_coding;locus_tag=MKAN_RS27790;old_locus_tag=MKAN_29105', 'ID=gene-MKAN_RS27790', 'Dbxref=GeneID:29696757', 'Name=MKAN_RS27790', 'gbkey=Gene', 'gene_biotype=protein_coding', 'locus_tag=MKAN_RS27790', 'old_locus_tag=MKAN_29105', None, -1423, 950, 6007]]


with open('MKAN_Plas_TA_Intra_4-24-19.csv','w+') as outfile:
    csv.writer(outfile).writerows(Intra_TAs_Plas)

with open('MKAN_Plas_TA_Inter_4-24-19.csv','w+') as outfile:
    csv.writer(outfile).writerows(Inter_TAs_Plas)



Intra_TAs_Plas_II = pd.read_csv('MKAN_Plas_TA_Intra_4-24-19a.csv',index_col=False)
Intra_TAs_Plas_II['TA_frm_Genestrt'] = Intra_TAs_Plas_II['TA_hit'] - Intra_TAs_Plas_II['start']

Intra_TAs_Plas_II[0:3]

Intra_TAs_Plas_II = pd.read_csv('MKAN_Plas_TA_Intra_4-24-19a.csv',index_col=False, names=['seqname','source','feature','start','end','score','strand','frame','attribute','gene1','Dbxref','Name','gbkey','gene_biotype','locus_tag1','old_locus_tag_1','Extra','TA_hit'])

Intra_TAs_Plas_II[0:3] 

Intra_TAs_Plas_II['Gene_Len'] = Intra_TAs_Plas_II['end'] - Intra_TAs_Plas_II['start']
Intra_TAs_Plas_II['TA_frm_Genestrt'] = Intra_TAs_Plas_II['TA_hit'] - Intra_TAs_Plas_II['start']
Intra_TAs_Plas_II['TA_pct_of_Gene'] = Intra_TAs_Plas_II['TA_frm_Genestrt'] / Intra_TAs_Plas_II['Gene_Len']
Intra_TAs_Plas_II[0:3]      

Intra_TAs_Plas_II['end'] = Intra_TAs_Plas_II['end'].astype(int)

Intra_TAs_Plas_II[0:3]



Inter_TAs_Plas_II = pd.read_csv('MKAN_Plas_TA_Inter_4-24-19.csv',index_col=False, names=['seqname1','source1','feature1','start1','end1','score1','strand1','frame1','attribute1','gene1','Dbxref1','Name1','gbkey1','gene_biotype1','locus_tag1','old_locus_tag_1','Extra1','EndDiff_1','StartDiff_1','seqname2','source2','feature2','start2','end2','score2','strand2','frame2','attribute2','gene2','Dbxref2','Name2','gbkey2','gene_biotype2','locus_tag2','old_locus_tag_2','Extra2','EndDiff_2','StartDiff_2','TA_hit'])

list(Inter_TAs_Plas_II.columns.values)
#['seqname1', 'source1', 'feature1', 'start1', 'end1', 'score1', 'strand1', 'frame1', 'attribute1', 'gene1', 'Dbxref1', 'Name1', 'gbkey1', 'gene_biotype1', 'locus_tag1', 'old_locus_tag_1', 'Extra1', 'EndDiff_1', 'StartDiff_1', 'seqname2', 'source2', 'feature2', 'start2', 'end2', 'score2', 'strand2', 'frame2', 'attribute2', 'gene2', 'Dbxref2', 'Name2', 'gbkey2', 'gene_biotype2', 'locus_tag2', 'old_locus_tag_2', 'Extra2', 'EndDiff_2', 'StartDiff_2', 'TA_hit']

Intra_TAs_Plas_III = Intra_TAs_Plas_II[['TA_hit','TA_frm_Genestrt','gene1','old_locus_tag_1', 'start', 'end','strand']]


Intra_TAs_Plas_III = Intra_TAs_Plas_II[['TA_hit','TA_frm_Genestrt','gene1','old_locus_tag_1', 'start', 'end','strand']]
Intra_TAs_Plas_III.rename(columns={'TA_hit' : 'pos' ,'TA_frm_Genestrt' : 'rel_pos_on_gene1', 'start': 'start1', 'end': 'end1','strand' : 'or1'}, inplace=True)


Inter_TAs_Plas_III = Inter_TAs_II[['TA_hit', 'EndDiff_1','gene1', 'old_locus_tag_1', 'start1', 'end1', 'strand1', 'gene2', 'old_locus_tag_2', 'start2', 'end2','strand2']]

Inter_TAs_Plas_IV.rename(columns={'TA_hit' : 'pos' ,'EndDiff_1' : 'rel_pos_on_gene1','strand_1' : 'or1','strand2' : 'or2'}, inplace=True)

Intra_TAs_Plas_III = Intra_TAs_Plas_II[['TA_hit','TA_frm_Genestrt','gene1','old_locus_tag_1', 'start', 'end','strand']]
Intra_TAs_Plas_III.rename(columns={'TA_hit' : 'pos' ,'TA_frm_Genestrt' : 'rel_pos_on_gene1', 'start': 'start1', 'end': 'end1','strand' : 'or1'}, inplace=True)

Inter_TAs_Plas_III = Inter_TAs_Plas_II[['TA_hit', 'EndDiff_1','gene1', 'old_locus_tag_1', 'start1', 'end1', 'strand1', 'gene2', 'old_locus_tag_2', 'start2', 'end2','strand2']]
Inter_TAs_Plas_III.rename(columns={'TA_hit' : 'pos' ,'EndDiff_1' : 'rel_pos_on_gene1','strand_1' : 'or1','strand2' : 'or2'}, inplace=True)

Inter_TAs_Plas_III[0:3]

Inter_TAs_Plas_II = pd.read_csv('MKAN_Plas_TA_Inter_4-24-19.csv',index_col=False, names=['TA_hit','seqname1','source1','feature1','start1','end1','score1','strand1','frame1','attribute1','gene1','Dbxref1','Name1','gbkey1','gene_biotype1','locus_tag1','old_locus_tag_1','Extra1','EndDiff_1','StartDiff_1','seqname2','source2','feature2','start2','end2','score2','strand2','frame2','attribute2','gene2','Dbxref2','Name2','gbkey2','gene_biotype2','locus_tag2','old_locus_tag_2','Extra2','EndDiff_2','StartDiff_2'])

Inter_TAs_Plas_II[0:5]

Inter_TAs_Plas_III = Inter_TAs_Plas_II[['TA_hit', 'EndDiff_1','gene1', 'old_locus_tag_1', 'start1', 'end1', 'strand1', 'gene2', 'old_locus_tag_2', 'start2', 'end2','strand2']]
Inter_TAs_Plas_III.rename(columns={'TA_hit' : 'pos' ,'EndDiff_1' : 'rel_pos_on_gene1','strand_1' : 'or1','strand2' : 'or2'}, inplace=True)
Inter_TAs_Plas_III[0:5]

Inter_TAs_Plas_III[0:19] 

Inter_TAs_Plas_II[8]

Inter_TAs_Plas_II[8:9]
     
Inter_TAs_Plas_II = pd.read_csv('MKAN_Plas_TA_Inter_4-24-19.csv',index_col=False, names=['seqname1','source1','feature1','start1','end1','score1','strand1','frame1','attribute1','gene1','Dbxref1','Name1','gbkey1','gene_biotype1','locus_tag1','old_locus_tag_1','Extra1','EndDiff_1','StartDiff_1','seqname2','source2','feature2','start2','end2','score2','strand2','frame2','attribute2','gene2','Dbxref2','Name2','gbkey2','gene_biotype2','locus_tag2','old_locus_tag_2','Extra2','EndDiff_2','StartDiff_2','TA_hit'])

#Inter_TAs_Plas_II[8:9]

Inter_TAs_Plas_II = pd.read_csv('MKAN_Plas_TA_Inter_4-24-19.csv',index_col=False, names=['seqname1','source1','feature1','start1','end1','score1','strand1','frame1','attribute1','gene1','Dbxref1','Name1','gbkey1','gene_biotype1','locus_tag1','old_locus_tag_1','Extra1','EndDiff_1','StartDiff_1','seqname2','source2','start2','end2','score2','strand2','frame2','attribute2','gene2','Dbxref2','Name2','gbkey2','gene_biotype2','locus_tag2','old_locus_tag_2','Extra2','EndDiff_2','StartDiff_2','TA_hit'])

 
Inter_TAs_Plas_III = Inter_TAs_Plas_II[['TA_hit', 'EndDiff_1','gene1', 'old_locus_tag_1', 'start1', 'end1', 'strand1', 'gene2', 'old_locus_tag_2', 'start2', 'end2','strand2']]
Inter_TAs_Plas_III.rename(columns={'TA_hit' : 'pos' ,'EndDiff_1' : 'rel_pos_on_gene1','strand_1' : 'or1','strand2' : 'or2'}, inplace=True)
#Inter_TAs_Plas_III[0:4]
 
Intra_TAs_Plas_III = Intra_TAs_Plas_II[Intra_TAs_Plas_II.source.notnull()]
Inter_TAs_Plas_III = Inter_TAs_Plas_II[Inter_TAs_Plas_II.source1.notnull()]

len(Intra_TAs_Plas_III)
#375
len(Inter_TAs_Plas_III)
#415
Inter_TAs_Plas_III[0:5]

Intra_TAs_Plas_III['TA_hit'] = Intra_TAs_Plas_III['TA_hit'].astype(int)

Intra_TAs_Plas_III = Intra_TAs_Plas_II[Intra_TAs_Plas_II.source.notnull()]
Inter_TAs_Plas_III = Inter_TAs_Plas_II[Inter_TAs_Plas_II.source1.notnull()]
Intra_TAs_Plas_IV = Intra_TAs_Plas_III[['TA_hit','TA_frm_Genestrt','gene1','old_locus_tag_1', 'start', 'end','strand']]
Intra_TAs_Plas_IV.rename(columns={'TA_hit' : 'pos' ,'TA_frm_Genestrt' : 'rel_pos_on_gene1', 'start': 'start1', 'end': 'end1','strand' : 'or1'}, inplace=True)

Inter_TAs_Plas_IV = Inter_TAs_Plas_III[['TA_hit', 'EndDiff_1','gene1', 'old_locus_tag_1', 'start1', 'end1', 'strand1', 'gene2', 'old_locus_tag_2', 'start2', 'end2','strand2']]
Inter_TAs_Plas_IV.rename(columns={'TA_hit' : 'pos' ,'EndDiff_1' : 'rel_pos_on_gene1','strand_1' : 'or1','strand2' : 'or2'}, inplace=True)

Inter_TAs_Plas_IV[0:4]
 
Intra_TAs_Plas_IV[0:4] 

Intra_TAs_Plas_IV['gene1'] = Intra_TAs_Plas_IV['gene1'].map(lambda x: x.lstrip('ID=gene-'))

Inter_TAs_Plas_IV['gene1'] = Intra_TAs_Plas_IV['gene1'].map(lambda x: x.lstrip('ID=gene-'))

Intra_TAs_Plas_IV['old_locus_tag_1'] = Intra_TAs_Plas_IV['old_locus_tag_1'].astype(str).map(lambda x: x.lstrip('old_locus_tag='))

Intra_TAs_Plas_IV[0:4]

Intra_TAs_Plas_IV['start1'] = Intra_TAs_Plas_IV['start1'].astype(int)
Intra_TAs_Plas_IV['end1'] = Intra_TAs_Plas_IV['end1'].astype(int)
Inter_TAs_Plas_IV[0:4]
 
                                                  ^
Intra_TAs_Plas_IV.to_csv('MKAN_Plas_TA_Intra_4-26-19.csv',header=True, index=False)
Inter_TAs_Plas_IV.to_csv('MKAN_Plas_TA_Inter_4-26-19.csv',header=True, index=False)

##############################

Inter_TAs_Plas_IV = pd.read_csv('MKAN_Plas_TA_Inter_4-26-19.csv',index_col=False)
Intra_TAs_Plas_IV = pd.read_csv('MKAN_Plas_TA_Intra_4-26-19.csv',index_col=False)

MKAN_Plas_TA = pd.read_csv('190313_Plas_Tnhits_Zaid.csv',index_col=False, skiprows=1, names=['Plas_Tnhits'])
MKAN_Plas_TA_list = MKAN_Plas_TA['Plas_Tnhits'].tolist()

Intra_TAs_Plas_IV[['pos','rel_pos_on_gene1']] = Intra_TAs_Plas_IV[['pos','rel_pos_on_gene1']].astype(int)

Intra_pos = Intra_TAs_Plas_IV['pos'].tolist()

Inter_TAs_Plas=[]
for line in Not_Intra_pos:
    Inter_TAs_Plas.append(TA_select_intergenic_2(line, MKAN_PlasGene_Select))





Intra_Plas_Table = pd.read_csv('MKAN_Plas_TA_Inter_FromNotIntras_4-29-19a.csv',index_col=False)

Intra_Plas_Table_NoNaN = Intra_Plas_Table[Intra_Plas_Table.feature_1.notnull()]

Intra_Plas_Table_NoNaN_II = Intra_Plas_Table_NoNaN[['TA_hit', 'endDiff_1','gene1', 'old_locus_tag_1', 'start1', 'end1', 'strand1', 'gene2', 'old_locus_tag_2', 'start2', 'end2','strand2']]
Intra_Plas_Table_NoNaN_II.rename(columns={'TA_hit' : 'pos' ,'EndDiff_1' : 'rel_pos_on_gene1','strand_1' : 'or1','strand2' : 'or2'}, inplace=True)


Intra_Plas_Table_NoNaN_II[['Tahit','endDiff_1','start_1','end_1','start_2','end_2']] = Intra_Plas_Table_NoNaN_II[['Tahit','endDiff_1','start_1','end_1','start_2','end_2']].astype(int)

Intra_Plas_Table_NoNaN_II['locus_tag_1']= Intra_Plas_Table_NoNaN_II['locus_tag_1'].astype(str).map(lambda x: x.lstrip('locus_tag='))
Intra_Plas_Table_NoNaN_II['locus_tag_2']= Intra_Plas_Table_NoNaN_II['locus_tag_2'].astype(str).map(lambda x: x.lstrip('locus_tag='))
Intra_Plas_Table_NoNaN_II['old_locus_tag_1']= Intra_Plas_Table_NoNaN_II['old_locus_tag_1'].astype(str).map(lambda x: x.lstrip('old_locus_tag='))
Intra_Plas_Table_NoNaN_II['old_locus_tag_2']= Intra_Plas_Table_NoNaN_II['old_locus_tag_2'].astype(str).map(lambda x: x.lstrip('old_locus_tag='))

Intra_Plas_Table_NoNaN_II = Intra_Plas_Table_NoNaN_II.rename(columns={'Tahit' : 'pos', 'endDiff_1':'rel_pos_on_gene1', 'locus_tag_1':'gene1' , 'start_1': 'start','end_1': 'end', 'strand_1':'or1','start_2':'start2', 'end_2':'end2', 'locus_tag_2':'gene2', 'strand_2':'or2'})

Plas_result= pd.concat([Intra_TAs_Plas_IV, Intra_Plas_Table_NoNaN_II], ignore_index=True)
Plas_result= result[['pos','rel_pos_on_gene1','gene1','old_locus_tag_1','start1','end1','or1','gene2','old_locus_tag_2','start2','end2','or2']]




MKAN_Chr_Table['gene1']= MKAN_Chr_Table['gene1'].astype(str).map(lambda x: x.rstrip('.pseudogene'))


MKAN_Chr_Table_II = MKAN_Chr_Table


MKAN_Chr_Table_II['old_locus_tag_1'] = MKAN_Chr_Table_II['old_locus_tag_1'].str.replace(r'ID=.*','')

MKAN_Chr_Table_II['old_locus_tag_1'] = MKAN_Chr_Table_II['old_locus_tag_1'].str.replace(r';product=.*','')
MKAN_Chr_Table_II['old_locus_tag_2'] = MKAN_Chr_Table_II['old_locus_tag_2'].str.replace(r';product=.*','')
MKAN_Chr_Table_II['old_locus_tag_2'] = MKAN_Chr_Table_II['old_locus_tag_2'].str.replace(r';pseudo=.*','')
MKAN_Chr_Table_II['old_locus_tag_1'] = MKAN_Chr_Table_II['old_locus_tag_1'].str.replace(r';pseudo=.*','')
