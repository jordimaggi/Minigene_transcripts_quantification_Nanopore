import subprocess
import datetime
import os
import glob
import sys
import pandas as pd
import openpyxl
import re

with open("sample_IDs.txt", "r") as f:
	samples = [x.rstrip() for x in f]

c = 0

for sample in samples:
    print('Processing results for {}'.format(sample))
    gff3 = '/Reference_sequence/RHO_minigene_{}_construct.gff3'.format(sample, sample)
    regions = pd.read_csv(gff3, sep='\t')
    regions = regions.drop(['##gff-version 3.1.26', 'Unnamed: 1', 'Unnamed: 5', 'Unnamed: 7'], axis=1)
    regions = regions.rename(columns={'Unnamed: 2': 'type', 'Unnamed: 3': 'start', 'Unnamed: 4': 'end', 'Unnamed: 6': 'strand', 'Unnamed: 8': 'annotations'})
    for i in range(regions.shape[0]):
        regions.loc[i,'annotations'] = regions.iloc[i].at['annotations'].split('ID=')[1].split(';')[0]

    exons = regions[regions['type'] == 'exon'].reset_index(drop=True)
    exons['length'] = exons['end'] - exons['start']
    exons_dict = {'start': [], 'end': [], 'name': [], 'length': []}
    for i in range(exons.shape[0]):
        exons_dict['start'].append(exons.iloc[i].at['start'])
        exons_dict['end'].append(exons.iloc[i].at['end'])
        exons_dict['name'].append(exons.iloc[i].at['annotations'])
        exons_dict['length'].append(exons.iloc[i].at['length'])

    SJs = {'donor': [], 'acceptor': [], 'junction': []}
    for i in range(exons.shape[0]):
            if len(SJs["donor"]) == 0:
                SJs['donor'].append(exons.iloc[i].at['end'])
            else:
                SJs['acceptor'].append(exons.iloc[i].at['start'])
                SJs['junction'].append(exons.iloc[i-1].at['annotations'] + '-' + exons.iloc[i].at['annotations'])
                if i+1 < exons.shape[0]:
                    SJs['donor'].append(exons.iloc[i].at['end'])
            if exons.iloc[i].at['annotations'] == 'RHO_ex3':
                exons.loc[i, 'length'] = 38
            elif exons.iloc[i].at['annotations'] == 'RHO_ex5':
                exons.loc[i, 'length'] = 81
            else:
                exons.loc[i, 'length'] = exons.iloc[i].at['end'] - exons.iloc[i].at['start']
    SJ_WT = pd.DataFrame.from_dict(SJs)
    SJs_new = {'donor': [], 'acceptor': [], 'sj': [], 'junction': [], 'WT_count': [], 'MT_count': []}
    for i in range(SJ_WT.shape[0]):
        donor = SJ_WT.iloc[i].at['donor']
        don = SJ_WT.iloc[i].at['junction'].split('-')[0]
        for i2 in range(SJ_WT.shape[0]):
            acceptor = SJ_WT.iloc[i2].at['acceptor']
            acc = SJ_WT.iloc[i2].at['junction'].split('-')[1]
            if acceptor > donor:
                SJs_new['donor'].append(donor)
                SJs_new['acceptor'].append(acceptor)
                SJs_new['junction'].append(don + '-' + acc)
                SJs_new['WT_count'].append(0)
                SJs_new['MT_count'].append(0)
                SJs_new['sj'].append('(' + str(donor) + ', ' + str(acceptor) + ')')

    SJ_all = pd.DataFrame.from_dict(SJs_new)


    file_WT = '/{}/NP_{}_WT_cDNA_NanoSplicer.hdf5.csv'.format(sample, sample)
    df_WT = pd.read_csv(file_WT, sep=',')
    df_WT = df_WT.drop(['Unnamed: 0', 'chrID'], axis=1)
    df_WT = df_WT.rename(columns={'loc': 'sj', 'Unnamed: 6': 'strand'})

    WT_sj = df_WT.loc[df_WT['JAQ'] >= 1].reset_index(drop=True)
    trans = {'SJs': [], 'mean_JAQ': []}
    old_ID = ''
    for i in range(WT_sj.shape[0]):
        ID = WT_sj.iloc[i].at['id']
        if ID == old_ID:
            continue
        sj_list = df_WT.loc[df_WT['id']==ID]['sj'].tolist()
        sj_s = "; ".join(sj_list)
        jaq = df_WT.loc[df_WT['id']==ID]['JAQ'].mean()
        trans['SJs'].append(sj_s)
        trans['mean_JAQ'].append(jaq)
        old_ID = ID

    trans_df = pd.DataFrame.from_dict(trans)
    trans_df = trans_df[trans_df['SJs'].str.contains(re.escape('(' + str(SJ_all.iloc[0].at['donor']) + ', ')) & trans_df['SJs'].str.contains(re.escape(str(SJ_all.iloc[-1].at['acceptor']) + ')' ))]
    read_number = trans_df.shape[0]

    filt_trans_WT = trans_df.value_counts(['SJs']).rename_axis(['SJs']).reset_index(name='counts')
    for i in range(filt_trans_WT.shape[0]):
        sj = filt_trans_WT.iloc[i].at['SJs']
        jaq = trans_df.loc[trans_df['SJs']==sj]['mean_JAQ'].mean()
        filt_trans_WT.loc[i, 'mean_JAQ'] = round(jaq, 2)
    filt_trans_WT = filt_trans_WT.loc[filt_trans_WT['counts'] >= 0.005*read_number].reset_index(drop=True)
    filt_trans_WT['percentage'] = round(filt_trans_WT['counts'] / read_number * 100, 2)
    filt_trans_WT['transcript'] = ''
    filt_trans_WT['exons'] = ''
    filt_trans_WT['length'] = 0
    for i in range(filt_trans_WT.shape[0]):
        exs_dict = {'start': [], 'end': []}
        junc = filt_trans_WT.iloc[i].at['SJs'].split('; ')
        for el in junc:
            start = el.split(', ')[1].split(')')[0]
            end = el.split(', ')[0].split('(')[1]
            if filt_trans_WT.loc[i, 'exons'] == '':
                filt_trans_WT.loc[i, 'exons'] += '({}, '.format(exons_dict['start'][0]) + end + '); (' + start + ', '
                exs_dict['start'].append(exons_dict['start'][0])
                exs_dict['end'].append(int(end))
                exs_dict['start'].append(int(start))
                if el == junc[-1]:
                    filt_trans_WT.loc[i, 'exons'] += str(exons.iloc[-1].at['end']) + ')'
                    exs_dict['end'].append(exons.iloc[-1].at['end'])
            else:
                filt_trans_WT.loc[i, 'exons'] += end + '); (' + start + ', '
                exs_dict['start'].append(int(start))
                exs_dict['end'].append(int(end))
                if el == junc[-1]:
                    filt_trans_WT.loc[i, 'exons'] += str(exons.iloc[-1].at['end']) + ')'
                    exs_dict['end'].append(exons.iloc[-1].at['end'])
        for i2 in range(len(exs_dict['start'])):
            st = exs_dict['start'][i2]
            en = exs_dict['end'][i2]
            leng = en -st
            filt_trans_WT.loc[i, 'length'] += leng
            if exs_dict['start'][i2] in exons_dict['start']:
                index = exons_dict['start'].index(exs_dict['start'][i2])
                if exs_dict['end'][i2] == exons_dict['end'][index]:
                    if filt_trans_WT.loc[i, 'transcript'] == '':
                        filt_trans_WT.loc[i, 'transcript'] += exons_dict['name'][index]
                    else:
                        filt_trans_WT.loc[i, 'transcript'] += '-' + exons_dict['name'][index]
                else:
                    if filt_trans_WT.loc[i, 'transcript'] == '':
                        filt_trans_WT.loc[i, 'transcript'] += 'unclear'
                    else:
                        filt_trans_WT.loc[i, 'transcript'] += '-unclear'
            else:
                if filt_trans_WT.loc[i, 'transcript'] == '':
                    filt_trans_WT.loc[i, 'transcript'] += 'unclear'
                else:
                    filt_trans_WT.loc[i, 'transcript'] += '-unclear'


    try:
        filt_trans_WT = filt_trans_WT[['transcript', 'length', 'counts', 'percentage', 'mean_JAQ', 'exons', 'SJs']]
    except:
        print(filt_trans_WT)

    gff3_dict = {'seqid': [], 'source': [], 'type': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []}
    seqid = gff3.split('/')[-1].split('.')[0]
    source = '.'

    for i in range(filt_trans_WT.shape[0]):
        #Specify transcript-level data
        gff3_dict['seqid'].append(seqid)
        gff3_dict['source'].append(source)
        gff3_dict['type'].append('mRNA')
        gff3_dict['start'].append(exons.iloc[0].at['start'])
        gff3_dict['end'].append(exons.iloc[-1].at['end'])
        gff3_dict['score'].append(filt_trans_WT.iloc[i].at['mean_JAQ'])
        gff3_dict['strand'].append('+')
        gff3_dict['phase'].append('.')
        if filt_trans_WT.iloc[i].at['transcript'].split('-') != exons.annotations.to_list():
            gff3_dict['attributes'].append('ID=transcript'+str(i+1)+';'+'Name='+filt_trans_WT.iloc[i].at['transcript']+';Note='+str(filt_trans_WT.iloc[i].at['percentage'])+';color=#ff0000')
        else:
            gff3_dict['attributes'].append('ID=transcript'+str(i+1)+';'+'Name='+filt_trans_WT.iloc[i].at['transcript']+';Note='+str(filt_trans_WT.iloc[i].at['percentage'])+';color=#00b300')


        #Add exons-level data
        exs = filt_trans_WT.iloc[i].at['exons'].split('; ')
        for i2 in range(len(exs)):
            start = exs[i2].split(', ')[0].split('(')[1]
            end = exs[i2].split(', ')[1].split(')')[0]
            gff3_dict['seqid'].append(seqid)
            gff3_dict['source'].append(source)
            gff3_dict['type'].append('exon')
            gff3_dict['start'].append(start)
            gff3_dict['end'].append(end)
            gff3_dict['score'].append(filt_trans_WT.iloc[i].at['mean_JAQ'])
            gff3_dict['strand'].append('+')
            gff3_dict['phase'].append('.')
            try:
                gff3_dict['attributes'].append('ID='+filt_trans_WT.iloc[i].at['transcript'].split('-')[i2]+';Name='+filt_trans_WT.iloc[i].at['transcript'].split('-')[i2]+';Parent=transcript'+str(i+1))
            except:
                print(filt_trans_WT.iloc[i].at['transcript'])
                gff3_dict['attributes'].append('ID='+str(i2))


    try:
        gff3_WT = pd.DataFrame.from_dict(gff3_dict)
    except:
        print(gff3_dict)
    file_line = """##gff-version 3.1.26\
    {}"""
    with open('/{}/{}_WT_transcripts.gff3'.format(sample, sample), 'w') as fp:
        fp.write(file_line.format(gff3_WT.to_csv(sep='\t', index=False)))


    file_MT = '/{}/NP_{}_MT_cDNA_NanoSplicer.hdf5.csv'.format(sample, sample)
    df_MT = pd.read_csv(file_MT, sep=',')
    df_MT = df_MT.drop(['Unnamed: 0', 'chrID'], axis=1)
    df_MT = df_MT.rename(columns={'loc': 'sj', 'Unnamed: 6': 'strand'})
    MT_sj = df_MT.loc[df_MT['JAQ'] >= 1].reset_index(drop=True)
    trans = {'SJs': [], 'mean_JAQ': []}
    old_ID = ''
    for i in range(MT_sj.shape[0]):
        ID = MT_sj.iloc[i].at['id']
        if ID == old_ID:
            continue
        sj_list = df_MT.loc[df_MT['id']==ID]['sj'].tolist()
        sj_s = "; ".join(sj_list)
        jaq = df_MT.loc[df_MT['id']==ID]['JAQ'].mean()
        trans['SJs'].append(sj_s)
        trans['mean_JAQ'].append(jaq)
        old_ID = ID

    trans_df = pd.DataFrame.from_dict(trans)
    trans_df = trans_df[trans_df['SJs'].str.contains(re.escape('(' + str(SJ_all.iloc[0].at['donor']) + ', ')) & trans_df['SJs'].str.contains(re.escape(str(SJ_all.iloc[-1].at['acceptor']) + ')' ))]
    read_number = trans_df.shape[0]

    filt_trans_MT = trans_df.value_counts(['SJs']).rename_axis(['SJs']).reset_index(name='counts')
    for i in range(filt_trans_MT.shape[0]):
        sj = filt_trans_MT.iloc[i].at['SJs']
        jaq = trans_df.loc[trans_df['SJs']==sj]['mean_JAQ'].mean()
        filt_trans_MT.loc[i, 'mean_JAQ'] = round(jaq, 2)
    filt_trans_MT = filt_trans_MT.loc[filt_trans_MT['counts'] >= 0.005*read_number].reset_index(drop=True)
    filt_trans_MT['percentage'] = round(filt_trans_MT['counts'] / read_number * 100, 2)
    filt_trans_MT['transcript'] = ''
    filt_trans_MT['exons'] = ''
    filt_trans_MT['length'] = 0
    for i in range(filt_trans_MT.shape[0]):
        exs_dict = {'start': [], 'end': []}
        junc = filt_trans_MT.iloc[i].at['SJs'].split('; ')
        for el in junc:
            start = el.split(', ')[1].split(')')[0]
            end = el.split(', ')[0].split('(')[1]
            if filt_trans_MT.loc[i, 'exons'] == '':
                filt_trans_MT.loc[i, 'exons'] += '({}, '.format(exons_dict['start'][0]) + end + '); (' + start + ', '
                exs_dict['start'].append(exons_dict['start'][0])
                exs_dict['end'].append(int(end))
                exs_dict['start'].append(int(start))
                if el == junc[-1]:
                    filt_trans_MT.loc[i, 'exons'] += str(exons.iloc[-1].at['end']) + ')'
                    exs_dict['end'].append(exons.iloc[-1].at['end'])
            else:
                filt_trans_MT.loc[i, 'exons'] += end + '); (' + start + ', '
                exs_dict['start'].append(int(start))
                exs_dict['end'].append(int(end))
                if el == junc[-1]:
                    filt_trans_MT.loc[i, 'exons'] += str(exons.iloc[-1].at['end']) + ')'
                    exs_dict['end'].append(exons.iloc[-1].at['end'])
        for i2 in range(len(exs_dict['start'])):
            st = exs_dict['start'][i2]
            en = exs_dict['end'][i2]
            leng = en -st
            filt_trans_MT.loc[i, 'length'] += leng
            if exs_dict['start'][i2] in exons_dict['start']:
                index = exons_dict['start'].index(exs_dict['start'][i2])
                if exs_dict['end'][i2] == exons_dict['end'][index]:
                    if filt_trans_MT.loc[i, 'transcript'] == '':
                        filt_trans_MT.loc[i, 'transcript'] += exons_dict['name'][index]
                    else:
                        filt_trans_MT.loc[i, 'transcript'] += '-' + exons_dict['name'][index]
                else:
                    if filt_trans_MT.loc[i, 'transcript'] == '':
                        filt_trans_MT.loc[i, 'transcript'] += 'unclear'
                    else:
                        filt_trans_MT.loc[i, 'transcript'] += '-unclear'
            else:
                if filt_trans_MT.loc[i, 'transcript'] == '':
                    filt_trans_MT.loc[i, 'transcript'] += 'unclear'
                else:
                    filt_trans_MT.loc[i, 'transcript'] += '-unclear'


    filt_trans_MT = filt_trans_MT[['transcript', 'length', 'counts', 'percentage', 'mean_JAQ', 'exons', 'SJs']]

    gff3_dict = {'seqid': [], 'source': [], 'type': [], 'start': [], 'end': [], 'score': [], 'strand': [], 'phase': [], 'attributes': []}
    seqid = gff3.split('/')[-1].split('.')[0]
    source = '.'

    for i in range(filt_trans_MT.shape[0]):
        #Specify transcript-level data
        gff3_dict['seqid'].append(seqid)
        gff3_dict['source'].append(source)
        gff3_dict['type'].append('mRNA')
        gff3_dict['start'].append(exons.iloc[0].at['start'])
        gff3_dict['end'].append(exons.iloc[-1].at['end'])
        gff3_dict['score'].append(filt_trans_MT.iloc[i].at['mean_JAQ'])
        gff3_dict['strand'].append('+')
        gff3_dict['phase'].append('.')
        if filt_trans_MT.iloc[i].at['transcript'].split('-') != exons.annotations.to_list():
            gff3_dict['attributes'].append('ID=transcript'+str(i+1)+';'+'Name='+filt_trans_MT.iloc[i].at['transcript']+';Note='+str(filt_trans_MT.iloc[i].at['percentage'])+';color=#ff0000')
        else:
            gff3_dict['attributes'].append('ID=transcript'+str(i+1)+';'+'Name='+filt_trans_MT.iloc[i].at['transcript']+';Note='+str(filt_trans_MT.iloc[i].at['percentage'])+';color=#00b300')

        #Add exons-level data
        exs = filt_trans_MT.iloc[i].at['exons'].split('; ')
        for i2 in range(len(exs)):
            start = exs[i2].split(', ')[0].split('(')[1]
            end = exs[i2].split(', ')[1].split(')')[0]
            gff3_dict['seqid'].append(seqid)
            gff3_dict['source'].append(source)
            gff3_dict['type'].append('exon')
            gff3_dict['start'].append(start)
            gff3_dict['end'].append(end)
            gff3_dict['score'].append(filt_trans_MT.iloc[i].at['mean_JAQ'])
            gff3_dict['strand'].append('+')
            gff3_dict['phase'].append('.')
            try:
                gff3_dict['attributes'].append('ID='+filt_trans_MT.iloc[i].at['transcript'].split('-')[i2]+';Name='+filt_trans_MT.iloc[i].at['transcript'].split('-')[i2]+';Parent=transcript'+str(i+1))
            except:
                print(filt_trans_MT.iloc[i].at['transcript'])
                gff3_dict['attributes'].append('ID='+str(i2))


    gff3_MT = pd.DataFrame.from_dict(gff3_dict)
    file_line = """##gff-version 3.1.26\
    {}"""
    with open('/{}/{}_MT_transcripts.gff3'.format(sample, sample), 'w') as fp:
        fp.write(file_line.format(gff3_MT.to_csv(sep='\t', index=False)))


    with pd.ExcelWriter('/{}/NanoSplicer_{}_quantifications.xlsx'.format(sample, sample)) as writer:
        filt_trans_WT.to_excel(writer, sheet_name='WT_transcripts', index=False)
        filt_trans_MT.to_excel(writer, sheet_name='MT_transcripts', index=False)

    c += 1
    print('Finished processing results for {}'.format(sample))
