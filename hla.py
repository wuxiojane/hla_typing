# -*- coding: utf-8 -*-
# !/usr/bin/python3

# @Date    : 2021-12-07
# @Author  : TTG
# @description：


import sys
import os
import re
import subprocess
import argparse
from multiprocessing import Process
import shutil
import pandas as pd

class Params:
    """
    Class for top-level system parameters
    """

    def __init__(self):

        usage = """
            python3 hla.py  { -1 FASTQ1 -2 FASTQ2 | -S BAM } [-N SAMPLENAME ] [-O OUTPUT_PATH]  [options]*
            """
        self.parser = argparse.ArgumentParser(usage=usage, prog='hla.py')
        self.parser.add_argument('-1', '--fastq1',
                                 dest='fastq1',
                                 default=None,
                                 help='Fastq file,fastq1')
        self.parser.add_argument('-2', '--fastq2',
                                 dest='fastq2',
                                 default=None,
                                 help='Fastq file,fastq2')
        self.parser.add_argument('-S', '--bam',
                                 dest='bam',
                                 default=None,
                                 help='Alignment files after genome alignment, both mapped and unmapped')
        self.parser.add_argument('-N', '--sampleName',
                                 dest='sampleName',
                                 default='sample',
                                 help='Sample name')
        self.parser.add_argument('-G', '--genome',
                                 dest='genome',
                                 default='hg19',
                                 help='Genome version [hg19/hg38], default:hg19')
        self.parser.add_argument('-T', '--threads',
                                 dest='threads',
                                 default='8',
                                 help='threads number')
        self.parser.add_argument('-D', '--hlahdDir',
                                 dest='hlahdDir',
                                 default='/Juno/software/hlahd.1.4.0/',
                                 help='hlahd DirPath')
        self.parser.add_argument('-O', '--outputPath',
                                 dest='outputPath',
                                 default='./',
                                 help='output path')                         
                                 

        args = self.parser.parse_args()
        if args.fastq1:
            self.FASTQ1 = args.fastq1
        else:
            self.FASTQ1 = 'False'
        if args.fastq2:
            self.FASTQ2 = args.fastq2
        else:
            self.FASTQ2 = 'False'
        if args.bam:
            self.BAM = args.bam
        else:
            self.BAM = 'False'
        if args.sampleName:
            self.SAMPLENAME = args.sampleName
        else:
            self.SAMPLENAME = 'sample'
        if args.genome:
            self.GENOME = args.genome
        else:
            self.GENOME = 'hg19'
        if args.threads:
            self.THREADS = args.threads
        else:
            self.THREADS = '8'
        if args.hlahdDir:
            self.HLAHDDIR = args.hlahdDir
        else:
            self.HLAHDDIR = '/software/hlahd.1.4.0/'  
        if args.outputPath:
            self.OUTPUTPATH = args.outputPath

        if not args.fastq1 and not args.bam:
            subprocess.check_call('python3 %s -h' % (sys.argv[0]), shell=True)
            print('\nThere is something wrong in the comand!')
            sys.exit(1)

    def check(self):
        logfile = self.OUTPUTPATH + '/' + self.SAMPLENAME + '/' + self.SAMPLENAME + '.HLA_typing.log'
        f = open(logfile, 'w')
        f.write('fastq1 :   ' + self.FASTQ1 + '\n')
        f.write('fastq2 :   ' + self.FASTQ2 + '\n')
        f.write('bam    :   ' + self.BAM + '\n')
        f.write('threads    :   ' + self.THREADS + '\n')
        f.write('reference  :   ' + self.GENOME + '\n')
        f.write('sampleName  :  ' + self.SAMPLENAME + '\n')
        f.write('outputPath  :  ' + self.OUTPUTPATH + self.SAMPLENAME + '\n')
        f.write('HLA typing results  :  ' + self.OUTPUTPATH + '/' + self.SAMPLENAME + '/' + self.SAMPLENAME + '.HLA.typing.result.txt' + '\n')
        f.close()


def fromfq(threads, index, fq1, fq2, outsam, outbam):
    args = 'bowtie2 -p %s -x %s -1 %s -2 %s -S %s' % (threads, index, fq1, fq2, outsam)
    subprocess.check_call(args, shell=True)
    args2 = 'samtools view -h -b -F 4 %s > %s' % (outsam, outbam)
    subprocess.check_call(args2, shell=True)


def frombam(inputbam, mappedbam, unmappedbam, mergebam):
    args = 'samtools view -h -b %s chr6:28477797-33480577 > %s' % (inputbam, mappedbam)
    subprocess.check_call(args, shell=True)
    args2 = 'samtools view -h -b -f 4 %s > %s' % (inputbam, unmappedbam)
    subprocess.check_call(args2, shell=True)
    args3 = 'samtools merge %s %s %s' % (mergebam, mappedbam, unmappedbam)
    subprocess.check_call(args3, shell=True)


def sortbam(inputbam, outputbam):
    args = 'samtools sort -n %s -o %s' % (inputbam, outputbam)
    subprocess.check_call(args, shell=True)


def bamtofq(inputbam, outfq1, outfq2):
    args = 'samtools fastq %s -1 %s -2 %s -0 /dev/null -s /dev/null -n' % (inputbam, outfq1, outfq2)
    subprocess.check_call(args, shell=True)


def hlahd(threads, hlahddir, fq1, fq2, sample, outpath):
    freq_data = hlahddir + '/freq_data'
    split_d = hlahddir + '/HLA_gene.split.3.32.0.txt'
    dicd = hlahddir + '/dictionary'
    args = 'hlahd.sh -t %s -f %s %s %s %s %s %s %s' % (threads, freq_data,fq1, fq2, split_d, dicd, sample, outpath)
    subprocess.check_call(args, shell=True)


def removeTempFolder(tempdir):
    shutil.rmtree(tempdir)


def extract(resDir, hlahddir, sample):
    genelis = ['A', 'B', 'C', 'DRB1', 'DQB1']
    maplist = resDir + '/' + sample + '/maplist/maplist.txt'
    for i in range(len(genelis)):
        gene = genelis[i]
        maplist_L = resDir + '/' + sample + '/maplist/maplist_' + gene + '.txt'
        est_fi = resDir + '/' + sample + '/result/' + sample + '_' + gene + '.est.txt'
        read_fi = resDir + '/' + sample + '/result/' + sample + '_' + gene + '.read.txt'
        log_fi = resDir + '/' + sample + '/log/' + sample + '_' + gene + '.log'
        freq_fi = hlahddir +'/freq_data/' + gene + '_count.txt'
        tmp_fi = resDir + '/' + sample + '/result/' +sample + '.HLA.typing.tmp.txt'
        final_fi = resDir + '/' + sample + '.HLA.typing.result.txt'
        args = 'hla_est -L 100 --ml 50 -m 100 --hth 4.0 %s %s %s -o %s --pread %s > %s' % (gene, maplist_L, maplist, est_fi, read_fi, log_fi)
        subprocess.check_call(args, shell=True)
        args2 = 'pick_up_allele %s %s -f %s >> %s' % (gene, est_fi, freq_fi, tmp_fi)
        subprocess.check_call(args2, shell=True)
    f = open(tmp_fi, 'r')
    dat = f.readlines()
    f.close()
    bb = {}
    for line in dat:
        line = line.strip()
        line = re.sub(r'HLA.{2,5}\*','',line)
        kk = line.split('\t')[0] + "*"
        vv1 = line.split('\t')[1].rsplit(':',1)[0]
        vv2 = line.split('\t')[2].rsplit(':',1)[0]
        if vv2 =='-':
            vv2 = vv1
    vv = [vv1,vv2]
    bb[kk]=vv
    pd.DataFrame(bb).to_csv(final_fi,index=False,sep= '\t')

params = Params()

sampleDir = params.OUTPUTPATH + '/' + params.SAMPLENAME
tmpDir = params.OUTPUTPATH + '/' + params.SAMPLENAME + '/tmp'
bamDir = params.OUTPUTPATH + '/' + params.SAMPLENAME + '/bam'
fqDir = params.OUTPUTPATH + '/' + params.SAMPLENAME + '/fastq'
hlaDir = params.OUTPUTPATH + '/' + params.SAMPLENAME + '/HLA'
os.makedirs(sampleDir, exist_ok=True)
os.makedirs(tmpDir, exist_ok=True)
os.makedirs(bamDir, exist_ok=True)
os.makedirs(fqDir, exist_ok=True)
os.makedirs(hlaDir, exist_ok=True)

params.check()

INDEX = params.HLAHDDIR + '/data/hla_gene'
outfq1 = fqDir + '/' + params.SAMPLENAME + '.hla.R1.fq'
outfq2 = fqDir + '/' + params.SAMPLENAME + '.hla.R2.fq'

if params.FASTQ1 != 'False':
    FASTQ1 = params.FASTQ1
    FASTQ2 = params.FASTQ2
    outsam = tmpDir + '/' + params.SAMPLENAME + '.hlamap.sam'
    outbam = tmpDir + '/' + params.SAMPLENAME + '.hla.mapped.bam'
    sortedbam = bamDir + '/' + params.SAMPLENAME + '.hla.mapped.sorted.bam'

    rfromfq = Process(target=fromfq, args=(params.THREADS, INDEX, FASTQ1, FASTQ2, outsam, outbam))
    if os.path.exists(outbam):
        pass
    else:
        try:
            rfromfq.start()
            rfromfq.join()
        except:
            print('from fastq to bam ,failed')
            sys.exit(1)
    rsortbam = Process(target=sortbam, args=(outbam, sortedbam))
    if os.path.exists(sortedbam):
        pass
    else:
        try:
            rsortbam.start()
            rsortbam.join()
        except:
            print('sort bam ,failed')
            sys.exit(1)

    rbamtofq = Process(target=bamtofq, args=(sortedbam, outfq1, outfq2))
    if os.path.exists(outfq1) and os.path.exists(outfq1):
        pass
    else:
        try:
            rbamtofq.start()
            rbamtofq.join()
        except:
            print('bamToFastq ,failed')
            sys.exit(1)

if params.BAM != 'False':
    mappedbam = bamDir + '/' + params.SAMPLENAME + '.mapped.bam'
    unmappedbam = bamDir + '/' + params.SAMPLENAME + '.unmapped.bam'
    mergebam = tmpDir + '/' + params.SAMPLENAME + '.merge.bam'
    sortedbam_ = bamDir + '/' + params.SAMPLENAME + '.merge.sorted.bam'


    rfrombam = Process(target=frombam, args=(params.BAM, mappedbam, unmappedbam, mergebam))
    if os.path.exists(unmappedbam) and os.path.exists(mappedbam):
        pass
    else:
        try:
            rfrombam.start()
            rfrombam.join()
        except:
            print('from bam ,failed')
            sys.exit(1)

    rsortbam = Process(target=sortbam, args=(mergebam, sortedbam_))
    if os.path.exists(sortedbam_):
        pass
    else:
        try:
            rsortbam.start()
            rsortbam.join()
        except:
            print('sort bam ,failed')
            sys.exit(1)

    rbamtofq = Process(target=bamtofq, args=(sortedbam_, outfq1, outfq2))
    if os.path.exists(outfq1) and os.path.exists(outfq1):
        pass
    else:
        try:
            rbamtofq.start()
            rbamtofq.join()
        except:
            print('bamToFastq ,failed')
            sys.exit(1)

rhlahd = Process(target=hlahd, args=(params.THREADS, params.HLAHDDIR, outfq1, outfq2, params.SAMPLENAME, hlaDir))
try:
    rhlahd.start()
    rhlahd.join()
except:
    print('hlahd typing,failed')
    sys.exit(1)
rextract = Process(target=extract, args=(hlaDir, params.HLAHDDIR, params.SAMPLENAME))
try:
    rextract.start()
    rextract.join()
except:
    print('hlahd extract results,failed')
    sys.exit(1)

removeTempFolder(tmpDir)
