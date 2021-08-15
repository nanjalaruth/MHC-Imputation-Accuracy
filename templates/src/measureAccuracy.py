#!/users/nanje/miniconda3/bin/python

#-*- coding: utf-8 -*-
import sys, os

HLA_names = ["A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1"]

def measureAccuracy(answerfile, predictfile, genes, outfile=None, __asSTDOUT = False, __only4digits=False):


    if isinstance(genes, list):
        if len(genes) == 0 or (len(genes) == 1 and genes[0] == "all"):
            genes = "A B C DRB1 DPA1 DPB1 DQA1 DQB1".split()
    elif isinstance(genes, str):
        if genes == 'all':
            genes = "A B C DRB1 DPA1 DPB1 DQA1 DQB1".split()


    # Accuracy dictionary
    __RETURN__ = {'2D' : {_hla: None for _hla in HLA_names},
                  '4D' : {_hla: None for _hla in HLA_names}}

    __asFileWrite = outfile and not __asSTDOUT

    if __asFileWrite:
        fo = open(outfile, 'w')


    for gene in genes:
        correct2d=0
        total2d=0
        correct4d=0
        total4d=0

        answers2d={}
        answers4d={}
        with open(answerfile) as fin:
            for l in fin:
                c=l.split()
                ID=c[0]+' '+c[1]
                if c[2] == gene and len(c) > 3:
                    answers2d[ID]=c[3].split(',')
                    answers4d[ID]=c[4].split(',')


        with open(predictfile) as fin:
            for l in fin:
                c=l.split()
                ID=c[0]+' '+c[1]
                if c[2] == gene:
                    predict2d=c[3].split(',')
                    predict4d=c[4].split(',')
                    if ID in answers2d:
                        # (correct, total)=compare_and_score(predict2d, answers2d[ID])
                        (correct, total)=compare_and_score2(predict2d, answers2d[ID])
                        correct2d+=correct
                        total2d+=total
                    if ID in answers4d:
                        # (correct, total)=compare_and_score(predict4d, answers4d[ID])
                        (correct, total)=compare_and_score2(predict4d, answers4d[ID])
                        correct4d+=correct
                        total4d+=total


        if not __only4digits:

            try:
                __RETURN__['2D'][gene] = float(correct2d)/total2d
            except ZeroDivisionError:
                print("[measureAccuracy.py::Error] No 2-digit alleles for HLA_{}.".format(gene))
                continue

        try:
            __RETURN__['4D'][gene] = float(correct4d)/total4d
        except ZeroDivisionError:
            print("[measureAccuracy.py::Error] No 4-digit alleles for HLA_{}.".format(gene))
            continue


        if __asSTDOUT:

            if not __only4digits:
                sys.stdout.write("%s\t2D\t%.5f\n"%(gene, float(correct2d)/total2d))
            sys.stdout.write("%s\t4D\t%.5f\n"%(gene, float(correct4d)/total4d))

        if __asFileWrite:

            if not __only4digits:
                fo.write("%s\t2D\t%.5f\n"%(gene, float(correct2d)/total2d))
            fo.write("%s\t4D\t%.5f\n"%(gene, float(correct4d)/total4d))

    if __asFileWrite:
        fo.close()


    return outfile




def compare_and_score(predict, answer):
    ## return value is correct, total
    if answer[0] == '':
        answer[0]='answer1'
    if answer[1] == '':
        answer[1]='answer2'
    if predict[0] == '':
        predict[0]='predict1'
    if predict[1] == '':
        predict[1]='predict2'

    correct=max((answer[0]==predict[0])+(answer[1]==predict[1]), \
            (answer[0]==predict[1])+(answer[1]==predict[0]))
    return(correct,2)


def compare_and_score2(predict, answer):

    # Made by Wanson.

    if answer[0] == '' or answer[1] == '':
        return (0, 0)

    if predict[0] == '':
        predict[0]='predict1'
    if predict[1] == '':
        predict[1]='predict2'

    correct=max((answer[0]==predict[0])+(answer[1]==predict[1]), \
            (answer[0]==predict[1])+(answer[1]==predict[0]))

    return(correct, 2)




if __name__ == "__main__":

    """
    < measureAccuracy.py >
    
    Module to get accuracy(%) of given '*.alleles' file.
    
    
    """

    answerfile = sys.argv[1]
    predictfile = sys.argv[2]
    genes = sys.argv[3:]

    measureAccuracy(answerfile, predictfile, genes, __asSTDOUT=True)

    ### < Manual Testing >

    # # 1958BC / T1DGC_REF
    #
    # #_2_Plain
    # measureAccuracy(answerfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/answer/1958BC_IC.maf0.001.alleles',
    #                 predictfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_2_Plain/1958BC_T1DGC_REF.Plain.MHC.HLA_IMPUTATION_OUT.alleles',
    #                 genes='all', __only4digits=True, outfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_2_Plain/1958BC_T1DGC_REF.Plain.MHC.HLA_IMPUTATION_OUT.alleles.fixed.accuracy')
    #
    # # _3_MM
    # measureAccuracy(answerfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/answer/1958BC_IC.maf0.001.alleles',
    #                 predictfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_3_MM/1958BC_T1DGC_REF.MM.MHC.HLA_IMPUTATION_OUT.alleles',
    #                 genes='all', __only4digits=True, outfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_3_MM/1958BC_T1DGC_REF.MM.MHC.HLA_IMPUTATION_OUT.alleles.fixed.accuracy')
    #
    # # _4_AGM_HapMap_Map
    # measureAccuracy(answerfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/answer/1958BC_IC.maf0.001.alleles',
    #                 predictfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_4_AGM_HapMap_Map/1958BC_T1DGC_REF.AGM_HapMap_Map.MHC.HLA_IMPUTATION_OUT.alleles',
    #                 genes='all', __only4digits=True, outfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_4_AGM_HapMap_Map/1958BC_T1DGC_REF.AGM_HapMap_Map.MHC.HLA_IMPUTATION_OUT.alleles.fixed.accuracy')
    #
    # # _5_AGM
    # measureAccuracy(answerfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/answer/1958BC_IC.maf0.001.alleles',
    #                 predictfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_5_AGM/1958BC_T1DGC_REF.AGM.MHC.HLA_IMPUTATION_OUT.alleles',
    #                 genes='all', __only4digits=True, outfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_5_AGM/1958BC_T1DGC_REF.AGM.MHC.HLA_IMPUTATION_OUT.alleles.fixed.accuracy')
    #
    # # _6_MM_AGM
    # measureAccuracy(answerfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/answer/1958BC_IC.maf0.001.alleles',
    #                 predictfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_6_MM_AGM/1958BC_T1DGC_REF.MM.AGM.MHC.HLA_IMPUTATION_OUT.alleles',
    #                 genes='all', __only4digits=True, outfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/1958BC_T1DGC_REF/trial1/_6_MM_AGM/1958BC_T1DGC_REF.MM.AGM.MHC.HLA_IMPUTATION_OUT.alleles.fixed.accuracy')
    #
    #
    # # HM_CEU / T1DGC_REF
    # measureAccuracy(answerfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/data/answer/HM_CEU_REF.bgl.phased.FIDadj.alleles.answer',
    #                 predictfile='/Users/wansun/Dropbox/_Sync_MyLaptop/Projects/CookHLA/Figure/AccuracyTables/HM_CEU_T1DGC_REF/trial3/HM_CEU_T1DGC_REF.MM.AGM.MHC.HLA_IMPUTATION_OUT.alleles',
    #                 genes='all', __asSTDOUT=True, __only4digits=True)

