#!/usr/bin/env python3

from Bio import SeqIO
from Bio import Align
from Bio.SeqRecord import SeqRecord
from Bio.SeqRecord import Seq
import pandas as pd
import argparse


def assign(samplefas, CSStable):

    CSSpos = []
    CSS_df = pd.read_csv(CSStable, dtype=object)
    
    # Get zero-based polymorphic sites from CSStable
    for pos in CSS_df['POS']:
        site = int(pos)
        pysite = site-1
        CSSpos.append(pysite)
    
    # Extract subsequence of defining positions for sample
    samplefullseq = SeqIO.read(samplefas, 'fasta')
    sampleSNPs = []
    for val in CSSpos:
        sampleSNPs.append(samplefullseq[val])
    sampleCSS = SeqRecord(Seq(''.join(sampleSNPs)), samplefullseq.id, '', '')
    nCount = sampleCSS.count('N')

    # Define aligner and scoring parameters
    aligner = Align.PairwiseAligner(mismatch_score = -8, open_gap_score = -30,
                                    extend_gap_score = -50,  wildcard = 'N')

    # Drop non-clade columns
    CSS_df.drop('POS', axis=1, inplace=True)
    CSS_df.drop('Monomorphic', axis=1, inplace=True)
    CSS_df.drop('Numclades', axis=1, inplace=True)
    
    # Calculate number of clade specific SNPS (i.e rows in table)
    CSStested = CSS_df.shape[0]

    # Create Seq objects from each clade CSS and score alignment with sample 
    clade_dfs = [CSS_df[[col]] for col in CSS_df.columns]
    clade_scores = []
    for x in clade_dfs:
        clade = x.columns[-1]
        bases = ''.join(x.iloc[:,0])
        fasta = SeqRecord(Seq(bases), clade, '', '')
        score = aligner.score(sampleCSS, fasta)
        alignment = aligner.align(sampleCSS, fasta)
        count = alignment[0].counts()
        #print(alignment[0])
        matches = round((count.identities/CSStested*100), 3)
        mismatches = round((((count.mismatches)-(nCount))/CSStested*100), 3)
        nocov = round((nCount/CSStested*100), 3)
        gaps = round(((count.gaps)/CSStested*100), 3)

        clade_scores.append({'Sample':samplefullseq.id, 'flag':'', 'group':fasta.id,
                            'CSSTested':CSStested, 'Score':score, 'matches':matches,
                            'mismatches':mismatches, 'nocoverage':nocov, 'anomolous':gaps})
    
    data = pd.DataFrame(clade_scores)
    sortedScore_df = data.sort_values('Score', ascending=False).head(10)
    
    sortedScore_df.loc[(sortedScore_df['Score'] <= 2500), 'flag'] = 'Outlier'
    sortedScore_df.loc[(sortedScore_df['Score'] > 2800), 'flag'] = 'Match'
    sortedScore_df.loc[(sortedScore_df['Score'] > 2500) & (sortedScore_df['Score'] <= 2800),
                       'flag'] = 'Tentative'
    #sortedScore_df.loc[(sortedScore_df['flag'] == 'Outlier'), 'group'] = 'NoMatch'
    print(sortedScore_df)

    
    # Create output compatable with btb-seq
    # Sample,flag,group,CSSTested,matches,mismatches,noCoverage,anomalous
    """
    'Sample' is samplefullseq.id
    'flag' determined from how close score is to numder of sites TODO - set threshold e.g. 95%?
    'group' top clade from sortedScore_df
    'CSSTested' = CSS_df.shape[0]
    'matches' = identity counts: https://biopython.org/docs/latest/api/Bio.Align.html#Bio.Align.Alignment.counts
    'mismatches' = mismatch counts - as above 
    'noCoverage' = sampleCSS.count('N')
    'anomalous' = gap counts - as above
    """
    
    sortedScore_df.head(1).to_csv("{}_cladematch.csv".format(samplefullseq.id), index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('samplefas')
    parser.add_argument('CSStable')

    args = parser.parse_args()

    assign(**vars(args))
