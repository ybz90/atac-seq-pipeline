import os
import subprocess


# Input location of fastq reads for alignment
# e.g. fastq_loc = '/mnt/promoter4/Bin_Raw_Sequence_Data/2016_07_08_DIP_B/PE-07082016-DIPB'
fastq_loc = 'FASTQ_INPUT_FOLDER'

# Output location for nodup bam alignments
nodup_bam_loc = 'NODUP_BAM_LOC'
if not os.path.exists(os.path.expanduser(nodup_bam_loc)):
    os.makedirs(os.path.expanduser(nodup_bam_loc))

# Define prefixes to use and corresponding fastq index
fastq_prefix = {
    'e11.5FB_1': 'HH167',
}


'''
################################
# ALIGN FASTQ READS WITH BOWTIE2
################################

Set the requisite parameters indicating the location of the fastq.gz
reads and desired output location of bowtie2 alignments. Also, define
new prefixes for each experiment and specify the corresponding fastq
index prefix.
'''

# Shared silencer2 bowtie2 indices
mm10_bw2_index = '/mnt/silencer2/share/bowtie2_indexes/mm10'

for curr_bam in fastq_prefix.keys():

    curr_fastq = fastq_prefix[curr_bam]

    cmd_bowtie = pipe(['bowtie2 -p %d ' % 6 +
                       '-t -X2000 --no-mixed --no-discordant -x %s ' % mm10_bw2_index +
                       '-1 %s/%s_R1.fastq.bz2 ' % (fastq_loc, curr_fastq) +
                       '-2 %s/%s_R2.fastq.bz2 ' % (fastq_loc, curr_fastq),
                       'samtools view -bS - > %s/%s.bam ' % (nodup_bam_loc, curr_bam)
                       ])
    print cmd_bowtie
    subprocess.call(cmd_bowtie, shell=True)


'''
Remove duplicate reads from the paired alignments and generate index.
'''

for curr_bam in fastq_prefix.keys():
    cmd_rmdup = pipe([
        'samtools sort %s/%s.bam' % (nodup_bam_loc, curr_bam),
        'samtools rmdup - %s/%s_sorted_nodup.bam' % (
            nodup_bam_loc, curr_bam)
    ])
    cmd_index = 'samtools_index %s/%s_sorted_nodup.bam %s/%s_sorted_nodup.bai' % (
        nodup_bam_loc, curr_bam, nodup_bam_loc, curr_bam)
    print '; '.join([cmd_rmdup, cmd_index])
    subprocess.call('; '.join([cmd_rmdup, cmd_index]), shell=True)
