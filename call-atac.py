import numpy as np
import os
import subprocess
from multiprocessing import Pool, Queue, Process, Manager, cpu_count

# Set global parameters

# Chromosome sizes for clipping and bigwig generation
mm_chr_sz = '~/annotation/chromosomes/mm10/mm10.chrom.sizes'

# Input location for nodup bam alignments
nodup_bam_loc = 'NODUP_BAM_LOC'

# Specify experiment condition
# location further processing and analysis for experiment
expt_prefix = 'EXPT_PREFIX'
pipe_out_loc = 'EXPT_OUT_LOC'



'''
#################################
# LOAD MODULES AND INIT FUNCTIONS
#################################
'''

#
def pipe(cmd_list):
    return ' | '.join(cmd_list)

# Function to run shell command -- crude approach to parallelization;
# argument is entire shell command, so multi_arg is list of commands to run
def run_bash(command):
    print command
    #subprocess.call(command, shell=True)
    # Set analysis output directory in run_bash before every command
    subprocess.call('cd %s; %s' %(pipe_out_loc, command), shell=True)
    return None

#
def get_numpeaks(file_loc):
    with open(file_loc) as in_file:
        for num, line in enumerate(in_file):
            pass
    return num

# Remove control_lambda and treat_pileups, if no other bdgcmp operations
# will be performed
def remove_bdg(outdirs):
    cmd_rm = ['rm -f %s/*.bdg ' %outdir for outdir in outdirs]
    return '; '.join(cmd_rm)

# Remove mitochondrial reads and adjust reads for Tn5 insertion
def cmd_procalign(in_bam_loc, expt_prefix, in_bam_rep):

    prefix = expt_prefix + '_' + str(in_bam_rep)
    prefix2 = expt_prefix + '.' + str(in_bam_rep)

    '''Faster method might be to awk the header to get all non-chrM headers,
    and then just use samtools view chr*.'''

    cmd_remove_mito = pipe([
        'samtools view -h %s/%s_sorted_nodup.bam' % (in_bam_loc, prefix),
        'awk -F "\\t" \'$3 != "chrM" { print $0 }\'',
        'samtools view -bS - > ' +
        '%s.sorted.nodup.non-chrM.bam' % prefix2,
    ])
    # print cmd_remove_mito

    cmd_bam2bed = pipe([
        'bamToBed -i ' +
        '%s.sorted.nodup.non-chrM.bam' % prefix2,
        'awk -F "\\t" \'BEGIN {OFS=FS} ' +
        '$6 == "+" {$2 = $2 + 4} $6 == "-" {$3 = $3 - 5}; 1\'',
        'gzip -c > ' +
        '%s.sorted.nodup.non-chrM.tn5.bed.gz ' % prefix2
    ])
    # print cmd_bam2bed

    return '; '.join([cmd_remove_mito, cmd_bam2bed])

# Generate flagstat summary
def cmd_flagstat(expt_prefix):
    flagstat = ['samtools flagstat ' + '%s.%s.sorted.nodup.non-chrM.bam' % (expt_prefix, rep_no) +
                ' > flagstat%s.txt ' % rep_no for rep_no in ['1', '2']]
    return '; '.join(flagstat)

# Get counts
def cmd_readcounts(flagstat_loc):
    # Load flagstat counts; if pooled, add them together for scale_factor
    counts = {}
    for curr_no in [1,2]:
        with open('%s/flagstat%d.txt' % (flagstat_loc, curr_no)) as flagstat:
            num_mapped = flagstat.readlines()[4].split(' ')[0]
        num_scaled = int(num_mapped) / 1000000.
        #print num_mapped, num_scaled
        counts[curr_no] = num_scaled
    # for pooled obvi, will be R1 + R2
    return counts

# Call narrowPeaks, and generate treat_pileup and control_lambda bdg's
def cmd_callpeak(prefix, t_files, outdir, p_val, shift_val, ext_val):

    command = 'macs2 callpeak ' + \
        '-t %s ' %t_files + \
        '-f BED -n %s ' % prefix + \
        '--outdir %s ' % outdir + \
        '-g mm -p %s --nomodel --shift -%s --extsize %s ' % (p_val, shift_val, ext_val) + \
        '--keep-dup all -B --SPMR --call-summits '

    return command

# Fold-enrichment
def cmd_FE(mm_chr_sz, outdir, prefix):

    cmd_bdgcmp = 'macs2 bdgcmp ' + \
            '-t %s/%s_treat_pileup.bdg ' %(outdir, prefix) + \
            '-c %s/%s_control_lambda.bdg ' %(outdir, prefix) + \
            '--outdir %s -o %s.FE.bdg ' % (outdir, prefix) + \
            '-m FE'
    #print cmd_bdgcmp
    cmd_trimsort = pipe([
        'sort -k1,1 -k2,2n %s/%s.FE.bdg' % (outdir, prefix),
        'slopBed -i stdin -g %s -b 0' % mm_chr_sz,
        'bedClip stdin %s %s/%s.FE.trim.bdg ' % (mm_chr_sz, outdir, prefix)
    ])
    #print cmd_trimsort
    cmd_bw = 'bedGraphToBigWig %s/%s.FE.trim.bdg %s %s/%s.FE.bw ' % (outdir, prefix, mm_chr_sz, outdir, prefix)
    #print cmd_bw

    return '; '.join([cmd_bdgcmp, cmd_trimsort, cmd_bw])

# P-value
def cmd_ppois(mm_chr_sz, outdir, prefix, scale_factor):

    cmd_bdgcmp = 'macs2 bdgcmp ' + \
            '-t %s/%s_treat_pileup.bdg ' %(outdir, prefix) + \
            '-c %s/%s_control_lambda.bdg ' %(outdir, prefix) + \
            '--outdir %s -o %s.pval.bdg ' % (outdir, prefix) + \
            '-m ppois -S %s' % str(scale_factor)
    #print cmd_bdgcmp
    cmd_trimsort = pipe([
        'sort -k1,1 -k2,2n %s/%s.pval.bdg' % (outdir, prefix),
        'slopBed -i stdin -g %s -b 0' % mm_chr_sz,
        'bedClip stdin %s %s/%s.pval.trim.bdg ' % (mm_chr_sz, outdir, prefix)
    ])
    #print cmd_trimsort
    cmd_bw = 'bedGraphToBigWig %s/%s.pval.trim.bdg %s %s/%s.pval.bw ' % (outdir, prefix, mm_chr_sz, outdir, prefix)
    #print cmd_bw

    return '; '.join([cmd_bdgcmp, cmd_trimsort, cmd_bw])


'''
################################
# READ-ALIGNMENT POST-PROCESSING
################################

First, remove mitochondrial reads by checking header or for ChrM.

Convert bam to bed s.t. paired reads are treated as individual
single reads, so both ends (cut sites) will pileup instead of just
the 5' end of a read pair. (-bampe will retain mate information.)

Adjust the ends to account for actual cut site due to Tn5 insertion
position: add 4 to '+' strands and subtract 5 from '-' strands.
'''

# Post-process read alignments for both replicates
procalign_cmds = [cmd_procalign(
    nodup_bam_loc, expt_prefix, rep_num) for rep_num in [1, 2]]
# print '\n'.join(procalign_cmds) + '\n'

if __name__ == '__main__':
    pool = Pool(processes=2)
    pool.map(run_bash, procalign_cmds)
    pool.close()

# Generate flagstat summary statistics for both reps, load into counts dict
# print cmd_flagstat(expt_prefix) + '\n'
run_bash(cmd_flagstat(expt_prefix))

counts = cmd_readcounts(os.path.expanduser(pipe_out_loc))
print counts, '\n'


'''
##########################################
# PEAK CALLING AND SIGNAL TRACK GENERATION
##########################################

Call peaks on individual and pooled replicates. Generate FE and
-log10(p-val) signal tracks.
'''

# Params for peak calling
shift_val = '37'
ext_val = '73'
p_val = '1e-2'

# File locations of post-processed alignments for both replicates
R1_bedgz = '%s.%d.sorted.nodup.non-chrM.tn5.bed.gz' % (expt_prefix, 1)
R2_bedgz = '%s.%d.sorted.nodup.non-chrM.tn5.bed.gz' % (expt_prefix, 2)
pooled_bedgz = ' '.join([R1_bedgz, R2_bedgz])

# Call peaks on both replicates and pooled reads
callpeak1 = cmd_callpeak('%s.%d' % (expt_prefix, 1),
                         R1_bedgz, 'R1', p_val, shift_val, ext_val)
callpeak2 = cmd_callpeak('%s.%d' % (expt_prefix, 2),
                         R2_bedgz, 'R2', p_val, shift_val, ext_val)
callpeakp = cmd_callpeak('%s.%s' % (expt_prefix, 'pooled'),
                         pooled_bedgz, 'pooled', p_val, shift_val, ext_val)
# print '\n'.join([callpeak1, callpeak2, callpeakp]) + '\n'

if __name__ == '__main__':
    pool = Pool(processes=3)
    pool.map(run_bash, [callpeak1, callpeak2, callpeakp])
    pool.close()


'''
Generate comparative bedgraph from treat_pileup and control_lambda,
trim the coordinates outside of chromosome size limits (due to MACS bug),
then convert bedgraph to bigWig. Delete the intermediate files.
'''

# Generate FE signal tracks
FE1 = cmd_FE(mm_chr_sz, 'R1', '%s.%d' % (expt_prefix, 1))
FE2 = cmd_FE(mm_chr_sz, 'R2', '%s.%d' % (expt_prefix, 2))
FEp = cmd_FE(mm_chr_sz, 'pooled', '%s.%s' % (expt_prefix, 'pooled'))
# print '\n'.join([FE1, FE2, FEp]) + '\n'

if __name__ == '__main__':
    pool = Pool(processes=3)
    pool.map(run_bash, [FE1, FE2, FEp])
    pool.close()

'''
Since --SPMR was used during peak calling, to generate the p-value track,
we must apply a scaling factor of the # of reads per million of the input
BED files. (Normally, we would take the lower of this value for treatment
and control, but there is no control for ATAC-seq.)

Use flagstat to count the number of total mapped reads in the non-chrM bam
file. For the pooled peaks, add the # of reads for each replicate. (An
alternate way to calculate this value is to simply count the number of lines
in the compressed BED file; zcat * | wc -l)
'''

# Generate p-value signal tracks
pv1 = cmd_ppois(mm_chr_sz, 'R1', '%s.%d' % (expt_prefix, 1), counts[1])
pv2 = cmd_ppois(mm_chr_sz, 'R2', '%s.%d' % (expt_prefix, 2), counts[2])
pvp = cmd_ppois(mm_chr_sz, 'pooled', '%s.%s' %
                (expt_prefix, 'pooled'), counts[1] + counts[2])
# print '\n'.join([pv1, pv2, pvp]) + '\n'

if __name__ == '__main__':
    pool = Pool(processes=3)
    pool.map(run_bash, [pv1, pv2, pvp])
    pool.close()


# Remove temporary files e.g. bedgraphs
print remove_bdg(['R1','R2','pooled', 'pseudoreps'])
# run_bash(remove_bdg(['R1','R2','pooled', 'pseudoreps']))
