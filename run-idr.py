import numpy as np
import os
import subprocess
import sys
from multiprocessing import Pool, Queue, Process, Manager, cpu_count

# Set global parameters

if len(sys.argv) == 4:

    # Chromosome sizes for clipping and bigwig generation
    mm_chr_sz = sys.argv[3]

    # Specify experiment condition
    # location further processing and analysis for experiment
    expt_prefix = sys.argv[1]
    pipe_out_loc = sys.argv[2]

# If single argument '-m', use the following manually input
# parameters
#elif len(sys.argv) == 2 and sys.argv[-1] == '-m':
elif len(sys.argv) > 1 and sys.argv[1] == '-m':

    # Chromosome sizes for clipping and bigwig generation
    mm_chr_sz = '~/annotation/chromosomes/mm10/mm10.chrom.sizes'

    # Specify experiment condition
    # location further processing and analysis for experiment
    expt_prefix = ''
    pipe_out_loc = ''

# If no args or incorrect number input, break and print help
else:
    print 'Error: Missing arguments. Exiting.\n\n' + \
          'Usage: python run-idr.py [options] <expt prefix> <output dir> <genome sz>\n\n' + \
          'Options: \n' + \
          '  -m   : Use parameters manually specified in script.\n' + \
          '         This option will ignore following arguments.'
    exit()

# Print verbose params
print 'Running script: \'run-idr.py\'. \n\n' + \
      '- Experiment prefix      : %s\n' % expt_prefix + \
      '- Peak call output dir   : %s\n' % pipe_out_loc + \
      '- Genome size            : %s\n' % mm_chr_sz


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
    print(command)
    #subprocess.call(command, shell=True)
    # Set analysis output directory in run_bash before every command
    subprocess.call('cd %s; %s' %(pipe_out_loc, command), shell=True)
    return None

#
def get_numpeaks(file_loc):
    with open(file_loc) as in_file:
        for num, line in enumerate(in_file):
            pass
    return num + 1

# Remove control_lambda and treat_pileups, if no other bdgcmp operations
# will be performed
def remove_bdg(outdirs):
    cmd_rm = ['rm -f %s/*.bdg ' %outdir for outdir in outdirs]
    return '; '.join(cmd_rm)

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

# Get pseudoreplicates for R1, R2
def pseudoreps(prefix, n_lines):

    # Split the non-chrM.tn5.bed.gz alignment for current replicate
    # Cat and gzip each split as pseudorep 1 and 2 respectively, remove temp files
    cmd_split = pipe(['gzip -dc %s.sorted.nodup.non-chrM.tn5.bed.gz' %prefix,
               'shuf',
               'split -a 2 -d -l %d - temp_split_%s' %(n_lines, prefix) ])
    cmd_PR1 = pipe(['cat temp_split_%s00' %prefix,
                   'gzip -c > pseudoreps/%s.PR1.bed.gz' %prefix])
    cmd_PR2 = pipe(['cat temp_split_%s01' %prefix,
                   'gzip -c > pseudoreps/%s.PR2.bed.gz' %prefix])

    return '; '.join([cmd_split, cmd_PR1, cmd_PR2, 'rm temp_split_%s* ' %prefix])

# Combine into "pooled" pseudoreplicates
def pooled_pseudoreps(expt_prefix):

    cmd_pooled_PR1 = pipe([
            'gzip -dc pseudoreps/%s.1.PR1.bed.gz pseudoreps/%s.2.PR1.bed.gz' %(expt_prefix, expt_prefix),
            'gzip -c > pseudoreps/%s.pooled.PR1.bed.gz' %expt_prefix
        ])
    cmd_pooled_PR2 = pipe([
            'gzip -dc pseudoreps/%s.1.PR2.bed.gz pseudoreps/%s.2.PR2.bed.gz' %(expt_prefix, expt_prefix),
            'gzip -c > pseudoreps/%s.pooled.PR2.bed.gz' %expt_prefix
        ])
    return '; '.join([cmd_pooled_PR1, cmd_pooled_PR2])

#
def run_IDR(oracle_peaks, rep_peaks, idr_thrsh, idr_prefix):

    # Run IDR on (pseudo)replicate peak sets v pooled oracle set
    #--use-best-multisummit-IDR
    idr_cmd = 'idr --plot --use-best-multisummit-IDR ' + \
    '--soft-idr-threshold %2.2f ' %idr_thrsh + \
    '--rank signal.value ' + \
    '--output-file IDR/%s.idr_values ' %idr_prefix + \
    '--log-output-file IDR/%s.log ' %idr_prefix + \
    '--peak-list %s ' %oracle_peaks + \
    '--samples %s ' %rep_peaks

    # Get peaks that fit IDR threshold
    log_thrsh = int(round(-125*np.log2(idr_thrsh)))
    idr_peaks = pipe([
            'cat IDR/%s.idr_values' %idr_prefix,
            'awk \'BEGIN{OFS="\\t"} $5>=%d {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10}\'' %log_thrsh,
            'sort -k7n,7n > IDR/%s.narrowPeak' %idr_prefix
        ])

    return '; '.join([idr_cmd, idr_peaks])




'''
##########################################
# GENERATE PSEUDOREPLICATES AND CALL PEAKS
##########################################

Split each replicate's post-processed alignments (non-chrM.tn5.bed.gz) and shuffle, then split
into an equal number of lines, (num_reads + 1)/2 to generate pseudoreplicates of each true replicate.
Combine the pseudoreplicates together to create the pooled pseudoreplicates.

Then, call peaks on all six pseudoreplicates:
- replicate 1: PR1, PR2
- replicate 2: PR1, PR
- pooled: PR1, PR2
'''

counts = cmd_readcounts(os.path.expanduser(pipe_out_loc))
print(counts, '\n')

run_bash('if ! [ -d pseudoreps ]; then mkdir pseudoreps; fi')

# Generate pseudoreplicates
cmd_pseudo_R1 = pseudoreps('%s.%d' % (
    expt_prefix, 1), (int(counts[1] * 1000000) + 1) / 2)
cmd_pseudo_R2 = pseudoreps('%s.%d' % (
    expt_prefix, 2), (int(counts[2] * 1000000) + 1) / 2)
cmd_pseudo_pooled = pooled_pseudoreps(expt_prefix)
# print '\n'.join([cmd_pseudo_R1, cmd_pseudo_R2, cmd_pseudo_pooled]) + '\n'

if __name__ == '__main__':
    pool = Pool(processes=2)
    pool.map(run_bash, [cmd_pseudo_R1, cmd_pseudo_R2])
    pool.close()
run_bash(cmd_pseudo_pooled)


''' '''
# Params for peak calling
shift_val = '37'
ext_val = '73'
p_val = '1e-2'

# Call peaks on each pseudoreplicate
callpeak_pseudo_cmds = []
for curr_rep in ['1', '2', 'pooled']:
    PR1 = 'pseudoreps/%s.%s.PR1.bed.gz ' % (expt_prefix, curr_rep)
    PR2 = 'pseudoreps/%s.%s.PR2.bed.gz ' % (expt_prefix, curr_rep)
    callpeak_pseudo_cmds.append(
        cmd_callpeak('%s.%s.PR1' % (expt_prefix, curr_rep),
                     PR1, 'pseudoreps', p_val, shift_val, ext_val)
    )
    callpeak_pseudo_cmds.append(
        cmd_callpeak('%s.%s.PR2' % (expt_prefix, curr_rep),
                     PR2, 'pseudoreps', p_val, shift_val, ext_val)
    )
# print '\n'.join(callpeak_pseudo_cmds) + '\n'

if __name__ == '__main__':
    pool = Pool(processes=6)
    pool.map(run_bash, callpeak_pseudo_cmds)
    pool.close()


'''
###############################################
# RUN IDR AND CALCULATE REPRODUCIBILITY METRICS
###############################################

The following four IDR analyses will be performed, using the indicated peak sets
as samples" and "peak-list" respectively:

1. Replicate 1, self-reproducibility
   Rep1.PR1, Rep1.PR2; vs Rep1

2. Replicate 2, self-reproducibility
   Rep1.PR2, Rep2.PR2; vs Rep2

3. True replicate reproducibility
   Rep1, Rep2; vs Pooled

4. Pooled pseudoreplicate reproducibility
   Pooled.PR1, Pooled.PR2; vs Pooled

After filtering by IDR, parse the output peaks by log-transformed IDR score to
select retain only the peaks that pass the IDR cutoff threshold, and cut the
additional columns beyond standard narrowPeak format.
'''

run_bash('if ! [ -d IDR ]; then mkdir IDR; fi')

# Set IDR threshold, default = 0.05
idr_thrsh = 0.05

# Run IDR filtering on all four tests
idr_cmds = []

# Replicate 1, self-reproducibility
R1_oracle_peaks = 'R1/%s.1_peaks.narrowPeak' % expt_prefix
R1_rep_peaks = ' '.join([
    'pseudoreps/%s.1.PR1_peaks.narrowPeak' % expt_prefix,
    'pseudoreps/%s.1.PR2_peaks.narrowPeak' % expt_prefix
])
R1_idr_prefix = 'rep1.self-pseudorep'
idr_cmds.append(run_IDR(R1_oracle_peaks, R1_rep_peaks,
                        idr_thrsh, R1_idr_prefix))
# Replicate 2, self-reproducibility
R2_oracle_peaks = 'R2/%s.2_peaks.narrowPeak' % expt_prefix
R2_rep_peaks = ' '.join([
    'pseudoreps/%s.2.PR1_peaks.narrowPeak' % expt_prefix,
    'pseudoreps/%s.2.PR2_peaks.narrowPeak' % expt_prefix
])
R2_idr_prefix = 'rep2.self-pseudorep'
idr_cmds.append(run_IDR(R2_oracle_peaks, R2_rep_peaks,
                        idr_thrsh, R2_idr_prefix))
# True replicate reproducibility
TR_oracle_peaks = 'pooled/%s.pooled_peaks.narrowPeak' % expt_prefix
TR_rep_peaks = ' '.join([
    'R1/%s.1_peaks.narrowPeak' % expt_prefix,
    'R2/%s.2_peaks.narrowPeak' % expt_prefix
])
TR_idr_prefix = 'pooled.truerep'
idr_cmds.append(run_IDR(TR_oracle_peaks, TR_rep_peaks,
                        idr_thrsh, TR_idr_prefix))
# Pooled pseudoreplicate reproducibility
PR_oracle_peaks = 'pooled/%s.pooled_peaks.narrowPeak' % expt_prefix
PR_rep_peaks = ' '.join([
    'pseudoreps/%s.pooled.PR1_peaks.narrowPeak' % expt_prefix,
    'pseudoreps/%s.pooled.PR2_peaks.narrowPeak' % expt_prefix
])
PR_idr_prefix = 'pooled.pseudorep'
idr_cmds.append(run_IDR(PR_oracle_peaks, PR_rep_peaks,
                        idr_thrsh, PR_idr_prefix))

# print('\n'.join(idr_cmds) + '\n')

if __name__ == '__main__':
    pool = Pool(processes=4)
    pool.map(run_bash, idr_cmds)
    pool.close()


'''
Using the IDR filtered final peak sets, calculate reproducibility metrics and
choose an 'optimal' set for further analysis. The 'conservative' peak set is the
IDR filtered true replicates set; however, if the number of peaks in the IDR
pooled pseudoreplicates is greater, this will be chosen as the optimal set.

There are two reproducibility metrics: rescue ratio, and self-consistency-ratio.
Rescue ratio is the ratio between the number of peaks in IDR filtered TR and PR
sets. Self consistency ratio is the ratio between the number of peaks in the IDR
filtered R1 and R2 sets.  If BOTH ratios are <2, then the reproducibility test
criteria are passed.
'''

# Get numpeaks yo
n_TR = get_numpeaks(os.path.expanduser(
    '%s/IDR/pooled.truerep.narrowPeak' % pipe_out_loc))
n_PR = get_numpeaks(os.path.expanduser(
    '%s/IDR/pooled.pseudorep.narrowPeak' % pipe_out_loc))
n_R1 = get_numpeaks(os.path.expanduser(
    '%s/IDR/rep1.self-pseudorep.narrowPeak' % pipe_out_loc))
n_R2 = get_numpeaks(os.path.expanduser(
    '%s/IDR/rep2.self-pseudorep.narrowPeak' % pipe_out_loc))

# Choose optimal set
if n_TR >= n_PR:
    opt_result = 'Optimal set is IDR True Replicates.'
else:
    opt_result = 'Optimal set is IDR Pooled Pseudoreplicates'

# Calculate reproducibility metrics
rescue_ratio = float(max(n_PR, n_TR)) / float(min(n_PR, n_TR))
self_consistency_ratio = float(max(n_R1, n_R2)) / float(min(n_R1, n_R2))
if rescue_ratio > 2 and self_consistency_ratio > 2:
    test_result = 'Reproducibility Test: Fail'
elif rescue_ratio > 2 or self_consistency_ratio > 2:
    test_result = 'Reproducibility Test: Borderline'
else:
    test_result = 'Reproducibility Test: Pass'
ratio_results = 'Rescue ratio: %s, Self-Consistency Ratio: %s' % (
    rescue_ratio, self_consistency_ratio)

# Save results output to IDR folder
print('\n'.join([opt_result, test_result, ratio_results]))
with open(os.path.expanduser('%s/IDR/IDR_results.txt' % pipe_out_loc), 'w') as idr_results:
    idr_results.write('\n'.join(['n_TR, n_PR, n_R1, n_R2', '%d, %d, %d, %d' % (
        n_TR, n_PR, n_R1, n_R2), opt_result, test_result, ratio_results]))


# Remove temporary files e.g. bedgraphs
# print(remove_bdg(['R1','R2','pooled', 'pseudoreps']))
run_bash(remove_bdg(['R1','R2','pooled', 'pseudoreps']))
