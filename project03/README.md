# Introduction
This project implements a Gibbs Sampler, a stochastic (randomized) algorithm used to identify conserved DNA patterns (motifs) within a set of sequences. The goal is to discover biological signals that are hidden within genomic noise. It does this by iteratively refining motif position predictions by building probabilistic models from sequence alignments.

The Evolution of the Project:
Initially, this project focused on identifying the binding site for the p53 tumor suppressor protein using human ChIP-seq data and then later refining the data using peak calling, but due to unsuccessful attempts, we switched to a bacterial dataset to verify the accuracy of our algorithm

# BAM to FASTA Workflow  
```
# Get peaks from bam file
macs2 callpeak -t SRR9090854.subsampled_5pct.bam -f BAMPE -n p53_study -g hs -q 0.01

# Take peak summits and expand start and end positions by 50bp in each direction
bedtools slop -i p53_study_summits.bed -g <List of chromosomes and corresponding sizes in txt file> -b 50 > p53_study_peaks100.bed

# Convert chromosome names in bed file to agree with NCBI chromosome naming convention in GRCh38.p14 FASTA
awk 'BEGIN {OFS="\t"} NR==FNR {map[$1]=$2; next} \
{chr=$1; gsub(/^chr/, "", chr); if (chr in map) {$1=map[chr]; print}}' \
bed/chr_mapping.txt p53_study_peaks_100.bed > output.bed

# Convert bed file to FASTA based on NCBI GRCh38.p14 FASTA
bedtools getfasta -fi <Human Genome FASTA file> 
```
# Pseudocode
Pseudocode for processing the Bacterial Dataset

```
1. INITIALIZE ENVIRONMENT
   - Load NumPy, Seqlogo, and parser functions
   - Load utility functions for k-mer extraction and other helper functions

2. PRE-PROCESS DATA
   - Open bacterial genomic FASTA and GFF files.
   - For each entry marked as a CDS:
     - Extract the 50bp sequence immediately upstream of the promoter.
     - Check if sequence contains "AGGAGG".
     - Add matching sequences to a master list.

3. INITIALIZE GIBBS SAMPLER
    - For every sequence in the list:
    - Randomly select a 10bp window motif as a starting point.
    - Define background noise to ignore.

4. SAMPLING LOOP for 50,000 Iterations
   - Pick one sequence from the list at random.
   - Temporarily remove its motif from the current model.
   - Build a Position Frequency Matrix (PFM) from all other sequences.
   - Convert PFM to a Weight Matrix (PWM) with pseudocounts.
   - Generate all possible 10bp k-mers in the picked sequence.
   - Score each window against the PWM.
   - Probabilistically select a new window based on its score using the log2 transform.
   - Update the master list with the new window.

5. OUTPUT & VISUALIZATION
   - Take the final PFM and transpose it.
   - Convert counts into probabilities for normalization.
   - Scale the height of letters based on Information Content.
   - Render and display the final Sequence Logo.
```

# Successes
- Successfully implemented and validated our Gibbs Sampler by switching from the p-53 dataset to the bacterial dataset
- Overcame long run times and debugging difficulties by implementing smaller text sequences and fewer iterations to verify code behavior before scaling up to larger datasets
- Was able to implement extraction from raw sequences (FASTA/GFF) and was able to use peak calling (MACS2) for p-53 dataset, although decided not to go with this due to A-T biases
- Correctly implemented complex data transformations, such as converting logarithmic PWM scores into linear weights to enable probabilistic motif selection.
- Successfully navigated significant dependency issues with the seqlogo package and generated a meaningful seqlogo result that matches our expectations

# Struggles
- Handling large datasets and determining how to approach testing and building models due to extremely long run times.
- Understanding data transformations
- Seqlogo dependencies
- Faced significant issues installing and configuring MACS2 and its associated dependencies within the local environment
- Struggled to define peak parameters that were specific enough to capture binding sites without being so broad that they introduced excessive genomic noise.
- Formatting data after peak calling
- Encountered persistent AT-rich biases in the p53 dataset, where the model converged on high-frequency A/T regions rather than the biological RRRCWWGYYY motif.
- challenging to maintain the state of the data through various mathematical transformations, particularly when reversing scores to weights for the sampling step.

# Personal Reflections
## Group Leader
Fardina Tabassum- It was a bit confusing at first to figure out which files to work with and what sort of information each file has, and what it is used for, but after reviewing the lecture and consulting with my group members, I was able to learn about them. It was also hard to understand how many times to iterate so that the program would not be computationally intensive. We also had issues displaying the seqlogo. It was my first time working with a large dataset, so I found it a bit challenging as the run times were very long, and it took a lot of time to debug the code. We combatted this by using text sequences and fewer iterations to check if the code was working first. Working with MACS2 was a bit frustrating for me to install at first, but my teammates helped me figure it out. We made multiple attempts using the p-53 dataset to generate motifs at first, but we were unable to resolve seeing AT biases in our reads. We then decided to switch over to the bacteria dataset, and we found out that our sampler does indeed work so that was a relief. Overall, my team was very accommodating to my lack of understanding and helped me with a lot of the code clarification. They were able to dedicate a good amount of time to this project, trying different code logic to get the dataset to work. We were able to meet and collaborate engagingly so I am thankful for having such a dedicated team. This project really tested my understanding of algorithms and made me realize how important it is to understand the biological context as well. It also gave me an insight into how researchers work when they have unknown datasets and how frustrating it is, as we don't know what result to expect. 

## Other member
Nicholas Bottomley- This assignment came with many new hurdles that we had not previously had much experience with. The major one being handling large datasets. It became apparent very quickly that it was not feasible to test and build our model if we were to run it on the full dataset every time that we wanted to check how it was behaving. To overcome this, we took a subset of the main dataset, so that we could see how the model would behave in a more reasonable amount of time. Part of building the model was also keeping track of the state of data and what transformations it had undergone in case those transformations had to be reversed, which was the case when converting the pwm-derived scores into weights for determining a new motif to use in the pfm. Once we had a functioning algorithm, we still faced the challenge of how to handle such a large dataset, much of which contained duplicative information. This is where we were introduced to the method of peak calling. Peak calling tools are very useful for this exact purpose because they can simplify broad sets of genomic data into just the regions with the most coverage. There were some difficulties getting MACS2 installed, but it was successfully installed and run after some troubleshooting. It was at this point that we encountered our last major hurdle, trying to get good peaks which would contain the p53 binding sites we were looking for without being overly broad, as that would defeat the purpose of finding the peaks. An unforeseen issue that came with this peak calling was AT-rich regions around the peaks, which would lead the model to converge to a strong A or T bias, clearly not reflecting the motif that we were searching for. Outside of the computational realm, we also faced issues with properly installing and enabling dependencies for the seqlogo package, with most errors seeming to arise from problems with ghostscript permissions and fonts that are required for creation of a sequence logo.

Tien Nguyen- In this project, I implemented the Gibbs Motif Finder algorithm, which challenged me to understand the basic principles of motif discovery and how the provided functions worked together. Although it was difficult at first, the process helped strengthen my problem-solving and algorithmic thinking skills. After being given extra time, I explored ChIP-seq analysis using MACS2, bedtools, and MEME to identify DNA motifs from sequencing data. This experience improved my confidence in using bioinformatics tools and helped me connect computational methods with biological interpretation.

# Generative AI Appendix
Anthropic. (2025). Claude Sonnet 4.6 [Large language model]. https://www.anthropic.com  
Claude was used to understand the biological context behind the datasets. Additionally, Claude was used to help with install and use of file conversion packages, such as macs2 and bedtools.
