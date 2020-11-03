# Set up a conda environment with python 3.6 https://uoa-eresearch.github.io/eresearch-cookbook/recipe/2014/11/20/conda/

conda create -n SomaticHaplotype python=3.6 anaconda

# you may need to initialize conda and source your .bashrc or .bash_profile (etc.) or restart the shell
#conda init [your shell]

# install pysam https://pysam.readthedocs.io/en/latest/
# pysam packages requires python 3.6 or below
conda install -n SomaticHaplotype -c bioconda pysam

# install pyVCF https://pyvcf.readthedocs.io/en/latest/
conda install -n SomaticHaplotype -c bioconda pyvcf

# install igraph https://igraph.org/python/
conda install -n SomaticHaplotype -c conda-forge python-igraph

# install samtools and bcftools for ancestry data variant calling  
# https://samtools.github.io/bcftools/
conda install -n SomaticHaplotype -c bioconda samtools=1.9
conda install -n SomaticHaplotype -c bioconda bcftools

# if you will use Beagle/IBD tools, then install a compatible Java version 8 distribution
conda install -n SomaticHaplotype -c anaconda openjdk

# at start of each session, activate SomaticHaplotype environment
conda activate SomaticHaplotype

# to exit SomaticHaplotype (or any) environment, close your terminal, or 
#conda deactivate

# to delete SomaticHaplotype environment
#conda remove -n SomaticHaplotype --all

################################################################################
# create a test bam based on a short genomic region and index it
# -h includes header and reads; -b compresses output to bam format
#samtools view -h -b your.bam your_region > your_destination/your.new.bam
#samtools index your_destination/your.new.bam
################################################################################
