#######################################
# Download beagle, Refined IBD, merge #
#######################################

tools_dir=tools
mkdir -p $tools_dir

# Get beagle.21Sep19.ec3.jar
wget https://faculty.washington.edu/browning/beagle/beagle.21Sep19.ec3.jar -P $tools_dir

# Get Refined IBD
wget http://faculty.washington.edu/browning/refined-ibd/refined-ibd.16May19.ad5.jar -P $tools_dir

# Get merge IBD
wget https://faculty.washington.edu/browning/refined-ibd/merge-ibd-segments.16May19.ad5.jar -P $tools_dir

wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples_v3.20130502.ALL.panel
