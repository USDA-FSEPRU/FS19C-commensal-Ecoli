# Pre-Methods

Written summary of methods performed prior to analyzing complete set of FS19C sequence data to practice using tools (bbmap). Includes lab notes of how methods performed. Methods were copied from `OneNote/FS19C/FS19C lab notebook` entries

## 17Aug2020 - 14Sept2020: Review bbmap scripts on Daniel's 10 E. coli isolates

**17Aug2020** Talked with Daniel on the phone who answered my questions:
1. Do you know where the DNA for the 10 in-house E. coli isolates are stored? In one of the -20C freezers near large autoclave room in Vijay's lab. DNA are stored on a rack.

2. Is the Qubit readings for the DNA written down somewhere? Unknown

3. Where are the stocks for the isolates? Vijay already answered

4. Where's the long read sequence data for the 10 E. coli isolates stored at? On your account on Ceres?
  * Daniel already has short-read sequence data. He ran MiSeq specifically for the 10 isolates.
  * Find the assembled sequences on SciNet
  * The original sequence data on Q drive: David Alt --> March 2020, DNielsen or Nielsen
	* In I drive (department drive): Daniel's folder --> Ecoliassemblyinfo: spreadsheet that tells you about assemblies of the 10 isolates and the 4 challenge isolates from FS19. Tells you about colicins, replicons, # contigs
    * I drive --> Nielsen  --> Work Laptop --> FS19 --> `EcoliAssemblyDocuments.xlsx`
    * xlsx info: assembly, careful # contigs, isolate # contigs, core peptide, core immunity genes, bottromycin, replicons, other, # polishing iterations

5. Are these the 10?
```
1803
2397
2880
2936
2938
429
5523
6673
6677
6681
```		
  * May take some hunting
  * Email: daniel.nielsen@uspto.gov

6. Is it your GitHub or wherever that can tell me what you've done so far on the 10 isolates? Already answered above

**27Aug2020**
1. I found Daniel's list of isolates on Q drive and only found three of the 10 isolates that Vijay showed me.
  * Q drive > DAlt > Illumina MiSeq Output > 03-02-20 DNielsen

2. Downloaded the 10 isolates to local Mac desktop

3. Will try practicing Jules' BBMap pipeline later

**10Sept2020**
1. Transferred *.fastq.gz files from local desktop to scinet folder
```
scp *.gz scinet:~/fsepru_kmou/EcoliGenomesPractice/
```

2. Downloaded latest Bbmap tools from: https://sourceforge.net/projects/bbmap/

3. Untar bbmap file

4. Transferred bbmap directory to scinet EcoliGenomesPractice folder

**11Sept2020** Emailed Jules about how to get bbmap started. He sent me slurm template and directions. I need to parse through.
His response:
```
Keep in mind that there are many ways to accomplish this, but what follows is how I like to do this kind of task.
 
Here are the general steps I like to take:
1)     Make a list of sample names
2)     make a template script that has the word “REPLACE” where you would like to swap in the sample names
3)     use sed to replace the word REPLACE in your template with your sample name
 
See the attached template script as an example of what I have done in the past.
The file ‘samples.txt’ should be a text file that contains one sample name per line.
 
Here is the bash code that will generate a slurm script for each of your samples.  (the double quotes here are important, single quotes wont work).
 
while read line
do
cat SRAassemblyPipeline.SLURM_TEMPLATE | sed “s/REPLACE/$line/g” > “$line”.slurm
done < samples.txt
```

**14Sept2020**
1. Made text list of sample names -- do I keep fastq.gz or just the part before "fastq.gz"

2. In Dnielsen directory, on terminal, write:
```
ls >> samples.txt
```

3. Opened samples.txt and manually changed names of samples to:
```
1803
2397
2880
2936
2938
429
5523
6673
6677
6681
```

4. I had asked Jules the following questions:
  * When I'm making the list of samples, how much of the sample's file name should I keep? For example, the file `429_S6_L001_R1_001.fastq.gz` - should I list this as 429 and ignore everything else after that? 
	* I was reading through the template script and for line 22, when it says `in1=REPLACE_1.fastq.gz`, does this mean it will find, for example, a file called `429_S6_L001_R1_001.fastq.gz` and it's counterpart `429_S6_L001_R2_001.fastq.gz` for `in2=REPLACE_2.fastq.gz` and replace "REPLACE" with 429 from samples.txt?
	* On line 19, what does source do? What should I put instead? Working directory? I looked through my bbmap directory and didn't find templates directory or `adapt_polish.sh`.
	* Once I have all that figured out, then I can make the bash script you wrote at the end of the last email, save as whatever_batch.sh, run sbatch whatever_batch.sh, and it'll do its magic?

5. His responses:
```
(DONE)1)      I tend to remove all of the extra information from the illumina output files, but that choice is up to you.  I guess I forgot to include that my first step is to rename the R1 and R2 files to only contain the sample ID information, and not any of the S*_L00*_ etc info.  Again this is just how I do stuff, there is no ‘correct’ way to do this.
(DONE) 2)      Whenever you see the word “REPLACE” in my template script it will be replaced by exactly the same word that is contained in each line of your samples.txt file.  It is up to you to make sure that replacement points to a legitimate file.
(DONE) 3)      adapt_polish.sh is a little bash script that contains a function I wrote to help automate polishing steps with pilon.  You don’t necessarily need to do this, it’s just an extra step I included.  The “source” command essentially executes that script which will make the function available for use later in the script.  For now I recommend commenting these lines out and just focus on getting the QC and assembly steps to work.  Once you have that we can add other extra things back in.
4)     “while read XXX” is how you initiate  a loop in bash that will do their actions once for each line in a file.  At the very end of this loop you can see a “ done < samples.txt”  this tells bash to use the samples.txt file to supply the lines that the loop will read.  So this loop will output a slurm script for each sample in your samples.txt file.  You can then submit those slurm scripts to the queue.
https://linuxize.com/post/bash-while-loop/#read-a-file-line-by-line
```

6. Completed the following that Jules described above and ran `1803.slurm`.

7. Slurm script failed because it couldn't find reformat.sh script. I asked Jules how to tell where to locate script. He said the following:
```
You have downloaded bbtools but you have not informed your system where to find these scripts.  You need to tell your system where they are located.
There are two ways to accomplish this.
1)      The bad way:
				a.      Whenever you call these commands, type out the full path to the location of each.
				b.      For example: /home/Kathy.mou/software/bbmap/reformat.sh in1=… 
2)      The better way:
				a.      You need to add the location of the bbmap scripts to your “path”
				b.      A path is a list of directories where the computer will look for executables
You add things to your path by editing your “.bash_profile” file.  This is a hidden file located in your home directory.
I have attached my .bash_profile file so you can see an example of how to do this.
Basically PATH is a bash variable of colon delimited paths, so you need to add the location of the bbtools scripts to this variable, separated by a colon.
```

8. Need to parse through Jules' script

## 15Sept2020 - 2020: Practice bbmap scripts on 4 FS19C E. coli isolates to determine whether to run with 100X or 250X coverage assemblies

**15Sept2020**
1. Got email from Darrell that NovaSeq data is available. I am downloading to local computer.

2. On SciNet, I ended up moving bbmap from EcoliGenomesPractice to home directory
```
cd ~/bbmap
du -sh #look at how big folder is
vi ~/.bash #then press tab twice to see what other .bash* files there are
vi ~/.bash_profile
```

3. Edited `.bash_profile` in vi by adding ":$HOME/bbmap" to this line:
```
PATH=$PATH:$HOME/bin:$HOME/bbmap
source ~/.bash_profile
echo $PATH
#Output: /usr/lib64/qt-3.3/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/kathy.mou/bin:/home/kathy.mou/bin:/home/kathy.mou/bbmap
cd fsepru_kmou/EcoliGenomesPractice/
sbatch 1803.slurm
```

4. Made sure I corrected the path to bbmap scripts within 1803.slurm

5. It couldn't run `spades.py` because I didn't have the executable available.

6. Moved bbmap to directory "software" in home directory

7. Download Spades binaries for Mac from: http://cab.spbu.ru/files/release3.13.0/manual.html#sec2

8. Transfer `Spades-3.13.0-Darwin` folder to SciNet:
```
scp -r Spades* scinet:~/software
```

9. Repeated editing bash_profile to include spades
```
PATH=$PATH:$HOME/bin:$HOME/software/SPAdes-3.13.0-Darwin/bin
```

10. I ran source and echo

11. Tried to run `spades.py` but still no luck. Will try running sample SLURM job following this: https://bioinformaticsworkbook.org/dataAnalysis/GenomeAssembly/Assemblers/spades.html#gsc.tab=0

**16Sept2020**
1. Jules helped to transfer NovaSeq data from Q drive to my fsep/kmou directory on Ceres.

2. Untar file
```
tar -xvf 1_33298_01_1-A01-1-428RN3A_HTFVy_1801.tar
```

3. Make `samples.txt` locally and upload to Ceres

4. Read through Jules' SRAAssemblyPipeline scripts and compare with `assemblyPipeline.sh` and other scripts, bbmap text files
  * `Reformat.sh`
    * Question for Jules about sbt=500000000 parameter
  * `Bbduk.sh`
    * Leave as is
  * `Bbmerge.sh`
    * Leave as is
  * `Clumpify.sh`
    * Leave as is
  * `Tadpole.sh`
    * Leave as is
5. Question about interleaved=t

6. Question about `spades.py` t, m, k settings

7. Emailed Jules my questions

8. Downloaded spades v 3.14.1 for Linux: https://cab.spbu.ru/files/release3.14.1/manual.html

9. Uploaded to Ceres ~/software directory. Untar file. Updated `.bash_profile` to include correct path and spades version
```
PATH=$PATH:$HOME/bin:$HOME/software/SPAdes-3.14.1-Linux/bin
```

10. Needed to source `.bash_profile` and `echo $PATH`

11. Ran spades.py in 1803.slurm and seems to be working
```spades.py --12 1803_ecct.fq.gz -o 1803_spades_out --only-assembler -t 16 -m 32 -k 25,55,95,125
```
**17Sept2020**
1. Analyzed spades from stdout.5099458.ceres19-compute-62
```
* Assembled contigs are in /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/contigs.fasta
 * Assembled scaffolds are in /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/scaffolds.fasta
 * Paths in the assembly graph corresponding to the contigs are in /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/contigs.paths
 * Paths in the assembly graph corresponding to the scaffolds are in /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/scaffolds.paths
 * Assembly graph is in /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/assembly_graph.fastg
 * Assembly graph in GFA format is in /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/assembly_graph_with_scaffolds.gfa
======= SPAdes pipeline finished.
SPAdes log can be found here: /lustre/project/fsepru/kmou/EcoliGenomesPractice/1803_spades_out/spades.log
```

2. Question about setting sbt parameter - see Jules' email. Also emailed David Alt

**7Oct2020**
1. Because there's still general issues logging in to Ceres (this was posted on SCINet basecamp message board since yesterday Oct. 6) so I asked Jules if he could transfer novaseq data to Atlas in fsepru project directory. He completed that with no issues.

2. Made kathy.mou directory in fsepru directory, transferred novaseq file to my directory

3. Untar novaseq file using `tar -xvf <file>`

4. Transferred `samples.txt` (only lists 4 samples) from local computer to Atlas: `/project/fsepru/kathy.mou/lane1`

5. Downloaded latest Bbmap tools from: https://sourceforge.net/projects/bbmap/

6. Transferred bbmap directory to Atlas in home directory under software directory. Untar file

7. Scp FS19C SLURM script to lane1 folder - script set up for 100X coverage

8. Copy the four fastq files, SLURM script, `samples.txt` to `FS19C_4Samples`

9. Found out Ceres is back online so I will pause what I'm doing on Atlas and go back to Ceres.

10.
```
scp SRAassemblyPipeline.FS19C.SLURM_TEMPLATE ceres:/project/fsepru/kmou/FS19C/
```

11. Scp samples.txt (only 4 samples) from local computer to ceres
```
1-A01-1-428RN3A_S8_L001_R1_001.fastq.gz
1-A01-1-428RN3A_S8_L001_R2_001.fastq.gz
1-B08-20-427FEC_S27_L001_R1_001.fastq.gz
1-B08-20-427FEC_S27_L001_R2_001.fastq.gz
1-H10-94-439FED_S101_L001_R1_001.fastq.gz
1-H10-94-439FED_S101_L001_R2_001.fastq.gz
1-H12-96-441FEC_S103_L001_R1_001.fastq.gz
1-H12-96-441FEC_S103_L001_R2_001.fastq.gz
```

12. Copy replace.batch from `EColiGenomesPractice` directory to `FS19C_4Samples` directory. Create SLURM scripts for four files
```
bash replace.batch
```

13. Run at 100X coverage
```
sbatch 1-428RN3A.slurm
Submitted batch job 5138118
```
  * Job completed
  * `Covstats.text`
```
sbatch 20-427FEC.slurm
Submitted batch job 5138138
```
  * Job completed WITH WARNINGS… Ask Jules??
```
  === Error correction and assembling warnings:
	 * 0:00:26.748   161M / 11G   WARN    General                 (pair_info_count.cpp       : 341)   Unable to estimate insert size for paired library #0
	 * 0:00:26.748   161M / 11G   WARN    General                 (pair_info_count.cpp       : 347)   None of paired reads aligned properly. Please, check orientation of your read pairs.
	 * 0:00:26.750   161M / 11G   WARN    General                 (repeat_resolving.cpp      :  63)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
```

```
sbatch 94-439FED.slurm
Submitted batch job 5138140
```
  * Job completed
  * `Covstats.txt`

```
sbatch 96-441FEC.slurm
Submitted batch job 5138141
```
  * Job failed
  * File is empty (in stdout file)

14. Move output files to a different folder `FS19_4Samples100`

15. Rerun replace.batch to generate new slurm scripts with 250X coverage (5500000*250=1375000000 )

16. Run at 250X coverage
```
Sbatch 1-428RN3A.slurm
Submitted batch job 5138154
```
* Job completed

```
sbatch 20-427FEC.slurm
Submitted batch job 5138157
```
  * Job completed WITH WARNINGS!
```
  === Error correction and assembling warnings:
		 * 0:00:31.862   157M / 11G   WARN    General                 (pair_info_count.cpp       : 341)   Unable to estimate insert size for paired library #0
		 * 0:00:31.863   157M / 11G   WARN    General                 (pair_info_count.cpp       : 347)   None of paired reads aligned properly. Please, check orientation of your read pairs.
		 * 0:00:31.864   157M / 11G   WARN    General                 (repeat_resolving.cpp      :  63)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
```		

```
sbatch 94-439FED.slurm
Submitted batch job 5138158
```
  * Job completed

```
sbatch 96-441FEC.slurm
Submitted batch job 5138159
```
  * Job failed
  * File is empty (in stdout file)

17. Compile results of assemblies by making table, email Jules

| Sample File | 100X | Any missing output? Any extra output? | 250X | Any missing output? Any extra output? |
| --- | ---- | ---------- | -- | ------|
| 1-428RN3A |job completed | NA| job completed| No |
| 20-427FEC |	job completed with WARNINGS |	Missing: contigs.paths, scaffolds.fasta, scaffolds.paths<br /> 	Extra: warnings.log		| job completed with WARNINGS |	Missing: contigs.paths, scaffolds.fasta, scaffolds.paths<br />Extra: warnings.log
| 94-439FED	| job completed |	No |	job completed | No|
| 96-441FEC	| job failed	| Empty directory	| job failed |	Empty directory |

18. Are the warnings the same for #20 regardless of coverage? Yes

19. Contents in `FS19C_4Samples100X/1-428RN3A_spades_out/`
```
assembly_graph_after_simplification.gfa
assembly_graph.fastg
assembly_graph_with_scaffolds.gfa
before_rr.fasta
contigs.fasta
contigs.paths
dataset.info - nothing much there (just path name for reads)
input_dataset.yaml
K125
```
* Within the K files are:
```
assembly_graph_after_simplification.gfa
assembly_graph.fastg
assembly_graph_with_scaffolds.gfa
before_rr.fasta
configs/ # directory
		bgc_mode.info          detail_info_printer.info  meta_mode.info      simplification.info
		careful_mda_mode.info  distance_estimation.info  moleculo_mode.info  toy.info
		careful_mode.info      isolate_mode.info         pe_params.info      tsa.info
		config.info            large_genome_mode.info    plasmid_mode.info
		construction.info      mda_mode.info             rna_mode.info
final_contigs.fasta
final_contigs.paths
final.lib_data
path_extend/ # empty directory
scaffolds.fasta
scaffolds.paths
```
* Other K directories that have the same output as K125:
```
K25
K55
K95
misc/broken_scaffolds.fasta
params.txt
pipeline_state/ # directory
	stage_0_before_start  stage_10_as_finish  stage_11_bs  stage_12_terminate  stage_1_preprocess_start  stage_2_preprocess_12  stage_3_preprocess_finish  stage_4_as_start  stage_5_k25  stage_6_k55  stage_7_k95  stage_8_k125  stage_9_copy_files
run_spades.sh
run_spades.yaml
scaffolds.fasta
scaffolds.paths
spades.log
split_input/ # directory
	1-428RN3A_ecct_1.fastq  1-428RN3A_ecct_2.fastq
tmp directory/ # empty directory
```

20. Read Spades manual for what output files mean: https://cab.spbu.ru/files/release3.13.0/manual.html#sec3.5
```
<output_dir>/scaffolds.fasta # contains resulting scaffolds (recommended for use as resulting sequences)
<output_dir>/contigs.fasta # contains resulting contigs
<output_dir>/assembly_graph.gfa # contains SPAdes assembly graph and scaffolds paths in GFA 1.0 format
<output_dir>/assembly_graph.fastg # contains SPAdes assembly graph in FASTG format
<output_dir>/contigs.paths # contains paths in the assembly graph corresponding to contigs.fasta (see details below)
<output_dir>/scaffolds.paths # contains paths in the assembly graph corresponding to scaffolds.fasta (see details below)
```
21. For later: Read Jules' email, troubleshoot. Read through NARMS resources, standard E. coli genome size: 6MB

**19Oct2020 and 20Oct2020**

1. Troubleshoot warning messages for 20-427FEC
```
Unable to estimate insert size for paired library (pair_info_count.cpp)
Non of paired reads aligned properly. Check orientation of read pairs (pair_info_count.cpp)
Insert size was not estimated for any of the paired libraries, repeat resolution module will not run (repeat_resolving.cpp)
```
2. Searched `pair_info_count.cpp` in `spades.log` and stdout (`vi stdout.5138138.ceres14-compute-65`) with vi and saw the same following logs:
100X
```
 0:00:23.973   165M / 11G   INFO    General                 (pair_info_count.cpp       : 322)   Min edge length for estimation: 87821
  0:00:23.973   165M / 11G   INFO    General                 (pair_info_count.cpp       : 333)   Estimating insert size for library #0
  0:00:23.974   165M / 11G   INFO    General                 (pair_info_count.cpp       : 190)   Estimating insert size (takes a while)
  0:00:24.071   438M / 11G   INFO    General                 (pair_info_count.cpp       :  39)   Selecting usual mapper
  0:00:25.525   434M / 11G   INFO    General                 (sequence_mapper_notifier.h:  94)   Total 1807994 reads processed
  0:00:26.742   434M / 11G   INFO    General                 (pair_info_count.cpp       : 208)   Edge pairs: 4521
  0:00:26.742   434M / 11G   INFO    General                 (pair_info_count.cpp       : 212)   0 paired reads (0% of all) aligned to long edges
  0:00:26.748   161M / 11G   WARN    General                 (pair_info_count.cpp       : 341)   Unable to estimate insert size for paired library #0
  0:00:26.748   161M / 11G   WARN    General                 (pair_info_count.cpp       : 347)   None of paired reads aligned properly. Please, check orientation of your read pairs.
 0:00:26.749   161M / 11G   INFO   StageManager             (stage.cpp                 : 166)   STAGE == Distance Estimation
  0:00:26.749   161M / 11G   INFO    General                 (distance_estimation.cpp   : 183)   Clearing raw paired index
  0:00:26.749   161M / 11G   INFO   StageManager             (stage.cpp                 : 166)   STAGE == Repeat Resolving
  0:00:26.750   161M / 11G   WARN    General                 (repeat_resolving.cpp      :  63)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
```

250X (`spades.log` and `vi stdout.5138157.ceres14-compute-65`)
```
0:00:28.421   164M / 11G   INFO    General                 (pair_info_count.cpp       : 322)   Min edge length for estimation: 87821
  0:00:28.421   164M / 11G   INFO    General                 (pair_info_count.cpp       : 333)   Estimating insert size for library #0
  0:00:28.422   164M / 11G   INFO    General                 (pair_info_count.cpp       : 190)   Estimating insert size (takes a while)
  0:00:28.522   436M / 11G   INFO    General                 (pair_info_count.cpp       :  39)   Selecting usual mapper
  0:00:30.638   430M / 11G   INFO    General                 (sequence_mapper_notifier.h:  94)   Total 2915760 reads processed
  0:00:31.856   430M / 11G   INFO    General                 (pair_info_count.cpp       : 208)   Edge pairs: 5059
  0:00:31.856   430M / 11G   INFO    General                 (pair_info_count.cpp       : 212)   0 paired reads (0% of all) aligned to long edges
  0:00:31.862   157M / 11G   WARN    General                 (pair_info_count.cpp       : 341)   Unable to estimate insert size for paired library #0
  0:00:31.863   157M / 11G   WARN    General                 (pair_info_count.cpp       : 347)   None of paired reads aligned properly. Please, check orientation of your read pairs.
0:00:31.863   157M / 11G   INFO   StageManager             (stage.cpp                 : 166)   STAGE == Distance Estimation
  0:00:31.863   157M / 11G   INFO    General                 (distance_estimation.cpp   : 183)   Clearing raw paired index
  0:00:31.864   157M / 11G   INFO   StageManager             (stage.cpp                 : 166)   STAGE == Repeat Resolving
  0:00:31.864   157M / 11G   WARN    General                 (repeat_resolving.cpp      :  63)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
```

3. Compared contigs.fasta sizes between 20-427FEC (5.5M) and 1-428RN3A (5.0M). Similar…

4. What is long edge? Longer than a length threshold LT (https://academic.oup.com/bioinformatics/article/35/14/i61/5529115)
Repeat resolution: help resolve repeats; there are algorithms for resolving repeats (https://academic.oup.com/bioinformatics/article/30/12/i293/384610)

5. Google error messages
  * https://github.com/ablab/spades/issues/315
  * https://github.com/ablab/spades/issues/403

6. Ask Jules what to do about corrupt reads

7. Compare 100X to 250X with `covstats.txt`

8. 1-428RN3A @ 100X (find "covstats):

```
vi stderr.5138118.ceres14-compute-65
	------------------   Results   ------------------
	Genome:                 1
	Key Length:             13
	Max Indel:              200
	Minimum Score Ratio:    0.56
	Mapping Mode:           normal
	Reads Used:             3728938 (549990749 bases)
	Mapping:                26.085 seconds.
	Reads/sec:              142952.33
	kBases/sec:             21084.41


	Pairing data:           pct pairs       num pairs       pct bases          num bases

	mated pairs:             98.1920%         1830760        98.1660%          539903804
	bad pairs:                1.2401%           23122         1.2611%            6935878
	insert size avg:          365.39


	Read 1 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.8881%         1862383        99.8873%          274684879
	unambiguous:             99.7431%         1859679        99.7554%          274322245
	ambiguous:                0.1450%            2704         0.1319%             362634
	low-Q discards:           0.0000%               0         0.0000%                  0

	perfect best site:       86.9527%         1621207        86.8843%          238927421
	semiperfect site:        87.1917%         1625662        87.1270%          239594706
	rescued:                  0.1736%            3236

	Match Rate:                   NA               NA        99.7535%          274024033
	Error Rate:              12.6750%          236131         0.2316%             636317
	Sub Rate:                12.6316%          235323         0.2243%             616157
	Del Rate:                 0.0560%            1043         0.0059%              16330
	Ins Rate:                 0.0661%            1231         0.0014%               3830
	N Rate:                   0.3154%            5876         0.0149%              40859


	Read 2 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.4581%         1854365        99.4534%          273492646
	unambiguous:             99.3135%         1851669        99.3218%          273130814
	ambiguous:                0.1446%            2696         0.1316%             361832
	low-Q discards:           0.0002%               3         0.0001%                180

	perfect best site:       83.3223%         1553518        83.1513%          228662684
	semiperfect site:        83.5518%         1557797        83.3844%          229303760
	rescued:                  0.1821%            3395

	Match Rate:                   NA               NA        99.5713%          272350996
	Error Rate:              16.0007%          296793         0.4130%            1129652
	Sub Rate:                15.9601%          296040         0.4002%            1094740
	Del Rate:                 0.0612%            1136         0.0113%              30854
	Ins Rate:                 0.0680%            1261         0.0015%               4058
	N Rate:                   0.2718%            5041         0.0157%              42852

	Reads:                                  3728938
	Mapped reads:                           3746597
	Mapped bases:                           556335989
	Ref scaffolds:                          266
	Ref bases:                              5101802

	Percent mapped:                         100.474
	Percent proper pairs:                   98.832
	Average coverage:                       109.047
	Average coverage with deletions:        108.297
	Standard deviation:                     48.568
	Percent scaffolds with any coverage:    97.74
	Percent of reference bases covered:     99.98

	Total time:             30.656 seconds.
```

9. 1-428RN3A @ 250X:

```
vi stderr.5138154.ceres14-compute-8
	------------------   Results   ------------------

	Genome:                 1
	Key Length:             13
	Max Indel:              200
	Minimum Score Ratio:    0.56
	Mapping Mode:           normal
	Reads Used:             7981132 (1177132132 bases)

	Mapping:                56.199 seconds.
	Reads/sec:              142014.70
	kBases/sec:             20945.66


	Pairing data:           pct pairs       num pairs       pct bases          num bases

	mated pairs:             98.1891%         3918302        98.1627%         1155504968
	bad pairs:                1.2451%           49687         1.2662%           14904394
	insert size avg:          366.77


	Read 1 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.8895%         3986156        99.8885%          587908797
	unambiguous:             99.7331%         3979915        99.7460%          587070056
	ambiguous:                0.1564%            6241         0.1425%             838741
	low-Q discards:           0.0000%               0         0.0000%                  0

	perfect best site:       86.9527%         3469903        86.8838%          511367722
	semiperfect site:        87.1876%         3479280        87.1225%          512772510
	rescued:                  0.1733%            6914

	Match Rate:                   NA               NA        99.7524%          586487523
	Error Rate:              12.6741%          505369         0.2319%            1363657
	Sub Rate:                12.6317%          503678         0.2247%            1321251
	Del Rate:                 0.0613%            2443         0.0059%              34693
	Ins Rate:                 0.0600%            2392         0.0013%               7713
	N Rate:                   0.3206%           12783         0.0157%              92310


	Read 2 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.4590%         3968978        99.4540%          585353736
	unambiguous:             99.3012%         3962678        99.3099%          584505702
	ambiguous:                0.1579%            6300         0.1441%             848034
	low-Q discards:           0.0001%               3         0.0000%                180

	perfect best site:       83.3115%         3324602        83.1393%          489330778
	semiperfect site:        83.5399%         3333714        83.3712%          490695668
	rescued:                  0.1801%            7186

	Match Rate:                   NA               NA        99.5700%          582899760
	Error Rate:              16.0080%          635528         0.4131%            2418120
	Sub Rate:                15.9673%          633914         0.4007%            2346032
	Del Rate:                 0.0683%            2710         0.0108%              63464
	Ins Rate:                 0.0625%            2482         0.0015%               8624
	N Rate:                   0.2841%           11279         0.0170%              99320

	Reads:                                  7981132
	Mapped reads:                           8019821
	Mapped bases:                           1190841534
	Ref scaffolds:                          276
	Ref bases:                              5103674

	Percent mapped:                         100.485
	Percent proper pairs:                   98.842
	Average coverage:                       233.330
	Average coverage with deletions:        231.721
	Standard deviation:                     102.626
	Percent scaffolds with any coverage:    98.19
	Percent of reference bases covered:     99.98

	Total time:             60.012 seconds.
```

10. Compare 1-428RN3A results side by side and found out that:
  * Same genome, key length, max indel, minimum score ratio, mapping mode stats. Reads used 100X = 3728938 (549990749 bases), 250X = 7981132 (1177132132 bases)
	* 250X takes twice as long to map. Reads/sec comparable. kBases/sec slightly higher for 100X.
	* Pairing data
    * Mated pairs: 100X slightly higher percentage for pairs and bases.
		*  Bad pairs: similar (250X slightly higher at 1.2451% vs 100X at 1.2401%)
		*  Insert size average: comparable
		* Read 1 data - comparable
		* Read 2 data - comparable
		* Ref scaffolds: slightly more at 250X (276) than 100X (266)
		* Ref bases: slightly higher at 250X (5103674) vs 100X (51001802)
		* Percent scaffolds with any coverage: slightly higher for 250X (98.19) vs 100X (97.74)
		* Total time: 30.656s (100X), 60.012s (250X)

11. Compare 94-439FED @100X:

```
vi stderr.5138140.ceres14-compute-8
	------------------   Results   ------------------

	Genome:                 1
	Key Length:             13
	Max Indel:              200
	Minimum Score Ratio:    0.56
	Mapping Mode:           normal
	Reads Used:             3688760 (549991013 bases)

	Mapping:                32.931 seconds.
	Reads/sec:              112015.43
	kBases/sec:             16701.41


	Pairing data:           pct pairs       num pairs       pct bases          num bases

	mated pairs:             97.6830%         1801646        97.6711%          537182330
	bad pairs:                1.7163%           31655         1.7260%            9492690
	insert size avg:          368.61


	Read 1 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.8565%         1841734        99.8561%          274598917
	unambiguous:             99.7539%         1839841        99.7606%          274336126
	ambiguous:                0.1026%            1893         0.0956%             262791
	low-Q discards:           0.0000%               0         0.0000%                  0

	perfect best site:       84.7542%         1563190        84.7211%          232978435
	semiperfect site:        85.0598%         1568826        85.0284%          233823495
	rescued:                  0.1793%            3307

	Match Rate:                   NA               NA        99.6902%          273766609
	Error Rate:              14.7904%          272534         0.2869%             787768
	Sub Rate:                14.7533%          271851         0.2794%             767220
	Del Rate:                 0.0629%            1159         0.0067%              18465
	Ins Rate:                 0.0366%             675         0.0008%               2083
	N Rate:                   0.3847%            7088         0.0229%              63005


	Read 2 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.4318%         1833901        99.4300%          273428808
	unambiguous:             99.3301%         1832025        99.3352%          273168125
	ambiguous:                0.1017%            1876         0.0948%             260683
	low-Q discards:           0.0000%               0         0.0000%                  0

	perfect best site:       80.1897%         1479003        80.1171%          220319224
	semiperfect site:        80.4937%         1484610        80.4228%          221159830
	rescued:                  0.1972%            3637

	Match Rate:                   NA               NA        99.5147%          272131088
	Error Rate:              19.0716%          349892         0.4610%            1260539
	Sub Rate:                19.0345%          349212         0.4491%            1228197
	Del Rate:                 0.0700%            1285         0.0108%              29441
	Ins Rate:                 0.0407%             746         0.0011%               2901
	N Rate:                   0.3541%            6496         0.0244%              66622

	Reads:                                  3688760
	Mapped reads:                           3718435
	Mapped bases:                           558122804
	Ref scaffolds:                          327
	Ref bases:                              5571353

	Percent mapped:                         100.804
	Percent proper pairs:                   98.648
	Average coverage:                       100.177
	Average coverage with deletions:        99.483
	Standard deviation:                     61.052
	Percent scaffolds with any coverage:    99.39
	Percent of reference bases covered:     99.99

	Total time:             36.755 seconds.
```

12. 94-439FED @250X:
```
vi stderr.5138158.ceres19-compute-55
	 ------------------   Results   ------------------
	Genome:                 1
	Key Length:             13
	Max Indel:              200
	Minimum Score Ratio:    0.56
	Mapping Mode:           normal
	Reads Used:             6848350 (1021069332 bases)

	Mapping:                40.818 seconds.
	Reads/sec:              167776.91
	kBases/sec:             25015.06

	Pairing data:           pct pairs       num pairs       pct bases          num bases

	mated pairs:             97.6805%         3344752        97.6682%          997260144
	bad pairs:                1.7242%           59041         1.7341%           17705866
	insert size avg:          369.70

	Read 1 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.8561%         3419248        99.8556%          509795498
	unambiguous:             99.7504%         3415627        99.7571%          509292481
	ambiguous:                0.1057%            3621         0.0985%             503017
	low-Q discards:           0.0000%               0         0.0000%                  0

	perfect best site:       84.7664%         2902549        84.7335%          432592159
	semiperfect site:        85.0674%         2912857        85.0362%          434137555
	rescued:                  0.1776%            6081

	Match Rate:                   NA               NA        99.6899%          508254618
	Error Rate:              14.7827%          505695         0.2874%            1465350
	Sub Rate:                14.7458%          504430         0.2788%            1421306
	Del Rate:                 0.0638%            2181         0.0079%              40076
	Ins Rate:                 0.0372%            1271         0.0008%               3968
	N Rate:                   0.3824%           13080         0.0227%             115606

	Read 2 data:            pct reads       num reads       pct bases          num bases

	mapped:                  99.4366%         3404883        99.4345%          507649541
	unambiguous:             99.3310%         3401266        99.3359%          507146405
	ambiguous:                0.1056%            3617         0.0986%             503136
	low-Q discards:           0.0000%               0         0.0000%                  0

	perfect best site:       80.1783%         2745445        80.1056%          408968444
	semiperfect site:        80.4804%         2755789        80.4093%          410519187
	rescued:                  0.1972%            6751

	Match Rate:                   NA               NA        99.5143%          505241380
	Error Rate:              19.0836%          650029         0.4611%            2341060
	Sub Rate:                19.0470%          648784         0.4487%            2278158
	Del Rate:                 0.0700%            2383         0.0114%              57731
	Ins Rate:                 0.0406%            1383         0.0010%               5171
	N Rate:                   0.3562%           12134         0.0246%             124832

	Reads:                                  6848350
	Mapped reads:                           6905093
	Mapped bases:                           1036410249
	Ref scaffolds:                          334
	Ref bases:                              5572991

	Percent mapped:                         100.829
	Percent proper pairs:                   98.666
	Average coverage:                       185.970
	Average coverage with deletions:        184.683
	Standard deviation:                     112.917
	Percent scaffolds with any coverage:    99.40
	Percent of reference bases covered:     99.99

	Total time:             45.197 seconds.
```

13. Compare 94-439FED results side by side and found out that:
  * Same genome, key length, max indel, minimum score ratio, mapping mode stats. Reads used 100X = 3688760 (549991013 bases), 250X = 6848350 (1021069332 bases)
  * 250X took 8 seconds longer to map. Reads/sec and kBases/sec are higher for 250X.
  * Pairing data
	 * Mated pairs: comparable between the two
   * Bad pairs: similar (250X slightly higher at 1.7242% vs 100X at 1.7163%)
   * Insert size average: comparable
	* Read 1 data: comparable
	* Read 2 data: comparable
	* Ref scaffolds: slightly more at 250X (334) than 100X (327)
	* Ref bases: slightly higher at 250X (5572991) vs 100X (5571353)
	* Percent scaffolds with any coverage: comparable
	* Total time: 36.755s (100X), 45.197s (250X)

14. Based on conclusions for comparing coverage results for 1-428RN3A and 94-439FED, I will run at 250X coverage for all the samples and for ones that didn't work, troubleshoot those and ask Jules

15. Looked up quast: https://github.com/ablab/quast
  * Would I need to run this tool or is the template script Jules sent me "evaluation" section good enough for genome assembly evaluation?

16. Emailed Jules a bunch of questions on:
  * Not getting paired reads because they may be corrupt- what to do
  * How to evaluate genome assemblies (bbamp.sh, other scripts in the slurm template)

**21Oct2020**
Jules response and my attempts to troubleshoot with 100X directory

```
1. Corrupted PE reads
I would go back to the original raw data and do some investigation
2. Whats the file sizes for the R1 and the R2 files? Are they similar
-rw-rw-r--.  1 kathy.mou       proj-fsepru 311M Oct  7 15:02 1-428RN3A_1.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru 317M Oct  7 15:02 1-428RN3A_2.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru  386 Oct  7 15:02 96-441FEC_1.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru  377 Oct  7 15:02 96-441FEC_2.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru 267M Oct  7 15:02 94-439FED_1.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru 273M Oct  7 15:02 94-439FED_2.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru 232M Oct  7 15:01 20-427FEC_2.fastq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru 232M Oct  7 15:01 20-427FEC_1.fastq.gz
Sizes of R1 and R2 are slightly different for all samples except for 20 (exact same sizes)
Run “repair.sh” from the bbtools suite that can help diagnose issues.
bbsplitpairs.sh in=20-427FEC.fq.gz out=fixed.20-427FEC.fq.gz outs=singletons.20-427FEC.fq.gz minlen=70 fint
Submitted batch job 5176991
If that didn't work, try
repair.sh in=20-427FEC.fq.gz out=fixed.20-427FEC.fq.gz outs=singletons.20-427FEC.fq.gz repair
The scaffolds are just contigs that have been stitched together using paired end data. That is they know how far R1 and R2 should be so if there are read pairs that span the gap between contigs they will scaffold this gap together using the PE info.
Because your reads aren’t paired they cant do this and only the contigs file is generated
2. Assessing assembly quality
	a. The main things to look at are how many contigs you have, The N50 (which will give you an idea of how much of your assembly is in long contigs), and the total length of the assembly. These things you determine without calling bbmap. Check out “stats.sh” from the bbtools suite. (Or quast etc)
	b. An important thing to look at from the “Evaluation” section is the file output by the bbmap call, the “covstats” file.  Bbmap is a read mapping tool,  it maps your short reads back to a “reference”, in this case the reference is your draft assembly. The covstats file is the statistics for how your reads map back to your assembly.  Using this file you can tell if there are certain contigs that have low coverage, or are only partially covered.  If there are some contigs with strange coverage, you can omit them.   
3. Those extra scripts.
	a. The Rscript line calls an Rscript I wrote and the covstats file output by bbmap to remove contigs from the assembly that have low coverage relative to the length weighted average coverage of the longest 10 contigs.  Basically if the chromosome is at 100x, contigs need to have at least 30x coverage to be kept.  I can provide you with this script if you are interested.
	b. The Rscript just outputs a file of contig names to keep.  This is used by “filterbyname.sh” to keep only the good contigs
	c. The adapt_polish.sh line provides a function that is a wrapper for pilon, where you can specify the number of pilon iterations you want to run, and it will terminate early if there are no further changes to be made.  I can provide this to you as well if you’d like.  
		a. Pilon: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0112963
 
How about you try these tips out, don’t worry about the extra scripts for now.  Get to the point where you make a decision about how much data to include using a few genomes (100x vs 250x) and then we can. Meet up and discuss things and I can show you how to make use of the extra scripts.
```

**22Oct2020**
1. Looked at stdout for job 5176991 and saw no output. Looked at stderr and got this
```
java -ea -Xmx200m -cp /home/kathy.mou/software/bbmap/current/ jgi.SplitPairsAndSingles in=20-427FEC.fq.gz out=fixed.20-427FEC.fq.gz outs=singletons.20-427FEC.fq.gz minlen=70 fint
Executing jgi.SplitPairsAndSingles [in=20-427FEC.fq.gz, out=fixed.20-427FEC.fq.gz, outs=singletons.20-427FEC.fq.gz, minlen=70, fint]
Paired input disabled; running in FixInterleaving mode
Started output stream.
Input:                          5957666 reads           893649900 bases.
Result:                         5957666 reads (100.00%)         893649900 bases (100.00%)
Pairs:                          0 reads (0.00%)         0 bases (0.00%)
Singletons:                     5957666 reads (100.00%)         893649900 bases (100.00%)
Time:                           9.201 seconds.
Reads Processed:       5957k    647.50k reads/sec
Bases Processed:        893m    97.13m bases/sec
```

2. The only output I have is `singletons.20-427FEC.fq.gz`. Try running with `repair.sh` instead
```
repair.sh in=20-427FEC.fq.gz out=nonint_fixed.20-427FEC.fq.gz outs=nonint_singletons.20-427FEC.fq.gz repair
Submitted batch job 5181381bbmap.sh
```

3. This seemed to produce both outputs (`nonint_fixed.20-427FEC.fq.gz`, `nonint_singletons.20-427FEC.fq.gz`)

4. Looked at `nonint_fixed.20-427FEC.fq.gz` (garbled letters and symbols like `20-427FEC.fq.gz`), `nonint_singletons.20-427FEC.fq.gz`

5. Renamed `20-427FEC.fq.gz` to `old_20-427FEC.fq.gz`

6. Edited `20-427FEC.slurm` using the new fixed 20-427FEC file (renamed `nonint_fixed.20-427FEC.fq.gz` to `20-427FEC_temp.fq.gz`). Run `bbduk.sh` down to `bbmap.sh`

7. I edited reformat.sh to include corrected sbt parameter: 6000000 * 100 = 600000000
```
reformat.sh in=20-427FEC_temp.fq.gz sbt=600000000 out=20-427FEC_subsamp.fq.gz interleaved=t
rm 20-427FEC_temp.fq.gz; ln -s 20-427FEC_subsamp.fq.gz 20-427FEC_temp.fq.gz
```

8. Make new folder `20-427FEC_old` to move all the files generated from `20-427FEC.slurm` pipeline the first time around.

9. Run `20-427FEC.slurm` on bash
```
Submitted batch job 5181664
```

10. Job failed because I still had `20-427FEC_spades_out/` folder. Moved that to `20-427FEC_old` directory.

11. Modify `20-427FEC.slurm` to run `repair.sh`, output `20-427FEC_temp.fq.gz` and run all the way to `bbmap.sh`. Run slurm job.
```
Submitted batch job 5181817
```

12. Job failed.
```
vi stdout.5181817.ceres18-compute-24.
```

13. Got same error correction and assembling warning messages (couldn't estimate insert size for paired library. None of paired reads aligned properly. Insert size not estimated for any of the paired libraries). I will try running `bbsplitpairs.sh` again

```
bbsplitpairs.sh in=old_20-427FEC.fq.gz out=fixedint.20-427FEC.fq.gz outs=intsingletons.20-427FEC.fq.gz minlen=70 fint
Submitted batch job 5182016
ls -alth
-rw-rw-r--.  1 kathy.mou       proj-fsepru   20 Oct 22 17:01 fixedint.20-427FEC.fq.gz
-rw-rw-r--.  1 kathy.mou       proj-fsepru 253M Oct 22 17:01 intsingletons.20-427FEC.fq.gz
```

14. 'Good' pairs (small file size) - `fixedint.20-427FEC.fq.gz`
	'Good' singletons (same file size as `old_20-427FEC.fq.gz`)- `intsingletons.20-427FEC.fq.gz`

15. Need to ask Jules what I did wrong with `repair.sh` or `bbsplitpair.sh`

16. In the mean time, try `stats.sh`
  * The main things to look at are how many contigs you have, The N50 (which will give you an idea of how much of your assembly is in long contigs), and the total length of the assembly. These things you determine without calling bbmap. Check out `stats.sh` from the bbtools suite. (Or quast etc)

17. 1-428RN3A 100X: Created `evalassembly1-428RN3A.slurm` with this command to run:
```
stats.sh in=scaffolds.fasta
```

18. Moved slurm template to `1-428RN3A_spades_out` directory and ran batch job on slurm
```
Submitted batch job 5182243
```

19. Examined output

```
vi stdout.5182243.ceres18-compute-24
		A       C       G       T       N       IUPAC   Other   GC      GC_stdev
		0.2472  0.2548  0.2515  0.2465  0.0001  0.0000  0.0000  0.5063  0.0811

		Main genome scaffold total:             266
		Main genome contig total:               270
		Main genome scaffold sequence total:    5.102 MB
		Main genome contig sequence total:      5.101 MB        0.006% gap
		Main genome scaffold N/L50:             9/144.648 KB
		Main genome contig N/L50:               10/144.648 KB
		Main genome scaffold N/L90:             37/24.344 KB
		Main genome contig N/L90:               39/24.315 KB
		Max scaffold length:                    624.593 KB
		Max contig length:                      624.593 KB
		Number of scaffolds > 50 KB:            26
		% main genome in scaffolds > 50 KB:     82.63%

		Minimum         Number          Number          Total           Total           Scaffold
		Scaffold        of              of              Scaffold        Contig          Contig
		Length          Scaffolds       Contigs         Length          Length          Coverage
		--------        --------------  --------------  --------------  --------------  --------
		    All                    266             270       5,101,802       5,101,492    99.99%
		    100                    266             270       5,101,802       5,101,492    99.99%
		    250                    214             218       5,091,714       5,091,404    99.99%
		    500                    164             168       5,074,158       5,073,848    99.99%
		   1 KB                    123             127       5,046,376       5,046,066    99.99%
		 2.5 KB                     80              84       4,974,532       4,974,222    99.99%
		   5 KB                     65              69       4,921,799       4,921,489    99.99%
		  10 KB                     53              57       4,838,375       4,838,065    99.99%
		  25 KB                     36              40       4,587,442       4,587,132    99.99%
		  50 KB                     26              28       4,215,872       4,215,762   100.00%
		 100 KB                     16              17       3,425,646       3,425,546   100.00%
		 250 KB                      5               6       1,902,738       1,902,638    99.99%
		 500 KB                      1               1         624,593         624,593   100.00%
```

20. What is N50
  * http://www.metagenomics.wiki/pdf/definition/assembly/n50
```
N50 is the shortest contig length needed to cover 50% of the genome. 
Meaning
			-> Half of the genome sequence is covered by contigs larger than or equal the N50 contig size. 
			-> The sum of the lengths of all contigs of size N50 or longer contain at least 50 percent of the total genome sequence.
```

21. What is L50
  * https://www.molecularecologist.com/2017/03/whats-n50/

22. 1-428RN3A 250X: Copied `evalassembly1-428RN3A.slurm` to `1-428RN3A_spades_out/` directory in `FS19C_4Samples250X` directory

23. Ran batch job
```
Submitted batch job 5185424
```

24. Examined output
```
vi stdout.5185424.ceres14-compute-54
		A       C       G       T       N       IUPAC   Other   GC      GC_stdev
		0.2463  0.2537  0.2525  0.2475  0.0001  0.0000  0.0000  0.5062  0.0809

		Main genome scaffold total:             276
		Main genome contig total:               280
		Main genome scaffold sequence total:    5.104 MB
		Main genome contig sequence total:      5.103 MB        0.008% gap
		Main genome scaffold N/L50:             10/133.539 KB
		Main genome contig N/L50:               11/133.539 KB
		Main genome scaffold N/L90:             39/24.344 KB
		Main genome contig N/L90:               41/23.661 KB
		Max scaffold length:                    624.593 KB
		Max contig length:                      624.593 KB
		Number of scaffolds > 50 KB:            26
		% main genome in scaffolds > 50 KB:     80.67%

		Minimum         Number          Number          Total           Total           Scaffold
		Scaffold        of              of              Scaffold        Contig          Contig
		Length          Scaffolds       Contigs         Length          Length          Coverage
		--------        --------------  --------------  --------------  --------------  --------
		    All                    276             280       5,103,674       5,103,274    99.99%
		    100                    276             280       5,103,674       5,103,274    99.99%
		    250                    218             222       5,092,326       5,091,926    99.99%
		    500                    166             170       5,074,388       5,073,988    99.99%
		   1 KB                    126             130       5,047,527       5,047,127    99.99%
		 2.5 KB                     81              85       4,972,451       4,972,051    99.99%
		   5 KB                     67              71       4,922,488       4,922,088    99.99%
		  10 KB                     56              60       4,845,791       4,845,391    99.99%
		  25 KB                     38              42       4,578,560       4,578,160    99.99%
		  50 KB                     26              27       4,117,150       4,117,050   100.00%
		 100 KB                     16              17       3,320,408       3,320,308   100.00%
		 250 KB                      4               5       1,646,796       1,646,696    99.99%
		 500 KB                      1               1         624,593         624,593   100.00%
```

25. 94-439FED 100X: Created `evalassembly94-439FED.slurm` with this command to run:
```
stats.sh in=scaffolds.fasta
```

26. Moved slurm template to `94-439FED_spades_out` directory and ran batch job on slurm
```
Submitted batch job 5185420
```

27. Examined output
```
vi stdout.5185420.ceres14-compute-54
		A       C       G       T       N       IUPAC   Other   GC      GC_stdev
		0.2469  0.2531  0.2538  0.2462  0.0001  0.0000  0.0000  0.5069  0.0687

		Main genome scaffold total:             327
		Main genome contig total:               335
		Main genome scaffold sequence total:    5.571 MB
		Main genome contig sequence total:      5.571 MB        0.014% gap
		Main genome scaffold N/L50:             16/93.761 KB
		Main genome contig N/L50:               16/93.761 KB
		Main genome scaffold N/L90:             62/16.697 KB
		Main genome contig N/L90:               64/13.684 KB
		Max scaffold length:                    378.381 KB
		Max contig length:                      378.381 KB
		Number of scaffolds > 50 KB:            34
		% main genome in scaffolds > 50 KB:     74.03%

		Minimum         Number          Number          Total           Total           Scaffold
		Scaffold        of              of              Scaffold        Contig          Contig
		Length          Scaffolds       Contigs         Length          Length          Coverage
		--------        --------------  --------------  --------------  --------------  --------
		    All                    327             335       5,571,353       5,570,553    99.99%
		    100                    327             335       5,571,353       5,570,553    99.99%
		    250                    300             308       5,566,114       5,565,314    99.99%
		    500                    260             268       5,551,064       5,550,264    99.99%
		   1 KB                    188             196       5,500,541       5,499,741    99.99%
		 2.5 KB                    117             125       5,377,812       5,377,012    99.99%
		   5 KB                     92              96       5,284,197       5,283,797    99.99%
		  10 KB                     73              77       5,151,930       5,151,530    99.99%
		  25 KB                     52              53       4,807,516       4,807,416   100.00%
		  50 KB                     34              35       4,124,649       4,124,549   100.00%
		 100 KB                     14              15       2,616,370       2,616,270   100.00%
		 250 KB                      1               1         378,381         378,381   100.00%
```

28. 94-439FED 250X: Copied `evalassembly94-439FED.slurm` to `FS19C_4Samples250X/94-439FED_spades_out/` directory

29. Ran batch job
```
Submitted batch job 5185753
```

30. Examined output
```
vi stdout.5185753.ceres19-compute-55
		A       C       G       T       N       IUPAC   Other   GC      GC_stdev
		0.2471  0.2541  0.2528  0.2460  0.0001  0.0000  0.0000  0.5069  0.0702

		Main genome scaffold total:             334
		Main genome contig total:               342
		Main genome scaffold sequence total:    5.573 MB
		Main genome contig sequence total:      5.572 MB        0.014% gap
		Main genome scaffold N/L50:             16/93.761 KB
		Main genome contig N/L50:               16/93.761 KB
		Main genome scaffold N/L90:             62/15.462 KB
		Main genome contig N/L90:               64/13.684 KB
		Max scaffold length:                    378.381 KB
		Max contig length:                      378.381 KB
		Number of scaffolds > 50 KB:            34
		% main genome in scaffolds > 50 KB:     74.20%

		Minimum         Number          Number          Total           Total           Scaffold
		Scaffold        of              of              Scaffold        Contig          Contig
		Length          Scaffolds       Contigs         Length          Length          Coverage
		--------        --------------  --------------  --------------  --------------  --------
		    All                    334             342       5,572,991       5,572,191    99.99%
		    100                    334             342       5,572,991       5,572,191    99.99%
		    250                    303             311       5,566,940       5,566,140    99.99%
		    500                    264             272       5,552,181       5,551,381    99.99%
		   1 KB                    187             195       5,497,903       5,497,103    99.99%
		 2.5 KB                    116             124       5,375,773       5,374,973    99.99%
		   5 KB                     93              97       5,290,519       5,290,119    99.99%
		  10 KB                     72              76       5,142,966       5,142,566    99.99%
		  25 KB                     52              53       4,818,161       4,818,061   100.00%
		  50 KB                     34              35       4,135,277       4,135,177   100.00%
		 100 KB                     14              15       2,626,942       2,626,842   100.00%
		 250 KB                      1               1         378,381         378,381   100.00%
```

31. Need to ask Jules what threshold is for determining "there is a large discrepancy and thus will use this coverage" vs either one is fine
```
An important thing to look at from the “Evaluation” section is the file output by the bbmap call, the “covstats” file.  Bbmap is a read mapping tool,  it maps your short reads back to a “reference”, in this case the reference is your draft assembly. The covstats file is the statistics for how your reads map back to your assembly.  Using this file you can tell if there are certain contigs that have low coverage, or are only partially covered.  If there are some contigs with strange coverage, you can omit them.
```

32. Download `covstats.txt` files to local and open in text editor
  * 1-428RN3A
    * 100X
		* 250X
	* 94-439FED
		* 100X
		* 250X

33. Try running the three strains (2 that worked, 1 that is not supposed to work) at 6Mb coverage to see how that affects stats for genome assembly evaluation

34. Make new directory 5Mbgenome and transfer all files to that.

35. Modify slurm template for all four files (250X) at the `reformat.sh` step to change 100X to 250X. Make new copy of slurm, with `SampleID.6Mb.slurm`

36. Copy/move slurm template files and `SampleID_temp.fq.gz` to `FS19C_4Samples250X` directory or `FS19C_4Samples100X` directory

37. Edit slurm template files:
  * Comment out 1st `reformat.sh` command
	* `reformat.sh`:
	 * 100X: sbt=600000000
	 * 250X: sbt 600000000*250=1500000000
	* `spades.py`:
	 * 100X: `-o SampleID_6Mb_100X_spades_out`
	 * 250X: `-o SampleID_6Mb_250X_spades_out`
	* `bbmap.sh`:
	 * 100X: `ref=SampleID_6Mb_100X_spades_out/scaffolds.fasta, covstats=SampleID_6Mb_100X_covstats.txt`
	 * 250X: `ref=SampleID_6Mb_250X_spades_out/scaffolds.fasta, covstats=SampleID_6Mb_250X_covstats.txt`

38. Completed slurm template file editing for 100X, 250X of the following, with notes:
	* 1-428
	* 20-427FEC
	 * 100X (I used the repaired `temp.fq.gz` file, not the original)
	 * 250X
	  * Try `repair.sh` with rest of script
		* `repair.sh in=20-427FEC.fq.gz out=20-427FEC_temp.fq.gz outs=nonint_singletons.20-427FEC.fq.gz repair`
	* 94-439FED
  * 96-441FEC

**27Oct2020** Run jobs on slurm
1. 1-428RN3A 100X
```
Submitted batch job 5203616
```
  * Looked at stderr, assembly was successful

2. 1-428RN3A250X
```
Submitted batch job 5203622
```
  * Looked at stderr, assembly was successful

3. 20-427FEC 100X: I used the repaired `temp.fq.gz` file, not the original
```
Submitted batch job 5203618
```
  * Got the same warnings as before even with new repair.sh
```
    	=== Error correction and assembling warnings:
			 * 0:00:33.273   159M / 11G   WARN    General                 (pair_info_count.cpp       : 341)   Unable to estimate insert size for paired library #0
			 * 0:00:33.274   159M / 11G   WARN    General                 (pair_info_count.cpp       : 347)   None of paired reads aligned properly. Please, check orientation of your read pairs.
			 * 0:00:33.275   159M / 11G   WARN    General                 (repeat_resolving.cpp      :  63)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
```			

4. 20-427FEC 250X
```
Submitted batch job 5203623
```
  * Got the same warnings as before even with new `repair.sh`
```
		=== Error correction and assembling warnings:
		 * 0:00:32.863   162M / 11G   WARN    General                 (pair_info_count.cpp       : 341)   Unable to estimate insert size for paired library #0
		 * 0:00:32.864   162M / 11G   WARN    General                 (pair_info_count.cpp       : 347)   None of paired reads aligned properly. Please, check orientation of your read pairs.
		 * 0:00:32.865   162M / 11G   WARN    General                 (repeat_resolving.cpp      :  63)   Insert size was not estimated for any of the paired libraries, repeat resolution module will not run.
```

5. 94-439FED 100X
		a. 100X
			i. Submitted batch job 5203619
			ii. Looked at stderr, assembly was successful
		b. 250X
			i. Submitted batch job 5203626
			ii. Looked at stderr, assembly was successful
	4. 96-441FEC
		a. 100X
			i. Submitted batch job 5203621
			ii. Looked at stderr, assembly unsuccessful
		b. 250X
			i. Submitted batch job 5203628
			ii. Looked at stderr, assembly unsuccessful
2. Run stats.sh for samples 1 and 94
	1. Make copy of slurm files, rename as SampleID.6Mb.stats.slurm, comment all scripts, add stats.sh in=scaffolds.fasta, run on slurm
		a. 1-428RN3A
			i. 100X
				1) Submitted batch job 5204491
			ii. 250X
				1) Submitted batch job 5204494
		b. 94-439FED
			i. 100X
				1) Submitted batch job 5204492
			ii. 250X
				1) Submitted batch job 5204493
