# FCS-GX
This file outlines how to run FCS-GX from the source code. 

### Building requirements

- GCC >= 12.2.0
- cmake >= 3.1
- python >= 3.8

### Execution requirements

- A host with sufficient RAM to hold the database and accessory files (approximately 470 GiB). A host with 32-64 CPUs and 512 GiB of RAM is sufficient for execution. Trying to run on a server without sufficient memory will result in extremely long run times (as much as a 10000x difference in performance). Optionally, sufficient disk space to save a local copy of the database files to avoid subsequent downloads from NCBI's FTP site.  
- A genome assembly in FASTA format.
- The tax-id of the organism.
- For downloading the database, `rclone` or `aria2c` installed and available from $PATH.

### Build
- From the repository root, run `make`.
- After a successful build, run the following for the help menu. `./dist/run_gx --help`

   
### Verify functionality by using a small *test-only* database.
- Create a folder in your shared memory space. This is where you will save the GX database.  
    ```mkdir /dev/shm/gxdb```
    
- Download the test-only GX database to your shared memory for testing purposes.  
    ```./scripts/sync_files.py get --mft=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest  --dir=/dev/shm/gxdb```

- Retrieve this test fasta file.  
    ```curl -LO https://zenodo.org/records/10932013/files/FCS_combo_test.fa```  

- Make a directory for output files.  
    ```mkdir gx_out```

- Run GX.  
    ```./dist/run_gx --fasta=FCS_combo_test.fa --tax-id=4932 --gx-db=/dev/shm/gxdb/test-only --out-dir=./gx_out/```  

- A successful run should produce these two files inside of gx_out directory:
    ```
    FCS_combo_test.4932.fcs_gx_report.txt
    FCS_combo_test.4932.taxonomy.rpt
    ```

### For normal runs, use the complete *all* database:
- Download the database to the local disk, and then copy it to a RAM-backed tmpfs or ramfs location (we will use `/dev/shm/gxdb` in examples below).
   ```
   ./scripts/sync_files.py get --mft=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir=/path/to/local/disk/gxdb
   ./scripts/sync_files.py get --mft=/path/to/local/disk/gxdb/all.manifest --dir=/dev/shm/gxdb
   ```

- Alternatively, you may skip the disk and download to tmpfs directly.
   ```
   ./scripts/sync_files.py get --mft=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest  --dir=/dev/shm/gxdb
   ```

- Alternatively, you can `vmtouch` downloaded database on local disk and use it directly instead of tmpfs.
   ```
   vmtouch -m1000G -v -t /path/to/local/disk/gxdb
   vmtouch -m1000G -v    /path/to/local/disk/gxdb  # execute again to verify that the database files are in RAM and have not been swapped-out.
   ```

- Run GX.  
``` ./dist/run_gx --fasta FCS_combo_test.fa --tax-id 4932 --gx-db /dev/shm/gxdb/all --out-dir ./gx_out/```  


### Environment Variables 

You can run FCS-GX with environment variables to control the number of CPUs used (GX_NUM_CORES) or to exclude alignments to particular tax-ids (GX_ALIGN_EXCLUDE_TAXA). Multiple tax-ids may be provided as a comma-separated list. Note this only works for bottom-level tax-ids explicitly in the database, e.g. setting GX_ALIGN_EXCLUDE_TAXA=33208 will not exclude all metazoan hits. 

For example, to run a genome with 8 CPUs and excluding alignments to Toxoplasma gondii, set the variables as follows prior to running GX:   

```
export GX_NUM_CORES=8
export GX_ALIGN_EXCLUDE_TAXA=5811
```

### Useful GX subcommands
The sequences used to build the gx database are listed in the file *all.seq_info.tsv.gz* within the gxdb folder. From there, you can select the sequences of your choice, and then generate the fasta files using the *gx get-fasta* subcommand:   
 
```
./dist/gx get-fasta --db /dev/shm/gxdb/all.gxi --input 3col.txt  --output out.txt
```

The input file, which is provided by the user, is a tab delimited, 3 column file in the following format, along with the header:
```
cat 3col.txt 
##[["GX locs",1,1]]
NC_060925.1     .       .
```

To get the fasta for a specific set of coordinates, format your input file with the start and end coordinates in the 2nd and 3rd column, respectively:
```
##[["GX locs",1,1]]
NC_060925.1     1       200
```
### FCS-GX wiki
Please see the [FCS-GX wiki](https://github.com/ncbi/fcs/wiki/FCS-GX-quickstart) for more details on input requirements, outputs, and troubleshooting.

## Releases

Please see [release history](https://github.com/ncbi/fcs/releases) 

## Citing FCS-GX
Astashyn A, Tvedte ES, Sweeney D, Sapojnikov V, Bouk N, Joukov V, Mozes E, Strope PK, Sylla PM, Wagner L, Bidwell SL, Brown LC, Clark K, Davis EW, Smith-White B, Hlavina W, Pruitt KD, Schneider VA, Murphy TD. Rapid and sensitive detection of genome contamination at scale with FCS-GX. Genome Biol. 2024 Feb 26;25(1):60. doi: 10.1186/s13059-024-03198-7. PMID: 38409096; PMCID: PMC10898089.

[Read the FCS-GX paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-024-03198-7)

## Licensing
The NCBI FCS tool suite software authored by NCBI is a "United States
Government Work" under the terms of the United States Copyright
Act. It was written as part of the authors' official duties as United
States Government employees and thus cannot be copyrighted. This
software is freely available to the public for use. The National
Library of Medicine and the U.S. Government have not placed any
restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy
and reliability of the software and data, the NLM and the
U.S. Government do not and cannot warrant the performance or results
that may be obtained by using this software or data. The NLM and the
U.S. Government disclaim all warranties, express or implied, including
warranties of performance, merchantability or fitness for any
particular purpose.

Please cite NCBI in any work or product based on this material.

## Funding
This work was supported by the National Center for Biotechnology Information of the National Library of Medicine (NLM), National Institutes of Health.

FCS is part of the [NIH Comparative Genomics Resource (CGR)](https://www.ncbi.nlm.nih.gov/comparative-genomics-resource/), an NLM project to establish an ecosystem to facilitate reliable comparative genomics analyses for all eukaryotic organisms.

## Contact

Please create an [issue](https://github.com/ncbi/fcs-gx/issues) if you encounter any problems.
For all other questions or comments, please contact us at refseq-support@nlm.nih.gov
