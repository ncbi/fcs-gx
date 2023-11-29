# FCS-GX
This file outlines how to run FCS-GX from the source code. 

### Building requirements

- GCC >= 7.3
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
    ```curl -LO https://github.com/ncbi/fcs/raw/main/examples/fcsgx_test.fa.gz```  

- Make a directory for output files.  
    ```mkdir gx_out```

- Run GX.  
    ```./dist/run_gx --fasta=fcsgx_test.fa.gz --tax-id=6973 --gx-db=/dev/shm/gxdb/test-only --out-dir=./gx_out/```  

- A successful run should produce these two files inside of gx_out directory:
    ```
    fcsgx_test.fa.6973.fcs_gx_report.txt
    fcsgx_test.fa.6973.taxonomy.rpt
    ```

### For normal runs, use the complete *all* database:
- Download the database to the local disk, and then copy it to a RAM-backed tmpfs or ramfs location (we will use `/dev/shm/gxdb` in examples below).
   ```
   ./scripts/sync_files.py sync-in-place --mft=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest --dir=/path/to/local/disk/gxdb
   ./scripts/sync_files.py sync-in-place --mft=/path/to/local/disk/gxdb/all.manifest --dir=/dev/shm/gxdb
   ```

- Alternatively, you may skip the disk and download to tmpfs directly.
   ```
   ./scripts/sync_files.py  sync-in-place --mft=https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/latest/all.manifest  --dir=/dev/shm/gxdb
   ```

- Alternatively, you can `vmtouch` downloaded database on local disk and use it directly instead of tmpfs.
   ```
   vmtouch -m1000G -v -t /path/to/local/disk/gxdb
   vmtouch -m1000G -v    /path/to/local/disk/gxdb  # execute again to verify that the database files are in RAM and have not been swapped-out.
   ```

- Run GX.  
``` ./dist/run_gx --fasta fcsgx_test.fa.gz --tax-id 6973 --gx-db /dev/shm/gxdb/all --out-dir ./gx_out/```  


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
### Output files
Please see the [GitHub documentation](https://github.com/ncbi/fcs/wiki/FCS-GX) for more details on the output files.

## Releases

Please see [release history](https://github.com/ncbi/fcs/releases) 

## Contact

Please create an [issue](https://github.com/ncbi/fcs-gx/issues) if you encounter any problems.
For all other questions or comments, please contact us at refseq-support@nlm.nih.gov
