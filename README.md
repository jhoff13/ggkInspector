# ggkInspector 

ggkInspector is a python program to be run on Biotite for catching errors in ggKbase import files in a user friendly manner. Once all of the flagged errors are fixed the file will be prepared for import.
</br></br>
Follow the ggkhelp page for more info: https://ggkbase-help.berkeley.edu/overview/data-import/. </br>
Import tsv template can be downloaded here: https://ggkbase-help.berkeley.edu/wp-content/uploads/2015/11/ggkbase_bulk_upload_template.tsv
## Running the cmd:
```
python /home/jhoff/scripts/ggkInspector.py -t /path/to/ggkimport.tsv
```
## Tips for a succesful import:
ggKbase has a particular formatting scheme that trips users up in the same places. Here are some things to consider as you start tbe metagenomic's pipeline.</br>
- **Slugs**: Slugs are the project names that appear in the URL and therefore cannot have any disruptive charecters (ie: whitespaces, '.', '/', '*').
- The Slug must be present in the assembly fasta headers and the subsequent downstream files. It is therefore easiest to make sure this is consistent before starting the pipeline before generating the downstream files.
- Changing fasta headers using *sed*. If you need to make changes to fasta headers here is a cmd you can run to do that. However you should be very careful as this will alter the file. Echo the cmds and test run on a copy fasta first.
```
for i in path/to/dir/*.fasta; do echo sed -i 's/current header/fixed header/g' $i; done | bash
```
- You can validate that that the target string is present with *grep*. 
```
for i in path/to/dir/*.fasta; do grep "target string" $i; done | wc -l
```
- Read length and coverage must be present in the assembly fasta headers which need to be artifically added to certain data (curated genomes, NCBI genomes, etc.). You can add this to the end of the header: read_length_0 read_count_1.

## User Feedback:
Please report errors using this google form: https://forms.gle/dFPX7L9C2fryY9LS6 </br>
