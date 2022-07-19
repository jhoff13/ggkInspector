#Import Libraries
import pandas as pd
pd.set_option('display.max_colwidth', None)
from Bio import SeqIO
import os, sys, glob, argparse, re

#Import Variables:

parser = argparse.ArgumentParser(description='Check the viabilitiy of the tsv file for ggkBase import.\n\
                                                Follow the formatting from this sheet: https://ggkbase-help.berkeley.edu/wp-content/uploads/2015/11/ggkbase_bulk_upload_template.tsv\n\
                                                Pass in a path to the .tsv file you wish to import -t')
parser.add_argument('-t','--tsv',type=str,help='File path to tsv file for ggkBase import')
parser.add_argument('-c','--comprehensive',type=bool,default=False,help='Boolean value to use heuristics to save time [default==False]')
parser.add_argument('-p','--path_check',type=bool,default=False,help='Boolean value to skip to checking read and assembly paths [default==False]')
parser.add_argument('-v','--verbose',type=bool,default=True,help='Boolean value to display messages [default==True]')
parser.add_argument('-l','--log',type=bool,default=False,help='Boolean value to write import_logger.log [default==False]')

args = parser.parse_args()
input_DF = pd.read_csv(args.tsv,'\t')
    

#logging functions and DF-----------------------------------

first_run = True
CWD = os.getcwd()
def log_writter(string,log_BOOL=args.log):
    """Initializes 'import_logger.log' file in working directory and appends."""
    if log_BOOL:
        global first_run
        #initialize log file:
        if first_run:
            output_log = open('{}/import_logger.log'.format(CWD),'w')
            output_log.write('ggkInspector log'+'\n\n')
            first_run = False
        #append to log:
        output_log = open('{}/import_logger.log'.format(CWD),'a')
        output_log.write(string)
        output_log.close()

def recorder(input_string, disburse='both'):
    """Handles message recording based on verbosity"""
    if disburse == 'both':
        if args.verbose:
            print(input_string)
        if args.log:
            log_writter(input_string)
    if disburse == 'display' and args.verbose:
        print(input_string)
    if disburse == 'log' and args.log:
        log_writter(input_string)

Error_DF, k = pd.DataFrame(columns={'Error_Type','Error_Source'}), 0 # logs errors for appending to log
def error_recorder(error_type,error_source):
    global Error_DF, k
    Error_DF = Error_DF.append(pd.DataFrame({'Error_Source':error_source,'Error_Type':error_type},index=[k])).reset_index(drop=True)
    k+=1

def terminator(complete=False):
    """Decides to terminate program if errors have accumalated"""
    if Error_DF.shape[0] > 0:
        recorder('\n------> Exiting Program -- Fix Errors Below and re-run <------\n')
        recorder(Error_DF,'display')
        sys.exit()
    elif complete==True:
        recorder('\n------> No Errors Detected! Exiting Program <------\n')
        recorder('                  Have a great day!')
        sys.exit()
    
#Workflow Functions----------------------------------


def empty_value_checker(DF):
    """Check for empty values in columns"""
    recorder('-'*25+'Checking for Empty Values'+'-'*25+'\n')
    name_error = False
    bool_DF = pd.isnull(DF).any().reset_index().rename({'index':'column_name',0:'boolean_value'},axis=1)
    for i, row in bool_DF.iterrows():
        #print(i)
        if row.boolean_value:
            error_recorder('Missing value','{}th row of column: {}'.format(i,row.column_name))
            recorder('Missing value in column: "{}"'.format(row.column_name),'display')
            name_error = True
    if name_error == False:
        recorder('No empty values!')

        
def slug_checker(DF):
    slug_series=DF['slug* - project names as it should appear in the URL']
    """Verify no invalid characters in slug"""
    recorder('\n'+'-'*25+'Verifying Slug Names'+'-'*25+'\n')
    slug_error = False
    for invalid_character in ['/','\.','\*']:
        if DF['slug* - project names as it should appear in the URL'].str.contains(invalid_character).any():
            error_value = 'Value Error: "{}" is not a valid character in a slug'.format(invalid_character)
            error_sources = (DF[DF['slug* - project names as it should appear in the URL'].str.contains(invalid_character)]\
                    ['slug* - project names as it should appear in the URL'])
            error_recorder(error_value,error_sources)
            recorder(error_value)
            recorder(error_sources)
            slug_error = True
    if slug_error == False:
        recorder('Valid slug names!')
            
def read_path_checker(DF):
    """Check Read Paths to see if valid"""
    recorder('\n'+'-'*25+'Verifying Read Paths'+'-'*25+'\n')
    read_path_error = False
    read_cols = ['read_file_1* -- path to read file *P.E.1','read_file_2(*) -- path to read file *P.E.2']
    for col in read_cols:
        for cell in DF[col]:
            if os.path.isfile(cell) == False:
                error_recorder('Read Path Error','Path: "{}", Column: "{}"'.format(cell,col))
                recorder('Read Path Error: "{}"'.format(cell))
                read_path_error = True
    if read_path_error == False:            
        recorder('Read Paths Valid!')

def wildcard_fixer(raw_path):
    """Accounts for paths without a wildcard, returns with * if necessary"""
    if bool(re.search('\*$',raw_path)) == False:
        raw_path = raw_path+'*' 
        #print('Path Corrected: {}'.format(raw_path))
    return raw_path
        
def assembly_path_checker(DF,min_files=15):
    """Check Assembly Paths, proper files included and not empty"""
    recorder('\n'+'-'*25+'Checking Assemblies'+'-'*25+'\n')
    assembly_path_error = False
    for assembly_path in DF['assembly basename path*']:
        used_assembly_path = wildcard_fixer(assembly_path) #dont want to report mutated string            
        list_of_ass_paths = glob.glob(used_assembly_path)
        if len(list_of_ass_paths) == 0:
            recorder('Assembly Path Error: "{}"'.format(assembly_path))
            error_recorder('Assembly Path Error','Path: "{}"'.format(assembly_path))
            assembly_path_error = True
        elif len(list_of_ass_paths) <= min_files:
            recorder('Missing Files: {}'.format(assembly_path))
            error_recorder('Missing Files in Assembly','Path: "{}"'.format(assembly_path))
            assembly_path_error = True
        elif len(list_of_ass_paths) >= min_files: #Check if files full 
            for file in list_of_ass_paths:
                if os.path.getsize(glob.glob(file)[0]) == 0:
                    recorder('Empty Files: {}'.format(file))
                    error_recorder('Empty Files in Assembly','Path: "{}"'.format(file))
                    assembly_path_error = True
    if assembly_path_error == False:            
        recorder('Assemblies Paths Valid!')

def fasta_grabber(Assembly_path):
    """Returns the fasta file out of the assembly paths independent of file extension naming (assumes fasta will have the least extensions). Helper function for fasta_peeper"""
    assembly_path_DF = pd.DataFrame(glob.glob(Assembly_path)).rename({0:'path'},axis=1)
    assembly_path_DF['file'] = assembly_path_DF.path.apply(lambda x: os.path.basename(x))
    assembly_path_DF['extension_number'] = assembly_path_DF.file.str.split('.').apply(lambda x: len(x))
    fasta_DF = assembly_path_DF[assembly_path_DF.extension_number == assembly_path_DF.extension_number.min()]
    returned_fasta_path = fasta_DF.path.iloc[0]
    if fasta_DF.shape[0]>1:
        recorder('Detecting multiple fasta files: {}'.format(', '.join(list(fasta_DF.file))))
    if fasta_DF.shape[0]==0:
        recorder('No Fasta Detected!') #Buff this out - ask for user input
        returned_fasta_path = False
    #recorder('Checking the following fasta file for slug matching: {}'.format(returned_fasta_path)) #Verbose filter
    return returned_fasta_path

def sed_presciption(DF):
    """Scripts sed cmd for Slug Mismatch Errors. Helper function for fasta_peeper"""
    output_cmd_loc = '{}/ggkInspector_sed.cmd'.format(CWD)
    output_cmd = open(output_cmd_loc,'w')
    recorder('\n---> Writing sed cmds to fix mismatch error. Run the following command in terminal: "bash {}"\nDouble check cmd before running. This may not resolve the issue for all the files.'.format(output_cmd_loc))
    #TO DO: Fix current_header to pull Basename - run Seq.IO and find shared commonalities 
    for i,row in DF.iterrows():
        cmd = "for i in {}; do sed -i 's/{}/{}/g' $i; done"\
                    .format(row.assembly_path,row.current_header,row.slug) #TO DO: SWAP current_header with BASENAME
        #print(cmd)  
        output_cmd.write(cmd+'\n')
    output_cmd.close()
    
    
def fasta_peeper(DF):
    """Check that fasta files header matches slug and reports coverage - check only first header. If there is a Slug Mismatch Error scripts sed cmd"""
    assembly_error, sed_remedy = False, False
    temp_DF = DF[['slug* - project names as it should appear in the URL','assembly basename path*']]\
                .rename({'slug* - project names as it should appear in the URL':'slug','assembly basename path*':'assembly_path'},axis=1)
    temp_DF['assembly_path'] = temp_DF.assembly_path.apply(lambda x: wildcard_fixer(x))
    for i,row in temp_DF.iterrows():
        fasta_path = fasta_grabber(row.assembly_path) 
        temp_DF.loc[i,'fasta_loc'] = fasta_path
        first_header = next(SeqIO.parse(fasta_path, "fasta")).description
        if row.slug in first_header:
            #print("Fasta: {}, Header: {}".format(fasta_path,first_header)) # Might want to include with verbosity filter
            pass
        elif row.slug not in first_header:
            recorder('Slug Mismatch Error: The slug is not detected in the header of {}'.format(fasta_path))
            error_recorder('Slug Mismatch Error: The slug is not detected in the fasta header.','Slug: "{}", First header: "{}", Fasta Path: "{}"'.format(fasta_path,first_header,fasta_path))
            temp_DF.loc[i,'current_header'] = first_header
            assembly_error, sed_remedy = True, True
            
        if "read_length" in first_header and "read_count" in first_header:
            #print('Mapping info included') # Might want to include with verbosity filter
            pass
        elif "read_length" not in first_header and "read_count" not in first_header:
            print('Read Mapping Formatting Error') 
            recorder('Read Mapping Error: The read mapping info is not detected in the fasta header: "{}"')
            error_recorder('Read Mapping Error: The read mapping info is not detected in the fasta header.','Fasta Path: "{}"'.format(fasta_path,first_header,fasta_path))
            assembly_error = True
    
    if sed_remedy:
        sed_presciption(temp_DF.dropna())
        
    if assembly_error == False:            
        recorder('Slugs names match header!')
        
        
        
#Task Manager---------------------------------------

canon_column_name = ['name*', 'slug* - project names as it should appear in the URL',
       'project_group* (url / slug of existing project group or name of new project group)',
       'description', 'date_collected*', 'location*', 'sequencing_facility*',
       'read_length', 'read_file_1* -- path to read file *P.E.1',
       'read_file_2(*) -- path to read file *P.E.2', 'read_processing*',
       'total_read_bp', 'assembly_type*', 'total_assembled_bp',
       'assembly basename path*']

def main(args):
    '''Main entry point to program'''
    #Assert column names:
    assert (input_DF.columns == canon_column_name).all(), 'Error: Invalid Column Names. Follow this sheet: https://ggkbase-help.berkeley.edu/wp-content/uploads/2015/11/ggkbase_bulk_upload_template.tsv'
    
    #Naming Convention Check:
    empty_value_checker(input_DF)
    slug_checker(input_DF)
    terminator()
    
    #Check Paths:
    read_path_checker(input_DF)
    assembly_path_checker(input_DF,12) #make an input arg?
    terminator()
    
    #Verify Slug Header:
    fasta_peeper(input_DF) 
    terminator(complete=True)
    
if __name__ == '__main__':
    main(args)
