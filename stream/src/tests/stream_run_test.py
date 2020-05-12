import subprocess
import os 
from pathlib import Path
import filecmp
import tarfile
import tempfile
from shutil import rmtree
import sys
import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np

_root = os.path.abspath(os.path.dirname(__file__))

def stream_test_Nestorowa_2016():

	workdir = os.path.join(_root,'datasets/Nestorowa_2016/')

	temp_folder=tempfile.gettempdir()

	tar = tarfile.open(workdir+'output/stream_result.tar.gz')
	tar.extractall(path=temp_folder)
	tar.close()
	ref_temp_folder=os.path.join(temp_folder,'stream_result')

											 
	print(workdir+'data_Nestorowa.tsv.gz')
	#subprocess.call(['stream','-m', workdir+'data_Nestorowa.tsv.gz', '-l', workdir+'cell_label.tsv.gz', '-c', workdir+'cell_label_color.tsv.gz', '-g', workdir+'gene_list.tsv.gz', '--TG,', '--DE', '--LG','-o',workdir+'stream_result/'])

	input_file=os.path.join(workdir,'data_Nestorowa.tsv.gz')
	label_file=os.path.join(workdir,'cell_label.tsv.gz')
	label_color_file=os.path.join(workdir,'cell_label_color.tsv.gz')
	gene_list_file = os.path.join(workdir,'gene_list.tsv.gz')
	comp_temp_folder = os.path.join(temp_folder,'stream_result_comp')

	stream_cmd='stream -m {0} -l {1} -c {2} -g {3} -o {4} --dr_method mlle --nb_pct 0.1 --clustering ap --TG --DE --LG --n_jobs 4'.format(input_file,label_file,label_color_file,gene_list_file,comp_temp_folder)
	print(stream_cmd)
	code = subprocess.call(stream_cmd,shell=True) 
	if(code!=0):
		print('STREAM running failed.')
		sys.exit(1)

	print(ref_temp_folder)
	print(comp_temp_folder)

	pathlist = Path(ref_temp_folder)
	for path in pathlist.glob('**/*'):
	    if path.is_file() and (not path.name.startswith('.')):
	        file = os.path.relpath(str(path),ref_temp_folder)
	        if(file.endswith('pdf')):
	            if(os.path.getsize(os.path.join(comp_temp_folder,file))>0):
	                print('The file %s passed' %file)
	            else:
	                print('Error! The file %s is not matched' %file)
	                sys.exit(1)
	        else:
	            checklist = list()
	            df_ref = pd.read_csv(os.path.join(ref_temp_folder,file),sep='\t')
	            df_comp = pd.read_csv(os.path.join(comp_temp_folder,file),sep='\t')
	            for c in df_ref.columns:
	                if(is_numeric_dtype(df_ref[c])):
	                    checklist.append(all(np.isclose(df_ref[c],df_comp[c])))
	                else:
	                    checklist.append(all(df_ref[c]==df_comp[c]))
	            if(all(checklist)):
	                print('The file %s passed' %file)
	            else:
	                print('Error! The file %s is not matched' %file)
	                sys.exit(1)

	print('Successful!')

	rmtree(comp_temp_folder,ignore_errors=True)
	rmtree(ref_temp_folder,ignore_errors=True)


def main():
	stream_test_Nestorowa_2016()
	

if __name__ == "__main__":
    main()