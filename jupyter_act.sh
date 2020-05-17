jpy=$1  #/home/ye/anaconda3/envs/r-base/bin/jupyter-notebook
host="10.100.110.103"
nohup $jpy --ip $host  --no-browser --port 8888 1>nb.log 2>&1 &
