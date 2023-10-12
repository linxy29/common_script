## only one command is enough
wget -c --user="jwkho@hku.hk" --password="HoLab2013!"  --mirror ftp://human.big.ac.cn/HRA000183

#!/bin/bash
#target=HRA000183
#for i in `cat HRA000183.txt|tr  "\t" "/"`;do 
#echo ${i};
#wget -c --user="jwkho@hku.hk" --password="HoLab2013!"  --mirror ftp://human.big.ac.cn/${target}/${i}
#done