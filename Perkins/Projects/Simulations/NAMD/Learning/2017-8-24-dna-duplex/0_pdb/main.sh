AMBERHOME="/usr/lib/amber16"
NABHOME="${AMBERHOME}/dat"
PATH="$PATH:$AMBERHOME/bin"

nab -o nab.out generate.nab
# generate the pdb and clean up after...
./nab.out ; rm nab.out generate.c
