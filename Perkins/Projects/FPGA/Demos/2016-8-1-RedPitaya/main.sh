#main.sh: automates creationg of build enviromnent, following instructions from:
#  https://github.com/RedPitaya/RedPitaya/blob/master/README.md, asusming 
# that 'Software requirements' section has been cmopleted

# Args:
# $1: Build directory. Default, assumes home
# $2: Vivado Binary Path, assumes 2016.2 by default
WorkingDir=${1:-~}
VivadoBin=${2:-"/opt/Xilinx/Vivado/2016.2/bin/"}
#  Currently, the configuration file assumes that 'patrick' and 'perkins'
# are the two users needed. If you have a different username, add it to the
# config file "schroot_config" to the relevant csv-separated lists

# note: you *must* have downloaded Xilinx (the 'setup' script can be run to 
# create the directories needed
# generic dependencies
sudo apt-get install make curl xz-utils
# U-Boot build dependencies
sudo apt-get install libssl-dev device-tree-compiler u-boot-tools
#Vivado requires a gmake executable which does not exist on Ubuntu. It is necessary to create a symbolic link to the regular make executable.
sudo ln -s /usr/bin/make /usr/bin/gmake
# remember where we ran this from, to get the config file
Dir=$PWD
# build in the home directory
cd $WorkingDir
git clone https://github.com/RedPitaya/RedPitaya.git 
cd RedPitaya
# set up environmnetal variables
sudo bash settings.sh
# add vivado to path... just in case
export PATH=$PATH:$VivadoBin
# get ready for downloads
mkdir -p dl
export DL=$PWD/dl
# download everything
wget http://downloads.redpitaya.com/ubuntu/redpitaya_ubuntu-latest.tar.gz
sudo chown root:root redpitaya_ubuntu-latest.tar.gz
sudo chmod 664 redpitaya_ubuntu-latest.tar.gz
# create the configuraton file we need in this directory 
MagicFile="/etc/schroot/chroot.d/red-pitaya-ubuntu.conf"
ConfigFile="schroot_config"
sudo cp $Dir/$ConfigFile $MagicFile
sudo chmod 777 $MagicFile
# build everything
sudo make -f Makefile.x86
schroot -c red-pitaya-ubuntu <<- EOL_CHROOT
make -f Makefile.arm CROSS_COMPILE="" REVISION=$GIT_COMMIT_SHORT
EOL_CHROOT
sudo make zip
