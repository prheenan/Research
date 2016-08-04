#setup.sh: automates setting up build enviromnent, following instructions from:
#  https://github.com/RedPitaya/RedPitaya/blob/master/README.md, asusming 
# that 'Software requirements' section has been cmopleted

#  Currently, the configuration file assumes that 'patrick' and 'perkins'
# are the two users needed. If you have a different username, add it to the
# config file "schroot_config" to the relevant csv-separated lists

# note: you *must* have downloaded Xilinx (the 'setup' script can be run to 
# create the directories needed
# generic dependencies
sudo apt-get install make curl xz-utils
# U-Boot build dependencies
sudo apt-get install libssl-dev device-tree-compiler u-boot-tools
#Vivado requires a gmake executable which does not exist on Ubuntu.
#It is necessary to create a symbolic link to the regular make executable.
sudo ln -s /usr/bin/make /usr/bin/gmake
# for ARM emulation, if we are using a x86 CPU, we need the fllowing package
# (Correspondence with Iztok.Jeras@redpitaya.com, prh, 8/3/2016)
sudo apt-get install qemu qemu-user qemu-user-static
# download schroot for emulation
sudo apt-get install schroot
# remember where we ran this from, to get the config file
Dir=`pwd`
# build in the home directory
cd ~
git clone https://github.com/RedPitaya/RedPitaya.git 
cd RedPitaya
# get ready for downloads
sudo mkdir -p dl
export DL=$PWD/dl
# download everything, change the permissions how we want them to be changed
wget http://downloads.redpitaya.com/ubuntu/redpitaya_ubuntu-latest.tar.gz
sudo chown root:root redpitaya_ubuntu-latest.tar.gz
sudo chmod 664 redpitaya_ubuntu-latest.tar.gz
