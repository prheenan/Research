# main.sh: after (1) running setup.sh, (2) installing vivado, be in the same
# directory as main.sh and run

# this essentially copies the configuration file, then runs the
# entire build process from the ground up 

Dir=`pwd`
cd ~
cd RedPitaya
# source the settings
. settings.sh
# create the 'magic' configuration file we want
MagicFile="/etc/schroot/chroot.d/red-pitaya-ubuntu.conf"
ConfigFile="schroot_config"
sudo cp $Dir/$ConfigFile $MagicFile
# create the configuraton file we need in this directory 
# build everything
sudo make -f Makefile.x86
schroot -c red-pitaya-ubuntu <<- EOL_CHROOT
make -f Makefile.arm CROSS_COMPILE="" REVISION=$GIT_COMMIT_SHORT
EOL_CHROOT
sudo make zip
