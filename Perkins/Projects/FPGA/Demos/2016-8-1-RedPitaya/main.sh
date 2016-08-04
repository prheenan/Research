Dir=`pwd`
cd ~
cd RedPitaya
# source the settings
. settings.sh
# create the configuraton file we need in this directory 
# build everything
sudo make -f Makefile.x86
schroot -c red-pitaya-ubuntu <<- EOL_CHROOT
make -f Makefile.arm CROSS_COMPILE="" REVISION=$GIT_COMMIT_SHORT
EOL_CHROOT
sudo make zip
