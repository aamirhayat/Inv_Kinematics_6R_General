# Bash file to make the executable of Inverse kinematics code for Generalized 6R Robot
# Color output in bash https://stackoverflow.com/questions/5947742/how-to-change-the-output-color-of-echo-in-linux

echo -e "\e[1;36mInput to this code are DH, Rotation Matrix and Position listed in /IK_CODE/input.dat\e[0m"
echo -e "\e[0;34ma1 a2 ... a6\e[0m"
echo -e "\e[0;34md1 d2 ... d6\e[0m"
echo -e "\e[0;34malpha1 alpha2 ... alpha6\e[0m"
echo -e "\e[0;34mR11 R12 R13\e[0m"
echo -e "\e[0;34mR21 R22 R23\e[0m"
echo -e "\e[0;34mR31 R32 R33\e[0m"
echo -e "\e[0;34mT1  T2  T3\e[0m"
# By A. A. Hayat
echo -e "\e[0;37m****** General IKin by Dr. Dinesh Manocha*****\e[0m"
echo -e "\e[1;33m****** Run .sh by by A. A. Hayat ***********\e[0m"
echo -e "Make Sure that the address of folders EISPACK and IK_CODE is changed in ikin_runme.sh"
echo "*******************************************************"

read -p "Press ENTER to continue"

# Rename the file locations of EISPACK and IK_CODE as in your system

# Delete the .obj file from the folder location first
cd /home/iit/6R_G/EISPACK && rm -r *.o ; cd -
cd /home/iit/6R_G/IK_CODE  && rm -r *.o ; cd -
echo ".obj files deleted from the EISPACK and IK_CODE"
echo "*******************************************************"
read -p "Press ENTER to continue"
# Make file .obj of the fortran files
cd /home/iit/6R_G/EISPACK && gfortran -c * ; cd -
cd /home/iit/6R_G/IK_CODE && gfortran -c * ; cd -
# alias short='~/home/iit/6R_G/IK_CODE' & 
# Output executable file name is General6R_ikin_AH and it can be changed 
cd /home/iit/6R_G/IK_CODE && gfortran -o General6R_ikin_AH \
/home/iit/6R_G/IK_CODE/main.o \
/home/iit/6R_G/IK_CODE/matrix.o \
/home/iit/6R_G/IK_CODE/read.o \
/home/iit/6R_G/IK_CODE/read_ip.o \
/home/iit/6R_G/IK_CODE/eigen1.o \
/home/iit/6R_G/IK_CODE/improve.o \
/home/iit/6R_G/IK_CODE/reduce.o \
/home/iit/6R_G/IK_CODE/degen.o \
/home/iit/6R_G/IK_CODE/setup.o \
/home/iit/6R_G/IK_CODE/setup2.o \
/home/iit/6R_G/IK_CODE/check_degen.o \
/home/iit/6R_G/IK_CODE/gauss.o \
/home/iit/6R_G/EISPACK/dgeco.o \
/home/iit/6R_G/EISPACK/dgedi.o \
/home/iit/6R_G/EISPACK/rg.o \
/home/iit/6R_G/EISPACK/rgg.o \
/home/iit/6R_G/EISPACK/hybrj.o \
/home/iit/6R_G/EISPACK/dsvdc.o \
/home/iit/6R_G/EISPACK/dgesl.o ; cd -

echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
echo "General6R_ikin_AH is the executable | Run it from terminal using ./General6R_ikin_AH"
echo "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"

# Delete the .obj file from the folder location first
cd /home/iit/6R_G/EISPACK && rm -r *.o ; cd -
cd /home/iit/6R_G/IK_CODE  && rm -r *.o ; cd -
# cd /home/iit/6R_G/IK_CODE
# exec bash
cp /home/iit/6R_G/IK_CODE/General6R_ikin_AH /home/iit/6R_G/Z_Exe_only



