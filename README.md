# Inv_Kinematics_6R_General
This project is to get the generalized inverse kinematics solution

## Folder list
  * EISPACK : Contains fortran code
  * IK_CODE : contains c and fortran files
  * Z_exe_only : The .exe file generated with the input.dat and z_output.dat files (Copy the generated executable in this folder)
  * README.odt : Help file to generate the executable

Versions Required for above project:
* gcc
* fortran 77

  ## Terminal Command
```
./ikin_runme.sh
```
   * Ensure to change the folder location in the ikin_runme.sh file

  ### Output
    ~/IK_CODE/General6R_ikin_AH
```
./General6R_ikin_AH
```

### Input to the executable
0.18 0.6 0.12 0 0 0  <br />
0.4 0 0 0.62 0 0.115   <br />
90 180 -90 90 -90 0  <br />

-0.8086   0  0.5883  <br />
0.0   1.0   0.0  <br />
-0.5883   0.0   -0.8086  <br />
0.08   0.1   1.2  <br />

10 20 60 40 50 60  <br />


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  <br />
% Link lengths  <br />
% joint offset  <br />
% Twist anlge  <br />
The above DH parameters are taken from the RoboAnalyzer  and the visualization can be done there for the 8 IKin solutions
Paper reference is: <br />

### Output to the inverse kinematics
* From General6R_ikin_AH
![Output Image](https://github.com/aamirhayat/Inv_Kinematics_6R_General/blob/master/output.png)
* From RoboAnalyzer
![Output Image](https://github.com/aamirhayat/Inv_Kinematics_6R_General/blob/master/image_RA_Ikin.JPG)
