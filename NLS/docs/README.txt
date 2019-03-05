-----------------------------
NLS++ Usage Instructions
-----------------------------

1. Open the NLS++ project and application in Visual Studio

- Navigate to tests/prj/vc7/
- Double click tnls.vcproj 
- This will open the Solution 'tnls' which contains 2 projects: nls and tnls
- If you are running a versions of VS other than 2007, you will get a message that 
  the solution needs to be converted, agree to the conversion
- Disable incremental linking on the tnls project by right clicking on the tnls 
  project and select 'Properties’, then select Linker->General->Enable Incremental
  Linking -> “No”

2. Compile and Run the solution either from Visual Studio (a) or from the command line (b)

(a) Compile and Run the solution from Visual Studio

- Add the model and algorithm files you want to run to tests/prj/vc7
- Note:  the algorithm and model files for the unidimensional function 
  are already in the directory)
- Right click on the tnls project and select 'Properties'
- Click 'Debugging' under 'Configuration Properties'
- Enter the names of the model and algorithm files in the 'Command Arguments' box
  (e.g.: model_unidimensional.inp alg_gdcm.inp)
- Click 'Apply' then 'OK'
- Right click on the tnls project and select 'Build Order'
- Make sure that nls is listed before tnls. If not click the 'Dependencies' tab,
  then select 'tnls' under the Projects menu, and click the check box next to 'nls'
- Click 'OK'
- Run in Debug mode by clicking the green play button at the top of the screen
- The output files(s) will be written to the tests/prj/vc7 directory
- The executable file used to run the project from the command line will be written 
  to bin/vc7/tnls.exe
- Note: Additional directories and files were automatically created when the 
  program is compiled, namely bin/, libd/, objd/ and a .user file in tests/prj/vc7
  that saves your configuration settings (e.g. command arguments, build order, etc.)
  
(b) Compile the solution in Visual Studio and Run from the Command Line

- Right click on the tnls project and select 'Build Order'
- Make sure that nls is listed before tnls. If not click the 'Dependencies' tab,
  then select 'tnls' under the Projects menu, and click the check box next to 'nls'
- Click 'OK'
- Right click on the tnls project and select 'Clean Solution'
- Right click on the tnls project and select 'Build Solution'
- The executable file will be written to bin/vc7/tnls.exe
- Add the model and algorithm files you want to run to bin/vc7
- Open a command prompt
- Navigate to bin/vc7
- Run the command: tnls.exe [model_file.inp] [algorithm_file.inp]

  
