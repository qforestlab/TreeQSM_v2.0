# Setup
TreeQSM v2.0 setup instruction manual (For Windows10)
<br> <br>
**1.	Setup WSL2 (Windows Subsystem for Linux 2)**
<br> <br>
a.	Download WSL2 by following the instructions from: 
https://learn.microsoft.com/en-us/windows/wsl/install
<br>
WSL default installation will allow you to run Ubuntu on a Windows machine. Keep in mind that as long as WSL is running, a portion of system RAM will be reserved (you can see this as Vmmem in task manager processes)
<br> <br>
b.	Launch WSL by typing in ‘WSL’ to the windows search bar. You will be prompted to create a username and password that is completely independent of your windows user account. Follow these steps for more information on setting up a WSL account: https://learn.microsoft.com/en-us/windows/wsl/setup/environment#set-up-your-linux-username-and-password
<br> <br>
c.	Set up directories in the Linux subsystem for TreeQSM. Navigate to your user folder and create a TreeQSM/v2 directory by running:
<br> <br>
cd ~
<br>
sudo mkdir TreeQSM
<br> <br>
**2.	X11 Forwarding (required to activate Matlab through the gui)**
<br>
There is a great guide to have X11 forwarding work with WSL2 and automatically adapt to WSL IP address changes:
https://thomasward.com/wsl2-x11/
<br>
Recommendation: under the section “Configure Windows Firewall” try using the second method listed as manually adding the rule through a GUI instead of running the WSL2 script.
<br><br>
**3.	Setup Matlab:**
<br><br>
a.	 In WSL, run:
<br><br>
sudo apt-get install –y unzip
<br>
wget https://www.mathworks.com/mpm/glnxa64/mpm
<br>
chmod +x mpm
<br>
sudo ./mpm install --release=R2022a --destination=/usr/local/MATLAB/R2021a --products MATLAB Parallel_Computing_Toolbox Statistics_and_Machine_Learning_Toolbox
<br><br>
This will download and install the matlab package manager, and then will install matlab with the two required packages to run TreeQSM v2. Feel free to change the destination location if desired, but remember where it is because we will need the matlab path later. If you run into download issues, check that your vpn is not on. If you run into folder permission issues, run the sudo command in front of ./mpm.
<br><br>
b.	Download ResearchLicense.dat from the Q-ForestLab github or teams Laser Scanning Group Software_Code channel. Place it somewhere easy to find on WSL, like in your TreeQSM folder.
<br><br>
c.	Next, run:
<br><br>
cd ~
<br>
cd */usr/local/MATLAB/R2022a/bin*
<br>
sudo ./activate_matlab.sh
<br><br>
Note that in the previous code block, the *italicized* pathway needs to match your matlab installation path.
<br><br>
d.	The MathWorks Software Activation gui should be open now. It will look like this:
![image](https://github.com/user-attachments/assets/005ee861-ae02-43d0-9bd8-e72b73cdecf3)
<br>
Select “Activate manually without the Internet”
<br>
Click on the “Browse…” button, navigate to, and select the ReseachLicense.dat file you previously downloaded. Hit next, and the license should be activated. If you run into a file permission error, it’s most likely because you ran ./activate_matlab.sh without sudo.
<br><br>
**4.	Setup QSM scripts**
<br><br>
a.	A complete folder setup can be found in the github or in the Laser Scanning Group teams Software_Code channel. Look for v2.7z and download.
<br><br>
b.	You can extract the files to your TreeQSM folder in your Linux subsystem.
<br><br>
c.	There are a few file paths in the QSM python and batch scripts that will need to be double checked before running:
<br><br>
i.	In batch_run_model.py, two file paths need to match your location for runCylinderModel.py (in your v2 folder) and the location of your matlab run script (found where you installed matlab). See below:
![image](https://github.com/user-attachments/assets/85aed8fb-616d-4d13-9a7e-fc3f3ac11b80)
<br> <br>
ii.	In runCylinderModel.py, one file path marked below needs to match the location of your matlab run script. See below:
![image](https://github.com/user-attachments/assets/33e2049a-3018-4673-b2d8-b857f8d176a2)
<br><br>
iii.	In runCylinderModel.m (different than runCylinderModel.py), two file paths need to match the location of your v2 file that contains all the v2 QSM scripts. Note that the second file path ends with a ‘/’ and the first does not. 
![image](https://github.com/user-attachments/assets/3a52c9dc-8a00-4665-b7a1-3e51c420d2af)
<br><br>
Setup is complete.
<br><br>


# Run QSMs Locally
1.	Put any point clouds to process in v2/pointcloud
<br><br>
2.	Run: 
<br><br>
python batch_run_model.py i
<br><br>
Where ‘i’ is the amount of cores you want to dedicate to parallel processing. When in doubt just test it with 1. The QSMs should be placed in v2/results, alond with .tracker files. If batch_run_model.py says it is starting processing on pointclouds but no processing occurs, it’s most likely because .tracker files for those trees already exist in v2/results, and removing those .tracker files will fix this problem.
<br><br>
3.	Next, to optimize based on generated QSMs, run:
<br><br>
python optimiseCylinderModel.py -c pointcloud/ -m results/
<br><br>
4.	The first generated QSM from optimal parameters for each point cloud can now be found in opt_1mod. Note that v2 .mat files for QSM output are structured differently than later TreeQSM versions, but can still be summarized and graphed using ITSMe.


# Run QSMs on Cave013
1. First check if you have an account, if not ask someone to create one for you.
<br><br>
2. You can use the matlab installation in Zane’s folder, so matlab doesn’t need to be installed. 
<br><br>
3. Copy Zane’s TreeQSM folder and remove all the .txt, .mat and .tracker files in the pointcloud, results, opt, opt_param, and opt_1mod files.
<br><br>
4. Follow step 4: Setup QSM scripts, from the setup section above.
<br><br>
5. Copy your pointclouds from your local machine to your “pointcloud” folder on CAVE013. This can be done remotely with the following command (recommended to read about scp file transfering):
<br><br>
scp -r -P 2225 /C:/wherever/you/keep/your/pointclouds <username>@cave013.ugent.be:/Stor2/<username>/TreeQSM/v2/
<br><br>
Note that in the previous code block, you must substitute your own file pathways and username
<br><br>
6. To connect to cave013 remotely through any terminal, make sure you have the university VPN enabled, then use:
<br><br>
ssh -p 2225 <username>@cave013.ugent.be
<br><br>
7. You can start a new 'screen' to ensure your session will continue when you disconnect from Cave013. **Keep in mind your screens will run indefinitely until you shut them down, potentially reserving computing resources**:
<br><br>
a. Creating a screen:
<br>
screen –S ‘screenname’
<br><br>
Detaching from a screen:
<br>
crtl+a+d
<br><br>
Display which screens are running:
<br>
screen –ls
<br><br>
Reattach to a screen:
<br>
screen –r ‘screenname’
<br><br>
Shutdown a screen:
<br>
screen –X –S ‘screenname’ quit
<br><br>
8. To copy the results to your own computer again, close the session with the Cave013 and run the following:
<br><br>
scp -r -P 2225 <username>@cave013.ugent.be:/Stor2/<username>/TreeQSM/v2/results C:/yourfilepath/TreeQSM/v2/results
<br><br>
You can download extra folders such as “opt”, or “results” by tweaking the file path in the previous line of code
<br><br>
9. When the data is transferred back to your local machine, you can run the quality control and image generation locally.
<br><br>
OTHER STUFF TO KNOW:
<br>
•	When you want to run a new batch of QSMs, you can store all the old point clouds/QSMs in the storage folder.
<br>
•	Every point cloud gets a tracker, and *once it’s there it will skip this point cloud when processing*. Remove trackers if you want to reprocess a pointcloud.
<br>
•	Don’t forget to shut down your screens
<br>
•	It's good QSM practice to save all results and optimization parameters, especially if it's a dataset for publication
<br><br>

# Further questions?
Contact Zane Cooper
<br>
Zane.Cooper@Ugent.be
