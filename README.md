**Data files & variables**


-   ***data_Exp1.mat*** contains data from 20 participants from Exp. 1: Choice-reaction time task with color judgments and difficulty judgments 


-   ID: subject ID
-   Task: 'Color' = Color judgment task (Training); 'Difficulty' = Difficulty judgment task
-   sColCoh1 and sColCoh2: signed color coherence of S1 and S2 (only in difficulty judgment task), neg = yellow; pos = blue
-   Color1 and Color2: color dominance of S1 and S2 (only in difficulty judgment task), 0 = yellow; 1 = blue
-   Choice: In Color judgment task: 0 = yellow; 1 = blue; In Difficulty judgment task: 0 = S2; 1 = S1
-   Accuracy: 0 = error; 1 = correct
-   RT: reaction time in seconds


-   ***data_Exp2_RT.mat*** contains data from 3 participants from Exp. 2: Choice reaction time task of difficulty judgments with known vs. unknown color


-   Same variable names/coding as in data_Exp1.mat
-   Cond: 1 = unknown color; 2 = known color (also see variable CondLabel)
-   sDiffCoh: signed difference in color coherence = abs(sColCoh1)-abs(sColCoh2)


-   ***data_Exp2_VD.mat* **contains data from 3 participants from Exp. 2: Variable duration task of difficulty judgments with known vs. unknown color


-   Same variable names/coding as in data_Exp1.mat and data_Exp2_RT.mat
-   stimDur: stimulus duration in ms


**Reproduce figures in paper**


-   Each folder called 'Fig*' contains a main script called run_fig*.m that recreates the corresponding figure in the paper
-   For figures showing model fits, the script will by default load saved model fits that are located in the subfolders called 'fits'
-   Some folders also contain a function called fit_difficulty*.m. These are optional functions that aren't called in the main scripts, but can be executed independently to re-fit the models (see below).


**Refit models**


-   Execute fit_difficulty*.m functions to create new fits. The new fits will be saved in the subfolder 'fits' with a different name (e.g., fits_ID1_new.mat)
-   If you want to plot the new fits, make sure you change the files that are being loaded in the main script run_fig*.m
-   Before refitting the models: Check fitting parameters carefully! E.g., the number of iterations (= 1 by default for faster execution; 10 in the paper), starting values & upper/lower bound of each parameter etc.
-   ***The model fits may take several days to run (even for a single iteration).*** Change code to execute parallel loops according to set-up.


**Additional functions**


-   Most scripts will require some of the custom functions in the 'functions' folder. 
-   The subfolder 'Fit_Difficulty_ChoiceRT' contains functions used for fitting the 4 models of difficulty choice-RT (see main Figures 3 and 6)
-   The subfolder 'Fit_Color_ChoiceRT' contains functions used for fitting color choice-RT for Exp. 1 (see main Figure 1b)