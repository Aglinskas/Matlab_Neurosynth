# Matlab_Neurosynth
Runs meta-analysis on neurosynth database in Matlab

.. work still in progress. 


What it does: 
###- Part 1 of the code generates a list of words;
iterate by changing to words to be included and exluded(1) and pressing CMD+ENTER till you're satisfied with the 'FINAL LIST'(2)

<img width="812" alt="napkin 03-02-17 4 39 00 pm" src="https://cloud.githubusercontent.com/assets/15108226/22597335/7762ec0a-ea2f-11e6-9633-0fcd2180f69a.png">

### - Part 2 of the code does the following; 
- Searches the neurosynth _features_ database for studies that 'load' (as defined by neurosynth) on those words

![image](https://cloud.githubusercontent.com/assets/15108226/22597895/c8c8493a-ea31-11e6-85cd-320e8dd96847.png)


![screen shot 2017-02-03 at 16 44 03](https://cloud.githubusercontent.com/assets/15108226/22597961/0be27e66-ea32-11e6-9818-957678933d5f.jpg)


- takes the voxels reported by those studies and plots them onto an empty brain by adding +1 every time that coordinate gets reported; 
for example if exact coordinate [42 -46 -22] get reported three times, it has value of 3; [42 -46 -23] get reported 1 time, it gets value of 1; 
- smooths result with kernel chosen.
- exact file gets exported as meta.nii
-  smoothed file gets exported as smeta.nii

<img width="1098" alt="napkin 03-02-17 5 05 06 pm" src="https://cloud.githubusercontent.com/assets/15108226/22598204/f1b6aaca-ea32-11e6-80e4-dfa86ac20487.png">



##How to run: 
*click for video*
[![Running the script](https://cloud.githubusercontent.com/assets/15108226/22597159/c0a849a6-ea2e-11e6-932b-a5fac84f3351.jpg)](https://www.youtube.com/watch?v=3z9eThlg45w "Quick walkthrough")

