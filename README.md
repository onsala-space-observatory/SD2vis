# SD2VIS

## Installation

Steps to install the `SD2vis` task into `casa`

 1. Clone the git repository into a directory of your choice
 (e.g., $HOME/.casa/NordicTools)

``` shell
cd $HOME/.casa/NordicTools
git clone <repository url>
cd SD2vis
buildmytasks
```
 2. Edit the file `$HOME/.casa/init.py`. Add the line:

``` shell
execfile('$HOME/.casa/NordicTools/SD2vis/mytasks.py')
```

That's it! You should be able to run the new task in CASA! Just doing:

``` shell
tget SD2vis
```

inside `casa` should load the task. To get help, just type `help SD2vis`
