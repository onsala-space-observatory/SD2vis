# SD2VIS

## Installation

Steps to install the `SD2vis` task into `casa`

 1. Clone the git repository into a directory of your choice
 (e.g., $HOME/.casa/NordicTools)

``` shell
cd $HOME/.casa/NordicTools
git clone <repository url>
cd SD2vis
buildmytasks --module fakeobs SD2vis.xml
```
 2. Inside `casa` add the folder to your `PYTHONPATH`:

``` python
CASA <1>: sys.path.insert(0, <path to SD2vis folder>)
CASA <2>: from SD2vis.gotasks.SD2vis import SD2vis
CASA <3>: inp(SD2vis)

```

That's it! Enjoy!
