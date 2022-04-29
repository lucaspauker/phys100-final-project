# phys100-final-project

## HOW TO RUN THE CODE

### CONDA
#### Set up conda env for the first time
```
conda env create -f env.yml
```
#### Update conda env when you pull code
```
conda env update --file env.yml --prune
```
#### Do this before committing
```
conda env export --no-builds > env.yml
```
#### Activate environment
Let's make sure we are using the same conda environment
```
conda activate phys100py3
```

### GITHUB
#### Pull code (do this when you want most updated code)
Also, make sure you pull master before pushing new code.
```
git pull origin master
```
#### Push your new code
The idea is to create a new branch, commit to the branch, push the branch, then merge it into the repo.
This is good since we can make sure that there are no breaking changes and is good practice in general.
```
git checkout -b mybranch
git add FILE1 FILE2 ...
git commit -m "What your changes are"
git push origin mybranch
```
