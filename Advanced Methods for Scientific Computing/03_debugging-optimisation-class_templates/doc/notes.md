# Debugging

In order debug, we need to deactivate all optimizations as they can cause change in code.\
Thus, we will have for compilation
```g++ -O0 -ggdb <filename>```\
Now, we can enter `gdb <filename>` which opens an interactive shell. This will automatically open the first error.\
By typing `up`, we will get a recursive trace rollback of the error.\
`Ctrl+x` and `Ctrl+d` are used to get out.

# Class templates
We you have templates, it is mandatory to put definitions in header files.\
Since this can create large and messy files, we can split it; in order to have just a single header file to be included,
it is a good idea to put an `#INCLUDE` at the end of the main header file.