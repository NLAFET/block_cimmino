TODO
=====

## C interface ##

It's easy to implement a C-interface for `ABCD`. Two files are included in `test/src/`, `abcd_c.h` and `c_interface.cc`. The C-interface is a structure called `abcd_solver` and it replicates part (only public) members of the C++ class `abcd`.

### What's done? ###

We have already defined an example in the aforementioned files, the test file is also included `c_test.c` and compiles into an executable `test_abcd_c_interface`. The `icntl` array is the only interfaced member right now, the structure looks like:

    struct abcd_solver
    {
      int *icntl;
    };
    
The initializing function is called `new_solver()` and returns a pointer to an allocated structure:

    struct abcd_solver* new_solver()
    {
      struct abcd_solver *solver = (struct abcd_solver *)malloc(sizeof(*solver));

      obj = new abcd();
      
      solver->icntl = &obj->icntl[0];
      // TODO: put the rest of associations here!

      return solver;
    }
    
In this function we create an instance of the C++ class `abcd` and
associate the pointer to the first entry in the C++ `icntl` vector to
the C-array with the same name. Here is a simple working example in C

    #include "abcd_c.h"
    //....
    typedef struct abcd_solver abcd;
    int i;

    abcd *obj = new_solver();

    for(i = 0; i < 20; i++)
        printf("%d   %d\n", i, obj->icntl[i]);
        
Notice that `abcd` in this case is not related to the C++ class but to the C-structure as the C++ header is not included.

### What's left? ###

Almost everything! 

- Add the rest of associations
- Add the function call to run the solver (it should call the function `bc()`.
