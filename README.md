# psCRN-simulator
A compiler and simulator for the verification of psCRN protocols, as introduced and detailed in [Tai A., Condon A. (2019) Error-Free Stable Computation with Polymer-Supplemented Chemical Reaction Networks. In: Thachuk C., Liu Y. (eds) DNA Computing and Molecular Programming. DNA 2019. Lecture Notes in Computer Science, vol 11648. Springer, Cham](https://link.springer.com/chapter/10.1007/978-3-030-26807-7_11)

# Introduction to psCRNs
The psCRN model is a mathematical model describing how monomers capable of forming polymers interact in a well-mixed test tubes. As many psCRN protocols rely heavily on program counters, the number of reactions, reactants, or products can easily number in the hundreds or thousands. Therefore, we also designed a simple programming language to facilitate reasoning about them, and wrote this simulator to take advantage of that.

# Understanding the psCRN coding language
The psCRN coding language relies on atomic instructions, each which corresponds to a set of reactions. A psCRN program is a set of these instructions, each one which must finish executing before moving onto the next. Therefore, each instruction in the program is associated with an instruction number *i*, and a unique "program counter species" *L_i*. In our paper, an example set of atomic instructions we gave were:

```
instruction: inc(sigma)
    L_{i} + B -> L_{i}* + {sigma}
    L_{i}* + A_{sigma} -> L_{i+1} + I_{sigma}

instruction: dec(sigma)
    L_{i} + I_{sigma} -> L_{i}* + A_{sigma} 
    L_{i}* + {sigma} -> L_{i+1} + B

instruction: jump-if-empty(sigma, k)
    L_{i} + bot_{sigma} -> L_{k} + bot_{sigma}
    L_{i} + I_{sigma} -> L_{i+1} + I_{sigma}

instruction: goto(k)
    L_{i} + B -> L_{k} + B

instruction: create(sigma)
    L_{i} + B -> L_{i+1} + {sigma}

instruction: destroy(sigma)
    L_{i} + {sigma} -> L_{i+1} + B
```

Using only these instructions we can build more complicated functions, for example:

```
function: add-to(sigma, sigma')
i:  goto(i.1)
i.1:jump-if-empty(sigma, i.6)
i.2:    dec(sigma)
i.3:    inc(sigma')
i.4:    inc(copy)
i.5:    goto(i.1)
i.6:jump-if-empty(copy, i.10)
i.7:    dec(copy)
i.8:    inc(sigma)
i.9:    goto(i.6)
i.10:goto(i+1)
```

Then of course, we can use these functions and instructions together to build a complete psCRN program, without ever writing down a reaction:

  ```
  1: add-to(X,X')
  2: add-to(X',X'')
  3: jump-if-empty(X', 7)
  4: dec(X')
  5: add-to(X'', Y_int)
  6: goto(3)
  7: halt
  ```

# The application
Our application consists of two components, 1) a simple compiler and 2) the simulator proper.

The compiler is optional; it converts a high-level psCRN program containing written functions and instructions, into a low-level set of reactions that constitute the psCRN protocol. To this end, for it to be recognized by the compiler, the psCRN program must contain three sections:

1. the label "polymers:", followed by a list of polymer species, separated by commas.
2. the set of high-level functions the program requires. The main program itself should be listed as a function, and labeled "main:". If input-detection is desired, add a function labeled "restart:". Other functions should be labeled with "function: func_name", where "func_name" is how that function is called in the program. All lines of the function body should be numbered as described in the paper. Use *i* for all non-restart functions; the restart function should be labeled with *R* instead.
3. the set of low-level instructions the program requires. Each should be labeled with "instruction: instr_name", where "inst_name" is how that instruction is called within the program. The lines after should be the set of reactions corresponding to that instruction, where each reaction takes up a separate line and follows the format "A + B -> C + D." If the instruction contains variables, enclose the variable component with braces.

The data folder contains more concrete working examples.
To run the compiler, run in the command line:

```python src/compile.py input_file```

To run the simulator, in the command line, run:

 ```python src/run.py input_file output_name```
 
 Where the input is simply a file containing the set of reactions in the format "A + B -> C + D", separated line by line. The compiler always generates a valid input for the simulator. Output name determines how the data from the simulation will be saved. The simulator can be customized through the run.ini file in the main directory.
