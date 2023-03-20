# learner-strips
Repository for learning (first-order) STRIPS models from graph-based representation of state spaces, either complete or partial.
The graph-based representation can be obtained from observable traces of (black box states) nodes and
action labels; for example, images and action labels. The assumptions being that different nodes correspond
to different planning states, and different planning states correspond to different nodes.

The learning of first-order STRIPS models is cast as a combinatorial optimization task where a simplest (first-order
STRIPS) model that fits a given training data is sought. The space of models is finitely bounded by setting values for
hypeparameters that bound the space (e.g. max arity of predicates/actions, max number of actions, max number of 
predicates, max number of objects, etc).

Two encodings for learning are presented. The original encoding using SAT and that is no longer developed, and
a somewhat less efficient but much more compact and flexible encoding using Answer Set Programming (ASP). 
The later one is being actively modified/extended.

This work is based on the following publications:
1. Blai Bonet and Hector Geffner.
*Learning First-Order Symbolic Representations for Planning from the Structure of the State Space.*
Proc. 24th European Conf. on Artificial Intelligence (ECAI). Santiago de Compostela, Spain. 2020. Pages 2322-2329.
2. Ivan D. Rodriguez, Blai Bonet, Javier Romero, and Hector Geffner.
*Learning First-Order Representations for Planning from Black-Box States: New Results.*
Proc. 18th Int. Conf. on Principles of Knowledge Representation and Reasoning (KR). Hanoi, Vietnam. 2021. Pages 539-548.

The repository depends on the following submodules:
* ``dfa-sat`` that include support for reading/writing for labeled directed graphs.
* ``sat-engine`` that provides C++ tools for writing SAT theories.

The repository is organized with the following folders:
* ``dfas`` that contains the input graphs from which a model is learned,
* ``sat`` that contains the SAT-based implementation of the learner, and
* ``asp`` that contains the ASP-based implementation of the learner.

## Inputs

The input for learning is a number of labeled directed graphs that encode (opaque) state spaces: each node in
the graph stands for a state and directed edges by state transitions associated to *ground actions*. The edges
are tagged with high-level action labels. These labels do not indicate arity nor object information. These labels
are not required, yet a graph without transition labels is considered to be a completely different input from 
the same graph with transition labels, and thus two completely different models can be expected.

A labeled directed graph is described as a text file in the following format:
1. The first line contains the reserved word ``dfa`` followed by two integers: the number of nodes in graph, and the second is expected to be ``-1``.
2. The second line specifies the set of labels; first, an integer that tells the number of labels, followed by a space-separated string labels.
3. The third line contains ``1 0``.
4. The following lines, one per each node in order. Each line begins with the number of outgoing edges at the node, then for each edge, the label of the edge followed by the index of the node pointed by the edge.

For example, the following file describes a directed graph with 22 nodes
and labels ``PUTDOWN``, ``PICK``, ``UNSTACK``, and ``STACK`` that corresponds
to the state space of a Blocksworld problem with 4 operators and 3 blocks.
Such a description is found in ```dfas/blocks4ops_3.dfa```.

```
dfa 22 -1
4 PUTDOWN PICK UNSTACK STACK
1 0
3 PICK 1 PICK 2 PICK 3
3 PUTDOWN 0 STACK 8 STACK 9
3 PUTDOWN 0 STACK 6 STACK 7
3 PUTDOWN 0 STACK 4 STACK 5
2 UNSTACK 3 PICK 15
2 UNSTACK 3 PICK 14
2 UNSTACK 2 PICK 13
2 UNSTACK 2 PICK 12
2 UNSTACK 1 PICK 11
2 UNSTACK 1 PICK 10
2 PUTDOWN 9 STACK 21
2 PUTDOWN 8 STACK 20
2 PUTDOWN 7 STACK 19
2 PUTDOWN 6 STACK 18
2 PUTDOWN 5 STACK 17
2 PUTDOWN 4 STACK 16
1 UNSTACK 15
1 UNSTACK 14
1 UNSTACK 13
1 UNSTACK 12
1 UNSTACK 11
1 UNSTACK 10
```

### Examples

The folder ```dfas``` contains several DFAs (labeled and directed graphs) from different domains.
These can be used directly by the SAT-based approach for synthesis. For using the ASP-based approach,
these must be processed to generate a ```.lp``` description of the corresponding DFA.


## SAT-based approach

The SAT-based approach mainly consists of a C++ program that generates a SAT theory from given inputs and
values for the hyperparameters. Although this program can be called directly and then the SAT solver, we
provide python scripts that implement the solution pipeline.

The source for the C++ program called ```strips``` is in ```sat/src```. To construct it, a simple ```make```
in the folder should suffice, yet some update to ```makefile``` may be required. The program requires the
```boost``` library and the submodules ```dfa-sat``` and ```sat-engine``` that can be installed with the
command ```git sumodule update --init```.

Directly executing ```src/strips``` gives a description of the different options and usages provided.
One such usage allows for the generation of ```.dot``` (graphical) description of an input DFA.

In the first paper described above, we solve a model learning problem by solving many (hundreds or thousands)
SAT theories using a computation cluster, and then selecting the solutions that verify over more bigger
instances and that are simpler in terms of number of actions and predicates, their arities, etc. For this,
we make use of the python script described below.


### Generation of graphs from DFAs

To generate a ```.pdf``` file that graphically depicts a DFA, execute:

```
./strips --dump-ts-dot --output <filename>.dot <filename>.dfa
dot -Tpdf -O <filename>.dot
```

The first command uses ```strips``` to generate a ```.dot``` file from the ```.dfa``` file, while
the second generates the ```.pdf``` file. For this, you need to have installed ```graphviz```.


### Python scripts

The folder ```sat/scripts``` contains the python script ```experiment.py``` that implements a complete
pipeline for solving and verifying. The script reads a *record* from a *benchmarks* file, construct
the SAT theory according to the values of the hyperparameters specified in the records, and, if successful,
verify the obtained STRIPS model on a number of input graphs, also specified in the record.
It is assumed that the SAT solver ```glucose``` is in the path.

Several benchmarks files are provided in the folder ```sat/scripts/benchmarks```. Each such file contains
one record per line. As an example, let us consider the file ```sat/scripts/benchmarks/benchmarks_grid4ops_1r.txt```
that contains the single record:

```
1 grid4ops_n4853  4 2 2 2 2  2 1 1  6 2 1  0 4 0 0 grid4ops_3x4   verify 8 grid4ops_4x4 grid4ops_5x6
```

The first integer (0 or 1) is for debugging purposes and indicates whether the record provides a solution
that yields a model that verifies. The second field is the name of the "experiment". The next integer
tells how many action schemas to synthesize together with the maximum arity for each one; in the example,
we instruct to do synthesis of 4 action schemas, each of maximum arity 2. The next integer tells how many
first-order predicates to include in the model together with their maximum arities; in the example, we
instruct to use 2 predicates of maximum arity 1. The following 3 integers indicate the total number
of different *first-order atoms*, and the number of *static* unary nd binary predicates; in the example,
the record tells to use 7 meta feature, 2 static unary predicates, and 1 binary predicate.

The record then contains one DFA accompanied by 4 integers: number of different ground atoms,
number of objects, and lower and upper bounds on the number of true atoms per state. Typically,
the all these integers are equal to 0 except the number of objects which must be a positive
integer. Zero values for the others tell the solver to construct a theory that is not constraint
on such dimensions. The name of the DFA appears after the 4 integers. At least one DFA should be
provided, but more than one is also supported.

Finally, the record specified the DFAs on which to *verify* an obtained model. In the example, if some
model is found, it should be verified in ```grid4ops_4x4``` and ```grid4ops_5x6``` using at most 8 objects.

In order to use the script, execute the following command from the ```sat/scripts``` folder:
```
python experiment.py --remove_dir --dfas_path ../../dfas --verbose 1 benchmarks/benchmarks_grid4ops_1r.txt 0
```
This tell the script to use record number 0 from the provided file, together with the options
```--remove_dir``` that remove any existing folder with the name ```grid4ops_n4853```,
```--dfas_path``` to indicate where to find the DFAs, and ```--verbose 1``` to obtain interesting
information. However, all the resulting files together with a log are stored in the folder ```grid4ops_n4853```.
Running ```experiment.py``` with the option ```--help``` provides information about additional flags.


## ASP-based approach

### Graphs

The DFAs that serve as inputs to the SAT-solver are preprocessed into ``.lp`` files that are understood by  ``clingo``.
The preprocessing is straightforward and consists of describing the graph using the atoms:
* ``node/1`` to describe the nodes in the graph with non-negative integer indices; e.g., ``node(0)``.
* ``labelname/2`` assigns symbolic names to indices for labels; e.g., ``labelname(`,"NEWTOWER")``.
* ``edge/1`` describes the (directed) edges in the graph where the argument is a pair of indices; e.g., ``edge((0,1))``.
* ``tlabel/2`` that assign labels to edges (e.g., ``tlabel((0,1),2)``). It is assumed that the pair is also declared as an edge and that the label index is also declared with ``labelname/2``.

Making a ``.lp`` file from a ``.dfa`` file is done with the script ``make_lp_from_dfa.py`` in the folder ``asp/scripts``. 
It takes two positional arguments, the first the path to the ``.dfa`` file and the second the path to the (new) ``.lp`` file.
For example, the following generates ``blocks3ops_5.lp`` in the current folder when executed in ``asp/scripts``:

```
$ python make_lp_from_dfa.py ../../dfas/blocks3ops_5.dfa .
```

### Solver

The ASP-based solver consists of a python driver program called ```incremental_solver.py``` located in
the ```scripts/``` folder together with ```.lp``` Clingo programs located in the ```clingo/``` folder.
There are different versions of the ```.lp``` programs (some of them currently under development) and
the current "official" and default release called ```kr21```.

The driver accepts a number of arguments that affects where the inputs/outputs are located,
options that affect the behavior of the driver, and options that affect the behavior of the solver.
As in the SAT-based approach, the "experiments" are described in ```.txt``` files, one experiment
per line, that are located in the ```benchmarks/``` folder. The options for the solver are:

```
usage: incremental_solver.py [-h] [--debug_level DEBUG_LEVEL]
                             [--inverse_actions] [--label_partitioning]
                             [--heuristics] [--no_invariants] [--no_optimize]
                             [--opt_val {1,2,3}]
                             [--incremental num-graphs max-depth]
                             [--version VERSION] [--add_lp ADD_LP]
                             [--add_flag ADD_FLAG] [--threads THREADS]
                             [--sat_prepro SAT_PREPRO] [--onlyver] [--skipver]
                             [--results RESULTS] [--graph_path SAMPLE_PATH]
                             [--partial_graph_path PARTIAL_SAMPLE_PATH]
                             [--solver_path SOLVER_PATH] [--remove_dir]
                             [--mem_bound MEM_BOUND] [--time_bound TIME_BOUND]
                             [--time_bound_ver TIME_BOUND_VER]
                             benchmarks record

positional arguments:
  benchmarks            Filename of file containing benchmarks
  record                Record index into benchmarks file

optional arguments:
  -h, --help            show this help message and exit
  --debug_level DEBUG_LEVEL
                        Set debug level (default=0)

preprocessing of input graphs:
  --inverse_actions     Identify inverse actions in input graph
  --label_partitioning  Identify L1/L2 label partitioning in input graph

options for solver:
  --heuristics          Apply heuristics when solving
  --no_invariants       Don't enforce invariants when solving
  --no_optimize         Don't do optimization when solving
  --opt_val {1,2,3}     Set method for choosing state valuation, only for 'mf'
                        version (default=1)
  --incremental num-graphs max-depth
                        Set options for incremental learning (default=None)
  --version VERSION     Set solver version (default=kr21)

options for Clingo:
  --add_lp ADD_LP       Add additional .lp file for solver (possibly multiple
                        times)
  --add_flag ADD_FLAG   Add additional flag for solver (possibly multiple
                        times)
  --threads THREADS     Set number of threads for solver (default=6)
  --sat_prepro SAT_PREPRO
                        Set SAT preprocessing option (default=0)

verification:
  --onlyver             Only do verification
  --skipver             Skip verification step

paths:
  --results RESULTS     Path to results folders (default='')
  --graph_path SAMPLE_PATH
                        Path to graphs (default='../graphs/full')
  --partial_graph_path PARTIAL_SAMPLE_PATH
                        Path to partial graphs (default='../graphs/partial')
  --solver_path SOLVER_PATH
                        Path to solver files (default='../clingo')
  --remove_dir          Discard existing files in results folder (if exists)

bounds:
  --mem_bound MEM_BOUND
                        Set memory bound for solver in MBs (default=None)
  --time_bound TIME_BOUND
                        Set time bound in seconds for synthesis (0 means no
                        bound, default=0)
  --time_bound_ver TIME_BOUND_VER
                        Set time bound in seconds for verification (0 means no
                        bound, default=0)
```

For the default behavior, the solver produces results that correspond (or are close)
to the ones reported in the KR2021 paper.

For example, to solve the blocks3ops problem and store the results in the folder
```results/```, it is enough to do the following:

```python incremental_solver.lp --remove_dir --results results benchmarks/blocks3ops_1r.txt 0```

This instruct the solver to run the experiment described in the first line (record number 0)
in the file ```benchmarks/blocks3ops_1r.txt```. The flag ```--remove_dir``` instruct the
solver to first erase the output folder for the experiment (if it exists).

The flag ```--version``` instruct the solver to read the ```.lp``` files from the corresponding
subfolder in ```clingo/```. The default behavior is to use the files in ```clingo/kr21/```.
Other options enable/disable the use of all the files in this folder (e.g., ```--no_optimize```,
```--no_invariants```, and ```--heuristics```). The ```--heuristics``` option does not yield
conclusive results, sometimes it improves performance and other times it degrades it.

