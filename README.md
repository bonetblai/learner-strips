# learner-strips
Repo for learning (first-order) STRIPS models from graph-based representation of state spaces, either complete or partial.
The graph-based representation can be obtained from observable traces of (black box states) noes and
action labels; for example, images and action labels. The assumptions being that different nodes correspond
to different planning states, and different planning states correspond to different nodes.

The learning of first-order STRIPS models is cast as a combinatorial optimization task where a simplest (first-order
STRIPS) model that fits a given training data is sought. The space of models is finitely bounded by setting values for
hypeparameters that bound the space (e.g. max arity of predicates/actions, max actions, max predicates, max objects, etc).

Two encodings for learning are presented. The original encoding using SAT and that is no longer developed, and
a somewhat less efficient but much more compact and flexible encoding using ASP. The later one is being actively
modified/extended.

This work is based on the following publications:
1. Blai Bonet and Hector Geffner.
*Learning First-Order Symbolic Representations for Planning from the Structure of the State Space.*
Proc. 24th European Conf. on Artificial Intelligence (ECAI). Santiago de Compostela, Spain. 2020. Pages 2322-2329.
2. Ivan D. Rodriguez, Blai Bonet, Javier Romero, and Hector Geffner.
*Learning First-Order Representations for Planning from Black-Box States: New Results.*
Proc. 18th Int. Conf. on Principles of Knowledge Representation and Reasoning (KR). Hanoi, Vietnam. 2021. Pages 539-548.

The repo depends on the following submodules:
1. ``dfa-sat`` that include support for reading/writing for labeled directed graphs.
2. ``sat-engine`` that provides C++ tools for writing SAT theories.

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


