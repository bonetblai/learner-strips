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
1. First line contains the reserved word ``dfa`` followed by two integers: the number of nodes in graph, and the second is expected to be ``-1``.
2. Second line specifies the set of labels; first, an integer that tells the number of labels, followed by a space-separated string labels.
3. Third lines contains ``1 0``.
4. Following lines, one per each node in order. Each line begins with the number of outgoing edges at the node, then for each edge, the label of the edge followed by the index of the node pointed by the edge.

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


## SAT-Based Approach

The SAT-based approach mainly consists of a C++ program that generates a SAT theory from given inputs and
values for the hyperparameters. Although this program can be called directly and then the SAT solver, we
provide python scripts that implement the solution pipeline.

The source for the C++ program called ```strips``` is in ```sat/src```. To construct it, a simple ```make```
in the folder should suffice, yet some update to ```makefile``` may be required. The program requires the
```boost``` library and the submodules ```dfa-sat``` and ```sat-engine``` that can be installed with the
command ```git sumodule update --init```.

Directly executing ```src/strips``` gives a description of the different options and usages provided.
One such usage allows for the generation of ```.dot``` (graphical) description of an input DFA.


### Generation of graphs for DFAs

To generate a ```.dot``` file that graphically depicts a DFA as a ```.pdf``` file, execute:

```
./strips --dump-ts-dot --output <filename>.dot <filename>.dfa
dot -Tpdf -O <filename>.dot
```

### Scripts



The learning task is expressed as a SAT-

A SAT theory whose models (if any) are in correspondence with the STRIPS models compatible with the input graph for a given set of *hyperparameters*
is obtained with the executable ``src/strips``. The argumens are:

1. The number of action schemas ``<num-action-schemas>``.
2. The max arity of each action schema, in order: ``<max-arity-first-action> ... <max-arity-last-action>``.
3. The number of atom schemas: ``<num-atom-schemas>``.
4. The max arity of each atom schema, in order: ``<max-arity-first-atom> ... <max-arity-last-atom>``.
5. The number of ''first-order atoms'' (called meta-features) in the STRIPS model: ``<num-meta-features>``.
6. The number of *static* unary and binary predicates: ``<num-static-unary-predicates> <num-static-binary-predicates>``.
7. For each input graph, the number of ''features'', number of objects, lower and upper bounds of number of feature per state, and prefix for graph name: ``[<num-features> <num-objects> <LB-on-sum-features-per-state> <UB-on-sum-features-per-state> <prefix>]``.

Typically, ``<num-features>``, ``<LB-on-sum-features-per-state>``, and ``<UB-on-sum-features-per-state>`` are zero.
The executable ``src/strips`` allows for a number of options, use ``--help`` to see a full description.

## Scripts to Run Experiments

Scripts in folder ``aws/`` can be used to run experiments without the need to run ``src/strips`` directly.
