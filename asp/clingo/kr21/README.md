
* ```base2.lp``` contains the basic solver which it is enough to find a
complete model that solves an input graph.

* ```optimize.lp``` instruct the solver to find the simplest model determined
by the optimization criteria that solves an input graph.

* ```heuristics.lp``` define "heuristics" to improve performance. However,
these heuristics do not always translate into better performace. It is still
an open question how to improve the heuristics. By default, the solver does
not use the heuristics. The flag ```--heuristics``` enforces them.

* The other programs add additional constraints, some of them redundant, with
the idea of improving performace. Automatically learned invariants that are
then enforced is done with ```invariants4a.lp``` while the ```invariants4.lp```
is not currently used.
To disable the use of invariants, use option ```--no_invariants```.

* Finally, ```constraints_blai.lp``` and ```constraints_javier.lp``` define
additional constraints and symmetry-reduction rules with the idea of 
improving performace.

