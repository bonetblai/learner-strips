%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

% Allocate max_predicates predicates of which the last max_static are
% static predicates. Some predicates may be left unused; these will be
% the last ones in each block of fluent or static predicates.

% max_precs and max_effects bound the number of preconditions and effects
% for each action.

#const max_predicates=5. % total number of predicates
#const max_static=2.     % these last number of predicates are static
#const max_precs=6.      % for each action
#const max_effects=6.    % for each action
#const opt_unequal=1.    % inequality constraints for action arguments
#const opt_synthesis=1.  % enable synthesis of model
#const opt_prec=2.       % 2 enforce positive prec, 1 prefer positive prec, else disabled

% Constants needed and set at call time:
%#const num_objects >= 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Marked and relevant nodes

#defined disabled/1.
#defined filename/2.
#defined partial/2.

marked(I,S) :- node(I,S), not disabled(I), not partial(I,File) : filename(I,File).

% To bypass instance I altogether       : assert disabled(I)
% To consider "some" nodes of instance I: remove disabled(I), assert partial(I,File), filename(I,File), and add marked(I,S) for target nodes

% relevant(I,S): marked or child of marked
relevant(I,S) :- marked(I,S).
relevant(I,S) :- marked(I,S1), tlabel(I,(S1,S),B).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicates and static predicates (p_static), objects and actions

pred(1..max_predicates).
p_static(max_predicates-max_static+1..max_predicates).
:- p_static(P), P = 1..max_predicates-max_static.
object(I,1..num_objects) :- instance(I).
action(1..N) :- num_actions(N), opt_synthesis = 1.
#defined num_actions/1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Atoms: ungrounded and grounded. Predicate's arity is either 0, 1, or 2.
%
% Ungrounded atoms are pairs (P,(V1 V2)) where P is predicate, and
% V1 and V2 are action arguments >= 1. For nullary predicates,
% V1 = V2 = 0, and V1 = V2 for unary.
%
% Grounded atoms are pairs (P,(O1,O2)) where O1 and O2 are objects.
% As before, if P is nullary O1 = O2 = 0, and O1 = O2 if P is unary.
%
% Actions on the other hand are of arity <= 3. Therefore, any ungrounded
% atom (P,(V1,V2)) in either precondition or effect has 1 <= V1, V2 <= 3,
% or V1 = V2 = 0 if action is nullary.
%
% We thus define tuples of arguments, tuples of objects, and function
% that maps pair (V1,V2) and obj tuple (grounded action argument)
% (O1,O2,O3) into pair (O4,O5).

% Tuples of variables for effects and preconditions

argtuple((0,0)).
argtuple((V1,V2)) :- V1 = 1..3, V2 = 1..3.

% Tuples of objects for atoms

objtuple(I,( 0, 0   ),2) :- instance(I).
objtuple(I,(O1,O2   ),2) :- object(I,O1), object(I,O2).
objtuple(I,(O1,O2,O3),3) :- object(I,O1), object(I,O2), object(I,O3).

% Map atoms arguments (argtuple) and action arguments (objtuples)
% into tuples of 2 objects (by selection)

% map(T,OO1,OO2): from input tuple OO1, produce output OO2 from refs in T

map((0,0),(O1,O2,O3),( 0, 0)) :- objtuple(I,(O1,O2,O3),3). % for nullary predicates
map((1,1),(O1,O2,O3),(O1,O1)) :- objtuple(I,(O1,O2,O3),3).
map((1,2),(O1,O2,O3),(O1,O2)) :- objtuple(I,(O1,O2,O3),3).
map((1,3),(O1,O2,O3),(O1,O3)) :- objtuple(I,(O1,O2,O3),3).
map((2,1),(O1,O2,O3),(O2,O1)) :- objtuple(I,(O1,O2,O3),3).
map((2,2),(O1,O2,O3),(O2,O2)) :- objtuple(I,(O1,O2,O3),3).
map((2,3),(O1,O2,O3),(O2,O3)) :- objtuple(I,(O1,O2,O3),3).
map((3,1),(O1,O2,O3),(O3,O1)) :- objtuple(I,(O1,O2,O3),3).
map((3,2),(O1,O2,O3),(O3,O2)) :- objtuple(I,(O1,O2,O3),3).
map((3,3),(O1,O2,O3),(O3,O3)) :- objtuple(I,(O1,O2,O3),3).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actions: atoms, variables, arities, and labels

% a_atom(J,A,M): action A uses atom M = (P,(V1,V2)) as prec (J=1) or effect (J=2)
%   a_atom(A,M): action A uses atom M = (P,(V1,V2))
%    a_var(A,V): action A uses var V
%  a_arity(A,N): action A has arity N where N in {0..3}

a_atom(1,A,M) :- prec(A,M,X).
a_atom(2,A,M) :- eff(A,M,X).
a_atom(A,M)   :- a_atom(J,A,M), J = 1..2.

a_var(A,V) :- a_atom(A,(_,(V,_))).
a_var(A,V) :- a_atom(A,(_,(_,V))).
:- action(A), a_var(A,V), V > 1, not a_var(A,V-1).
:- action(A), not a_var(A,V) : V = 0..3. % redundant but small

a_arity(A,N)     :- a_var(A,N), not a_var(A,N+1).
a_arity_geq(A,M) :- a_arity(A,N), M = 1..3, N = M..3.

% Action labels and consistent labels in each instance
% Actions are counted in python: num_actions.
% There are as many Actions as there are Labels for Actions.
% Same Label = Same Action
% label(L) There is the Label L.
% label(A,L) Action A has Label L. (bijective)
% labelname(I,B,L) in Instance I local Action B has Label L
% labelmatch(I,A,B) in Instance I local Action B is global Action A.

label(L) :- label(A,L).
label(L) :- labelname(I,A,L).
:- labelname(I,A,L1), labelname(I,A,L2), L1 != L2.

% If synthesis, assign schema labels from transition labels
label(A,L) :- labelname(I,A,L),
              not label(B,L) : action(B), B != A;
              not labelname(J,A,L2) : label(L2), L2 != L, J = 1..I-1;
              opt_synthesis = 1.
{ label(A,L) : label(L) } = 1 :- action(A), opt_synthesis = 1.
:- label(A,L1), label(A,L2), L1 != L2.

% Establish correspondence between schema and transition labels
labelmatch(I,A,B) :- label(A,L), labelname(I,B,L).
:- labelname(I,B,L), not label(A,L) : action(A).
:- labelname(I,B,L), not labelmatch(I,A,B) : label(A,L).

% inverse agregation
#defined inverse/3.
inverse_canidate(A,B) :- inverse(I,C,D), labelmatch(I,A,C), labelmatch(I,B,D).
inverse_compatible(I,A,B) :- inverse_canidate(A,B), inverse(I,C,D), labelmatch(I,A,C), labelmatch(I,B,D).
inverse_compatible(I,A,B) :- inverse_canidate(A,B), instance(I),  { labelmatch(I,A,C) } = 0.
inverse(A,B) :- inverse_canidate(A,B), inverse_compatible(I,A,B) : instance(I).
safe_inverse(A,B) :- inverse(A,B), inverse(B,A), A != B, {inverse(A,_)}=1, {inverse(B,_)}=1, {inverse(_,A)}=1, {inverse(_,B)}=1.

action_comp(I,A,B) :- action(A), labelmatch(I,A,B).
action_comp(I,A,0) :- action(A), labelmatch(I,A,_).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actions: preconditions and effects

% prec(A,(P,T),0): atom (P,T) is negative prec of action A
% prec(A,(P,T),1): atom (P,T) is positive prec of action A
% if action A is already synthesized, use frozen(A)
#defined frozen/1.
#defined act_inverse/1.
#defined partitioning/0.

{ prec(A,(P,T),0..1) : pred(P), argtuple(T) } max_precs :- action(A), not partitioning, not frozen(A), not act_inverse(A), opt_synthesis = 1, opt_prec != 2.
{ prec(A,(P,T),   1) : pred(P), argtuple(T) } max_precs :- action(A), not partitioning, not frozen(A), not act_inverse(A), opt_synthesis = 1, opt_prec  = 2.
:- prec(A,M,0), opt_prec = 2.        % rule out negative precs
:- prec(A,M,0), prec(A,M,1).         % rule out contradictory precs
:- prec(_,(P,(0,0)),_), p_static(P). % CHECK: not nullary prec on static predicates

% eff(A,(P,T),0): atom (P,T) is negative effect of action A
% eff(A,(P,T),1): atom (P,T) is positive effect of action A

% at least one effect (assumes there are no noop actions, i.e. self loops in graphs)
1 { eff(A,(P,T),0..1) : pred(P), argtuple(T), not p_static(P) } max_effects :- action(A), not partitioning, not frozen(A), not act_inverse(A), opt_synthesis = 1.
:- eff(A,M,0), eff(A,M,1).           % rule out contradictory effects

% No noop effects and no noop actions
:- prec(A,M,X), eff(A,M,X).

% Inverses impose heavy restrictions:
:-  eff(A,(P,T),0),                    not prec(A,(P,T),  1),           inverse(A,B).
:-  eff(A,(P,_),X),                       { eff(B,(P,_),1-X)} = 0,      inverse(A,B).
:-  eff(A,(P,T),X),                       { eff(B,(P,T),1-X)} = 0, safe_inverse(A,B).
:- prec(B,(P,T),X), not eff(B,(P,T),1-X), {prec(A,(P,_),  X)} = 0,      inverse(A,B).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicates

% p_arity(P,N): predicate P has arity N where N in {0..3}
%    p_used(P): predicate P is used in encoding

% CHECK: Method 1: Choose arity for predicates (Is it better?)
{ p_arity(P,0..2) } 1 :- pred(P), opt_synthesis = 1.
p_used(P)             :- p_arity(P,_).
:- p_arity(P,0), a_atom(_,(P,(V1, _))), V1 > 0.
:- p_arity(P,1), a_atom(_,(P,(V1,V2))), V1 != V2.
:- p_arity(P,K), a_atom(_,(P,( 0, 0))),  K > 0.

%*
% CHECK: Alternatively, from definition of actions infer arity of predicates
p_arity(P,2) :- a_atom(_,(P,(V1,V2))), V1 != V2.
p_arity(P,1) :- p_reflexive(P), not p_arity(P,2).
p_arity(P,0) :- a_atom(_,(P,(0,0))).
p_used(P)    :- p_arity(P,_).

p_reflexive(P) :- a_atom(_,(P,(V,V))), V > 0.
:- p_arity(P,0), p_arity(P,2).
:- p_arity(P,0), p_reflexive(P).
:- p_arity(P,1), p_arity(P,2).   % redundant
:- p_arity(P,0), p_static(P).    % CHECK: can nullary predicates be static?
*%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% States: valuations of ground atoms, either static or dynamic

% CHECK: Method 1: Choose true atoms and set the rest to false
{ val(I,(P,( 0, 0)),     1) }     :- pred(P),     p_static(P), p_arity(P,0), objtuple(I,( 0, 0),2),                                            opt_val = 1.
{ val(I,(P,( 0, 0)),S,   1) }     :- pred(P), not p_static(P), p_arity(P,0), objtuple(I,( 0, 0),2), relevant(I,S),                             opt_val = 1.
{ val(I,(P,( O, O)),     1) }     :- pred(P),     p_static(P), p_arity(P,1), objtuple(I,( O, O),2),                O != 0,                     opt_val = 1.
{ val(I,(P,( O, O)),S,   1) }     :- pred(P), not p_static(P), p_arity(P,1), objtuple(I,( O, O),2), relevant(I,S), O != 0,                     opt_val = 1.
{ val(I,(P,(O1,O2)),     1) }     :- pred(P),     p_static(P), p_arity(P,2), objtuple(I,(O1,O2),2),                O1 != 0,                    opt_val = 1.
{ val(I,(P,(O1,O2)),S,   1) }     :- pred(P), not p_static(P), p_arity(P,2), objtuple(I,(O1,O2),2), relevant(I,S), O1 != 0,                    opt_val = 1.
val(I,(P,(O1,O2)),  0)            :- pred(P),     p_static(P), p_arity(P,_), objtuple(I,(O1,O2),2),                not val(I,(P,(O1,O2)),  1), opt_val = 1.
val(I,(P,(O1,O2)),S,0)            :- pred(P), not p_static(P), p_arity(P,_), objtuple(I,(O1,O2),2), relevant(I,S), not val(I,(P,(O1,O2)),S,1), opt_val = 1.

% CHECK: Method 2: First, set false atoms of non-matching arity, then choose full valuation for the rest
val(I,(P,(O1,O2)),  0)            :- pred(P),     p_static(P), p_arity(P,0), objtuple(I,(O1,O2),2),                O1 != 0,                    opt_val = 2.
val(I,(P,(O1,O2)),S,0)            :- pred(P), not p_static(P), p_arity(P,0), objtuple(I,(O1,O2),2), relevant(I,S), O1 != 0,                    opt_val = 2.
val(I,(P,( 0, 0)),  0)            :- pred(P),     p_static(P), p_arity(P,N), objtuple(I,( 0, 0),2),                N = 1..2,                   opt_val = 2.
val(I,(P,( 0, 0)),S,0)            :- pred(P), not p_static(P), p_arity(P,N), objtuple(I,( 0, 0),2), relevant(I,S), N = 1..2,                   opt_val = 2.
val(I,(P,(O1,O2)),  0)            :- pred(P),     p_static(P), p_arity(P,1), objtuple(I,(O1,O2),2),                O1 != O2,                   opt_val = 2.
val(I,(P,(O1,O2)),S,0)            :- pred(P), not p_static(P), p_arity(P,1), objtuple(I,(O1,O2),2), relevant(I,S), O1 != O2,                   opt_val = 2.

{ val(I,(P,( 0, 0)),  0..1) } = 1 :- pred(P),     p_static(P), p_arity(P,0), objtuple(I,( 0, 0),2),                                            opt_val = 2.
{ val(I,(P,( 0, 0)),S,0..1) } = 1 :- pred(P), not p_static(P), p_arity(P,0), objtuple(I,( 0, 0),2), relevant(I,S),                             opt_val = 2.
{ val(I,(P,( O, O)),  0..1) } = 1 :- pred(P),     p_static(P), p_arity(P,1), objtuple(I,( O, O),2),                O != 0,                     opt_val = 2.
{ val(I,(P,( O, O)),S,0..1) } = 1 :- pred(P), not p_static(P), p_arity(P,1), objtuple(I,( O, O),2), relevant(I,S), O != 0,                     opt_val = 2.
{ val(I,(P,(O1,O2)),  0..1) } = 1 :- pred(P),     p_static(P), p_arity(P,2), objtuple(I,(O1,O2),2),                O1 != 0,                    opt_val = 2.
{ val(I,(P,(O1,O2)),S,0..1) } = 1 :- pred(P), not p_static(P), p_arity(P,2), objtuple(I,(O1,O2),2), relevant(I,S), O1 != 0,                    opt_val = 2.

% CHECK: Method 3: Choose full valuation for all atoms, then constraint allowed values
{ val(I,(P,OO),  0..1) } = 1      :- pred(P),     p_static(P), objtuple(I,OO,2),                                                               opt_val = 3.
{ val(I,(P,OO),S,0..1) } = 1      :- pred(P), not p_static(P), objtuple(I,OO,2), relevant(I,S),                                                opt_val = 3.

% Forced false atoms by predicate arity
p_arity_1_or_2(P) :- p_arity(P,N), N = 1..2,                                                                                                   opt_val = 3.

% If P hasn't arity 2, atoms P(O1,O2) where O1 != O2 are false
:- val(I,(P,(O1,O2)),S,1), not p_arity(P,2), O1 != O2,                                                                                         opt_val = 3.

% If P hasn't arity 1 or 2, atoms P(O,O) where O != 0 are false
:- val(I,(P,( O, O)),S,1), not p_arity_1_or_2(P), O != 0,                                                                                      opt_val = 3.

% If P hasn't arity 0, P(0,0) is false
:- val(I,(P,( 0, 0)),  1), not p_arity(P,0),                                                                                                   opt_val = 3.
:- val(I,(P,( 0, 0)),S,1), not p_arity(P,0),                                                                                                   opt_val = 3.

% Previous constraints imply all P(O1,O2) are false if P has no arity


% Different nodes are different states (for non-noisy nodes)
:- relevant(I,S1), relevant(I,S2), S1 < S2, val(I,(P,T),S2,V) : val(I,(P,T),S1,V); not noisy(S1), not noisy(S2).

% Directive to avoid warnings when noisy/1 isn't defined
#defined noisy/1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postulate and check action argument inequalities for grounding
% (enabled iff opt_unequal=1)

% unequal(A,(i,j)): no ground A(OO) for OO[i] == OO[j]

{ unequal(A,(1,2)) }                   :- a_arity_geq(A,2), not frozen(A), opt_unequal = 1, opt_synthesis = 1.
{ unequal(A,(1,3)); unequal(A,(2,3)) } :-     a_arity(A,3), not frozen(A), opt_unequal = 1, opt_synthesis = 1.

nappl(I,A,OO1)  :- a_arity_geq(A,2), objtuple(I,OO1,3), unequal(A,(1,2)), OO1 = (O1,O1,O3).
nappl(I,A,OO1)  :-     a_arity(A,3), objtuple(I,OO1,3), unequal(A,(1,3)), OO1 = (O1,O2,O1).
nappl(I,A,OO1)  :-     a_arity(A,3), objtuple(I,OO1,3), unequal(A,(2,3)), OO1 = (O1,O2,O2).

% Inequalities forced by schemas
diff(A,V1,V2) :- eff(A,(P,( V,V1)),X), eff(A,(P,( V,V2)),1-X), V1 < V2, opt_unequal = 1.
diff(A,V1,V2) :- eff(A,(P,(V1, V)),X), eff(A,(P,(V2, V)),1-X), V1 < V2, opt_unequal = 1.
diff(A,V1,V2) :- eff(A,(P,(V1,V1)),X), eff(A,(P,(V2,V2)),1-X), V1 < V2, opt_unequal = 1.
:- diff(A,V1,V2), not unequal(A,(V1,V2)), opt_unequal = 1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-applicability and applicability of grounded actions

% Ground action A(OO1) isn't applicable if it contains a false static precondition
nappl(I,A,OO1  ) :- val(I,(P,OO2),  V), map(T,OO1,OO2), prec(A,(P,T),1-V),     p_static(P).

% Ground action A(OO1) isn't applicable at S if it contains a false precondition in S
nappl(I,A,OO1,S) :- val(I,(P,OO2),S,V), map(T,OO1,OO2), prec(A,(P,T),1-V), not p_static(P).

% Ground action A(OO1) isn't applicable at S if A(OO1) isn't applicable (statically)
nappl(I,A,OO1,S) :- nappl(I,A,OO1), marked(I,S).

% Ground action A(OO1) is applicable at S if not non-applicable.
% For nullary actions, OO1 = (1,1,1), for unary OO1 = (O1,1,1), for binary
% OO1 = (O1,O2,1), and OO1 = (O1,O2,O3) for ternary actions.

appl(I,A,OO1,S) :- action(A), objtuple(I,OO1,3), marked(I,S), not nappl(I,A,OO1,S),
                   OO1 = ( 1, 1, 1).
appl(I,A,OO1,S) :- action(A), objtuple(I,OO1,3), marked(I,S), not nappl(I,A,OO1,S),
                   OO1 = (O1, 1, 1), a_arity_geq(A,1), O1 != 1.
appl(I,A,OO1,S) :- action(A), objtuple(I,OO1,3), marked(I,S), not nappl(I,A,OO1,S),
                   OO1 = (O1,O2, 1), a_arity_geq(A,2), O2 != 1.
appl(I,A,OO1,S) :- action(A), objtuple(I,OO1,3), marked(I,S), not nappl(I,A,OO1,S),
                   OO1 = (O1,O2,O3), a_arity_geq(A,3), O3 != 1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Successors

% Create missing unknown labels
tlabel(I,(S1,S2),0) :- {tlabel(I,(S1,S2),B) : labelmatch(I,_,B)} = 0, edge(I,(S1,S2)).

% If A(OO) is applicable at S1, choose a next node among those give by edges
{ next(I,A,OO,S1,S2) : tlabel(I,(S1,S2),B), action_comp(I,A,B) } = 1 :- appl(I,A,OO,S1), marked(I,S1).

% If action A maps S1 into S2, there is instance of A assigned to edge
:- marked(I,S1), tlabel(I,(S1,S2),B), { next(I,A,OO,S1,S2) : action_comp(I,A,B) } = 0.

%For inverses also the oposite direction must be present
:- next(I,A, _,S1,S2),  { next(I,B, _,S2,S1) } = 0,      inverse(A,B).
:- next(I,A,OO,S1,S2),  { next(I,B,OO,S2,S1) } = 0, safe_inverse(A,B).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check action effects on successor states (factored test)

% A maps S1 into S2 with O1 in V1, O2 in V2, and O3 in V3
pnext(I,0,A, 0,(O1,O2,O3),S1,S2) :- tlabel(I,(S1,S2),B), action(A), S2 >= 0, marked(I,S1), labelmatch(I,A,B), objtuple((O1,O2,O3),3). % for nullary predicates, object selection is irrelevant.
pnext(I,0,A, 0,(O1,O2,O3),S1,S2) :- next(I,A,(O1,O2,O3),S1,S2),     S2 >= 0. % for nullary predicates
pnext(I,1,A,O1,(O1,O2,O3),S1,S2) :- next(I,A,(O1,O2,O3),S1,S2),     S2 >= 0.
pnext(I,2,A,O2,(O1,O2,O3),S1,S2) :- next(I,A,(O1,O2,O3),S1,S2),     S2 >= 0.
pnext(I,3,A,O3,(O1,O2,O3),S1,S2) :- next(I,A,(O1,O2,O3),S1,S2),     S2 >= 0.

% If A goes from S1 to S2 with O1 in V1 and O2 in V2, the effects must hold in S2
:- eff(A,(P,(V1,V2)),X), pnext(I,V1,A,O1,OO,S1,S2), pnext(I,V2,A,O2,OO,S1,S2), val(I,(P,(O1,O2)),S2,1-X).

% For inverse all effects must be effective
:-   eff(A,(P,(V1,V2)),X),   pnext(V1,A,O1,OO,S1,S2), pnext(V2,A,O2,OO,S1,S2), val((P,(O1,O2)),S1,  X), inverse(A,B).

% an effect should be visible at one node at least
eff_used(A,(P,(V1,V2))) :- pnext(I,V1,A,O1,OO,S1,S2), pnext(I,V2,A,O2,OO,S1,S2), val(I,(P,(O1,O2)),S1,1-X), marked(I,S1), eff(A,(P,(V1,V2)),X).
:- eff(A,(P,(V1,V2)),X), not eff_used(A,(P,(V1,V2))).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertia

% There cannot be a non-caused change of a fluent
:- marked(I,S1), edge(I,(S1,S2)),    val(I,K,S1,V), val(I,K,S2,1-V), { caused(I,A,S1,S2,OO,K) } = 0.
:- marked(I,S1), next(I,A,OO,S1,S2), val(I,K,S1,V), val(I,K,S2,1-V), { caused(I,A,S1,S2,OO,K) } = 0.


% If A goes from S1 to S2 with O1 in V1 and O2 in V2, P(O1,O2) is caused if P(V1,V2) is effect
caused(I,A,S1,S2,OO,(P,(O1,O2))) :- affects(A,(P,(V1,V2))), pnext(I,V1,A,O1,OO,S1,S2), pnext(I,V2,A,O2,OO,S1,S2).
affects(A,(P,T)) :- eff(A,(P,T),X).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional constraints (some redundant)

% No more changes than the max of effects
:- marked(I,S1), edge(I,(S1,S2)), max_effects+1 #count{ M : val(I,M,S1,V), val(I,M,S2,1-V) }.

% If A has X outgoing edges in S1, it should be applicable X times
:- marked(I,S1), action(A), labelmatch(I,A,B), X = { tlabel(I,(S1,S2),B) }, { appl(I,A,OO,S1) } < X.
:- marked(I,S1), X = { tlabel(I,(S1,S2),C) : labelmatch(I,_,C) }, { appl(I,B,OO,S1) : action(B) } = Y, {tlabel(I,(S1,S2),0)} > Y - X.
:- { appl(I,A,OO,S1) : next(I,A,OO,S1,S2) } = 0 , action(A), tlabel(I,(S1,S2),B), labelmatch(I,A,B).
:- { appl(I,A,OO,S1) : next(I,A,OO,S1,S2), action(A) } = 0 , tlabel(I,(S1,S2),0).

% If there are X edges, there should be X applicable actions
:- action(A), labelmatch(I,A,B), X = { tlabel(I,(S1,S2),B) : marked(I,S1) }, { appl(I,A,OO,S) } < X.
:- instance(I), X = { tlabel(I,(S1,S2),A) : marked(I,S1) }, { appl(I,B,OO,S) } < X.

% Define root if no root node is specified
 min(I,S) :- marked(I,S), S <= S1 : marked(I,S1).
root(I,S) :- min(I,S), not root(I,N) : marked(I,N), N != S.

:- tlabel(I,(S1,S2),_), not edge(I,(S1,S2)).

% Indirect label distribution
#defined edge_num/3.
#defined edge_map/3.
:- edge_num(I,(S1,S2),J), next(I,A,_,S1,S2), edge_map(I,C,M), labelmatch(I,A,C), M > 1, N != M-1 : edge_num(I,(S3,S4),K), next(I,B,_,S3,S4), edge_map(I,D,N), labelmatch(I,B,D), J >= K.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display

#show action/1.
#show p_arity/2.
#show p_static/1.
#show p_used/1.

#show a_arity/2.
#show unequal/2.
#show label/2.
#show labelmatch/3.
#show labelname/3.
#show prec/3.
#show eff/3.

#show object/2.
#show val(I,K,S) : val(I,K,S,1).
#show val(I,K  ) : val(I,K,  1).

#show marked/2.
#show relevant/2.

