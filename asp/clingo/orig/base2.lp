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

% Constants needed and set at call time:
%#const num_objects >= 1


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicates and static predicates (p_static), objects and actions

pred(1..max_predicates).
p_static(max_predicates-max_static+1..max_predicates).
:- p_static(P), P = 1..max_predicates-max_static.
object(1..num_objects).
action(1..N) :- num_actions(N), opt_synthesis = 1.
#defined num_actions/1.

%%%action(L) :- tlabel(T,L). % BLAI: old


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

objtuple(( 0, 0   ),2).% for nullary predicates
objtuple((O1,O2   ),2) :- object(O1), object(O2).
objtuple((O1,O2,O3),3) :- object(O1), object(O2), object(O3).

% Map atoms arguments (argtuple) and action arguments (objtuples)
% into tuples of 2 objects (by selection)

% map(T,OO1,OO2): from input tuple OO1, produce output OO2 from refs in T

map((0,0),(O1,O2,O3),( 0, 0)) :- objtuple((O1,O2,O3),3). % for nullary predicates
map((1,1),(O1,O2,O3),(O1,O1)) :- objtuple((O1,O2,O3),3).
map((1,2),(O1,O2,O3),(O1,O2)) :- objtuple((O1,O2,O3),3).
map((1,3),(O1,O2,O3),(O1,O3)) :- objtuple((O1,O2,O3),3).
map((2,1),(O1,O2,O3),(O2,O1)) :- objtuple((O1,O2,O3),3).
map((2,2),(O1,O2,O3),(O2,O2)) :- objtuple((O1,O2,O3),3).
map((2,3),(O1,O2,O3),(O2,O3)) :- objtuple((O1,O2,O3),3).
map((3,1),(O1,O2,O3),(O3,O1)) :- objtuple((O1,O2,O3),3).
map((3,2),(O1,O2,O3),(O3,O2)) :- objtuple((O1,O2,O3),3).
map((3,3),(O1,O2,O3),(O3,O3)) :- objtuple((O1,O2,O3),3).


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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Actions: preconditions and effects

% prec(A,(P,T),0): atom (P,T) is negative prec of action A
% prec(A,(P,T),1): atom (P,T) is positive prec of action A

{ prec(A,(P,T),0..1) : pred(P), argtuple(T) } max_precs :- action(A), opt_synthesis = 1.
:- prec(A,M,0), prec(A,M,1).         % rule out contradictory precs
:- prec(_,(P,(0,0)),_), p_static(P). % CHECK: not nullary prec on static predicates

% eff(A,(P,T),0): atom (P,T) is negative effect of action A
% eff(A,(P,T),1): atom (P,T) is positive effect of action A

% at least one effect (assumes there are no noop actions, i.e. self loops in graphs)
1 { eff(A,(P,T),0..1) : pred(P), argtuple(T), not p_static(P) } max_effects :- action(A), opt_synthesis = 1.
:- eff(A,M,0), eff(A,M,1).           % rule out contradictory effects

% No noop effects and no noop actions
:- prec(A,M,X), eff(A,M,X).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Predicates

% p_arity(P,N): predicate P has arity N where N in {0..3}
%    p_used(P): predicate P is used in encoding

p_arity(P,2) :- a_atom(_,(P,(V1,V2))), V1 != V2.
p_arity(P,1) :- p_reflexive(P), not p_arity(P,2).
p_arity(P,0) :- a_atom(_,(P,(0,0))).
p_used(P)    :- p_arity(P,_).

p_reflexive(P) :- a_atom(_,(P,(V,V))), V > 0.
:- p_arity(P,0), p_arity(P,2).
:- p_arity(P,0), p_reflexive(P).
:- p_arity(P,1), p_arity(P,2).   % redundant
:- p_arity(P,0), p_static(P).    % CHECK: can nullary predicates be static?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% States: valuations of ground atoms, either static or dynamic

{ val((P,OO),  0..1) } = 1 :- pred(P),     p_static(P), objtuple(OO,2).
{ val((P,OO),S,0..1) } = 1 :- pred(P), not p_static(P), objtuple(OO,2), node(S).

% Forced false atoms by predicate arity
p_arity_1_or_2(P) :- p_arity(P,N), N = 1..2.

% If P hasn't arity 2, atoms P(O1,O2) where O1 != O2 are false
:- val((P,(O1,O2)),S,1), not p_arity(P,2), O1 != O2.

% If P hasn't arity 1 or 2, atoms P(O,O) where O != 0 are false
:- val((P,(O,O)),S,1), not p_arity_1_or_2(P), O != 0.

% If P hasn't arity 0, P(0,0) is false
:- val((P,(0,0)),1),   not p_arity(P,0).
:- val((P,(0,0)),S,1), not p_arity(P,0).

% Previous constraints imply all P(O1,O2) are false if P has no arity


% Different nodes are different states (for non-noisy nodes)
:- node(S1), node(S2), S1 < S2, val((P,T),S2,V) : val((P,T),S1,V); not noisy(S1), not noisy(S2).

% Directive to avoid warnings when noisy/1 isn't defined
#defined noisy/1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postulate and check action argument inequalities for grounding
% (enabled iff opt_unequal=1)

% unequal(A,(i,j)): no ground A(OO) for OO[i] == OO[j]

{ unequal(A,(1,2)) }                   :- a_arity_geq(A,2), opt_unequal = 1.
{ unequal(A,(1,3)); unequal(A,(2,3)) } :-     a_arity(A,3), opt_unequal = 1.

nappl(A,OO1) :- a_arity_geq(A,2), objtuple(OO1,3), unequal(A,(1,2)), OO1 = (O1,O1,O3).
nappl(A,OO1) :-     a_arity(A,3), objtuple(OO1,3), unequal(A,(1,3)), OO1 = (O1,O2,O1).
nappl(A,OO1) :-     a_arity(A,3), objtuple(OO1,3), unequal(A,(2,3)), OO1 = (O1,O2,O2).

% Inequalities forced by schemas
diff(A,V1,V2) :- eff(A,(P,( V,V1)),X), eff(A,(P,( V,V2)),1-X), V1 < V2, opt_unequal = 1.
diff(A,V1,V2) :- eff(A,(P,(V1, V)),X), eff(A,(P,(V2, V)),1-X), V1 < V2, opt_unequal = 1.
diff(A,V1,V2) :- eff(A,(P,(V1,V1)),X), eff(A,(P,(V2,V2)),1-X), V1 < V2, opt_unequal = 1.
:- diff(A,V1,V2), not unequal(A,(V1,V2)), opt_unequal = 1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-applicability and applicability of grounded actions

% Ground action A(OO1) isn't applicable if it contains a false static precondition
nappl(A,OO1  ) :- val((P,OO2),  V), map(T,OO1,OO2), prec(A,(P,T),1-V),     p_static(P).

% Ground action A(OO1) isn't applicable at S if it contains a false precondition in S
nappl(A,OO1,S) :- val((P,OO2),S,V), map(T,OO1,OO2), prec(A,(P,T),1-V), not p_static(P).

% Ground action A(OO1) isn't applicable at S if A(OO1) isn't applicable (statically)
nappl(A,OO1,S) :- nappl(A,OO1), node(S).

% Ground action A(OO1) is applicable at S if not non-applicable.
% For nullary actions, OO1 = (1,1,1), for unary OO1 = (O1,1,1), for binary
% OO1 = (O1,O2,1), and OO1 = (O1,O2,O3) for ternary actions.

appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = ( 1, 1, 1).
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = (O1, 1, 1), a_arity_geq(A,1), O1 != 1.
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = (O1,O2, 1), a_arity_geq(A,2), O2 != 1.
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = (O1,O2,O3), a_arity_geq(A,3), O3 != 1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Successors

% If A(OO) is applicable at S1, choose a next node among those give by edges
{ next(A,OO,S1,S2) : tlabel((S1,S2),A) } = 1 :- appl(A,OO,S1).

% If action A maps S1 into S2, there is instance of A assigned to edge
:- tlabel((S1,S2),A), { next(A,OO,S1,S2) } != 1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check action effects on successor states (factored test)

% A maps S1 into S2 with O1 in V1, O2 in V2, and O3 in V3
pnext(0,A, 0,S1,S2) :-        tlabel((S1,S2),A), S2 >= 0. % for nullary predicates
pnext(1,A,O1,S1,S2) :- next(A,(O1,O2,O3),S1,S2), S2 >= 0.
pnext(2,A,O2,S1,S2) :- next(A,(O1,O2,O3),S1,S2), S2 >= 0.
pnext(3,A,O3,S1,S2) :- next(A,(O1,O2,O3),S1,S2), S2 >= 0.

% If A goes from S1 to S2 with O1 in V1 and O2 in V2, the effects must hold in S2
:- eff(A,(P,(V1,V2)),X), pnext(V1,A,O1,S1,S2), pnext(V2,A,O2,S1,S2), val((P,(O1,O2)),S2,1-X).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inertia

% There cannot be a non-caused change of a fluent
:- edge((S1,S2)), val(K,S1,V), val(K,S2,1-V), not caused(S1,S2,K).

% If A goes from S1 to S2 with O1 in V1 and O2 in V2, P(O1,O2) is caused if P(V1,V2) is effect
caused(S1,S2,(P,(O1,O2))) :- affects(A,(P,(V1,V2))), pnext(V1,A,O1,S1,S2), pnext(V2,A,O2,S1,S2).
affects(A,(P,T)) :- eff(A,(P,T),X).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Additional constraints (some redundant)

% No more changes than the max of effects
:- edge((S1,S2)), max_effects+1 #count{ M : val(M,S1,V), val(M,S2,1-V) }.

% If A has X outgoing edges in S1, it should be applicable X times
:- node(S1), action(A), X = { tlabel((S1,S2),A) }, { appl(A,OO,S1) } != X.

% If there are X edges, there should be X applicable actions
:- action(A), X = { tlabel(E,A) }, { appl(A,OO,S) } != X.
:- X = { tlabel(E,A) }, { appl(B,OO,S) } != X.

% Define root(0) by default if no root node specified
root(0) :- not root(N) : node(N), N != 0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display

#show action/1.
#show p_arity/2.
#show p_static/1.
#show p_used/1.

#show a_arity/2.
#show unequal/2.
#show labelname/2.
#show prec/3.
#show eff/3.

#show object/1.
#show val(K,S) : val(K,S,1).
#show val(K  ) : val(K,  1).

