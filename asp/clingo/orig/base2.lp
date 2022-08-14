%
% Constants
%

#const max_predicates=5.
#const static=2.
#const max_effects=6.    % for each action
#const max_precs=6.      % for each action


%
% Actions, predicates, static predicates and objects
%

action(L) :- tlabel(T,L).

pred(1..max_predicates).

p_static(max_predicates-static+1..max_predicates).

object(1..num_objects).


% Tuples of objects for effects and preconditions
%

argtuple((0,0)).
argtuple((V1,V2)) :- V1=1..3, V2=1..3.


%
% Generate effects of actions
%

% at least one effect? (yes if there are no 0-length loops)
% excludes static predicates
1 { eff(A,(P,T),0..1) : pred(P), argtuple(T), not p_static(P) } max_effects :- action(A), opt_synthesis=1.

:- eff(A,M,0), eff(A,M,1).


%
% Generate preconditions of actions
%

{ prec(A,(P,T),0..1) : pred(P), argtuple(T) } max_precs :- action(A), opt_synthesis=1.

:- prec(A,M,0), prec(A,M,1).
:- prec(_,(P,(0,0)),_), p_static(P).

% no effects to make preconditions persist
:- prec(A,M,V), eff(A,M,V).


%
% Variables used by A
%

cond(0,A,M,X) :- prec(A,M,X).
cond(1,A,M,X) :-  eff(A,M,X).

used(A,V) :- cond(_,A,(_,(V,_)),_).
used(A,V) :- cond(_,A,(_,(_,V)),_).
%
:- action(A), not used(A,V) : V=0..3. % redundant but small
%
% if used(A,V) then also used(A,V-1) (for V>1)
:- action(A), used(A,V), V>1, not used(A,V-1).


%
% Action and predicate arities
%

% a_arity(A,V) where V \in {0..3}
a_arity(A,V) :- used(A,V), not used(A,V+1). 

a_arity_geq(A,V) :- V=1..3, Ar = V..3, a_arity(A,Ar).

p_arity(P,2) :- cond(_,_,(P,(V1,V2)),_), V1 != V2.
%    
p_reflexive(P) :- cond(_,_,(P,(V,V)),_), V > 0.
p_arity(P,1) :- p_reflexive(P), not p_arity(P,2).
%
p_arity(P,0) :- cond(_,_,(P,(0,0)),_).
:- p_arity(P,0),   p_arity(P,2).
:- p_arity(P,0), p_reflexive(P).
%
p_used(P) :- p_arity(P,_).


%
% Tuples of objects for atoms
%

objtuple(( 0, 0   ),2).                            % for nullary predicates
objtuple((O1,O2   ),2) :- object(O1), object(O2).
objtuple((O1,O2,O3),3) :- object(O1), object(O2), object(O3).

%
% Possible valuations of ground atoms
%

{ val((P,OO),  0..1) } = 1 :- pred(P),     p_static(P), objtuple(OO,2).
{ val((P,OO),S,0..1) } = 1 :- pred(P), not p_static(P), objtuple(OO,2), node(S).


%
% Map 2-argtuples and 3-objtuples to 2-objtuples
%

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


%
% Variable inequalities
%

{ unequal(A,(1,2)) }                   :- a_arity_geq(A,2), unequal=1.
{ unequal(A,(1,3)); unequal(A,(2,3)) } :-     a_arity(A,3), unequal=1.

nappl(A,OO1) :- a_arity_geq(A,2), objtuple(OO1,3), unequal(A,(1,2)), OO1 = (O1,O1,O3).
nappl(A,OO1) :-     a_arity(A,3), objtuple(OO1,3), unequal(A,(1,3)), OO1 = (O1,O2,O1).
nappl(A,OO1) :-     a_arity(A,3), objtuple(OO1,3), unequal(A,(2,3)), OO1 = (O1,O2,O2).

swap(A, V1, V2) :- eff(A,(P,(V,V1)),X), eff(A,(P,(V,V2)),1-X), V1 < V2, unequal=1.
swap(A, V1, V2) :- eff(A,(P,(V1,V)),X), eff(A,(P,(V2,V)),1-X), V1 < V2, unequal=1.
swap(A, V1, V2) :- eff(A,(P,(V1,V1)),X), eff(A,(P,(V2,V2)),1-X), V1 < V2, unequal=1.

:- swap(A, V1, V2), not unequal(A,(V1,V2)), unequal=1.

%
% Check preconditions
%

% static
nappl(A,OO1  ) :- val((P,OO2),  V), map(T,OO1,OO2), prec(A,(P,T),1-V),     p_static(P).
nappl(A,OO1,S) :- nappl(A,OO1), node(S).
% dynamic
nappl(A,OO1,S) :- val((P,OO2),S,V), map(T,OO1,OO2), prec(A,(P,T),1-V), not p_static(P).
%
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = ( 1, 1, 1).
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = (O1, 1, 1), a_arity_geq(A,1), O1 != 1.
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = (O1,O2, 1), a_arity_geq(A,2), O2 != 1.
appl(A,OO1,S) :- action(A), objtuple(OO1,3), node(S), not nappl(A,OO1,S),
                 OO1 = (O1,O2,O3), a_arity_geq(A,3), O3 != 1.

% if A(OO) is applicable at S1, choose a next node
{ next(A,OO,S1,S2) : tlabel((S1,S2),A) } = 1 :- appl(A,OO,S1).

% if action A goes from S1 to S2,
% then there must be an instance of A assigned to that edge
:- tlabel((S1,S2),A), { next(A,OO,S1,S2) } != 1.


%
% Check effects
%

% projected next on the arguments of the action
pnext(0,A, 0,S1,S2) :-        tlabel((S1,S2),A), S2 >= 0. % for nullary predicates
pnext(1,A,O1,S1,S2) :- next(A,(O1,O2,O3),S1,S2), S2 >= 0.
pnext(2,A,O2,S1,S2) :- next(A,(O1,O2,O3),S1,S2), S2 >= 0.
pnext(3,A,O3,S1,S2) :- next(A,(O1,O2,O3),S1,S2), S2 >= 0.

% if A goes from S1 to S2 with O1 in V1 and O2 in V2,
% then the effects must hold in S2
:- eff(A,(P,(V1,V2)),X), pnext(V1,A,O1,S1,S2), pnext(V2,A,O2,S1,S2),
   val(  (P,(O1,O2)),S2,1-X).

% inertia

% there cannot be a non-caused change of a fluent
:- edge((S1,S2)), val(K,S1,V), val(K,S2,1-V), not caused(S1,S2,K).

% projection on eff/2
eff(A,(P,T)) :- eff(A,(P,T),V).

% If A goes from S1 to S2 with O1 in V1 and O2 in V2,
% then if P(V1,V2) is an effect P(O1,O2) is caused
caused(S1,S2,(P,(O1,O2))) :- eff(A,(P,(V1,V2))), pnext(V1,A,O1,S1,S2), pnext(V2,A,O2,S1,S2).


%
% Check that different nodes are different states
%

:- node(S1), node(S2), S1 < S2, val((P,T), S2, V) : val((P,T), S1, V);
   not noisy(S1), not noisy(S2).

#defined noisy/1.

%
% Display
%
#show a_arity/2.
#show p_arity/2.
#show p_static/1.
#show p_used/1.

#show eff/3.
#show prec/3.
#show labelname/2.
#show val(K,S) : val(K,S,1).
#show val(K  ) : val(K,  1).
#show object/1.
#show unequal/2.
%
% Heuristics
%
%*
#heuristic    a_arity(A,Ar). [1,level]
#heuristic    p_arity(A,Ar). [1,level]
#heuristic       eff(A,M,V). [1,level]
#heuristic      prec(A,M,V). [2,level]
#heuristic       val(K,S,V). [1,level]
#heuristic next(A,OO,S1,S2). [1,level]
*%


%
% Additional constraints
%

% no more changes than the max of effects
:- edge((S1,S2)), max_effects+1 #count{K : val(K,S1,V), val(K,S2,1-V) }.

% if A has X outgoing edges in S1, then it should be applicable at least X times
:- node(S1), action(A), X = { tlabel((S1,S2),A) }, { appl(A,OO,S1) } != X.

% if there are X edges, then there should be X applicable actions
:- X = { tlabel(E,L) }, { appl(A,B,C) } != X.

% if P has not arity 2 then the atoms P(O1,O2) where O1 != O2 are false
:- val((P,(O1,O2)),S,1), not p_arity(P,2), O1 != O2.

% if P has not arity 1 or 2 then the atoms P(O,O) where O != 0 are false
p_arity_1_or_2(P) :- p_arity(P,Ar), Ar=1..2.
:- val((P,(O,O)),S,1), not p_arity_1_or_2(P), O != 0.

% if P has not arity 0 then P(0,0) is false
:- val((P,(0,0)),1),   not p_arity(P,0).
:- val((P,(0,0)),S,1), not p_arity(P,0).    

% By the previous constraints, 
% if P has no arity then all P(O1,O2) are false

% define root(0) by default if no root node specified
root(0) :- not root(N) : node(N), N != 0.

