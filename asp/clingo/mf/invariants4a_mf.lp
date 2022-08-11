%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

% Max number of invariants to synthetize
#const max_num_invariants=1.

% Max true atoms per state.
#const max_true_atoms_per_state=0.

% Must have binary predicates in encoding.
#const opt_must_have_binary_predicates=0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose number and type of invariants to use

  invariant(1)   :-                 max_num_invariants > 0.
{ invariant(N) } :- invariant(N-1), N <= max_num_invariants.

% Schemas for invariants
%
% Unary predicates can only enter as type-1 predicates.
% Binary predicates can enter as type-2 or type-3 predicates.
% If there are binary predicates, each invariant must contain at least one such predicate.

  { inv_schema(N,P,1)                                             } :- invariant(N), pred(P), p_arity(P,1), max_num_invariants > 0.
1 { inv_schema(N,P,2..3) : pred(P), p_arity(P,2), not p_static(P) } :- invariant(N), some_binary_predicate, max_num_invariants > 0.

inv_non_empty(N) :- invariant(N), pred(P), inv_schema(N,P,1..3).
some_binary_predicate :- pred(P), p_arity(P,2), not p_static(P).
:- not some_binary_predicate, opt_must_have_binary_predicates != 0.

% Invariants must be different
:- inv_schema(N1,P,X) : inv_schema(N2,P,X); invariant(N1), invariant(N2), N1 != N2.

% Each non-static binary predicate must appear in some type-1 invariant
inv_used_pred(P) :- pred(P), p_arity(P,2), not p_static(P), invariant(N), inv_schema(N,P,2..3).
:- pred(P), p_arity(P,2), not p_static(P), not inv_used_pred(P), max_num_invariants > 0.

% Display
#show num_invariants(N) : invariant(N), not invariant(N+1).
#show inv_schema/3.
#show some_binary_predicate/0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enforce invariants

% An invariant defines a set of true atoms for a given object O and node S.
% The invariant is valid only if the set contains at most one atom.

% For inv_schema(N,P,1) and P static, set contains P(O) if P(O) is true
% For inv_schema(N,P,1) and P non-static, set contains P(O) if P(O) is true in S
% For inv_schema(N,P,2) and P non-static, set contains atoms in { P(O,O2) : O2 } true in S
% For inv_schema(N,P,3) and P non-static, set contains atoms in { P(O2,O) : O2 } true in S

{ %%val(I,(P,(O2,O2)),S,1) : object(I,O2), inv_schema(N,P,0), not p_static(P);   % unary { P(O2) : O2 }
  val(I,(P, (O,O)),  1) :               inv_schema(N,P,1),     p_static(P);   % (static) unary { P(O) }
  val(I,(P, (O,O)),S,1) :               inv_schema(N,P,1), not p_static(P);   % unary { P(O) }
  val(I,(P,(O,O2)),S,1) : object(I,O2), inv_schema(N,P,2), not p_static(P);   % binary { P(O,O2) : O2 }
  val(I,(P,(O2,O)),S,1) : object(I,O2), inv_schema(N,P,3), not p_static(P)    % binary { P(O2,O) : O2 }
} 1 :- invariant(N), inv_non_empty(N), object(I,O), relevant(I,S), max_num_invariants > 0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bound number of true atoms per state

{ bound_on_true_atoms_per_state(I,1..max_true_atoms_per_state) } = 1 :- instance(I), max_true_atoms_per_state > 0.
{ val(I,(P,OO),S,1) : pred(P), not p_static(P), objtuple(I,OO,2) } N :- relevant(I,S), bound_on_true_atoms_per_state(I,N), max_true_atoms_per_state != 0.

% Display
#show bound_on_true_atoms_per_state/2.

