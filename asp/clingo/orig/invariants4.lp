%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constants

#const max_num_invariants=0.
#const max_true_atoms_per_state=0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose number and type of invariants to use

#heuristic invariant(N) : N=1..max_num_invariants. [N, true]
  invariant(1)   :-                 max_num_invariants != 0.
{ invariant(N) } :- invariant(N-1), N <= max_num_invariants.

% Schemas for invariants
% { inv_schema(N, P, 0) : not p_static(P);
%   inv_schema(N, P, 1)
% } :- invariant(N), pred(P), p_arity(P, 1), max_num_invariants != 0.
  { inv_schema(N, P, 1) } :- invariant(N), pred(P), p_arity(P, 1), max_num_invariants != 0.
1 { inv_schema(N, P, 2..3) : pred(P), p_arity(P, 2), not p_static(P) } :- invariant(N), max_num_invariants != 0.

% Different invariants
:- inv_schema(N1, P, X) : inv_schema(N2, P, X); invariant(N1), invariant(N2), N1 != N2.

% Show statements
#show num_invariants(N) : invariant(N), not invariant(N+1).
#show inv_schema/3.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Enforce invariants

{ %% val((P,(O2, O2)), S, 1) : object(O2), inv_schema(N, P, 0), not p_static(P);   % unary { P(O2) : O2 }
  val((P,  (O, O)),    1) :             inv_schema(N, P, 1),     p_static(P);   % (static) unary { P(O) }
  val((P,  (O, O)), S, 1) :             inv_schema(N, P, 1), not p_static(P);   % unary { P(O) }
  val((P, (O, O2)), S, 1) : object(O2), inv_schema(N, P, 2), not p_static(P);   % binary { P(O,O2) : O2 }
  val((P, (O2, O)), S, 1) : object(O2), inv_schema(N, P, 3), not p_static(P)    % binary { P(O2,O) : O2 }
} = 1 :- invariant(N), object(O), node(S), max_num_invariants != 0.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bound number of true atoms per state

%#heuristic bound_true_atoms_per_state(N) : N=1..max_true_atoms_per_state. [100-N, true]
{ bound_true_atoms_per_state(1..max_true_atoms_per_state) } = 1 :- max_true_atoms_per_state != 0.
{ val((P, OO), S, 1) : pred(P), not p_static(P), objtuple(OO, 2) } N :- node(S), bound_true_atoms_per_state(N), max_true_atoms_per_state != 0.

% Show statements
#show bound_true_atoms_per_state/1.

