%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization (ordered by priority)

% Prefer models with minimum sum of actions' cardinalities
#minimize { (1+N)**N@5, A : a_arity(A, N) }.

% Prefer models with minimum sum of (non-static) predicates' cardinalities
#minimize { (1+N)**N@4, P : p_arity(P, N), not p_static(P) }.

% Prefer models with minimum sum of (static) predicates' cardinalities
#minimize { (1+N)**N@3, P : p_arity(P, N), p_static(P) }.

% Prefer models with minimum number of true atoms per state
#minimize {        N@2, N : bound_on_true_atoms_per_state(N) }.

% Prefer positive preconditions
#minimize {        1@1, prec(A,(P,T),0) : prec(A,(P,T),0) , opt_prec  = 1}.
