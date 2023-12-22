%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Optimization (ordered by priority)

% Prefer models with minimum sum of actions' cardinalities
#minimize { (1+N)**N@60,     A : a_arity(A,N) }.

% Prefer models with minimum sum of (non-static) predicates' cardinalities
#minimize { (1+N)**N@50,     P : p_arity(P,N), not p_static(P) }.

% Prefer models with minimum sum of (static) predicates' cardinalities
%#minimize { (1+N)**N@40,     P : p_arity(P,N), p_static(P) }. % (post paper)

% Prefer models with minimum number of preconditions (post paper)
#minimize {        1@30, A,M,V : prec(A,M,V) }.

% Prefer models with minimum number of effects (post paper)
#minimize {        1@20, A,M,V : eff(A,M,V) }.

% Prefer models with minimum number of true atoms per state
#minimize {        N@10,   I,N : bound_on_true_atoms_per_state(I,N) }.

% Prefer positive preconditions
#minimize {        1@1, prec(A,(P,T),0) : prec(A,(P,T),0) , opt_prec  = 1}.
