%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Previsouly in base

%*
#heuristic    a_arity(A,Ar). [1,level]
#heuristic    p_arity(A,Ar). [1,level]
#heuristic       eff(A,M,V). [1,level]
#heuristic      prec(A,M,V). [2,level]
#heuristic       val(K,S,V). [1,level]
#heuristic next(A,OO,S1,S2). [1,level]
*%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Previsouly in constraints_javier

#heuristic  eff(A,M,1-V) : not prec(A,M,V), opt_req_prec = 2. [1,false]
#heuristic prec(A,M,1-V) :      eff(A,M,V), opt_req_prec = 2. [1, true]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Previsouly in invariants4a

#heuristic invariant(N) : N=1..max_num_invariants. [N, true]
%#heuristic bound_on_true_atoms_per_state(N) : N=1..max_true_atoms_per_state. [100-N, true]

