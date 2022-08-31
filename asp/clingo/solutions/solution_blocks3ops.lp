% Predicates
pred(1).      % ontable(?b)
p_arity(1,1).
pred(2).      % clear(?b,?b) / below(?b1,?b2) "?b2 on ?b1"
p_arity(2,2).

% NEWTOWER(?x,?y): 2p + 3e / Unstack ?x from ?y
action(1).
a_arity(1,2).
prec(1,(1,(1,1)),0) :- opt_req_prec > 0. % -ontable(?x)
prec(1,(2,(1,1)),1). % clear(?x)
prec(1,(2,(2,2)),0) :- opt_req_prec > 0. % -clear(?y)
prec(1,(2,(2,1)),1). % below(?y,?x)
eff(1,(1,(1,1)),1).  % ontable(?x)
eff(1,(2,(2,1)),0).  % -below(?y,?x)
eff(1,(2,(2,2)),1).  % clear(?y)
unequal(1,(1,2)).

% STACK(?x,?y): 3p + 3e / Stack ?x on ?y
action(2).
a_arity(2,2).
prec(2,(1,(1,1)),1). % ontable(?x)
prec(2,(2,(1,1)),1). % clear(?x)
prec(2,(2,(2,2)),1). % clear(?y)
prec(2,(2,(2,1)),0) :- opt_req_prec > 0. % -below(?y,?x)
eff(2,(1,(1,1)),0).  % -ontable(?x)
eff(2,(2,(2,1)),1).  % below(?y,?x)
eff(2,(2,(2,2)),0).  % -clear(?x)
unequal(2,(1,2)).

% MOVE(?x,?y,?z): 3p + 4e / Move ?y from ?z to ?x
action(3).
a_arity(3,3).
prec(3,(2,(3,3)),0) :- opt_req_prec > 0. % -clear(?z)
prec(3,(2,(3,2)),1). % below(?z,?y)
prec(3,(2,(2,2)),1). % clear(?y)
prec(3,(2,(1,2)),0) :- opt_req_prec > 0. % -below(?x,?y)
prec(3,(2,(1,1)),1). % clear(?x)
eff(3,(2,(3,3)),1).  % clear(?z)
eff(3,(2,(3,2)),0).  % -below(?z,?y)
eff(3,(2,(1,2)),1).  % below(?x,?y)
eff(3,(2,(1,1)),0).  % -clear(?x)
unequal(3,(1,2)).
unequal(3,(1,3)).
unequal(3,(2,3)).

