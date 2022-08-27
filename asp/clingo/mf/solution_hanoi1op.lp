% Predicates
pred(1).      % clear(?d)
p_arity(1,1).
pred(2).      % on(?d1,?d2)
p_arity(2,2).
pred(3).      % BIGGER(?d1,?d2)
p_arity(3,2).
p_static(3).

% MOVE(?d,?to,?fr): 4p + 4e
action(1).
a_arity(1,3).
%%prec(1,(3,(3,1)),1). % BIGGER(?to,?d)
prec(1,(3,(2,1)),1). % BIGGER(?fr,?d)
prec(1,(1,(1,1)),1). % clear(?d)
%%prec(1,(1,(3,3)),0). % -clear(?to)
prec(1,(1,(2,2)),1). % clear(?fr)
prec(1,(2,(1,3)),1). % on(?d,?to)
%%prec(1,(2,(3,1)),0). % -on(?to,?d)
%%prec(1,(2,(1,2)),0). % -on(?d,?fr)
eff(1,(1,(3,3)),1).  % clear(?to)
eff(1,(1,(2,2)),0).  % -clear(?fr)
eff(1,(2,(1,2)),1).  % on(?d,?fr)
eff(1,(2,(1,3)),0).  % -on(?d,?to)
unequal(1,(1,2)).
unequal(1,(1,3)).
unequal(1,(2,3)).

