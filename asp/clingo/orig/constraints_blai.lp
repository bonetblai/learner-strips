% Break symmetries for predicates
:- p_arity(P, X), p_arity(P+1, Y), X > Y, P=1..max_predicates-static-1.
:- not p_used(P), p_used(P+1), P=1..max_predicates-static-1.
:- p_arity(P, X), p_arity(P+1, Y), X > Y, P=max_predicates-static+1..max_predicates-1.
:- not p_used(P), p_used(P+1), P=max_predicates-static+1..max_predicates-1.

% Non-static predicates must differ at some state
:- pred(P), not p_static(P), p_arity(P, 1), root(R), val((P,   (O, O)), S, X) : node(S), object(O), val((P, (O, O)), R, X), S != R.
:- pred(P), not p_static(P), p_arity(P, 2), root(R), val((P, (O1, O2)), S, X) : node(S), object(O1), object(O2), val((P, (O1, O2)), R, X), S != R.

