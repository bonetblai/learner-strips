%
% Additional constraints
%

%
% Options:
% * req_prec={1..2} (2 requires --heuristic=Domain)
% * psym={1}
% * osym={1..3}
% * asym={1..3}
%

% every effect changes some precondition
:- eff(A,M,V), not prec(A,M,1-V),           req_prec=1.
#heuristic  eff(A,M,1-V) : not prec(A,M,V), req_prec=2. [1,false]
#heuristic prec(A,M,1-V) :      eff(A,M,V), req_prec=2. [1, true]

% symmetry breaking on predicates
%
p_next(P,P+1) :- pred(P), pred(P+1), not p_static(P), not p_static(P+1), psym=1.
p_next(P,P+1) :- pred(P), pred(P+1),     p_static(P),     p_static(P+1), psym=1.
%    
:- p_next(P,PP), p_arity(P,A), p_arity(PP,AA), A > AA,                   psym=1.
%
:- p_next(P,PP), not p_used(P), p_used(PP),                              psym=1.
%
a_uses_p(A,P) :-  eff(A,(P,T),V),                                        psym=1.
a_uses_p(A,P) :- prec(A,(P,T),V),                                        psym=1.
%
% action A cannot use PP if 
%   the previous actions and A itself 
%   have not used a previous P (that has the same arity)
:- p_next(P,PP), p_arity(P,Ar), p_arity(PP,Ar), a_uses_p(A,PP),          psym=1,
   #false: a_uses_p(AA,P), AA<=A.

% symmetry breaking on objects
:- object(O), object(O+1), M = max_predicates, root(R), osym=1,
   1 #sum{ -((1+num_objects)**((M-P)*2 + 1)),P,OO : val((P,  (O,OO)),R,1);
           -((1+num_objects)**((M-P)*2 + 1)),P,OO : val((P,  (O,OO)),  1);
           -((1+num_objects)**((M-P)*2    )),P,OO : val((P,  (OO,O)),R,1); 
           -((1+num_objects)**((M-P)*2    )),P,OO : val((P,  (OO,O)),  1); 
             (1+num_objects)**((M-P)*2 + 1 ),P,OO : val((P,(O+1,OO)),R,1);
             (1+num_objects)**((M-P)*2 + 1 ),P,OO : val((P,(O+1,OO)),  1);
             (1+num_objects)**((M-P)*2     ),P,OO : val((P,(OO,O+1)),R,1);
             (1+num_objects)**((M-P)*2     ),P,OO : val((P,(OO,O+1)),  1)}.

% simpler symmetry breaking on objects
:- object(O), object(O+1), root(R), osym=2,
   1 #sum{ -1,P,OO,1 : val((P,  (O,OO)),R,1);
           -1,P,OO,2 : val((P,  (OO,O)),R,1);
           -1,P,OO,1 : val((P,  (O,OO)),  1);
           -1,P,OO,2 : val((P,  (OO,O)),  1);
            1,P,OO,1 : val((P,(O+1,OO)),R,1);
            1,P,OO,2 : val((P,(OO,O+1)),R,1);
            1,P,OO,1 : val((P,(O+1,OO)),  1);
            1,P,OO,2 : val((P,(OO,O+1)),  1)}.

% symmetry breaking on objects separating the levels
%
% O occurs less or equal times than O+1 on the first argument of P
lesseq(O,O+1,P,1) :- object(O), object(O+1), pred(P), root(R), osym=3,
    0 #sum{ -1,OO : val((P,  (O,OO)),  1);
            -1,OO : val((P,  (O,OO)),R,1);
             1,OO : val((P,(O+1,OO)),  1);
             1,OO : val((P,(O+1,OO)),R,1)}.
%
% O occurs less or equal times than O+1 on the second argument of P
lesseq(O,O+1,P,2) :- object(O), object(O+1), pred(P), root(R), osym=3,
    0 #sum{ -1,OO : val((P,  (OO,O)),  1);
            -1,OO : val((P,  (OO,O)),R,1);
             1,OO : val((P,(OO,O+1)),  1);
             1,OO : val((P,(OO,O+1)),R,1)}.
%
% O occurs less times than O+1 on the first argument of P
  less(O,O+1,P,1) :- object(O), object(O+1), pred(P), root(R), osym=3,
    1 #sum{ -1,OO : val((P,  (O,OO)),  1);
            -1,OO : val((P,  (O,OO)),R,1);
             1,OO : val((P,(O+1,OO)),  1);
             1,OO : val((P,(O+1,OO)),R,1)}.
%
% O occurs less times than O+1 on the second argument of P
  less(O,O+1,P,2) :- object(O), object(O+1), pred(P), root(R), osym=3,
    1 #sum{ -1,OO : val((P,  (OO,O)),  1);
            -1,OO : val((P,  (OO,O)),R,1);
             1,OO : val((P,(OO,O+1)),  1);
             1,OO : val((P,(OO,O+1)),R,1)}.
%
% it cannot be that O occurs less than O+1 on argument I of P, 
% and it occurs less or equal on the previous pairs (PP,II)
:- object(O), object(O+1), osym=3,
   less(O,O+1,P,I), lesseq(O,O+1,PP,II) : pred(PP), II=1..2, (PP,II) < (P,I).


% symmetry breaking on action variables
:- V=1..2, Vs=3, action(A), M=max_predicates, asym=1, 
   V=1 : sreversible(A);
   1 #sum{ -((1+Vs)**((M-P)*2 + 1)),C,P,VV : cond(C,A,(P,( V,  VV)));
           -((1+Vs)**((M-P)*2    )),C,P,VV : cond(C,A,(P,( VV,  V)));
             (1+Vs)**((M-P)*2 + 1 ),C,P,VV : cond(C,A,(P,(V+1, VV)));
             (1+Vs)**((M-P)*2     ),C,P,VV : cond(C,A,(P,( VV,V+1)))}.

% simpler symmetry breaking on action variables
:- V=1..2, action(A), asym=2, 
   V=1 : sreversible(A);
   1 #sum{ -1,C,P,VV,1 : cond(C,A,(P,(  V, VV)));
           -1,C,P,VV,2 : cond(C,A,(P,( VV,  V)));
            1,C,P,VV,1 : cond(C,A,(P,(V+1, VV)));
            1,C,P,VV,2 : cond(C,A,(P,( VV,V+1)))}.

% symmetry breaking on action variables separating the levels
%
% V occurs less or equal times than V+1 on the first argument of P
% (in a condition of type C of action A)
lesseq(V,V+1,A,C,P,1) :- V=1..2, action(A), C=1..2, pred(P), asym=3,
    V=1 : sreversible(A);
    0 #sum{ -1,VV : cond(C,A,(P,(  V,VV)));
             1,VV : cond(C,A,(P,(V+1,VV)))}.
%
% V occurs less or equal times than V+1 on the second argument of P
lesseq(V,V+1,A,C,P,2) :- V=1..2, action(A), C=1..2, pred(P), asym=3,
    V=1 : sreversible(A);
    0 #sum{ -1,VV : cond(C,A,(P,(VV,  V)));
             1,VV : cond(C,A,(P,(VV,V+1)))}.
%
% V occurs less times than V+1 on the first argument of P
less(V,V+1,A,C,P,1) :- V=1..2, action(A), C=1..2, pred(P), asym=3,
    V=1 : sreversible(A);
    1 #sum{ -1,VV : cond(C,A,(P,(  V,VV)));
             1,VV : cond(C,A,(P,(V+1,VV)))}.
%
% V occurs less times than V+1 on the second argument of P
less(V,V+1,A,C,P,2) :- V=1..2, action(A), C=1..2, pred(P), asym=3,
    V=1 : sreversible(A);
    1 #sum{ -1,VV : cond(C,A,(P,(VV,  V)));
             1,VV : cond(C,A,(P,(VV,V+1)))}.
%
% it cannot be that O occurs less than O+1 on argument I of P, 
% and it occurs less or equal on the previous pairs (PP,II)
:- V=1..2, action(A), asym=3, 
   V=1 : sreversible(A);
     less(V,V+1,A, C, P, I), 
   lesseq(V,V+1,A,CC,PP,II) : CC=1..2, pred(PP), II=1..2, (CC,PP,II) < (C,P,I).

% this clingo directive avoids warnings when sreversible/1 is not defined
#defined sreversible/1.

% auxiliary atoms
cond(1,A,M) :-  eff(A,M,V), asym=1..3.
cond(2,A,M) :- prec(A,M,V), asym=1..3.
