%
% Additional constraints
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Options:

% * opt_req_prec={1..2} (2 requires --heuristic=Domain) (effects change preconditions)
% * opt_psym={0..1}                                     (symmetry breaking on predicates)
% * opt_osym={0..3}                                     (symmetry breaking on objects)
% * opt_asym={0..3}                                     (symmetry breaking on actions)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Every effect changes some precondition

:- eff(A,M,V), not prec(A,M,1-V),           opt_req_prec=1.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on predicates

p_next(P,P+1) :- pred(P), pred(P+1), not p_static(P), not p_static(P+1), opt_psym=1.
p_next(P,P+1) :- pred(P), pred(P+1),     p_static(P),     p_static(P+1), opt_psym=1.

% Order predicates by use and arity
:- p_next(P,PP), p_arity(P,A), p_arity(PP,AA), A > AA,                   opt_psym=1.
:- p_next(P,PP), not p_used(P), p_used(PP),                              opt_psym=1.

% Order predicates/actions by predicate/arity use
a_uses_p(A,P) :-  eff(A,(P,T),V),                                        opt_psym=1.
a_uses_p(A,P) :- prec(A,(P,T),V),                                        opt_psym=1.

% Action A cannot use PP if the previous actions and A itself
% have not used a previous P (that has the same arity)
:- p_next(P,PP), p_arity(P,Ar), p_arity(PP,Ar), a_uses_p(A,PP),          opt_psym=1,
   #false: a_uses_p(AA,P), AA<=A.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on objects: (partial) order objects by appearance in denotations of root

% Option opt_osym=1
:- object(O), object(O+1), M=max_predicates, root(R), opt_osym=1,
   1 #sum{ -((1+num_objects)**((M-P)*2 + 1)),P,OO : val((P,  (O,OO)),R,1);
           -((1+num_objects)**((M-P)*2 + 1)),P,OO : val((P,  (O,OO)),  1);
           -((1+num_objects)**((M-P)*2    )),P,OO : val((P,  (OO,O)),R,1); 
           -((1+num_objects)**((M-P)*2    )),P,OO : val((P,  (OO,O)),  1); 
             (1+num_objects)**((M-P)*2 + 1 ),P,OO : val((P,(O+1,OO)),R,1);
             (1+num_objects)**((M-P)*2 + 1 ),P,OO : val((P,(O+1,OO)),  1);
             (1+num_objects)**((M-P)*2     ),P,OO : val((P,(OO,O+1)),R,1);
             (1+num_objects)**((M-P)*2     ),P,OO : val((P,(OO,O+1)),  1)}.

% Option opt_osym=2: simplest?
:- object(O), object(O+1), root(R), opt_osym=2,
   1 #sum{ -1,P,OO,1 : val((P,  (O,OO)),R,1);
           -1,P,OO,2 : val((P,  (OO,O)),R,1);
           -1,P,OO,1 : val((P,  (O,OO)),  1);
           -1,P,OO,2 : val((P,  (OO,O)),  1);
            1,P,OO,1 : val((P,(O+1,OO)),R,1);
            1,P,OO,2 : val((P,(OO,O+1)),R,1);
            1,P,OO,1 : val((P,(O+1,OO)),  1);
            1,P,OO,2 : val((P,(OO,O+1)),  1)}.

% Option opt_osym=3: separating the levels

% O occurs less or equal times than O+1 on the first argument of P
lesseq(O,O+1,P,1) :- object(O), object(O+1), pred(P), root(R), opt_osym=3,
    0 #sum{ -1,OO : val((P,  (O,OO)),  1);
            -1,OO : val((P,  (O,OO)),R,1);
             1,OO : val((P,(O+1,OO)),  1);
             1,OO : val((P,(O+1,OO)),R,1)}.

% O occurs less or equal times than O+1 on the second argument of P
lesseq(O,O+1,P,2) :- object(O), object(O+1), pred(P), root(R), opt_osym=3,
    0 #sum{ -1,OO : val((P,  (OO,O)),  1);
            -1,OO : val((P,  (OO,O)),R,1);
             1,OO : val((P,(OO,O+1)),  1);
             1,OO : val((P,(OO,O+1)),R,1)}.

% O occurs less times than O+1 on the first argument of P
  less(O,O+1,P,1) :- object(O), object(O+1), pred(P), root(R), opt_osym=3,
    1 #sum{ -1,OO : val((P,  (O,OO)),  1);
            -1,OO : val((P,  (O,OO)),R,1);
             1,OO : val((P,(O+1,OO)),  1);
             1,OO : val((P,(O+1,OO)),R,1)}.

% O occurs less times than O+1 on the second argument of P
  less(O,O+1,P,2) :- object(O), object(O+1), pred(P), root(R), opt_osym=3,
    1 #sum{ -1,OO : val((P,  (OO,O)),  1);
            -1,OO : val((P,  (OO,O)),R,1);
             1,OO : val((P,(OO,O+1)),  1);
             1,OO : val((P,(OO,O+1)),R,1)}.

% It cannot be that O occurs less than O+1 on argument I of P, and
% it occurs less or equal on the previous pairs (PP,II)
:- object(O), object(O+1), opt_osym=3,
   less(O,O+1,P,I), lesseq(O,O+1,PP,II) : pred(PP), II=1..2, (PP,II) < (P,I).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Symmetry breaking on action variables

% Directive avoids warnings when sreversible/1 is not defined
#defined sreversible/1.

% Auxiliary atoms
cond(1,A,M) :-  eff(A,M,V), opt_asym=1..3.
cond(2,A,M) :- prec(A,M,V), opt_asym=1..3.

% Option opt_asym=1
:- V=1..2, Vs=3, action(A), M=max_predicates, opt_asym=1, 
   V=1 : sreversible(A);
   1 #sum{ -((1+Vs)**((M-P)*2 + 1)),C,P,VV : cond(C,A,(P,( V,  VV)));
           -((1+Vs)**((M-P)*2    )),C,P,VV : cond(C,A,(P,( VV,  V)));
             (1+Vs)**((M-P)*2 + 1 ),C,P,VV : cond(C,A,(P,(V+1, VV)));
             (1+Vs)**((M-P)*2     ),C,P,VV : cond(C,A,(P,( VV,V+1)))}.

% Option opt_asym=2: simplest?
:- V=1..2, action(A), opt_asym=2, 
   %V=1 : sreversible(A);
   1 #sum{ -1,C,P,VV,1 : cond(C,A,(P,(  V, VV)));
           -1,C,P,VV,2 : cond(C,A,(P,( VV,  V)));
            1,C,P,VV,1 : cond(C,A,(P,(V+1, VV)));
            1,C,P,VV,2 : cond(C,A,(P,( VV,V+1)))}.

% Option opt_asym=3: separating the levels

% V occurs less or equal times than V+1 on the first argument of P (in a condition of type C of action A)
lesseq(V,V+1,A,C,P,1) :- V=1..2, action(A), C=1..2, pred(P), opt_asym=3,
    V=1 : sreversible(A);
    0 #sum{ -1,VV : cond(C,A,(P,(  V,VV)));
             1,VV : cond(C,A,(P,(V+1,VV)))}.

% V occurs less or equal times than V+1 on the second argument of P
lesseq(V,V+1,A,C,P,2) :- V=1..2, action(A), C=1..2, pred(P), opt_asym=3,
    V=1 : sreversible(A);
    0 #sum{ -1,VV : cond(C,A,(P,(VV,  V)));
             1,VV : cond(C,A,(P,(VV,V+1)))}.

% V occurs less times than V+1 on the first argument of P
less(V,V+1,A,C,P,1) :- V=1..2, action(A), C=1..2, pred(P), opt_asym=3,
    V=1 : sreversible(A);
    1 #sum{ -1,VV : cond(C,A,(P,(  V,VV)));
             1,VV : cond(C,A,(P,(V+1,VV)))}.

% V occurs less times than V+1 on the second argument of P
less(V,V+1,A,C,P,2) :- V=1..2, action(A), C=1..2, pred(P), opt_asym=3,
    V=1 : sreversible(A);
    1 #sum{ -1,VV : cond(C,A,(P,(VV,  V)));
             1,VV : cond(C,A,(P,(VV,V+1)))}.

% It cannot be that O occurs less than O+1 on argument I of P, and
% it occurs less or equal on the previous pairs (PP,II)
:- V=1..2, action(A), opt_asym=3, 
   V=1 : sreversible(A);
     less(V,V+1,A, C, P, I), 
   lesseq(V,V+1,A,CC,PP,II) : CC=1..2, pred(PP), II=1..2, (CC,PP,II) < (C,P,I).

